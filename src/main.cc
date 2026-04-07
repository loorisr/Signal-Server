/****************************************************************************\
*  Signal Server: Radio propagation simulator by Alex Farrant QCVS, 2E0TDW   *
******************************************************************************
*    SPLAT! Project started in 1997 by John A. Magliacane, KD2BD             *
******************************************************************************
*         Please consult the SPLAT! documentation for a complete list of     *
*         individuals who have contributed to this project.                  *
******************************************************************************
*                                                                            *
*  This program is free software; you can redistribute it and/or modify it   *
*  under the terms of the GNU General Public License as published by the     *
*  Free Software Foundation; either version 2 of the License or any later    *
*  version.                                                                  *
*                                                                            *
*  This program is distributed in the hope that it will useful, but WITHOUT  *
*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or     *
*  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License     *
*  for more details.                                                         *
*                                                                            *
\****************************************************************************/

/* GDAL must come before any header that defines MAX */
#include <gdal.h>
#include <cpl_conv.h>
#include <ogr_srs_api.h>

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <errno.h>
#include <limits.h>

#include "main.hh"
#include "common.hh"
#include "inputs.hh"
#include "outputs.hh"
#include "models/itwom3.0.hh"
#include "models/los.hh"
#include "models/pel.hh"
#include "image.hh"

#include <chrono>
#include <thread>

#include <spdlog/spdlog.h>

int ippd = 1200;
int MAX_DISTANCE_DEGRES = 3; // max distance : 3° so around 300 km
int ARRAYSIZE = (MAX_DISTANCE_DEGRES * ippd) + 10;

char DEM_path[255], *color_file = NULL;
std::string output_filename;

double max_range = 0.0,  dpp, ppd, samples_per_radian,
    fzone_clearance = 0.6, clutter, lat, lon, txh, tercon, terdic,
    north, east, south, west, dBm, loss, field_strength,
    min_north = 90, max_north = -90, min_lon = 180.0, max_lon = -180.0,
    min_lat = 90.0, max_lat = -90.0,
    delta=0, rxGain=0, antenna_rotation,
    antenna_downtilt,antenna_dt_direction, cropLat=-70, cropLon=0,cropLonNeg=0;

int mpi, max_elevation = -32768, min_elevation = 32768,
    contour_threshold, debug = 0,
    height, width;

std::atomic<int> cnt_point_to_point_ITM{0};
std::atomic<int> cnt_point_to_point{0};
std::atomic<int> cnt_computeLoss{0};
std::atomic<int> cnt_PlotPropPath{0};
std::atomic<int> cnt_ReadPath{0};
std::atomic<int> cnt_PlotPropagation{0};

bool got_elevation_pattern = false, got_azimuth_pattern = false, dbm = false;
bool geotiff = false;
bool ngs = false;
bool to_stdout = false, cropping = true;
int knifeedge = 0;
int pmenv = 1;

__thread double *elev;
__thread struct path path;
site tx_site;
site rx_site;
double         **dem_data   = nullptr;
unsigned char **dem_signal = nullptr;
int dem_min_lat   = 0;
int dem_min_lon   = 0;
int dem_width_px  = 0;
int dem_height_px = 0;

struct LR LR;
struct region region;

/* Convert (lat, lon) to flat-array pixel coordinates.
 * x increases northward from the south edge of the loaded area.
 * y increases westward  from the east  edge of the loaded area.
 * Returns true and sets x_out/y_out when the point is inside the array. */
static bool find_dem_xy(double lat, double lon, int &x_out, int &y_out)
{
    if (!dem_data) return false;
    int x = (int)rint(ppd  * (lat - dem_min_lat));
    int y = (int)rint(ppd * (lon - dem_min_lon));
    if (x < 0 || x >= dem_height_px || y < 0 || y >= dem_width_px) return false;
    x_out = x;
    y_out = y;
    return true;
}

void PutSignal(double lat, double lon, unsigned char signal)
{
    int x, y;
    if (find_dem_xy(lat, lon, x, y))
        //dem_signal[x][y] = MAX(signal, GetSignal(lat, lon));
        dem_signal[x][y] = signal;
}

unsigned char GetSignal(double lat, double lon)
{
    int x, y;
    if (!find_dem_xy(lat, lon, x, y)) return 0;
    return dem_signal[x][y];
}

double GetElevation(struct site location)
{
    int x, y;
    if (!find_dem_xy(location.lat, location.lon, x, y)) return -5000.0;
    return (double)dem_data[x][y];
}


double ElevationAngle(struct site source, struct site destination)
{
    /* This function returns the angle of elevation (in degrees)
       of the destination as seen from the source location.
       A positive result represents an angle of elevation (uptilt),
       while a negative result represents an angle of depression
       (downtilt), as referenced to a normal to the center of
       the earth. */

    double a, b, dx;

    a = GetElevation(destination) + destination.alt + EARTHRADIUS;
    b = GetElevation(source) + source.alt + EARTHRADIUS;

    dx = Distance(source, destination);

    /* Apply the Law of Cosines */

    return ((180.0 *
         (acos(((b * b) + (dx * dx) - (a * a)) / (2.0 * b * dx))) /
         PI) - 90.0);
}

/* This function generates a sequence of latitude and
       longitude positions between source and destination
       locations along a great circle path, and stores
       elevation and distance information for points
       along that path in the "path" structure. 
*/
void ReadPath(struct site source, struct site destination)
{
    if (debug) cnt_ReadPath++;
    int c;
    double azimuth, distance, lat1, lon1, beta, den, num,
        lat2, lon2, total_distance, dx, dy, path_length,
        m_per_sample;
    struct site tempsite;

    lat1 = source.lat * DEG2RAD;
    lon1 = source.lon * DEG2RAD;
    lat2 = destination.lat * DEG2RAD;
    lon2 = destination.lon * DEG2RAD;

    azimuth = Azimuth(source, destination) * DEG2RAD;
    total_distance = Distance(source, destination);  

    dx = samples_per_radian * (lon1 - lon2);
    dy = samples_per_radian * (lat1 - lat2);
    path_length = sqrt((dx * dx) + (dy * dy));
    m_per_sample = total_distance / path_length;

    for (distance = 0.0, c = 0;
         (total_distance != 0.0 && distance <= total_distance
          && c < ARRAYSIZE); c++, distance = m_per_sample * (double)c) {

        beta = distance / EARTHRADIUS; 
        lat2 = asin(sin(lat1) * cos(beta) + cos(azimuth) * sin(beta) * cos(lat1));
        num = cos(beta) - (sin(lat1) * sin(lat2));
        den = cos(lat1) * cos(lat2);

        if (azimuth == 0.0 && (beta > HALFPI - lat1))
            lon2 = lon1 + PI;

        else if (azimuth == HALFPI && (beta > HALFPI + lat1))
            lon2 = lon1 + PI;

        else if (fabs(num / den) > 1.0)
            lon2 = lon1;

        else {
            if ((PI - azimuth) >= 0.0)
                lon2 = lon1 - arccos(num, den);
            else
                lon2 = lon1 + arccos(num, den);
        }

        if (lon2 < -PI)
            lon2 += TWOPI;
        else if (lon2 > PI)
            lon2 -= TWOPI;

        lat2 = lat2 * RAD2DEG;
        lon2 = lon2 * RAD2DEG;

        path.lat[c] = lat2;
        path.lon[c] = lon2;
        tempsite.lat = lat2;
        tempsite.lon = lon2;
        path.elevation[c] = GetElevation(tempsite);
        path.distance[c] = distance;
    }

    /* Make sure exact destination point is recorded at path.length-1 */
    /* Check if really useful */

    path.lat[c] = destination.lat;
    path.lon[c] = destination.lon;
    path.elevation[c] = GetElevation(destination);
    path.distance[c] = total_distance;
    c++;

    path.length = c;
}

double ElevationAngle2(struct site source, struct site destination, double er)
{
    /* This function returns the angle of elevation (in degrees)
       of the destination as seen from the source location, UNLESS
       the path between the sites is obstructed, in which case, the
       elevation angle to the first obstruction is returned instead.
       "er" represents the earth radius. */

    int x;
    char block = 0;
    double source_alt, destination_alt, cos_xmtr_angle,
        cos_test_angle, test_alt, elevation, distance,
        source_alt2, first_obstruction_angle = 0.0;
    struct path temp;

    temp = path;

    ReadPath(source, destination);

    distance = Distance(source, destination); 
    source_alt = er + source.alt + GetElevation(source);
    destination_alt = er + destination.alt + GetElevation(destination);
    source_alt2 = source_alt * source_alt;

    /* Calculate the cosine of the elevation angle of the
       destination (receiver) as seen by the source (transmitter). */

    cos_xmtr_angle =
        ((source_alt2) + (distance * distance) -
         (destination_alt * destination_alt)) / (2.0 * source_alt *
                             distance);

    /* Test all points in between source and destination locations to
       see if the angle to a topographic feature generates a higher
       elevation angle than that produced by the destination.  Begin
       at the source since we're interested in identifying the FIRST
       obstruction along the path between source and destination. */

    for (x = 2, block = 0; x < path.length && block == 0; x++) {
        distance = path.distance[x];

        test_alt =
            EARTHRADIUS + (path.elevation[x] ==
                   0.0 ? path.elevation[x] : path.elevation[x] +
                   clutter);

        cos_test_angle =
            ((source_alt2) + (distance * distance) -
             (test_alt * test_alt)) / (2.0 * source_alt * distance);

        /* Compare these two angles to determine if
           an obstruction exists.  Since we're comparing
           the cosines of these angles rather than
           the angles themselves, the sense of the
           following "if" statement is reversed from
           what it would be if the angles themselves
           were compared. */

        if (cos_xmtr_angle >= cos_test_angle) {
            block = 1;
            first_obstruction_angle =
                ((acos(cos_test_angle)) * RAD2DEG) - 90.0;
        }
    }

    if (block)
        elevation = first_obstruction_angle;

    else
        elevation = ((acos(cos_xmtr_angle)) * RAD2DEG) - 90.0;

    path = temp;

    return elevation;
}

void ObstructionAnalysis(struct site xmtr, struct site rcvr, double f,
             FILE *outfile)
{
    /* Perform an obstruction analysis along the
       path between receiver and transmitter. */

    int x;
    struct site site_x;
    double h_r, h_t, h_x, h_r_orig, cos_tx_angle, cos_test_angle,
        cos_tx_angle_f1, cos_tx_angle_fpt6, d_tx, d_x,
        h_r_f1, h_r_fpt6, h_f, h_los, lambda = 0.0;
    char string[255], string_fpt6[255], string_f1[255];

    ReadPath(xmtr, rcvr);
    h_r = GetElevation(rcvr) + rcvr.alt + EARTHRADIUS;
    h_r_f1 = h_r;
    h_r_fpt6 = h_r;
    h_r_orig = h_r;
    h_t = GetElevation(xmtr) + xmtr.alt + EARTHRADIUS;
    d_tx = Distance(rcvr, xmtr); 
    cos_tx_angle =
        ((h_r * h_r) + (d_tx * d_tx) - (h_t * h_t)) / (2.0 * h_r * d_tx);
    cos_tx_angle_f1 = cos_tx_angle;
    cos_tx_angle_fpt6 = cos_tx_angle;

    if (f)
        lambda = 299792458.0 / (f * 1e6);

    if (clutter > 0.0) {
        fprintf(outfile, "Terrain has been raised by");
        fprintf(outfile, " %.2f meters", clutter);
        fprintf(outfile, " to account for ground clutter.\n\n");
    }

    /* At each point along the path calculate the cosine
       of a sort of "inverse elevation angle" at the receiver.
       From the antenna, 0 deg. looks at the ground, and 90 deg.
       is parallel to the ground.

       Start at the receiver.  If this is the lowest antenna,
       then terrain obstructions will be nearest to it.  (Plus,
       that's the way ppa!'s original los() did it.)

       Calculate cosines only.  That's sufficient to compare
       angles and it saves the extra computational burden of
       acos().  However, note the inverted comparison: if
       acos(A) > acos(B), then B > A. */

    for (x = path.length - 1; x > 0; x--) {
        site_x.lat = path.lat[x];
        site_x.lon = path.lon[x];
        site_x.alt = 0.0;

        h_x = GetElevation(site_x) + EARTHRADIUS + clutter;
        d_x = Distance(rcvr, site_x);

        /* Deal with the LOS path first. */

        cos_test_angle =
            ((h_r * h_r) + (d_x * d_x) -
             (h_x * h_x)) / (2.0 * h_r * d_x);

        if (cos_tx_angle > cos_test_angle) {
            if (h_r == h_r_orig)
                fprintf(outfile,
                    "Between RX and TX, obstructions were detected :\n\n");

            if (site_x.lat >= 0.0) {
                fprintf(outfile,
                    "   %8.4f N,%9.4f W, %5.2f kilometers, %6.2f meters AMSL\n",
                    site_x.lat, site_x.lon,
                    d_x / 1000.0,
                    h_x - EARTHRADIUS);
            }

            else {
                fprintf(outfile,
                    "   %8.4f S,%9.4f W, %5.2f kilometers, %6.2f meters AMSL\n",
                    -site_x.lat, site_x.lon,
                    d_x / 1000.0,
                    h_x - EARTHRADIUS);
            }
        }

        while (cos_tx_angle > cos_test_angle) {
            h_r += 1;
            cos_test_angle =
                ((h_r * h_r) + (d_x * d_x) -
                 (h_x * h_x)) / (2.0 * h_r * d_x);
            cos_tx_angle =
                ((h_r * h_r) + (d_tx * d_tx) -
                 (h_t * h_t)) / (2.0 * h_r * d_tx);
        }

        if (f) {
            /* Now clear the first Fresnel zone... */

            cos_tx_angle_f1 =
                ((h_r_f1 * h_r_f1) + (d_tx * d_tx) -
                 (h_t * h_t)) / (2.0 * h_r_f1 * d_tx);
            h_los =
                sqrt(h_r_f1 * h_r_f1 + d_x * d_x -
                 2 * h_r_f1 * d_x * cos_tx_angle_f1);
            h_f = h_los - sqrt(lambda * d_x * (d_tx - d_x) / d_tx);

            while (h_f < h_x) {
                h_r_f1 += 1;
                cos_tx_angle_f1 =
                    ((h_r_f1 * h_r_f1) + (d_tx * d_tx) -
                     (h_t * h_t)) / (2.0 * h_r_f1 * d_tx);
                h_los =
                    sqrt(h_r_f1 * h_r_f1 + d_x * d_x -
                     2 * h_r_f1 * d_x * cos_tx_angle_f1);
                h_f =
                    h_los -
                    sqrt(lambda * d_x * (d_tx - d_x) / d_tx);
            }

            /* and clear the 60% F1 zone. */

            cos_tx_angle_fpt6 =
                ((h_r_fpt6 * h_r_fpt6) + (d_tx * d_tx) -
                 (h_t * h_t)) / (2.0 * h_r_fpt6 * d_tx);
            h_los =
                sqrt(h_r_fpt6 * h_r_fpt6 + d_x * d_x -
                 2 * h_r_fpt6 * d_x * cos_tx_angle_fpt6);
            h_f =
                h_los -
                fzone_clearance * sqrt(lambda * d_x * (d_tx - d_x) /
                           d_tx);

            while (h_f < h_x) {
                h_r_fpt6 += 1;
                cos_tx_angle_fpt6 =
                    ((h_r_fpt6 * h_r_fpt6) + (d_tx * d_tx) -
                     (h_t * h_t)) / (2.0 * h_r_fpt6 * d_tx);
                h_los =
                    sqrt(h_r_fpt6 * h_r_fpt6 + d_x * d_x -
                     2 * h_r_fpt6 * d_x *
                     cos_tx_angle_fpt6);
                h_f =
                    h_los -
                    fzone_clearance * sqrt(lambda * d_x *
                               (d_tx - d_x) / d_tx);
            }
        }
    }

    if (h_r > h_r_orig) {
        snprintf(string, 150,
             "\nAntenna at Rx must be raised to at least %.2f meters AGL\nto clear all obstructions detected.\n",
             h_r - GetElevation(rcvr) - EARTHRADIUS);
    }

    else
        snprintf(string, 150,
             "\nNo obstructions to LOS path due to terrain were detected\n");

    if (f) {
        if (h_r_fpt6 > h_r_orig) {
            snprintf(string_fpt6, 150,
                 "\nAntenna at Rx must be raised to at least %.2f meters AGL\nto clear %.0f%c of the first Fresnel zone.\n",
                 h_r_fpt6 - GetElevation(rcvr) - EARTHRADIUS,
                 fzone_clearance * 100.0, 37);
        }

        else
            snprintf(string_fpt6, 150,
                 "\n%.0f%c of the first Fresnel zone is clear.\n",
                 fzone_clearance * 100.0, 37);

        if (h_r_f1 > h_r_orig) {
            snprintf(string_f1, 150,
                 "\nAntenna at Rx must be raised to at least %.2f meters AGL\nto clear the first Fresnel zone.\n",
                 h_r_f1 - GetElevation(rcvr) - EARTHRADIUS);
        }

        else
            snprintf(string_f1, 150,
                 "\nThe first Fresnel zone is clear.\n");
    }

    fprintf(outfile, "%s", string);

    if (f) {
        fprintf(outfile, "%s", string_f1);
        fprintf(outfile, "%s", string_fpt6);
    }

}

void free_dem(void)
{
    for (int i = 0; i < dem_height_px; i++) {
        delete [] dem_data[i];
        delete [] dem_signal[i];
    }
    delete [] dem_data;   dem_data   = nullptr;
    delete [] dem_signal; dem_signal = nullptr;
    dem_width_px = dem_height_px = 0;
}

void free_elev(void) {
  delete [] elev;
}

void free_path(void)
{
    delete [] path.lat;
    delete [] path.lon;
    delete [] path.elevation;
    delete [] path.distance;
}

void alloc_elev(void)
{
  elev  = new double[ARRAYSIZE + 10];
}

/* Allocate the flat DEM arrays for a known bounding box.
 * Called from LoadTopoData() after tiles_lat/tiles_lon are computed. */
void alloc_dem(int min_lat, int min_lon, int tiles_lat, int tiles_lon)
{
    if (dem_data) free_dem();

    dem_min_lat   = min_lat;
    dem_min_lon   = min_lon;
    dem_height_px = tiles_lat * ippd;
    dem_width_px  = tiles_lon * ippd;

    dem_data   = new double         *[dem_height_px];
    dem_signal = new unsigned char *[dem_height_px];
    for (int i = 0; i < dem_height_px; i++) {
        dem_data[i]   = new double        [dem_width_px]();
        dem_signal[i] = new unsigned char[dem_width_px]();
    }
}

void alloc_path(void)
{
    path.lat = new double[ARRAYSIZE + 1];
    path.lon = new double[ARRAYSIZE + 1];
    path.elevation = new double[ARRAYSIZE + 1];
    path.distance = new double[ARRAYSIZE + 1];
}

void do_allocs(void)
{
    alloc_elev();
    alloc_path();
    /* dem is allocated later by alloc_dem() once the bounding box is known */
}

void write_geotiff_from_canvas(const uint8_t *canvas, int img_width, int img_height, const char *filename)
{
    char tif_file[300];
    double ulx, uly, lrx, lry;

    /* Build .tif output path */
    size_t len = strlen(filename);
    if (len > 4 && strcmp(filename + len - 4, ".ppm") == 0)
        snprintf(tif_file, sizeof(tif_file), "%.*s.tif", (int)(len - 4), filename);
    else
        snprintf(tif_file, sizeof(tif_file), "%s.tif", filename);

    /* Compute geographic bounds */
    if (cropping) {
        ulx = tx_site.lon - cropLon;
        uly = tx_site.lat + cropLat;
        lrx = tx_site.lon + cropLon;
        lry = tx_site.lat - cropLat;
    } else {
        ulx = west;
        uly = max_north;
        lrx = east;
        lry = min_north;
    }
    
    GDALDriverH drv = GDALGetDriverByName("GTiff");
    if (drv == NULL) {
        spdlog::error("write_geotiff_from_canvas: GTiff GDAL driver not available");
        return;
    }

    char *create_opts[] = {
        (char *)"COMPRESS=LZW",
        (char *)"PREDICTOR=2",
        (char *)"TILED=YES",
        (char *)"BLOCKXSIZE=256",
        (char *)"BLOCKYSIZE=256",
        NULL
    };
    GDALDatasetH ds = GDALCreate(drv, tif_file, img_width, img_height, 4, GDT_Byte, create_opts);
    if (ds == NULL) {
        spdlog::error("write_geotiff_from_canvas: failed to create {}", tif_file);
        return;
    }

    double gt[6] = {
        ulx,
        (lrx - ulx) / img_width,
        0.0,
        uly,
        0.0,
        (lry - uly) / img_height   /* negative → north-up raster */
    };
    GDALSetGeoTransform(ds, gt);

    OGRSpatialReferenceH srs = OSRNewSpatialReference(NULL);
    OSRImportFromEPSG(srs, 4326);
    char *wkt = NULL;
    OSRExportToWkt(srs, &wkt);
    GDALSetProjection(ds, wkt);
    CPLFree(wkt);
    OSRDestroySpatialReference(srs);

    GDALSetRasterColorInterpretation(GDALGetRasterBand(ds, 4), GCI_AlphaBand);

    /* Build RGBA buffer: white pixels become transparent, others opaque */
    int npixels = img_width * img_height;
    std::vector<uint8_t> rgba(npixels * 4);
    for (int i = 0; i < npixels; ++i) {
        uint8_t r = canvas[i * RGB_SIZE + 0];
        uint8_t g = canvas[i * RGB_SIZE + 1];
        uint8_t b = canvas[i * RGB_SIZE + 2];
        rgba[i * 4 + 0] = r;
        rgba[i * 4 + 1] = g;
        rgba[i * 4 + 2] = b;
        rgba[i * 4 + 3] = (r == 255 && g == 255 && b == 255) ? 0 : 255;
    }

    /* Write interleaved RGBA buffer to 4 separate GDAL bands */
    int bandMap[4] = {1, 2, 3, 4};
    CPLErr err = (CPLErr)GDALDatasetRasterIO(ds, GF_Write,
        0, 0, img_width, img_height,
        rgba.data(), img_width, img_height, GDT_Byte,
        4, bandMap,
        4,              /* nPixelSpace */
        img_width * 4,  /* nLineSpace  */
        1);             /* nBandSpace  */
    if (err != CE_None)
        spdlog::error("write_geotiff_from_canvas: RasterIO write failed for {}", tif_file);
    else
        spdlog::info("GeoTIFF written: {}", tif_file);

    GDALClose(ds);
}

int main(int argc, char *argv[])
{
    auto start_time = std::chrono::steady_clock::now();

    int x, y, z = 0, ppa = 0, normalise = 0,
      result,
      number_threads = std::max(4u, (std::thread::hardware_concurrency() / 2) * 2);

    PropModel prop_model;

    double rxlat, rxlon;

    char mapfile[255], antenna_file[255];
    char *az_filename, *el_filename = NULL;

    double altitudeLR = 0.0;

    spdlog::info("Version {}.{} ({} {})", VER_MAJ, VER_MIN, GIT_BRANCH, GIT_COMMIT_HASH);
    spdlog::info("    Compile date: {} {}", __DATE__, __TIME__);
    spdlog::info("");

    if (argc == 1) {
        spdlog::info("License: GNU General Public License (GPL) version 2");
        spdlog::info("");
        spdlog::info("Radio propagation simulator by Alex Farrant QCVS, 2E0TDW");
        spdlog::info("Based upon SPLAT! by John Magliacane, KD2BD");
        spdlog::info("Some feature enhancements/additions by Aaron A. Collins, N9OZB");
        spdlog::info("Additional improvements and multithreading fixes by P. McDonnell, W3AXL");
        spdlog::info("");
        spdlog::info("Usage: signalserver [data options] [input options] [antenna options] [output options] -o outputfile");
        spdlog::info("");
        spdlog::info("Data:");
        spdlog::info("     -dem Directory containing Copernicus DEM GeoTIFF COG tiles");
        spdlog::info("                  (Copernicus_DSM_COG_30_N##_00_?###_00_DEM.tif for 1200 ppd,");
        spdlog::info("                   Copernicus_DSM_COG_10_N##_00_?###_00_DEM.tif for 3600 ppd)");
        spdlog::info("     -color File to pre-load .scf/.lcf/.dcf for Signal/Loss/dBm color palette");
        spdlog::info("Input:");
        spdlog::info("     -lat Tx Latitude (decimal degrees) -70/+70");
        spdlog::info("     -lon Tx Longitude (decimal degrees) -180/+180");
        spdlog::info("     -rla (Optional) Rx Latitude for PPA (decimal degrees) -70/+70");
        spdlog::info("     -rlo (Optional) Rx Longitude for PPA (decimal degrees) -180/+180");
        spdlog::info("     -f Tx Frequency (MHz) 20MHz to 100GHz (LOS after 20GHz)");
        spdlog::info("     -erp Tx Total Effective Radiated Power in Watts (dBd) inc Tx+Rx gain. 2.14dBi = 0dBd");
        spdlog::info("     -gc Random ground clutter (meters)");
        spdlog::info("     -te Terrain code 1-6 (optional - 1. Water, 2. Marsh, 3. Farmland,");
        spdlog::info("          4. Mountain, 5. Desert, 6. Urban");
        spdlog::info("     -terdic Terrain dielectric value 2-80 (optional)");
        spdlog::info("     -tercon Terrain conductivity 0.01-0.0001 (optional)");
        spdlog::info("     -cl Climate code 1-7 (optional - 1. Equatorial 2. Continental subtropical");
        spdlog::info("          3. Maritime subtropical 4. Desert 5. Continental temperate");
        spdlog::info("          6. Maritime temperate (Land) 7. Maritime temperate (Sea)");
        spdlog::info("     -rel Reliability for ITM model (% of 'time') 1 to 99 (optional, default 50%)");
        spdlog::info("     -conf Confidence for ITM model (% of 'situations') 1 to 99 (optional, default 50%)");
        spdlog::info("     -number_threads Number of number_threads to divide the plot rectangle into (must be even and > 4)");
        spdlog::info("     -hd Use HD mode (30m), per defaut 90m");
        spdlog::info("Output:");
        spdlog::info("     -o basename (Output file basename - required, min 5 chars)");
        spdlog::info("     -dbm Plot Rxd signal power instead of field strength in dBuV/m");
        spdlog::info("     -rt Rx Threshold (dB / dBm / dBuV/m)");
        spdlog::info("     -R Radius (kilometers)");
        spdlog::info("     -pm Propagation model. 1: ITM, 2: LOS, 3: Hata, 4: ECC33,");
        spdlog::info("          5: SUI, 6: COST-Hata, 7: FSPL, 8: ITWOM, 9: Ericsson,");
        spdlog::info("          10: Plane earth, 11: Egli VHF/UHF, 12: Soil");
        spdlog::info("     -pe Propagation model mode: 1=Urban,2=Suburban,3=Rural");
        spdlog::info("     -ked Knife edge diffraction (Already on for ITM)");
        spdlog::info("     -geotiff Output a geotiff file");
        spdlog::info("Antenna:");
        spdlog::info("     -ant (antenna pattern file basename+path for .az and .el files)");
        spdlog::info("     -txh Tx Height (above ground)");
        spdlog::info("     -rxh Rx Height(s) (optional. Default=0.1)");
        spdlog::info("     -rxg Rx gain dBd (optional for PPA text report)");
        spdlog::info("     -hp Horizontal Polarisation (default=vertical)");
        spdlog::info("     -rot  (  0.0 - 359.0 degrees, default 0.0) Antenna Pattern Rotation");
        spdlog::info("     -dt   ( -10.0 - 90.0 degrees, default 0.0) Antenna Downtilt");
        spdlog::info("     -dtdir ( 0.0 - 359.0 degrees, default 0.0) Antenna Downtilt Direction");
        spdlog::info("Debugging:");
        spdlog::info("     -t Terrain greyscale background");
        spdlog::info("     -dbg Verbose debug messages");
        spdlog::info("     -ng Normalise Path Profile graph");

        return 1;
    }

    GDALAllRegister();

    do_allocs();

    y = argc - 1;
    dbm = false;
    DEM_path[0] = 0;
    mapfile[0] = 0;
    clutter = 0.0;
    color_file = NULL;
    path.length = 0;
    fzone_clearance = 0.6;
    contour_threshold = 0;
    max_range = 1.0;
    prop_model = ITM_LR;
    lat = 0;
    lon = 0;
    txh = 0;
    ngs = true;			// no greyscale background

    sscanf("0.1", "%lf", &altitudeLR);

    // Defaults
    LR.eps_dielect = 15.0;	// Farmland
    LR.sgm_conductivity = 0.005;	// Farmland
    LR.eno_ns_surfref = 301.0;
    LR.frq_mhz = 19.0;	// Deliberately too low
    LR.radio_climate = 5;	// continental
    LR.pol = 1;		// vert
    LR.conf = 0.50;
    LR.rel = 0.50;
    LR.erp = 0.0;		// will default to Path Loss

    antenna_rotation = -1;  // unique defaults to test usage
    antenna_downtilt = 99.0; // don't mess with them!
    antenna_dt_direction = -1;
    antenna_file[0] = '\0';

    tx_site.lat = 91.0;
    tx_site.lon = 181.0;
    rx_site.lat = 91.0;
    rx_site.lon = 181.0;

    /* Scan for command line arguments */

    for (x = 1; x <= y; x++) {

        if (strcmp(argv[x], "-R") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0] && argv[z][0] != '-') {
                sscanf(argv[z], "%lf", &max_range);

            }
        }

        if (strcmp(argv[x], "-gc") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0] && argv[z][0] != '-') {
                sscanf(argv[z], "%lf", &clutter);

                if (clutter < 0.0)
                    clutter = 0.0;

                
            }
        }

        if (strcmp(argv[x], "-ant") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0] && argv[z][0] != '-') {
                strncpy(antenna_file, argv[z], 253);
            }
        }

        if (strcmp(argv[x], "-rot") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0] && argv[z][0] != '-') {
                sscanf(argv[z], "%lf", &antenna_rotation);

                if (antenna_rotation < 0.0)
                    antenna_rotation = 0.0;
                if (antenna_rotation > 359.0)
                    antenna_rotation = 0.0;
            }
        }

        if (strcmp(argv[x], "-dt") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0]) {	/* A minus argument is legal here */
                sscanf(argv[z], "%lf", &antenna_downtilt);
                if (antenna_downtilt < -10.0)
                    antenna_downtilt = -10.0;
                if (antenna_downtilt > 90.0)
                    antenna_downtilt = 90.0;
            }
        }

        if (strcmp(argv[x], "-dtdir") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0] && argv[z][0] != '-') {
                sscanf(argv[z], "%lf", &antenna_dt_direction);

                if (antenna_dt_direction < 0.0)
                    antenna_dt_direction = 0.0;
                if (antenna_dt_direction > 359.0)
                    antenna_dt_direction = 0.0;
            }
        }

        if (strcmp(argv[x], "-o") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0] && argv[z][0] != '-') {
                // sanity check length
                if(strlen(argv[z]) < 5){
                    spdlog::error("Output name is too short. Must be at least 5 chars");
                    exit(1);
                }

                /* Antenna pattern files have the same basic name as the output file
                 * but with a different extension. If they exist, load them now */
                size_t base_len = strlen(argv[z]);
                if(base_len >= sizeof(mapfile)){
                    spdlog::error("Output name too long (max {} chars)", sizeof(mapfile)-2);
                    exit(1);
                }
                // Copy base name into mapfile and tx_site structures
                strncpy(mapfile, argv[z], sizeof(mapfile)-1);
                mapfile[sizeof(mapfile)-1] = '\0';
                output_filename = argv[z];

                const char *az_base = (antenna_file[0] != '\0') ? antenna_file : argv[z];
                const char *el_base = az_base; // same logic
                size_t az_needed = strlen(az_base) + strlen(AZ_FILE_SUFFIX) + 1;
                size_t el_needed = strlen(el_base) + strlen(EL_FILE_SUFFIX) + 1;

                az_filename = (char*)calloc(az_needed, sizeof(char));
                if(az_filename == NULL) return ENOMEM;
                el_filename = (char*)calloc(el_needed, sizeof(char));
                if(el_filename == NULL){
                    free(az_filename);
                    return ENOMEM;
                }
                snprintf(az_filename, az_needed, "%s%s", az_base, AZ_FILE_SUFFIX);
                snprintf(el_filename, el_needed, "%s%s", el_base, EL_FILE_SUFFIX);

                if( (result = LoadPAT(az_filename,el_filename)) != 0 ){
                    spdlog::error("Permissions error reading antenna pattern file");
                    free(az_filename);
                    free(el_filename);
                    exit(result);
                }
                free(az_filename);
                free(el_filename);
            } else if (z <= y && argv[z][0] && argv[z][0] == '-' && argv[z][1] == '\0' ) {
                /* Handle writing image data to stdout */
                to_stdout = true;
                mapfile[0] = '\0';
                output_filename.clear();
                spdlog::error("Writing data to stdout");
            }
        }

        if (strcmp(argv[x], "-rt") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0])	/* A minus argument is legal here */
                sscanf(argv[z], "%d", &contour_threshold);
        }

        if (strcmp(argv[x], "-hd") == 0) {
            spdlog::info("    hd mode");
            free_elev();
            free_path();
            free_dem();
            ippd = 3600;
            ARRAYSIZE = (MAX_DISTANCE_DEGRES * ippd) + 10;
            do_allocs();
            spdlog::info("    Built for {} ppd", ippd);
        }

        if (strcmp(argv[x], "-t") == 0) {
            ngs = false;	// greyscale background
        }

        if (strcmp(argv[x], "-dbm") == 0)
            dbm = true;

        if (strcmp(argv[x], "-geotiff") == 0) {
            geotiff = true;
        }

        if (strcmp(argv[x], "-dem") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0] && argv[z][0] != '-') {
                strncpy(DEM_path, argv[z], 253);
                DEM_path[253] = '\0';
            }
        }
        
        if (strcmp(argv[x], "-lat") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0]) {
                tx_site.lat = atof(argv[z]);
            }
        }
        if (strcmp(argv[x], "-lon") == 0) {
            z = x + 1;
            if (z <= y && argv[z][0]) {
                tx_site.lon = atof(argv[z]);
            }
        }
        //Switch to Path Profile Mode if Rx co-ords specified
        if (strcmp(argv[x], "-rla") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0]) {
                ppa = 1;
                rx_site.lat = atof(argv[z]);

            }
        }
        if (strcmp(argv[x], "-rlo") == 0) {
            z = x + 1;
            if (z <= y && argv[z][0]) {
                rx_site.lon = atof(argv[z]);
            }
        }

        if (strcmp(argv[x], "-txh") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0] && argv[z][0] != '-') {
                sscanf(argv[z], "%f", &tx_site.alt);

            }
        }

        if (strcmp(argv[x], "-rxh") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0] && argv[z][0] != '-') {
                sscanf(argv[z], "%lf", &altitudeLR);
                sscanf(argv[z], "%f", &rx_site.alt);
            }
        }

        if (strcmp(argv[x], "-rxg") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0] && argv[z][0] != '-') {
                sscanf(argv[z], "%lf", &rxGain);
            }
        }

        if (strcmp(argv[x], "-f") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0] && argv[z][0] != '-') {
                sscanf(argv[z], "%lf", &LR.frq_mhz);
            }
        }

        if (strcmp(argv[x], "-erp") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0] && argv[z][0] != '-') {
                sscanf(argv[z], "%lf", &LR.erp);
            }
        }

        if (strcmp(argv[x], "-cl") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0] && argv[z][0] != '-') {

                sscanf(argv[z], "%d", &LR.radio_climate);

            }
        }
        if (strcmp(argv[x], "-te") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0] && argv[z][0] != '-') {
                int ter;
                sscanf(argv[z], "%d", &ter);

                switch (ter) {
                case 1:	// Water
                    terdic = 80;
                    tercon = 0.010;
                    break;

                case 2:	// Marsh
                    terdic = 12;
                    tercon = 0.007;
                    break;

                case 3:	// Farmland
                    terdic = 15;
                    tercon = 0.005;
                    break;

                case 4:	//Mountain
                    terdic = 13;
                    tercon = 0.002;
                    break;
                case 5:	//Desert
                    terdic = 13;
                    tercon = 0.002;
                    break;
                case 6:	//Urban
                    terdic = 5;
                    tercon = 0.001;
                    break;
                }
                LR.eps_dielect = terdic;
                LR.sgm_conductivity = tercon;

            }
        }

        if (strcmp(argv[x], "-terdic") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0] && argv[z][0] != '-') {

                sscanf(argv[z], "%lf", &terdic);

                LR.eps_dielect = terdic;

            }
        }

        if (strcmp(argv[x], "-tercon") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0] && argv[z][0] != '-') {

                sscanf(argv[z], "%lf", &tercon);

                LR.sgm_conductivity = tercon;

            }
        }

        if (strcmp(argv[x], "-hp") == 0) {
            // Horizontal polarisation (0)
            LR.pol = 0;
        }

        if (strcmp(argv[x], "-dbg") == 0) {
            debug = 1;
        }

        /*Prop model */

        if (strcmp(argv[x], "-pm") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0]) {
                int temp;
                sscanf(argv[z], "%d", &temp);
                prop_model = (PropModel)temp;
            }
        }
        // Prop model variant eg. urban/suburban
        if (strcmp(argv[x], "-pe") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0]) {
                sscanf(argv[z], "%d", &pmenv);
            }
        }
        //Knife edge diffraction
        if (strcmp(argv[x], "-ked") == 0) {
            z = x + 1;
            knifeedge = 1;
        }

        //Normalise Path Profile chart
        if (strcmp(argv[x], "-ng") == 0) {
            z = x + 1;
            normalise = 1;
        }

        // Reliability % for ITM model
        if (strcmp(argv[x], "-rel") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0]) {
                sscanf(argv[z], "%lf", &LR.rel);
                LR.rel=LR.rel/100;
            }
        }
        // Confidence % for ITM model
        if (strcmp(argv[x], "-conf") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0]) {
                sscanf(argv[z], "%lf", &LR.conf);
                LR.conf=LR.conf/100;
            }
        }
        // LossColors for the -scf, -dcf and -lcf, depending on mode
        if (strcmp(argv[x], "-color") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0]) {
                color_file = (char*) calloc(PATH_MAX+1, sizeof(char));
                if (color_file == NULL)
                    return ENOMEM;
                strncpy(color_file, argv[z], 253);
                color_file[253] = '\0';
            }
        }

        // number_threads to divide plot by
        if (strcmp(argv[x], "-number_threads") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0]) {
                sscanf(argv[z], "%d", &number_threads);
            }
        }
    }

    if (debug) {
        spdlog::set_level(spdlog::level::debug);
        spdlog::debug("Debug logging enabled");
    } else {
        spdlog::set_level(spdlog::level::info);
    }

    /* ERROR DETECTION */
    if (tx_site.lat > 90 || tx_site.lat < -90) {
        spdlog::error("Either the lat was missing or out of range!");
        exit(EINVAL);

    }
    if (tx_site.lon > 180.0 || tx_site.lon < -180.0) {
        spdlog::error("Either the lon was missing or out of range! (expected -180 to +180)");
        exit(EINVAL);

    }
    if (LR.frq_mhz < 20 || LR.frq_mhz > 100000) {
        spdlog::error("Either the Frequency was missing or out of range!");
        exit(EINVAL);
    }
    if (LR.erp > 500000000) {
        spdlog::error("Power was out of range!");
        exit(EINVAL);

    }
    if (LR.eps_dielect > 80 || LR.eps_dielect < 0.1) {
        spdlog::error("Ground Dielectric value out of range!");
        exit(EINVAL);

    }
    if (LR.sgm_conductivity > 0.01 || LR.sgm_conductivity < 0.000001) {
        spdlog::error("Ground conductivity out of range!");
        exit(EINVAL);

    }

    if (tx_site.alt < 0 || tx_site.alt > 60000) {
        spdlog::error("Tx altitude above ground was too high: {}",
            tx_site.alt);
        exit(EINVAL);
    }
    if (altitudeLR < 0 || altitudeLR > 60000) {
        spdlog::error("Rx altitude above ground was too high!");
        exit(EINVAL);
    }

    if (ippd < 300 || ippd > 10000) {
        spdlog::error("resolution out of range!");
        exit(EINVAL);
    }

    if (contour_threshold < -200 || contour_threshold > 240) {
        spdlog::error("Receiver threshold out of range (-200 / +240)");
        exit(EINVAL);
    }
    if (prop_model > 2 && prop_model < 7 && LR.frq_mhz < 150) {
        spdlog::error("Frequency too low for Propagation model");
        exit(EINVAL);
    }

    if (to_stdout == true && ppa != 0) {
        spdlog::error("Cannot write to stdout in ppa mode");
        exit(EINVAL);
    }

    /* Ensure a trailing '/' is present in DEM_path */

    if (DEM_path[0]) {
        x = strlen(DEM_path);

        if (DEM_path[x - 1] != '/' && x != 0) {
            spdlog::debug("Appending / to Copernicus directory");
            DEM_path[x] = '/';
            DEM_path[x + 1] = 0;
        }
    }

    if (number_threads % 2 != 0 || number_threads < 4) {
        spdlog::error("Number of number_threads must be even and greater than 4");
    }

    spdlog::info("-------------------------------- Plot Information --------------------------------");
    spdlog::info("    TX site parameters: {:.6f}N, {:.6f}W, {:.0f} m AGL", tx_site.lat, tx_site.lon, tx_site.alt);
    spdlog::info("    Plot parameters: {:.2f}-km radius, resolution of {} ppd", max_range, ippd);
    spdlog::info("    Model parameters: {} MHz at {} W EIRP (dBd), {}% confidence", LR.frq_mhz, LR.erp, (uint8_t)(LR.conf * 100));
    spdlog::info("    Map number_threads: {}", number_threads);
    spdlog::info("");
    spdlog::info("    Directories:");
    spdlog::info("        DEM: {}", DEM_path);
    spdlog::info(VERT_SEP);

    /**
     * Calculate the required data bounds to the nearest whole degree using WGS 84 approximation
     * https://en.wikipedia.org/wiki/Geographic_coordinate_system#Length_of_a_degree
    */

    // Get latitude in radians
    double tx_lat_rad = tx_site.lat * DEG2RAD;

    // Find the distance in lat and lon per degree using the above referenced formulas
    double m_per_deg_lon = (111412.84 * cos(tx_lat_rad)) - (93.5 * cos(3 * tx_lat_rad)) + (0.118 * cos(5 * tx_lat_rad));
    double m_per_deg_lat = 111132.92 - (559.82 * cos(2 * tx_lat_rad)) + (1.175 * cos(4 * tx_lat_rad)) - (0.0023 * cos(6 * tx_lat_rad));

    // Calculate angular distance from the above numbers
    double dist_deg_lon = (max_range * 1000.0) / m_per_deg_lon;
    double dist_deg_lat = (max_range * 1000.0) / m_per_deg_lat;

    spdlog::debug("Radius of {:.3f} km is approx {:.6f} deg EW and {:.6f} deg NS", max_range, dist_deg_lon, dist_deg_lat);

    // Calculate our plot bounds based on these numbers
    min_lon = tx_site.lon - dist_deg_lon;
    max_lon = tx_site.lon + dist_deg_lon;
    min_lat = tx_site.lat - dist_deg_lat;
    max_lat = tx_site.lat + dist_deg_lat;

    // If doing P2P analysis, we need to make sure the RX site is within our whole degree bounds as well, so data is loaded
    // TODO: update this so it makes sense with the new approach

    if (ppa == 1) {
        rxlat = rx_site.lat;
        rxlon = rx_site.lon;

        if (rxlat < min_lat)
            min_lat = rxlat;

        if (rxlat > max_lat)
            max_lat = rxlat;

        if (LonDiff(rxlon, min_lon) < 0.0)
            min_lon = rxlon;

        if (LonDiff(rxlon, max_lon) >= 0.0)
            max_lon = rxlon;

        spdlog::debug("RX site location expanded plot bounds to {:.6f}N {:.6f}E to {:.6f}N {:.6f}E", min_lat, min_lon, max_lat, max_lon);
    }

    bbox plot_bounds;
    plot_bounds.lower_right = {min_lat, min_lon};
    plot_bounds.upper_left = {max_lat, max_lon};

    spdlog::debug("Calculated plot boundaries: {:.6f}N {:.6f}E to {:.6f}N {:.6f}E", 
        plot_bounds.lower_right.lat, 
        plot_bounds.lower_right.lon, 
        plot_bounds.upper_left.lat, 
        plot_bounds.upper_left.lon
    );

    /* Load the required tiles */
    // DEM first

    if( (result = LoadTopoData(plot_bounds)) != 0 ){
        // This only fails on errors loading DEM tiles
        spdlog::error("Error loading topo data");
        return result;
    }

    ppd=(double)ippd;
    samples_per_radian = ppd * (180.0 / PI);

    width = (unsigned)(ippd * (max_lon - min_lon));
    height = (unsigned)(ippd * (max_north - min_north));

    dpp = 1 / ppd;
    mpi = ippd-1; 

    if (ppa == 0) {
        if (prop_model == LOS) {  // Model 2 = LOS
            cropping = false; // TODO: File is written in DoLOS() so this needs moving to PlotPropagation() to allow styling, cropping etc
            PlotLOSMap(tx_site, altitudeLR, number_threads);
            DoLOS(mapfile);
        } else {
            // 90% of effort here
            PlotPropagationRadius(tx_site, max_range, altitudeLR, prop_model, (uint8_t)number_threads);
            spdlog::debug("Finished PlotPropagationRadius()");

            if (cropping) {
                // CROPPING Factor determined in propPathLoss().
                // cropLon is the circle radius in pixels at it's widest (east/west) 
                cropLon*=dpp; // pixels to degrees
                max_north=cropLat; // degrees
                min_lon = tx_site.lon - cropLon; // western crop boundary
                cropLat-=tx_site.lat; // angle from tx to edge

                spdlog::debug("Cropping 1: min_lon: {:.4f} cropLat: {:.4f} cropLon: {:.4f} longitude: {:.4f} dpp {:.5f}",min_lon,cropLat,cropLon,tx_site.lon,dpp);

                width=(int)((cropLon*ppd)*2);
                height=(int)((cropLat*ppd)*2);

                spdlog::debug("Cropping 2: min_lon: {:.4f} cropLat: {:.4f} cropLon: {:.4f} longitude: {:.4f} width {:d}",min_lon,cropLat,cropLon,tx_site.lon,width);

                if (width > 3600 * 10 || cropLon < 0) {
                    spdlog::error("FATAL BOUNDS! min_lon: {:.4f} cropLat: {:.4f} cropLon: {:.7f} longitude: {:.5f}",min_lon,cropLat,cropLon,tx_site.lon);
                    return 0;
                }
            }

            // Write bitmap
            if (LR.erp == 0.0)
                DoPathLoss(mapfile);
            else if (dbm)
                DoRxdPwr((to_stdout == true ? NULL : mapfile));
            else
                    if ((result = DoSigStr(mapfile)) != 0)
                    return result;
        }


        if (cropping) {
            spdlog::info("Area boundaries:{:.6f} | {:.6f} | {:.6f} | {:.6f} ", tx_site.lat+cropLat, tx_site.lon+cropLon, tx_site.lat-cropLat,tx_site.lon-cropLon);
        } else {
            spdlog::info("Area boundaries:{:.6f} | {:.6f} | {:.6f} | {:.6f} ",max_north,east,min_north,west);
        }

    } else {
        PlotPath(tx_site, rx_site);
        PathReport(tx_site, rx_site, output_filename.c_str(), 0, prop_model, rxGain);
        // Order flipped for benefit of graph. Makes no difference to data.
        SeriesData(rx_site, tx_site, output_filename.c_str(), 1, normalise);
    }
    fflush(stderr);

    auto end_time = std::chrono::steady_clock::now();
    double elapsed_s = std::chrono::duration<double>(end_time - start_time).count();
    fprintf(stderr, "Execution time: %.3f seconds\n", elapsed_s);

    if (debug) {
        fprintf(stderr, "--- Function call counters ---\n");
        fprintf(stderr, "  point_to_point_ITM : %d\n", cnt_point_to_point_ITM.load());
        fprintf(stderr, "  point_to_point     : %d\n", cnt_point_to_point.load());
        fprintf(stderr, "  computeLoss        : %d\n", cnt_computeLoss.load());
        fprintf(stderr, "  PlotPropPath       : %d\n", cnt_PlotPropPath.load());
        fprintf(stderr, "  ReadPath           : %d\n", cnt_ReadPath.load());
        fprintf(stderr, "  PlotPropagation    : %d\n", cnt_PlotPropagation.load());
    }

    return 0;
}
