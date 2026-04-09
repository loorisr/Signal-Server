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
#include "cmdline.hh"
#include "inputs.hh"
#include "outputs.hh"
#include "models/itwom3.0.hh"
#include "models/los.hh"
#include "models/pel.hh"

#include <algorithm>
#include <chrono>
#include <thread>

#include <spdlog/spdlog.h>

int ippd = 1200;
int MAX_DISTANCE_DEGRES = 3; // max distance : 3° so around 300 km
int ARRAYSIZE = (MAX_DISTANCE_DEGRES * ippd) + 10;

char DEM_path[255];
std::string color_palette = "heat";
std::string output_filename;

double max_range = 0.0, dpp, ppd, samples_per_radian,
    fzone_clearance = 0.6, clutter, tercon, terdic,
    north, east, south, west, dBm, loss, field_strength,
    min_north = 90, max_north = -90, min_lon = 180.0, max_lon = -180.0,
    rxGain=0, antenna_rotation,
    antenna_downtilt, antenna_dt_direction;

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
bool cropping = true;
int knifeedge = 0;
int pmenv = 1;

__thread double *elev;
__thread struct path path;
site tx_site;
site rx_site;
double         **dem_data   = nullptr;
int **dem_signal = nullptr;
int dem_min_lat   = 0;
int dem_min_lon   = 0;
int dem_width_px  = 0;
int dem_height_px = 0;

struct LR LR;

/* Convert (lat, lon) to flat-array pixel coordinates.
 * x increases northward from the south edge of the loaded area.
 * y increases westward  from the east  edge of the loaded area.
 * Returns true and sets x_out/y_out when the point is inside the array. */
bool find_dem_xy(double lat, double lon, int &x_out, int &y_out)
{
    if (!dem_data) return false;
    int x = (int)rint(ppd * (lat - dem_min_lat));
    int y = (int)rint(ppd * (lon - dem_min_lon));
    if (x < 0 || x >= dem_height_px || y < 0 || y >= dem_width_px) return false;
    x_out = x;
    y_out = y;
    return true;
}

void PutSignal(double lat, double lon, int signal)
{
    int x, y;
    if (find_dem_xy(lat, lon, x, y))
        dem_signal[x][y] = signal;
}

int GetSignal(double lat, double lon)
{
    int x, y;
    if (!find_dem_xy(lat, lon, x, y)) return 0;
    return dem_signal[x][y];
}

double GetElevation(struct site location)
{
    int x, y;
    if (!find_dem_xy(location.lat, location.lon, x, y)) return -5000.0;
    return dem_data[x][y];
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

    double cos_angle = std::clamp(((b * b) + (dx * dx) - (a * a)) / (2.0 * b * dx), -1.0, 1.0);
    return (acos(cos_angle) * RAD2DEG) - 90.0;
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
        lat2 = asin(std::clamp(sin(lat1) * cos(beta) + cos(azimuth) * sin(beta) * cos(lat1), -1.0, 1.0));
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
        path.elevation[c] = GetElevation({lat2, lon2, 0.0f});
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
                (acos(std::clamp(cos_test_angle, -1.0, 1.0)) * RAD2DEG) - 90.0;
        }
    }

    if (block)
        elevation = first_obstruction_angle;

    else
        elevation = (acos(std::clamp(cos_xmtr_angle, -1.0, 1.0)) * RAD2DEG) - 90.0;

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
    dem_signal = new int *[dem_height_px];
    for (int i = 0; i < dem_height_px; i++) {
        dem_data[i]   = new double        [dem_width_px]();
        dem_signal[i] = new int[dem_width_px]();
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

void write_geotiff_rgba(const uint8_t *rgba, int img_width, int img_height, const char *filename)
{
    char tif_file[300];
    double ulx, uly, lrx, lry;

    /* Build .tif output path */
    snprintf(tif_file, sizeof(tif_file), "%s.tif", filename);

    /* Compute geographic bounds */
    ulx = min_lon;
    uly = max_north;
    lrx = max_lon;
    lry = min_north;

    GDALDriverH drv = GDALGetDriverByName("GTiff");
    if (drv == NULL) {
        spdlog::error("write_geotiff_rgba: GTiff GDAL driver not available");
        return;
    }

    char *create_opts[] = {
        (char *)"COMPRESS=LZW",
        (char *)"PREDICTOR=2",
        (char *)"TILED=YES",
        (char *)"BLOCKXSIZE=256",
        (char *)"BLOCKYSIZE=256",
        (char *)"INTERLEAVE=PIXEL",
        (char *)"NUM_THREADS=ALL_CPUS",
        NULL
    };
    GDALDatasetH ds = GDALCreate(drv, tif_file, img_width, img_height, 4, GDT_Byte, create_opts);
    if (ds == NULL) {
        spdlog::error("write_geotiff_rgba: failed to create {}", tif_file);
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

    /* Write interleaved RGBA buffer directly to 4 separate GDAL bands */
    int bandMap[4] = {1, 2, 3, 4};
    CPLErr err = (CPLErr)GDALDatasetRasterIO(ds, GF_Write,
        0, 0, img_width, img_height,
        (void *)rgba, img_width, img_height, GDT_Byte,
        4, bandMap,
        4,              /* nPixelSpace */
        img_width * 4,  /* nLineSpace  */
        1);             /* nBandSpace  */
    if (err != CE_None)
        spdlog::error("write_geotiff_rgba: RasterIO write failed for {}", tif_file);
    else
        spdlog::info("GeoTIFF written: {}", tif_file);

    GDALClose(ds);
}

int main(int argc, char *argv[])
{
    auto start_time = std::chrono::steady_clock::now();

    int result;

    CmdlineArgs args = parse_cmdline(argc, argv);
    int ppa             = args.ppa;
    int normalise       = args.normalise;
    int number_threads  = args.number_threads;
    double altitudeLR   = args.altitudeLR;
    PropModel prop_model = args.prop_model;
    char *mapfile       = args.mapfile;


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
    double dist_deg_lon = max_range / m_per_deg_lon;
    double dist_deg_lat = max_range / m_per_deg_lat;

    spdlog::debug("Radius of {:.3f} m is approx {:.6f} deg EW and {:.6f} deg NS", max_range, dist_deg_lon, dist_deg_lat);

    // Calculate our plot bounds based on these numbers
    min_lon = tx_site.lon - dist_deg_lon;
    max_lon = tx_site.lon + dist_deg_lon;
    double min_lat = tx_site.lat - dist_deg_lat;
    double max_lat = tx_site.lat + dist_deg_lat;

    // If doing P2P analysis, we need to make sure the RX site is within our whole degree bounds as well, so data is loaded
    // TODO: update this so it makes sense with the new approach

    if (ppa == 1) {
        if (rx_site.lat < min_lat)
            min_lat = rx_site.lat;

        if (rx_site.lat > max_lat)
            max_lat = rx_site.lat;

        if (LonDiff(rx_site.lon, min_lon) < 0.0)
            min_lon = rx_site.lon;

        if (LonDiff(rx_site.lon, max_lon) >= 0.0)
            max_lon = rx_site.lon;

        spdlog::debug("RX site location expanded plot bounds to {:.6f}N {:.6f}E to {:.6f}N {:.6f}E", min_lat, min_lon, max_lat, max_lon);
    }

    bbox plot_bounds = getCircularBoundingBox({tx_site.lat, tx_site.lon}, max_range);

    spdlog::debug("Calculated plot boundaries: {:.6f}N {:.6f}E to {:.6f}N {:.6f}E", 
        plot_bounds.lower_left.lat, 
        plot_bounds.lower_left.lon, 
        plot_bounds.upper_right.lat, 
        plot_bounds.upper_right.lon
    );

    /* Load the required DEM tiles */
    if( (result = LoadTopoData(plot_bounds)) != 0 ){
        // This only fails on errors loading DEM tiles
        spdlog::error("Error loading topo data");
        return result;
    }

    ppd=(double)ippd;
    samples_per_radian = ppd * RAD2DEG;

    width = (unsigned)(ippd * (max_lon - min_lon));
    height = (unsigned)(ippd * (max_north - min_north));

    dpp = 1.0 / ppd;
    mpi = ippd-1; 

    if (ppa == 0) {
        if (prop_model == LOS) {  // Model 2 = LOS
            PlotLOSMap(tx_site, max_range, altitudeLR, number_threads);
            DoLOS(mapfile);
        } else {
            // 90% of effort here
            PlotPropagationRadius(tx_site, max_range, altitudeLR, prop_model, (uint8_t)number_threads);
            spdlog::debug("Finished PlotPropagationRadius()");

            if (cropping) {
                spdlog::debug("Cropping 1: N: {:.4f} S: {:.4f} E: {:.4f} W: {:.4f} dpp {:.5f}",plot_bounds.upper_right.lat, plot_bounds.lower_left.lat, plot_bounds.upper_right.lon, plot_bounds.lower_left.lon, ppd);

                max_north = plot_bounds.upper_right.lat;
                min_north = plot_bounds.lower_left.lat;
                max_lon = plot_bounds.upper_right.lon;
                min_lon = plot_bounds.lower_left.lon;
                width = (unsigned)(ippd * (max_lon - min_lon));
                height = (unsigned)(ippd * (max_north - min_north));
                spdlog::info("width/height: {}/{}", width, height);
            }

            // Write bitmap
            if (LR.erp == 0.0)
                DoPathLoss(mapfile);
            else if (dbm)
                DoRxdPwr(mapfile);
            else if ((result = DoSigStr(mapfile)) != 0)
                return result;
        }


        spdlog::info("Area boundaries:{:.6f} | {:.6f} | {:.6f} | {:.6f} ",max_north,max_lon,min_north,min_lon);

    } else {
        PlotPath(tx_site, rx_site);
        PathReport(tx_site, rx_site, output_filename.c_str(), 0, prop_model, rxGain);
        // Order flipped for benefit of graph. Makes no difference to data.
        SeriesData(rx_site, tx_site, output_filename.c_str(), 1, normalise);
    }

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
