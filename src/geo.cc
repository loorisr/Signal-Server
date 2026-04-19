#include <algorithm>
#include <math.h>

#include <spdlog/spdlog.h>

#include "geo.hh"
#include "main.hh"
#include "common.hh"


// Compute the wrapped longitude difference between two coordinates.
float LonDiff(float lon1, float lon2)
{
    float diff = lon1 - lon2;
    if (diff <= -180.0) return diff + 360.0;
    if (diff >= 180.0)  return diff - 360.0;
    return diff;
}

// Compute great-circle distance between two sites.
float Distance(struct site site1, struct site site2)
{
    float lat1 = site1.lat * DEG2RAD;
    float lon1 = site1.lon * DEG2RAD;
    float lat2 = site2.lat * DEG2RAD;
    float lon2 = site2.lon * DEG2RAD;

    float dot = sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2);
    dot = std::clamp(dot, -1.0f, 1.0f);  /* guard against floating-point overshoot */

    return EARTHRADIUS * acos(dot);
}

// Compute the forward azimuth from the source site to the destination site.
float Azimuth(struct site source, struct site destination)
{
    const float src_lat  = source.lat      * DEG2RAD;
    const float src_lon  = source.lon      * DEG2RAD;
    const float dest_lat = destination.lat * DEG2RAD;
    const float dest_lon = destination.lon * DEG2RAD;
    const float dlon     = dest_lon - src_lon;

    float az = atan2(sin(dlon) * cos(dest_lat),
                      cos(src_lat) * sin(dest_lat) - sin(src_lat) * cos(dest_lat) * cos(dlon));
    if (az < 0.0) az += TWOPI;
    return az * RAD2DEG;
}

// Return the WGS84 Earth radius at the given latitude.
float earthRadius(float lat)
{
    // Convert latitude to rad
    float lat_rad = lat * DEG2RAD;

    float An = WGS84_a * WGS84_a * cos(lat_rad);
    float Bn = WGS84_b * WGS84_b * sin(lat_rad);
    float Ad = WGS84_a * cos(lat_rad);
    float Bd = WGS84_b * sin(lat_rad);

    return sqrt((An*An + Bn*Bn) / (Ad*Ad + Bd*Bd));
}

// Project a point from a center coordinate by distance and bearing.
coord getPointAtDistance(const coord center, float distance, float bearing_rad)
{
    coord endCoords;

    // 1. Pre-calculate radians and Earth radius ratio
    const float start_lat_rad = center.lat * DEG2RAD;
    const float start_lon_rad = center.lon * DEG2RAD;
    const float dR            = distance / earthRadius(center.lat);

    // 2. Pre-calculate common trig values to avoid redundant calls
    const float sin_lat = sin(start_lat_rad);
    const float cos_lat = cos(start_lat_rad);
    const float sin_dR  = sin(dR);
    const float cos_dR  = cos(dR);
    const float sin_brg = sin(bearing_rad);
    const float cos_brg = cos(bearing_rad);

    // 3. Calculate Latitude
    const float sin_end_lat = sin_lat * cos_dR + cos_lat * sin_dR * cos_brg;
    const float end_lat_rad = asin(sin_end_lat);
    
    // 4. Calculate Longitude
    // Optimization: Reuse sin_end_lat to avoid another sin() call
    const float end_lon_rad = start_lon_rad + atan2(
        sin_brg * sin_dR * cos_lat,
        cos_dR - sin_lat * sin_end_lat
    );

    endCoords.lat = end_lat_rad * RAD2DEG;
    endCoords.lon = end_lon_rad * RAD2DEG;

    return endCoords;
}

// Build a latitude/longitude bounding box that encloses a circular area.
bbox getCircularBoundingBox(coord center, float radius)
{
    // Result bbox
    bbox result;

    // Convert input degrees to rads
    float lat_rad = center.lat * DEG2RAD;
    float lon_rad = center.lon * DEG2RAD;

    // Get earth's radius at the specified latitude (km)
    float e_rad = earthRadius(center.lat);

    // Get parallel radius at latitude (km)
    float p_rad = e_rad * cos(lat_rad);

    // Calculate bounds (radians)
    float latMin = lat_rad - (radius / e_rad);
    float latMax = lat_rad + (radius / e_rad);
    float lonMin = lon_rad - (radius / p_rad);
    float lonMax = lon_rad + (radius / p_rad);

    result.lower_left = { (float)(latMin * RAD2DEG), (float)(lonMin * RAD2DEG) };
    result.upper_right = { (float)(latMax * RAD2DEG), (float)(lonMax * RAD2DEG) };

    return result;
}

/* DEM access */

// Translate geographic coordinates into DEM pixel coordinates.
bool find_dem_xy(float lat, float lon, int &x_out, int &y_out)
{
    if (!dem_data) return false;
    int x = (int)rint(ppd * (lat - dem_min_lat));
    int y = (int)rint(ppd * (lon - dem_min_lon));
    if (x < 0 || x >= dem_height_px || y < 0 || y >= dem_width_px) return false;
    x_out = x;
    y_out = y;
    return true;
}

// Store the computed signal value at the DEM cell for a location.
void PutSignal(float lat, float lon, int signal)
{
    int x, y;
    if (find_dem_xy(lat, lon, x, y))
        dem_signal[x][y] = signal;
}

// Read the stored signal value at the DEM cell for a location.
int GetSignal(float lat, float lon)
{
    int x, y;
    if (!find_dem_xy(lat, lon, x, y)) return 0;
    return dem_signal[x][y];
}

// Read the terrain elevation for a site from the loaded DEM.
float GetElevation(struct site location)
{
    int x, y;
    if (!find_dem_xy(location.lat, location.lon, x, y)) return -5000.0;
    return dem_data[x][y];
}

/* Path / elevation geometry */

float ElevationAngle(struct site source, struct site destination)
{
    /* This function returns the angle of elevation (in degrees)
       of the destination as seen from the source location.
       A positive result represents an angle of elevation (uptilt),
       while a negative result represents an angle of depression
       (downtilt), as referenced to a normal to the center of
       the earth. */

    float a, b, dx;

    a = GetElevation(destination) + destination.alt + EARTHRADIUS;
    b = GetElevation(source) + source.alt + EARTHRADIUS;

    dx = Distance(source, destination);

    /* Apply the Law of Cosines */

    float cos_angle = std::clamp(((b * b) + (dx * dx) - (a * a)) / (2.0f * b * dx), -1.0f, 1.0f);
    return (acos(cos_angle) * RAD2DEG) - 90.0;
}

// Allocate the per-thread path buffers for a new sample count.
static void alloc_path(int size)
{
    path.lat = new float[size];
    path.lon = new float[size];
    path.elevation = new float[size];
    path.distance = new float[size];
    path_allocated = size;
}

// Generate the sampled great-circle path and cache terrain samples along it.
void ReadPath(struct site source, struct site destination)
{
    if (debug) cnt_ReadPath++;

    float azimuth, distance, lat2, lon2, beta,
        total_distance, dx, dy, path_length, m_per_sample;

    const float lat1 = source.lat * DEG2RAD;
    const float lon1 = source.lon * DEG2RAD;

    azimuth        = Azimuth(source, destination) * DEG2RAD;
    total_distance = Distance(source, destination);

    dx          = ppd * RAD2DEG * (lon1 - destination.lon * DEG2RAD);
    dy          = ppd * RAD2DEG * (lat1 - destination.lat * DEG2RAD);
    path_length = sqrt((dx * dx) + (dy * dy));
    m_per_sample = total_distance / path_length;

    /* +2: one extra for the destination point written after the loop, one for safety */
    int needed = (int)ceil(path_length) + 2;
    if (needed > path_allocated) {
        free_path();
        alloc_path(needed);
    }
    /* elev needs path.length + 2 slots (indices 0..path.length+1) */
    if (needed + 2 > elev_allocated) {
        delete [] elev;
        elev = new float[needed + 2];
        elev_allocated = needed + 2;
    }

    /* Hoist loop-invariant trig */
    const float sin_lat1 = sin(lat1);
    const float cos_lat1 = cos(lat1);
    const float cos_az   = cos(azimuth);
    const float sin_az   = sin(azimuth);

    int c = 0;
    if (total_distance != 0.0) {
        for (distance = 0.0; distance <= total_distance; c++, distance = m_per_sample * c) {
            beta = distance / EARTHRADIUS;
            const float cos_beta  = cos(beta);
            const float sin_beta  = sin(beta);
            const float sin_lat2  = std::clamp(sin_lat1 * cos_beta + cos_az * sin_beta * cos_lat1, -1.0f, 1.0f);
            lat2 = asin(sin_lat2) * RAD2DEG;
            lon2 = (lon1 + atan2(sin_az * sin_beta * cos_lat1, cos_beta - sin_lat1 * sin_lat2)) * RAD2DEG;

            path.lat[c]       = lat2;
            path.lon[c]       = lon2;
            path.elevation[c] = GetElevation({lat2, lon2, 0.0});
            path.distance[c]  = distance;
        }
    }

    /* Make sure exact destination point is recorded at path.length-1 */
    path.lat[c]       = destination.lat;
    path.lon[c]       = destination.lon;
    path.elevation[c] = GetElevation(destination);
    path.distance[c]  = total_distance;

    path.length = c + 1;
}

// Return the destination elevation angle or the first obstruction angle.
float ElevationAngle2(struct site source, struct site destination, float er)
{
    /* This function returns the angle of elevation (in degrees)
       of the destination as seen from the source location, UNLESS
       the path between the sites is obstructed, in which case, the
       elevation angle to the first obstruction is returned instead.
       "er" represents the earth radius. */

    int x;
    char block = 0;
    float source_alt, destination_alt, cos_xmtr_angle,
        cos_test_angle, test_alt, elevation, distance,
        source_alt2, first_obstruction_angle = 0.0;

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
                (acos(std::clamp(cos_test_angle, -1.0f, 1.0f)) * RAD2DEG) - 90.0;
        }
    }

    if (block)
        elevation = first_obstruction_angle;

    else
        elevation = (acos(std::clamp(cos_xmtr_angle, -1.0f, 1.0f)) * RAD2DEG) - 90.0;

    return elevation;
}

// Report terrain obstructions and Fresnel-zone clearance along a path.
void ObstructionAnalysis(struct site xmtr, struct site rcvr, float f,
             FILE *outfile)
{
    /* Perform an obstruction analysis along the
       path between receiver and transmitter. */

    int x;
    float h_r, h_t, h_x, h_r_orig, cos_tx_angle, cos_test_angle,
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
        h_x = path.elevation[x] + EARTHRADIUS + clutter;
        d_x = d_tx - path.distance[x];

        /* Deal with the LOS path first. */

        cos_test_angle =
            ((h_r * h_r) + (d_x * d_x) -
             (h_x * h_x)) / (2.0 * h_r * d_x);

        if (cos_tx_angle > cos_test_angle) {
            if (h_r == h_r_orig)
                fprintf(outfile,
                    "Between RX and TX, obstructions were detected :\n\n");

            if (path.lat[x] >= 0.0) {
                fprintf(outfile,
                    "   %8.4f N,%9.4f W, %5.2f kilometers, %6.2f meters AMSL\n",
                    path.lat[x], path.lon[x],
                    d_x / 1000.0,
                    h_x - EARTHRADIUS);
            } else {
                fprintf(outfile,
                    "   %8.4f S,%9.4f W, %5.2f kilometers, %6.2f meters AMSL\n",
                    -path.lat[x], path.lon[x],
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

/* Memory management */

// Release the loaded DEM and signal grids.
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

// Release the per-thread elevation buffer.
void free_elev(void) {
  delete [] elev;
  elev_allocated = 0;
}

// Release the per-thread path buffers.
void free_path(void)
{
    delete [] path.lat;
    delete [] path.lon;
    delete [] path.elevation;
    delete [] path.distance;
    path_allocated = 0;
}


/* Allocate the flat DEM arrays for a known bounding box.
 * Called from LoadTopoData() after tiles_lat/tiles_lon are computed. */
// Allocate the flat DEM arrays for a known bounding box.
void alloc_dem(int min_lat, int min_lon, int tiles_lat, int tiles_lon)
{
    if (dem_data) free_dem();

    dem_min_lat   = min_lat;
    dem_min_lon   = min_lon;
    dem_height_px = tiles_lat * ppd;
    dem_width_px  = tiles_lon * ppd;

    dem_data   = new float         *[dem_height_px];
    dem_signal = new int *[dem_height_px];
    for (int i = 0; i < dem_height_px; i++) {
        dem_data[i]   = new float        [dem_width_px]();
        dem_signal[i] = new int[dem_width_px]();
    }
}

// TODO: temporary test — delete after validation
void test_Azimuth(void)
{
    struct { const char *name; site a; site b; float expected; } cases[] = {
        // Due North: az should be 0
        {"N",  {45.0, 0.0, 0.0}, {46.0, 0.0, 0.0}, 0.0},
        {"E",  {0.0, 5.0, 0.0}, {0.0, 7.0, 0.0}, 90.0},
        {"S",  {45.0, 5.0, 0.0}, {44.0, 5.0, 0.0}, 180.0},
        {"W",  {45.0, 7.0, 0.0}, {45.0, 5.0, 0.0}, 270.0},
        {"NE", {45.7, 6.4, 0.0}, {45.8, 6.5, 0.0}, 45.0},
    };

    spdlog::info("--- Azimuth unit tests ---");
    bool all_ok = true;
    for (auto &c : cases) {
        float got = Azimuth(c.a, c.b);
        float err = fabs(got - c.expected);
        bool ok = err < 1.0; // 1 degree tolerance
        spdlog::info("  {}: expected {:.1f}, got {:.2f} {}",
                     c.name, c.expected, got, ok ? "OK" : "FAIL");
        if (!ok) all_ok = false;
    }
    if (!all_ok)
        spdlog::error("Azimuth test FAILED");
    else
        spdlog::info("Azimuth test PASSED");
    spdlog::info("--------------------------");
}
