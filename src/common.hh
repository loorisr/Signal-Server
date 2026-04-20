#ifndef _COMMON_HH_
#define _COMMON_HH_

#include <atomic>
#include <cstdint>
#include <string>

#define GAMMA 		2.5
#define ONE_OVER_GAMMA 	(1.0 / GAMMA)


#ifndef PI
  #define PI		3.141592653589793f
#endif

#ifndef TWOPI
  #define TWOPI		6.283185307179586f
#endif

#ifndef HALFPI
  #define HALFPI	1.570796326794896f
#endif

// Multiplier to convert decimal degrees to radians
#define DEG2RAD		1.74532925199e-02f
#define RAD2DEG		57.2957795130823f
// Radius of the earth, in meters
#define	EARTHRADIUS	6371000.0f
#define FOUR_THIRDS	1.3333333333333f
#define KM_PER_DEG_LAT 110.754f

#define	FOUR_THIRDS_EARTH (FOUR_THIRDS * EARTHRADIUS)
//#define MAX(x,y)((x)>(y)?(x):(y))


struct site {
	float lat;
	float lon;
	float alt;
};

struct path {
	float *lat;
	float *lon;
	float *elevation;
	float *distance;
	int length;
};

struct LR {
	float eps_dielect;
	float sgm_conductivity;
	float eno_ns_surfref;
	float frq_mhz;
	float conf;
	float rel;
	float erp;
	int radio_climate;
	int pol;
	float antenna_pattern[361][1001];
};

struct coord {
    float lat;
    float lon;
};

struct bbox {
    coord lower_left;   /* SW corner: min lat, min lon */
    coord upper_right;  /* NE corner: max lat, max lon */
};

enum PropagationMode {
    PROP_MODE_NONE = 0,
    PROP_MODE_LOS = 1 << 0,
    PROP_MODE_1_HRZN = 1 << 1,
    PROP_MODE_2_HRZN = 1 << 2,
    PROP_MODE_PEAK = 1 << 3,
    PROP_MODE_DIFFRACTION = 1 << 4,
    PROP_MODE_TROPOSCATTER = 1 << 5,
};

inline PropagationMode operator|(PropagationMode lhs, PropagationMode rhs)
{
    return static_cast<PropagationMode>(static_cast<int>(lhs) | static_cast<int>(rhs));
}

inline PropagationMode &operator|=(PropagationMode &lhs, PropagationMode rhs)
{
    lhs = lhs | rhs;
    return lhs;
}


extern float min_lat;
extern float max_lat;
extern float min_lon;
extern float max_lon;
extern int max_elevation;
extern int min_elevation;
extern int contour_threshold;
extern bool knifeedge;
extern int pmenv;
extern int number_threads;
extern int dh_n;
extern int fast_dh_stride;

extern bool ppa;
extern int normalise;
extern std::string mapfile;

extern float max_range;
extern int ppd;
extern float fzone_clearance;
extern float clutter;
extern float dBm;
extern float loss;
extern float field_strength;
extern float altitudeLR;

extern std::string DEM_path;
extern std::string color_palette;
extern std::string output_filename;

extern bool got_elevation_pattern;
extern bool got_azimuth_pattern;
extern bool dbm;
extern bool geotiff;
extern bool ngs;

extern struct site tx_site;
extern struct site rx_site;

extern float antenna_rotation;
extern float antenna_downtilt;
extern float antenna_dt_direction;
extern float rxGain;
extern float tercon;
extern float terdic;

/* Flat DEM arrays: dem_data[x][y], dem_signal[x][y]
 * x: 0 = global south edge, increases northward  (rows = ppd * tiles_lat)
 * y: 0 = global west  edge, increases eastward   (cols = ppd * tiles_lon)
 * Allocated by alloc_dem() after LoadTopoData knows the bounding box. */
extern float         **dem_data;
extern int **dem_signal;
/* Geographic origin of the flat arrays (south-west corner, integer degrees) */
extern int dem_min_lat;   /* southernmost tile_lat */
extern int dem_min_lon;   /* westernmost  tile_lon */
extern int dem_width_px;  /* tiles_lon * ppd */
extern int dem_height_px; /* tiles_lat * ppd */
extern thread_local float *elev;
extern thread_local int elev_allocated;
extern thread_local struct path path;
extern thread_local int path_allocated;
extern struct LR LR;

extern int debug;

extern std::atomic<int> cnt_point_to_point_ITM;
extern std::atomic<int> cnt_point_to_point;
extern std::atomic<int> cnt_computeLoss;
extern std::atomic<int> cnt_PlotPropPath;
extern std::atomic<int> cnt_ReadPath;
extern std::atomic<int> cnt_PlotPropagation;

#endif /* _COMMON_HH_ */
