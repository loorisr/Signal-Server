#ifndef _COMMON_HH_
#define _COMMON_HH_

#include <atomic>
#include <cstdint>
#include <string>

#define GAMMA 		2.5
#define ONE_OVER_GAMMA 	(1.0 / GAMMA)


#ifndef PI
  #define PI		3.141592653589793
#endif

#ifndef TWOPI
  #define TWOPI		6.283185307179586
#endif

#ifndef HALFPI
  #define HALFPI	1.570796326794896
#endif

// Multiplier to convert decimal degrees to radians
#define DEG2RAD		1.74532925199e-02
#define RAD2DEG		57.2957795130823
// Radius of the earth, in meters
#define	EARTHRADIUS	6371000.0
#define FOUR_THIRDS	1.3333333333333
#define KM_PER_DEG_LAT 110.754

#define	FOUR_THIRDS_EARTH (FOUR_THIRDS * EARTHRADIUS)
//#define MAX(x,y)((x)>(y)?(x):(y))


struct site {
	double lat;
	double lon;
	float alt;
};

struct path {
	double *lat;
	double *lon;
	double *elevation;
	double *distance;
	int length;
};

struct LR {
	double eps_dielect;
	double sgm_conductivity;
	double eno_ns_surfref;
	double frq_mhz;
	double conf;
	double rel;
	double erp;
	int radio_climate;
	int pol;
	float antenna_pattern[361][1001];
};

struct region {
	unsigned char color[128][3];
	int level[128];
	int levels;
};

struct coord {
    double lat;
    double lon;
};

struct bbox {
    coord lower_right;
    coord upper_left;
};

extern int MAX_DISTANCE_DEGRES;
extern int ARRAYSIZE;

extern double min_north;
extern double max_north;
extern double min_lon;
extern double max_lon;
extern double min_lat;
extern double max_lat;
extern int ippd;
extern int mpi;
extern int max_elevation;
extern int min_elevation;
extern int contour_threshold;
extern int width;
extern int height;
extern int knifeedge;
extern int pmenv;

extern double north;
extern double east;
extern double south;
extern double west;
extern double max_range;
extern double dpp;
extern double ppd;
extern double samples_per_radian;
extern double fzone_clearance;
extern double clutter;
extern double dBm;
extern double loss;
extern double field_strength;
extern __thread double *elev;
extern double delta;

extern char string[];
extern char DEM_path[];
extern std::string output_filename;

extern bool got_elevation_pattern;
extern bool got_azimuth_pattern;
extern bool dbm;
extern bool geotiff;
extern bool ngs;

/* Flat DEM arrays: dem_data[x][y], dem_signal[x][y]
 * x: 0 = global south edge, increases northward  (rows = ippd * tiles_lat)
 * y: 0 = global east  edge, increases westward   (cols = ippd * tiles_lon)
 * Allocated by alloc_dem() after LoadTopoData knows the bounding box. */
extern double         **dem_data;
extern int **dem_signal;
/* Geographic origin of the flat arrays (south-west corner, integer degrees) */
extern int dem_min_lat;   /* southernmost tile_lat */
extern int dem_min_lon;   /* westernmost  tile_lon */
extern int dem_width_px;  /* tiles_lon * ippd */
extern int dem_height_px; /* tiles_lat * ippd */
extern __thread struct path path;
extern struct LR LR;
extern struct region region;

extern int debug;

extern std::atomic<int> cnt_point_to_point_ITM;
extern std::atomic<int> cnt_point_to_point;
extern std::atomic<int> cnt_computeLoss;
extern std::atomic<int> cnt_PlotPropPath;
extern std::atomic<int> cnt_ReadPath;
extern std::atomic<int> cnt_PlotPropagation;

#endif /* _COMMON_HH_ */
