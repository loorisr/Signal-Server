/**
 * @file geo.hh
 * @ingroup geo
 * @file geo.cc
 * @ingroup geo
 * 
*/

#ifndef __GEO_HH_
#define __GEO_HH_

#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include "common.hh"

/// WGS84 semi-axes
#define WGS84_a 6378137.0
#define WGS84_b 6356752.3

double LonDiff(double lon1, double lon2);
double Distance(struct site site1, struct site site2);
double Azimuth(struct site source, struct site destination);

/// @brief Calculate the approximate radius of the earth at a given latitude, using the WGS84 model
/// @cite http://en.wikipedia.org/wiki/Earth_radius
/// @param lat latitude in degrees
/// @return earth radius in km
double earthRadius(double lat);

/// @brief Get the lat & lon of a point at a certain distance and heading from a given point
/// @cite Adapted from https://www.movable-type.co.uk/scripts/latlong.html and https://stackoverflow.com/a/7835325
/// @param start_lat starting latitude in degrees
/// @param start_lon starting longitude in degrees
/// @param distance distance in km
/// @param bearing bearing in radians
/// @return coodinates of the resultant point in decimal degrees
coord getPointAtDistance(coord center, double distance, double bearing_rad);

/// @brief Get the bounding box for a circle at a given lat/lon and radius
/// @param center center coordinates in decimal degrees
/// @param radius radius in km
/// @return bounding box
bbox getCircularBoundingBox(coord center, double radius);

/* DEM access */
bool find_dem_xy(double lat, double lon, int &x_out, int &y_out);
double GetElevation(struct site location);
void PutSignal(double lat, double lon, int signal);
int GetSignal(double lat, double lon);

/* Path / elevation geometry */
double ElevationAngle(struct site source, struct site destination);
double ElevationAngle2(struct site source, struct site destination, double er);
void ReadPath(struct site source, struct site destination);
void ObstructionAnalysis(struct site xmtr, struct site rcvr, double f, FILE *outfile);

/* Memory management */
void free_elev(void);
void free_path(void);
void alloc_dem(int min_lat, int min_lon, int tiles_lat, int tiles_lon);
void free_dem(void);

#endif