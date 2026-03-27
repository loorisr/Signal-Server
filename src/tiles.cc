#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include "tiles.hh"
#include "common.hh"

#define MAX_LINE 50000

/* Computes the distance between two long/lat points */
double haversine_formula(double th1, double ph1, double th2, double ph2)
{
	#define TO_RAD (3.1415926536 / 180)
	int R = 6371;
	double dx, dy, dz;
	ph1 -= ph2;
	ph1 *= TO_RAD, th1 *= TO_RAD, th2 *= TO_RAD;
	dz = sin(th1) - sin(th2);
	dx = cos(ph1) * cos(th1) - cos(th2);
	dy = sin(ph1) * cos(th1);
	return asin(sqrt(dx * dx + dy * dy + dz * dz) / 2) * 2 * R;
}

/*
 * tile_destroy
 * This function simply destroys any data associated with a tile
 */
void tile_destroy(tile_t* tile){
	if (tile->data != NULL)
		free(tile->data);
}

