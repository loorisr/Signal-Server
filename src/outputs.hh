#ifndef _OUTPUT_HH_
#define _OUTPUT_HH_

#include <stdint.h>

#include "models/los.hh"

void write_geotiff_rgba(const uint8_t *rgba, int img_width, int img_height, const char *filename);

void DoPathLoss(char *filename);

int DoSigStr(char *filename);

void DoRxdPwr(char *filename);

void DoLOS(char *filename);

void PathReport(struct site source, struct site destination, const char *name, char graph_it, PropModel propmodel, double rxGain);

void SeriesData(struct site source, struct site destination, const char *name, unsigned char fresnel_plot, unsigned char normalised);

#endif /* _OUTPUT_HH_ */
