#ifndef _OUTPUT_HH_
#define _OUTPUT_HH_

#include "models/los.hh"

void DoPathLoss(char *filename);

int DoSigStr(char *filename);

void DoRxdPwr(char *filename);

void DoLOS(char *filename);

void PathReport(struct site source, struct site destination, const char *name, char graph_it, PropModel propmodel, double rxGain);

void SeriesData(struct site source, struct site destination, const char *name, unsigned char fresnel_plot, unsigned char normalised);

#endif /* _OUTPUT_HH_ */
