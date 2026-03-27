#ifndef _OUTPUT_HH_
#define _OUTPUT_HH_

#include "models/los.hh"

void DoPathLoss(char *filename, unsigned char ngs, struct site *xmtr);

int DoSigStr(char *filename, unsigned char ngs, struct site *xmtr);

void DoRxdPwr(char *filename, unsigned char ngs, struct site *xmtr);

void DoLOS(char *filename, unsigned char ngs, struct site *xmtr);

void PathReport(struct site source, struct site destination, char *name, char graph_it, PropModel propmodel, int pmenv,
                double rxGain);

void SeriesData(struct site source, struct site destination, char *name, unsigned char fresnel_plot, unsigned char normalised);

#endif /* _OUTPUT_HH_ */
