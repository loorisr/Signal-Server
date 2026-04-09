#ifndef _INPUTS_HH_
#define _INPUTS_HH_

#include "common.hh"

int LoadCopernicus(int tile_lat, int tile_lon);
int LoadPAT(char *az_filename, char *el_filename);
int LoadTopoData(bbox region);

#endif /* _INPUTS_HH_ */
