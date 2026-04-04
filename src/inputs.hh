#ifndef _INPUTS_HH_
#define _INPUTS_HH_

#include "common.hh"
extern char scf_file[255];

int LoadSDF_SDF(char *name, int winfiles);
int LoadCopernicus(int tile_lat, int tile_lon);
int LoadPAT(char *az_filename, char *el_filename);
int LoadSignalColors(struct site xmtr);
int LoadLossColors(struct site xmtr);
int LoadDBMColors(struct site xmtr);
int LoadTopoData(bbox region);
static const char AZ_FILE_SUFFIX[] = ".az";
static const char EL_FILE_SUFFIX[] = ".el";

#endif /* _INPUTS_HH_ */
