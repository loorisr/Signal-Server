#ifndef _CMDLINE_HH_
#define _CMDLINE_HH_

#include "models/los.hh"

struct CmdlineArgs {
    int ppa;
    int normalise;
    int number_threads;
    double altitudeLR;
    PropModel prop_model;
    char mapfile[255];
};

/* Parse command-line arguments, initialise globals, validate inputs, and print
 * a summary.  Returns a CmdlineArgs struct with values that main() needs after
 * the parse phase.  Calls exit() on fatal errors (bad arguments, etc.). */
CmdlineArgs parse_cmdline(int argc, char *argv[]);

#endif /* _CMDLINE_HH_ */
