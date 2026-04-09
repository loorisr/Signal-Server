#ifndef _CMDLINE_HH_
#define _CMDLINE_HH_

#include "models/los.hh"

/* Parse command-line arguments, initialise globals, validate inputs, and print
 * a summary.  Calls exit() on fatal errors (bad arguments, etc.). */
void parse_cmdline(int argc, char *argv[]);

#endif /* _CMDLINE_HH_ */
