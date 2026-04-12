#ifndef _ITWOM30_HH_
#define _ITWOM30_HH_

#include "../common.hh"

void point_to_point_ITM(double tht_m, double rht_m, double eps_dielect,
			double sgm_conductivity, double eno_ns_surfref,
			double frq_mhz, int radio_climate, int pol,
			double conf, double rel, double &dbloss, PropagationMode &mode,
			int &errnum);
void point_to_point(double tht_m, double rht_m, double eps_dielect,
		    double sgm_conductivity, double eno_ns_surfref,
		    double frq_mhz, int radio_climate, int pol, double conf,
		    double rel, double &dbloss, PropagationMode &mode, int &errnum);

/* Precomputed context for repeated ITM calls where only rht_m and elev[] vary.
 * Initialise once with ITM_ctx_init(), then call point_to_point_ITM_fast()
 * for every receiver point in an area coverage run. */
struct ITM_ctx {
	double tht_m;       /* transmitter AGL height (m) */
	double eno;         /* surface refractivity (N-units) */
	double wn;          /* wave number: frq_mhz / 47.7 */
	double zgndreal;    /* precomputed ground impedance (real) */
	double zgndimag;    /* precomputed ground impedance (imag) */
	double zc;          /* qerfi(conf) */
	double zr;          /* qerfi(rel) */
	double fs_fmhz;     /* 32.45 + 20*log10(frq_mhz) — constant part of free-space loss */
	int klim;           /* radio climate */
	int mdvar;          /* mode variance (12) */
};

ITM_ctx ITM_ctx_init(double tht_m, double eps_dielect, double sgm_conductivity,
		     double eno_ns_surfref, double frq_mhz, int radio_climate,
		     int pol, double conf, double rel);

void point_to_point_ITM_fast(const ITM_ctx &ctx, double rht_m, double elev[],
			     double &dbloss, PropagationMode &mode, int &errnum);

#endif /* _ITWOM30_HH_ */
