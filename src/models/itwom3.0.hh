#ifndef _ITWOM30_HH_
#define _ITWOM30_HH_

#include "../common.hh"

void point_to_point_ITM(float tht_m, float rht_m, float eps_dielect,
			float sgm_conductivity, float eno_ns_surfref,
			float frq_mhz, int radio_climate, int pol,
			float conf, float rel, float &dbloss, PropagationMode &mode,
			int &errnum);
void point_to_point(float tht_m, float rht_m, float eps_dielect,
		    float sgm_conductivity, float eno_ns_surfref,
		    float frq_mhz, int radio_climate, int pol, float conf,
		    float rel, float &dbloss, PropagationMode &mode, int &errnum);

/* Precomputed context for repeated ITM calls where only rht_m and elev[] vary.
 * Initialise once with ITM_ctx_init(), then call point_to_point_ITM_fast()
 * for every receiver point in an area coverage run. */
struct ITM_ctx {
	float tht_m;       /* transmitter AGL height (m) */
	float eno;         /* surface refractivity (N-units) */
	float wn;          /* wave number: frq_mhz / 47.7 */
	float zgndreal;    /* precomputed ground impedance (real) */
	float zgndimag;    /* precomputed ground impedance (imag) */
	float zc;          /* qerfi(conf) */
	float zr;          /* qerfi(rel) */
	float fs_fmhz;     /* 32.45 + 20*log10(frq_mhz) — constant part of free-space loss */
	int klim;           /* radio climate */
	int mdvar;          /* mode variance (12) */
};

ITM_ctx ITM_ctx_init(float tht_m, float eps_dielect, float sgm_conductivity,
		     float eno_ns_surfref, float frq_mhz, int radio_climate,
		     int pol, float conf, float rel);

void point_to_point_ITM_fast(const ITM_ctx &ctx, float rht_m, float elev[],
			     float &dbloss, PropagationMode &mode, int &errnum);

#endif /* _ITWOM30_HH_ */
