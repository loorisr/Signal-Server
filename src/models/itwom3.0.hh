#ifndef _ITWOM30_HH_
#define _ITWOM30_HH_

#include <vector>

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

/* Per-ray prefix sums of the terrain profile.  The profile of a radial is
 * a growing prefix of the same elev[] array, so values that qlrpfl currently
 * recomputes with O(n) scans can be answered in O(1):
 *
 *   p0[k] = sum_{i=2..k} elev[i]
 *   p1[k] = sum_{i=2..k} i * elev[i]
 *
 * Both arrays must be valid for indices up to the last profile point. */
struct ITM_ray_ctx {
	std::vector<double> p0;
	std::vector<double> p1;

	/* Build the prefix sums for elev[2..last_index]. */
	void build(const float elev[], int last_index);
};

void point_to_point_ITM_fast(const ITM_ctx &ctx, float rht_m, float elev[],
			     float &dbloss, PropagationMode &mode, int &errnum,
			     float *cached_dh, bool reuse_cached_dh,
			     const ITM_ray_ctx *ray = nullptr);

/* Horizon state produced by hzns().  Setting an override makes the next
 * hzns() call reuse these values instead of scanning the profile. */
struct ITM_horizon {
	float dl0, dl1;
	float the0, the1;
};

void ITM_horizon_override(const ITM_horizon *hz); /* nullptr disables */
ITM_horizon ITM_last_horizon();

#endif /* _ITWOM30_HH_ */
