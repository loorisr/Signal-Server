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

/* Per-ray data for the -fast mode.  The ITM profile of receiver point y is
 * the prefix elev[0..y+1] of the same array for the whole radial, so the
 * following can be built once per ray instead of once per pixel:
 *
 *   dtx[i]   : distance from TX of profile point i, accumulated exactly as
 *              hzns() does (float d += xi), so it is bit-identical;
 *   tx_c[i]  : (z_i - z_tx) / dtx[i], the TX elevation angle of point i
 *              without the earth-curvature term, and its prefix argmax:
 *              the TX horizon of any receiver point in O(1);
 *   hull     : upper convex hull of the curvature-corrected profile, from
 *              which the RX horizon is a tangent query in O(log n);
 *   p0, p1   : prefix sums of the profile, for O(1) zsys and z1sq1().
 *
 * The exact (default) mode does not use this context. */
struct ITM_ray_ctx {
	bool fast_horizons = false;  /* -fast: use hzns_fast() and prefix sums */
	float tx_height = 0.0f;      /* TX antenna height (m) */
	float eno_ns = 301.0f;       /* surface refractivity, for qc_ref */
	int last_i = 0;              /* largest valid point index i (pfl[2+i]) */
	double qc_ref = 0.0;         /* reference curvature term (gme/2) */

	std::vector<float> dtx;
	std::vector<float> tx_c;     /* (z_i - z_tx) / dtx[i] */
	std::vector<int> tx_best;    /* prefix argmax of tx_c - qc_ref*dtx */
	std::vector<double> hx, hg;  /* RX hull points: x_i, z_i - qc_ref*x_i^2 */
	std::vector<double> p0, p1;  /* profile prefix sums */

	/* Incremental upper convex hull of points 1..added. */
	struct RxHull {
		std::vector<int> v;
		int added = 0;
		void reset();
		void add_upto(int last, const double *x, const double *g);
		int query(double D, double W, const double *x, const double *g) const;
	};
	/* [0]: receiver points in ray order, [1]: -fast look-ahead points.
	 * Each sequence of profile lengths is monotone. */
	mutable RxHull hull[2];
	mutable int stream = 0;

	/* Build from elev[2..last_index] (elev[1] = sample spacing). */
	void build(const float elev[], int last_index);
};

void point_to_point_ITM_fast(const ITM_ctx &ctx, float rht_m, float elev[],
			     float &dbloss, PropagationMode &mode, int &errnum,
			     float *cached_dh, bool reuse_cached_dh,
			     const ITM_ray_ctx *ray = nullptr);

float ITM_delta_h_fast(const ITM_ctx &ctx, float rht_m, float elev[],
		       const ITM_ray_ctx *ray = nullptr);

#endif /* _ITWOM30_HH_ */
