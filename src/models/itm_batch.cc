/* -fast: ITM evaluation in batches of pixels, with lrprop() vectorised.
 *
 * Kept in its own translation unit on purpose: adding code to itwom3.0.cc
 * changes the compiler's inlining decisions there, and the exact mode must
 * stay bit-identical to the original code.
 *
 * Per pixel: prepare (scalar: profile analysis, horizons, delta-h,
 * effective heights), lrprop (vectorised over 8 pixels in structure-of-
 * arrays form for the common line-of-sight-region case, scalar otherwise),
 * finish (scalar: mode, avar, free-space term).  Results differ from
 * point_to_point_ITM_fast() by float rounding only (vector libm). */

#include <math.h>
#include <complex>
#include <algorithm>
#include <vector>
#include <cfloat>

#include "../common.hh"
#include "itwom3.0.hh"
#include "itm_types.hh"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define FORTRAN_DIM(x, y) (((x) - (y) > 0.0) ? (x) - (y) : 0.0)

using namespace std;

/* Functions of itwom3.0.cc */
void hzns(float pfl[], prop_type &prop);
void hzns_fast(float pfl[], prop_type &prop, const ITM_ray_ctx &ray);
void z1sq1(float z[], const float &x1, const float &x2, float &z0, float &zn);
float d1thx(float pfl[], const float &x1, const float &x2);
float calculate_delta_h_adjustable(float pfl[], const float &x1, const float &x2, int N);
void lrprop(float d, prop_type &prop, propa_type &propa);
float avar(float zzt, float zzl, float zzc, prop_type &prop, propv_type &propv);

static void z1sq1_prefixed(const float z[], const float &x1,
                           const float &x2, float &z0, float &zn,
                           const ITM_ray_ctx &ray)
{
    float xn, xa, xb, x, a, b;
    xn = z[0];
    xa = int(MAX(x1 / z[1], 0.0f));
    xb = xn - int(FORTRAN_DIM(xn, x2 / z[1]));
    if (xb <= xa) {
        xa = FORTRAN_DIM(xa, 1.0f);
        xb = xn - FORTRAN_DIM(xn, xb + 1.0f);
    }

    const int ja = (int)xa;
    const int jb = (int)xb;
    const int n  = jb - ja;
    xa = xb - xa;
    x  = -0.5f * xa;
    xb += x;

    const float p_first = z[ja + 2];
    const float p_last  = z[jb + 2];

    double interior0 = 0.0;
    double interior1 = 0.0;
    if (n > 0) {
        interior0 = ray.p0[jb + 1] - ray.p0[ja + 2];
        interior1 = (ray.p1[jb + 1] - ray.p1[ja + 2])
                  - (double)(ja + 2) * interior0;
    }

    const double a_sum = 0.5 * ((double)p_first + p_last) + interior0;
    /* b_sum initial term is 0.5*(p_first - p_last)*x; the loop then adds
     * p[i]*(x+i) for i=1..n-1, hence b_sum = x*(a_sum - p_last) + interior1. */
    const double b_sum = x * (a_sum - (double)p_last) + interior1;

    a = (float)a_sum / xa;
    b = (float)b_sum * 12.0f / ((xa * xa + 2.0f) * xa);
    z0 = a - b * xb;
    zn = a + b * (xn - xb);
}

static void qlrpfl_no_lrprop(float pfl[], int klimx, int mdvarx, prop_type & prop, propa_type & /*propa*/, propv_type & propv,
            float *cached_dh, bool reuse_cached_dh, const ITM_ray_ctx *ray)
{
	float xl[2], q, za, zb, temp;

	prop.dist = pfl[0] * pfl[1];
	const unsigned int np = (unsigned int)pfl[0];

	if (ray && ray->fast_horizons)
		hzns_fast(pfl, prop, *ray);
	else
		hzns(pfl, prop);

	xl[0] = MIN(15.0 * prop.hg[0], 0.1 * prop.dl[0]);
	xl[1] = prop.dist - MIN(15.0 * prop.hg[1], 0.1 * prop.dl[1]);

	if (reuse_cached_dh && cached_dh != NULL) {
		prop.dh = *cached_dh;
	} else {
		if (dh_n < 0)
			prop.dh = d1thx(pfl, xl[0], xl[1]);
		else
			prop.dh = calculate_delta_h_adjustable(pfl, xl[0], xl[1], dh_n);

		if (cached_dh != NULL)
			*cached_dh = prop.dh;
	}

	if (prop.dl[0] + prop.dl[1] > 1.5 * prop.dist) {
		if (ray && ray->fast_horizons)
			z1sq1_prefixed(pfl, xl[0], xl[1], za, zb, *ray);
		else
			z1sq1(pfl, xl[0], xl[1], za, zb);
		prop.he[0] = prop.hg[0] + FORTRAN_DIM(pfl[2], za);
		prop.he[1] = prop.hg[1] + FORTRAN_DIM(pfl[np + 2], zb);

		prop.dl[0] = sqrt(2.0 * prop.he[0] / prop.gme) * exp(-0.07 * sqrt(prop. dh / MAX(prop. he[0], 5.0)));
		prop.dl[1] = sqrt(2.0 * prop.he[1] / prop.gme) * exp(-0.07 * sqrt(prop. dh / MAX(prop. he[1], 5.0)));

		q = prop.dl[0] + prop.dl[1];

		if (q <= prop.dist) {	/* if there is a rounded horizon, or two obstructions, in the path */
			/* q=pow(prop.dist/q,2.0); */
			temp = prop.dist / q;
			q = temp * temp;

			prop.he[0] *= q;	/* tx effective height set to be path dist/distance between obstacles */
			prop.dl[0] = sqrt(2.0 * prop.he[0] / prop.gme) * exp(-0.07 * sqrt(prop.dh / MAX(prop.he[0], 5.0)));
			prop.he[1] *= q;	/* tx effective height set to be path dist/distance between obstacles */
			prop.dl[1] = sqrt(2.0 * prop.he[1] / prop.gme) * exp(-0.07 * sqrt(prop.dh / MAX(prop.he[1], 5.0)));
		}

		/* original empirical adjustment?  uses delta-h to adjust grazing angles */
		q = sqrt(2.0 * prop.he[0] / prop.gme);
		prop.the[0] = (0.65 * prop.dh * (q / prop.dl[0] - 1.0) - 2.0 * prop.he[0]) / q;
		q = sqrt(2.0 * prop.he[1] / prop.gme);
		prop.the[1] = (0.65 * prop.dh * (q / prop.dl[1] - 1.0) - 2.0 * prop.he[1]) / q;
	}

	else {
		if (ray && ray->fast_horizons) {
			z1sq1_prefixed(pfl, xl[0], 0.9 * prop.dl[0], za, q, *ray);
			z1sq1_prefixed(pfl, prop.dist - 0.9 * prop.dl[1], xl[1], q, zb, *ray);
		} else {
			z1sq1(pfl, xl[0], 0.9 * prop.dl[0], za, q);
			z1sq1(pfl, prop.dist - 0.9 * prop.dl[1], xl[1], q, zb);
		}
		prop.he[0] = prop.hg[0] + FORTRAN_DIM(pfl[2], za);
		prop.he[1] = prop.hg[1] + FORTRAN_DIM(pfl[np + 2], zb);
	}

	prop.mdp = -1;
	propv.lvar = MAX(propv.lvar, 3);

	if (mdvarx >= 0) {
		propv.mdvar = mdvarx;
		propv.lvar = MAX(propv.lvar, 4);
	}

	if (klimx > 0) {
		propv.klim = klimx;
		propv.lvar = 5;
	}

}

static void itm_prepare(const ITM_ctx &ctx, float rht_m, float elev[],
			float *cached_dh, bool reuse_cached_dh, const ITM_ray_ctx *ray,
			prop_type &prop, propa_type &propa, propv_type &propv)
{
	float zsys = 0;
	unsigned int ja, jb, i, np;

	prop.hg[0]  = ctx.tht_m;
	prop.hg[1]  = rht_m;
	propv.klim  = ctx.klim;
	prop.kwx    = 0;
	propv.lvar  = 5;
	prop.mdp    = -1;

	np = (unsigned int)elev[0];
	ja = (unsigned int)(3.0 + 0.1 * elev[0]);
	jb = np - ja + 6;

	if (ray && ray->fast_horizons) {
		/* -fast: prefix sums, sum(elev[ja-1 .. jb-1]) */
		zsys = (float)(ray->p0[jb - 1] - ray->p0[ja - 2]);
	} else {
		for (i = ja - 1; i < jb; ++i)
			zsys += elev[i];
	}
	zsys /= (jb - ja + 1);

	/* Inline the constant-optimised portion of qlrps():
	 *   wn and zgnd are precomputed; only ens and gme depend on zsys. */
	prop.wn       = ctx.wn;
	prop.ens      = ctx.eno;
	if (zsys != 0.0)
		prop.ens *= exp(-zsys / 9460.0);
	prop.gme      = 157e-9 * (1.0 - 0.04665 * exp(prop.ens / 179.3));
	prop.zgndreal = ctx.zgndreal;
	prop.zgndimag = ctx.zgndimag;

	propv.mdvar = ctx.mdvar;
	qlrpfl_no_lrprop(elev, propv.klim, propv.mdvar, prop, propa, propv,
	                 cached_dh, reuse_cached_dh, ray);
}

static void itm_finish(const ITM_ctx &ctx, prop_type &prop, propa_type &propa,
		       propv_type &propv, float &dbloss, PropagationMode &mode,
		       int &errnum)
{
	float q, fs;

	fs = ctx.fs_fmhz + 20.0f * log10(prop.dist) - 60.0f;
	q  = prop.dist - propa.dla;

	if (int(q) < 0.0)
		mode = PROP_MODE_LOS;
	else if (int(q) == 0.0)
		mode = PROP_MODE_1_HRZN;
	else if (int(q) > 0.0)
		mode = PROP_MODE_2_HRZN;

	if (mode != PROP_MODE_LOS) {
		if (prop.dist <= propa.dlsa || prop.dist <= propa.dx)
			mode |= PROP_MODE_DIFFRACTION;
		else if (prop.dist > propa.dx)
			mode |= PROP_MODE_TROPOSCATTER;
	}

	dbloss = avar(ctx.zr, 0.0, ctx.zc, prop, propv) + fs;
	errnum = prop.kwx;
}

/* ------------------------------------------------------------------------
 * lrprop() for the line-of-sight region (dist < dlsa), one pixel per SIMD
 * lane.  Straight transcription of lrprop/adiff/fht/aknfe/alos with
 * prop.mdp < 0 (initialisation for each pixel), the branches turned into
 * selects and the complex arithmetic written out.  Pixels beyond dlsa
 * (troposcatter region, rare) use the scalar lrprop().
 * ------------------------------------------------------------------------ */

namespace {

constexpr int W = 8;   /* pixels per batch (one AVX2 register of floats) */

static inline __attribute__((always_inline)) float v_aknfe(float v2)
{
	const float a = 6.02f + 9.11f * sqrtf(MAX(v2, 0.0f)) - 1.27f * v2;
	const float b = 12.953f + 10.0f * log10f(MAX(v2, 1e-30f));
	return (v2 < 5.76f) ? a : b;
}

static inline __attribute__((always_inline)) float v_fht(float x, float pk)
{
	const float lx = log10f(MAX(x, 1e-30f));
	const float w = -logf(pk);
	const float fa = -117.0f + ((x > 1.0f) ? 40.0f * lx : 0.0f);
	const float fb = 2.5e-5f * x * x / pk - 8.686f * w - 15.0f;
	const float small = (pk < 1.0e-5f || x * w * w * w > 5495.0f) ? fa : fb;
	float big = 0.05751f * x - 10.0f * lx;
	const float w2 = 0.0134f * x * expf(-0.005f * x);
	big = (x < 2000.0f) ? (1.0f - w2) * big + w2 * (40.0f * lx - 117.0f) : big;
	return (x < 200.0f) ? small : big;
}

static inline __attribute__((always_inline)) float v_h0f(float r, float et)
{
	/* a[], b[] tables of h0f() for it = 1..5, selected without lookups */
	const float fit = floorf(MAX(et, 0.0f));          /* (int)et for et >= 0 */
	const float itf = MIN(MAX(fit, 1.0f), 5.0f);
	const float q = (fit <= 0.0f || fit >= 5.0f) ? 0.0f : et - fit;

	const float a0 = (itf == 1.0f) ? 25.0f : (itf == 2.0f) ? 80.0f : (itf == 3.0f) ? 177.0f : (itf == 4.0f) ? 395.0f : 705.0f;
	const float b0 = (itf == 1.0f) ? 24.0f : (itf == 2.0f) ? 45.0f : (itf == 3.0f) ? 68.0f : (itf == 4.0f) ? 80.0f : 105.0f;
	const float a1 = (itf == 1.0f) ? 80.0f : (itf == 2.0f) ? 177.0f : (itf == 3.0f) ? 395.0f : 705.0f;
	const float b1 = (itf == 1.0f) ? 45.0f : (itf == 2.0f) ? 68.0f : (itf == 3.0f) ? 80.0f : 105.0f;
	const float t = 1.0f / r;
	const float x = t * t;
	const float h = 4.343f * logf((a0 * x + b0) * x + 1.0f);
	const float h1 = (1.0f - q) * h + q * 4.343f * logf((a1 * x + b1) * x + 1.0f);
	return (q != 0.0f) ? h1 : h;
}

static inline __attribute__((always_inline)) float v_ahd(float td)
{
	const float a = (td <= 10e3f) ? 133.4f : (td <= 70e3f) ? 104.6f : 71.8f;
	const float b = (td <= 10e3f) ? 0.332e-3f : (td <= 70e3f) ? 0.212e-3f : 0.157e-3f;
	const float c = (td <= 10e3f) ? -4.343f : (td <= 70e3f) ? -1.086f : 2.171f;
	return a + b * td + c * logf(MAX(td, 1e-30f));
}

/* Inputs of one batch of W pixels (structure of arrays). */
struct LaneIn {
	float hg0[W], hg1[W], he0[W], he1[W], dl0[W], dl1[W], th0[W], th1[W];
	float dh[W], gme[W], dist[W], ens[W], rch0[W];
};

struct LaneOut {
	float dls0[W], dls1[W], dlsa[W], dla[W], tha[W], emd[W], aed[W];
	float ak1[W], ak2[W], ael[W];          /* line-of-sight region */
	float ems[W], aes[W], dx[W];           /* beyond the horizon */
	float aref[W];
};

/* Part of lrprop() common to both regions: horizons and the diffraction
 * line (adiff at d3 and d4). */
static inline __attribute__((always_inline)) void v_common(float hg0, float hg1, float he0, float he1, float dl0, float dl1,
		     float th0, float th1, float gme, float dh, float wn, float qkz,
		     float &dls0, float &dls1, float &dlsa, float &dla, float &tha,
		     float &emd, float &aed, float &xae)
{
	dls0 = sqrtf(2.0f * he0 / gme);
	dls1 = sqrtf(2.0f * he1 / gme);
	dlsa = dls0 + dls1;
	dla = dl0 + dl1;
	tha = MAX(th0 + th1, -dla * gme);

	/* adiff(0): coefficients */
	float q = hg0 * hg1;
	const float qk = he0 * he1 - q;
	q += 10.0f;                                /* prop.mdp < 0 */
	const float wd1 = sqrtf(1.0f + qk / q);
	const float xd1 = dla + tha / gme;
	q = (1.0f - 0.8f * expf(-dlsa / 50e3f)) * dh;
	q *= 0.78f * expf(-powf(q / 16.0f, 0.25f));
	const float afo = MIN(15.0f, 2.171f * logf(1.0f + 4.77e-4f * hg0 * hg1 * wn * q));

	float aht = 20.0f, xht = 0.0f;
	{
		const float a = 0.5f * (dl0 * dl0) / he0;
		const float wa = cbrtf(a * wn);
		const float pk = qkz / wa;
		const float qq = (1.607f - pk) * 151.0f * wa * dl0 / a;
		xht += qq;
		aht += v_fht(qq, pk);
	}
	{
		const float a = 0.5f * (dl1 * dl1) / he1;
		const float wa = cbrtf(a * wn);
		const float pk = qkz / wa;
		const float qq = (1.607f - pk) * 151.0f * wa * dl1 / a;
		xht += qq;
		aht += v_fht(qq, pk);
	}

	xae = 1.0f / cbrtf(wn * gme * gme);
	const float d3 = MAX(dlsa, 1.3787f * xae + dla);
	const float d4 = d3 + 2.7574f * xae;

	/* adiff(d3), adiff(d4) */
	float ad[2];
	const float dd[2] = { d3, d4 };
	for (int k = 0; k < 2; ++k) {
		const float d = dd[k];
		const float th = tha + d * gme;
		const float ds = d - dla;
		float qa = 0.0795775f * wn * ds * th * th;
		const float v = v_aknfe(qa * dl0 / (ds + dl0)) + v_aknfe(qa * dl1 / (ds + dl1));
		const float a = ds / th;
		const float wa = cbrtf(a * wn);
		const float pk = qkz / wa;
		qa = (1.607f - pk) * 151.0f * wa * th + xht;
		const float ar = 0.05751f * qa - 4.343f * logf(qa) - aht;
		qa = (wd1 + xd1 / d) * MIN((1.0f - 0.8f * expf(-d / 50e3f)) * dh * wn, 6283.2f);
		const float wd = 25.1f / (25.1f + sqrtf(qa));
		ad[k] = ar * wd + (1.0f - wd) * v + afo;
	}
	emd = (ad[1] - ad[0]) / (d4 - d3);
	aed = ad[0] - emd * d3;
}

/* Line-of-sight region (dist < dlsa): alos() fit. */
void lrprop_los_lanes(const LaneIn &in, LaneOut &out, float wn, float zr, float zi)
{
	const float qkz = 1.0f / sqrtf(zr * zr + zi * zi);   /* 1/|zgnd| */

#pragma omp simd
	for (int l = 0; l < W; ++l) {
		const float he0 = in.he0[l], he1 = in.he1[l];
		const float dh = in.dh[l], dist = in.dist[l];
		float dls0, dls1, dlsa, dla, tha, emd, aed, xae;
		v_common(in.hg0[l], in.hg1[l], he0, he1, in.dl0[l], in.dl1[l], in.th0[l], in.th1[l],
		         in.gme[l], dh, wn, qkz, dls0, dls1, dlsa, dla, tha, emd, aed, xae);

		/* line-of-sight coefficients (lrprop, dist < dlsa) */
		const float wls = 0.021f / (0.021f + wn * dh / MAX(10e3f, dlsa));
		const float d2 = dlsa;
		const float a2 = aed + d2 * emd;
		float d0 = 1.908f * wn * he0 * he1;
		float d1;
		if (aed >= 0.0f) {
			d0 = MIN(d0, 0.5f * dla);
			d1 = d0 + 0.25f * (dla - d0);
		} else {
			d1 = MAX(-aed / emd, 0.25f * dla);
		}

		/* alos(d1), alos(d0) */
		float al[2];
		const float dx[2] = { d1, d0 };
		for (int k = 0; k < 2; ++k) {
			const float d = dx[k];
			float qq = (1.0f - 0.8f * expf(-d / 50e3f)) * dh;
			const float s = 0.78f * qq * expf(-powf(qq / 16.0f, 0.25f));
			qq = he0 + he1;
			const float sps = qq / sqrtf(d * d + qq * qq);
			/* r = (sps - z) / (sps + z) * exp(-min(10, wn*s*sps)) */
			const float nr = sps - zr, ni = -zi, dr = sps + zr, di = zi;
			const float den = dr * dr + di * di;
			const float e = expf(-MIN(10.0f, wn * s * sps));
			float rr = (nr * dr + ni * di) / den * e;
			float ri = (ni * dr - nr * di) / den * e;
			const float r2 = rr * rr + ri * ri;
			if (r2 < 0.25f || r2 < sps) {
				const float f = sqrtf(sps / r2);
				rr *= f;
				ri *= f;
			}
			const float v = emd * d + aed;
			qq = wn * he0 * he1 * 2.0f / d;
			if (qq > 1.57f)
				qq = 3.14f - 2.4649f / qq;
			const float cr = cosf(qq) + rr, ci = -sinf(qq) + ri;
			al[k] = (-4.343f * logf(cr * cr + ci * ci) - v) * wls + v;
		}
		const float a1 = al[0], a0 = al[1];

		float ak1, ak2;
		const float ak1_plain = (a2 - a1) / (d2 - d1);
		if (d0 < d1) {
			const float ql = logf(d2 / d0);
			const float ak2a = MAX(0.0f, ((d2 - d0) * (a1 - a0) - (d1 - d0) * (a2 - a0)) /
			                             ((d2 - d0) * logf(d1 / d0) - (d1 - d0) * ql));
			const bool wq = aed >= 0.0f || ak2a > 0.0f;
			if (wq) {
				const float ak1w = (a2 - a0 - ak2a * ql) / (d2 - d0);
				if (ak1w < 0.0f) {
					ak2 = MAX(a2 - a0, 0.0f) / ql;
					ak1 = (ak2 == 0.0f) ? emd : 0.0f;
				} else {
					ak1 = ak1w;
					ak2 = ak2a;
				}
			} else {
				ak2 = 0.0f;
				ak1 = (ak1_plain <= 0.0f) ? emd : ak1_plain;
			}
		} else {
			ak2 = 0.0f;
			ak1 = (ak1_plain <= 0.0f) ? emd : ak1_plain;
		}
		const float ael = a2 - ak1 * d2 - ak2 * logf(d2);
		float aref = (dist > 0.0f) ? ael + ak1 * dist + ak2 * logf(dist) : 0.0f;
		aref = MAX(aref, 0.0f);

		out.dls0[l] = dls0;
		out.dls1[l] = dls1;
		out.dlsa[l] = dlsa;
		out.dla[l] = dla;
		out.tha[l] = tha;
		out.emd[l] = emd;
		out.aed[l] = aed;
		out.ak1[l] = ak1;
		out.ak2[l] = ak2;
		out.ael[l] = ael;
		out.aref[l] = aref;
	}
}

/* One ascat() evaluation at distance d (prop.mdp < 0 sequence: ascat(0),
 * then ascat(d6), then ascat(d5); h0s carries the state between calls). */
static inline __attribute__((always_inline)) float v_ascat(float d, float &h0s, float ad, float rr, float etq,
		     float th0, float th1, float he0, float he1, float gme,
		     float tha, float ens, float wn)
{
	/* Branch-free form of ascat(): both cases are computed and selected. */
	const float th = th0 + th1 + d * gme;
	float r2 = 2.0f * wn * th;
	const float r1 = r2 * he0;
	r2 *= he1;

	float ss = (d - ad) / (d + ad);
	float q = rr / ss;
	ss = MAX(0.1f, ss);
	q = MIN(MAX(0.1f, q), 10.0f);
	const float z0 = (d - ad) * (d + ad) * th * 0.25f / d;
	float temp = MIN(1.7f, z0 / 8.0e3f);
	temp = temp * temp * temp * temp * temp * temp;
	const float et = (etq * expf(-temp) + 1.0f) * z0 / 1.7556e3f;
	const float ett = MAX(et, 1.0f);
	const float r1s = MAX(r1, 1e-30f), r2s = MAX(r2, 1e-30f);
	float h0n = (v_h0f(r1s, ett) + v_h0f(r2s, ett)) * 0.5f;
	h0n += MIN(h0n, (1.38f - logf(ett)) * logf(ss) * logf(q) * 0.49f);
	h0n = MAX(h0n, 0.0f);
	const float t2 = (1.0f + 1.4142f / r1s) * (1.0f + 1.4142f / r2s);
	const float h0e = et * h0n + (1.0f - et) * 4.343f *
	                  logf(MAX((t2 * t2) * (r1 + r2) / (r1 + r2 + 2.8284f), 1e-30f));
	h0n = (et < 1.0f) ? h0e : h0n;
	h0n = (h0n > 15.0f && h0s >= 0.0f) ? h0s : h0n;

	/* h0s > 15: reuse the previous h0; otherwise r1, r2 < 0.2 makes
	 * ascat() return 1001 without updating h0s. */
	const float h0 = (h0s > 15.0f) ? h0s : h0n;

	const float tb = tha + d * gme;
	const float v = v_ahd(tb * d) + 4.343f * logf(47.7f * wn * (tb * tb * tb * tb)) -
	                0.1f * (ens - 301.0f) * expf(-tb * d / 40e3f) + h0;
	const float ret = (h0s <= 15.0f && r1 < 0.2f && r2 < 0.2f) ? 1001.0f : v;
	h0s = (h0s <= 15.0f && r1 < 0.2f && r2 < 0.2f) ? h0s : h0;
	return ret;
}

/* Beyond the smooth-earth horizon (dist >= dlsa): diffraction and
 * troposcatter lines. */
void lrprop_scat_lanes(const LaneIn &in, LaneOut &out, float wn, float zr, float zi)
{
	const float qkz = 1.0f / sqrtf(zr * zr + zi * zi);

#pragma omp simd
	for (int l = 0; l < W; ++l) {
		const float he0 = in.he0[l], he1 = in.he1[l], gme = in.gme[l];
		const float th0 = in.th0[l], th1 = in.th1[l], ens = in.ens[l];
		const float dist = in.dist[l];
		float dls0, dls1, dlsa, dla, tha, emd, aed, xae;
		v_common(in.hg0[l], in.hg1[l], he0, he1, in.dl0[l], in.dl1[l], th0, th1,
		         gme, in.dh[l], wn, qkz, dls0, dls1, dlsa, dla, tha, emd, aed, xae);

		/* ascat(0) */
		const float ad0 = in.dl0[l] - in.dl1[l];
		const float rr0 = he1 / in.rch0[l];
		const float ad = (ad0 < 0.0f) ? -ad0 : ad0;
		const float rr = (ad0 < 0.0f) ? 1.0f / rr0 : rr0;
		const float etq = (5.67e-6f * ens - 2.32e-3f) * ens + 0.031f;
		float h0s = -15.0f;

		const float d5 = dla + 200e3f;
		const float d6 = d5 + 200e3f;
		const float a6 = v_ascat(d6, h0s, ad, rr, etq, th0, th1, he0, he1, gme, tha, ens, wn);
		const float a5 = v_ascat(d5, h0s, ad, rr, etq, th0, th1, he0, he1, gme, tha, ens, wn);

		const bool ok = (a5 < 1000.0f);
		const float ems_ok = (a6 - a5) / 200e3f;
		const float dx_ok = MAX(dlsa, MAX(dla + 0.3f * xae * logf(47.7f * wn),
		                                  (a5 - aed - ems_ok * d5) / (emd - ems_ok)));
		const float ems = ok ? ems_ok : emd;
		const float dx = ok ? dx_ok : 10.e6f;
		const float aes = ok ? (emd - ems_ok) * dx_ok + aed : aed;
		float aref = (dist > dx) ? aes + ems * dist : aed + emd * dist;
		aref = MAX(aref, 0.0f);

		out.dls0[l] = dls0;
		out.dls1[l] = dls1;
		out.dlsa[l] = dlsa;
		out.dla[l] = dla;
		out.tha[l] = tha;
		out.emd[l] = emd;
		out.aed[l] = aed;
		out.ems[l] = ems;
		out.aes[l] = aes;
		out.dx[l] = dx;
		out.aref[l] = aref;
	}
}

/* The validity flags (prop.kwx) set by lrprop(): cheap, kept scalar. */
void lrprop_kwx(prop_type &prop, const propa_type &propa)
{
	complex<float> prop_zgnd(prop.zgndreal, prop.zgndimag);

	if (prop.wn < 0.838 || prop.wn > 210.0)
		prop.kwx = MAX(prop.kwx, 1);
	for (int j = 0; j < 2; j++)
		if (prop.hg[j] < 1.0 || prop.hg[j] > 1000.0)
			prop.kwx = MAX(prop.kwx, 1);
	for (int j = 0; j < 2; j++)
		if (abs(prop.the[j]) > 200e-3
		    || prop.dl[j] < 0.1 * propa.dls[j]
		    || prop.dl[j] > 3.0 * propa.dls[j])
			prop.kwx = MAX(prop.kwx, 3);
	if (prop.ens < 250.0 || prop.ens > 400.0 || prop.gme < 75e-9
	    || prop.gme > 250e-9
	    || prop_zgnd.real() <= abs(prop_zgnd.imag())
	    || prop.wn < 0.419 || prop.wn > 420.0)
		prop.kwx = 4;
	for (int j = 0; j < 2; j++)
		if (prop.hg[j] < 0.5 || prop.hg[j] > 3000.0)
			prop.kwx = 4;

	const float dmin = abs(prop.he[0] - prop.he[1]) / 200e-3;
	if (prop.dist > 0.0) {
		if (prop.dist > 1000e3)
			prop.kwx = MAX(prop.kwx, 1);
		if (prop.dist < dmin)
			prop.kwx = MAX(prop.kwx, 3);
		if (prop.dist < 1e3 || prop.dist > 2000e3)
			prop.kwx = 4;
	}
}

} // namespace

/* ------------------------------------------------------------------------ */

struct ITM_batch::Impl {
	const ITM_ctx *ctx = nullptr;
	std::vector<prop_type> P;
	std::vector<propa_type> A;
	std::vector<propv_type> V;
	std::vector<float> loss;
	std::vector<PropagationMode> mode;
	std::vector<int> idx_los, idx_scat;
};

ITM_batch::ITM_batch() : impl(new Impl) {}
ITM_batch::~ITM_batch() { delete impl; }

void ITM_batch::clear()
{
	impl->P.clear();
	impl->A.clear();
	impl->V.clear();
	impl->loss.clear();
	impl->mode.clear();
}

int ITM_batch::size() const { return (int)impl->P.size(); }

void ITM_batch::add(const ITM_ctx &ctx, float rht_m, float elev[], float *cached_dh,
		    bool reuse_cached_dh, const ITM_ray_ctx *ray)
{
	if (debug) cnt_point_to_point_ITM++;
	impl->ctx = &ctx;
	impl->P.emplace_back();
	impl->A.emplace_back();
	impl->V.emplace_back();
	itm_prepare(ctx, rht_m, elev, cached_dh, reuse_cached_dh, ray,
	            impl->P.back(), impl->A.back(), impl->V.back());
}

void ITM_batch::run(bool vectorised)
{
	const int n = size();
	impl->loss.resize((size_t)n);
	impl->mode.resize((size_t)n);
	if (n == 0)
		return;

	if (!vectorised) {
		for (int i = 0; i < n; ++i)
			lrprop(0.0, impl->P[(size_t)i], impl->A[(size_t)i]);
	} else {
		/* Split the pixels by region (the smooth-earth horizon distance
		 * dlsa is cheap to compute) so that each batch of W lanes runs a
		 * single code path. */
		std::vector<int> &los = impl->idx_los, &scat = impl->idx_scat;
		los.clear();
		scat.clear();
		for (int i = 0; i < n; ++i) {
			const prop_type &p = impl->P[(size_t)i];
			const float dlsa = sqrtf(2.0f * p.he[0] / p.gme) + sqrtf(2.0f * p.he[1] / p.gme);
			(p.dist < dlsa ? los : scat).push_back(i);
		}

		const prop_type &p0 = impl->P[0];
		LaneIn in;
		LaneOut out;
		for (int pass = 0; pass < 2; ++pass) {
			const std::vector<int> &idx = pass == 0 ? los : scat;
			const int m_all = (int)idx.size();
			for (int c = 0; c < m_all; c += W) {
				const int m = MIN(W, m_all - c);
				for (int l = 0; l < W; ++l) {
					const prop_type &p = impl->P[(size_t)idx[(size_t)(c + (l < m ? l : 0))]];
					in.hg0[l] = p.hg[0]; in.hg1[l] = p.hg[1];
					in.he0[l] = p.he[0]; in.he1[l] = p.he[1];
					in.dl0[l] = p.dl[0]; in.dl1[l] = p.dl[1];
					in.th0[l] = p.the[0]; in.th1[l] = p.the[1];
					in.dh[l] = p.dh; in.gme[l] = p.gme; in.dist[l] = p.dist;
					in.ens[l] = p.ens; in.rch0[l] = p.rch[0];
				}
				if (pass == 0)
					lrprop_los_lanes(in, out, p0.wn, p0.zgndreal, p0.zgndimag);
				else
					lrprop_scat_lanes(in, out, p0.wn, p0.zgndreal, p0.zgndimag);

				for (int l = 0; l < m; ++l) {
					const int i = idx[(size_t)(c + l)];
					prop_type &p = impl->P[(size_t)i];
					propa_type &a = impl->A[(size_t)i];
					/* float rounding can put a pixel on the other side of
					 * dlsa than the split above: use the scalar path */
					if ((pass == 0) != (p.dist < out.dlsa[l])) {
						lrprop(0.0, p, a);
						continue;
					}
					a.dls[0] = out.dls0[l];
					a.dls[1] = out.dls1[l];
					a.dlsa = out.dlsa[l];
					a.dla = out.dla[l];
					a.tha = out.tha[l];
					a.emd = out.emd[l];
					a.aed = out.aed[l];
					if (pass == 0) {
						a.ak1 = out.ak1[l];
						a.ak2 = out.ak2[l];
						a.ael = out.ael[l];
					} else {
						a.ems = out.ems[l];
						a.aes = out.aes[l];
						a.dx = out.dx[l];
					}
					p.aref = out.aref[l];
					lrprop_kwx(p, a);
				}
			}
		}
	}

	for (int i = 0; i < n; ++i) {
		int errnum = 0;
		itm_finish(*impl->ctx, impl->P[(size_t)i], impl->A[(size_t)i], impl->V[(size_t)i],
		           impl->loss[(size_t)i], impl->mode[(size_t)i], errnum);
	}
}

float ITM_batch::loss(int i) const { return impl->loss[(size_t)i]; }
PropagationMode ITM_batch::mode(int i) const { return impl->mode[(size_t)i]; }
