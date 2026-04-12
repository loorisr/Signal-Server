
/******************************************************************************
* ITWOM version 3.0a, January 20, 2011  File: itwom3.0a.cpp                  *
* Provenance:  Further test version of itwom2.0m re adj to Hrzn rangefactors  *
* 1. This file is based on a thorough debugging, completion, and update of    *
* the ITM, based on an original, public domain version of this file obtained  *
* from: ftp://flattop.its.bldrdoc.gov/itm/ITMDLL.cpp prior to May, 2007. C++  *
* routines for this program are taken from a translation of the FORTRAN code  *
* written by U.S. Department of Commerce NTIA/ITS Institute for               *
* Telecommunication Sciences Irregular Terrain Model (ITM) (Longley-Rice).    *
* 2. The Linux version of this file incorporates improvements suggested by a  *
* study of changes made to file itm.cpp by J. D. McDonald to remove Microsoft *
* Windows dll-isms and to debug an ambguity in overloaded calls.              *
* 3. The Linux version of this file also incorporates improvements suggested  *
* by a study of further modifications made to itm.cpp by John A. Magliacane   *
* to remove unused variables, unneeded #includes, and to replace pow()        *
* statements with explicit multiplications to improve execution speed and     *
* accuracy.                                                                   *
* 4. On August 19, 2007 this file was modified by Sid Shumate to include      *
* changes and updates included in version 7.0 of ITMDLL.cpp, which was        *
* released by the NTIA/ITS on June 26, 2007. With correction set SS1 and      *
* SS2: itm71.cpp.                                                             *
* 5. On Feb. 5, 2008 this file became v.1.0 of the ITWOM with the addition,   *
* by Sid Shumate, of multiple corrections, the replacement of subroutines     *
* lrprop and alos with lrprop2 and alos2, and the addition of subroutine      *
* saalos to incorporate Radiative Transfer Engine (RTE) computations in the   *
* line of sight range.                                                        *
* Update 8 Jun 2010 to modify alos to match 2010 series of IEEE-BTS           *
* newsletter articles                                                          *
* Update June 12, 2010 to z version to change test outputs                    *
* Offshoot start date June 23, 2010 to start itwom2.0 dual version for FCC.   *
* Update to 2.0b July 25 to correct if statement errors in adiff2 re two peak *
* calculations starting at line 525                                            *
* Development to 2.0c 8 Aug 2010 after modifying saalos and adiff for full    *
* addition of saalos treatment to post obstruction calculations and            *
* debugging.                                                                   *
* Modified to make 1st obs loss=5.8 only, no clutter loss considered          *
*                                                                              *
* Commented out unused variables and calculations to eliminate gcc warnings   *
*    (-Wunused-but-set-variable)  -- John A. Magliacane -- July 25, 2013      *
******************************************************************************/

#include <math.h>
#include <complex>
#include <assert.h>
#include <string.h>
#include <algorithm>
#include <math.h>
#include <immintrin.h>

#include "../common.hh"

#define THIRD (1.0/3.0)
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define FORTRAN_DIM(x, y) (((x) - (y) > 0.0) ? (x) - (y) : 0.0)

using namespace std;

/* Structure to store a complex number (real + imaginary parts). */
struct tcomplex {
	double tcreal;
	double tcimag;
};

/*
 * prop_type : Main parameters of the propagation model.
 *   aref        : Reference attenuation (dB)
 *   dist        : Total distance between transmitter and receiver (m)
 *   hg[2]       : Antenna heights above ground (m) [0]=TX, [1]=RX
 *   rch[2]      : Site heights above sea level (m)
 *   wn          : Wave number (normalized frequency)
 *   dh          : Terrain roughness (delta-h, m)
 *   dhd          : Delta-h for diffraction (optional)
 *   ens         : Surface atmospheric refractivity (N-units)
 *   encc        : Vegetation/clutter refractivity (N-units)
 *   cch         : Vegetation cover / clutter height (m)
 *   cd          : Clutter density
 *   gme         : Effective Earth curvature (m^-1)
 *   zgndreal    : Real part of normalized surface impedance
 *   zgndimag    : Imaginary part of normalized surface impedance
 *   he[2]       : Effective antenna heights (m) [0]=TX, [1]=RX
 *   dl[2]       : Distances to the horizon (m) [0]=TX, [1]=RX
 *   the[2]      : Grazing angles at the horizon (rad) [0]=TX, [1]=RX
 *   tiw         : Terrain profile interval width (m)
 *   ght         : TX antenna height above ground (m)
 *   ghr         : RX antenna height above ground (m)
 *   rph         : Reflection point height in profile (m)
 *   hht         : TX obstacle peak height (m)
 *   hhr         : RX obstacle peak height (m)
 *   tgh         : TX height above local ground (m)
 *   tsgh        : TX site height above global ground (m)
 *   thera       : Receiver approach slope angle (rad)
 *   thenr       : Angle of last profile segment on RX side (rad)
 *   rpl         : Reflection point index in profile
 *   kwx         : Error/warning code (0=OK)
 *   mdp         : Calculation mode (-1=point-to-point, 0=area in progress, 1=area init)
 *   ptx         : Transmitter polarization (0=H, 1=V, 2=circular)
 *   los         : Line-of-sight indicator (1=LOS, 0=NLOS)
 */
struct prop_type {
	double aref;
	double dist;
	double hg[2];
	double rch[2];
	double wn;
	double dh;
	double dhd;
	double ens;
	double encc;
	double cch;
	double cd;
	double gme;
	double zgndreal;
	double zgndimag;
	double he[2];
	double dl[2];
	double the[2];
	double tiw;
	double ght;
	double ghr;
	double rph;
	double hht;
	double hhr;
	double tgh;
	double tsgh;
	double thera;
	double thenr;
	int rpl;
	int kwx;
	int mdp;
	int ptx;
	int los;
};

/*
 * propv_type : Statistical variability parameters.
 *   sgc    : Composite standard deviation of the distribution (dB)
 *   lvar   : Required update level (0=none, ..., 5=all)
 *   mdvar  : Variability mode (0-3, ±10, ±20)
 *   klim   : Radio climate index (1=Equatorial ... 7=Maritime Temperate)
 */
struct propv_type {
	double sgc;
	int lvar;
	int mdvar;
	int klim;
};

/*
 * propa_type : Intermediate coefficients computed by lrprop/lrprop2.
 *   dlsa        : Radio line-of-sight distance (smooth earth, m)
 *   dx          : Diffraction/troposcatter transition distance (m)
 *   ael         : LOS model y-intercept (dB)
 *   ak1         : Linear slope of LOS model (dB/m)
 *   ak2         : Logarithmic coefficient of LOS model
 *   aed         : Diffraction attenuation at reference distance (dB)
 *   emd         : Slope of diffraction model (dB/m)
 *   aes         : Troposcatter attenuation at transition distance (dB)
 *   ems         : Slope of troposcatter model (dB/m)
 *   dls[2]      : Smooth-earth line-of-sight distances [0]=TX, [1]=RX (m)
 *   dla         : Sum of TX + RX horizon distances (m)
 *   tha         : Total TX+RX grazing angle (rad)
 */
struct propa_type {
	double dlsa;
	double dx;
	double ael;
	double ak1;
	double ak2;
	double aed;
	double emd;
	double aes;
	double ems;
	double dls[2];
	double dla;
	double tha;
};

/*
 * aknfe - Knife-edge diffraction attenuation
 *
 * Computes the attenuation in dB from the Fresnel parameter v² (v squared).
 * Uses two empirical formulas depending on the value of v²:
 *   - If v² < 5.76  : polynomial formula
 *   - Otherwise     : logarithmic (asymptotic) formula
 *
 * Parameter:
 *   fresnel_v_squared (former: v2) : square of the Fresnel-Kirchhoff parameter
 * Returns: attenuation in dB
 */
double aknfe(const double &fresnel_v_squared) /* former parameter name: v2 */
{
	double attenuation_dB; /* former: a */

	if (fresnel_v_squared < 5.76)
		attenuation_dB = 6.02 + 9.11 * sqrt(fresnel_v_squared) - 1.27 * fresnel_v_squared;
	else
		attenuation_dB = 12.953 + 10 * log10(fresnel_v_squared);

	return attenuation_dB;
}

/*
 * fht - Height-gain function for spherical-Earth diffraction
 *
 * Computes the height gain (dB) used in smooth-Earth diffraction calculations.
 * Two regimes based on the value of x (normalized distance):
 *   - x < 200.0  : near field, computation via surface-wave attenuation
 *   - x >= 200.0 : far field, asymptotic formula with interpolation
 *
 * Parameters:
 *   norm_distance (former: x)      : dimensionless normalized distance
 *   surface_impedance (former: pk) : normalized surface impedance (|zgnd|^-1)
 * Returns: height gain in dB
 */
double fht(const double &norm_distance, const double &surface_impedance) /* former: x, pk */
{
	double log_inv_impedance; /* former: w */
	double height_gain_dB;   /* former: fhtv */

	if (norm_distance < 200.0) {
		log_inv_impedance = -log(surface_impedance);

		if (surface_impedance < 1.0e-5 || norm_distance * log_inv_impedance * log_inv_impedance * log_inv_impedance > 5495.0) {
			height_gain_dB = -117.0;

			if (norm_distance > 1.0)
				height_gain_dB = 40.0 * log10(norm_distance) + height_gain_dB;
		} else
			height_gain_dB = 2.5e-5 * norm_distance * norm_distance / surface_impedance - 8.686 * log_inv_impedance - 15.0;
	}

	else {
		height_gain_dB = 0.05751 * norm_distance - 10.0 * log10(norm_distance);

		if (norm_distance < 2000.0) {
			log_inv_impedance = 0.0134 * norm_distance * exp(-0.005 * norm_distance);
			height_gain_dB = (1.0 - log_inv_impedance) * height_gain_dB + log_inv_impedance * (40.0 * log10(norm_distance) - 117.0);
		}
	}
	return height_gain_dB;
}

/*
 * h0f - Height correction factor for tropospheric scatter (H0 factor)
 *
 * Computes a corrective H0 factor (dB) depending on the path ratio r and
 * the terrain/troposphere parameter et (effective turbulence index).
 * Uses 5-entry empirical lookup tables for interpolation.
 *
 * Parameters:
 *   path_ratio (former: r)        : ratio of scattering distances
 *   turbulence_index (former: et) : effective atmospheric turbulence index
 * Returns: H0 factor in dB
 */
double h0f(double path_ratio, double turbulence_index) /* former: r, et */
{
	/* Empirical coefficients a[] and b[] for 5 turbulence levels */
	double a[5] = { 25.0, 80.0, 177.0, 395.0, 705.0 };
	double b[5] = { 24.0, 45.0, 68.0, 80.0, 105.0 };
	double interp_frac, inv_r_squared; /* former: q, x */
	double h0_value, temp;             /* former: h0fv, temp */
	int table_index;                   /* former: it */

	table_index = (int)turbulence_index;

	if (table_index <= 0) {
		table_index = 1;
		interp_frac = 0.0;
	}

	else if (table_index >= 5) {
		table_index = 5;
		interp_frac = 0.0;
	}

	else
		interp_frac = turbulence_index - table_index;

	/* x=pow(1.0/r,2.0); */
	temp = 1.0 / path_ratio;
	inv_r_squared = temp * temp;

	h0_value = 4.343 * log((a[table_index - 1] * inv_r_squared + b[table_index - 1]) * inv_r_squared + 1.0);

	if (interp_frac != 0.0)
		h0_value =
		    (1.0 - interp_frac) * h0_value + interp_frac * 4.343 * log((a[table_index] * inv_r_squared + b[table_index]) * inv_r_squared +
					       1.0);

	return h0_value;
}

/*
 * ahd - Tropospheric scatter attenuation as a function of distance (Above-horizon diffraction)
 *
 * Computes the basic attenuation (dB) for tropospheric scatter at long range.
 * Three distance regimes covered by distinct empirical coefficients:
 *   - td <= 10 km
 *   - 10 km < td <= 70 km
 *   - td > 70 km
 *
 * Parameter:
 *   scatter_distance (former: td) : scattering distance (m)
 * Returns: attenuation in dB
 */
double ahd(double scatter_distance) /* former: td */
{
	int regime_index;  /* former: i */
	double a[3] = { 133.4, 104.6, 71.8 };
	double b[3] = { 0.332e-3, 0.212e-3, 0.157e-3 };
	double c[3] = { -4.343, -1.086, 2.171 };

	if (scatter_distance <= 10e3)
		regime_index = 0;

	else if (scatter_distance <= 70e3)
		regime_index = 1;

	else
		regime_index = 2;

	return a[regime_index] + b[regime_index] * scatter_distance + c[regime_index] * log(scatter_distance);
}

/*
 * abq_alos - Squared magnitude of a complex number (Absolute value squared)
 *
 * Computes |r|² = re² + im² for a double-precision complex number.
 * Used in the reflection coefficient calculation in alos/alos2.
 *
 * Parameter:
 *   reflection_coeff (former: r) : complex number (reflection coefficient)
 * Returns: squared magnitude (dimensionless)
 */
double abq_alos(complex < double >reflection_coeff) /* former: r */
{
	return reflection_coeff.real() * reflection_coeff.real() + reflection_coeff.imag() * reflection_coeff.imag();
}

/*
 * saalos - Vegetation cover attenuation in line-of-sight (Surface/Above-Average Line-of-Sight)
 *
 * Computes the additional attenuation (dB) due to vegetation cover (clutter) on a
 * line-of-sight or near-LOS path (ITWOM RTE).
 * Takes into account:
 *   - the geometry of penetration into the canopy (refraction at the air/canopy interface)
 *   - transmitter polarization (prop.ptx)
 *   - TX height relative to the canopy
 *
 * Parameters:
 *   path_distance (former: d)  : path distance (m)
 *   prop                       : propagation parameters
 *   propa                      : intermediate coefficients (unused here)
 * Returns: vegetation/clutter attenuation in dB
 */
double saalos(double path_distance, prop_type & prop, propa_type & /*propa*/) /* former: d */
{
	/* Geometric variables for refraction at the air/canopy interface */
	double surface_refractivity;         /* former: ensa  — surface refractivity (1 + N*1e-6) */
	double canopy_refractivity;          /* former: encca — canopy refractivity */
	double interp_q;                     /* former: q     — general intermediate variable */
	double path_dist_iter;               /* former: dp    — current path distance (iteration) */
	double horiz_dist_to_canopy;         /* former: dx    — horizontal distance to canopy */
	double earth_curve_height;           /* former: hc    — elevation due to Earth curvature */
	double tx_effective_height;          /* former: hone  — effective TX height relative to RX site */
	double incidence_angle;              /* former: tic   — total incidence angle at interface */
	double sin_incidence;                /* former: stic  */
	double cos_incidence;                /* former: ctic  */
	double refracted_angle;             /* former: ttc   — refracted angle inside canopy */
	double cos_refracted;               /* former: cttc  */
	double transmission_coeff_sq;       /* former: rsp   — power transmission coefficient */
	double reflection_coeff_sq;         /* former: tsp   — power reflection coefficient */
	double slant_path_in_canopy;        /* former: crpc  — slant path length inside canopy */
	double entry_angle_complement;      /* former: ssnps — complement of entry angle */
	double entry_dist_to_canopy;        /* former: d1a   — horizontal distance to canopy entry point */
	double reflection_coeff_temp;       /* former: rsp (reused locally as Fresnel q) */
	double clutter_attenuation;         /* former: arte  — computed attenuation (dB) */
	double freq_factor;                 /* former: zi    — normalized frequency distance factor */
	double path_dist_km;                /* former: pdk   — distance in km */
	double tip_angle, ucrpc_length;     /* former: tip, ucrpc */
	double tx_rx_height_diff;           /* former: tvsr  — TX-RX height difference */
	double saalos_result = 0.0;         /* former: saalosv */

	interp_q = 0.0;

	if (path_distance == 0.0) {
		transmission_coeff_sq = 1.0;
		reflection_coeff_sq = 0.0;
		entry_dist_to_canopy = 50.0;
		saalos_result = 0.0;
	} else if (prop.hg[1] > prop.cch) {
		/* RX antenna above canopy: no vegetation attenuation */
		saalos_result = 0.0;
	} else {
		path_dist_iter = path_distance;
		path_dist_km = path_dist_iter / 1000.0;
		transmission_coeff_sq = 1.0;
		reflection_coeff_sq = 0.0;
		entry_dist_to_canopy = path_dist_iter;

		/* hone = TX height relative to ground level at RX site */
		tx_effective_height = prop.tgh + prop.tsgh - (prop.rch[1] - prop.hg[1]);

		if (prop.tgh > prop.cch) { /* TX above canopy */
			surface_refractivity = 1 + prop.ens * 0.000001;
			canopy_refractivity  = 1 + prop.encc * 0.000001;
			path_dist_iter = path_distance;

			/* Iterations to find the canopy entry point */
			for (int j = 0; j < 5; ++j) {
				double tde = path_dist_iter / 6378137.0;           /* angle sous-tendu (rad) */
				earth_curve_height = (prop.cch + 6378137.0) * (1 - cos(tde));
				horiz_dist_to_canopy = (prop.cch + 6378137.0) * sin(tde);
				ucrpc_length =
				    sqrt((tx_effective_height - prop.cch + earth_curve_height) * (tx_effective_height -
							prop.cch + earth_curve_height) +
					 (horiz_dist_to_canopy * horiz_dist_to_canopy));
				double ctip = (tx_effective_height - prop.cch + earth_curve_height) / ucrpc_length;
				tip_angle = acos(ctip);
				incidence_angle = tip_angle + tde;
				incidence_angle = MAX(0.0, incidence_angle);
				sin_incidence = sin(incidence_angle);
				double sta = (surface_refractivity / canopy_refractivity) * sin_incidence;
				refracted_angle = asin(sta);
				cos_refracted = sqrt(1 - (sin(refracted_angle)) * (sin(refracted_angle)));
				slant_path_in_canopy = (prop.cch - prop.hg[1]) / cos_refracted;
				if (slant_path_in_canopy >= path_dist_iter) {
					slant_path_in_canopy = path_dist_iter - 1 / path_dist_iter;
				}

				entry_angle_complement = (3.1415926535897 / 2) - incidence_angle;
				entry_dist_to_canopy = (slant_path_in_canopy * sin(refracted_angle)) / (1 - 1 / 6378137.0);
				path_dist_iter = path_distance - entry_dist_to_canopy;
			}

			cos_incidence = cos(incidence_angle);

			/*
			 * If the ucrpc path touches the canopy before the end,
			 * the entry point is shifted toward the transmitter.
			 */
			if (entry_angle_complement <= 0.0) {
				entry_dist_to_canopy = MIN(0.1 * path_distance, 600.0);
				slant_path_in_canopy = entry_dist_to_canopy;
				/* TX height redefined slightly above canopy */
				tx_effective_height = prop.cch + 1;
				reflection_coeff_sq = .997;
				transmission_coeff_sq = 1 - reflection_coeff_sq;
			} else {
				/* Compute reflection coefficient based on polarization */
				if (prop.ptx >= 1) { /* vertical or circular polarization */
					interp_q = ((surface_refractivity * cos_refracted -
					      canopy_refractivity * cos_incidence) / (surface_refractivity * cos_refracted +
								       canopy_refractivity * cos_incidence));
					reflection_coeff_sq = interp_q * interp_q;
					transmission_coeff_sq = 1 - reflection_coeff_sq;

					if (prop.ptx == 2) { /* circular polarization */
						interp_q = ((surface_refractivity * cos_incidence -
						      canopy_refractivity * cos_refracted) / (surface_refractivity *
								       cos_incidence +
								       canopy_refractivity *
								       cos_refracted));
						reflection_coeff_sq =
						    ((surface_refractivity * cos_refracted -
						      canopy_refractivity * cos_incidence) / (surface_refractivity *
								       cos_refracted +
								       canopy_refractivity *
								       cos_incidence));
						reflection_coeff_sq = (interp_q * interp_q + reflection_coeff_sq * reflection_coeff_sq) / 2;
						transmission_coeff_sq = 1 - reflection_coeff_sq;
					}
				} else { /* horizontal polarization (ptx == 0) */
					interp_q = ((surface_refractivity * cos_incidence -
					      canopy_refractivity * cos_refracted) / (surface_refractivity * cos_incidence +
								       canopy_refractivity * cos_refracted));
					reflection_coeff_sq = interp_q * interp_q;
					transmission_coeff_sq = 1 - reflection_coeff_sq;
				}
			}

			/* tvsr = TX height difference above RX antenna */
			tx_rx_height_diff = MAX(0.0, prop.tgh + prop.tsgh - prop.rch[1]);

			if (entry_dist_to_canopy < 50.0) {
				clutter_attenuation = 0.0195 * slant_path_in_canopy - 20 * log10(transmission_coeff_sq);
			}

			else {
				if (entry_dist_to_canopy < 225.0) {
					if (tx_rx_height_diff > 1000.0) {
						interp_q = entry_dist_to_canopy * (0.03 * exp(-0.14 * path_dist_km));
					} else {
						interp_q = entry_dist_to_canopy * (0.07 * exp(-0.17 * path_dist_km));
					}

					clutter_attenuation =
					    interp_q + (0.7 * path_dist_km -
						 MAX(0.01,
						       log10(prop.wn * 47.7) -
						       2)) * (prop.hg[1] / tx_effective_height);
				}

				else {
					interp_q = 0.00055 * (path_dist_km) +
					    log10(path_dist_km) * (0.041 -
							  0.0017 * sqrt(tx_effective_height) +
							  0.019);

					clutter_attenuation =
					    entry_dist_to_canopy * interp_q -
					    (18 * log10(reflection_coeff_sq)) /
					    (exp(tx_effective_height / 37.5));

					freq_factor = 1.5 * sqrt(tx_effective_height - prop.cch);

					if (path_dist_km > freq_factor) {
						interp_q = (path_dist_km -
						     freq_factor) * 10.2 *
						    ((sqrt
						      (MAX
						       (0.01,
							log10(prop.wn * 47.7) -
							2.0))) / (100 - freq_factor));
					} else {
						interp_q = ((freq_factor -
						      path_dist_km) / freq_factor) * (-20.0 *
								    MAX(0.01,
									  log10
									  (prop.
									   wn *
									   47.7)
									  -
									  2.0))
						    / sqrt(tx_effective_height);
					}
					clutter_attenuation = clutter_attenuation + interp_q;
				}
			}
		} else { /* TX below or at canopy level */
			interp_q = (prop.cch - prop.tgh) * (2.06943 -
					     1.56184 * exp(1 /
							   prop.cch -
							   prop.tgh));
			interp_q = interp_q + (17.98 -
				 0.84224 * (prop.cch -
					    prop.tgh)) * exp(-0.00000061 * path_distance);
			clutter_attenuation = interp_q + 1.34795 * 20 * log10(path_distance + 1.0);
			clutter_attenuation =
			    clutter_attenuation -
			    (MAX(0.01, log10(prop.wn * 47.7) - 2)) *
			    (prop.hg[1] / prop.tgh);
		}
		saalos_result = clutter_attenuation;
	}
	return saalos_result;
}

/*
 * adiff - Diffraction attenuation (original ITM model)
 *
 * Computes the diffraction attenuation (dB) beyond the horizon over smooth Earth.
 * Must be called first with d==0.0 to initialize the static coefficients
 * (wd1, xd1, afo, qk, aht, xht), then with d>0 for each distance.
 *
 * Parameters:
 *   path_dist (former: d) : path distance (m), 0 for initialization
 *   prop                  : propagation parameters
 *   propa                 : intermediate coefficients
 * Returns: diffraction attenuation in dB
 */
double adiff(double path_dist, prop_type & prop, propa_type & propa) /* former: d */
{
	complex < double >prop_zgnd(prop.zgndreal, prop.zgndimag);
	/* Thread-local static coefficients, initialized at d==0 */
	static __thread double wd1, xd1, afo, qk, aht, xht;
	double slant_dist;          /* former: a  — normalized slant distance */
	double interp_q;            /* former: q  */
	double surface_impedance;   /* former: pk — local surface impedance */
	double beyond_horizon_dist; /* former: ds — distance beyond the horizon */
	double total_horizon_angle; /* former: th — total grazing angle */
	double aperture_param;      /* former: wa — diffraction aperture parameter */
	double round_earth_loss;    /* former: ar — rounded-Earth attenuation */
	double knife_weight;        /* former: wd — knife-edge/rounded weight */
	double diffraction_loss;    /* former: adiffv */

	if (path_dist == 0) {
		/* --- Coefficient initialization phase --- */
		interp_q = prop.hg[0] * prop.hg[1];
		qk = prop.he[0] * prop.he[1] - interp_q;

		if (prop.mdp < 0.0)
			interp_q += 10.0;

		wd1 = sqrt(1.0 + qk / interp_q);
		xd1 = propa.dla + propa.tha / prop.gme;
		interp_q = (1.0 - 0.8 * exp(-propa.dlsa / 50e3)) * prop.dh;
		interp_q *= 0.78 * exp(-pow(interp_q / 16.0, 0.25));
		afo =
		    MIN(15.0,
			  2.171 * log(1.0 +
				      4.77e-4 * prop.hg[0] * prop.hg[1] *
				      prop.wn * interp_q));
		qk = 1.0 / abs(prop_zgnd);
		aht = 20.0;
		xht = 0.0;

		/* Compute height gain for TX and RX */
		for (int j = 0; j < 2; ++j) {
			/* a=0.5*pow(prop.dl[j],2.0)/prop.he[j]; */
			slant_dist = 0.5 * (prop.dl[j] * prop.dl[j]) / prop.he[j];
			aperture_param = pow(slant_dist * prop.wn, THIRD);
			surface_impedance = qk / aperture_param;
			interp_q = (1.607 - surface_impedance) * 151.0 * aperture_param * prop.dl[j] / slant_dist;
			xht += interp_q;
			aht += fht(interp_q, surface_impedance);
		}

		diffraction_loss = 0.0;
	}

	else {
		/* --- Compute attenuation for distance path_dist --- */
		total_horizon_angle = propa.tha + path_dist * prop.gme;
		beyond_horizon_dist = path_dist - propa.dla;
		/* q=0.0795775*prop.wn*ds*pow(th,2.0); */
		interp_q = 0.0795775 * prop.wn * beyond_horizon_dist * total_horizon_angle * total_horizon_angle;
		diffraction_loss =
		    aknfe(interp_q * prop.dl[0] / (beyond_horizon_dist + prop.dl[0])) +
		    aknfe(interp_q * prop.dl[1] / (beyond_horizon_dist + prop.dl[1]));
		slant_dist = beyond_horizon_dist / total_horizon_angle;
		aperture_param = pow(slant_dist * prop.wn, THIRD);
		surface_impedance = qk / aperture_param;
		interp_q = (1.607 - surface_impedance) * 151.0 * aperture_param * total_horizon_angle + xht;
		round_earth_loss = 0.05751 * interp_q - 4.343 * log(interp_q) - aht;
		interp_q = (wd1 +
		     xd1 / path_dist) * MIN(((1.0 - 0.8 * exp(-path_dist / 50e3)) * prop.dh *
				       prop.wn), 6283.2);
		knife_weight = 25.1 / (25.1 + sqrt(interp_q));
		diffraction_loss = round_earth_loss * knife_weight + (1.0 - knife_weight) * diffraction_loss + afo;
	}

	return diffraction_loss;
}

/*
 * adiff2 - Extended diffraction attenuation (ITWOM model)
 *
 * Improved version of adiff() for the ITWOM model.
 * Handles 1 or 2 obstacles, vegetation (via saalos), and various
 * receiver grazing angle cases.
 * Must be called first with d==0.0 to initialize the static coefficients.
 *
 * Parameters:
 *   path_dist (former: d) : path distance (m), 0 for initialization
 *   prop                  : propagation parameters
 *   propa                 : intermediate coefficients
 * Returns: diffraction attenuation in dB
 */
double adiff2(double path_dist, prop_type & prop, propa_type & propa) /* former: d */
{
	complex < double >prop_zgnd(prop.zgndreal, prop.zgndimag);

	/* Thread-local static coefficients, initialized at path_dist==0 */
	static __thread double wd1, xd1, qk, aht, xht;
	/* Slant distances between TX, obstacles and RX */
	static __thread double
	    dist_tx_obs_top,    /* former: toh  — TX obstacle height above LOS */
	    dist_tx_obs_top2,   /* former: toho — TX obstacle height (variant) */
	    dist_rx_obs_top,    /* former: roh  — RX obstacle height above LOS */
	    dist_rx_obs_top2,   /* former: roho — RX obstacle height (variant) */
	    slant_tx_obs,       /* former: dto  */
	    slant_tx_obs2,      /* former: dto1 */
	    slant_txrx_obs,     /* former: dtro */
	    slant_rxtx_obs,     /* former: drto */
	    slant_rx_obs,       /* former: dro  */
	    slant_rx_obs2,      /* former: dro2 */
	    slant_total,        /* former: dtr  */
	    slant_hilltop1,     /* former: dhh1 */
	    slant_hilltop2,     /* former: dhh2 */
	    /* dhec, */
	    slant_tx_obs_base,  /* former: dtof  — path to base of TX obstacle */
	    slant_tx_obs2_base, /* former: dto1f */
	    slant_rx_obs_base,  /* former: drof  — path to base of RX obstacle */
	    slant_rx_obs2_base; /* former: dro2f */

	double slant_dist;          /* former: a  */
	double interp_q;            /* former: q  */
	double surface_impedance;   /* former: pk */
	double rx_dist_beyond;      /* former: rd — RX distance beyond obstacle */
	double beyond_horizon_dist; /* former: ds */
	double beyond_horizon_dist_limited; /* former: dsl */
	/* dfdh, */
	double total_horizon_angle; /* former: th */
	double aperture_param;      /* former: wa */
	/* ar, wd, sf1, */
	double foliage_scatter2;    /* former: sf2 — foliage scatter factor 2 obstacles */
	/* ec, */
	double fresnel_param;       /* former: vv */
	double ke_phase_diff;       /* former: kedr */
	double ke_amplitude;        /* former: arp */
	double scatter_delay;       /* former: sdr */
	double phase_diff_rad;      /* former: pd  */
	double scatter_amplitude;   /* former: srp */
	double ke_voltage;          /* former: kem */
	double combined_scatter_dB; /* former: csd */
	double scatter_voltage;     /* former: sdl */
	double diffraction_loss2 = 0.0; /* former: adiffv2 */
	double clutter_loss = 0.0;      /* former: closs */

	/* sf1=1.0; */ /* average foliage scatter factor 1 obstacle */
	foliage_scatter2 = 1.0; /* average foliage scatter factor 2 obstacles */

	/* dfdh=prop.dh; */
	/* ec=0.5*prop.gme; */

	/* adiff2 must be called first with path_dist==0.0 to initialize coefficients */
	if (path_dist == 0) {
		interp_q = prop.hg[0] * prop.hg[1];
		qk = prop.he[0] * prop.he[1] - interp_q;
		/* dhec=2.73; */

		if (prop.mdp < 0.0)
			interp_q += 10.0;

		/* Coefficients for the 4-ray model over rounded Earth */
		wd1 = sqrt(1.0 + qk / interp_q);
		xd1 = propa.dla + propa.tha / prop.gme;
		interp_q = (1.0 - 0.8 * exp(-propa.dlsa / 50e3)) * prop.dh;
		interp_q *= 0.78 * exp(-pow(interp_q / 16.0, 0.25));
		qk = 1.0 / abs(prop_zgnd);
		aht = 20.0;
		xht = 0.0;
		slant_dist = 0.5 * (prop.dl[0] * prop.dl[0]) / prop.he[0];
		aperture_param = pow(slant_dist * prop.wn, THIRD);
		surface_impedance = qk / aperture_param;
		interp_q = (1.607 - surface_impedance) * 151.0 * aperture_param * prop.dl[0] / slant_dist;
		xht = interp_q;
		aht += fht(interp_q, surface_impedance);

		if ((int (prop.dl[1]) == 0.0)||(prop.the[1] > 0.2)) {
			xht += xht;
			aht += (aht - 20.0);
		}

		else {
			slant_dist = 0.5 * (prop.dl[1] * prop.dl[1]) / prop.he[1];
			aperture_param = pow(slant_dist * prop.wn, THIRD);
			surface_impedance = qk / aperture_param;
			interp_q = (1.607 - surface_impedance) * 151.0 * aperture_param * prop.dl[1] / slant_dist;
			xht += interp_q;
			aht += fht(interp_q, surface_impedance);
		}
		diffraction_loss2 = 0.0;
	}

	else {
		total_horizon_angle = propa.tha + path_dist * prop.gme;

		beyond_horizon_dist_limited = MAX(path_dist - propa.dla, 0.0);
		beyond_horizon_dist = path_dist - propa.dla;
		slant_dist = beyond_horizon_dist / total_horizon_angle;
		aperture_param = pow(slant_dist * prop.wn, THIRD);
		surface_impedance = qk / aperture_param;

		/* Compute obstacle heights relative to the LOS */
		dist_tx_obs_top =
		    prop.hht - (prop.rch[0] -
				prop.dl[0] * ((prop.rch[1] - prop.rch[0]) /
					      prop.dist));
		dist_rx_obs_top =
		    prop.hhr - (prop.rch[0] -
				(prop.dist -
				 prop.dl[1]) * ((prop.rch[1] -
						 prop.rch[0]) / prop.dist));
		dist_tx_obs_top2 =
		    prop.hht - (prop.rch[0] -
				(prop.dl[0] +
				 beyond_horizon_dist_limited) * ((prop.hhr - prop.rch[0]) / (prop.dist -
							     prop.dl[1])));
		dist_rx_obs_top2 =
		    prop.hhr - (prop.hht -
				beyond_horizon_dist_limited * ((prop.rch[1] - prop.hht) / beyond_horizon_dist_limited));

		/* Slant distances between the various points */
		slant_tx_obs = sqrt(prop.dl[0] * prop.dl[0] + dist_tx_obs_top * dist_tx_obs_top);
		slant_tx_obs += prop.gme * prop.dl[0];
		slant_tx_obs2 = sqrt(prop.dl[0] * prop.dl[0] + dist_tx_obs_top2 * dist_tx_obs_top2);
		slant_tx_obs2 += prop.gme * prop.dl[0];
		slant_txrx_obs =
		    sqrt((prop.dl[0] + beyond_horizon_dist_limited) * (prop.dl[0] + beyond_horizon_dist_limited) +
			 prop.hhr * prop.hhr);
		slant_txrx_obs += prop.gme * (prop.dl[0] + beyond_horizon_dist_limited);
		slant_rxtx_obs =
		    sqrt((prop.dl[1] + beyond_horizon_dist_limited) * (prop.dl[1] + beyond_horizon_dist_limited) +
			 prop.hht * prop.hht);
		slant_rxtx_obs += prop.gme * (prop.dl[1] + beyond_horizon_dist_limited);
		slant_rx_obs = sqrt(prop.dl[1] * prop.dl[1] + dist_rx_obs_top * dist_rx_obs_top);
		slant_rx_obs += prop.gme * (prop.dl[1]);
		slant_rx_obs2 = sqrt(prop.dl[1] * prop.dl[1] + dist_rx_obs_top2 * dist_rx_obs_top2);
		slant_rx_obs2 += prop.gme * (prop.dl[1]);
		slant_total =
		    sqrt(prop.dist * prop.dist +
			 (prop.rch[0] - prop.rch[1]) * (prop.rch[0] -
							prop.rch[1]));
		slant_total += prop.gme * prop.dist;
		slant_hilltop1 =
		    sqrt((prop.dist - propa.dla) * (prop.dist - propa.dla) +
			 dist_tx_obs_top2 * dist_tx_obs_top2);
		slant_hilltop1 += prop.gme * (prop.dist - propa.dla);
		slant_hilltop2 =
		    sqrt((prop.dist - propa.dla) * (prop.dist - propa.dla) +
			 dist_rx_obs_top2 * dist_rx_obs_top2);
		slant_hilltop2 += prop.gme * (prop.dist - propa.dla);

		/* Distances to the obstacle base (below vegetation cover) */
		slant_tx_obs_base =
		    sqrt(prop.dl[0] * prop.dl[0] +
			 (dist_tx_obs_top - prop.cch) * (dist_tx_obs_top - prop.cch));
		slant_tx_obs_base += prop.gme * prop.dl[0];
		slant_tx_obs2_base =
		    sqrt(prop.dl[0] * prop.dl[0] +
			 (dist_tx_obs_top2 - prop.cch) * (dist_tx_obs_top2 - prop.cch));
		slant_tx_obs2_base += prop.gme * prop.dl[0];
		slant_rx_obs_base =
		    sqrt(prop.dl[1] * prop.dl[1] +
			 (dist_rx_obs_top - prop.cch) * (dist_rx_obs_top - prop.cch));
		slant_rx_obs_base += prop.gme * (prop.dl[1]);
		slant_rx_obs2_base =
		    sqrt(prop.dl[1] * prop.dl[1] +
			 (dist_rx_obs_top2 - prop.cch) * (dist_rx_obs_top2 - prop.cch));
		slant_rx_obs2_base += prop.gme * (prop.dl[1]);

		/* Prepare saalos parameters for post-obstacle reception */
		prop.tgh = prop.cch + 1.0;
		prop.tsgh = prop.hhr;
		rx_dist_beyond = prop.dl[1];

		/* Diffraction calculation for 2 obstacles */
		if (int (beyond_horizon_dist) > 0) { /* 2 obstacles */
			if (int (prop.dl[1]) > 0.0) { /* RX beyond the 2nd peak */
				/* Rounded-Earth attenuation */
				interp_q = (1.607 - surface_impedance) * 151.0 * aperture_param * total_horizon_angle + xht;
				/* ar=0.05751*q-10*log10(q)-aht; */

				/* Knife-edge/rounded weighting */
				interp_q = (1.0 - 0.8 * exp(-path_dist / 50e3)) * prop.dh;
				interp_q = (wd1 + xd1 / path_dist) * MIN((interp_q * prop.wn), 6283.2);
				/* wd=25.1/(25.1+sqrt(q)); */

				interp_q = 0.6365 * prop.wn;

				if (prop.the[1] < 0.2) { /* RX grazing angle < 0.2 rad */
					/* Knife-edge attenuation for 2 obstacles */

					if (prop.hht < 3400) { /* TX peak below vegetation limit */
						fresnel_param = interp_q * abs(slant_tx_obs2 + slant_hilltop1 - slant_txrx_obs);
						diffraction_loss2 = -18.0 + foliage_scatter2 * aknfe(fresnel_param);
					} else {
						fresnel_param = interp_q * abs(slant_tx_obs2 + slant_hilltop1 - slant_txrx_obs);
						diffraction_loss2 = aknfe(fresnel_param);
					}

					if (prop.hhr < 3400) {
						fresnel_param = interp_q * abs(slant_rx_obs2 + slant_hilltop2 - slant_rxtx_obs);
						diffraction_loss2 += (-18.0 + foliage_scatter2 * aknfe(fresnel_param));
					} else {
						fresnel_param = interp_q * abs(slant_rx_obs2 + slant_hilltop2 - slant_rxtx_obs);
						diffraction_loss2 += aknfe(fresnel_param);
					}
					/* Add vegetation attenuation */
					clutter_loss = saalos(rx_dist_beyond, prop, propa);
					diffraction_loss2 += MIN(22.0, clutter_loss);

				} else { /* RX too close to the 2nd obstacle */

					/* Knife-edge for the 1st obstacle */
					if (prop.hht < 3400) {
						fresnel_param = interp_q * abs(slant_tx_obs2 + slant_hilltop1 - slant_txrx_obs);
						diffraction_loss2 = -18.0 + foliage_scatter2 * aknfe(fresnel_param);
					} else {
						fresnel_param = interp_q * abs(slant_tx_obs2 + slant_hilltop1 - slant_txrx_obs);
						diffraction_loss2 = aknfe(fresnel_param);
					}

					/* Clutter path loss beyond the 2nd peak */
					if (prop.the[1] < 1.22) {
						rx_dist_beyond = prop.dl[1];

						if (prop.the[1] > 0.6) { /* through vegetation going downhill */
							prop.tgh = prop.cch;
						} else { /* near vegetation, RX in vegetation on slope */
							fresnel_param = 0.6365 * prop.wn *
							    abs(slant_rx_obs2 + slant_hilltop2 - slant_rxtx_obs);
						}
						diffraction_loss2 += aknfe(fresnel_param);
						clutter_loss = saalos(rx_dist_beyond, prop, propa);
						diffraction_loss2 += MIN(clutter_loss, 22.0);
					} else { /* RX very close to a cliff or building */
						diffraction_loss2 = 5.8 + 25.0;
					}
				}
			} else { /* RX at the top of the 2nd obstacle */
				fresnel_param = 0.6365 * prop.wn * abs(slant_tx_obs + slant_rx_obs - slant_total);
				diffraction_loss2 = 5.8 + aknfe(fresnel_param);
			}
		} else { /* Single obstacle */

			if (int (prop.dl[1]) > 0.0) { /* RX beyond the 1st peak */

				if (prop.the[1] < 0.2) { /* grazing angle < 0.2 rad */
					fresnel_param = 0.6365 * prop.wn * abs(slant_tx_obs + slant_rx_obs - slant_total);

					if (prop.hht < 3400) {
						scatter_voltage = 18.0;
						scatter_voltage = pow(10, (-scatter_voltage / 20));
						/* Knife-edge phase shift relative to direct LOS */
						ke_phase_diff =
						    0.159155 * prop.wn *
						    abs(slant_tx_obs + slant_rx_obs - slant_total);
						ke_amplitude = abs(ke_phase_diff - (int (ke_phase_diff)));
						ke_voltage = aknfe(fresnel_param);
						ke_voltage = pow(10, (-ke_voltage / 20));
						/* Scatter path phase shift relative to direct LOS */
						scatter_delay =
						    0.5 +
						    0.159155 * prop.wn *
						    abs(slant_tx_obs_base + slant_rx_obs_base - slant_total);
						scatter_amplitude = abs(scatter_delay - (int (scatter_delay)));
						/* Phase difference between scatter and knife-edge (rad) */
						phase_diff_rad = 6.283185307 * abs(scatter_amplitude - ke_amplitude);
						/* Clamp pd between 0 and pi, correct for 3rd and 4th quadrant */
						if (phase_diff_rad >= 3.141592654) {
							phase_diff_rad = 6.283185307 - phase_diff_rad;
							combined_scatter_dB =
							    abq_alos(complex <double>(scatter_voltage, 0) +
								     complex <double>(ke_voltage * -cos(phase_diff_rad),
									       ke_voltage * -sin(phase_diff_rad)));
						} else {
							combined_scatter_dB =
							    abq_alos(complex <double>(scatter_voltage, 0) +
								     complex <double>(ke_voltage * cos(phase_diff_rad),
									       ke_voltage * sin(phase_diff_rad)));
						}
						/* combined_scatter_dB=MAX(combined_scatter_dB,0.0009); // limit to 30.45 dB */
						diffraction_loss2 = -3.71 - 10 * log10(combined_scatter_dB);
					} else {
						diffraction_loss2 = aknfe(fresnel_param);
					}
					/* Add vegetation attenuation */
					clutter_loss = saalos(rx_dist_beyond, prop, propa);
					diffraction_loss2 += MIN(clutter_loss, 22.0);

				} else { /* grazing angle too high */

					if (prop.the[1] < 1.22) {
						rx_dist_beyond = prop.dl[1];

						if (prop.the[1] > 0.6) { /* through vegetation going downhill */
							prop.tgh = prop.cch;
						} else { /* on slope above vegetation */
							fresnel_param = 0.6365 * prop.wn *
							    abs(slant_tx_obs + slant_rx_obs - slant_total);
							diffraction_loss2 = aknfe(fresnel_param);
						}
						clutter_loss = saalos(rx_dist_beyond, prop, propa);
						diffraction_loss2 += MIN(22.0, clutter_loss);
					} else { /* RX very close to a cliff or building */
						diffraction_loss2 = 5.8 + 25.0;
					}
				}
			} else { /* RX at the top of the 1st obstacle */
				diffraction_loss2 = 5.8;
			}
		}
	}
	return diffraction_loss2;
}

/*
 * ascat - Tropospheric scatter attenuation
 *
 * Computes the attenuation (dB) due to tropospheric scatter beyond the horizon.
 * Must be called first with d==0.0 to initialize the static coefficients.
 *
 * Parameters:
 *   path_dist (former: d) : path distance (m), 0 for initialization
 *   prop                  : propagation parameters
 *   propa                 : intermediate coefficients
 * Returns: troposcatter attenuation in dB (or 1001.0 if not applicable)
 */
double ascat(double path_dist, prop_type & prop, propa_type & propa) /* former: d */
{
	/* Thread-local static coefficients, initialized at path_dist==0 */
	static __thread double
	    horizon_asymmetry,    /* former: ad  — horizon distance asymmetry */
	    horizon_ratio,        /* former: rr  — effective height ratio */
	    turbulence_factor,    /* former: etq — empirical turbulence factor */
	    h0_prev;              /* former: h0s — H0 from previous iteration */

	double h0_current;        /* former: h0  — current H0 factor */
	double radius1, radius2;  /* former: r1, r2 — TX and RX Fresnel radii */
	double scatter_height;    /* former: z0  — scatter volume height */
	double distance_ratio;    /* former: ss  — distance ratio */
	double height_ratio;      /* former: et  — turbulence/height ratio */
	double height_ratio_ceil; /* former: ett — ett = max(et, 1.0) */
	double horizon_angle;     /* former: th  — total grazing angle */
	double interp_q;          /* former: q   */
	double scatter_loss;      /* former: ascatv */
	double temp;

	if (path_dist == 0.0) {
		horizon_asymmetry = prop.dl[0] - prop.dl[1];
		horizon_ratio = prop.he[1] / prop.rch[0];

		if (horizon_asymmetry < 0.0) {
			horizon_asymmetry = -horizon_asymmetry;
			horizon_ratio = 1.0 / horizon_ratio;
		}

		turbulence_factor = (5.67e-6 * prop.ens - 2.32e-3) * prop.ens + 0.031;
		h0_prev = -15.0;
		scatter_loss = 0.0;
	}

	else {
		if (h0_prev > 15.0)
			h0_current = h0_prev;
		else {
			horizon_angle = prop.the[0] + prop.the[1] + path_dist * prop.gme;
			radius2 = 2.0 * prop.wn * horizon_angle;
			radius1 = radius2 * prop.he[0];
			radius2 *= prop.he[1];

			if (radius1 < 0.2 && radius2 < 0.2)
				return 1001.0; /* early exit: no scatter possible */

			distance_ratio = (path_dist - horizon_asymmetry) / (path_dist + horizon_asymmetry);
			interp_q = horizon_ratio / distance_ratio;
			distance_ratio = MAX(0.1, distance_ratio);
			interp_q = MIN(MAX(0.1, interp_q), 10.0);
			scatter_height = (path_dist - horizon_asymmetry) * (path_dist + horizon_asymmetry) * horizon_angle * 0.25 / path_dist;

			temp = MIN(1.7, scatter_height / 8.0e3);
			temp = temp * temp * temp * temp * temp * temp;
			height_ratio = (turbulence_factor * exp(-temp) + 1.0) * scatter_height / 1.7556e3;

			height_ratio_ceil = MAX(height_ratio, 1.0);
			h0_current = (h0f(radius1, height_ratio_ceil) + h0f(radius2, height_ratio_ceil)) * 0.5;
			h0_current +=
			    MIN(h0_current,
				  (1.38 - log(height_ratio_ceil)) * log(distance_ratio) * log(interp_q) * 0.49);
			h0_current = FORTRAN_DIM(h0_current, 0.0);

			if (height_ratio < 1.0) {
				temp = ((1.0 + 1.4142 / radius1) * (1.0 + 1.4142 / radius2));
				h0_current = height_ratio * h0_current + (1.0 -
						height_ratio) * 4.343 * log((temp * temp) * (radius1 + radius2)
							  / (radius1 + radius2 + 2.8284));
			}

			if (h0_current > 15.0 && h0_prev >= 0.0)
				h0_current = h0_prev;
		}

		h0_prev = h0_current;
		horizon_angle = propa.tha + path_dist * prop.gme;
		scatter_loss =
		    ahd(horizon_angle * path_dist) +
		    4.343 * log(47.7 * prop.wn * (horizon_angle * horizon_angle * horizon_angle * horizon_angle)) -
		    0.1 * (prop.ens - 301.0) * exp(-horizon_angle * path_dist / 40e3) + h0_current;
	}

	return scatter_loss;
}

/*
 * qerfi - Inverse of the complementary normal distribution function (Inverse Q-function / probit)
 *
 * Computes the inverse of Q(x) = 1 - Φ(x) (complementary CDF of the
 * standard normal distribution). Implemented via Hastings polynomial approximation.
 *
 * Parameter:
 *   cumulative_prob (former: q) : cumulative probability (0 < q < 1)
 * Returns: value z such that Q(z) = q
 */
double qerfi(double cumulative_prob) /* former: q */
{
    static const double c0 = 2.515516698;
    static const double c1 = 0.802853;
    static const double c2 = 0.010328;
    static const double d1 = 1.432788;
    static const double d2 = 0.189269;
    static const double d3 = 0.001308;

    double centered_prob = 0.5 - cumulative_prob;              /* ancien : x */
    double approx_param  = sqrt(-2.0 * log(MAX(0.5 - fabs(centered_prob), 1e-6))); /* ancien : t */
    double result_value  = approx_param - ((c2 * approx_param + c1) * approx_param + c0) /
                           (((d3 * approx_param + d2) * approx_param + d1) * approx_param + 1.0); /* ancien : v */

    return centered_prob < 0.0 ? -result_value : result_value;
}

/*
 * qlrps - Radio-physical propagation parameter initialization
 *         (Quasi-smooth earth Long-Range Propagation Setup)
 *
 * Computes and stores in prop the parameters derived from physical conditions:
 *   - wave number wn (from frequency)
 *   - surface refractivity ens (corrected with elevation zsys)
 *   - effective Earth curvature gme (depends on ens)
 *   - normalized surface impedance zgnd (depends on eps, sgm, wn, pol)
 *
 * Parameters:
 *   freq_mhz (former: fmhz)           : frequency in MHz
 *   sys_elevation (former: zsys)      : mean system elevation (m)
 *   surf_refractivity (former: en0)   : sea-level surface refractivity (N-units)
 *   polarization (former: ipol)       : polarization (0=H, otherwise=V)
 *   relative_permittivity (former: eps) : relative permittivity of ground
 *   ground_conductivity (former: sgm) : ground conductivity (S/m)
 *   prop                              : parameter structure (modified on output)
 */
void qlrps(double freq_mhz, double sys_elevation, double surf_refractivity,
           int polarization, double relative_permittivity,
           double ground_conductivity, prop_type & prop)
/* former parameters: fmhz, zsys, en0, ipol, eps, sgm */
{
	double earth_curvature_constant = 157e-9; /* former: gma */

	prop.wn  = freq_mhz / 47.7;
	prop.ens = surf_refractivity;

	if (sys_elevation != 0.0)
		prop.ens *= exp(-sys_elevation / 9460.0);

	prop.gme = earth_curvature_constant * (1.0 - 0.04665 * exp(prop.ens / 179.3));

	complex < double >impedance_factor, prop_zgnd(prop.zgndreal, prop.zgndimag); /* former: zq, prop_zgnd */
	impedance_factor = complex < double >(relative_permittivity, 376.62 * ground_conductivity / prop.wn);
	prop_zgnd = sqrt(impedance_factor - 1.0);

	if (polarization != 0.0)
		prop_zgnd = prop_zgnd / impedance_factor;

	prop.zgndreal = prop_zgnd.real();
	prop.zgndimag = prop_zgnd.imag();
}

/*
 * alos - Line-of-sight attenuation (original ITM model)
 *        (Above Line-of-Sight attenuation)
 *
 * Computes the attenuation (dB) on a line-of-sight path, accounting for
 * ground reflection (2-ray model) and terrain roughness.
 * Must be called first with d==0.0 to initialize the static coefficient wls.
 *
 * Parameters:
 *   path_dist (former: d) : path distance (m), 0 for initialization
 *   prop                  : propagation parameters
 *   propa                 : intermediate coefficients
 * Returns: LOS attenuation in dB
 */
double alos(double path_dist, prop_type & prop, propa_type & propa) /* former: d */
{
	complex < double >prop_zgnd(prop.zgndreal, prop.zgndimag);
	static __thread double terrain_weight; /* former: wls — terrain roughness weight */
	complex < double >reflection_coeff;   /* former: r   */
	double terrain_roughness;             /* former: s   */
	double sin_grazing_angle;             /* former: sps */
	double interp_q;                      /* former: q   */
	double los_attenuation;               /* former: alosv */

	if (path_dist == 0.0) {
		terrain_weight =
		    0.021 / (0.021 +
			     prop.wn * prop.dh / MAX(10e3, propa.dlsa));
		los_attenuation = 0.0;
	}

	else {
		interp_q = (1.0 - 0.8 * exp(-path_dist / 50e3)) * prop.dh;
		terrain_roughness = 0.78 * interp_q * exp(-pow(interp_q / 16.0, 0.25));
		interp_q = prop.he[0] + prop.he[1];
		sin_grazing_angle = interp_q / sqrt(path_dist * path_dist + interp_q * interp_q);
		reflection_coeff = (sin_grazing_angle - prop_zgnd) / (sin_grazing_angle +
				prop_zgnd) * exp(-MIN(10.0,
							prop.wn * terrain_roughness *
							sin_grazing_angle));
		interp_q = abq_alos(reflection_coeff);

		if (interp_q < 0.25 || interp_q < sin_grazing_angle)
			reflection_coeff = reflection_coeff * sqrt(sin_grazing_angle / interp_q);

		los_attenuation = propa.emd * path_dist + propa.aed;
		interp_q = prop.wn * prop.he[0] * prop.he[1] * 2.0 / path_dist;

		if (interp_q > 1.57)
			interp_q = 3.14 - 2.4649 / interp_q;

		los_attenuation =
		    (-4.343 *
		     log(abq_alos(complex < double >(cos(interp_q), -sin(interp_q)) + reflection_coeff)) -
		     los_attenuation) * terrain_weight + los_attenuation;
	}
	return los_attenuation;
}

/*
 * alos2 - Line-of-sight attenuation (improved ITWOM model)
 *         (Above Line-of-Sight attenuation version 2)
 *
 * Improved version of alos() for the ITWOM model.
 * Incorporates vegetation cover correction (saalos) when the RX antenna
 * is below the canopy and terrain angles are low.
 * Also handles reflection correction in point-to-point mode (mdp < 0).
 *
 * Parameters:
 *   path_dist (former: d) : path distance (m), 0 for initialization
 *   prop                  : propagation parameters
 *   propa                 : intermediate coefficients
 * Returns: LOS attenuation in dB (capped at 22 dB)
 */
double alos2(double path_dist, prop_type & prop, propa_type & propa) /* former: d */
{
	complex < double >prop_zgnd(prop.zgndreal, prop.zgndimag);
	complex < double >reflection_coeff;  /* former: r     */
	double clutter_dist_tx;              /* former: cd    — TX distance inside canopy */
	double clutter_dist_rx;              /* former: cr    — RX distance inside canopy */
	double reflect_point_dist;           /* former: dr    — distance to reflection point */
	double rx_antenna_height;            /* former: hr    — effective RX height */
	double rx_ground_height;             /* former: hrg   — RX height above ground */
	double tx_antenna_height;            /* former: ht    — effective TX height */
	double tx_ground_height;             /* former: htg   — TX height above ground */
	double reflect_point_height;         /* former: hrp   — reflection point height */
	double path_phase_factor;            /* former: re    — phase factor (squared magnitude) */
	double terrain_roughness;            /* former: s     */
	double sin_grazing_angle;            /* former: sps   */
	double interp_q;                     /* former: q     */
	double total_path_dist;              /* former: pd    — total distance */
	double reflect_earth_curve;          /* former: drh   — height due to Earth curvature at reflection point */
	double los_attenuation;              /* former: alosv */

	clutter_dist_tx = 0.0;
	clutter_dist_rx = 0.0;
	tx_ground_height = prop.hg[0];
	rx_ground_height = prop.hg[1];
	tx_antenna_height = prop.ght;
	rx_antenna_height = prop.ghr;
	/* rp=prop.rpl; */
	reflect_point_height = prop.rph;
	total_path_dist = prop.dist;

	if (path_dist == 0.0) {
		los_attenuation = 0.0;
	}

	else {
		interp_q = prop.he[0] + prop.he[1];
		sin_grazing_angle = interp_q / sqrt(total_path_dist * total_path_dist + interp_q * interp_q);
		interp_q = (1.0 - 0.8 * exp(-total_path_dist / 50e3)) * prop.dh;

		if (prop.mdp < 0) { /* point-to-point mode */
			reflect_point_dist = total_path_dist / (1 + rx_ground_height / tx_ground_height);

			if (reflect_point_dist < (0.5 * total_path_dist)) {
				reflect_earth_curve =
				    6378137.0 - sqrt(-(0.5 * total_path_dist) * (0.5 * total_path_dist) +
						     6378137.0 * 6378137.0 +
						     (0.5 * total_path_dist -
						      reflect_point_dist) * (0.5 * total_path_dist - reflect_point_dist));
			} else {
				reflect_earth_curve =
				    6378137.0 - sqrt(-(0.5 * total_path_dist) * (0.5 * total_path_dist) +
						     6378137.0 * 6378137.0 +
						     (reflect_point_dist - 0.5 * total_path_dist) * (reflect_point_dist - 0.5 * total_path_dist));
			}

			/* If far from TX and RX below canopy */
			if ((sin_grazing_angle < 0.05) && (prop.cch > rx_ground_height) && (prop.dist < prop.dl[0])) {
				clutter_dist_tx = MAX(0.01, total_path_dist * (prop.cch - rx_ground_height) / (tx_ground_height - rx_ground_height));
				clutter_dist_rx = MAX(0.01, total_path_dist - reflect_point_dist + reflect_point_dist * (prop.cch - reflect_earth_curve) / tx_ground_height);
				interp_q = ((1.0 - 0.8 * exp(-total_path_dist / 50e3)) * prop.dh *
				     (MIN(-20 * log10(clutter_dist_tx / clutter_dist_rx), 1.0)));
			}
		}

		terrain_roughness = 0.78 * interp_q * exp(-pow(interp_q / 16.0, 0.25));
		interp_q = exp(-MIN(10.0, prop.wn * terrain_roughness * sin_grazing_angle));
		reflection_coeff = interp_q * (sin_grazing_angle - prop_zgnd) / (sin_grazing_angle + prop_zgnd);
		interp_q = abq_alos(reflection_coeff);
		interp_q = MIN(interp_q, 1.0);

		if (interp_q < 0.25 || interp_q < sin_grazing_angle) {
			reflection_coeff = reflection_coeff * sqrt(sin_grazing_angle / interp_q);
		}

		interp_q = prop.wn * prop.he[0] * prop.he[1] / (total_path_dist * 3.1415926535897);

		if (prop.mdp < 0) { /* point-to-point mode: use actual antenna heights */
			interp_q = prop.wn * ((tx_antenna_height - reflect_point_height) * (rx_antenna_height - reflect_point_height)) /
			           (total_path_dist * 3.1415926535897);
		}
		interp_q -= floor(interp_q);

		if (interp_q < 0.5) {
			interp_q *= 3.1415926535897;
		} else {
			interp_q = (1 - interp_q) * 3.1415926535897;
		}

		/* Note: positive sign before sin() — complex conjugate removed */
		path_phase_factor = abq_alos(complex < double >(cos(interp_q), sin(interp_q)) + reflection_coeff);
		los_attenuation = -10 * log10(path_phase_factor);
		prop.tgh  = prop.hg[0];                         /* TX height AGL */
		prop.tsgh = prop.rch[0] - prop.hg[0];           /* TX site elevation (AMSL) */

		/* If RX below canopy and terrain angles are low: add saalos */
		if ((prop.hg[1] < prop.cch) && (prop.thera < 0.785) && (prop.thenr < 0.785)) {
			if (sin_grazing_angle < 0.05) {
				los_attenuation = los_attenuation + saalos(total_path_dist, prop, propa);
			} else {
				los_attenuation = saalos(total_path_dist, prop, propa);
			}
		}
	}
	los_attenuation = MIN(22.0, los_attenuation);
	return los_attenuation;
}

/*
 * qlra - Antenna parameter initialization for area mode
 *        (Quasi-smooth earth Long-Range Antenna setup)
 *
 * Computes effective antenna heights he[], horizon distances dl[],
 * and grazing angles the[] based on site type (kst).
 * Updates variability flags in propv.
 *
 * Parameters:
 *   antenna_site_type (former: kst[])  : site type [0]=TX, [1]=RX
 *                                        (0=flat terrain, 1=sea shore, 2=hilltop)
 *   climate_index (former: klimx)      : climate index (1-7), 0=unchanged
 *   variability_mode (former: mdvarx)  : variability mode, <0=unchanged
 *   prop                               : propagation parameters (modified)
 *   propv                              : variability parameters (modified)
 */
void qlra(int antenna_site_type[], int climate_index, int variability_mode,
          prop_type & prop, propv_type & propv)
/* former parameters: kst, klimx, mdvarx */
{
	double interp_q; /* former: q */

	for (int j = 0; j < 2; ++j) {
		if (antenna_site_type[j] <= 0)
			prop.he[j] = prop.hg[j];
		else {
			interp_q = 4.0;

			if (antenna_site_type[j] != 1)
				interp_q = 9.0;

			if (prop.hg[j] < 5.0)
				interp_q *= sin(0.3141593 * prop.hg[j]);

			prop.he[j] =
			    prop.hg[j] + (1.0 +
					  interp_q) * exp(-MIN(20.0,
							    2.0 * prop.hg[j] /
							    MAX(1e-3, prop.dh)));
		}

		interp_q = sqrt(2.0 * prop.he[j] / prop.gme);
		prop.dl[j] =
		    interp_q * exp(-0.07 * sqrt(prop.dh / MAX(prop.he[j], 5.0)));
		prop.the[j] =
		    (0.65 * prop.dh * (interp_q / prop.dl[j] - 1.0) -
		     2.0 * prop.he[j]) / interp_q;
	}

	prop.mdp = 1;
	propv.lvar = MAX(propv.lvar, 3);

	if (variability_mode >= 0) {
		propv.mdvar = variability_mode;
		propv.lvar = MAX(propv.lvar, 4);
	}

	if (climate_index > 0) {
		propv.klim = climate_index;
		propv.lvar = 5;
	}
}

/*
 * lrprop - Longley-Rice propagation computation (original ITM model)
 *          (Long-Range Propagation)
 *
 * Computes the reference attenuation prop.aref (dB) for a given distance d,
 * combining the LOS (alos), diffraction (adiff) and troposcatter (ascat) models.
 * Must be called first with d != 0 when prop.mdp != 0 to initialize
 * the coefficients, then with prop.dist for each distance to compute.
 *
 * Parameters:
 *   path_dist (former: d) : path distance (m)
 *   prop                  : propagation parameters (modified)
 *   propa                 : intermediate coefficients (modified)
 */
void lrprop(double path_dist, prop_type & prop, propa_type & propa) /* former: d */
{
	/* PaulM_lrprop used for ITM */
	static __thread bool los_initialized, scatter_initialized;   /* former: wlos, wscat */
	static __thread double min_valid_dist, airy_earth_radius;    /* former: dmin, xae */
	complex < double >prop_zgnd(prop.zgndreal, prop.zgndimag);
	/* Calculation points for LOS and scatter curve fitting */
	double a0, a1, a2, a3, a4, a5, a6;
	double d0, d1, d2, d3, d4, d5, d6;
	bool curve_valid; /* ancien : wq */
	double interp_q;  /* ancien : q  */
	int j;

	if (prop.mdp != 0) {
		/* Compute smooth-Earth line-of-sight distances */
		for (j = 0; j < 2; j++)
			propa.dls[j] = sqrt(2.0 * prop.he[j] / prop.gme);

		propa.dlsa = propa.dls[0] + propa.dls[1];
		propa.dla  = prop.dl[0] + prop.dl[1];
		propa.tha  = MAX(prop.the[0] + prop.the[1], -propa.dla * prop.gme);
		los_initialized    = false;
		scatter_initialized = false;

		/* Check for out-of-range parameters */
		if (prop.wn < 0.838 || prop.wn > 210.0)
			prop.kwx = MAX(prop.kwx, 1);

		for (j = 0; j < 2; j++)
			if (prop.hg[j] < 1.0 || prop.hg[j] > 1000.0)
				prop.kwx = MAX(prop.kwx, 1);

		for (j = 0; j < 2; j++)
			if (abs(prop.the[j]) > 200e-3
			    || prop.dl[j] < 0.1 * propa.dls[j]
			    || prop.dl[j] > 3.0 * propa.dls[j])
				prop.kwx = MAX(prop.kwx, 3);

		if (prop.ens < 250.0 || prop.ens > 400.0 || prop.gme < 75e-9
		    || prop.gme > 250e-9
		    || prop_zgnd.real() <= abs(prop_zgnd.imag())
		    || prop.wn < 0.419 || prop.wn > 420.0)
			prop.kwx = 4;

		for (j = 0; j < 2; j++)
			if (prop.hg[j] < 0.5 || prop.hg[j] > 3000.0)
				prop.kwx = 4;

		min_valid_dist = abs(prop.he[0] - prop.he[1]) / 200e-3;
		interp_q = adiff(0.0, prop, propa);
		/* xae=pow(prop.wn*pow(prop.gme,2.),-THIRD); */
		airy_earth_radius = pow(prop.wn * (prop.gme * prop.gme), -THIRD);
		d3 = MAX(propa.dlsa, 1.3787 * airy_earth_radius + propa.dla);
		d4 = d3 + 2.7574 * airy_earth_radius;
		a3 = adiff(d3, prop, propa);
		a4 = adiff(d4, prop, propa);
		propa.emd = (a4 - a3) / (d4 - d3);
		propa.aed = a3 - propa.emd * d3;
	}

	if (prop.mdp >= 0) {
		prop.mdp  = 0;
		prop.dist = path_dist;
	}

	if (prop.dist > 0.0) {
		if (prop.dist > 1000e3)
			prop.kwx = MAX(prop.kwx, 1);

		if (prop.dist < min_valid_dist)
			prop.kwx = MAX(prop.kwx, 3);

		if (prop.dist < 1e3 || prop.dist > 2000e3)
			prop.kwx = 4;
	}

	/* LOS zone */
	if (prop.dist < propa.dlsa) {
		if (!los_initialized) {
			interp_q = alos(0.0, prop, propa);
			d2 = propa.dlsa;
			a2 = propa.aed + d2 * propa.emd;
			d0 = 1.908 * prop.wn * prop.he[0] * prop.he[1];

			if (propa.aed >= 0.0) {
				d0 = MIN(d0, 0.5 * propa.dla);
				d1 = d0 + 0.25 * (propa.dla - d0);
			} else
				d1 = MAX(-propa.aed / propa.emd, 0.25 * propa.dla);

			a1 = alos(d1, prop, propa);
			curve_valid = false;

			if (d0 < d1) {
				a0 = alos(d0, prop, propa);
				interp_q = log(d2 / d0);
				propa.ak2 =
				    MAX(0.0,
					  ((d2 - d0) * (a1 - a0) -
					   (d1 - d0) * (a2 - a0)) / ((d2 - d0) * log(d1 / d0) - (d1 - d0) * interp_q));
				curve_valid = propa.aed >= 0.0 || propa.ak2 > 0.0;

				if (curve_valid) {
					propa.ak1 = (a2 - a0 - propa.ak2 * interp_q) / (d2 - d0);

					if (propa.ak1 < 0.0) {
						propa.ak1 = 0.0;
						propa.ak2 = FORTRAN_DIM(a2, a0) / interp_q;

						if (propa.ak2 == 0.0)
							propa.ak1 = propa.emd;
					}
				} else {
					propa.ak2 = 0.0;
					propa.ak1 = (a2 - a1) / (d2 - d1);

					if (propa.ak1 <= 0.0)
						propa.ak1 = propa.emd;
				}
			} else {
				propa.ak1 = (a2 - a1) / (d2 - d1);
				propa.ak2 = 0.0;

				if (propa.ak1 <= 0.0)
					propa.ak1 = propa.emd;
			}

			propa.ael     = a2 - propa.ak1 * d2 - propa.ak2 * log(d2);
			los_initialized = true;
		}

		if (prop.dist > 0.0)
			prop.aref =
			    propa.ael + propa.ak1 * prop.dist +
			    propa.ak2 * log(prop.dist);
	}

	/* Troposcatter zone */
	if (prop.dist <= 0.0 || prop.dist >= propa.dlsa) {
		if (!scatter_initialized) {
			interp_q = ascat(0.0, prop, propa);
			d5 = propa.dla + 200e3;
			d6 = d5 + 200e3;
			a6 = ascat(d6, prop, propa);
			a5 = ascat(d5, prop, propa);

			if (a5 < 1000.0) {
				propa.ems = (a6 - a5) / 200e3;
				propa.dx =
				    MAX(propa.dlsa,
					  MAX(propa.dla +
						0.3 * airy_earth_radius * log(47.7 * prop.wn),
						(a5 - propa.aed -
						 propa.ems * d5) / (propa.emd - propa.ems)));
				propa.aes =
				    (propa.emd - propa.ems) * propa.dx +
				    propa.aed;
			} else {
				propa.ems = propa.emd;
				propa.aes = propa.aed;
				propa.dx  = 10.e6;
			}

			scatter_initialized = true;
		}

		if (prop.dist > propa.dx)
			prop.aref = propa.aes + propa.ems * prop.dist;
		else
			prop.aref = propa.aed + propa.emd * prop.dist;
	}

	prop.aref = MAX(prop.aref, 0.0);
}

/*
 * lrprop2 - Longley-Rice propagation computation (improved ITWOM model)
 *           (Long-Range Propagation version 2)
 *
 * ITWOM version of lrprop(). Handles area and point-to-point modes,
 * using alos2() and adiff2() instead of alos() and adiff().
 * In point-to-point mode (iw > 0), chooses between diffraction and troposcatter
 * based on the computed attenuation.
 *
 * Parameters:
 *   path_dist (former: d) : path distance (m)
 *   prop                  : propagation parameters (modified)
 *   propa                 : intermediate coefficients (modified)
 */
void lrprop2(double path_dist, prop_type & prop, propa_type & propa) /* former: d */
{
	/* ITWOM_lrprop2 */
	static __thread bool los_initialized, scatter_initialized;    /* former: wlos, wscat */
	static __thread double min_valid_dist, airy_earth_radius;     /* former: dmin, xae */
	complex < double >prop_zgnd(prop.zgndreal, prop.zgndimag);
	double current_dist;    /* former: pd1 — current distance */
	double a0, a1, a2, a3, a4, a5, a6;
	double d0, d1, d2, d3, d4, d5, d6;
	double interval_width;  /* former: iw  — profile interval width */
	bool curve_valid;       /* former: wq  */
	double interp_q;        /* former: q   */
	int j;

	interval_width = prop.tiw;
	current_dist   = prop.dist;
	propa.dx       = 2000000.0;

	if (prop.mdp != 0) { /* non-area mode in progress */
		for (j = 0; j < 2; j++)
			propa.dls[j] = sqrt(2.0 * prop.he[j] / prop.gme);

		propa.dlsa = propa.dls[0] + propa.dls[1];
		propa.dlsa = MIN(propa.dlsa, 1000000.0);
		propa.dla  = prop.dl[0] + prop.dl[1];
		propa.tha  = MAX(prop.the[0] + prop.the[1], -propa.dla * prop.gme);
		los_initialized     = false;
		scatter_initialized = false;

		/* Check for out-of-range parameters */
		if (prop.wn < 0.838 || prop.wn > 210.0)
			prop.kwx = MAX(prop.kwx, 1);

		for (j = 0; j < 2; j++)
			if (prop.hg[j] < 1.0 || prop.hg[j] > 1000.0)
				prop.kwx = MAX(prop.kwx, 1);

		if (abs(prop.the[0]) > 200e-3)
			prop.kwx = MAX(prop.kwx, 3);

		if (abs(prop.the[1]) > 1.220)
			prop.kwx = MAX(prop.kwx, 3);

		if (prop.ens < 250.0 || prop.ens > 400.0 || prop.gme < 75e-9
		    || prop.gme > 250e-9
		    || prop_zgnd.real() <= abs(prop_zgnd.imag())
		    || prop.wn < 0.419 || prop.wn > 420.0)
			prop.kwx = 4;

		for (j = 0; j < 2; j++)
			if (prop.hg[j] < 0.5 || prop.hg[j] > 3000.0)
				prop.kwx = 4;

		min_valid_dist    = abs(prop.he[0] - prop.he[1]) / 200e-3;
		interp_q          = adiff2(0.0, prop, propa);
		airy_earth_radius = pow(prop.wn * (prop.gme * prop.gme), -THIRD);
		d3 = MAX(propa.dlsa, 1.3787 * airy_earth_radius + propa.dla);
		d4 = d3 + 2.7574 * airy_earth_radius;
		a3 = adiff2(d3, prop, propa);
		a4 = adiff2(d4, prop, propa);
		propa.emd = (a4 - a3) / (d4 - d3);
		propa.aed = a3 - propa.emd * d3;
	}

	if (prop.mdp >= 0) { /* area mode initialization */
		prop.mdp  = 0;
		prop.dist = path_dist;
	}

	if (prop.dist > 0.0) {
		if (prop.dist > 1000e3) /* > 1000 km: warning */
			prop.kwx = MAX(prop.kwx, 1);

		if (prop.dist < min_valid_dist)
			prop.kwx = MAX(prop.kwx, 3);

		if (prop.dist < 1e3 || prop.dist > 2000e3)
			prop.kwx = 4;
	}

	/* LOS zone */
	if (prop.dist < propa.dlsa) {

		if (interval_width <= 0.0) { /* area mode */

			if (!los_initialized) {
				interp_q = alos2(0.0, prop, propa);
				d2 = propa.dlsa;
				a2 = propa.aed + d2 * propa.emd;
				d0 = 1.908 * prop.wn * prop.he[0] * prop.he[1];

				if (propa.aed > 0.0) {
					prop.aref = propa.aed + propa.emd * prop.dist;
				} else {
					if (propa.aed == 0.0) {
						d0 = MIN(d0, 0.5 * propa.dla);
						d1 = d0 + 0.25 * (propa.dla - d0);
					} else { /* aed < 0 */
						d1 = MAX(-propa.aed / propa.emd, 0.25 * propa.dla);
					}
					a1 = alos2(d1, prop, propa);
					curve_valid = false;

					if (d0 < d1) {
						a0 = alos2(d0, prop, propa);
						a2 = MIN(a2, alos2(d2, prop, propa));
						interp_q = log(d2 / d0);
						propa.ak2 =
						    MAX(0.0,
							  ((d2 - d0) * (a1 - a0) -
							   (d1 - d0) * (a2 - a0)) /
							  ((d2 - d0) * log(d1 / d0) -
							   (d1 - d0) * interp_q));
						curve_valid = propa.aed >= 0.0 || propa.ak2 > 0.0;

						if (curve_valid) {
							propa.ak1 =
							    (a2 - a0 - propa.ak2 * interp_q) /
							    (d2 - d0);

							if (propa.ak1 < 0.0) {
								propa.ak1 = 0.0;
								propa.ak2 = FORTRAN_DIM(a2, a0) / interp_q;

								if (propa.ak2 == 0.0)
									propa.ak1 = propa.emd;
							}
						}
					}

					if (!curve_valid) {
						propa.ak1 = FORTRAN_DIM(a2, a1) / (d2 - d1);
						propa.ak2 = 0.0;

						if (propa.ak1 == 0.0)
							propa.ak1 = propa.emd;
					}
					propa.ael     = a2 - propa.ak1 * d2 - propa.ak2 * log(d2);
					los_initialized = true;
				}
			}
		} else { /* ITWOM point-to-point mode */

			if (!los_initialized) {
				interp_q      = alos2(0.0, prop, propa); /* initialize coefficients */
				los_initialized = true;
			}

			if (prop.los == 1) { /* line of sight */
				prop.aref = alos2(current_dist, prop, propa);
			} else {
				if (int (prop.dist - prop.dl[0]) == 0) { /* at the horizon */
					prop.aref = 5.8 + alos2(current_dist, prop, propa);
				} else if (int (prop.dist - prop.dl[0]) > 0.0) { /* beyond the horizon */
					interp_q = adiff2(0.0, prop, propa);
					prop.aref = adiff2(current_dist, prop, propa);
				} else {
					prop.aref = 1.0;
				}
			}
		}
	}

	/* Troposcatter zone */
	if (prop.dist <= 0.0 || prop.dist >= propa.dlsa) {
		if (interval_width == 0.0) { /* area mode */
			if (!scatter_initialized) {
				interp_q = ascat(0.0, prop, propa);
				d5 = propa.dla + 200e3;
				d6 = d5 + 200e3;
				a6 = ascat(d6, prop, propa);
				a5 = ascat(d5, prop, propa);

				if (a5 < 1000.0) {
					propa.ems = (a6 - a5) / 200e3;
					propa.dx =
					    MAX(propa.dlsa,
						  MAX(propa.dla +
							0.3 * airy_earth_radius * log(47.7 * prop.wn),
							(a5 - propa.aed - propa.ems * d5) /
							(propa.emd - propa.ems)));

					propa.aes =
					    (propa.emd - propa.ems) * propa.dx +
					    propa.aed;
				} else {
					propa.ems = propa.emd;
					propa.aes = propa.aed;
					propa.dx  = 10000000;
				}
				scatter_initialized = true;
			}

			if (prop.dist > propa.dx) {
				prop.aref = propa.aes + propa.ems * prop.dist;
			} else {
				prop.aref = propa.aed + propa.emd * prop.dist;
			}
		} else { /* ITWOM point-to-point mode: choose the best model */

			if (!scatter_initialized) {
				d5 = 0.0;
				d6 = 0.0;
				interp_q = ascat(0.0, prop, propa);
				a6 = ascat(current_dist, prop, propa);
				interp_q = adiff2(0.0, prop, propa);
				a5 = adiff2(current_dist, prop, propa);

				if (a5 <= a6) {
					propa.dx  = 10000000;
					prop.aref = a5; /* diffraction better or equal */
				} else {
					propa.dx  = propa.dlsa;
					prop.aref = a6; /* troposcatter better */
				}
				scatter_initialized = true;
			}
		}
	}
	prop.aref = MAX(prop.aref, 0.0);
}

/*
 * curve - Empirical variability curve
 *
 * Computes a variability value as a function of the effective reference distance de,
 * using a sigmoid curve weighted by a normalized distance term.
 *
 * Parameters:
 *   c1, c2        : sigmoid curve coefficients
 *   x1, x2, x3   : shape and position parameters
 *   effective_dist (former: de) : effective reference distance (m)
 * Returns: variability value (dB)
 */
double curve(double const &c1, double const &c2, double const &x1,
	     double const &x2, double const &x3, double const &effective_dist) /* former: de */
{
	/* return (c1+c2/(1.0+pow((de-x2)/x3,2.0)))*pow(de/x1,2.0)/(1.0+pow(de/x1,2.0)); */
	double sigmoid_term;       /* former: temp1 — sigmoid term */
	double distance_norm_term; /* former: temp2 — normalized distance term */

	sigmoid_term       = (effective_dist - x2) / x3;
	distance_norm_term = effective_dist / x1;

	sigmoid_term       *= sigmoid_term;
	distance_norm_term *= distance_norm_term;

	return (c1 + c2 / (1.0 + sigmoid_term)) * distance_norm_term / (1.0 + distance_norm_term);
}

/*
 * avar - Statistical attenuation variability computation (Attenuation Variability)
 *
 * Computes the total attenuation (dB) including statistical variability for
 * time (zzt), location (zzl) and composite (zzc) confidences.
 * Initializes variability parameters as needed (propv.lvar > 0).
 *
 * Parameters:
 *   time_z (former: zzt)       : normal deviate for time variability
 *   location_z (former: zzl)   : normal deviate for location variability
 *   composite_z (former: zzc)  : normal deviate for composite variability
 *   prop                       : propagation parameters
 *   propv                      : variability parameters
 * Returns: total attenuation in dB
 */
double avar(double time_z, double location_z, double composite_z,
            prop_type & prop, propv_type & propv)
/* former: zzt, zzl, zzc */
{
	static __thread int variability_mode_decoded; /* former: kdv */
	static __thread double
	    effective_dist_threshold,  /* former: dexa */
	    effective_dist,            /* former: de */
	    median_variability,        /* former: vmd */
	    variance_short,            /* former: vs0 */
	    sigma_location,            /* former: sgl */
	    sigma_time_minus,          /* former: sgtm */
	    sigma_time_plus,           /* former: sgtp */
	    sigma_time_decay,          /* former: sgtd */
	    threshold_time_decay,      /* former: tgtd */
	    freq_factor_minus,         /* former: gm */
	    freq_factor_plus,          /* former: gp */
	    cv1, cv2, yv1, yv2, yv3,
	    csm1, csm2, ysm1, ysm2, ysm3,
	    csp1, csp2, ysp1, ysp2, ysp3,
	    csd1, zd,
	    cfm1, cfm2, cfm3,
	    cfp1, cfp2, cfp3;

	/* Empirical arrays indexed by climate type (klim-1) */
	double bv1[7] = { -9.67, -0.62, 1.26, -9.21, -0.62, -0.39, 3.15 };
	double bv2[7] = { 12.7, 9.19, 15.5, 9.05, 9.19, 2.86, 857.9 };
	double xv1[7] = { 144.9e3, 228.9e3, 262.6e3, 84.1e3, 228.9e3, 141.7e3, 2222.e3 };
	double xv2[7] = { 190.3e3, 205.2e3, 185.2e3, 101.1e3, 205.2e3, 315.9e3, 164.8e3 };
	double xv3[7] = { 133.8e3, 143.6e3, 99.8e3, 98.6e3, 143.6e3, 167.4e3, 116.3e3 };
	double bsm1[7] = { 2.13, 2.66, 6.11, 1.98, 2.68, 6.86, 8.51 };
	double bsm2[7] = { 159.5, 7.67, 6.65, 13.11, 7.16, 10.38, 169.8 };
	double xsm1[7] = { 762.2e3, 100.4e3, 138.2e3, 139.1e3, 93.7e3, 187.8e3, 609.8e3 };
	double xsm2[7] = { 123.6e3, 172.5e3, 242.2e3, 132.7e3, 186.8e3, 169.6e3, 119.9e3 };
	double xsm3[7] = { 94.5e3, 136.4e3, 178.6e3, 193.5e3, 133.5e3, 108.9e3, 106.6e3 };
	double bsp1[7] = { 2.11, 6.87, 10.08, 3.68, 4.75, 8.58, 8.43 };
	double bsp2[7] = { 102.3, 15.53, 9.60, 159.3, 8.12, 13.97, 8.19 };
	double xsp1[7] = { 636.9e3, 138.7e3, 165.3e3, 464.4e3, 93.2e3, 216.0e3, 136.2e3 };
	double xsp2[7] = { 134.8e3, 143.7e3, 225.7e3, 93.1e3, 135.9e3, 152.0e3, 188.5e3 };
	double xsp3[7] = { 95.6e3, 98.6e3, 129.7e3, 94.2e3, 113.4e3, 122.7e3, 122.9e3 };
	double bsd1[7] = { 1.224, 0.801, 1.380, 1.000, 1.224, 1.518, 1.518 };
	double bzd1[7] = { 1.282, 2.161, 1.282, 20., 1.282, 1.282, 1.282 };
	double bfm1[7] = { 1.0, 1.0, 1.0, 1.0, 0.92, 1.0, 1.0 };
	double bfm2[7] = { 0.0, 0.0, 0.0, 0.0, 0.25, 0.0, 0.0 };
	double bfm3[7] = { 0.0, 0.0, 0.0, 0.0, 1.77, 0.0, 0.0 };
	double bfp1[7] = { 1.0, 0.93, 1.0, 0.93, 0.93, 1.0, 1.0 };
	double bfp2[7] = { 0.0, 0.31, 0.0, 0.19, 0.31, 0.0, 0.0 };
	double bfp3[7] = { 0.0, 2.00, 0.0, 1.79, 2.00, 0.0, 0.0 };

	static __thread bool suppress_short_term, suppress_location; /* former: ws, w1 */
	double time_threshold = 7.8, location_threshold = 24.0;      /* former: rt, rl */
	double total_attenuation;  /* former: avarv */
	double interp_q;           /* former: q     */
	double variance_total;     /* former: vs    */
	double zt, zl, zc;
	double sigma_time_used;    /* former: sgt   */
	double variability_offset; /* former: yr    */
	double temp1, temp2;
	int climate_index = propv.klim - 1; /* ancien : temp_klim */

	if (propv.lvar > 0) {
		switch (propv.lvar) {
		default:
			if (propv.klim <= 0 || propv.klim > 7) {
				propv.klim = 5;
				climate_index = 4;
				prop.kwx = MAX(prop.kwx, 2);
			}

			/* Load variability coefficients for this climate */
			cv1 = bv1[climate_index]; cv2 = bv2[climate_index];
			yv1 = xv1[climate_index]; yv2 = xv2[climate_index]; yv3 = xv3[climate_index];
			csm1 = bsm1[climate_index]; csm2 = bsm2[climate_index];
			ysm1 = xsm1[climate_index]; ysm2 = xsm2[climate_index]; ysm3 = xsm3[climate_index];
			csp1 = bsp1[climate_index]; csp2 = bsp2[climate_index];
			ysp1 = xsp1[climate_index]; ysp2 = xsp2[climate_index]; ysp3 = xsp3[climate_index];
			csd1 = bsd1[climate_index]; zd   = bzd1[climate_index];
			cfm1 = bfm1[climate_index]; cfm2 = bfm2[climate_index]; cfm3 = bfm3[climate_index];
			cfp1 = bfp1[climate_index]; cfp2 = bfp2[climate_index]; cfp3 = bfp3[climate_index];
			[[fallthrough]];

		case 4:
			variability_mode_decoded = propv.mdvar;
			suppress_short_term = variability_mode_decoded >= 20;

			if (suppress_short_term)
				variability_mode_decoded -= 20;

			suppress_location = variability_mode_decoded >= 10;

			if (suppress_location)
				variability_mode_decoded -= 10;

			if (variability_mode_decoded < 0 || variability_mode_decoded > 3) {
				variability_mode_decoded = 0;
				prop.kwx = MAX(prop.kwx, 2);
			}
			[[fallthrough]];

		case 3:
			interp_q = log(0.133 * prop.wn);
			freq_factor_minus = cfm1 + cfm2 / ((cfm3 * interp_q * cfm3 * interp_q) + 1.0);
			freq_factor_plus  = cfp1 + cfp2 / ((cfp3 * interp_q * cfp3 * interp_q) + 1.0);
			[[fallthrough]];

		case 2:
			effective_dist_threshold =
			    sqrt(18e6 * prop.he[0]) + sqrt(18e6 * prop.he[1]) +
			    pow((575.7e12 / prop.wn), THIRD);
			[[fallthrough]];

		case 1:
			if (prop.dist < effective_dist_threshold)
				effective_dist = 130e3 * prop.dist / effective_dist_threshold;
			else
				effective_dist = 130e3 + prop.dist - effective_dist_threshold;
		}

		median_variability = curve(cv1, cv2, yv1, yv2, yv3, effective_dist);
		sigma_time_minus   = curve(csm1, csm2, ysm1, ysm2, ysm3, effective_dist) * freq_factor_minus;
		sigma_time_plus    = curve(csp1, csp2, ysp1, ysp2, ysp3, effective_dist) * freq_factor_plus;
		sigma_time_decay   = sigma_time_plus * csd1;
		threshold_time_decay = (sigma_time_plus - sigma_time_decay) * zd;

		if (suppress_location)
			sigma_location = 0.0;
		else {
			interp_q = (1.0 - 0.8 * exp(-prop.dist / 50e3)) * prop.dh * prop.wn;
			sigma_location = 10.0 * interp_q / (interp_q + 13.0);
		}

		if (suppress_short_term)
			variance_short = 0.0;
		else {
			/* vs0=pow(5.0+3.0*exp(-de/100e3),2.0); */
			temp1 = (5.0 + 3.0 * exp(-effective_dist / 100e3));
			variance_short = temp1 * temp1;
		}

		propv.lvar = 0;
	}

	zt = time_z;
	zl = location_z;
	zc = composite_z;

	/* Decode variability mode */
	switch (variability_mode_decoded) {
	case 0:
		zt = zc;
		zl = zc;
		break;
	case 1:
		zl = zc;
		break;
	case 2:
		zl = zt;
	}

	if (fabs(zt) > 3.1 || fabs(zl) > 3.1 || fabs(zc) > 3.1)
		prop.kwx = MAX(prop.kwx, 1);

	/* Select time variability slope */
	if (zt < 0.0)
		sigma_time_used = sigma_time_minus;
	else if (zt <= zd)
		sigma_time_used = sigma_time_plus;
	else
		sigma_time_used = sigma_time_decay + threshold_time_decay / zt;

	temp1 = sigma_time_used * zt;
	temp2 = sigma_location * zl;

	variance_total = variance_short +
	                 (temp1 * temp1) / (time_threshold + zc * zc) +
	                 (temp2 * temp2) / (location_threshold + zc * zc);

	/* Compute variability offset by mode */
	if (variability_mode_decoded == 0) {
		variability_offset = 0.0;
		propv.sgc = sqrt(sigma_time_used * sigma_time_used + sigma_location * sigma_location + variance_total);
	} else if (variability_mode_decoded == 1) {
		variability_offset = sigma_time_used * zt;
		propv.sgc = sqrt(sigma_location * sigma_location + variance_total);
	} else if (variability_mode_decoded == 2) {
		variability_offset = sqrt(sigma_time_used * sigma_time_used + sigma_location * sigma_location) * zt;
		propv.sgc = sqrt(variance_total);
	} else {
		variability_offset = sigma_time_used * zt + sigma_location * zl;
		propv.sgc = sqrt(variance_total);
	}

	total_attenuation = prop.aref - median_variability - variability_offset - propv.sgc * zc;

	/* Compress negative variability */
	if (total_attenuation < 0.0)
		total_attenuation = total_attenuation * (29.0 - total_attenuation) / (29.0 - 10.0 * total_attenuation);

	return total_attenuation;
}

/*
 * hzns - Radio horizon computation (ITM model)
 *        (Horizon finder)
 *
 * Determines the grazing angles (prop.the[]) and horizon distances
 * (prop.dl[]) for TX and RX from the terrain profile.
 * Uses a maximum-angle search from each endpoint.
 *
 * Parameters:
 *   terrain_profile (former: pfl[]) : terrain profile [0]=num points, [1]=spacing(m), [2..n+2]=elevations(m)
 *   prop                            : propagation parameters (modified)
 */
void hzns(double terrain_profile[], prop_type & prop) /* former: pfl */
{
    const int num_points = static_cast<int>(terrain_profile[0]); /* former: np */
    if (num_points < 2) {
        prop.dl[0] = prop.dist;
        prop.dl[1] = prop.dist;
        prop.the[0] = 0.0;
        prop.the[1] = 0.0;
        return;
    }

    const double step_dist    = terrain_profile[1];        /* former: xi — point spacing (m) */
    const double total_dist   = prop.dist;                 /* former: dist */
    const double tx_elev      = terrain_profile[2] + prop.hg[0]; /* former: za — TX antenna elevation */
    const double rx_elev      = terrain_profile[num_points + 2] + prop.hg[1]; /* former: zb — RX antenna elevation */
    const double earth_curve_factor = 0.5 * prop.gme;     /* former: qc — Earth curvature factor */
    const double horizon_angle_init = earth_curve_factor * total_dist; /* former: q_dist */

    double tx_horizon_angle = (rx_elev - tx_elev) / total_dist - horizon_angle_init; /* former: the0 */
    double rx_horizon_angle = -((rx_elev - tx_elev) / total_dist) - horizon_angle_init; /* former: the1 */
    double tx_horizon_dist  = total_dist; /* former: dl0 */
    double rx_horizon_dist  = total_dist; /* former: dl1 */

    double dist_from_tx = 0.0;  /* former: sa */
    double dist_from_rx = total_dist; /* former: sb */
    bool no_horizon_found = true; /* former: wq */

    const double* profile_ptr = &terrain_profile[3]; /* pointer to elevations */

    for (int i = 1; i < num_points; ++i) {
        dist_from_tx += step_dist;
        dist_from_rx -= step_dist;

        const double elev_val = *profile_ptr++; /* current elevation */

        /* Check TX-side horizon */
        double q = elev_val - (earth_curve_factor * dist_from_tx + tx_horizon_angle) * dist_from_tx - tx_elev;
        if (q > 0.0) {
            tx_horizon_angle += q / dist_from_tx;
            tx_horizon_dist   = dist_from_tx;
            no_horizon_found  = false;
        }

        /* Check RX-side horizon (only if a TX horizon has been found) */
        if (!no_horizon_found) {
            q = elev_val - (earth_curve_factor * dist_from_rx + rx_horizon_angle) * dist_from_rx - rx_elev;
            if (q > 0.0) {
                rx_horizon_angle += q / dist_from_rx;
                rx_horizon_dist   = dist_from_rx;
            }
        }
    }

    prop.the[0] = tx_horizon_angle;
    prop.the[1] = rx_horizon_angle;
    prop.dl[0]  = tx_horizon_dist;
    prop.dl[1]  = rx_horizon_dist;
}

/*
 * hzns2 - Extended radio horizon computation (ITWOM model)
 *         (Horizon finder version 2)
 *
 * Improved version of hzns() for ITWOM.
 * Also computes: prop.hht (TX obstacle peak height), prop.hhr (RX),
 * prop.los (line-of-sight indicator), prop.rpl (reflection point index),
 * prop.rph (reflection point height).
 *
 * Parameters:
 *   terrain_profile (former: pfl[]) : terrain profile
 *   prop                            : propagation parameters (modified)
 *   propa                           : intermediate coefficients (unused)
 */
void hzns2(double terrain_profile[], prop_type & prop, propa_type & /*propa*/) /* former: pfl */
{
	bool no_tx_horizon_found; /* former: wq */
	int num_points, reflect_point_index, i, j; /* former: np, rp, i, j */
	double step_dist, tx_elev, rx_elev;        /* former: xi, za, zb */
	double earth_curve_factor, interp_q;       /* former: qc, q */
	double dist_from_rx, dist_from_tx;         /* former: sb, sa */
	double reflect_dist, inter_obstacle_gap;   /* former: dr, dshh */

	num_points     = (int)terrain_profile[0];
	step_dist      = terrain_profile[1];
	tx_elev        = terrain_profile[2] + prop.hg[0];
	rx_elev        = terrain_profile[num_points + 2] + prop.hg[1];
	prop.tiw       = step_dist;
	prop.ght       = tx_elev;
	prop.ghr       = rx_elev;
	earth_curve_factor = 0.5 * prop.gme;
	interp_q       = earth_curve_factor * prop.dist;
	prop.the[1]    = atan((rx_elev - tx_elev) / prop.dist);
	prop.the[0]    = (prop.the[1]) - interp_q;
	prop.the[1]    = -prop.the[1] - interp_q;
	prop.dl[0]     = prop.dist;
	prop.dl[1]     = prop.dist;
	prop.hht       = 0.0;
	prop.hhr       = 0.0;
	prop.los       = 1; /* initial assumption: line of sight */

	if (num_points >= 2) {
		dist_from_tx = 0.0;
		dist_from_rx = prop.dist;
		no_tx_horizon_found = true;

		/* Search for TX horizon (from TX toward RX) */
		for (j = 1; j < num_points; j++) {
			dist_from_tx += step_dist;
			interp_q = terrain_profile[j + 2] - (earth_curve_factor * dist_from_tx + prop.the[0]) * dist_from_tx - tx_elev;

			if (interp_q > 0.0) {
				prop.los = 0; /* NLOS detected */
				prop.the[0] += interp_q / dist_from_tx;
				prop.dl[0]   = dist_from_tx;
				prop.the[0]  = MIN(prop.the[0], 1.569);
				prop.hht     = terrain_profile[j + 2];
				no_tx_horizon_found = false;
			}
		}

		if (!no_tx_horizon_found) {
			/* Search for RX horizon (from RX toward TX) */
			for (i = 1; i < num_points; i++) {
				dist_from_rx -= step_dist;
				interp_q = terrain_profile[num_points + 2 - i] - (earth_curve_factor * (prop.dist - dist_from_rx) +
						       prop.the[1]) *
				    (prop.dist - dist_from_rx) - rx_elev;
				if (interp_q > 0.0) {
					prop.the[1] += interp_q / (prop.dist - dist_from_rx);
					prop.the[1]  = MIN(prop.the[1], 1.57);
					prop.the[1]  = MAX(prop.the[1], -1.568);
					prop.hhr     = terrain_profile[num_points + 2 - i];
					prop.dl[1]   = MAX(0.0, prop.dist - dist_from_rx);
				}
			}
			prop.the[0] =
			    atan((prop.hht - tx_elev) / prop.dl[0]) -
			    0.5 * prop.gme * prop.dl[0];
			prop.the[1] =
			    atan((prop.hhr - rx_elev) / prop.dl[1]) -
			    0.5 * prop.gme * prop.dl[1];
		}
	}

	/* Compute reflection point distance and its index */
	if ((prop.dl[1]) < (prop.dist)) {
		inter_obstacle_gap = prop.dist - prop.dl[0] - prop.dl[1];

		if (int (inter_obstacle_gap) == 0) { /* single obstacle */
			reflect_dist = prop.dl[1] / (1 + rx_elev / prop.hht);
		} else { /* two obstacles */
			reflect_dist = prop.dl[1] / (1 + rx_elev / prop.hhr);
		}
	} else { /* line of sight */
		reflect_dist = (prop.dist) / (1 + rx_elev / tx_elev);
	}

	reflect_point_index = 2 + (int)(floor(0.5 + reflect_dist / step_dist));
	prop.rpl = reflect_point_index;
	prop.rph = terrain_profile[reflect_point_index];
}

/*
 * z1sq1 - Linear least-squares fit on terrain profile (ITM)
 *         (Linear least-squares fit on terrain profile, version 1)
 *
 * Computes via least squares the start elevation z0 and end elevation zn of a line
 * fitted over section [x1, x2] of profile z[]. Used by qlrpfl to estimate
 * effective antenna heights.
 *
 * Parameters:
 *   terrain_profile (former: z[])       : terrain profile
 *   range_start (former: x1)            : section start (normalized distance)
 *   range_end (former: x2)              : section end
 *   slope_intercept_start (former: z0)  : fitted elevation at start (output)
 *   slope_intercept_end (former: zn)    : fitted elevation at end (output)
 */
void z1sq1(double terrain_profile[], const double &range_start, const double &range_end,
           double &slope_intercept_start, double &slope_intercept_end)
/* former: z, x1, x2, z0, zn */
{
	/* Used only with ITM 1.2.2 */
	double total_points, idx_start, idx_end; /* former: xn, xa, xb */
	double running_x, sum_elev, sum_weighted; /* former: x, a, b */
	int num_pts, start_idx, end_idx;          /* former: n, ja, jb */

	total_points = terrain_profile[0];
	idx_start    = int (MAX(range_start / terrain_profile[1], 0.0));
	idx_end      = total_points - int (FORTRAN_DIM(total_points, range_end / terrain_profile[1]));

	if (idx_end <= idx_start) {
		idx_start = FORTRAN_DIM(idx_start, 1.0);
		idx_end   = total_points - FORTRAN_DIM(total_points, idx_end + 1.0);
	}

	start_idx = (int)idx_start;
	end_idx   = (int)idx_end;
	num_pts   = end_idx - start_idx;
	idx_start = idx_end - idx_start;
	running_x = -0.5 * idx_start;
	idx_end  += running_x;

	sum_elev    = 0.5 * (terrain_profile[start_idx + 2] + terrain_profile[end_idx + 2]);
	sum_weighted = 0.5 * (terrain_profile[start_idx + 2] - terrain_profile[end_idx + 2]) * running_x;

	for (int i = 2; i <= num_pts; ++i) {
		++start_idx;
		running_x    += 1.0;
		sum_elev     += terrain_profile[start_idx + 2];
		sum_weighted += terrain_profile[start_idx + 2] * running_x;
	}

	sum_elev     /= idx_start;
	sum_weighted  = sum_weighted * 12.0 / ((idx_start * idx_start + 2.0) * idx_start);
	slope_intercept_start = sum_elev - sum_weighted * idx_end;
	slope_intercept_end   = sum_elev + sum_weighted * (total_points - idx_end);
}

/*
 * z1sq2 - Linear least-squares fit on terrain profile (ITWOM)
 *         (Linear least-squares fit on terrain profile, version 2)
 *
 * Corrected version of z1sq1() for ITWOM. Computes fitted elevations z0 and zn
 * with a slightly different algorithm (even point count, revised weighting).
 *
 * Parameters: same as z1sq1()
 */
void z1sq2(double terrain_profile[], const double &range_start, const double &range_end,
           double &slope_intercept_start, double &slope_intercept_end)
/* former: z, x1, x2, z0, zn */
{
	/* corrected for use with ITWOM */
	double total_points, idx_start, idx_end;     /* former: xn, xa, xb */
	double running_x, sum_elev, sum_weighted;    /* former: x, a, b */
	double sum_x_squared;                        /* former: bn */
	int num_pts, start_idx, end_idx;             /* former: n, ja, jb */

	total_points = terrain_profile[0];
	idx_start    = int (MAX(range_start / terrain_profile[1], 0.0));
	idx_end      = total_points - int (FORTRAN_DIM(total_points, range_end / terrain_profile[1]));

	if (idx_end <= idx_start) {
		idx_start = FORTRAN_DIM(idx_start, 1.0);
		idx_end   = total_points - FORTRAN_DIM(total_points, idx_end + 1.0);
	}

	start_idx      = (int)idx_start;
	end_idx        = (int)idx_end;
	idx_start      = (2 * int ((idx_end - idx_start) / 2)) - 1;
	running_x      = -0.5 * (idx_start + 1);
	idx_end       += running_x;
	start_idx      = end_idx - 1 - (int)idx_start;
	num_pts        = end_idx - start_idx;
	sum_elev       = (terrain_profile[start_idx + 2] + terrain_profile[end_idx + 2]);
	sum_weighted   = (terrain_profile[start_idx + 2] - terrain_profile[end_idx + 2]) * running_x;
	sum_x_squared  = 2 * (running_x * running_x);

	for (int i = 2; i <= num_pts; ++i) {
		++start_idx;
		running_x      += 1.0;
		sum_x_squared  += (running_x * running_x);
		sum_elev       += terrain_profile[start_idx + 2];
		sum_weighted   += terrain_profile[start_idx + 2] * running_x;
	}

	sum_elev     /= (idx_start + 2);
	sum_weighted  = sum_weighted / sum_x_squared;
	slope_intercept_start = sum_elev - (sum_weighted * idx_end);
	slope_intercept_end   = sum_elev + (sum_weighted * (total_points - idx_end));
}

/*
 * qtile - k-th order statistic (k-th order statistic / quickselect)
 *
 * Finds the k-th largest element of array a[] using nth_element.
 * The array is partially sorted in place.
 * Used to compute roughness percentiles in d1thx/d1thx2.
 *
 * WARNING: array a[] is modified in place (partial sort).
 *
 * Parameters:
 *   total_count (former: nn) : number of elements (indices 0..nn)
 *   values (former: a[])     : array to process (modified in place)
 *   rank_index (former: ir)  : desired rank (0-based)
 * Returns: value at rank rank_index (in descending order)
 */
double qtile(const int &total_count, double values[], const int &rank_index) /* former: nn, a, ir */
{
	int k = MAX(0, MIN(rank_index, total_count));
    std::nth_element(values, values + k, values + total_count + 1, std::greater<double>());
    return values[k];
}

/*
 * get_two_qtiles - Simultaneous computation of two order statistics (Two order statistics)
 *
 * Efficiently computes the difference between the ka-th and kb-th largest
 * elements of array a[] (ka < kb). Optimizes by making two nth_element calls
 * over reduced ranges, avoiding a full sort.
 *
 * Parameters:
 *   values (former: a[])    : input array (modified in place)
 *   count (former: n)       : number of elements
 *   rank_lower (former: ka) : rank of first element (descending order)
 *   rank_upper (former: kb) : rank of second element (descending order, > ka)
 * Returns: values[ka] - values[kb] (difference between the two quantiles)
 */
double get_two_qtiles(double values[], int count, int rank_lower, int rank_upper) /* former: a, n, ka, kb */
{
    /* First call: find the quantile at the larger rank (rank_upper) */
    std::nth_element(values, values + rank_upper, values + count, std::greater<double>());
    double upper_quantile = values[rank_upper]; /* former: res2 */

    /* Second call: search only in the left portion [0, rank_upper] */
    std::nth_element(values, values + rank_lower, values + rank_upper, std::greater<double>());
    double lower_quantile = values[rank_lower]; /* former: res1 */

    return lower_quantile - upper_quantile;
}

/*
 * qerf - Complementary normal distribution function (Q-function / complementary error function)
 *
 * Computes Q(z) = 1 - Φ(z) via Horner polynomial approximation.
 * Returns the probability that X > z for X ~ N(0,1).
 *
 * Parameter:
 *   normal_deviate (former: z) : standardized normal deviate
 * Returns: Q(z) ∈ [0, 1]
 */
double qerf(const double &normal_deviate) /* former: z */
{
	double b1 = 0.319381530, b2 = -0.356563782, b3 = 1.781477937;
	double b4 = -1.821255987, b5 = 1.330274429;
	double rp = 4.317008, rrt2pi = 0.398942280;
	double t, x, qerfv;

	x = normal_deviate;
	t = fabs(x);

	if (t >= 10.0)
		qerfv = 0.0;
	else {
		t     = rp / (t + rp);
		qerfv = exp(-0.5 * x * x) * rrt2pi *
		        ((((b5 * t + b4) * t + b3) * t + b2) * t + b1) * t;
	}

	if (x < 0.0)
		qerfv = 1.0 - qerfv;

	return qerfv;
}

/*
 * d1thx - Terrain roughness delta-h computation (ITM model)
 *         (Terrain height difference - delta h)
 *
 * Computes the interquartile terrain roughness (m) between distances x1 and x2
 * by interpolating profile pfl[], fitting a least-squares line, and computing
 * the 90th-10th percentile difference (ITM version, z1sq1).
 *
 * Parameters:
 *   terrain_profile (former: pfl[]) : terrain profile
 *   dist_start (former: x1)        : start distance (m)
 *   dist_end (former: x2)          : end distance (m)
 * Returns: delta-h roughness in meters
 */
double d1thx(double terrain_profile[], const double &dist_start, const double &dist_end) /* former: pfl, x1, x2 */
{
	int num_points, lower_rank, upper_rank, sample_count, profile_idx, j; /* former: np, ka, kb, n, k, j */
	double terrain_roughness, total_samples, start_frac, end_frac;         /* former: d1thxv, sn, xa, xb */
	double *resampled_profile; /* ancien : s */

	num_points = (int)terrain_profile[0];
	start_frac = dist_start / terrain_profile[1];
	end_frac   = dist_end   / terrain_profile[1];
	terrain_roughness = 0.0;

	if (end_frac - start_frac < 2.0) /* profile too short, return immediately */
		return terrain_roughness;

	lower_rank   = (int)(0.1 * (end_frac - start_frac + 8.0));
	lower_rank   = MIN(MAX(4, lower_rank), 25);
	sample_count = 10 * lower_rank - 5;
	upper_rank   = sample_count - lower_rank + 1;
	total_samples = sample_count - 1;
	assert((resampled_profile = new double[sample_count + 2]) != 0);
	resampled_profile[0] = total_samples;
	resampled_profile[1] = 1.0;
	end_frac   = (end_frac - start_frac) / total_samples;
	profile_idx = (int)(start_frac + 1.0);
	start_frac -= (double)profile_idx;

	/* Resample profile by linear interpolation */
	for (j = 0; j < sample_count; j++) {
		while (start_frac > 0.0 && profile_idx < num_points) {
			start_frac -= 1.0;
			++profile_idx;
		}
		resampled_profile[j + 2] = terrain_profile[profile_idx + 2] + (terrain_profile[profile_idx + 2] - terrain_profile[profile_idx + 1]) * start_frac;
		start_frac += end_frac;
	}

	/* Linear fit and detrending */
	z1sq1(resampled_profile, 0.0, total_samples, start_frac, end_frac);
	end_frac = (end_frac - start_frac) / total_samples;

	for (j = 0; j < sample_count; j++) {
		resampled_profile[j + 2] -= start_frac;
		start_frac += end_frac;
	}

	/* Compute interquartile roughness (percentile difference) */
	terrain_roughness = get_two_qtiles(resampled_profile + 2, sample_count - 1, lower_rank - 1, upper_rank - 1);
	terrain_roughness /= 1.0 - 0.8 * exp(-(dist_end - dist_start) / 50.0e3);
	delete[] resampled_profile;

	return terrain_roughness;
}

/*
 * d1thx2 - Terrain roughness delta-h computation (ITWOM model)
 *          (Terrain height difference - delta h, version 2)
 *
 * ITWOM version of d1thx(). Identical in logic but uses z1sq2()
 * for the linear fit and supports finer maximum terrain resolution
 * (adaptive kmx).
 *
 * Parameters:
 *   terrain_profile (former: pfl[]) : terrain profile
 *   dist_start (former: x1)        : start distance (m)
 *   dist_end (former: x2)          : end distance (m)
 *   propa                           : unused
 * Returns: delta-h roughness in meters
 */
double d1thx2(double terrain_profile[], const double &dist_start, const double &dist_end,
              propa_type & /*propa*/) /* former: pfl, x1, x2 */
{
	int num_points, lower_rank, upper_rank, sample_count, profile_idx, max_rank, j; /* former: np, ka, kb, n, k, kmx, j */
	double terrain_roughness, total_samples, start_frac, end_frac, interp_frac;      /* former: d1thx2v, sn, xa, xb, xc */
	double *resampled_profile; /* ancien : s */

	num_points = (int)terrain_profile[0];
	start_frac = dist_start / terrain_profile[1];
	end_frac   = dist_end   / terrain_profile[1];
	terrain_roughness = 0.0;

	if (end_frac - start_frac < 2.0) /* profile too short */
		return terrain_roughness;

	lower_rank   = (int)(0.1 * (end_frac - start_frac + 8.0));
	max_rank     = MAX(25, (int)(83350 / (terrain_profile[1])));
	lower_rank   = MIN(MAX(4, lower_rank), max_rank);
	sample_count = 10 * lower_rank - 5;
	upper_rank   = sample_count - lower_rank + 1;
	total_samples = sample_count - 1;
	assert((resampled_profile = new double[sample_count + 2]) != 0);
	resampled_profile[0] = total_samples;
	resampled_profile[1] = 1.0;
	end_frac    = (end_frac - start_frac) / total_samples;
	profile_idx = (int (start_frac + 1.0));
	interp_frac = start_frac - (double (profile_idx));

	/* Resample profile by linear interpolation */
	for (j = 0; j < sample_count; j++) {
		while (interp_frac > 0.0 && profile_idx < num_points) {
			interp_frac -= 1.0;
			++profile_idx;
		}
		resampled_profile[j + 2] = terrain_profile[profile_idx + 2] + (terrain_profile[profile_idx + 2] - terrain_profile[profile_idx + 1]) * interp_frac;
		interp_frac += end_frac;
	}

	/* Linear fit (ITWOM version) and detrending */
	z1sq2(resampled_profile, 0.0, total_samples, start_frac, end_frac);
	end_frac = (end_frac - start_frac) / total_samples;

	for (j = 0; j < sample_count; j++) {
		resampled_profile[j + 2] -= start_frac;
		start_frac += end_frac;
	}

	terrain_roughness = qtile(sample_count - 1, resampled_profile + 2, lower_rank - 1) -
	                    qtile(sample_count - 1, resampled_profile + 2, upper_rank - 1);
	terrain_roughness /= 1.0 - 0.8 * exp(-(dist_end - dist_start) / 50.0e3);
	delete[] resampled_profile;
	return terrain_roughness;
}

/*
 * qlrpfl - Terrain profile setup for ITM model
 *          (Quasi-smooth earth Long-Range Propagation profile setup)
 *
 * Computes from the terrain profile all quantities needed by lrprop():
 *   - horizons (hzns), roughness (d1thx)
 *   - effective heights (he[]), horizon distances (dl[]), angles (the[])
 * Then initializes lrprop with d=0.
 *
 * Parameters:
 *   terrain_profile (former: pfl[])    : terrain profile
 *   climate_index (former: klimx)      : radio climate index
 *   variability_mode (former: mdvarx)  : variability mode
 *   prop, propa, propv                 : parameter structures (modified)
 */
void qlrpfl(double terrain_profile[], int climate_index, int variability_mode,
            prop_type & prop, propa_type & propa, propv_type & propv)
/* former: pfl, klimx, mdvarx */
{
	int num_points, j;           /* former: np, j */
	double window_limits[2];     /* former: xl — fitting window limits */
	double interp_q, tx_ground, rx_ground, scale_factor; /* former: q, za, zb, temp */

	prop.dist = terrain_profile[0] * terrain_profile[1];
	num_points = (int)terrain_profile[0];
	hzns(terrain_profile, prop);

	for (j = 0; j < 2; j++)
		window_limits[j] = MIN(15.0 * prop.hg[j], 0.1 * prop.dl[j]);

	window_limits[1] = prop.dist - window_limits[1];
	prop.dh = d1thx(terrain_profile, window_limits[0], window_limits[1]);

	if (prop.dl[0] + prop.dl[1] > 1.5 * prop.dist) {
		/* LOS path or multiple close obstructions */
		z1sq1(terrain_profile, window_limits[0], window_limits[1], tx_ground, rx_ground);
		prop.he[0] = prop.hg[0] + FORTRAN_DIM(terrain_profile[2], tx_ground);
		prop.he[1] = prop.hg[1] + FORTRAN_DIM(terrain_profile[num_points + 2], rx_ground);

		for (j = 0; j < 2; j++)
			prop.dl[j] =
			    sqrt(2.0 * prop.he[j] / prop.gme) * exp(-0.07 *
								    sqrt(prop.dh / MAX(prop.he[j], 5.0)));

		interp_q = prop.dl[0] + prop.dl[1];

		if (interp_q <= prop.dist) { /* rounded horizon or two obstacles */
			scale_factor = prop.dist / interp_q;
			interp_q     = scale_factor * scale_factor;

			for (j = 0; j < 2; j++) {
				prop.he[j] *= interp_q; /* adjusted effective height */
				prop.dl[j] =
				    sqrt(2.0 * prop.he[j] / prop.gme) *
				    exp(-0.07 * sqrt(prop.dh / MAX(prop.he[j], 5.0)));
			}
		}

		for (j = 0; j < 2; j++) { /* empirical grazing angle adjustment */
			interp_q = sqrt(2.0 * prop.he[j] / prop.gme);
			prop.the[j] =
			    (0.65 * prop.dh * (interp_q / prop.dl[j] - 1.0) -
			     2.0 * prop.he[j]) / interp_q;
		}
	}

	else {
		/* Separate TX and RX fit on the slope leading to the obstacle */
		z1sq1(terrain_profile, window_limits[0], 0.9 * prop.dl[0], tx_ground, interp_q);
		z1sq1(terrain_profile, prop.dist - 0.9 * prop.dl[1], window_limits[1], interp_q, rx_ground);
		prop.he[0] = prop.hg[0] + FORTRAN_DIM(terrain_profile[2], tx_ground);
		prop.he[1] = prop.hg[1] + FORTRAN_DIM(terrain_profile[num_points + 2], rx_ground);
	}

	prop.mdp  = -1;
	propv.lvar = MAX(propv.lvar, 3);

	if (variability_mode >= 0) {
		propv.mdvar = variability_mode;
		propv.lvar  = MAX(propv.lvar, 4);
	}

	if (climate_index > 0) {
		propv.klim = climate_index;
		propv.lvar = 5;
	}

	lrprop(0.0, prop, propa);
}

/*
 * qlrpfl2 - Terrain profile setup for ITWOM model
 *           (Quasi-smooth earth Long-Range Propagation profile setup, version 2)
 *
 * ITWOM version of qlrpfl(). Uses hzns2(), d1thx2(), z1sq2() and lrprop2().
 * In ITWOM mode (no np or spacing > 150 m), also computes the receiver
 * approach angles (prop.thera, prop.thenr) for saalos.
 *
 * Parameters: same as qlrpfl()
 */
void qlrpfl2(double terrain_profile[], int climate_index, int variability_mode,
             prop_type & prop, propa_type & propa, propv_type & propv)
/* former: pfl, klimx, mdvarx */
{
	int num_points, j;           /* former: np, j */
	double window_limits[2];     /* former: xl */
	double total_horizon_dist;   /* former: dlb — sum of TX+RX horizon distances */
	double interp_q = 1.0;       /* former: q */
	double tx_ground, rx_ground; /* former: za, zb */
	double scale_factor;         /* former: temp */
	double approach_dist;        /* former: rad — distance for RX approach angle */
	double approach_elev1, approach_elev2; /* former: rae1, rae2 — elevations for angle computation */

	prop.dist  = terrain_profile[0] * terrain_profile[1];
	num_points = (int)terrain_profile[0];
	hzns2(terrain_profile, prop, propa);
	total_horizon_dist = prop.dl[0] + prop.dl[1];
	prop.rch[0] = prop.hg[0] + terrain_profile[2];
	prop.rch[1] = prop.hg[1] + terrain_profile[num_points + 2];

	for (j = 0; j < 2; j++)
		window_limits[j] = MIN(15.0 * prop.hg[j], 0.1 * prop.dl[j]);

	window_limits[1] = prop.dist - window_limits[1];
	prop.dh = d1thx2(terrain_profile, window_limits[0], window_limits[1], propa);

	if ((num_points < 1) || (terrain_profile[1] > 150.0)) {
		/* Trans-horizon or large grid spacing mode */
		if (total_horizon_dist < 1.5 * prop.dist) {
			/* Trans-horizon: diffraction over mutual horizon or obstruction(s) */
			z1sq2(terrain_profile, window_limits[0], 0.9 * prop.dl[0], tx_ground, interp_q);
			z1sq2(terrain_profile, prop.dist - 0.9 * prop.dl[1], window_limits[1], interp_q, rx_ground);
			prop.he[0] = prop.hg[0] + FORTRAN_DIM(terrain_profile[2], tx_ground);
			prop.he[1] = prop.hg[1] + FORTRAN_DIM(terrain_profile[num_points + 2], rx_ground);
		} else {
			/* LOS path */
			z1sq2(terrain_profile, window_limits[0], window_limits[1], tx_ground, rx_ground);
			prop.he[0] = prop.hg[0] + FORTRAN_DIM(terrain_profile[2], tx_ground);
			prop.he[1] = prop.hg[1] + FORTRAN_DIM(terrain_profile[num_points + 2], rx_ground);

			for (j = 0; j < 2; j++)
				prop.dl[j] =
				    sqrt(2.0 * prop.he[j] / prop.gme) *
				    exp(-0.07 * sqrt(prop.dh / MAX(prop.he[j], 5.0)));

			/* Obstruction correction — ITM FORTRAN legacy, rarely active */
			if ((prop.dl[0] + prop.dl[1]) <= prop.dist) {
				scale_factor = prop.dist / (prop.dl[0] + prop.dl[1]);
				interp_q     = scale_factor * scale_factor;
			}

			for (j = 0; j < 2; j++) {
				prop.he[j] *= interp_q;
				prop.dl[j] =
				    sqrt(2.0 * prop.he[j] / prop.gme) *
				    exp(-0.07 * sqrt(prop.dh / MAX(prop.he[j], 5.0)));
			}

			/* Recompute grazing angles */
			for (j = 0; j < 2; j++) {
				interp_q = sqrt(2.0 * prop.he[j] / prop.gme);
				prop.the[j] =
				    (0.65 * prop.dh * (interp_q / prop.dl[j] - 1.0) -
				     2.0 * prop.he[j]) / interp_q;
			}
		}
	}

	else { /* ITWOM mode: compute receiver approach angles for saalos */
		prop.he[0] = prop.hg[0] + (terrain_profile[2]);
		prop.he[1] = prop.hg[1] + (terrain_profile[num_points + 2]);

		approach_dist = (prop.dist - 500.0);

		if (prop.dist > 550.0) {
			/* Approach slope angle over the last 500 meters */
			z1sq2(terrain_profile, approach_dist, prop.dist, approach_elev1, approach_elev2);
		} else {
			approach_elev1 = 0.0;
			approach_elev2 = 0.0;
		}

		prop.thera = atan(abs(approach_elev2 - approach_elev1) / prop.dist);

		if (approach_elev2 < approach_elev1) {
			prop.thera = -prop.thera; /* downward slope toward receiver */
		}

		/* Angle of the last profile segment on the RX side */
		prop.thenr =
		    atan(MAX(0.0, (terrain_profile[num_points + 2] - terrain_profile[num_points + 1])) / terrain_profile[1]);
	}

	prop.mdp  = -1;
	propv.lvar = MAX(propv.lvar, 3);

	if (variability_mode >= 0) {
		propv.mdvar = variability_mode;
		propv.lvar  = MAX(propv.lvar, 4);
	}

	if (climate_index > 0) {
		propv.klim = climate_index;
		propv.lvar = 5;
	}

	lrprop2(0.0, prop, propa);
}

//***************************************************************************************
//* Point-To-Point Mode Calculations
//***************************************************************************************

/*
 * point_to_point_ITM - Point-to-point propagation computation (classic ITM model)
 *
 * Computes the propagation loss (dbloss, dB) between two points using the
 * original ITM (Longley-Rice) model with:
 *   - propagation over irregular terrain (profile elev[])
 *   - radio-physical parameters (frequency, climate, polarization, ground)
 *   - confidence (conf) and reliability (rel) statistics
 *
 * Parameters:
 *   tx_height_m (former: tht_m)                    : TX antenna height (m AGL)
 *   rx_height_m (former: rht_m)                    : RX antenna height (m AGL)
 *   dielectric_const (former: eps_dielect)         : relative permittivity of ground
 *   conductivity (former: sgm_conductivity)        : ground conductivity (S/m)
 *   surface_refractivity (former: eno_ns_surfref)  : surface refractivity (N-units)
 *   freq_mhz (former: frq_mhz)                    : frequency (MHz)
 *   radio_climate                                   : radio climate index (1-7)
 *   polarization (former: pol)                     : polarization (0=H, 1=V)
 *   confidence (former: conf)                      : confidence level (0.01-0.99)
 *   reliability (former: rel)                      : reliability level (0.01-0.99)
 *   dbloss                                          : computed propagation loss (dB, output)
 *   mode                                            : detected propagation mode (output)
 *   errnum                                          : error code (output, 0=OK)
 */
void point_to_point_ITM(double tx_height_m, double rx_height_m,
                        double dielectric_const, double conductivity,
                        double surface_refractivity, double freq_mhz,
                        int radio_climate, int polarization,
                        double confidence, double reliability,
                        double &dbloss, PropagationMode &mode, int &errnum)
/* former: tht_m, rht_m, eps_dielect, sgm_conductivity, eno_ns_surfref, frq_mhz, pol, conf, rel */

/******************************************************************************

Note that point_to_point has become point_to_point_ITM for use as the old ITM

	pol:
		0-Horizontal, 1-Vertical

	radio_climate:
		1-Equatorial, 2-Continental Subtropical,
		3-Maritime Tropical, 4-Desert, 5-Continental Temperate,
		6-Maritime Temperate, Over Land, 7-Maritime Temperate,
		Over Sea

	conf, rel: .01 to .99

	elev[]: [num points - 1], [delta dist(meters)],
	        [height(meters) point 1], ..., [height(meters) point n]

	errnum: 0- No Error.
		1- Warning: Some parameters are nearly out of range.
		            Results should be used with caution.
		2- Note: Default parameters have been substituted for
		         impossible ones.
		3- Warning: A combination of parameters is out of range.
			    Results are probably invalid.
		Other-  Warning: Some parameters are out of range.
			Results are probably invalid.

*****************************************************************************/
{
	if (debug) cnt_point_to_point_ITM++;
	prop_type prop;
	propv_type propv;
	propa_type propa;
	double mean_elevation = 0;  /* former: zsys */
	double confidence_z, reliability_z; /* former: zc, zr */
	double refractivity, interp_q;      /* former: eno, q */
	long profile_start, profile_end, i, num_profile_pts; /* former: ja, jb, i, np */
	/* double dkm, xkm; */
	double free_space_loss; /* former: fs */

	prop.hg[0] = tx_height_m;
	prop.hg[1] = rx_height_m;
	propv.klim = radio_climate;
	prop.kwx   = 0;
	propv.lvar = 5;
	prop.mdp   = -1;
	confidence_z = qerfi(confidence);
	reliability_z = qerfi(reliability);
	num_profile_pts = (long)elev[0];
	/* dkm=(elev[1]*elev[0])/1000.0; */
	/* xkm=elev[1]/1000.0; */
	refractivity = surface_refractivity;

	/* Compute mean profile elevation to correct refractivity */
	profile_start = (long)(3.0 + 0.1 * elev[0]);
	profile_end   = num_profile_pts - profile_start + 6;

	for (i = profile_start - 1; i < profile_end; ++i)
		mean_elevation += elev[i];

	mean_elevation /= (profile_end - profile_start + 1);
	interp_q = refractivity;

	propv.mdvar = 12;
	qlrps(freq_mhz, mean_elevation, interp_q, polarization,
	      dielectric_const, conductivity, prop);
	qlrpfl(elev, propv.klim, propv.mdvar, prop, propa, propv);
	free_space_loss = 32.45 + 20.0 * log10(freq_mhz) + 20.0 * log10(prop.dist / 1000.0);
	interp_q = prop.dist - propa.dla;

	/* Determine propagation mode */
	if (int(interp_q) < 0.0)
		mode = PROP_MODE_LOS;
	else if (int(interp_q) == 0.0)
		mode = PROP_MODE_1_HRZN;
	else if (int(interp_q) > 0.0)
		mode = PROP_MODE_2_HRZN;

	if (mode != PROP_MODE_LOS) {
		if (prop.dist <= propa.dlsa || prop.dist <= propa.dx)
			mode |= PROP_MODE_DIFFRACTION;
		else if (prop.dist > propa.dx)
			mode |= PROP_MODE_TROPOSCATTER;
	}

	dbloss = avar(reliability_z, 0.0, confidence_z, prop, propv) + free_space_loss;
	errnum = prop.kwx;
}

/*
 * point_to_point - Point-to-point propagation computation (ITWOM model)
 *
 * ITWOM version of point_to_point_ITM(). Uses qlrpfl2() and lrprop2().
 * Accounts for vegetation cover (cch, encc, cd) and supports circular
 * polarization (pol=2). Vegetation parameters are fixed to predefined
 * values compatible with ITU-R P.1546-2 / FCC.
 *
 * Parameters: same as point_to_point_ITM() except:
 *   polarization : 0=H, 1=V, 2=circular
 */
void point_to_point(double tx_height_m, double rx_height_m,
                    double dielectric_const, double conductivity,
                    double surface_refractivity, double freq_mhz,
                    int radio_climate, int polarization,
                    double confidence, double reliability,
                    double &dbloss, PropagationMode &mode, int &errnum)
/* former: tht_m, rht_m, eps_dielect, sgm_conductivity, eno_ns_surfref, frq_mhz, pol, conf, rel */

/******************************************************************************

	Note that point_to_point_two has become point_to_point
	for drop-in interface to splat.cpp.
	The new variable inputs,
	double enc_ncc_clcref,
	double clutter_height,
	double clutter_density,
	double delta_h_diff, and
	int mode_var)
	have been given fixed values below.

	pol:
		0-Horizontal, 1-Vertical, 2-Circular

	radio_climate:
		1-Equatorial, 2-Continental Subtropical,
		3-Maritime Tropical, 4-Desert, 5-Continental Temperate,
		6-Maritime Temperate, Over Land, 7-Maritime Temperate,
		Over Sea

	conf, rel: .01 to .99

	elev[]: [num points - 1], [delta dist(meters)],
	        [height(meters) point 1], ..., [height(meters) point n]

	clutter_height  	25.2 meters for compatibility with ITU-R P.1546-2.

	clutter_density 	1.0 for compatibility with ITU-R P.1546-2.

	delta_h_diff		optional delta h for beyond line of sight. 90 m. average.
				setting to 0.0 will default to use of original internal
				use of delta-h for beyond line-of-sight range.

	mode_var		set to 12; or to 1 for FCC ILLR;  see documentation

	enc_ncc_clcref 		clutter refractivity; 1000 N-units to match ITU-R P.1546-2

	eno=eno_ns_surfref	atmospheric refractivity at sea level; 301 N-units nominal
				(ranges from 250 for dry, hot day to 450 on hot, humid day]
				(stabilizes near 301 in cold, clear weather)

	errnum: 0- No Error.
		1- Warning: Some parameters are nearly out of range.
		            Results should be used with caution.
		2- Note: Default parameters have been substituted for
		         impossible ones.
		3- Warning: A combination of parameters is out of range.
			    Results are probably invalid.
		Other-  Warning: Some parameters are out of range.
			Results are probably invalid.

*****************************************************************************/
{
	if (debug) cnt_point_to_point++;
	prop_type prop;
	propv_type propv;
	propa_type propa;

	double mean_elevation = 0;  /* former: zsys */
	double confidence_z, reliability_z; /* former: zc, zr */
	double refractivity, interp_q;      /* former: eno, q */
	long profile_start, profile_end, i, num_profile_pts; /* former: ja, jb, i, np */
	/* double dkm, xkm; */
	double slant_dist, free_space_loss; /* former: tpd, fs */

	prop.hg[0] = tx_height_m;
	prop.hg[1] = rx_height_m;
	propv.klim = radio_climate;
	prop.kwx   = 0;
	propv.lvar = 5;
	prop.mdp   = -1;
	prop.ptx   = polarization;
	prop.thera = 0.0;
	prop.thenr = 0.0;
	confidence_z = qerfi(confidence);
	reliability_z = qerfi(reliability);
	num_profile_pts = (long)elev[0];
	/* dkm=(elev[1]*elev[0])/1000.0; */
	/* xkm=elev[1]/1000.0; */
	refractivity = surface_refractivity;

	/* Predefined values for the base version without additional inputs */
	prop.encc = 1000.00; /* vegetation canopy refractivity (N-units) */
	prop.cch  = 22.5;    /* vegetation canopy height (m) — ILLR calibration;
	                        use 25.3 for ITU-P1546-2 calibration */
	prop.cd   = 1.00;    /* vegetation canopy density */
	int variability_mode_preset = 1; /* former: mode_var — mode 1 for FCC compatibility */
	prop.dhd  = 0.0;     /* delta-h for diffraction (preset) */

	profile_start = (long)(3.0 + 0.1 * elev[0]);
	profile_end   = num_profile_pts - profile_start + 6;

	for (i = profile_start - 1; i < profile_end; ++i)
		mean_elevation += elev[i];

	mean_elevation /= (profile_end - profile_start + 1);
	interp_q = refractivity;

	propv.mdvar = variability_mode_preset;
	qlrps(freq_mhz, mean_elevation, interp_q, polarization,
	      dielectric_const, conductivity, prop);
	qlrpfl2(elev, propv.klim, propv.mdvar, prop, propa, propv);

	/* Slant distance for free-space loss computation */
	slant_dist =
	    sqrt((prop.he[0] - prop.he[1]) * (prop.he[0] - prop.he[1]) +
		 (prop.dist) * (prop.dist));
	free_space_loss = 32.45 + 20.0 * log10(freq_mhz) + 20.0 * log10(slant_dist / 1000.0);
	interp_q = prop.dist - propa.dla;

	/* Determine propagation mode */
	if (int(interp_q) < 0.0)
		mode = PROP_MODE_LOS;
	else if (int(interp_q) == 0.0)
		mode = PROP_MODE_1_HRZN;
	else if (int(interp_q) > 0.0)
		mode = PROP_MODE_2_HRZN;

	if (mode != PROP_MODE_LOS) {
		if (prop.dist <= propa.dlsa || prop.dist <= propa.dx) {
			if (int(prop.dl[1]) == 0.0)
				mode |= PROP_MODE_PEAK;
			else
				mode |= PROP_MODE_DIFFRACTION;
		} else if (prop.dist > propa.dx) {
			mode |= PROP_MODE_TROPOSCATTER;
		}
	}

	dbloss = avar(reliability_z, 0.0, confidence_z, prop, propv) + free_space_loss;
	errnum = prop.kwx;
}
