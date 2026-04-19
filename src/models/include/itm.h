#include <complex>
#include <math.h>
#include <algorithm>
#include <stdlib.h>
#include <vector>

using namespace std;

// Export the DLL functions as "C" and not C++
//#define DLLEXPORT extern "C" __declspec(dllexport)
#define DLLEXPORT 
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define DIM(x, y) (((x) > (y)) ? (x - y) : (0))
#ifndef PI
#define PI                                      3.1415926535897932384
#endif
#define SQRT2                                   sqrt(2)
#define a_0__meter                              6370e3
#define a_9000__meter                           9000e3
#define THIRD                                   1.0 / 3.0

#define MODE__P2P                               0
#define MODE__AREA                              1


/////////////////////////////
// Data Structures

struct IntermediateValues
{
    float theta_hzn[2];        // Terminal horizon angles
    float d_hzn__meter[2];     // Terminal horizon distances, in meters
    float h_e__meter[2];       // Terminal effective heights, in meters
    float N_s;                 // Surface refractivity, in N-Units
    float delta_h__meter;      // Terrain irregularity parameter, in meters
    float A_ref__db;           // Reference attenuation, in dB
    float A_fs__db;            // Free space basic transmission loss, in dB
    float d__km;               // Path distance, in km
    int mode;                   // Mode of propagation value
};

/////////////////////////////
// Main ITM Functions

DLLEXPORT int ITM_P2P_TLS(const float h_tx__meter, const float h_rx__meter, const float pfl[], const int climate, const float N_0, const float f__mhz,
    const int pol, const float epsilon, const float sigma, const int mdvar, const float time, const float location, const float situation,
    float *A__db, long *warnings);
DLLEXPORT int ITM_P2P_TLS_Ex(const float h_tx__meter, const float h_rx__meter, const float pfl[], const int climate, const float N_0, const float f__mhz,
    const int pol, const float epsilon, const float sigma, const int mdvar, const float time, const float location, const float situation,
    float *A__db, long *warnings, IntermediateValues *interValues);
DLLEXPORT int ITM_P2P_CR(const float h_tx__meter, const float h_rx__meter, const float pfl[], const int climate, const float N_0, const float f__mhz,
    const int pol, const float epsilon, const float sigma, const int mdvar, const float confidence, const float reliability,
    float *A__db, long *warnings);
DLLEXPORT int ITM_P2P_CR_Ex(const float h_tx__meter, const float h_rx__meter, const float pfl[], const int climate, const float N_0, const float f__mhz,
    const int pol, const float epsilon, const float sigma, const int mdvar, const float confidence, const float reliability,
    float *A__db, long *warnings, IntermediateValues *interValues);
DLLEXPORT int ITM_AREA_TLS(const float h_tx__meter, const float h_rx__meter, const int tx_site_criteria, const int rx_site_criteria, const float d__km,
    const float delta_h__meter, const int climate, const float N_0, float f__mhz, const int pol, const float epsilon, const float sigma,
    const int mdvar, const float time, const float location, const float situation, float *A__db, long *warnings);
DLLEXPORT int ITM_AREA_TLS_Ex(const float h_tx__meter, const float h_rx__meter, const int tx_site_criteria, const int rx_site_criteria, const float d__km,
    const float delta_h__meter, const int climate, const float N_0, const float f__mhz, const int pol, const float epsilon, const float sigma,
    const int mdvar, const float time, const float location, const float situation, float *A__db, long *warnings, IntermediateValues *interValues);
DLLEXPORT int ITM_AREA_CR(const float h_tx__meter, const float h_rx__meter, const int tx_site_criteria, const int rx_site_criteria, const float d__km,
    const float delta_h__meter, const int climate, const float N_0, const float f__mhz, const int pol, const float epsilon, const float sigma,
    const int mdvar, const float confidence, const float reliability, float *A__db, long *warnings);
DLLEXPORT int ITM_AREA_CR_Ex(const float h_tx__meter, const float h_rx__meter, const int tx_site_criteria, const int rx_site_criteria, const float d__km,
    const float delta_h__meter, const int climate, const float N_0, const float f__mhz, const int pol, const float epsilon, const float sigma,
    const int mdvar, const float confidence, const float reliability, float *A__db, long *warnings, IntermediateValues *interValues);

/////////////////////////////
// ITM Helper Functions

DLLEXPORT float ComputeDeltaH(const float pfl[], const float d_start__meter, const float d_end__meter);
DLLEXPORT float DiffractionLoss(const float d__meter, const float d_hzn__meter[2], const float h_e__meter[2], const complex<float> Z_g,
    const float a_e__meter, const float delta_h__meter, const float h__meter[2], const int mode, const float theta_los, const float d_sML__meter, const float f__mhz);
DLLEXPORT float FFunction(const float td);
DLLEXPORT void FindHorizons(const float pfl[], const float a_e__meter, const float h__meter[2], float theta_hzn[2], float d_hzn__meter[2]);
DLLEXPORT float FreeSpaceLoss(const float d__meter, const float f__mhz);
DLLEXPORT float FresnelIntegral(const float v2);
DLLEXPORT float H0Function(const float r, float eta_s);
DLLEXPORT float HeightFunction(const float x__km, const float K);
DLLEXPORT void InitializeArea(const int site_criteria[2], const float gamma_e, const float delta_h__meter,
    const float h__meter[2], float h_e__meter[2], float d_hzn__meter[2], float theta_hzn[2]);
DLLEXPORT void InitializePointToPoint(const float f__mhz, const float h_sys__meter, const float N_0, const int pol, const float epsilon, 
    const float sigma, complex<float> *Z_g, float *gamma_e, float *N_s);
DLLEXPORT float InverseComplementaryCumulativeDistributionFunction(const float q);
DLLEXPORT float KnifeEdgeDiffraction(const float d__meter, const float f__mhz, const float a_e__meter, const float theta_los, const float d_hzn__meter[2]);
DLLEXPORT void LinearLeastSquaresFit(const float pfl[], const float d_start, const float d_end, float *fit_y1, float *fit_y2);
DLLEXPORT float LineOfSightLoss(const float d__meter, const float h_e__meter[2], const complex<float> Z_g, const float delta_h__meter,
    const float M_d, const float A_d0, const float d_sML__meter, const float f__mhz);
DLLEXPORT int LongleyRice(const float theta_hzn[2], const float f__mhz, const complex<float> Z_g, const float d_hzn__meter[2], const float h_e__meter[2], 
    const float gamma_e, const float N_s, const float delta_h__meter, const float h__meter[2], const float d__meter, const int mode, float *A_ref__db, 
    long *warnings, int *propmode);
DLLEXPORT void QuickPfl(const float pfl[], const float gamma_e, const float h__meter[2], float theta_hzn[2], float d_hzn__meter[2], 
    float h_e__meter[2], float *delta_h__meter, float *d__meter);
DLLEXPORT float SigmaHFunction(const float delta_h__meter);
DLLEXPORT float SmoothEarthDiffraction(const float d__meter, const float f__mhz, const float a_e__meter, const float theta_los, 
    const float d_hzn__meter[2], const float h_e__meter[2], const complex<float> Z_g);
DLLEXPORT float TerrainRoughness(const float d__meter, const float delta_h__meter);
DLLEXPORT float TroposcatterLoss(const float d__meter, const float theta_hzn[2], const float d_hzn__meter[2], const float h_e__meter[2], 
    const float a_e__meter, const float N_s, const float f__mhz, const float theta_los, float *h0);
DLLEXPORT int ValidateInputs(const float h_tx__meter, const float h_rx__meter, const int climate, const float time,
    const float location, const float situation, const float N_0, const float f__mhz, const int pol,
    const float epsilon, const float sigma, const int mdvar, long *warnings);
DLLEXPORT float Variability(const float time, const float location, const float situation, const float h_e__meter[2], const float delta_h__meter,
    const float f__mhz, const float d__meter, const float A_ref__db, const int climate, const int mdvar, long *warnings);
