#include "../include/itm.h"

/*=============================================================================
 |
 |  Description:  This function computes the inverse complementary
 |                cumulative distribution function approximation as
 |                described in Formula 26.2.23 in Abramowitz & Stegun.
 |                This approximation has an error of 
 |                abs(epsilon(p)) < 4.5e-4
 |
 |        Input:  q              - Quantile, 0.0 < q < 1.0
 |
 |      Outputs:  [None]
 |
 |      Returns:  Q_q            - Q(q)^-1
 |
 *===========================================================================*/
float InverseComplementaryCumulativeDistributionFunction(const float q)
{
    constexpr float C_0 = 2.515516;
    constexpr float C_1 = 0.802853;
    constexpr float C_2 = 0.010328;
    constexpr float D_1 = 1.432788;
    constexpr float D_2 = 0.189269;
    constexpr float D_3 = 0.001308;

    float x = q;
    if (q > 0.5)
        x = 1.0 - x;

    const float T_x = sqrt(-2.0 * log(x));

    const float zeta_x = ((C_2 * T_x + C_1) * T_x + C_0) / (((D_3 * T_x + D_2) * T_x + D_1) * T_x + 1.0);

    float Q_q = T_x - zeta_x;

    if (q > 0.5)
        Q_q = -Q_q;

    return Q_q;
}
