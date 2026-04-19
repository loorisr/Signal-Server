#include "../include/itm.h"

/*=============================================================================
 |
 |  Description:  Free space basic transmission loss equation
 |
 |        Input:  d__meter       - Path distance, in meters
 |                f__mhz         - Frequency, in MHz
 |
 |      Outputs:  [None]
 |
 |      Returns:  A_fs__db       - Free space basic transmission loss, in dB
 |
 *===========================================================================*/
float FreeSpaceLoss(const float d__meter, const float f__mhz)
{
    return 32.45 + 20.0 * log10(f__mhz) + 20.0 * log10(d__meter / 1000.0);
}