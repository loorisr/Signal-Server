#include <algorithm>
#include "geo.hh"

double arccos(double x, double y)
{
    if (y == 0.0) return 0.0;
    double result = acos(x / y);
    return y < 0.0 ? PI + result : result;
}

double LonDiff(double lon1, double lon2)
{
    double diff = lon1 - lon2;
    if (diff <= -180.0) return diff + 360.0;
    if (diff >= 180.0)  return diff - 360.0;
    return diff;
}

double Distance(struct site site1, struct site site2)
{
    double lat1 = site1.lat * DEG2RAD;
    double lon1 = site1.lon * DEG2RAD;
    double lat2 = site2.lat * DEG2RAD;
    double lon2 = site2.lon * DEG2RAD;

    double dot = sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2);
    dot = std::clamp(dot, -1.0, 1.0);  /* guard against floating-point overshoot */

    return EARTHRADIUS * acos(dot);
}

double Azimuth(struct site source, struct site destination)
{
    double dest_lat, dest_lon, src_lat, src_lon,
        beta, azimuth, diff, num, den, fraction;

    dest_lat = destination.lat * DEG2RAD;
    dest_lon = destination.lon * DEG2RAD;

    src_lat = source.lat * DEG2RAD;
    src_lon = source.lon * DEG2RAD;

    beta =
        acos(sin(src_lat) * sin(dest_lat) +
         cos(src_lat) * cos(dest_lat) * cos(src_lon - dest_lon));

    num = sin(dest_lat) - (sin(src_lat) * cos(beta));
    den = cos(src_lat) * sin(beta);

    if (den == 0.0)
        return 0.0;  /* source at pole or source == destination */

    fraction = std::clamp(num / den, -1.0, 1.0);
    azimuth = acos(fraction);

    diff = dest_lon - src_lon;
    if (diff <= -PI) diff += TWOPI;
    if (diff >=  PI) diff -= TWOPI;

    if (diff > 0.0)
        azimuth = TWOPI - azimuth;

    return (azimuth * RAD2DEG);
}

double earthRadius(double lat)
{
    // Convert latitude to rad
    double lat_rad = lat * DEG2RAD;

    double An = WGS84_a * WGS84_a * cos(lat_rad);
    double Bn = WGS84_b * WGS84_b * sin(lat_rad);
    double Ad = WGS84_a * cos(lat_rad);
    double Bd = WGS84_b * sin(lat_rad);

    return sqrt((An*An + Bn*Bn) / (Ad*Ad + Bd*Bd));
}

coord getPointAtDistance(const coord center, double distance, double bearing_rad)
{
    coord endCoords;

    // 1. Pre-calculate radians and Earth radius ratio
    const double start_lat_rad = center.lat * DEG2RAD;
    const double start_lon_rad = center.lon * DEG2RAD;
    const double dR            = distance / earthRadius(center.lat);

    // 2. Pre-calculate common trig values to avoid redundant calls
    const double sin_lat = sin(start_lat_rad);
    const double cos_lat = cos(start_lat_rad);
    const double sin_dR  = sin(dR);
    const double cos_dR  = cos(dR);
    const double sin_brg = sin(bearing_rad);
    const double cos_brg = cos(bearing_rad);

    // 3. Calculate Latitude
    const double sin_end_lat = sin_lat * cos_dR + cos_lat * sin_dR * cos_brg;
    const double end_lat_rad = asin(sin_end_lat);
    
    // 4. Calculate Longitude
    // Optimization: Reuse sin_end_lat to avoid another sin() call
    const double end_lon_rad = start_lon_rad + atan2(
        sin_brg * sin_dR * cos_lat,
        cos_dR - sin_lat * sin_end_lat
    );

    endCoords.lat = end_lat_rad * RAD2DEG;
    endCoords.lon = end_lon_rad * RAD2DEG;

    return endCoords;
}

bbox getCircularBoundingBox(coord center, double radius)
{
    // Result bbox
    bbox result;

    // Convert input degrees to rads
    double lat_rad = center.lat * DEG2RAD;
    double lon_rad = center.lon * DEG2RAD;

    // Get earth's radius at the specified latitude (km)
    double e_rad = earthRadius(center.lat);

    // Get parallel radius at latitude (km)
    double p_rad = e_rad * cos(lat_rad);

    // Calculate bounds (radians)
    double latMin = lat_rad - (radius / e_rad);
    double latMax = lat_rad + (radius / e_rad);
    double lonMin = lon_rad - (radius / p_rad);
    double lonMax = lon_rad + (radius / p_rad);
    
    result.lower_left = { latMin * RAD2DEG, lonMin * RAD2DEG };
    result.upper_right = { latMax * RAD2DEG, lonMax * RAD2DEG };

    return result;
}