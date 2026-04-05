#include "geo.hh"

double earthRadius(double lat)
{
    // Convert latitude to rad
    double lat_rad = lat * DEG2RAD;

    double An = WGS84_a * WGS84_a * cos(lat_rad);
    double Bn = WGS84_b * WGS84_b * sin(lat_rad);
    double Ad = WGS84_a * cos(lat_rad);
    double Bd = WGS84_b * sin(lat_rad);

    return double(sqrt( (An*An + Bn*Bn) / (Ad*Ad + Bd*Bd) ) / (double)1000);
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
    
    result.lower_right = { latMin * RAD2DEG, lonMin * RAD2DEG };
    result.upper_left = { latMax * RAD2DEG, lonMax * RAD2DEG };

    return result;
}