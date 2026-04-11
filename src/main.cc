/****************************************************************************\
*  Signal Server: Radio propagation simulator by Alex Farrant QCVS, 2E0TDW   *
******************************************************************************
*    SPLAT! Project started in 1997 by John A. Magliacane, KD2BD             *
******************************************************************************
*         Please consult the SPLAT! documentation for a complete list of     *
*         individuals who have contributed to this project.                  *
******************************************************************************
*                                                                            *
*  This program is free software; you can redistribute it and/or modify it   *
*  under the terms of the GNU General Public License as published by the     *
*  Free Software Foundation; either version 2 of the License or any later    *
*  version.                                                                  *
*                                                                            *
*  This program is distributed in the hope that it will useful, but WITHOUT  *
*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or     *
*  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License     *
*  for more details.                                                         *
*                                                                            *
\****************************************************************************/

/* GDAL must come before any header that defines MAX */
#include <gdal.h>
#include <cpl_conv.h>
#include <ogr_srs_api.h>

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <errno.h>
#include <limits.h>

#include "main.hh"
#include "common.hh"
#include "cmdline.hh"
#include "inputs.hh"
#include "outputs.hh"
#include "models/itwom3.0.hh"
#include "models/los.hh"
#include "models/pel.hh"

#include <algorithm>
#include <chrono>
#include <thread>

#include <spdlog/spdlog.h>

int ppd = 1200;
int MAX_DISTANCE_DEGRES = 3; // max distance : 3° so around 300 km
int ARRAYSIZE = (MAX_DISTANCE_DEGRES * ppd) + 10;

char DEM_path[255];
std::string color_palette = "heat";
std::string output_filename;

double max_range = 0.0,
    fzone_clearance = 0.6, clutter, tercon, terdic,
    north, east, south, west, dBm, loss, field_strength,
    min_north = 90, max_north = -90, min_lon = 180.0, max_lon = -180.0,
    rxGain=0, antenna_rotation,
    antenna_downtilt, antenna_dt_direction;

int max_elevation = -32768, min_elevation = 32768,
    contour_threshold, debug = 0;

std::atomic<int> cnt_point_to_point_ITM{0};
std::atomic<int> cnt_point_to_point{0};
std::atomic<int> cnt_computeLoss{0};
std::atomic<int> cnt_PlotPropPath{0};
std::atomic<int> cnt_ReadPath{0};
std::atomic<int> cnt_PlotPropagation{0};

bool got_elevation_pattern = false, got_azimuth_pattern = false, dbm = false;
bool geotiff = false;
bool ngs = false;
bool cropping = true;
int knifeedge = 0;
int pmenv = 1;
int number_threads = 4;
double altitudeLR = 1;
PropModel prop_model = ITM_LR;
bool ppa = false;
int normalise = 0;
char mapfile[255] = "";

__thread double *elev;
__thread struct path path;
site tx_site;
site rx_site;
double         **dem_data   = nullptr;
int **dem_signal = nullptr;
int dem_min_lat   = 0;
int dem_min_lon   = 0;
int dem_width_px  = 0;
int dem_height_px = 0;

struct LR LR;

int main(int argc, char *argv[])
{
    auto start_time = std::chrono::steady_clock::now();

    parse_cmdline(argc, argv);

    if (ppa) {

    } else {
        bbox geo_bounds = getCircularBoundingBox({tx_site.lat, tx_site.lon}, max_range);
    }
    


    // Get latitude in radians
    double tx_lat_rad = tx_site.lat * DEG2RAD;

    // Find the distance in lat and lon per degree using the above referenced formulas
    double m_per_deg_lon = (111412.84 * cos(tx_lat_rad)) - (93.5 * cos(3 * tx_lat_rad)) + (0.118 * cos(5 * tx_lat_rad));
    double m_per_deg_lat = 111132.92 - (559.82 * cos(2 * tx_lat_rad)) + (1.175 * cos(4 * tx_lat_rad)) - (0.0023 * cos(6 * tx_lat_rad));

    // Calculate angular distance from the above numbers
    double dist_deg_lon = max_range / m_per_deg_lon;
    double dist_deg_lat = max_range / m_per_deg_lat;

    spdlog::debug("Radius of {:.3f} m is approx {:.6f} deg EW and {:.6f} deg NS", max_range, dist_deg_lon, dist_deg_lat);

    // Calculate our plot bounds based on these numbers
    min_lon = tx_site.lon - dist_deg_lon;
    max_lon = tx_site.lon + dist_deg_lon;
    double min_lat = tx_site.lat - dist_deg_lat;
    double max_lat = tx_site.lat + dist_deg_lat;

    // If doing P2P analysis, we need to make sure the RX site is within our whole degree bounds as well, so data is loaded
    // TODO: update this so it makes sense with the new approach

    if (ppa) {
        if (rx_site.lat < min_lat)
            min_lat = rx_site.lat;

        if (rx_site.lat > max_lat)
            max_lat = rx_site.lat;

        if (LonDiff(rx_site.lon, min_lon) < 0.0)
            min_lon = rx_site.lon;

        if (LonDiff(rx_site.lon, max_lon) >= 0.0)
            max_lon = rx_site.lon;

        spdlog::debug("RX site location expanded plot bounds to {:.6f}N {:.6f}E to {:.6f}N {:.6f}E", min_lat, min_lon, max_lat, max_lon);
    }

    bbox plot_bounds = getCircularBoundingBox({tx_site.lat, tx_site.lon}, max_range);

    spdlog::debug("Calculated plot boundaries: {:.6f}N {:.6f}E to {:.6f}N {:.6f}E", 
        plot_bounds.lower_left.lat, 
        plot_bounds.lower_left.lon, 
        plot_bounds.upper_right.lat, 
        plot_bounds.upper_right.lon
    );

    /* Load the required DEM tiles */
    if( (LoadTopoData(plot_bounds)) != 0 ){
        // This only fails on errors loading DEM tiles
        spdlog::error("Error loading topo data");
    }

    if (!ppa) {
        PlotPropagationRadius(tx_site);
        spdlog::debug("Finished PlotPropagationRadius()");

        if (cropping) {
            spdlog::debug("Cropping 1: N: {:.4f} S: {:.4f} E: {:.4f} W: {:.4f}",plot_bounds.upper_right.lat, plot_bounds.lower_left.lat, plot_bounds.upper_right.lon, plot_bounds.lower_left.lon);
            max_north = plot_bounds.upper_right.lat;
            min_north = plot_bounds.lower_left.lat;
            max_lon = plot_bounds.upper_right.lon;
            min_lon = plot_bounds.lower_left.lon;
        }

        // Write image
        if (prop_model == LOS)
            DoLOS(mapfile);
        else if (LR.erp == 0.0)
            DoPathLoss(mapfile);
        else if (dbm)
            DoRxdPwr(mapfile);
        else DoSigStr(mapfile);

        spdlog::info("Area boundaries:{:.6f} | {:.6f} | {:.6f} | {:.6f} ",max_north,max_lon,min_north,min_lon);

    } else {
        PlotPath(tx_site, rx_site);
        PathReport(tx_site, rx_site, output_filename.c_str(), 0, prop_model, rxGain);
        // Order flipped for benefit of graph. Makes no difference to data.
        SeriesData(rx_site, tx_site, output_filename.c_str(), 1, normalise);
    }

    auto end_time = std::chrono::steady_clock::now();
    double elapsed_s = std::chrono::duration<double>(end_time - start_time).count();
    fprintf(stderr, "Execution time: %.3f seconds\n", elapsed_s);

    if (debug) {
        fprintf(stderr, "--- Function call counters ---\n");
        fprintf(stderr, "  point_to_point_ITM : %d\n", cnt_point_to_point_ITM.load());
        fprintf(stderr, "  point_to_point     : %d\n", cnt_point_to_point.load());
        fprintf(stderr, "  computeLoss        : %d\n", cnt_computeLoss.load());
        fprintf(stderr, "  PlotPropPath       : %d\n", cnt_PlotPropPath.load());
        fprintf(stderr, "  ReadPath           : %d\n", cnt_ReadPath.load());
        fprintf(stderr, "  PlotPropagation    : %d\n", cnt_PlotPropagation.load());
    }

    return 0;
}
