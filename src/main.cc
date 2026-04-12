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

std::string DEM_path;
std::string color_palette = "heat";
std::string output_filename;

double max_range = 0.0,
    fzone_clearance = 0.6, clutter, tercon, terdic,
    dBm, loss, field_strength,
    min_lat = 90, max_lat = -90, min_lon = 180.0, max_lon = -180.0,
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
bool knifeedge = false;
int pmenv = 1;
int number_threads = 4;
double altitudeLR = 1;
PropModel prop_model = ITM_LR;
bool ppa = false;
int normalise = 0;
std::string mapfile;

thread_local double *elev;
thread_local int elev_allocated = 0;
thread_local struct path path;
thread_local int path_allocated = 0;
site tx_site;
site rx_site;
double         **dem_data   = nullptr;
int **dem_signal = nullptr;
int dem_min_lat   = 0;
int dem_min_lon   = 0;
int dem_width_px  = 0;
int dem_height_px = 0;

struct LR LR;

// Initialize the run, load terrain data, and dispatch the selected outputs.
int main(int argc, char *argv[])
{
    auto start_time = std::chrono::steady_clock::now();

    parse_cmdline(argc, argv);
    bbox geo_bounds;

    if (ppa) { // TODO: handle +-90 and +/- 180° wrap 
        min_lat = MIN(rx_site.lat, tx_site.lat);
        max_lat = MAX(rx_site.lat, tx_site.lat);
        min_lon = MIN(rx_site.lon, tx_site.lon);
        max_lon = MAX(rx_site.lon, tx_site.lon);
        geo_bounds.lower_left = {min_lat, min_lon};
        geo_bounds.upper_right = {max_lat, max_lon};
    } else {
        geo_bounds = getCircularBoundingBox({tx_site.lat, tx_site.lon}, max_range);
    }
    
    spdlog::info("Calculated boundaries: {:.6f}N {:.6f}E to {:.6f}N {:.6f}E", 
        geo_bounds.lower_left.lat, 
        geo_bounds.lower_left.lon, 
        geo_bounds.upper_right.lat, 
        geo_bounds.upper_right.lon
    );
    

    /* Load the required DEM tiles */
    if (LoadTopoData(geo_bounds) != 0) {
        spdlog::error("Error loading topo data");
        return 1;
    }

    //for plots
    min_lon = geo_bounds.lower_left.lon;
    max_lon = geo_bounds.upper_right.lon;
    min_lat = geo_bounds.lower_left.lat;
    max_lat = geo_bounds.upper_right.lat;

    if (!ppa) {
        PlotPropagationRadius(tx_site);
        spdlog::debug("Finished PlotPropagationRadius()");

        // Write image
        if (prop_model == LOS)
            DoLOS(mapfile.c_str());
        else if (LR.erp == 0.0)
            DoPathLoss(mapfile.c_str());
        else if (dbm)
            DoRxdPwr(mapfile.c_str());
        else DoSigStr(mapfile.c_str());

        spdlog::debug("Area boundaries:{:.6f} | {:.6f} | {:.6f} | {:.6f} ",geo_bounds.upper_right.lat, geo_bounds.lower_left.lat, geo_bounds.upper_right.lon, geo_bounds.lower_left.lon);

    } else {
        PlotPath(tx_site, rx_site);
        PathReport(tx_site, rx_site, output_filename.c_str(), 0, prop_model, rxGain);
        // Order flipped for benefit of graph. Makes no difference to data.
        SeriesData(rx_site, tx_site, output_filename.c_str(), 1, normalise);
    }

    auto end_time = std::chrono::steady_clock::now();
    double elapsed_s = std::chrono::duration<double>(end_time - start_time).count();
    spdlog::info("Execution time: {:.3f} seconds", elapsed_s);

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
