#ifndef _LOS_HH_
#define _LOS_HH_

#include <stdio.h>
#include <stdint.h>
#include <atomic>
#include <thread>
#include <future>
#include <ctime>

#include "../common.hh"

// Propagation models supported
enum PropModel {
    ITM_P2P = 0,
    ITM_LR = 1,
    LOS = 2,
    HATA = 3,
    ECC33 = 4,
    SUI = 5,
    COST231_HATA = 6,
    ITU_R = 7,
    ITWOM_3 = 8,
    ERICSSON = 9,
    PLANE_EARTH = 10,
    ELGI_V_U = 11,
    SOIL = 12,
};

// Rectangular bounding box propagation range
struct PropagationRange {
    double min_lon, max_lon, min_north, max_north;
    double altitude;
    bool eastwest, los;
    site source;
    PropModel prop_model;
    int knifeedge, pmenv;
};

// Angular propagation area
struct PropagationRadius {
    double start_angle_rad, stop_angle_rad;
    double radius;
    double altitude;
    bool los;
    site source;
    PropModel prop_model;
    int knifeedge, pmenv, points;
};

// Struct for storing thread progress
struct progress_t {
    int id;
    std::atomic<unsigned int> count {} ;
    std::atomic<unsigned int> total {} ;
};

void PlotLOSPath(struct site source, struct site destination);

void PlotPropPath(struct site source, struct site destination, PropModel propmodel, int knifeedge, int pmenv);

void PlotLOSMap(struct site source, double altitude, uint8_t number_threads);

/// @brief Plot propagation using a center point and circular radius. This plots around a circle instead of a rectangular bounding box and is theoretically more efficient.
/// @param source source transmitter
/// @param range maximum plot rage in miles or km
/// @param altitude altitude in ft or m
/// @param plot_filename output plot filename
/// @param prop_model propagation model to use
/// @param number_threads number_threads to split the plot circle into (must be a multiple of 2 or 3)
void PlotPropagationRadius(struct site source, double range, double altitude, PropModel prop_model, int knifeedge, int pmenv, 
                            uint8_t number_threads);

void PlotPath(struct site source, struct site destination);

double computeLoss(PropModel model, double tx_alt, double rx_alt, double rx_terrain_alt,
                   double dkm, int pmenv, char *strmode, int &errnum);

#endif /* _LOS_HH_ */
