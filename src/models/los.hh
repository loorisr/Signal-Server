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
    ITM_NTIA = 13,
};

// Angular propagation area
struct PropagationRadius {
    double start_angle_rad, stop_angle_rad;
    site source;
    int points;
};

extern PropModel prop_model;

// Struct for storing thread progress
struct progress_t {
    int id;
    std::atomic<unsigned int> count {} ;
    std::atomic<unsigned int> total {} ;
};

// Plot line-of-sight signal along the path between two sites.
void PlotLOSPath(struct site source, struct site destination);

// Plot propagation loss or signal along the path between two sites.
void PlotPropPath(struct site source, struct site destination);

// Plot a radial propagation surface centered on the source site.
void PlotPropagationRadius(struct site source);

// Check whether the direct path between two sites is blocked.
void PlotPath(struct site source, struct site destination);

// Compute propagation loss for the selected model and path geometry.
double computeLoss(PropModel model, double tx_alt, double rx_alt, double rx_terrain_alt, double dm, PropagationMode &mode, int &errnum);

#endif /* _LOS_HH_ */
