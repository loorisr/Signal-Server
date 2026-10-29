#include <algorithm>
#include <cstdlib>
#include <stdio.h>
#include <math.h>
#include <atomic>
#include "../main.hh"
#include "los.hh"
#include "cost.hh"
#include "ecc33.hh"
#include "ericsson.hh"
#include "fspl.hh"
#include "hata.hh"
#include "itwom3.0.hh"
#include "include/itm.h"
#include "include/Enums.h"
#include "sui.hh"
#include "pel.hh"
#include "egli.hh"
#include "soil.hh"
#include "../geo.hh"
#include <spdlog/spdlog.h>
#include <vector>
#include <limits.h>

//#define MIN(a,b) (((a)<(b))?(a):(b))
//#define MAX(a,b) (((a)>(b))?(a):(b))

namespace {
	//bool ***processed; //unused

    // Storage for processing threads
    std::vector<std::thread> threads;

    // Storage for processing thread futures
    std::vector<std::future<void *>> futures;

    // Thread progress vector
    std::vector<progress_t> thread_progress;

    // Global ray angles, built once so a pixel can be assigned deterministically
    // to the same ray regardless of how work is split across threads.
    std::vector<float> ray_angles;

    // owner_map[x * dem_width_px + y] = global index of the first ray that
    // reaches this pixel in the sequential ray order.  This makes the
    // multi-threaded result identical to the single-threaded one.
    std::vector<int32_t> owner_map;

    // Ray currently being processed by this worker (thread-local).
    thread_local int current_ray_index = -1;

    // Dynamic work queue: threads pull the next global ray index.  Ownership
    // is already fixed by owner_map, so the order does not affect the output.
    std::atomic<int> next_ray_index{0};
    int total_ray_count = 0;

	// Check whether this pixel belongs to the ray currently being processed.
	// Also returns the DEM pixel coordinates so the caller does not have to
	// recompute them when storing the result.
	bool can_process(float lat, float lon, int &x_out, int &y_out)
	{
		int x = (int)rint(ppd  * (lat - dem_min_lat));
		int y = (int)rint(ppd * (lon - dem_min_lon));

		if (x < 0 || x >= dem_height_px || y < 0 || y >= dem_width_px)
			return false;

		x_out = x;
		y_out = y;

		/* The owner map gives the first ray that reaches this pixel.  A ray
		 * can cross the same pixel more than once (adjacent path samples
		 * round to the same cell); keep only the first crossing by
		 * invalidating the entry after it has been claimed. */
		int32_t &owner = owner_map[(size_t)x * dem_width_px + y];
		if (owner == current_ray_index) {
			owner = -2;
			return true;
		}
		return false;
	}

	/* Compute the pixels crossed by rays [first, last) into per-ray lists.
	 * Runs on worker threads; the lists are merged in global ray order by
	 * build_owner_map() so the ownership is deterministic. */
	void collect_ray_pixels(struct site source, int first, int last,
	                        std::vector<std::vector<int32_t>> &ray_pixels)
	{
		for (int r = first; r < last; ++r) {
			coord point = getPointAtDistance({source.lat, source.lon},
			                                 max_range, ray_angles[r]);
			site edge;
			edge.lat = point.lat;
			edge.lon = point.lon;
			edge.alt = altitudeLR;

			ReadPath(source, edge);

			std::vector<int32_t> &pixels = ray_pixels[r];
			pixels.clear();
			pixels.reserve((size_t)MAX(path.length - 2, 0));

			for (int y = 2; y < path.length - 1; ++y) {
				int x  = (int)rint(ppd * (path.lat[y] - dem_min_lat));
				int yy = (int)rint(ppd * (path.lon[y] - dem_min_lon));
				if (x < 0 || x >= dem_height_px || yy < 0 || yy >= dem_width_px)
					continue;
				pixels.push_back((int32_t)((size_t)x * dem_width_px + yy));
			}
		}
	}

	/* Assign every pixel to the first ray that reaches it in the global ray
	 * order.  The path sampling is parallel, but the merge is sequential and
	 * ordered, so the result is identical to the old single-threaded pass. */
	void build_owner_map(struct site source, int total_rays)
	{
		owner_map.assign((size_t)dem_height_px * dem_width_px, -1);

		std::vector<std::vector<int32_t>> ray_pixels((size_t)total_rays);

		const int workers = MAX(1, number_threads);
		const int base = total_rays / workers;
		const int remainder = total_rays % workers;

		std::vector<std::future<void>> futures;
		futures.reserve((size_t)workers);

		int first = 0;
		for (int w = 0; w < workers; ++w) {
			const int count = base + (w < remainder ? 1 : 0);
			if (count > 0)
				futures.push_back(std::async(std::launch::async,
					collect_ray_pixels, source, first, first + count,
					std::ref(ray_pixels)));
			first += count;
		}

		for (auto &f : futures)
			f.get();

		/* Merge in global ray order: the first ray wins. */
		for (int r = 0; r < total_rays; ++r) {
			for (int32_t idx : ray_pixels[(size_t)r]) {
				if (owner_map[(size_t)idx] < 0)
					owner_map[(size_t)idx] = r;
			}
		}
	}

    // Compute one angular slice of the propagation radius.
    void* radiusPropagation(progress_t &progress, void *parameters)
    {
        // Create a prop radius from our parameters
        PropagationRadius *r = (PropagationRadius*)parameters;
        // Thread buffer allocation

		/* elev and path are allocated dynamically by ReadPath on first call */

        spdlog::debug("Starting radius worker {}, {:d} ppd", progress.id, ppd);

        progress.total.store(total_ray_count);
        progress.count.store(0);

        site edge;
        edge.alt = altitudeLR;

        for (;;) {
            const int ray = next_ray_index.fetch_add(1, std::memory_order_relaxed);
            if (ray >= total_ray_count)
                break;

            const float rad = ray_angles[ray];

            // Get coordinates of point on circle
            coord point = getPointAtDistance({r->source.lat, r->source.lon}, max_range, rad);
            // Create site for prop path
            edge.lat = point.lat;
            edge.lon = point.lon;

            // Plot
            if (prop_model == LOS) {
                current_ray_index = -1;
                PlotLOSPath(r->source, edge);
            } else {
                current_ray_index = ray;
                PlotPropPath(r->source, edge);
            }

            progress.count++;
        }

        current_ray_index = -1;

        // Free the buffers we made earlier
        free_elev();
        free_path();

        return NULL;
    }
}

/*
 * Acute Angle from Rx point to an obstacle of height (opp) and
 * distance (adj)
 */
// Return the acute angle from a receiver point to an obstacle.
static float incidenceAngle(float opp, float adj)
{
	return atan2(opp, adj) * RAD2DEG;
}

/*
 * Knife edge diffraction:
 * This is based upon a recognised formula like Huygens, but trades
 * thoroughness for increased speed which adds a proportional diffraction
 * effect to obstacles.
 */
// Estimate knife-edge diffraction loss from the terrain profile.
static float ked(float freq, float rxh, float dm)
{
	float obh, obd, rxobaoi = 0, d;

	obh = 0;		// Obstacle height
	obd = 0;		// Obstacle distance

	// walk along path
	for (int n = 2; n < (dm / elev[1]); n++) {

		d = (n - 2) * elev[1];	// no of points * delta = km

		//Find dip(s)
		if (elev[n] < obh) {

			// Angle from Rx point to obstacle
			rxobaoi = incidenceAngle((obh - (elev[n] + rxh)), d - obd);
		} else {
			// Line of sight or higher
			rxobaoi = 0;
		}

		//note the highest point
		if (elev[n] > obh) {
			obh = elev[n];
			obd = d;
		}

	}

	if (rxobaoi >= 0) {
		return (rxobaoi / (300 / freq))+3;	// Diffraction angle divided by wavelength (m)
	} else {
		return 1;
	}
}

// to complete
// Mark line-of-sight points along the source-to-destination path.
void PlotLOSPath(struct site source, struct site destination)
{
    /* This function analyzes the path between the source and
       destination locations. It determines which points along
       the path have line-of-sight visibility to the source.
    */

    int x, y;
    char block = 0;
    float cos_test_angle, cos_rcvr_angle;
    float distance, tx_alt, distance_test, tx_alt2, test_alt, test_alt2, dest_alt, dest_alt2;

    ReadPath(source, destination);

    tx_alt = EARTHRADIUS + source.alt + GetElevation(source);
    tx_alt2 = tx_alt * tx_alt;

    distance = 0.0;

    for (y = 1; (y < (path.length - 1)) && distance <= max_range; y++) {
        distance = path.distance[y];
        dest_alt = EARTHRADIUS + destination.alt + path.elevation[y];
        dest_alt2 = dest_alt * dest_alt;
        
        cos_rcvr_angle = ((tx_alt2) + (distance * distance) - (dest_alt2)) / (2.0 * tx_alt * distance);
        cos_rcvr_angle = MIN(1.0, MAX(-1.0, cos_rcvr_angle));

        for (x = 2, block = 0; (x < y && block == 0); x++) {
            distance_test = path.distance[x];
            test_alt = EARTHRADIUS + path.elevation[x] + clutter;
            test_alt2 = test_alt * test_alt;

            /* Calculate the cosine of the elevation between transmitter and this test point. */
            cos_test_angle = ((tx_alt2) + (distance_test * distance_test) - (test_alt2)) / (2.0 * tx_alt * distance_test);
            cos_test_angle = std::clamp(cos_test_angle, -1.0f, 1.0f);

            if (cos_rcvr_angle >= cos_test_angle) { //block
                block = 1;
            } else {
                PutSignal(path.lat[x], path.lon[x], 200);
            }
        }
    }
}

// Select and run the configured propagation model for the current path.
float computeLoss(PropModel model, float tx_alt, float rx_alt, float rx_terrain_alt,
                   float dm, PropagationMode &mode, int &errnum)
{
    if (debug) cnt_computeLoss++;
    float loss = 0.0;
    float dkm = dm / 1000.0;

    switch (model) {
    case ITM_LR:
        point_to_point_ITM(tx_alt, rx_alt,
                           LR.eps_dielect, LR.sgm_conductivity, LR.eno_ns_surfref,
                           LR.frq_mhz, LR.radio_climate, LR.pol, LR.conf, LR.rel,
                           loss, mode, errnum);
        break;
    case ITM_NTIA:
        {
            long warnings = 0;

            errnum = ITM_P2P_CR(tx_alt, rx_alt, elev, LR.radio_climate, LR.eno_ns_surfref,
                                LR.frq_mhz, LR.pol, LR.eps_dielect, LR.sgm_conductivity,
                                MDVAR__SINGLE_MESSAGE_MODE, LR.conf * 100.0, LR.rel * 100.0,
                                &loss, &warnings);
            mode = PROP_MODE_NONE;
        }
        break;
    case HATA:
        loss = HATApathLoss(LR.frq_mhz, tx_alt, rx_terrain_alt, dkm, pmenv);
        break;
    case ECC33:
        loss = ECC33pathLoss(LR.frq_mhz, tx_alt, rx_terrain_alt, dkm, pmenv);
        break;
    case SUI:
        loss = SUIpathLoss(LR.frq_mhz, tx_alt, rx_terrain_alt, dkm, pmenv);
        break;
    case COST231_HATA:
        loss = COST231pathLoss(LR.frq_mhz, tx_alt, rx_terrain_alt, dkm, pmenv);
        break;
    case ITU_R:
        loss = FSPLpathLoss(LR.frq_mhz, dkm, false);
        break;
    case ITWOM_3:
        point_to_point(tx_alt, rx_alt,
                       LR.eps_dielect, LR.sgm_conductivity, LR.eno_ns_surfref,
                       LR.frq_mhz, LR.radio_climate, LR.pol, LR.conf, LR.rel,
                       loss, mode, errnum);
        break;
    case ERICSSON:
        loss = EricssonpathLoss(LR.frq_mhz, tx_alt, rx_terrain_alt, dkm, pmenv);
        break;
    case PLANE_EARTH:
        loss = PlaneEarthLoss(dkm, tx_alt, rx_terrain_alt);
        break;
    case ELGI_V_U:
        loss = EgliPathLoss(LR.frq_mhz, tx_alt, rx_terrain_alt, dkm);
        break;
    case SOIL:
        loss = SoilPathLoss(LR.frq_mhz, dkm, LR.eps_dielect);
        break;
    default:
        spdlog::warn("Defaulting to ITM propagation model");
        point_to_point_ITM(tx_alt, rx_alt,
                           LR.eps_dielect, LR.sgm_conductivity, LR.eno_ns_surfref,
                           LR.frq_mhz, LR.radio_climate, LR.pol, LR.conf, LR.rel,
                           loss, mode, errnum);
        break;
    }

    return loss;
}

/**
 * Calculate propagation for the points on a line between two coordinates
 *
 * @param source - the source site object
 * @param destination - the destination site object
*/
// Plot the full propagation result for every point along a path.
void PlotPropPath(
    struct site source,
    struct site destination
)
{
	if (debug) cnt_PlotPropPath++;
	int x, y, errnum, azimuth, px, py;
	char block = 0;
	PropagationMode mode = PROP_MODE_NONE;
	float loss, pattern = 0.0,
	    xmtr_alt, dest_alt, xmtr_alt2, dest_alt2,
	    cos_rcvr_angle, cos_test_angle = 0.0, test_alt,
	    elevation = 0.0, distance = 0.0, distance_test,
	    field_strength = 0.0, dBm, diffloss;

	ReadPath(source, destination);

	azimuth = rint(Azimuth(source, destination)); // for antenna pattern

	elev[1] = path.distance[2] - path.distance[1]; // distance between samples
    // add clutter
	for (x = 1; x < path.length - 1; x++)
		elev[x + 2] = clutter + path.elevation[x];

	/* Copy ending points without clutter */

	elev[2] = path.elevation[0];
	elev[path.length + 1] = path.elevation[path.length - 1];

	/* The ITM profile for receiver point y is the prefix elev[2..y+1].
	 * Build the prefix sums once per radial so qlrpfl() no longer rescans
	 * the profile for zsys and for each z1sq1() window. */
	thread_local ITM_ray_ctx ray_ctx;
	ray_ctx.build(elev, path.length + 1);

	/* Constant parts of the dBm / field-strength conversions. */
	const float dbm_const = (LR.erp > 0.0f)
		? 30.0f + 10.0f * log10f(LR.erp) + 2.14f : 0.0f;
	const float field_const = (LR.erp > 0.0f)
		? 139.4f + 20.0f * log10f(LR.frq_mhz) + 10.0f * log10f(LR.erp / 1000.0f)
		: 0.0f;


	/* Since the only energy the Longley-Rice model considers
	   reaching the destination is based on what is scattered
	   or deflected from the first obstruction along the path,
	   we first need to find the location and elevation angle
	   of that first obstruction (if it exists).  This is done
	   using a 4/3rds Earth radius to match the model used by
	   Longley-Rice.  This information is required for properly
	   integrating the antenna's elevation pattern into the
	   calculation for overall path loss. */

	xmtr_alt = FOUR_THIRDS_EARTH + source.alt + path.elevation[0];
	xmtr_alt2 = xmtr_alt * xmtr_alt;

	/* Init the precomputed ITM context once for this path; avoids repeating
	 * constant-parameter work (qerfi, complex sqrt, etc.) on every point. */
	const bool fast_itm = (prop_model == ITM_LR);
	ITM_ctx itm_ctx;
	if (fast_itm)
		itm_ctx = ITM_ctx_init(source.alt, LR.eps_dielect, LR.sgm_conductivity,
		                       LR.eno_ns_surfref, LR.frq_mhz, LR.radio_climate,
		                       LR.pol, LR.conf, LR.rel);
	float segment_start_dh = 0.0f;
	float segment_end_dh = 0.0f;
	int segment_start_y = -1;
	int segment_end_y = -1;
	bool interp_segment_active = false;

	/* Quadratic interpolation needs the anchor before the current segment. */
	float previous_start_dh = 0.0f;
	int previous_start_y = -1;
	bool previous_start_valid = false;

	for (y = 2; (y < (path.length - 1) && distance <= max_range);  y++) {
		/* Process this point only if it
		   has not already been processed. */
		if (can_process(path.lat[y], path.lon[y], px, py)) {

			distance = path.distance[y];
			dest_alt = FOUR_THIRDS_EARTH + destination.alt + path.elevation[y];
			dest_alt2 = dest_alt * dest_alt;

			/* Calculate the cosine of the elevation of
			   the receiver as seen by the transmitter. */

			cos_rcvr_angle = ((xmtr_alt2) + (distance * distance) - (dest_alt2)) / (2.0 * xmtr_alt * distance);
            cos_rcvr_angle = MIN(1.0, MAX(-1.0, cos_rcvr_angle));

			if (got_elevation_pattern) {
				/* Determine the elevation angle to the first obstruction
				   along the path IF elevation pattern data is available */
				for (x = 2, block = 0; (x < y && block == 0); x++) {
					distance_test = path.distance[x];

					test_alt = FOUR_THIRDS_EARTH + path.elevation[x] + clutter;

					/* Calculate the cosine of the elevation
					   angle of the terrain (test point)
					   as seen by the transmitter. */

					cos_test_angle = ((xmtr_alt2) + (distance_test * distance_test) - (test_alt * test_alt)) / (2.0 * xmtr_alt * distance_test);
                    cos_test_angle = MIN(1.0, MAX(-1.0, cos_test_angle));    

					/* Compare these two angles to determine if
					   an obstruction exists.  Since we're comparing
					   the cosines of these angles rather than
					   the angles themselves, the sense of the
					   following "if" statement is reversed from
					   what it would be if the angles themselves
					   were compared. */

					if (cos_rcvr_angle >= cos_test_angle)
						block = 1;
				}

				if (block)
					elevation = (acos(std::clamp(cos_test_angle, -1.0f, 1.0f)) * RAD2DEG) - 90.0;
				else
					elevation = (acos(std::clamp(cos_rcvr_angle, -1.0f, 1.0f)) * RAD2DEG) - 90.0;
			}

			/* Determine attenuation for each point along the
			   path using a prop model starting at y=2 (number_of_points = 1), the
			   shortest distance terrain can play a role in
			   path loss. */

			elev[0] = y - 1;	/* (number of points - 1) */

			if (path.elevation[y] < 1) {
				path.elevation[y] = 1;
			}

			if (fast_itm) {
				if (debug) cnt_computeLoss++;
				if (fast_dh_stride > 0 && y > fast_dh_stride) {
					const bool start_new_segment = (!interp_segment_active || y >= segment_end_y);

					if (start_new_segment) {
						const float old_start_dh = segment_start_dh;
						const int old_start_y = segment_start_y;
						const bool old_valid = (old_start_y >= 0);

						/* Compute actual loss and delta-h at this segment start. */
						point_to_point_ITM_fast(itm_ctx, destination.alt, elev, loss, mode, errnum,
						                        &segment_start_dh, false, &ray_ctx);

						const int next_y = MIN(y + fast_dh_stride, path.length - 2);
						segment_start_y = y;
						segment_end_y = next_y;

						if (fast_quadratic && old_valid) {
							previous_start_dh = old_start_dh;
							previous_start_y = old_start_y;
							previous_start_valid = true;
						} else {
							previous_start_valid = false;
						}

						if (next_y > y) {
							/* Look ahead once to get delta-h at the next full-compute point. */
							const float saved_elev0 = elev[0];
							float tmp_loss = 0.0f;
							PropagationMode tmp_mode = PROP_MODE_NONE;
							int tmp_errnum = 0;

							elev[0] = next_y - 1;
							point_to_point_ITM_fast(itm_ctx, destination.alt, elev, tmp_loss, tmp_mode, tmp_errnum,
							                        &segment_end_dh, false, &ray_ctx);
							elev[0] = saved_elev0;
							interp_segment_active = true;
						} else {
							interp_segment_active = false;
						}
					} else {
						float interp_dh;
						if (fast_quadratic && previous_start_valid) {
							const float x  = path.distance[y];
							const float x0 = path.distance[previous_start_y];
							const float x1 = path.distance[segment_start_y];
							const float x2 = path.distance[segment_end_y];

							const float l0 = ((x - x1) * (x - x2)) /
							                 ((x0 - x1) * (x0 - x2));
							const float l1 = ((x - x0) * (x - x2)) /
							                 ((x1 - x0) * (x1 - x2));
							const float l2 = ((x - x0) * (x - x1)) /
							                 ((x2 - x0) * (x2 - x1));
							interp_dh = previous_start_dh * l0
							          + segment_start_dh * l1
							          + segment_end_dh * l2;

							const float lo = MIN(previous_start_dh,
							                     MIN(segment_start_dh, segment_end_dh));
							const float hi = MAX(previous_start_dh,
							                     MAX(segment_start_dh, segment_end_dh));
							interp_dh = std::clamp(interp_dh, lo, hi);
						} else {
							const float denom = float(segment_end_y - segment_start_y);
							const float t = (denom > 0.0f)
								? (float(y - segment_start_y) / denom) : 0.0f;
							interp_dh = segment_start_dh
							          + t * (segment_end_dh - segment_start_dh);
						}

						float dh_for_point = interp_dh;
						point_to_point_ITM_fast(itm_ctx, destination.alt, elev, loss, mode, errnum,
						                        &dh_for_point, true, &ray_ctx);
					}
				} else {
					interp_segment_active = false;
					point_to_point_ITM_fast(itm_ctx, destination.alt, elev, loss, mode, errnum,
					                        NULL, false, &ray_ctx);
				}
			} else {
				loss = computeLoss(prop_model, source.alt, destination.alt,
				                   path.elevation[y] + destination.alt, distance,
				                   mode, errnum);
			}

			if (knifeedge && prop_model > 1) {
				diffloss = ked(LR.frq_mhz, destination.alt, distance);
				loss += diffloss;	
			}
			//Key stage. Link dB for p2p is returned as 'loss'.


			/* Integrate the antenna's radiation
			   pattern into the overall path loss. */

			x = (int)rint(10.0 * (10.0 - elevation));

			if (x >= 0 && x <= 1000) {
				pattern = LR.antenna_pattern[azimuth][x]; /* already in dB */
				loss -= pattern;
			}

			if (LR.erp != 0.0) {
				if (dbm) {
					/* dBm is based on EIRP (ERP + 2.14).  The pow/log10
					 * pair of the original formula is constant minus loss. */
					dBm = dbm_const - loss;
					PutSignalXY(px, py, dBm);
				}
				else {
					field_strength = field_const - loss;
					PutSignalXY(px, py, field_strength);
				}
			}

			else {
				PutSignalXY(px, py, loss);
			}
		}
	}
}

// Divide a circular plot area into worker slices and render each one.
void PlotPropagationRadius(struct site source)
{
    // Get plot type string
	if (debug) {
        char plotType[32];
        if (LR.erp == 0.0)
            snprintf(plotType, sizeof(plotType), "path loss");
        else if (dbm)
            snprintf(plotType, sizeof(plotType), "signal power level");
        else
            snprintf(plotType, sizeof(plotType), "field strength");
        // Print debug
        spdlog::debug("Plotting {} contours out to a radius of {:.2f} km with Rx antenna(s) at {:.2f} m AGL",
                plotType,
                max_range,
                altitudeLR
        );
	}

    // Optional clutter debug print
	if (clutter > 0.0)
        spdlog::debug("Using {:.2f} meters of ground clutter", clutter);

    // TX site location print
    spdlog::debug("TX site location: {:.6f}N {:.6f}W at {:.2f} m AGL", source.lat, source.lon, source.alt);

    // Get bounding box of plot
    bbox bounds = getCircularBoundingBox( {source.lat, source.lon}, max_range);

    // Calculate plot width & height in degrees
    float plot_width = bounds.upper_right.lon - bounds.lower_left.lon;
    float plot_height = bounds.upper_right.lat - bounds.lower_left.lat; 

    // Ramanujan
    int circle_pixels = (int)ceil(ppd * M_PI / 2.0 * (3.0 * (plot_width + plot_height) - sqrt((3.0 * plot_width + plot_height) * (3.0 * plot_height + plot_width))));

    /* Build the global ray sequence once.  Using the same sequence for any
     * thread count makes the ownership map (and therefore the output)
     * independent of the number of workers. */
    ray_angles.resize(circle_pixels);
    /* Match the historical single-thread angle sequence exactly: the old code
     * used 360 * DEG2RAD (a float) divided by the ray count. */
    const float full_circle_rad = (float)(360.0 * DEG2RAD);
    const float ray_step = full_circle_rad / (float)circle_pixels;
    /* volatile prevents -ffast-math from rewriting the recurrence as i*step,
     * which would round differently from the historical loop. */
    volatile float rad = 0.0f;
    for (int i = 0; i < circle_pixels; i++) {
        ray_angles[i] = rad;
        rad = rad + ray_step;
    }

    // Split the global rays into number_threads contiguous slices
    std::vector<PropagationRadius> radii((size_t)number_threads);
    for (int i = 0; i < number_threads; i++)
        radii[(size_t)i].source = source;

    spdlog::debug("With {} threads and {} total points, rays are distributed dynamically ({:.2f} degree step)",
                  number_threads, circle_pixels, ray_step * RAD2DEG);

    // Size our progress vector appropriately
    thread_progress = std::vector<progress_t>(number_threads);

    /* Pre-compute which ray owns each pixel.  Doing this once reproduces the
     * old single-threaded first-claim order exactly, while allowing the
     * expensive ITM pass to run in parallel. */
    if (prop_model != LOS)
        build_owner_map(source, circle_pixels);

    next_ray_index.store(0, std::memory_order_relaxed);
    total_ray_count = circle_pixels;

    // Start the workers; each one pulls the next ray from the shared queue.
    for (size_t i = 0; i < radii.size(); i++) {
        thread_progress[i].id = (int)i;
        spdlog::debug("Starting calc worker {}", i);
        futures.push_back( std::async( std::launch::async, radiusPropagation, std::ref(thread_progress[i]), &radii[i] ) );
    }

    // Wait for futures to finish
    spdlog::debug("Waiting for threads to finish...");
    for (auto& f : futures)
        f.get();
    futures.clear();
}

// Check whether any terrain point blocks the direct path.
void PlotPath(struct site source, struct site destination)
{
	/* This function analyzes the path between the source and
	   destination locations. */

	char block;
	int x, y;
	float cos_xmtr_angle, cos_test_angle, test_alt;
	float distance, rx_alt, tx_alt;

	ReadPath(source, destination);

	for (y = 1; y < path.length; y++) {
		distance = path.distance[y];
		tx_alt = EARTHRADIUS + source.alt + path.elevation[0];
		rx_alt = EARTHRADIUS + destination.alt + path.elevation[y];

		/* Calculate the cosine of the elevation of the
		   transmitter as seen at the temp rx point. */

		cos_xmtr_angle =
		    ((rx_alt * rx_alt) + (distance * distance) -
		     (tx_alt * tx_alt)) / (2.0 * rx_alt * distance);

		for (x = y - 1, block = 0; x >= 0 && block == 0; x--) {
			distance = path.distance[y] - path.distance[x];
			if (distance == 0.0)
				continue;
			test_alt =
			    EARTHRADIUS + (path.elevation[x] == 0.0
			                   ? path.elevation[x]
			                   : path.elevation[x] + clutter);

			cos_test_angle =
			    ((rx_alt * rx_alt) + (distance * distance) -
			     (test_alt * test_alt)) / (2.0 * rx_alt * distance);

			/* Compare these two angles to determine if
			   an obstruction exists.  Since we're comparing
			   the cosines of these angles rather than
			   the angles themselves, the following "if"
			   statement is reversed from what it would
			   be if the actual angles were compared. */

			if (cos_xmtr_angle >= cos_test_angle)
				block = 1;
		}
	}
}
