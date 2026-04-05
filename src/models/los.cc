#include <stdio.h>
#include <math.h>
#include "../main.hh"
#include "los.hh"
#include "cost.hh"
#include "ecc33.hh"
#include "ericsson.hh"
#include "fspl.hh"
#include "hata.hh"
#include "itwom3.0.hh"
#include "sui.hh"
#include "pel.hh"
#include "egli.hh"
#include "soil.hh"
#include "../geo.hh"
#include <mutex>
#include <spdlog/spdlog.h>
#include <vector>
#include <limits.h>


#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

namespace {
	//bool ***processed; //unused
	bool has_init_processed = false;

    // Storage for processing threads
    std::vector<std::thread> threads;

    // Storage for processing thread futures
    std::vector<std::future<void *>> futures;

    // Mutex for processed vector
    std::mutex maskMutex;

    // Thread progress vector
    std::vector<progress_t> thread_progress;

    // 2D processed-flag array matching the flat DEM: processedPoints[x][y]
    std::vector<std::vector<bool>> processedPoints;

	void init_processed()
	{
        processedPoints.assign(dem_height_px, std::vector<bool>(dem_width_px, false));
        spdlog::debug("Initialized processedPoints vector of side {}x{}", dem_height_px, dem_width_px);
		has_init_processed = true;
	}

	bool can_process(double lat, double lon)
	{
		int x = (int)rint(ppd  * (lat - dem_min_lat));
		int y = (int)rint(ppd * (lon - dem_min_lon));
		bool rtn = false;

		if (x < 0 || x >= dem_height_px || y < 0 || y >= dem_width_px)
			return false;

		if (!processedPoints[x][y]) {
			maskMutex.lock();
			if (!processedPoints[x][y]) {
				rtn = true;
				processedPoints[x][y] = true;
			}
			maskMutex.unlock();
		}
		return rtn;
	}

    /**
     * Calulate a propagation for a specific range
     * 
     * @param *parameters parameters object for the propagation range
    */
	void* rangePropagation(progress_t &progress, void *parameters)
	{
        // Create propagationRange object based on our parameters
		PropagationRange *v = (PropagationRange*)parameters;

		alloc_elev();
		alloc_path();

        // Check if we're plotting a single line
        if (v->min_north == v->max_north && v->min_lon == v->max_lon) {
            spdlog::warn("Propagation plot range is a single point!");
        }

        // If our min & max lon coords are the same, it's a vertical line
        bool vertical = (v->min_lon == v->max_lon) ? true : false;

        // Calculate total number of points we are going to process
        unsigned int totalPoints = vertical ? (int)((v->max_north - v->min_north) / dpp) : (int)((v->max_lon - v->min_lon) / dpp);
        progress.total.store(totalPoints);

        // Init the count
        progress.count.store(0);

        spdlog::debug("Starting rangePropagation for {} range {:.6f}N {:.6f}E to {:.6f}N {:.6f}E, {} points at {:.8f} dpp [Segment {}]",
            vertical ? "vertical" : "horizontal",
            v->min_north, v->min_lon, v->max_north, v->max_lon, progress.total.load(), dpp, progress.id);

        // Init our varaibles for tracking position over the loop
        double lat = v->min_north;
        double lon = v->min_lon;
        int y = 0;
        // Iterate
		do {
			if (lon > 180.0)
				lon -= 360.0;

			site edge;
			edge.lat = lat;
			edge.lon = lon;
			edge.alt = v->altitude;

            //spdlog::debug("Plotting propagation path to {:.6f}N {:.6f}W", edge.lat, edge.lon);

			if(v->los)
				PlotLOSPath(v->source, edge);
			else
				PlotPropPath(v->source, edge, v->prop_model,
					v->knifeedge, v->pmenv);
            // Increment our counters
			++y;
            progress.count++;
            // Incremenet our lat/lon as needed
			if(vertical) {
                lat = (double)v->min_north + (dpp * (double)y);
            } else {
			    lon = (double)v->min_lon + (dpp * (double)y);
            }

        } while ( vertical ? (lat < (double)v->max_north) : (lon <= (double)v->max_lon) );

        free_elev();
        free_path();
		
		return NULL;
	}

    void* radiusPropagation(progress_t &progress, void *parameters)
    {
        // Create a prop radius from our parameters
        PropagationRadius *r = (PropagationRadius*)parameters;
        // Thread buffer allocation
        
		alloc_elev();
		alloc_path();

        // Check if our start & stop angles are the same
        if (r->start_angle_rad == r->stop_angle_rad)
        {
            spdlog::warn("Start & stop angles are the same, this radius segment will be a single line");
        }

        spdlog::debug("Starting radiusPropagation for range {:.2f} to {:.2f}, {} points, {:.8f} dpp",
            r->start_angle_rad * RAD2DEG, r->stop_angle_rad * RAD2DEG, r->points, dpp);

        // Get the amount in radians to increment per iteration
        double rps = (r->stop_angle_rad - r->start_angle_rad) / r->points;

        progress.total.store(r->points);

        progress.count.store(0);
        // Iterate
        double rad = r->start_angle_rad;
        site edge;
        edge.alt = r->altitude;
        for (int i = 0; i < r->points; i++)
        {
            // Get coordinates of point on circle
            coord point = getPointAtDistance({r->source.lat, r->source.lon}, r->radius, rad);
            // Create site for prop path
            edge.lat = point.lat;
            edge.lon = point.lon;

            // Plot
            if (r->los)
                PlotLOSPath(r->source, edge);
            else
                PlotPropPath(r->source, edge, r->prop_model, r->knifeedge, r->pmenv);

            // Increment
            rad += rps;
            progress.count++;
        }

        // Double check we covered the whole range
        if (rad < (r->stop_angle_rad - 0.01f))
        {
            spdlog::warn("Only got to {:.2f} degrees when we expected to get to {:.2f} degrees", rad * RAD2DEG, r->stop_angle_rad * RAD2DEG);
        }

        // Free the buffers we made earlier
        free_elev();
        free_path();

        return NULL;
    }

	/// @brief Wait for all threads in our threads array to finish
	void finishThreads()
	{
        for (auto& th : threads)
            th.join();
	}

    /// @brief Wait for the progress accumulators to finish, then finish out any running threads
    void finishProgress()
    {
        unsigned int total_points = 0;
        unsigned int points_processed = 0;
        
        // Calculate the expected total progress count
        for (const auto& p : thread_progress)
        {
            if (p.total <= 0) { std::this_thread::sleep_for( std::chrono::milliseconds(2) ); }
            total_points += p.total;
        }

        spdlog::debug("{} total points to process", total_points);
        
        // Await progress completion
        while (points_processed < total_points)
        {
            std::this_thread::sleep_for( std::chrono::milliseconds(500) );

            // Reset count for this check
            points_processed = 0;

            for (const auto& p : thread_progress)
            {
                points_processed += p.count;
            }

            // Update print
            spdlog::info("[{: 3d}%] Processing {}/{} points", int(points_processed * 100 / total_points), points_processed, total_points);
        }
    }
}

/*
 * Acute Angle from Rx point to an obstacle of height (opp) and
 * distance (adj)
 */
static double incidenceAngle(double opp, double adj)
{
	return atan2(opp, adj) * 180 / PI;
}

/*
 * Knife edge diffraction:
 * This is based upon a recognised formula like Huygens, but trades
 * thoroughness for increased speed which adds a proportional diffraction
 * effect to obstacles.
 */
static double ked(double freq, double rxh, double dkm)
{
	double obh, obd, rxobaoi = 0, d;

	obh = 0;		// Obstacle height
	obd = 0;		// Obstacle distance

	dkm = dkm * 1000;	// KM to metres

	// walk along path
	for (int n = 2; n < (dkm / elev[1]); n++) {

		d = (n - 2) * elev[1];	// no of points * delta = km

		//Find dip(s)
		if (elev[n] < obh) {

			// Angle from Rx point to obstacle
			rxobaoi =
			    incidenceAngle((obh - (elev[n] + rxh)), d - obd);
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

void PlotLOSPath(struct site source, struct site destination)
{
    /* This function analyzes the path between the source and
       destination locations. It determines which points along

       the path have line-of-sight visibility to the source.
       Points with line-of-sight visibility to the source are
           stored by setting bit 1 in the mask[][] array, which are
           displayed in green when PPM maps are later generated by ss. */

    bool bStop;
    int x, iCounter;
    double cos_angle, cos_test_angle, cos_horizon_angle, cos_limit_angle, rx_alt2;
    double distance, rx_alt, tx_alt, limit_alt, distance2, tx_alt2, test_alt, test_alt2, limit_alt2;

    ReadPath(source, destination);

    distance = 0.0;
    tx_alt = 0.0;
    distance2 = 0.0;
    tx_alt2 = 0.0;
    rx_alt2 = 0.0;
    test_alt = 0.0;
    test_alt2 = 0.0;
    limit_alt2 = 0.0;
    cos_horizon_angle = 1.0;
    bStop = false;
    iCounter = 0;

    /* altitude limit of 10000 meters */
    limit_alt = EARTHRADIUS + 10000.0;
    limit_alt2 = limit_alt * limit_alt;

    tx_alt = EARTHRADIUS + source.alt + path.elevation[0];
    tx_alt2 = tx_alt * tx_alt;

    for (x = 0; (bStop == false) && (x < (path.length - 1)) && (path.distance[x] <= max_range * 1000.0); x++) {

        if (x > 0) {
            distance = path.distance[x];
            distance2 = distance * distance;

            rx_alt = EARTHRADIUS + destination.alt + path.elevation[x];
            rx_alt2 = rx_alt * rx_alt;

            /* Calculate the cosine of the elevation between
               transmitter and receiver. */

            cos_angle = (distance2 + tx_alt2 - rx_alt2) / (2.0 * distance * tx_alt);

            cos_angle = MIN(1.0, MAX(-1.0, cos_angle));


            test_alt = EARTHRADIUS + (path.elevation[x] == 0.0 ? path.elevation[x] : path.elevation[x] + clutter);
            test_alt2 = test_alt * test_alt;

            /* Calculate the cosine of the elevation between
               transmitter and this test point. */

            cos_test_angle = (distance2 + tx_alt2 - test_alt2) / (2.0 * distance * tx_alt);
        }
        else {
            cos_angle = -1.0;
            cos_test_angle = 1.0;
        }

        if (cos_test_angle < cos_horizon_angle) {
            cos_horizon_angle = cos_test_angle;
        }

        if ((x > 0) && (cos_horizon_angle < 0.0)) {
            if (iCounter > 10) {
                /* Check for Mount Everest in line-of-sight visibility */

                /* Calculate the cosine of the elevation between
                   transmitter and altitude limit. */

                cos_limit_angle = (distance2 + tx_alt2 - limit_alt2) / (2.0 * distance * tx_alt);

                cos_limit_angle = MIN(1.0, MAX(-1.0, cos_limit_angle));

                if (cos_limit_angle > cos_horizon_angle) {
                    bStop = true;
                }

                iCounter = 0;
            }
            else {
                iCounter++;
            }
        }
    }
}

double computeLoss(PropModel model, double tx_alt, double rx_alt, double rx_terrain_alt,
                   double dkm, int pmenv, char *strmode, int &errnum)
{
    if (debug) cnt_computeLoss++;
    double loss = 0.0;

    switch (model) {
    case ITM_LR:
        point_to_point_ITM(tx_alt, rx_alt,
                           LR.eps_dielect, LR.sgm_conductivity, LR.eno_ns_surfref,
                           LR.frq_mhz, LR.radio_climate, LR.pol, LR.conf, LR.rel,
                           loss, strmode, errnum);
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
                       loss, strmode, errnum);
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
                           loss, strmode, errnum);
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
void PlotPropPath(
    struct site source,
    struct site destination,
    PropModel prop_model,
	int knifeedge,
    int pmenv
)
{
	if (debug) cnt_PlotPropPath++;
	int x, y, ifs, ofs, errnum;
	char block = 0, strmode[100];
	double loss, azimuth, pattern = 0.0,
	    xmtr_alt, dest_alt, xmtr_alt2, dest_alt2,
	    cos_rcvr_angle, cos_test_angle = 0.0, test_alt,
	    elevation = 0.0, distance = 0.0,
	    field_strength = 0.0, rxp, dBm, diffloss;
	struct site temp;
	float dkm;

	ReadPath(source, destination);


	for (x = 1; x < path.length - 1; x++)
		elev[x + 2] = (path.elevation[x] == 0.0 ? path.elevation[x] : (clutter + path.elevation[x]));

	/* Copy ending points without clutter */

	elev[2] = path.elevation[0];

	elev[path.length + 1] = path.elevation[path.length - 1];

	/* Since the only energy the Longley-Rice model considers
	   reaching the destination is based on what is scattered
	   or deflected from the first obstruction along the path,
	   we first need to find the location and elevation angle
	   of that first obstruction (if it exists).  This is done
	   using a 4/3rds Earth radius to match the model used by
	   Longley-Rice.  This information is required for properly
	   integrating the antenna's elevation pattern into the
	   calculation for overall path loss. */
	//if(debug)
	//	fprintf(stderr,"four_thirds_earth %.1f source.alt %.1f path.elevation[0] %.1f\n",four_thirds_earth,source.alt,path.elevation[0]);
	for (y = 2; (y < (path.length - 1) && path.distance[y] <= max_range * 1000.0);
	     y++) {
		/* Process this point only if it
		   has not already been processed. */

		if (can_process(path.lat[y], path.lon[y])) {

			distance = path.distance[y];
			xmtr_alt = FOUR_THIRDS_EARTH + source.alt + path.elevation[0];
			dest_alt = FOUR_THIRDS_EARTH + destination.alt + path.elevation[y];
			dest_alt2 = dest_alt * dest_alt;
			xmtr_alt2 = xmtr_alt * xmtr_alt;

			/* Calculate the cosine of the elevation of
			   the receiver as seen by the transmitter. */

			cos_rcvr_angle =
			    ((xmtr_alt2) + (distance * distance) -
			     (dest_alt2)) / (2.0 * xmtr_alt * distance);

			
            cos_rcvr_angle = MIN(1.0, MAX(-1.0, cos_rcvr_angle));

			if (got_elevation_pattern) {
				/* Determine the elevation angle to the first obstruction
				   along the path IF elevation pattern data is available
				   or an output (.ano) file has been designated. */

				for (x = 2, block = 0; (x < y && block == 0);
				     x++) {
					distance = path.distance[x];

					test_alt =
					    FOUR_THIRDS_EARTH +
					    (path.elevation[x] ==
					     0.0 ? path.elevation[x] : path.
					     elevation[x] + clutter);

					/* Calculate the cosine of the elevation
					   angle of the terrain (test point)
					   as seen by the transmitter. */

					cos_test_angle =
					    ((xmtr_alt2) +
					     (distance * distance) -
					     (test_alt * test_alt)) / (2.0 *
								       xmtr_alt
								       *
								       distance);


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
					elevation =
					    ((acos(cos_test_angle)) * RAD2DEG) -
					    90.0;
				else
					elevation =
					    ((acos(cos_rcvr_angle)) * RAD2DEG) -
					    90.0;
			}

			/* Determine attenuation for each point along the
			   path using a prop model starting at y=2 (number_of_points = 1), the
			   shortest distance terrain can play a role in
			   path loss. */

			elev[0] = y - 1;	/* (number of points - 1) */

			/* Distance between elevation samples (meters) */

			elev[1] = path.distance[y] - path.distance[y - 1];

			if (path.elevation[y] < 1) {
				path.elevation[y] = 1;
			}

			dkm = (elev[1] * elev[0]) / 1000;	// km

			loss = computeLoss(prop_model, source.alt, destination.alt,
			                   path.elevation[y] + destination.alt, dkm, pmenv,
			                   strmode, errnum);

			if (knifeedge == 1 && prop_model > 1) {
				diffloss =
				    ked(LR.frq_mhz, destination.alt, dkm);
				loss += (diffloss);	// ;)
			}
			//Key stage. Link dB for p2p is returned as 'loss'.

			temp.lat = path.lat[y];
			temp.lon = path.lon[y];

			azimuth = (Azimuth(source, temp));

			/* Integrate the antenna's radiation
			   pattern into the overall path loss. */

			x = (int)rint(10.0 * (10.0 - elevation));

			if (x >= 0 && x <= 1000) {
				azimuth = rint(azimuth);

				pattern =
				    (double)LR.antenna_pattern[(int)azimuth][x];

				if (pattern != 0.0) {
					pattern = 20.0 * log10(pattern);
					loss -= pattern;
				}
			}

			if (LR.erp != 0.0) {
				if (dbm) {
					/* dBm is based on EIRP (ERP + 2.14) */

					rxp =
					    LR.erp /
					    (pow(10.0, (loss - 2.14) / 10.0));

					dBm = 10.0 * (log10(rxp * 1000.0));


					/* Scale roughly between 0 and 255 */

					ifs = 200 + (int)rint(dBm);

					if (ifs < 0)
						ifs = 0;

					if (ifs > 255)
						ifs = 255;

					ofs =
					    GetSignal(path.lat[y], path.lon[y]);

					if (ofs > ifs)
						ifs = ofs;

					PutSignal(path.lat[y], path.lon[y],
						  (unsigned char)ifs);

				}

				else {
					field_strength =
					    (139.4 +
					     (20.0 * log10(LR.frq_mhz)) -
					     loss) +
					    (10.0 * log10(LR.erp / 1000.0));

					ifs = 100 + (int)rint(field_strength);

					if (ifs < 0)
						ifs = 0;

					if (ifs > 255)
						ifs = 255;

					ofs =
					    GetSignal(path.lat[y], path.lon[y]);

					if (ofs > ifs)
						ifs = ofs;

					PutSignal(path.lat[y], path.lon[y],
						  (unsigned char)ifs);

				}
			}

			else {
				if (loss > 255)
					ifs = 255;
				else
					ifs = (int)rint(loss);
				
				ofs = GetSignal(path.lat[y], path.lon[y]);

				if (ofs < ifs && ofs != 0)
					ifs = ofs;

				PutSignal(path.lat[y], path.lon[y],
					  (unsigned char)ifs);
			}
		}
	}

	if(path.lat[y]>cropLat)
		cropLat=path.lat[y];

	
	if(y>cropLon)
		cropLon=y;

	//if(cropLon>180)
	//	cropLon-=360;
}

void PlotLOSMap(struct site source, double altitude,
		uint8_t number_threads)
{
	/* This function performs a 360 degree sweep around the
	   transmitter site (source location), and plots the
	   line-of-sight coverage of the transmitter on the ss
	   generated topographic map based on a receiver located
	   at the specified altitude (in feet AGL).  Results
	   are stored in memory, and written out in the form
	   of a topographic map when the WritePPM() function
	   is later invoked. */

	// Four sections start here
	// Process north edge east/west, east edge north/south,
	// south edge east/west, west edge north/south
	double range_min_lon[] = {min_lon, min_lon, min_lon, max_lon};
	double range_min_north[] = {max_north, min_north, min_north, min_north};
	double range_max_lon[] = {max_lon, min_lon, max_lon, max_lon};
	double range_max_north[] = {max_north, max_north, min_north, max_north};
	PropagationRange *r = new PropagationRange[number_threads];

    // Size our progress vector appropriately
    thread_progress = std::vector<progress_t>(number_threads);

	for(int i = 0; i < number_threads; ++i) {
        r[i].los = true;

		r[i].eastwest = (range_min_lon[i] == range_max_lon[i] ? false : true);
		r[i].min_lon = range_min_lon[i];
		r[i].max_lon = range_max_lon[i];
		r[i].min_north = range_min_north[i];
		r[i].max_north = range_max_north[i];

		r[i].altitude = altitude;
		r[i].source = source;

        // Set the segment id
        thread_progress[i].id = i;

		futures.push_back( std::async( std::launch::async, rangePropagation, std::ref(thread_progress[i]), &r[i] ) );
	}

	finishThreads();

	delete[] r;
}
void PlotPropagationRadius(struct site source, double range, 
                            double altitude, PropModel prop_model, int knifeedge, int pmenv,
                            uint8_t number_threads)
{

    // Ensure number_threads is a logical value
    if ((number_threads % 2 != 0) && (number_threads % 3 != 0))
    {
        spdlog::error("Segment number must be an multiple of either 2 or 3!");
        exit(1);
    }

    // Get plot type string
    char plotType[32];
	if (LR.erp == 0.0 && debug)
		sprintf(plotType, "path loss");
	else {
		if (debug) {
			if (dbm)
				sprintf(plotType, "signal power level");
			else
				sprintf(plotType, "field strength");
		}
	}
    // Print debug
	spdlog::debug("Plotting {} contours of \"{}\" out to a radius of {:.2f} km with Rx antenna(s) at {:.2f} m AGL",
            plotType,
			source.name,
			range,
			altitude
    );

    // Optional clutter debug print
	if (clutter > 0.0)
        spdlog::debug("Using {:.2f} meters of ground clutter", clutter);

    // TX site location print
    spdlog::debug("TX site location: {:.6f}N {:.6f}W at {:.2f} m AGL", source.lat, source.lon, source.alt);

    // Get bounding box of plot
    bbox bounds = getCircularBoundingBox( {source.lat, source.lon}, range);

    // Calculate plot width & height in degrees
    double plot_width = bounds.upper_left.lon - bounds.lower_right.lon;
    double plot_height = bounds.upper_left.lat - bounds.lower_right.lat; 

    // Calculate the radius of our circle, in pixels
    //double radius_px = (plot_width / 2.0) * ppd;

    // Calculate the radius of our circle, in pixels
    //double radius_px = (plot_width + plot_height)/2.0/2.0 * ppd;

    // Calculate the total number of pixels/points in our plot circle using the midpoint circle algorithm
    // Borrowed from https://math.stackexchange.com/a/167310
    // We use the upper bound to ensure we don't miss any points
    //int circle_pixels = (int)ceil(radius_px * 2.0 * PI);
    // Ramanujan
    int circle_pixels = (int)ceil(ppd * M_PI / 2.0 * (3.0 * (plot_width + plot_height) - sqrt((3.0 * plot_width + plot_height) * (3.0 * plot_height + plot_width))));

    // Calculate the size of each angular degree section, in rads
    double section_size_rad = 360.0 / number_threads * DEG2RAD;

    // Calculate the number of points/pixels in each segment
    int section_pixels = (int)(circle_pixels / number_threads);

    // Create our ranges
    std::vector<PropagationRadius> radii;

    // Iterate through our number_threads
    for (int i = 0; i < number_threads; i++)
    {
        // Create a new radius
        PropagationRadius propRadius;
        // Populate static data
        propRadius.source = source;
        propRadius.radius = range;
        propRadius.points = section_pixels;
        propRadius.altitude = altitude;
        propRadius.prop_model = prop_model;
        propRadius.knifeedge = knifeedge;
        propRadius.pmenv = pmenv;
        // We're not doing LOS
        propRadius.los = false;
        // Calculate start and stop angles
        propRadius.start_angle_rad = i * section_size_rad;
        propRadius.stop_angle_rad = (i + 1) * section_size_rad;

        // Add to list
        radii.push_back(propRadius);
    }

	spdlog::debug("With {} threads and {} total points, our circular area will be divided into {:.2f}-degree (or {} point) segments", 
                    number_threads, circle_pixels, section_size_rad * RAD2DEG, section_pixels);

    // Make sure we didn't do anythng wrong
    if (radii.size() != number_threads) {
        spdlog::error("Our vector of radii ({}) does not match expected segment count {}", radii.size(), number_threads);
        exit(1);
    }

    // Size our progress vector appropriately
    thread_progress = std::vector<progress_t>(number_threads);

    // Init our vector for storing processing progress
    if (!has_init_processed) {
        init_processed();
    }

    // Iterate over the final list of ranges
    for (size_t i = 0; i < radii.size(); i++) {
        // Set the segment id
        thread_progress[i].id = i;
        // Start a thread
        spdlog::debug("Starting calc thread for radius segment {:.2f} to {:.2f}", radii[i].start_angle_rad * RAD2DEG, radii[i].stop_angle_rad * RAD2DEG);
        futures.push_back( std::async( std::launch::async, radiusPropagation, std::ref(thread_progress[i]), &radii[i] ) );

    }

    // Wait for futures to finish
    spdlog::debug("Waiting for threads to finish...");
    //finishProgress(); //indicate the progress but it is slow
    for (auto& f : futures)
        f.get();
    futures.clear();

    // Clean up our radii
	for(size_t i = 0; i < radii.size(); i++){
		radii.erase(radii.begin() + i);
	}
}

void PlotPath(struct site source, struct site destination)
{
	/* This function analyzes the path between the source and
	   destination locations.  It determines which points along
	   the path have line-of-sight visibility to the source.
	   Points along with path having line-of-sight visibility
	   to the source at an AGL altitude equal to that of the
	   destination location are stored by setting bit 1 in the
	   mask[][] array, which are displayed in green when PPM
	   maps are later generated by SPLAT!. */

	char block;
	int x, y;
	register double cos_xmtr_angle, cos_test_angle, test_alt;
	double distance, rx_alt, tx_alt;

	ReadPath(source, destination);

	for (y = 0; y < path.length; y++) {
		/* Test this point only if it hasn't been already
		   tested and found to be free of obstructions. */

			distance = path.distance[y];
			tx_alt = EARTHRADIUS + source.alt + path.elevation[0];
			rx_alt = EARTHRADIUS + destination.alt + path.elevation[y];

			/* Calculate the cosine of the elevation of the
			   transmitter as seen at the temp rx point. */

			cos_xmtr_angle =
			    ((rx_alt * rx_alt) + (distance * distance) -
			     (tx_alt * tx_alt)) / (2.0 * rx_alt * distance);

			for (x = y, block = 0; x >= 0 && block == 0; x--) {
				distance = path.distance[y] - path.distance[x];
				test_alt =
				    EARTHRADIUS + (path.elevation[x] ==
						   0.0 ? path.
						   elevation[x] : path.
						   elevation[x] + clutter);

				cos_test_angle =
				    ((rx_alt * rx_alt) + (distance * distance) -
				     (test_alt * test_alt)) / (2.0 * rx_alt *
							       distance);

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
