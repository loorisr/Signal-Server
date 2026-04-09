#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <algorithm>
#include <climits>
#include <optional>
#include <vector>

#include "tinycolormap.hpp"

#include <spdlog/spdlog.h>

#include "common.hh"
#include "main.hh"
#include "inputs.hh"
#include "models/los.hh"

static tinycolormap::ColormapType get_colormap()
{
	if (color_palette == "parula")    return tinycolormap::ColormapType::Parula;
	if (color_palette == "magma")     return tinycolormap::ColormapType::Magma;
	if (color_palette == "inferno")   return tinycolormap::ColormapType::Inferno;
	if (color_palette == "plasma")    return tinycolormap::ColormapType::Plasma;
	if (color_palette == "viridis")   return tinycolormap::ColormapType::Viridis;
	if (color_palette == "cividis")   return tinycolormap::ColormapType::Cividis;
	if (color_palette == "heat")      return tinycolormap::ColormapType::Heat;
	if (color_palette == "hot")       return tinycolormap::ColormapType::Hot;
	if (color_palette == "turbo")     return tinycolormap::ColormapType::Turbo;
	if (color_palette == "hsv")       return tinycolormap::ColormapType::HSV;
	if (color_palette == "cubehelix") return tinycolormap::ColormapType::Cubehelix;
	if (color_palette == "github")    return tinycolormap::ColormapType::Github;
	if (color_palette == "gray")      return tinycolormap::ColormapType::Gray;
	spdlog::warn("Unknown color palette '{}', falling back to heat", color_palette);
	return tinycolormap::ColormapType::Heat;
}

template<typename ClassifyFn>
static void render_geotiff(const char *filename,
                            tinycolormap::ColormapType colormap,
                            double scale_min, double scale_max, bool reverse,
                            ClassifyFn classify)
{
	const double conversion = 255.0 / pow((double)(max_elevation - min_elevation), ONE_OVER_GAMMA);
	const double scale_range = scale_max - scale_min;

	spdlog::debug("Writing \"{}\" ({} x {} pixmap image)...", filename, width, height);

	if (!geotiff || filename == NULL) return;

	std::vector<uint8_t> rgba((size_t)width * height * 4, 0);
	uint8_t *p = rgba.data();

	for (int y = 0; y < (int)height; y++) {
		double lat = max_north - dpp * y;
		for (int x = 0; x < (int)width; x++) {
			double lon = min_lon + dpp * x;
			uint8_t r = 0, g = 0, b = 0, a = 0;
			int x0, y0;
			if (find_dem_xy(lat, lon, x0, y0)) {
				auto sig = classify(dem_signal[x0][y0]);
				if (!sig) {
					if (!ngs) {
						/* Grayscale terrain */
						unsigned terrain = (unsigned)(0.5 + pow((double)(dem_data[x0][y0] - min_elevation), ONE_OVER_GAMMA) * conversion);
						r = g = b = (uint8_t)terrain;
						a = 255;
					}
				} else {
					double t = std::clamp((*sig - scale_min) / scale_range, 0.0, 1.0);
					if (reverse) t = 1.0 - t;
					auto c = tinycolormap::GetColor(t, colormap);
					r = c.ri(); g = c.gi(); b = c.bi(); a = 255;
				}
			}
			*p++ = r; *p++ = g; *p++ = b; *p++ = a;
		}
	}
	write_geotiff_rgba(rgba.data(), width, height, filename);
}

void DoPathLoss(char *filename)
{
	render_geotiff(filename, get_colormap(), 80.0, 230.0, true,
		[](int sig) -> std::optional<double> {
			if (sig == 0 || (contour_threshold != 0 && sig > abs(contour_threshold)))
				return std::nullopt;
			return sig;
		});
}

int DoSigStr(char *filename)
{
	render_geotiff(filename, get_colormap(), 108.0, 228.0, false,
		[](int sig) -> std::optional<double> {
			if (contour_threshold != 0 && sig < contour_threshold)
				return std::nullopt;
			return sig;
		});
	return 0;
}

void DoRxdPwr(char *filename)
{
	render_geotiff(filename, get_colormap(),
		(double)contour_threshold, 0.0, false,
		[](int sig) -> std::optional<double> {
			if (contour_threshold != 0 && sig < contour_threshold)
				return std::nullopt;
			return sig;
		});
}

void DoLOS(char *filename)
{
	render_geotiff(filename, get_colormap(), 0.0, 1.0, /*reverse=*/false,
		[](int) -> std::optional<double> { return std::nullopt; });
}

void PathReport(struct site source, struct site destination, const char *name,
		char graph_it, PropModel propmodel, double rxGain)
{
	/* This function writes a PPA Path Report (name.txt) to
	   the filesystem.  If (graph_it == 1), then gnuplot is invoked
	   to generate an appropriate output file indicating the Longley-Rice
	   model loss between the source and destination locations.
	   "filename" is the name assigned to the output file generated
	   by gnuplot.  The filename extension is used to set gnuplot's
	   terminal setting and output file type.  If no extension is
	   found, .png is assumed. */

	int x, y, z, errnum;
	char basename[255], term[30], ext[15], report_name[80], block = 0;
	PropagationMode mode = PROP_MODE_NONE;
	double maxloss = -100000.0, minloss = 100000.0, angle1, angle2,
	    azimuth, pattern = 1.0, patterndB = 0.0,
	    total_loss = 0.0, cos_xmtr_angle, cos_test_angle = 0.0,
	    source_alt, test_alt, dest_alt, source_alt2, dest_alt2,
	    distance, elevation,
	    free_space_loss = 0.0, eirp =
	    0.0, voltage, rxp, power_density, dkm;
	FILE *fd = NULL, *fd2 = NULL;

	snprintf(report_name, 80, "%s.txt", name);

	fd2 = fopen(report_name, "w");
	if (fd2 == NULL) {
		spdlog::error("PathReport: cannot open {} for writing: {}", report_name, strerror(errno));
		return;
	}

	fprintf(fd2, "\n\t\t--==[ Path Profile Analysis ]==--\n\n");
	fprintf(fd2, "Transmitter site: Tx\n");

	fprintf(fd2, "Site location: %.4f, %.4f\n", source.lat, source.lon);

	fprintf(fd2, "Ground elevation: %.2f meters AMSL\n", GetElevation(source));
	fprintf(fd2,
		"Antenna height: %.2f meters AGL / %.2f meters AMSL\n",
		source.alt, source.alt + GetElevation(source));

	const double link_distance = Distance(source, destination);

	azimuth = Azimuth(source, destination);
	angle1 = ElevationAngle(source, destination);
	angle2 = ElevationAngle2(source, destination, EARTHRADIUS);

	if (got_azimuth_pattern || got_elevation_pattern) {
		x = (int)rint(10.0 * (10.0 - angle2));

		if (x >= 0 && x <= 1000)
			pattern =
			    (double)LR.antenna_pattern[(int)rint(azimuth)][x];

		patterndB = 20.0 * log10(pattern);
	}

	fprintf(fd2, "Distance to Rx: %.2f kilometers\n",
		link_distance / 1000.0);

	fprintf(fd2, "Azimuth to Rx: %.2f degrees grid\n",
		azimuth);


	fprintf(fd2, "Downtilt angle to Rx: %+.4f degrees\n",
		angle1);



	/* Receiver */

	fprintf(fd2, "\nReceiver site: Rx\n");

	fprintf(fd2, "Site location: %.4f, %.4f\n", destination.lat, destination.lon);

	fprintf(fd2, "Ground elevation: %.2f meters AMSL\n", GetElevation(destination));
	fprintf(fd2,
		"Antenna height: %.2f meters AGL / %.2f meters AMSL\n",
		destination.alt, destination.alt + GetElevation(destination));

	fprintf(fd2, "Distance to Rx: %.2f kilometers\n",
		link_distance / 1000.0);

	azimuth = Azimuth(destination, source);

	angle1 = ElevationAngle(destination, source);
	angle2 = ElevationAngle2(destination, source, EARTHRADIUS);

	fprintf(fd2, "Azimuth to Tx: %.2f degrees grid\n", azimuth);


	fprintf(fd2, "Downtilt angle to Rx: %+.4f degrees\n",
		angle1);

	if (LR.frq_mhz > 0.0) {
		fprintf(fd2, "\n\nPropagation model: ");

		switch (propmodel) {
		case ITM_P2P:
		case ITM_LR:
			fprintf(fd2, "Irregular Terrain Model\n");
			break;
		case ITM_NTIA:
			fprintf(fd2, "NTIA Irregular Terrain Model\n");
			break;
		case LOS:
			fprintf(fd2, "Line of sight\n");
			break;
		case HATA:
			fprintf(fd2, "Okumura-Hata\n");
			break;
		case ECC33:
			fprintf(fd2, "ECC33 (ITU-R P.529)\n");
			break;
		case SUI:
			fprintf(fd2, "Stanford University Interim\n");
			break;
		case COST231_HATA:
			fprintf(fd2, "COST231-Hata\n");
			break;
		case ITU_R:
			fprintf(fd2, "Free space path loss (ITU-R.525)\n");
			break;
		case ITWOM_3:
			fprintf(fd2, "ITWOM 3.0\n");
			break;
		case ERICSSON:
			fprintf(fd2, "Ericsson\n");
			break;
		case PLANE_EARTH:
			fprintf(fd2, "Plane Earth\n");
			break;
		case ELGI_V_U:
			fprintf(fd2, "Egli VHF/UHF\n");
			break;
		case SOIL:
			fprintf(fd2, "Soil\n");
			break;
		}

		fprintf(fd2, "Model sub-type: ");

		switch (pmenv) {
		case 1:
			fprintf(fd2, "City / Conservative\n");
			break;
		case 2:
			fprintf(fd2, "Suburban / Average\n");
			break;
		case 3:
			fprintf(fd2, "Rural / Optimistic\n");
			break;
		}
		fprintf(fd2, "Earth's Dielectric Constant: %.3lf\n",
			LR.eps_dielect);
		fprintf(fd2, "Earth's Conductivity: %.3lf Siemens/meter\n",
			LR.sgm_conductivity);
		fprintf(fd2,
			"Atmospheric Bending Constant (N-units): %.3lf ppm\n",
			LR.eno_ns_surfref);
		fprintf(fd2, "Frequency: %.3lf MHz\n", LR.frq_mhz);
		fprintf(fd2, "Radio Climate: %d (", LR.radio_climate);

		switch (LR.radio_climate) {
		case 1:
			fprintf(fd2, "Equatorial");
			break;

		case 2:
			fprintf(fd2, "Continental Subtropical");
			break;

		case 3:
			fprintf(fd2, "Maritime Subtropical");
			break;

		case 4:
			fprintf(fd2, "Desert");
			break;

		case 5:
			fprintf(fd2, "Continental Temperate");
			break;

		case 6:
			fprintf(fd2, "Maritime Temperate, Over Land");
			break;

		case 7:
			fprintf(fd2, "Maritime Temperate, Over Sea");
			break;

		default:
			fprintf(fd2, "Unknown");
		}

		fprintf(fd2, ")\nPolarisation: %d (", LR.pol);

		if (LR.pol == 0)
			fprintf(fd2, "Horizontal");

		if (LR.pol == 1)
			fprintf(fd2, "Vertical");

		fprintf(fd2, ")\nFraction of Situations: %.1lf%c\n",
			LR.conf * 100.0, 37);
		fprintf(fd2, "Fraction of Time: %.1lf%c\n", LR.rel * 100.0, 37);

		if (LR.erp != 0.0) {
			fprintf(fd2, "\nReceiver gain: %.1f dBd / %.1f dBi\n", rxGain, rxGain+2.14);
			fprintf(fd2, "Transmitter ERP plus Receiver gain: ");

			if (LR.erp < 1.0)
				fprintf(fd2, "%.1lf milliwatts",
					1000.0 * LR.erp);

			if (LR.erp >= 1.0 && LR.erp < 10.0)
				fprintf(fd2, "%.1lf Watts", LR.erp);

			if (LR.erp >= 10.0 && LR.erp < 10.0e3)
				fprintf(fd2, "%.0lf Watts", LR.erp);

			if (LR.erp >= 10.0e3)
				fprintf(fd2, "%.3lf kilowatts", LR.erp / 1.0e3);

			dBm = 10.0 * (log10(LR.erp * 1000.0));
			fprintf(fd2, " (%+.2f dBm)\n", dBm);
			fprintf(fd2, "Transmitter ERP minus Receiver gain: %.2f dBm\n", dBm-rxGain);

			/* EIRP = ERP + 2.14 dB */

			fprintf(fd2, "Transmitter EIRP plus Receiver gain: ");

			eirp = LR.erp * 1.636816521;

			if (eirp < 1.0)
				fprintf(fd2, "%.1lf milliwatts", 1000.0 * eirp);

			if (eirp >= 1.0 && eirp < 10.0)
				fprintf(fd2, "%.1lf Watts", eirp);

			if (eirp >= 10.0 && eirp < 10.0e3)
				fprintf(fd2, "%.0lf Watts", eirp);

			if (eirp >= 10.0e3)
				fprintf(fd2, "%.3lf kilowatts", eirp / 1.0e3);

			dBm = 10.0 * (log10(eirp * 1000.0));
			fprintf(fd2, " (%+.2f dBm)\n", dBm);

			// Rx gain
			fprintf(fd2, "Transmitter EIRP minus Receiver gain: %.2f dBm\n", dBm-rxGain);
		}

		fprintf(fd2, "\nSummary for the link between Rx and Tx:\n\n");

		if (patterndB != 0.0)
		        fprintf(fd2, "Rx antenna pattern towards Tx: %.3f (%.2f dB)\n",
				pattern,
				patterndB);

		ReadPath(source, destination);	/* source=TX, destination=RX */

		/* Copy elevations plus clutter along
		   path into the elev[] array. */

		for (x = 1; x < path.length - 1; x++)
			elev[x + 2] =
			    (path.elevation[x] ==
			     0.0 ? path.elevation[x] : (clutter + path.elevation[x]));

		/* Copy ending points without clutter */

		elev[2] = path.elevation[0];
		elev[path.length + 1] = path.elevation[path.length - 1];

		azimuth = rint(Azimuth(source, destination));

		for (y = 2; y < (path.length - 1); y++) {	/* path.length-1 avoids LR error */
			distance = path.distance[y];

			source_alt = FOUR_THIRDS_EARTH + source.alt + path.elevation[0];
			dest_alt = FOUR_THIRDS_EARTH + destination.alt +
			    path.elevation[y];
			dest_alt2 = dest_alt * dest_alt;
			source_alt2 = source_alt * source_alt;

			/* Calculate the cosine of the elevation of
			   the receiver as seen by the transmitter. */

			cos_xmtr_angle =
			    ((source_alt2) + (distance * distance) -
			     (dest_alt2)) / (2.0 * source_alt * distance);

			if (got_elevation_pattern) {
				/* If an antenna elevation pattern is available, the
				   following code determines the elevation angle to
				   the first obstruction along the path. */

				for (x = 2, block = 0; x < y && block == 0; x++) {
					distance =
					    path.distance[y] - path.distance[x];
					test_alt =
					    FOUR_THIRDS_EARTH +
					    path.elevation[x];

					/* Calculate the cosine of the elevation
					   angle of the terrain (test point)
					   as seen by the transmitter. */

					cos_test_angle =
					    ((source_alt2) +
					     (distance * distance) -
					     (test_alt * test_alt)) / (2.0 *
								       source_alt
								       *
								       distance);

					/* Compare these two angles to determine if
					   an obstruction exists.  Since we're comparing
					   the cosines of these angles rather than
					   the angles themselves, the sense of the
					   following "if" statement is reversed from
					   what it would be if the angles themselves
					   were compared. */

					if (cos_xmtr_angle >= cos_test_angle)
						block = 1;
				}

				/* At this point, we have the elevation angle
				   to the first obstruction (if it exists). */
			}

			/* Determine path loss for each point along the
			   path using Longley-Rice's point_to_point mode
			   starting at x=2 (number_of_points = 1), the
			   shortest distance terrain can play a role in
			   path loss. */

			elev[0] = y - 1;	/* (number of points - 1) */

			/* Distance between elevation samples (meters) */

			elev[1] = path.distance[y] - path.distance[y - 1];

			dkm = (elev[1] * elev[0]) / 1000;	// km

			loss = computeLoss(propmodel, source.alt, destination.alt,
			                   path.elevation[y] + destination.alt, dkm,
			                   mode, errnum);

			if (block)
				elevation =
				    (acos(std::clamp(cos_test_angle, -1.0, 1.0)) * RAD2DEG) - 90.0;
			else
				elevation =
				    (acos(std::clamp(cos_xmtr_angle, -1.0, 1.0)) * RAD2DEG) - 90.0;

			/* Integrate the antenna's radiation
			   pattern into the overall path loss. */

			x = (int)rint(10.0 * (10.0 - elevation));

			if (x >= 0 && x <= 1000) {
				pattern =
				    (double)LR.antenna_pattern[(int)azimuth][x];

				if (pattern != 0.0){
					patterndB = 20.0 * log10(pattern);
				}else{
					patterndB = 0.0;
				}
			}

			else
				patterndB = 0.0;

			total_loss = loss - patterndB;

			if (total_loss > maxloss)
				maxloss = total_loss;

			if (total_loss < minloss)
				minloss = total_loss;
		}

		distance = link_distance;

		if (distance != 0.0) {
			free_space_loss =
			    36.6 + (20.0 * log10(LR.frq_mhz)) +
			    (20.0 * log10(distance));
			fprintf(fd2, "Free space path loss: %.2f dB\n",
				free_space_loss);
		}

		fprintf(fd2, "Computed path loss: %.2f dB\n", loss);


        if ((loss * 1.5) < free_space_loss) {
			fprintf(fd2, "Model error! Computed loss of %.1f dB is less than free space loss of %.1f dB. Check your inputs for model %d\n", loss, free_space_loss, propmodel);
			spdlog::error("Model error! Computed loss of {:.1f} dB is less than free space loss of {:.1f} dB. Check your inputs for model {}", loss, free_space_loss, static_cast<int>(propmodel));
			fclose(fd2);
			return;
        }

		if (free_space_loss != 0.0)
			fprintf(fd2,
				"Attenuation due to terrain shielding: %.2f dB\n",
				loss - free_space_loss);

		if (patterndB != 0.0)
		        fprintf(fd2,"Total path loss including Tx antenna pattern: %.2f dB\n",
				total_loss);

		if (LR.erp != 0.0) {
			field_strength =
			    (139.4 + (20.0 * log10(LR.frq_mhz)) - total_loss) +
			    (10.0 * log10(LR.erp / 1000.0));

			/* dBm is referenced to EIRP */

			rxp = eirp / (pow(10.0, (total_loss / 10.0)));
			dBm = 10.0 * (log10(rxp * 1000.0));
			power_density =
			    (eirp /
			     (pow
			      (10.0, (total_loss - free_space_loss) / 10.0)));
			/* divide by 4*PI*distance_in_meters squared */
			power_density /= (4.0 * PI * distance * distance *
					  2589988.11);

			fprintf(fd2, "Field strength at Rx: %.2f dBuV/meter\n",
				field_strength);
			fprintf(fd2, "Signal power level at Rx: %+.2f dBm\n",
				dBm);
			fprintf(fd2,
				"Signal power density at Rx: %+.2f dBW per square meter\n",
				10.0 * log10(power_density));
			voltage =
			    1.0e6 * sqrt(50.0 *
					 (eirp /
					  (pow
					   (10.0,
					    (total_loss - 2.14) / 10.0))));
			fprintf(fd2,
				"Voltage across 50 ohm dipole at Rx: %.2f uV (%.2f dBuV)\n",
				voltage,
				20.0 * log10(voltage));

			voltage =
			    1.0e6 * sqrt(75.0 *
					 (eirp /
					  (pow
					   (10.0,
					    (total_loss - 2.14) / 10.0))));
			fprintf(fd2,
				"Voltage across 75 ohm dipole at Rx: %.2f uV (%.2f dBuV)\n",
				voltage,
				20.0 * log10(voltage));
		}

		if (propmodel == ITM_LR || propmodel == ITM_P2P || propmodel == ITM_NTIA) {
			fprintf(fd2, "Longley-Rice model error number: %d",
				errnum);

			switch (errnum) {
			case 0:
				fprintf(fd2, " (No error)\n");
				break;

			case 1:
				fprintf(fd2,
					"\n  Warning: Some parameters are nearly out of range.\n");
				fprintf(fd2,
					"  Results should be used with caution.\n");
				break;

			case 2:
				fprintf(fd2,
					"\n  Note: Default parameters have been substituted for impossible ones.\n");
				break;

			case 3:
				fprintf(fd2,
					"\n  Warning: A combination of parameters is out of range for this model.\n");
				fprintf(fd2,
					"  Results should be used with caution.\n");
				break;

			default:
				fprintf(fd2,
					"\n  Warning: Some parameters are out of range for this model.\n");
				fprintf(fd2,
					"  Results should be used with caution.\n");
			}
		}

	}

	ObstructionAnalysis(source, destination, LR.frq_mhz, fd2);
	fclose(fd2);

	spdlog::debug("Path loss: {:.1f} dB, Received Power: {:.1f} dBm, Field strength {:.1f} dBuV", loss, dBm, field_strength);

	/* Skip plotting the graph if ONLY a path-loss report is needed. */

	if (graph_it) {
		if (name[0] == '.') {
			/* Default filename and output file type */

			strncpy(basename, "profile", 8);
			strncpy(term, "png", 4);
			strncpy(ext, "png", 4);
		}

		else {
			/* Extract extension and terminal type from "name" */

			ext[0] = 0;
			y = strlen(name);
			strncpy(basename, name, 254);

			for (x = y - 1; x > 0 && name[x] != '.'; x--) ;

			if (x > 0) {	/* Extension found */
				for (z = x + 1; z <= y && (z - (x + 1)) < 10;
				     z++) {
					ext[z - (x + 1)] = tolower(name[z]);
					term[z - (x + 1)] = name[z];
				}

				ext[z - (x + 1)] = 0;	/* Ensure an ending 0 */
				term[z - (x + 1)] = 0;
				basename[x] = 0;
			}
		}

		if (ext[0] == 0) {	/* No extension -- Default is png */
			strncpy(term, "png", 4);
			strncpy(ext, "png", 4);
		}

		/* Either .ps or .postscript may be used
		   as an extension for postscript output. */

		if (strncmp(term, "postscript", 10) == 0)
			strncpy(ext, "ps", 3);

		else if (strncmp(ext, "ps", 2) == 0)
			strncpy(term, "postscript enhanced color", 26);

		fd = fopen("ppa.gp", "w");
		if (fd == NULL) {
			spdlog::error("PathReport: cannot open ppa.gp for writing: {}", strerror(errno));
			fclose(fd2);
			return;
		}

		fprintf(fd, "set grid\n");
		fprintf(fd, "set yrange [%2.3f to %2.3f]\n", minloss, maxloss);
		fprintf(fd, "set encoding iso_8859_1\n");
		fprintf(fd, "set term %s\n", term);
		fprintf(fd,
			"set title \"Path Loss Profile Along Path Between Rx and Tx (%.2f%c azimuth)\"\n",
			Azimuth(destination,
							       source), 176);

		fprintf(fd,
			"set xlabel \"Distance Between Rx and Tx (%.2f kilometers)\"\n",
			Distance(destination, source));

		if (got_azimuth_pattern || got_elevation_pattern)
			fprintf(fd,
				"set ylabel \"Total Path Loss (including TX antenna pattern) (dB)");
		else
			fprintf(fd, "set ylabel \"Longley-Rice Path Loss (dB)");

		fprintf(fd, "\"\nset output \"%s.%s\"\n", basename, ext);
		fprintf(fd,
			"plot \"profile.gp\" title \"Path Loss\" with lines\n");

		fclose(fd);

		x = system("gnuplot ppa.gp");
	}

}

void SeriesData(struct site source, struct site destination, const char *name,
		unsigned char fresnel_plot, unsigned char normalised)
{
	int x;
	char profilename[255], referencename[255], cluttername[255],
	    curvaturename[255], fresnelname[255], fresnel60name[255];
	double a, b, c, height = 0.0, refangle, cangle, maxheight =
	    -100000.0, minheight = 100000.0, lambda = 0.0, f_zone =
	    0.0, fpt6_zone = 0.0, nm = 0.0, nb = 0.0, ed = 0.0, es = 0.0, r =
	    0.0, d = 0.0, d1 = 0.0, terrain, azimuth, distance;
	struct site remote;
	FILE *fd = NULL, *fd1 = NULL, *fd2 = NULL, *fd3 = NULL, *fd4 =
	    NULL, *fd5 = NULL;

	ReadPath(destination, source);
	azimuth = Azimuth(destination, source);
	distance = Distance(destination, source);
	refangle = ElevationAngle(destination, source);
	b = GetElevation(destination) + destination.alt + EARTHRADIUS;

	spdlog::debug("SeriesData: az = {}, dist = {}, ref = {}, b = {}", azimuth, distance, refangle, b);
	
	if (fresnel_plot) {
		lambda = 299792458.0 / (LR.frq_mhz * 1e6);
		d = path.distance[path.length - 1];
	}

	if (normalised) {
		ed = GetElevation(destination);
		es = GetElevation(source);
		nb = -destination.alt - ed;
		nm = (-source.alt - es - nb) / (path.distance[path.length - 1]);
	}

	snprintf(profilename,   sizeof(profilename),   "%s_profile",   name);
	snprintf(referencename, sizeof(referencename), "%s_reference", name);
	snprintf(cluttername,   sizeof(cluttername),   "%s_clutter",   name);
	snprintf(curvaturename, sizeof(curvaturename), "%s_curvature", name);
	snprintf(fresnelname,   sizeof(fresnelname),   "%s_fresnel",   name);
	snprintf(fresnel60name, sizeof(fresnel60name), "%s_fresnel60", name);

	fd = fopen(profilename, "wb");
	if (fd == NULL) { spdlog::error("SeriesData: cannot open {}", profilename); return; }
	if (clutter > 0.0) {
		fd1 = fopen(cluttername, "wb");
		if (fd1 == NULL) { spdlog::error("SeriesData: cannot open {}", cluttername); fclose(fd); return; }
	}
	fd2 = fopen(referencename, "wb");
	if (fd2 == NULL) { spdlog::error("SeriesData: cannot open {}", referencename); fclose(fd); if (fd1) fclose(fd1); return; }
	fd5 = fopen(curvaturename, "wb");
	if (fd5 == NULL) { spdlog::error("SeriesData: cannot open {}", curvaturename); fclose(fd); if (fd1) fclose(fd1); fclose(fd2); return; }

	if ((LR.frq_mhz >= 20.0) && (LR.frq_mhz <= 100000.0) && fresnel_plot) {
		fd3 = fopen(fresnelname, "wb");
		if (fd3 == NULL) { spdlog::error("SeriesData: cannot open {}", fresnelname); fclose(fd); if (fd1) fclose(fd1); fclose(fd2); fclose(fd5); return; }
		fd4 = fopen(fresnel60name, "wb");
		if (fd4 == NULL) { spdlog::error("SeriesData: cannot open {}", fresnel60name); fclose(fd); if (fd1) fclose(fd1); fclose(fd2); fclose(fd5); fclose(fd3); return; }
	}

	for (x = 0; x < path.length - 1; x++) {
		remote.lat = path.lat[x];
		remote.lon = path.lon[x];
		remote.alt = 0.0;
		terrain = GetElevation(remote);
		if (x == 0)
			terrain += destination.alt;	/* RX antenna spike */

		a = terrain + EARTHRADIUS;
		cangle = Distance(destination, remote) / EARTHRADIUS;
		c = b * sin(refangle * DEG2RAD + HALFPI) / sin(HALFPI -
							       refangle *
							       DEG2RAD -
							       cangle);
		height = a - c;

		/* Per Fink and Christiansen, Electronics
		 * Engineers' Handbook, 1989:
		 *
		 *   H = sqrt(lamba * d1 * (d - d1)/d)
		 *
		 * where H is the distance from the LOS
		 * path to the first Fresnel zone boundary.
		 */

		if ((LR.frq_mhz >= 20.0) && (LR.frq_mhz <= 100000.0)
		    && fresnel_plot) {
			d1 = path.distance[x];
			f_zone = -1.0 * sqrt(lambda * d1 * (d - d1) / d);
			fpt6_zone = f_zone * fzone_clearance;
		}

		if (normalised) {
			r = -(nm * path.distance[x]) - nb;
			height += r;

			if ((LR.frq_mhz >= 20.0) && (LR.frq_mhz <= 100000.0)
			    && fresnel_plot) {
				f_zone += r;
				fpt6_zone += r;
			}
		}

		else
			r = 0.0;

		if (height > 0) {
			fprintf(fd, "%.3f %.3f\n", path.distance[x], height);
		}

		if (fd1 != NULL && x > 0 && x < path.length - 2)
			fprintf(fd1, "%.3f %.3f\n", path.distance[x],
				(terrain == 0.0 ? height : (height + clutter)));

		fprintf(fd2, "%.3f %.3f\n", path.distance[x], r);
		fprintf(fd5, "%.3f %.3f\n", path.distance[x], height - terrain);

		if ((LR.frq_mhz >= 20.0) && (LR.frq_mhz <= 100000.0)
		    && fresnel_plot) {
			fprintf(fd3, "%.3f %.3f\n", path.distance[x], f_zone);
			fprintf(fd4, "%.3f %.3f\n", path.distance[x], fpt6_zone);

			if (f_zone < minheight)
				minheight = f_zone;
		}

		if ((height + clutter) > maxheight)
			maxheight = height + clutter;

		if (height < minheight)
			minheight = height;

		if (r > maxheight)
			maxheight = r;

	}			// End of loop

	if (normalised)
		r = -(nm * path.distance[path.length - 1]) - nb;
	else
		r = 0.0;

	fprintf(fd, "%.3f %.3f", path.distance[path.length - 1], r);
	fprintf(fd2, "%.3f %.3f", path.distance[path.length - 1], r);

	if ((LR.frq_mhz >= 20.0) && (LR.frq_mhz <= 100000.0) && fresnel_plot) {
		fprintf(fd3, "%.3f %.3f", path.distance[path.length - 1], r);
		fprintf(fd4, "%.3f %.3f", path.distance[path.length - 1], r);
	}

	if (r > maxheight)
		maxheight = r;

	if (r < minheight)
		minheight = r;

	fclose(fd);

	if (fd1 != NULL)
		fclose(fd1);

	fclose(fd2);
	fclose(fd5);

	if ((LR.frq_mhz >= 20.0) && (LR.frq_mhz <= 100000.0) && fresnel_plot) {
		fclose(fd3);
		fclose(fd4);
	}

}
