/* GDAL must come before any header that defines MAX, because some GDAL
 * versions guard both MIN and MAX under a single #if !defined(MAX) block. */
#include <gdal.h>

#include <errno.h>
#include <limits.h>
#include <math.h>
#include <spdlog/spdlog.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "common.hh"
#include "main.hh"

extern char *color_file;

extern double antenna_rotation, antenna_downtilt, antenna_dt_direction;

int LoadPAT(char *az_filename, char *el_filename)
{
	/* This function reads and processes antenna pattern (.az
		 and .el) files that may correspond in name to previously
		 loaded ss .lrp files or may be user-supplied by cmdline.  */

	int a, b, w, x, y, z, last_index, next_index, span;
	char string[255], *pointer = NULL;
	float az, xx, elevation, amplitude, rotation, valid1, valid2, delta, azimuth[361], azimuth_pattern[361], el_pattern[10001],
			elevation_pattern[361][1001], slant_angle[361], tilt, mechanical_tilt = 0.0, tilt_azimuth, tilt_increment, sum;
	FILE *fd = NULL;
	unsigned char read_count[10001];

	rotation = 0.0;

	got_azimuth_pattern = false;
	got_elevation_pattern = false;

	/* Load .az antenna pattern file */

	if (az_filename != NULL && (fd = fopen(az_filename, "r")) == NULL && errno != ENOENT)
		/* Any error other than file not existing is an error */
		return errno;

	if (fd != NULL) {
		spdlog::debug("Antenna Pattern Azimuth File = [{}]", az_filename);

		/* Clear azimuth pattern array */
		for (x = 0; x <= 360; x++) {
			azimuth[x] = 0.0;
			read_count[x] = 0;
		}

		/* Read azimuth pattern rotation
			 in degrees measured clockwise
			 from true North. */

		if (fgets(string, 254, fd) == NULL) {
			// fprintf(stderr,"Azimuth read error\n");
			// exit(0);
		}
		pointer = strchr(string, ';');

		if (pointer != NULL) *pointer = 0;

		if (antenna_rotation != -1)  // If cmdline override
			rotation = (float)antenna_rotation;
		else
			sscanf(string, "%f", &rotation);

		spdlog::debug("Antenna Pattern Rotation = {}", rotation);

		/* Read azimuth (degrees) and corresponding
			 normalized field radiation pattern amplitude
			 (0.0 to 1.0) until EOF is reached. */

		if (fgets(string, 254, fd) == NULL) {
			// fprintf(stderr,"Azimuth read error\n");
			// exit(0);
		}
		pointer = strchr(string, ';');

		if (pointer != NULL) *pointer = 0;

		sscanf(string, "%f %f", &az, &amplitude);

		do {
			x = (int)rintf(az);

			if (x >= 0 && x <= 360 && fd != NULL) {
				azimuth[x] += amplitude;
				read_count[x]++;
			}

			if (fgets(string, 254, fd) == NULL)
				break;
			pointer = strchr(string, ';');

			if (pointer != NULL) *pointer = 0;

			sscanf(string, "%f %f", &az, &amplitude);

		} while (1);

		fclose(fd);
		fd = NULL;

		/* Handle 0=360 degree ambiguity */

		if ((read_count[0] == 0) && (read_count[360] != 0)) {
			read_count[0] = read_count[360];
			azimuth[0] = azimuth[360];
		}

		if ((read_count[0] != 0) && (read_count[360] == 0)) {
			read_count[360] = read_count[0];
			azimuth[360] = azimuth[0];
		}

		/* Average pattern values in case more than
			 one was read for each degree of azimuth. */

		for (x = 0; x <= 360; x++) {
			if (read_count[x] > 1) azimuth[x] /= (float)read_count[x];
		}

		/* Interpolate missing azimuths
			 to completely fill the array */

		last_index = -1;
		next_index = -1;

		for (x = 0; x <= 360; x++) {
			if (read_count[x] != 0) {
				if (last_index == -1)
					last_index = x;
				else
					next_index = x;
			}

			if (last_index != -1 && next_index != -1) {
				valid1 = azimuth[last_index];
				valid2 = azimuth[next_index];

				span = next_index - last_index;
				delta = (valid2 - valid1) / (float)span;

				for (y = last_index + 1; y < next_index; y++) azimuth[y] = azimuth[y - 1] + delta;

				last_index = y;
				next_index = -1;
			}
		}

		/* Perform azimuth pattern rotation
			 and load azimuth_pattern[361] with
			 azimuth pattern data in its final form. */

		for (x = 0; x < 360; x++) {
			y = x + (int)rintf(rotation);

			if (y >= 360) y -= 360;

			azimuth_pattern[y] = azimuth[x];
		}

		azimuth_pattern[360] = azimuth_pattern[0];

		got_azimuth_pattern = true;
	}

	/* Read and process .el file */

	if (el_filename != NULL && (fd = fopen(el_filename, "r")) == NULL && errno != ENOENT)
		/* Any error other than file not existing is an error */
		return errno;

	if (fd != NULL) {
		spdlog::debug("Antenna Pattern Elevation File = [{}]", el_filename);

		/* Clear azimuth pattern array */

		for (x = 0; x <= 10000; x++) {
			el_pattern[x] = 0.0;
			read_count[x] = 0;
		}

		/* Read mechanical tilt (degrees) and
			 tilt azimuth in degrees measured
			 clockwise from true North. */

		if (fgets(string, 254, fd) == NULL) {
			// fprintf(stderr,"Tilt read error\n");
			// exit(0);
		}
		pointer = strchr(string, ';');

		if (pointer != NULL) *pointer = 0;

		sscanf(string, "%f %f", &mechanical_tilt, &tilt_azimuth);

		if (antenna_downtilt != 99.0) {    // If Cmdline override
			if (antenna_dt_direction == -1)  // dt_dir not specified
				tilt_azimuth = rotation;       // use rotation value
			mechanical_tilt = (float)antenna_downtilt;
		}

		if (antenna_dt_direction != -1)  // If Cmdline override
			tilt_azimuth = (float)antenna_dt_direction;

		spdlog::debug("Antenna Pattern Mechamical Downtilt = {}", mechanical_tilt);
		spdlog::debug("Antenna Pattern Mechanical Downtilt Direction = {}", tilt_azimuth);

		/* Read elevation (degrees) and corresponding
			 normalized field radiation pattern amplitude
			 (0.0 to 1.0) until EOF is reached. */

		if (fgets(string, 254, fd) == NULL) {
			// fprintf(stderr,"Ant elevation read error\n");
			// exit(0);
		}
		pointer = strchr(string, ';');

		if (pointer != NULL) *pointer = 0;

		sscanf(string, "%f %f", &elevation, &amplitude);

		while (1) {
			/* Read in normalized radiated field values
				 for every 0.01 degrees of elevation between
				 -10.0 and +90.0 degrees */

			x = (int)rintf(100.0 * (elevation + 10.0));

			if (x >= 0 && x <= 10000) {
				el_pattern[x] += amplitude;
				read_count[x]++;
			}

			if (fgets(string, 254, fd) == NULL)
				break;
			pointer = strchr(string, ';');
			if (pointer != NULL) *pointer = 0;

			sscanf(string, "%f %f", &elevation, &amplitude);
		}

		fclose(fd);

		/* Average the field values in case more than
			 one was read for each 0.01 degrees of elevation. */

		for (x = 0; x <= 10000; x++) {
			if (read_count[x] > 1) el_pattern[x] /= (float)read_count[x];
		}

		/* Interpolate between missing elevations (if
			 any) to completely fill the array and provide
			 radiated field values for every 0.01 degrees of
			 elevation. */

		last_index = -1;
		next_index = -1;

		for (x = 0; x <= 10000; x++) {
			if (read_count[x] != 0) {
				if (last_index == -1)
					last_index = x;
				else
					next_index = x;
			}

			if (last_index != -1 && next_index != -1) {
				valid1 = el_pattern[last_index];
				valid2 = el_pattern[next_index];

				span = next_index - last_index;
				delta = (valid2 - valid1) / (float)span;

				for (y = last_index + 1; y < next_index; y++) el_pattern[y] = el_pattern[y - 1] + delta;

				last_index = y;
				next_index = -1;
			}
		}

		/* Fill slant_angle[] array with offset angles based
			 on the antenna's mechanical beam tilt (if any)
			 and tilt direction (azimuth). */

		if (mechanical_tilt == 0.0) {
			for (x = 0; x <= 360; x++) slant_angle[x] = 0.0;
		}

		else {
			tilt_increment = mechanical_tilt / 90.0;

			for (x = 0; x <= 360; x++) {
				xx = (float)x;
				y = (int)rintf(tilt_azimuth + xx);

				while (y >= 360) y -= 360;

				while (y < 0) y += 360;

				if (x <= 180) slant_angle[y] = -(tilt_increment * (90.0 - xx));

				if (x > 180) slant_angle[y] = -(tilt_increment * (xx - 270.0));
			}
		}

		slant_angle[360] = slant_angle[0]; /* 360 degree wrap-around */

		for (w = 0; w <= 360; w++) {
			tilt = slant_angle[w];

			/** Convert tilt angle to
							an array index offset **/

			y = (int)rintf(100.0 * tilt);

			/* Copy shifted el_pattern[10001] field
				 values into elevation_pattern[361][1001]
				 at the corresponding azimuth, downsampling
				 (averaging) along the way in chunks of 10. */

			for (x = y, z = 0; z <= 1000; x += 10, z++) {
				for (sum = 0.0, a = 0; a < 10; a++) {
					b = a + x;

					if (b >= 0 && b <= 10000) sum += el_pattern[b];
					if (b < 0) sum += el_pattern[0];
					if (b > 10000) sum += el_pattern[10000];
				}

				elevation_pattern[w][z] = sum / 10.0;
			}
		}

		got_elevation_pattern = true;

		for (x = 0; x <= 360; x++) {
			for (y = 0; y <= 1000; y++) {
				if (got_elevation_pattern)
					elevation = elevation_pattern[x][y];
				else
					elevation = 1.0;

				if (got_azimuth_pattern)
					az = azimuth_pattern[x];
				else
					az = 1.0;

				LR.antenna_pattern[x][y] = az * elevation;
			}
		}
	}
	return 0;
}

int LoadSignalColors(struct site xmtr)
{
	int x, y, ok, val[4];
	char filename[255], string[80], *pointer = NULL, *s;
	FILE *fd = NULL;

	if (color_file != NULL && color_file[0] != 0)
		for (x = 0; color_file[x] != '.' && color_file[x] != 0 && x < 250; x++) filename[x] = color_file[x];
	else
		for (x = 0; xmtr.filename[x] != '.' && xmtr.filename[x] != 0 && x < 250; x++) filename[x] = xmtr.filename[x];

	filename[x] = '.';
	filename[x + 1] = 's';
	filename[x + 2] = 'c';
	filename[x + 3] = 'f';
	filename[x + 4] = 0;

	/* Default values */

	region.level[0] = 128;
	region.color[0][0] = 255;
	region.color[0][1] = 0;
	region.color[0][2] = 0;

	region.level[1] = 118;
	region.color[1][0] = 255;
	region.color[1][1] = 165;
	region.color[1][2] = 0;

	region.level[2] = 108;
	region.color[2][0] = 255;
	region.color[2][1] = 206;
	region.color[2][2] = 0;

	region.level[3] = 98;
	region.color[3][0] = 255;
	region.color[3][1] = 255;
	region.color[3][2] = 0;

	region.level[4] = 88;
	region.color[4][0] = 184;
	region.color[4][1] = 255;
	region.color[4][2] = 0;

	region.level[5] = 78;
	region.color[5][0] = 0;
	region.color[5][1] = 255;
	region.color[5][2] = 0;

	region.level[6] = 68;
	region.color[6][0] = 0;
	region.color[6][1] = 208;
	region.color[6][2] = 0;

	region.level[7] = 58;
	region.color[7][0] = 0;
	region.color[7][1] = 196;
	region.color[7][2] = 196;

	region.level[8] = 48;
	region.color[8][0] = 0;
	region.color[8][1] = 148;
	region.color[8][2] = 255;

	region.level[9] = 38;
	region.color[9][0] = 80;
	region.color[9][1] = 80;
	region.color[9][2] = 255;

	region.level[10] = 28;
	region.color[10][0] = 0;
	region.color[10][1] = 38;
	region.color[10][2] = 255;

	region.level[11] = 18;
	region.color[11][0] = 142;
	region.color[11][1] = 63;
	region.color[11][2] = 255;

	region.level[12] = 8;
	region.color[12][0] = 140;
	region.color[12][1] = 0;
	region.color[12][2] = 128;

	region.levels = 13;

	/* Don't save if we don't have an output file */
	if ((fd = fopen(filename, "r")) == NULL && xmtr.filename[0] == '\0') return 0;

	if (fd == NULL) {
		if ((fd = fopen(filename, "w")) == NULL) return errno;

		for (x = 0; x < region.levels; x++)
			fprintf(fd, "%3d: %3d, %3d, %3d\n", region.level[x], region.color[x][0], region.color[x][1], region.color[x][2]);

		fclose(fd);
	}

	else {
		x = 0;
		s = fgets(string, 80, fd);

		while (x < 128 && s != NULL) {
			pointer = strchr(string, ';');

			if (pointer != NULL) *pointer = 0;

			ok = sscanf(string, "%d: %d, %d, %d", &val[0], &val[1], &val[2], &val[3]);

			if (ok == 4) {
				spdlog::debug("LoadSignalColors() {}: {}, {}, {}", val[0], val[1], val[2], val[3]);

				for (y = 0; y < 4; y++) {
					if (val[y] > 255) val[y] = 255;

					if (val[y] < 0) val[y] = 0;
				}

				region.level[x] = val[0];
				region.color[x][0] = val[1];
				region.color[x][1] = val[2];
				region.color[x][2] = val[3];
				x++;
			}

			s = fgets(string, 80, fd);
		}

		fclose(fd);
		region.levels = x;
	}
	return 0;
}

int LoadLossColors(struct site xmtr)
{
	int x, y, ok, val[4];
	char filename[255], string[80], *pointer = NULL, *s;
	FILE *fd = NULL;

	if (color_file != NULL && color_file[0] != 0)
		for (x = 0; color_file[x] != '.' && color_file[x] != 0 && x < 250; x++) filename[x] = color_file[x];
	else
		for (x = 0; xmtr.filename[x] != '.' && xmtr.filename[x] != 0 && x < 250; x++) filename[x] = xmtr.filename[x];

	filename[x] = '.';
	filename[x + 1] = 'l';
	filename[x + 2] = 'c';
	filename[x + 3] = 'f';
	filename[x + 4] = 0;

	/* Default values */

	region.level[0] = 80;
	region.color[0][0] = 255;
	region.color[0][1] = 0;
	region.color[0][2] = 0;

	region.level[1] = 90;
	region.color[1][0] = 255;
	region.color[1][1] = 128;
	region.color[1][2] = 0;

	region.level[2] = 100;
	region.color[2][0] = 255;
	region.color[2][1] = 165;
	region.color[2][2] = 0;

	region.level[3] = 110;
	region.color[3][0] = 255;
	region.color[3][1] = 206;
	region.color[3][2] = 0;

	region.level[4] = 120;
	region.color[4][0] = 255;
	region.color[4][1] = 255;
	region.color[4][2] = 0;

	region.level[5] = 130;
	region.color[5][0] = 184;
	region.color[5][1] = 255;
	region.color[5][2] = 0;

	region.level[6] = 140;
	region.color[6][0] = 0;
	region.color[6][1] = 255;
	region.color[6][2] = 0;

	region.level[7] = 150;
	region.color[7][0] = 0;
	region.color[7][1] = 208;
	region.color[7][2] = 0;

	region.level[8] = 160;
	region.color[8][0] = 0;
	region.color[8][1] = 196;
	region.color[8][2] = 196;

	region.level[9] = 170;
	region.color[9][0] = 0;
	region.color[9][1] = 148;
	region.color[9][2] = 255;

	region.level[10] = 180;
	region.color[10][0] = 80;
	region.color[10][1] = 80;
	region.color[10][2] = 255;

	region.level[11] = 190;
	region.color[11][0] = 0;
	region.color[11][1] = 38;
	region.color[11][2] = 255;

	region.level[12] = 200;
	region.color[12][0] = 142;
	region.color[12][1] = 63;
	region.color[12][2] = 255;

	region.level[13] = 210;
	region.color[13][0] = 196;
	region.color[13][1] = 54;
	region.color[13][2] = 255;

	region.level[14] = 220;
	region.color[14][0] = 255;
	region.color[14][1] = 0;
	region.color[14][2] = 255;

	region.level[15] = 230;
	region.color[15][0] = 255;
	region.color[15][1] = 194;
	region.color[15][2] = 204;

	region.levels = 16;
	/*	region.levels = 120; // 240dB max PL */

	/*	for(int i=0; i<region.levels;i++){
					region.level[i] = i*2;
					region.color[i][0] = i*2;
					region.color[i][1] = i*2;
					region.color[i][2] = i*2;
			}
	*/
	/* Don't save if we don't have an output file */
	if ((fd = fopen(filename, "r")) == NULL && xmtr.filename[0] == '\0') return 0;

	if (fd == NULL) {
		if ((fd = fopen(filename, "w")) == NULL) return errno;

		for (x = 0; x < region.levels; x++)
			fprintf(fd, "%3d: %3d, %3d, %3d\n", region.level[x], region.color[x][0], region.color[x][1], region.color[x][2]);

		fclose(fd);

		spdlog::error("loadLossColors: fopen fail: {}", filename);
	}

	else {
		x = 0;
		s = fgets(string, 80, fd);

		while (x < 128 && s != NULL) {
			pointer = strchr(string, ';');

			if (pointer != NULL) *pointer = 0;

			ok = sscanf(string, "%d: %d, %d, %d", &val[0], &val[1], &val[2], &val[3]);

			if (ok == 4) {
				spdlog::debug("LoadLossColors() {}: {}, {}, {}", val[0], val[1], val[2], val[3]);

				for (y = 0; y < 4; y++) {
					if (val[y] > 255) val[y] = 255;

					if (val[y] < 0) val[y] = 0;
				}

				region.level[x] = val[0];
				region.color[x][0] = val[1];
				region.color[x][1] = val[2];
				region.color[x][2] = val[3];
				x++;
			}

			s = fgets(string, 80, fd);
		}

		fclose(fd);
		region.levels = x;
	}
	return 0;
}

int LoadDBMColors(struct site xmtr)
{
	int x, y, ok, val[4];
	char filename[255], string[80], *pointer = NULL, *s;
	FILE *fd = NULL;

	if (color_file != NULL && color_file[0] != 0)
		for (x = 0; color_file[x] != '.' && color_file[x] != 0 && x < 250; x++) filename[x] = color_file[x];
	else
		for (x = 0; xmtr.filename[x] != '.' && xmtr.filename[x] != 0 && x < 250; x++) filename[x] = xmtr.filename[x];

	filename[x] = '.';
	filename[x + 1] = 'd';
	filename[x + 2] = 'c';
	filename[x + 3] = 'f';
	filename[x + 4] = 0;

	/* Default values */

	region.level[0] = 0;
	region.color[0][0] = 255;
	region.color[0][1] = 0;
	region.color[0][2] = 0;

	region.level[1] = -10;
	region.color[1][0] = 255;
	region.color[1][1] = 128;
	region.color[1][2] = 0;

	region.level[2] = -20;
	region.color[2][0] = 255;
	region.color[2][1] = 165;
	region.color[2][2] = 0;

	region.level[3] = -30;
	region.color[3][0] = 255;
	region.color[3][1] = 206;
	region.color[3][2] = 0;

	region.level[4] = -40;
	region.color[4][0] = 255;
	region.color[4][1] = 255;
	region.color[4][2] = 0;

	region.level[5] = -50;
	region.color[5][0] = 184;
	region.color[5][1] = 255;
	region.color[5][2] = 0;

	region.level[6] = -60;
	region.color[6][0] = 0;
	region.color[6][1] = 255;
	region.color[6][2] = 0;

	region.level[7] = -70;
	region.color[7][0] = 0;
	region.color[7][1] = 208;
	region.color[7][2] = 0;

	region.level[8] = -80;
	region.color[8][0] = 0;
	region.color[8][1] = 196;
	region.color[8][2] = 196;

	region.level[9] = -90;
	region.color[9][0] = 0;
	region.color[9][1] = 148;
	region.color[9][2] = 255;

	region.level[10] = -100;
	region.color[10][0] = 80;
	region.color[10][1] = 80;
	region.color[10][2] = 255;

	region.level[11] = -110;
	region.color[11][0] = 0;
	region.color[11][1] = 38;
	region.color[11][2] = 255;

	region.level[12] = -120;
	region.color[12][0] = 142;
	region.color[12][1] = 63;
	region.color[12][2] = 255;

	region.level[13] = -130;
	region.color[13][0] = 196;
	region.color[13][1] = 54;
	region.color[13][2] = 255;

	region.level[14] = -140;
	region.color[14][0] = 255;
	region.color[14][1] = 0;
	region.color[14][2] = 255;

	region.level[15] = -150;
	region.color[15][0] = 255;
	region.color[15][1] = 194;
	region.color[15][2] = 204;

	region.levels = 16;

	/* Don't save if we don't have an output file */
	if ((fd = fopen(filename, "r")) == NULL && xmtr.filename[0] == '\0') return 0;

	if (fd == NULL) {
		if ((fd = fopen(filename, "w")) == NULL) return errno;

		for (x = 0; x < region.levels; x++)
			fprintf(fd, "%+4d: %3d, %3d, %3d\n", region.level[x], region.color[x][0], region.color[x][1], region.color[x][2]);

		fclose(fd);
	}

	else {
		x = 0;
		s = fgets(string, 80, fd);

		while (x < 128 && s != NULL) {
			pointer = strchr(string, ';');

			if (pointer != NULL) *pointer = 0;

			ok = sscanf(string, "%d: %d, %d, %d", &val[0], &val[1], &val[2], &val[3]);

			if (ok == 4) {
				spdlog::debug("LoadDBMColors() {}: {}, {}, {}", val[0], val[1], val[2], val[3]);

				if (val[0] < -200) val[0] = -200;

				if (val[0] > +40) val[0] = +40;

				region.level[x] = val[0];

				for (y = 1; y < 4; y++) {
					if (val[y] > 255) val[y] = 255;

					if (val[y] < 0) val[y] = 0;
				}

				region.color[x][0] = val[1];
				region.color[x][1] = val[2];
				region.color[x][2] = val[3];
				x++;
			}

			s = fgets(string, 80, fd);
		}

		fclose(fd);
		region.levels = x;
	}
	return 0;
}

/**
 * Load a single Copernicus DSM GeoTIFF COG tile into a free dem[] page.
 *
 * @param tile_lat  Southern latitude of the tile (e.g. 44 for 44N–45N)
 * @param tile_lon  East-positive western edge of tile (e.g. 4 for 4E–5E,
 *                  -73 for 73W–72W). min_lon = tile_lon, max_lon = tile_lon + 1.
 *
 * Filename built from DEM_path:
 *   ippd==3600 → Copernicus_DSM_COG_10_N##_00_?###_00_DEM.tif
 *   otherwise  → Copernicus_DSM_COG_30_N##_00_?###_00_DEM.tif
 *
 * Returns 1 on success, 0 if already loaded, negative errno on error.
 */
static std::mutex copernicus_mutex;

int LoadCopernicus(int tile_lat, int tile_lon)
{
    /* Build the Copernicus filename. */
    char ew = (tile_lon >= 0) ? 'E' : 'W';
    int  lon_abs = abs(tile_lon);
    char ns = (tile_lat >= 0) ? 'N' : 'S';
    int  lat_abs = abs(tile_lat);
    const char *res_str = (ippd == 3600) ? "10" : "30";

    char filename[64];
    snprintf(filename, sizeof(filename),
             "Copernicus_DSM_COG_%s_%c%02d_00_%c%03d_00_DEM.tif",
             res_str, ns, lat_abs, ew, lon_abs);

    /* Try current working directory first, then DEM_path */
    char path_plus_name[PATH_MAX];
    strncpy(path_plus_name, filename, sizeof(path_plus_name) - 1);
    path_plus_name[sizeof(path_plus_name) - 1] = '\0';

    GDALDatasetH ds = GDALOpen(path_plus_name, GA_ReadOnly);
    if (ds == NULL) {
        snprintf(path_plus_name, sizeof(path_plus_name), "%s%s", DEM_path, filename);
        ds = GDALOpen(path_plus_name, GA_ReadOnly);
    }
    if (ds == NULL) {
        spdlog::debug("LoadCopernicus: file not found: {}", filename);
        return -ENOENT;
    }

    spdlog::debug("LoadCopernicus: loading \"{}\"...", path_plus_name);

    GDALRasterBandH band = GDALGetRasterBand(ds, 1);
    int nodata_valid = 0;
    double nodata_val = GDALGetRasterNoDataValue(band, &nodata_valid);

    int src_x = GDALGetRasterXSize(ds);
    int src_y = GDALGetRasterYSize(ds);

    float *buf = new float[ippd * ippd];
    CPLErr err = (CPLErr)GDALRasterIO(band, GF_Read,
                                      0, 0, src_x, src_y,
                                      buf, ippd, ippd,
                                      GDT_Float32, 0, 0);
    GDALClose(ds);

    if (err != CE_None) {
        spdlog::error("LoadCopernicus: RasterIO failed for {}", path_plus_name);
        delete[] buf;
        return -EIO;
    }

    /* Compute global pixel offset for this tile.
     * x increases northward: tile row 0 (south edge) maps to gx_base.
     * y increases westward:  tile col 0 (west  edge) maps to gy_base + (ippd-1).
     *
     * GeoTIFF layout: row 0 = north, col 0 = west.
     * Flat array:
     *   x = gx_base + (ippd-1-r)   — row r=0 (north) → highest x in tile
     *   y = gy_base + (ippd-1-c)   — col c=0 (west)  → highest y in tile
     */
    int gx_base = (tile_lat - dem_min_lat) * ippd;
    int gy_base = (tile_lon - dem_min_lon) * ippd;

    int tile_min_el = 32768, tile_max_el = -32768;

    for (int r = 0; r < ippd; r++) {
        int gx = gx_base + (ippd - 1 - r);
        for (int c = 0; c < ippd; c++) {
            int gy = gy_base + c;
            float fval = buf[r * ippd + c];
            short val = (nodata_valid && fval == (float)nodata_val) ? 0 : (short)roundf(fval);

            dem_data[gx][gy]   = val;
            dem_signal[gx][gy] = 0;

            if (val > tile_max_el) tile_max_el = val;
            if (val < tile_min_el) tile_min_el = val;
        }
    }

    delete[] buf;

    {
        std::lock_guard<std::mutex> lock(copernicus_mutex);

        if (tile_min_el < min_elevation) min_elevation = tile_min_el;
        if (tile_max_el > max_elevation) max_elevation = tile_max_el;

        float f_max_north = (float)(tile_lat + 1);
        float f_min_north = (float)tile_lat;
        float f_max_lon   = (float)(tile_lon + 1);
        float f_min_lon   = (float)tile_lon;

        if (max_north == -90 || f_max_north > max_north) max_north = f_max_north;
        if (min_north ==  90 || f_min_north < min_north) min_north = f_min_north;
        if (f_max_lon > max_lon) max_lon = f_max_lon;
        if (f_min_lon < min_lon) min_lon = f_min_lon;
    }

    spdlog::info("LoadCopernicus: loaded {} (el {}/{}m, bounds {:.0f}N {:.0f}E → {:.0f}N {:.0f}E)",
                 filename, tile_min_el, tile_max_el,
                 (float)tile_lat, (float)tile_lon,
                 (float)(tile_lat + 1), (float)(tile_lon + 1));

    return 1;
}

/**
 * Load the required Topo data for the given lat/lon box
*/
int LoadTopoData(bbox region)
{
    spdlog::info("Loading topo data for boundaries: ({:.6f}N, {:.6f}E) to ({:.6f}N, {:.6f}E)",
        region.lower_right.lat,
        region.lower_right.lon,
        region.upper_left.lat,
        region.upper_left.lon
    );

    // TODO: we don't handle loading data around 0/360 W

    // Get the nearest whole number lat/lons based on coords
    int r_min_lat = floor(region.lower_right.lat);
    int r_max_lat = ceil(region.upper_left.lat);
    int r_min_lon = floor(region.lower_right.lon);
    int r_max_lon = ceil(region.upper_left.lon);

    int tiles_lat = r_max_lat - r_min_lat;
    int tiles_lon = r_max_lon - r_min_lon;

    // Sanity check
    if (!tiles_lat || !tiles_lon) {
        spdlog::error("Our plot area gave us {} x {} tiles which is invalid!", tiles_lat, tiles_lon);
        exit(1);
    }

    /* Allocate flat DEM arrays now that the bounding box is known */
    alloc_dem(r_min_lat, r_min_lon, tiles_lat, tiles_lon);

    // Load the data
    for (int x = 0; x < tiles_lon; x++) {
        for (int y = 0; y < tiles_lat; y++) {
            int tile_lon = r_min_lon + x;
            int tile_lat = r_min_lat + y;
            spdlog::debug("Loading topo for tile {}N {}E to {}N {}E", tile_lat, tile_lon, tile_lat + 1, tile_lon + 1);
            int success = LoadCopernicus(tile_lat, tile_lon);
            if (success < 0 && success != -ENOENT) {
                return -success;
            }
        }
    }

	return 0;
}

