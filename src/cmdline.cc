/* GDAL must come before any header that defines MAX */
#include <gdal.h>

#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <thread>

#include <spdlog/spdlog.h>

#include "cmdline.hh"
#include "main.hh"
#include "common.hh"
#include "inputs.hh"
#include "models/los.hh"

// Parse command-line options, apply defaults, and load optional antenna files.
void parse_cmdline(int argc, char *argv[])
{
    int x, y, z = 0;

    ppa            = false;
    normalise      = 0;
    altitudeLR     = 1;
    prop_model     = ITM_LR;
    number_threads = std::thread::hardware_concurrency();
    mapfile.clear();

    char antenna_file[255];
    char *az_filename = NULL, *el_filename = NULL;

    spdlog::info("Version {} ({} {})", VERSION, GIT_BRANCH, GIT_COMMIT_HASH);
    spdlog::info("    Compile date: {} {}", __DATE__, __TIME__);
    spdlog::info("");

    if (argc == 1) {
        spdlog::info("License: GNU General Public License (GPL) version 2");
        spdlog::info("");
        spdlog::info("Radio propagation simulator by Alex Farrant QCVS, 2E0TDW");
        spdlog::info("Based upon SPLAT! by John Magliacane, KD2BD");
        spdlog::info("Some feature enhancements/additions by Aaron A. Collins, N9OZB");
        spdlog::info("Additional improvements and multithreading fixes by P. McDonnell, W3AXL");
        spdlog::info("Additional improvements and cleanup by Loris, F4ITL");
        spdlog::info("");
        spdlog::info("Usage: signalserver [data options] [input options] [antenna options] [output options] -o outputfile");
        spdlog::info("");
        spdlog::info("Data:");
        spdlog::info("     -dem Directory containing Copernicus DEM GeoTIFF COG tiles");
        spdlog::info("                  (Copernicus_DSM_COG_30_N##_00_?###_00_DEM.tif for 1200 ppd,");
        spdlog::info("                   Copernicus_DSM_COG_10_N##_00_?###_00_DEM.tif for 3600 ppd)");
        spdlog::info("     -color Color palette: heat (default), jet, turbo, viridis, magma, plasma, inferno, hot, parula, gray, hsv, cubehelix, cividis, github");
        spdlog::info("Input:");
        spdlog::info("     -lat Tx Latitude (decimal degrees) -70/+70");
        spdlog::info("     -lon Tx Longitude (decimal degrees) -180/+180");
        spdlog::info("     -rla (Optional) Rx Latitude for PPA (decimal degrees) -70/+70");
        spdlog::info("     -rlo (Optional) Rx Longitude for PPA (decimal degrees) -180/+180");
        spdlog::info("     -f Tx Frequency (MHz) 20MHz to 100GHz (LOS after 20GHz)");
        spdlog::info("     -erp Tx Total Effective Radiated Power in Watts (dBd) inc Tx+Rx gain. 2.14 dBi = 0 dBd");
        spdlog::info("     -gc Random ground clutter (meters)");
        spdlog::info("     -te Terrain code 1-6 (optional - 1. Water, 2. Marsh, 3. Farmland,");
        spdlog::info("          4. Mountain, 5. Desert, 6. Urban");
        spdlog::info("     -terdic Terrain dielectric value 2-80 (optional)");
        spdlog::info("     -tercon Terrain conductivity 0.01-0.0001 (optional)");
        spdlog::info("     -cl Climate code 1-7 (optional - 1. Equatorial 2. Continental subtropical");
        spdlog::info("          3. Maritime subtropical 4. Desert 5. Continental temperate");
        spdlog::info("          6. Maritime temperate (Land) 7. Maritime temperate (Sea)");
        spdlog::info("     -rel Reliability for ITM model (% of 'time') 1 to 99 (optional, default 50%)");
        spdlog::info("     -conf Confidence for ITM model (% of 'situations') 1 to 99 (optional, default 50%)");
        spdlog::info("     -number_threads Number of worker threads to divide the plot rectangle into (must be even and >= 4)");
        spdlog::info("     -hd Use HD mode (30m), per defaut 90m");
        spdlog::info("Output:");
        spdlog::info("     -o basename (Output file basename - required, min 5 chars)");
        spdlog::info("     -dbm Plot Rxd signal power instead of field strength in dBuV/m");
        spdlog::info("     -rt Rx Threshold (dB / dBm / dBuV/m)");
        spdlog::info("     -R Radius (kilometers)");
        spdlog::info("     -pm Propagation model. 1: ITM, 2: LOS, 3: Hata, 4: ECC33,");
        spdlog::info("          5: SUI, 6: COST-Hata, 7: FSPL, 8: ITWOM, 9: Ericsson,");
        spdlog::info("          10: Plane earth, 11: Egli VHF/UHF, 12: Soil, 13: NTIA ITM");
        spdlog::info("     -pe Propagation model mode: 1=Urban,2=Suburban,3=Rural");
        spdlog::info("     -ked Knife edge diffraction (Already on for ITM)");
        spdlog::info("     -geotiff Output a geotiff file");
        spdlog::info("Antenna:");
        spdlog::info("     -ant (antenna pattern file basename+path for .az and .el files)");
        spdlog::info("     -txh Tx Height (above ground)");
        spdlog::info("     -rxh Rx Height(s) (optional. Default=1)");
        spdlog::info("     -rxg Rx gain dBd (optional for PPA text report)");
        spdlog::info("     -hp Horizontal Polarisation (default=vertical)");
        spdlog::info("     -rot  (  0.0 - 359.0 degrees, default 0.0) Antenna Pattern Rotation");
        spdlog::info("     -dt   ( -10.0 - 90.0 degrees, default 0.0) Antenna Downtilt");
        spdlog::info("     -dtdir ( 0.0 - 359.0 degrees, default 0.0) Antenna Downtilt Direction");
        spdlog::info("Debugging:");
        spdlog::info("     -t Terrain greyscale background");
        spdlog::info("     -dbg Verbose debug messages");
        spdlog::info("     -ng Normalise Path Profile graph");

        exit(1);
    }

    GDALAllRegister();

    y = argc - 1;
    dbm = false;
    DEM_path.clear();
    clutter = 0.0;
    path.length = 0;
    fzone_clearance = 0.6;
    contour_threshold = 0;
    max_range = 1000.0; // meters
    ngs = true;			// no greyscale background

    // Defaults
    LR.eps_dielect = 15.0;	// Farmland
    LR.sgm_conductivity = 0.005;	// Farmland
    LR.eno_ns_surfref = 301.0;
    LR.frq_mhz = 19.0;	// Deliberately too low
    LR.radio_climate = 5;	// continental
    LR.pol = 1;		// vert
    LR.conf = 0.50;
    LR.rel = 0.50;
    LR.erp = 0.0;		// will default to Path Loss

    antenna_rotation = -1;  // unique defaults to test usage
    antenna_downtilt = 99.0; // don't mess with them!
    antenna_dt_direction = -1;
    antenna_file[0] = '\0';

    tx_site.lat = 91.0;
    tx_site.lon = 181.0;
    rx_site.lat = 91.0;
    rx_site.lon = 181.0;

    /* Scan for command line arguments */

    for (x = 1; x <= y; x++) {

        if (strcmp(argv[x], "-R") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0] && argv[z][0] != '-') {
                sscanf(argv[z], "%lf", &max_range);
                max_range *= 1000.0;
            }
        }

        if (strcmp(argv[x], "-gc") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0] && argv[z][0] != '-') {
                sscanf(argv[z], "%lf", &clutter);

                if (clutter < 0.0)
                    clutter = 0.0;


            }
        }

        if (strcmp(argv[x], "-ant") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0] && argv[z][0] != '-') {
                strncpy(antenna_file, argv[z], 253);
            }
        }

        if (strcmp(argv[x], "-rot") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0] && argv[z][0] != '-') {
                sscanf(argv[z], "%lf", &antenna_rotation);

                if (antenna_rotation < 0.0)
                    antenna_rotation = 0.0;
                if (antenna_rotation > 359.0)
                    antenna_rotation = 0.0;
            }
        }

        if (strcmp(argv[x], "-dt") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0]) {	/* A minus argument is legal here */
                sscanf(argv[z], "%lf", &antenna_downtilt);
                if (antenna_downtilt < -10.0)
                    antenna_downtilt = -10.0;
                if (antenna_downtilt > 90.0)
                    antenna_downtilt = 90.0;
            }
        }

        if (strcmp(argv[x], "-dtdir") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0] && argv[z][0] != '-') {
                sscanf(argv[z], "%lf", &antenna_dt_direction);

                if (antenna_dt_direction < 0.0)
                    antenna_dt_direction = 0.0;
                if (antenna_dt_direction > 359.0)
                    antenna_dt_direction = 0.0;
            }
        }

        if (strcmp(argv[x], "-o") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0] && argv[z][0] != '-') {
                // sanity check length
                if(strlen(argv[z]) < 5){
                    spdlog::error("Output name is too short. Must be at least 5 chars");
                    exit(1);
                }

                mapfile = argv[z];
                output_filename = argv[z];
            } else if (z <= y && argv[z][0] && argv[z][0] == '-' && argv[z][1] == '\0' ) {
                /* default file name */
                output_filename = "output";
            }
        }

        if (strcmp(argv[x], "-rt") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0])	/* A minus argument is legal here */
                sscanf(argv[z], "%d", &contour_threshold);
        }

        if (strcmp(argv[x], "-hd") == 0) {
            spdlog::info("    hd mode");
            ppd = 3600;

            spdlog::info("    Built for {} ppd", ppd);
        }

        if (strcmp(argv[x], "-t") == 0) {
            ngs = false;	// greyscale background
        }

        if (strcmp(argv[x], "-dbm") == 0)
            dbm = true;

        if (strcmp(argv[x], "-geotiff") == 0) {
            geotiff = true;
        }

        if (strcmp(argv[x], "-dem") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0] && argv[z][0] != '-')
                DEM_path = argv[z];
        }

        if (strcmp(argv[x], "-lat") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0]) {
                tx_site.lat = atof(argv[z]);
            }
        }
        if (strcmp(argv[x], "-lon") == 0) {
            z = x + 1;
            if (z <= y && argv[z][0]) {
                tx_site.lon = atof(argv[z]);
            }
        }
        //Switch to Path Profile Mode if Rx co-ords specified
        if (strcmp(argv[x], "-rla") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0]) {
                ppa = true;
                rx_site.lat = atof(argv[z]);

            }
        }
        if (strcmp(argv[x], "-rlo") == 0) {
            z = x + 1;
            if (z <= y && argv[z][0]) {
                rx_site.lon = atof(argv[z]);
            }
        }

        if (strcmp(argv[x], "-txh") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0] && argv[z][0] != '-') {
                sscanf(argv[z], "%lf", &tx_site.alt);

            }
        }

        if (strcmp(argv[x], "-rxh") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0] && argv[z][0] != '-') {
                sscanf(argv[z], "%lf", &altitudeLR);
                sscanf(argv[z], "%lf", &rx_site.alt);
            }
        }

        if (strcmp(argv[x], "-rxg") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0] && argv[z][0] != '-') {
                sscanf(argv[z], "%lf", &rxGain);
            }
        }

        if (strcmp(argv[x], "-f") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0] && argv[z][0] != '-') {
                sscanf(argv[z], "%lf", &LR.frq_mhz);
            }
        }

        if (strcmp(argv[x], "-erp") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0] && argv[z][0] != '-') {
                sscanf(argv[z], "%lf", &LR.erp);
            }
        }

        if (strcmp(argv[x], "-cl") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0] && argv[z][0] != '-') {

                sscanf(argv[z], "%d", &LR.radio_climate);

            }
        }
        if (strcmp(argv[x], "-te") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0] && argv[z][0] != '-') {
                int ter;
                sscanf(argv[z], "%d", &ter);

                switch (ter) {
                case 1:	// Water
                    terdic = 80;
                    tercon = 0.010;
                    break;

                case 2:	// Marsh
                    terdic = 12;
                    tercon = 0.007;
                    break;

                case 3:	// Farmland
                    terdic = 15;
                    tercon = 0.005;
                    break;

                case 4:	//Mountain
                    terdic = 13;
                    tercon = 0.002;
                    break;
                case 5:	//Desert
                    terdic = 13;
                    tercon = 0.002;
                    break;
                case 6:	//Urban
                    terdic = 5;
                    tercon = 0.001;
                    break;
                }
                LR.eps_dielect = terdic;
                LR.sgm_conductivity = tercon;

            }
        }

        if (strcmp(argv[x], "-terdic") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0] && argv[z][0] != '-') {

                sscanf(argv[z], "%lf", &terdic);

                LR.eps_dielect = terdic;

            }
        }

        if (strcmp(argv[x], "-tercon") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0] && argv[z][0] != '-') {

                sscanf(argv[z], "%lf", &tercon);

                LR.sgm_conductivity = tercon;

            }
        }

        if (strcmp(argv[x], "-hp") == 0) {
            // Horizontal polarisation (0)
            LR.pol = 0;
        }

        if (strcmp(argv[x], "-dbg") == 0) {
            debug = 1;
        }

        /*Prop model */

        if (strcmp(argv[x], "-pm") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0]) {
                int temp;
                sscanf(argv[z], "%d", &temp);
                prop_model = (PropModel)temp;
            }
        }
        // Prop model variant eg. urban/suburban
        if (strcmp(argv[x], "-pe") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0]) {
                sscanf(argv[z], "%d", &pmenv);
            }
        }
        //Knife edge diffraction
        if (strcmp(argv[x], "-ked") == 0) {
            knifeedge = true;
        }

        //Normalise Path Profile chart
        if (strcmp(argv[x], "-ng") == 0) {
            normalise = 1;
        }

        // Reliability % for ITM model
        if (strcmp(argv[x], "-rel") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0]) {
                sscanf(argv[z], "%lf", &LR.rel);
                LR.rel=LR.rel/100;
            }
        }
        // Confidence % for ITM model
        if (strcmp(argv[x], "-conf") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0]) {
                sscanf(argv[z], "%lf", &LR.conf);
                LR.conf=LR.conf/100;
            }
        }
        if (strcmp(argv[x], "-color") == 0) {
            z = x + 1;
            if (z <= y && argv[z][0] && argv[z][0] != '-')
                color_palette = argv[z];
        }

        // number_threads to divide plot by
        if (strcmp(argv[x], "-number_threads") == 0) {
            z = x + 1;

            if (z <= y && argv[z][0]) {
                sscanf(argv[z], "%d", &number_threads);
            }
        }
    }

    if (debug) {
        spdlog::set_level(spdlog::level::debug);
        spdlog::debug("Debug logging enabled");
    } else {
        spdlog::set_level(spdlog::level::info);
    }

    /* Load antenna pattern files now that all arguments have been parsed.
     * Use -ant basename if provided, otherwise fall back to the output basename. */
    if (!output_filename.empty()) {
        const char *az_base = (antenna_file[0] != '\0') ? antenna_file : output_filename.c_str();
        size_t az_needed = strlen(az_base) + 3 + 1;

        az_filename = (char*)calloc(az_needed, sizeof(char));
        if (az_filename == NULL) exit(ENOMEM);
        el_filename = (char*)calloc(az_needed, sizeof(char));
        if (el_filename == NULL) { free(az_filename); exit(ENOMEM); }

        snprintf(az_filename, az_needed, "%s.az", az_base);
        snprintf(el_filename, az_needed, "%s.el", az_base);

        int result;
        if ((result = LoadPAT(az_filename, el_filename)) != 0) {
            spdlog::error("Permissions error reading antenna pattern file");
            free(az_filename);
            free(el_filename);
            exit(result);
        }
        free(az_filename);
        free(el_filename);
    }

    /* ERROR DETECTION */
    if (tx_site.lat > 90 || tx_site.lat < -90) {
        spdlog::error("Either the lat was missing or out of range!");
        exit(EINVAL);

    }
    if (tx_site.lon > 180.0 || tx_site.lon < -180.0) {
        spdlog::error("Either the lon was missing or out of range! (expected -180 to +180)");
        exit(EINVAL);

    }
    if (LR.frq_mhz < 20 || LR.frq_mhz > 100000) {
        spdlog::error("Either the Frequency was missing or out of range!");
        exit(EINVAL);
    }
    if (LR.erp > 500000000) {
        spdlog::error("Power was out of range!");
        exit(EINVAL);

    }
    if (LR.eps_dielect > 80 || LR.eps_dielect < 0.1) {
        spdlog::error("Ground Dielectric value out of range!");
        exit(EINVAL);

    }
    if (LR.sgm_conductivity > 0.01 || LR.sgm_conductivity < 0.000001) {
        spdlog::error("Ground conductivity out of range!");
        exit(EINVAL);

    }

    if (tx_site.alt < 0 || tx_site.alt > 60000) {
        spdlog::error("Tx altitude above ground was too high: {}",
            tx_site.alt);
        exit(EINVAL);
    }
    if (altitudeLR < 0 || altitudeLR > 60000) {
        spdlog::error("Rx altitude above ground was too high!");
        exit(EINVAL);
    }

    if (ppd < 300 || ppd > 10000) {
        spdlog::error("resolution out of range!");
        exit(EINVAL);
    }

    if (contour_threshold < -200 || contour_threshold > 240) {
        spdlog::error("Receiver threshold out of range (-200 / +240)");
        exit(EINVAL);
    }
    if (prop_model > 2 && prop_model < 7 && LR.frq_mhz < 150) {
        spdlog::error("Frequency too low for Propagation model");
        exit(EINVAL);
    }

    /* Ensure a trailing '/' is present in DEM_path */

    if (!DEM_path.empty() && DEM_path.back() != '/') {
        spdlog::debug("Appending / to Copernicus directory");
        DEM_path += '/';
    }

    if (number_threads < 1) {
        spdlog::error("Number of worker threads must be >= 1");
        exit(EINVAL);
    }

    spdlog::info("-------------------------------- Plot Information --------------------------------");
    spdlog::info("    TX site parameters: {:.6f}N, {:.6f}E, {:.0f} m AGL", tx_site.lat, tx_site.lon, tx_site.alt);
    spdlog::info("    Plot parameters: {:.2f}-km radius, resolution of {} ppd", max_range/1000, ppd);
    spdlog::info("    Model parameters: {} MHz at {} W EIRP (dBd), {}% confidence", LR.frq_mhz, LR.erp, (uint8_t)(LR.conf * 100));
    spdlog::info("    Worker threads: {}", number_threads);
    spdlog::info("");
    spdlog::info("    Directories:");
    spdlog::info("        DEM: {}", DEM_path);
    spdlog::info(VERT_SEP);
}
