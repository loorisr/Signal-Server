# Signal Server

This is a fork of [W3AXL Signal-Server](https://github.com/W3AXL/Signal-Server).

Thanks to all contributors of SPLAT! and Signal Server!

Developped to be used with [this API](https://github.com/loorisr/splat-API)

It is optimized for speed for area calculations.

Improvements:
- speed: **6 times faster** than the original code (1.15s vs 6.98s for a 100 km radius area calculation).
- use the total number of cores of a computer. Can be overrided by `-number_threads X`
- better compilation flags
- `-geotiff` for geotiff image output
- `-hd` for HD mode (30m resolution)
- use [Copernicus DEM](https://dataspace.copernicus.eu/explore-data/data-collections/copernicus-contributing-missions/collections-description/COP-DEM) (GLO-30 or GLO-90) DEM (more accurate than SRTM) https://registry.opendata.aws/copernicus-dem/
- simplification: remove of useless vars, functions
- No LIDAR mode
- No more sdf for DEM, no more strm2sdf
- better code lisibility, comments
- several color palettes with -color
- add NTIA official ITM code
- generate geotiff only (no more heavy ppm files)
- `-dh` to select Delta-H computation method for ITM (default: ITM legacy code; or `calculate_delta_h_adjustable` with N interpolated points, no interpolation, recommended >=50, allows for 30-40% faster calculation)
- `-fast [N]` to speed up ITM by recalculating Delta-H every N points after the first N points from TX on each radial, with linear interpolation between points (default when enabled: `N=10`)


# Signal Server
Multi-threaded radio propagation simulator based upon SPLAT! by Alex Farrant QCVS, 2E0TDW and P. McDonnell, W3AXL. 

SPLAT! Project started in 1997 by John A. Magliacane, KD2BD

Some additional features and fixes by Aaron A. Collins, N9OZB
                                      Tom Hayward, KD7LXL
                                      Darcy Buskermolen, VA7DBI
                                      McDonnell, W3AXL

This server application will generate RF coverage predictions, producing either 2D profile plots (Point-to-Point) or 360 degree polar plots in WGS-84 projection. 

## Requirements
* C++14-conformant C++ compiler (GCC,G++ / clang)
* Build environment for C++ (linker, C++ Standard Library and so forth) 
* CMake v3.5 or newer
* library pthread: POSIX threads library
* library dl: Open and close a shared object - POSIX conform
* library z: zlib is a general-purpose lossless data-compression library
* Multicore CPU (optional but recomended)
* ~2GB Memory
* Copernicus terrain tile(s)

Signal Server is a very resource intensive multicore application. Only publish it for common use if you know what you are doing and you are advised to wrap it with another script to perform input validation.

## File extensions and types used by signalserver:
```
.az SPLAT! antenna pattern azimuth data file
.el SPLAT! antenna pattern elevation data file
```

## Install dependencies
Assuming Debian / Ubuntu, this will fetch the core libraries needed to build it as well as an image library for manipulating outputs.
```
sudo apt-get install g++ cmake spdlog
```

## Installation
Change into the source directory to build the binaries.
```
mkdir build
cd build
cmake ../src
make
```

## Parameters
```
License: GNU General Public License (GPL) version 2

Radio propagation simulator by Alex Farrant QCVS, 2E0TDW
Based upon SPLAT! by John Magliacane, KD2BD
Some feature enhancements/additions by Aaron A. Collins, N9OZB
Additional improvements and multithreading fixes by P. McDonnell, W3AXL
Additional improvements and cleanup by Loris, F4ITL

Usage: signalserver [data options] [input options] [antenna options] [output options] -o outputfile

Data:
     -dem Directory containing Copernicus DEM GeoTIFF COG tiles
                  (Copernicus_DSM_COG_30_N##_00_?###_00_DEM.tif for 1200 ppd,
                   Copernicus_DSM_COG_10_N##_00_?###_00_DEM.tif for 3600 ppd)
     -color Color palette: heat (default), jet, turbo, viridis, magma, plasma, inferno, hot, parula, gray, hsv, cubehelix, cividis, github
Input:
     -lat Tx Latitude (decimal degrees) -70/+70
     -lon Tx Longitude (decimal degrees) -180/+180
     -rla (Optional) Rx Latitude for PPA (decimal degrees) -70/+70
     -rlo (Optional) Rx Longitude for PPA (decimal degrees) -180/+180
     -f Tx Frequency (MHz) 20MHz to 100GHz (LOS after 20GHz)
     -erp Tx Total Effective Radiated Power in Watts (dBd) inc Tx+Rx gain. 2.14 dBi = 0 dBd
     -gc Random ground clutter in meters (optional, default 0)
     -te Terrain code 1-6 (optional - 1. Water, 2. Marsh, 3. Farmland,
          4. Mountain, 5. Desert, 6. Urban
     -terdic Terrain dielectric value 2-80 (optional, default 15)
     -tercon Terrain conductivity 0.01-0.0001 (optional, default 0.005)
     -cl Climate code 1-7 (optional - 1. Equatorial 2. Continental subtropical
          3. Maritime subtropical 4. Desert 5. Continental temperate
          6. Maritime temperate (Land) 7. Maritime temperate (Sea)
          (optional, default 5. Continental temperate)
     -rel Reliability for ITM model (% of 'time') 1 to 99 (optional, default 50%)
     -conf Confidence for ITM model (% of 'situations') 1 to 99 (optional, default 50%)
     -number_threads Number of worker threads to divide the plot rectangle into (default = CPU core numbers)
     -hd Use HD mode (1'=30m), (optional, default 3'=90m)
     -dh Delta-H method for ITM: N points for calculate_delta_h_adjustable (0=no interpolation). Default: itm legacy code.
     -fast [N] Recompute ITM delta-h every N points after the first N points from TX on each radial (ITM only), using linear interpolation in between. Default N=10
Output:
     -o basename (Output file basename - required, min 5 chars)
     -dbm Plot Rxd signal power (in dBm)instead of field strength (in dBuV/m)
     -rt Rx Threshold (dBm / dBuV/m)
     -R Radius (kilometers)
     -pm Propagation model. 1: ITM, 2: LOS, 3: Hata, 4: ECC33,
          5: SUI, 6: COST-Hata, 7: FSPL, 8: ITWOM, 9: Ericsson,
          10: Plane earth, 11: Egli VHF/UHF, 12: Soil, 13: NTIA ITM
     -pe Propagation model mode: 1=Urban,2=Suburban,3=Rural
     -ked Knife edge diffraction (Already on for ITM)
     -geotiff Output a geotiff file
Antenna:
     -ant (antenna pattern file basename+path for .az and .el files)
     -txh Tx Height (above ground)
     -rxh Rx Height(s) (optional. Default=1)
     -rxg Rx gain dBd (optional for PPA text report)
     -hp Horizontal Polarisation (default=vertical)
     -rot  (  0.0 - 359.0 degrees, default 0.0) Antenna Pattern Rotation
     -dt   ( -10.0 - 90.0 degrees, default 0.0) Antenna Downtilt
     -dtdir ( 0.0 - 359.0 degrees, default 0.0) Antenna Downtilt Direction
Debugging:
     -t Terrain greyscale background
     -dbg Verbose debug messages
     -ng Normalise Path Profile graph
```

### REFERENCE DATA
Signal server is designed for most of the environments and climates on Planet Earth but Polar region support is limited above extreme latitudes. (Svalbard is ok).

It can run with or without terrain data and can even be used to simulate radiation of other EM emissions like light.

#### -dem
##### Directory containing Digital Elevation Models (DEM)
Use Copernicus DEM GeoTIFF COG tiles. The parser expects filenames such as `Copernicus_DSM_COG_30_N##_00_?###_00_DEM.tif` for 1200 ppd and `Copernicus_DSM_COG_10_N##_00_?###_00_DEM.tif` for 3600 ppd.


### Antenna radiation pattern(s)
Antenna  pattern  data is read from a pair of files having the same base name as the output file (-o), but with  .az  and .el  extensions  for azimuth and elevation patterns.  This name can be overridden with the command line option -ant followed by a base filename or path and base filename for the .az and .el files.  The program will first search for command line specified antenna files with -ant, then it will default to searching for antenna files having the same name as the output file.  
```
045.0
0       0.8950590
1       0.8966406
2       0.8981447
3       0.8995795
4       0.9009535
5       0.9022749
```
The  first  line of the .az file specifies the direction measured clockwise in degrees from  True  North.  This is followed by azimuth headings (0 to 360 degrees) and their associated normalized field patterns (0.000 to 1.000) separated by whitespace.  This direction can be overridden by the command line option -rot followed by a positive direction value from 0 to 360 degrees.  This value will override any direction specified in the first line of the .az file.  If not specified on the command line, the program will default to using the direction from the first line of the .az file.

The  structure of SPLAT! elevation pattern files is slightly different. The first line of the .el file specifies the amount of mechanical  beamtilt  applied  to  the  antenna.   A downward tilt is expressed as a positive angle, while an upward tilt is expressed as a negative angle.  This data is followed by the azimuthal direction of the tilt, 0 to 360 degrees, separated by whitespace.  These can be overridden by the command line options -dt for downtilt and -dtdir for downtilt azimuthal direction.  The options, like -rot above, will override any values embedded in the .el file if specified on the command line.

The remainder of the file consists of elevation angles and their radiation pattern (0.000 to 1.000) values separated by whitespace.  Elevation angles must  be  specified  over  a -10.0  to  +90.0  degree  range.  In  this  example,  the  antenna  is  tilted down 2 degrees towards an azimuth of 045.0 degrees.
```
2.0    045.0
-10.0   0.172
-9.5    0.109
-9.0    0.115
-8.5    0.155
-8.0    0.157
```
