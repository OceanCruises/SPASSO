![cropped-campagnebioswot2018-scaled-8](https://user-images.githubusercontent.com/48068153/229452531-b29daedb-bbdb-498f-9d7e-f339390a5c9d.jpg)

**SPASSO** (**S**oftware **P**ackage for **A**daptive **S**atellite-based **S**ampling for **O**ceanographic cruises is a Python code designed to compute daily maps of ocean state, based on available satellite data, to help guiding oceanographic cruise sampling. It can be used during oceanographic campaigns (near-real-time mode) or pre- and post-cruise to prepare sampling strategy and analyze observations. SPASSO code includes an independant Lagrangian diagnostic package (LAMTA) available <a href="https://github.com/rousseletL/lamtaLR" target="_blank">here</a>.

SPASSO website: https://spasso.mio.osupytheas.fr/

**SPASSO** automatically and independently download daily satellite-derived data such as surface height, velocity, temperature, salinity and chlorophyll-a
concentration (freely distributed on [CMEMS](https://data.marine.copernicus.eu/products) with support from CNES). Data are mapped over the studied region and diagnostics for lagrangian analysis are computed.
SPASSO outputs (figures and data files) can be directly send to user’s email and can be gathered in
a daily bulletin (.tex file).

# Installation
## Requirements
SPASSO is a package designed for **Python >= 3.9**, tested on Linux and Mac systems. It also requires a LaTeX compiler (`latexmk` or `pdflatex`).

Required python packages:\
`termcolor`, `pandas`, `netCDF4`, `motuclient`, `matplotlib`, `cmocean`, `scipy`, `xarray`, `basemap` (`pyproj`, `pyshp`, `geos`), `basemap-data-hires`, `tabulate`, `simplekml`, `requests`, `importlib-metadate`, `pylatex`, `copernicusmarine`

> [!TIP]
> python packages can be installed using the `pip3 install <package_name>` or `python -m pip install <package_name>` commands. [How to install and use pip3](https://www.activestate.com/resources/quick-reads/how-to-install-and-use-pip3/#:~:text=Pip3%20is%20the%20official%20package,in%20the%20Python%20standard%20library.).

## Download  
There are different options too clone the SPASSO repository from GitHub:
- Above the list of files, click **<> Code**, and then **Download ZIP**.
- Open a terminal, change the working directory to the location of the local cloned directory
```
git clone https://github.com/OceanCruises/SPASSO.git 
```

## Directory tree
- `Cruises/`: contains the different cruise directories.
    - `CruiseName#1/` : contains the configuration file and all SPASSO outputs organized in the corresponding directories.
- `Data/`: contains data downloaded from CMEMS.
    - `BATHY/`: must contain a NETCDF  le including global bathymetry that can
be downloaded from NOAA.
- [`Doc/`](Doc/): contains user manual and useful documents such as references.
- `src/`: contains all source code:
    - `Bulletin.py`: code preparing a tex and pdf bulletin file including all
maps computed with SPASSO.
    - `colormaps.py`: contains some useful colormaps.
    - `Diagnostics.py`: code computing Eulerian and Lagrangian diagnostics.
    - `Fields.py`: code to download, load and create netcdf for all the fields
used by SPASSO such as satellite, diagnostics and bathymetry.
    - `Functions.py`: additional useful functions (ex: climatologies).
    - `GlobalVars.py`: code loading all the global variables such as directory paths, dates, diagnostic parameters,
emailing parameters, and bulletin parameters.
    - `Library.py`: code for SPASSO run functions to print messages,
write in a log fille, send email, copy files, clean directories, execute shell requests..
.
    - `PlotField.py`: plot maps.
- `Spasso.py`: main code launching SPASSO software.

## How to run SPASSO
1. Open a terminal
2. Define the name of the Cruise:
```
cr=WMedSeaExample
```
3. Change working directory to the source code repository:
```
cd spasso2.0/src/
```
4. Launch SPASSO:
```
clear; python3 Spasso.py $cr
```

> More detailed information on how to set up a new SPASSO configuration can be found in the [user manual](Doc/usermanual.pdf). Section 5 of the manual is a Tutorial on the western Mediterranean Sea.

# SPASSO manuscript and code

The manuscript detailing the first release of SPASSO, version 2.0, is currently in revision for JTECH:

*Louise Rousselet, Francesco d’Ovidio, Lloyd Izard, Alice Della Penna, Anne Petrenko, et al.. A
Software Package for an Adaptive Satellite-based Sampling for Oceanographic cruises (SPASSOv2.0):
tracking fine scale features for physical and biogeochemical studies. 2024. hal-04705438* [preprint](https://hal.science/hal-04705438v1/file/software.pdf)

# SPASSO fundings and supports
SPASSO development has been supported by the following organisations:
- LOCEAN (Laboratoire d'Océanographie et du Climat: Expérimentations et Approches Numériques), Sorbonne University
- MIO (Mediterranean Institute of Oceanography), Aix-Marseille University
- CNES (Centre National d'Études Spatiales)
