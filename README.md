![cropped-campagnebioswot2018-scaled-8](https://user-images.githubusercontent.com/48068153/229452531-b29daedb-bbdb-498f-9d7e-f339390a5c9d.jpg)

**SPASSO** (**S**oftware **P**ackage for **A**daptive **S**atellite-based **S**ampling for **O**ceanographic cruises is a Python code designed to compute daily maps of ocean state, based on available satellite data, to help guiding oceanographic cruise sampling. It can be used during oceanographic campaigns (near-real-time mode) or pre- and post-cruise to prepare sampling strategy and analyze observations. SPASSO code includes an independant Lagrangian diagnostic package (LAMTA) available <a href="https://github.com/rousseletL/lamtaLR" target="_blank">here</a>.

SPASSO website: https://spasso.mio.osupytheas.fr/

# Installation
## Requirements
SPASSO is a package designed for **Python >= 3.9**, tested on Linux and Mac systems. It also requires a LaTeX compiler (`latexmk` or `pdflatex`).

Required python packages:\
`termcolor`, `pandas`, `netCDF4`, `motuclient`, `matplotlib`, `cmocean`, `scipy`, `xarray`, `basemap` (`pyproj`, `pyshp`, `geos`), `basemap-data-hires`, `tabulate`, `simplekml`, `requests`, `importlib-metadate`, `pylatex`, `copernicusmarine`

> [!TIP]
> python packages can be installed using the `pip3 install <package_name>` command.

## Download  
There are different options too clone the SPASSO repository from GitHub:
- Above the list of files, click **<> Code**, and then **Download ZIP**.
- Open a terminal, change the working directory to the location of the local cloned directory\
```git clone https://github.com/OceanCruises/SPASSO.git 
```

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

The manuscript detailing the first release of SPASSO, version 2.0, is currently in preparation:

*Rousselet, L., d'Ovidio, F., Nencioli, F., ...(in revision) A Software Package for an Adaptive Satellite-based Sampling for Oceanographic cruises (SPASSOv2.0): tracking fine-scale features for biophysical studies.*

# SPASSO fundings and supports
SPASSO development has been supported by the following organisations:
- LOCEAN (Laboratoire d'Océanographie et du Climat: Expérimentations et Approches Numériques), Sorbonne University
- MIO (Mediterranean Institute of Oceanography), Aix-Marseille University
- CNES (Centre National d'Études Spatiales)
