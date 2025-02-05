"""
Main code for SPASSO run. 
Launch all SPASSO commands.
"""

import subprocess
import sys

def check_and_install_requirements(requirements_file="../requirements.txt"):
    """Check if required packages are installed. If not, install them automatically."""
    try:
        with open(requirements_file, "r") as f:
            packages = [line.strip() for line in f.readlines() if line.strip() and not line.startswith("#")]

        missing_packages = []
        for package in packages:
            package_name = package.split("==")[0]  # Extract package name
            try:
                __import__(package_name)  # Try importing package
            except ImportError:
                missing_packages.append(package)

        if missing_packages:
            print(f"Installing missing packages: {', '.join(missing_packages)}")
            subprocess.check_call([sys.executable, "-m", "pip", "install", *missing_packages])
        else:
            print("All required packages are installed.")

    except FileNotFoundError:
        print(f"Error: `{requirements_file}` not found! Ensure the file exists.")


# Run package check before executing the main script
check_and_install_requirements()

import sys
import GlobalVars, Diagnostics
import Library
import Fields
import Bulletin
import PlotField
import Functions
from Diagnostics import Launch

def spasso():
    ################# Welcome to SPASSO
    Library.welcome_message()
    Library.printMainMessage('Initialization')
    cruise = Library.choose_cruise()
    
    GlobalVars.init_dataDir(cruise)
    GlobalVars.init_prodDate()
    Library.Listproducts()

############################### GET DATA AND PLOT
    Library.printMainMessage("Get data and make figures")
    Library.tic()
    #get products and options
    products = [str(x) for x in GlobalVars.config.get('products','products').split(',')]
    opt = [str(x) for x in GlobalVars.config.get('plot_options','options').split(',')]
    mode = GlobalVars.config.get('cruises','mode')
    outmode = GlobalVars.config.get('cruises','outmode')
    
    for pr in products:
        nprod = GlobalVars.config.get('products',pr+'prod')
        #download data
        Library.printMessage("Downloading "+nprod)
        eval('Fields.'+nprod+'.download()')
        if mode == 'DT':
            date = GlobalVars.all_dates['datec_'+pr.lower()]
            if outmode == 'clim':
                Functions.climatology(nprod,date)
                
        #plot field
        Library.printMessage("Ploting "+nprod)
        PlotField.PlotField.Plot(cruise,pr,opt,type=None)
    Library.toc('satellite figures')
    
################################ DIAGNOSTICS
    Library.printMainMessage("Computing diagnostics")

    ### Eulerian computation
    Library.tic()
    if GlobalVars.Eul['diag'] is not None:
        if outmode == 'clim':
            Diagnostics.Launch(cruise,'eulerian',clim='on')
        else:
            Diagnostics.Launch(cruise,'eulerian')
    Library.toc('Eulerian diagnostics')
    ### Lagrangian computation
    if GlobalVars.Lag['diag'] is not None:
        if outmode == 'clim':
            Diagnostics.Launch(cruise,'lagrangian',clim='on')   
        else:
            Diagnostics.Launch(cruise,'lagrangian')
    
    ############################### COPYING FIGURE AND SAVED
    
    Library.printMainMessage("COPYING AND ZIPPING DATA")
    Library.cleantmp()
    Library.copyfiles()
    
    ############################### CREATE BULLETIN
    if GlobalVars.Bull['authors']!=None:
        Library.printMainMessage("CREATING BULLETIN")
        Bulletin.create(cruise)
    
    # ############################### SENDING EMAIL
    if GlobalVars.Email['sender']!=None:
        Library.printMainMessage("SENDING EMAIL")
        Library.send_email(cruise)

    return 0

############################### PROGRAM STARTS HERE
if __name__ == '__main__':

################# Usage

    args = sys.argv

    if (len(args) != 2):
        Library.usage()
        sys.exit(2)

    ################# Initilization from 'config.ini' file

    cruise = args[1] # get config.ini from defined cruise
    GlobalVars.configIni(cruise)
    Library.clean_wrk()
    GlobalVars.FigParam()
    GlobalVars.DiagParam()
    GlobalVars.EmailParam()
    GlobalVars.BulletinParam()
    GlobalVars.LibrariesPaths()
    
    ################# Get program starting time
    GlobalVars.init_date()

    ################# lauching the program with: spasso()
    if (spasso() == 0):
        Library.exit_program(0)
    else:
        Library.exit_program(1)

############################### end of program  ###############################
