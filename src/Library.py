#!/usr/bin/python3
"""
List of required functions used in SPASSO run:
    - Printing information on screen:
        - Done, printMessage, printInfo, printMainMessage, welcomeMessage
    - Create and write a Log file: Logfile, ListProducts
    - Define cruise
    - exit program with error message
    - Check if file already exit: ExistingFile
    - Clean working directory: clean_wrk
    - Create a dictionnary with all parameters for a specific field: GetVars
    - Send email: send_email
    - Copying files from working directory to specific directories (Figures/, 
    Processed, Bulletin/): copy_files
    - Create kml output files: make_kml
    
"""
from termcolor import cprint
import os
from os.path import basename
import glob
import sys
sys.stdout.reconfigure(encoding='utf-8')
import datetime
import time
import smtplib
import numpy as np
from netCDF4 import Dataset
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.application import MIMEApplication
from email.mime.base import MIMEBase
from email import encoders
import requests
from simplekml import (Kml, OverlayXY, ScreenXY, Units, RotationXY,
                       AltitudeMode, Camera)
import matplotlib.pyplot as plt
import shutil
from pathlib import Path

import GlobalVars,Fields

def usage():
    cprint("\t Should be called with two arguments, such as: python3 Spasso.py WMedSeaExample", 'red', attrs=['bold'])

def Done(*args):
    if not args:
        cprint("\t Done", 'green', attrs=['bold'])
        Logfile("\n \t Done \n")
    else:
        for arg in args:
            cprint("\t"+arg, 'green', attrs=['bold'])
            Logfile("\n \t "+arg+" \n")

def Logfile(txt):
    if os.path.isfile(GlobalVars.Dir['dir_wrk']+'spassoLog.txt'):
        f = open(GlobalVars.Dir['dir_wrk']+'spassoLog.txt','a')
        if isinstance(txt, str):
            f.write(txt+'\n')
        elif isinstance(txt, dict):
            for key, value in txt.items(): 
                     f.write('\t\t%s : %s\n' % (key, value))        
        f.close()
    else:
        f = open(GlobalVars.Dir['dir_wrk']+'spassoLog.txt','w')
        f.write(txt+'\n')
        f.close()
    
def printMessage(message):
    txt = '\t## ' + message
    cprint(txt,'cyan',attrs=['bold'])
    Logfile(txt)
    
def printInfo(message):
    txt = '\t '+ message
    cprint(txt, 'yellow', attrs=['bold'])

def Listproducts():
    dic = dict(GlobalVars.config.items('products'))
    Logfile('List of products:')
    Logfile(dic)

    
def choose_cruise():
    cruises = os.listdir("../Cruises")
    found_campaigns = "\tThe following campaigns (repositories) were found by SPASSO:\n" \
        + '\t\t' + str(cruises) + '\n'
    Logfile(found_campaigns)
    user_selection_config = GlobalVars.config.get('cruises', 'cruise')
    current_campaign = "\tSpasso is actually working on the " + user_selection_config + " cruise.\n"
    cprint(current_campaign, 'yellow', attrs=['bold'])
    Logfile(current_campaign)
    return user_selection_config

## current date
def get_current_date():
	today = datetime.datetime.now()
	return today.strftime("%Y/%m/%d %H:%M")

def printMainMessage(message):
    txt = '\n##################################################\n'\
        + '##\t' + message + '\n' \
            + '##################################################\n'
    #cprint('\n\t##################################################', 'blue', attrs=['bold'])
    #cprint('\t## ' + message, 'blue', attrs=['bold'])
    #cprint('\t##################################################\n', 'blue', attrs=['bold'])
    cprint(txt,'blue', attrs=['bold'])
    Logfile(txt)

def welcome_message():
    welcome_message = "\n##########################################################\n\
##\t\t\t\t\t\t\t\n\
##\t\t\t\t\t\t\t\n\
##\t SPASSO software V2.0\t\n\
##\tContact: " + str(GlobalVars.config.get('email', 'sender_mail')) + "\t\t\n\
##\t\t\t\t\t\t\t\n\
##########################################################\n"
    
    cprint(welcome_message, 'blue', attrs=['bold'])
    Logfile(welcome_message)

def exit_program(error_value):
    """Handles SPASSO exit with appropriate messages and error codes."""
    
    if error_value == 1:
        cprint("\n Program exit with errors.\n", 'red', attrs=['bold'])
        sys.exit(1)
    elif error_value == 2:
        cprint("\t Correct data.", 'yellow', attrs=['bold'])
        sys.exit(2)
    else:
        txt = f"""
##########################################################
##                                                    
##    Successfully ended program on: 
##    {get_current_date()}
##    Thanks for using SPASSO. 
##                                                    
##########################################################
        """
        cprint(txt, 'blue', attrs=['bold'])
        Logfile(txt)

        # Move log file into Logs/
        log_filename = f"spassoLog{GlobalVars.all_dates['today']}.txt"
        log_source = Path(GlobalVars.Dir['dir_wrk']) / "spassoLog.txt"
        log_dest = Path(GlobalVars.Dir['dir_logs']) / log_filename

        # Ensure Logs directory exists
        log_dest.parent.mkdir(parents=True, exist_ok=True)

        try:
            shutil.move(log_source, log_dest)
            print(f"Log file moved to: {log_dest}")
        except Exception as e:
            print(f"Error moving log file: {e}")

        sys.exit(0)

def ExistingFile(fname,date):
    fn = str(fname)
    ff = glob.glob(fn)
    exf = False if not ff else os.path.exists(ff[0])
    if exf == True:
        Done(date+' file downloaded')
        ff = fn
    elif exf == False:
        if len(ff)>1:
            Done('Several files available: '+date+' '+GlobalVars.all_dates['d0']+' used.')
            files_ascending = sorted(ff, key=lambda t: os.stat(t).st_mtime)
            if GlobalVars.all_dates['d0']=='d0':
                ff = files_ascending[0]
            elif GlobalVars.all_dates['d0']=='d1':
                ff = files_ascending[1]
        else:
            cprint('\tCould not find '+date+' file to load !','red',attrs=['bold'])
	
    return exf,ff

def clean_wrk():
    wrk_dir = GlobalVars.Dir['dir_wrk']

    # Check if directory exists
    if os.path.exists(wrk_dir):
        print(f"🧹 Cleaning work directory: {wrk_dir}")
        # Remove all files in directory
        for file in os.listdir(wrk_dir):
            file_path = os.path.join(wrk_dir, file)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.remove(file_path)  # Remove file
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)  # Remove folder
            except Exception as e:
                print(f" Error deleting {file_path}: {e}")
    else:
        print(f"Work directory does not exist, nothing to clean: {wrk_dir}")

def execute_req(req):
    Logfile("\tTrying to execute:\n")
    Logfile("\t" + str(req) + "\n")   
    try:
        os.system(req)
    except:
        Logfile("\tERROR in:\n" + req + "\nExiting program.")
        cprint("\tERROR in:\n" + req + "\nExiting program.", 'red',attrs=['bold'])
        exit()

def UpdateVars(var,**kwargs):
    for key, value in kwargs.items():
        up = {key: value}
        var.update(up)
        
def GetVars(data):
    if 'cls' in data:
        user   = GlobalVars.config.get('userpwd', 'userCLS')
        pwd    = GlobalVars.config.get('userpwd', 'pwdCLS')
        name   = GlobalVars.config.get('products', data+'_name')
        arc    = GlobalVars.config.get('products', data+'_arc')
    elif 'h8' in data:
        user   = GlobalVars.config.get('userpwd', 'userH8')
        pwd    = GlobalVars.config.get('userpwd', 'pwdH8')
    elif 'sen3' in data:
        user   = GlobalVars.config.get('userpwd', 'userSEN3')
        pwd    = GlobalVars.config.get('userpwd', 'pwdSEN3')
        id0    = [str(x) for x in GlobalVars.config.get('products', data+'_id0').split(',')]
        id1    = [str(x) for x in GlobalVars.config.get('products', data+'_id1').split(',')]
    else:
        user   = GlobalVars.config.get('userpwd', 'userCMEMS')
        pwd    = GlobalVars.config.get('userpwd', 'pwdCMEMS')
        

    direct     = GlobalVars.repositories['dir_'+data]
    idd        = GlobalVars.config.get('products', data+'_id')    
    prod       = GlobalVars.config.get('products', data+'prod')
    year       = GlobalVars.all_dates['year_'+data]
    month      = GlobalVars.all_dates['month_'+data]
    day        = GlobalVars.all_dates['day_'+data]
    date       = GlobalVars.all_dates['date_'+data]
    datef      = GlobalVars.all_dates['datef_'+data]
    datec      = GlobalVars.all_dates['datec_'+data]
    dir_wrk    = GlobalVars.Dir['dir_wrk']
    logf       = GlobalVars.Dir['dir_wrk'] + 'spassoLog.txt'
    mode       = GlobalVars.config.get('cruises', 'mode')
    Lon        = [float(x) for x in GlobalVars.config.get('cruise_param','Lon').split(',')]
    Lat        = [float(x) for x in GlobalVars.config.get('cruise_param','Lat').split(',')]
    
    var = {"direct":direct,
           "id":idd,
           "user":user,
           "pwd":pwd,
           "prod":prod,
           "year":year,
           "month":month,
           "day":day,
           "date":date,
           "datef":datef,
           "datec":datec,
           "dir_wrk":dir_wrk,
           "logf":logf,
           "mode":mode,
           "Lon":Lon,
           "Lat":Lat,
        }

    if 'name' in locals(): var['name'] = name
    if 'arc' in locals(): var['arc'] = arc
    if 'id0' in locals(): var['id0'] = id0
    if 'id1' in locals(): var['id1'] = id1
    
    return var

def send_email(cruise):
    msg = MIMEMultipart()
    msg["Subject"] = " "+cruise+" Spasso files"
    msg["From"] = GlobalVars.Email['sender']
    msg["To"] = ",".join(GlobalVars.Email['receiver'])
    # write the text/plain part
    text = """\
    Dear Spasso user,
    
    Please find attached the daily Figures computed by SPASSO for the """+cruise+""" cruise.
    
    *** This email was automatically generated by Python3 from @satellite machine ***
    
    """

    # convert to MIMEText objects and add them to the MIMEMultipart message
    msg.attach(MIMEText(text, "plain"))
    
    #attach figures and tex
    files = []
    if 'tar' in GlobalVars.Email['attach']:
        files += glob.glob(GlobalVars.Dir['dir_wrk']+'*.tar.gz')
    if 'tex' in GlobalVars.Email['attach']:
        files += glob.glob(GlobalVars.Dir['dir_wrk']+'*bulletin*tex')
    for f in files:
        with open(f, "rb") as fil:
            part = MIMEApplication(fil.read(),Name=basename(f))
        msg.attach(part)

    #attach pdf
    if 'pdf' in GlobalVars.Email['attach']:
        pdfpath = glob.glob(GlobalVars.Dir['dir_wrk']+'*bulletin*.pdf')
        pdfname = os.path.basename(pdfpath[0])
        binpdf = open(pdfpath[0],'rb')
        payload = MIMEBase('application', 'octate-stream', Name=pdfname)
        payload.set_payload((binpdf).read())
        encoders.encode_base64(payload)
        msg.attach(payload)

    # send your email
    with smtplib.SMTP(GlobalVars.Email['smtp'],GlobalVars.Email['port']) as server:
        server.connect(GlobalVars.Email['smtp'],GlobalVars.Email['port'])
        server.ehlo()
        server.starttls()
        server.ehlo()
        if GlobalVars.Email['login']!='':
            server.login(GlobalVars.Email['login'],GlobalVars.Email['password'])
        server.sendmail(GlobalVars.Email['sender'],GlobalVars.Email['receiver'],msg.as_string())
    txt = 'Email was sent to:'+','.join(GlobalVars.Email['receiver'])
    Done(txt)
    return
 
def copyfiles():
    """Copies relevant files to Figures and Processed directories and zips figures."""
    os.chdir(GlobalVars.Dir['dir_wrk'])

    # Copy PNG files to Figures directory
    try:
        figures_dir = Path("../Figures")
        figures_dir.mkdir(exist_ok=True)
        
        for file in Path(GlobalVars.Dir['dir_wrk']).glob("*.png"):
            print(file)
            shutil.copy(file, figures_dir)
        
        Done("Copy to Figures/ done.")
    except Exception as e:
        print(f"Error copying PNG files: {e}")

    # Handle KMZ files
    try:
        kmz_files = list(Path(GlobalVars.Dir['dir_wrk']).glob("*.kmz"))
        if kmz_files:
            kmz_dir = Path(GlobalVars.Dir['dir_fig']) / "kmz"
            kmz_dir.mkdir(parents=True, exist_ok=True)

            # Copy daily KMZ files
            for file in Path(GlobalVars.Dir['dir_wrk']).glob("Figures_oftheday_*.kmz"):
                dest_file = Path(GlobalVars.Dir['dir_wrk']) / f"{GlobalVars.all_dates['today']}_{file.name}"
                shutil.copy(file, dest_file)

            # Move all KMZ files to kmz directory
            for file in kmz_files:
                shutil.copy(file, kmz_dir)

            Done("Copy to Figures/kmz/ done.")
    except Exception as e:
        print(f"Error handling KMZ files: {e}")

    # Copy NetCDF files to Processed directory
    try:
        processed_dir = Path("../Processed")
        processed_dir.mkdir(exist_ok=True)
        
        for file in Path(GlobalVars.Dir['dir_wrk']).glob("*.nc"):
            shutil.copy(file, processed_dir)

        Done("Copy to Processed/ done.")
    except Exception as e:
        print(f"Error copying NetCDF files: {e}")

    # Zip all PNG files
    try:
        zip_filename = f"{GlobalVars.all_dates['ref']}_Figures.tar.gz"
        shutil.make_archive(zip_filename.replace(".tar.gz", ""), "gztar", root_dir=GlobalVars.Dir['dir_wrk'], base_dir=".", verbose=True)
        Done(f"Figures are zipped in {zip_filename}.")
    except Exception as e:
        print(f"Error zipping Figures: {e}")

def cleantmp():
    """Deletes temporary PNG files in the working directory."""
    
    try:
        figs = list(Path(GlobalVars.Dir['dir_wrk']).glob("*tmp*png"))
        print(figs)
        for ff in figs:
            print(ff)
            ff.unlink()
        Done(f"Deleted {len(figs)} temporary PNG files.")
    except Exception as e:
        print(f"Error deleting temporary PNG files: {e}")

#Homemade version of matlab tic and toc functions
def tic():
    global startTime_for_tictoc
    startTime_for_tictoc = time.time()

def toc(comment):
    if 'startTime_for_tictoc' in globals():
        txt = "Elapsed time to compute "+comment+" is " + str(time.time() - startTime_for_tictoc) + " seconds.\n\n"
    else:
        txt = "Toc: start time not set"
    Logfile(txt)

#### Functions to save kml outputs
def make_kml(llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat,
             figs, colorbar=None, **kw):
    """TODO: LatLon bbox, list of figs, optional colorbar figure,
    and several simplekml kw..."""

    kml = Kml()
    altitude = kw.pop('altitude', 2e6)
    roll = kw.pop('roll', 0)
    tilt = kw.pop('tilt', 0)
    altitudemode = kw.pop('altitudemode', AltitudeMode.relativetoground)
    camera = Camera(latitude=np.mean([urcrnrlat, llcrnrlat]),
                    longitude=np.mean([urcrnrlon, llcrnrlon]),
                    altitude=altitude, roll=roll, tilt=tilt,
                    altitudemode=altitudemode)

    kml.document.camera = camera
    draworder = 0
    names = kw['name']
    ii = 0
    for fig in figs:  # NOTE: Overlays are limited to the same bbox.
        draworder += 1
        ground = kml.newgroundoverlay(name='GroundOverlay')
        ground.draworder = draworder
        ground.visibility = kw.pop('visibility', 1)
        ground.name = names[ii]
        ground.color = kw.pop('color', '9effffff')
        ground.atomauthor = kw.pop('author', 'ocefpaf')
        ground.latlonbox.rotation = kw.pop('rotation', 0)
        ground.description = kw.pop('description', 'Matplotlib figure')
        ground.gxaltitudemode = kw.pop('gxaltitudemode',
                                       'clampToSeaFloor')
        ground.icon.href = fig
        ground.latlonbox.east = llcrnrlon
        ground.latlonbox.south = llcrnrlat
        ground.latlonbox.north = urcrnrlat
        ground.latlonbox.west = urcrnrlon
        ii = ii+1

    if colorbar:  # Options for colorbar are hard-coded (to avoid a big mess).
        screen = kml.newscreenoverlay(name='ScreenOverlay')
        screen.icon.href = colorbar
        screen.overlayxy = OverlayXY(x=0, y=0,
                                     xunits=Units.fraction,
                                     yunits=Units.fraction)
        screen.screenxy = ScreenXY(x=0.015, y=0.075,
                                   xunits=Units.fraction,
                                   yunits=Units.fraction)
        screen.rotationXY = RotationXY(x=0.5, y=0.5,
                                       xunits=Units.fraction,
                                       yunits=Units.fraction)
        screen.size.x = 0
        screen.size.y = 0
        screen.size.xunits = Units.fraction
        screen.size.yunits = Units.fraction
        screen.visibility = 1

    kmzfile = kw.pop('kmzfile', 'overlay.kmz')
    kml.savekmz(kmzfile)

def gearth_fig(llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat, pixels=1024):
    """Return a Matplotlib `fig` and `ax` handles for a Google-Earth Image."""
    aspect = np.cos(np.mean([llcrnrlat, urcrnrlat]) * np.pi/180.0)
    xsize = np.ptp([urcrnrlon, llcrnrlon]) * aspect
    ysize = np.ptp([urcrnrlat, llcrnrlat])
    aspect = ysize / xsize

    if aspect > 1.0:
        figsize = (10.0 / aspect, 10.0)
    else:
        figsize = (10.0, 10.0 * aspect)

    if False:
        plt.ioff()  # Make `True` to prevent the KML components from poping-up.
    figtmp = plt.figure(figsize=figsize,
                     frameon=False,
                     dpi=pixels//10)
    # KML friendly image.  If using basemap try: `fix_aspect=False`.
    axtmp = figtmp.add_axes([0, 0, 1, 1])
    axtmp.set_xlim(llcrnrlon, urcrnrlon)
    axtmp.set_ylim(llcrnrlat, urcrnrlat)
    return figtmp, axtmp
