#!/usr/bin/python3

from termcolor import cprint
import os
from os.path import basename
import glob
import sys
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

import GlobalVars,Fields

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

    if error_value == 1:
        cprint("\nProgram exit with errors.\n", 'red', attrs=['bold'])
        sys.exit(1)
    elif error_value == 2:
        cprint("\tCorrect data.", 'yellow', attrs=['bold'])
        sys.exit(2)
    else:
        txt = "\n##########################################################\n\
##\t\t\t\t\t\t\t\n\
##\t Successfully end program on: \n\
##\t "+get_current_date()+"\n\
##\t Thanks for using SPASSO.\n\
##\t\t\t\t\t\t\t\n\
##########################################################\n"
        cprint(txt, 'blue', attrs=['bold'])
        Logfile(txt)
        #move log file into Logs/
        logf = "spassoLog"+GlobalVars.all_dates['today']+".txt"
        req = "mv "+GlobalVars.Dir['dir_wrk']+"spassoLog.txt "\
            + GlobalVars.Dir['dir_logs']+logf
        os.system(req)
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
    os.system("rm "+GlobalVars.Dir['dir_wrk']+"*.*")

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
    os.chdir(GlobalVars.Dir['dir_wrk'])
    execute_req("cp *.png ../Figures")
    Done('Copy in Figures/ done.')
    if glob.glob(GlobalVars.Dir['dir_wrk']+'*.kmz'):
        dirk = GlobalVars.Dir['dir_fig']+'kmz/'
        isexists = os.path.exists(dirk)
        if not isexists:
            os.makedirs(dirk)
        listf = glob.glob(GlobalVars.Dir['dir_wrk']+'Figures_oftheday_*.kmz')
        for ff in listf:
            req = "cp "+ff+" "+GlobalVars.Dir['dir_wrk']+"/"+GlobalVars.all_dates['today']+"_"+os.path.basename(ff)
            execute_req(req)
        execute_req("cp *.kmz "+dirk)
        Done('Copy in Figures/kmz/ done.')

    execute_req("cp *.nc ../Processed")
    Done('Copy in Processed/ done.')
    execute_req("tar -czf "+GlobalVars.all_dates['ref']+"_Figures.tar.gz *.png")
    Done('Figures are zipped in '+GlobalVars.all_dates['ref']+'_Figures.tar.gz .')
    return

def cleantmp():
        #delete tmp png files
        figs = glob.glob(GlobalVars.Dir['dir_wrk']+'/*tmp*png')
        for ff in figs:
            os.remove(ff)

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


def get_SEN3json(datef,id0,id1):
    """
    Get products id and full name for sentinel3 data
    product date depends on satellite passing time
    so an adjustement window over the product time is required (dday1, dday2)
    """
    
    def get_url(dday1,dday2,id0,id1,**kwargs):
        url = "https://catalogue.dataspace.copernicus.eu/odata/v1/Products?$filter=contains(Name,'S3') "\
        +"and contains(Name,'_OL_2_WFR____') and contains(Name,'179_"+id0+"_') and contains(Name,'_"+id1+"_MAR_O_NR_003.SEN3')"\
            +" and ContentDate/Start gt "+dday1+" and ContentDate/Start lt "+dday2
        
        if 'Namenb' in kwargs:
            url = "https://catalogue.dataspace.copernicus.eu/odata/v1/Products?$filter=contains(Name,'S3') "\
            +"and contains(Name,'_OL_2_WFR____') and contains(Name,'"+kwargs['Namenb']+"_"+id0+"_') and contains(Name,'_"+id1+"_MAR_O_NR_003.SEN3')"\
                +" and ContentDate/Start gt "+dday1+" and ContentDate/Start lt "+dday2
    
        r = requests.get(url)
        data = r.json()
        value = data['value']
        return value
    
    datet = datetime.datetime.strptime(datef,'%Y%m%d').date()
    refd = datetime.date(2023,3,22)
    refh = datetime.datetime.strptime('0930','%H%M')
    deltad = (datet-refd).days
    hd1 = refh - (datetime.timedelta(minutes=1)*deltad) #one minute delay every day
    hd2 = hd1 + datetime.timedelta(minutes=45)
    shd1 = datetime.datetime.strftime(hd1,'%H:%M:%S')
    shd2 = datetime.datetime.strftime(hd2,'%H:%M:%S')
    dday1 = datetime.datetime.strftime(datet,'%Y-%m-%d')+'T'+shd1+'.000Z'
    dday2 = datetime.datetime.strftime(datet,'%Y-%m-%d')+'T'+shd2+'.000Z'
    
    value = get_url(dday1,dday2,id0,id1)
    
    if value:
        idp = value[0]['Id']
        name = value[0]['Name']
    else:
        hd2 = hd2 + datetime.timedelta(minutes=45)
        shd2 = datetime.datetime.strftime(hd2,'%H:%M:%S')
        dday2 = datetime.datetime.strftime(datet,'%Y-%m-%d')+'T'+shd2+'.000Z'
        value = get_url(dday1,dday2,id0,id1)
        if value:
            idp = value[0]['Id']
            name = value[0]['Name']
        else:
            value = get_url(dday1,dday2,id0,id1,Namenb='180')
            if value:
                idp = value[0]['Id']
                name = value[0]['Name']
            else:
                idp,name = None,None
    
    return idp,name

def CreateSEN3(var):
    for nf in range(len(var['date'])):
        for id0 in var['id0']:
            for id1 in var['id1']:
                exf,ff=ExistingFile(var['direct']+"/"+var['date'][nf]+"_"+id0+"_geo_coordinates_"+id1+".nc",var['date'][nf])
                if ff:
                    file = Dataset(var['direct']+"/"+var['date'][nf]+"_"+id0+"_geo_coordinates_"+id1+".nc")
                    lon = file.variables['longitude'][:,:]
                    lat = file.variables['latitude'][:,:]
                    file = Dataset(var['direct']+"/"+var['date'][nf]+"_"+id0+"_chl_oc4me_"+id1+".nc")
                    chl = file.variables['CHL_OC4ME'][:]
                    #collate variables
                    if 'chln' not in locals():
                        lonn = np.empty(np.shape(lon))
                        latn = np.empty(np.shape(lat))
                        chln = np.empty(np.shape(chl))
                        
                    #check dimensions
                    if np.shape(chl)[0]!=np.shape(chln)[0]:
                        dif = np.abs(np.shape(chl)[0]-np.shape(chln)[0])
                        lonn = np.dstack((lonn,lon[:-dif,:]))
                        latn = np.dstack((latn,lat[:-dif,:]))
                        chln = np.dstack((chln,chl[:-dif,:]))
                    else:
                        lonn = np.dstack((lonn,lon))
                        latn = np.dstack((latn,lat))
                        chln = np.dstack((chln,chl))
                else:
                    cprint('\tCould not find '+var['date'][nf]+' file to load !','red',attrs=['bold'])
    #createnc and save
        if 'chln' in locals():
            chln = chln[:,:,1:]
            lonn = lonn[:,:,1:]
            latn = latn[:,:,1:]
            fname = var['dir_wrk']+'/'+var['date'][nf]+'_'+var['prod']+'.nc'
            Fields.Sentinel3_CHL(fname).createnc3Dll(lonn,latn,chln,'Sentinel-3 '+var['date'][nf]+' CHL')
    return 

def get_SEN3token(user,pwd):
    URL = "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token"

    payload = {
        'grant_type': 'password',
        'client_id': 'cdse-public',
        'username': user,
        'password': pwd
    }
     
    r = requests.post(URL,headers={"Content-Type":"application/x-www-form-urlencoded"},
        data=payload)
    token = r.json()['access_token']
    return token

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
