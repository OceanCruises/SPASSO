import configparser  # pour recuperer variable dans fichier config.ini
from termcolor import cprint
import sys
import datetime
import os

import Library

##########################################
def init_prodDate():
    global all_dates
    global current_date_for_bulletin    
    
    Library.printMessage('Defining dates for satellite products')
    mode = config.get('cruises','mode')
    Library.Done("SPASSO run in mode: "+mode)
    
    #date of the day
    current_dt = datetime.datetime.strptime(starting_date, "%Y/%m/%d %H:%M")
    txt = "\tProgram started on: " + str(current_dt)
    cprint(txt, 'yellow', attrs=['bold'])
    Library.Logfile(txt)
    
    #init
    prod = [str(x) for x in config.get('products','products').split(',')]
    all_dates = {}
    
    if mode == 'NRT':
        
        day0_ago = current_dt
        day1_ago = datetime.datetime.today() - datetime.timedelta(days=1) # current_dt not working ?
        day2_ago = datetime.datetime.today() - datetime.timedelta(days=2)
        day3_ago = datetime.datetime.today() - datetime.timedelta(days=3)
        day4_ago = datetime.datetime.today() - datetime.timedelta(days=4)
        day5_ago = datetime.datetime.today() - datetime.timedelta(days=5)
        day6_ago = datetime.datetime.today() - datetime.timedelta(days=6)
        
        for pr in prod:
            ago = config.get('products',pr+'_date')
            all_dates['date_' + pr.lower()]  = [eval('day'+ago+'.strftime("%Y%m%d")')]
            all_dates['year_' + pr.lower()]  = [eval('day'+ago+'.strftime("%Y")')]
            all_dates['month_' + pr.lower()] = [eval('day'+ago+'.strftime("%m")')]
            all_dates['day_' + pr.lower()] = [eval('day'+ago+'.strftime("%d")')]
            all_dates['datef_' + pr.lower()]  = [eval('day'+ago+'.strftime("%Y-%m-%d")')]
        all_dates['today'] = current_dt.strftime('%Y%m%d')
        all_dates['ref'] = current_dt.strftime('%Y%m%d')
        all_dates['d0'] = str(config.get('cruises','d0'))

    elif mode == 'DT':
        dtmode = config.get('cruises','dtmode')
        refdate = [str(x) for x in config.get('cruises','refdate').split(',')]
        refdate = [datetime.datetime.strptime(x,"%Y-%m-%d") for x in refdate]
        if dtmode == 'range':
            refdate = [refdate[0]+datetime.timedelta(days=x) \
                       for x in range(((refdate[1]+datetime.timedelta(days=1))-refdate[0]).days)]
        
        for pr in prod:
            if 'CLS' in pr:
                all_dates['date_' + pr.lower()] = [x.strftime('%Y-%m-%d') for x in refdate]
            else:
                all_dates['date_' + pr.lower()] = [x.strftime('%Y%m%d') for x in refdate]
            all_dates['year_' + pr.lower()] = [x.strftime('%Y') for x in refdate]
            all_dates['month_' + pr.lower()] = [x.strftime('%m') for x in refdate]
            all_dates['day_' + pr.lower()] = [x.strftime('%d') for x in refdate]
            all_dates['datef_' + pr.lower()] = [x.strftime('%Y-%m-%d') for x in refdate]
        all_dates['today'] = current_dt.strftime('%Y%m%d')
        all_dates['ref'] = refdate[0].strftime('%Y%m%d')
        all_dates['d0'] = str(config.get('cruises','d0'))

    Library.Logfile(all_dates)
    Library.Done()
    
##########################################
def init_dataDir(cruise):
    Library.printMessage('Creating paths to data directories for ' + str(cruise))
    
    global repositories

    directories(cruise)
    repositories = {}
    
    for key in config['products']:
        if 'path' in key:
            dirn = 'dir_' + key[5:]
            if 'phy' in key:
                path = Dir['dir_data']+ 'ALTI/' + config['products'][key]
            elif 'sst' in key:
                path = Dir['dir_data']+ 'SST/' + config['products'][key]
            elif 'sss' in key:
                path = Dir['dir_data']+ 'SSS/' + config['products'][key]
            elif 'chl' in key:
                path = Dir['dir_data']+ 'CHL/' + config['products'][key]
            elif 'medsea' in key:
                path = Dir['dir_data']+ 'MEDSEA/' + config['products'][key]
            
            if 'cls' in key:
                name = config['products']['name_'+key[5:]]
                path += '/' + name
            repositories[dirn] = path

    #check if repositories exist else create
    for rep in repositories:
        isexists = os.path.exists(repositories[rep])
        if not isexists:
            txt = '\t\t '+ repositories[rep] +'-> '+str(isexists)
            os.makedirs(repositories[rep])
            Library.Logfile(txt)
 
    Library.Logfile(repositories)
    Library.Done()

#######
def directories(cruise):
    global Dir

    main_path = os.path.abspath(os.path.join(os.getcwd(), os.pardir))+'/'
    cruise_path = main_path + "Cruises/" + cruise + "/" 
    
    Dir = {
        'main_path' : main_path,
        'cruise_path' : cruise_path,
		'dir_data' : main_path + 'Data/',
		'dir_wrk' : cruise_path + 'Wrk/',
        'dir_bull' : cruise_path + 'Bulletin/',
		'dir_src' : main_path + 'src/',
		'dir_proc' : cruise_path + 'Processed/',
		'dir_fig' : cruise_path + 'Figures/',
		'dir_logs' : cruise_path + 'Logs/',
        'dir_bathy' : main_path + 'Data/BATHY/'
	}
    
    #check if repositories exist else create
    for rep in Dir:
        isexists = os.path.exists(Dir[rep])
        if not isexists:
            os.makedirs(Dir[rep])

def DiagParam():
    global Lag,Eul

    item = config.items('Lagrangian')
    diag = [str(x) for x in item[0][1].split(',')]
    if diag[0] == 'None': diag = None      
    Lag = {
        'diag' : diag,
        'products' : [str(x) for x in item[1][1].split(',')],
        'mode' : str(item[2][1]),
		'dayv' : item[3][1],
		'loni' : [float(x) for x in item[4][1].split(',')],
        'lati' : [float(x) for x in item[5][1].split(',')],
		'numdays' : int(item[6][1]),
		'delta0' : float(item[7][1]),
        'PeriodicBC' : item[8][1],
        'method' : str(item[9][1]),
        'numstep' : int(item[10][1]),
        'bathyfile' : str(item[11][1]),
        'bathylvl' : float(item[12][1]),
        'sstprod' : str(item[13][1]),
        'sstadvd' : int(item[14][1]),
	}
    
    item = config.items('Eulerian')
    diag = [str(x) for x in item[0][1].split(',')]
    if diag[0] == 'None': diag = None
    Eul = {
        'diag' : diag,
        'products' : [str(x) for x in item[1][1].split(',')],
		'dayv' : item[2][1],
		'loni' : [float(x) for x in item[3][1].split(',')],
        'lati' : [float(x) for x in item[4][1].split(',')],
		'delta0' : float(item[5][1]),
		'UVunit' : item[6][1]
	}
    
def EmailParam():
    global Email

    item = config.items('email')
    if item[3][1]:
        port = int(item[3][1])
    else:
        port = 0
    send = str(item[0][1])
    if send == 'None': send = None
        
    Email = {
        'sender' : send,
        'receiver' : [str(x) for x in item[1][1].split(',')],
        'smtp' : str(item[2][1]),
		'port' : port,
		'login' : str(item[4][1]),
        'password' : str(item[5][1]),
        'attach' : [str(x) for x in item[6][1].split(',')]
	}

def BulletinParam():
    global Bull

    item = config.items('bulletin')
    auth = [str(x) for x in item[0][1].split(',')]
    if auth[0] == 'None': auth = None
    
    Bull = {
        'authors' : auth,
        'acknow' : str(item[1][1]),
	}
    if len(item)>2:
        Bull.update({'swotCOn': [str(x) for x in item[2][1].split(',')],
                     'swotCOlon' : [float(x) for x in item[3][1].split(',')],
                     'swotCOlat' : [float(x) for x in item[4][1].split(',')]})
    
def LibrariesPaths():
    global Lib

    item = config.items('library')      
    Lib = {
        'motulib' : str(item[0][1]),
        'latexcompiler' : str(item[1][1]),
	}

##
def configIni(cruise):
    global config
    
    directories(cruise)
    configFile = Dir['cruise_path'] + 'config_' + cruise + '.ini'
    
    try:
        Library.Logfile("\tReading configuration file: " + str(configFile))
        f = open(configFile, "r")
        Library.Logfile("\tOK.")
    except:
        Library.Logfile("\t\tCould not find configuration file.")
        sys.exit(1)

    config = configparser.ConfigParser()
    config.read(configFile)
    f.close()

##
def init_date():
    global starting_date
    starting_date = Library.get_current_date()
    
def LATEXini(cruise):
    global latexini
    
    directories(cruise)
    configFile = Dir['cruise_path'] + 'LATEXtools_' + cruise + '.ini'
    
    try:
        Library.printMessage("\tReading LATEXtools initialization file: " + str(configFile))
        f = open(configFile, "r")
    except:
        sys.exit(1)

    latexini = configparser.ConfigParser()
    latexini.read(configFile)
    f.close()