#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Defined functions that can be useful in SPASSO configuration:
    - climatology: compute a climatology from satellite date.
    - SWOT_passing_time: find passing time of SWOT during the CalVal orbit 
    
Created on Thu Jun 16 14:25:16 2022

@author: lrousselet
"""
import GlobalVars, Library, Fields
import glob, os
import numpy as np
import pandas as pd
import xarray as xr
import datetime
import shutil

def climatology(prod,date):
    Library.printMessage("Computing climatology for "+prod)
    
    dir_wrk         = GlobalVars.Dir['dir_wrk']
    fname           = glob.glob(dir_wrk+'/*'+prod+'*.nc')
    var = 0
    u   = 0
    v   = 0
    it = 0
    for file in fname:
        it += 1
        #load file
        field = eval('Fields.'+prod+'(file).loadnc()')
        tmpvar = np.ma.filled(field['var'],np.nan)
        var += tmpvar
        if 'u' and 'v' in field:
            u += field['u']
            v += field['v']

    var = var/it
    lon = field['lon']
    lat = field['lat']
    if 'u' and 'v' in field:
        u = u/it
        v = v/it
        
    title = prod + 'Climatology between '+date[0]+' and '+date[1]
    ncfile = dir_wrk + date[0] + '-' + date[-1] + '_' + prod + '_clim.nc'
    eval('Fields.'+prod+'(ncfile).createnc(lon,lat,var,title,vu=u,vv=v)')
    
    # remove daily nc file
    [os.remove(file) for file in fname]
    
    # old #copy nc file in Processed
    # req = "cp "+ncfile+" "+GlobalVars.Dir['dir_proc']
    # Library.execute_req(req)
    
    ncfile_dest = os.path.join(GlobalVars.Dir['dir_proc'], os.path.basename(ncfile))
    shutil.copy(ncfile, ncfile_dest)  # Works on both Windows & Linux
    print(f"Copied {ncfile} to {ncfile_dest}")
    
    return

def SWOT_passing_time(crossover):
    def find_cotime(point_A, latitude, longitude, time):
        # Convert the coordinates to radians
        latitude = np.radians(latitude)
        longitude = np.radians(longitude)
        point_A = (np.radians(point_A[0]), np.radians(point_A[1]))

        # Define the radius of the Earth in meters
        R = 6371e3

        # Compute the haversine distance between point A and each point in the array
        a = np.sin((latitude - point_A[0]) / 2)**2 + np.cos(point_A[0]) * np.cos(latitude) * np.sin((longitude - point_A[1]) / 2)**2
        c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
        distances = R * c

        # Find the index of the closest point
        #print(distances.size,time.size)
        min_index = distances.argmin()
        if (distances.min()/1e3)>90:
            return np.nan*np.ones(176)
        else:
            # Get the closest point
            #closest_point = (np.degrees(latitude[min_index]), np.degrees(longitude[min_index]))
            time0=time.values[min_index]
            dt=86400-562.464 #back-shift 562 seconds every day
            freq_base = int(dt)
            freq_remain = (dt - freq_base)
            a = pd.date_range(time0,periods=176,freq=f"{freq_base}s{freq_remain * 1000:.0f}L")
            datef = a.strftime('%Y-%m-%d %H:%M:%S')
            return datef
    
    #load calval orbit
    coord=xr.open_dataset(GlobalVars.Dir['dir_data']+'SWOT/calval_orbit.nc')
    lat=coord['latitude'].values.flatten()
    lon=coord['longitude'].values.flatten()
    lat0=np.r_[0,lat]
    msk=np.diff(lat0)>0
    
    #initiate panda frame
    df = pd.DataFrame({})
    
    #loop over points
    for i in range(0,len(crossover['name'])):
        coordinate=crossover['coordinate'][i]
        region_name=crossover['name'][i]

        #Ascending
        asc=find_cotime(coordinate,lat[msk],lon[msk],coord['time'][msk])
        
        #Decending
        dsc=find_cotime(coordinate,lat[~msk],lon[~msk],coord['time'][~msk])
      
        #Print only time in the near future
        nbday = 5
        now = (datetime.datetime.strptime(GlobalVars.all_dates['today'],'%Y%m%d')).strftime('%Y-%m-%d')
        delay = (datetime.datetime.strptime(GlobalVars.all_dates['today'],'%Y%m%d')+datetime.timedelta(days=nbday)).strftime('%Y-%m-%d')
        asct,dsct = asc,dsc
        if not isinstance(asc,np.ndarray):
            mask = (asc > now) & (asc <= delay)
            asct = asct[mask]
            if len(asct)>nbday: asct = asct[0:nbday]
            df[region_name+' Asc'] = asct
        if not isinstance(dsc,np.ndarray):
            mask = (dsc > now) & (dsc <= delay)
            dsct = dsct[mask]
            if len(dsct)>nbday: dsct = dsct[0:nbday]
            df[region_name+' Desc'] = dsct

        #pd.DataFrame({'Ascending':asc,'Descending':dsc}).to_csv('%s_crossover_time.csv'%GlobalVars.Dir['dir_proc']+region_name)
    
    return df
