#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Define satellite and diagnostic fields.
    - Download satellite data from remote server.
    - Functions to load and create netcdf files for each fields.
    - Parameter for each fields are defined here (variable names and units, 
                                                  colorbars).

Created on Thu Jun  2 13:54:17 2022

@author: lrousselet
"""
import GlobalVars, Library
from netCDF4 import Dataset
import cmocean as cm_oc
import cblind as cb
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap
import numpy as np
import datetime
import pandas as pd
import glob
import copernicusmarine

class Load():
    """
    Class to load satellite data from netcdf file created by SPASSO.
    """
    def loadnc(self):
        file = Dataset(self.fname)
        self.lon = file.variables[self.lon_name][:]
        self.lat = file.variables[self.lat_name][:]
        if self.var_name == ' ':
            self.var = np.full((len(self.lat),len(self.lon)),np.nan)
        else:
            if isinstance(self.var_name,str):
                if self.lon.ndim==3:
                    self.var = file.variables[self.var_name][:,:,:]
                else:
                    self.var = file.variables[self.var_name][0,:,:]   
            elif isinstance(self.var_name,tuple):
                if len(self.var_name)==2:
                    self.var = file.variables[self.var_name[0]][0,:,:],file.variables[self.var_name[1]][0,:,:]
                elif len(self.var_name)==3:
                    self.var = file.variables[self.var_name[0]][0,:,:],file.variables[self.var_name[1]][0,:,:],file.variables[self.var_name[2]][0,:,:]
                elif len(self.var_name)==4:
                    self.var = file.variables[self.var_name[0]][0,:,:],file.variables[self.var_name[1]][0,:,:],file.variables[self.var_name[2]][0,:,:],file.variables[self.var_name[3]][0,:,:]
            #if sst data
            if hasattr(self,'K'):
                self.var = self.var - 273.15
            #if chl data in log10 base
            if hasattr(self,'log10'):
                self.var = 10**self.var
        #if u and v data
        if hasattr(self,'u_name'):
            self.u = file.variables[self.u_name][0,:,:]
            self.v = file.variables[self.v_name][0,:,:]
            if hasattr(self,'dimz'):
                self.u,self.v = self.u[int(self.dimz),:,:],self.v[int(self.dimz),:,:]
            if hasattr(self,'tit'):
                field = {'lon':self.lon,'lat':self.lat,'var':self.var,'u':self.u,'v':self.v,'cm':self.cmap,'title':self.tit}
            else:
                field = {'lon':self.lon,'lat':self.lat,'var':self.var,'u':self.u,'v':self.v,'cm':self.cmap}
        else:
            if hasattr(self,'tit'):
                field = {'lon':self.lon,'lat':self.lat,'var':self.var,'cm':self.cmap,'title':self.tit}
            else:
                field = {'lon':self.lon,'lat':self.lat,'var':self.var,'cm':self.cmap}
        if hasattr(self,'colnorm'):
            field['colnorm'] = self.colnorm
           
        file.close()
        return field
    
    def LoadLag(self,numdays,product,**kwargs):
        """
        Function to load velocity field from netcdf file for Lagrangian computation.
        
        Output: dictioonnary with all velocity fields loaded to compute Lagrangian
        trajectories. Numdays is the number of daily velocity fields to load.
        """
        RT=6371e3
        self.u_all = []
        self.v_all = []
        self.dates = []
        self.day2 = self.date
        self.day1 = datetime.date.fromordinal(datetime.date.toordinal(self.day2)-numdays)
        all_days = pd.date_range(self.day1,self.day2)
        
        for i in range(len(all_days)):
            dayv = datetime.datetime.strftime(all_days[i],"%Y%m%d")
            paths = self.data_dir+"/*"+dayv+"_*.nc"
            exf,ff = Library.ExistingFile(paths,dayv)
            if exf == True:
                filei = glob.glob(paths)[0]
            elif exf == False and len(ff)<1:
                #download data
                eval(product+'.download(date=dayv)')
                filei = glob.glob(self.data_dir+"/*"+dayv+"_*.nc")[0]
            file = Dataset(filei)
            self.dates.append(datetime.date.toordinal(all_days[i]))
            self.lon =file.variables[self.lon_name][:]
            self.lat =file.variables[self.lat_name][:]
            self.u =file.variables[self.u_name][0,:,:].T #time,lon,lat
            self.v =file.variables[self.v_name][0,:,:].T
            [X,Y]=np.meshgrid(self.lon,self.lat)
            self.u_dd=(self.u/(RT*np.cos(Y.T/180*np.pi))*180/np.pi)*24*60*60 ## convert from m/s to degree/day
            self.v_dd=(self.v*180/np.pi/RT)*24*60*60 ## convert from m/s to degree/day
            self.u_dd = self.u_dd.filled(np.nan)
            self.v_dd = self.v_dd.filled(np.nan)
            self.u_all.append(self.u_dd)
            self.v_all.append(self.v_dd)
  
        self.u_all=np.array(self.u_all)
        self.v_all=np.array(self.v_all)
        self.dates=np.array(self.dates)

        field = {'lon':self.lon,'lat':self.lat,'u':self.u_all,'v':self.v_all,'dates':self.dates}
        return field 

class Create():
    """
    Class to create a netcdf file with satellite data cropped on the 
    domain defined in the config.ini file.
    """
    def createnc(self,lon,lat,vvar,title,vu=None,vv=None,vvar2=None,vvar3=None):
        file = Dataset(self.fname,mode='w',format='NETCDF4_CLASSIC') 
        file.createDimension('lon', len(lon))
        file.createDimension('lat', len(lat))
        file.createDimension(self.d3_name, 1)
        file.title = title
        latitude = file.createVariable(self.lat_name, np.float32, ('lat',))
        latitude.units = 'degrees_north'
        longitude = file.createVariable(self.lon_name, np.float32, ('lon',))
        longitude.units = 'degrees_east'
        longitude[:] = lon
        latitude[:] = lat
        if isinstance(self.var_name,str):
            if self.var_name == ' ':
                var = file.createVariable('None', np.float32,(self.d3_name,'lat','lon'))
            else:
                var = file.createVariable(self.var_name, np.float32,(self.d3_name,'lat','lon'))
                var.units = self.var_units
            if hasattr(self,'K'):
                var[:,:,:] = vvar + 273.15
                #if chl data in log10 base
            elif hasattr(self,'log10'):
                var[:,:,:] = np.log10(vvar)
            else:
                var[:,:,:] = vvar
        elif isinstance(self.var_name,tuple):
            var = file.createVariable(self.var_name[0], np.float32,(self.d3_name,'lat','lon'))
            var.units = self.var_units[0]
            var[:,:,:] = vvar
            if vvar2 is not None:
                var2 = file.createVariable(self.var_name[1], np.float32,(self.d3_name,'lat','lon'))
                var2.units = self.var_units[1]
                var2[:,:,:] = vvar2
            if vvar3 is not None:
                var3 = file.createVariable(self.var_name[2], np.float32,(self.d3_name,'lat','lon'))
                var3.units = self.var_units[2]
                var3[:,:,:] = vvar3
        if hasattr(self,'u_name'):
            u = file.createVariable(self.u_name, np.float32,(self.d3_name,'lat','lon'))
            u.units = self.u_units
            v = file.createVariable(self.v_name, np.float32,(self.d3_name,'lat','lon'))
            v.units = self.v_units
            u[:,:,:] = vu
            v[:,:,:] = vv
            
        file.close();
        txt = self.fname + ' created.'
        Library.Logfile(txt)

    def createnc3Dll(self,lon,lat,vvar,title):
        file = Dataset(self.fname,mode='w',format='NETCDF4_CLASSIC') 
        file.createDimension('lon', np.shape(lon)[0])
        file.createDimension('lat', np.shape(lat)[1])
        file.createDimension(self.d3_name, np.shape(lon)[2])
        file.title = title
        latitude = file.createVariable(self.lat_name, np.float32, ('lon','lat',self.d3_name))
        latitude.units = 'degrees_north'
        longitude = file.createVariable(self.lon_name, np.float32, ('lon','lat',self.d3_name))
        longitude.units = 'degrees_east'
        longitude[:,:,:] = lon
        latitude[:,:,:] = lat
        if isinstance(self.var_name,str):
            var = file.createVariable(self.var_name, np.float32,('lon','lat',self.d3_name))
            var.units = self.var_units
            var[:,:,:] = vvar
            
        file.close();
        txt = self.fname + ' created.'
        Library.Logfile(txt)

##################################################################
## COPERNICUS DATA (CMEMS)
##################################################################
class Copernicus_PHY(Load,Create): 
    def __init__(self,fname,**kwargs):
        self.fname = fname
        data       = GlobalVars.config.get('products', 'phy_data')
        var        = Library.GetVars(data)
        self.data_dir   = var['direct']
        self.lon_name='longitude'
        self.lat_name='latitude'
        self.d3_name = 'time'
        self.var_name = 'adt'
        self.var_units = 'm'
        self.u_name='ugos'
        self.v_name='vgos'
        self.u_units = 'm/s'
        self.v_units = 'm/s'
        self.cmap = cb.cbmap("cb.pregunta_r")
        if 'dayv' in kwargs: 
            self.date = datetime.datetime.strptime(kwargs['dayv'],"%Y-%m-%d")
        else: 
            self.date = datetime.datetime.strptime(var['date'][0],"%Y%m%d")
    
    def download(**kwargs):
        # setting global variables to local for a shorter req
        data       = GlobalVars.config.get('products', 'phy_data')
        var        = Library.GetVars(data)
        
        if ('date' in kwargs):
            ddate = [kwargs['date']]
            year = [datetime.datetime.strftime(datetime.datetime.strptime(ddate[0],'%Y%m%d'),'%Y')]
            mo = [datetime.datetime.strftime(datetime.datetime.strptime(ddate[0],'%Y%m%d'),'%m')]
        else:
            ddate = var['date']
            year = var['year']
            mo = var['month']
        
        for nf in range(len(ddate)):
            # dowloading data in /DATA
            date_range = "*/" + year[nf] + "/" + mo[nf] + "/*_allsat_" + data + \
                "_*_" + ddate[nf] + "_*.nc"
            query_metadata = copernicusmarine.get(dataset_id=var['id'],username=var['user'],
                                             password=var['pwd'],output_directory=var['direct'],
                                             filter=date_range,no_directories='Y',overwrite='Y')
            for file in query_metadata.files:
                get_files = file.file_path
            exf,ff = Library.ExistingFile(get_files,ddate[nf])

            if len(kwargs)==0 or 'cp' in kwargs:
                # copying data in /Wrk
                req_cp = "cp " + ff +' '+ var['dir_wrk'] + ddate[nf] + "_" + var['prod'] + '.nc'
                Library.execute_req(req_cp)
            
        return

class Copernicus_PHYTOT(Load,Create):
    def __init__(self,fname,**kwargs):
        self.fname = fname
        data = GlobalVars.config.get('products', 'phytot_data')
        var        = Library.GetVars(data)
        self.data_dir   = var['direct']
        self.lon_name='longitude'
        self.lat_name='latitude'
        self.d3_name = 'depth'
        self.var_name = ' '
        self.u_name = 'uo'
        self.v_name = 'vo'
        self.u_units = 'm/s'
        self.v_units = 'm/s'
        self.dimz = '0'
        self.cmap = cb.cbmap("cb.pregunta_r")
        if 'dayv' in kwargs: 
            self.date = datetime.datetime.strptime(kwargs['dayv'],"%Y-%m-%d")
        else: 
            self.date = datetime.datetime.strptime(var['date'][0],"%Y%m%d")
    
    def download():
        # setting global variables to local for a shorter req
        data       = GlobalVars.config.get('products', 'phytot_data')
        var = Library.GetVars(data)
        
        for nf in range(len(var['date'])):
            # dowloading data in /DATA
            date_range = "*/" + var['year'][nf] + "/" + var['month'][nf] + \
                "/dataset-uv-nrt-daily_" + var['date'][nf] + "*_*.nc"
            query_metadata = copernicusmarine.get(dataset_id=var['id'],username=var['user'],
                                             password=var['pwd'],output_directory=var['direct'],
                                             filter=date_range,overwrite='Y')
            
            for file in query_metadata.files:
                get_files = file.file_path
            exf,ff = Library.ExistingFile(get_files,var['date'][nf])
            
            # copying data in /Wrk
            if ff:
                req_cp = "cp " + ff +' '+ var['dir_wrk'] + var['date'][nf] + "_" + var['prod'] + '.nc'
                Library.execute_req(req_cp)
            
        return

class Copernicus_PHYEURO(Load,Create): 
    def __init__(self,fname,**kwargs):
        self.fname = fname
        data       = GlobalVars.config.get('products', 'phyeuro_data')
        var        = Library.GetVars(data)
        self.data_dir   = var['direct']
        self.lon_name='longitude'
        self.lat_name='latitude'
        self.d3_name = 'time'
        self.var_name = 'adt'
        self.var_units = 'm'
        self.u_name='ugos'
        self.v_name='vgos'
        self.u_units = 'm/s'
        self.v_units = 'm/s'
        self.cmap = cb.cbmap("cb.pregunta_r")
        if 'dayv' in kwargs: 
            self.date = datetime.datetime.strptime(kwargs['dayv'],"%Y-%m-%d")
        else: 
            self.date = datetime.datetime.strptime(var['date'][0],"%Y%m%d")
    
    def download(**kwargs):
        # setting global variables to local for a shorter req
        data       = GlobalVars.config.get('products', 'phyeuro_data')
        var        = Library.GetVars(data)
        
        if ('date' in kwargs):
            ddate = [kwargs['date']]
            year = [datetime.datetime.strftime(datetime.datetime.strptime(ddate[0],'%Y%m%d'),'%Y')]
            mo = [datetime.datetime.strftime(datetime.datetime.strptime(ddate[0],'%Y%m%d'),'%m')]
        else:
            ddate = var['date']
            year = var['year']
            mo = var['month']
        
        for nf in range(len(ddate)):
            # dowloading data in /DATA
            date_range = "*/" + year[nf] + "/" + mo[nf] + "/nrt_europe_allsat"+ \
            "_*_" + ddate[nf] + "_*.nc"
            query_metadata = copernicusmarine.get(dataset_id=var['id'],username=var['user'],
                                             password=var['pwd'],output_directory=var['direct'],
                                             filter=date_range,overwrite='Y')
            
            for file in query_metadata.files:
                get_files = file.file_path
            exf,ff = Library.ExistingFile(get_files,ddate[nf])
            
            if len(kwargs)==0 or 'cp' in kwargs:
                # copying data in /Wrk
                req_cp = "cp " + ff +' '+ var['dir_wrk'] + ddate[nf] + "_" + var['prod'] + '.nc'
                Library.execute_req(req_cp)
        return

class Copernicus_PHY_WIND(Load,Create): 
    def __init__(self,fname,**kwargs):
        self.fname = fname
        data       = GlobalVars.config.get('products', 'phy_wind_data')
        var        = Library.GetVars(data)
        self.data_dir   = var['direct']
        self.lon_name='lon'
        self.lat_name='lat'
        self.d3_name = 'time'
        self.var_name = 'wind_curl'
        self.var_units = 's-1'
        self.u_name='eastward_wind'
        self.v_name='northward_wind'
        self.u_units = 'm/s'
        self.v_units = 'm/s'
        self.cmap = plt.get_cmap('BrBG')
    
    def download(**kwargs):
        # setting global variables to local for a shorter req
        data       = GlobalVars.config.get('products', 'phy_wind_data')
        var        = Library.GetVars(data)
        
        if ('date' in kwargs):
            ddate = [kwargs['date']]
            year = [datetime.datetime.strftime(datetime.datetime.strptime(ddate[0],'%Y%m%d'),'%Y')]
            mo = [datetime.datetime.strftime(datetime.datetime.strptime(ddate[0],'%Y%m%d'),'%m')]
        else:
            ddate = var['date']
            year = [var['year']]
            mo = [var['month']]
        
        #closest hour
        tmpdate = datetime.datetime.strptime(ddate,'%Y%m%d%H%M')
        if tmpdate.minute >= 30:
            tmpdate.replace(second=0, microsecond=0, minute=0, hour=tmpdate.hour+1)
        else:
            tmpdate.replace(second=0, microsecond=0, minute=0)
        ddate = [datetime.datetime.strftime(tmpdate,'%Y%m%d%H')]
        
        for nf in range(len(ddate)):
            # dowloading data in /DATA
            date_range = "*/" + year[nf] + "/" + mo[nf] +\
                "/cmems_obs-wind_glo_phy_nrt_l4_0.125deg_PT1H_" + ddate[nf] + "_*.nc"
            query_metadata = copernicusmarine.get(dataset_id=var['id'],username=var['user'],
                                             password=var['pwd'],output_directory=var['direct'],
                                             filter=date_range,overwrite='Y')
            
            for file in query_metadata.files:
                get_files = file.file_path
            exf,ff = Library.ExistingFile(get_files,ddate[nf])
            
            if len(kwargs)==0 or 'cp' in kwargs:
                # copying data in /Wrk
                req_cp = "cp " + ff +' '+ var['dir_wrk'] + ddate[nf] + "_" + var['prod'] + '.nc'
                Library.execute_req(req_cp)
            
        return
    
class Copernicus_SST_L4(Load,Create):
    def __init__(self,fname,**kwargs):
        self.fname = fname
        self.lon_name='lon'
        self.lat_name='lat'
        self.d3_name = 'time'
        self.var_name = 'analysed_sst'
        self.var_units = 'Kelvin'
        self.K = True
        self.cmap = cb.cbmap("cb.iris_r")
    
    def download():       
        # setting global variables to local for a shorter req
        data       = GlobalVars.config.get('products', 'sst_l4_data')
        var        = Library.GetVars(data)
        
        for nf in range(len(var['date'])):
            # dowloading data in /DATA
            date_range = "*/" + var['year'][nf] + "/" + var['month'][nf] + "/" + var['date'][nf] + \
                "120000-UKMO-L4_GHRSST-SSTfnd-OSTIA-GLOB-v02.0-fv02.0.nc"
            query_metadata = copernicusmarine.get(dataset_id=var['id'],username=var['user'],
                                             password=var['pwd'],output_directory=var['direct'],
                                             filter=date_range,overwrite='Y')
            
            for file in query_metadata.files:
                get_files = file.file_path
            exf,ff = Library.ExistingFile(get_files,var['date'][nf])
            # copying data in /Wrk and rename
            if ff:
                req_cp = "cp " + ff +' '+var['dir_wrk'] + var['date'][nf] + "_" + var['prod'] + '.nc'
                Library.execute_req(req_cp)
            
        return

class Copernicus_SST_BAL_L4(Load,Create):
    def __init__(self,fname,**kwargs):
        self.fname = fname
        self.lon_name='lon'
        self.lat_name='lat'
        self.d3_name = 'time'
        self.var_name = 'analysed_sst'
        self.var_units = 'Kelvin'
        self.K = True
        self.cmap = cb.cbmap("cb.iris_r")
    
    def download():       
        # setting global variables to local for a shorter req
        data       = GlobalVars.config.get('products', 'sst_ball4_data')
        var        = Library.GetVars(data)
        
        for nf in range(len(var['date'])):
            # dowloading data in /DATA
            date_range = "*/" + var['year'][nf] + "/" + var['month'][nf] + "/" + var['date'][nf] + \
                "000000-DMI-L4_GHRSST-SSTfnd-DMI_OI-NSEABALTIC-v02.0-fv01.0.nc"
            query_metadata = copernicusmarine.get(dataset_id=var['id'],username=var['user'],
                                             password=var['pwd'],output_directory=var['direct'],
                                             filter=date_range,overwrite='Y')
            
            for file in query_metadata.files:
                get_files = file.file_path
            exf,ff = Library.ExistingFile(get_files,var['date'][nf])
            # copying data in /Wrk and rename
            if ff:
                req_cp = "cp " + ff +' '+var['dir_wrk'] + var['date'][nf] + "_" + var['prod'] + '.nc'
                Library.execute_req(req_cp)
            
        return
    
class Copernicus_SSS_L4(Load,Create):
    def __init__(self,fname,**kwargs):
        self.fname = fname
        data = GlobalVars.config.get('products', 'sss_l4_data')
        var        = Library.GetVars(data)
        self.data_dir   = var['direct']
        self.lon_name='lon'
        self.lat_name='lat'
        self.d3_name = 'time'
        self.var_name = 'sos'
        self.var_units = 'psu'
        self.cmap = cm_oc.cm.haline
    
    def download():       
        # setting global variables to local for a shorter req
        data       = GlobalVars.config.get('products', 'sss_l4_data')
        var        = Library.GetVars(data)
        
        for nf in range(len(var['date'])):
            #data available weekly / find closest available date
            # dowloading data in /DATA
            date_range = "*/" + var['year'][nf] + "/" + var['month'][nf] + "/" + \
                "dataset-sss-ssd-nrt-daily_"+var['date'][nf]+"T*.nc"
            query_metadata = copernicusmarine.get(dataset_id=var['id'],username=var['user'],
                                             password=var['pwd'],output_directory=var['direct'],
                                             filter=date_range,overwrite='Y')
            
            for file in query_metadata.files:
                get_files = file.file_path
            exf,ff = Library.ExistingFile(get_files,var['date'][nf])
            # copying data in /Wrk and rename
            if ff:
                req_cp = "cp " + ff +' '+var['dir_wrk'] + var['date'][nf] + "_" + var['prod'] + '.nc'
                Library.execute_req(req_cp)
            
        return
    
class Copernicus_CHL_L3(Load,Create):
    def __init__(self,fname,**kwargs):
        self.fname = fname
        self.lon_name='lon'
        self.lat_name='lat'
        self.d3_name = 'time'
        self.var_name = 'CHL'
        self.var_units = 'mg/m3'
        self.cmap = 'Greens'
        self.colnorm = 'PowerNorm'
    
    def download():
        # setting global variables to local for a shorter req
        data       = GlobalVars.config.get('products', 'chl_l3_data')
        var        = Library.GetVars(data)

        for nf in range(len(var['date'])):  
            # dowloading data in /DATA
            date_range = "*/" + var['year'][nf] + "/" + var['month'][nf] + "/" + var['date'][nf] + \
                "_cmems_obs-oc_glo_bgc-plankton_nrt_l3-multi-4km_P1D.nc"
            query_metadata = copernicusmarine.get(dataset_id=var['id'],username=var['user'],
                                             password=var['pwd'],output_directory=var['direct'],
                                             filter=date_range,overwrite='Y')
            
            for file in query_metadata.files:
                get_files = file.file_path
            exf,ff = Library.ExistingFile(get_files,var['date'][nf])
            # copying data in /Wrk
            if ff:
                req_cp = "cp " + ff+' '+var['dir_wrk'] + var['date'][nf] + "_" + var['prod'] + '.nc'
                Library.execute_req(req_cp)
            
        return
    
class Copernicus_CHL_L4(Load,Create):
    def __init__(self,fname,**kwargs):
        self.fname = fname
        self.lon_name='lon'
        self.lat_name='lat'
        self.d3_name = 'time'
        self.var_name = 'CHL'
        self.var_units = 'mg/m3'
        self.cmap = 'Greens'
        self.colnorm = 'PowerNorm'

    def download():
        # setting global variables to local for a shorter req
        data       = GlobalVars.config.get('products', 'chl_l4_data')
        var        = Library.GetVars(data)

        for nf in range(len(var['date'])):  
            # dowloading data in /DATA
            date_range = "*/" + var['year'][nf] + "/" + var['month'][nf] + "/" + var['date'][nf] + \
                "_cmems_obs-oc_glo_bgc-plankton_nrt_l4-gapfree-multi-4km_P1D.nc"
            query_metadata = copernicusmarine.get(dataset_id=var['id'],username=var['user'],
                                             password=var['pwd'],output_directory=var['direct'],
                                             filter=date_range,overwrite='Y')
            
            for file in query_metadata.files:
                get_files = file.file_path
            exf,ff = Library.ExistingFile(get_files,var['date'][nf])
            # copying data in /Wrk
            if ff:
                req_cp = "cp " + ff +' '+var['dir_wrk'] + var['date'][nf] + "_" + var['prod'] + '.nc'
                Library.execute_req(req_cp)
            
        return
    
class Copernicus_CHL_L4_DT(Load,Create):
    def __init__(self,fname,**kwargs):
        self.fname = fname
        self.lon_name='lon'
        self.lat_name='lat'
        self.d3_name = 'time'
        self.var_name = 'CHL'
        self.var_units = 'mg/m3'
        self.cmap = 'Greens'
        self.colnorm = 'PowerNorm'
 
    def download():
        # setting global variables to local for a shorter req
        data       = GlobalVars.config.get('products', 'chl_l4dt_data')
        var        = Library.GetVars(data)

        for nf in range(len(var['date'])):  
            # dowloading data in /DATA
            date_range = "*/" + var['year'][nf] + "/" + var['month'][nf] + "/" + var['date'][nf] + \
                "_cmems_obs-oc_glo_bgc-plankton_myint_l4-gapfree-multi-4km_P1D.nc"
            query_metadata = copernicusmarine.get(dataset_id=var['id'],username=var['user'],
                                             password=var['pwd'],output_directory=var['direct'],
                                             filter=date_range,overwrite='Y')
            
            for file in query_metadata.files:
                get_files = file.file_path
            exf,ff = Library.ExistingFile(get_files,var['date'][nf])

            # copying data in /Wrk
            if ff:
                req_cp = "cp " + ff+' '+var['dir_wrk'] + var['date'][nf] + "_" + var['prod'] + '.nc'
                Library.execute_req(req_cp)
            
        return

class Copernicus_CHL_BAL(Load,Create):
    def __init__(self,fname,**kwargs):
        self.fname = fname
        self.lon_name='lon'
        self.lat_name='lat'
        self.d3_name = 'time'
        self.var_name = 'CHL'
        self.var_units = 'mg/m3'
        self.cmap = 'Greens'
        self.colnorm = 'PowerNorm'
 
    def download():
        # setting global variables to local for a shorter req
        data       = GlobalVars.config.get('products', 'chl_bal_data')
        var        = Library.GetVars(data)

        for nf in range(len(var['date'])):  
            # dowloading data in /DATA
            date_range = "*/" + var['year'][nf] + "/" + var['month'][nf] + "/" + var['date'][nf] + \
                "_cmems_obs-oc_bal_bgc-plankton_nrt_l3-olci-300m_P1D.nc"
            query_metadata = copernicusmarine.get(dataset_id=var['id'],username=var['user'],
                                             password=var['pwd'],output_directory=var['direct'],
                                             filter=date_range,overwrite='Y')
            
            for file in query_metadata.files:
                get_files = file.file_path
            exf,ff = Library.ExistingFile(get_files,var['date'][nf])
            # copying data in /Wrk
            if ff:
                req_cp = "cp " + ff+' '+var['dir_wrk'] + var['date'][nf] + "_" + var['prod'] + '.nc'
                Library.execute_req(req_cp)
            
        return
    
class Copernicus_MEDSEA_WAVF(Load,Create):
    def __init__(self,fname,**kwargs):
        self.fname = fname
        self.lon_name='longitude'
        self.lat_name='latitude'
        self.d3_name = 'time'
        self.var_name = 'VHM0','VMDR','VTM01_WW','VHM0_WW' #Spectral significant wave height (Hm0), wind wave from direction, wind wave mean period, wind wave significant height
        self.var_units = 'm','degree','s','m'
        self.u_name='VSDX' #Stokes drift U
        self.v_name='VSDY' #Stokes drift V
        self.u_units = 'm/s'
        self.v_units = 'm/s'
        self.cmap = plt.get_cmap('Blues'),plt.get_cmap('gnuplot'),plt.get_cmap('Purples'),'PuRd'
        self.tit = 'Wave_Height','WWave_from_dir','WWave_mean_per','WWave_height'

    def download():
        # setting global variables to local for a shorter req
        data       = GlobalVars.config.get('products', 'medsea_wavf_data')
        var        = Library.GetVars(data)

        for nf in range(len(var['date'])):  
            # dowloading data in /DATA
            date_range = "*/" + var['year'][nf] + "/" + var['month'][nf] + "/" + var['date'][nf] + \
                "12_h-HCMR--WAVE-MEDWAM4-MEDATL-*.nc"
            query_metadata = copernicusmarine.get(dataset_id=var['id'],username=var['user'],
                                             password=var['pwd'],output_directory=var['direct'],
                                             filter=date_range,overwrite='Y')
            
            for file in query_metadata.files:
                get_files = file.file_path
            exf,ff = Library.ExistingFile(get_files,var['date'][nf])
            # copying data in /Wrk
            if ff:
                req_cp = "cp " + ff + ' '+var['dir_wrk'] + var['date'][nf] + "_" + var['prod'] + '.nc'
                Library.execute_req(req_cp)
            
        return

##################################################################
## CLS DATA 
##################################################################
class CLS_PHY(Load,Create):
    def __init__(self,fname,**kwargs):
        data       = GlobalVars.config.get('products', 'phy_cls_data')
        var        = Library.GetVars(data)
        self.data_dir   = var['direct']
        self.fname = fname
        self.lon_name='lon'
        self.lat_name='lat'
        self.d3_name = 'time'
        self.var_name = ' '
        self.var_units = 'None'
        self.u_name='u'
        self.v_name='v'
        self.u_units = 'm/s'
        self.v_units = 'm/s'
        self.cmap = plt.get_cmap('PRGn')
        if 'dayv' in kwargs: 
            self.date = datetime.datetime.strptime(kwargs['dayv'],"%Y-%m-%d")
        else: 
            self.date = datetime.datetime.strptime(var['datec'][0],"%Y%m%d")
 
    def download(**kwargs):
        # setting global variables to local for a shorter req
        data       = GlobalVars.config.get('products', 'phy_cls_data')
        var        = Library.GetVars(data)
        
        if ('date' in kwargs):
            ddate = [datetime.datetime.strftime(datetime.datetime.strptime(kwargs['date'],'%Y%m%d'),'%Y-%m-%d')]
            ddatec = [kwargs['date']]
        else:
            ddate = var['date']
            ddatec = var['datec']
        
        for nf in range(len(ddate)):
            # dowloading data in /DATA
            req_wget = GlobalVars.Lib['motulib']+"motuclient -q -u "+var['user']+" -p "+var['pwd']+" -m https://motu-"+var['arc']+"datastore.cls.fr/motu-web/Motu -s "\
                +var['id']+" -d "+var['name']+" -x "+str(var['Lon'][0])+" -X "+str(var['Lon'][1])+" -y "+str(var['Lat'][0])+" -Y "+str(var['Lat'][1])\
                    +" -t "+ddate[nf]+" -T "+ddate[nf]+ " --outputWritten netcdf -v"\
                        +" surface_eastward_geostrophic_sea_water_velocity -v surface_northward_geostrophic_sea_water_velocity -o "\
                            +var['direct']+"/ -f "+ddatec[nf]+"_"+var['prod']+".nc"
                
            Library.execute_req(req_wget)
            exf,ff = Library.ExistingFile(var['direct']+"/"+ddatec[nf]+"_"+var['prod']+".nc",ddate[nf])

            # copying data in /Wrk
            if ff:
                req_cp = "cp " + ff +' '+ var['dir_wrk']
                Library.execute_req(req_cp)
            # copy global data in data folder
            req_wget = GlobalVars.Lib['motulib']+"motuclient -q -u "+var['user']+" -p "+var['pwd']+" -m https://motu-"+var['arc']+"datastore.cls.fr/motu-web/Motu -s "\
                +var['id']+" -d "+var['name']+" -x 0 -X 359.98 -y -79 -Y 80"\
                    +" -t "+ddate[nf]+" -T "+ddate[nf]+ " --outputWritten netcdf -v"\
                        +" surface_eastward_geostrophic_sea_water_velocity -v surface_northward_geostrophic_sea_water_velocity -o "\
                            +var['direct']+"/ -f "+ddatec[nf]+"_"+var['prod']+".nc"
                
            Library.execute_req(req_wget)
            
        return

class CLS_CHL(Load,Create):
    def __init__(self,fname,**kwargs):
        self.fname = fname
        self.lon_name='longitude'
        self.lat_name='latitude'
        self.d3_name = 'time'
        self.var_name = 'chl_a'
        self.var_units = 'log10(mg/m3)'
        self.cmap = 'YlGnBu_r'
        self.log10 = True
        self.colnorm = 'PowerNorm'
    
    def download():
        # setting global variables to local for a shorter req
        data       = GlobalVars.config.get('products', 'chl_cls_data')
        var        = Library.GetVars(data)

        for nf in range(len(var['date'])):  
            # dowloading data in /DATA
            req_wget = GlobalVars.Lib['motulib']+"motuclient -q -u "+var['user']+" -p "+var['pwd']+" -m https://motu-"+var['arc']+"datastore.cls.fr/motu-web/Motu -s "\
                +var['id']+" -d "+var['name']+" -x "+str(var['Lon'][0])+" -X "+str(var['Lon'][1])+" -y "+str(var['Lat'][0])+" -Y "+str(var['Lat'][1])\
                    +" -t "+var['date'][nf]+" -T "+var['date'][nf]+ " --outputWritten netcdf -o "+var['direct']\
                        +"/ -f "+var['datec'][nf]+"_"+var['prod']+".nc"
                
            Library.execute_req(req_wget)
            exf,ff=Library.ExistingFile(var['direct']+"/"+var['datec'][nf]+"_"+var['prod']+".nc",var['date'][nf])

            # copying data in /Wrk
            if ff:
                req_cp = "cp " + ff +' '+ var['dir_wrk']
                Library.execute_req(req_cp)
            
        return
    
class CLS_CHL_5d(Load,Create):
    def __init__(self,fname,**kwargs):
        self.fname = fname
        self.lon_name='longitude'
        self.lat_name='latitude'
        self.d3_name = 'time'
        self.var_name = 'chl_a'
        self.var_units = 'log10(mg/m3)'
        self.cmap = 'YlGnBu_r'
        self.log10 = True
        self.colnorm = 'PowerNorm'
    
    def download():
        # setting global variables to local for a shorter req
        data       = GlobalVars.config.get('products', 'chl_cls5_data')
        var        = Library.GetVars(data)

        for nf in range(len(var['date'])):  
            # dowloading data in /DATA
            req_wget = GlobalVars.Lib['motulib']+"motuclient -q -u "+var['user']+" -p "+var['pwd']+" -m https://motu-"+var['arc']+"datastore.cls.fr/motu-web/Motu -s "\
                +var['id']+" -d "+var['name']+" -x "+str(var['Lon'][0])+" -X "+str(var['Lon'][1])+" -y "+str(var['Lat'][0])+" -Y "+str(var['Lat'][1])\
                    +" -t "+var['date'][nf]+" -T "+var['date'][nf]+ " --outputWritten netcdf -o "+var['direct']\
                        +"/ -f "+var['datec'][nf]+"_"+var['prod']+".nc"
                
            Library.execute_req(req_wget)
            exf,ff=Library.ExistingFile(var['direct']+"/"+var['datec'][nf]+"_"+var['prod']+".nc",var['date'][nf])

            # copying data in /Wrk
            if ff:
                req_cp = "cp " + ff +' '+ var['dir_wrk']
                Library.execute_req(req_cp)
            
        return

class CLS_CHL_10d(Load,Create):
    def __init__(self,fname,**kwargs):
        self.fname = fname
        self.lon_name='longitude'
        self.lat_name='latitude'
        self.d3_name = 'time'
        self.var_name = 'chl_a'
        self.var_units = 'log10(mg/m3)'
        self.cmap = 'YlGnBu_r'
        self.log10 = True
        self.colnorm = 'PowerNorm'

    def download():
        # setting global variables to local for a shorter req
        data       = GlobalVars.config.get('products', 'chl_cls10_data')
        var        = Library.GetVars(data)

        for nf in range(len(var['date'])):  
            # dowloading data in /DATA
            req_wget = GlobalVars.Lib['motulib']+"motuclient -q -u "+var['user']+" -p "+var['pwd']+" -m https://motu-"+var['arc']+"datastore.cls.fr/motu-web/Motu -s "\
                +var['id']+" -d "+var['name']+" -x "+str(var['Lon'][0])+" -X "+str(var['Lon'][1])+" -y "+str(var['Lat'][0])+" -Y "+str(var['Lat'][1])\
                    +" -t "+var['date'][nf]+" -T "+var['date'][nf]+ " --outputWritten netcdf -o "+var['direct']\
                        +"/ -f "+var['datec'][nf]+"_"+var['prod']+".nc"
                
            Library.execute_req(req_wget)
            exf,ff=Library.ExistingFile(var['direct']+"/"+"/"+var['datec'][nf]+"_"+var['prod']+".nc",var['date'][nf])

            # copying data in /Wrk
            if ff:
                req_cp = "cp " + ff +' '+ var['dir_wrk']
                Library.execute_req(req_cp)
            
        return

class CLS_SST(Load,Create):
    def __init__(self,fname,**kwargs):
        self.fname = fname
        self.lon_name='NbLongitudes'
        self.lat_name='NbLatitudes'
        self.d3_name = 'time'
        self.var_name = 'Grid_0001'
        self.var_units = 'degreesC'
        self.cmap = cm_oc.cm.thermal
 
    def download():
        # setting global variables to local for a shorter req
        data       = GlobalVars.config.get('products', 'sst_cls_data')
        var        = Library.GetVars(data)
        
        for nf in range(len(var['date'])):  
            # dowloading data in /DATA
            req_wget = GlobalVars.Lib['motulib']+"motuclient -q -u "+var['user']+" -p "+var['pwd']+" -m https://motu-"+var['arc']+"datastore.cls.fr/motu-web/Motu -s "\
                +var['id']+" -d "+var['name']+" -x "+str(var['Lon'][0])+" -X "+str(var['Lon'][1])+" -y "+str(var['Lat'][0])+" -Y "+str(var['Lat'][1])\
                    +" -t "+var['date'][nf]+" -T "+var['date'][nf]+ " --outputWritten netcdf -o "+var['direct']\
                        +"/ -f "+var['datec'][nf]+"_"+var['prod']+".nc"
                
            Library.execute_req(req_wget)
            exf,ff=Library.ExistingFile(var['direct']+"/"+var['datec'][nf]+"_"+var['prod']+".nc",var['date'][nf])

            # copying data in /Wrk
            if ff:
                req_cp = "cp " + ff +' '+ var['dir_wrk']
                Library.execute_req(req_cp)
            
        return
        
class CLS_SST_5d(Load,Create):
    def __init__(self,fname,**kwargs):
        self.fname = fname
        self.lon_name='NbLongitudes'
        self.lat_name='NbLatitudes'
        self.d3_name = 'time'
        self.var_name = 'Grid_0001'
        self.var_units = 'degreesC'
        self.cmap = cm_oc.cm.thermal
 
    def download():
        # setting global variables to local for a shorter req
        data       = GlobalVars.config.get('products', 'sst_cls5_data')
        var        = Library.GetVars(data)
        
        for nf in range(len(var['date'])):  
            # dowloading data in /DATA
            req_wget = GlobalVars.Lib['motulib']+"motuclient -q -u "+var['user']+" -p "+var['pwd']+" -m https://motu-"+var['arc']+"datastore.cls.fr/motu-web/Motu -s "\
                +var['id']+" -d "+var['name']+" -x "+str(var['Lon'][0])+" -X "+str(var['Lon'][1])+" -y "+str(var['Lat'][0])+" -Y "+str(var['Lat'][1])\
                    +" -t "+var['date'][nf]+" -T "+var['date'][nf]+ " --outputWritten netcdf -o "+var['direct']\
                        +"/ -f "+var['datec'][nf]+"_"+var['prod']+".nc"
                
            Library.execute_req(req_wget)
            exf,ff=Library.ExistingFile(var['direct']+"/"+var['datec'][nf]+"_"+var['prod']+".nc",var['date'][nf])

            # copying data in /Wrk
            if ff:
                req_cp = "cp " + ff +' '+ var['dir_wrk']
                Library.execute_req(req_cp)
            
        return
    
class CLS_SST_7d(Load,Create):
    def __init__(self,fname,**kwargs):
        self.fname = fname
        self.lon_name='NbLongitudes'
        self.lat_name='NbLatitudes'
        self.d3_name = 'time'
        self.var_name = 'Grid_0001'
        self.var_units = 'degreesC'
        self.cmap = cm_oc.cm.thermal
  
    def download():
        # setting global variables to local for a shorter req
        data       = GlobalVars.config.get('products', 'sst_cls7_data')
        var        = Library.GetVars(data)

        for nf in range(len(var['date'])):  
            # dowloading data in /DATA
            req_wget = GlobalVars.Lib['motulib']+"motuclient -q -u "+var['user']+" -p "+var['pwd']+" -m https://motu-"+var['arc']+"datastore.cls.fr/motu-web/Motu -s "\
                +var['id']+" -d "+var['name']+" -x "+str(var['Lon'][0])+" -X "+str(var['Lon'][1])+" -y "+str(var['Lat'][0])+" -Y "+str(var['Lat'][1])\
                    +" -t "+var['date'][nf]+" -T "+var['date'][nf]+ " --outputWritten netcdf -o "+var['direct']\
                        +"/ -f "+var['datec'][nf]+"_"+var['prod']+".nc"
                
            Library.execute_req(req_wget)
            exf,ff=Library.ExistingFile(var['direct']+"/"+var['datec'][nf]+"_"+var['prod']+".nc",var['date'][nf])

            # copying data in /Wrk
            if ff:
                req_cp = "cp " + ff +' '+ var['dir_wrk']
                Library.execute_req(req_cp)
            
        return
##################################################################
## SENTINEL PRODUCT 
##################################################################   
class Sentinel3_CHL(Load,Create):
    def __init__(self,fname,**kwargs):
        self.fname = fname
        self.lon_name='longitude'
        self.lat_name='latitude'
        self.d3_name = 'polygon'
        self.var_name = 'CHL_OC4ME'
        self.var_units = 'log10(mg/m3)'
        self.cmap = 'YlGnBu_r'
        self.log10 = True
        self.colnorm = 'PowerNorm'
  
    def download():
        # setting global variables to local for a shorter req
        data       = GlobalVars.config.get('products', 'chl_sen3_data')
        var        = Library.GetVars(data)

        for nf in range(len(var['date'])):
            for id0 in var['id0']:
                for id1 in var['id1']:
                    idp,name = Library.get_SEN3json(var['date'][nf],id0,id1)
                    if idp!=None:
                        # dowloading data in /DATA
                        token = Library.get_SEN3token(var['user'],var['pwd'])
                    
                        req_wget = "wget -nv -a "+var['logf']+" --header 'Authorization: Bearer "+token+"' "\
                        +"'https://catalogue.dataspace.copernicus.eu/odata/v1/Products("+idp\
                            +")/Nodes("+name+")/Nodes(geo_coordinates.nc)/$value' "\
                                +"-O "+var['direct']+"/"+var['date'][nf]+"_"+id0+"_geo_coordinates_"+id1+".nc --no-check-certificate"
                        Library.execute_req(req_wget)
                        req_wget = "wget -nv -a "+var['logf']+" --header 'Authorization: Bearer "+token+"' "\
                        +"'https://catalogue.dataspace.copernicus.eu/odata/v1/Products("+idp\
                            +")/Nodes("+name+")/Nodes(chl_oc4me.nc)/$value' "\
                                +"-O "+var['direct']+"/"+var['date'][nf]+"_"+id0+"_chl_oc4me_"+id1+".nc --no-check-certificate"
                        Library.execute_req(req_wget)            
            
            #Create unique netcdf file containing all polygons
            Library.CreateSEN3(var)
        return

##################################################################
## JAXA Himawari data
## https://www.eorc.jaxa.jp/ptree/
##################################################################  
class H8_SST_daily(Load,Create):
    def __init__(self,fname,**kwargs):
        self.fname = fname
        self.lon_name='lon'
        self.lat_name='lat'
        self.d3_name = 'time'
        self.var_name = 'sea_surface_temperature'
        self.var_units = 'Kelvin'
        self.K = True
        self.cmap = cm_oc.cm.thermal

    def download():
        # setting global variables to local for a shorter req
        data       = GlobalVars.config.get('products', 'sst_h8d_data')
        var        = Library.GetVars(data)

        for nf in range(len(var['date'])):  
            # dowloading data in /DATA
            req_wget = "wget -r --mirror -nd --directory-prefix=" + var['direct'] + \
                " -nv --no-proxy --user=" + var['user'] + " --password=" + var['pwd'] + \
                    " -a " + var['logf'] + " ftp://" + var['path'] + "/" + var['year'][nf] + var['month'][nf] +\
                        "/" + var['day'][nf]+ "/" + var['date'][nf] + \
                        "000000-JAXA-L3C_GHRSST-SSTskin-H09_AHI_NRT-v2.1_daily-v02.0-fv01.0.nc"

            Library.execute_req(req_wget)
            exf,ff = Library.ExistingFile(var['direct'] + "/" + var['date'][nf] + "000000-JAXA-L3C_GHRSST-SSTskin-H09_AHI_NRT-v2.1_daily-v02.0-fv01.0.nc",var['date'][nf])

            # copying data in /Wrk
            if ff:
                req_cp = "cp " + ff + ' '+var['dir_wrk'] + var['date'][nf] + "_" + var['prod'] + '.nc'
                Library.execute_req(req_cp)
            
        return
##################################################################
## EULERIAN AND LAGRANGIAN DIAGNOSTICS 
##################################################################   
class FTLE(Load,Create):
    def __init__(self,fname,**kwargs):
        self.fname = fname
        self.lon_name='lons'
        self.lat_name='lats'
        self.d3_name = 'time'
        self.var_name = 'ftle_lyap'
        self.var_units = 'day^{-1}'
        #customized colormap
        thresh = float(GlobalVars.config.get('plot_param','ftlethresh'))
        fmax = [float(x) for x in GlobalVars.config.get('plot_param','ftlemax').split(',')]
        if np.isnan(thresh):
            self.cmap = 'binary'
        else:
            oce = cm.get_cmap('binary', 256)
            newcolors = oce(np.linspace(0, 1, 256))
            whi = np.array([1, 1, 1, 1])
            nb = int(256/(fmax[0]/thresh))
            newcolors[:nb, :] = whi
            newcmp = ListedColormap(newcolors)
            self.cmap = newcmp   
        self.tit = 'Finite Time Lyapunov Exponent'
  
class LLADV(Load,Create):
    def __init__(self,fname,**kwargs):
        self.fname = fname
        self.lon_name='lons'
        self.lat_name='lats'
        self.d3_name = 'time'
        self.var_name = 'lonf','latf'
        self.var_units = '\Delta Lon [$^\circ$]','\Delta Lat [$^\circ$]'
        self.cmap = cm_oc.cm.curl,cm_oc.cm.delta
        self.tit = 'LonAdv','LatAdv'
    
class OWTRAJ(Load,Create):
    def __init__(self,fname,**kwargs):
        self.fname = fname
        self.lon_name='lons'
        self.lat_name='lats'
        self.d3_name = 'time'
        self.var_name = 'owdisp'
        self.var_units = ''
        self.cmap = cm_oc.cm.deep
        self.tit = 'Retention parameter'
    
class TIMEFROMBATHY(Load,Create):
    def __init__(self,fname,**kwargs):
        self.fname = fname
        self.lon_name='lon'
        self.lat_name='lat'
        self.d3_name = 'time'
        self.var_name = 'timfb','latfb','lonfb'
        self.var_units = 'days','Lat [$^\circ$]','Lon [$^\circ$]'
        self.cmap = 'turbo','turbo','turbo'
        bathy = str(abs(int(GlobalVars.Lag['bathylvl'])))
        self.tit = 'Timefrombathy_'+bathy+'m','Latfrombathy_'+bathy+'m','Lonfrombathy_'+bathy+'m'
 
class SSTADV(Load,Create):
    def __init__(self,fname,**kwargs):
        self.fname = fname
        self.lon_name='lon'
        self.lat_name='lat'
        self.d3_name = 'time'
        self.var_name = 'sstadv'
        self.var_units = 'degreesC'
        self.cmap = 'inferno'
        self.tit = 'Tracer advection'
    
class OW(Load,Create):
    def __init__(self,fname,**kwargs):
        self.fname = fname
        self.lon_name='lon'
        self.lat_name='lat'
        self.d3_name = 'time'
        self.var_name = 'ow'
        self.var_units = 'd$^{-2}$'
        self.cmap = cb.cbmap("cb.bird_r")
        self.tit = 'Okubo-Weiss parameter'
    
class KE(Load,Create):
    def __init__(self,fname,**kwargs):
        self.fname = fname
        self.lon_name='lon'
        self.lat_name='lat'
        self.d3_name = 'time'
        self.var_name = 'KE'
        self.var_units = 'cm$^{2}$.s$^{-2}$'
        self.cmap = 'turbo'
        self.tit = 'Kinetic energy'

class dSST(Load,Create):
    def __init__(self,fname,**kwargs):
        self.fname = fname
        self.lon_name='lon'
        self.lat_name='lat'
        self.d3_name = 'time'
        self.var_name = 'SSTgrad','SSTgradir'
        self.var_units = 'degreesC/pixel',' '
        self.cmap = 'afmhot_r','PRGn'
        self.tit = 'SSTgrad','SSTgradir'
    
##################################################################
## OTHER FIELDS
##################################################################   
class ETOPO():
    def loadnc(fname,**kwargs):
        file = Dataset(fname)
        lon = file.variables['lon'][:]
        lat = file.variables['lat'][:]
        if 'rlon' in kwargs:
            rlon = kwargs['rlon']
            indx = [index for index, item in enumerate(lon) if (item > rlon[0] and item < rlon[1])]
        else:
            indx = range(0,len(lon))
        if 'rlat' in kwargs:
            rlat = kwargs['rlat']
            indy = [index for index, item in enumerate(lat) if (item > rlat[0] and item < rlat[1])]
        else:
            indy = range(0,len(lat))
            

        z = file.variables['z'][indy,indx]
        lon = lon[indx]
        lat = lat[indy]
        file.close()
        field = {'lon':lon,'lat':lat,'z':z}
        return field
