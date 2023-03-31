#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 11:41:58 2022

@author: lrousselet
"""
import GlobalVars
import Library
import Fields

import matplotlib
matplotlib.use('Agg')
import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
import matplotlib 
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'

import glob, os
from mpl_toolkits.basemap import Basemap

class PlotField():
    def Plot(cruise,Field, *args, **kwargs):
        def moorings(mymap):
            lon = [float(x) for x in GlobalVars.config.get('moorings','lon').split(',')]
            lat = [float(x) for x in GlobalVars.config.get('moorings','lat').split(',')]
            x,y = mymap(lon,lat)
            mymap.plot(x,y,'x',color='r',zorder=1)
            return
        
        def stations(mymap):
            lon_stat = [float(x) for x in GlobalVars.config.get('stations','coordlon').split(',')]
            lat_stat = [float(x) for x in GlobalVars.config.get('stations','coordlat').split(',')]
            x_stat,y_stat = mymap(lon_stat,lat_stat)
            mymap.scatter(x_stat,y_stat,s=8,color='r',zorder=1)
            return
        
        def waypoints(mymap):
            lon_wayp = [float(x) for x in GlobalVars.config.get('waypoints','waylon').split(',')]
            lat_wayp = [float(x) for x in GlobalVars.config.get('waypoints','waylat').split(',')]
            x_wayp,y_wayp = mymap(lon_wayp,lat_wayp)
            mymap.plot(x_wayp,y_wayp,'-',color='r',zorder=1)
            return
        
        def cities(mymap):
            clon = [float(x) for x in GlobalVars.config.get('cities','clon').split(',')]
            clat = [float(x) for x in GlobalVars.config.get('cities','clat').split(',')]
            #cname = [str(x) for x in GlobalVars.config.get('cities','cname').split(',')]
            x_c,y_c = mymap(clon,clat)
            mymap.plot(x_c,y_c,'x',color='r',zorder=1,markersize=20)
            #plt.annotate(cname,(x_c,y_c))
            return

        def PeriodicLon(Lon,lon):
            minlon, maxlon = np.min(lon),np.max(lon)
            if Lon[1] < Lon[0]: 
                Lonp = Lon[1] + 360
                if Lonp > maxlon and minlon < 0:
                    lon[lon<0] += 360
                    sl = np.argsort(lon)
                else:
                    sl = None
            else:
                sl = None
                Lonp = Lon[1]
            if minlon<0 and maxlon==180.0:
                lon[lon<0] += 360
            if minlon>Lon[1]:
                lon -= 360
            return sl,Lonp
            
        #get param
        GlobalVars.configIni(cruise)
        GlobalVars.directories(cruise)
        Lon = [float(x) for x in GlobalVars.config.get('cruise_param','Lon').split(',')]
        Lat = [float(x) for x in GlobalVars.config.get('cruise_param','Lat').split(',')]
        fmin = [float(x) for x in GlobalVars.config.get('plot_param',Field+'min').split(',')]
        fmax = [float(x) for x in GlobalVars.config.get('plot_param',Field+'max').split(',')]
        funit = [str(x) for x in GlobalVars.config.get('plot_param',Field+'unit').split(',')]
        par = [eval(x) for x in GlobalVars.config.get('plot_param','parallels').split(';')]
        mer = [eval(x) for x in GlobalVars.config.get('plot_param','meridians').split(';')]
        
        if kwargs['type']== 'Lagrangian' or kwargs['type']== 'Eulerian':
            fprod = [Field]
        elif 'LATEXtools' in kwargs and kwargs['LATEXtools']==True:
            fprod = [str(x) for x in GlobalVars.latexini.get('products',Field+'prod').split(',')]
        else:
            fprod = [str(x) for x in GlobalVars.config.get('products',Field+'prod').split(',')]
        
        # how many plots ?
        nb_domain = GlobalVars.config.getint('cruise_param','nb_domain')
        nc = 0

        for nf in fprod:
            fname = glob.glob(GlobalVars.Dir['dir_wrk']+'/*'+nf+'*.nc')
            for file in fname:
                #load file
                field = eval('Fields.'+nf+'(file).loadnc()')
                lon = field['lon']
                lat = field['lat']
                var = field['var']
                cm = field['cm']
                if 'u' and 'v' in field:
                    u = field['u']
                    v = field['v']
                #add a complement to fig title
                if 'title' in field:
                    tit = field['title']
                    if type(tit) is not tuple: tit = [tit]
                # if multiple variable in a netcdf
                if isinstance(var,tuple):
                    varnb = len(var)
                else: 
                    varnb = 1
                    var = [var]
                    cm = [cm]
                #check for Periodic longitude boundary
                sl,Lonp = PeriodicLon(Lon,lon)
                if sl is not None:
                    lon = lon[sl]
                    for nv in range(0,varnb): var[nv] = var[nv][:,sl]
                    if 'u' and 'v' in field: u,v = u[:,sl],v[:,sl]
                count = 0
                #start loop on figures
                nv0 = 0
                for ii in range(0,nb_domain):
                    iLL = ii + count
                    count = count + 1
                    for nv in range(nv0,nv0+varnb):
                        #count figure
                        nc = nc + 1
                        #define domain
                        lontmp = Lon[iLL:(iLL+1)+1]
                        sl,Lonp = PeriodicLon(lontmp,lon)
                        mymap=Basemap(projection='merc',llcrnrlat=Lat[iLL],urcrnrlat=Lat[iLL+1],llcrnrlon=Lon[iLL],urcrnrlon=Lonp,resolution='h')
                        fig=plt.figure()
                        if len(np.shape(lon))==2:
                            long,latg = lon,lat
                        else:
                            long, latg = np.meshgrid(lon, lat)
                        (x,y)=mymap(long,latg)
                        cax1=mymap.pcolormesh(x,y,np.squeeze(var[nv-nv0]),cmap=cm[nv-nv0],zorder=-1,vmin=fmin[ii+nv],vmax=fmax[ii+nv])
                        if np.isnan(var[nv-nv0]).all()==False:
                            cbar1=fig.colorbar(cax1, orientation='vertical',shrink=0.5)
                            cbar1.ax.set_ylabel(funit[nv-nv0])
                        if 'u' and 'v' in field:
                            uvscale = [float(x) for x in GlobalVars.config.get('plot_param',Field+'uv').split(',')]
                            step = [int(x) for x in GlobalVars.config.get('plot_param',Field+'uvstep').split(',')]
                            Q = mymap.quiver(x[::step[1],::step[0]],y[::step[1],::step[0]],u[::step[1],::step[0]],v[::step[1],::step[0]],linewidth=0.05,scale=uvscale[ii])
                            plt.quiverkey(Q, 1.25, 1.05, 0.5, '0.5 m/s', labelpos='S')
                        mymap.drawcoastlines()
                        mymap.fillcontinents(color='0.83',lake_color='0.83',zorder=2)
                        mymap.drawparallels(par[ii],labels=[1,0,0,0],fontsize=10)
                        mymap.drawmeridians(mer[ii],labels=[0,0,0,1],fontsize=10)
                        if 'tit' in locals():
                            titlefig = os.path.basename(file)[0:8]+' '+tit[nv-nv0]
                            if len(tit)>1:
                                filefig = file[:-3] +'_'+tit[nv-nv0]+'_'+str(ii)+'.png'
                            else:
                                filefig = file[:-3] +'_'+str(ii)+'.png'
                        else:
                            titlefig = os.path.basename(file)[:-3]
                            filefig = file[:-3] +'_'+str(ii)+'.png'
                        plt.title(titlefig)
                        # plot options
                        for arg in args:
                            for opt in arg:
                                eval(opt + "(mymap)")
                        plt.savefig(filefig,bbox_inches='tight',dpi=150)
                        txt = '\t\t Figure: ' + filefig
                        Library.Logfile(txt)
                        plt.close(fig)
                    nv0 = nv
                    
                    
            nfig = str(nc) + ' figure(s) done.'
            Library.Done(nfig)
