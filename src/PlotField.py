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
import matplotlib.colors as colors
import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
import matplotlib 
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'

import glob, os
from netCDF4 import Dataset
from datetime import timedelta, datetime
from mpl_toolkits.basemap import Basemap

#PlotField.PlotField.Plot(cruise,pr,opt,type=None)

class PlotField():
    def Outputs(outp,pltarg):
        '''
        Call sub functions to save figures in other file format.
        
        Parameters
        ----------
        outp : output type
              Contains a list of output type. Must match the name of children 
              definition below.
        pltarg : Dictionnary with plot arguments
                Dictionnary with all plot arguments needed to create temporary 
                     figures to be saved in "outp" format.

        Returns
        -------
        None.
        
        Definitions
        -------
        kml(): create kmz file for each figure for Google earth vizionning 
        '''
        def kml(pltarg):
            '''
            Parameters
            ----------
            pltarg : Dictionnary with all plot arguments needed for kml creation
            '''
            lon_stat = [float(x) for x in GlobalVars.config.get('stations','coordlon').split(',')]
            lat_stat = [float(x) for x in GlobalVars.config.get('stations','coordlat').split(',')]
            stat = [str(x) for x in GlobalVars.config.get('stations','name').split(',')]
            #overlay 1
            figkml, axkml = Library.gearth_fig(llcrnrlon=pltarg['llon'],
                                  llcrnrlat=pltarg['llat'],
                                  urcrnrlon=pltarg['ulon'],
                                  urcrnrlat=pltarg['ulat'])
            if pltarg['var'].ndim==3:
                for kk in range(np.shape(pltarg['var'])[2]):
                    axkml.pcolormesh(pltarg['lon'][:,:,kk],pltarg['lat'][:,:,kk],
                                     pltarg['var'][:,:,kk],cmap=pltarg['cmap'],
                                          vmin=pltarg['min'],vmax=pltarg['max'])
            else:
                if pltarg['lon'].ndim==3:
                    axkml.pcolormesh(np.squeeze(pltarg['lon']),np.squeeze(pltarg['lat']),pltarg['var'],cmap=pltarg['cmap'],
                                      vmin=pltarg['min'],vmax=pltarg['max'])
                else:
                    axkml.pcolormesh(pltarg['lon'],pltarg['lat'],pltarg['var'],cmap=pltarg['cmap'],
                                  vmin=pltarg['min'],vmax=pltarg['max'])
            axkml.scatter(lon_stat,lat_stat,s=5,color='k',zorder=1)
            axkml.plot(lon_stat,lat_stat,'-',color='k',zorder=1)
            if pltarg['ii'] > 0:
                lons = [x - (10/60) for x in lon_stat]
                lats = [x - (5/60) for x in lat_stat]
                for i, txt in enumerate(stat):
                    plt.annotate(txt,(lons[i],lats[i]),size=12)
            axkml.set_axis_off()
            fn = pltarg['file'][:-1]+'_tmp1.png'
            figkml.savefig(fn,transparent=True,format='png')
            plt.close(figkml)
            figs=[fn]
            names=[pltarg['title']]
            
            #optional overlay 2
            if 'u' in pltarg:
                step = pltarg['step']
                figkml2, axkml2 = Library.gearth_fig(llcrnrlon=pltarg['llon'],
                                      llcrnrlat=pltarg['llat'],
                                      urcrnrlon=pltarg['ulon'],
                                      urcrnrlat=pltarg['ulat'])
                axkml2.quiver(pltarg['lon'][::step[1],::step[0]],pltarg['lat'][::step[1],::step[0]],
                              pltarg['u'][::step[1],::step[0]],pltarg['v'][::step[1],::step[0]],
                              linewidth=0.05,scale=pltarg['uvscale'])
                axkml2.set_axis_off()
                fn = pltarg['file'][:-1]+'_tmp2.png'
                figkml2.savefig(fn,transparent=True,format='png')
                plt.close(figkml2)
                figs.append(fn)
                names.append('UV '+pltarg['title'])
            
            #make kml
            Library.make_kml(llcrnrlon=pltarg['llon'], llcrnrlat=pltarg['llat'],
                      urcrnrlon=pltarg['ulon'], urcrnrlat=pltarg['ulat'],
                      figs=figs,kmzfile=pltarg['file']+'kmz',name=names)
            txt = os.path.basename(pltarg['file'])+'kmz created.'
            Library.Done(txt)
                
            return
        
        def Mergekml(pltarg):
            '''
            Parameters
            ----------
            pltarg : Dictionnary with all plot arguments needed for kml creation
            '''
            count=0
            for ii in range(0,pltarg['nd']):
                iLL = ii + count
                figs = glob.glob(GlobalVars.Dir['dir_wrk']+'/*'+str(ii)+'_tmp*.png')
                names=[os.path.basename(x)[:-9] for x in figs]
                #make kml
                Library.make_kml(llcrnrlon=pltarg['Lon'][iLL],llcrnrlat=pltarg['Lat'][iLL],
                          urcrnrlon=pltarg['Lon'][iLL+1],urcrnrlat=pltarg['Lat'][iLL+1],
                          figs=figs,kmzfile=GlobalVars.Dir['dir_wrk']+'Figures_oftheday_'+str(ii)+'.kmz',
                          name=names)
                count = count+1
                
            return 
        
        for out in outp:
            eval(out+"(pltarg)")
            
        return
    
    def Plot(cruise,Field, *args, **kwargs):
        def MedZONEX(mymap,ii,nc,**kwargs):
            col =['r','k']
            lin = ['dashed','solid']
            wid = [0.5,0.5]
            rep = str(GlobalVars.config.get('MedZONEX','direct'))
            files = sorted(glob.glob(GlobalVars.Dir['dir_data']+rep+'/*.shp'))
            co=0
            for ff in files:
                res = mymap.readshapefile(ff[:-4],\
                                        os.path.basename(ff[:-4]),color=col[co],linewidth=wid[co])
                seg = res[-1]
                seg.set_linestyle(lin[co])
                co = co+1
            
            return
        
        def AlgEEZ(mymap,ii,nc,**kwargs):
            rep = str(GlobalVars.config.get('AlgEEZ','direct'))
            listf = glob.glob(GlobalVars.Dir['dir_data']+rep+'/*.shp')
            for ff in listf:
                res = mymap.readshapefile(ff[:-4],ff[:-4],color='k',linewidth=0.5)
                seg = res[-1]
                seg.set_linestyle('--')

            return

        def SMODEb(mymap,ii,nc,**kwargs):
            """ Add polygon to show S-MODE IOP2 operations area."""
                
            coord = [[38.342, -126.25], [37.707, -123.99], [37.75, -123.354], [37.00, -122.92], [36.337, -124.36], [36.00, -124.16],[35.60, -125.515]]
            coord.append(coord[0]) #repeat the first point to create a 'closed loop'
            
            ys, xs = zip(*coord) #create lists of x and y values
            x,y = mymap(xs,ys)
            mymap.plot(x,y,color='r',linewidth=1,zorder=3)
            # mark a known place to help us geo-locate ourselves
            SF_lon=-(122+25/60)
            SF_lat= 37+60/60
            xSF,ySF = mymap(SF_lon,SF_lat)
            txSF,tySF = mymap(SF_lon-5/60, SF_lat+5/60)
            mymap.plot(xSF,ySF, 'o',color='k',markersize=3, zorder=3)
            plt.annotate('San Francisco', xy=(txSF,tySF))
            return
        
        def moorings(mymap,ii,nc,**kwargs):
            lon = [float(x) for x in GlobalVars.config.get('moorings','lon').split(',')]
            lat = [float(x) for x in GlobalVars.config.get('moorings','lat').split(',')]
            x,y = mymap(lon,lat)
            mymap.plot(x,y,'x',color='r',zorder=1)
            return
        
        def Locations(mymap,ii,nc,**kwargs):
            lon = [float(x) for x in GlobalVars.config.get('Locations','lon').split(',')]
            lat = [float(x) for x in GlobalVars.config.get('Locations','lat').split(',')]
            col = [str(x) for x in GlobalVars.config.get('Locations','col').split(',')]
            x,y = mymap(lon,lat)
            for jj in range(0,len(col)):
                mymap.plot(x[jj],y[jj],'x',color=col[jj],zorder=1)
            return
        
        def Transects(mymap,ii,nc,**kwargs):
            file = [str(x) for x in GlobalVars.config.get('Transects','file').split(',')]
            col = [str(x) for x in GlobalVars.config.get('Transects','col').split(',')]
            for jj in range(len(file)):
                x,y = [],[]
                filei = glob.glob(GlobalVars.Dir['dir_data']+'TRANSECTS/'+
                              GlobalVars.Param['cruise']+'/'+file[jj])[0]
                with open(filei,'r') as f:
                    for line in f:
                        xtmp, ytmp, ztmp = line.split()
                        x.append(float(xtmp))
                        y.append(float(ytmp))
                xx,yy = mymap(x,y)
                mymap.plot(xx,yy,'-',color=col[jj],zorder=1,linewidth=1.5,label=file[jj])
        
            return
        
        def LIONmooring(mymap,ii,nc,**kwargs):
            lon = [float(x) for x in GlobalVars.config.get('LIONmooring','lon').split(',')]
            lat = [float(x) for x in GlobalVars.config.get('LIONmooring','lat').split(',')]
            x,y = mymap(lon,lat)
            mymap.plot(x,y,'x',color='r',zorder=1)
            return
        
        def MooseT02(mymap,ii,nc,**kwargs):
            lon = [float(x) for x in GlobalVars.config.get('MooseT02','lon').split(',')]
            lat = [float(x) for x in GlobalVars.config.get('MooseT02','lat').split(',')]
            x,y = mymap(lon,lat)
            mymap.plot(x,y,'--',color='r',zorder=1,linewidth=0.8)
            return
        
        def delmar(mymap,ii,nc,**kwargs):
            col =['k','r','orange','m','b','k','k']
            lin = ['dashed','solid','solid','solid','solid','dashed','dotted']
            wid = [0.5,0.5,0.5,0.5,0.5,0.5,0.5]
            rep = str(GlobalVars.config.get('delmar','direct'))
            listf = sorted(glob.glob(GlobalVars.Dir['dir_data']+rep))
            co=0
            for i in listf:
                if not i.startswith('.'):
                    file = glob.glob(i+'/SHAPE/*.shp')
                    res = mymap.readshapefile(file[0][:-4],\
                                        os.path.basename(file[0][:-4]),color=col[co],linewidth=wid[co])
                    seg = res[-1]
                    seg.set_linestyle(lin[co])
                    co = co+1
            return
        
        def ZEE(mymap,ii,nc,**kwargs):
            files = [str(x) for x in GlobalVars.config.get('ZEE','files').split(',')]
            for f in files:
                fname = GlobalVars.Dir['dir_data']+'ZEE/'+f+'.txt'
                with open(fname,'r') as data:
                    x = []
                    y = []
                    for line in data:
                        p = line.split(',')
                        x.append(float(p[0]))
                        y.append(float(p[1]))
                    xx,yy = mymap(x,y)
                    mymap.plot(xx,yy,'-',color='k',zorder=1,linewidth=0.5);
            return
        
        def SWOTswath(mymap,ii,nc,**kwargs):
            file = str(GlobalVars.config.get('SWOTswath','file'))
            swath = sio.loadmat(GlobalVars.Dir['dir_data']+'SWOT/'+file)
            lon,lat = swath['lon'],swath['lat']
            for i in range(0,len(lon)-1):
                if abs(lon[i]-lon[i+1])<10:
                    x,y = mymap(lon[i:(i+1)+1],lat[i:(i+1)+1])
                    mymap.plot(x,y,'-',color='gold',zorder=1,linewidth=0.5);
            return
            
        def glider(mymap,ii,nc,**kwargs):
            file = [str(x) for x in GlobalVars.config.get('glider','file').split(',')]
            col = ['r','tab:orange']
            for ii in range(len(file)):
                filei = glob.glob(GlobalVars.Dir['dir_data']+'GLIDER/'+
                              GlobalVars.Param['cruise']+'/'+file[ii])[0]
                fn = Dataset(filei)
                lon = fn.variables[str(GlobalVars.config.get('glider','varnLon'))][:]
                lat = fn.variables[str(GlobalVars.config.get('glider','varnLat'))][:]
                time = fn.variables[str(GlobalVars.config.get('glider','varnTime'))][:]
                #remove NaN values if any
                time = [x for x in time if ~np.isnan(x)]
                lon = [x for x in lon if ~np.isnan(x)]
                lat = [x for x in lat if ~np.isnan(x)]
                #check for daily trajectory
                start = datetime(1,1,1,0,0,0)
                delta = [(timedelta(x)-timedelta(days=367)) for x in time]# Create a time delta object from the number of days
                dates = [start+dd for dd in delta]
                dday = GlobalVars.all_dates['ref_all']
                if 'dayv' in kwargs:
                    dayv = datetime.strptime(kwargs['dayv'],'%Y-%m-%d')
                else:
                    dayv = datetime.strptime(dday[nc],'%Y%m%d')
                index = [i for i,dt in enumerate(dates) if 
                             (dt <= dayv+timedelta(days=1) and dt >= dayv)]
                lon = [lon[i] for i in index]
                lat = [lat[i] for i in index]
                # plot glider for this day
                x,y = mymap(lon,lat)
                mymap.plot(x,y,'-',color=col[ii],zorder=1,linewidth=1.5,label=file[ii])
                plt.legend(loc='lower right',fontsize='small')
            return
        
        def stations(mymap,ii,nc,**kwargs):
            lon_stat = [float(x) for x in GlobalVars.config.get('stations','coordlon').split(',')]
            lat_stat = [float(x) for x in GlobalVars.config.get('stations','coordlat').split(',')]
            stat = [str(x) for x in GlobalVars.config.get('stations','name').split(',')]
            x_stat,y_stat = mymap(lon_stat,lat_stat)
            mymap.scatter(x_stat,y_stat,s=4,color='k',zorder=1)
            mymap.plot(x_stat,y_stat,'-',color='k',zorder=1,linewidth=0.5)
            lons = [x - (30/60) for x in lon_stat]
            lats = [x - (5/60) for x in lat_stat]
            tx,ty = mymap(lons,lats)
            if ii >= 0:
                for i, txt in enumerate(stat):
                    plt.annotate(txt,(tx[i],ty[i]),size=7,color='r')
            return
        
        def waypoints(mymap,ii,nc,**kwargs):
            lon_wayp = [float(x) for x in GlobalVars.config.get('waypoints','waylon').split(',')]
            lat_wayp = [float(x) for x in GlobalVars.config.get('waypoints','waylat').split(',')]
            x_wayp,y_wayp = mymap(lon_wayp,lat_wayp)
            mymap.plot(x_wayp,y_wayp,'-',color='r',zorder=1)
            return
        
        def cities(mymap,ii,nc,**kwargs):
            clon = [float(x) for x in GlobalVars.config.get('cities','clon').split(',')]
            clat = [float(x) for x in GlobalVars.config.get('cities','clat').split(',')]
            x_c,y_c = mymap(clon,clat)
            mymap.plot(x_c,y_c,'x',color='r',zorder=1,markersize=20)
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
        if not GlobalVars.config.get('plot_param',Field+'min'):
            fmin,fmax =[],[]
        else:
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
            if 'dayv' in kwargs:
                dayv = datetime.strftime(
                    datetime.strptime(kwargs['dayv'],'%Y-%m-%d')
                    ,'%Y%m%d')
                fname = glob.glob(GlobalVars.Dir['dir_wrk']+'/'+dayv+'*'+nf+'*.nc')
            else:
                fname = glob.glob(GlobalVars.Dir['dir_wrk']+'/*'+nf+'*.nc')
            fname.sort(key=os.path.getmtime)
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
                        if fmin==[]:
                            vmin,vmax = np.nanmin(var[nv-nv0]),np.nanmax(var[nv-nv0])
                        else:
                            vmin,vmax = fmin[ii+nv],fmax[ii+nv]
                            
                        if lon.ndim==1:
                            long, latg = np.meshgrid(lon, lat)
                            (x,y)=mymap(long,latg)
                            if 'colnorm' in field:
                                if field['colnorm']=='PowerNorm':
                                    cax1=mymap.pcolormesh(x,y,np.squeeze(var[nv-nv0]),norm=colors.PowerNorm(gamma=0.5,vmin=vmin,vmax=vmax),cmap=cm[nv-nv0],zorder=-1)
                            else:
                                cax1=mymap.pcolormesh(x,y,np.squeeze(var[nv-nv0]),cmap=cm[nv-nv0],zorder=-1,vmin=vmin,vmax=vmax)   
                        elif lon.ndim==3:
                            long,latg = lon,lat
                            for nP in range(np.shape(lon)[2]):
                                tmplong,tmplatg = lon[:,:,nP],lat[:,:,nP]
                                (x,y)=mymap(tmplong,tmplatg)
                                tmpvar = var[nv-nv0][:,:,nP]
                                if 'colnorm' in field:
                                    if field['colnorm']=='PowerNorm':
                                        cax1=mymap.pcolormesh(x,y,np.squeeze(tmpvar),norm=colors.PowerNorm(gamma=0.5,vmin=vmin,vmax=vmax),cmap=cm[nv-nv0],zorder=-1)
                                else:
                                    cax1=mymap.pcolormesh(x,y,np.squeeze(tmpvar),cmap=cm[nv-nv0],zorder=-1,vmin=vmin,vmax=vmax) 
                                   
                        Library.printInfo('min: '+str(np.nanmin(var[nv-nv0]))+', max: '+str(np.nanmax(var[nv-nv0])))
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
                                if len(opt)>0:
                                    if 'dayv' in kwargs:
                                        dayv = kwargs['dayv']
                                        eval(opt + "(mymap,ii,nc-1,dayv=dayv)")
                                    else:
                                        eval(opt + "(mymap,ii,nc-1)")
                        plt.savefig(filefig,bbox_inches='tight',dpi=GlobalVars.Fig['dpi'])
                        txt = '\t\t Figure: ' + filefig
                        Library.Logfile(txt)
                        plt.close(fig)
                        #extra outputs
                        out = [str(x) for x in GlobalVars.config.get('plot_options','outopt').split(',')]
                        if out:
                            pltarg = {'lon':long,'lat':latg,'var':np.squeeze(var[nv-nv0]),
                                      'cmap':cm[nv-nv0],'min':vmin,
                                      'max':vmax,'llon':Lon[iLL],
                                      'ulon':Lonp,'llat':Lat[iLL],'ulat':Lat[iLL+1],
                                      'file':filefig[:-3],'title':titlefig,
                                      'nd':nb_domain,'Lon':Lon,'Lonp':Lonp,'Lat':Lat,
                                      'ii':ii}
                            if 'u' and 'v' in field:
                                pltarg['u'] = u
                                pltarg['v'] = v
                                pltarg['uvscale'] = uvscale[ii]
                                pltarg['step'] = step
                            PlotField.Outputs(out,pltarg)

                    nv0 = nv
            nfig = str(nc) + ' figure(s) done.'
            Library.Done(nfig)
