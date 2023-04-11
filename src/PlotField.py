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
            #overlay 1
            figkml, axkml = Library.gearth_fig(llcrnrlon=pltarg['llon'],
                                  llcrnrlat=pltarg['llat'],
                                  urcrnrlon=pltarg['ulon'],
                                  urcrnrlat=pltarg['ulat'])
            axkml.pcolormesh(pltarg['lon'],pltarg['lat'],pltarg['var'],cmap=pltarg['cmap'],
                                  vmin=pltarg['min'],vmax=pltarg['max'])
            axkml.set_axis_off()
            fn = GlobalVars.Dir['dir_wrk']+'overlaytmp1.png'
            figkml.savefig(fn,transparent=False,format='png')
            figs=[fn]
            
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
                fn = GlobalVars.Dir['dir_wrk']+'overlaytmp2.png'
                figkml2.savefig(fn,transparent=False,format='png')
                figs.append(fn)
            
            #make kml
            Library.make_kml(llcrnrlon=pltarg['llon'], llcrnrlat=pltarg['llat'],
                      urcrnrlon=pltarg['ulon'], urcrnrlat=pltarg['ulat'],
                      figs=figs,kmzfile=pltarg['file']+'kmz',name=pltarg['title'])
            txt = os.path.basename(pltarg['file'])+'kmz created.'
            Library.Done(txt)
            
            #delete tmp png files
            for ff in figs:
                os.remove(ff)
                
            return 
        
        for out in outp:
            eval(out+"(pltarg)")
            
        return
    
    def Plot(cruise,Field, *args, **kwargs):
        def MedZONEX(mymap):
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
        
        def AlgEEZ(mymap):
            rep = str(GlobalVars.config.get('AlgEEZ','direct'))
            listf = glob.glob(GlobalVars.Dir['dir_data']+rep+'/*.shp')
            for ff in listf:
                res = mymap.readshapefile(ff[:-4],ff[:-4],color='k',linewidth=0.5)
                seg = res[-1]
                seg.set_linestyle('--')

            return

        def SMODEb(mymap):
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
        
        def moorings(mymap):
            lon = [float(x) for x in GlobalVars.config.get('moorings','lon').split(',')]
            lat = [float(x) for x in GlobalVars.config.get('moorings','lat').split(',')]
            x,y = mymap(lon,lat)
            mymap.plot(x,y,'x',color='r',zorder=1)
            return
        
        def LIONmooring(mymap):
            lon = [float(x) for x in GlobalVars.config.get('LIONmooring','lon').split(',')]
            lat = [float(x) for x in GlobalVars.config.get('LIONmooring','lat').split(',')]
            x,y = mymap(lon,lat)
            mymap.plot(x,y,'x',color='r',zorder=1)
            return
        
        def MooseT02(mymap):
            lon = [float(x) for x in GlobalVars.config.get('MooseT02','lon').split(',')]
            lat = [float(x) for x in GlobalVars.config.get('MooseT02','lat').split(',')]
            x,y = mymap(lon,lat)
            mymap.plot(x,y,'--',color='r',zorder=1,linewidth=0.8)
            return
        
        def delmar(mymap):
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
        
        def ZEE(mymap):
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
        
        def SWOTswath(mymap):
            file = str(GlobalVars.config.get('SWOTswath','file'))
            swath = sio.loadmat(GlobalVars.Dir['dir_data']+'SWOT/'+file)
            lon,lat = swath['lon'],swath['lat']
            for i in range(0,len(lon)-1):
                if abs(lon[i]-lon[i+1])<10:
                    x,y = mymap(lon[i:(i+1)+1],lat[i:(i+1)+1])
                    mymap.plot(x,y,'-',color='gold',zorder=1,linewidth=0.5);
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
                        #extra outputs
                        out = [str(x) for x in GlobalVars.config.get('plot_options','outopt').split(',')]
                        if out:
                            pltarg = {'lon':long,'lat':latg,'var':np.squeeze(var[nv-nv0]),
                                      'cmap':cm[nv-nv0],'min':fmin[ii+nv],
                                      'max':fmax[ii+nv],'llon':Lon[iLL],
                                      'ulon':Lonp,'llat':Lat[iLL],'ulat':Lat[iLL+1],
                                      'file':filefig[:-3],'title':titlefig}
                            if 'u' and 'v' in field:
                                pltarg['u'] = u
                                pltarg['v'] = v
                                pltarg['uvscale'] = uvscale[ii]
                                pltarg['step'] = step
                            PlotField.Outputs(out,pltarg)

                    nv0 = nv
                    
                    
            nfig = str(nc) + ' figure(s) done.'
            Library.Done(nfig)
