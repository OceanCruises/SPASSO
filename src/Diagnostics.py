"""
Lagrangian and Eulerian diagnostics computation.

Three classes are defined:
    - Lagrangian: includes code to compute particle trajectories and derived
    diagnostics
    - Eulerian: includes code to compute Eulerian diagnostics from velocity 
    field
    - ParticleSet: includes code to initialize numerical particle before
    advection

@author: lrousselet
"""

from scipy.interpolate import interp2d,interpn,dfitpack
import numpy as np
import datetime as dt
import warnings
import sys
import glob
# Spasso spec
import GlobalVars, Fields, PlotField, Library


class Lagrangian():
    
    def diag(self,diag=None,method=None,f=None,**kwargs):
        """ Initialize and launch the diagnostics requested and method set to advect the Lagrangian particles previsouly set by ParticleSet
        
        :param diag: list of the requested Lagrangian diagnostics. 
        diag = ['FTLE'] or ['LLADV']. 
        For multiple diagnostics: diag = ['FTLE','LLADV']
        
        :param method: method for particle advection. Default is set to Rnga-kutta 4 'rk4flat'.
        
        :param f: function to get particle new position
        
        :output: outputs are concatenated in list 'out' starting with particle
        trajectories as a first dict. Then outputs are in the same order as
        listed in diag.
        """
        if method != None:
            if 'Library' in sys.modules.keys(): Library.tic()
            out = []    
            if 'numstep' in kwargs:
                Nstep = kwargs['numstep']
            else:
                Nstep = 4 #default value
                warntxt = "Warning: 'numstep' is not defined -> using default value (4)"
                warnings.warn(warntxt)
                if 'Library' in sys.modules.keys(): Library.Logfile(warntxt)
            self.Nstep = Nstep
            
            if method == 'rk1flat':
                trjf = self.rk1flat(f,Nstep,**kwargs)
            elif method == 'rk4flat':
                trjf = self.rk4flat(f,Nstep,**kwargs)
            else:
                warntxt = "Warning: 'method' is not defined -> using default value (rk4flat)"
                warnings.warn(warntxt)
                if 'Library' in sys.modules.keys(): Library.Logfile(warntxt)
                trjf = self.rk4flat(f,Nstep,**kwargs)
            out.append(trjf)
            if 'Library' in sys.modules.keys(): Library.toc('Lagrangian trajectories')
            
        if diag != None:
            for i in diag:
                if i == 'LLADV': dd = self.LLADV(trjf,**kwargs)
                if i == 'SSTADV': dd = self.SSTADV(trjf,**kwargs)
                if i == 'FTLE': dd = self.FTLE(trjf,**kwargs)
                if i == 'OWTRAJ': dd = self.OWTRAJ(trjf,**kwargs)
                if i == 'TIMEFROMBATHY': dd = self.TIMEFROMBATHY(trjf,**kwargs)
                out.append(dd)
                
        return out
    
    def backonsphere(self,x,y):
        """
        Returns
        -------
        New particle x,y positions when they reach domain boundary.
        Used when Periodic boundary is set to True and coordinates set to
        spherical.
        """
        x = np.array(x).astype('float64')
        y = np.array(y).astype('float64')
        
        #zonal
        pos90p = y > 90
        pos90m = y < -90
        if pos90p.any:
            y[pos90p] = 180 - y[pos90p]
            x[pos90p] = x[pos90p] + 180
        if pos90m.any:
            y[pos90m] = -180 - y[pos90p]
            x[pos90m] = x[pos90m] + 180
        
        #meridional
        lon = self.lon
        xm = x < lon[0]
        xp = x > lon[-1]
        if xm.any:
            x[xm] += lon[-1] - lon[0]
        if xp.any:
            x[xp] -= lon[-1] - lon[0]       
        return x,y

    def interpf(self,t,x,y,**kwargs):
        """
        Velocity field interpolation at particle position using interpn function.

        Parameters
        ----------
        t : time
        x : particle x position
        y : particle y position
        **kwargs : coordinates

        Returns
        -------
        u_in : u interpolated value at particle position.
        v_in : v interpolated value at particle position.

        """
        udim, vdim = np.asarray(self.u_nonan), np.asarray(self.v_nonan)
        if np.size(t)!=np.size(x):
            t = np.tile(t,len(x))
        if udim.ndim == 2 and vdim.ndim == 2:
            if 'coordinates' in kwargs and kwargs['coordinates']=='spherical':
                xn,yn = Lagrangian.backonsphere(self,x,y)
                new_grid = list(zip(xn,yn))
            else:
                new_grid = list(zip(x,y))
            if self.lon.ndim == 2:
                lon,lat = self.lon[:,0],self.lat[0,:]
            else:
                lon,lat = self.lon,self.lat
            u_in = interpn((lon,lat),self.u_nonan,new_grid,bounds_error=False,fill_value=np.nan)
            v_in = interpn((lon,lat),self.v_nonan,new_grid,bounds_error=False,fill_value=np.nan)
        elif udim.ndim == 3 and udim.ndim == 3:
            if 'coordinates' in kwargs and kwargs['coordinates']=='spherical':
                xn,yn = Lagrangian.backonsphere(self,x,y)
                new_grid = list(zip(t,xn,yn))
            else:
                new_grid = list(zip(t,x,y))
            u_in = interpn((self.dates,self.lon,self.lat),self.u_nonan,new_grid,bounds_error=False,fill_value=np.nan)
            v_in = interpn((self.dates,self.lon,self.lat),self.v_nonan,new_grid,bounds_error=False,fill_value=np.nan)
        return u_in,v_in
    
    def interp2d_pairs(*args,**kwargs):
        """ Same interface as interp2d but the returned interpolant will evaluate 
        its inputs as pairs of values.
        """
        # Internal function, that evaluates pairs of values, output has the same shape as input
        def interpolant(x,y,f):
            x,y = np.asarray(x), np.asarray(y)
            return (dfitpack.bispeu(f.tck[0], f.tck[1], f.tck[2], f.tck[3], f.tck[4], x.ravel(), y.ravel())[0]).reshape(x.shape)
        # Wrapping the scipy interp2 function to call out interpolant instead
        return lambda x,y: interpolant(x,y,interp2d(*args,**kwargs))

    def rk1flatstep(self,t,x,y,f,h): 
        xp,yp=f(self,t,x,y) #with mathematical formalism, this is d(pts)/dt, or pts', that is, ptsp
        x_n=x+xp*h
        y_n=y+yp*h
        return x_n,y_n

    def rk4flatstep(self,t,x,y,f,h,**kwargs): 
        k1 = h * np.asarray(f(self,t,x,y,**kwargs))
        k2 = h * np.asarray(f(self,(t+h/2),(x+k1[0]/2),(y+k1[1]/2),**kwargs))
        k3 = h * np.asarray(f(self,(t+h/2),(x+k2[0]/2),(y+k2[1]/2),**kwargs))
        k4 = h * np.asarray(f(self,(t+h), (x+k3[0]), (y+k3[1]),**kwargs))
        k = (k1+2*k2+2*k3+k4)/6
        x_n = x + k[0]
        y_n = y + k[1]
        return x_n,y_n

    def rk1flat(self,f,Nstep,**kwargs):
        """
        First order Runga-Kutta method

        Parameters
        ----------
        f : interpolation function (default is interpf)
        Nstep : number of step

        Returns
        -------
        trjf : dictionnary including particle trajectories: trjx, trjy
                                     particle initial positions: lons, lats
                                     particle final positions: lonf, latf
        """
        t_v = self.pt
        x = self.px
        y = self.py
        trjx = []
        trjy = []
        trjx.append(x)
        trjy.append(y)
        t=t_v[0]
        h=(t_v[1]-t_v[0])/(Nstep*self.numdays)
        xn=x
        yn=y
        for i in range(Nstep*self.numdays):
            t=t+h
            xn,yn = Lagrangian.rk1flatstep(self,t,xn,yn,f,h)
            if 'coordinates' in kwargs and kwargs['coordinates']=='spherical':
                xn,yn = Lagrangian.backonsphere(self,xn,yn)
            trjx.append(xn)
            trjy.append(yn)
        trjf = {'lons':self.lons,'lats':self.lats,'trjx':trjx,'trjy':trjy,'lonf':trjx[-1:],'latf':trjy[-1:]}
        return trjf

    def rk4flat(self,f,Nstep,**kwargs):
        """
        Fourth order Runga-Kutta method

        Parameters
        ----------
        f : interpolation function (default is interpf)
        Nstep : number of step

        Returns
        -------
        trjf : dictionnary including particle trajectories: trjx, trjy, trjt
                                     particle initial positions: lons, lats
                                     particle final positions: lonf, latf
        """
        t_v = self.pt
        x = self.px
        y = self.py
        trjx,trjy,trjt = ([] for i in range(3))
#            if (key == 'noise'):
#                noise=value    #Not yet implemented;     
        if np.size(np.shape(t_v))==1:
            h = (t_v[1]-t_v[0])/(Nstep*self.numdays)
        else:
            h = (t_v[1,0]-t_v[0,0])/(Nstep*self.numdays)
        xn = x
        yn = y
        trjx.append(x)
        trjy.append(y)
        trjt.append(t_v[0])
        t = t_v[0]-h
        for i in range(Nstep*self.numdays):
            t = t+h
            xn,yn = Lagrangian.rk4flatstep(self,t,xn,yn,f,h,**kwargs)
            if 'coordinates' in kwargs and kwargs['coordinates']=='spherical':
                xn,yn = Lagrangian.backonsphere(self,xn,yn)
            trjx.append(xn)
            trjy.append(yn)
            trjt.append(t)
        
        trjf = {'lons':self.lons,'lats':self.lats,'trjx':trjx,'trjy':trjy,'trjt':trjt,'lonf':trjx[-1:],'latf':trjy[-1:]}
        return trjf

    def LLADV(self,trjf,**kwargs):
        """ Compute Lon/Lat advections
        :param trj: dictionnary with particle trajectories from advection 
        
        :output lladv: lons/lats are longitudes and latitudes for mapping; 
        lonf_map and latf_map are longitude and latitude advections respectively 
        formatted for mapping.
        """
        if 'Library' in sys.modules.keys():
            Library.tic()
            
        lonf = trjf['lonf']
        latf = trjf['latf']
        lons = trjf['lons']
        lats = trjf['lats']
        if 'dayv' in kwargs:
            dayv = kwargs['dayv']
        else:
            print("Missing 'dayv' argument", file=sys.stderr)
        
        if 'PeriodicBC' in kwargs and kwargs['PeriodicBC']==True:
            if self.loni[1] < self.loni[0]:
                lons0 = lons
                lons0[lons0<0] += 360
                lonf[0][lonf[0]<0]  += 360
        else:
            lons0 = lons                

        [Xs,Ys]=np.meshgrid(lons0,lats)
        lonf_map=np.reshape(lonf,(np.shape(Xs)[0],np.shape(Xs)[1]))
        latf_map=np.reshape(latf,(np.shape(Xs)[0],np.shape(Xs)[1]))
        lonf_map = Xs - lonf_map
        latf_map = Ys - latf_map
        
        if 'PeriodicBC' in kwargs and kwargs['PeriodicBC']==True:
            if self.loni[1] < self.loni[0]:
                lonf_map[lonf_map>180] -= 360
        
        lladv = {'lons':lons,'lats':lats,'lonf_map':lonf_map,'latf_map':latf_map}
        
        ### saving ###
        if 'output' in kwargs and kwargs['output']=='netcdf':
            if 'GlobalVars' in sys.modules.keys():
                prod = kwargs['product']
                date = dt.datetime.strftime(dt.datetime.strptime(dayv,'%Y-%m-%d'),'%Y%m%d')
                fname = GlobalVars.Dir['dir_wrk']+date+'_LLADV_'+prod+'.nc'
                title = 'LON/LAT ADVECTION '+date+' (computed with '+prod+')'
                Fields.LLADV(fname).createnc(lons,lats,lonf_map,title,vvar2=latf_map)
            else:
                warnings.warn("Warning: Use Save.py to save your data")
             
        return lladv
    
    def SSTADV(self,trjf,**kwargs):
        """ Compute SST advection from lon/lat advection
        :param lladv: lon/lat advection (returned from 'LLADV')
        
        :output sstadv: lons/lats are longitudes and latitudes for mapping; lonf_map and latf_map are longitude and latitude advections respectively formatted for mapping.
        """
        trjx = trjf['trjx']
        trjy = trjf['trjy']
        lons = trjf['lons']
        lats = trjf['lats']
        
        if 'dayv' in kwargs:
            dayv = kwargs['dayv']
        else:
            print("Missing 'dayv' argument", file=sys.stderr)
        
        if 'PeriodicBC' in kwargs and kwargs['PeriodicBC']==True:
            if self.loni[1] < self.loni[0]:
                lons0 = lons
                lons0[lons0<0] += 360
        else:
            lons0 = lons                

        #compute lon/lat adv at (t0 - n)
        if 'daysst' in kwargs:
            day = kwargs['daysst']
        else:
            if 'GlobalVars' in sys.modules.keys():
                day = GlobalVars.Lag['sstadvd']
            else:
                day = 3
                warnings.warn("Warning: 'daysst' is not defined -> using default value (3)")

        [Xs,Ys]=np.meshgrid(lons0,lats)
        lonf = np.asarray(trjx)
        latf = np.asarray(trjy)
        lonf,latf = lonf[day,:],latf[day,:]
        lonf_map=np.reshape(lonf,(np.shape(Xs)[0],np.shape(Xs)[1]))
        latf_map=np.reshape(latf,(np.shape(Xs)[0],np.shape(Xs)[1]))
        lonf_map = lonf_map.flatten()
        latf_map = latf_map.flatten()
        
        if 'PeriodicBC' in kwargs and kwargs['PeriodicBC']==True:
            if self.loni[1] < self.loni[0]:
                lonf_map[lonf_map>180] -= 360

        #load sst map at (t0 - n)
        if 'sstfield' in kwargs:
            field = kwargs['sstfield']
            lon = field['lon']
            lat = field['lat']
            var = field['sst']
        elif 'GlobalVars' in sys.modules.keys():
            nprod = GlobalVars.config.get('products',GlobalVars.Lag['sstprod']+'prod')
            data = GlobalVars.config.get('products',GlobalVars.Lag['sstprod']+'_data')
            var = Library.GetVars(data)
            tmpd = var['date'][0]
            fname = glob.glob(GlobalVars.Dir['dir_wrk']+'/'+tmpd+'_'+nprod+'.nc')
            if fname:
                field = eval('Fields.'+nprod+'(fname[0]).loadnc()')
                lon = field['lon']
                lat = field['lat']
                var = field['var']
            else:
                lon,lat,var = [],[],[]
                warntxt = "Warning: No SST file, sst advection is empty."
                warnings.warn(warntxt)
                Library.Logfile(warntxt)
        else:
            warnings.warn('Missing SST field.')

        if isinstance(var,list):
            sst_map = np.zeros((np.shape(Xs)[0],np.shape(Xs)[1]))
        else:
            # Create the interpolant (same interface as interp2d)
            f = Lagrangian.interp2d_pairs(lon,lat,var,kind='linear')
            # Evaluate the interpolant on each pairs of x and y values
            sst = f(lonf_map,latf_map)
            sst_map = np.reshape(sst,(np.shape(Xs)[0],np.shape(Xs)[1]))
         
        sstadv = {'lons':lons,'lats':lats,'sstadv':sst_map}
        
        ### saving ###
        if 'output' in kwargs and kwargs['output']=='netcdf':
            if 'GlobalVars' in sys.modules.keys():
                prod = kwargs['product']
                date = dt.datetime.strftime(dt.datetime.strptime(dayv,'%Y-%m-%d'),'%Y%m%d')
                fname = GlobalVars.Dir['dir_wrk']+date+'_SSTADV_'+prod+'.nc'
                title = prod + 'SST ADVECTION '+date
                Fields.SSTADV(fname).createnc(lons,lats,sst_map,title)
            else:
                warnings.warn("Warning: Use Save.py to save your data")
        return sstadv
    
    def FTLE(self,trjf,**kwargs):
        """
        Compute Finite Time Lyapunov Exponents from initial and final positions
        of computed particle trajectories.

        Parameters
        ----------
        trjf : dictionnary including particle trajectories, initial and final
        positions.

        Returns
        -------
        ftle : dictionnary with lon/lat grid (lons/lats) and gridded values of FTLE.

        """
        lons = trjf['lons']
        lats = trjf['lats']
        lonf = trjf['lonf']
        latf = trjf['latf']
        if 'numdays' in kwargs:
            numdays = int(kwargs['numdays'])
        else:
            print("Missing 'numdays' argument to compute FTLE", file=sys.stderr)
        if 'dayv' in kwargs:
            dayv = kwargs['dayv']
        else:
            print("Missing 'dayv' argument", file=sys.stderr)
        
        if 'PeriodicBC' in kwargs and kwargs['PeriodicBC']==True:
            if self.loni[1] < self.loni[0]:
                lons0 = lons
                lons0[lons0<0] += 360
                lonf[0][lonf[0]<0]  += 360
        else:
            lons0 = lons
                
        [Xs,Ys]=np.meshgrid(lons0,lats)
        lon0_map=Xs
        lat0_map=Ys
        lonf_map=np.reshape(lonf,(np.shape(Xs)[0],np.shape(Xs)[1]))
        latf_map=np.reshape(latf,(np.shape(Xs)[0],np.shape(Xs)[1]))
        ### Gradients of final positions ###
        [d_lonyf,d_lonxf]=np.gradient(lonf_map)
        [d_latyf,d_latxf]=np.gradient(latf_map)
        d_lonxf=d_lonxf*np.cos(Ys/180*np.pi)
        ### Gradients of initial positions ###
        [d_lony0,d_lonx0]=np.gradient(lon0_map)
        [d_laty0,d_latx0]=np.gradient(lat0_map)
        d_lonx0=d_lonx0*np.cos(Ys/180*np.pi)
        ### Final separation ###
        Xgradf=(d_lonxf**2)+(d_latxf**2)
        Ygradf=(d_lonyf**2)+(d_latyf**2)
        XYgradf=[]
        XYgradf.append(Xgradf)
        XYgradf.append(Ygradf)
        final_separation=np.max(XYgradf,0)
        ### Initial separation ###
        initial_separation=(d_lonx0**2)+(d_latx0**2)
        ### FTLE ###
        ftle_lyap=(np.log(final_separation/initial_separation)/(numdays))/2
        ftle = {'lons':lons,'lats':lats,'ftle':ftle_lyap}
        
        ### saving ###
        if 'output' in kwargs and kwargs['output']=='netcdf':
            if 'GlobalVars' in sys.modules.keys():
                prod = kwargs['product']
                date = dt.datetime.strftime(dt.datetime.strptime(dayv,'%Y-%m-%d'),'%Y%m%d')
                fname = GlobalVars.Dir['dir_wrk']+date+'_FTLE_'+prod+'.nc'
                title = prod + 'FTLE '+date
                Fields.FTLE(fname).createnc(lons,lats,ftle_lyap,title)
            else:
                warnings.warn("Warning: Use Save.py to save your data")
        return ftle
    
    def OWTRAJ(self,trjf,**kwargs):
        """
        Compute Okubo-Weiss parameter along each particle trajectories.
        Identify consecutive negative values of Okubo-Weiss and compute a gridded
        mean time (i.e. retention parameter) of consecutive negative Okubo-Weiss 
        values over the domain.

        Parameters
        ----------
        trjf : dictionnary including particle trajectories, initial and final
        positions.
        
        Returns
        -------
        owdisp : dictionnary including lon/lat grid (lons/lats) and gridded value
        of the retention parameter.

        """
        if 'Library' in sys.modules.keys():
            Library.tic() 
        
        t = np.ravel(np.array(trjf['trjt'])[::self.Nstep,:])
        tv = np.array(trjf['trjt'])[::self.Nstep,0]
        x = np.ravel(np.array(trjf['trjx'])[::self.Nstep,:])
        y = np.ravel(np.array(trjf['trjy'])[::self.Nstep,:])
        sz = np.shape(np.array(trjf['trjx'])[::self.Nstep,:])
        lons = trjf['lons']
        lats = trjf['lats']

        
        if 'dayv' in kwargs:
            dayv = kwargs['dayv']
        else:
            print("Missing 'dayv' argument (format: %Y-%m-%d) to compute OWTRAJ", file=sys.stderr)
            
        if 'ds' in kwargs:
            ds = kwargs['ds']
        else:
            print("Missing 'ds' argument to compute OWTRAJ", file=sys.stderr)
        
        if 'PeriodicBC' in kwargs and kwargs['PeriodicBC']==True:
            if self.loni[1] < self.loni[0]:
                lons0 = lons
                lons0[lons0<0] += 360
                x[x<0] += 360
        else:
            lons0 = lons
            
        [Xs,Ys]=np.meshgrid(lons0,lats)
        RT=6371e5
        dUdx,dUdy,dVdx,dVdy,U,V = ([] for i in range(6))
            
        #convert t in days since beginning of integration
        dayv0 = dt.datetime.toordinal(dt.datetime.strptime(dayv,'%Y-%m-%d').date())
        tv -= dayv0
        convx = np.pi/180.*(RT*np.cos(y/180*np.pi))
        convy = 1/180*np.pi*RT
        velx = Lagrangian.interpf(self,t,x+ds,y,**kwargs)
        vely = Lagrangian.interpf(self,t,x,y+ds,**kwargs)
        velmx = Lagrangian.interpf(self,t,x-ds,y,**kwargs)
        velmy = Lagrangian.interpf(self,t,x,y-ds,**kwargs)
        dsx = 2*ds*convx
        dsy = 2*ds*convy
        dUdx = (velx[0]-velmx[0])/dsx
        dUdy = (vely[0]-velmy[0])/dsy
        dVdx = (velx[1]-velmx[1])/dsx
        dVdy = (vely[1]-velmy[1])/dsy     
        
        sn = dUdx-dVdy
        ss = dUdy+dVdx
        vor = -dUdy+dVdx
        ow = (sn**2)+(ss**2)-(vor**2);
        owm = np.reshape(ow,sz)

        posexit = []
        for ct in range(sz[1]):
            tmp = np.sign(owm[:,ct])
            maxv = np.nanmax(tmp)
            exit1 = np.argmax(tmp)
            exit = tv[exit1]
            if (maxv<0): exit = tv[-1]
            posexit.append(exit)  
        
        owd=np.reshape(posexit,(np.shape(Xs)[0],np.shape(Xs)[1]))
        owdisp = {'lons':lons,'lats':lats,'owdisp':owd}
        
        ### saving ###
        if 'output' in kwargs and kwargs['output']=='netcdf':
            if 'GlobalVars' in sys.modules.keys():
                prod = kwargs['product']
                date = dt.datetime.strftime(dt.datetime.strptime(dayv,'%Y-%m-%d'),'%Y%m%d')
                fname = GlobalVars.Dir['dir_wrk']+date+'_OWTRAJ_'+prod+'.nc'
                title = prod + 'OWTRAJ '+date
                Fields.OWTRAJ(fname).createnc(lons,lats,owd,title)
            else:
                warnings.warn("Warning: Use Save.py to save your data")
                
        if 'Library' in sys.modules.keys():
            Library.toc('OWTRAJ')   
        return owdisp
    
    def TIMEFROMBATHY(self,trjf,**kwargs):
        """
        Computed gridded mean time since particle last contact over a specified
        bathymetry value.

        Parameters
        ----------
        trjf : dictionnary including particle trajectories, initial and final
        positions.

        Returns
        -------
        timfbathy : dictionnary including lon/lat griid, gridded time from 
        bathymetry (timfb), gridded intial longitude (lonfb)/latitude (latfb) 
        of the particle when they were over the specified bathymetry value.

        """
        lons = trjf['lons']
        lats = trjf['lats']
        trjx = np.asarray(trjf['trjx'])
        trjy = np.asarray(trjf['trjy'])
                
        if 'dayv' in kwargs:
            dayv = kwargs['dayv']
        else:
            print("Missing 'dayv' argument (format: %Y-%m-%d) to compute OWTRAJ", file=sys.stderr)
             
        if 'PeriodicBC' in kwargs and kwargs['PeriodicBC']==True:
            if self.loni[1] < self.loni[0]:
                trjx[trjx<0] += 360
                
        #extract grid from bathy file in subdomain
        if 'bathyfield' in kwargs:
            field = kwargs['bathyfield']
            lon = field['lon']
            lat = field['lat']
            z = field['z']
            trjx = trjx[1::kwargs['numstep']]
            trjy = trjy[1::kwargs['numstep']]
            bathylvl = kwargs['bathylvl']
        elif 'GlobalVars' in sys.modules.keys():
            trjx = trjx[1::GlobalVars.Lag['numstep']]
            trjy = trjy[1::GlobalVars.Lag['numstep']]
            bathylvl = GlobalVars.Lag['bathylvl']
            file = GlobalVars.Dir['dir_bathy']+GlobalVars.Lag['bathyfile']
            field = Fields.ETOPO.loadnc(file,rlon=GlobalVars.Lag['loni'],rlat=GlobalVars.Lag['lati'])
            lon = field['lon']
            lat = field['lat']
            z = field['z']
        else:
            print("Missing bathymetry field: bathyfield = field", file=sys.stderr)

        # Create the interpolant (same interface as interp2d)
        f = Lagrangian.interp2d_pairs(lon,lat,z,kind='linear')
        # Evaluate the interpolant on each pairs of x and y values
        trjd = f(trjx,trjy)
        
        touched = []
        nottouched = []
        for ct in range(0,np.shape(trjx)[1]):
            touch = [i for i,v in enumerate(trjd[:,ct]) if v > bathylvl]
            if not touch: 
                touch=0 #not touched
                nottouched.append(ct)
            else: 
                touch = min(touch)
            touched.append(touch)
        
        touchedlat = []
        touchedlon = []

        for ct in range(0,len(touched)):
            touchedlat.append(trjy[touched[ct],ct])
            touchedlon.append(trjx[touched[ct],ct])

        touched,touchedlat,touchedlon = np.asarray(touched,dtype=float),np.asarray(touchedlat,dtype=float),np.asarray(touchedlon,dtype=float)
        touched[nottouched] = np.nan
        touchedlat[nottouched] = np.nan
        touchedlon[nottouched] = np.nan
        
        touched = touched.reshape((len(lats),len(lons)))
        touchedlat = touchedlat.reshape((len(lats),len(lons)))
        touchedlon = touchedlon.reshape((len(lats),len(lons)))
        
        timfbathy = {'lons':lons,'lats':lats,'timfb':touched,'latfb':touchedlat,'lonfb':touchedlon}
        
        ### saving ###
        if 'output' in kwargs and kwargs['output']=='netcdf':
            if 'GlobalVars' in sys.modules.keys():
                prod = kwargs['product']
                date = dt.datetime.strftime(dt.datetime.strptime(dayv,'%Y-%m-%d'),'%Y%m%d')
                fname = GlobalVars.Dir['dir_wrk']+date+'_TIMEFROMBATHY_'+prod+'.nc'
                title = prod + 'Time from Bathy '+date
                Fields.TIMEFROMBATHY(fname).createnc(lons,lats,touched,title,vvar2=touchedlat,vvar3=touchedlon)
            else:
                warnings.warn("Warning: Use Save.py to save your data")
        return timfbathy
    

class ParticleSet(Lagrangian):
    """
    Set particle initialization parameters.
        - from_input: particles are initialized on the position (px, py) and 
        time defined by the user
        - from_grid: particles initialized on a grid defined by longitude min,
        longitude max et delta step.
    
    Default parameter to compute Lagrangian diagnostics from satellite velocity
    field is from_grid.
    """
    def __init__(self,pt=None,px=None,py=None,lons=None,lats=None,numdays=None,loni=None,lati=None,delta0=None,dayv=None,**kwargs):
        self.pt = pt
        self.px = px
        self.py = py
        self.lons = lons
        self.lats = lats
        self.loni = loni
        self.lati = lati
        if 'fieldset' in kwargs:
            ff = kwargs.get('fieldset')
            self.lon = ff['lon']
            self.lat = ff['lat']
            self.u = ff['u']
            self.v = ff['v']
            self.dates = ff['dates']
            self.u_nonan=np.where(np.isnan(self.u),0,self.u)
            self.v_nonan=np.where(np.isnan(self.v),0,self.v)
            ParticleSet.check_dimensions(self)

        if 'PeriodicBC' in kwargs and kwargs['PeriodicBC']==True:
            dlon = np.array(np.diff(self.lon))
            self.lon = np.hstack([self.lon[0]-dlon[0]/2,self.lon,self.lon[-1]+dlon[-1]/2])    
            if self.u.ndim == 2 and self.u.ndim == 2:
                self.u = np.vstack([self.u[0,:],self.u,self.u[-1,:]])
                self.v = np.vstack([self.v[0,:],self.v,self.v[-1,:]])
                self.u_nonan = np.vstack([self.u_nonan[0,:],self.u_nonan,self.u_nonan[-1,:]])
                self.v_nonan = np.vstack([self.v_nonan[0,:],self.v_nonan,self.v_nonan[-1,:]])
            elif self.u.ndim == 3 and self.v.ndim == 3:
                self.u = [np.vstack([self.u[i,0,:],self.u[i,:,:],self.u[i,-1,:]]) for i in range(self.u.shape[0])]
                self.v = [np.vstack([self.v[i,0,:],self.v[i,:,:],self.v[i,-1,:]]) for i in range(self.v.shape[0])]
                self.u_nonan = [np.vstack([self.u_nonan[i,0,:],self.u_nonan[i,:,:],self.u_nonan[i,-1,:]]) for i in range(self.u_nonan.shape[0])]
                self.v_nonan = [np.vstack([self.v_nonan[i,0,:],self.v_nonan[i,:,:],self.v_nonan[i,-1,:]]) for i in range(self.v_nonan.shape[0])]

        if numdays is None:
            self.numdays = 1
        else:
            self.numdays = numdays
            
        return
    
    @classmethod
    def from_input(cls,pt,px,py,**kwargs):
        px = px
        py = py
        pt = pt
        lons,lats = px,py
        return cls(pt,px,py,lons,lats,**kwargs)

    @classmethod
    def from_grid(cls,numdays,loni,lati,delta0,dayv,**kwargs):
        if 'PeriodicBC' in kwargs and kwargs['PeriodicBC']==True:
            if loni[0] > loni[1]:
                bnd = ParticleSet.get_boundary(**kwargs)
                xi = loni[1] + (bnd[1]-bnd[0])
                px = np.arange(loni[0],xi,delta0)
                py = np.arange(lati[0],lati[1],delta0)
        elif 'PeriodicBC' in kwargs and kwargs['PeriodicBC']==False:
            print(kwargs['PeriodicBC'])
            if(((loni[0])==(loni[1]))&((lati[0])!=(lati[1]))):
                px = loni[0]
                py = np.arange(lati[0],lati[1],delta0)
            elif(((loni[0])!=(loni[1]))&((lati[0])==(lati[1]))):
                py = lati[0]
                px = np.arange(loni[0],loni[1],delta0)
            elif(((loni[0])==(loni[1]))&((lati[0])==(lati[1]))):
                px = loni[0]
                py = lati[0]
            else:
                px = np.arange(loni[0],loni[1],delta0)
                py = np.arange(lati[0],lati[1],delta0)
        else:
            if(((loni[0])==(loni[1]))&((lati[0])!=(lati[1]))):
                px = loni[0]
                py = np.arange(lati[0],lati[1],delta0)
            elif(((loni[0])!=(loni[1]))&((lati[0])==(lati[1]))):
                py = lati[0]
                px = np.arange(loni[0],loni[1],delta0)
            elif(((loni[0])==(loni[1]))&((lati[0])==(lati[1]))):
                px = loni[0]
                py = lati[0]
            else:
                px = np.arange(loni[0],loni[1],delta0)
                py = np.arange(lati[0],lati[1],delta0)

        lons,lats = px,py
        [X,Y] = np.meshgrid(px,py)
        px = np.ravel(X)
        py = np.ravel(Y)
        if dayv == None:
            pt0 = np.array([numdays] * len(px))
            ptf = np.array([0] * len(px))
        else:
            day2 = dt.datetime.strptime(dayv,'%Y-%m-%d').date()
            day2j = dt.datetime.toordinal(day2)
            if 'mode' in kwargs:
                if kwargs['mode'] == 'backward':
                    day1j = day2j-numdays
                elif kwargs['mode'] == 'forward':
                    day1j = day2j+numdays
            else:
                #default is backward
                day1j = day2j-numdays       
            pt0 = np.array([day2j] * len(px))
            ptf = np.array([day1j] * len(px))                  
        pt = np.array([pt0,ptf])

        return cls(pt,px,py,lons,lats,numdays,loni,lati,**kwargs)

    def get_boundary(**kwargs):
        try:
            if 'xmin' and 'xmax' in kwargs:
                xmin, xmax = kwargs['xmin'], kwargs['xmax']
            elif 'ymin' and 'ymax' in kwargs:
                ymin, ymax = kwargs['ymin'], kwargs['ymax']
            elif 'fieldset' in kwargs:
                ff = kwargs.get('fieldset')
                xmin, xmax = ff['lon'][0], ff['lon'][-1]
                ymin, ymax = ff['lat'][0], ff['lat'][-1]
        except Exception as e:
            print(e)
            
        return [xmin,xmax,ymin,ymax]
    
    def check_dimensions(self):
        #lon,lat,dates should be one dimensional
        if self.lon.ndim > 1:
            if self.lon.ndim == 2 and self.lon.ndim == self.lat.ndim:
                self.lon,self.lat = self.lon[:,0],self.lat[0,:]
                if np.all(np.diff(self.lon))==0 and np.all(np.diff(self.lat))==0:
                    self.lon,self.lat = self.lon[0,:],self.lat[:,0]
            elif self.lon.ndim == 3:
                self.lon,self.lat = self.lon[0,:,0],self.lat[0,0,:]
                if np.all(np.diff(self.lon))==0 and np.all(np.diff(self.lat))==0:
                    self.lon,self.lat = self.lon[:,0,0],self.lat[0,:,0]
                    if np.all(np.diff(self.lon))==0 and np.all(np.diff(self.lat))==0:
                        self.lon,self.lat = self.lon[0,0,:],self.lat[:,0,0]
        return

class Eulerian():
    """
    Compute Eulerian diagnostic from velocity field input.
    """
    def __init__(self,fieldset=None,dayv=None):
        self.RT = 6371e5
        if fieldset!=None:
            self.lon = fieldset['lon']
            self.lat = fieldset['lat']
            if 'u' in fieldset or 'v' in fieldset:
                self.u = fieldset['u']
                self.v = fieldset['v']
                if isinstance(self.u, np.ma.MaskedArray):
                    self.u = self.u.filled(np.nan)
                    self.v = self.v.filled(np.nan)
            elif 'var' in fieldset:
                self.var = fieldset['var']
                if hasattr(self,'u'):
                    delattr(self,'u')
                    delattr(self,'v')
            if 'dates' in fieldset:
                self.dates = fieldset['dates']
        else:
            print("Missing 'fieldset' argument. Cannot compute Eulerian diags.", file=sys.stderr)
            
        if dayv!=None:
            self.dayv = dayv
        else:
            self.dayv = self.dates[-1]
            print("Missing 'dayv' argument to compute Eulerian diag. Default value is used (i.e. last date of field)", file=sys.stderr)
        return
    
    def diag(self,diag=None,**kwargs):
        out = []
        dout = []
        if diag!=None:
            for i in diag:
                if i == 'KE' and hasattr(self,'u'): dd = self.KE(**kwargs), dout.append(i)
                if i == 'OW' and hasattr(self,'u'): dd = self.OW(**kwargs), dout.append(i)
                if i == 'dSST' and hasattr(self,'var'): dd = self.dSST(**kwargs), dout.append(i)
                if 'dd' in dir(): 
                    out.append(dd)
                    del dd
        return out, dout
    
    def KE(self,**kwargs):
        if self.u.ndim == 2:
            U = self.u
            V = self.v
        elif self.u.ndim == 3:
                day = dt.datetime.toordinal(dt.datetime.strptime(self.dayv,'%Y-%m-%d').date())
                idd = np.where(self.dates == day)
                U = np.squeeze(self.u[idd,:,:])
                V = np.squeeze(self.v[idd,:,:])
        if 'UVunit' in kwargs:
            if kwargs['UVunit']=='m/s':
                Ucms,Vcms = U*1e2,V*1e2
            if kwargs['UVunit']=='cm/s':
                Ucms,Vcms = U,V

        if 'delta' and 'lon' and 'lat' in kwargs:
            delta0 = kwargs['delta']
            loni = kwargs['lon']
            lati = kwargs['lat']
            lon = np.arange(loni[0],loni[1],delta0/2)
            lat = np.arange(lati[0],lati[1],delta0/2)
            [X,Y] = np.meshgrid(lon,lat)
            lon0 = np.array(self.lon)
            lat0 = np.array(self.lat)
            u_nonan = np.where(np.isnan(Ucms),0,Ucms)
            v_nonan = np.where(np.isnan(Vcms),0,Vcms)

            if np.shape(u_nonan) == (lat0.size,lon0.size):
                fu = interp2d(lon0,lat0,u_nonan, kind='cubic')
                fv = interp2d(lon0,lat0,v_nonan, kind='cubic')
            elif np.shape(u_nonan) == (lon0.size,lat0.size):
                fu = interp2d(lon0,lat0,u_nonan.T, kind='cubic')
                fv = interp2d(lon0,lat0,v_nonan.T, kind='cubic')
            Ucms = fu(lon, lat)
            Vcms = fv(lon, lat)
        else:
            print("Missing 'delta', 'lon' and 'lat' arguments to compute Eulerian diag.", file=sys.stderr)

        E = (Ucms**2) + (Vcms**2)
        
        ### saving ###
        if 'output' in kwargs and kwargs['output']=='netcdf':
            if 'GlobalVars' in sys.modules.keys():
                prod = kwargs['product']
                date = dt.datetime.strftime(dt.datetime.strptime(self.dayv,'%Y-%m-%d'),'%Y%m%d')
                fname = GlobalVars.Dir['dir_wrk']+date+'_KE_'+prod+'.nc'
                title = prod + 'KE '+date
                Fields.KE(fname).createnc(lon,lat,E,title)
            else:
                warnings.warn("Warning: Use Save.py to save your data")
                
        KE = {'lon':X,'lat':Y,'KE':E}
        return KE
    
    def OW(self,**kwargs):
        if self.u.ndim == 2:
            U = self.u
            V = self.v
        elif self.u.ndim == 3:
                day = dt.datetime.toordinal(dt.datetime.strptime(self.dayv,'%Y-%m-%d').date())
                idd = np.where(self.dates == day)
                U = np.squeeze(self.u[idd,:,:])
                V = np.squeeze(self.v[idd,:,:])
        [X,Y]=np.meshgrid(self.lon,self.lat)
        if np.shape(X)!=np.shape(U): [Y,X]=np.meshgrid(self.lat,self.lon)
        if 'UVunit' in kwargs:
            if kwargs['UVunit']=='m/s':
                Ucms,Vcms = U*1e2,V*1e2
            if kwargs['UVunit']=='cm/s':
                Ucms,Vcms = U,V
        if 'delta' and 'lon' and 'lat' in kwargs:
            delta0 = kwargs['delta']
            loni = kwargs['lon']
            lati = kwargs['lat']
            lon = np.arange(loni[0],loni[1],delta0/2)
            lat = np.arange(lati[0],lati[1],delta0/2)
            [X,Y] = np.meshgrid(lon,lat)
            lon0 = np.array(self.lon)
            lat0 = np.array(self.lat)
            u_nonan = np.where(np.isnan(Ucms),0,Ucms)
            v_nonan = np.where(np.isnan(Vcms),0,Vcms)
            if np.shape(u_nonan) == (lat0.size,lon0.size):
                fu = interp2d(lon0,lat0,u_nonan, kind='cubic')
                fv = interp2d(lon0,lat0,v_nonan, kind='cubic')
            elif np.shape(u_nonan) == (lon0.size,lat0.size):
                fu = interp2d(lon0,lat0,u_nonan.T, kind='cubic')
                fv = interp2d(lon0,lat0,v_nonan.T, kind='cubic')
            Ucms = fu(lon, lat)
            Vcms = fv(lon, lat)
        else:
            print("Missing 'delta', 'lon' and 'lat' arguments to compute Eulerian diag.", file=sys.stderr)
        
        dUdy,dUdx = np.gradient(Ucms)
        dVdy,dVdx = np.gradient(Vcms)
        
        tmp,Dx = np.gradient(X/180*np.pi*self.RT*np.cos(Y*np.pi/180))
        Dy,tmp = np.gradient(Y/180*np.pi*self.RT)
        
        dUdx = dUdx/Dx
        dVdx = dVdx/Dx
        dUdy = dUdy/Dy
        dVdy = dVdy/Dy 
        
        sn = dUdx - dVdy
        ss = dUdy + dVdx
        vor = -dUdy + dVdx
        ow = sn**2 + ss**2 - vor**2
        ow = ow*(60*60*24)**2 #daily
        
        ### saving ###
        if 'output' in kwargs and kwargs['output']=='netcdf':
            if 'GlobalVars' in sys.modules.keys():
                prod = kwargs['product']
                date = dt.datetime.strftime(dt.datetime.strptime(self.dayv,'%Y-%m-%d'),'%Y%m%d')
                fname = GlobalVars.Dir['dir_wrk']+date+'_OW_'+prod+'.nc'
                title = prod + 'OW '+date
                Fields.OW(fname).createnc(lon,lat,ow,title)
            else:
                warnings.warn("Warning: Use Save.py to save your data")

        OW = {'lon':X,'lat':Y,'sn':sn,'ss':ss,'vor':vor,'ow':ow}
        return OW
    
    def dSST(self,**kwargs):
        dSSTdy,dSSTdx = np.gradient(self.var)
        
        SSTgrad = np.sqrt(dSSTdx**2 + dSSTdy**2)
        SSTgradir = np.arctan2(dSSTdy,dSSTdx)
        
        ### saving ###
        if 'output' in kwargs and kwargs['output']=='netcdf':
            if 'GlobalVars' in sys.modules.keys():
                prod = kwargs['product']
                date = dt.datetime.strftime(dt.datetime.strptime(self.dayv,'%Y-%m-%d'),'%Y%m%d')
                fname = GlobalVars.Dir['dir_wrk']+date+'_dSST_'+prod+'.nc'
                title = prod + 'dSST '+date
                Fields.dSST(fname).createnc(self.lon,self.lat,SSTgrad,title,vvar2=SSTgradir)
            else:
                warnings.warn("Warning: Use Save.py to save your data")

        dSST = {'lon':self.lon,'lat':self.lat,'SSTgrad':SSTgrad,'SSTgradir':SSTgradir}
        return dSST

#Spasso spec   
def Launch(cruise,approach):
    opt = [str(x) for x in GlobalVars.config.get('plot_options','options').split(',')]
    if approach == 'eulerian':
        for pr in GlobalVars.Eul['products']:
            if pr:
                nprod = GlobalVars.config.get('products',pr+'prod')
                var = Library.GetVars(GlobalVars.config.get('products',pr+'_data'))
                if GlobalVars.Eul['dayv']=='default':
                    dayv = var['datef']
                    tmp = var['date']
                else:
                    dayv = [str(x) for x in GlobalVars.Eul['dayv'].split(',')]
                    tmp = []
                    for nf in dayv:
                        tmp.append(var['date'][var['datef'].index(nf)])
                for nf in range(len(dayv)):
                    Library.printMessage("Computing Eulerian for "+nprod+" "+dayv[nf])
                    path = GlobalVars.Dir['dir_wrk']+'/'+tmp[nf]+'*'+nprod+'*.nc'
                    exf,ff = Library.ExistingFile(path,tmp[nf])
                    if exf == True:
                        fname = glob.glob(path)[0]
                    elif exf == False:
                        #download data
                        eval("Fields."+nprod+".download(date=tmp[nf],cp='yes')")
                        fname = glob.glob(path)[0]
                    field = eval('Fields.'+nprod+'(fname).loadnc()')
                    out,dout = Eulerian(field,dayv[nf]).diag(diag=GlobalVars.Eul['diag'],
                                                              UVunit=GlobalVars.Eul['UVunit'],delta=GlobalVars.Eul['delta0'],
                                                              lon=GlobalVars.Eul['loni'],lat=GlobalVars.Eul['lati'],
                                                              output='netcdf',product=nprod)
                    #plot field
                    if out:
                        for di in range(0,len(out)):
                            Library.printMessage("Ploting "+dout[di]+" for "+nprod)
                            PlotField.PlotField.Plot(cruise,dout[di],opt,type='Eulerian')
                    else:
                        Library.Done('None')
    elif approach == 'lagrangian':
        for pr in GlobalVars.Lag['products']:
            if pr:
                nprod = GlobalVars.config.get('products',pr+'prod')
                var = Library.GetVars(GlobalVars.config.get('products',pr+'_data'))
                if GlobalVars.Lag['dayv']=='default':
                    dayv = var['datef']
                    tmp = var['date']
                else:
                    dayv = [str(x) for x in GlobalVars.Lag['dayv'].split(',')]
                    tmp = []
                    for nf in dayv:
                        tmp.append(var['date'][var['datef'].index(nf)])
                for nf in range(len(dayv)):
                    Library.printMessage("Computing lagrangian for "+nprod+" "+dayv[nf])
                    
                    path = GlobalVars.Dir['dir_wrk']+'/'+tmp[nf]+'*'+nprod+'*.nc'
                    exf,ff = Library.ExistingFile(path,tmp[nf])
                    if exf == True:
                        fname = glob.glob(path)[0]
                    elif exf == False:
                        #download data
                        eval("Fields."+nprod+".download(date=tmp[nf],cp='yes')")
                        fname = glob.glob(path)[0]
                    field = eval('Fields.'+nprod+'(fname,dayv=dayv[nf]).LoadLag(GlobalVars.Lag["numdays"],product=nprod)')
                    pset = ParticleSet.from_grid(GlobalVars.Lag['numdays'],GlobalVars.Lag['loni'],
                                              GlobalVars.Lag['lati'],GlobalVars.Lag['delta0'],
                                              dayv[nf],fieldset=field,mode=GlobalVars.Lag['mode'],
                                              PeriodicBC=GlobalVars.Lag['PeriodicBC'])
                    out = pset.diag(diag=GlobalVars.Lag['diag'],method=GlobalVars.Lag['method'],
                                    f=Lagrangian.interpf,numstep=GlobalVars.Lag['numstep'],
                                    coordinates='spherical',numdays=GlobalVars.Lag['numdays'],
                                    dayv=dayv[nf],ds=1/6,output='netcdf',product=nprod,PeriodicBC=GlobalVars.Lag['PeriodicBC'])
        
                    #plot field
                    if out:
                        for di in range(0,len(out)-1):
                            Library.printMessage("Ploting "+GlobalVars.Lag['diag'][di]+" for "+nprod)
                            PlotField.PlotField.Plot(cruise,GlobalVars.Lag['diag'][di],opt,type='Lagrangian')                  
                    else:
                        Library.Done('None')
    return
