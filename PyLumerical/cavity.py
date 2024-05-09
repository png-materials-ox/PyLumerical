# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

class Cavity:
    """ Define the geometry of a Fabry-Perot optical microcavity
    
    Author       : Gareth Sion jones
    Date         : May, 2024
    Version      : 0.1
    Organisation : University of Oxford
    """
    
    def __init__(self, roc=2, wlen=637e-09, q=2, ncav=2.1, rough=False, fdtd=None):
        """
        
        Constructor method

        Params:
            roc (float, optional): 
                Radius of curvature of the cavity mirror. Default is 2.
            wlen (float, optional): 
                Wavelength of light in meters. Default is 637e-09.
            q (int, optional): 
                Longitudinal mode index of the cavity. Default is 2.
            ncav (float, optional): 
                Refractive index of the cavity material. Default is 2.1.
            rough (bool, optional): 
                Flag indicating whether the cavity mirrors are rough. Default is False.
            fdtd (object, optional):
                Lumerical FDTD object for simulation. Default is None.
        
        Returns:

        Raises:
            IOError: If no Lumerical FDTD object is defined in the constructor.
        """
        self.roc = roc
        self.wlen = wlen
        self.q = q
        self.ncav = ncav
        self.rough = rough
        self.wlen_effective = self.wlen/self.ncav
        
        self.L_depth = .5 * self.wlen 
        self.feature_width = 2*np.sqrt(self.roc**2-(self.roc-self.L_depth)**2)
        self.xy_span = 1.5*self.feature_width
        self.xy_span_bleed = self.xy_span + 300e-09
        
        self.fdtd = fdtd
        
        if not self.fdtd:
            raise IOError ("A lumerical FDTD object must be defined in the constructor")
            
    def cavity_length(self):
        """Calculate the length of the optical cavity.

        Params:
    
        Returns:
            float: 
                Length of the optical cavity.

        """
        return .5*(self.q * self.wlen_effective) # Optical cavity length
            
    def beam_waist(self):
        """Calculate the cavity beam waist.

        Params:
    
        Returns:
            float: 
                beam waist.

        """
        
        L = self.cavity_length()
        w0 = np.sqrt(self.wlen_effective*L/ np.pi*np.sqrt(self.roc/L-1))
        return w0
    
    def rayleigh_range(self):
        """Calculate the cavity rayleigh range.

        Params:
    
        Returns:
            float: 
                rayleigh range.

        """
        w0 = self.beam_waist()
        return np.pi*w0**2/self.wlen_effective
    
    def beam_waist_at_feature(self):
        """Calculate the cavity beam waist at the feature.

        Params:
    
        Returns:
            float: 
                beam waist at feature.

        """
        L = self.cavity_length()
        w0 = self.beam_waist() 
        rayleigh_range = self.rayleigh_range()
        wz = w0*np.sqrt(1+(L / rayleigh_range())**2)
        return wz
    
    def guoy_shift(self):
        """Calculate the cavity Guoy shift.

        Params:
    
        Returns:
            float: 
                Guoy shift.

        """
        return np.arctan(self.cavity_length() / self.rayleigh_range())    
            
    def convex_feature(self, resolution=512, plotting=False):
        """Calculate the convex feature of the cavity surface and optionally plot it.

        Params:
            resolution (int, optional): 
                Resolution of the grid. Default is 512.
            plotting (bool, optional): 
                Flag indicating whether to plot the feature. Default is False.

        Returns:
            dict: Dictionary containing grid data and mirror separation.

        """
       
        feature_halfwidth = self.feature_width/2;
        L = self.cavity_length()
           
        # Create grid for plotting
        x = np.linspace(-self.xy_span/2, self.xy_span/2, resolution) # x and y both range from negative to positive.
        y = x
        X,Y = np.meshgrid(x, y)
        Z = np.ones((resolution, resolution), dtype=float)

        # Calculate surface
        for a in range(0, len(x)):
            for b in range(0, len(y)):
                R_buff= np.sqrt(x[a]**2 + y[b]**2)
                if R_buff < feature_halfwidth:
                    Z[a,b]= -np.sqrt(self.roc**2 - R_buff**2) + np.sqrt(self.roc**2 - feature_halfwidth**2)
                else:
                    Z[a,b]= 0
        if plotting:
            plt.pcolor(X,Y,Z)
            plt.show()
            
        depth = Z.max() - Z.min()
        mir_sep = L-depth

        return {'X':X, 'Y':Y, 'Z':Z, 'x':x, 'y':y, 'mir_sep':mir_sep}  

        
    def planar_dbr(self, num=10, n1=2.21, n2=1.42, nsub=1.45, Lsub=27e-09):
        """Generate a planar distributed Bragg reflector (DBR) structure.

        Params:
            num (int, optional): 
                Number of DBR layers. Default is 10.
            n1 (float, optional): 
                Refractive index of the first material. Default is 2.21.
            n2 (float, optional): 
                Refractive index of the second material. Default is 1.42.
            nsub (float, optional): 
                Refractive index of the substrate. Default is 1.45.
            Lsub (float, optional): 
                Thickness of the substrate in meters. Default is 27e-09.

        Returns:

        """
        thickness_1 = self.wlen/(4*n1);                      
        thickness_2 = self.wlen/(4*n2);                      
                   
        dbr_thickness = thickness_1 + thickness_2

        #Generate planar upper mirror
        self.fdtd.addrect()
        self.fdtd.set("x max",self.xy_span_bleed/2)
        self.fdtd.set("x min",-self.xy_span_bleed/2)
        self.fdtd.set("y max",self.xy_span_bleed/2)
        self.fdtd.set("y min",-self.xy_span_bleed/2)
        self.fdtd.set("index",n2)
        self.fdtd.set("z min",0)
        self.fdtd.set("z max",thickness_2)
        self.fdtd.set("name","hi layer template")
        self.fdtd.set('render type','wireframe')

        self.fdtd.addrect()
        self.fdtd.set("x max",self.xy_span_bleed/2)
        self.fdtd.set("x min",-self.xy_span_bleed/2)
        self.fdtd.set("y max",self.xy_span_bleed/2)
        self.fdtd.set("y min",-self.xy_span_bleed/2)
        self.fdtd.set("index",n1)
        self.fdtd.set("z min",thickness_2)
        self.fdtd.set("z max",dbr_thickness)
        self.fdtd.set("name","lo layer template")
        self.fdtd.set('render type','wireframe')

        #DBR code stacks object called layer_profile
        for i in range(1, num):
            self.fdtd.select("hi layer template")
            self.fdtd.copy(0,0,(i-1)*dbr_thickness)
            self.fdtd.set("name", "hi layer")
            self.fdtd.addtogroup('planar mirror')
            self.fdtd.select("lo layer template")
            self.fdtd.copy(0,0,(i-1)*dbr_thickness)
            self.fdtd.set("name", "lo layer")
            self.fdtd.addtogroup('planar mirror')

        self.fdtd.select("hi layer template");
        self.fdtd.delete()
        self.fdtd.select("lo layer template");
        self.fdtd.delete()

        # Create substrate extending to Region_max
        #groupscope('::model');
        self.fdtd.addrect()
        self.fdtd.set("x max",self.xy_span_bleed/2)
        self.fdtd.set("x min",-self.xy_span_bleed/2)
        self.fdtd.set("y max",self.xy_span_bleed/2)
        self.fdtd.set("y min",-self.xy_span_bleed/2)
        self.fdtd.set("index",nsub)
        self.fdtd.set("z min",(num-1)*dbr_thickness)
        self.fdtd.set("z max",((num-1)*dbr_thickness) + Lsub) 
        self.fdtd.set("name","substrate")
        self.fdtd.set('render type','wireframe')
        self.fdtd.addtogroup('planar mirror')
    
    def convex_dbr(self, num=10, n1=2.21, n2=1.42, nsub=1.45, Lsub=27e-09):
        """Generate a convex distributed Bragg reflector (DBR) structure.

        Params:
            num (int, optional): 
                Number of DBR layers. Default is 10.
            n1 (float, optional): 
                Refractive index of the first material. Default is 2.21.
            n2 (float, optional): 
                Refractive index of the second material. Default is 1.42.
            nsub (float, optional): 
                Refractive index of the substrate. Default is 1.45.
            Lsub (float, optional): 
                Thickness of the substrate in meters. Default is 27e-09.

        Returns:

        """
        
        thickness_1 = self.wlen/(4*n1);                      
        thickness_2 = self.wlen/(4*n2);                      
                   
        dbr_thickness = thickness_1 + thickness_2
        
        feat = self.convex_feature()
        # Create upper surface with thickness depth 
        self.fdtd.addimport()
        self.fdtd.importsurface2(feat['Z'], feat['x'], feat['y'], 1) 
        self.fdtd.set("name","hi layer template")
        self.fdtd.set('render type','wireframe')
        self.fdtd.set("index",n2)
        self.fdtd.set("z max", -feat['mir_sep'])
        self.fdtd.set("z min", -self.cavity_length() - thickness_2)
        
        self.fdtd.addimport()
        self.fdtd.importsurface2(feat['Z'], feat['x'], feat['y'], 1) 
        self.fdtd.set("name","lo layer template")
        self.fdtd.set('render type','wireframe')
        self.fdtd.set("index",n1)
        self.fdtd.set("z max", -feat['mir_sep'] - thickness_2)
        self.fdtd.set("z min", -self.cavity_length() - dbr_thickness)
        
        for i in range(num):
            self.fdtd.select("hi layer template")
            self.fdtd.copy(0,0,(i-1)*(-dbr_thickness))
            self.fdtd.set("index",n2)
            self.fdtd.set("name","hi layer")
            self.fdtd.addtogroup('featured mirror')
            self.fdtd.select("lo layer template")
            self.fdtd.copy(0,0,(i-1)*(-dbr_thickness))
            self.fdtd.set("index",n1)
            self.fdtd.set("name","low layer")
            self.fdtd.addtogroup('featured mirror')

        
        self.fdtd.select("hi layer template")
        self.fdtd.delete()
        self.fdtd.select("lo layer template")
        self.fdtd.delete()
        
        # Create substrate extending to Region_min
        self.fdtd.addimport()
        self.fdtd.importsurface2(feat['Z'], feat['x'], feat['y'], 1)
        self.fdtd.set("name","substrate")
        self.fdtd.set('render type','wireframe')
        self.fdtd.set("index", nsub)
        self.fdtd.set("z max", -feat['mir_sep'] - ((num-1)*dbr_thickness))
        self.fdtd.set("z min", -self.cavity_length() - ((num-1)*dbr_thickness) - Lsub)
        self.fdtd.addtogroup('featured mirror')

    
    def cavity_medium(self):
        """Create a dielectric medium representing the cavity.

        This method adds a rectangular dielectric medium object to the simulation domain,
        representing the optical cavity.
        
        Params:
            
        Returns:

        """
        self.fdtd.addrect()
        self.fdtd.set('name','dielectric medium')
        self.fdtd.set('index',self.ncav)
        self.fdtd.set("x max",self.xy_span_bleed/2)
        self.fdtd.set("x min",-self.xy_span_bleed/2)
        self.fdtd.set("y max",self.xy_span_bleed/2)
        self.fdtd.set("y min",-self.xy_span_bleed/2)
        self.fdtd.set('z max',0)
        self.fdtd.set('z min',-self.cavity_length()) 
        
    def build_cavity(self, num_planar=10, num_feat=10, n1=2.21, n2=42, nsub=1.45, Lsub=27e-09, resolution=512):
        """Build the complete optical cavity structure.

        This method constructs the optical cavity by adding planar DBR layers, the cavity medium, and convex DBR layers.

        Params:
            num_planar (int, optional): 
                Number of planar DBR layers. Default is 10.
            num_feat (int, optional): 
                Number of convex DBR layers. Default is 10.
            n1 (float, optional): 
                Refractive index of the first material. Default is 2.21.
            n2 (float, optional): 
                Refractive index of the second material. Default is 42.
            nsub (float, optional): 
                Refractive index of the substrate. Default is 1.45.
            Lsub (float, optional): 
                Thickness of the substrate in meters. Default is 27e-09.
            resolution (int, optional): 
                Resolution of the grid. Default is 512.
            
        Returns:

        """
        self.planar_dbr(num=num_planar, n1=n1, n2=n2, nsub=nsub, Lsub=Lsub)
        self.cavity_medium()
        self.convex_dbr(num=num_feat, n1=n1, n2=n2, nsub=nsub, Lsub=Lsub)    
        
        self.fdtd.addstructuregroup()
        self.fdtd.set('name','structure');
        self.fdtd.select('planar mirror');
        self.fdtd.addtogroup('structure');
        self.fdtd.select('featured mirror');
        self.fdtd.addtogroup('structure');
        self.fdtd.select('dielectric medium');
        self.fdtd.set('override mesh order from material database',1);
        self.fdtd.set('mesh order',3);
        self.fdtd.addtogroup('structure');
        
        self.fdtd.select('structure::featured mirror::hi layer');
        
        self.fdtd.select('structure');
        self.fdtd.set('first axis','y');
        self.fdtd.set('rotation 1',180);
        self.fdtd.set('x',0);
        self.fdtd.set('y',0);
        self.fdtd.set('z',0);
        
        self.fdtd.select('structure::planar mirror');   
        self.fdtd.set('x',0);
        self.fdtd.set('y',0);
        self.fdtd.set('z',0);
        
        self.fdtd.select('structure::featured mirror');
        self.fdtd.set('x',0);
        self.fdtd.set('y',0);
        self.fdtd.set('z',0);
        