# -*- coding: utf-8 -*-

class Simulation:
    """A class representing a simulation environment for FDTD simulations.
    
    Author       : Gareth Sion jones
    Date         : May, 2024
    Version      : 0.1
    Organisation : University of Oxford
    """
    
    def __init__(self, fdtd=None, xy_span_bleed=1e-06, runtime=1000e-15, meshacc=3, z_min=0, z_max=1e-06):
        """
        Constructor method

        Params:
            fdtd (object, optional): 
                Lumerical FDTD object for simulation. Default is None.
            xy_span_bleed (float, optional): 
                Span of the simulation region in the xy plane with bleed. Default is 1e-06.
            runtime (float, optional): 
                Total runtime of the simulation in seconds. Default is 1000e-15.
            meshacc (int, optional): 
                Mesh accuracy for simulation. Default is 3.
            z_min (float, optional): 
                Minimum z-coordinate of the simulation region. Default is 0.
            z_max (float, optional): 
                Maximum z-coordinate of the simulation region. Default is 1e-06.

        Returns:

        Raises:
            IOError: 
                If no Lumerical FDTD object is defined in the constructor.

        """
        self.xy_span_bleed = xy_span_bleed
        self.runtime = runtime
        self.meshacc = meshacc
        self.x_min = -.5*self.xy_span_bleed;
        self.x_max = .5*self.xy_span_bleed;
        self.y_min = -.5*self.xy_span_bleed;
        self.y_max = .5*self.xy_span_bleed;
        self.z_min = z_min
        self.z_max = z_max
        
        self.fdtd = fdtd
        if not self.fdtd:
            raise IOError ("A lumerical FDTD object must be defined in the constructor")
    
    def fdtd_region(self, x_min_bc="PML", y_min_bc="PML", z_min_bc="PML", 
                    dt_stab=0.99, fdtd_layers=8, min_layers=8, max_layers=64,
                    autoshutoff=1e-05):
        """Configure FDTD simulation region parameters.

        This method sets up the FDTD simulation region with specified boundary conditions,
        simulation time, mesh accuracy, and PML settings.

        Params:
            x_min_bc (str, optional): 
                Boundary condition at the minimum x-coordinate. Default is "PML".
            y_min_bc (str, optional):
                Boundary condition at the minimum y-coordinate. Default is "PML".
            z_min_bc (str, optional): 
                Boundary condition at the minimum z-coordinate. Default is "PML".
            dt_stab (float, optional): 
                Stability factor for time step. Default is 0.99.
            fdtd_layers (int, optional): 
                Number of PML layers. Default is 8.
            min_layers (int, optional): 
                Minimum number of PML layers. Default is 8.
            max_layers (int, optional): 
                Maximum number of PML layers. Default is 64.
            autoshutoff (float, optional): 
                Automatic shutoff threshold. Default is 1e-05.
        
        Returns:

        """
        self.fdtd.addfdtd()
        self.fdtd.set("x max",self.x_max);
        self.fdtd.set("x min",self.x_min);
        self.fdtd.set("y max",self.y_max);
        self.fdtd.set("y min",self.y_min);
        self.fdtd.set("z min",self.z_min);
        self.fdtd.set("z max",self.z_max); 
        self.fdtd.set("x min bc",x_min_bc);
        self.fdtd.set("y min bc",y_min_bc); 
        self.fdtd.set("z min bc",z_min_bc); 
        self.fdtd.set('simulation time',self.runtime);
        self.fdtd.set('mesh accuracy',self.meshacc);
        self.fdtd.set('dt stability factor',dt_stab);
        self.fdtd.set('pml layers',fdtd_layers);
        self.fdtd.set('pml min layers',min_layers);
        self.fdtd.set('pml max layers',max_layers);
        self.fdtd.set("auto shutoff min",autoshutoff);
        
    def add_mesh(self, name='mesh', x=0, y=0, zmin=0, zmax=1e-06,
                 xspan=1e-06, yspan=1e-06, dx=0.1e-06, dy=0.1e-06, dz=0.1e-06,
                 based_on_struct=False, struct=""):
        """Add an overriding mesh region to the simulation.

        Params:
            name (str, optional): 
                Name of the mesh. Default is 'mesh'.
            x (float, optional):
                x-coordinate of the mesh origin. Default is 0.
            y (float, optional): 
                y-coordinate of the mesh origin. Default is 0.
            zmin (float, optional): 
                Minimum z-coordinate of the mesh region. Default is 0.
            zmax (float, optional): 
                Maximum z-coordinate of the mesh region. Default is 1e-06.
            xspan (float, optional): 
                x-span of the mesh region. Default is 1e-06.
            yspan (float, optional): 
                y-span of the mesh region. Default is 1e-06.
            dx (float, optional): 
                Spatial step size in the x-direction. Default is 0.1e-06.
            dy (float, optional): 
                Spatial step size in the y-direction. Default is 0.1e-06.
            dz (float, optional): 
                Spatial step size in the z-direction. Default is 0.1e-06.
            based_on_struct (bool, optional):
                Choose whether to base the mesh geometry on an existing structure
            struct (string, optional):
                structure on which to base the mesh, if based_on_struct is True
                
        Returns:
            

        """
        self.fdtd.addmesh()
        self.fdtd.set("name",name)
        self.fdtd.set("override x mesh", 1)
        self.fdtd.set("override y mesh", 1)
        self.fdtd.set("override z mesh", 1)
        self.fdtd.set("dx", dx)
        self.fdtd.set("dy", dy)
        self.fdtd.set("dz", dz)
        
        if based_on_struct:
            self.fdtd.set("based on a structure", 1);
            self.fdtd.set("structure", struct);
        else:
            self.fdtd.set("x",x)
            self.fdtd.set("x span",xspan)
            self.fdtd.set("y",y)
            self.fdtd.set("y span",yspan)
            self.fdtd.set("z min",zmin)
            self.fdtd.set("z max",zmax)
    