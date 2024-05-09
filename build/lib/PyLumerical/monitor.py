# -*- coding: utf-8 -*-
import numpy as np

class Monitor:
    """A class for creating monitors in FDTD simulations.

    Author       : Gareth Sion jones
    Date         : May, 2024
    Version      : 0.1
    Organisation : University of Oxford

    """
    def __init__(self, fdtd=None, monitor_pt=100, monitor_span=10e-09, wlen=637e-09):
        """
        Constructor method
        
        Params:
            fdtd (object, optional): 
                Lumerical FDTD object for simulation. Default is None.
            monitor_pt (int, optional): 
                Number of frequency points to monitor. Default is 100.
            monitor_span (float, optional): 
                Wavelength span to monitor in meters. Default is 10e-09.
            wlen (float, optional): 
                Center wavelength to monitor in meters. Default is 637e-09.
    
        Returns:
    
        Raises:
            IOError: 
                If no Lumerical FDTD object is defined in the constructor.
        """
        
        self.fdtd = fdtd
        if not self.fdtd:
            raise IOError ("A lumerical FDTD object must be defined in the constructor")
            
        self.monitor_pt = monitor_pt
        self.monitor_span = monitor_span
        self.wlen = wlen
        
        self.fdtd.setglobalmonitor('use source limits',0);                     #The user will specify
        self.fdtd.setglobalmonitor('frequency points', self.monitor_pt);
        self.fdtd.setglobalmonitor('wavelength span', self.monitor_span);
        self.fdtd.setglobalmonitor('wavelength center', self.wlen);
    
    def power_monitor(self, name='power_monitor', montype="2D Z-normal", plane="xy",
                      x_min=0, x_max=1e-06, y_min=0, y_max=1e-06, z_min=0, z_max=1e-06, 
                      x=1e-06, y=1e-06, z=1e-06,
                      apod="Start", apod_center=0e-15, apod_time_width=0e-15):
        """
            Add a power monitor to the simulation.

            This method adds a power monitor to the simulation with specified parameters.

            Params:
                name (str, optional): 
                    Name of the power monitor. Default is 'power_monitor'.
                montype (str, optional): 
                    Type of the power monitor. Default is "2D Z-normal".
                plane (str, optional): 
                    Plane of the power monitor ('xy', 'xz', or 'yz'). Default is "xy".
                x_min (float, optional): 
                    Minimum x-coordinate of the monitor region. Default is 0.
                x_max (float, optional): 
                    Maximum x-coordinate of the monitor region. Default is 1e-06.
                y_min (float, optional): 
                    Minimum y-coordinate of the monitor region. Default is 0.
                y_max (float, optional): 
                    Maximum y-coordinate of the monitor region. Default is 1e-06.
                z_min (float, optional):
                    Minimum z-coordinate of the monitor region. Default is 0.
                z_max (float, optional): 
                    Maximum z-coordinate of the monitor region. Default is 1e-06.
                x (float, optional):
                    x-coordinate of the monitor point. Default is 1e-06.
                y (float, optional): 
                    y-coordinate of the monitor point. Default is 1e-06.
                z (float, optional): 
                    z-coordinate of the monitor point. Default is 1e-06.
                apod (str, optional): 
                    Apodization type for the monitor. Default is "Start".
                apod_center (float, optional): 
                    Apodization center in seconds. Default is 0e-15.
                apod_time_width (float, optional): 
                    Apodization time width in seconds. Default is 0e-15.
                    
            Returns:
                
        """
        self.fdtd.addpower()
        self.fdtd.set('name', name)
        self.fdtd.set('monitor type', montype)
        
        if plane=="xy":
            self.fdtd.set("x max", x_max)
            self.fdtd.set("x min", x_min)
            self.fdtd.set("y max", y_max)
            self.fdtd.set("y min", y_min)
            self.fdtd.set('z',z)
        
        if plane=="xz":
            self.fdtd.set("x max", x_max)
            self.fdtd.set("x min", x_min)
            self.fdtd.set("y", y)
            self.fdtd.set("z min", z_min)
            self.fdtd.set('z max', z_max)
        
        if plane=="yz":
            self.fdtd.set("x", x)
            self.fdtd.set("y min", y_min)
            self.fdtd.set("y max", y_max)
            self.fdtd.set("z min", z_min)
            self.fdtd.set('z max', z_max)
        
        self.fdtd.set('Output Px',1)
        self.fdtd.set('Output Py',1)
        self.fdtd.set('Output Pz',1)
        self.fdtd.set('Output Hx',1)
        self.fdtd.set('Output Hy',1)
        self.fdtd.set('Output Hz',1)
        self.fdtd.set('apodization',apod)
        self.fdtd.set('apodization center', apod_center)
        self.fdtd.set('apodization time width', apod_time_width)
        
    def index_monitor(self, name="n", monitor_type="2D Y-normal",
                       x_min=0, x_max=1e-06, z_min=0, z_max=1e-06):
        """
            Add an index monitor to the simulation.
    
            This method adds an index monitor to the simulation with specified parameters.
    
            Params:
                name (str, optional): 
                    Name of the index monitor. Default is "n".
                monitor_type (str, optional): 
                    Type of the index monitor. Default is "2D Y-normal".
                x_min (float, optional): 
                    Minimum x-coordinate of the monitor region. Default is 0.
                x_max (float, optional): 
                    Maximum x-coordinate of the monitor region. Default is 1e-06.
                z_min (float, optional): 
                    Minimum z-coordinate of the monitor region. Default is 0.
                z_max (float, optional): 
                    Maximum z-coordinate of the monitor region. Default is 1e-06.
            
            Returns:

        """
        self.fdtd.addindex()
        self.fdtd.set('monitor type', monitor_type)
        self.fdtd.set('name', name)
        self.fdtd.set('y',0)
        self.fdtd.set("x min", x_min)
        self.fdtd.set("x max", x_max)
        self.fdtd.set("z min", z_min)
        self.fdtd.set("z max", z_max)
        self.fdtd.set("spatial interpolation", 'NEAREST MESH CELL')
        
    def Q_monitor(self, Qmonitor_zspan=0, Qmonitor_zlayer=0, t_sample=10, dipole_shift=0):
        """
            Add Q monitors to the simulation.
    
            This method adds time monitors to the simulation for extracting information 
            about quality factor
    
            Params:
                Qmonitor_zspan (float, optional): 
                    Half span in the z-direction for the Q monitors in nanometers. Default is 0.
                Qmonitor_zlayer (int, optional): 
                    Number of layers in the z-direction for the Q monitors. Default is 0.
                t_sample (int, optional): 
                    Minimum sampling per cycle for the monitors. Default is 10.
                dipole_shift (float, optional): 
                    Shift along the z-axis of the dipole position in meters. Default is 0.
            
            Returns:

        """
        # specify the half span in 3D
        x_span = 0.1e-6;
        y_span = 0.1e-6;
        z_span = Qmonitor_zspan*1e-9;
        
        # error checking on the inputs
        nx = 2
        ny = 2
        nz = Qmonitor_zlayer
        
        if ((x_span==0) or (nx < 1)):
            nx = 1 
            
        if ((y_span==0) or (ny < 1)):
            ny = 1 
            
        if ((z_span==0) or (nz < 1)):
            nz = 1 
        
        if nx == 1:
            x_span = 0
            
        if ny == 1:
            y_span = 0
            
        if nz == 1:
            z_span = 0
        
        # define position vectors for monitors
        xpos = np.linspace(0.03*1e-6, 0.03*1e-6+x_span, nx);
        ypos = np.linspace(0.03*1e-6, 0.03*1e-6+y_span, ny);
        zpos = np.linspace(0, z_span, nz);
        
        # add the monitors
        self.fdtd.addanalysisgroup()
        self.fdtd.set('name', 'Q monitors')
        self.fdtd.set('x',0)
        self.fdtd.set('y',0)
        self.fdtd.set('z',0)
        self.fdtd.set('use relative coordinates',0)
        self.fdtd.groupscope("Q monitors")
        
        mon_counter = 0;
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    mon_counter = mon_counter + 1
                    self.fdtd.addtime()
                    self.fdtd.set("name","t_h_"+ str(mon_counter))
                    self.fdtd.set("x",xpos[i])
                    self.fdtd.set("y",ypos[j])
                    self.fdtd.set("z",dipole_shift-zpos[k])
                    self.fdtd.set("min sampling per cycle", t_sample)

        self.fdtd.groupscope("::model")