# -*- coding: utf-8 -*

class AnalysisObjects:
    
    def __init__(self, fdtd=None, mon=None):
        self.fdtd = fdtd
        # self.sim = simulation.Simulation(fdtd=fdtd, xy_span_bleed=xy_span, z_min=0, z_max=1e-06)
        self.mon = mon
        
    def farfield(self, xy_span=1e-06, z_span=1e-06, x=0, y=0, z=0, theta_max=90,
                 N_theta=180, Nphi=91):
        
        self.fdtd.unselectall()
        self.fdtd.addobject('scat_ff_radiation')
        self.fdtd.select('ff_radiation_plot') 
        self.fdtd.set('name','farfield')
        self.fdtd.set('x span', xy_span)
        self.fdtd.set('y span', xy_span)
        self.fdtd.set('z span', z_span)
        self.fdtd.set('x', x);
        self.fdtd.set('y', y);
        self.fdtd.set('z', z) # Need to work out how to find the direct centre of the overall structure

        with open('farfield_script.txt', 'r', encoding="utf8") as f:
            script = f.read().rstrip()
    
        self.fdtd.set('analysis script', script)
        self.fdtd.addanalysisresult('T')
        self.fdtd.addanalysisresult('Purcell')
        self.fdtd.addanalysisresult('P_vs_theta')
        self.fdtd.addanalysisprop("near field points per wavelength", 0, 3)
        self.fdtd.addanalysisprop("N phi per half", 0, Nphi)
        self.fdtd.addanalysisprop("project all wavelengths", 0, 1)
        self.fdtd.addanalysisprop('include_z_min_boundary', 0, 1);
        self.fdtd.set('theta max', theta_max)
        self.fdtd.set('N theta', N_theta)
        self.fdtd.unselectall()

    def Qfactor(self):
        
        self.fdtd.unselectall()
        self.fdtd.select('source')
        fstart = self.fdtd.get('frequency start')
        fstop = self.fdtd.get('frequency stop')
        
        self.fdtd.addobject('Qanalysis')
        self.fdtd.select('Qanalysis')
        self.fdtd.set('x span', 0.1e-06)
        self.fdtd.set('y span', 0.1e-06)
        self.fdtd.set('f min', fstart)
        self.fdtd.set('f max', fstop)
        self.fdtd.set('t start', 100e-15)
        self.fdtd.addanalysisresult('f_spectrum')
        self.fdtd.unselectall()

    def mode_volume_3D(self, xy_span=1e-06, z_span=1e-06, x=0, y=0, z=0):
        self.fdtd.addobject('mode_volume')
        self.fdtd.select('mode_volume')
        self.fdtd.set('x span', xy_span)
        self.fdtd.set('y span', xy_span)
        self.fdtd.set('z span', z_span)
        self.fdtd.set('x', x)
        self.fdtd.set('y', y)
        self.fdtd.set('z', z)
        self.fdtd.unselectall
        
    def mode_volume_2D(self, xy_span_pml=1e-06, apod="Start", apod_center=0, apod_start_w=100e-09,
                       z_min=0, z_max=1e-06):
        # self.fdtd.unselectall
        self.fdtd.addanalysisgroup()
        self.fdtd.set('name', 'mode volume 2D')
        
        self.mon.power_monitor(name='xy_middle', montype="2D Z-normal", plane="xy",
                              x_min=-.5*xy_span_pml, x_max=.5*xy_span_pml, 
                              y_min=-.5*xy_span_pml, y_max=.5*xy_span_pml, 
                              z=0,
                              apod="Start", apod_center=apod_center, apod_time_width=apod_start_w)        
        
        self.mon.power_monitor(name='xz_middle', montype="2D Y-normal", plane="xz",
                              x_min=-.5*xy_span_pml, x_max=.5*xy_span_pml, 
                              z_min=z_min, z_max=z_max, 
                              y=0,
                              apod="Start", apod_center=apod_center, apod_time_width=apod_start_w)    
        
        self.mon.power_monitor(name='yz_middle', montype="2D X-normal", plane="yz",
                              y_min=-.5*xy_span_pml, y_max=.5*xy_span_pml, 
                              z_min=z_min, z_max=z_max, 
                              x=0,
                              apod="Start", apod_center=apod_center, apod_time_width=apod_start_w)    
        
        self.fdtd.select('xy_middle')
        self.fdtd.addtogroup('mode volume 2D')        
        
        self.fdtd.select('xz_middle')
        self.fdtd.addtogroup('mode volume 2D')        
        
        self.fdtd.select('yz_middle')
        self.fdtd.addtogroup('mode volume 2D')   
        
        with open('modeVol_2D.txt', 'r', encoding="utf8") as f:
            script = f.read().rstrip()
    
        self.fdtd.set('analysis script', script)
        # self.fdtd.addanalysisresult('T')
        # self.fdtd.addanalysisresult('Purcell')
        # self.fdtd.addanalysisresult('P_vs_theta')
        # self.fdtd.addanalysisprop("near field points per wavelength", 0, 3)
        # self.fdtd.addanalysisprop("N phi per half", 0, Nphi)
        # self.fdtd.addanalysisprop("project all wavelengths", 0, 1)
        # self.fdtd.addanalysisprop('include_z_min_boundary', 0, 1);
        # self.fdtd.set('theta max', theta_max)
        # self.fdtd.set('N theta', N_theta)
        # self.fdtd.unselectall()
        
        self.fdtd.unselectall()
    
