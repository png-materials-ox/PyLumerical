# -*- coding: utf-8 -*-
"""
Created on Wed May 15 12:21:44 2024

@author: mans4209
"""

class AnalysisObjects:
    
    def __init__(self, fdtd=None):
        self.fdtd = fdtd
        
    def farfield(self, xy_span=1e-06, z_span=1e-06, x=0, y=0, z=0, theta_max=90,
                 N_theta=180, Nphi=91):
        
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

    def Qfactor(self):
        
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

    def mode_volume_3D(self, xy_span=1e-06, z_span=1e-06, x=0, y=0, z=0):
        self.fdtd.addobject('mode_volume')
        self.fdtd.select('mode_volume')
        self.fdtd.set('x span', xy_span)
        self.fdtd.set('y span', xy_span)
        self.fdtd.set('z span', z_span)
        self.fdtd.set('x', x)
        self.fdtd.set('y', y)
        self.fdtd.set('z', z)