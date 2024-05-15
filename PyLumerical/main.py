# -*- coding: utf-8 -*-
"""
Created on Wed May  8 14:35:52 2024

@author: mans4209
"""

import numpy as np 
import imp
# import os
import cavity
import simulation
import source 
import monitor
import analysisobjects

# os.add_dll_directory("C:\\Program Files\\Lumerical\\v232\\api\\python\\")
lumapi = imp.load_source("lumapi","C:\\Program Files\\Lumerical\\v232\\api\\python\\lumapi.py")

save = False
filename = "testrun.fsp"

fdtd = lumapi.FDTD()

wlen = .637e-06
roc = 2e-06
q = 4
ncav = 1

# BUILD CAVITY
##############################################################################
cav = cavity.Cavity(roc=roc, wlen=wlen, q=q, ncav=ncav, rough=False, fdtd=fdtd)
cav.build_cavity(num_planar=10, num_feat=10, n1=2.21, n2=1.42, nsub=1.45, Lsub=20e-09, resolution=512)

# ADD FDTD REGION 
##############################################################################
fdtd.select('structure::planar mirror::substrate');
Region_min = -fdtd.get('z max');
fdtd.select('structure::featured mirror::substrate');
Region_max = abs(fdtd.get('z min') - .5*wlen);

sim = simulation.Simulation(fdtd=fdtd, xy_span_bleed=cav.xy_span_bleed, 
                            runtime=3000e-15, meshacc=3, z_min=Region_min, 
                            z_max=Region_max)

sim.fdtd_region(x_min_bc="Anti-Symmetric", y_min_bc="Symmetric", z_min_bc="PML", 
                    dt_stab=0.99, fdtd_layers=8, min_layers=8, max_layers=64,
                    autoshutoff=1e-05)

# ADD MESH OVERRIDE
sim.add_mesh(name='planar mesh', dx=0.01e-06, dy=0.01e-06, dz=0.01e-06, 
             based_on_struct=True, struct="planar mirror")

sim.add_mesh(name='feat mesh', dx=0.01e-06, dy=0.01e-06, dz=0.01e-06, 
             based_on_struct=True, struct="featured mirror")

# ADD SOURCE
##############################################################################
src = source.Source(fdtd=fdtd, wlen=wlen)
src.dipole(theta=90, shift=0, emission_width=100e-09)

# ADD MONITORS
##############################################################################
pml_thickness = 0e-09
xy_span_pml = cav.xy_span + pml_thickness

mon = monitor.Monitor(fdtd=fdtd)

apod_center = 100e-15
apod_start_w = 10e-15

# mon.Q_monitor(Qmonitor_zspan=10, Qmonitor_zlayer=1, t_sample=10, dipole_shift=0)

mon.index_monitor(name="n", monitor_type="2D Y-normal",
                        x_min=-.5*xy_span_pml, x_max=.5*xy_span_pml, 
                        z_min=Region_min-.5*pml_thickness, z_max=Region_max + .5*pml_thickness)

mon.power_monitor(name='xy_exoplanar', montype="2D Z-normal", plane="xy",
                      x_min=-.5*xy_span_pml, x_max=.5*xy_span_pml, 
                      y_min=-.5*xy_span_pml, y_max=.5*xy_span_pml, 
                      z=Region_max,
                      apod="Start", apod_center=apod_center, apod_time_width=apod_start_w)

mon.power_monitor(name='xy_mid', montype="2D Z-normal", plane="xy",
                      x_min=-.5*xy_span_pml, x_max=.5*xy_span_pml, 
                      y_min=-.5*xy_span_pml, y_max=.5*xy_span_pml, 
                      z=0,
                      apod="Start", apod_center=apod_center, apod_time_width=apod_start_w)

mon.power_monitor(name='xy_exofeatured', montype="2D Z-normal", plane="xy",
                      x_min=-.5*xy_span_pml, x_max=.5*xy_span_pml, 
                      y_min=-.5*xy_span_pml, y_max=.5*xy_span_pml, 
                      z=Region_min,
                      apod="Start", apod_center=apod_center, apod_time_width=apod_start_w)


mon.power_monitor(name='xz_middle', montype="2D Y-normal", plane="xz",
                      x_min=-.5*xy_span_pml, x_max=.5*xy_span_pml, 
                      z_min=Region_min-pml_thickness, z_max=Region_max+pml_thickness, 
                      y=0,
                      apod="Start", apod_center=apod_center, apod_time_width=apod_start_w)    

mon.power_monitor(name='xz_edge', montype="2D Y-normal", plane="xz",
                      x_min=-.5*xy_span_pml, x_max=.5*xy_span_pml, 
                      z_min=Region_min-pml_thickness, z_max=Region_max+pml_thickness,  
                      y = .5*cav.xy_span,
                      apod="Start", apod_center=apod_center, apod_time_width=apod_start_w)    

mon.power_monitor(name='yz_middle', montype="2D X-normal", plane="yz",
                      y_min=-.5*xy_span_pml, y_max=.5*xy_span_pml, 
                      z_min=Region_min-pml_thickness, z_max=Region_max+pml_thickness, 
                      x=0,
                      apod="Start", apod_center=apod_center, apod_time_width=apod_start_w)    

mon.power_monitor(name='yz_edge', montype="2D X-normal", plane="yz",
                      y_min=-.5*xy_span_pml, y_max=.5*xy_span_pml, 
                      z_min=Region_min-pml_thickness, z_max=Region_max+pml_thickness, 
                      x=.5*cav.xy_span,
                      apod="Start", apod_center=apod_center, apod_time_width=apod_start_w)  

objects = analysisobjects.AnalysisObjects(fdtd=fdtd)

xy_span = cav.xy_span
z_span = abs(Region_min) + abs(Region_max)

objects.farfield(xy_span=xy_span, z_span=z_span, x=0, y=0, z=0.55e-06, theta_max=90,
                 N_theta=180, Nphi=91)

objects.Qfactor()

objects.mode_volume_3D(xy_span=xy_span, z_span=z_span, x=0, y=0, z=0.55e-06)

fdtd.select('FDTD');
fdtd.set('x min bc','PML');
fdtd.set('y min bc','PML');  

if save:
    fdtd.save(filename)
    
# fdtd.run()
