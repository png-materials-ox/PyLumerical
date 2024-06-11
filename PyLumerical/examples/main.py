# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 17:46:36 2024

@author: mans4209
"""

# -*- coding: utf-8 -*-
"""
Created on Wed May  8 14:35:52 2024

@author: mans4209
"""

import numpy as np 
import imp
import os
from PyLumerical import cavity, simulation, source, monitor, analysisobjects

wlen = .637e-06
roc = 8e-06
qlist =  np.arange(4, 13, 1)

ncav = 2.41

for q in qlist:
    
    # os.add_dll_directory("C:\\Program Files\\Lumerical\\v232\\api\\python\\")
    lumapi = imp.load_source("lumapi","C:\\Program Files\\Lumerical\\v232\\api\\python\\lumapi.py")

    save = False

    fdtd = lumapi.FDTD()

    filename = "20240605_monolithic_diamond_cavity_roc_q{:d}_Lsub20nm_theta90//".format(q)
    
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
    # sim.add_mesh(name='planar mesh', dx=0.01e-06, dy=0.01e-06, dz=0.01e-06, 
    #              based_on_struct=True, struct="planar mirror")
    
    # sim.add_mesh(name='feat mesh', dx=0.01e-06, dy=0.01e-06, dz=0.01e-06, 
    #              based_on_struct=True, struct="featured mirror")
    
    # ADD SOURCE
    ##############################################################################
    src = source.Source(fdtd=fdtd, wlen=wlen)
    src.dipole(theta=90, shift=0, emission_width=100e-09)
    
    # ADD MONITORS
    ##############################################################################
    pml_thickness = 0e-09
    xy_span_pml = cav.xy_span + pml_thickness
    
    mon = monitor.Monitor(fdtd=fdtd)
    
    # mon.Q_monitor(Qmonitor_zspan=10, Qmonitor_zlayer=1, t_sample=10, dipole_shift=0)
    
    # mon.index_monitor(name="n", monitor_type="2D Y-normal",
    #                         x_min=-.5*xy_span_pml, x_max=.5*xy_span_pml, 
    #                         z_min=Region_min-.5*pml_thickness, z_max=Region_max + .5*pml_thickness)
    
    
    objects = analysisobjects.AnalysisObjects(fdtd=fdtd, mon=mon)
    
    xy_span = cav.xy_span
    z_span = abs(Region_min) + abs(Region_max)
    
    
    objects = analysisobjects.AnalysisObjects(fdtd=fdtd, mon=mon)
    
    xy_span = cav.xy_span
    z_span = abs(Region_min) + abs(Region_max)
    
    
    objects.Qfactor(x=0.1e-06, y=0.1e-06, z=0.1e-06)
    
    # fdtd.unselectall()
    
    z_ff = .5*(abs(Region_max) - abs(Region_min))
    objects.farfield(xy_span=xy_span, z_span=z_span, x=0, y=0, z=z_ff, theta_max=90,
                      N_theta=180, Nphi=91)
    
    # objects.mode_volume_3D(xy_span=xy_span, z_span=z_span, x=0, y=0, z=0.55e-06)
    
    objects.mode_volume_2D(xy_span_pml=xy_span_pml, apod="Start", apod_center=100e-15, apod_start_w=10e-15, z_min=Region_min - pml_thickness, z_max=Region_max + pml_thickness)
    
    # fdtd.select('FDTD');
    # fdtd.set('x min bc','PML');
    # fdtd.set('y min bc','PML');  
    
    if save:
        if not os.path.exists(filename):
            os.mkdir(filename)
            
        try:
            fdtd.save(filename + 'fsp')
        except:
            FileExistsError("Couldn't write file")
        # fdtd.savedata(filename + 'lnp')
        
    # fdtd.run()
