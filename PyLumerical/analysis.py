# -*- coding: utf-8 -*-

import numpy as np
import imp
# import os
import scipy.constants as sc
import cavityanalysis
import lumericalanalysis
# import analysisbuilder

# os.add_dll_directory("C:\\Program Files\\Lumerical\\v232\\api\\python\\")
lumapi = imp.load_source("lumapi","C:\\Program Files\\Lumerical\\v232\\api\\python\\lumapi.py")

fdtd = lumapi.FDTD()

pdir = 'C:\\Users\\mans4209\\Documents\\LumericalFiles\\Philippa\\test\\'

fdtd.load(pdir + 'fsp.fsp')

analysis = lumericalanalysis.LumericalAnalysis(fdtd=fdtd)

resonances = analysis.resonances()
# w = analysis.omega()
w0 = 2*np.pi*resonances['f0'].values
Q = resonances['Q'].values

eigenfreqs = analysis.complex_eigenfrequency(omega0=w0, Q=Q) 

f_spectrum = analysis.f_spectrum(plotting=True, saveplot=False)

Efields = analysis.electric_field_magnitude()

# farfield = analysis.farfield_analysis()

V = analysis.mode_volume_2D_QNM(dipole_shift=0)

# Farfield Analysis


# # Qfactor analysis
# fdtd.runanalysis('Qanalysis')
# resonances = fdtd.getresult('Qanalysis', 'resonances')
# spectrum = fdtd.getresult('Qanalysis', 'spectrum')
# f_spectrum = fdtd.getresult('Qanalysis', 'f_spectrum')
# Q = fdtd.getresult('Qanalysis', 'Q')

# wres = resonances[:,0]


# monitor_names = {'m0':'n', 'm1':'xy_exoplanar', 'm2':'xy_middle',
#                  'm3':'xy_exofeatured', 'm4':'xz_middle','m5':'xz_edge',
#                  'm6':'yz_middle', 'm7':'yz_edge'}

# builder = analysisbuilder.AnalysisBuilder(fdtd=fdtd, pdir=pdir, monitor_names=monitor_names, min_wlen=630e-09, max_wlen=638e-09)

# fdtd.load(pdir + 'fsp.fsp')
# fdtd.loaddata(pdir + 'lnp.ldf')



# xy_span = fdtd.getnamed("structure::dielectric mediu", "x span")
# wlen = fdtd.getnamed("source", "center wavelength")

# x = fdtd.getdata(m2,"x") 
# y = fdtd.getdata(m2,"y")
# z = fdtd.getdata(m4,"z") 

# x_pts = np.size(x)     # Points in the axis
# z_pts = np.size(z)

# x_res = xy_span/x_pts # resolution of simulation region
# midpoint = np.floor(x_pts/2)+1

# f = np.squeeze(fdtd.getdata(m2,"f")) #Gets frequency data and removes singleton dimensions
# lam = sc.c/f

# ### Define the freq ROI for the peak interested ####

# max_wlen = 638.0*1e-09
# min_wlen = 630.0*1e-09

# max_wlen = fdtd.find(lam, max_wlen)
# min_wlen = fdtd.find(lam, min_wlen)
# w_range_max = (2*np.pi*sc.c)/min_wlen
# w_range_min = (2*np.pi*sc.c)/max_wlen
# f_range_max = sc.c/min_wlen
# f_range_min = sc.c/max_wlen

# cavana = cavityanalysis.CavityAnalysis(builder=builder, fdtd=fdtd)
# cavana.Qfactor()


# ## Decide whether to run the convergence test
# Q_test=0;

# fullrun=0;
# savedatafiles=1;
# savemonitor=1;

# lam_max=638.0; #nm
# lam_min=630.0; #nm


# ###################################################
# ### Define the freq ROI for the peak interested ####
# #lamoff = 0.2;
# lam_max_i=find(ax_lam,lam_max*1e-9);
# lam_min_i=find(ax_lam,lam_min*1e-9);
# w_range_max=2*pi*c/lam_min/1e-9;
# w_range_min=2*pi*c/lam_max/1e-9;
# f_range_max=c/lam_min/1e-9;
# f_range_min=c/lam_max/1e-9;
# ###################################################