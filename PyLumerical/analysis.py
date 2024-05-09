# -*- coding: utf-8 -*-

import numpy as np
import imp
# import os

# os.add_dll_directory("C:\\Program Files\\Lumerical\\v232\\api\\python\\")
lumapi = imp.load_source("lumapi","C:\\Program Files\\Lumerical\\v232\\api\\python\\lumapi.py")

fdtd = lumapi.FDTD()

pdir = 'Y:\\from_lumerical\\from_sda3\\lumerical\\Gareth\\characterising_cavity_designs\\RoC4um_CavWidth3128nm_CavLength1080nm_Nf11Np16_MeshAcc3_q8_nmedium2p41runtime=8e-12_bareCavity\\'
varname = pdir + 'lnp.ldf'    

fdtd.loaddata(varname)

# Monitor names
m0 = 'n'
m1 = 'xy_exoplanar'
m2 = 'xy_middle'
m3 = 'xy_exofeatured'
m4 = 'xz_middle'
m5 = 'xz_edge'
m6 = 'yz_middle'
m7 = 'yz_edge'



# ## Decide whether to run the convergence test
# Q_test=0;

# fullrun=0;
# savedatafiles=1;
# savemonitor=1;

# lam_max=638.0; #nm
# lam_min=630.0; #nm


# #### Simulation info ####
# ax_x = getdata(m2,"x"); #Gets data from monitor m2 called x
# ax_y = getdata(m2,"y");
# ax_z = getdata(m4,"z"); #Gets data from monitor m4 called z
# ax_x_pt=size(ax_x);     # Points in the axis
# ax_z_pt=size(ax_z);
# ax_x_res=XY_span/ax_x_pt; 	
# i_xmid=floor(ax_x_pt(1)/2)+1;
# ax_f= pinch(getdata(m2,"f")); #Gets frequency data and removes singleton dimensions
# ax_lam= c/ax_f; #Calculates wavelengths = c / f
# lam_emt= lambda_res;    #Emitter wavelength

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