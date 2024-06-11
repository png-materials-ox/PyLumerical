# -*- coding: utf-8 -*-
"""
Created on Sun May 26 16:30:57 2024

@author: mans4209
"""
import matplotlib.pyplot as plt
from  matplotlib import patches
from matplotlib.figure import Figure
from matplotlib import rcParams
import numpy as np
import pandas as pd
import scipy.constants as sc
import progressbar

class LumericalAnalysis:
    
    def __init__(self, fdtd=None):
        self.fdtd = fdtd
        if not self.fdtd:
            raise IOError ("A lumerical FDTD object must be defined in the constructor")
    
    def electric_field_magnitude(self) -> pd.DataFrame:
        t = self.fdtd.getdata('Qanalysis::t1', 't').flatten()
        Ex = np.squeeze(self.fdtd.getdata('Qanalysis::t1', 'Ex'))
        Ey = np.squeeze(self.fdtd.getdata('Qanalysis::t1', 'Ey'))
        Ez = np.squeeze(self.fdtd.getdata('Qanalysis::t1', 'Ez'))
        
        E = np.sqrt(Ex**2 + Ey**2 + Ez**2)
        return pd.DataFrame(data={'t':t, 'Ex':Ex, 'Ey':Ey, 'Ez':Ez, 'E':E})    
    
    def f_spectrum(self, plotting=False, saveplot=False) -> list:
        if not self.fdtd.havedata('Qanalysis'):
            self.fdtd.runanalysis
            
        # wlen = self.fdtd.getresult('Qanalysis', 'w').flatten()
        f_spectrum = self.fdtd.getresult('Qanalysis', 'f_spectrum').flatten()
        
        if plotting:
            plt.plot(f_spectrum)
            plt.ylabel('Response (AU)')
            if saveplot:
                plt.savefig('f_spectrum.jpg', bbox_inches="tight")
            plt.show()
            
        return f_spectrum
    
    def omega(self):
        if not self.fdtd.havedata('Qanalysis'):
            self.fdtd.runanalysis     
        return self.fdtd.getresult('Qanalysis', 'w').flatten()
        
    def resonances(self) -> pd.DataFrame:
        if not self.fdtd.havedata('Qanalysis'):
            self.fdtd.runanalysis     
        
        res = self.fdtd.getresult('Qanalysis', 'resonances')
        df = pd.DataFrame(data={'f0':res[:,0], 'decay':res[:,1], 'Q':res[:,2],'amp':res[:,3], 'phase':res[:,4], 'error':res[:,5]})
        return df
    
    def complex_eigenfrequency(self, omega0=1, Q=100):
        omega_imag = 2*Q
        return omega0 - 1j * omega_imag
    
    # def transfer_function(self, eigenfreqs=[], w=[], wr=[], Q=[]):
    #     r = 1
    #     b = 0
        
    #     return (r)/(1j* (w - wr) - (wr/(2*Q)))
        
        
    
    def farfield_analysis(self, NA=0.9):
        # start_unix_time = now;
        if not self.fdtd.havedata('farfield'):
            self.fdtd.runanalysis
            
                
        ##############################################
        # Plot results from single simulation
        
        ## Scaling factors for plots
        um = 1e6
        nm = 1e9
        
        ##############################################
        # Get field from z-normal monitor (farfield::z2)
        x = self.fdtd.getdata("farfield::z2", "x")
        y = self.fdtd.getdata("farfield::z2", "y")
        E2 = self.fdtd.getelectric("farfield::z2")
        self.fdtd.image(x*um, y*um, self.fdtd.pinch(E2, 4, 1), "x (um)", "y (um)", "")
        
        ##############################################
        # Get results from "farfield" analysis group
        ## Total transmission
        T = self.fdtd.getresult("farfield", "T")
        fig = plt.figure()
        plt.plot(T['lambda']*nm, T['T']) 
        plt.xlabel('$\lambda$ (nm)')
        plt.ylabel('Total transmission')
        plt.show()
        
        ## Poynting vector in far field
        S = self.fdtd.getresult("farfield", "farfield")
        self.fdtd.visualize(S)
        
        ## P_vs_theta
        P_vs_theta = self.fdtd.getresult("farfield", "P_vs_theta")
        fig = plt.figure()
        plt.plot(P_vs_theta['theta_degrees'], self.fdtd.pinch(P_vs_theta['P'], 2, 1))
        plt.xlabel('theta (degrees)')
        plt.ylabel('Normalized Power')
        plt.show()
        
        ## Purcell factor
        Purcell = self.fdtd.getresult("farfield", "Purcell")
        fig = plt.figure()
        plt.plot(Purcell['lambda']*nm, Purcell['purcell'])
        plt.xlabel('lambda (nm)')
        plt.ylabel('Purcell factor')
        plt.show()
        
        ## Power in angular cone
        # create 2D mesh of theta
        cone_angle = np.arcsin(NA)*180/np.pi
        Theta = self.fdtd.meshgridx(P_vs_theta['theta_radians'],P_vs_theta['f'])
        # integrate over specified cone
        T34 = 0.5*2*np.pi*np.real(self.fdtd.integrate(P_vs_theta['P']*self.fdtd.sin(Theta)*(Theta<=cone_angle*np.pi/180),1,P_vs_theta['theta_radians']))
        T34 = T34.flatten()
        
        # plot final results
        fig = plt.figure()
        plt.plot(P_vs_theta['lambda']*nm,T34, label="Normalized power")
        plt.plot(P_vs_theta['lambda']*nm,T34/Purcell['purcell'], label="Optical extraction efficiency")
        plt.xlabel('wavelength (nm)')
        plt.legend()
        plt.show()
        
        ff = self.fdtd.getresult('farfield', 'farfield')
        f = self.fdtd.getdata("farfield::x2", "f").flatten()
        #f = getdata('xy_mid','f');
        wlen = P_vs_theta['lambda'].flatten()
        Fp = Purcell['purcell'].flatten()
        theta = P_vs_theta['theta_degrees'].flatten()
        P_v_theta = self.fdtd.pinch(P_vs_theta['P'], 2, 1).flatten()
        Twlen = T['lambda'].flatten()*nm
        trans = T['T'].flatten()
        sp = self.fdtd.sourcepower(f).flatten()
        dp = self.fdtd.dipolepower(f).flatten()
        
        fig = plt.figure()
        plt.plot(wlen, trans*(self.fdtd.sourcepower(f)/self.fdtd.dipolepower(f)))
        plt.show()
        
        #vtksave(fname + ".vtu",ff);
        return {'f':f, 'sp':sp, 'dp':dp, 'wlen':wlen, 'Fp':Fp, 'cone_angle':cone_angle, 
                'T34':T34, 'theta':theta, 'P_v_theta':P_v_theta, 
                'Twlen':Twlen, 'trans':trans}
        
    
    def mode_volume_2D_QNM(self, dipole_shift=0, resonances=None):
        
        if resonances is None:
            resonances = self.resonances()
            
        wr = 2 * np.pi * resonances.f0.values
        Q = resonances.Q.values
        f_spectrum = self.f_spectrum(plotting=False, saveplot=False)
        
        w_qnm = wr - 1j*(1/(2*Q))
        
        x = self.fdtd.getdata('mode volume 2D::xy_middle',"x")   
        y = self.fdtd.getdata('mode volume 2D::xy_middle',"y")
        z = self.fdtd.getdata('mode volume 2D::xz_middle',"z")   
        
        x_pt = np.shape(x)           # Points in the axis
        z_pt = np.shape(z)           # Points in the z axis            
        
        f = self.fdtd.getresult('mode volume 2D::xy_middle', 'f').flatten()
        fres_indices = self.fdtd.find(f, resonances.f0.values).flatten()
        fres_indices = [int(f) for f in fres_indices]
        
        wres = 2 * np.pi * f[fres_indices]
        w = 2*np.pi*f
        
        midpoint_z_ind = self.fdtd.round(.5*self.fdtd.size(x , 1))    # find the z midpoint index value
        midpoint_x_ind = int(self.fdtd.floor(x_pt[0]/2)+1)
        
        # Integration ranges (0 -> rmax, zmin -> zmax)
        z_int_range = z; 
        r_int_range = x[0:midpoint_x_ind]  
            
        # Get electric fields
        E_xz = self.fdtd.pinch(self.fdtd.getelectric('mode volume 2D::xz_middle'))
        E_yz = self.fdtd.pinch(self.fdtd.getelectric('mode volume 2D::yz_middle'))
        
        # Get refractive indices
        n_xz = self.fdtd.pinch(self.fdtd.getdata('mode volume 2D::n',"index_z"))
        
        nv_zpos_ind = int(self.fdtd.find(z, dipole_shift))            # Position of the NV
        n_nv_xz = np.real(n_xz[midpoint_x_ind, nv_zpos_ind]) # Refractive index at the position of the NV
        n_nv_yz = np.real(n_xz[midpoint_x_ind, nv_zpos_ind])                                   
        
        # Predefine arrays for looping
        Vol_abs_xz = np.ones(len(fres_indices))
        Vol_lam_xz = np.ones(len(fres_indices))
        
        Vol_abs_yz = np.ones(len(fres_indices))
        Vol_lam_yz = np.ones(len(fres_indices))
        
        with progressbar.ProgressBar(max_value=len(fres_indices)) as bar:
        
            for i in range(0, len(fres_indices)-1):
                
                n = n_xz[0:midpoint_x_ind,0:z_pt[0]]
                
                wm = w_qnm[i]
                wn = w_qnm[i+1]
                
                delta_wqnm = wm * np.real(n**2) - wn*np.real(n**2)
                
                
                Em_xz_res = self.fdtd.pinch(E_xz, 3, fres_indices[i])
                Em_xz_mid = Em_xz_res[0:midpoint_x_ind, 0:z_pt[0]]                     # half of the x-z cut
                
                En_xz_res = self.fdtd.pinch(E_xz, 3, fres_indices[i+1])
                En_xz_mid = En_xz_res[0:midpoint_x_ind,0:z_pt[0]]  
                
                E_qnm_xz = Em_xz_mid * delta_wqnm * En_xz_mid
                eps_E_xz = np.real(n**2)*Em_xz_mid
                
                
                eps_E_xz_at_nv = eps_E_xz[-1, nv_zpos_ind] # THIS IS WRONG!!!!!!!!
                eps_E_xz_max = max(eps_E_xz.flatten())
                
                V0_xz = 2*np.pi*E_qnm_xz
                Vol_raw1_xz = self.fdtd.integrate(V0_xz,2,z_int_range)
                # Vol_raw2_xz = 1e18*abs(self.fdtd.integrate(Vol_raw1_xz*r_int_range,1,r_int_range))
                Vol_raw2_xz = abs(self.fdtd.integrate(Vol_raw1_xz*r_int_range,1,r_int_range))
                
                
                ############ TESTING #################
                
                eps = np.real(n**2)  # use real part of esp
                sz = self.fdtd.size(E_qnm_xz).flatten()
                sz = [int(s) for s in sz]
                
                
                E_qnm_test = np.sqrt(Em_xz_mid) * delta_wqnm * np.sqrt(En_xz_mid)
                
                mode = eps*E_qnm_test
                
                mode[0:, 1:] = mode[0:, 1:] / (sc.epsilon_0 * max(mode[0:,1:].flatten()))
                
                V1 = self.fdtd.integrate(mode,2,z_int_range)
                V2 = self.fdtd.integrate(V1, 1, r_int_range)
                
                ############         #################
                
                
                
                n_i = 2.41 # CHANGE
                
                Vol_abs_xz[i] = Vol_raw2_xz/(eps_E_xz_max)                        # unit:um^3, for air-like mode
                #Vol_lam_xz[i] = Vol_abs_xz[i]/(wres[i]**3 * 1e18)   
                Vol_lam_xz[i] = Vol_abs_xz[i] * (n_i / (sc.c/(2*np.pi*wres[i])))**3 * 1e-18
                
                
                # REPEAT for YZ
                Em_yz_res = self.fdtd.pinch(E_yz, 3, fres_indices[i])
                Em_yz_mid = Em_yz_res[0:midpoint_x_ind, 0:z_pt[0]]                     # half of the x-z cut
                
                En_yz_res = self.fdtd.pinch(E_yz, 3, fres_indices[i+1])
                En_yz_mid = En_yz_res[0:midpoint_x_ind,0:z_pt[0]]  
                
                E_qnm_yz = Em_yz_mid * delta_wqnm * En_yz_mid
                eps_E_yz = np.real(n**2)*Em_yz_mid
                
                
                eps_E_yz_at_nv = eps_E_yz[-1, nv_zpos_ind] # THIS IS WRONG!!!!!!!!
                eps_E_yz_max = max(eps_E_yz.flatten())
                
                V0_yz = 2*np.pi*E_qnm_yz
                Vol_raw1_yz = self.fdtd.integrate(V0_yz,2,z_int_range)
                # Vol_raw2_yz = 1e18*abs(self.fdtd.integrate(Vol_raw1_yz*r_int_range,1,r_int_range))
                Vol_raw2_yz = abs(self.fdtd.integrate(Vol_raw1_yz*r_int_range,1,r_int_range))
                
                Vol_abs_yz[i] = Vol_raw2_yz/(eps_E_yz_max)
                # Vol_lam_yz[i] = Vol_abs_yz[i]/(wres[i]**3)   
                Vol_lam_yz[i] = Vol_abs_yz[i] * (n_i / (sc.c/(2*np.pi*wres[i])))**3 * 1e-18
                
                
                
                
                
                # eps = np.real(n^2)  # use real part of esp
                # sz = np.size(E2)
                
                # mode = eps*E2
                # for (i=1:sz(4)) {  # normalize each frequency  
                # mode(1:sz(1),1:sz(2),1:sz(3),i) =     mode(1:sz(1),1:sz(2),1:sz(3),i) /
                #                                   max(mode(1:sz(1),1:sz(2),1:sz(3),i));
                # }
                # V = integrate2(mode,1:3,x,y,z);
                
                
                bar.update(i)
         
            
        ## Average mode volume
        Vol_abs_avg = ((Vol_abs_xz + Vol_abs_yz)/2);
        Vol_lam_avg = (Vol_lam_xz + Vol_lam_yz)/2;
        n_nv_avg = (n_nv_xz + n_nv_yz)/2;
        In2_nv_avg = (eps_E_xz_at_nv + eps_E_yz_at_nv)/2;
        In2_max_avg = (eps_E_xz_max + eps_E_yz_max)/2;             # unit:lam^3    
        
        return pd.DataFrame(data={'Vol_abs_xz':Vol_abs_xz, 
                                  'Vol_lam_xz':Vol_lam_yz,
                                  'Vol_abs_yz':Vol_abs_xz,
                                  'Vol_lam_yz':Vol_lam_yz,
                                  'Vol_abs_avg':Vol_abs_avg,
                                  'Vol_lam_avg':Vol_lam_avg
                                  })
    
    def mode_volume_2D(self):
        if not self.fdtd.havedata('mode volume 2D'):
            self.fdtd.runanalysis
    
    
        Vol_abs_xz = self.fdtd.getresult('mode volume 2D', 'Vol_abs_xz').flatten()
        Vol_lam_xz = self.fdtd.getresult('mode volume 2D', 'Vol_lam_xz').flatten()
        Vol_abs_yz = self.fdtd.getresult('mode volume 2D', 'Vol_abs_yz').flatten()
        Vol_lam_yz = self.fdtd.getresult('mode volume 2D', 'Vol_lam_yz').flatten()
        Vol_abs_avg = self.fdtd.getresult('mode volume 2D', 'Vol_abs_avg').flatten()
        Vol_lam_avg = self.fdtd.getresult('mode volume 2D', 'Vol_lam_avg').flatten()
        
        return pd.DataFrame(data={'Vol_abs_xz':Vol_abs_xz, 
                                      'Vol_lam_xz':Vol_lam_yz,
                                      'Vol_abs_yz':Vol_abs_xz,
                                      'Vol_lam_yz':Vol_lam_yz,
                                      'Vol_abs_avg':Vol_abs_avg,
                                      'Vol_lam_avg':Vol_lam_avg
                                      })