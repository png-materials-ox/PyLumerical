# -*- coding: utf-8 -*-
"""
Created on Sun May 26 16:30:57 2024

@author: mans4209
"""
import matplotlib.pyplot as plt
from  matplotlib import patches
from matplotlib.figure import Figure
from matplotlib.colors import LogNorm
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
        """
        Retrieves and computes the magnitude of the electric field from FDTD simulation data.

        This function extracts time-domain electric field components (Ex, Ey, Ez) from the 
        simulation data at the monitor named 'Qanalysis::t1'. It then calculates the total 
        electric field magnitude as:

            E = sqrt(Ex^2 + Ey^2 + Ez^2)

        Returns:
            pd.DataFrame: A DataFrame containing the time values ('t') and the corresponding 
            electric field components ('Ex', 'Ey', 'Ez') along with the computed field magnitude ('E').
        """
        t = self.fdtd.getdata('Qanalysis::t1', 't').flatten()
        Ex = np.squeeze(self.fdtd.getdata('Qanalysis::t1', 'Ex'))
        Ey = np.squeeze(self.fdtd.getdata('Qanalysis::t1', 'Ey'))
        Ez = np.squeeze(self.fdtd.getdata('Qanalysis::t1', 'Ez'))
        
        E = np.sqrt(Ex**2 + Ey**2 + Ez**2)
        return pd.DataFrame(data={'t':t, 'Ex':Ex, 'Ey':Ey, 'Ez':Ez, 'E':E})    
    
    def f_spectrum(self, plotting=False, saveplot=False) -> list:
        """
        Retrieves the frequency-domain spectrum from FDTD simulation data.

        This function extracts the frequency-domain response (`f_spectrum`) from the 
        analysis group 'Qanalysis'. If the analysis data is not available, it ensures 
        that the analysis has been run.

        Args:
            plotting (bool, optional): If True, plots the frequency spectrum.
            saveplot (bool, optional): If True and `plotting` is enabled, saves the plot as 'f_spectrum.jpg'.

        Returns:
            list: A flattened list containing the frequency-domain spectrum data.

        Notes:
            - The function assumes that the analysis has already been performed in the simulation.
            - If `plotting` is enabled, the spectrum is visualized using Matplotlib.
            - If `saveplot` is enabled, the plot is saved to 'f_spectrum.jpg' before being displayed.
        """
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
        """
        Retrieves the angular frequency values from FDTD simulation data.

        This function extracts the angular frequency (`w`) from the analysis group 'Qanalysis'. 
        If the analysis data is not available, it ensures that the analysis has been run.

        Returns:
            np.ndarray: A flattened array containing the angular frequency values.

        Notes:
            - The function checks if 'Qanalysis' data is available; if not, it triggers the analysis.
            - The returned values correspond to the angular frequency ω (rad/s).
        """
        if not self.fdtd.havedata('Qanalysis'):
            self.fdtd.runanalysis     
        return self.fdtd.getresult('Qanalysis', 'w').flatten()
        
    def resonances(self) -> pd.DataFrame:
        """
        Retrieves and organizes resonance parameters from FDTD simulation data.

        This function extracts resonance data from the analysis group 'Qanalysis'. If the 
        analysis data is not available, it ensures that the analysis has been run. The 
        extracted resonance parameters include central frequency (f0), decay rate, quality 
        factor (Q), amplitude, phase, and error.

        Returns:
            pd.DataFrame: A DataFrame containing the following resonance parameters:
                - 'f0': Central resonance frequency.
                - 'decay': Decay rate of the resonance.
                - 'Q': Quality factor of the resonance.
                - 'amp': Amplitude of the resonance mode.
                - 'phase': Phase of the resonance mode.
                - 'error': Estimation error of the resonance parameters.

        Notes:
            - The function checks if 'Qanalysis' data is available; if not, it triggers the analysis.
            - The returned DataFrame provides a structured view of resonance characteristics.
        """
        if not self.fdtd.havedata('Qanalysis'):
            self.fdtd.runanalysis     
        
        res = self.fdtd.getresult('Qanalysis', 'resonances')
        df = pd.DataFrame(data={'f0':res[:,0], 'decay':res[:,1], 'Q':res[:,2],'amp':res[:,3], 'phase':res[:,4], 'error':res[:,5]})
        return df
    
    def complex_eigenfrequency(self, omega0=1, Q=100):
        """
        Computes the complex eigenfrequency of a resonance mode.

        This function calculates the complex eigenfrequency as:

            ω_complex = ω0 - 1j * (2 * Q)

        where:
            - ω0 (omega0) is the central resonance frequency (default: 1).
            - Q is the quality factor of the resonance (default: 100).
            - The imaginary part represents the damping rate.

        Args:
            omega0 (float, optional): Central resonance frequency. Default is 1.
            Q (float, optional): Quality factor of the resonance. Default is 100.

        Returns:
            complex: The computed complex eigenfrequency.

        Notes:
            - The imaginary part is defined as -2Q, which models energy dissipation.
            - A higher Q leads to a smaller imaginary component, indicating lower losses.
        """
        omega_imag = 2*Q
        return omega0 - 1j * omega_imag
    
    # def transfer_function(self, eigenfreqs=[], w=[], wr=[], Q=[]):
    #     r = 1
    #     b = 0
        
    #     return (r)/(1j* (w - wr) - (wr/(2*Q)))
        
        
    def farfield_analysis(self, NA=0.9):
        """
        Performs a far-field optical analysis using FDTD simulation data.

        This function extracts and visualizes various far-field optical properties, including 
        total transmission, Poynting vector distribution, angular power distribution, Purcell 
        factor, and power within a specified numerical aperture (NA). If the 'farfield' data 
        is not available, it ensures that the analysis is performed.

        Args:
            NA (float, optional): Numerical aperture used to compute the power within an angular 
                cone. Default is 0.9.

        Returns:
            dict: A dictionary containing the following computed far-field parameters:
                - 'f': Frequency values from the far-field analysis.
                - 'sp': Source power spectrum.
                - 'dp': Dipole power spectrum.
                - 'wlen': Wavelength values (nm).
                - 'Fp': Purcell factor as a function of wavelength.
                - 'cone_angle': Angular cone corresponding to the given NA (degrees).
                - 'T34': Normalized power within the angular cone.
                - 'theta': Angular distribution of power (degrees).
                - 'P_v_theta': Power distribution as a function of angle.
                - 'Twlen': Wavelength values for total transmission (nm).
                - 'trans': Total transmission spectrum.

        Notes:
            - The function checks if 'farfield' data exists before extracting results.
            - It visualizes key results, including transmission, power vs. angle, and Purcell factor.
            - The total power within an angular cone is computed by integrating over a 
            specified angle derived from the numerical aperture (NA).
            - Uses Lumerical FDTD functions for field extraction, integration, and visualization.
        """

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
        """
        Computes the 2D mode volume of quasi-normal modes (QNM) using FDTD simulation data.

        This function extracts mode volumes along the xz- and yz-planes from a 2D mode volume 
        analysis and integrates the electric field energy density over space to obtain absolute 
        and wavelength-normalized mode volumes. It supports dipole position adjustments and 
        processes results from an FDTD simulation.

        Args:
            dipole_shift (float, optional): The shift in the dipole position along the z-axis.
                Default is 0.
            resonances (pd.DataFrame, optional): A DataFrame containing resonance information
                (frequencies, Q-factors, etc.). If not provided, it is obtained via self.resonances().

        Returns:
            pd.DataFrame: A DataFrame containing the computed mode volumes, including:
                - 'Vol_abs_xz': Absolute mode volume in the xz-plane (μm³).
                - 'Vol_lam_xz': Wavelength-normalized mode volume in the xz-plane (λ³).
                - 'Vol_abs_yz': Absolute mode volume in the yz-plane (μm³).
                - 'Vol_lam_yz': Wavelength-normalized mode volume in the yz-plane (λ³).
                - 'Vol_abs_avg': Averaged absolute mode volume (μm³).
                - 'Vol_lam_avg': Averaged wavelength-normalized mode volume (λ³).

        Notes:
            - The function retrieves the mode field data and refractive index from FDTD.
            - It calculates mode volume by integrating over the electric field energy density.
            - The numerical aperture (NA) is used to define an integration range.
            - Uses the quasi-normal mode frequencies derived from the resonance data.
            - A progress bar is displayed during computations for multiple resonances.
        """
        
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
        """
        Retrieves the 2D mode volume results from an FDTD simulation.

        This function extracts precomputed mode volume data along the xz- and yz-planes 
        from an existing FDTD analysis. It ensures that the data is available and, if not, 
        triggers the necessary analysis before extracting the results.

        Returns:
            pd.DataFrame: A DataFrame containing mode volume results, including:
                - 'Vol_abs_xz': Absolute mode volume in the xz-plane (μm³).
                - 'Vol_lam_xz': Wavelength-normalized mode volume in the xz-plane (λ³).
                - 'Vol_abs_yz': Absolute mode volume in the yz-plane (μm³).
                - 'Vol_lam_yz': Wavelength-normalized mode volume in the yz-plane (λ³).
                - 'Vol_abs_avg': Averaged absolute mode volume (μm³).
                - 'Vol_lam_avg': Averaged wavelength-normalized mode volume (λ³).

        Notes:
            - If the 'mode volume 2D' dataset is unavailable, the function runs the necessary FDTD analysis.
            - Extracts results from the FDTD solver and organizes them into a structured DataFrame.
            - The wavelength-normalized mode volumes are given in units of λ³.
        """
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
    
    def get_field_profiles(self, monitor):
        """
        Retrieves the electromagnetic field profiles from a specified FDTD monitor.

        This function extracts the spatial and frequency-dependent electric field (E) 
        and polarization (P) components from a given monitor in the FDTD simulation.

        Args:
            monitor (str): The name of the FDTD monitor from which to extract field data.

        Returns:
            dict: A dictionary containing the extracted field components and spatial coordinates:
                - 'f': Frequency values.
                - 'x': X-coordinates of the monitor grid.
                - 'y': Y-coordinates of the monitor grid.
                - 'z': Z-coordinates of the monitor grid.
                - 'Ex': X-component of the electric field.
                - 'Ey': Y-component of the electric field.
                - 'Ez': Z-component of the electric field.
                - 'Px': X-component of the polarization field.
                - 'Py': Y-component of the polarization field.
                - 'Pz': Z-component of the polarization field.

        Notes:
            - The function assumes that the specified monitor has the required field components.
            - Ensure that the monitor name corresponds to a valid monitor in the FDTD simulation.
        """
        f = self.fdtd.getresult(monitor, 'f')
        x = self.fdtd.getresult(monitor, 'x')
        y = self.fdtd.getresult(monitor, 'y')
        z = self.fdtd.getresult(monitor, 'z')
        Ex = self.fdtd.getresult(monitor, 'Ex')
        Ey = self.fdtd.getresult(monitor, 'Ey')
        Ez = self.fdtd.getresult(monitor, 'Ez')
        Px = self.fdtd.getresult(monitor, 'Px')
        Py = self.fdtd.getresult(monitor, 'Py')
        Pz = self.fdtd.getresult(monitor, 'Pz')
        
        return {'f':f, 'x':x, 'y':y, 'z':z, 'Ex':Ex, 'Ey':Ey, 'Ez':Ez, 'Px':Px, 'Py':Py, 'Pz':Pz}
    
    def get_single_field_profile(self, monitor, wlen=637e-09):
        """
        Extracts a single field profile at a specified wavelength from an FDTD monitor.

        This function retrieves the electric field components from a given monitor and
        computes the total field magnitude at the closest available wavelength.

        Args:
            monitor (str): The name of the FDTD monitor from which to extract the field data.
            wlen (float, optional): The target wavelength (in meters) for extracting the field profile.
                                    Default is 637 nm (commonly used for NV centers).

        Returns:
            dict: A dictionary containing spatial coordinates and the field magnitude:
                - 'x': 1D array of x-coordinates (in micrometers).
                - 'y': 1D array of y-coordinates (in micrometers).
                - 'X': 2D meshgrid of x-coordinates.
                - 'Y': 2D meshgrid of y-coordinates.
                - 'Z': 2D array of the total electric field magnitude at the specified wavelength.

        Notes:
            - The function retrieves all frequency components from the monitor and finds 
            the closest match to the desired wavelength.
            - The electric field magnitude is computed as sqrt(Ex² + Ey² + Ez²).
            - The spatial coordinates are converted from meters to micrometers for visualization.
        """
        field_profiles = self.get_field_profiles(monitor)
        f = field_profiles['f'].flatten()
        wavelength = sc.c / f
        
        # TODO CHANGE FOR GENERALITY
        x = field_profiles['x'].flatten()*1e06
        y = field_profiles['z'].flatten()*1e06
        X, Y = np.meshgrid(y, x)
        
        Ex = np.squeeze(field_profiles['Ex'])
        Ey = np.squeeze(field_profiles['Ey'])
        Ez = np.squeeze(field_profiles['Ez'])
        Z = np.sqrt(Ex**2 + Ey**2 + Ez**2)
        
        Z = np.real(np.sqrt(Ex**2 + Ey**2 + Ez**2))
        
        wlen_idx = np.abs(wavelength-wlen).argmin()
        
        Z = np.squeeze(Z)[:,:,wlen_idx]
        
        return {'x':x, 'y':y, 'X':X, 'Y':Y, 'Z':Z}
    
    def plot_field_profile(self, monitor, wlen=637e-09, savefig=False, 
                           title='', filename=''):
        """
        Plots the field profile at a specified wavelength from a given FDTD monitor.

        This function generates a 2D color plot of the electric field magnitude at the
        specified wavelength. The plot is created based on the electric field components
        (Ex, Ey, Ez) from the FDTD simulation results.

        Args:
            monitor (str): The name of the FDTD monitor from which to extract the field data.
            wlen (float, optional): The target wavelength (in meters) for extracting the field profile.
                                    Default is 637 nm (commonly used for NV centers).
            savefig (bool, optional): If True, saves the plot as an image file. Default is False.
            title (str, optional): The title of the plot. Default is an empty string.
            filename (str, optional): The filename to save the plot if `savefig` is True. Default is an empty string.

        Returns:
            None: The function displays the plot and optionally saves it.

        Notes:
            - The function retrieves the field components (Ex, Ey, Ez) from the monitor and computes
            the total electric field magnitude as sqrt(Ex² + Ey² + Ez²).
            - The spatial coordinates (x, y) are converted to micrometers (µm) for visualization.
            - The closest available wavelength from the monitor data is selected based on the provided `wlen`.
            - The color plot uses a linear colormap for visualization, with the field magnitude displayed.

        Example:
            plot_field_profile('monitor1', wlen=637e-09, savefig=True, 
                            title='Field Profile at 637 nm', filename='field_profile.png')
        """

        field_profiles = self.get_field_profiles(monitor)
        f = field_profiles['f'].flatten()
        wavelength = sc.c / f
        
        # TODO CHANGE FOR GENERALITY
        X, Y = np.meshgrid(field_profiles['z'].flatten()*1e06, field_profiles['x'].flatten()*1e06)
        
        Ex = np.squeeze(field_profiles['Ex'])
        Ey = np.squeeze(field_profiles['Ey'])
        Ez = np.squeeze(field_profiles['Ez'])
        Z = np.sqrt(Ex**2 + Ey**2 + Ez**2)
        
        Z = np.real(np.sqrt(Ex**2 + Ey**2 + Ez**2))
        
        wlen_idx = np.abs(wavelength-wlen).argmin()
        
        Z = np.squeeze(Z)[:,:,wlen_idx]
        
        
        fig, ax = plt.subplots()
        
        # c = ax.pcolor(Yyz, Xyz, Zyz, shading='auto',
        #                norm=LogNorm(vmin=Zyz.min(), vmax=Zyz.max()), cmap='afmhot')
        c = ax.pcolor(Y, X, Z, shading='auto')
        fig.colorbar(c, ax=ax)
        ax.set_title(title)
        ax.set_xlabel('y ($\mu$m)')
        ax.set_ylabel('z ($\mu$m)')
        if savefig:
            plt.savefig(filename, bbox_inches="tight")
        plt.show()
                