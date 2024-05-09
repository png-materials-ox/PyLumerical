# -*- coding: utf-8 -*-
class Source:
    """A class representing a light source in FDTD simulations.

    Author       : Gareth Sion jones
    Date         : May, 2024
    Version      : 0.1
    Organisation : University of Oxford

    """
    
    def __init__(self, fdtd=None, wlen=637e-06):
        """
        Constructor method
        Params:
            fdtd (object, optional):
                Lumerical FDTD object for simulation. Default is None.
            wlen (float, optional): 
                    Wavelength of the light source in meters. Default is 637e-06.
                
        Returns:
            
        Raises:
            IOError: 
                If no Lumerical FDTD object is defined in the constructor.
        """
        self.wlen = wlen
        self.fdtd = fdtd
        if not self.fdtd:
            raise IOError ("A lumerical FDTD object must be defined in the constructor")
    
    def dipole(self, theta=90, shift=0, emission_width=100e-09):
        """Add a dipole source to the simulation.

        This method adds a dipole source to the simulation with specified parameters.

        Params:
            theta (float, optional): 
                Angle of the dipole orientation in degrees. Default is 90.
            shift (float, optional): 
                Shift along the z-axis from the origin in meters. Default is 0.
            emission_width (float, optional): 
                Wavelength span of the emission in meters. Default is 100e-09.

        """
        self.fdtd.adddipole()
        self.fdtd.set('theta',90)
        self.fdtd.set('x',0)
        self.fdtd.set('y',0)
        self.fdtd.set('z',shift)
        self.fdtd.set('center wavelength',self.wlen)
        self.fdtd.set('wavelength span',emission_width)
        self.fdtd.set('optimize for short pulse',0)

