Examples
========

3D Fabry-Perot Open-Access Optical Microcavity
----------------------------------------------

Building the Simulation File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The full example code for building the geometry can be found :doc:`here <example1_geometry>`

First, start by importing the necessary python libraries

.. code-block:: python
    
    import numpy as np 
    import imp
    # import os
    import cavity
    import simulation
    import source 
    import monitor

'imp' is a standard python library which allows us to import a local file. We need to link the main script to the 
Lumerical python API file which is stored locally on our PC within the Lumerical folder. For example, the lumapi
may be found at the following address (depending on distribution).

.. code-block:: python
    
    # os.add_dll_directory("C:\\Program Files\\Lumerical\\v232\\api\\python\\")
    lumapi = imp.load_source("lumapi","C:\\Program Files\\Lumerical\\v232\\api\\python\\lumapi.py")

Above, the **os.add_dll_directory** is only used if there is any difficulty adding the file.

We can choose whether or not to save the file, and asign a filename if necessary:

.. code-block:: python
    
    save = False
    filename = "testrun.fsp"

The FDTD object from the lumapi now needs to be instantiated:

.. code-block:: python
        
    fdtd = lumapi.FDTD()

With this, we can access the necessary Lumerical methods. Now, we can begin building the simulation geometry. 
We'll start by providing some necessary constants:

.. code-block:: python
    
    wlen = .637e-06 # centre wavelength
    roc = 2e-06     # cavity radius of curvature
    q = 4           # cavity longitudinal mode index
    ncav = 1        # refractive index of cavity medium (1 if open-access)

The Cavity() constructor now needs to be instantiated:

.. code-block:: python
    
    cav = cavity.Cavity(roc=roc, wlen=wlen, q=q, ncav=ncav, rough=False, fdtd=fdtd)
    cav.build_cavity(num_planar=10, num_feat=10, n1=2.21, n2=1.42, nsub=1.45, Lsub=20e-09, resolution=512)

This will build a plano-convex cavity with a specified number of mirrors, indices, and substrate parameters. 
It may be necessary, at the time being, to modify this specific specific method in the cavity.py library,
since this has not yet been fully developed.

To determine the dimensions of needed for the FDTD simulation region, we need to find the minimum and maximum 
z-dimensions of the cavity geometry:

.. code-block:: python
    
    fdtd.select('structure::planar mirror::substrate')
    Region_min = -fdtd.get('z max')
    fdtd.select('structure::featured mirror::substrate')
    Region_max = abs(fdtd.get('z min') - .5*wlen)

Now we can create the FDTD region. Since the cavity geometry in this case is cyllindrically symmetric, we can 
cut the simulation time by exploiting the symmetry, with Symmetric and Anti-Symmetric boundary conditions.

.. code-block:: python
    
    sim = simulation.Simulation(fdtd=fdtd, xy_span_bleed=cav.xy_span_bleed, 
                            runtime=3000e-15, meshacc=3, z_min=Region_min, 
                            z_max=Region_max)

    sim.fdtd_region(x_min_bc="Symmetric", y_min_bc="Anti-Symmetric", z_min_bc="PML", 
                        dt_stab=0.99, fdtd_layers=8, min_layers=8, max_layers=64,
                        autoshutoff=1e-05)

Since the DBR mirror have small features, we need to capture the behaviour here more carefully. We can do this 
with a seperate mesh for each mirror:

.. code-block:: python 

    sim.add_mesh(name='planar mesh', dx=0.01e-06, dy=0.01e-06, dz=0.01e-06, 
                based_on_struct=True, struct="planar mirror")

    sim.add_mesh(name='feat mesh', dx=0.01e-06, dy=0.01e-06, dz=0.01e-06, 
                based_on_struct=True, struct="featured mirror")

We can now add a source in the cavity. In this example we will use a dipole source:

.. code-block:: python
    
    src = source.Source(fdtd=fdtd, wlen=wlen)
    src.dipole(theta=90, shift=0, emission_width=100e-09)

.. code-block:: python
    
    pml_thickness = 0e-09
    xy_span_pml = cav.xy_span + pml_thickness

    mon = monitor.Monitor(fdtd=fdtd)

    apod_center = 0e-15
    apod_start_w = 100e-15

.. code-block:: python
    
    mon.Q_monitor(Qmonitor_zspan=10, Qmonitor_zlayer=1, t_sample=10, dipole_shift=0)


.. code-block:: python
    
    mon.index_monitor(name="n", monitor_type="2D Y-normal",
                        x_min=-.5*xy_span_pml, x_max=.5*xy_span_pml, 
                        z_min=Region_min-.5*pml_thickness, z_max=Region_max + .5*pml_thickness)

.. code-block:: python
    
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

.. code-block:: python
    
    if save:
        fdtd.save(filename)
    
    fdtd.run()

Analysis
~~~~~~~~

.. note::

   This section is still under heavy development. This analysis has been done well in Lumerical's in-built 
   scripting language, but at present, I am porting this over to python in the most natural way I can. As a 
   result, this will take a bit of time to get fully finished.

As before, we start by importing the necessary libraries. this time, we need to import the cavityanalysis class, 
which contains the required methods for calculating Quality Factor, Purcell Factor, Mode Volume, etc.

.. code-block:: python
    
    import numpy as np
    import imp
    # import os
    import scipy.constants as sc
    import cavityanalysis

The FDTD object must be instantiated again:

.. code-block:: python
    
    # os.add_dll_directory("C:\\Program Files\\Lumerical\\v232\\api\\python\\")
    lumapi = imp.load_source("lumapi","C:\\Program Files\\Lumerical\\v232\\api\\python\\lumapi.py")

    fdtd = lumapi.FDTD()


Now, the lumerical simulation file (.fsp) and data file (.ldf) must be loaded. These contain all of the simulation 
structures and variables.

.. code-block:: python
    
    pdir = '< DIRECTORY CONTAINING SIMULATION .FSP FILE >'

    fdtd.load(pdir + 'fsp.fsp')
    fdtd.loaddata(pdir + 'lnp.ldf')

We will need to extract information from the simulation monitors. We can get the monitor names using the FDTD 
object, but since we already know the monitor names (since we defined them), we can make this easier by just assinging 
the monitor variable names directly:

.. code-block:: python
    
    # Monitor names
    m0 = 'n'
    m1 = 'xy_exoplanar'
    m2 = 'xy_middle'
    m3 = 'xy_exofeatured'
    m4 = 'xz_middle'
    m5 = 'xz_edge'
    m6 = 'yz_middle'
    m7 = 'yz_edge'

We can acquire some required variables directly from the simulation data, and use these to define further necessary
variables:

.. code-block:: python
    
    xy_span = fdtd.getnamed("structure::dielectric mediu", "x span")
    wlen = fdtd.getnamed("source", "center wavelength")

    x = fdtd.getdata(m2,"x") 
    y = fdtd.getdata(m2,"y")
    z = fdtd.getdata(m4,"z") 

    x_pts = np.size(x)     # Points in the axis
    z_pts = np.size(z)

    x_res = xy_span/x_pts # resolution of simulation region
    midpoint = np.floor(x_pts/2)+1

    f = np.squeeze(fdtd.getdata(m2,"f")) #Gets frequency data and removes singleton dimensions
    lam = sc.c/f

The first thing we will calculate is the quality factor. to do this, we must first define a range of wavelengths 
over which to search for the peak centre frequency and bandwidth:

.. code-block:: python
    
    ### Define the freq ROI for the peak interested ####

    max_wlen = 638.0*1e-09
    min_wlen = 630.0*1e-09

    max_wlen = fdtd.find(lam, max_wlen)
    min_wlen = fdtd.find(lam, min_wlen)
    w_range_max = (2*np.pi*sc.c)/min_wlen
    w_range_min = (2*np.pi*sc.c)/max_wlen
    f_range_max = sc.c/min_wlen
    f_range_min = sc.c/max_wlen

Now, we can call the cavityanalysis class, and instantiate the QFactor method:

.. code-block:: python
    
    cavana = cavityanalysis.CavityAnalysis(fdtd=fdtd)
    cavana.Qfactor()