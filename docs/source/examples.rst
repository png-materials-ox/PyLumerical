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
    import os
    from PyLumerical import cavity, simulation, source, monitor, analysisobjects

'imp' is a standard python library which allows us to import a local file. We need to link the main script to the 
Lumerical python API file which is stored locally on our PC within the Lumerical folder. For example, the lumapi
may be found at the following address (depending on distribution).

.. code-block:: python
    
    # os.add_dll_directory("C:\\Program Files\\Lumerical\\v232\\api\\python\\")
    lumapi = imp.load_source("lumapi","C:\\Program Files\\Lumerical\\v232\\api\\python\\lumapi.py")

Above, the **os.add_dll_directory** is only used if there is any difficulty adding the file.

We can choose whether or not to save the file, and asign a filename if necessary:

.. code-block:: python
    
    save = True
    filename = "DIRECTORY_NAME"


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

    sim.fdtd_region(x_min_bc="Anti-Symmetric", y_min_bc="Symmetric", z_min_bc="PML", 
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

We may want to add a number of different monitors. For example, we can add a refractive index monitor 
and a power monitor. To set the geometry, we first need to determine the xy span. Then, we may wish to 
set the apodisation. This is a cutoff for when to start / stop collecting field data. We may need this if, for 
example, we want to remove the impact of the initial emitter pulse from our measurements. We set the apodisation 
such that we start collecting field data after the pulse has rung down, and with a filter window width (see 
Lumericals documentation for more information).

.. code-block:: python
    
    pml_thickness = 0e-09
    xy_span_pml = cav.xy_span + pml_thickness

    mon = monitor.Monitor(fdtd=fdtd)

    apod_center = 100e-15
    apod_start_w = 10e-15
    
    mon.index_monitor(name="n", monitor_type="2D Y-normal",
                        x_min=-.5*xy_span_pml, x_max=.5*xy_span_pml, 
                        z_min=Region_min-.5*pml_thickness, z_max=Region_max + .5*pml_thickness)
    
    mon.power_monitor(name='xy_power', montype="2D Z-normal", plane="xy",
                      x_min=-.5*xy_span_pml, x_max=.5*xy_span_pml, 
                      y_min=-.5*xy_span_pml, y_max=.5*xy_span_pml, 
                      z=Region_max,
                      apod="Start", apod_center=apod_center, apod_time_width=apod_start_w)

For cavity simulations, we usually want to calculate the Quality factor, mode volume, purcell factor, and 
perhaps collection efficiency. We can start by adding special structures known as 'analysis groups'. These are 
created in the AnalysisObjects class. Here, we will look at farfield collection efficiency, quality factor,
and mode volume. Once the simulation has run, we will run these analyses and extract the desired data.

.. code-block:: python 

    objects = analysisobjects.AnalysisObjects(fdtd=fdtd)

    z_span = abs(Region_min) + abs(Region_max)

    objects.farfield(xy_span=xy_span, z_span=z_span, x=0, y=0, z=0.55e-06, theta_max=90,
                    N_theta=180, Nphi=91)

    objects.Qfactor()

    objects.mode_volume_3D(xy_span=xy_span, z_span=z_span, x=0, y=0, z=0.55e-06)

Before running any simulation, it is generally a good idea to check if the index data makes sense visually.
So, it can be useful to temporarily change the boundary conditions back to all PML by uncommenting the associated
variables, and commenting out the **fdtd.run()** line. We also want to save our simulation file to a well named 
directory.

.. code-block:: python
    
    #fdtd.select('FDTD')
    #fdtd.set('x min bc','PML')
    #fdtd.set('y min bc','PML')

    if save:
        if not os.path.exists(filename):
            os.mkdir(filename)
            
        else:
            try:
                fdtd.save(filename + 'fsp')
            except:
                FileExistsError("Couldn't write file")

Analysis
~~~~~~~~

.. note::

   This section is still under heavy development. This analysis has been done well in Lumerical's in-built 
   scripting language, but at present, I am porting this over to python in the most natural way I can. As a 
   result, this will take a bit of time to get fully finished.

This package contains bespoke custom scripts alongside invokations of Lumerical's in-built analysis groups. 
For the present example, we combine both custom and in-built methods. Please be aware that the 'imp' library has been deprecated for modern python versions, so you may need to use a different library. For more details, see the 'usage' page.

.. code-block:: python

    import numpy as np
    import imp
    import scipy.constants as sc
    from PyLumerical import cavityanalysis, lumericalanalysis
    import matplotlib.pyplot as plt

Start by loading the lumerical simulation file which has previously been run and contains simulation data.

.. code-block:: python

    # os.add_dll_directory("C:\\Program Files\\Lumerical\\v232\\api\\python\\")
    lumapi = imp.load_source("lumapi","C:\\Program Files\\Lumerical\\v232\\api\\python\\lumapi.py")

    fdtd = lumapi.FDTD()

    pdir = "<DIRECTORY CONTAINING SIMULATION FILE WITH DATA>"

    fdtd.load(pdir + 'fsp.fsp')

The Lumerical analysis class can now be instantiated. This class contains the methods necessary for running and analysing the quality factor, purcell factor, mode volume, 
collection efficiency, etc, of the microcavity.

.. code-block:: python

    analysis = lumericalanalysis.LumericalAnalysis(fdtd=fdtd)

The resonant frequencies, quality factors, etc of each peak:

.. code-block:: python 

    resonances = analysis.resonances()
    w0 = 2*np.pi*resonances['f0'].values
    Q = resonances['Q'].values

    eigenfreqs = analysis.complex_eigenfrequency(omega0=w0, Q=Q) 
    f_spectrum = analysis.f_spectrum(plotting=True, saveplot=False)

To get the electric field profiles:

.. code-block:: python

    Efields = analysis.electric_field_magnitude()
    farfield = analysis.farfield_analysis(NA=0.9)