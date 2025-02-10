Usage
=====

Installation
------------

PyLumerical can be installed directly as-is by downloading using pip:

.. code-block:: console

    $ pip install PyLumerical

At present however, since the project is still under development, it is best to download it from source from 
github.

.. code-block:: console

    $ git clone https://github.com/png-materials-ox/PyLumerical.git

This is a private repository, and so at present only those with permission may download this. When the 
repository has been downloaded, it can be used in one of two ways. First, the software package can be directly 
installed into the python namespace by navigating to the directory and running the setup.py module:

.. code-block:: console

    $ cd PyLumerical
    $ python setup.py install

Now, the package can be imported in your python code by:

.. code-block:: python

    import PyLumerical

Alternatively, the software can be used in the local 'source' folder, with a main calling script. This might be 
helpful when wanting to make some intitial additions to the existing code, e.g. adding extra functionality, but 
is not appropriate for long term use. To do this, navigate to the source directory, and make a main file:

.. code-block:: console 

    $ cd PyLumerical
    $ touch main.py

In the main file, call the necessary libraries, e.g.:

.. code-block:: python 

    import cavity
    import simulation
    import source
    import monitor

To use PyLumerical, a user must make a call the Lumerical API file stored locally in the Lumerical directory. For python distributions up to and including Python 3.11, this can be done using Python's in-built **imp** library:

.. code-block:: python

    import imp
    lumapi = imp.load_source("lumapi","C:\\Program Files\\Lumerical\\v232\\api\\python\\lumapi.py")

However, in Python 3.12 **imp** was deprecated, and for distributions higher than this, it is now necessary to use python's in-built **importlib** library:

**Windows**

.. code-block:: python

    import importlib.util
    spec_win = importlib.util.spec_from_file_location('lumapi', 'C:\\Program Files\\Lumerical\\v242\\api\\python\\lumapi.py')
    lumapi = importlib.util.module_from_spec(spec_win) #windows
    spec_win.loader.exec_module(lumapi)

**Linux**

.. code-block:: python

    import importlib.util
    spec_lin = importlib.util.spec_from_file_location('lumapi', "/opt/lumerical/v242/api/python/lumapi.py")
    lumapi = importlib.util.module_from_spec(spec_lin)
    spec_lin.loader.exec_module(lumapi)

Now with the lumapi object, and FDTD session can be created by:

.. code-block:: python

    fdtd = lumapi.FDTD()