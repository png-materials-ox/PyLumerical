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