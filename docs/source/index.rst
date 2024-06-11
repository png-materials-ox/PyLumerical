.. PyLumerical documentation master file, created by
   sphinx-quickstart on Thu May  9 13:49:19 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PyLumerical's documentation!
=======================================

**PyLumerical** is a Python library for generating structures for Ansys's Lumerical FDTD software, and for 
performing subsequent post-process analysis. The primary focus of this libray is to investigate the opticsal 
response of microcavities; in particuar Fabry-Perot structures. However, the library can be extended to a range 
of optical structures.  

This library has been developed at Oxford University's Photonic Nanomaterials Group (PNG) 
(https://www-png.materials.ox.ac.uk) by Gareth Siôn Jones. It provides an interface to Lumerical's FDTD software, 
and at times uses it's in-built analysis objects, for example, for calculating Q-Factor. This library also includes 
bespoke analysis objects for calculating mode volume and Purcell enhancement.


.. note::

   This project is under active development.

.. toctree::
   :maxdepth: 4
   :caption: Contents:

   usage     
   
   examples
   
   theory

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
