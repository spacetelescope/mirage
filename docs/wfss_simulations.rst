.. _wfss_data:

Simulating WFSS data with Mirage
================================

Mirage can be used to simulate Wide Field Slitless Spectroscopy (WFSS) mode data for NIRCam and NIRISS, using the wfss_simulator.py module. To produce these simulations, Mirage constructs one or more imaging mode :ref:`seed images <seed_images>`, along with an associated segmentation map. These seed images and segmentation map are then passed to the disperser software, which is in the NIRCAM_Gsim<ADD LINK HERE> package. The disperser then takes the seed images, along with an optional input file containing object spectra, and disperses the signal across the detector in the same manner as the NIRCam and NIRISS grisms.<ADD LINK TO GRISM PAGES>

This mode generally uses more input compared to the creation of imaging mode data, as spectral information on the input sources may be provided. There are two methods that spectral information



yaml file: must specify wfss mode and grism_source_image = True. The appropriate grism must be specified, as well as a crossing filter



Inputs
------

There are three types of inputs that can be used to create WFSS data. The first is the same :ref:`yaml parameter file <input_yaml_file_parameters>`that is used when creating imaging mode data. Along with the yaml files, the appropriate :ref:`ascii source catalogs <catalogs>` must be provided. The third input, which is optional, is an hdf5<ADD LINK HERE> file that contains the spectra for some or all of the targets that are in the source catalogs. Below we describe how to use these inputs to create WFSS data.


.. _single_yaml:

Single yaml file
++++++++++++++++

In the simplest case a single :ref:`yaml parameter file <input_yaml_file_parameters> is provided, along with a source catalog containing target magnitudes in a single filter. In this case, Mirage converts the provided magnitude values to flux densities, and the disperser assumes a flat continuum spanning the entire wavelength range of interest.

If the source catalog contains magnitudes in multiple filters, Mirage will, for each source, linearly interpolate the source magnitudes in order to construct a continuum spectrum. If the provided magnitudes do not cover the entire wavelength range necessary for the dispersion, then Mirage will optionally extrapolate the continuum spectrum to cover the full wavelength range. (FOR THIS TO WORK YOU NEED TO CALL CATALOGS_TO_HDF5. INTEGRATE THIS WITH WFSS_SIMULATOR SO THAT IT IS INVISIBLE TO THE USER)

.. tip::
    In this case where only ascii source catalogs are provided, source spectra will not contain any features (emission, absorption), but will rather be smooth continuum spectra. In order to simulate sources with emission or absorption features, this information must be added via the hdf5 file described below.


Multiple yaml files
+++++++++++++++++++

Another way to produce data with smooth continuum spectra is to provide multiple yaml files, where each yaml file will produce a seed image through a different filter. In this case, the mutiple seed images will be used to calculate source flux densities, rather than these calculations being done using the source catalogs as input.

THIS MODE WILL BE RARE NOW THAT WE HAVE THE ABILITY TO CALCULATE THIS INFO FROM A SOURCE CATALOG


Yaml file plus SED file
+++++++++++++++++++++++

In order to create simulated data with more realistic spectra, users can provide an optional hdf5 file that contains spectra for some or all of the targets listed in the source catalogs. The spectra in this file can have any shape, including emission and absorption features. The spectrum for each source in this file is contained in a "dataset". A dataset is a dictionary containing "wavelengths" and "fluxes" keys. The values associated with each of these fields is a list of floating point numbers. Units can be specified by adding a string as a dataset attribute. If no units are provided, Mirage assumes that wavelengths are in units of microns, and flux densisites are in units of F_lambda (erg/sec/cm^2/A).


POINT TO EXAMPLE OF HOW TO CREATE THIS FILE. DESCRIBE WHAT IS IN THIS FILE


Input Catalogs
++++++++++++++

ascii source catalog

SED file


Create the data
---------------

Use wfss_simulator