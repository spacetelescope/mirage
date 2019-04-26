.. _wfss_data:

Simulating WFSS data with Mirage
================================

Mirage can be used to simulate Wide Field Slitless Spectroscopy (WFSS) mode data for NIRCam and NIRISS, using the *wfss_simulator.py* module. To produce these simulations, *Mirage* constructs one or more imaging mode :ref:`seed images <seed_images>`, along with an associated segmentation map. These seed images and segmentation map are then passed to the disperser software, which is in the `NIRCAM_Gsim <https://github.com/npirzkal/NIRCAM_Gsim>`_ package. The disperser then takes the seed images, along with an optional input file containing object spectra, and disperses the signal across the detector in the same manner as the `NIRCam <https://jwst-docs.stsci.edu/display/JTI/NIRCam+Grisms>`_ and `NIRISS grisms <https://jwst-docs.stsci.edu/display/JTI/NIRISS+Wide+Field+Slitless+Spectroscopy>`_.

This mode generally uses more input compared to the creation of imaging mode data, as spectral information on the input sources may be provided. There are two methods that spectral information



yaml file: must specify wfss mode and grism_source_image = True. The appropriate grism must be specified, as well as a crossing filter



Inputs
------

There are three types of inputs that can be used to create WFSS data. The first is the same :ref:`yaml parameter file <example_yaml>` that is used when creating imaging mode data. Along with the yaml files, the appropriate :ref:`ascii source catalogs <catalogs>` must be provided. The third input, which is optional, is an `hdf5 <https://www.h5py.org/>`_ file that contains the spectra for some or all of the targets that are in the source catalogs. Below we describe how to use these inputs to create WFSS data.

Note that when using yaml files as inputs, at least one of these files must have the :ref:`mode <mode>` set to "wfss", :ref:`grism_source_image <grism_source_image>` set to True, and a grism placed in the appropriate :ref:`filter <filter>` or :ref:`pupil <pupil>` entry. This will be done automatically for the appropriate files during the yaml generation process if :ref:`generating yaml files from an APT proposal <from_apt>`.


.. _single_yaml:

Single yaml file
++++++++++++++++

In the simplest case a single :ref:`yaml parameter file <example_yaml>` is provided to the WFSS simulator module, along with a source catalog containing target magnitudes in a single filter. In this case, Mirage converts the provided magnitude values to flux densities, and the disperser assumes a flat continuum spanning the entire wavelength range of the filter/grism combination.

If the source catalog contains magnitudes in multiple filters, *Mirage* will, for each source, linearly interpolate the source magnitudes in order to construct a continuum spectrum. If the provided magnitudes do not cover the entire wavelength range necessary for the dispersion, then *Mirage* will optionally extrapolate the continuum spectrum to cover the full wavelength range.

.. tip::
    In this case where only ascii source catalogs are provided, source spectra will not contain any features (emission, absorption), but will rather be smooth continuum spectra. In order to simulate sources with emission or absorption features, this information must be added via the :ref:`hdf5 file <yaml_plus_sed>` described below.


Multiple yaml files
+++++++++++++++++++

Another way to produce data with smooth continuum spectra is to provide multiple yaml files, where each yaml file will produce a seed image through a different filter. In this case, the mutiple seed images will be used to calculate source flux densities, rather than these calculations being done using the source catalogs as input. One of these yaml files must specify WFSS mode with the requested grism, as described above. The other yaml files should specify imaging mode. As this method produces the same output as that when a :ref:`single yaml file <single_yaml>` and source catalog with multiple magnitude columns is provided, but with more calculations taking a longer time, we recommend against using this strategy.

.. _yaml_plus_sed:

Yaml file plus SED file
+++++++++++++++++++++++

In order to create simulated data with more realistic spectra, users can provide an optional `hdf5 <https://www.h5py.org/>`_ file that contains spectra for some or all of the targets listed in the source catalogs. The spectra in this file can have any shape, and can include emission and absorption features. The spectrum for each source is contained in a "dataset", which is read in as a python dictionary containing "wavelengths" and "fluxes" keys. The values associated with each of these fields is a list of floating point numbers. Units can be specified by adding a string as a dataset attribute. If no units are provided, Mirage assumes that wavelengths are in units of microns, and flux densisites are in units of F_lambda (erg/sec/cm:sup:`2` /A). *Mirage* contains functions that can be used to `create hdf5 files <https://github.com/spacetelescope/mirage/blob/master/mirage/catalogs/hdf5_catalog.py>`_ with their target SEDs.

.. tip::

    The `NIRISS WFSS example notebook <https://github.com/spacetelescope/mirage/blob/master/examples/NIRISS_WFSS_data_creation_example.ipynb>`_ shows examples of how to create your own hdf5 catalog file, and how to create WFSS data using the methods described above.



