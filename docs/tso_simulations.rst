.. _tso_data:

Simulating TSO data with Mirage
===============================

Mirage can be used to simulate both grism and imaging time series data for NIRCam. Grism Time Series Observations (TSO) are constructed using the *grism_tso_simulator.py* module. Imaging TSO data are produced using the same *catalog_seed_image.py* module that is used for standard imaging observations.

The easiest way to see how to create both types of TSO data is by reading through the `TSO example notebook <https://github.com/spacetelescope/mirage/blob/master/examples/NIRCam_TSO_examples.ipynb>`_.

Yaml files
----------

Generation of the appropriate input yaml files :ref:`from an APT file <from_apt>` is fully supported, as with other observing modes. Note that TSO observations often include a target acquisition (TA) exposure at the beginning of the observation. Mirage will generate a yaml file corresponding to the TA exposure in addition to the TSO exposure. Further, when simulating Grism Time Series observations, the two shortwave detectors with the save field of view as the lower half of the longwave detector will collect imaging time series data of the scene that is dispersed in the long wave detector. Mirage also generates yaml files for these accompanying observations.

User Inputs
-----------
As shown in the `TSO example notebook <https://github.com/spacetelescope/mirage/blob/master/examples/NIRCam_TSO_examples.ipynb>`_, there are several required user inputs that are unique to TSO data. These inputs are listed in the subsections below. See the notebook for details on how to create each of these inputs.

Mirage uses two TSO-specific types of source catalogs. Descriptions of the :ref:`**GrismTSOCatalog** <grism_tso_cat>` and :ref:`**ImagingTSOCatalog** <imaging_tso_cat>` are given on the :ref:`catalogs <catalogs>` page.

Note that Mirage relies on the `Batman <https://www.cfa.harvard.edu/~lkreidberg/batman/>`_ package to generate and work with lightcurves. Many of the input quantities required in the Grism TSO catalog are Batman-specific.

Both Grism and Imaging TSO data
+++++++++++++++++++++++++++++++

- **Background source catalog:** Point source, galaxy, or extended source :ref:`catalog <catalogs>` containing any background sources you wish to include in the simulation.

Grism TSO data
++++++++++++++

- **SED file of TSO target:** This is the spectrum of the unocculted star associated with the TSO object.
- **Transmission spectrum of planet:** Wavelength-dependent effective radius of the planet, in units of the stellar radius.
- **Grism TSO source catalog:** Mirage-specific source catalog containing all information relevant to the source

Imaging TSO data
++++++++++++++++

- **Lightcurve file:** Tabulated list of the parent body's flux versus time
- **Imaging TSO source catalog:** Mirage-specific source catalog containing all information relevant to the source

