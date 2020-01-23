1.3
===

Installation
------------

setup.py has been modified to support installation via pip and Pypi. Installation documentation has been updated to describe the new process.


Gain Values
-----------

Update observation_generator.py, wfss_simulator.py, grism_tso_simualtor.py to use the mean gain value stored in utils/constants.py rather than the values in the gain reference file when translating the dispersed seed image from units of e-/sec to ADU/sec.

Flat Field
----------

Seed images are now multiplied by the flat field reference file rather than the pixel area map reference file in order to get the surface brightnesses correct. Or more simply, since one of the steps in the JWST calibration pipeline is to divide by the flat field, we must multiply by the flat field when creating the data. See #430. For imaging/time series modes, the flat field is multiplied into the seed image. For NIRCam WFSS mode, the flat field is multiplied in to the dispsersed seed image. For NIRISS WFSS, the flat field is multiplied in to the seed image prior to dispersing. This is because the flat field reference file contains both pixel-to-pixel differences in respsonse, as well as images of the occulting spots, which are in the optical train. Ideally the occulting spots would be multiplied into the seed image prior to dispersing, and then the pixel-to-pixel flat would be multiplied into the seed image after dispersing. Unfortunately these two effects are mixed in the flat field reference file and cannot be separated. This will have some implications for calibrated data products.


Backgrounds
-----------

Update the calculation of background signals to better match the values calcualted by the ETC. Values generally are within 10% of those from the ETC, although there are some filters/pointings/levels where the values differ by up to 20%. #430

For the purposes of calculating the background signal in NIRISS WFSS simulations, the system throughput is simply set to 80% of the throughput in the imaging mode (with appropriate filter).

Fixed bug in the Grism TSO simulator where the background signal was being added twice.

Changed the code so that the Grism TSO simulator works in the case where no background source catalogs are provided. In this case, a dummy background point source catalog is generated, as the calculation and addition of background is done using the background sources.

calculate_background function moved into backgrounds.py so that it can be more easily used by modules other than catalog_seed_image.

Besancon Model Query
--------------------

Code relating to the production of Besancon model source catalogs has been updated to reflect the new workflow for querying and retrieving data. This is now a 2-step process. Users must create an account on the Besancon model website. Queries can then be submitting using the `catalogs.create_catalog.besancon` function. The user must then wait for an email which contains a link to download the resulting catalog. Conversion of this catalog to Mirage-format can then proceed. See the `Catalog_Generation_Tools.ipynb` notebook for details.

Non-sidereal
------------

Segmentation map addition bug corrected. Example notebook and input yaml file updated.


1.2.2
=====

Versioning
----------

Update package versioning to be done with setuptools-scm rather than relic.


1.2.1
=====

TSO Modes
---------

- Updated documentation on readthedocs with information on TSO mode work


1.2
===

TSO Modes
---------

- Add the ability to simulate both grism and imaging time series observations for NIRCam. Example notebook included.


1.1.5
=====

PSF Selection
-------------

- Fix bug in PSF library selection code for observations using one of NIRCam's filters present in the pupil wheel. The bug was preventing the correct library file from being found. (#420)


1.1.4
=====

WCS keywords
------------

- Correct the input RA and Dec used to calculate the values of the PC matrix. Remove the calculation of CRVAL1,2 from set_telescope_pointing.py since it is already done in observation_generator.py (#419)


1.1.3
=====

Yaml Generator
--------------

- Update generator to produce yaml files only for the detectors used with a given aperture. e.g. SUB400P with the NIRCam B module only uses NIRCam B1 and B5 detectors. With this update,
yaml files will only be produced for B1 and B5, whereas previously yaml files were generated for all 5 B module detectors. This change only affects NIRCam.


1.1.2
=====

WFSS
----

- Update functionality for rescaling input spectra to desired magnitude in given instrument/filter. Rescaling is now done via synphot's renormalize() function in the prpoper photon-weighted units. (#412)

Catalogs
--------

- Change photometric system in catalog output from 2MASS query from ABmag to Vegamag (#415)

Seed Image
----------

- Remove filter substring from seed image output file name in the case of FGS simulations (#415)


1.1.1
=====

WFSS
----

- Update background scaling calcultions. NIRISS scales pre-existing background image. NIRCam creates image from jwst_background-provided date or level [#399]
