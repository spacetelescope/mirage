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
