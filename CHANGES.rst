Not Yet Tagged
==============

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
