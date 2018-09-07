Creation of Seed Images
=======================

A "seed image" is a noiseless image that contains signal only from simulated astronomical sources, created during the first stage of Mirage. In addition to containing no noise, there are no cosmic rays nor dark current or background signal. The seed image can be used as a "truth image" for the purposes of quality checks on the final output data.

Mirage generates seed images in one of two ways:

1. :ref:`From source catalogs <source_catalogs>`
2. :ref:`Extraction from a fits file containing a distortion-free image <mosaic_input>` (e.g. HUDF, GOODS, CANDELS, etc)

.. _source_catalogs:

Seed images from source catalogs
--------------------------------

To create a seed image from user-specified source catalogs, the *catalog_seed_image.py* module is used. There are two types of observations supported by this method, sidereal and non-sidereal, and several catalogs to go with each. The :ref:`Catalogs <catalogs>` page shows examples of each type of source catalog that can be used by Mirage.

For each type of target, Mirage creates or retrieves a small stamp image of the target based on the user-supplied parameters in the catalog. This stamp image is then placed into the seed image at the x,y or RA, Dec location specified by the user.




.. _mosaic_input:

Seed images from distortion-free images
---------------------------------------
