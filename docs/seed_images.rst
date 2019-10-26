.. _seed_images:

Creation of Seed Images
=======================

A "seed image" is a noiseless image that contains signal only from simulated astronomical sources, created during the first stage of Mirage. In addition to containing no noise, there are no cosmic rays nor dark current. The seed image can be used as a "truth image" for the purposes of quality checks on the final output data.

Mirage generates seed images in one of two ways:

1. :ref:`From source catalogs <source_catalogs>`
2. :ref:`Extraction from a fits file containing a distortion-free image <mosaic_input>` (e.g. HUDF, GOODS, CANDELS, etc)

.. _source_catalogs:

Seed images from source catalogs
--------------------------------

To create a seed image from user-specified source catalogs, the *catalog_seed_image.py* module is used. There are two types of observations supported by this method, sidereal and non-sidereal, and several catalogs to go with each. The :ref:`Catalogs <catalogs>` page shows examples of each type of source catalog that can be used by Mirage.

For each type of target, Mirage creates or retrieves a small stamp image of the target based on the user-supplied parameters in the catalog. This stamp image is then placed into the seed image at the x,y or RA, Dec location specified by the user. Source locations on the detector take into account instrument distortion. If source locations in the catalogs are specified using RA and Dec, `Mirage` will use the distortion model from a distortion reference file if supplied in the input yaml file. If no reference file is given, `Mirage` will fall back to using distortion information in the Science Instrument Aperture File (`SIAF <https://jwst-docs.stsci.edu/display/JTI/JWST+Field+of+View>`_).

For point sources, `Mirage` selects the appropriate PSF for each source from a PSF library. This PSF is scaled to the requested brightness, and then placed in the correct location on the detector.

Galaxies in `Mirage` are represented by 2 dimensional Seric profiles. These are created using astropy's `Sersic2D model <http://docs.astropy.org/en/stable/api/astropy.modeling.functional_models.Sersic2D.html>`_, scaled to the requested brightness, and placed in the seed image.

"External sources" is the term used for stamp images that are read in from supplied fits files. In this case, the user-suppied fits file is read in and the image is scaled and placed in the seed image at the requested location. See the :ref:`Catalogs <catalogs>` page for more details.

After simulated astronomical sources are added to the seed image, background signals can be added using the `bkgrnd` entry in the input yaml file. Users can request a constant background in units of ADU/pixel/second, or a general background level. For the latter, `Mirage` uses the `jwst_backgrounds <https://github.com/spacetelescope/jwst_backgrounds>`_ package to add zodiacal and thermal background signal. In this case, "low", "medium", or "high" are the three values that can be used. The definition of these terms mirrors those used in `APT <http://www.stsci.edu/hst/proposing/apt>`_ and the `JWST ETC <https://jwst.etc.stsci.edu/>`_.

For a given pointing, `jwst_backgrounds` will calculate the total expected background level (zodiacal plus thermal) for each day over the course of a year. A "low" background is defined as the background level corresponding to the 10th percentile of this distribution of background levels. "Medium" is 50th percentile, and "high" is 90th percentile.

Users can also provide a custom background through the use of an input fits file referenced in the "zodiacal" entry of the input yaml file. In this case, the contents of the file are read in and added to the seed image.

Finally, `Mirage` uses the pixel area map to mimic the effects of distortion on the simulated sources. While the locations of the simulated sources within the seed image accurately account for distortion, the brightness of the sources are not affected. The seed image (minus any point sources) is multiplied by the pixel area map in order to adjust the extended sources’ brightness to account for the relative differences in pixel areas across the detector.

.. _mosaic_input:

Seed images from distortion-free images
---------------------------------------

Another way to create a seed image is through the use of a large field-of-view image from a FITS file. This file should contain a north-up, distortion-free count rate image at an arbitrary pixel scale with a proper world coordinate system (WCS) in the header.

In this case, the portion of the image corresponding to the requested RA, Dec of the simulation is extracted from the input, and some of the JWST calibration pipeline functionality is used to resample it onto the appropriate JWST instrument pixel grid while also introducing the appropriate distortion for the requested detector. This functionality is a new implementation of Astrodrizzle’s `blotting function <https://drizzlepac.readthedocs.io/en/deployment/ablot.html>`_, and works in essentially the same way. For more details on Astrodrizzle, see the `DrizzlePac website <http://drizzlepac.stsci.edu/>`_. The simulator then saves this blotted image in the seed image format used by subsequent steps of MIRAGE.

In order to create a seed image from this input file, the simulator uses the **fits_seed_image.py**, function. For convenience, this function accepts the same yaml input file as the other parts of the simulator. To see an example of how to create a seed image from an input mosaic file, see the `Simulated_data_from_mosaic_image.ipynb notebook <https://github.com/spacetelescope/mirage/blob/master/examples/Simulated_data_from_mosaic_image.ipynb>`_ in the `examples directory <https://github.com/spacetelescope/mirage/tree/master/examples>`_ of the Mirage repository.



