.. _ghosts:

Addition of Ghosts to Mirage Simulations
========================================

NIRCam
++++++

The addition of ghosts is not currently supported.


FGS
+++

The addition of ghosts is not currently supported.


NIRISS
++++++

Optical ghosts in NIRISS have been observed and characterized in ground testing data. A detailed description of the ghosts is given in the `Ghosts section <https://jwst-docs.stsci.edu/near-infrared-imager-and-slitless-spectrograph/niriss-instrumentation/niriss-gr150-grisms#NIRISSGR150Grisms-Ghosts>`_ of the NIRISS GR150 Grisms Jdox page.


Controlling the Addition of Ghosts
----------------------------------

Optical ghosts can be added to simulated NIRISS imaging and WFSS observations. Their addition is optional, and is controlled by the **sim_Signals:add_ghosts** and **sim_Signals:PSFConvolveGhosts** keywords in the yaml input files.

If **add_ghosts** is True, Mirage will calculate, for each point source in your source catalog, the location and magnitude of the optical ghost associated with the optical elements being used. Within Mirage, once the ghost locations and magnitudes are known, the ghosts are treated as :ref:`extended sources <source_catalogs>`, which rely upon fits files containing stamp images. The stamp images are read in, scaled to the requested brightness, and added to the seed image. The location and brightness of the ghost are calculated using the "gap summary" file in the **config** directory of the Mirage repository. This file specifies the observed offsets, as well as count rate of a ghost for a given filter, pupil, and source location and count rate on the detector.

.. tip::

    The current "gap summary" file allows for the addition of ghost sources for observations that use F090W, F115W, F140M, F150W, and F200W. For other filters, the ghost positions and/or magnitudes have not been characterized, and therefore ghosts will not be added. If you wish to run Mirage with a "gap summary" file other than that in the config directory, you can either replace the file in the config directory with your own, or update the value of NIRISS_GHOST_GAP_FILE in the constants.py file of the repository to point to your file. This may be useful for testing new ghost position calibration results.

Mirage currently carries one fits file in its collection of :ref:`reference files <reference_files>` that provides an image of an optical ghost from a point source on the NIRISS detector. This image was created from ground testing data. By default, this file is used when adding ghosts to the data. **NOTE that this stamp image applies to point sources only. By default, no ghosts are added for galaxy or extended sources.**

Users may add ghosts associated with galaxy or extended sources (or override the default stamp image for point sources) by supplying their own stamp images and specifying which stamp images are to be associated with each source within the source catalogs. An example of this is shown in the :ref:`Specifying Ghost Stamp Images <specifying_ghost_stamp_images>` section below.

If the **PSFConvolveGhosts** keyword is True, then the stamp image will be convolved with the instrumental PSF prior to adding the ghost image to the data. By default, this is set to False. Users can control its value, along with that of **add_ghosts**, in the call to the :ref:`yaml_generator <add_ghosts>`.

Specifying Ghost Stamp Images
-----------------------------

In this section we give a brief example of creating a source catalog containing a list of stamp images to be used in the creation of ghost sources. By default, a source catalog will not be created with a column for these stamp images. In such a case, Mirage will fall back to using the default stamp image within the collection of reference files for point sources only. The example below creates a source catalog for extended sources with ghost stamp images. This will allow Mirage to create a ghost source for each extended source. The *niriss_ghost_stamp* keyword can enable this behavior when creating any type of source catalog (e.g. PointSourceCatalog, GalaxyCatalog, MovingPointSourceCatalog, etc).

In this example we have a catalog with three sources. Stamp images of the sources themselves are contained in the files *source1.fits*, *source2.fits*, and *source3.fits*. Stamp images of the corresponding ghosts are in *ghost1.fits*, *ghost2.fits*, and *ghost3.fits*.

::

	from mirage.catalogs import catalog_generator as cg

	ra = [12.2, 12.2001, 12.2002]
	dec = [23.1, 23.1001, 23.1002]
	source_stamp_files = ['source1.fits', 'source2.fits', 'source3.fits']
	position_angle = [0., 0., 0.]
	magnitudes = [15.5, 16.3, 17.4]
	ghost_stamp_files = ['ghost1.fits', 'ghost2.fits', 'ghost3.fits']

	cat = cg.ExtendedCatalog(ra=ra, dec=dec, filenames=source_stamp_files,
							 position_angle=position_angle,
							 niriss_ghost_stamp=ghost_stamp_files)
	cat.add_magnitude_column(magnitudes, instrument='niriss', filter_name='f380m')
	cat.save('extended_sources_with_ghosts.cat')

Display the table:

::

	cat.table

	index x_or_RA y_or_Dec   filename   niriss_f380m_magnitude niriss_ghost_stamp pos_angle
    int64 float64 float64     str12            float64               str11         float64
    ----- ------- -------- ------------ ---------------------- ------------------ ---------
    1    12.2     23.1 source1.fits                   15.5        ghost1.fits       0.0
    2 12.2001  23.1001 source2.fits                   16.3        ghost2.fits       0.0
    3 12.2002  23.1002 source3.fits                   17.4        ghost3.fits       0.0
