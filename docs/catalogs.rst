.. _catalogs:

Source Catalog File Formats
===========================

Mirage accepts 7 different types of ascii source catalogs, all of which can be generated using `Mirage's` :ref:`catalog generation <catalog_generation>` functionality. These include:

1. :ref:`Point sources <point_source>`
2. :ref:`Galaxies (2d Sersic profiles) <galaxies>`
3. :ref:`Extended targets <extended_obj>`
4. :ref:`Non-sidereal source <nonsidereal>`
5. :ref:`Moving point sources <moving_point_source>`
6. :ref:`Moving Sersic sources <moving_sersic>`
7. :ref:`Moving extended sources <moving_extended>`
8. :ref:`Grism time series sources <grism_tso_cat>`
9. :ref:`Imaging time series sources <imaging_tso_cat>`

.. tip::
    For information on the optional hdf5 catalog that can be used when simulating WFSS data, see the :ref:`WFSS simulation page<wfss_data>`.

The type of source catalog(s) you need depends on the scene you are simulating. The table below shows the possible catalogs for sidereal and non-sidereal observations. For sidereal observations, the point source, galaxy, and extended source catalogs can be used to create static sources in the scene. The three types of moving target catalogs can be used to create sources that trail across the field of view during the observation. For non-sidereal observations, the non-sidereal source catalog is used to create a static source in the scene. This is the non-sidereal target which JWST is tracking. The other types of catalogs can then all be used to create trailed sources. In this case, the inverse of the non-sidereal source's velocity is applied to all targets in the other catalogs. In this way, Mirage will create trailed background stars and galaxies, as well as trailed non-sidereal targets that have velocities other than that of the main target.

+---------------------+------------------------+----------------------------------+
| Type of observation | Static source catalogs | Trailed (moving) source catalogs |
+=====================+========================+==================================+
|      Sidereal       |  - Point source        |       - Moving point source      |
|                     |  - Galaxies            |       - Moving Sersic sources    |
|                     |  - Extended targets    |       - Moving extended sources  |
+---------------------+------------------------+----------------------------------+
|    Non-sidereal     |  - Non-sidereal source |       - Point source             |
|                     |                        |       - Galaxies                 |
|                     |                        |       - Extended targets         |
|                     |                        |       - Moving point source      |
|                     |                        |       - Moving Sersic sources    |
|                     |                        |       - Moving extended sources  |
+---------------------+------------------------+----------------------------------+
|     Grism TSO       |  - Point source        |       - Not supported            |
|                     |  - Galaxies            |                                  |
|                     |  - Extended targets    |                                  |
|                     |  - Grism Time Series   |                                  |
+---------------------+------------------------+----------------------------------+
|    Imaging TSO      |  - Point source        |       - Not supported            |
|                     |  - Galaxies            |                                  |
|                     |  - Extended targets    |                                  |
|                     |  - Imaging Time Series |                                  |
+---------------------+------------------------+----------------------------------+


Common formatting details
-------------------------
Mirage scans the top 4 lines of each catalog for certain phrases that can be used to characterize the inputs. These phrases are used to specify the units of source locations or velocities, as well as the magnitude system to use. Multiple phrases can be used in a single catalog, but only one phrase per line is allowed. Examples are shown in the catalogs below.

The locations of sources can be specified in RA, Dec or in (x,y) pixel locations on the detector. If you wish to provide positions in units of (x, y) detector pixels, then the string ‘position_pixels’ must be added after the # in one of the top 4 lines of the file.

.. tip::

    RA and Dec values can be given in decimal degrees, colon-separated values (HH:MM:SS and DD:MM:SS), or in more conventional string formats, **but all sources in a given catalog must have the same format**:

    53.08864      -27.83999

    03:32:21.273  -27:50:23.983

    03h32m21.273s -27d50m23.983s


Mirage uses AB magnitudes as the default for input sources. However, you can change the magnitude system by specifying an alternative in one of the top 4 lines. The three acceptible options are **vegamag**, **stmag**, and **abmag**. **All sources in a given catalog must be in the same magnitude system.**

For moving targets (both those that are moving across the field of view, as well as non-sidereal targets), the default unit for velocity is arcseconds per hour. If you wish to instead use pixels per hour, then **velocity_pixels** must be added to one of the 4 top lines of the catalog.

.. _point_source:

Point Sources
-------------
Point sources are specified using a catalog that includes the locations of the sources in RA and Dec (or x,y pixel locations on the detector) and the corresponding magnitudes through the filter specified by the user. Currently the simulator supports the use of ABMAG [Oke, 1983]_, STMAG [Stone, 1996]_ , and VEGAMAG () systems, with ABMAG as the default.

An example point source catalog is shown below with the positions given in RA and Dec.

::

	#
	#
	#
	# abmag
	#
	# Magnitudes are converted from input flux densities.
	x_or_RA          y_or_Dec      nircam_f200w_magnitude
	53.0886395   -27.8399952              20.0
	53.0985009   -27.8398137              19.2

Mirage looks for the exact column names shown above when reading point source catalogs. Changing these column names will cause the simulator to fail when attempting to read in the file.


.. [Oke, 1983] `ApJ 266, 713 <https://ui.adsabs.harvard.edu/#abs/1983ApJ...266..713O/abstract>`_
.. [Stone, 1996] `ApJS 107, 423 <https://ui.adsabs.harvard.edu/#abs/1996ApJS..107..423S/abstract>`_

.. _galaxies:

Galaxies (aka 2D Sersic profiles)
---------------------------------

Below is an example of a galaxy source catalog. In this case, galaxy positions are given in RA, Dec decimal degrees, and the half light radii are in units of arcseconds. The half light radius can also be specified in units of pixels. In that case, you must add **radius_pixels** after the # in one of the top four lines.

Position angle values are defined as degrees east of north.

The simulator software looks for the exact column names shown below when reading these catalogs. Changing these column names will cause the simulator to fail when attempting to read in the file.

::

	#
	#
	#
	# abmag
	#
	# Magnitudes are converted from input flux densities.
	x_or_RA         y_or_Dec     radius    ellipticity    pos_angle       sersic_index      niriss_f200w_magnitude
	53.05           -27.83        0.17        0.46         104.35              3.3                 18.06
	53.10           -27.83        0.73        0.01         195.50              2.7                 16.86

.. _extended_obj:

Extended Objects
----------------

The extended object catalog lists files containing stamp images to be added to the seed image. For example, a source such as a nebula or spiral galaxy that cannot be simulated via a 2-dimensional Sersic profile can be added by placing an image of the source in a fits file. `Mirage` will then read in, scale, and add this image to the seed image.

It is assumed that the fits file contains an array in the 0th or 1st extension. The array can be any size. If it is larger than the field of view of the simulated data, then it is cropped by placing the center of the extended stamp image at the specified x,y or RA, Dec location on the detector, and cropping any areas that fall outside of the detector.

Each row of this catalog contains the name of a FITS file containing the image to use, along with the RA, Dec (or x,y) position of the source, the position angle to use, and the source’s magnitude. The position angle is the angle in degrees East of North of the stamp image. For example, if you have an image of a spiral galaxy extracted from a dataset where the angle from North to the y-axis of the array is 330 degrees, then 330 should go into the pos_angle column of the catalog (assuming you wish to keep the real orientation of the object in the simulated data). Mirange will combine this pos_angle value with the roll angle of the telescope in the simulated scene to properly rotate the stamp image. In the case where you have a stamp image of a generic spiral galaxy that you would like to have at various position angles in the simulated data, then you can modify the pos_angle value appropriately.

For stamp images where it may not make sense to specify a magnitude (such as a galaxy cluster), it is possible to specify ‘None’ as the magnitude. In this case the code assumes that the data contained in the fits file is in units of ADU per second, and will not rescale the data before adding to the seed image. However, the user can also adjust the signal rate of all extended sources through the use of the extendedScale field in the input yaml file. This is a multiplicative factor to apply to the data in the fits file prior to adding the source to the seed image.

::

	#
	#
	#
	#
	# Columns 1 and 2 can be either x,y positions on the detector aperture (e.g.
	# 0,0 is lower left corner of the full frame of the subarray used for the
	# output) or RA,Dec location of the center of the source. If they are x,y
	# positions, make the top line of the file '# position_pixel'
	#
	#
	#
	x_or_RA        y_or_Dec       pos_angle      nircam_f200w_magnitude       filename
	359.65          0.0006           20                 16.000             ring_nebula.fits


.. _nonsidereal:

Non-sidereal Source
-------------------

This catalog is used when creating non-sidereal simulated exposures. In this case, all targets other than the non-sidereal target will then trail through the field of view during the observation. This mode is meant to simulate observations of solar system targets with non-sidereal velocities. This catalog should contain only one entry, with RA, Dec or x, y position, as well as velocity values (arcsec/hour or pixels/hour) and object magnitude.

::

	#
	#
	#
	# abmag
	#
	# radius can also be in units of pixels or arcseconds. Put 'radius_pixels' at top of file
	# to specify radii in pixels.
	# position angle is given in degrees counterclockwise.
	# An "object" value containing 'point' will be interpreted as a point source.
	# Anything containing "sersic" will create a 2D sersic profile.
	# Any other value will be interpreted as an extended source.
	# x_or_RA_velocity is the proper motion of the target in units of arcsec (or pixels) per hour
	# Y_or_Dec_velocity is the proper motion of the target in units of arcsec (or pixels) per hour
	# if the units are pixels per hour, include 'velocity pixels' in line 2 above.
	object       x_or_RA    y_or_Dec   x_or_RA_velocity    y_or_Dec_velocity     nircam_f200w_magnitude
	pointSource  53.101      -27.801       2103840.              0.0                       17.

.. _moving_point_source:

Moving Point Sources
--------------------

The moving point source catalog contains a list of point sources to move through the field of view during the integration. Similar to the static point source catalog, the position of each object (at the beginning of the integration) in RA, Dec or x,y must be provided, along with the object's magnitude in the filter used for the simulation. In addition, the velocity of the object must be specified. This is done in units of delta RA, delta Dec (arcsec/hour), or delta x, delta y (pixels/hour). ‘velocity_pixels’ must be placed in one of the top 4 lines of the file if the provided velocities are in units of pixels per hour rather than arcseconds per hour.

Below is an example catalog:

::

	#
	#abmag
	#
	#
	# List of point sources to create as moving targets (KBOs, asteroids, etc)
	# position can be x,y or RA,Dec. If x,y, put the phrase 'position_pixels' in one
	# of the top 4 lines of the file.
	# Velocity can be in units of pix/hour or arcsec/hour.
	# If using pix/hour, place 'velocity_pixels' in the second line of the file.
	# Note that if using velocities of pix/hour, the results will not be
	# strictly correct because in reality distortion will cause object's
	# velocities to vary in pixels/hour. Velocities in arcsec/hour will be
	# constant.
	x_or_RA    y_or_Dec   nircam_f200w_magnitude  x_or_RA_velocity   y_or_Dec_velocity
	53.0985    -27.8015       14                        180                 180

.. _moving_sersic:

Moving 2D Sersic Objects
------------------------

This option may be useful for simulating moving moons around a primary target that is being tracked. Similar to the static galaxy inputs, each moving target in this catalog must have an initial position in RA, Dec or x,y specified, along with a radius in arcseconds or pixels, ellipticity, position angle, Sersic index, and magnitude. In addition, velocities in the RA, Dec or x,y directions must be specified in units of arcseconds or pixels per hour.

::

	#
	#
	#
	#abmag
	# Columns 1 and 2 can be either x,y positions on the detector aperture (e.g.
	# 0,0 is lower left corner of the full frame of the subarray used for the
	# output) or RA,Dec location of the center of the source. If they are x,y
	# positions, make the top line of the file '# position_pixels'
	#
	# radius is the half-light radius in pixels or arcseconds. If in pixels
	# make the second line of the file '# radius_pixels
	#
	# pos_angle is the position angle of the semimajor axis, in degrees.
	# 0 causes the semi-major axis to be horizontal.
	x_or_RA   y_or_Dec  radius  ellipticity  pos_angle  sersic_index  nircam_f200w_magnitude  x_or_RA_velocity  y_or_Dec_velocity
	354.765   0.00064    1.0       0.25         20          2.0            16.000                  -0.5              -0.02


.. _moving_extended:

Moving Extended Sources
-----------------------

Similar to the catalog of static extended targets, this catalog contains a fits filename for each source containing the stamp image to use for the object, along with an initial position in RA, Dec or x,y, the object's magnitude, and position angle (of the array as read in from the fits file). In addition, velocities in the RA, Dec (arcsec/hour) or x,y directions (pixels/hour) must be specified.

::

	#
	#
	#
	#abmag
	# List of stamp image files to read in and use to create moving targets.
	# This is the method to use in order to create moving targets of
	# extended sources, like planets, moons, etc.
	# position can be x,y or RA,Dec. Velocity can be in units of pix/hour or arcsec/hour.
	# If using pix/hour, place 'velocity_pixels' in one of the top 4 lines.
	# Note that if using velocities of pix/hour, the results will not be
	# strictly correct because in reality distortion will cause object's
	# velocities to vary in pixels/sec. Velocities in arcsec/hour will be
	# constant.
	filename            x_or_RA    y_or_Dec   nircam_f200w_magnitude   pos_angle    x_or_RA_velocity   y_or_Dec_velocity
	ring_nebula.fits    0.007       0.003             12.0               0.0             -0.5               -0.02


.. _grism_tso_cat:

Grism Time Series Sources
-------------------------

When creating Grism TSO data, this catalog should contain the information on the TSO source parent body only. Other (background) sources should be placed in a catalog appropriate to their source type (point source, galaxy, extended). The index number for the TSO source should be set to 99999, to help ensure it will take presidence over the background sources when constructing the segmentation map. See the `TSO example notebook <https://github.com/spacetelescope/mirage/blob/master/examples/NIRCam_TSO_examples.ipynb>`_ for an example of how to use Mirage's catalog generation functionality to create a Grism TSO source catalog. As with the other catalog types, the source location can be given in RA, Dec or pixel x, y.

Most of the other columns in the catalog are specific to the `Batman <https://www.cfa.harvard.edu/~lkreidberg/batman/>`_ package, which Mirage uses when creating TSO data. See the Batman documentation for lists of the possible limb darkening models and associated coefficients. **Make sure that all time related entries use the same units.** Start_time and End_time are for the lightcurve relative to the start time of the exposure. Code development was done using a Start_time of 0.0 (i.e. the beginning of the exposure) and an End_time set to the maximum time the user wishes to simulate. If the End_time of the lightcurve comes before the end of the exposure, Mirage will automatically extend the lightcurve to fully encompass the exposure.

The Transmission_spectrum column should hold the name of an ascii file containing the transmission curve of the source. This is the Wavelength-dependent effective radius of the planet, in units of the stellar radius. Again, see the `TSO example notebook <https://github.com/spacetelescope/mirage/blob/master/examples/NIRCam_TSO_examples.ipynb>`_ for an example of this file.

Finally, this catalog contains magnitude columns similar to those in other catalog types. Note that along with this catalog, Grism TSO simulations require a file containing the SED of the star. If the SED in that file is given in flux density units, then the magnitudes in this catalog will be ignored. If the SED is normalized or scaled in some way, then the SED will be renormalized to the magnitude in this catalog corresponding to the filter of the TSO observation.

::

    # position_RA_Dec
    # vegamag
    #
    #
    index   x_or_RA     y_or_Dec   Semimajor_axis_in_stellar_radii  Orbital_inclination_deg  Eccentricity  Longitude_of_periastron  Limb_darkening_model  Limb_darkening_coeffs  Time_units  Start_time  End_time  Time_of_inferior_conjunction  Orbital_period      Transmission_spectrum      nircam_f444w_magnitude  nircam_f322w2_magnitude
    99999 66.37090333 -30.60044722             9.37                        83.3                   0.0               90.0                 nonlinear        "0.5, 0.1, 0.1, -0.1"     second       0.0       580.0             280.0                   3162.24     ./transmission_spectrum.txt             9.0                     9.05


.. _imaging_tso_cat:

Imaging Time Series Sources
---------------------------

When creating Imaging TSO data, this catalog should contain the information on the TSO source parent body only. Other (background) sources should be placed in a catalog appropriate to their source type (point source, galaxy, extended). The index number for the TSO source should be set to 99999, to help ensure it will take presidence over the background sources when constructing the segmentation map. See the `TSO example notebook <https://github.com/spacetelescope/mirage/blob/master/examples/NIRCam_TSO_examples.ipynb>`_ for an example of how to use Mirage's catalog generation functionality to create an Imaging TSO source catalog. As with the other catalog types, the source location can be given in RA, Dec or pixel x, y.

The lightcurve_file column should contain the name of an hdf5 file that contains a dataset with wavelength and flux arrays that describe the transit lightcurve at the wavelength range corresponding to the filter used in the observation. Again, see the `TSO example notebook <https://github.com/spacetelescope/mirage/blob/master/examples/NIRCam_TSO_examples.ipynb>`_ for an example of how to create this file.

Finally, this catalog contains magnitude columns similar to those in other catalog types.

::

    # position_RA_Dec
    # vegamag
    #
    #
    index    x_or_RA     y_or_Dec         lightcurve_file       nircam_f182m_magnitude  nircam_f210m_magnitude   nircam_f470n_magnitude
    99999  66.37090333 -30.60044722 ./example_lightcurve.hdf5          10.0                       9.5                    9.0
