.. _catalog_generation:

Tools for Creating Source Catalogs
==================================

To help users create source catalogs in the proper format, `Mirage` contains the `catalog generator <https://github.com/spacetelescope/mirage/blob/master/mirage/catalogs/catalog_generator.py>`_ and `create catalog <https://github.com/spacetelescope/mirage/blob/master/mirage/catalogs/create_catalog.py>`_ modules. *catalog_generator.py* contains classes for each of the types of catalogs accepted by `Mirage`. The *create_catalog.py* module contains functions that will query several astronomical databases and attempt to construct realistic point source and extragalactic catalogs for a given pointing. There is also a function for automated catalog creation from an input APT file. See the `Catalog Generation Tools notebook <https://github.com/spacetelescope/mirage/blob/master/examples/Catalog_Generation_Tools.ipynb>`_ in the `Mirage` repository for examples of all of this functionality.

.. _catalog_generator:

Catalog_generator
-----------------

The classes in this module can be used to simplify `Mirage` catalog creation from lists of source coordinates and magnitudes. Generally the user instantiates the class for the desired type of catalog with a list of RA and Dec or detector x, y source positions. Magnitude columns can then be added to the catalog by providing a list of magnitudes as well as the instrument and filter name associated with the magnitudes. Catalogs may contain multiple magnitude columns covering multiple instruments and filters. Examples of how to use the catalog generator functionality are given in the `Catalog Generation Tools notebook <https://github.com/spacetelescope/mirage/blob/master/examples/Catalog_Generation_Tools.ipynb>`_ notebook. Below, we show examples for the most common catalogs: point sources and galaxies.

Point Source Catalog
++++++++++++++++++++

The examples below create point source catalogs containing 10 sources. By supplying source locations via the `ra` and `dec` keywords, `Mirage` assumes that the locations are Right Ascention and Declination values given in decimal degrees. (Will it work with string representations?) If you wish to specify source locations in units of row and column indexes on the detector, then the `x` and `y` keywords should be used, as seen in the second example below.

Columns containing AB magnitudes in the NIRCam F200W and F212N filters, as well as the NIRISS F090W filter are added to the catalog using the `add_magnitude_column` method. Note that the instrument and filter names associated with the magnitude values must be supplied in these calls. The `mag_sys` keyword is optional and allows uses to specify the magnitude system. Allowed values are "abmag", "stmag", and "vegamag". The default is "abmag".

Once the catalog has been created and at least one magnitude column has been added, the table can be printed to the screen using the `table` method.

The `save` method will save the table in an ascii file in the appropriate format for Mirage to use.

::

    import numpy as np
    from mirage.catalogs import catalog_generator

    ra_list = np.random.random(10) + 53.5
    dec_list = np.random.random(10) - 67.2

    nrc_f200w_mag = np.random.random(10) + 16.
    nrc_f212n_mag = np.random.random(10) + 19.
    nis_f090w_mag = np.random.random(10) + 15.5

    ptsrc = catalog_generator.PointSourceCatalog(ra=ra_list, dec=dec_list)
    ptsrc.add_magnitude_column(nrc_f200w_mag, instrument='nircam', filter='F200W', mag_sys='abmag')
    ptsrc.add_magnitude_column(nrc_f212n_mag, instrument='nircam', filter='F212N', mag_sys='abmag')
    ptsrc.add_magnitude_column(nis_f090w_mag, instrument='niriss', filter='F090W', mag_sys='abmag')
    ptsrc.save('point_sources.cat')
    ptsrc.table()

::

    x_list = np.random.random(10) * 2048
    y_list = np.random.random(10) * 2048

    ptsrc = catalog_generator.PointSourceCatalog(x=x_list, y=y_list)
    ptsrc.add_magnitude_column(nrc_f200w_mag, instrument='nircam', filter='F200W', mag_sys='abmag')
    ptsrc.add_magnitude_column(nrc_f212n_mag, instrument='nircam', filter='F212N', mag_sys='abmag')
    ptsrc.add_magnitude_column(nis_f090w_mag, instrument='niriss', filter='F090W', mag_sys='abmag')
    ptsrc.save('point_sources_xy.cat')

Galaxy Catalog
++++++++++++++

The example below creates a galaxy catalog. The main difference compared to the point source catalog is that the user must provide values for the ellipticity, sersic index, and position angle when instantiating the class.

::

    import numpy as np
    from mirage.catalogs import catalog_generator

    ra_list = np.random.random(10) + 53.5
    dec_list = np.random.random(10) - 67.2
    ellipticity = np.random.random(10) * 0.75
    sersic_index = np.random.random(10) * 4.
    position_angle = np.random.random(10) * 359.

    nrc_f200w_mag = np.random.random(10) + 16.
    nrc_f212n_mag = np.random.random(10) + 19.
    nis_f090w_mag = np.random.random(10) + 15.5

    gal = catalog_generator.GalaxyCatalog(ra=ra_list, dec=dec_list, ellipticity=ellipticity,
                                          sersic_index=sersic_index, position_angle=position_angle)
    gal.add_magnitude_column(nrc_f200w_mag, instrument='nircam', filter='F200W', mag_sys='abmag')
    gal.add_magnitude_column(nrc_f212n_mag, instrument='nircam', filter='F212N', mag_sys='abmag')
    gal.add_magnitude_column(nis_f090w_mag, instrument='niriss', filter='F090W', mag_sys='abmag')
    gal.save('galaxies.cat')
    gal.table()



.. _create_catalogs:

Create_catalogs
---------------

The functions in this module use astroquery (ADD LINK) to search astronomical databases and retrieve source lists for a given pointing. In this way, a user can quickly generate reasonably realistic catalogs of point sources and galaxies for a given pointing.

.. _foreground_stars:

Foreground Stars
++++++++++++++++

A realistic list of foreground stars is compiled by querying the 2MASS, WISE, and XXX catalogs using the given pointing. Using the retrieved magnitudes in the various bands associated with these surveys, `Mirage` converts these to magnitude values in the requested NIRCam or NIRISS filters. Note that these queries return stars only down to about V=16 (CONFIRM). For dimmer stars, `Mirage` queries the Besancon model (reference here). This process is described in the :ref:`Background Stars <background_stars>` section below.

::

    from mirage.catalogs import create_catalogs

    ra = 53.1
    dec = -67.8

    something = create_catalogs.queryXXXXX(ra=ra, dec=dec)
    print(something)

.. _background_stars:

Background Stars
++++++++++++++++

To obtain a list of stars dimmer than those returned in the :ref:`Foreground Stars <foreground_stars>` search, `Mirage` uses `astroquery` to query the Besancon model of stars in the Milky Way. This query returns a representative sample of stars for a given pointing, including a realistic stellar density and realistic luminosity distribution. However these stars are not actual stars in the sky.

::

    from mirage.catalogs import create_catalogs

    ra = 53.1
    dec = -67.8

    something = create_catalogs.besancon(ra=ra, dec=dec)
    print(something)




Background Galaxies
+++++++++++++++++++

For a given pointing, Mirage can generate a catalog containing a representative sample of background galaxies. Similar to the background star functionality described above, Mirage will generate a catalog containing a pseudo-realistic density of galaxies across the field at reasonable magnitudes. To accomplish this, `Miage` queries the **(DOUBLE CHECK)** 3D HST (add link) catalog and extracts an appropriate subarray.....blah blah