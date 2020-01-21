.. _catalog_generation:

Tools for Creating Source Catalogs
==================================

To help users create source catalogs in the proper format, `Mirage` contains the `catalog generator <https://github.com/spacetelescope/mirage/blob/master/mirage/catalogs/catalog_generator.py>`_ and `create catalog <https://github.com/spacetelescope/mirage/blob/master/mirage/catalogs/create_catalog.py>`_ modules. *catalog_generator.py* contains classes for each of the types of catalogs accepted by `Mirage`. The *create_catalog.py* module contains functions that will query several astronomical databases and attempt to construct realistic point source and extragalactic catalogs for a given pointing. There is also a function for automated catalog creation from an input APT file. See the `Catalog Generation Tools notebook <https://github.com/spacetelescope/mirage/blob/master/examples/Catalog_Generation_Tools.ipynb>`_ in the `Mirage` repository for examples of all of this functionality.

.. _catalog_generator:

Catalog_generator - Create catalogs from a list of sources
----------------------------------------------------------

The classes in this module can be used to simplify `Mirage` catalog creation from lists of source coordinates and magnitudes. Generally the user instantiates the class for the desired type of catalog with a list of RA and Dec or detector x, y source positions. Magnitude columns can then be added to the catalog by providing a list of magnitudes as well as the instrument and filter name associated with the magnitudes. Catalogs may contain multiple magnitude columns covering multiple instruments and filters. Examples of how to use the catalog generator functionality are given in the `Catalog Generation Tools notebook <https://github.com/spacetelescope/mirage/blob/master/examples/Catalog_Generation_Tools.ipynb>`_ notebook. Below, we show examples for the most common catalogs: point sources and galaxies.

Point Source Catalog
++++++++++++++++++++

The examples below create point source catalogs containing 10 sources. By supplying source locations via the `ra` and `dec` keywords, `Mirage` assumes that the locations are Right Ascention and Declination values given in decimal degrees. If you wish to specify source locations in units of row and column indexes on the detector, then the `x` and `y` keywords should be used, as seen in the second example below.

Columns containing AB magnitudes in the NIRCam F200W and F212N filters, as well as the NIRISS F090W filter are added to the catalog using the `add_magnitude_column` method. Note that the instrument and filter names associated with the magnitude values must be supplied in these calls. The `mag_sys` keyword is optional and allows uses to specify the magnitude system. Allowed values are "abmag", "stmag", and "vegamag". The default is "abmag". **Note that all magnitude columns in a given catalog MUST be in the same magnitude system.**

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

Create_catalogs - create catalogs using online astronomical databases
---------------------------------------------------------------------

The functions in this module use `astroquery <https://astroquery.readthedocs.io/en/latest/>`_ to search astronomical databases and retrieve source lists for a given pointing. In this way, a user can quickly generate reasonably realistic catalogs of point sources and galaxies for a given pointing.

The **get_all_catalogs** function takes the RA and Dec of a particular pointing along with the width in arcseconds of the area for which to produce the catalog, and queries multiple databases to produce a point source catalog. An example call to create a 120 x 120 arcsecond catalog is shown below. The resulting point source catalog can then be placed in the :ref:`pointSource <pointsource>` entry of the :ref:`yaml input file <example_yaml>`. The *besancon_catalog_file* in the command below is the result from a query of the `Besancon model <https://model.obs-besancon.fr/modele_home.php>`_. Details of how to query the model and download the result are shown in the :ref:`Background Stars <background_stars>` section below.

::

    from mirage.catalogs import create_catalogs

    ra = 80.4  # degrees
    dec = -69.8  # degrees
    box_width = 120  # arcseconds
    filter_list = ['F444W', 'F480M']
    cat, mag_column_names = create_catalog.get_all_catalogs(ra, dec, box_width, besancon_catalog_file='besancon.cat',
                                                            instrument='NIRCAM', filters=filter_list
                                                            )


.. _foreground_stars:

Foreground Stars
++++++++++++++++

A realistic list of foreground stars is compiled by querying the `2MASS <https://astroquery.readthedocs.io/en/latest/irsa/irsa.html>`_, `WISE <https://astroquery.readthedocs.io/en/latest/irsa/irsa.html>`_, and `GAIA <https://astroquery.readthedocs.io/en/latest/gaia/gaia.html>`_ catalogs using the given pointing. Using the retrieved magnitudes in the various bands associated with these surveys, Mirage converts these to magnitude values in the requested NIRCam or NIRISS filters. Note that these queries return stars only down to about V=16. For dimmer stars, you can query the `Besancon model <https://model.obs-besancon.fr/modele_home.php>`_. This process is described in the :ref:`Background Stars <background_stars>` section below.


.. _background_stars:

Background Stars
++++++++++++++++

To obtain a list of stars dimmer than those returned in the :ref:`Foreground Stars <foreground_stars>` search, Mirage uses queries the `Besancon model <https://model.obs-besancon.fr/modele_home.php>`_ of stars in the Milky Way. This query returns a **representative sample** (in terms of luminosity distribution) of stars for a given pointing, including a realistic stellar density and realistic luminosity distribution. Note that these stars are not actual stars in the sky. Due to the way in which the model is queried and results are returned, the use of a Besancon-derived catalog is a two-step process. First, you must create an account on the `Besancon model page <https://model.obs-besancon.fr/modele_home.php>`_. Once the account is activated, you can query the model using Mirage's wrapper function, as shown below. A more complete example of this is given in the `Example uses of Mirage catalog generators <https://github.com/spacetelescope/mirage/blob/master/examples/Catalog_Generation_Tools.ipynb>`_ jupyter notebook.

::

    from mirage.catalogs import create_catalog
    ra = 224.2  # degrees
    dec = -65.54  # degrees
    box_width = 200  # arcseconds
    create_catalog.besancon(ra, dec, box_width, username='hilbert', kmag_limits=(17, 30))

Once the query is complete, you will receive an email with a link to download the resulting ascii table. With the saved table in hand, you can then transform the source magnitudes from JHK to the JWST filters of interest, and combine the catalog with query results from GAIA/2MASS/WISE. This combined catalog can then be used as input to a Mirage simulation. Again, see the `Example uses of Mirage catalog generators <https://github.com/spacetelescope/mirage/blob/master/examples/Catalog_Generation_Tools.ipynb>`_ jupyter notebook for an example.


Background Galaxies
+++++++++++++++++++

For a given pointing, Mirage can also generate a catalog containing a **representative sample** of background galaxies. Similar to the Besancon query described above, Mirage will generate a catalog containing a realistic density of galaxies across the field at reasonable magnitudes. To accomplish this, Mirage queries the `GOODS-S catalog from 3DHST <https://3dhst.research.yale.edu/Data.php>`_ and extracts an appropriate number of galaxies to populate the catalog at a reasonable density. Currently this function will fail if the user requests a catalog with an area larger than the GOODS-S field: 606,909 arcsec :sup:`2`. An example is shown below. The resulting file can then be placed in the :ref:`galaxyListFile <galaxylistfile>` entry of the :ref:`yaml input file <example_yaml>`.