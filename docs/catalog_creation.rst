.. _create_catalogs:

Creation of Source Catalogs
===========================

To help users create source catalogs in the proper format, `Mirage` contains the *catalog_generator.py* and *create_catalog.py* modules. *catalog_generator.py* contains classes for each of the :ref:`types of catalogs <catalogs>` accepted by `Mirage`. The *create_catalog.py* module contains functions that will query several astronomical databases and attempt to construct realistic point source and extragalactic catalogs for a given pointing. There is also a function for automated catalog creation from an input APT file. See the XXXXXX notebook in the `Mirage` repository for examples of all of this functionality.

.. _catalog_generator:

Catalog_generator
-----------------

The classes in this module can be used to simplify `Mirage` catalog creation from lists of source coordinates and magnitudes. Generally the user instantiates the class for the desired type of catalog with a list of RA and Dec or detector x, y source positions. Magnitude columns can then be added to the catalog by providing a list of magnitudes as well as the instrument and filter name associated with the magnitudes. Catalogs may contain multiple magnitude columns covering multiple instruments and filters.


.. _create_catalogs:

Create_catalogs
---------------

The functions in this module use astroquery (ADD LINK) to search astronomical databases and retrieve source lists for a given pointing.