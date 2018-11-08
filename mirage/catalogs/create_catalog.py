#! /usr/bin/env python

"""
Tools for generating Mirage-compatible catalogs from surveys
"""

import astropy.units as u
from astroquery.irsa import Irsa

from mirage.catalogs.catalog_generator import PointSourceCatalog, GalaxyCatalog


def get_2MASS_ptsrc_catalog(ra, dec, box_width):
    """Query the 2MASS All-Sky Point Source Catalog in a square region around the RA and Dec
    provided. Box width must be in units of arcseconds

    Parameters
    ----------
    something
    """
    # Don't artificially limit how many sources are returned
    Irsa.ROW_LIMIT = 1000000

    ra_dec_string = "{}  {}".format(ra, dec)
    query_table = Irsa.query_region(ra_dec_string, catalog='fp_psc', spatial='Box',
                                    width=box_width * u.arcsec)

    # Exclude any entries with missing RA or Dec values
    radec_mask = filter_bad_ra_dec(query_table)
    cat = PointSourceCatalog(ra=query_table['ra'].data.data[radec_mask],
                             dec=query_table['dec'].data.data[radec_mask])

    # Add the J, H, and K magnitudes as they may be useful for magnitude conversions later
    # Add the values that have had fill_values applied. The fill_value is 1e20.
    for key in ['j_m', 'h_m', 'k_m']:
        data = query_table[key].filled().data
        cat.add_magnitude_column(data, instrument='2MASS', filter_name=key)

    return cat


#def from_luminosity_function(self, luminosity_function):
#    more customizable


#def galactic_plane(self):
#    create representative scene looking into galactic galactic_plane


#def galactic_bulge(self):
#    looking into galactic_bulge


#def out_of_plane(self):
#    looking out of the plane of galaxy


def filter_bad_ra_dec(table_data):
    """Use the column masks to find which entries have bad RA or Dec values.
    These will be excluded from the Mirage catalog

    Parameters
    ----------
    something
    """
    ra_data = table_data['ra'].data.data
    ra_mask = ~table_data['ra'].data.mask
    dec_data = table_data['dec'].data.data
    dec_mask = ~table_data['dec'].data.mask
    position_mask = ra_mask & dec_mask
    return position_mask
