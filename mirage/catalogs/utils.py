#! /usr/bin/env python

"""This module contains utility functions related to source catalogs
"""
import os

from astropy.io import ascii
import numpy as np
from mirage.utils.constants import IMAGING_ALLOWED_CATALOGS, WFSS_ALLOWED_CATALOGS, \
                                   TS_IMAGING_ALLOWED_CATALOGS, TS_GRISM_ALLOWED_CATALOGS


def catalog_index_check(catalogs):
    """Check to see if there are any overlaps in the index
    values of the provided catalogs. Search only by looking
    at the minimum and maximum index values, rather than
    searching through all the individual values.

    Parameters
    ----------
    catalogs : list
        List of catalog filenames

    Returns
    -------
    overlaps : bool
        True if there are overlapping index values

    max_index : int
        Maximum index value across all input catalogs
    """
    min_indexes = []
    max_indexes = []
    for catalog in catalogs:
        cat = ascii.read(catalog)
        if 'index' in cat.colnames:
            # Collect the min and max index values from each catalog
            min_val = np.min(cat['index'])
            min_indexes.append(min_val)
            max_indexes.append(np.max(cat['index']))
            if min_val == 0:
                raise ValueError(("{} has a source with an index value of 0. This is not allowed. "
                                  "Zero is reserved for pixels in the segmentation map with no source flux.".format(catalog)))
            if min_val < 0:
                raise ValueError("{} has sources with negative index numbers. This is not allowed.".format(catalog))

        else:
            raise ValueError("{} does not have an 'index' column. Cannot compare to other catalogs.".format(catalog))

    # Reorder the min_indexes to be monotonically increasing.
    # Apply the same sorting to the max_indexes and the catalog list
    min_indexes = np.array(min_indexes)
    max_indexes = np.array(max_indexes)
    sorter = np.argsort(min_indexes)
    min_indexes = min_indexes[sorter]
    max_indexes = max_indexes[sorter]

    # Create an array of min and max indexes and make sure the list is
    # monotonically increasing
    indexes = [[mn, mx] for mn, mx in zip(min_indexes, max_indexes)]
    indexes = np.array([item for sublist in indexes for item in sublist])

    # Check that each element of the array is less than the preceding element
    overlaps = ~np.all(indexes[1:] >= indexes[:-1])

    # Also return the maximum index value
    max_index = np.max(indexes)

    return overlaps, max_index


def determine_used_cats(obs_mode, cat_dict):
    """Return a list of the source catalogs that will be used by Mirage,
    based on the observation mode

    Parameters
    ----------
    obs_mode : str
        e.g. 'imaging', 'wfss'

    cat_dict : dict
        Dictionary containing catalog names. Keys should match those
        in the catalog entries of the Mirage yaml input file

    Returns
    -------
    cats : list
        List of catalogs that will be used for the given observing mode
    """
    if obs_mode in ['imaging', 'pom', 'ami', 'coron']:
        possible_cats = [cat_dict[entry] for entry in IMAGING_ALLOWED_CATALOGS]
    elif obs_mode == 'wfss':
        possible_cats = [cat_dict[entry] for entry in WFSS_ALLOWED_CATALOGS]
    elif obs_mode == 'ts_imaging':
        possible_cats = [cat_dict[entry] for entry in TS_IMAGING_ALLOWED_CATALOGS]
    elif obs_mode == 'ts_grism':
        possible_cats = [cat_dict[entry] for entry in TS_GRISM_ALLOWED_CATALOGS]

    # Remove any catalogs that are set to None
    cats = [ele for ele in possible_cats if str(ele).lower() != 'none']
    return cats


def get_nonsidereal_catalog_name(cat_dict, target_name, instrument_name):
    """Given a dictionary or nested dictionary of source catalogs,
    return the name of the non-sidereal catalog for the given target
    name.

    Parameters
    ----------
    cat_dict : dict
        Dictionary of source catalogs as input by the user

    target_name : str
        Target name to search for

    instrument_name : str
        Name of instrument used to observe ``target_name``

    Returns
    -------
    cat_file : str
        Name of source catalog
    """
    target_cat = cat_dict[target_name]
    if 'moving_target_to_track' in target_cat.keys():
        if not os.path.isfile(target_cat['moving_target_to_track']):
            raise ValueError(("{} is listed as the non-sidereal target catalog for "
                              "{}, but it appears this file does not exist.".format(target_cat['moving_target_to_track'],
                                                                                    target_name)))
        else:
            cat_file = target_cat['moving_target_to_track']
    else:
        if instrument_name not in target_cat.keys():
            raise ValueError(("Catalog dictionary does not contain a 'moving_target_to_track' catalog "
                              "for {}. Unable to proceed.".format(target_name)))
        else:
            if 'moving_target_to_track' in target_cat[instrument_name].keys():
                if not os.path.isfile(target_cat[instrument_name]['moving_target_to_track']):
                    raise ValueError(("{} is listed as the non-sidereal target catalog for "
                                      "{}, but it appears this file does not exist."
                                      .format(target_cat['moving_target_to_track'], target_name)))
                else:
                    cat_file = target_cat[instrument_name]['moving_target_to_track']
            else:
                raise ValueError(("Catalog dictionary does not contain a 'moving_target_to_track' catalog "
                                  "for {}. Unable to proceed.".format(target_name)))
    return cat_file


def read_nonsidereal_catalog(filename):
    """Read in a Mirage formatted non-sidereal source catalog

    Paramters
    ---------
    filename : str
        Name of ascii catalog

    Returns
    -------
    catalog_table : astropy.table.Table
        Catalog contents

    pixelflag : bool
        True if the source position is in units of detector (x, y)
        False if RA, Dec

    pixelvelflag : bool
        True if the source velocity is given in pixels/hour.
        False if arcsec/hour
    """
    catalog_table = ascii.read(filename, comment='#')

    # Check to see whether the position is in x,y or ra,dec
    pixelflag = False
    try:
        if 'position_pixels' in catalog_table.meta['comments'][0:4]:
            pixelflag = True
    except:
        pass

    # If present, check whether the velocity entries are pix/sec
    # or arcsec/sec.
    pixelvelflag = False
    try:
        if 'velocity_pixels' in catalog_table.meta['comments'][0:4]:
            pixelvelflag = True
    except:
        pass
    return catalog_table, pixelflag, pixelvelflag
