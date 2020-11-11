#! /usr/bin/env python

"""This module contains utility functions related to source catalogs
"""
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
    catalogs = catalogs[sorter]

    # Create an array of min and max indexes and make sure the list is
    # monotonically increasing
    indexes = [[mn, mx] for mn, mx in zip(min_indexes, max_indexes)]
    indexes = [item for sublist in indexes for item in sublist]

    # Check that each element of the array is less than the preceding element
    overlaps = ~np.all(indexes[1:] >= indexes[:-1])
    return overlaps


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
    if obs_mode == 'imaging':
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
