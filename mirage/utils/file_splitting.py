#! /usr/bin/env python

"""Utility for determining how to split up seed images and dark files
that are too large. This is meant to roughly mimic the way that the
DMS pipeline will do the same job. Files that are spllit will have
``segNNN`` added to their names, where ``NNN`` is a number.
"""
import numpy as np
from mirage.utils.constants import FILE_SPLITTING_LIMIT


def find_file_splits(xdim, ydim, groups, integrations, pixel_limit=FILE_SPLITTING_LIMIT):
    """Determine the frame and/or integration numbers where a file
    should be split in order to keep the file size reasonable.

    Parameters
    ----------
    xdim : int
        Number of columns in the aperture

    ydim : int
        Number of rows in the aperture

    groups : int
        Number of groups (or frames) in an integration

    integrations : int
        Number of integrations in exposure

    pixel_limit : int
        Proxy for file size limit. Number of pixel read outs to treat
        as the upper limit to be contained in a single file.

    Returns
    -------
    split : bool
        Whether or not the exposure needs to be split into multiple
        files

    group_list : numpy.ndarray
        1D array listing the beginning group number of each file split

    integration_list : numpy.ndarray
        1d array listing the beginning integration number of each file
        split
    """
    pix_per_group = ydim * xdim
    pix_per_int = groups * pix_per_group
    observation = pix_per_int * integrations

    # Default = no splitting
    split = False
    group_list = np.array([0, groups+1])
    integration_list = np.array([0, integrations+1])

    # Check for splitting between groups first
    # i.e. splitting within an integration
    if pix_per_int > pixel_limit:
        split = True
        delta_group = np.int(pixel_limit / pix_per_group)
        group_list = np.arange(0, groups, delta_group)
        group_list = np.append(group_list, groups)
        print('Splitting within each integration:')
        print('group_list: ', group_list)
        integration_list = np.arange(integrations + 1)
        print('integration_list: ', integration_list)
    elif observation > pixel_limit:
        split = True
        print('Splitting by integration:')
        group_list = np.array([0, groups])
        delta_int = np.int(pixel_limit / pix_per_int)
        integration_list = np.arange(0, integrations, delta_int)
        integration_list = np.append(integration_list, integrations)
        print('group_list: ', group_list)
        print('integration_list: ', integration_list)

    return split, group_list, integration_list
