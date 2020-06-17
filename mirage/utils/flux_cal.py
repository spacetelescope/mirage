"""This module contains functions for working with flux calibration
information, such as filter-based zeropoints, photflam, photfnu, and
pivot wavelegth values
"""
from astropy.io import ascii
from astropy.table import Column
import copy
import numpy as np


def add_detector_to_zeropoints(detector, zeropoint_table):
    """Manually add detector dependence to the zeropoint table for
    NIRCam and NIRISS simualtions. This is being done as a placeholder
    for the future, where we expect zeropoints to be detector-dependent.

    Parameters
    ----------
    detector : str
        Name of detector to add to the table

    zeropoint_table : astropy.table.Table
        Table of filter zeropoint information

    Returns
    -------
    base_table : astropy.table.Table
        Copy of ``zeropoint_table`` with Detector column added
    """
    # Add "Detector" to the list of column names
    base_table = copy.deepcopy(zeropoint_table)
    num_entries = len(zeropoint_table)
    det_column = Column(np.repeat(detector, num_entries), name="Detector")
    base_table.add_column(det_column, index=0)
    return base_table


def fluxcal_info(fluxcal_file, instrument, filter_value, pupil_value, detector, module):
    """Retrive basic flux calibration information from the ascii file in
    the repository

    Parameters
    ----------
    params : dict
        Nested dictionary containing the input yaml file parameters

    Returns
    -------
    detector : str
        Name of the detector (e.g. NRCA1)
    """
    zpts = read_zeropoint_file(fluxcal_file)

    # In the future we expect zeropoints to be detector dependent, as they
    # currently are for FGS. So if we are working with NIRCAM or NIRISS,
    # manually add a Detector key to the dictionary as a placeholder.
    if instrument.lower() in ["nircam", "niriss"]:
        zps = add_detector_to_zeropoints(detector, zpts)
    else:
        zps = copy.deepcopy(zpts)

    # Get the photflambda and photfnu values that go with
    # the filter
    if instrument.lower() == 'nircam':
        # WLP8 and WLM8 have the same throughput, so the zeropoint file
        # contains only the entry for WLP8. If the user gave WLM8, then
        # be sure to look for the corresponding WLP8 entry.
        matching_pupil = pupil_value
        if matching_pupil == 'WLM8':
            matching_pupil = 'WLP8'

        mtch = ((zps['Detector'] == detector) &
                (zps['Filter'] == filter_value) &
                (zps['Pupil'] == matching_pupil) &
                (zps['Module'] == module))

    elif instrument.lower() in ['niriss', 'fgs']:
        matching_filter = filter_value
        if filter_value.upper() in ['CLEAR', 'CLEARP']:
            matching_filter = pupil_value

        mtch = ((zps['Detector'] == detector) &
                (zps['Filter'] == matching_filter) &
                (zps['Module'] == module))

    # Make sure the requested filter/pupil is allowed
    if not any(mtch):
        raise ValueError(("WARNING: requested filter and pupil values of {} and {} are not in the list of "
                          "possible options.".format(filter_value, pupil_value)))

    vegazeropoint = zps['VEGAMAG'][mtch][0]
    photflam = zps['PHOTFLAM'][mtch][0]
    photfnu = zps['PHOTFNU'][mtch][0]
    pivot = zps['Pivot_wave'][mtch][0]

    return vegazeropoint, photflam, photfnu, pivot


def read_zeropoint_file(filename):
    """Read in the ascii table containing all of the flux calibration
    information

    Parameters
    ----------
    filename : str
        Name of ascii file

    Returns
    -------
    flux_table : astropy.table.Table
        Table of flux calibration information
    """
    flux_table = ascii.read(filename)
    return flux_table
