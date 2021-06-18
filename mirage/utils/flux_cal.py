"""This module contains functions for working with flux calibration
information, such as filter-based zeropoints, photflam, photfnu, and
pivot wavelegth values
"""
from astropy.io import ascii
from astropy.table import Column
import copy
import numpy as np
import scipy.special as sp

from mirage.utils.utils import standardize_filters


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
    fluxcal_file : str
        Name of file containing the zeropoint information for all filters

    instrument : str
        Instrument name

    filter_value : str
        Name of the optical element in filter wheel

    pupil_value : str
        Name of optical element in the pupil wheel

    detector : str
        Detector name

    module : str
        Module name

    Returns
    -------
    vegazeropoint : float
        Zeropoint of the given filter in vegamags

    photflam : float
        Photflam value for the given filter

    photfnu : float
        Photfnu value for the give filter

    pivot : float
        Pivot wavelength in microns for the given filter
    """
    zpts = read_zeropoint_file(fluxcal_file)

    # In the future we expect zeropoints to be detector dependent, as they
    # currently are for FGS. So if we are working with NIRCAM or NIRISS,
    # manually add a Detector key to the dictionary as a placeholder.
    if instrument.lower() == "niriss":
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

        # For entries that include the grism, substitute CLEAR for the GRISM
        if matching_pupil in ['GRISMR', 'GRISMC', 'GDHS0', 'GDHS60']:
            matching_pupil = 'CLEAR'

        mtch = ((zps['Detector'] == detector) &
                (zps['Filter'] == filter_value) &
                (zps['Pupil'] == matching_pupil) &
                (zps['Module'] == module))

    elif instrument.lower() in ['niriss', 'fgs']:
        matching_filter = filter_value.upper()
        if filter_value.upper() in ['CLEAR', 'CLEARP', 'GR150R', 'GR150C']:
            matching_filter = pupil_value.upper()

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


def mag_col_name_to_filter_pupil(colname):
    """Given a magnitude column name from a source catalog
    find the filter and pupil name

    Parameters
    ----------
    colname : str
        Column name (e.g 'nircam_f090w_magnitude', 'nircam_f090w_clear_magnitude')

    Returns
    -------
    filter_name : str
        Name of filter (e.g. 'f090w')

    pupil_name : str
        Name of pupil (e.g. 'clear')
    """
    mag_parts = colname.split('_')
    instrument = mag_parts[0].lower()
    if len(mag_parts) == 4:
        filt_pup_str = '{}/{}'.format(mag_parts[1], mag_parts[2])
    elif len(mag_parts) == 3:
        filt_pup_str = mag_parts[1]

    std_filt_pup_str = standardize_filters(instrument, [filt_pup_str])
    std_parts = std_filt_pup_str[0].split('/')
    if len(std_parts) == 2:
        filter_name = std_parts[0].lower()
        pupil_name = std_parts[1].lower()
    elif len(std_parts) == 1:
        filter_name = std_parts[0]
        pupil_name = 'NONE'
    return filter_name, pupil_name


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


def sersic_total_signal(effective_radius, sersic_index):
    """Calculate the total signal (out to infinity) associated with a 2D
    Sersic. Equation taken from JWST ETC. See here for more info:
    http://ned.ipac.caltech.edu/level5/March05/Graham/Graham2.html

    Parameters
    ----------
    effective_radius : float
        Radius that contains half the flux. R_e

    sersic_index : float
        Sersic index

    Returns
    -------
    sersic_total : float
        Total signal associated with the 2D Sersic profile
    """
    b_n = sp.gammaincinv(2 * sersic_index, 0.5)
    sersic_total = effective_radius**2 * 2 * np.pi * sersic_index * np.exp(b_n)/(b_n**(2 * sersic_index)) * sp.gamma(2 * sersic_index)
    return sersic_total


def sersic_fractional_radius(effective_radius, sersic_index, fraction_of_total, ellipticity):
    """For a 2D Sersic profile, calculate the semi-major and semi-minor axes
    lengths that contain a specified fraction of the total signal.
    See here for more info:
    http://ned.ipac.caltech.edu/level5/March05/Graham/Graham2.html

    Parameters
    ----------
    effective_radius : float
        Radius that contains half the flux. R_e

    sersic_index : float
        Sersic index

    fraction_of_total : float
        Fraction of the total signal desired within the calculated
        semi-major and semi-minor axes

    ellipticity : float
        Ellipticity of the 2D Sersic profile

    Returns
    -------
    radius : float
        Effective radius that encompasses the requested signal, assuming a circular profile

    semi_major : float
        Semi-major axis size that encompasses the requested signal

    semi-minor : float
        Semi-minor axis size that encompasses the requested signal
    """
    if fraction_of_total > 1.0:
        raise ValueError("fraction_of_total must be <= 1")

    b_n = sp.gammaincinv(2 * sersic_index, 0.5)
    x = sp.gammaincinv(2*sersic_index, fraction_of_total)
    sersic_total = effective_radius**2 * 2 * np.pi * sersic_index * np.exp(b_n)/(b_n**(2 * sersic_index)) * sp.gammainc(2 * sersic_index, x)
    radius = effective_radius * (x / b_n)**sersic_index
    semi_major = np.sqrt(radius**2 /  (1. - ellipticity))
    semi_minor = semi_major * (1. - ellipticity)
    return radius, semi_major, semi_minor