"""This module contains functions for working with flux calibration
information, such as filter-based zeropoints, photflam, photfnu, and
pivot wavelegth values
"""
from astropy.io import ascii
from astropy.table import Column
import copy
import numpy as np
import scipy.special as sp


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


def fluxcal_info(params, usefilt, detector, module):
    """Retrive basic flux calibration information from the ascii file in
    the repository

    Parameters
    ----------
    params : dict
        Nested dictionary containing the input yaml file parameters

    usefile : str
        Either 'filter' or 'pupil', corresponding to which yaml
        parameter contains the optical element to use for flux cal

    Returns
    -------
    detector : str
        Name of the detector (e.g. NRCA1)
    """
    zpts = read_zeropoint_file(params['Reffiles']['flux_cal'])

    # In the future we expect zeropoints to be detector dependent, as they
    # currently are for FGS. So if we are working with NIRCAM or NIRISS,
    # manually add a Detector key to the dictionary as a placeholder.
    if params["Inst"]["instrument"].lower() in ["nircam", "niriss"]:
        zps = add_detector_to_zeropoints(detector, zpts)
    else:
        zps = copy.deepcopy(zpts)

    # Make sure the requested filter is allowed
    if params['Readout'][usefilt] not in zps['Filter']:
        raise ValueError(("WARNING: requested filter {} is not in the list of "
                          "possible filters.".format(params['Readout'][usefilt])))

    # Get the photflambda and photfnu values that go with
    # the filter
    mtch = ((zps['Detector'] == detector) &
            (zps['Filter'] == params['Readout'][usefilt]) &
            (zps['Module'] == module))
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