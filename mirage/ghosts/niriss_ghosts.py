#! /usr/bin/env python

"""This module contains functions related to the addition of ghosts to NIRISS
simulations. Code originally written by Takahiro Morishita.
"""
import logging
import numpy as np
from astropy.io import fits,ascii

from mirage.utils.constants import DEFAULT_NIRISS_PTSRC_GHOST_FILE


def determine_ghost_stamp_filename(row, source_type):
    """For a given type of source (point_source, galaxy, extended), determine
    the filename to use for the ghost stamp image. This may be provided in the
    source catalog, or there are default values to fall back to (for point sources)

    Parameters
    ----------
    row : astropy.table.row.Row
        Row from a source catalog table

    source_type : str
        Type of source contained in ```row```. Can be 'point_source', 'galaxy',
        or 'extended'

    Returns
    -------
    ghost_filename : str
        Name of a fits file containing a stamp image to use for the ghost source
        associated with the source in ```row```.
    """
    if source_type.lower() == 'point_source':
        default = DEFAULT_NIRISS_PTSRC_GHOST_FILE
    elif source_type.lower() == 'galaxy':
        default = None
    elif source_type.lower() == 'extended':
        default = None
    else:
        raise ValueError('Invalid source_type. Unable to continue.')

    if 'niriss_ghost_stamp' in row.colnames:
        if row['niriss_ghost_stamp'] is not None and row['niriss_ghost_stamp'].lower() != 'none':
            ghost_stamp_filename = row['niriss_ghost_stamp']
        else:
            if default is not None:
                ghost_stamp_filename = default
            else:
                raise ValueError(('Attempting to add a ghost source corresponding to a {} source '
                                  'but no niriss_ghost_stamp image file is given in the source catalog, '
                                  'and the default ghost stamp image file is None. Unable to create ghost source.'
                                  .format(source_type)))
    else:
        if default is not None:
            ghost_stamp_filename = default
        else:
            raise ValueError(('Attempting to add a ghost source corresponding to a {} source '
                              'source, but no niriss_ghost_stamp column is present in the source catalog, '
                              'and the default ghost stamp image file is None. Unable to create ghost source.'
                              .fomat(source_type)))
    return ghost_stamp_filename


def get_gap(filter_name, pupil_name, gap_file):
    """Get GAP coordinates corresponding to input filter. GAP is based on CV3 data.
    We currently do not know fractional flux, tab_gap['frac_50'], i.e. there may be positional dependence too.

    Parameters
    ----------
    filter_name : str
        Optical element in the filter wheel

    pupil_name : str
        Optical element in the pupil wheel

    gap_file : str
        Name of ASCII file contianing ghost location offsets

    Returns
    -------
    xgap : float
        X-coordinate on the detector of the ghost image

    ygap : float
        Y-coordinate on the detector of the ghost image

    frac : float
        Fraction of ghost flux relative to the source flux.

    """
    logger = logging.getLogger('mirage.ghosts.niriss_ghosts.get_gap')
    tab_gap = ascii.read(gap_file)
    iix = np.where((tab_gap['filt'] == filter_name.upper()) & (tab_gap['pupil'] == pupil_name.upper()))

    if len(iix[0]) > 0:
        xgap, ygap = tab_gap['gapx_50'][iix[0][0]], tab_gap['gapy_50'][iix[0][0]]
        frac = tab_gap['frac_50'][iix[0][0]]
    elif len(iix[0]) == 0:
        logger.info('Filter {} not found in the ghost gap file {}.'.format(filter_name, gap_file))
        return np.nan, np.nan, np.nan

    return xgap, ygap, frac


def get_ghost(x, y, flux, filter_name, pupil_name, gap_file, shift=0):
    """
    Calculates expected ghost positions given position of a source.

    Parameters
    ----------
    x : float or numpy.array
        X-coordinate of real source

    y : float or numpy.array
        Y-coordinate of real source

    flux : float or numpy.array
        fluxes for sources (array)

    filter_name : str
        Optical element in the filter wheel

    pupil_name : str
        Optical element in the pupil wheel

    gap_file : str
        Name of ASCII file contianing ghost location offsets

    shift : int
        may be used, because photutils is 0-based, while ds9 is not.

    Returns
    -------
    xgs : float or numpy.array
        X-coordinates of ghosts

    ygs : float or numpy.array
        Y-coordinates of ghosts

    flux_gs : float or numpy.array
        fluxes of ghosts
    """
    xgap, ygap, frac = get_gap(filter_name, pupil_name, gap_file)

    xgap += shift
    ygap += shift

    xgs = 2 * xgap - x
    ygs = 2 * ygap - y

    try:
        flux_gs = flux * frac / 100.
    except:
        flux_gs = np.zeros(len(xgap), 'float')

    return xgs, ygs, flux_gs
