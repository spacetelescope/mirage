#! /usr/bin/env python

"""This module contains functions related to the addition of ghosts to NIRISS
simulations. Code originally written by Takahiro Morishita.
"""
import logging
import numpy as np
from astropy.io import fits, ascii
from astropy.table import Table

from mirage.utils.constants import DEFAULT_NIRISS_PTSRC_GHOST_FILE
from mirage.utils.flux_cal import fluxcal_info
from mirage.utils.utils import countrate_to_magnitude, magnitude_to_countrate


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
    logger = logging.getLogger('mirage.ghosts.niriss_ghosts.determine_ghost_stamp_filename')
    if source_type.lower() == 'point_source':
        default = DEFAULT_NIRISS_PTSRC_GHOST_FILE
    elif source_type.lower() == 'galaxies':
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
                ghost_stamp_filename = None
                logger.info(('No niriss_ghost_stamp filename for this source in the source catalog, and the default '
                             'file is set to None. Skipping ghost addition for this source.'))
    else:
        if default is not None:
            ghost_stamp_filename = default
        else:
            ghost_stamp_filename = None
            logger.info(('No niriss_ghost_stamp column in source catalog, and the default file is set to None. '
                         'Skipping ghost addition for this source.'))
    return ghost_stamp_filename


def get_gap(filter_name, pupil_name, gap_file, log_skipped_filters=True):
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

    log_skipped_filters : bool
        If True, any filter magnitudes that are not converted because
        the filter is not in the gap summary file will be logged.

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
        if log_skipped_filters:
            logger.info('Filter/Pupil {}/{} not found in the ghost gap file {}. Unable to add ghosts.'.format(filter_name, pupil_name, gap_file))
        return np.nan, np.nan, np.nan

    return xgap, ygap, frac


def get_ghost(x, y, flux, filter_name, pupil_name, gap_file, shift=0, log_skipped_filters=True):
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

    log_skipped_filters : bool
        Passed into get_gap. If True, any filter magnitudes that are not converted because
        the filter is not in the gap summary file will be logged.

    Returns
    -------
    xgs : float or numpy.array
        X-coordinates of ghosts

    ygs : float or numpy.array
        Y-coordinates of ghosts

    flux_gs : float or numpy.array
        fluxes of ghosts
    """
    xgap, ygap, frac = get_gap(filter_name, pupil_name, gap_file, log_skipped_filters=log_skipped_filters)

    xgap += shift
    ygap += shift

    xgs = 2 * xgap - x
    ygs = 2 * ygap - y

    try:
        flux_gs = flux * frac / 100.
    except:
        flux_gs = np.zeros(len(xgap), 'float')

    return xgs, ygs, flux_gs


def source_mags_to_ghost_mags(row, flux_cal_file, magnitude_system, gap_summary_file, filter_value_ghost, log_skipped_filters=False):
    """Works only for NIRISS. Given a row from a source catalog, create a ghost source
    catalog containing the magnitudes of the ghost associated with the source in all
    filters contained in the original catalog.

    Parameters
    ----------
    row : astropy.table.Row
        Row containing a source from a source catalog

    flux_cal_file : str
        Name of file containing the flux calibration information
        for all filters

    magnitude_system : str
        Magnitude system of the magnitudes in the input catalog

    gap_summary_file : str
        Name of file that controls ghost locations relative to sources

    log_skipped_filters : bool
        If False, notices of skipped magnitude translations will not be logged.
        This is convenient for catalogs with lots of sources, because otherwise
        there can be many repeats of the message.

    filter_value_ghost : str
        CLEAR or GR150, to select frac50 of ghosts from gap_summary_file.

    Returns
    -------
    ghost_row : astropy.table.Table
        Single row table containing the ghost magnitudes in all filters

    skipped_non_niriss_cols : bool
        If True, non NIRISS magnitude columns were found and skipped

    """
    logger = logging.getLogger('mirage.ghosts.niriss_ghosts.source_mags_to_ghost_mags')
    mag_cols = [key for key in row.colnames if 'magnitude' in key]
    ghost_row = Table()
    skipped_non_niriss_cols = False
    for mag_col in mag_cols:

        # Ghost magnitudes can currently be calcuated only for NIRISS
        if 'nircam' in mag_col.lower() or 'fgs' in mag_col.lower():
            skipped_non_niriss_cols = True
            continue

        col_parts = mag_col.split('_')
        if len(col_parts) == 3:
            filt = col_parts[1]
        else:
            raise ValueError('Unsupported magnitude column name for ghost conversion: {}'.format(mag_col))

        wave = int(filt[1:4])
        if wave > 200:
            filter_value = filt
            pupil_value = 'CLEARP'
        else:
            filter_value = 'CLEAR'
            pupil_value = filt

        # Get basic flux calibration information for the filter
        vegazeropoint, photflam, photfnu, pivot = \
            fluxcal_info(flux_cal_file, 'niriss', filter_value, pupil_value, 'NIS', 'N')

        # Convert source magnitude to count rate
        mag = float(row[mag_col])
        countrate = magnitude_to_countrate('niriss', filt, magnitude_system, mag,
                                            photfnu=photfnu, photflam=photflam,
                                            vegamag_zeropoint=vegazeropoint)

        # Get the count rate associated with the ghost
        if filter_value_ghost[0] == 'G':
            filter_value_ghost = 'GR150'
        _, _, ghost_countrate = get_ghost(1024, 1024, countrate, filter_value_ghost, pupil_value, gap_summary_file,
                                          log_skipped_filters=log_skipped_filters)

        # Convert count rate to magnitude
        ghost_mag = countrate_to_magnitude('niriss', filt, magnitude_system, ghost_countrate,
                                            photfnu=photfnu, photflam=photflam,
                                            vegamag_zeropoint=vegazeropoint)
        ghost_row[mag_col] = [ghost_mag]
    return ghost_row, skipped_non_niriss_cols

