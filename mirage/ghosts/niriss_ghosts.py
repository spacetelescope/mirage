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

