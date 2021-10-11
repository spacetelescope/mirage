"""SIAF interface module to support accessing SIAF information.

This module provides ``mirage`` with functions to interface SIAF content via the pysiaf module.

Authors
-------

    - Johannes Sahlmann
    - Bryan Hilbert

Use
---

    This module can be imported and used with

    ::

        from mirage.utils import siaf_interface

"""
import os
import logging
import numpy as np

import pysiaf
from pysiaf import iando

from mirage.logging import logging_functions
from ..utils import rotations
from ..utils import set_telescope_pointing_separated as set_telescope_pointing
from mirage.utils.constants import LOG_CONFIG_FILENAME, STANDARD_LOGFILE_NAME


classdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../'))
log_config_file = os.path.join(classdir, 'logging', LOG_CONFIG_FILENAME)
logging_functions.create_logger(log_config_file, STANDARD_LOGFILE_NAME)


def aperture_ra_dec(siaf_instance, aperture_name, ra, dec, telescope_roll, output_apertures):
    """For a given aperture with a known RA, Dec, and telescope roll angle,
    calculate the RA, Dec values at the reference location in a list of
    other apertures.

    Parameters
    ----------
    siaf_instance : pysiaf.Siaf
        Instance of SIAF for a single instrument

    aperture_name : str
        Aperture name (e.g. "NRCA1_FULL")

    ra : float
        RA value of pointing in degrees

    dec : float
        Dec value of pointing in degrees

    telescope_roll : float
        PA_V3, Position angle of the telescope in degrees

    output_apertures : list
        List of aperture names to calculate RA, Dec for

    Returns
    -------
    aperture_pointing : dict
        Dictionary with output_apertures as keys. Values are (ra, dec)
        tuples
    """
    local_roll, att_matrix, fullframesize, subarray_boundaries = get_siaf_information(siaf_instance,
                                                                                      aperture_name,
                                                                                      ra, dec,
                                                                                      telescope_roll)
    aperture_pointing = {}
    for aperture in output_apertures:
        siaf_out = siaf_instance[aperture]
        out_ra, out_dec = pysiaf.rotations.pointing(att_matrix, siaf_out.V2Ref, siaf_out.V3Ref)
        aperture_pointing[aperture] = (out_ra, out_dec)
    return aperture_pointing


def aperture_xy_to_radec(x, y, instrument, aperture, fiducial_ra, fiducial_dec, pav3):
    """For a given aperture and roll angle, translate a given detector
    (x, y) location to RA, Dec

    Parameters
    ----------
    x : float
        X-coordinate within ```aperture```

    y : float
        Y-coordinate within ```aperture```

    instrument : str
        Name of JWST instrument (e.g. 'nircam')

    aperture : str
        Name of aperture (e.g. 'NRCA1_FULL')

    fiducial_ra : float
        Right ascention value at the reference location of the aperture,
        in decimal degrees

    fiducial_dec : float
        Declination value at the reference location of the aperture,
        in decimal degrees

    pav3 : float
        Telescope roll angle, in degrees

    Returns
    -------
    ra : float
        RA corresponding to (x, y)

    dec : float
        Dec corresponding to (x, y)
    """
    instrument_siaf = siaf_interface.get_instance(instrument)
    siaf = instrument_siaf[aperture]
    local_roll, attitude_matrix, ffsize, \
            subarray_bounds = get_siaf_information(instrument, aperture, fiducial_ra,
                                                   fiducial_dec, pav3)
    loc_v2, loc_v3 = siaf.sci_to_tel(x + 1, y + 1)
    ra, dec = pysiaf.utils.rotations.pointing(attitude_matrix, loc_v2, loc_v3)
    return ra, dec


def get_instance(instrument):
    """Return an instance of a pysiaf.Siaf object for the given instrument

    Parameters
    ----------
    instrument : str
        Name of instrument

    Returns
    -------
    siaf : pysiaf.Siaf
        Siaf object for the requested instrument
    """
    siaf = pysiaf.Siaf(instrument)
    return siaf


def get_siaf_information(siaf_instance, aperture_name, ra, dec, telescope_roll, v2_arcsec=None,
                         v3_arcsec=None, verbose=False):
    """Use pysiaf to get aperture information.

    Parameters
    ----------
    siaf_instance : pysiaf.Siaf
        Instance of SIAF for a single instrument

    aperture_name : str
        Aperture name (e.g. "NRCA1_FULL")

    ra : float
        RA value of pointing in degrees

    dec : float
        Dec value of pointing in degrees

    telescope_roll : float
        PA_V3, Position angle of the telescope in degrees

    v2_arcsec : float
        The V2 value in arcseconds of the reference location for the
        instrument aperture

    v3_arcsecc : float
        The V3 value in arcseconds of the reference location for the
        instrument aperture

    verbose : bool
        Print extra information to the screen

    Returns
    -------
    local_roll : float
        Local roll angle at the reference location of the aperture

    att_matrix : matrix
        Attitude matrix used to relate RA, Dec, local roll angle to V2, V3

    fullframesize : int
        Number of columns in the  given aperture

    subarray_boundaries : list
        List of full-frame coordinates corresponding to the minimum and maximum
        values of x and y in the given aperture
    """
    logger = logging.getLogger('mirage.utils.siaf_interface.get_siaf_information')

    # Select the correct aperture
    siaf = siaf_instance[aperture_name]

    if v2_arcsec is None:
        v2_arcsec = siaf.V2Ref
    if v3_arcsec is None:
        v3_arcsec = siaf.V3Ref

    local_roll = set_telescope_pointing.compute_local_roll(telescope_roll,
                                                           ra, dec, v2_arcsec, v3_arcsec)

    # Create attitude_matrix
    att_matrix = rotations.attitude(v2_arcsec, v3_arcsec, ra, dec, local_roll)

    # Get full frame size
    fullframesize = siaf.XDetSize

    # Subarray boundaries in full frame coordinates
    try:
        xcorner, ycorner = sci_subarray_corners(siaf_instance.instrument, aperture_name, siaf=siaf_instance)
        subarray_boundaries = [xcorner[0], ycorner[0], xcorner[1], ycorner[1]]
    except (RuntimeError, TypeError) as e:  # e.g. NIRSpec NRS_FULL_MSA aperture
        if verbose:
            logger.info('get_siaf_information raised error:\n{}\nIgnoring it.'.format(e))
        subarray_boundaries = [0, 0, 0, 0]
    return local_roll, att_matrix, fullframesize, subarray_boundaries


def sci_subarray_corners(instrument, aperture_name, siaf=None, verbose=False):
    """Return the two opposing aperture corners in the SIAF Science frame of the full-frame SCA.

    This function serves as interface between the SIAF information accessible via the pysiaf package
    and the subarray information formatted for use by mirage.

    Parameters
    ----------
    instrument : str
        JWST instrument name with correct capitalization
    aperture_name : str
        SIAF aperture name
    siaf : pysiaf.Siaf
        SIAF instance for a single instrument
    verbose : bool
        Verbose output on/off

    Returns
    -------
    x_sci, y_sci : tuple of numpy arrays
        Subarray corner coordinates

    """
    logger = logging.getLogger('mirage.utils.get_siaf_information.sci_subarray_corners')

    # get SIAF
    if siaf is None:
        siaf = get_instance(instrument)

    # get master aperture names
    siaf_detector_layout = iando.read.read_siaf_detector_layout()
    master_aperture_names = siaf_detector_layout['AperName'].data

    # read pysiaf aperture definition file containing DetRef and SciRef values
    siaf_aperture_definitions = iando.read.read_siaf_aperture_definitions(instrument)

    # aperture object
    aperture = siaf[aperture_name]

    # aperture corners in SIAF detector coordinates
    x_det, y_det = aperture.corners('det', rederive=True)

    # determine parent aperture, i.e. the underlying full frame SCA aperture
    index = siaf_aperture_definitions['AperName'].tolist().index(aperture_name)
    aperture._parent_apertures = siaf_aperture_definitions['parent_apertures'][index]

    # If multiuple apertures are listed as parents keep only the first
    if ';' in aperture._parent_apertures:
        logger.info('Multiple parent apertures: {}'.format(aperture._parent_apertures))
        aperture._parent_apertures = aperture._parent_apertures.split(';')[0]

    if aperture_name in master_aperture_names:
        # if master aperture, use it directly to transform to science frame
        x_sci, y_sci = aperture.det_to_sci(x_det, y_det)
    elif aperture._parent_apertures is not None:
        # use parent aperture for transformation
        if verbose:
            logger.info('Using parent {} for {}'.format(aperture._parent_apertures, aperture_name))
        x_sci, y_sci = siaf[aperture._parent_apertures].det_to_sci(x_det, y_det)
        aperture = siaf[aperture._parent_apertures]

    if instrument.lower() == 'nircam':
        if aperture.DetSciParity == 1:
            corner_index = np.array([1, 3])
        elif aperture.DetSciParity == -1:
            # NIRCam will always fall in here, except in case of non-dms orientation
            corner_index = np.array([0, 2])
        x_corner = x_sci[corner_index]
        y_corner = y_sci[corner_index]
    elif instrument.lower() == 'niriss':
        x_corner_index = np.array([0, 2])
        y_corner_index = np.array([0, 2])
        if aperture_name == 'NIS_CEN_OSS':
            x_corner_index = np.array([1, 3])
            y_corner_index = np.array([3, 1])
        x_corner = x_sci[x_corner_index]
        y_corner = y_sci[y_corner_index]
        if aperture_name in ['NIS_SUBSTRIP96', 'NIS_SUBSTRIP256']:
            x_corner = [1, 2048]
            if aperture_name == 'NIS_SUBSTRIP96':
                y_corner = [1803, 1898]
            if aperture_name == 'NIS_SUBSTRIP256':
                y_corner = [1793, 2048]
    elif instrument.lower() == 'fgs':
        x_corner_index = np.array([0, 2])
        y_corner_index = np.array([0, 2])
        if aperture_name == 'FGS1_FULL_OSS':
            x_corner_index = np.array([1, 3])
            y_corner_index = np.array([3, 1])
        if aperture_name == 'FGS2_FULL_OSS':
            x_corner_index = np.array([1, 3])
            y_corner_index = np.array([1, 3])
        x_corner = x_sci[x_corner_index]
        y_corner = y_sci[y_corner_index]
    else:
        raise NotImplementedError(("Instrument {} not supported for SIAF subarray corners"
                                   .format(instrument)))

    # account for mirage conventions (e.g. 0-based indexing)
    # we also want integer values as these will be indexes
    x_corner = np.array([np.ceil(x_corner[0]) - 1, np.floor(x_corner[1]) - 1])
    y_corner = np.array([np.ceil(y_corner[0]) - 1, np.floor(y_corner[1]) - 1])
    return x_corner.astype(int), y_corner.astype(int)
