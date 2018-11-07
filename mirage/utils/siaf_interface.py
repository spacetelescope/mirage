"""SIAF interface module to support accessing SIAF information.

This module provides ``mirage`` with functions to interface SIAF content via the pysiaf module.

Authors
-------

    - Johannes Sahlmann

Use
---

    This module can be imported and used with

    ::

        from mirage.utils import siaf_interface

"""
import os
import numpy as np

import pysiaf
from pysiaf import iando
from pysiaf.constants import JWST_DELIVERY_DATA_ROOT
from ..utils import rotations
from ..utils import set_telescope_pointing_separated as set_telescope_pointing


def get_siaf_information(instrument, aperture_name, ra, dec, telescope_roll, v2_arcsec=None, v3_arcsec=None, verbose=False):
    """Use pysiaf to get aperture information.

    Parameters
    ----------
    instrument : str
        Instrument name.

    aperture_name : str
        Aperture name (e.g. "NRCA1_FULL")
    """
    # Temporary fix to access good NIRCam distortion coefficients which
    # which are not yet in the PRD
    if instrument.lower() == 'nircam':
        print("NOTE: Using pre-delivery SIAF data for {}".format(aperture_name))
        1/0
        if instrument == 'NIRCAM':
            instrument = 'NIRCam'
        pre_delivery_dir = os.path.join(JWST_DELIVERY_DATA_ROOT, 'NIRCam')
        siaf = pysiaf.Siaf(instrument, basepath=pre_delivery_dir)[aperture_name]
    else:
        siaf = pysiaf.Siaf(instrument)[aperture_name]

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
        xcorner, ycorner = sci_subarray_corners(instrument, aperture_name)
        subarray_boundaries = [xcorner[0], ycorner[0], xcorner[1], ycorner[1]]
    except (RuntimeError, TypeError) as e: # e.g. NIRSpec NRS_FULL_MSA aperture
        if verbose:
            print('get_siaf_information raised error:\n{}\nIgnoring it.'.format(e))
        subarray_boundaries = [0, 0, 0, 0]
    return siaf, local_roll, att_matrix, fullframesize, subarray_boundaries


def sci_subarray_corners(instrument, aperture_name, verbose=False):
    """Return the two opposing aperture corners in the SIAF Science frame of the full-frame SCA.

    This function serves as interface between the SIAF information accessible via the pysiaf package
    and the subarray information formatted for use by mirage.

    Parameters
    ----------
    instrument : str
        JWST instrument name with correct capitalization
    aperture_name : str
        SIAF aperture name
    verbose : bool
        Verbose output on/off

    Returns
    -------
    x_sci, y_sci : tuple of numpy arrays
        Subarray corner coordinates

    """
    # get SIAF
    siaf = pysiaf.Siaf(instrument)

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

    if aperture_name in master_aperture_names:
        # if master aperture, use it directly to transform to science frame
        x_sci, y_sci = aperture.det_to_sci(x_det, y_det)
    elif aperture._parent_apertures is not None:
        # use parent aperture for transformation
        if verbose:
            print('Using parent {} for {}'.format(aperture._parent_apertures, aperture_name))
        x_sci, y_sci = siaf[aperture._parent_apertures].det_to_sci(x_det, y_det)
        aperture = siaf[aperture._parent_apertures]

    if instrument == 'NIRCam':
        if aperture.DetSciParity == 1:
            corner_index = np.array([1, 3])
        elif aperture.DetSciParity == -1:
            # NIRCam will always fall in here, except in case of non-dms orientation
            corner_index = np.array([0, 2])
        x_corner = x_sci[corner_index]
        y_corner = y_sci[corner_index]
    elif instrument == 'NIRISS':
        x_corner_index = np.array([0, 2])
        y_corner_index = np.array([0, 2])
        if aperture_name == 'NIS_CEN_OSS':
            x_corner_index = np.array([1, 3])
            y_corner_index = np.array([3, 1])
        x_corner = x_sci[x_corner_index]
        y_corner = y_sci[y_corner_index]
        if aperture_name in ['NIS_SUBSTRIP96', 'NIS_SUBSTRIP256']:
            x_corner = [1, 2048]
            y_corner = [1, 2048]
    elif instrument == 'FGS':
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
    return x_corner.astype(np.int), y_corner.astype(np.int)
