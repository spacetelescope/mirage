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

import numpy as np

import pysiaf
from pysiaf import iando
from ..utils import rotations
from ..utils import set_telescope_pointing_separated as set_telescope_pointing


def get_siaf_information(instrument, aperture, ra, dec, telescope_roll):
    """Use pysiaf to get aperture information

    Parameters:
    -----------
    instrument : str
        Instrument name.

    aperture : str
        Aperture name (e.g. "NRCA1_FULL")

    Returns:
    --------
    None
    """
    siaf = pysiaf.Siaf(instrument)[aperture]
    local_roll = set_telescope_pointing.compute_local_roll(telescope_roll,
                                                           ra, dec, siaf.V2Ref,
                                                           siaf.V3Ref)
    # Create attitude_matrix
    att_matrix = rotations.attitude(siaf.V2Ref, siaf.V3Ref, ra, dec, local_roll)

    # Get full frame size
    fullframesize = siaf.XDetSize

    # Subarray boundaries in full frame coordinates
    xcorner, ycorner = sci_subarray_corners(instrument, aperture)
    subarray_boundaries = [xcorner[0], ycorner[0], xcorner[1], ycorner[1]]
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
            corner_index = np.array([0, 2])
    else:
        raise NotImplementedError

    x_corner = x_sci[corner_index]
    y_corner = y_sci[corner_index]

    # account for mirage conventions (e.g. 0-based indexing)
    # we also want integer values as these will be indexes
    x_corner = np.array([np.ceil(x_corner[0]) - 1, np.floor(x_corner[1]) - 1])
    y_corner = np.array([np.ceil(y_corner[0]) - 1, np.floor(y_corner[1]) - 1])
    return x_corner.astype(np.int), y_corner.astype(np.int)
