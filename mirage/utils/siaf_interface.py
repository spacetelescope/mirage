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

    # account for mirage conventions (e.g. 0-based indexing)
    x_sci -= 1.5
    y_sci -= 1.5
    if instrument == 'NIRCam':
        if aperture.DetSciParity == 1:
            corner_index = np.array([1, 3])
        elif aperture.DetSciParity == -1:
            corner_index = np.array([0, 2])
    else:
        raise NotImplementedError

    return x_sci[corner_index], y_sci[corner_index]
