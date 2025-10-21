"""Unit tests for flux calibration/zeropoints

Authors
-------
    - Bryan Hilbert

Use
---
    Ensure you have pytest installed. Then, simply run pytest in any
    parent directory of mirage/tests/:
    >>> pytest
"""
from astropy.table import Table
import importlib.resources as resources
import numpy as np
import os

from mirage.utils import flux_cal
from mirage.utils.constants import MODULE_PATH

CONFIG_DIR = os.path.join(MODULE_PATH, 'config')


def test_add_detector_to_zeropoints():
    """Test addition of column to table
    """
    detector = 'NRCA1'
    tab = Table()
    tab['index'] = np.arange(5)
    tab['information'] = [1.2, 2.3, 3.4, 4.5, 5.6]

    updated_tab = flux_cal.add_detector_to_zeropoints(detector, tab)
    assert np.all(updated_tab['Detector'].data == np.array([detector] * 5))
    assert np.all(updated_tab['index'].data == tab['index'].data)
    assert np.all(updated_tab['information'].data == tab['information'].data)


def test_fluxcal_info():
    """Test that zeropoint information for the exposure is correctly retrieved
    """
    params = {'Inst': {"instrument": 'NIRCAM'},
              'Readout': {'filter': 'F200W', 'pupil': 'CLEAR'},
              'Reffiles': {'flux_cal': os.path.join(CONFIG_DIR, 'NIRCam_zeropoints.list')}
              }

    detector = 'NRCA1'
    module = 'A'
    vegazp, photflam, photfnu, pivot = flux_cal.fluxcal_info(params['Reffiles']['flux_cal'], 'NIRCAM',
                                                             params['Readout']['filter'],
                                                             params['Readout']['pupil'], detector, module)

    assert vegazp == 25.6
    assert photflam == 3.3443932204899414e-21
    assert photfnu == 4.4058412e-31
    assert pivot == 1.988
