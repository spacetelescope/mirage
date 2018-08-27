"""Test the functions provided by siaf_interface.

Authors
-------
    - Johannes Sahlmann

Use
---
    >>> pytest -s test_siaf_interface.py


"""
from pysiaf import iando

from mirage.utils import siaf_interface


def test_sci_subarray_corners():
    """Unit test for siaf_interface.sci_subarray_corners."""
    instrument = 'NIRCam'
    siaf_detector_layout = iando.read.read_siaf_detector_layout()
    master_aperture_names = siaf_detector_layout['AperName'].data
    for aperture_name in master_aperture_names:
        if 'NRC' not in aperture_name:
            continue
        x_sci, y_sci = siaf_interface.sci_subarray_corners(instrument, aperture_name)
        assert len(x_sci) == 2
        assert len(y_sci) == 2
