"""System test of mirage/NIRCam that creates simulations based on a .yaml file.

Authors
-------
    - Johannes Sahlmann

Use
---
    >>> pytest -s test_nircam_imaging.py


"""

import os
import pytest

from mirage.scripts import imaging_simulator as im

# os.environ['MIRAGE_DATA'] = ''
os.environ['TEST_DATA'] = os.path.join(os.path.dirname(__file__), 'test_data/NIRCam')

@pytest.mark.xfail
def test_nircam_imaging():
    m = im.ImgSim()
    m.paramfile = os.path.join(os.path.dirname(__file__), 'test_data/NIRCam/nircam_imaging_example.yaml')
    m.create()
