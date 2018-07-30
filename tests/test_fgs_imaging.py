"""System test of mirage/FGS that creates simulations based on a .yaml file.

Authors
-------
    - Johannes Sahlmann

Use
---
    >>> pytest -s test_fgs_imaging.py


"""

import os
import pytest

from mirage import imaging_simulator as im

# os.environ['MIRAGE_DATA'] = ''
os.environ['TEST_DATA'] = os.path.join(os.path.dirname(__file__), 'test_data/FGS')

@pytest.mark.xfail
def test_fgs_imaging():
    m = im.ImgSim()
    m.paramfile = os.path.join(os.path.dirname(__file__), 'test_data/FGS/fgs_imaging_example.yaml')
    m.create()

