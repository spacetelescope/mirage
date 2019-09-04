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

os.environ['TEST_FGS_DATA'] = os.path.join(os.path.dirname(__file__), 'test_data/FGS')

# Determine if tests are being run on Travis
ON_TRAVIS = 'travis' in os.path.expanduser('~')


@pytest.mark.skipif(ON_TRAVIS,
                    reason="Cannot access mirage data in the central storage directory from Travis CI.")
def test_fgs_imaging():
    m = im.ImgSim(offline=True)
    m.paramfile = os.path.join(os.path.dirname(__file__), 'test_data/FGS/fgs_imaging_example.yaml')
    m.create()
