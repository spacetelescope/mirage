"""System test of mirage/NIRISS that creates simulations based on a .yaml file.

Authors
-------
    - Johannes Sahlmann

Use
---
    >>> pytest -s test_niriss_imaging.py


"""

import os
import pytest

from mirage import imaging_simulator as im

os.environ['TEST_NIRISS_DATA'] = os.path.join(os.path.dirname(__file__), 'test_data/NIRISS')

# Determine if tests are being run on Github Actions CI
ON_GITHUB = '/home/runner' in os.path.expanduser('~')


@pytest.mark.skipif(ON_GITHUB,
                   reason="Cannot access mirage data in the central storage directory from Github Actions CI.")
def test_niriss_imaging():
    m = im.ImgSim(offline=True)
    m.paramfile = os.path.join(os.path.dirname(__file__), 'test_data/NIRISS/niriss_imaging_example.yaml')
    m.create()
