"""System test of mirage/NIRISS that creates SOSS simulations based on a .yaml file.

Authors
-------
    - Joe Filippazzo

Use
---
    >>> pytest -s test_niriss_soss.py


"""

import os
import pytest

from mirage import soss_simulator as ss

os.environ['TEST_NIRISS_DATA'] = os.path.join(os.path.dirname(__file__), 'test_data/NIRISS')

# Determine if tests are being run on Travis
ON_TRAVIS =  'travis' in os.path.expanduser('~')


@pytest.mark.skipif(ON_TRAVIS, reason="Cannot access mirage data in the central storage directory from Travis CI.")
def test_niriss_soss_substrip256_clear():
    m = ss.SossSim(offline=True)
    m.paramfile = os.path.join(os.path.dirname(__file__), 'test_data/NIRISS/niriss_soss_substrip256_clear.yaml')
    m.create()

@pytest.mark.skipif(ON_TRAVIS, reason="Cannot access mirage data in the central storage directory from Travis CI.")
def test_niriss_soss_substrip256_f277w():
    m = ss.SossSim(offline=True)
    m.paramfile = os.path.join(os.path.dirname(__file__), 'test_data/NIRISS/niriss_soss_substrip256_f277w.yaml')
    m.create()

@pytest.mark.skipif(ON_TRAVIS, reason="Cannot access mirage data in the central storage directory from Travis CI.")
def test_niriss_soss_substrip96_clear():
    m = ss.SossSim(offline=True)
    m.paramfile = os.path.join(os.path.dirname(__file__), 'test_data/NIRISS/niriss_soss_substrip96_clear.yaml')
    m.create()

@pytest.mark.skipif(ON_TRAVIS, reason="Cannot access mirage data in the central storage directory from Travis CI.")
def test_niriss_soss_substrip96_f277w():
    m = ss.SossSim(offline=True)
    m.paramfile = os.path.join(os.path.dirname(__file__), 'test_data/NIRISS/niriss_soss_substrip96_f277w.yaml')
    m.create()

@pytest.mark.skipif(ON_TRAVIS, reason="Cannot access mirage data in the central storage directory from Travis CI.")
def test_niriss_soss_fullframe_clear():
    m = ss.SossSim(offline=True)
    m.paramfile = os.path.join(os.path.dirname(__file__), 'test_data/NIRISS/niriss_soss_fullframe_clear.yaml')
    m.create()

@pytest.mark.skipif(ON_TRAVIS, reason="Cannot access mirage data in the central storage directory from Travis CI.")
def test_niriss_soss_fullframe_f277w():
    m = ss.SossSim(offline=True)
    m.paramfile = os.path.join(os.path.dirname(__file__), 'test_data/NIRISS/niriss_soss_fullframe_f277w.yaml')
    m.create()

