"""Test the functions provided by mirage.yaml.yaml_generator

Authors
-------
    - Lauren Chambers

Use
---
    >>> pytest test_yaml_generator.py
"""
import os

import numpy as np
import pytest

from mirage.yaml.yaml_generator import SimInput

# Define directory and file locations
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

# Reset the MIRAGE_DATA env variable to be a real location so yaml_generator
# doesn't croak
os.environ["MIRAGE_DATA"] = __location__


def test_get_psf_path():
    """Test that the get_psf_path method of the yaml_generator.SimInput class
    is working as expected.
    """
    # Make an instance of the SimInput class
    input_xml = os.path.join(__location__, 'test_data', 'NIRCam', '1144-OTE-10.xml')
    pointing_file = os.path.join(__location__, 'test_data', 'NIRCam', '1144-OTE-10.pointing')
    yam = SimInput(input_xml, pointing_file, offline=True)

    # Fake the activity IDs and instrument
    n_activities = 101
    act_ids = [np.base_repr(i, 36).zfill(2)  for i in range(n_activities)]
    yam.info = {}
    yam.info['act_id'] = act_ids
    yam.info['Instrument'] = ['NIRCam'] * n_activities

    # Test for a default path
    paths_out = yam.get_psf_path()
    assert len(paths_out) == n_activities,\
        'Default PSF path not properly provided.'
    assert 'nircam/gridded_psf_library' in paths_out[0],\
        'Default PSF path not properly provided.'
    np.testing.assert_array_equal(paths_out, paths_out[0],
                                  err_msg='Default PSF path not properly defined.')

    # Test for a single path
    yam.psf_paths = __location__
    paths_out = yam.get_psf_path()
    assert paths_out == [__location__] * n_activities,\
        'Single PSF path not properly assigned.'

    # Test for a list of paths of incorrect length
    with pytest.raises(ValueError) as e:
        yam.psf_paths = [__location__] * 99
        yam.get_psf_path()
        assert 'Please provide the psf_paths in the form of a list of strings ' \
               'with a length equal to the number of activities in the APT ' \
               'program (101), not equal to 99.' in e, \
            'Failed to reject psf_path of incorrect length.'

    # Test for a list of paths of correct length
    list_101_paths = [__location__] * 50 + [os.path.dirname(__location__)] * 51
    yam.psf_paths = list_101_paths
    paths_out = yam.get_psf_path()
    assert paths_out == sorted(list_101_paths),\
        'List of PSF paths not properly assigned.'

    # Test for a completely invalid path
    with pytest.raises(TypeError) as e:
        yam.psf_paths = 3.054
        yam.get_psf_path()
        assert 'Please provide the psf_paths in the form of a list or ' \
               'string, not float' in e, \
            'Failed to reject psf_path of incorrect type.'
