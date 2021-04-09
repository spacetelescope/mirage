"""Test the functions provided by mirage.psf.deployments

Authors
-------
    - Lauren Chambers

Use
---
    >>> pytest test_deployments.py
"""
import glob
import os
import shutil
import sys

import numpy as np
import pytest
import webbpsf

from .utils import parametrized_data
from mirage.psf.deployments import generate_random_ote_deployment, load_ote_from_deployment_yaml
from mirage.utils.utils import ensure_dir_exists

# Define directory and file locations
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
TEMP_TEST_DIRECTORY = os.path.join(__location__, 'temp_data', 'test_deployments')

# Load parametrized data
PARAMETRIZED_DATA = parametrized_data()['test_deployments']

# Determine if tests are being run on Github Actions CI
ON_GITHUB = '/home/runner' in os.path.expanduser('~')

# Determine the version of python used. For python 3.8 and above
# webbpsf is installed via pip, which means the data files will not
# be accessible and any test that relies on webbpsf should be skipped
python_version = sys.version[0:3]
testable_versions = ['3.6', '3.7']
skip_versions = ['3.8', '3.9']
if python_version in skip_versions:
    SKIP_WEBBPSF = True
else:
    SKIP_WEBBPSF = False


@pytest.fixture(scope="module")
def test_directory(test_dir=TEMP_TEST_DIRECTORY):
    """Create a test directory.

    Parameters
    ----------
    test_dir : str
        Path to directory used for testing

    Yields
    -------
    test_dir : str
        Path to directory used for testing
    """
    # Create directory and yield its name
    ensure_dir_exists(test_dir)  # creates directory with default mode=511
    yield test_dir

    # Remove directory
    if os.path.isdir(test_dir):
        shutil.rmtree(test_dir)


@pytest.fixture(scope="function")
def remove_yamls_and_fits(test_directory):
    """Deletes YAML and FITS files produced by an OTE deployment function.

    Parameters
    ----------
    test_directory : str
        Path to directory used for testing
    """

    yaml_search_path = os.path.join(test_directory, 'deployment_errors*.yaml')
    yaml_files = glob.glob(yaml_search_path)
    opd_search_path = os.path.join(test_directory, 'OPD*.fits')
    opd_files = glob.glob(opd_search_path)
    yield

    print("Removing OPD yamls and FITS files")
    for file in yaml_files + opd_files:
        os.remove(file)


@pytest.mark.skipif((ON_GITHUB and SKIP_WEBBPSF), reason='Webbpsf data files cannot be downloaded via pip')
def test_generate_random_ote_deployment(test_directory, remove_yamls_and_fits):
    """Test the creation of a WebbPSF adjustable OTE object representing a
    perturbed OTE mirror state by randomly generating mirror deployment errors

    Parameters
    ----------
    test_directory : str
        Path to directory used for testing
    remove_yamls_and_fits
        Pytest fixture to clear files after test is run
    """
    # Make the OTE randomly
    ote, segment_tilts, ote_opd_with_tilts = generate_random_ote_deployment(
        test_directory, reduction_factor=0.2, save=True
    )

    # Ensure the OPD yaml and fits files were saved
    yaml_search_path = os.path.join(test_directory, 'deployment_errors*.yaml')
    yaml_files = glob.glob(yaml_search_path)
    opd_search_path = os.path.join(test_directory, 'OPD*.fits')
    opd_files = glob.glob(opd_search_path)
    assert len(yaml_files) >=2, 'Failed to generate randomly deployed OTE'
    assert len(opd_files) >=2, 'Failed to generate randomly deployed OTE'

    # Ensure the resultant OTE is as expected
    assert isinstance(ote, webbpsf.opds.OTE_Linear_Model_WSS), 'Invalid format for OTE.'
    assert ote.shape == (1024, 1024), 'Invalid format for OTE.'
    assert ote.segment_state.shape == (19, 6), 'Invalid format for OTE.'
    check_removed_tip_tilt = ote.segment_state[:-1, :3] == 0.
    assert check_removed_tip_tilt.all(), 'Segment tilts were not removed.'
    check_perturbed = ote.segment_state[:-1, 3:] != 0.
    assert check_perturbed.all(), 'Segment tilts were not applied.'

    # Ensure the segment tilts array is as expected
    assert segment_tilts.shape == (19, 2), 'Invalid format for segment tilts.'
    check_perturbed_tilts = segment_tilts != 0.
    assert check_perturbed_tilts.all(), 'Segment tilts were not applied.'

    # Ensure the OTE array with tilts preserved is as expected
    assert ote_opd_with_tilts.shape == (1024, 1024), 'Invalid format for OTE.'
    check_different_opd = ote_opd_with_tilts == ote.opd
    assert not check_different_opd.all(), 'Segment tilts were not removed.'


@pytest.mark.skipif((ON_GITHUB and SKIP_WEBBPSF), reason='Webbpsf data files cannot be downloaded via pip')
def test_load_ote_from_deployment_yaml(test_directory, remove_yamls_and_fits):
    """Test the creation of a WebbPSF adjustable OTE object representing a
    perturbed OTE mirror state by loading from a YAML file.

    Parameters
    ----------
    test_directory : str
        Path to directory used for testing
    remove_yamls_and_fits
        Pytest fixture to clear files after test is run
    """
    # Load known data for comparison
    test_data = PARAMETRIZED_DATA['test_load_ote_from_deployment_yaml']

    # Load the OTE from a YAML file
    deployments_file = os.path.join(__location__, 'test_data', 'deployment_errors_test.yaml')
    ote, segment_tilts, ote_opd_with_tilts = load_ote_from_deployment_yaml(
        deployments_file, test_directory, save=True
    )

    # Ensure the OPD yaml and fits files were saved
    yaml_search_path = os.path.join(test_directory, 'deployment_errors*.yaml')
    yaml_files = glob.glob(yaml_search_path)
    opd_search_path = os.path.join(test_directory, 'OPD*.fits')
    opd_files = glob.glob(opd_search_path)

    #assert len(yaml_files) >= 2, 'Failed to generate randomly deployed OTE'
    assert len(opd_files) >= 2, 'Failed to generate randomly deployed OTE'

    # Ensure the resultant OTE is as expected
    assert isinstance(ote, webbpsf.opds.OTE_Linear_Model_WSS), 'Invalid format for OTE.'
    assert ote.shape == (1024, 1024), 'Invalid format for OTE.'
    assert ote.segment_state.shape == (19, 6), 'Invalid format for OTE.'
    np.testing.assert_array_equal(ote.segment_state[:-1, :3], 0.,
                                  err_msg='Segment tilts were not removed.')
    assert not np.array_equiv(ote.segment_state[:-1, 3:], 0.), \
        'Segment tilts were not applied.'

    assert np.isclose(np.mean(ote.opd), -2.1281188008526972e-08), \
        'OTE was not loaded correctly from file.'
    np.testing.assert_allclose(ote.segment_state[:-1, 3:],
                               np.array(test_data['segment_state_no_tilts']),
                               err_msg="OTE was not loaded correctly from file.")


    # Ensure the segment tilts array is as expected
    assert segment_tilts.shape == (19, 2), 'Invalid format for segment tilts.'
    assert not np.array_equiv(segment_tilts, 0.), 'Segment tilts were not applied.'
    np.testing.assert_allclose(segment_tilts,
                               test_data['segment_tilts'],
                               err_msg="OTE was not loaded correctly from file.")

    # Ensure the OTE array with tilts preserved is as expected
    assert ote_opd_with_tilts.shape == (1024, 1024), 'Invalid format for OTE.'
    assert not np.array_equal(ote_opd_with_tilts, ote.opd), 'Segment tilts were not removed.'
    assert np.isclose(np.mean(ote_opd_with_tilts), 1.7097442327942417e-05), \
        'OTE was not loaded correctly from file.'
