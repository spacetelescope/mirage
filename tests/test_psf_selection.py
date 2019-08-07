"""Test the functions provided by mirage.psf.psf_selection

Authors
-------
    - Lauren Chambers

Use
---
    >>> pytest test_psf_selection.py
"""
import os
import shutil

from astropy.io import fits
import photutils
import pytest

from mirage.psf.psf_selection import get_library_file, _load_itm_library
from mirage.utils.utils import ensure_dir_exists

# Define directory and file locations
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
TEMP_TEST_DIRECTORY = os.path.join(__location__, 'temp_data', 'test_psf_selection')


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


@pytest.fixture(scope="module")
def test_segment_psf_library_file(test_directory):
    """Download and save locally a segment PSF library file from Box for testing

    Yields
    -------
    test_lib_filename : str
        Path to test segment PSF library used for testing
    """
    # Download file and yield its name
    test_lib_filename = os.path.join(test_directory,
                                     'nircam_nrca3_f212n_fovp1024_samp1_npsf1_seg12.fits')
    url = "https://stsci.box.com/shared/static/4c0em1yhb1qsvrpku7j0jw1tztucvkih.fits"
    with fits.open(url) as hdulist:
        hdulist.writeto(test_lib_filename)
        yield test_lib_filename


@pytest.fixture(scope="module")
def itm_file(test_directory):
    """Download and save locally an ITM image file from Box for testing

    Yields
    -------
    itm_filename : str
        Path to ITM image file used for testing
    """
    # Download file and yield its name
    itm_filename = os.path.join(test_directory,
                                     'itm_file.fits')
    url = "https://stsci.box.com/shared/static/lsvb6acrlubwnlzz7wkkpbjgm01bbzgb.fits"
    with fits.open(url) as hdulist:
        hdulist.writeto(itm_filename)
        yield itm_filename


def test_get_segment_psf_library_file(test_segment_psf_library_file):
    """Test the identification of a segment PSF library file

    Parameters
    ----------
    test_segment_psf_library_file : str
        Path to test segment PSF library used for testing
    """
    instrument = "NIRCam"
    detector = "NRCA3"
    filt = "F212N"
    pupil = "CLEAR"
    wfe = ''
    wfe_group = 0
    library_path = os.path.dirname(test_segment_psf_library_file)
    segment_id = 12

    # Ensure the right file is found given the right parameters
    match_file = get_library_file(instrument, detector, filt, pupil, wfe, wfe_group,
                     library_path, wings=False, segment_id=segment_id)
    assert match_file == test_segment_psf_library_file, \
        "Failed to match segment PSF library"

    # Ensure an error is raised if no file exists that matches
    with pytest.raises(ValueError) as e:
        detector = "NRCB5"
        segment_id = 11
        get_library_file(instrument, detector, filt, pupil, wfe, wfe_group,
                         library_path, wings=False, segment_id=segment_id)
        assert "No PSF library file found matching requested parameters" in e


def test_load_itm_library(itm_file):
    """Test the loading of a GriddedPSFModel object from an ITM FITS file

    Parameters
    ----------
    itm_filename : str
        Path to ITM image file used for testing
    """
    lib_model = _load_itm_library(itm_file)

    assert isinstance(lib_model, photutils.psf.models.GriddedPSFModel), \
        'ITM PSF library not created correctly'
    assert lib_model.grid_xypos == [(1023.0, 1023.0)], \
        'ITM PSF library not created correctly'
    assert lib_model.oversampling == 1, \
        'ITM PSF library not created correctly'
    assert lib_model.data.shape == (1, 2048, 2048), \
        'ITM PSF library not created correctly'
