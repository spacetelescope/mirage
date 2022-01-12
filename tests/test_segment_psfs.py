"""Test the functions provided by mirage.psf.segment_psfs

Authors
-------
    - Lauren Chambers

Use
---
    >>> pytest test_segment_psfs.py
"""
from ast import literal_eval
import glob
import os
import shutil
import sys

from astropy.io import fits
import numpy as np
import photutils
import pytest
from webbpsf.utils import to_griddedpsfmodel

from .utils import parametrized_data
from mirage.psf.deployments import generate_random_ote_deployment
from mirage.psf.psf_selection import get_library_file
from mirage.psf.segment_psfs import (get_segment_library_list, get_segment_offset,
                                     get_gridded_segment_psf_library_list, generate_segment_psfs)
from mirage.utils.utils import ensure_dir_exists

# Define directory and file locations
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
TEMP_TEST_DIRECTORY = os.path.join(__location__, 'temp_data', 'test_segment_psfs')
TEST_DATA_DIR = os.path.expandvars("$MIRAGE_DATA/test_data/test_segment_psfs")

# Load parametrized data
PARAMETRIZED_DATA = parametrized_data()['test_segment_psfs']
SUCCESS_GENERATE_SEGMENT_PSF = None

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

# Define default inputs
INSTRUMENT = 'NIRCam'
DETECTOR = 'NRCA3'
FILTER = 'F212N'


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
def test_library_file(test_directory):
    """Download and save locally a segment PSF library file from Box for testing

    Yields
    -------
    test_lib_filename : str
        Path to test library used for testing
    """
    # Download file and yield its name
    test_lib_filename = os.path.join(test_directory, 'test_library',
                                     'nircam_nrca3_f212n_fovp1024_samp1_npsf1_seg12.fits')
    ensure_dir_exists(os.path.dirname(test_lib_filename))
    url = "https://stsci.box.com/shared/static/4c0em1yhb1qsvrpku7j0jw1tztucvkih.fits"
    with fits.open(url) as hdulist:
        hdulist.writeto(test_lib_filename)
        yield test_lib_filename


@pytest.mark.skipif((ON_GITHUB and SKIP_WEBBPSF), reason='Webbpsf data files cannot be downloaded via pip')
def test_generate_segment_psfs(test_directory):
    """Test the generation of 18 segment PSFs with a randomly-perturbed
    mirror state.

    Parameters
    ----------
    test_directory : str
        Path to directory used for testing
    """
    test_data = PARAMETRIZED_DATA['test_generate_segment_psfs']

    # Randomly generate OTE state with deployment errors
    ote, segment_tilts, ote_opd_with_tilts = generate_random_ote_deployment(
        test_directory, reduction_factor=0.2
    )
    # Ensure the OPD yaml and fits files were saved
    yaml_search_path = os.path.join(test_directory, 'deployment_errors*.yaml')
    opd_search_path = os.path.join(test_directory, 'OPD*.fits')
    assert len(glob.glob(yaml_search_path)) >= 2, 'Failed to generate randomly deployed OTE'
    assert len(glob.glob(opd_search_path)) >= 2, 'Failed to generate randomly deployed OTE'

    # Ensure it doesn't work with mismatched filter & detector
    with pytest.raises(ValueError) as e:
        generate_segment_psfs(ote, segment_tilts, test_directory, filters=['F212N'],
                              detectors='NRCA5', fov_pixels=10, overwrite=False)
        assert "No matching filters and detectors given" in e, \
            'Failed to catch mismatched filter and detector'

    # Make 18 segment PSF library files
    generate_segment_psfs(ote, segment_tilts, test_directory, filters=['F212N'],
                          detectors='NRCA3', fov_pixels=101, overwrite=False)

    # Ensure that all 18 segment PSF library files were created
    fits_search_path = os.path.join(test_directory, 'nircam*.fits')
    library_names = test_data

    global SUCCESS_GENERATE_SEGMENT_PSF
    SUCCESS_GENERATE_SEGMENT_PSF = True

    for name in library_names:
        lib_success = os.path.join(test_directory, name) in glob.glob(fits_search_path)
        SUCCESS_GENERATE_SEGMENT_PSF = SUCCESS_GENERATE_SEGMENT_PSF & lib_success
        assert lib_success, 'Failed to create file: {}'.format(os.path.join(test_directory, name))


@pytest.mark.skipif(ON_GITHUB,
                   reason="Cannot access mirage data in the central storage directory from Github Actions CI.")
def test_get_segment_library_list_remote():
    """Test construction of segment PSF libraries using data on the central
    storage directory.
    """
    # Get a list of segment PSF library files
    library_list = get_segment_library_list(INSTRUMENT, DETECTOR, FILTER, TEST_DATA_DIR)

    # Ensure the right files were found
    library_names = PARAMETRIZED_DATA['test_get_segment_library_list']
    assert len(library_list) == 18, 'Did not find all 18 segment libraries'
    assert sorted(library_list) == sorted([os.path.join(TEST_DATA_DIR, name) for name in library_names]), \
        'Did not find all expected 18 segment libraries'


def test_get_segment_library_list_local(test_directory):
    """Test construction of segment PSF libraries using data generated by
    test_generate_segment_psfs().

    Parameters
    ----------
    test_directory : str
        Path to directory used for testing
    """
    # Skip if test_generate_segment_psfs failed
    if not SUCCESS_GENERATE_SEGMENT_PSF:
        pytest.skip("Cannot load library files generated in test_generate_segment_psfs.")

    # Get a list of segment PSF library files
    library_list = get_segment_library_list(INSTRUMENT, DETECTOR, FILTER, test_directory)

    # Ensure the right files were found
    library_names = PARAMETRIZED_DATA['test_generate_segment_psfs']
    assert len(library_list) == 18, 'Did not find all 18 segment libraries'
    assert sorted(library_list) == sorted([os.path.join(test_directory, name) for name in library_names]), \
        'Did not find all expected 18 segment libraries'


test_data = PARAMETRIZED_DATA['test_get_segment_offset_remote']
get_segment_offset_remote_parameters = []
for i, tuple_string in enumerate(test_data):
    offset_tuple =   literal_eval(tuple_string)
    get_segment_offset_remote_parameters.append([i+1, offset_tuple])
@pytest.mark.parametrize('segment_number, correct_offset', get_segment_offset_remote_parameters)
@pytest.mark.skipif(ON_GITHUB,
                   reason="Cannot access mirage data in the central storage directory from Github Actions CI.")
def test_get_segment_offset_remote(segment_number, correct_offset):
    """Test the extraction of segment offsets from segment PSF libraries on the
    central storage directory.

    Parameters
    ----------
    segment_number : int
        Segment ID, i.e 3
    correct_offset : tuple
        Expected x and y offset in arcseconds
    """
    library_path = TEST_DATA_DIR
    library_list = get_segment_library_list(INSTRUMENT, DETECTOR, FILTER, library_path)
    print('CPL:',segment_number, DETECTOR)
    x_arcsec, y_arcsec = get_segment_offset(segment_number, DETECTOR, library_list)
    assert (x_arcsec, y_arcsec) == correct_offset, 'Incorrect conversion of segment offsets'


def test_get_segment_offset_local_created(test_directory):
    """Test the extraction of segment offsets from segment PSF libraries
    generated by test_generate_segment_psfs().

    Parameters
    ----------
    test_directory : str
        Path to directory used for testing
    """
    # Skip if test_generate_segment_psfs failed
    if not SUCCESS_GENERATE_SEGMENT_PSF:
        pytest.skip("Cannot load library files generated in test_generate_segment_psfs.")
    library_path = test_directory
    library_list = get_segment_library_list(INSTRUMENT, DETECTOR, FILTER, library_path)

    # Just try to get one offset and make sure it doesn't croak
    x_arcsec, y_arcsec = get_segment_offset(1, DETECTOR, library_list)
    assert type(x_arcsec) == np.float64, 'Incorrect conversion of segment offsets'
    assert type(y_arcsec) == np.float64, 'Incorrect conversion of segment offsets'


def test_get_segment_offset_local_stored(test_library_file):
    """Test the extraction of segment offsets from a single segment PSF library
    in the temp_data directory.
    """
    # Get test library file
    test_library_dir = os.path.dirname(test_library_file)

    seg_id = 12
    segment_file = get_library_file(
        INSTRUMENT, DETECTOR, FILTER, 'CLEAR', '', 0, test_library_dir,
        segment_id=seg_id
    )

    # Set segment = 0 and expect error:
    with pytest.raises(AssertionError) as e:
        library_list = [segment_file]
        get_segment_offset(0, DETECTOR, library_list)
        assert "Uh-oh. The segment ID of the library does not match the requested segment." in e, \
            'Failed to catch mismatch between segment IDs.'

    # Put segment 12th on a faked library list and extract the segment offset
    library_list = [''] * 18
    library_list[11] = segment_file
    x_arcsec, y_arcsec = get_segment_offset(12, DETECTOR, library_list)
    assert np.all(np.isclose((x_arcsec, y_arcsec), (-7.412675266715986, -7.905276016530717), atol=1e-14)), \
        'Incorrect conversion of segment offsets'


@pytest.mark.skipif(ON_GITHUB,
                   reason="Cannot access mirage data in the central storage directory from Github Actions CI.")
def test_get_gridded_segment_psf_library_list_remote():
    """Test the loading of segment PSF libraries as a GriddedPSFModel on the
    central storage directory.
    """
    library_path = TEST_DATA_DIR

    libraries = get_gridded_segment_psf_library_list(INSTRUMENT, DETECTOR, FILTER,
                                         library_path)
    assert len(libraries) == 18, 'Did not find all 18 segment libraries'
    for i, lib_model in enumerate(libraries):
        assert isinstance(lib_model, photutils.psf.models.GriddedPSFModel), \
            'Segment PSF library not created correctly'
        assert lib_model.grid_xypos == [(1023.5, 1023.5)], \
            'Segment PSF library not created correctly'
        assert lib_model.oversampling == 1, \
            'Segment PSF library not created correctly'
        for k in ['segid', 'segname', 'xtilt', 'ytilt']:
            assert k in list(lib_model.meta.keys()), \
                'Segment PSF library not created correctly'
        assert lib_model.meta['segid'][0] == i + 1, \
            'Segment PSF library not created correctly'
        assert lib_model.data.shape == (1, 1024, 1024), \
            'Segment PSF library not created correctly'


def test_get_gridded_segment_psf_library_list_local():
    """Test the loading of a single segment PSF library as a GriddedPSFModel in
    the test_data directory.
    """
    # Skip if test_generate_segment_psfs failed
    if not SUCCESS_GENERATE_SEGMENT_PSF:
        pytest.skip("Cannot load library files generated in test_generate_segment_psfs.")

    library_path = TEMP_TEST_DIRECTORY

    libraries = get_gridded_segment_psf_library_list(INSTRUMENT, DETECTOR, FILTER,
                                         library_path)
    assert len(libraries) == 18, 'Did not find all 18 segment libraries'
    for i, lib_model in enumerate(libraries):
        assert isinstance(lib_model, photutils.psf.models.GriddedPSFModel), \
            'Segment PSF library not created correctly'
        assert lib_model.grid_xypos == [(1023.0, 1023.0)], \
            'Segment PSF library not created correctly'
        assert lib_model.oversampling == 1, \
            'Segment PSF library not created correctly'
        for k in ['segid', 'segname', 'xtilt', 'ytilt']:
            assert k in list(lib_model.meta.keys()), \
                'Segment PSF library not created correctly'
        assert lib_model.meta['segid'][0] == i + 1, \
            'Segment PSF library not created correctly'
        assert lib_model.data.shape == (1, 101, 101), \
            'Segment PSF library not created correctly'


def test_to_gridded_psfmodel(test_library_file):
    """Test that the example library file can be correctly loaded as a
    GriddedPSFModel using the webbpsf.utils.to_griddedpsfmodel function.

    Note that this is more a test of webbpsf than of MIRaGe, but it is
    required for MIRaGe to work!
    """
    with fits.open(test_library_file) as hdulist:
        lib_model = to_griddedpsfmodel(hdulist)

    assert isinstance(lib_model, photutils.psf.models.GriddedPSFModel), \
        'Segment PSF library not created correctly'
    assert lib_model.grid_xypos == [(1023.5, 1023.5)], \
        'Segment PSF library not created correctly'
    assert lib_model.oversampling == 1, \
        'Segment PSF library not created correctly'
    for k in ['segid', 'segname', 'xtilt', 'ytilt']:
        assert k in list(lib_model.meta.keys()), \
            'Segment PSF library not created correctly'
    assert lib_model.meta['segid'][0] == 12, \
        'Segment PSF library not created correctly'
    assert lib_model.data.shape == (1, 1024, 1024), \
        'Segment PSF library not created correctly'
