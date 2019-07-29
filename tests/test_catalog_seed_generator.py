#! /usr/bin/env python

"""Tests for ``catalog_seed_generator.py`` module

Authors
-------

    Bryan Hilbert

Use
---

    These tests can be run via the command line:

    ::

        pytest -s test_catalog_seed_generator.py
"""
from astropy.table import Table
import numpy as np
import os
import webbpsf

from mirage.seed_image import catalog_seed_image

# Determine if tests are being run on Travis
ON_TRAVIS = 'travis' in os.path.expanduser('~')
if ON_TRAVIS:
    env_path = os.path.expandvars('$CONDA_PREFIX')
    os.environ["WEBBPSF_PATH"] = os.path.join(env_path, 'share/webbpsf-data')


def create_dummy_psf_grid():
    """Use webbpsf to create a griddedPSFModel object"""
    nrc = webbpsf.NIRCam()
    nrc.filter = 'F200W'
    nrc.detector = 'NRCB1'
    grid = nrc.psf_grid(num_psfs=1, all_detectors=False, oversample=1, fov_pixels=301,
                        add_distortion=False, save=False)
    return grid


def test_overlap_coordinates_full_frame():
    """Test the function that calculates the coordinates for the
    overlap between two stamp images
    """
    # Create instance of Catalog_seed
    seed = catalog_seed_image.Catalog_seed(offline=True)

    # Create a point sources table
    tab = Table()
    tab['index'] = np.arange(10).astype(np.int) + 1
    tab['pixelx'] = [1000, 50, 1000, 1900, -0.5, 1033.587, 622.63, 1523.75, 1013.413, 1124.371]
    tab['pixely'] = [1000, 400, 25, 2000, -0.5, 1023.425, 1024.223, 1013.25, 223.575, 1822.72]
    tab['RA_degrees'] = [12.008729, 11.999914, 12.003422, 11.995729, 12.000036, 11.99918, 12.008729,
                         11.999914, 12.003422, 11.995729]
    tab['Dec_degrees'] = [-0.00885006, 1.7231344e-09, -1.5512744e-05, -4.9949034e-05, -0.006866415,
                          0.0068373946, -0.00885006, 1.7231344e-09, -1.5512744e-05, -4.9949034e-05]
    tab['magnitude'] = np.repeat(14, 10)
    tab['countrate_e/s'] = np.repeat(161630.72, 10)
    tab['counts_per_frame_e'] = np.repeat(1735391.9, 10)

    expected_results = []
    expected_results.append((850, 1151, 850, 1151, 0, 301, 0, 301))
    expected_results.append((0, 201, 250, 551, 100, 301, 0, 301))
    expected_results.append((850, 1151, 0, 176, 0, 301, 125, 301))
    expected_results.append((1750, 2048, 1850, 2048, 0, 298, 0, 198))

    stamp_dims = (301, 301)
    stamp_x = 150
    stamp_y = 150

    # Check the basic function
    aperture_dims = (2048, 2048)
    for index in range(4):
        results = seed.cropped_coords(tab[index]['pixelx'], tab[index]['pixely'], aperture_dims, stamp_x,
                                      stamp_y, stamp_dims, ignore_detector=False)
    assert results == expected_results[index]

    # Check the function that wraps around the most basic function
    seed.ffsize = 2048
    seed.subarray_bounds = [0, 0, 2047, 2047]
    for index in range(4):
        results = seed.create_psf_stamp_coords(tab[index]['pixelx'], tab[index]['pixely'], stamp_dims,
                                               stamp_x, stamp_y, coord_sys='full_frame',
                                               ignore_detector=False)
        ijkl = []
        for ele in results[4:]:
            ijkl.extend(list(ele))
        ijkl = tuple(ijkl)
        assert ijkl == expected_results[index]

    # Create dummy PSF library
    seed.psf_library = create_dummy_psf_grid()
    seed.psf_library_x_dim = seed.psf_library.data.shape[2] / seed.psf_library.oversampling
    seed.psf_library_y_dim = seed.psf_library.data.shape[1] / seed.psf_library.oversampling

    # Check the wrapped function, but call using 'aperture' rather than
    # 'full_frame'
    expected_k1l1 = [(850, 1151, 850, 1151, 0, 301, 0, 301),
                     (0, 201, 250, 551, 0, 201, 0, 301),
                     (850, 1151, 0, 176, 0, 301, 0, 176),
                     (1750, 2048, 1850, 2048, 0, 298, 0, 198),
                     (0, 150, 0, 150, 0, 150, 0, 150)]

    seed.coord_adjust['xoffset'] = 0  # Already defined, but let's be explicit
    seed.coord_adjust['yoffset'] = 0  # Already defined, but let's be explicit
    seed.output_dims = [2048, 2048]
    seed.psf_library_core_x_dim = 301
    seed.psf_library_core_y_dim = 301
    seed.params = {'simSignals': {}}
    seed.add_psf_wings = False
    seed.psf_library_oversamp = 1

    for index in range(5):
        psf, min_x, min_y, add_wings = seed.create_psf_stamp(tab[index]['pixelx'], tab[index]['pixely'],
                                                             stamp_dims[1], stamp_dims[0],
                                                             ignore_detector=False)
        stamp_x_loc = 301 // 2 - min_x
        stamp_y_loc = 301 // 2 - min_y
        updated_psf_dimensions = psf.shape

        xap, yap, xpts, ypts, (i1, i2), (j1, j2), (k1, k2), \
            (l1, l2) = seed.create_psf_stamp_coords(tab[index]['pixelx'], tab[index]['pixely'],
                                                    updated_psf_dimensions, stamp_x_loc, stamp_y_loc,
                                                    coord_sys='aperture')
        assert (i1, i2, j1, j2, k1, k2, l1, l2) == expected_k1l1[index]


def test_overlap_coordinates_subarray():
    """Test the function that calculates the coordinates for the
    overlap between two stamp images
    """
    seed = catalog_seed_image.Catalog_seed(offline=True)

    # Create a point sources table
    tab = Table()
    tab['index'] = np.arange(5).astype(np.int) + 1
    tab['pixelx'] = [200, 50, 200, 350, -0.5]
    tab['pixely'] = [200, 200, 25, 350, -0.5]
    tab['RA_degrees'] = [12.008729, 11.999914, 12.003422, 11.995729, 11.995729]
    tab['Dec_degrees'] = [-0.00885006, 1.7231344e-09, -1.5512744e-05, -4.9949034e-05, -4.9949034e-05]
    tab['magnitude'] = np.repeat(14, 5)
    tab['countrate_e/s'] = np.repeat(161630.72, 5)
    tab['counts_per_frame_e'] = np.repeat(1735391.9, 5)

    subarray_expected_results = []
    subarray_expected_results.append((50, 351, 50, 351, 0, 301, 0, 301))
    subarray_expected_results.append((0, 201, 50, 351, 100, 301, 0, 301))
    subarray_expected_results.append((50, 351, 0, 176, 0, 301, 125, 301))
    subarray_expected_results.append((200, 400, 200, 400, 0, 200, 0, 200))
    subarray_expected_results.append((0, 150, 0, 150, 151, 301, 151, 301))

    # Test the basic function using a subarray
    stamp_dims = (301, 301)
    stamp_x = 150
    stamp_y = 150
    aperture_dims = (400, 400)  # Assume 400x400 pixel subarray
    for index in range(5):
        results = seed.cropped_coords(tab[index]['pixelx'], tab[index]['pixely'], aperture_dims, stamp_x,
                                      stamp_y, stamp_dims, ignore_detector=False)
    assert results == subarray_expected_results[index]

    # Test the wrapper function
    seed.coord_adjust['xoffset'] = 0  # Already defined, but let's be explicit
    seed.coord_adjust['yoffset'] = 0  # Already defined, but let's be explicit
    seed.output_dims = [400, 400]  # Assume 400x400 pixel subarray
    seed.subarray_bounds = [1648, 1648, 2047, 2047]
    seed.ffsize = 2048
    seed.psf_library_core_x_dim = 301
    seed.psf_library_core_y_dim = 301
    seed.params = {'simSignals': {}}
    seed.add_psf_wings = False
    seed.psf_library_oversamp = 1
    for index in range(5):
        results = seed.create_psf_stamp_coords(tab[index]['pixelx'], tab[index]['pixely'], stamp_dims,
                                               stamp_x, stamp_y, coord_sys='aperture', ignore_detector=False)
        ijkl = []
        for ele in results[4:]:
            ijkl.extend(list(ele))
        ijkl = tuple(ijkl)
    assert ijkl == subarray_expected_results[index]

    # Dummy PSF library
    seed.psf_library = create_dummy_psf_grid()
    seed.psf_library_x_dim = seed.psf_library.data.shape[2] / seed.psf_library.oversampling
    seed.psf_library_y_dim = seed.psf_library.data.shape[1] / seed.psf_library.oversampling

    # Test the wrapper function using a modified PSF stamp image shape
    expected_k1l1 = [(50, 351, 50, 351, 0, 301, 0, 301),
                     (0, 201, 50, 351, 100, 301, 0, 301),
                     (50, 351, 0, 176, 0, 301, 125, 301),
                     (200, 400, 200, 400, 0, 200, 0, 200),
                     (0, 150, 0, 150, 151, 301, 151, 301)]
    for index in range(5):
        psf, min_x, min_y, add_wings = seed.create_psf_stamp(tab[index]['pixelx'], tab[index]['pixely'],
                                                             stamp_dims[1], stamp_dims[0],
                                                             ignore_detector=False)
        stamp_x_loc = 301 // 2 - min_x
        stamp_y_loc = 301 // 2 - min_y
        updated_psf_dimensions = psf.shape

        xap, yap, xpts, ypts, (i1, i2), (j1, j2), (k1, k2), \
            (l1, l2) = seed.create_psf_stamp_coords(tab[index]['pixelx'], tab[index]['pixely'],
                                                    updated_psf_dimensions, stamp_x_loc, stamp_y_loc,
                                                    coord_sys='aperture')
        assert (i1, i2, j1, j2, k1, k2, l1, l2) == expected_k1l1[index]
