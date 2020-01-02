'''Define unit tests that cover the addition of a TSO source to a seed image

Authors
-------
    - Bryan Hilbert

Use
---
    Ensure you have pytest installed. Then, simply run pytest in any
    parent directory of mirage/tests/:
    >>> pytest
'''
import astropy.units as u
import numpy as np
from synphot import units

from mirage.seed_image import segmentation_map as segmap
from mirage.seed_image import tso
from mirage.utils.constants import FLAMBDA_CGS_UNITS

def test_add_tso_sources():
    """Test the addition of the TSO source to the seed image
    """
    background_seed_image = np.zeros((10, 10))
    background_seed_image[1:4, 1:4] = 100

    background_seg_map = np.zeros((10, 10)).astype(np.int)
    background_seg_map[1:4, 1:4] = 2

    tso_seed = np.zeros((10, 10))
    tso_seed[3:6, 3:6] = 1000
    tso_seeds = [tso_seed]

    tso_seg = segmap.SegMap()
    tso_seg.ydim, tso_seg.xdim = (10, 10)
    tso_seg.initialize_map()
    tso_seg.segmap[3:6, 3:6] = 99999
    tso_segs = [tso_seg]

    light_curve = {}
    light_curve['times'] = np.arange(0, 10) * u.second
    light_curve['fluxes'] = np.repeat(1., len(light_curve['times'])) * FLAMBDA_CGS_UNITS
    lcs = [light_curve]

    frame_time = 2.
    total_frames = 5
    exposure_frames = 5
    frames_per_int = 5
    num_ints = 1
    resets = 1

    seed, seg = tso.add_tso_sources(background_seed_image, background_seg_map, tso_seeds, tso_segs, lcs, frame_time,
                                    total_frames, exposure_frames, frames_per_int, num_ints, resets,
                                    starting_time=0, starting_frame=0, samples_per_frametime=5)

    # Make sure the seed image is 4D and the segmentation map is 2D
    assert len(seed.shape) == 4
    assert len(seg.shape) == 2

    # Check the signal levels in the seed image frames
    assert np.all(seed[0, :, 2, 2] == np.array([200., 400., 600., 800., 1000.]))
    assert np.all(seed[0, :, 5, 5] == np.array([2000., 4000., 6000., 8000., 10000.]))

    # Check that both objects are in the segmentation map
    assert np.unique(seg[1:3, 1:3]) == [2]
    assert np.unique(seg[3:6, 3:6]) == [99999]


def test_check_lightcurve_time():
    """Check that the lightcurve time is compared to the total
    exposure time correctly
    """
    exposure_time = 1000.
    frame_time = 1.

    # Case where light curve is too short
    light_curve = {}
    light_curve['times'] = np.arange(0, 500) * u.second
    light_curve['fluxes'] = np.repeat(1., len(light_curve['times'])) * FLAMBDA_CGS_UNITS
    new_lc = tso.check_lightcurve_time(light_curve, exposure_time, frame_time)
    assert np.max(new_lc['times'].data) >= exposure_time

    # Case where lightcurve start time is < 0
    light_curve = {}
    light_curve['times'] = np.arange(-10, 500) * u.second
    light_curve['fluxes'] = np.repeat(1., len(light_curve['times'])) * FLAMBDA_CGS_UNITS
    new_lc = tso.check_lightcurve_time(light_curve, exposure_time, frame_time)
    assert np.min(new_lc['times'].data) == 0.

    # Case where lightcurve start time is > 0
    light_curve = {}
    light_curve['times'] = np.arange(30, 500) * u.second
    light_curve['fluxes'] = np.repeat(1., len(light_curve['times'])) * FLAMBDA_CGS_UNITS
    new_lc = tso.check_lightcurve_time(light_curve, exposure_time, frame_time)
    assert np.min(new_lc['times'].data) == 0.


def test_interpolate_lightcurve():
    """Test that a lightcurve is correctly interpolated to have the
    correct number of entries
    """
    samples_per_frame_time = 5
    frame_time = 1.

    # Create simple lightcurve
    light_curve = {}
    light_curve['times'] = np.arange(0, 10) * u.second
    light_curve['fluxes'] = np.repeat(1., len(light_curve['times'])) * FLAMBDA_CGS_UNITS

    # Interpolate
    interp = tso.interpolate_lightcurve(light_curve, samples_per_frame_time, frame_time)

    assert len(interp['times'].value) == (10 - 1) * frame_time * (samples_per_frame_time - 1)
    assert np.all(interp['times'].value[0:5] == [0, 0.25, 0.5, 0.75, 1.0])


def test_update_segmentation_map():
    """Test that the segmentation map is updated correctly when adding a
    new object
    """
    seg_map = np.zeros((10, 10))
    seg_map[3:6, 3:6] = 10

    new_obj = np.zeros((10, 10))
    new_obj[5:8, 5:8] = 999

    new_map = tso.update_segmentation_map(seg_map, new_obj)

    assert np.unique(new_map[3:5, 3:5]) == [10]
    assert np.unique(new_map[5:8, 5:8]) == [999]

