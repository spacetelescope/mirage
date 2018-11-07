import os

import astropy.convolution
from astropy.io import fits
import numpy as np
import pytest
import webbpsf

from mirage.psf.psf_library import CreatePSFLibrary


def test_all_filters_and_detectors():
    """Check that setting filters and detectors to all works"""

    # Case 1: Setting filters="all" and "detectors="all"
    inst1 = CreatePSFLibrary(instrument="FGS", filters="all", detectors="all",
                             num_psfs=1, save=False)
    grid1 = inst1.create_files()
    inst2 = CreatePSFLibrary(instrument="FGS", filters="FGS", detectors=["FGS1", "FGS2"],
                             num_psfs=1, save=False)
    grid2 = inst2.create_files()

    # Check outputs are the same
    assert np.array_equal(grid1[0][0].data, grid2[0][0].data)

    # Case 2: NIRCam should be able to pull all the appropriate detectors based on the filter type (SW vs LW)
    longfilt = "F250M"
    shortfilt = "F140M"
    inst3 = CreatePSFLibrary(instrument="NIRCam", filters=[shortfilt, longfilt],
                             detectors="all", num_psfs=1, add_distortion=False,
                             fov_pixels=1, oversample=2, save=False)
    grid1, grid2 = inst3.create_files()

    # Check that only and all the SW detectors are in the first file
    det_list = []
    for i in range(len(grid1[0].data)):
        det_list.append(grid1[0].header["DETNAME{}".format(i)])
    assert len(grid1[0].data) == len(CreatePSFLibrary.nrca_short_detectors)
    assert set(det_list) == set(CreatePSFLibrary.nrca_short_detectors)

    # Check that only and all the LW detectors are in the second file
    det_list = []
    for i in range(len(grid2[0].data)):
        det_list.append(grid2[0].header["DETNAME{}".format(i)])
    assert len(grid2[0].data) == len(CreatePSFLibrary.nrca_long_detectors)
    assert set(det_list) == set(CreatePSFLibrary.nrca_long_detectors)


def test_compare_to_calc_psf():
    """Check that the output grid has the expected PSFs in the right grid locations by comparing to calc_psf"""

    oversample = 2
    fov_pixels = 10

    # Create a PSF grid
    inst = CreatePSFLibrary(instrument="FGS", filters="FGS", detectors="FGS1", num_psfs=4,
                            oversample=oversample, fov_pixels=fov_pixels, save=False)
    grid = inst.create_files()

    # Pull one of the PSFs out of the grid
    psfnum = 1
    loc = grid[0][0].header["DET_YX{}".format(psfnum)]  # (y,x) location
    pos = grid[0][0].header["DET_JI{}".format(psfnum)]  # (j,i) position
    locy = int(loc.split()[0][1:-1])
    locx = int(loc.split()[1][:-1])
    posj = int(pos.split()[0][1:-1])
    posi = int(pos.split()[1][:-1])
    gridpsf = grid[0][0].data[0, posj, posi, :, :]

    # Using grid header data, create the expected same PSF via calc_psf + convolution
    fgs = webbpsf.FGS()
    fgs.detector = "FGS1"
    fgs.filter = "FGS"
    fgs.detector_position = (locx, locy)
    calcpsf = fgs.calc_psf(oversample=oversample, fov_pixels=fov_pixels)["OVERDIST"].data
    kernel = astropy.convolution.Box2DKernel(width=oversample)
    convpsf = astropy.convolution.convolve(calcpsf, kernel)

    # Compare to make sure they are in fact the same PSF
    assert np.array_equal(gridpsf, convpsf)


def test_nircam_errors():
    """Check that there are errors for incorrect value setting - particularly with NIRCam"""

    longfilt = "F250M"
    shortfilt = "F140M"
    longdet = "NRCB5"
    shortdet = "NRCA3"

    # Shouldn't error - applying SW to SW and LW to LW
    inst1 = CreatePSFLibrary(instrument="NIRCam", filters=longfilt, detectors=longdet,
                             add_distortion=False, num_psfs=1, fov_pixels=1, save=False)  # no error
    inst1.create_files()
    inst2 = CreatePSFLibrary(instrument="NIRCam", filters=shortfilt, detectors=shortdet,
                             add_distortion=False, num_psfs=1, fov_pixels=1, save=False)  # no error
    inst2.create_files()

    # Should error - Bad filter/detector combination (LW filt to SW det)
    with pytest.raises(ValueError) as excinfo:
        inst3 = CreatePSFLibrary(instrument="NIRCam", filters=longfilt, detectors=shortdet,
                                 add_distortion=False, num_psfs=1, fov_pixels=1, save=False)  # error
        inst3.create_files()
    assert "ValueError" in str(excinfo)

    # Should error - Bad filter/detector combination (SW filt to LW det)
    with pytest.raises(ValueError) as excinfo:
        inst4 = CreatePSFLibrary(instrument="NIRCam", filters=shortfilt, detectors=longdet,
                                 add_distortion=False, num_psfs=1, fov_pixels=1, save=False)  # error
        inst4.create_files()
    assert "ValueError" in str(excinfo)

    # Should error - Bad num_psfs entry (must be a square number)
    with pytest.raises(ValueError) as excinfo:
        inst5 = CreatePSFLibrary(instrument="NIRCam", filters=longfilt, detectors=longdet,
                                 add_distortion=False, num_psfs=3, fov_pixels=1, save=False)  # error
        inst5.create_files()
    assert "ValueError" in str(excinfo)


def test_one_psf():
    """Check that setting num_psfs = 1 produces the PSF in the expected location"""

    oversample = 2
    fov_pixels = 101

    # Create 2 cases with different locations: the default center and a set location
    inst1 = CreatePSFLibrary(instrument="NIRISS", filters="F090W", detectors="NIS",
                             num_psfs=1, add_distortion=True, oversample=oversample,
                             fov_pixels=fov_pixels, save=False)
    grid1 = inst1.create_files()
    inst2 = CreatePSFLibrary(instrument="NIRISS", filters="F090W", detectors="NIS",
                             num_psfs=1, add_distortion=True, oversample=oversample,
                             fov_pixels=fov_pixels, psf_location=(0, 10), save=False)
    grid2 = inst2.create_files()

    assert grid1[0][0].header["DET_YX0"] == "(1024, 1024)"  # the default is the center of the NIS aperture
    assert grid2[0][0].header["DET_YX0"] == "(0, 10)"  # (y,x)

    # Compare to the WebbPSF calc_psf output to make sure it's placing the PSF in the right location
    nis = webbpsf.NIRISS()
    nis.filter = "F090W"
    nis.detector_position = (10, 0)  # (x,y)
    calc = nis.calc_psf(add_distortion=True, oversample=oversample, fov_pixels=fov_pixels)
    kernel = astropy.convolution.Box2DKernel(width=oversample)
    convpsf = astropy.convolution.convolve(calc["OVERDIST"].data, kernel)

    assert np.array_equal(convpsf, grid2[0][0].data[0, 0, 0, :, :])


def test_saving(tmpdir):
    """Test saving files works properly"""

    # Create a temp directory to place file in
    path = str(tmpdir)

    # Test using default calc_psf values
    inst = CreatePSFLibrary(instrument="FGS", filters="FGS", detectors="FGS1", opd_type="predicted",
                            opd_number=3, num_psfs=4, oversample=2, fov_pixels=10,
                            save=True, fileloc=path, filename=None)
    grid = inst.create_files()

    # Check that the saved file matches the returned file (and thus that the save worked through properly)
    with fits.open(os.path.join(path, "fgs_fgs_fovp10_samp2_npsf4.fits")) as infile:
        assert infile[0].header == grid[0][0].header
        assert np.array_equal(infile[0].data, grid[0][0].data)

    # Remove temp directory
    tmpdir.remove()


def test_setting_inputs():
    """Test that no errors are raised when filters and detectors can be set different ways"""

    # Pass strings for filter and detector
    inst1 = CreatePSFLibrary(instrument="NIRCam", filters="F090W", detectors="NRCA4",
                             num_psfs=1, fov_pixels=1, oversample=2, save=False)
    grid1 = inst1.create_files()

    # Pass a list of filters
    inst2 = CreatePSFLibrary(instrument="NIRCam", filters=["F090W", "F140M"], detectors="NRCA4",
                             num_psfs=1, fov_pixels=1, oversample=2, save=False)
    grid2 = inst2.create_files()

    # Pass lists for filter and detector
    inst3 = CreatePSFLibrary(instrument="NIRCam", filters=["F090W"], detectors=["NRCA4"],
                             num_psfs=1, fov_pixels=1, oversample=2, save=False)
    grid3 = inst3.create_files()

    # Pass the name "longwave" to both filter and detector
    inst4 = CreatePSFLibrary(instrument="NIRCam", filters="longwave", detectors="longwave",
                             num_psfs=1, fov_pixels=1, oversample=2, save=False)
    grid4 = inst4.create_files()

    # # Pass the name "shortwave" to both filter and detector - This case takes 6min alone to run - excluding for time
    # inst5 = CreatePSFLibrary(instrument="NIRCam", filters="shortwave", detectors="shortwave",
    #                          num_psfs=1, fov_pixels=1, oversample=2, save=False)
    # grid5 = inst5.create_files()

    # Length is the number of filters, Shape of those objects follows (det, j, i, y, x)
    assert len(grid1) == 1
    assert grid1[0][0].data.shape == (1, 1, 1, 2, 2)
    assert len(grid2) == 2
    assert grid2[0][0].data.shape == (1, 1, 1, 2, 2)
    assert len(grid3) == 1
    assert grid3[0][0].data.shape == (1, 1, 1, 2, 2)
    assert len(grid4) == 16
    assert grid4[0][0].data.shape == (2, 1, 1, 2, 2)
    # assert len(grid5) == 13
    # assert grid5[0][0].data.shape == (8, 1, 1, 2, 2)

    # Check that passing strings and lists produces the same result
    assert np.array_equal(grid1[0][0].data, grid3[0][0].data)
