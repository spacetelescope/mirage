'''Define unit tests for calculation of single frame exposure times

Authors
-------
    - Bryan Hilbert

Use
---
    Ensure you have pytest installed. Then, simply run pytest in any
    parent directory of mirage/tests/:
    >>> pytest
'''

import numpy as np

from mirage.utils.utils import calc_frame_time


def test_fgs_cal_frametime():
    """These tests apply to the FGS CAL apertures only. ACQ1, ACQ2, ID,
    TRK, and FG are not yet supported
    """
    fgs_full = calc_frame_time('fgs', 'FGS_', 2048, 2048, 4)
    assert np.isclose(fgs_full, 10.73677, rtol=0., atol=1e-5)

    fgs_128 = calc_frame_time('fgs', 'FGS_', 128, 128, 1)
    assert np.isclose(fgs_128, 0.1820, rtol=0., atol=1e-5)

    fgs_32 = calc_frame_time('fgs', 'FGS_', 32, 32, 1)
    assert np.isclose(fgs_32, 0.01254, rtol=0., atol=1e-5)

    fgs_8 = calc_frame_time('fgs', 'FGS_', 8, 8, 1)
    assert np.isclose(fgs_8, 0.00126, rtol=0., atol=1e-5)


def test_nircam_frametime():
    """Test NIRCam exposure times
    """
    nrc_full = calc_frame_time('nircam', 'NRCA1_FULL', 2048, 2048, 4)
    assert np.isclose(nrc_full, 10.73677, rtol=0., atol=1e-5)

    nrc_640 = calc_frame_time('nircam', 'NRCA1_SUB640', 640, 640, 1)
    assert np.isclose(nrc_640, 4.18584, rtol=0., atol=1e-5)

    nrc_320 = calc_frame_time('nircam', 'NRCA1_SUB320', 320, 320, 1)
    assert np.isclose(nrc_320, 1.06904, rtol=0., atol=1e-5)

    nrc_160 = calc_frame_time('nircam', 'NRCA1_SUB160', 160, 160, 1)
    assert np.isclose(nrc_160, 0.27864, rtol=0., atol=1e-5)

    nrc_64 = calc_frame_time('nircam', 'NRCB4_SUB64P', 64, 64, 1)
    assert np.isclose(nrc_64, 0.05016, rtol=0., atol=1e-5)

    nrc_32 = calc_frame_time('nircam', 'NRC_SUB32TATS', 32, 32, 1)
    assert np.isclose(nrc_32, 0.01496, rtol=0., atol=1e-5)

    nrc_subgrism256_1 = calc_frame_time('nircam', 'NRC_SUBGRISM256', 2048, 256, 1)
    print(nrc_subgrism256_1, 5.31480)
    #assert np.isclose(nrc_subgrism256_1, 5.29420, rtol=0., atol=1e-5)

    nrc_subgrism256_4 = calc_frame_time('nircam', 'NRC_SUBGRISM256', 2048, 256, 4)
    print(nrc_subgrism256_4, 1.34669)
    #assert np.isclose(nrc_subgrism256_4, 1.34669, rtol=0., atol=1e-5)

    nrc_subgrism128_1 = calc_frame_time('nircam', 'NRC_SUBGRISM128', 2048, 128, 1)
    print(nrc_subgrism128_1, 2.67800)
    #assert np.isclose(nrc_subgrism128_1, 2.6574, rtol=0., atol=1e-5)

    nrc_subgrism128_4 = calc_frame_time('nircam', 'NRC_SUBGRISM128', 2048, 128, 4)
    assert np.isclose(nrc_subgrism128_4, 0.67597, rtol=0., atol=1e-5)

    nrc_subgrism64_1 = calc_frame_time('nircam', 'NRC_SUBGRISM64', 2048, 64, 1)
    assert np.isclose(nrc_subgrism64_1, 1.35960, rtol=0., atol=1e-5)

    nrc_subgrism64_4 = calc_frame_time('nircam', 'NRC_SUBGRISM64', 2048, 64, 4)
    assert np.isclose(nrc_subgrism64_4, 0.34061, rtol=0., atol=1e-5)


def test_niriss_frametime():
    """Test NIRISS exposure times
    """
    nis_full = calc_frame_time('niriss', 'NIS_CEN', 2048, 2048, 4)
    assert np.isclose(nis_full, 10.73677, rtol=0., atol=1e-5)

    nis_80 = calc_frame_time('niriss', 'NIS_SUB80', 80, 80, 1)
    assert np.isclose(nis_80, 0.07544, rtol=0., atol=1e-5)

    nis_64 = calc_frame_time('niriss', 'NIS_SUBTAAMI', 64, 64, 1)
    assert np.isclose(nis_64, 0.05016, rtol=0., atol=1e-5)

    nis_wfss64r = calc_frame_time('niriss', 'NIS_WFSS64R', 64, 2048, 4)
    assert np.isclose(nis_wfss64r, 0.34061, rtol=0., atol=1e-5)

    nis_wfss64c = calc_frame_time('niriss', 'NIS_WFSS64C', 2048, 64, 1)
    assert np.isclose(nis_wfss64c, 1.55800, rtol=0., atol=1e-5)

    nis_wfss128r = calc_frame_time('niriss', 'NIS_WFSS128R', 128, 2048, 4)
    assert np.isclose(nis_wfss128r, 0.67597, rtol=0., atol=1e-5)

    nis_wfss128c = calc_frame_time('niriss', 'NIS_WFSS128C', 2048, 128, 1)
    assert np.isclose(nis_wfss128c, 2.87000, rtol=0., atol=1e-5)
