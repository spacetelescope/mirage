#! /usr/bin/env python

"""Test the locations of sources in various apertures using RADec_To_XY
Note that the truth values in this case were calculated using the
distortion reference files, which are not part of the repo and not
available to Github Actions CI
Values were calculated using run_all_apertures.py in
the test_mirage_subarray_source_location subdirectory
"""

from astropy.table import Table
import numpy as np
import pysiaf
import pytest

from mirage.seed_image import catalog_seed_image
from mirage.utils import siaf_interface


def find_distortion_file(inst, aper):
    """Identify the appropriate distortion reference file for a given
    instrument/aperture combination
    Parameters
    ----------
    inst : str
        Instrument name
    aper : str
        Aperture name (e.g NRCA1_FULL)
    """



def test_locations():
    """Test RA, Dec to x,y translation. This function does the
    translation using pysiaf. This is the truth case we will
    compare against the version translated via distortion
    reference file, which is in ``reference_file_values``. We have to
    test like this because Github Actions CI does not have access to the
    distortion reference file itself.
    """
    # Apertures to skip for this testing. These apertures don't use the
    # reference files in the intended way.
    apertures_to_skip = ['NRCA5_GRISMC_WFSS', 'NRCA5_GRISMR_WFSS', 'NRCB5_GRISMC_WFSS',
                         'NRCB5_GRISMR_WFSS', 'NIS_SUBSTRIP96', 'NIS_SUBSTRIP256', 'J-FRAME']

    # RA, Dec of pointing
    c = catalog_seed_image.Catalog_seed()

    pointing_ra = 12.0
    pointing_dec = 12.0
    rotation = 0.0

    ra_list = []
    dec_list = []
    delta_pix = np.array([0, 10, 20, 40, 75, 180, 500, 1000])
    delta_deg = delta_pix * 0.031 / 3600.
    for delta in delta_deg:
        ra_entry = pointing_ra + delta
        dec_entry = pointing_dec + delta
        ra_list.append(ra_entry)
        dec_list.append(dec_entry)

    xy_vals = {}
    instrument_list = ['niriss', 'fgs']

    for instrument in instrument_list:
        siaf = siaf_interface.get_instance(instrument)
        for aperture in siaf.apernames:
            if 'WEDGE' in aperture or 'MASK' in aperture or aperture in apertures_to_skip:
                continue

            c.local_roll, c.attitude_matrix, c.ffsize, \
                c.subarray_bounds = siaf_interface.get_siaf_information(siaf,
                                                                        aperture,
                                                                        pointing_ra, pointing_dec,
                                                                        rotation)
            c.coord_transform = None
            c.siaf = siaf[aperture]
            c.use_intermediate_aperture = False

            xvals_from_siaf = []
            yvals_from_siaf = []
            for ra, dec in zip(ra_list, dec_list):
                x, y = c.RADecToXY_astrometric(ra, dec)
                ra_check, dec_check, ra_str, dec_str = c.XYToRADec(x, y)
                assert np.isclose(ra, ra_check, rtol=0, atol=8e-6)
                assert np.isclose(dec, dec_check, rtol=0, atol=8e-6)
                xvals_from_siaf.append(x)
                yvals_from_siaf.append(y)


@pytest.mark.parametrize("pav3,galaxy_pos_angs,expected_galaxy_pos_angs", [
    (0., [0., 30., -80.], [90.1101252253335, 120.60294962655526, 8.463848460849974]),
    (30., [0., 30., -80.], [59.46543695693042, 90.1798908477186, -21.71832856291518])
    ])
def test_gal_rotations(pav3, galaxy_pos_angs, expected_galaxy_pos_angs):
    """Test the rotations used to place galaxy and extended sources in the correct
    orientation
    """
    aperture = 'NRCB2_FULL'
    pointing_ra = 12.008818
    pointing_dec = 45.008818

    c = catalog_seed_image.Catalog_seed()
    c.use_intermediate_aperture = False
    inst_siaf = siaf_interface.get_instance('nircam')
    c.siaf = inst_siaf[aperture]
    c.coord_transform = None
    c.local_roll, c.attitude_matrix, c.ffsize, \
                c.subarray_bounds = siaf_interface.get_siaf_information(inst_siaf,
                                                                        aperture,
                                                                        pointing_ra, pointing_dec,
                                                                        pav3)
    pa_table = Table()
    pa_table['RA_degrees'] = [pointing_ra] * len(galaxy_pos_angs)
    pa_table['Dec_degrees'] = [pointing_dec] * len(galaxy_pos_angs)
    pa_table['pos_angle'] = galaxy_pos_angs
    pa_table['pixelx'] = [1024] * len(galaxy_pos_angs)
    pa_table['pixely'] = [1024] * len(galaxy_pos_angs)
    gal_x_pas = [c.calc_x_position_angle(row) for row in pa_table]

    #gal_x_pas = [c.calc_x_position_angle(pa) for pa in galaxy_pos_angs]

    print(gal_x_pas)
    print(expected_galaxy_pos_angs)


    assert np.all(np.isclose(gal_x_pas, expected_galaxy_pos_angs, atol=0.1))


@pytest.mark.parametrize("pav3,ext_pos_angs,expected_ext_pos_angs", [
    (0., [0., 30., 330., -30.], [359.9664210290788, 389.9664210290788, 689.9664210290788, 329.9664210290788]),
    (30., [0., 30., 330., -30.], [29.898089979824874, 59.898089979824874, 359.8980899798249, -0.10191002017512574])
    ])
def test_ext_rotations(pav3, ext_pos_angs, expected_ext_pos_angs):
    """Test the rotations used to place galaxy and extended sources in the correct
    orientation
    """
    aperture = 'NRCB2_FULL'
    pointing_ra = 12.008818
    pointing_dec = 45.008818

    c = catalog_seed_image.Catalog_seed()
    c.use_intermediate_aperture = False
    inst_siaf = siaf_interface.get_instance('nircam')
    c.siaf = inst_siaf[aperture]
    c.local_roll, c.attitude_matrix, c.ffsize, \
                c.subarray_bounds = siaf_interface.get_siaf_information(inst_siaf,
                                                                        aperture,
                                                                        pointing_ra, pointing_dec,
                                                                        pav3)

    ext_x_pas = [c.calc_x_position_angle_extended(pa) for pa in ext_pos_angs]

    print(ext_x_pas)
    print(expected_ext_pos_angs)

    assert np.all(np.isclose(ext_x_pas, expected_ext_pos_angs, atol=0.1))




