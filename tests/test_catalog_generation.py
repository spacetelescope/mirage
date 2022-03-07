#! /usr/bin/env python
"""Test the catalog generation functionality in Mirage

Authors
-------
    - Bryan Hilbert

Use
---
    >>> pytest -s test_catalog_generation.py


"""

import copy
import numpy as np
import os

from astropy.io import ascii
from astropy.table import Table
import pytest

from mirage.catalogs import catalog_generator
from mirage.catalogs import create_catalog, utils

TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), 'test_data/')

# Determine if tests are being run on Github actions CI
ON_GITHUB = '/home/runner' in os.path.expanduser('~')

if not ON_GITHUB:
    orig_mirage_data = os.environ['MIRAGE_DATA']
os.environ['MIRAGE_DATA'] = '/test/'


def test_ptsrc_catalog_creation():
    """Test the generation of a basic point source catalog
    """
    ra = np.zeros(10) + 80.0
    dec = np.zeros(10) - 69.8
    mags = np.zeros(10) + 19.

    ptsrc = catalog_generator.PointSourceCatalog(ra=ra, dec=dec)
    ptsrc.add_magnitude_column(mags, instrument='nircam', filter_name='f090w')
    ptsrc.add_magnitude_column(mags, instrument='nircam', filter_name='f200w')
    ptsrc.add_magnitude_column(mags, column_name='nircam_f090w_wlp8_magnitude')

    assert all(ptsrc.ra == ra)
    assert all(ptsrc.dec == dec)
    assert ptsrc.table.colnames == ['index', 'x_or_RA', 'y_or_Dec', 'nircam_f090w_clear_magnitude',
                                    'nircam_f200w_clear_magnitude', 'nircam_f090w_wlp8_magnitude']


def test_galaxy_catalog_creation():
    """Test the creation of a basic galaxy object catalog
    """
    indexes = np.arange(10, 20)
    ra = np.zeros(10) + 80.0
    dec = np.zeros(10) - 69.8
    mags = np.zeros(10) + 19.
    radius = np.zeros(10) + 0.5
    ellip = np.zeros(10) + 0.45
    posang = np.zeros(10) + 27.
    sersic = np.zeros(10) + 3.3

    output_file = os.path.join(TEST_DATA_DIR, 'catalog_generation/galaxy_test.cat')
    gal = catalog_generator.GalaxyCatalog(ra=ra, dec=dec, ellipticity=ellip, radius=radius,
                                          sersic_index=sersic, position_angle=posang, radius_units='arcsec',
                                          starting_index=10)
    gal.add_magnitude_column(mags, instrument='nircam', filter_name='f090w')
    gal.save(output_file)
    as_read_in = ascii.read(output_file)

    assert gal.table.colnames == ['index', 'x_or_RA', 'y_or_Dec', 'pos_angle', 'sersic_index',
                                  'ellipticity', 'radius', 'nircam_f090w_clear_magnitude']
    assert all(gal.radius == radius)
    assert all(gal.ellipticity == ellip)
    assert all(gal.position_angle == posang)
    assert all(gal.sersic_index == sersic)
    assert all(gal.table['index'].data == indexes)
    assert all(gal.table == as_read_in)
    os.remove(output_file)


def test_moving_point_source_catalog_generation():
    """Look specifically at the position and location columns
    in the case where an ephemeris file is provided.
    """
    mags = [15, 16, 17]
    ra_list = [10., 11., 12.]
    dec_list = [34., 35., 36.]

    ra_vel_list = [0.2, 0.3, 0.4]
    dec_vel_list = [1.0, 1.1, 1.2]

    ephemeris_list = ['e1.txt', 'e2.txt', 'e3.txt']

    mpt = catalog_generator.MovingPointSourceCatalog(ra=ra_list, dec=dec_list, ra_velocity=ra_vel_list, dec_velocity=dec_vel_list,
                                                 ephemeris_file=ephemeris_list)
    mpt.add_magnitude_column(mags, magnitude_system='abmag', instrument='nircam', filter_name='f444w')

    comparison = Table()
    comparison['index'] = [1, 2, 3]
    comparison['x_or_RA'] = [np.nan, np.nan, np.nan]
    comparison['y_or_Dec'] = [np.nan, np.nan, np.nan]
    comparison['x_or_RA_velocity'] = [np.nan, np.nan, np.nan]
    comparison['y_or_Dec_velocity'] = [np.nan, np.nan, np.nan]
    comparison['ephemeris_file'] = [ephemeris_list]
    comparison['nircam_f444w_clear_magnidute'] = mags

    nan_cols = ['x_or_RA', 'y_or_Dec', 'x_or_RA_velocity', 'y_or_Dec_velocity']
    for colname in nan_cols:
        print()
        assert np.all(np.isnan(comparison[colname]) == np.isnan(mpt.table[colname]))
    assert np.all(comparison['ephemeris_file'] == mpt.table['ephemeris_file'])

    mixed_ephemeris_list = ['e1.txt', 'None', 'ex3.txt']
    mpt = catalog_generator.MovingPointSourceCatalog(ra=ra_list, dec=dec_list, ra_velocity=ra_vel_list, dec_velocity=dec_vel_list,
                                                     ephemeris_file=mixed_ephemeris_list)
    mpt.add_magnitude_column(mags, magnitude_system='abmag', instrument='nircam', filter_name='f444w')

    comparison['x_or_RA'] = [np.nan, 11., np.nan]
    comparison['y_or_Dec'] = [np.nan, 35., np.nan]
    comparison['x_or_RA_velocity'] = [np.nan, 0.3, np.nan]
    comparison['y_or_Dec_velocity'] = [np.nan, 1.1, np.nan]
    comparison['ephemeris_file'] = mixed_ephemeris_list

    for colname in nan_cols:
        assert np.all(np.isnan(comparison[colname]) == np.isnan(mpt.table[colname]))
    assert np.all(comparison['ephemeris_file'] == mpt.table['ephemeris_file'])

    mpt = catalog_generator.MovingPointSourceCatalog(ra=ra_list, dec=dec_list, ra_velocity=ra_vel_list, dec_velocity=dec_vel_list)
    mpt.add_magnitude_column(mags, magnitude_system='abmag', instrument='nircam', filter_name='f444w')

    comparison['x_or_RA'] = ra_list
    comparison['y_or_Dec'] = dec_list
    comparison['x_or_RA_velocity'] = ra_vel_list
    comparison['y_or_Dec_velocity'] = dec_vel_list
    comparison['ephemeris_file'] = ['None'] * 3

    for colname in nan_cols + ['ephemeris_file']:
        assert np.all(comparison[colname] == mpt.table[colname])


#@pytest.mark.skip(reason='Travis and astroquery not getting along')
def test_2mass_catalog_generation():
    """Test the generation of a catalog from a 2MASS query
    """
    two_mass_mirage, two_mass_query_results = create_catalog.get_2MASS_ptsrc_catalog(80.4, -69.8, 120)
    comparison_file = os.path.join(TEST_DATA_DIR, 'catalog_generation/TwoMass_test.cat')
    comparison_table = ascii.read(comparison_file)
    assert two_mass_query_results.colnames == ['ra', 'dec', 'clon', 'clat', 'err_maj', 'err_min', 'err_ang',
                                               'designation', 'j_m', 'j_cmsig', 'j_msigcom', 'j_snr', 'h_m',
                                               'h_cmsig', 'h_msigcom', 'h_snr', 'k_m', 'k_cmsig', 'k_msigcom',
                                               'k_snr', 'ph_qual', 'rd_flg', 'bl_flg', 'cc_flg', 'ndet', 'gal_contam',
                                               'mp_flg', 'hemis', 'xdate', 'scan', 'glon', 'glat', 'a', 'dist_opt',
                                               'phi_opt', 'b_m_opt', 'vr_m_opt', 'nopt_mchs', 'ext_key', 'j_h', 'h_k', 'j_k']
    assert two_mass_mirage.table.colnames == ['index', 'x_or_RA', 'y_or_Dec', '2mass_j_m_magnitude',
                                              '2mass_h_m_magnitude', '2mass_k_m_magnitude']
    assert len(two_mass_mirage.ra) == len(comparison_table['x_or_RA'].data)


def test_catalog_combination():
    """Test the combination of two existing catalogs into one
    """
    test_file = os.path.join(TEST_DATA_DIR, 'catalog_generation/TwoMass_test.cat')
    input_cat = ascii.read(test_file)

    # Create catalog 1
    test_cat_1 = catalog_generator.PointSourceCatalog(ra=input_cat['x_or_RA'], dec=input_cat['y_or_Dec'])
    test_cat_1.add_magnitude_column(input_cat['2mass_j_m_magnitude'], magnitude_system='vegamag', instrument='nircam', filter_name='f200w')
    test_cat_1.add_magnitude_column(input_cat['2mass_h_m_magnitude'], magnitude_system='vegamag', instrument='nircam', filter_name='f090w')
    orig_length = len(test_cat_1.ra)
    orig_cat_1 = copy.deepcopy(test_cat_1)

    # Create catalog 2
    test_cat_2 = catalog_generator.PointSourceCatalog(ra=input_cat['x_or_RA'] + 1., dec=input_cat['y_or_Dec'] + 1.)
    test_cat_2.add_magnitude_column(input_cat['2mass_h_m_magnitude'], magnitude_system='vegamag', instrument='nircam', filter_name='f200w')
    test_cat_2.add_magnitude_column(input_cat['2mass_h_m_magnitude'], magnitude_system='vegamag', instrument='nircam', filter_name='f090w')

    # Combine
    test_cat_1.add_catalog(test_cat_2)


    cat2_length = len(test_cat_2.ra)

    assert all(test_cat_1.ra[(0-cat2_length):] == test_cat_2.ra)
    assert all(test_cat_1.dec[(0-cat2_length):] == test_cat_2.dec)
    assert all(test_cat_1.ra[0: orig_length] == orig_cat_1.ra)
    assert all(test_cat_1.dec[0: orig_length] == orig_cat_1.dec)


@pytest.mark.skip(reason="Bug with the Besancon model in astroquery.")
def test_besancon_generation():
    """Test the creation of a catalog from a Besancon query
    """
    b_cat, b_query = create_catalog.besancon(80.4, -69.8, 120, coords='ra_dec', email='hilbert@stsci.edu')
    assert len(b_cat.table) > 0
    assert b_cat.table.colnames == ['index', 'x_or_RA', 'y_or_Dec', 'besancon_v_magnitude',
                                    'besancon_j_magnitude', 'besancon_h_magnitude', 'besancon_k_magnitude',
                                    'besancon_l_magnitude']


@pytest.mark.skip(reason="Bug with the Besancon model in astroquery.")
def test_2mass_plus_besaoncon_convenience_function():
    """Test the convenience function that queries 2MASS and Besancon
    and combines the two catalogs
    """
    twomass, query = create_catalog.get_2MASS_ptsrc_catalog(80.4, -69.8, 120)
    both = create_catalog.twoMASS_plus_background(80.4, -69.8, 120, kmag_limits=(16, 25),
                                                  email='hilbert@stsci.edu')
    assert len(twomass.table) < len(both.table)


@pytest.mark.skip(reason="Bug with the Besancon model in astroquery.")
def test_for_proposal():
    """Test the creation of source catalogs from a proposal"""
    xml = os.path.join(TEST_DATA_DIR, 'NIRCam/targets_with_large_separation.xml')
    pointing = xml.replace('.xml', '.pointing')
    output_directory = os.path.join(TEST_DATA_DIR, 'catalog_generation')
    p_cat, g_cat, p_name, g_name, p_map, g_map = create_catalog.for_proposal(xml, pointing,
                                                                             point_source=True,
                                                                             extragalactic=True,
                                                                             catalog_splitting_threshold=0.12,
                                                                             email='hilbert@stsci.edu',
                                                                             out_dir=output_directory,
                                                                             save_catalogs=True,
                                                                             besancon_seed=2345,
                                                                             galaxy_seed=2346)

    ptsrc_map = {'001': 'ptsrc_for_targets_with_large_separation_observations_001.cat',
                 '002': 'ptsrc_for_targets_with_large_separation_observations_002.cat',
                 '003': 'ptsrc_for_targets_with_large_separation_observations_003.cat'}

    galaxy_map = {'001': 'galaxies_for_targets_with_large_separation_observations_001.cat',
                  '002': 'galaxies_for_targets_with_large_separation_observations_002.cat',
                  '003': 'galaxies_for_targets_with_large_separation_observations_003.cat'}

    ptsrc_name = []
    galaxy_name = []
    for key in ptsrc_map:
        ptsrc_name.append(os.path.join(output_directory, ptsrc_map[key]))
        galaxy_name.append(os.path.join(output_directory, galaxy_map[key]))

    assert ptsrc_map == p_map
    assert galaxy_map == g_map
    assert ptsrc_name == p_name
    assert galaxy_name == g_name

    truth_ptsrc = ['ptsrc_1.cat', 'ptsrc_2.cat', 'ptsrc_3.cat']
    truth_galaxy = ['galaxy_1.cat', 'galaxy_2.cat', 'galaxy_3.cat']

    for pcat_new, pcat in zip(p_cat, truth_ptsrc):
        truth_file = os.path.join(TEST_DATA_DIR, 'catalog_generation/{}'.format(pcat))
        truth = ascii.read(truth_file)

        for col in pcat_new.table.colnames:
            assert all(pcat_new.table[col].data == truth[col].data)

    for gcat_new, gcat in zip(g_cat, truth_galaxy):
        truth_file = os.path.join(TEST_DATA_DIR, 'catalog_generation/{}'.format(gcat))
        truth = ascii.read(truth_file)

        for col in gcat_new.table.colnames:
            print('DATA: {}'.format(col))
            print(gcat_new.table[col].data)
            print('')
            print(truth[col].data)
            if 'magnitude' in col:
                assert np.allclose(gcat_new.table[col].data, truth[col].data, rtol=0, atol=0.5)
            else:
                assert all(gcat_new.table[col].data == truth[col].data)

    # Remove the catalog files just produced
    for del_ptsrc, del_galaxy in zip(ptsrc_name, galaxy_name):
        print(os.path.join(output_directory, del_ptsrc), os.path.join(output_directory, del_galaxy))
        os.remove(os.path.join(output_directory, del_ptsrc))
        os.remove(os.path.join(output_directory, del_galaxy))
        cwd = os.getcwd()
        try:
            os.remove(os.path.join(cwd, 'observation_list.yaml'))
            os.remove(os.path.join(cwd, 'expand_for_detectors.csv'))
            os.remove(os.path.join(cwd, 'Observation_table_for_targets_with_large_separation.csv'))
        except:
            pass


def test_get_all_catalogs():
    """Test the wrapper that queries anc combines catalogs from all sources"""
    ra = 80.4
    dec = -69.8
    width = 120.
    ins = 'NIRCAM'
    filters = ['F150W', 'F356W', 'F444W', 'F480M']
    besancon_file = os.path.join(TEST_DATA_DIR, 'catalog_generation/besancon_example.cat')

    cat, headers = create_catalog.get_all_catalogs(ra, dec, width, instrument=ins, filters=filters,
                                                   besancon_catalog_file=besancon_file)
    comparison_file = os.path.join(TEST_DATA_DIR, 'catalog_generation/get_all_catalogs_truth.cat')
    comparison_data = ascii.read(comparison_file)

    # Don't check the detailed results, since any changes to Besancon/WISE/GAIA/2MASS queries
    # could cause the test to fail without there being a problem with Mirage. Instead check that
    # all the expected columns are present.
    for col in comparison_data.colnames:
        assert col in cat.table.colnames, f"{col} column not in the returned catalog."


#@pytest.mark.skip(reason="Repeated HTML response errors.")
def test_gaia_query():
    """Test the GAIA query and transformation into a Mirage-format catalog"""
    ra = 80.4
    dec = -69.8
    box_width = 200.
    cat, query, gaia_2mass_cross, gaia_wise_cross = create_catalog.get_gaia_ptsrc_catalog(ra, dec, box_width)
    assert cat.table.colnames == ['index', 'x_or_RA', 'y_or_Dec', 'gaia_phot_g_mean_mag_magnitude',
                                  'gaia_phot_bp_mean_mag_magnitude', 'gaia_phot_rp_mean_mag_magnitude']


def test_random_ra_dec_values():
    """Test the random RA, Dec value generator used when getting Besancon
    sources"""
    num_stars = 10
    ra_min = 1.0
    ra_max = 1.1
    dec_min = 40.
    dec_max = 40.1
    ra1, dec1 = create_catalog.generate_ra_dec(num_stars, ra_min, ra_max, dec_min, dec_max,
                                               seed=37465)
    ra2, dec2 = create_catalog.generate_ra_dec(num_stars, ra_min, ra_max, dec_min, dec_max,
                                               seed=37465)
    assert np.all(ra1 == ra2)
    assert np.all(dec1 == dec2)
    assert np.min(ra1) >= ra_min
    assert np.max(ra1) <= ra_max
    assert np.min(dec1) >= dec_min
    assert np.max(dec1) <= dec_max


def test_cat_from_file():
    """Test function for reading ascii catalogs into catalog objects"""
    data_path = os.path.join(TEST_DATA_DIR, 'catalog_generation/')
    catalogs = {'ptsrc_1.cat': ('point_source', catalog_generator.PointSourceCatalog),
                'galaxy_1.cat': ('galaxy', catalog_generator.GalaxyCatalog),
                'extended_test.cat': ('extended', catalog_generator.ExtendedCatalog),
                'moving_point_source_test.cat': ('moving_point_source', catalog_generator.MovingPointSourceCatalog),
                'moving_sersic_source_test.cat': ('moving_sersic', catalog_generator.MovingSersicCatalog),
                'moving_extended_source_test.cat': ('moving_extended', catalog_generator.MovingExtendedCatalog),
                'nonsidereal_source_test.cat': ('non_sidereal', catalog_generator.NonSiderealCatalog)}
    for cat_name in catalogs:
        cat_path = os.path.join(data_path, cat_name)
        cat_object = catalog_generator.cat_from_file(cat_path, catalogs[cat_name][0])
        assert isinstance(cat_object, catalogs[cat_name][1])

if not ON_GITHUB:
    os.environ['MIRAGE_DATA'] = orig_mirage_data


def test_transform_johnson_to_jwst():
    jcat = Table()
    jcat['V'] = [10, 10, 10, 10]
    jcat['J'] = [10, 11, 12, 13]
    jcat['H'] = [10, 12, 14, 16]
    jcat['K'] = [10, 13, 16, 19]
    jcat['Av'] = [0.05, 0.05, 0.05, 0.05]

    filters = ['nircam_f090w_clear_magnitude', 'niriss_f090w_magnitude', 'fgs_guider1_magnitude']
    transform = create_catalog.transform_johnson_to_jwst(jcat, filters)

    truth = copy.deepcopy(jcat)
    truth['nircam_f090w_clear_magnitude'] = [10.015026, 11.381326, 12.881326, 14.381326]
    truth['niriss_f090w_magnitude'] = [10.016299, 11.379, 12.879, 14.379]

    assert jcat.colnames == truth.colnames
    assert np.allclose(jcat['nircam_f090w_clear_magnitude'], truth['nircam_f090w_clear_magnitude'])
    assert np.allclose(jcat['niriss_f090w_magnitude'], truth['niriss_f090w_magnitude'])


def test_johnson_catalog_to_mirage_catalog():
    input_catalog = os.path.join(TEST_DATA_DIR, 'catalog_generation/besancon_example.cat')
    filters = {'nircam': ['F090W/CLEAR', 'F322W2/F323N'],
               'niriss': ['F090W', 'F277W'],
               'fgs': ['GUIDER1']
               }
    transformed = create_catalog.johnson_catalog_to_mirage_catalog(input_catalog, filters)

    in_cat = ascii.read(input_catalog)
    ra_vals = in_cat['RAJ2000'].data
    dec_vals = in_cat['DECJ2000'].data
    nrc_090 = [13.484183, 12.336812, 14.283301, 12.838769]
    nrc_323 = [13.58035, 12.420761, 14.363094, 12.93159]
    nis_090 = [13.483351, 12.335981, 14.282558, 12.837949]
    nis_277 = [13.572666, 12.4141, 14.356856, 12.924085]
    fgs_g1 = [13.490332, 12.341099, 14.287583, 12.844421]
    truth = catalog_generator.PointSourceCatalog(ra=ra_vals, dec=dec_vals)
    truth.add_magnitude_column(nrc_090, column_name='nircam_f090w_clear_magnitude')
    truth.add_magnitude_column(nrc_323, instrument='nircam', filter_name='F323N')
    truth.add_magnitude_column(nis_090, column_name='niriss_f090w_magnitude')
    truth.add_magnitude_column(nis_277, instrument='niriss', filter_name='F277W')
    truth.add_magnitude_column(fgs_g1, column_name='fgs_guider1_magnitude')

    for col in truth.table.colnames:
        assert np.allclose(truth.table[col].data, transformed.table[col].data)


def test_catalog_index_check():
    """Test the code that assures no overlap of source indexes
    """
    #Make point source catalog
    ra = np.zeros(10) + 80.0
    dec = np.zeros(10) - 69.8
    mags = np.zeros(10) + 19.

    ptsrc = catalog_generator.PointSourceCatalog(ra=ra, dec=dec, starting_index=1)
    ptsrc.add_magnitude_column(mags, instrument='nircam', filter_name='f090w')
    ptsrc_output_file = os.path.join(TEST_DATA_DIR, 'catalog_generation/ptsrc_test.cat')
    ptsrc.save(ptsrc_output_file)

    # Make galaxy catalog
    indexes = np.arange(10, 20)
    radius = np.zeros(10) + 0.5
    ellip = np.zeros(10) + 0.45
    posang = np.zeros(10) + 27.
    sersic = np.zeros(10) + 3.3

    gal_output_file = os.path.join(TEST_DATA_DIR, 'catalog_generation/galaxy_test.cat')
    gal = catalog_generator.GalaxyCatalog(ra=ra, dec=dec, ellipticity=ellip, radius=radius,
                                          sersic_index=sersic, position_angle=posang, radius_units='arcsec',
                                          starting_index=11)
    gal.add_magnitude_column(mags, instrument='nircam', filter_name='f090w')
    gal.save(gal_output_file)

    # Make moving point source catalog
    ra_vel_list = [0.2, 0.3, 0.4]
    dec_vel_list = [1.0, 1.1, 1.2]
    ephemeris_list = ['e1.txt', 'e2.txt', 'e3.txt']

    mpt = catalog_generator.MovingPointSourceCatalog(ra=ra[0:3], dec=dec[0:3], ra_velocity=ra_vel_list, dec_velocity=dec_vel_list,
                                                 ephemeris_file=ephemeris_list, starting_index=21)
    mpt.add_magnitude_column(mags[0:3], magnitude_system='abmag', instrument='nircam', filter_name='f444w')
    mpt_output_file = os.path.join(TEST_DATA_DIR, 'catalog_generation/moving_ptsrc_test.cat')
    mpt.save(mpt_output_file)


    # First a case where there is no overlap
    cat_list = [ptsrc_output_file, gal_output_file, mpt_output_file]
    overlap = utils.catalog_index_check(cat_list)
    assert overlap == (False, 23)

    # Now a case where there is overlap
    gal = catalog_generator.GalaxyCatalog(ra=ra, dec=dec, ellipticity=ellip, radius=radius,
                                          sersic_index=sersic, position_angle=posang, radius_units='arcsec',
                                          starting_index=4)
    gal.add_magnitude_column(mags, instrument='nircam', filter_name='f090w')
    gal.save(gal_output_file)

    overlap = utils.catalog_index_check(cat_list)
    assert overlap == (True, 23)

    # Overlap, but indexes between catalogs are out of order
    ptsrc = catalog_generator.PointSourceCatalog(ra=ra, dec=dec, starting_index=25)
    ptsrc.add_magnitude_column(mags, instrument='nircam', filter_name='f090w')
    ptsrc.save(ptsrc_output_file)

    gal = catalog_generator.GalaxyCatalog(ra=ra, dec=dec, ellipticity=ellip, radius=radius,
                                          sersic_index=sersic, position_angle=posang, radius_units='arcsec',
                                          starting_index=1)
    gal.add_magnitude_column(mags, instrument='nircam', filter_name='f090w')
    gal.save(gal_output_file)

    mpt = catalog_generator.MovingPointSourceCatalog(ra=ra[0:3], dec=dec[0:3], ra_velocity=ra_vel_list, dec_velocity=dec_vel_list,
                                                 ephemeris_file=ephemeris_list, starting_index=29)
    mpt.add_magnitude_column(mags[0:3], magnitude_system='abmag', instrument='nircam', filter_name='f444w')
    mpt_output_file = os.path.join(TEST_DATA_DIR, 'catalog_generation/moving_ptsrc_test.cat')
    mpt.save(mpt_output_file)

    overlap = utils.catalog_index_check(cat_list)
    assert overlap == (True, 34)

    # No overlap, but indexes between catalogs are out of order
    ptsrc = catalog_generator.PointSourceCatalog(ra=ra, dec=dec, starting_index=25)
    ptsrc.add_magnitude_column(mags, instrument='nircam', filter_name='f090w')
    ptsrc.save(ptsrc_output_file)

    gal = catalog_generator.GalaxyCatalog(ra=ra, dec=dec, ellipticity=ellip, radius=radius,
                                          sersic_index=sersic, position_angle=posang, radius_units='arcsec',
                                          starting_index=1)
    gal.add_magnitude_column(mags, instrument='nircam', filter_name='f090w')
    gal.save(gal_output_file)

    mpt = catalog_generator.MovingPointSourceCatalog(ra=ra[0:3], dec=dec[0:3], ra_velocity=ra_vel_list, dec_velocity=dec_vel_list,
                                                 ephemeris_file=ephemeris_list, starting_index=36)
    mpt.add_magnitude_column(mags[0:3], magnitude_system='abmag', instrument='nircam', filter_name='f444w')
    mpt_output_file = os.path.join(TEST_DATA_DIR, 'catalog_generation/moving_ptsrc_test.cat')
    mpt.save(mpt_output_file)

    overlap = utils.catalog_index_check(cat_list)
    assert overlap == (False, 38)

    # Clean up
    os.remove(ptsrc_output_file)
    os.remove(gal_output_file)
    os.remove(mpt_output_file)


def test_determine_used_cats():
    """Test the evaluation of which source catalogs are used for an obs mode
    """
    catalogs = {}
    catalogs['pointsource'] = 'ptsrc.cat'
    catalogs['galaxyListFile'] = 'galaxies.cat'
    catalogs['extended'] = 'extended.cat'
    catalogs['movingTargetList'] = 'moving_ptsrc.cat'
    catalogs['movingTargetSersic'] = 'moving_galaxies.cat'
    catalogs['movingTargetExtended'] = 'moving_extended.cat'
    catalogs['movingTargetToTrack'] = 'nonsidereal.cat'
    catalogs['tso_imaging_catalog'] = 'img_tso.cat'
    catalogs['tso_grism_catalog'] = 'grism_tso.cat'

    img_cats = utils.determine_used_cats('imaging', catalogs)
    img_truth = ['ptsrc.cat', 'galaxies.cat', 'extended.cat', 'moving_ptsrc.cat',
                 'moving_galaxies.cat', 'moving_extended.cat', 'nonsidereal.cat']
    assert img_cats == img_truth

    wfss_cats = utils.determine_used_cats('wfss', catalogs)
    wfss_truth = ['ptsrc.cat', 'galaxies.cat', 'extended.cat']
    assert wfss_cats == wfss_truth

    tsoimg_cats = utils.determine_used_cats('ts_imaging', catalogs)
    tsoimg_truth = ['ptsrc.cat', 'galaxies.cat', 'extended.cat', 'img_tso.cat']
    assert tsoimg_cats == tsoimg_truth

    tsogrism_cats = utils.determine_used_cats('ts_grism', catalogs)
    tsogrism_truth = ['ptsrc.cat', 'galaxies.cat', 'extended.cat', 'grism_tso.cat']
    assert tsogrism_cats == tsogrism_truth


def test_nonsidereal_extended_source_catalog():
    """Test the creation of non-sidereal source catalogs
    """
    filename = 'Jupiter_stamp.fits'
    ephemeris_file = 'horizons_results.txt'
    ns = catalog_generator.NonSiderealCatalog(object_type=['extended'], ephemeris_file=[ephemeris_file],
                                              filenames=filename)
    ns.add_magnitude_column([2.0], magnitude_system='abmag', instrument='nircam', filter_name='f360m')
    ns.add_magnitude_column([2.0], magnitude_system='abmag', instrument='nircam', filter_name='f405n')
    expected_cols = ['index', 'object', 'x_or_RA', 'y_or_Dec', 'x_or_RA_velocity', 'y_or_Dec_velocity',
                     'filename', 'pos_angle', 'nircam_f360m_clear_magnitude', 'ephemeris_file',
                     'nircam_f444w_f405n_magnitude']
    for col in expected_cols:
        assert col in ns.table.colnames
    assert ns.table['object'] == 'extended'


def test_nonsidereal_sersic_source_catalog():
    """Test the creation of non-sidereal source catalogs
    """
    ephemeris_file = 'horizons_results.txt'
    ns = catalog_generator.NonSiderealCatalog(object_type=['sersic'], ephemeris_file=[ephemeris_file], ellipticity=[0.02], radius=[0.3], sersic_index=[1.4],
                 position_angle=[49.])
    ns.add_magnitude_column([2.0], magnitude_system='abmag', instrument='nircam', filter_name='f360m')
    ns.add_magnitude_column([2.0], magnitude_system='abmag', instrument='nircam', filter_name='f405n')
    expected_cols = ['index', 'object', 'x_or_RA', 'y_or_Dec', 'x_or_RA_velocity', 'y_or_Dec_velocity',
                     'sersic_index', 'ellipticity', 'pos_angle', 'radius', 'nircam_f360m_clear_magnitude',
                     'ephemeris_file', 'nircam_f444w_f405n_magnitude']
    for col in expected_cols:
        assert col in ns.table.colnames
    assert ns.table['object'] == 'sersic'


def test_nonsidereal_point_source_catalog():
    """Test the creation of non-sidereal source catalogs
    """
    ephemeris_file = 'horizons_results.txt'
    ns = catalog_generator.NonSiderealCatalog(object_type=['pointSource'], ephemeris_file=[ephemeris_file])
    ns.add_magnitude_column([2.0], magnitude_system='abmag', instrument='nircam', filter_name='f360m')
    ns.add_magnitude_column([2.0], magnitude_system='abmag', instrument='nircam', filter_name='f405n')
    expected_cols = ['index', 'object', 'x_or_RA', 'y_or_Dec', 'x_or_RA_velocity', 'y_or_Dec_velocity',
                     'nircam_f360m_clear_magnitude', 'ephemeris_file', 'nircam_f444w_f405n_magnitude']
    for col in expected_cols:
        assert col in ns.table.colnames
    assert ns.table['object'] == 'pointSource'

