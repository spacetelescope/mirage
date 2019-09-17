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
import pytest

from mirage.catalogs import catalog_generator
from mirage.catalogs import create_catalog

TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), 'test_data/')

# Determine if tests are being run on Travis
ON_TRAVIS = 'travis' in os.path.expanduser('~')

if not ON_TRAVIS:
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

    assert all(ptsrc.ra == ra)
    assert all(ptsrc.dec == dec)
    assert ptsrc.table.colnames == ['index', 'x_or_RA', 'y_or_Dec', 'nircam_f090w_magnitude',
                                    'nircam_f200w_magnitude']


def test_galaxy_catalog_creation():
    """Test the creation of a basic galaxy object catalog
    """
    ra = np.zeros(10) + 80.0
    dec = np.zeros(10) - 69.8
    mags = np.zeros(10) + 19.
    radius = np.zeros(10) + 0.5
    ellip = np.zeros(10) + 0.45
    posang = np.zeros(10) + 27.
    sersic = np.zeros(10) + 3.3

    output_file = os.path.join(TEST_DATA_DIR, 'catalog_generation/galaxy_test.cat')
    gal = catalog_generator.GalaxyCatalog(ra=ra, dec=dec, ellipticity=ellip, radius=radius,
                                          sersic_index=sersic, position_angle=posang, radius_units='arcsec')
    gal.add_magnitude_column(mags, instrument='nircam', filter_name='f090w')
    gal.save(output_file)
    as_read_in = ascii.read(output_file)

    assert gal.table.colnames == ['index', 'x_or_RA', 'y_or_Dec', 'pos_angle', 'sersic_index',
                                  'ellipticity', 'radius', 'nircam_f090w_magnitude']
    assert all(gal.radius == radius)
    assert all(gal.ellipticity == ellip)
    assert all(gal.position_angle == posang)
    assert all(gal.sersic_index == sersic)
    assert all(gal.table == as_read_in)
    os.remove(output_file)


def test_2mass_catalog_generation():
    """Test the generation of a catalog from a 2MASS query
    """
    two_mass_mirage, two_mass_query_results = create_catalog.get_2MASS_ptsrc_catalog(80.4, -69.8, 120)
    comparison_file = os.path.join(TEST_DATA_DIR, 'catalog_generation/TwoMass_test.cat')
    comparison_table = ascii.read(comparison_file)
    assert all(two_mass_mirage.ra == comparison_table['x_or_RA'].data)


def test_catalog_combination():
    """Test the combination of two existing catalogs into one
    """
    two_mass, query_results = create_catalog.get_2MASS_ptsrc_catalog(80.4, -69.8, 120)
    fake_nrc_mags = np.zeros(len(two_mass.ra)) + 16.5
    two_mass.add_magnitude_column(fake_nrc_mags, instrument='nircam', filter_name='f480m',
                                  magnitude_system='vegamag')
    orig_two_mass = copy.deepcopy(two_mass)

    two_mass2, query_results2 = create_catalog.get_2MASS_ptsrc_catalog(70.4, -69.8, 120)
    fake_nrc_mags = np.zeros(len(two_mass2.ra)) + 12.12
    two_mass2.add_magnitude_column(fake_nrc_mags, instrument='niriss', filter_name='f200w',
                                   magnitude_system='vegamag')
    two_mass.add_catalog(two_mass2)

    two_mass2_length = len(two_mass2.ra)
    orig_two_mass_length = len(orig_two_mass.ra)
    assert all(two_mass.ra[(0-two_mass2_length):] == two_mass2.ra)
    assert all(two_mass.dec[(0-two_mass2_length):] == two_mass2.dec)
    assert all(two_mass.ra[0: orig_two_mass_length] == orig_two_mass.ra)
    assert all(two_mass.dec[0: orig_two_mass_length] == orig_two_mass.dec)

    for key in two_mass.magnitudes.keys():
        assert key in ['2mass_j_m_magnitude', '2mass_h_m_magnitude', '2mass_k_m_magnitude',
                       'nircam_f480m_magnitude', 'niriss_f200w_magnitude']


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


@pytest.mark.skip(reason="Bug with the Besancon model in astroquery.")
def test_get_all_catalogs():
    """Test the wrapper that queries anc combines catalogs from all sources"""
    ra = 80.4
    dec = -69.8
    width = 120.
    ins = 'NIRCAM'
    filters = ['F150W', 'F356W', 'F444W', 'F480M']

    cat, headers = create_catalog.get_all_catalogs(ra, dec, width, kmag_limits=(13, 29),
                                                   email='hilbert@stsci.edu', instrument=ins, filters=filters,
                                                   besancon_seed=1234)
    comparison_file = os.path.join(TEST_DATA_DIR, 'catalog_generation/get_all_catalogs.cat')
    comparison_data = ascii.read(comparison_file)

    # Note that if Besancon/WISE/GAIA/2MASS query results change, this will
    # fail without there being a problem with Mirage.
    for col in cat.table.colnames:
        try:
            assert all(cat.table[col].data == comparison_data[col].data), \
                "Retrieved catalog does not match expected."
        except TypeError:
            assert False, "Retrieved catalog does not match expected."


def test_gaia_query():
    """Test the GAIA query and transformation into a Mirage-format catalog"""
    ra = 80.4
    dec = -69.8
    box_width = 200.
    cat, query, gaia_2mass_cross, gaia_wise_cross = create_catalog.get_gaia_ptsrc_catalog(ra, dec, box_width)
    assert len(cat.table) == 1153
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

if not ON_TRAVIS:
    os.environ['MIRAGE_DATA'] = orig_mirage_data
