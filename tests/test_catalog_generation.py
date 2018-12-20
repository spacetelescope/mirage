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

from mirage.catalogs import catalog_generator
from mirage.catalogs import create_catalog

TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), 'catalogs')


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
    assert ptsrc.table.colnames == ['index', 'x_or_RA', 'y_or_Dec', 'nircam_f090w_magnitude', 'nircam_f200w_magnitude']


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

    output_file = os.path.join(TEST_DATA_DIR, 'galaxy_test.cat')
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


def test_2mass_catalog_generation():
    """Test the generation of a catalog from a 2MASS query
    """
    two_mass_mirage, two_mass_query_results = create_catalog.get_2MASS_ptsrc_catalog(80.4, -69.8, 120)
    comparison_file = os.path.join(TEST_DATA_DIR, 'TwoMass_test.cat')
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


def test_besancon_generation():
    """Test the creation of a catalog from a Besancon query
    """
    b_cat, b_query = create_catalog.besancon(80.4, -69.8, 120, coords='ra_dec', email='hilbert@stsci.edu')
    assert len(b_cat.table) > 0
    assert b_cat.table.colnames == ['index', 'x_or_RA', 'y_or_Dec', 'besancon_v_magnitude',
                                    'besancon_j_magnitude', 'besancon_h_magnitude', 'besancon_k_magnitude',
                                    'besancon_l_magnitude']


def test_2mass_plus_besaoncon_convenience_function():
    """Test the convenience function that queries 2MASS and Besancon
    and combines the two catalogs
    """
    twomass, query = create_catalog.get_2MASS_ptsrc_catalog(80.4, -69.8, 120)
    both = create_catalog.twoMASS_plus_background(80.4, -69.8, 120, kmag_limits=(16, 25),
                                                  email='hilbert@stsci.edu')
    assert len(twomass.table) < len(both.table)
