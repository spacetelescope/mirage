"""Test the functions provided by spectra_from_catalog.py

Authors
-------
    - Bryan Hilbert

Use
---
    >>> pytest -s test_spectra_from_catalog.py


"""
import os
import numpy as np

import astropy.units as u

from mirage.catalogs import spectra_from_catalog as spec
from mirage.catalogs import hdf5_catalog as hdf5
from mirage.utils.constants import FLAMBDA_CGS_UNITS, FNU_CGS_UNITS


TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), 'test_data/hdf5_catalogs')


def test_get_filter_info():
    nrc = spec.get_filter_info(['nircam_f444w_magnitude'], 'abmag')
    assert nrc == {'nircam_f444w_magnitude': (6.6928e-22 * FLAMBDA_CGS_UNITS, 4.3637e-31 * FNU_CGS_UNITS, 27.3004,
                                              4.4212 * u.micron)}

    nis = spec.get_filter_info(['niriss_f200w_magnitude'], 'vegamag')
    assert nis == {'niriss_f200w_magnitude': (2.173398e-20 * FLAMBDA_CGS_UNITS, 2.879494e-31 * FNU_CGS_UNITS,
                                              26.04898, 1.9930 * u.micron)}

    fgs = spec.get_filter_info(['fgs_magnitude'], 'stmag')
    assert fgs == {'fgs_magnitude': (1.593395e-21 * FLAMBDA_CGS_UNITS, 3.324736e-32 * FNU_CGS_UNITS, 33.39426,
                                     2.5011 * u.micron)}


def test_hdf5_file_input():
    """Case where an hdf5 file is input. One of these also includes a
    normalized spectrum"""
    catfile = os.path.join(TEST_DATA_DIR, 'point_sources.cat')
    output_hdf5 = os.path.join(TEST_DATA_DIR, 'all_spectra.hdf5')
    sed_file = os.path.join(TEST_DATA_DIR, 'sed_file_with_normalized_dataset.hdf5')
    sed_catalog = spec.make_all_spectra(catfile, input_spectra_file=sed_file,
                                        normalizing_mag_column='nircam_f444w_magnitude',
                                        output_filename=output_hdf5)

    comparison = hdf5.open(os.path.join(TEST_DATA_DIR, 'output_spec_from_hdf5_input_including_normalized.hdf5'))
    constructed = hdf5.open(sed_catalog)
    for key in comparison:
        assert key in constructed.keys()
        assert all(comparison[key]["wavelengths"].value == constructed[key]["wavelengths"].value)
        assert all(comparison[key]["fluxes"].value == constructed[key]["fluxes"].value)
        assert comparison[key]["wavelengths"].unit == constructed[key]["wavelengths"].unit
        assert comparison[key]["fluxes"].unit == constructed[key]["fluxes"].unit

    cat_base = catfile.split('.')[0]
    outbase = cat_base + '_with_flambda.cat'
    flambda_output_catalog = os.path.join(TEST_DATA_DIR, outbase)
    os.remove(flambda_output_catalog)
    os.remove(sed_catalog)


def test_manual_inputs():
    """Case where spectra are input manually along side ascii catalog"""
    # Test case where spectra are input manually
    catfile = os.path.join(TEST_DATA_DIR, 'point_sources.cat')
    output_hdf5 = os.path.join(TEST_DATA_DIR, 'all_spectra.hdf5')
    hdf5file = os.path.join(TEST_DATA_DIR, 'sed_file_with_normalized_dataset.hdf5')
    sed_dict = hdf5.open(hdf5file)
    sed_catalog = spec.make_all_spectra(catfile, input_spectra=sed_dict,
                                        normalizing_mag_column='nircam_f444w_magnitude',
                                        output_filename=output_hdf5)
    constructed = hdf5.open(sed_catalog)
    comparison = hdf5.open(os.path.join(TEST_DATA_DIR, 'output_spec_from_manual_input.hdf5'))
    for key in comparison:
        assert key in constructed.keys()
        assert all(comparison[key]["wavelengths"].value == constructed[key]["wavelengths"].value)
        assert all(comparison[key]["fluxes"].value == constructed[key]["fluxes"].value)
        assert comparison[key]["wavelengths"].unit == constructed[key]["wavelengths"].unit
        assert comparison[key]["fluxes"].unit == constructed[key]["fluxes"].unit

    cat_base = catfile.split('.')[0]
    outbase = cat_base + '_with_flambda.cat'
    flambda_output_catalog = os.path.join(TEST_DATA_DIR, outbase)
    os.remove(flambda_output_catalog)
    os.remove(sed_catalog)


def test_manual_plus_file_inputs():
    """Case where spectra are input via hdf5 file as well as manually"""
    catfile = os.path.join(TEST_DATA_DIR, 'point_sources.cat')
    sed_file = os.path.join(TEST_DATA_DIR, 'sed_file_with_normalized_dataset.hdf5')
    output_hdf5 = os.path.join(TEST_DATA_DIR, 'all_spectra.hdf5')
    manual_sed = {}
    manual_sed[7] = {"wavelengths": [0.9, 1.4, 1.9, 3.5, 5.1]*u.micron,
                     "fluxes": [1e-17, 1.1e-17, 1.5e-17, 1.4e-17, 1.1e-17] * FLAMBDA_CGS_UNITS}
    sed_catalog = spec.make_all_spectra(catfile, input_spectra=manual_sed, input_spectra_file=sed_file,
                                        normalizing_mag_column='nircam_f444w_magnitude',
                                        output_filename=output_hdf5)
    comparison = hdf5.open(os.path.join(TEST_DATA_DIR, 'output_spec_from_file_plus_manual_input.hdf5'))
    constructed = hdf5.open(output_hdf5)
    for key in comparison:
        assert key in constructed.keys()
        assert all(comparison[key]["wavelengths"].value == constructed[key]["wavelengths"].value)
        assert all(comparison[key]["fluxes"].value == constructed[key]["fluxes"].value)
        assert comparison[key]["wavelengths"].unit == constructed[key]["wavelengths"].unit
        assert comparison[key]["fluxes"].unit == constructed[key]["fluxes"].unit

    cat_base = catfile.split('.')[0]
    outbase = cat_base + '_with_flambda.cat'
    flambda_output_catalog = os.path.join(TEST_DATA_DIR, outbase)
    os.remove(flambda_output_catalog)
    os.remove(sed_catalog)


def test_multiple_mag_columns():
    """Case where ascii catalog with multiple magnitude columns is input"""
    catfile = os.path.join(TEST_DATA_DIR, 'point_sources.cat')
    output_hdf5 = os.path.join(TEST_DATA_DIR, 'all_spectra.hdf5')
    sed_catalog = spec.make_all_spectra(catfile, output_filename=output_hdf5)
    constructed = hdf5.open(sed_catalog)
    comparison = hdf5.open(os.path.join(TEST_DATA_DIR, 'output_spec_from_multiple_filter.hdf5'))
    for key in comparison:
        assert key in constructed.keys()
        assert all(comparison[key]["wavelengths"].value == constructed[key]["wavelengths"].value)
        assert all(comparison[key]["fluxes"].value == constructed[key]["fluxes"].value)
        assert comparison[key]["wavelengths"].unit == constructed[key]["wavelengths"].unit
        assert comparison[key]["fluxes"].unit == constructed[key]["fluxes"].unit

    cat_base = catfile.split('.')[0]
    outbase = cat_base + '_with_flambda.cat'
    flambda_output_catalog = os.path.join(TEST_DATA_DIR, outbase)
    os.remove(flambda_output_catalog)
    os.remove(sed_catalog)


def test_convert_to_flam():
    """Test basic conversion"""
    magnitudes = np.array([15., 20.])

    photflam = 1e-20
    photfnu = 1e-33
    zeropoint = 25.0
    pivot = 4.44 * u.micron
    params = (photflam, photfnu, zeropoint, pivot)
    mag_sys = 'abmag'
    flam = spec.convert_to_flam(magnitudes, params, mag_sys)
    truth = np.array([5.51426449e-17, 5.51426449e-19])
    assert np.allclose(flam, truth, atol=1e-25)


def test_single_mag_column():
    """Case where input ascii catalog contains only one magntude column.
    In this case extrapolation is necessary."""
    catfile = os.path.join(TEST_DATA_DIR, 'point_sources_one_filter.cat')
    output_hdf5 = os.path.join(TEST_DATA_DIR, 'all_spectra.hdf5')
    sed_catalog = spec.make_all_spectra(catfile, output_filename=output_hdf5)
    comparison = hdf5.open(os.path.join(TEST_DATA_DIR, 'output_spec_from_one_filter.hdf5'))
    constructed = hdf5.open(output_hdf5)
    for key in comparison:
        assert key in constructed.keys()
        assert all(comparison[key]["wavelengths"].value == constructed[key]["wavelengths"].value)
        assert all(comparison[key]["fluxes"].value == constructed[key]["fluxes"].value)
        assert comparison[key]["wavelengths"].unit == constructed[key]["wavelengths"].unit
        assert comparison[key]["fluxes"].unit == constructed[key]["fluxes"].unit

    cat_base = catfile.split('.')[0]
    outbase = cat_base + '_with_flambda.cat'
    flambda_output_catalog = os.path.join(TEST_DATA_DIR, outbase)
    os.remove(flambda_output_catalog)
    os.remove(sed_catalog)


def test_multiple_ascii_catalogs():
    """Case where multiple ascii catalogs are input"""
    catfile = os.path.join(TEST_DATA_DIR, 'point_sources.cat')
    galfile = os.path.join(TEST_DATA_DIR, 'galaxies.cat')
    catalogs = [catfile, galfile]
    output_hdf5 = os.path.join(TEST_DATA_DIR, 'all_spectra.hdf5')
    sed_catalog = spec.make_all_spectra(catalogs, output_filename=output_hdf5)
    comparison = hdf5.open(os.path.join(TEST_DATA_DIR, 'source_sed_file_from_point_sources.hdf5'))
    constructed = hdf5.open(output_hdf5)
    for key in comparison:
        assert key in constructed.keys()
        assert all(comparison[key]["wavelengths"].value == constructed[key]["wavelengths"].value)
        assert all(comparison[key]["fluxes"].value == constructed[key]["fluxes"].value)
        assert comparison[key]["wavelengths"].unit == constructed[key]["wavelengths"].unit
        assert comparison[key]["fluxes"].unit == constructed[key]["fluxes"].unit
    for catfile in catalogs:
        cat_base = catfile.split('.')[0]
        outbase = cat_base + '_with_flambda.cat'
        flambda_output_catalog = os.path.join(TEST_DATA_DIR, outbase)
        os.remove(flambda_output_catalog)
    os.remove(sed_catalog)