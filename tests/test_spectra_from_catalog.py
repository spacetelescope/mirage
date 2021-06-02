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
import copy

from astropy.io import ascii
from astropy.table import Table
import astropy.units as u
from synphot import SourceSpectrum, Observation, units
from synphot.spectrum import SpectralElement
from synphot.models import Empirical1D
from synphot.models import BlackBodyNorm1D

from mirage.catalogs import spectra_from_catalog as spec
from mirage.catalogs import hdf5_catalog as hdf5
from mirage.utils.constants import FLAMBDA_CGS_UNITS, FNU_CGS_UNITS, MEAN_GAIN_VALUES
from mirage.utils.utils import get_filter_throughput_file, magnitude_to_countrate


TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), 'test_data/hdf5_catalogs')


def test_get_filter_info():
    nrc = spec.get_filter_info(['nircam_f444w_magnitude'], 'abmag')
    assert nrc == {'nircam_f444w_magnitude': (5.584676852068291e-22 * FLAMBDA_CGS_UNITS, 3.641228911284274e-31 * FNU_CGS_UNITS, 27.496928286774455,
                                              4.421150698281195 * u.micron)}

    nis = spec.get_filter_info(['niriss_f200w_magnitude'], 'vegamag')
    assert nis == {'niriss_f200w_magnitude': (2.173398e-21 * FLAMBDA_CGS_UNITS, 2.879494e-31 * FNU_CGS_UNITS,
                                              26.04898, 1.9930 * u.micron)}

    fgs = spec.get_filter_info(['fgs_magnitude'], 'stmag')
    assert fgs == {'fgs_magnitude': (1.593395e-22 * FLAMBDA_CGS_UNITS, 3.324736e-32 * FNU_CGS_UNITS, 33.39426,
                                     2.5011 * u.micron)}


def test_hdf5_file_input():
    """Case where an hdf5 file is input. One of these also includes a
    normalized spectrum"""
    catfile = os.path.join(TEST_DATA_DIR, 'point_sources.cat')
    output_hdf5 = os.path.join(TEST_DATA_DIR, 'all_spectra.hdf5')
    sed_file = os.path.join(TEST_DATA_DIR, 'sed_file_with_normalized_dataset.hdf5')
    sed_catalog = spec.make_all_spectra(catfile, input_spectra_file=sed_file,
                                        normalizing_mag_column='nircam_f444w_magnitude',
                                        output_filename=output_hdf5, module='A')

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
                                        output_filename=output_hdf5, module='A')
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
                                        output_filename=output_hdf5, module='A')
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
    flam = spec.convert_to_flam('nircam', 'F000W', magnitudes, params, mag_sys)
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


def test_spectra_rescaling():
    """Test the functionality for rescaling input spectra to a given
    magnitude in a given filter
    """
    # JWST primary mirror area in cm^2. Needed for countrate check
    # at the end
    primary_area = 25.326 * (u.m * u.m)

    # Create spectrum: one source to be normalized
    # and the other should not be
    waves = np.arange(0.4, 5.6, 0.01)
    flux = np.repeat(1e-16, len(waves))
    flux2 = np.repeat(4.24242424242e-18, len(waves))  # arbitrary value
    spectra = {1: {"wavelengths": waves * u.micron,
                   "fluxes": flux * u.pct},
               2: {"wavelengths": waves * u.micron,
                   "fluxes": flux2 * units.FLAM}}

    # Create source catalog containing scaling info
    catalog = Table()
    catalog['index'] = [1, 2]
    catalog['nircam_f322w2_magnitude'] = [18.] * 2
    catalog['niriss_f090w_magnitude'] = [18.] * 2
    #catalog['fgs_magnitude'] = [18.] * 2

    # Instrument info
    instrument = ['nircam', 'niriss']  # , 'fgs']
    filter_name = ['F322W2', 'F090W']  # , 'N/A']
    pupil_name = ['CLEAR', 'none']
    module = ['B', 'N']  # , 'F']
    detector = ['NRCB5', 'NIS']  # , 'GUIDER1']

    # Magnitude systems of renormalization magnitudes
    mag_sys = ['vegamag', 'abmag', 'stmag']

    # Loop over options and test each
    for inst, filt, pup, mod, det in zip(instrument, filter_name, pupil_name, module, detector):

        # Extract the appropriate column from the catalog information
        magcol = [col for col in catalog.colnames if inst in col]
        sub_catalog = catalog['index', magcol[0]]

        # Filter throughput files
        filter_thru_file = get_filter_throughput_file(instrument=inst, filter_name=filt,
                                                      pupil_name=pup, nircam_module=mod,
                                                      fgs_detector=det)

        # Retrieve the correct gain value that goes with the fluxcal info
        if inst == 'nircam':
            gain = MEAN_GAIN_VALUES['nircam']['lwb']
        elif inst == 'niriss':
            gain = MEAN_GAIN_VALUES['niriss']
        elif inst == 'fgs':
            gain = MEAN_GAIN_VALUES['fgs'][det.lower()]

        # Create filter bandpass object, to be used in the final
        # comparison
        filter_tp = ascii.read(filter_thru_file)
        bp_waves = filter_tp['Wavelength_microns'].data * u.micron
        thru = filter_tp['Throughput'].data
        bandpass = SpectralElement(Empirical1D, points=bp_waves, lookup_table=thru) / gain

        # Check the renormalization in all photometric systems
        for magsys in mag_sys:
            rescaled_spectra = spec.rescale_normalized_spectra(spectra, sub_catalog, magsys, filter_thru_file, gain)

            # Calculate the countrate associated with the renormalized
            # spectrum through the requested filter
            for dataset in rescaled_spectra:
                if dataset == 1:
                    # This block is for the spectra that are rescaled
                    rescaled_spectrum = SourceSpectrum(Empirical1D, points=rescaled_spectra[dataset]['wavelengths'],
                                                       lookup_table=rescaled_spectra[dataset]['fluxes'])

                    obs = Observation(rescaled_spectrum, bandpass, binset=bandpass.waveset)
                    renorm_counts = obs.countrate(area=primary_area)

                    # Calculate the countrate associated with an object of
                    # matching magnitude
                    if inst != 'fgs':
                        mag_col = '{}_{}_magnitude'.format(inst.lower(), filt.lower())
                    else:
                        mag_col = 'fgs_magnitude'
                    filt_info = spec.get_filter_info([mag_col], magsys)
                    magnitude = catalog[mag_col][dataset - 1]
                    photflam, photfnu, zeropoint, pivot = filt_info[mag_col]


                    print(photflam, photfnu, zeropoint, pivot)

                    check_counts = magnitude_to_countrate(inst, filt, magsys, magnitude, photfnu=photfnu.value,
                                                          photflam=photflam.value, vegamag_zeropoint=zeropoint)

                    if magsys != 'vegamag':
                        # As long as the correct gain is used, AB mag and ST mag
                        # count rates agree very well
                        tolerance = 0.0005
                    else:
                        # Vegamag count rates for NIRISS have a sligtly larger
                        # disagreement. Zeropoints were derived assuming Vega = 0.02
                        # magnitudes. This offset has been added to the rescaling
                        # function, but may not be exact.
                        tolerance = 0.0055

                    # This dataset has been rescaled, so check that the
                    # countrate from the rescaled spectrum matches that from
                    # the magnitude it was rescaled to
                    if isinstance(check_counts, u.quantity.Quantity):
                        check_counts = check_counts.value
                    if isinstance(renorm_counts, u.quantity.Quantity):
                        renorm_counts = renorm_counts.value

                    print(inst, filt, magsys, renorm_counts, check_counts, renorm_counts / check_counts)
                    assert np.isclose(renorm_counts, check_counts, atol=0, rtol=tolerance), \
                        print('Failed assertion: ', inst, filt, magsys, renorm_counts, check_counts,
                              renorm_counts / check_counts)
                elif dataset == 2:
                    # Not rescaled. In this case Mirage ignores the magnitude
                    # value in the catalog, so we can't check against check_counts.
                    # Just make sure that the rescaling function did not
                    # change the spectrum at all
                    assert np.all(spectra[dataset]['fluxes'] == rescaled_spectra[dataset]['fluxes'])
