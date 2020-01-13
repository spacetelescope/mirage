"""Unit tests for background generation

Authors
-------
    - Bryan Hilbert

Use
---
    Ensure you have pytest installed. Then, simply run pytest in any
    parent directory of mirage/tests/:
    >>> pytest
"""
import astropy.units as u
import datetime
from jwst_backgrounds import jbt
import numpy as np
import os
import pkg_resources
import pysiaf

from mirage.utils import backgrounds, file_io
from mirage.utils.constants import MEAN_GAIN_VALUES

package_path = pkg_resources.resource_filename('mirage', '')
CONFIG_DIR = os.path.join(package_path, 'config')


def test_day_of_year_background_spectrum():
    """Make sure the function returns a spectrum of wavelength and
    flux values.
    """
    ra = 57.2
    dec = -27.6
    observation_date = '2020-12-10'
    waves, signals = backgrounds.day_of_year_background_spectrum(ra, dec, observation_date)

    assert isinstance(waves, np.ndarray)
    assert isinstance(signals, np.ndarray)
    assert len(waves.shape) == 1
    assert len(signals.shape) == 1

    background = jbt.background(ra, dec, 4.)
    obsdate = datetime.datetime.strptime(observation_date, '%Y-%m-%d')
    obs_dayofyear = obsdate.timetuple().tm_yday

    assert obs_dayofyear in background.bkg_data['calendar']


def test_specific_day_of_year_background_spectrum():
    """Test this function using specific inputs and compare to the ETC
    output values"""
    sw_gain = MEAN_GAIN_VALUES['nircam']['swa']
    lw_gain = MEAN_GAIN_VALUES['nircam']['lwa']
    lw_etc = 2.26 / lw_gain  # 2.26 e/s/pix divided by gain 2.19 e/ADU, FOR LWA
    sw_etc = 0.20  / sw_gain  # 0.20 e/s/pix divided by gain 2.44 e/ADU for SWA

    # Use the NIRISS Focus Field
    ra = 85.22458
    dec = -69.5225
    obs_date = '2021-10-04'

    lw_filter_file = os.path.join(CONFIG_DIR, 'F444W_nircam_plus_ote_throughput_moda_sorted.txt')
    #lw_photflam = 7.7190e-22  # FLAM in cgs
    #lw_pivot = 4.3849  # microns
    lw_siaf = pysiaf.Siaf('nircam')['NRCA5_FULL']
    # Here: etc is 1.03, mirage is 0.84. This may be due to a bug in the ETC.

    sw_filter_file = os.path.join(CONFIG_DIR, 'F090W_nircam_plus_ote_throughput_moda_sorted.txt')
    #sw_photflam = 3.3895e-20  # FLAM in cgs
    #sw_pivot = 0.9034  # microns
    sw_siaf = pysiaf.Siaf('nircam')['NRCA2_FULL']

    waves, signals = backgrounds.day_of_year_background_spectrum(ra, dec, obs_date)

    sw_bkgd = backgrounds.calculate_background(ra, dec, sw_filter_file, True, sw_gain, sw_siaf,
                                               back_wave=waves, back_sig=signals)
    lw_bkgd = backgrounds.calculate_background(ra, dec, lw_filter_file, True, lw_gain, lw_siaf,
                                               back_wave=waves, back_sig=signals)

    assert np.isclose(sw_bkgd, sw_etc, atol=0, rtol=0.15)
    assert np.isclose(lw_bkgd, lw_etc, atol=0, rtol=0.25)

def test_find_low_med_high():
    """Test function for finding pre-defined low, medium, and high
    percentile values from an array"""
    array = np.arange(10) + 1
    values = backgrounds.find_low_med_high(array)
    assert values == [1., 5., 9.]


def test_find_percentile():
    """Test the find_percentile function"""
    array = np.arange(10) + 1
    perc = backgrounds.find_percentile(array, 0.6)
    assert perc == 6

    perc = backgrounds.find_percentile(array, 0.2)
    assert perc == 2


def test_low_med_high_background_spectrum():
    """Test the calculation of a background spectrum based on low, medium,
    high values"""
    detector = 'NRCA1'
    module = 'A'
    params = {'Telescope': {'ra': 57.2, 'dec': -27.6},
              'Readout': {'filter': 'F200W', 'pupil': 'CLEAR'},
              'simSignals': {'bkgdrate': 'low'},
              'Inst': {'instrument': 'NIRCAM'},
              'Reffiles': {'flux_cal': os.path.join(CONFIG_DIR, 'NIRCam_zeropoints.list')}
              }

    waves, signals = backgrounds.low_med_high_background_spectrum(params, detector, module)

    assert isinstance(waves, np.ndarray)
    assert isinstance(signals, np.ndarray)
    assert len(waves.shape) == 1
    assert len(signals.shape) == 1


def test_low_medium_high_background_value():
    """Test that the proper integrated background value is calculated for
    a given filter and level"""
    ra = 57.2
    dec = -27.6

    filt_waves = np.array([1., 2., 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 4., 5.])
    filt_thru = np.array([0., 0., 0., 1., 1., 1., 1., 1., 0., 0., 0.])

    siaf = pysiaf.Siaf('nircam')['NRCA1_FULL']
    bkgd_high = backgrounds.low_medium_high_background_value(ra, dec, "high", filt_waves, filt_thru, siaf)
    bkgd_med = backgrounds.low_medium_high_background_value(ra, dec, "medium", filt_waves, filt_thru, siaf)
    bkgd_low = backgrounds.low_medium_high_background_value(ra, dec, "low", filt_waves, filt_thru, siaf)

    assert bkgd_high > bkgd_med
    assert bkgd_med > bkgd_low


def test_specific_low_medium_high_background_value():
    """Test specific cases of this function and compare to ETC outputs
    """
    # etc = {'low': 0.24, 'medium': 0.26, 'high': 0.27}  # MJy/sr from webform
    etc = {'low': 1.246, 'medium': 1.352, 'high': 1.402}  # e-/sec/pixel measured from output image

    # Use the NIRISS Focus Field
    ra = 85.22458
    dec = -69.5225

    siaf = pysiaf.Siaf('niriss')['NIS_CEN']
    filter_throughput_file = os.path.join(CONFIG_DIR, 'f150w_niriss_throughput1.txt')
    filter_waves, filter_thru = file_io.read_filter_throughput(filter_throughput_file)

    bkgd_high = backgrounds.low_medium_high_background_value(ra, dec, "high", filter_waves, filter_thru, siaf)
    bkgd_med = backgrounds.low_medium_high_background_value(ra, dec, "medium", filter_waves, filter_thru, siaf)
    bkgd_low = backgrounds.low_medium_high_background_value(ra, dec, "low", filter_waves, filter_thru, siaf)

    assert np.isclose(bkgd_high, etc['high'], atol=0., rtol=0.05)
    assert np.isclose(bkgd_med, etc['medium'], atol=0., rtol=0.05)
    assert np.isclose(bkgd_low, etc['low'], atol=0., rtol=0.05)


def test_niriss_background_scaling():
    """Test that the scaling factor for NIRISS is correctly calculated"""
    detector = 'NIS'
    module = 'N'
    params = {'Telescope': {'ra': 57.2, 'dec': -27.6},
              'Readout': {'filter': 'F090W', 'pupil': 'CLEAR'},
              'simSignals': {'bkgdrate': 'medium'},
              'Inst': {'instrument': 'NIRISS'},
              'Reffiles': {'flux_cal': os.path.join(CONFIG_DIR, 'niriss_zeropoints.list')}
              }

    # Medium-scaled background
    medium = backgrounds.niriss_background_scaling(params, detector, module)

    # Low-scaled background
    params['simSignals']['bkgdrate'] = 'low'
    low = backgrounds.niriss_background_scaling(params, detector, module)
    assert low < medium
