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
import datetime
from jwst_backgrounds import jbt
import numpy as np
import os
import pkg_resources

from mirage.utils import backgrounds

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

    bkgd_high = backgrounds.low_medium_high_background_value(ra, dec, "high", filt_waves, filt_thru)
    bkgd_med = backgrounds.low_medium_high_background_value(ra, dec, "medium", filt_waves, filt_thru)
    bkgd_low = backgrounds.low_medium_high_background_value(ra, dec, "low", filt_waves, filt_thru)

    assert bkgd_high > bkgd_med
    assert bkgd_med > bkgd_low


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
