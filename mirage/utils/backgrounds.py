#! /usr/bin env python
"""This module contains functions for calculating background signals using
jwst_backgrounds
"""
import datetime
import numpy as np

from jwst_backgrounds import jbt

from mirage.utils.flux_cal import fluxcal_info

# Percentiles corresponding to "low", "medium", and "high" as used in
# the ETC
LOW = 0.1
MEDIUM = 0.5
HIGH = 0.9


def day_of_year_background_spectrum(ra, dec, observation_date):
    """Call jwst_backgrounds in order to produce an estimate of background
    versus wavelength for a particular pointing on a particular day of
    year.

    Parameters
    ----------
    ra : float
        Right Ascention of pointing

    dec : float
        Declination of pointing

    observation_date : str
        Requested day of year of the observation i.e. '2021-10-31'

    Returns
    -------
    background_waves : numpy.ndarray
        1D array with the wavelength values in microns associated with
        ``background_signals``

    background_sigmals : numpy.ndarray
        1D array containing background values in MJy/str
    """
    # Generate background spectra for all days
    background = jbt.background(ra, dec, 4.)

    # Confirm that the target is observable of the requested day of year
    obsdate = datetime.datetime.strptime(observation_date, '%Y-%m-%d')
    obs_dayofyear = obsdate.timetuple().tm_yday
    if obs_dayofyear not in background.bkg_data['calendar']:
        raise ValueError(("ERROR: The requested RA, Dec is not observable on {}. Either "
                          "specify a different day, or set simSignals:use_dateobs_for_background "
                          "to False.".format(observation_date)))

    # Extraxct the spectrum for the requested day
    match = obs_dayofyear == background.bkg_data['calendar']
    background_waves = background.bkg_data['wave_array']
    background_signals = background.bkg_data['total_bg'][match, :][0]

    return background_waves, background_signals


def find_low_med_high(array):
    """Given an array of values, find the value corresponding to the
    Nth percentile.

    Parameters
    ----------
    array : numpy.ndarray
        1D array of values

    Returns
    -------
    levels : list
        [low, medium, high] percentile values of array
    """
    levels = []
    for value in [LOW, MEDIUM, HIGH]:
        levels.append(find_percentile(array, value))
    return levels


def find_percentile(array, percentile):
    """Find the value corresponding to the ``percentile``
    percentile of ``array``.

    Parameters
    ----------
    array : numpy.ndarray
        Array of values to be searched

    percentile : float
        Percentile to search for. For example, to find the 50th percentile
        value, use 0.5
    """
    x = np.sort(array)
    y = np.arange(1, len(x) + 1) / len(x)
    value = np.interp(percentile, y, x)
    return value


def low_med_high_background_spectrum(param_dict, detector, module):
    """Call jwst_backgrounds in order to produce an estimate of background
    versus wavelength for a particular pointing that corresponds to one of
    "low", "medium", or "high", as with the ETC.

    Parameters
    ----------
    param_dict : dict
        Dictionary of observation information from an input yaml file

    detector : str
        Name of detector, e.g. "NRCA1"

    module : str
        Name of module, e.g. "A"

    Returns
    -------
    background_waves : numpy.ndarray
        1D array with the wavelength values in microns associated with
        ``background_signals``

    background_spec : numpy.ndarray
        1D array containing background values in MJy/str

    """
    # Generate background spectra for all days
    background = jbt.background(param_dict['Telescope']['ra'], param_dict['Telescope']['dec'], 4.)

    # Determine which optical element wheel contains the filter we want
    # to use for background estimation
    if param_dict['Readout']['pupil'][0].upper() == 'F':
        usefilt = 'pupil'
    else:
        usefilt = 'filter'

    # Get basic flux calibration information
    vegazp, photflam, photfnu, pivot_wavelength = fluxcal_info(param_dict, usefilt, detector, module)

    # Extract the spectrum value across all days at the pivot wavelength
    wave_diff = np.abs(background.bkg_data['wave_array'] - pivot_wavelength)
    bkgd_wave = np.where(wave_diff == np.min(wave_diff))[0][0]
    bkgd_at_pivot = background.bkg_data['total_bg'][:, bkgd_wave]

    # Now sort and determine the low/medium/high levels
    low, medium, high = find_low_med_high(bkgd_at_pivot)

    # Find the value based on the level in the yaml file
    background_level = param_dict['simSignals']['bkgdrate'].lower()
    if background_level == "low":
        level_value = low
    elif background_level == "medium":
        level_value = medium
    elif background_level == "high":
        level_value = high
    else:
        raise ValueError(("ERROR: Unrecognized background value: {}. Must be low, mediumn, or high"
                          .format(param_dict['simSignals']['bkgdrate'])))

    # Find the day with the background at the pivot wavelength that
    # is closest to the value associated with the requested level
    diff = np.abs(bkgd_at_pivot - level_value)
    mindiff = np.where(diff == np.min(diff))[0][0]
    background_spec = background.bkg_data['total_bg'][mindiff, :]

    return background.bkg_data['wave_array'], background_spec


def low_medium_high_background_value(ra, dec, background_level, filter_waves, filter_throughput):
    """Calculate the integrated background flux density for a given filter,
    using the filter's throughput curve and the user-input background level
    (e.g. "medium")

    Parameters
    ----------
    ra : float
        Right ascention of the pointing. Units are degrees

    dec : float
        Declineation of the pointing. Units are degrees

    background_level : str
        "low", "medium", or "high", just as with the ETC

    filter_waves : numpy.ndarray
        1d array of wavelengths in microns to be used along with
        ``filter_throughput``

    filter_throughput : numpy.ndarray
        1d array of filter throughput values to convolve with the background
        spectrum. Normalized units. 1.0 = 100% transmission.

    Returns
    -------
    value : float
        Background value corresponding to ``background_level``, integrated
        over the filter bandpass. Background units are MJy/str.
    """
    # Get background information
    bg = jbt.background(ra, dec, 4.)
    back_wave = bg.bkg_data['wave_array']
    bsigs = np.zeros(len(bg.bkg_data['total_bg'][:, 0]))
    for i in range(len(bg.bkg_data['total_bg'][:, 0])):
        back_sig = bg.bkg_data['total_bg'][i, :]

        # Interpolate background to match filter wavelength grid
        bkgd_interp = np.interp(filter_waves, back_wave, back_sig)

        # Combine
        filt_bkgd = bkgd_interp * filter_throughput

        # Integrate
        bsigs[i] = np.trapz(filt_bkgd, x=filter_waves)

    # Now sort and determine the low/medium/high levels
    low, medium, high = find_low_med_high(bsigs)

    # Find the value based on the level in the yaml file
    background_level = background_level.lower()
    if background_level == "low":
        value = low
    elif background_level == "medium":
        value = medium
    elif background_level == "high":
        value = high
    else:
        raise ValueError(("ERROR: Unrecognized background value: {}. Must be low, mediumn, or high"
                          .format(background_level)))
    return value


def nircam_background_spectrum(parameters, detector, module):
    """Generate a background spectrum by calling jwst_backgrounds and
    returning wavelengths and flux density values based on observation
    date or low/medium/high

    Parameters
    ----------
    parameters : dict
        Nested dictionary containing parameters pertaining to background
        generation. Designed around dictionary from reading in a Mirage
        input yaml file

    detector : str
        Detector name (e.g. 'NRCA1')

    module : str
        Module name (e.g. 'A')

    Returns
    -------
    waves : numpy.ndarray
        1D array of wavelength values (in microns)

    fluxes : numpy.ndarray
        1D array of flux density values
    """
    if parameters['simSignals']['use_dateobs_for_background']:
        print("Generating background spectrum for observation date: {}"
              .format(parameters['Output']['date_obs']))
        waves, fluxes = day_of_year_background_spectrum(parameters['Telescope']['ra'],
                                                        parameters['Telescope']['dec'],
                                                        parameters['Output']['date_obs'])
    else:
        if isinstance(parameters['simSignals']['bkgdrate'], str):
            if parameters['simSignals']['bkgdrate'].lower() in ['low', 'medium', 'high']:
                print("Generating background spectrum based on requested level of: {}"
                      .format(parameters['simSignals']['bkgdrate']))
                waves, fluxes = low_med_high_background_spectrum(parameters, detector, module)
            else:
                raise ValueError("ERROR: Unrecognized background rate. Must be one of 'low', 'medium', 'high'")
        else:
            raise ValueError(("ERROR: WFSS background rates must be one of 'low', 'medium', 'high', "
                              "or use_dateobs_for_background must be True "))
    return waves, fluxes


def niriss_background_scaling(param_dict, detector, module):
    """Determine the scaling factor needed to translate the pre-existing
    NIRISS WFSS background image to the requested signal level, which is
    one of "low", "medium", or "high", as with the ETC.

    Parameters
    ----------
    param_dict : dict
        Dictionary of observation information from an input yaml file

    detector : str
        Name of detector, e.g. "NRCA1"

    module : str
        Name of module, e.g. "A"

    Returns
    -------
    ratio : float
        Ratio of the background value at the pivot wavelength that
        corresponds to the requested level, ratioed to that for the
        "medium" level, which is what was used to create the NIRISS
        WFSS background images.
    """
    # Generate background spectra for all days
    background = jbt.background(param_dict['Telescope']['ra'], param_dict['Telescope']['dec'], 4.)

    # Determine which optical element wheel contains the filter we want
    # to use for background estimation
    if param_dict['Readout']['pupil'][0].upper() == 'F':
        usefilt = 'pupil'
    else:
        usefilt = 'filter'

    # Get basic flux calibration information
    vegazp, photflam, photfnu, pivot_wavelength = fluxcal_info(param_dict, usefilt, detector, module)

    # Extract the spectrum value across all days at the pivot wavelength
    wave_diff = np.abs(background.bkg_data['wave_array'] - pivot_wavelength)
    bkgd_wave = np.where(wave_diff == np.min(wave_diff))[0][0]
    bkgd_at_pivot = background.bkg_data['total_bg'][:, bkgd_wave]

    # Now sort and determine the low/medium/high levels
    low, medium, high = find_low_med_high(bkgd_at_pivot)

    # Find the value based on the level in the yaml file
    background_level = param_dict['simSignals']['bkgdrate'].lower()
    if background_level == "low":
        level_value = low
    elif background_level == "medium":
        level_value = medium
    elif background_level == "high":
        level_value = high
    else:
        raise ValueError(("ERROR: Unrecognized background value: {}. Must be low, medium, or high"
                          .format(params_dict['simSignals']['bkgdrate'])))
    return level_value
