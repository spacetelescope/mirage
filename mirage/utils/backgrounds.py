

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
    background_file : str
        Name of file into which the background spectrum was saved
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

    # Output file to save background spectrum into
    background_file = 'background.txt'
    background = jbt.get_background(ra, dec, 4., plot_background=False, plot_bathtub=False,
                                    thisday=obs_dayofyear, write_background=True, write_bathtub=False,
                                    background_file=background_file)
    return background_file


def find_low_med_high(array):
    """Given an array of values, find the value corresponding to the
    Nth percentile.
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
    background_file : str
        Name of file into which the background spectrum was saved
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

    # Will we need to save the background to a file because NIRCam_Gsim
    # requires it?
    save_to_ascii_file(background.bkg_data['wave_array'], background_spec)


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
        raise ValueError(("ERROR: Unrecognized background value: {}. Must be low, mediumn, or high"
                          .format(params_dict['simSignals']['bkgdrate'])))

    # The pre-computed background images for NIRISS were created using the
    # "medium" setting. So the scaling factor is the ratio of the requested
    # level to medium
    ratio = level_value / medium
    return ratio
