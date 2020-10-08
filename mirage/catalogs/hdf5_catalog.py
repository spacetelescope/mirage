#! /usr/bin/env python

"""This module contains code for reading and writing source spectra into
hdf5 files. This format is needed as input to the disperser package
used to create dispersed seed images

Authors
-------

    - Bryan Hilbert

Use
---

    This script is intended to be executed as such:

    ::

        from mirage.catalogs import hdf5_catalog
        spectra_dict = hdf5_catalog.open('my_catalog_file.hdf5')

        hdf5_catalog.save(spectra_dict, 'my_updated_catalog_file.hdf5')
"""

import astropy.units as u
import h5py
import logging
import os

from mirage.logging import logging_functions
from mirage.utils.constants import FLAMBDA_CGS_UNITS, FLAMBDA_MKS_UNITS, FNU_CGS_UNITS, FNU_MKS_UNITS, \
                                   LOG_CONFIG_FILENAME, STANDARD_LOGFILE_NAME

classdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../'))
log_config_file = os.path.join(classdir, 'logging', LOG_CONFIG_FILENAME)
logging_functions.create_logger(log_config_file, STANDARD_LOGFILE_NAME)


def open(filename):
    """Read in contents of an hdf5 file

    Parameters
    ----------
    filename : str
        Name of file to be opened

    Returns
    -------
    contents : dict
        Dictionary containing the contents of the file
        Dictionary format:
        keys are the index numbers of the sources corresponding to the segmentation map
        Each value is a dictionary containing keys 'wavelengths' and 'fluxes'.
        'wavelengths' is an astropy.units Quantity composed of a list of wavelength values
        and a wavelength unit
        'fluxes' is an astropy.units Quantity composed of a list of flux values with flux unit
    """
    logger = logging.getLogger('mirage.catalogs.hdf5_catalog.open')
    contents = {}
    with h5py.File(filename, 'r') as file_obj:
        no_wave_units = False
        no_flux_units = False
        for key in file_obj.keys():
            dataset = file_obj[key]
            try:
                wave_units_string = dataset.attrs['wavelength_units']
            except KeyError:
                wave_units_string = 'micron'
                no_wave_units = True
            try:
                flux_units_string = dataset.attrs['flux_units']
            except KeyError:
                flux_units_string = 'flam_cgs'
                no_flux_units = True

            # Catch common errors
            if wave_units_string.lower() in ['microns', 'angstroms', 'nanometers']:
                wave_units_string = wave_units_string[0:-1]

            # Convert the unit strings into astropy.units Unit object
            wave_units = string_to_units(wave_units_string)
            flux_units = string_to_units(flux_units_string)

            # Get the data
            waves = dataset[0] * wave_units
            fluxes = dataset[1] * flux_units

            # Convert wavelengths to microns and flux values to f_lambda in cgs
            if wave_units != u.micron:
                if wave_units.is_equivalent(u.micron):
                    waves = waves.to(u.micron)
                elif wave_units.is_equivalent(u.Hz):
                    waves = waves.to(u.micron, equivalencies=u.spectral())
                else:
                    raise ValueError("Wavelength units of {} in dataset {} are not compatible with microns."
                                     .format(wave_units, key))
            if flux_units != FLAMBDA_CGS_UNITS:
                if flux_units.is_equivalent(FLAMBDA_CGS_UNITS):
                    fluxes = fluxes.to(FLAMBDA_CGS_UNITS)
                elif flux_units.is_equivalent(FNU_CGS_UNITS):
                    fluxes = fluxes.to(FLAMBDA_CGS_UNITS, u.spectral_density(waves))
                elif flux_units == u.pct:
                    pass
                else:
                    raise ValueError("Flux density units of {} in dataset {} are not compatible with f_lambda."
                                     .format(flux_units, key))

            contents[int(key)] = {'wavelengths': waves, 'fluxes': fluxes}
    if no_wave_units:
        logger.info("{}: No wavelength units provided. Assuming MIRCONS.".format(filename))
    if no_flux_units:
        logger.info("{}: No flux density units provided. Assuming F_lambda in CGS (erg/sec/cm^2/A)".format(filename))
    return contents


def open_tso(filename):
    """Read in an hdf5 file containing TSO data

    Parameters
    ----------
    filename : str
        Name of file to be opened

    Returns
    -------
    contents : dict
        Dictionary containing the contents of the file
        Dictionary format:
        keys are the index numbers of the sources corresponding to the segmentation map
        Each value is a dictionary containing keys 'times' and 'fluxes'.
        'times' is an astropy.units Quantity composed of a list of time values
        and a time unit
        'fluxes' is an astropy.units Quantity composed of a list of flux values with flux unit
    """
    logger = logging.getLogger('mirage.catalogs.hdf5_catalog.open_tso')
    contents = {}
    with h5py.File(filename, 'r') as file_obj:
        no_time_units = False
        no_flux_units = False
        for key in file_obj.keys():
            dataset = file_obj[key]
            try:
                time_units_string = dataset.attrs['time_units']
            except KeyError:
                time_units_string = 'second'
                no_time_units = True
            try:
                flux_units_string = dataset.attrs['flux_units']
            except KeyError:
                flux_units_string = 'normalized'
                no_flux_units = True

            # Catch common errors
            if time_units_string.lower() in ['seconds', 'minutes', 'hours', 'days']:
                time_units_string = time_units_string[0:-1]

            # Convert the unit strings into astropy.units Unit object
            time_units = string_to_units(time_units_string)
            flux_units = string_to_units(flux_units_string)

            # Get the data
            times = dataset[0] * time_units
            fluxes = dataset[1] * flux_units

            # Convert times to seconds
            if time_units != u.second:
                if time_units.is_equivalent(u.second):
                    times = times.to(u.second)
                else:
                    raise ValueError("Time units of {} in dataset {} are not compatible with seconds."
                                     .format(time_units, key))
            if flux_units != u.pct:
                raise ValueError("Flux units are assumed to be normalized. But units in dataset {} are {}."
                                 .format(key, flux_units))

            contents[int(key)] = {'times': times, 'fluxes': fluxes}
    if no_time_units:
        logger.info("{}: No time units provided. Assuming SECONDS.".format(filename))
    if no_flux_units:
        logger.info("{}: No flux density units provided. Assuming light curve is normalized.".format(filename))
    return contents


def save(contents, filename, wavelength_unit='', flux_unit=''):
    """Save a dictionary into an hdf5 file

    Paramters
    ---------
    contents : dict
        Dictionary of data. Dictionary format:
        keys are the index numbers of the sources corresponding to the segmentation map
        Each value is a dictionary containing keys 'wavelengths' and 'fluxes'.
        'wavelengths' is an astropy.units Quantity composed of a list of wavelength values
        and wavelength units
        'fluxes' is an astropy.units Quantity composed of a list of flux values with flux units

    filename : str
        Name of hdf5 file to produce
    """
    with h5py.File(filename, "w") as file_obj:
        for key in contents.keys():
            flux = contents[key]['fluxes']
            wavelength = contents[key]['wavelengths']

            # If units are astropy.units Units objects, change to strings
            # If wavelengths are not a Quantity, fall back onto wavelength_unit
            # and flux_unit
            if isinstance(wavelength, u.quantity.Quantity):
                wavelength_units = units_to_string(wavelength.unit)
                wavelength_values = wavelength.value
            else:
                wavelength_units = wavelength_unit
                wavelength_values = wavelength
            if isinstance(flux, u.quantity.Quantity):
                flux_units = units_to_string(flux.unit)
                flux_values = flux.value
            else:
                flux_units = flux_unit
                flux_values = flux

            dset = file_obj.create_dataset(str(key), data=[wavelength_values, flux_values], dtype='f',
                                           compression="gzip", compression_opts=9)

            # Set dataset units. Not currently inspected by mirage.
            if wavelength_units != '':
                dset.attrs[u'wavelength_units'] = wavelength_units
            if flux_units != '':
                dset.attrs[u'flux_units'] = flux_units


def save_tso(contents, filename, time_unit='', flux_unit=''):
    """Save a dictionary of TSO data into an hdf5 file

    Paramters
    ---------
    contents : dict
        Dictionary of data. Dictionary format:
        keys are the index numbers of the sources corresponding to the segmentation map
        Each value is a dictionary containing keys 'times' and 'fluxes'.
        'times' is an astropy.units Quantity composed of a list oftime values
        and time units
        'fluxes' is an astropy.units Quantity composed of a list of flux values with flux units

    filename : str
        Name of hdf5 file to produce

    time_unit : str
        In the case where the ``times`` entry of ``contents`` is a numpy array
        rather than an astropy Quantity, this unit will be used for the ``times``
        entry. (e.g. 'second')

    flux_unit : str
        In the case where the ``fluxes`` entry of ``contents`` is a numpy array
        rather than an astropy Quantity, this unit will be used for the ``fluxes``
        entry. (e.g. 'flam_cgs')
    """
    with h5py.File(filename, "w") as file_obj:
        for key in contents.keys():
            flux = contents[key]['fluxes']
            time = contents[key]['times']

            # If units are astropy.units Units objects, change to strings
            # If wavelengths are not a Quantity, fall back onto wavelength_unit
            # and flux_unit
            if isinstance(time, u.quantity.Quantity):
                time_units = units_to_string(time.unit)
                time_values = time.value
            else:
                time_units = time_unit
                time_values = time
            if isinstance(flux, u.quantity.Quantity):
                flux_units = units_to_string(flux.unit)
                flux_values = flux.value
            else:
                flux_units = flux_unit
                flux_values = flux

            dset = file_obj.create_dataset(str(key), data=[time_values, flux_values], dtype='f',
                                           compression="gzip", compression_opts=9)

            # Set dataset units. Not currently inspected by mirage.
            if time_units != '':
                dset.attrs[u'time_units'] = time_units
            if flux_units != '':
                dset.attrs[u'flux_units'] = flux_units


def string_to_units(unit_string):
    """Convert a string containing units to an astropy.units Quantity

    Parameters
    ----------
    unit_string : str
        String containing units (e.g. 'erg/sec/cm/cm/A')

    Returns
    -------
    units : astropy.units Quantity
    """
    logger = logging.getLogger('mirage.catalogs.hdf5_catalog.string_to_units')
    if unit_string.lower() in ['flam', "flam_cgs"]:
        return FLAMBDA_CGS_UNITS
    elif unit_string.lower() == 'flam_mks':
        return FLAMBDA_MKS_UNITS
    elif unit_string.lower() in ['fnu', 'fnu_cgs']:
        return FNU_CGS_UNITS
    elif unit_string.lower() == 'fnu_mks':
        return FNU_MKS_UNITS
    elif unit_string in ['normalized', 'NORMALIZED', 'Normalized']:
        return u.pct
    else:
        try:
            return u.Unit(unit_string)
        except ValueError as e:
            logger.error(e)


def units_to_string(unit):
    """Convert the units of an astropy.units Quantity to a string

    Parameters
    ----------
    quantity : astropy.units Quantity

    Returns
    -------
    unit_string : str
        String representation of the units in quantity
    """
    if unit == FLAMBDA_CGS_UNITS:
        return 'flam_cgs'
    elif unit == FLAMBDA_MKS_UNITS:
        return 'flam_mks'
    elif unit == FNU_CGS_UNITS:
        return 'fnu_cgs'
    elif unit == FNU_MKS_UNITS:
        return 'fnu_mks'
    elif unit == u.pct:
        return 'normalized'
    else:
        return unit.to_string()
