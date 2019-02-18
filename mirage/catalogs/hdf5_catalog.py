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
    contents = {}
    with h5py.File(filename, 'r') as file_obj:
        for key in file_obj.keys():
            dataset = file_obj[key]
            wave_units_string = dataset.attrs['wavelength_units']
            flux_units_string = dataset.attrs['flux_units']

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
                else:
                    raise ValueError("Wavelength units of {} in dataset {} are not compatible with microns."
                                     .format(wave_units, key))
            f_lambda_cgs = u.erg / u.second / u.cm / u.cm / u.angstrom
            if flux_units != f_lambda_cgs:
                if flux_units.is_equivalent(f_lambda_cgs):
                    fluxes = fluxes.to(f_lambda_cgs)
                else:
                    raise ValueError("Flux units of {} in dataset {} are not compatible with f_lambda."
                                     .format(flux_units, key))

            contents[int(key)] = {'wavelengths': waves, 'fluxes': fluxes}
    return contents


def save(contents, filename, wavelength_unit=None, flux_unit=None):
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
            else:
                wavelength_units = wavelength_unit
            if isinstance(flux, u.quantity.Quantity):
                flux_units = units_to_string(flux.unit)
            else:
                flux_units = flux_unit

            dset = file_obj.create_dataset(str(key), data=[wavelength, flux], dtype='f', compression="gzip",
                                           compression_opts=9)

            # Set dataset units. Not currently inspected by mirage.
            dset.attrs[u'wavelength_units'] = wavelength_units
            dset.attrs[u'flux_units'] = flux_units


def string_to_units(unit_string):
    """Convert a string containing units to an astropy.units Quantity

    Parameters
    ----------
    unit_string : str
        String containing units (e.g. 'erg/sec/cm/cm/A/A')

    Returns
    -------
    units : astropy.units Quantity
    """
    try:
        return u.Unit(unit_string)
    except ValueError as e:
        print(e)


def units_to_string(quantity):
    """Convert the units of an astropy.units Quantity to a string

    Parameters
    ----------
    quantity : astropy.units Quantity

    Returns
    -------
    unit_string : str
        String representation of the units in quantity
    """
    return quantity.unit.to_string()
