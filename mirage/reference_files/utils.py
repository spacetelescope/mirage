#! /usr/bin/env python

"""
This module contains general functions useful for dealing with Mirage reference files

Author
------

    - Bryan Hilbert


Use
---

    This module can be used as such:
    ::
        from mirage.reference_files import utils
        params = {'INSTRUME': 'NIRCAM', 'DETECTOR': 'NRCA1'}
        reffiles = utils.get_transmission_file(params)
"""
from glob import glob
import os


def get_transmission_file(parameter_dict):
    """Determine the best transmission image file to use for the given
    instrument/detector setup. This image will be located in either the
    GRISM_NIRCAM or GRISM_NIRISS repositories.

    Parameters
    ----------
    parameter_dict : dict
        Dictionary of basic metadata from the file to be processed by the
        returned reference files (e.g. INSTRUME, DETECTOR, etc)

    Returns
    -------
    transmission_filename : str
        Full path to the transmission file
    """
    datadir = os.environ.get('MIRAGE_DATA')

    if parameter_dict['INSTRUME'].lower() == 'nircam':
        dirname = os.path.join(datadir, parameter_dict['INSTRUME'].lower(), 'GRISM_NIRCAM')

        # Assume that detector names in the transmission file names will
        # use 'NRCA5' rather than 'NRCALONG'
        if 'LONG' in parameter_dict['DETECTOR']:
            parameter_dict['DETECTOR'] = parameter_dict['DETECTOR'].replace('LONG', '5')

        if parameter_dict['DETECTOR'] not in ['NRCA5', 'NRCB5', 'NRCALONG', 'NRCBLONG']:
            # For NIRCam SW, we use the same file for all detectors/filters/pupils
            transmission_filename = os.path.join(dirname, 'NIRCam_SW_transmission_image.fits')
        elif parameter_dict['PUPIL'] not in ['GRISMR', 'GRISMC']:
            # For LW imaging, we use the same file for all detectors/filters/pupils
            transmission_filename = os.path.join(dirname, 'NIRCam_LW_imaging_transmission_image.fits')
        else:
            transmission_filename = search_transmission_files(dirname, parameter_dict)

    # For NIRISS we search for a detector/filter/pupil-dependent file
    elif parameter_dict['INSTRUME'].lower() == 'niriss':
        dirname = os.path.join(datadir, parameter_dict['INSTRUME'].lower(), 'GRISM_NIRISS')
        transmission_filename = search_transmission_files(dirname, parameter_dict)

    elif parameter_dict['INSTRUME'].lower() == 'fgs':
        dirname = os.path.join(datadir, 'niriss', 'GRISM_NIRISS')
        transmission_filename = os.path.join(dirname, 'FGS_transmission_image.fits')

    return transmission_file


def search_transmission_files(directory, parameters):
    """Find the transmission file within the collection of transmission files
    that matches the provided detector, filter and pupil values

    Parameters
    ----------
    directory : str
        Name of directory containing the transmission file collection

    parameter_dict : dict
        Dictionary of basic metadata from the file to be processed by the
        returned reference files (e.g. INSTRUME, DETECTOR, etc)

    Returns
    -------
    match[0] : str
        Full path to the transmission file
    """
    files = glob(os.path.join(directory, '*transmission*fits'))

    filt = parameters['FILTER']
    pupil = parameters['PUPIL']
    detector = parameters['DETECTOR']

    # Search using filenames
    if parameters['INSTRUME'].lower() == 'nircam':
        match = [file for file in files if ((filt in file) and (pupil in file) and (detector in file))]
    elif parameters['INSTRUME'].lower() == 'niriss':
        match = [file for file in files if ((filt in file) and (pupil in file))]
    else:
        raise ValueError('Only NIRCam and NIRISS are supported in search_transmission_files()')

    if len(match) == 1:
        return match[0]
    elif len(match) == 0:
        raise ValueError("No matching transmission image files found: {}".format(parameters))
    elif len(match) > 1:
        raise ValueError("More than one matching transmission image file found: {}".format(match))
