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
    GRISM_NIRCAM or GRISM_NIRISS repositories. For FGS, return a transmission
    image that is the size of the full frame detector and populated with 1.0
    for all pixels.

    Parameters
    ----------
    parameter_dict : dict
        Dictionary of basic metadata from the file to be processed by the
        returned reference files (e.g. INSTRUME, DETECTOR, etc)

    Returns
    -------
    filename : str
        Full path to the transmission file
    """
    if parameter_dict['INSTRUME'].lower() == 'nircam':
        repo = 'GRISM_NIRCAM'
        # Assume that detector names in the transmission file names will
        # use 'NRCA5' rather than 'NRCALONG'
        if 'LONG' in parameter_dict['DETECTOR']:
            parameter_dict['DETECTOR'] = parameter_dict['DETECTOR'].replace('LONG', '5')
    elif parameter_dict['INSTRUME'].lower() == 'niriss':
        repo = 'GRISM_NIRISS'
    elif parameter_dict['INSTRUME'].lower() == 'fgs':
        return 'create_in_place.fits'

    dirname = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../{}'.format(repo)))
    transmission_files = glob(os.path.join(dirname, '*transmission*fits'))

    filt = parameter_dict['FILTER']
    pupil = parameter_dict['PUPIL']
    detector = parameter_dict['DETECTOR']

    # Search using filenames
    if parameter_dict['INSTRUME'].lower() == 'nircam':
        match = [file for file in transmission_files if ((filt in file) and (pupil in file) and (detector in file))]
    elif parameter_dict['INSTRUME'].lower() == 'niriss':
        match = [file for file in transmission_files if ((filt in file) and (pupil in file))]

    if len(match) == 1:
        return match[0]
    elif len(match) == 0:
        raise ValueError("No matching transmission image files found: {}".format(parameter_dict))
    elif len(match) > 1:
        raise ValueError("More than one matching transmission image file found: {}".format(match))
