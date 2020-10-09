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
        dirname = os.path.join(datadir, parameter_dict['INSTRUME'].lower(), 'GRISM_NIRCAM/V2')

        module = parameter_dict['DETECTOR'][3]

        # Assume that detector names in the transmission file names will
        # use 'NRCA5' rather than 'NRCALONG'
        #if 'LONG' in parameter_dict['DETECTOR']:
        #    parameter_dict['DETECTOR'] = parameter_dict['DETECTOR'].replace('LONG', '5')

        if parameter_dict['DETECTOR'] not in ['NRCA5', 'NRCB5', 'NRCALONG', 'NRCBLONG']:
            # For NIRCam SW, we use the same file for all detectors/filters/pupils
            transmission_filename = None
        elif parameter_dict['PUPIL'] not in ['GRISMR', 'GRISMC']:
            # For LW imaging, we use the same file for all detectors/filters/pupils
            transmission_filename = None
        else:
            if module == 'A':
                transmission_filename = os.path.join(dirname, 'NIRCAM_LW_POM_ModA.fits')
            elif module == 'B':
                transmission_filename = os.path.join(dirname, 'NIRCAM_LW_POM_ModB.fits')

    # For NIRISS we search for a detector/filter/pupil-dependent file
    elif parameter_dict['INSTRUME'].lower() == 'niriss':
        dirname = os.path.join(datadir, parameter_dict['INSTRUME'].lower(), 'GRISM_NIRISS/V2')
        filt = parameter_dict['FILTER'].upper()
        if filt not in ['GR150R', 'GR150C']:
            transmission_filename = None
        else:
            transmission_filename = os.path.join(dirname, 'jwst_niriss_cv3_pomtransmission_{}_{}.fits'.format(filt, parameter_dict['PUPIL']))

    # For FGS we don't need to worry about a transmission file
    elif parameter_dict['INSTRUME'].lower() == 'fgs':
        transmission_filename = None

    return transmission_filename
