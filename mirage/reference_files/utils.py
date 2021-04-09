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
import logging
import os

import grismconf


def find_wfss_config_filename(pathname, instrument, filtername, mode):
    """Construct the name of the WFSS configuration file given instrument parameters

    Parameters
    ----------
    pathname : str
        Path where the configuration files are located

    instrument : str
        Instrument name (e.g. 'nircam')

    filtername : str
        Name of the crossing filter (e.g. 'f444w')

    mode : str
        String containing dispersion direction ('R' or 'C') as well as
        module name (for NIRCam). e.g. for NIRCam - modA_R and
        for NIRISS GR150R

    Returns
    -------
    config : str
        Full path and filename of the appropriate configuation file
    """
    config = os.path.join(pathname,"{}_{}_{}.conf".format(instrument.upper(),
                                                          filtername.upper(),
                                                          mode))
    return config


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
    logger = logging.getLogger('mirage.reference_files.utils.get_transmission_file')
    datadir = os.environ.get('MIRAGE_DATA')

    if parameter_dict['INSTRUME'].lower() == 'nircam':
        if parameter_dict['DETECTOR'] not in ['NRCA5', 'NRCB5', 'NRCALONG', 'NRCBLONG']:
            # For NIRCam SW, we don't use a transmission file
            transmission_filename = None
        elif parameter_dict['PUPIL'] not in ['GRISMR', 'GRISMC']:
            # For LW imaging, also don't use a transmission file
            transmission_filename = None
        else:
            module = parameter_dict['DETECTOR'][3]
            instrument = parameter_dict['INSTRUME'].lower()

            # POM transmission file is not grism orientation dependent. The file is the
            # same for GRISMR as GRISMC, so we can just use R here
            dmode = 'mod{}_R'.format(module)
            filt = parameter_dict['FILTER'].upper()

            loc = os.path.join(datadir, "{}/GRISM_{}/current".format(instrument.lower(),
                                                                     instrument.upper()))
            configuration_file = find_wfss_config_filename(loc, 'nircam', filt, dmode)
            c = grismconf.Config(configuration_file)
            transmission_filename = c.POM

    # For NIRISS we search for a detector/filter/pupil-dependent file
    elif parameter_dict['INSTRUME'].lower() == 'niriss':
        if parameter_dict['FILTER'] not in ['GR150R', 'GR150C']:
            # Imaging mode - no transmission file necessary
            transmission_filename = None
        else:
            # POM transmission file is not grism orientation dependent. The file is the
            # same for GR150R as GR150C, so we can just use R here
            instrument = parameter_dict['INSTRUME'].lower()
            dmode = 'GR150R'
            filt = parameter_dict['PUPIL'].upper()

            loc = os.path.join(datadir, "{}/GRISM_{}/current".format(instrument.lower(),
                                                                     instrument.upper()))

            configuration_file = find_wfss_config_filename(loc, 'niriss', filt, dmode)
            c = grismconf.Config(configuration_file)
            transmission_filename = c.POM

    # For FGS we don't need to worry about a transmission file
    elif parameter_dict['INSTRUME'].lower() == 'fgs':
        transmission_filename = None

    logger.info('POM Transmission filename: {}'.format(transmission_filename))
    return transmission_filename
