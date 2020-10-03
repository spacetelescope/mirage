#! /usr/bin/env python

"""This module contains functions related to logging in Mirage
"""

import logging
from logging.config import dictConfig
import os
import shutil
import yaml


def create_logger(filename, output_log_file):
    """Using the provided yaml file, create a logger

    Parameters
    ----------
    filename : str
        Name of a yaml file containing details on the logger, handler, etc
        to create.

    output_log_file : str
        Name of the text log file to output the log to.
    """
    with open(filename) as fobj:
        log_info = yaml.safe_load(fobj)

    # Set the filename of the output log
    log_info['handlers']['file']['filename'] = output_log_file

    # Create the logger using the dictionary from the yaml file
    dictConfig(log_info)
    logger = logging.getLogger('mirage')


def get_output_dir(yaml_file):
    """Get the output directory name from self.paramfile

    Parameters
    ----------
    yaml_file : str
        Name of yaml file to be inspected

    Returns
    -------
    outdir : str
        Output directory specified in self.paramfile
    """
    params = read_yaml(yaml_file)
    outdir = params['Output']['directory']
    return outdir


def move_logfile_to_standard_location(yaml_file, input_log_file, yaml_outdir=None):
    """Copy the log file from the current working directory and standard
    name to a ```mirage_logs``` subdirectory below the simulation data
    output directory. Rename the log to match the input yaml file name
    with a suffix of `log`

    Parameters
    ----------
    yaml_file : str
        Name of input yaml file used for the simulation

    input_log_file : str
        Name of the log file to be copied

    yaml_outdir : str
        Name of the output directory containing the simulated data. This
        is the directory listed in the Output:directory entry of the yaml
        file. It is here as an optional parameter to save having to read
        the yaml file
    """
    final_logfile_name = yaml_file.replace('.yaml', '.log')
    if yaml_outdir is None:
        yaml_outdir = get_output_dir(yaml_file)
    final_logfile_dir = os.path.join(yaml_outdir, 'mirage_logs')
    if not os.path.exists(final_logfile_dir):
        os.makedirs(final_logfile_dir)
    shutil.copy2(input_log_file, os.path.join(final_logfile_dir, final_logfile_name))


def read_yaml(filename):
    """Read the contents of a yaml file into a nested dictionary

    Parameters
    ----------
    filename : str
        Name of yaml file to be read in

    Returns
    -------
    data : dict
        Nested dictionary of file contents
    """
    try:
        with open(filename, 'r') as f:
            data = yaml.load(f, Loader=yaml.SafeLoader)
    except FileNotFoundError as e:
            print(e)
    return data

