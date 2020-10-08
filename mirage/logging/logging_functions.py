#! /usr/bin/env python

"""This module contains functions related to logging in Mirage
"""

import datetime
from functools import wraps
import logging
from logging.config import dictConfig
import os
import shutil
import traceback
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


def create_standard_logfile_name(base_file, log_type='catalog_seed'):
    """Construct the name for the final log file.

    Parameters
    ----------
    base_file : str
        Name of the file on which the name of the log file will be based.
        For a ```log_type``` of "catalog_seed" or "fits_seed" this should
        be one of the yaml input files used by i.e. catalog_seed_image.
        For a ```log_type``` of "yaml_generator" this should be the xml
        file from APT describing the proposal.

    log_type : str
        Type of log file being created. Must be one of:
        'catalog_seed': for log files from catalog_seed_image, dark_prep,
                        obs_generator, imaging_simulator, wfss_simulator,
                        grism_tso_simulator
        'fits_seed': for log files from fits_seed_image
        'yaml_generator': for log files from yaml_generator, apt_reader

    Reutrns
    -------
    log_name : str
        Name of the log file, based on base_file and log type
    """
    log_type = log_type.lower()
    if log_type == 'catalog_seed':
        log_name = base_file.replace('.yaml', '.log')
    elif log_type == 'fits_seed':
        log_name = base_file.replace('.yaml', 'fits_seed_image.log')
    elif log_type == 'yaml_generator':
        time_stamp = datetime.datetime.now().strftime('%d-%b-%YT%H:%m:%S:%f')
        log_name = base_file.replace('.xml', '_{}.log'.format(time_stamp))
    return log_name

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


def log_fail(func):
    """Decorator to log crashes in the decorated code.

    Parameters
    ----------
    func : func
        The function to decorate.

    Returns
    -------
    wrapped : func
        The wrapped function.
    """

    @wraps(func)
    def wrapped(*args, **kwargs):

        try:
            # Run the function
            func(*args, **kwargs)
            #logging.info('Completed Successfully')

        except Exception:
            logging.critical(traceback.format_exc())
            logging.critical('CRASHED')
            raise

    return wrapped


def move_logfile_to_standard_location(base_file, input_log_file, yaml_outdir=None, log_type='catalog_seed'):
    """Copy the log file from the current working directory and standard
    name to a ```mirage_logs``` subdirectory below the simulation data
    output directory. Rename the log to match the input yaml file name
    with a suffix of `log`

    Parameters
    ----------
    base_file : str
        Name of the file on which the name of the log file will be based.
        For a ```log_type``` of "catalog_seed" or "fits_seed" this should
        be one of the yaml input files used by i.e. catalog_seed_image.
        For a ```log_type``` of "yaml_generator" this should be the xml
        file from APT describing the proposal.

    input_log_file : str
        Name of the log file to be copied

    yaml_outdir : str
        Name of the output directory containing the simulated data. This
        is the directory listed in the Output:directory entry of the yaml
        file. It is here as an optional parameter to save having to read
        the yaml file

    log_type : str
        Type of log file being created. Must be one of:
        'catalog_seed': for log files from catalog_seed_image, dark_prep,
                        obs_generator, imaging_simulator, wfss_simulator,
                        grism_tso_simulator
        'fits_seed': for log files from fits_seed_image
        'yaml_generator': for log files from yaml_generator, apt_reader
    """
    # Construct the log filename
    final_logfile_name = create_standard_logfile_name(os.path.basename(base_file), log_type=log_type)

    # Get the directory name in which to place the log file
    if yaml_outdir is None:
        yaml_outdir = get_output_dir(base_file)
    final_logfile_dir = os.path.join(yaml_outdir, 'mirage_logs')
    if not os.path.exists(final_logfile_dir):
        os.makedirs(final_logfile_dir)

    # Copy the log into the directory
    shutil.copy2(input_log_file, os.path.join(final_logfile_dir, final_logfile_name))


def read_yaml(filename):
    """Read the contents of a yaml file into a nested dictionary. This
    here to prevent circular import error

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
            data = yaml.load(f, Loader=yaml.FullLoader)
    except FileNotFoundError as e:
            print(e)
    return data

