#! /usr/bin/env python

"""This module contains functions related to logging in Mirage
"""

import yaml
import logging
from logging.config import dictConfig

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
