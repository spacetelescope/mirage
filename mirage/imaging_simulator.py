#! /usr/bin/env python

'''
To make the generation of imaging (including moving target)
simulated integrations easier, combine the 3 relevant stages
of the simulator (seed image generator, dark prep,
obervation generator) into a single script.

HISTORY:
13 November 2017 - created, Bryan Hilbert
13 July 2018 - updated for name change to Mirage, Bryan Hilbert
'''

import os
import argparse
import logging
import yaml
import shutil

from .seed_image import catalog_seed_image
from .dark import dark_prep
from .ramp_generator import obs_generator
from .logging import logging_functions
from .utils import read_fits
from .utils.constants import LOG_CONFIG_FILENAME, STANDARD_LOGFILE_NAME
from .utils.utils import ensure_dir_exists, expand_environment_variable


classpath = os.path.dirname(__file__)
log_config_file = os.path.join(classpath, 'logging', LOG_CONFIG_FILENAME)
logging_functions.create_logger(log_config_file, STANDARD_LOGFILE_NAME)


class ImgSim():
    """Class to hold a simulated exposure

    Parameters
    ----------
    paramfile : str
        Name of yaml file to be used as simulator input.
        For details  on the information contained in the
        yaml files, see:
        https://mirage-data-simulator.readthedocs.io/en/latest/example_yaml.html

    override_dark : str or list
        List (or single filename) of outputs from a prior run of ``dark_prep`` to
        use when creating simulation. If set, the call to ``dark_prep`` will be
        skipped and these darks will be used instead. If None, ``dark_prep`` will
        be called and new dark objects will be created.

    """
    def __init__(self, paramfile=None, override_dark=None, offline=False):
        self.env_var = 'MIRAGE_DATA'
        datadir = expand_environment_variable(self.env_var, offline=offline)

        self.paramfile = paramfile
        self.override_dark = override_dark
        self.offline = offline

    def create(self):
        # Initialize the log using dictionary from the yaml file
        self.logger = logging.getLogger('mirage.imaging_simulator')
        self.logger.info('\n\nRunning imaging_simulator....\n')
        self.logger.info('using parameter file: {}'.format(self.paramfile))

        # Create seed image
        cat = catalog_seed_image.Catalog_seed(offline=self.offline)
        cat.paramfile = self.paramfile
        cat.make_seed()

        # Create observation generator object
        obs = obs_generator.Observation(offline=self.offline)

        # Prepare dark current exposure if
        # needed.
        if self.override_dark is None:
            self.logger.info('Perform dark preparation:')
            d = dark_prep.DarkPrep(offline=self.offline)
            d.paramfile = self.paramfile
            d.prepare()

            if len(d.dark_files) == 1:
                obs.linDark = d.prepDark
            else:
                obs.linDark = d.dark_files
        else:
            self.logger.info('\n\noverride_dark has been set. Skipping dark_prep.')
            if isinstance(self.override_dark, str):
                self.read_dark_product(self.override_dark)
                obs.linDark = self.prepDark
            elif isinstance(self.override_dark, list):
                obs.linDark = self.override_dark

        # Combine into final observation
        obs.paramfile = self.paramfile
        if len(cat.seed_files) == 1:
            obs.seed = cat.seedimage
            obs.segmap = cat.seed_segmap
            obs.seedheader = cat.seedinfo
        else:
            obs.seed = cat.seed_files
            self.logger.info('USING SEED_FILES')
            self.logger.info('{}'.format(cat.seed_files))
        obs.create()

        # Make useful information class attributes
        self.seedimage = cat.seedimage
        self.seed_segmap = cat.seed_segmap
        self.seedinfo = cat.seedinfo
        self.linDark = obs.linDark

        self.logger.info('\nImaging simulator complete')
        logging_functions.move_logfile_to_standard_location(self.paramfile, STANDARD_LOGFILE_NAME)

    def get_output_dir(self):
        """Get the output directory name from self.paramfile

        Returns
        -------
        outdir : str
            Output directory specified in self.paramfile
        """
        self.read_param_file()
        outdir = self.params['Output']['directory']
        return outdir

    def read_dark_product(self, file):
        # Read in dark product that was produced
        # by dark_prep.py
        self.prepDark = read_fits.Read_fits()
        self.prepDark.file = file
        self.prepDark.read_astropy()

    def add_options(self, parser=None, usage=None):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage, description="Wrapper for the creation of WFSS simulated exposures.")
        parser.add_argument("paramfile", help='Name of simulator input yaml file')
        parser.add_argument("--override_dark", help="If supplied, skip the dark preparation step and use the supplied dark to make the exposure", default=None)
        return parser


if __name__ == '__main__':

    usagestring = 'USAGE: imaging_simualtor.py file1.yaml'

    obs = ImgSim()
    parser = obs.add_options(usage=usagestring)
    args = parser.parse_args(namespace=obs)
    obs.create()
