#! /usr/bin/env python

'''
To make the generation of imaging (including moving target)
simulated integrations easier, combine the 3 relevant stages
of the simulator (seed image generator, dark prep,
obervation generator) into a single script.

Inputs:
paramfile - Name of yaml file to be used as simulator input.
            For details  on the information contained in the
            yaml files, see the readme file associated with
            the mirage github repo:
            https://github.com/spacetelescope/mirage.git

override_dark - If you wish to use a dark current file that
                has already gone through the dark_prep step
                of the pipeline and wish to use that for the
                simulation, set override_dark equal to the
                dark's filename. The dark_prep step will then
                be skipped.

HISTORY:
13 November 2017 - created, Bryan Hilbert
13 July 2018 - updated for name change to Mirage, Bryan Hilbert
'''

import os
import argparse

from .seed_image import catalog_seed_image
from .dark import dark_prep
from .ramp_generator import obs_generator
from .utils import read_fits


class ImgSim():
    def __init__(self, paramfile=None, override_dark=None):
        self.env_var = 'MIRAGE_DATA'
        datadir = os.environ.get(self.env_var)
        if datadir is None:
            raise ValueError(("WARNING: {} environment variable is not set."
                              "This must be set to the base directory"
                              "containing the darks, cosmic ray, PSF, etc"
                              "input files needed for the simulation."
                              "These files must be downloaded separately"
                              "from the Mirage package.".format(self.env_var)))
        self.paramfile = paramfile
        self.override_dark = override_dark

    def create(self):
        # Create seed image
        cat = catalog_seed_image.Catalog_seed()
        cat.paramfile = self.paramfile
        cat.make_seed()

        # Create observation generator object
        obs = obs_generator.Observation()

        # Prepare dark current exposure if
        # needed.
        if self.override_dark is None:
            print('Perform dark preparation:')
            d = dark_prep.DarkPrep()
            d.paramfile = self.paramfile
            d.prepare()
            obs.linDark = d.prepDark
        else:
            self.read_dark_product(self.override_dark)
            obs.linDark = self.prepDark

        # Combine into final observation
        obs.seed = cat.seedimage
        obs.segmap = cat.seed_segmap
        obs.seedheader = cat.seedinfo
        obs.paramfile = self.paramfile
        obs.create()

        # Make useful information class attributes
        self.seedimage = cat.seedimage
        self.seed_segmap = cat.seed_segmap
        self.seedinfo = cat.seedinfo
        self.linDark = obs.linDark

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
