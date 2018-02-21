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
            the nircam_simulator github repo:
            https://github.com/spacetelescope/nircam_simulator.git

override_dark - If you wish to use a dark current file that
                has already gone through the dark_prep step
                of the pipeline and wish to use that for the
                simulation, set override_dark equal to the
                dark's filename. The dark_prep step will then
                be skipped.

HISTORY:
13 November 2017 - created, Bryan Hilbert
'''

import os
import sys
import argparse
from nircam_simulator.scripts import catalog_seed_image
from nircam_simulator.scripts import dark_prep
from nircam_simulator.scripts import obs_generator
from nircam_simulator.scripts import read_fits

class ImgSim():
    def __init__(self):
        # Set the NIRCAM_SIM_DATA environment variable if it is not
        # set already. This is for users at STScI.
        self.datadir = os.environ.get('NIRCAM_SIM_DATA')
        if self.datadir is None:
            local = os.path.exists('/ifs/jwst/wit/nircam/nircam_simulator_data')
            if local:
                ddir = '/ifs/jwst/wit/nircam/nircam_simulator_data'
                os.environ['NIRCAM_SIM_DATA'] = ddir
                self.datadir = ddir
            else:
                print("WARNING: NIRCAM_SIM_DATA environment")
                print("variable is not set, and it appears that")
                print("you do not have access to the STScI")
                print("directory where the data are located.")
                sys.exit()

        self.paramfile = None
        self.override_dark = None

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
