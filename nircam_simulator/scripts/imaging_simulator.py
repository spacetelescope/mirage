#! /usr/bin/env python

'''
Combine the three stages of the simulator
(seed image generator, dark prep, obervation generator)
into a single script, for convenience.
'''

import os
import sys
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
                os.environ['NIRCAM_SIM_DATA'] = ('/ifs/jwst/wit/nircam/'
                                                 'nircam_simulator_data')
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

        
    def read_dark_product(self,file):
        # Read in dark product that was produced
        # by dark_prep.py
        self.prepDark = read_fits.Read_fits()
        self.prepDark.file = file
        self.prepDark.read_astropy()

        
