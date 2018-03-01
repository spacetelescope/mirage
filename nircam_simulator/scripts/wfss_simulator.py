#! /usr/bin/env python

'''
To make the generation of WFSS simulated integrations easier,
combine the 4 stages of the simulator
(seed image generator, disperser, dark prep, obervation generator)
into a single script.

Inputs:
paramfiles - List of yaml filenames. These files should be 
             inputs to the simulator. For details on the 
             information contained in the yaml files, see
             the readme file associated with the nircam_simulator
             github repo: 
             https://github.com/spacetelescope/nircam_simulator.git

crossing_filter - Name of the crossing filter to be used in
                  conjunction with the grism. All longwave
                  channel filter names are valid entries

module - Name of the NIRCam module to use for the simulation.
         Can be 'A' or 'B'

direction - Dispersion direction. Can be along rows ('R') or 
            along columns ('C')

override_dark - If you wish to use a dark current file that 
                has already gone through the dark_prep step
                of the pipeline and wish to use that for the
                simulation, set override_dark equal to the
                dark's filename. The dark_prep step will then
                be skipped.

HISTORY:
15 November 2017 - created, Bryan Hilbert
'''

import os
import sys
import argparse
from nircam_simulator.scripts import catalog_seed_image
from nircam_simulator.scripts import dark_prep
from nircam_simulator.scripts import obs_generator
from nircam_simulator.scripts import read_fits
from nircam_simulator.scripts import yaml_update
from NIRCAM_Gsim.grism_seed_disperser import Grism_seed

nircam_filters = ['F322W2','F277W','F356W','F444W','F250M'
                  ,'F300M','F335M','F360M','F410M','F430M'
                  ,'F323N','F405N','F466N','F470N']


class WFSSSim():
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

        self.paramfiles = None
        self.override_dark = None
        self.crossing_filter = None
        self.module = None
        self.direction = None
        self.prepDark = None
        self.save_dispersed_seed = True
        self.disp_seed_filename = None
        self.extrapolate_SED = False

    def create(self):
        # Make sure inputs are correct
        self.check_inputs()
        
        # Loop over the yaml files and create
        # a direct seed image for each
        imseeds = []
        # Create imaging seed images
        for pfile in self.paramfiles:
            cat = catalog_seed_image.Catalog_seed()
            cat.paramfile = pfile
            cat.make_seed()
            imseeds.append(cat.seed_file)

        # Create dispersed seed image from
        # the direct images
        dmode = 'mod{}_{}'.format(self.module,self.direction)
        loc = os.path.join(self.datadir,"GRISM_NIRCAM/")
        background_file = ("{}_{}_back.fits"
                           .format(self.crossing_filter,dmode))
        disp_seed = Grism_seed(imseeds, self.crossing_filter,
                               dmode, config_path=loc,
                               extrapolate_SED=self.extrapolate_SED)
        disp_seed.observation()
        disp_seed.finalize(Back = background_file)

        # Save the dispersed seed image
        if self.save_dispersed_seed:
            from astropy.io import fits
            hh00 = fits.PrimaryHDU()
            hh11 = fits.ImageHDU(disp_seed.final)
            hhll = fits.HDUList([hh00,hh11])
            if self.disp_seed_filename is None:
                pdir, pf = os.path.split(self.paramfiles[0])
                dname = 'dispersed_seed_image_for_' + pf + '.fits'
                self.disp_seed_filename = os.path.join(pdir, dname)
            hhll.writeto(self.disp_seed_filename, overwrite=True)
            print(("Dispersed seed image saved to {}"
                   .format(self.disp_seed_filename)))
            
        # Prepare dark current exposure if
        # needed.
        if self.override_dark is None:
            d = dark_prep.DarkPrep()
            d.paramfile = self.paramfiles[0]
            d.prepare()
            obslindark = d.prepDark
        else:
            self.read_dark_product(self.override_dark)
            obslindark = self.prepDark

        # Using the first of the imaging seed image yaml
        # files as a base, adjust to create the yaml file
        # for the creation of the final dispersed
        # integration
        y = yaml_update.YamlUpdate()
        y.file = self.paramfiles[0]
        y.filter = self.crossing_filter
        y.pupil = 'GRISM' + self.direction
        y.outname = ("wfss_dispersed_{}_{}.yaml"
                     .format(dmode,self.crossing_filter))
        y.run()
        
        # Combine into final observation
        obs = obs_generator.Observation()
        obs.linDark = obslindark
        obs.seed = disp_seed.final
        obs.segmap = cat.seed_segmap
        obs.seedheader = cat.seedinfo
        obs.paramfile = y.outname
        obs.create()

        
    def read_dark_product(self,file):
        # Read in dark product that was produced
        # by dark_prep.py
        self.prepDark = read_fits.Read_fits()
        self.prepDark.file = file
        self.prepDark.read_astropy()


    def check_inputs(self):
        # Make sure input parameters are good
        if self.module not in ['A','B']:
            self.invalid('module',self.module)
        else:
            self.module = self.module.upper()

        if self.direction not in ['R','C']:
            self.invalid('direction',self.direction)
        else:
            self.direction = self.direction.upper()

        if self.crossing_filter not in nircam_filters:
            self.invalid('crossing_filter',self.crossing_filter)
        else:
            self.crossing_filter = self.crossing_filter.upper()
        
        if self.override_dark is not None:
            avail = os.path.isfile(self.override_dark)
            if avail == False:
                print(("WARNING: {} does not exist."
                       .format(self.override_dark)))
                sys.exit()

        if len(self.paramfiles) < 2:
            print("WARNING: self.paramfiles must be a list")
            print("of 2 or more yaml files.")
            sys.exit()
        

    def invalid(self,field,value):
        print(("WARNING: invalid value for {}: {}"
               .format(field,value)))
        sys.exit()

        
    def add_options(self,parser = None, usage = None):
        if parser is None:
            parser = argparse.ArgumentParser(usage = usage, description="Wrapper for the creation of WFSS simulated exposures.")
        parser.add_argument("paramfiles",help='List of files describing the input parameters and instrument settings to use. (YAML format).',nargs='+')
        parser.add_argument("--crossing_filter",help = "Name of crossing filter to use in conjunction with the grism.",default=None)
        parser.add_argument("--module",help = "NIRCam module to use for simulation. Use 'A' or 'B'",default=None)
        parser.add_argument("--direction",help = "Direction of dispersion (along rows or along columns). Use 'R' or 'C'",default=None)
        parser.add_argument("--override_dark",help="If supplied, skip the dark preparation step and use the supplied dark to make the exposure", default=None)
        parser.add_argument("--extrapolate_SED", help="If true, the SED created from the filter-averaged magnitudes will be extrapolated to fill the wavelngth range of the grism", action='store_true')
        return parser


if __name__ == '__main__':

    usagestring = 'USAGE: wfss_simualtor.py file1.yaml file2.yaml --crossing_filter F444W --direction R --module A' 

    obs = WFSSSim()
    parser = obs.add_options(usage = usagestring)
    args = parser.parse_args(namespace = obs)
    obs.create()
