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
             the readme file associated with the mirage
             github repo:
             https://github.com/spacetelescope/mirage.git

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
13 July 2018 - updated for name change to Mirage, Bryan Hilbert
'''

import os
import sys
import argparse

from numpy import nanmedian, isfinite
from astropy.io import fits
from NIRCAM_Gsim.grism_seed_disperser import Grism_seed

from .seed_image import catalog_seed_image
from .dark import dark_prep
from .ramp_generator import obs_generator
from .utils import read_fits
from .utils.utils import expand_environment_variable
from .yaml import yaml_update

nircam_filters = ['F322W2', 'F277W', 'F356W', 'F444W', 'F250M', 'F300M',
                  'F335M', 'F360M', 'F410M', 'F430M', 'F323N', 'F405N',
                  'F466N', 'F470N']


class WFSSSim():
    def __init__(self, offline=False):
        # Set the MIRAGE_DATA environment variable if it is not
        # set already. This is for users at STScI.
        self.env_var = 'MIRAGE_DATA'
        self.datadir = expand_environment_variable(self.env_var, offline=offline)

        self.paramfiles = None
        self.override_dark = None
        self.crossing_filter = None
        self.module = None
        self.direction = None
        self.prepDark = None
        self.save_dispersed_seed = True
        self.disp_seed_filename = None
        self.extrapolate_SED = False
        self.fullframe_apertures = ["NRCA5_FULL", "NRCB5_FULL"]
        self.offline = offline

    def create(self):
        # Make sure inputs are correct
        self.check_inputs()

        # Loop over the yaml files and create
        # a direct seed image for each
        imseeds = []
        # Create imaging seed images
        for pfile in self.paramfiles:
            cat = catalog_seed_image.Catalog_seed(offline=self.offline)
            cat.paramfile = pfile
            cat.make_seed()
            imseeds.append(cat.seed_file)

        # Create dispersed seed image from
        # the direct images
        dmode = 'mod{}_{}'.format(self.module,self.direction)
        loc = os.path.join(self.datadir,"nircam/GRISM_NIRCAM/")
        background_file = ("{}_{}_back.fits"
                           .format(self.crossing_filter,dmode))
        disp_seed = Grism_seed(imseeds, self.crossing_filter,
                               dmode, config_path=loc,
                               extrapolate_SED=self.extrapolate_SED)
        disp_seed.observation()
        disp_seed.finalize(Back = background_file)

        # Get gain map
        gainfile = cat.params['Reffiles']['gain']
        gain, gainheader = self.read_gain_file(gainfile)

        # Disperser output is always full frame. Crop to the
        # requested subarray if necessary
        if cat.params['Readout']['array_name'] not in self.fullframe_apertures:
            print("Subarray bounds: {}".format(cat.subarray_bounds))
            print("Dispersed seed image size: {}".format(disp_seed.final.shape))
            disp_seed.final = self.crop_to_subarray(disp_seed.final, cat.subarray_bounds)
            gain = self.crop_to_subarray(gain, cat.subarray_bounds)
            # Segmentation map will be centered in a frame that is larger
            # than full frame by a factor of sqrt(2), so crop appropriately
            segy, segx = cat.seed_segmap.shape
            dx = int((segx - 2048) / 2)
            dy = int((segy - 2048) / 2)
            segbounds = [cat.subarray_bounds[0] + dx, cat.subarray_bounds[1] + dy,
                         cat.subarray_bounds[2] + dx, cat.subarray_bounds[3] + dy]
            cat.seed_segmap = self.crop_to_subarray(cat.seed_segmap, segbounds)

        # Convert seed image to ADU/sec to be consistent
        # with other simulator outputs
        disp_seed.final /= gain

        # Update seed image header to reflect the
        # division by the gain
        cat.seedinfo['units'] = 'ADU/sec'

        # Save the dispersed seed image
        if self.save_dispersed_seed:
            hh00 = fits.PrimaryHDU()
            hh11 = fits.ImageHDU(disp_seed.final)
            hhll = fits.HDUList([hh00,hh11])
            hhll[0].header['units'] = 'ADU/sec'
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
            d = dark_prep.DarkPrep(offline=self.offline)
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
        obs = obs_generator.Observation(offline=self.offline)
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

    def read_param_file(self, file):
        """
        Read in yaml simulator parameter file

        Parameters:
        -----------
        file -- Name of a yaml file in the proper format
                for mirage

        Returns:
        --------
        Nested dictionary with the yaml file's contents
        """
        import yaml
        try:
            with open(file, 'r') as infile:
                data = yaml.load(infile)
        except:
            raise IOError("WARNING: unable to open {}".format(file))

    def read_gain_file(self, file):
        """
        Read in CRDS-formatted gain reference file

        Paramters:
        ----------
        file -- Name of gain reference file

        Returns:
        --------
        Detector gain map (2d numpy array)
        """
        try:
            with fits.open(file) as h:
                image = h[1].data
                header = h[0].header
        except:
            raise IOError("WARNING: Unable to open gain file: {}".format(file))

        mngain = nanmedian(image)

        # Set pixels with a gain value of 0 equal to mean
        image[image == 0] = mngain
        # Set any pixels with non-finite values equal to mean
        image[~isfinite(image)] = mngain
        return image, header

    def crop_to_subarray(self, data, bounds):
        """
        Crop the given full frame array down to the appropriate
        subarray size and location based on the requested subarray
        name.

        Parameters:
        -----------
        data -- 2d numpy array. Full frame image. (2048 x 2048)
        bounds -- 4-element list containing the full frame indices that
                  define the position of the subarray.
                  [xstart, ystart, xend, yend]

        Returns:
        --------
        Cropped 2d numpy array
        """
        yl, xl = data.shape
        valid = [False, False, False, False]
        valid = [(b>=0 and b<xl) for b in bounds[0:3:2]]
        validy = [(b>=0 and b<yl) for b in bounds[1:4:2]]
        valid.extend(validy)

        print(valid)
        print(bounds)
        print(yl, xl)

        if all(valid):
            return data[bounds[1]:bounds[3] + 1, bounds[0]:bounds[2] + 1]
        else:
            raise ValueError(("WARNING: subarray bounds are outside the "
                              "dimensions of the input array."))

    def read_subarr_defs(self, subfile):
        # read in the file that contains a list of subarray
        # names and positions on the detector
        try:
            subdict = ascii.read(subfile, data_start=1, header_start=0)
            return subdict
        except:
            raise RuntimeError(("Error: could not read in subarray definitions file "
                                "{}".format(subfile)))

    def get_subarr_bounds(self, subname, sdict):
        # find the bounds of the requested subarray
        if subname in sdict['AperName']:
            mtch = subname == sdict['AperName']
            bounds = [sdict['xstart'].data[mtch][0], sdict['ystart'].data[mtch][0],
                      sdict['xend'].data[mtch][0], sdict['yend'].data[mtch][0]]
            return bounds
        else:
            raise ValueError(("WARNING: {} is not a subarray aperture name present "
                              "in the subarray definition file.".format(subname)))

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
