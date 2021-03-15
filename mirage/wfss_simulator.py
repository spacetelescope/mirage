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

import copy
import os
import sys
import argparse
import logging
import yaml

import numpy as np
from astropy.io import fits
from NIRCAM_Gsim.grism_seed_disperser import Grism_seed
import pysiaf

from .catalogs import spectra_from_catalog
from .seed_image import catalog_seed_image
from .dark import dark_prep
from .logging import logging_functions
from .ramp_generator import obs_generator
from .utils import backgrounds, read_fits
from .utils.constants import CATALOG_YAML_ENTRIES, MEAN_GAIN_VALUES, \
                             LOG_CONFIG_FILENAME, STANDARD_LOGFILE_NAME
from .utils.utils import expand_environment_variable

from .yaml import yaml_update

NIRCAM_GRISM_CROSSING_FILTERS = ['F322W2', 'F277W', 'F356W', 'F444W', 'F250M', 'F300M',
                                 'F335M', 'F360M', 'F410M', 'F430M', 'F323N', 'F405N',
                                 'F466N', 'F470N']
NIRISS_GRISM_CROSSING_FILTERS = ['F200W', 'F150W', 'F140M', 'F158M', 'F115W', 'F090W']

classpath = os.path.dirname(__file__)
log_config_file = os.path.join(classpath, 'logging', LOG_CONFIG_FILENAME)
logging_functions.create_logger(log_config_file, STANDARD_LOGFILE_NAME)


class WFSSSim():
    def __init__(self, paramfiles, SED_file=None, SED_normalizing_catalog_column=None,
                 final_SED_file=None, SED_dict=None, save_dispersed_seed=True, source_stamps_file=None,
                 extrapolate_SED=True, override_dark=None, disp_seed_filename=None, offline=False,
                 create_continuum_seds=True):

        # Set the MIRAGE_DATA environment variable if it is not
        # set already. This is for users at STScI.
        self.env_var = 'MIRAGE_DATA'
        self.datadir = os.environ.get(self.env_var)
        if self.datadir is None:
            raise ValueError(("WARNING: {} environment variable is not set."
                              "This must be set to the base directory"
                              "containing the darks, cosmic ray, PSF, etc"
                              "input files needed for the simulation."
                              "These files must be downloaded separately"
                              "from the Mirage package.".format(self.env_var)))

        # Make sure the input param file(s) is a list
        self.paramfiles = paramfiles
        if isinstance(paramfiles, str):
            self.paramfiles = [self.paramfiles]

        # Set the user-input parameters
        self.create_continuum_seds = create_continuum_seds
        self.SED_file = SED_file
        self.SED_dict = SED_dict
        self.SED_normalizing_catalog_column = SED_normalizing_catalog_column
        self.final_SED_file = final_SED_file
        self.override_dark = override_dark
        self.save_dispersed_seed = save_dispersed_seed
        self.source_stamps_file = source_stamps_file
        self.disp_seed_filename = disp_seed_filename
        self.extrapolate_SED = extrapolate_SED
        self.fullframe_apertures = ["NRCA5_FULL", "NRCB5_FULL", "NIS_CEN"]
        self.offline = offline

        # Make sure the right combination of parameter files and SED file
        # are given
        self.param_checks()

        # Attempt to find the crossing filter and dispersion direction
        # from the input paramfiles. Adjust any imaging mode parameter
        # files to have the mode set to wfss. This will ensure the seed
        # images will be the proper (expanded) dimensions
        self.paramfiles = self.find_param_info()

        # Make sure inputs are correct
        self.check_inputs()

        # Determine whether Mirage will create a final hdf5 file that
        # contains all sources. If only one paramfile is provided, this
        # must be done.
        if len(self.paramfiles) == 1:
            self.create_continuum_seds = True
        else:
            # If spectra are provided, this step must be done
            if self.SED_file is not None or self.SED_dict is not None:
                self.create_continuum_seds = True
            else:
                # If multiple paramfiles are given and no input spectra,
                # then we leave the choice up to the user, as either case
                # is ok.
                pass

    def create(self):
        """MAIN FUNCTION"""
        self.logger = logging.getLogger('mirage.wfss_simulator')
        self.logger.info('\n\nRunning wfss_simulator....\n')
        self.logger.info('using parameter files: ')
        for pfile in self.paramfiles:
            self.logger.info('{}'.format(pfile))

        # Loop over the yaml files and create
        # a direct seed image for each
        imseeds = []
        ptsrc_seeds = []
        galaxy_seeds = []
        extended_seeds = []
        for pfile in self.paramfiles:
            self.logger.info('Running catalog_seed_image for {}'.format(pfile))
            cat = catalog_seed_image.Catalog_seed(offline=self.offline)
            cat.paramfile = pfile
            cat.make_seed()
            imseeds.append(cat.seed_file)
            ptsrc_seeds.append(cat.ptsrc_seed_filename)
            galaxy_seeds.append(cat.galaxy_seed_filename)
            extended_seeds.append(cat.extended_seed_filename)

            # If Mirage is going to produce an hdf5 file of spectra,
            # then we only need a single direct seed image. Note that
            # find_param_info() has reordered the list such that the
            # wfss mode yaml file will be examined first.
            if self.create_continuum_seds:
                # Be sure to include any Mirage-created source catalogs containing
                # ghost sources in the list of catalogs whose sources will be included
                # in the hdf5 file
                if len(cat.ghost_catalogs) > 0:
                    self.catalog_files.extend(cat.ghost_catalogs)
                break

        # Create hdf5 file with spectra of all sources if requested.
        if self.create_continuum_seds:
            det_name = cat.params['Readout']['array_name'].split('_')[0]
            self.SED_file = spectra_from_catalog.make_all_spectra(self.catalog_files, input_spectra=self.SED_dict,
                                                                  input_spectra_file=self.SED_file,
                                                                  extrapolate_SED=self.extrapolate_SED,
                                                                  output_filename=self.final_SED_file,
                                                                  normalizing_mag_column=self.SED_normalizing_catalog_column,
                                                                  module=self.module, detector=det_name)

        # Location of the configuration files needed for dispersion
        loc = os.path.join(self.datadir, "{}/GRISM_{}/current".format(self.instrument,
                                                                      self.instrument.upper()))
        if not os.path.isdir(loc):
            raise ValueError(("{} directory is not present. GRISM_NIRCAM and/or GRISM_NIRISS portion of Mirage reference "
                              "file collection is out of date or not set up correctly. See "
                              "https://mirage-data-simulator.readthedocs.io/en/latest/reference_files.html foe details.".format(loc)))
        self.logger.info("Retrieving WFSS config files from: {}".format(loc))

        # Determine the name of the background file to use, as well as the
        # orders to disperse.
        if self.instrument == 'nircam':
            dmode = 'mod{}_{}'.format(self.module, self.dispersion_direction)
        elif self.instrument == 'niriss':
            dmode = 'GR150{}'.format(self.dispersion_direction)

        # Find the 1d spectrum of the background, based on date or 'low', 'medium', 'high'
        back_wave, back_sig = backgrounds.get_1d_background_spectrum(self.params, self.detector, self.module)

        # Default to extracting all orders
        orders = None

        # Call the disperser separately for each type of object: point sources
        # galaxies, extended objects
        disp_seed = np.zeros((cat.ffsize, cat.ffsize))
        background_done = False

        for seed_files in [ptsrc_seeds, galaxy_seeds, extended_seeds]:
            if seed_files[0] is not None:
                dispersed_objtype_seed = Grism_seed(seed_files, self.crossing_filter,
                                                    dmode, config_path=loc, instrument=self.instrument.upper(),
                                                    extrapolate_SED=self.extrapolate_SED, SED_file=self.SED_file,
                                                    SBE_save=self.source_stamps_file)
                dispersed_objtype_seed.observation(orders=orders)
                dispersed_objtype_seed.disperse(orders=orders)
                # Only include the background in one of the object type seed images
                if not background_done:
                    background_image = dispersed_objtype_seed.disperse_background_1D([back_wave, back_sig])
                    scaling_factor = np.max(background_image)
                    dispersed_objtype_seed.finalize(BackLevel=scaling_factor, tofits=self.background_image_filename)
                    background_done = True
                else:
                    dispersed_objtype_seed.finalize()
                disp_seed += dispersed_objtype_seed.final

        # Disperser output is always full frame. Remove the signal from
        # the refrence pixels now since we know exactly where they are
        disp_seed[0:4, :] = 0.
        disp_seed[2044:, :] = 0.
        disp_seed[:, 0:4] = 0.
        disp_seed[:, 2044:] = 0.

        # Crop to the requested subarray if necessary
        if cat.params['Readout']['array_name'] not in self.fullframe_apertures:
            self.logger.info("Subarray bounds: {}".format(cat.subarray_bounds))
            self.logger.info("Dispersed seed image size: {}".format(disp_seed.shape))
            disp_seed = self.crop_to_subarray(disp_seed, cat.subarray_bounds)

        # Save the dispersed seed image if requested
        # Save in units of e/s, under the idea that this should be a
        # "perfect" noiseless view of the scene that does not depend on
        # detector effects, such as gain.
        if self.save_dispersed_seed:
            self.save_dispersed_seed_image(disp_seed)

        # Convert seed image to ADU/sec to be consistent
        # with other simulator outputs
        if self.instrument == 'niriss':
            gain = MEAN_GAIN_VALUES['niriss']
        elif self.instrument == 'nircam':
            gain = MEAN_GAIN_VALUES['nircam']['lw{}'.format(self.module.lower())]

        disp_seed /= gain

        # Update seed image header to reflect the
        # division by the gain
        cat.seedinfo['units'] = 'ADU/sec'

        # Prepare dark current exposure if
        # needed.
        if self.override_dark is None:
            d = dark_prep.DarkPrep(offline=self.offline)
            d.paramfile = self.wfss_yaml
            d.prepare()

            if len(d.dark_files) == 1:
                obslindark = d.prepDark
            else:
                obslindark = d.dark_files
        else:
            self.logger.info('\n\noverride_dark has been set. Skipping dark_prep.')
            if isinstance(self.override_dark, str):
                self.read_dark_product()
                obslindark = self.prepDark
            elif isinstance(self.override_dark, list):
                obslindark = self.override_dark

        # Combine into final observation
        obs = obs_generator.Observation(offline=self.offline)
        obs.linDark = obslindark
        obs.seed = disp_seed
        obs.segmap = cat.seed_segmap
        obs.seedheader = cat.seedinfo
        #obs.paramfile = y.outname
        obs.paramfile = self.wfss_yaml
        obs.create()

    def param_checks(self):
        """Check parameter file inputs"""
        #if ((len(self.paramfiles) < 2) and (self.SED_file is None)):
        #    raise ValueError(("WARNING: Only one parameter file provided and no SED file. More "
        #                      "yaml files or an SED file needed in order to disperse."))

        if ((len(self.paramfiles) > 1) and (self.SED_file is not None)):
            raise ValueError(("WARNING: When using an SED file, you must provide only one parameter file."))

    def read_dark_product(self):
        # Read in dark product that was produced
        # by dark_prep.py
        self.prepDark = read_fits.Read_fits()
        self.prepDark.file = self.override_dark
        self.prepDark.read_astropy()

    def check_inputs(self):
        """Make sure input parameters are acceptible"""

        # ###################Instrument Name##################
        if self.instrument not in ['nircam', 'niriss']:
            self.invalid('instrument', self.instrument)

        # ###################Module Name##################
        if self.instrument == 'nircam':
            if self.module not in ['A', 'B']:
                self.invalid('module', self.module)
            else:
                self.module = self.module.upper()

            # ###################Crossing Filter##################
            if self.crossing_filter not in NIRCAM_GRISM_CROSSING_FILTERS:
                self.invalid('crossing_filter', self.crossing_filter)

        elif self.instrument == 'niriss':
            if self.crossing_filter not in NIRISS_GRISM_CROSSING_FILTERS:
                self.invalid('crossing_filter', self.crossing_filter)

        # ###################Dispersion Direction##################
        if self.dispersion_direction not in ['R', 'C']:
            self.invalid('dispersion_direction', self.dispersion_direction)

        # ###################Dark File to Use##################
        if self.override_dark is not None:
            dark_list = self.override_dark
            if isinstance(self.override_dark, str):
                dark_list = [self.override_dark]
            for darkfile in dark_list:
                avail = os.path.isfile(darkfile)
                if not avail:
                    raise FileNotFoundError(("WARNING: {} does not exist."
                                             .format(darkfile)))

    def find_param_info(self):
        """Extract dispersion direction and crossing filter from the input
        param files"""
        yamls_to_disperse = []
        self.catalog_files = []
        wfss_files_found = 0
        for i, pfile in enumerate(self.paramfiles):
            with open(pfile, 'r') as infile:
                params = yaml.safe_load(infile)

            cats = [params['simSignals'][cattype] for cattype in CATALOG_YAML_ENTRIES if 'tso_' not in cattype]
            cats = [e for e in cats if e.lower() != 'none']
            self.catalog_files.extend(cats)

            if i == 0:
                self.instrument = params['Inst']['instrument'].lower()
                if self.instrument == 'niriss':
                    self.module = 'N'
                    self.detector = 'NIS'
                elif self.instrument == 'nircam':
                    self.module = params['Readout']['array_name'][3]
                    self.detector = params['Readout']['array_name'][0:5]

            if params['Inst']['mode'].lower() == 'wfss':
                self.wfss_yaml = copy.deepcopy(pfile)
                self.params = params

                # Only 1 input yaml file should be for wfss mode
                wfss_files_found += 1
                if wfss_files_found == 2:
                    raise ValueError("WARNING: only one of the parameter files can be WFSS mode.")
                filter_name = params['Readout']['filter']
                pupil_name = params['Readout']['pupil']
                dispname = ('{}_dispersed_seed_image.fits'.format(params['Output']['file'].split('.fits')[0]))
                self.default_dispersed_filename = os.path.join(params['Output']['directory'], dispname)

                bkgd_output_file = '{}_background_image.fits'.format(params['Output']['file'].split('.fits')[0])
                self.background_image_filename = os.path.join(params['Output']['directory'], bkgd_output_file)

                # In reality, the grism elements are in NIRCam's pupil wheel, and NIRISS's
                # filter wheel. But in the APT xml file, NIRISS grisms are in the pupil
                # wheel and the crossing filter is listed in the filter wheel. At that
                # point, NIRISS and NIRCam are consistent, so let's keep with this reversed
                # information
                if self.instrument == 'niriss':
                    self.crossing_filter = pupil_name.upper()
                    self.dispersion_direction = filter_name[-1].upper()
                elif self.instrument == 'nircam':
                    self.crossing_filter = filter_name.upper()
                    self.dispersion_direction = pupil_name[-1].upper()
                # Prepend the wfss yaml file to the list, so that it is first
                yamls_to_disperse.insert(0, pfile)

            elif params['Inst']['mode'].lower() in ['imaging', 'pom']:
                # If the other yaml files are for imaging mode, we need to update them to
                # be wfss mode so that the resulting seed images have the correct dimensions.
                # Save these modified yaml files to new files.
                params['Inst']['mode'] = 'wfss'
                params['Output']['grism_source_image'] = True
                outdir, basename = os.path.split(pfile)
                modified_file = os.path.join(outdir, 'tmp_update_to_wfss_mode_{}'.format(basename))
                with open(modified_file, 'w') as output:
                    yaml.dump(params, output, default_flow_style=False)
                yamls_to_disperse.append(modified_file)

        if wfss_files_found == 0:
            raise ValueError(("WARNING: No WFSS mode parameter files found. One of the parameter "
                              "files must be wfss mode in order to define grism and crossing filter."))
        return yamls_to_disperse

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
        except (FileNotFoundError, IOError) as e:
            self.logger.info(e)

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
        except (FileNotFoundError, IOError) as e:
            self.logger.info(e)

        mngain = np.nanmedian(image)

        # Set pixels with a gain value of 0 equal to mean
        image[image == 0] = mngain
        # Set any pixels with non-finite values equal to mean
        image[~np.isfinite(image)] = mngain
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
        valid = [(b >= 0 and b < xl) for b in bounds[0:3:2]]
        validy = [(b >= 0 and b < yl) for b in bounds[1:4:2]]
        valid.extend(validy)

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
        except (FileNotFoundError, IOError) as e:
            self.logger.info(e)

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

    def invalid(self, field, value):
        raise ValueError(("WARNING: invalid value for {}: {}"
                          .format(field, value)))

    def save_dispersed_seed_image(self, seed_image):
        """Save the dispersed seed image"""
        primary_hdu = fits.PrimaryHDU()
        image_hdu = fits.ImageHDU(seed_image)
        hdu_list = fits.HDUList([primary_hdu, image_hdu])
        hdu_list[0].header['units'] = 'e/sec'
        hdu_list[1].header['units'] = 'e/sec'
        if self.disp_seed_filename is None:
            self.disp_seed_filename = self.default_dispersed_filename
        hdu_list.writeto(self.disp_seed_filename, overwrite=True)
        self.logger.info(("Dispersed seed image saved to {}".format(self.disp_seed_filename)))

    def add_options(self, parser=None, usage=None):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage, description=("Wrapper for the creation of"
                                                                       " WFSS simulated exposures."))
        parser.add_argument("paramfiles", help=('List of files describing the input parameters and '
                                                'instrument settings to use. (YAML format).'), nargs='+')
        parser.add_argument("--crossing_filter", help=("Name of crossing filter to use in conjunction "
                                                       "with the grism."), default=None)
        parser.add_argument("--module", help="NIRCam module to use for simulation. Use 'A' or 'B'",
                            default=None)
        parser.add_argument("--direction", help=("Direction of dispersion (along rows or along columns). "
                                                 "Use 'R' or 'C'"), default=None)
        parser.add_argument("--override_dark", help=("If supplied, skip the dark preparation step and "
                                                     "use the supplied dark to make the exposure"),
                            default=None)
        parser.add_argument("--extrapolate_SED", help=("If true, the SED created from the filter-averaged "
                                                       "magnitudes will be extrapolated to fill the "
                                                       "wavelngth range of the grism"), action='store_true')
        return parser


if __name__ == '__main__':

    usagestring = ('USAGE: wfss_simualtor.py file1.yaml file2.yaml --crossing_filter F444W --direction R '
                   '--module A')

    obs = WFSSSim()
    parser = obs.add_options(usage=usagestring)
    args = parser.parse_args(namespace=obs)
    obs.create()
