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
import yaml

import numpy as np
from astropy.io import fits
from NIRCAM_Gsim.grism_seed_disperser import Grism_seed
import pysiaf
from scipy.stats import sigmaclip

from .catalogs import spectra_from_catalog
from .seed_image import catalog_seed_image
from .dark import dark_prep
from .ramp_generator import obs_generator
from .utils import backgrounds, read_fits
from .utils.flux_cal import fluxcal_info
from .utils.constants import CATALOG_YAML_ENTRIES, MEAN_GAIN_VALUES, NIRISS_GRISM_THROUGHPUT_FACTOR
from .utils.utils import expand_environment_variable, get_filter_throughput_file

from .yaml import yaml_update

NIRCAM_GRISM_CROSSING_FILTERS = ['F322W2', 'F277W', 'F356W', 'F444W', 'F250M', 'F300M',
                                 'F335M', 'F360M', 'F410M', 'F430M', 'F323N', 'F405N',
                                 'F466N', 'F470N']
NIRISS_GRISM_CROSSING_FILTERS = ['F200W', 'F150W', 'F140M', 'F158M', 'F115W', 'F090W']


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

        # Loop over the yaml files and create
        # a direct seed image for each
        imseeds = []
        ptsrc_seeds = []
        galaxy_seeds = []
        extended_seeds = []
        for pfile in self.paramfiles:
            print('Running catalog_seed_image for {}'.format(pfile))
            cat = catalog_seed_image.Catalog_seed(offline=self.offline)
            cat.paramfile = pfile
            cat.make_seed()
            imseeds.append(cat.seed_file)
            ptsrc_seeds.append(cat.ptsrc_seed_filename)
            galaxy_seeds.append(cat.galaxy_seed_filename)
            extended_seeds.append(cat.extended_seed_filename)

            # Save the flat field file associated with the wfss yaml so
            # that we can apply it to the dispersed seed image if requested
            if pfile == self.wfss_yaml:
                self.flatfield = cat.flatfield

            # If Mirage is going to produce an hdf5 file of spectra,
            # then we only need a single direct seed image. Note that
            # find_param_info() has reordered the list such that the
            # wfss mode yaml file will be examined first.
            if self.create_continuum_seds:
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
        loc = os.path.join(self.datadir, "{}/GRISM_{}/".format(self.instrument,
                                                               self.instrument.upper()))

        # Determine the name of the background file to use, as well as the
        # orders to disperse.
        if self.instrument == 'nircam':
            dmode = 'mod{}_{}'.format(self.module, self.dispersion_direction)
            if self.params['simSignals']['use_dateobs_for_background']:
                print("Generating background spectrum for observation date: {}".format(self.params['Output']['date_obs']))
                back_wave, back_sig = backgrounds.day_of_year_background_spectrum(self.params['Telescope']['ra'],
                                                                                  self.params['Telescope']['dec'],
                                                                                  self.params['Output']['date_obs'])
            else:
                if isinstance(self.params['simSignals']['bkgdrate'], str):
                    if self.params['simSignals']['bkgdrate'].lower() in ['low', 'medium', 'high']:
                        print("Generating background spectrum based on requested level of: {}".format(self.params['simSignals']['bkgdrate']))
                        back_wave, back_sig = backgrounds.low_med_high_background_spectrum(self.params, self.detector,
                                                                                           self.module)
                    else:
                        raise ValueError("ERROR: Unrecognized background rate. Must be one of 'low', 'medium', 'high'")
                else:
                    raise ValueError(("ERROR: WFSS background rates must be one of 'low', 'medium', 'high', "
                                      "or use_dateobs_for_background must be True "))

        elif self.instrument == 'niriss':
            dmode = 'GR150{}'.format(self.dispersion_direction)
            background_file = "{}_{}_medium_background.fits".format(self.crossing_filter.lower(),
                                                                    dmode.lower())

            if isinstance(self.params['simSignals']['bkgdrate'], str):
                if self.params['simSignals']['bkgdrate'].lower() in ['low', 'medium', 'high']:
                    #scaling_factor = backgrounds.niriss_background_scaling(self.params, self.detector, self.module)

                    usefilt = 'pupil'

                    siaf_instance = pysiaf.Siaf('niriss')[self.params['Readout']['array_name']]
                    vegazp, photflam, photfnu, pivot_wavelength = fluxcal_info(self.params, usefilt, self.detector, self.module)

                    if os.path.split(self.params['Reffiles']['filter_throughput'])[1] == 'placeholder.txt' or self.params['Reffiles']['filter_throughput'] == 'config':
                        filter_file = get_filter_throughput_file(self.instrument, self.params['Readout'][usefilt])
                    else:
                        filter_file = self.params['Reffiles']['filter_throughput']

                    scaling_factor = backgrounds.calculate_background(self.params['Telescope']['ra'],
                                                                      self.params['Telescope']['dec'],
                                                                      filter_file,
                                                                      self.params['simSignals']['use_dateobs_for_background'],
                                                                      MEAN_GAIN_VALUES['niriss'], siaf_instance,
                                                                      level=self.params['simSignals']['bkgdrate'])

                    # Having the grism in the beam reduces the throughput by 20%.
                    # Mulitply that into the scaling factor
                    scaling_factor *= NIRISS_GRISM_THROUGHPUT_FACTOR

                    # Translate from ADU/sec/pix to e-/sec/pix since that is
                    # what the disperser works with
                    scaling_factor *= MEAN_GAIN_VALUES['niriss']

                else:
                    raise ValueError("ERROR: Unrecognized background rate. String value must be one of 'low', 'medium', 'high'")
            elif np.isreal(self.params['simSignals']['bkgdrate']):
                # The bkgdrate entry in the input yaml file is described as
                # the desired signal in ADU/sec/pixel IN A DIRECT IMAGE
                # Since we want e-/sec/pixel here for the disperser, multiply
                # by the gain as well as the throughput factor for the grism.
                scaling_factor = self.params['simSignals']['bkgdrate'] * MEAN_GAIN_VALUES['niriss'] * NIRISS_GRISM_THROUGHPUT_FACTOR

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
                    if self.instrument == 'nircam':
                        background_image = dispersed_objtype_seed.disperse_background_1D([back_wave, back_sig])
                        dispersed_objtype_seed.finalize(Back=background_image, BackLevel=None)
                    else:
                        # BackLevel is used as such: background / max(background) * BackLevel
                        # So we need to either set BackLevel equal to the requested level
                        # NOT THE RATIO OF THAT TO MEDIUM, or we need to open the background
                        # file and multiply it by the ratio of the requested level to medium.
                        # The former isn't quite correct because it'll be scaling the maximum
                        # value in the image to "low" or "high", rather than the median
                        full_background_file = os.path.join(loc, background_file)
                        background_image = fits.getdata(full_background_file)

                        # Before scaling the background image by the scaling_factor
                        # we need to normalize by the sigma-clipped mean value. This is
                        # because the background files were produced and scaled to the
                        # ETC "medium" level at some arbirtrary pointing, but the
                        # "medium" level is pointing-dependent. Current background files
                        # are scaled such that the "medium" value from the ETC is the
                        # sigma-clipped mean value.
                        clip, lo, hi = sigmaclip(background_image, low=3, high=3)
                        background_mean = np.mean(clip)
                        background_image = background_image / background_mean * scaling_factor
                        dispersed_objtype_seed.finalize(Back=background_image, BackLevel=None)

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
            print("Subarray bounds: {}".format(cat.subarray_bounds))
            print("Dispersed seed image size: {}".format(disp_seed.shape))
            disp_seed = self.crop_to_subarray(disp_seed, cat.subarray_bounds)

            # Segmentation map will be centered in a frame that is larger
            # than full frame by a factor of sqrt(2), so crop appropriately
            print("Need to make this work for subarrays...")
            segy, segx = cat.seed_segmap.shape
            dx = int((segx - 2048) / 2)
            dy = int((segy - 2048) / 2)
            segbounds = [cat.subarray_bounds[0] + dx, cat.subarray_bounds[1] + dy,
                         cat.subarray_bounds[2] + dx, cat.subarray_bounds[3] + dy]
            cat.seed_segmap = self.crop_to_subarray(cat.seed_segmap, segbounds)

        # If it is a NIRCam observation, multiply the dispsersed seed
        # image by the flat field. For NIRISS, the flat was multiplied
        # in to the direct seed image, as they need to have their
        # occulting spots in the direct seed image.
        if self.instrument == 'nircam':
            disp_seed *= self.flatfield

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
            obslindark = d.prepDark
        else:
            self.read_dark_product()
            obslindark = self.darkPrep

        # Using the first of the imaging seed image yaml
        # files as a base, adjust to create the yaml file
        # for the creation of the final dispersed
        # integration

        # I think we won't need this anymore assuming that
        # one of the input yaml files is for wfss mode
        #y = yaml_update.YamlUpdate()
        #y.file = self.paramfiles[0]
        #if self.instrument == 'nircam':
        #    y.filter = self.crossing_filter
        #    y.pupil = 'GRISM' + self.direction
        #elif self.instrument == 'niriss':
        #    y.filter = 'GR150' + self.direction
        #    y.pupil = self.crossing_filter
        #y.outname = ("wfss_dispersed_{}_{}.yaml"
        #             .format(dmode, self.crossing_filter))
        #y.run()

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
            avail = os.path.isfile(self.override_dark)
            if not avail:
                raise FileNotFoundError(("WARNING: {} does not exist."
                                         .format(self.override_dark)))

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

    """
    def generate_background_spectrum(self, ra, dec, day_of_year=None):
        Call jwst_backgrounds in order to produce an estimate of background
        versus wavelength for a particular pointing

        background = jbt.background(ra, dec, 4.)
        # If the user wants a background signal from a particular day,
        # then extract that array here
        if self.params['simSignals']['use_dateobs_for_background']:
            obsdate = datetime.datetime.strptime(self.params['Output']['date_obs'], '%Y-%m-%d')
            obs_dayofyear = obsdate.timetuple().tm_yday
            if obs_dayofyear not in background.bkg_data['calendar']:
                raise ValueError(("ERROR: The requested RA, Dec is not observable on {}. Either "
                                  "specify a different day, or set simSignals:use_dateobs_for_background "
                                  "to False.".format(self.params['Output']['date_obs'])))
            #match = obs_dayofyear == background.bkg_data['calendar']
            #back_wave = background.bkg_data['wave_array']
            #back_sig = background.bkg_data['total_bg'][match, :]
            background_file = 'background.txt'
            background = jbt.get_background(ra, dec, 4., plot_background=False, plot_bathtub=False,
                                            thisday=day_of_year, write_background=True, write_bathtub=False,
                                            background_file=background_file)
        else:
            # If the user wants a background chosen from the three general
            # levels (low, medium, of high) then we need to examine the
            # returned background spectra for all days, determine from the
            # histogram (of background values at the pivot wavelength) what
            # "low", "medium", and "high" correspond to, and then pick out
            # the spectrum from a day that is close to the requested level
            vegazp, photflam, photfnu, pivot_wavelength = fluxcal_info(self.params, usefilt, detector, module)
            wave_diff = np.abs(background.bkg_data['wave_array'] - pivot_wavelength)
            bkgd_wave = np.where(wave_diff == np.min(wave_diff))[0][0]
            bkgd_at_pivot = background.bkg_data['total_bg'][:, bkgd_wave]

            # Now sort and determine the low/medium/high levels
            low, medium, high = find_low_med_high(bkgd_at_pivot)

            # Find the value based on the level in the yaml file
            background_level = self.params['simSignals']['bkgdrate'].lower()
            if background_level == "low":
                level_value = low
            elif background_level == "medium":
                level_value = medium
            elif background_level == "high":
                level_value = high
            else:
                raise ValueError(("ERROR: Unrecognized background value: {}. Must be low, mediumn, or high"
                                  .format(self.params['simSignals']['bkgdrate'])))

            # Find the day with the background at the pivot wavelength that
            # is closest to the value associated with the requested level
            diff = np.abs(bkgd_at_pivot - level_value)
            mindiff = np.where(diff == np.min(diff))[0][0]
            background_spec = background.bkg_data['total_bg'][mindiff, :]

        return background_file
    """

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
            print(e)

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
            print(e)

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
            print(e)

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
        print(("Dispersed seed image saved to {}".format(self.disp_seed_filename)))

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
