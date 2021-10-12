#! /urs/bin/env python

"""Code for creating simulations of grism time series observations.


1. Make 2 temporary copies of the input yaml file.
   Copy 1 will contain all catalogs except the catalog with the TSO source
       and have the mode set to 'wfss'
   Copy 2 will contain only the TSO source catalog with mode set to 'wfss'

   NOTE THAT STEPS 2 AND 3 CAN BE COMBINED INTO A SINGLE CALL TO WFSS_SIMULATOR
2. Run the catalog_seed_generator on Copy 1. Result will be a padded countrate image
   containing only the non-TSO sources, which have constant brightness with time
   (dare we think about moving targets? - No, not now. Save for a later PR)
3. Run the disperser on the seed image with non-TSO sources. Result will be 2d
   dispersed countrate image

4. Run the catalog_seed_generator on Copy 2. Result will be a padded 2d countrate
   image containing only the TSO source
5. Read in transmission curve
6. Use batman to generate lightcurves from the transmission spectrum for each of
   a grid of wavelengths. (What grid? Dispersed image is 10A per res. element, so
   use that resolution?)
7. Run the disperser using the original, unaltered stellar spectrum. Set 'cache=True'.
   Result will be the dispersed countrate image of the TSO source for times when
   the lightcurve is 1.0 everywhere
8. Create an output array to hold the frame-by-frame seed image. Probably easier if this
   array holds just the signal accumulated in a given frame for that frame, rather than
   the cumulative signal in that frame. Note that there will
   most likely be file-splitting happening at this point...
9. Determine which frames of the exposure will take place with the unaltered stellar
   spectrum. This will be all frames where the associated lightcurve is 1.0 everywhere.
   These frames will be simply copies of the outputs from step 7 plus step 4. To get
   around the issue of file splitting, easier to just keep a list of frame index
   numbers which this situation applies to.
10.For each of the remaining frames, run the disperser with the appropriate lightcurve
   (after running interp1d to turn it into a function). The frame seed will be this
   output plus that from step 4
11.As you go, translate into a cumulative frame by frame seed image, and save to a
   fits file and reset the array variable as you get to the limit of each segment
   part.
12.Run the dark current prep step
13.Run the observation generator

"""


import copy
import os
import sys
import argparse
import datetime
import logging
import yaml

from astropy.io import ascii, fits
import astropy.units as u
from astropy.units.quantity import Quantity
import batman
import numpy as np
from NIRCAM_Gsim.grism_seed_disperser import Grism_seed
import pkg_resources
import pysiaf
from scipy.interpolate import interp1d, interp2d

from mirage import wfss_simulator
from mirage.catalogs import catalog_generator, spectra_from_catalog
from mirage.seed_image import catalog_seed_image
from mirage.dark import dark_prep
from mirage.logging import logging_functions
from mirage.ramp_generator import obs_generator
from mirage.reference_files import crds_tools
from mirage.utils import read_fits
from mirage.utils.constants import CATALOG_YAML_ENTRIES, MEAN_GAIN_VALUES, \
                                   LOG_CONFIG_FILENAME, STANDARD_LOGFILE_NAME
from mirage.utils.file_splitting import find_file_splits, SplitFileMetaData
from mirage.utils import utils, file_io, backgrounds
from mirage.utils.timer import Timer
from mirage.yaml import yaml_update


classpath = os.path.dirname(__file__)
log_config_file = os.path.join(classpath, 'logging', LOG_CONFIG_FILENAME)
logging_functions.create_logger(log_config_file, STANDARD_LOGFILE_NAME)


class GrismTSO():
    def __init__(self, parameter_file, SED_file=None, SED_normalizing_catalog_column=None,
                 final_SED_file=None, save_dispersed_seed=True, source_stamps_file=None,
                 extrapolate_SED=True, override_dark=None, disp_seed_filename=None, orders=["+1", "+2"],
                 lightcurves=None, lightcurve_times=None, lightcurve_wavelengths=None):
        """
        Parameters
        ----------
        parameter_file : str
            Name of input yaml file

        SED_file : str
            Name of hdf5 file containing spectra of some or all sources in
            the simulated data

        SED_normalizing_catalog_column : str
            Name of the column in the ascii source catalogs to use when
            normalizing input spectra in the ``SED_file`` IF some or all
            of the spectra in that file contain normalized fluxes

        final_SED_file : str
            Name of hdf5 file produced by Mirage that contains spectra for
            all simulated sources. In the case where some sources have only
            broadband magnitudes but no entries in the ``SED_file``, Mirage
            will interpolate/extrapolate to create continuum spectra and add
            these to the spectra in the ``SED_file``

        save_dispersed_seed : bool
            Save the dispersed seed image to a fits file

        source_stamps_file : str

        extrapolate_SED : bool
            Allow the disperser software to extrapolate when creating
            continuum spectra from source catalog magnitudes. This is
            done if the wavelength range covered by the filter magnitudes
            does not cover the entire wavelength range of the grism+filter

        override_dark : list
            List of outputs from a prior run of ``dark_prep`` to use when creating
            simulation. If set, the call to ``dark_prep`` will be skipped and these
            darks will be used instead. If None, ``dark_prep`` will be called and
            new dark objects will be created.

        disp_seed_filename : str
            Name of the file to save the dispsersed seed image into

        orders : list
            List of spectral orders to create during dispersion. Default
            for NIRCam is ["+1", "+2"]

        lightcurves : numpy.darray
            2D array containing lightcurves. If this is provided, the call to
            batman will be skipped. Array dimensions should be lightcurves[times, wavelengths]
            where the time values match the units/scale of the Start_time and End_time
            entries in the Mirage-formatted TSO source catalog.

        lightcurve_times : numpy.array
            1D array of times associated with ```lightcurves```. Times should match
            the units/scale of the Start_time and End_time entries in the Mirage-formatted
            TSO source catalog.

        lightcurve_wavelengths : numpy.array
            1D array of wavelengths associated with ```lightcurves```. Wavelengths should
            be in the same units as that of the transmission spectrum referenced in the
            Transmission_spectrum column of the TSO source catalog.
        """

        # Use the MIRAGE_DATA environment variable
        env_var = 'MIRAGE_DATA'
        self.datadir = os.environ.get(env_var)
        if self.datadir is None:
            raise ValueError(("WARNING: {} environment variable is not set."
                              "This must be set to the base directory"
                              "containing the darks, cosmic ray, PSF, etc"
                              "input files needed for the simulation."
                              "These files must be downloaded separately"
                              "from the Mirage package.".format(env_var)))
        self.modpath = pkg_resources.resource_filename('mirage', '')

        # Set the user-input parameters
        self.paramfile = parameter_file
        self.SED_file = SED_file
        self.SED_normalizing_catalog_column = SED_normalizing_catalog_column
        self.final_SED_file = final_SED_file
        self.save_dispersed_seed = save_dispersed_seed
        self.source_stamps_file = source_stamps_file
        self.extrapolate_SED = extrapolate_SED
        self.disp_seed_filename = disp_seed_filename
        self.orders = orders
        self.fullframe_apertures = ["NRCA5_FULL", "NRCB5_FULL", "NIS_CEN"]
        self.override_dark = override_dark
        self.lightcurves = lightcurves
        self.lightcurve_times = lightcurve_times
        self.lightcurve_wavelengths = lightcurve_wavelengths

        # Make sure the right combination of parameter files and SED file
        # are given
        self.param_checks()

        # Initialize timer
        self.timer = Timer()

    def calculate_exposure_time(self):
        """Calculate the total exposure time of the observation being
        simulated. Include time for resets between integrations

        Returns
        -------
        exposure_time : float
            Exposure time for the total exposuure, including reset frames,
            in seconds
        """
        self.frametime = utils.calc_frame_time(self.instrument, self.aperture, self.seed_dimensions[1],
                                               self.seed_dimensions[0], self.namps)
        return self.frametime * self.total_frames

    def create(self):
        """MAIN FUNCTION"""
        # Initialize the log using dictionary from the yaml file
        self.logger = logging.getLogger('mirage.grism_tso_simulator')
        self.logger.info('\n\nRunning grism_tso_simulator....\n')
        self.logger.info('using parameter file: {}'.format(self.paramfile))

        # Get parameters necessary to create the TSO data
        orig_parameters = self.get_param_info()
        subarray_table = utils.read_subarray_definition_file(orig_parameters['Reffiles']['subarray_defs'])
        orig_parameters = utils.get_subarray_info(orig_parameters, subarray_table)
        orig_parameters = utils.read_pattern_check(orig_parameters)

        self.basename = os.path.join(orig_parameters['Output']['directory'],
                                     orig_parameters['Output']['file'][0:-5].split('/')[-1])

        # Determine file splitting information. First get some basic info
        # on the exposure
        self.numints = orig_parameters['Readout']['nint']
        self.numgroups = orig_parameters['Readout']['ngroup']
        self.numframes = orig_parameters['Readout']['nframe']
        self.numskips = orig_parameters['Readout']['nskip']
        self.namps = orig_parameters['Readout']['namp']
        self.numresets = orig_parameters['Readout']['resets_bet_ints']
        self.frames_per_group, self.frames_per_int, self.total_frames = utils.get_frame_count_info(self.numints,
                                                                                                   self.numgroups,
                                                                                                   self.numframes,
                                                                                                   self.numskips,
                                                                                                   self.numresets)

        # Get gain map for later unit conversion
        #gainfile = orig_parameters['Reffiles']['gain']
        #gain, gainheader = file_io.read_gain_file(gainfile)

        # Make 2 copies of the input parameter file, separating the TSO
        # source from the other sources
        self.split_param_file(orig_parameters)

        self.logger.info('Splitting background and TSO source into multiple yaml files.')
        self.logger.info('Running background sources through catalog_seed_image.')
        self.logger.info('background param file is: {}'.format(self.background_paramfile))

        # Stellar spectrum hdf5 file will be required, so no need to create one here.
        # Create hdf5 file with spectra of all sources if requested
        self.final_SED_file = spectra_from_catalog.make_all_spectra(self.catalog_files,
                                                                    input_spectra_file=self.SED_file,
                                                                    extrapolate_SED=self.extrapolate_SED,
                                                                    output_filename=self.final_SED_file,
                                                                    normalizing_mag_column=self.SED_normalizing_catalog_column)

        bkgd_waves, bkgd_fluxes = backgrounds.get_1d_background_spectrum(orig_parameters,
                                                                         self.detector, self.module)

        # Run the catalog_seed_generator on the non-TSO (background) sources. Even if
        # no source catalogs are given, we run using the dummy catalog created earlier,
        # because we need to add the 2D dispersed background at this point.
        self.logger.info('Running catalog_seed_generator on background sources')
        background_direct = catalog_seed_image.Catalog_seed()
        background_direct.paramfile = self.background_paramfile
        background_direct.make_seed()
        background_segmentation_map = background_direct.seed_segmap

        # Run the disperser on the background sources. Add the background
        # signal here as well
        self.logger.info('\n\nDispersing background sources\n\n')

        background_done = False
        background_seed_files = [background_direct.ptsrc_seed_filename,
                                 background_direct.galaxy_seed_filename,
                                 background_direct.extended_seed_filename]
        for seed_file in background_seed_files:
            if seed_file is not None:
                self.logger.info("Dispersing seed image: {}".format(seed_file))

                # Generate the name of a file to save the dispersed background into.
                # We can put this here because this function is called with add_background=True
                # only once in the grism_tso_simulator
                bkgd_output_file = '{}_background_image.fits'.format(orig_parameters['Output']['file'].split('.fits')[0])
                background_image_filename = os.path.join(orig_parameters['Output']['directory'], bkgd_output_file)
                disp = self.run_disperser(seed_file, orders=self.orders,
                                          add_background=not background_done,
                                          background_waves=bkgd_waves,
                                          background_fluxes=bkgd_fluxes,
                                          finalize=True, background_image_output=background_image_filename)
                if not background_done:
                    # Background is added at the first opportunity. At this
                    # point, create an array to hold the final combined
                    # dispersed background
                    background_done = True
                    background_dispersed = copy.deepcopy(disp.final)
                else:
                    background_dispersed += disp.final

        # Run the catalog_seed_generator on the TSO source
        self.logger.info('Running catalog_seed_generator on TSO source')
        tso_direct = catalog_seed_image.Catalog_seed()
        tso_direct.paramfile = self.tso_paramfile
        tso_direct.make_seed()
        tso_segmentation_map = tso_direct.seed_segmap
        outside_tso_source = tso_segmentation_map == 0

        # Add any background sources to the segmentation map
        tso_segmentation_map[outside_tso_source] = background_segmentation_map[outside_tso_source]

        # Dimensions are (y, x)
        self.seed_dimensions = tso_direct.nominal_dims

        # Read in the transmission spectrum that goes with the TSO source
        tso_params = utils.read_yaml(self.tso_paramfile)
        tso_catalog_file = tso_params['simSignals']['tso_grism_catalog']

        tso_catalog = ascii.read(tso_catalog_file)

        transmission_file = tso_catalog['Transmission_spectrum'].data
        transmission_spectrum = ascii.read(transmission_file[0])

        # Calculate the total exposure time, including resets, to check
        # against the times provided in the catalog file.
        total_exposure_time = self.calculate_exposure_time() * u.second

        # Check to be sure the start and end times provided in the catalog
        # are enough to cover the length of the exposure.
        tso_catalog = self.tso_catalog_check(tso_catalog, total_exposure_time)

        if self.lightcurves is None:
            # If the user does not provide a 2D array of lightcurves, use
            # batman to create lightcurves from the transmission spectrum
            self.logger.info("Creating 2D array of lightcurves using batman package")
            lightcurves, times = self.make_lightcurves(tso_catalog, self.frametime, transmission_spectrum)
        else:
            if self.lightcurve_times is None:
                raise ValueError(("User-provided lightcurves are present, but associated times are not (using "
                                  "the 'lightcurve_times' keyword. Unable to continue."))

            if self.lightcurve_wavelengths is None:
                raise ValueError(("User-provided lightcurves are present, but associated wavelengths are not (using "
                                  "the 'lightcurve_wavelengths' keyword. Unable to continue."))

            if len(self.lightcurves.shape) != 2:
                raise ValueError(("User-provided lightcurves needs to be a 2D numpy array with dimensions."))

            self.logger.info(("User-input 2D array of lightcurves, lightcurve times, and wavelengths "
                              "will be used to create the data."))

            # Calculate the times associated with all frames of the exposure
            times = self.make_frame_times(tso_catalog)

            # Units checks for lightcurve times and wavelengths
            if isinstance(self.lightcurve_times, Quantity):
                lightcurve_times = self.lightcurve_times.to(u.second)
            else:
                self.logger.info('No units associated with lightcurve_times. Assuming seconds.')
                lightcurve_times = self.lightcurve_times

            if isinstance(self.lightcurve_wavelengths, Quantity):
                lightcurve_wavelengths = self.lightcurve_wavelengths.to(u.micron)
            else:
                self.logger.info("No units associated with lightcurve_wavelengths. Assuming microns.")
                lightcurve_wavelengths = self.lightcurve_wavelengths

            # If the user has provided a 2D array of lightcurves, plus associated 1D arrays of
            # times and wavelengths, then interpolate those lightcurves onto the grid of frame
            # times and transmission spectrum wavelengths.
            lc_function = interp2d(lightcurve_wavelengths, lightcurve_times, self.lightcurves)
            lightcurves = lc_function(transmission_spectrum['Wavelength'], times)

        # Determine which frames of the exposure will take place with the unaltered stellar
        # spectrum. This will be all frames where the associated lightcurve is 1.0 everywhere.
        transit_frames, unaltered_frames = self.find_transit_frames(lightcurves)
        self.logger.info('Frame numbers containing the transit: {} - {}'.format(np.min(transit_frames), np.max(transit_frames)))

        # Run the disperser using the original, unaltered stellar spectrum. Set 'cache=True'
        self.logger.info('\n\nDispersing TSO source\n\n')
        grism_seed_object = self.run_disperser(tso_direct.seed_file, orders=self.orders,
                                               add_background=False, cache=True, finalize=True)

        # Crop dispersed seed images to correct final subarray size
        #no_transit_signal = grism_seed_object.final
        no_transit_signal = utils.crop_to_subarray(grism_seed_object.final, tso_direct.subarray_bounds)
        background_dispersed = utils.crop_to_subarray(background_dispersed, tso_direct.subarray_bounds)

        # Create a reference pixel mask, and crop to the requeted aperture
        full_maskimage = np.zeros((tso_direct.ffsize, tso_direct.ffsize), dtype=np.int)
        full_maskimage[4:tso_direct.ffsize-4, 4:tso_direct.ffsize-4] = 1.
        refpix_mask = self.crop_to_aperture(orig_parameters, tso_direct.subarray_bounds, full_maskimage)

        # Zero-out the signal in the reference pixels
        background_dispersed *= refpix_mask
        grism_seed_object.final *= full_maskimage

        # Save the dispersed seed images if requested
        if self.save_dispersed_seed:
            h_back = fits.PrimaryHDU(background_dispersed)
            h_back.header['EXTNAME'] = 'BACKGROUND_SOURCES'
            h_tso = fits.ImageHDU(grism_seed_object.final)
            h_tso.header['EXTNAME'] = 'TSO_SOURCE'
            hlist = fits.HDUList([h_back, h_tso])
            disp_filename = '{}_dispersed_seed_images.fits'.format(self.basename)
            hlist.writeto(disp_filename, overwrite=True)
            self.logger.info('\nDispersed seed images (background sources and TSO source) saved to {}.\n\n'
                             .format(disp_filename))

        # Calculate file splitting info
        self.file_splitting()

        # Prepare for creating output files
        segment_file_dir = orig_parameters['Output']['directory']
        if orig_parameters['Readout']['pupil'][0].upper() == 'F':
            usefilt = 'pupil'
        else:
            usefilt = 'filter'
        segment_file_base = orig_parameters['Output']['file'].replace('.fits', '_')
        segment_file_base = '{}_{}_'.format(segment_file_base, orig_parameters['Readout'][usefilt])
        segment_file_base = os.path.join(segment_file_dir, segment_file_base)

        # Loop over frames and integrations up to the size of the segment
        # file.
        ints_per_segment = self.int_segment_indexes[1:] - self.int_segment_indexes[:-1]
        groups_per_segment = self.grp_segment_indexes[1:] - self.grp_segment_indexes[:-1]

        total_frame_counter = 0
        previous_segment = 1
        segment_part_number = 0
        segment_starting_int_number = 0

        self.segment_part_number = 0
        self.segment_ints = 0
        self.segment_frames = 0
        self.segment_frame_start_number = 0
        self.segment_int_start_number = 0
        self.part_int_start_number = 0
        self.part_frame_start_number = 0

        # Get split files' metadata
        split_meta = SplitFileMetaData(self.int_segment_indexes, self.grp_segment_indexes,
                                       self.file_segment_indexes, self.group_segment_indexes_g,
                                       self.frames_per_int, self.frames_per_group, self.frametime)

        # List of all output seed files
        self.seed_files = []

        counter = 0
        for i, int_dim in enumerate(ints_per_segment):
            int_start = self.int_segment_indexes[i]
            int_end = self.int_segment_indexes[i+1]

            # Start timer
            self.timer.start()

            for j, grp_dim in enumerate(groups_per_segment):
                initial_frame = self.grp_segment_indexes[j]
                # int_dim and grp_dim are the number of integrations and
                # groups in the current segment PART
                self.logger.info("\n\nCurrent segment part contains: {} integrations and {} groups.".format(int_dim, grp_dim))
                self.logger.info("Creating frame by frame dispersed signal")
                segment_seed = np.zeros((int_dim, grp_dim, self.seed_dimensions[0], self.seed_dimensions[1]))

                for integ in np.arange(int_dim):
                    overall_integration_number = int_start + integ
                    previous_frame = np.zeros(self.seed_dimensions)

                    for frame in np.arange(grp_dim):
                        #print('TOTAL FRAME COUNTER: ', total_frame_counter)
                        #print('integ and frame: ', integ, frame)
                        # If a frame is from the part of the lightcurve
                        # with no transit, then the signal in the frame
                        # comes from no_transit_signal
                        if total_frame_counter in unaltered_frames:
                            frame_only_signal = (background_dispersed + no_transit_signal) * self.frametime
                        # If the frame is from a part of the lightcurve
                        # where the transit is happening, then call the
                        # cached disperser with the appropriate lightcurve
                        elif total_frame_counter in transit_frames:
                            #print("{} is during the transit".format(total_frame_counter))
                            frame_transmission = lightcurves[total_frame_counter, :]
                            trans_interp = interp1d(transmission_spectrum['Wavelength'], frame_transmission)

                            for order in self.orders:
                                grism_seed_object.this_one[order].disperse_all_from_cache(trans_interp)
                            # Here is where we call finalize on the TSO object
                            # This will update grism_seed_object.final to
                            # contain the correct signal
                            grism_seed_object.finalize(Back=None, BackLevel=None)
                            cropped_grism_seed_object = utils.crop_to_subarray(grism_seed_object.final, tso_direct.subarray_bounds)
                            frame_only_signal = (background_dispersed + cropped_grism_seed_object) * self.frametime

                        # Now add the signal from this frame to that in the
                        # previous frame in order to arrive at the total
                        # cumulative signal
                        segment_seed[integ, frame, :, :] = previous_frame + frame_only_signal
                        previous_frame = copy.deepcopy(segment_seed[integ, frame, :, :])
                        total_frame_counter += 1

                    # At the end of each integration, increment the
                    # total_frame_counter by the number of resets between
                    # integrations
                    total_frame_counter += self.numresets

                # Use the split files' metadata
                self.segment_number = split_meta.segment_number[counter]
                self.segment_ints = split_meta.segment_ints[counter]
                self.segment_frames = split_meta.segment_frames[counter]
                self.segment_part_number = split_meta.segment_part_number[counter]
                self.segment_frame_start_number = split_meta.segment_frame_start_number[counter]
                self.segment_int_start_number = split_meta.segment_int_start_number[counter]
                self.part_int_start_number = split_meta.part_int_start_number[counter]
                self.part_frame_start_number = split_meta.part_frame_start_number[counter]
                counter += 1

                self.logger.info('Overall integration number: {}'.format(overall_integration_number))
                segment_file_name = '{}seg{}_part{}_seed_image.fits'.format(segment_file_base,
                                                                            str(self.segment_number).zfill(3),
                                                                            str(self.segment_part_number).zfill(3))


                self.logger.info('Segment int and frame start numbers: {} {}'.format(self.segment_int_start_number, self.segment_frame_start_number))
                #print('Part int and frame start numbers (ints and frames within the segment): {} {}'.format(self.part_int_start_number, self.part_frame_start_number))

                # Disperser output is always full frame. Crop to the
                # requested subarray if necessary
                if orig_parameters['Readout']['array_name'] not in self.fullframe_apertures:
                    self.logger.info("Dispersed seed image size: {}".format(segment_seed.shape))

                    # segment_seed has already been cropped, so no need to crop here...
                    #segment_seed = utils.crop_to_subarray(segment_seed, tso_direct.subarray_bounds)
                    #gain = utils.crop_to_subarray(gain, tso_direct.subarray_bounds)

                # Segmentation map will be centered in a frame that is larger
                # than full frame by a factor of sqrt(2), so crop appropriately
                self.logger.info('Cropping segmentation map to appropriate aperture')
                segy, segx = tso_segmentation_map.shape
                dx = int((segx - tso_direct.nominal_dims[1]) / 2)
                dy = int((segy - tso_direct.nominal_dims[0]) / 2)
                segbounds = [tso_direct.subarray_bounds[0] + dx, tso_direct.subarray_bounds[1] + dy,
                             tso_direct.subarray_bounds[2] + dx, tso_direct.subarray_bounds[3] + dy]


                tso_segmentation_map = utils.crop_to_subarray(tso_segmentation_map, segbounds)

                # Convert seed image to ADU/sec to be consistent
                # with other simulator outputs
                gain = MEAN_GAIN_VALUES['nircam']['lw{}'.format(self.module.lower())]
                segment_seed /= gain

                # Update seed image header to reflect the
                # division by the gain
                tso_direct.seedinfo['units'] = 'ADU/sec'

                # Save the seed image. Save in units of ADU/sec
                self.logger.info('Saving seed image')
                tso_seed_header = fits.getheader(tso_direct.seed_file)
                self.save_seed(segment_seed, tso_segmentation_map, tso_seed_header, orig_parameters) #,
                               #segment_number, segment_part_number)

            # Stop the timer and record the elapsed time
            self.timer.stop(name='seg_{}'.format(str(i+1).zfill(4)))

            # If there is more than one segment, provide an estimate of processing time
            self.logger.info('\n\nSegment {} out of {} complete.'.format(i+1, len(ints_per_segment)))
            if len(ints_per_segment) > 1:
                time_per_segment = self.timer.sum(key_str='seg_') / (i+1)
                estimated_remaining_time = time_per_segment * (len(ints_per_segment) - (i+1)) * u.second
                time_remaining = np.around(estimated_remaining_time.to(u.minute).value, decimals=2)
                finish_time = datetime.datetime.now() + datetime.timedelta(minutes=time_remaining)
                self.logger.info(('\nEstimated time remaining in this exposure: {} minutes. '
                                  'Projected finish time: {}\n'.format(time_remaining, finish_time)))

        # Prepare dark current exposure if
        # needed.
        if not self.override_dark:
            self.logger.info('Running dark prep')
            d = dark_prep.DarkPrep()
            d.paramfile = self.paramfile
            d.prepare()
            use_darks = d.dark_files
        else:
            self.logger.info('\noverride_dark is set. Skipping call to dark_prep and using these files instead.')
            use_darks = self.override_dark

        # Combine into final observation
        self.logger.info('Running observation generator')
        obs = obs_generator.Observation()
        obs.linDark = use_darks
        obs.seed = self.seed_files
        obs.segmap = tso_segmentation_map
        obs.seedheader = tso_direct.seedinfo
        obs.paramfile = self.paramfile
        obs.create()

        self.logger.info('\nGrism TSO simulator complete')
        logging_functions.move_logfile_to_standard_location(self.paramfile, STANDARD_LOGFILE_NAME)

    def file_splitting(self):
        """Determine file splitting details based on calculated data
        volume
        """
        frames_per_group = self.frames_per_int / self.numgroups
        self.split_seed, self.grp_segment_indexes, self.int_segment_indexes = find_file_splits(self.seed_dimensions[1],
                                                                                               self.seed_dimensions[0],
                                                                                               self.frames_per_int,
                                                                                               self.numints,
                                                                                               frames_per_group=frames_per_group)

        # If the file needs to be split, check to see what the splitting
        # would be in the case of groups rather than frames. This will
        # help align the split files between the seed image and the dark
        # object later (which is split by groups).
        if self.split_seed:
            forced_ints_per_file = int(self.frames_per_int / self.numgroups) * (self.int_segment_indexes[1] - self.int_segment_indexes[0])
            split_seed_g, self.group_segment_indexes_g, self.file_segment_indexes = find_file_splits(self.seed_dimensions[1],
                                                                                                self.seed_dimensions[0],
                                                                                                self.numgroups,
                                                                                                self.numints,
                                                                                                force_delta_int=forced_ints_per_file)

            # In order to avoid the case of having a single integration
            # in the final file, which leads to rate rather than rateints
            # files in the pipeline, check to be sure that the integration
            # splitting indexes indicate the last split isn't a single
            # integration
            if len(self.file_segment_indexes) > 2:
                delta_int = self.file_segment_indexes[1:] - self.file_segment_indexes[0: -1]
                if delta_int[-1] == 1 and delta_int[0] != 1:
                    self.file_segment_indexes[-2] -= 1
                    self.logger.info('Adjusted to avoid single integration: ', self.file_segment_indexes)

            # More adjustments related to segment numbers. We need to compare
            # the integration indexes for the seed images vs those for the final
            # data and make sure that there are no segments in the final data that
            # have no corresponding integrations from the seed images
            # Example: integration_segment_indexes = [0, 7, 8], and
            # self.file_segment_indexes = [0, 6, 8] - due to applying the adjustment
            # above. In this case as you step through integration_segment_indexes,
            # you see that (skipping 0), 7 and 8 both fall in the 6-8 bin in
            # self.file_segment_indexes. Nothing falls into the 0-7 bin, which
            # corresponds to segment 1. In this case, we need to adjust
            # integration_segment_indexes to make sure that all segments have data
            # associated with them.
            segnum_check = []
            for intnum in self.int_segment_indexes[1:]:
                segnum_check.append(np.where(intnum <= self.file_segment_indexes)[0][0])
            maxseg = max(segnum_check)
            for i in range(1, maxseg + 1):
                if i not in segnum_check:
                    self.int_segment_indexes = copy.deepcopy(self.file_segment_indexes)

        else:
            self.file_segment_indexes = np.array([0, self.numints])
            self.group_segment_indexes_g = np.array([0, self.numgroups])

        self.total_seed_segments = len(self.file_segment_indexes) - 1
        self.total_seed_segments_and_parts = (len(self.int_segment_indexes) - 1) * (len(self.grp_segment_indexes) - 1)

    @staticmethod
    def find_transit_frames(lightcurve_collection):
        """
        """
        no_transit = []
        transit = []
        for row in range(lightcurve_collection.shape[0]):
            lc = lightcurve_collection[row, :]
            if np.all(lc == 1.0):
                no_transit.append(row)
            else:
                transit.append(row)
        return transit, no_transit

    def get_param_info(self):
        """Collect needed information out of the parameter file. Check
        parameter values for correctness

        Returns
        -------
        parameters : dict
            Nested dictionary of parameters in the input yaml file
        """
        self.catalog_files = []
        parameters = utils.read_yaml(self.paramfile)

        # Create dictionary to use when looking in CRDS for reference files
        crds_dict = crds_tools.dict_from_yaml(parameters)

        # Expand reference file entries to be full path names
        parameters = utils.full_paths(parameters, self.modpath, crds_dict)

        # Find which background source catalogs are present
        supported_bkgd_types = ['pointsource', 'galaxyListFile', 'extended']

        try:
            CATALOG_YAML_ENTRIES.remove('tso_grism_catalog')
        except ValueError:
            pass

        try:
            CATALOG_YAML_ENTRIES.remove('tso_imaging_catalog')
        except ValueError:
            pass

        cats = []
        for cattype in CATALOG_YAML_ENTRIES:
            if cattype in supported_bkgd_types:
                if parameters['simSignals'][cattype].lower() != 'none':
                    cats.append(parameters['simSignals'][cattype])
            else:
                if parameters['simSignals'][cattype].lower() != 'none':
                    self.logger.info(parameters['simSignals'][cattype].lower(), type(parameters['simSignals'][cattype]))
                    raise ValueError('{} catalog: {} is unsupported in grism TSO mode.'.format(cattype, parameters['simSignals'][cattype]))

        self.catalog_files.extend(cats)

        if len(self.catalog_files) == 0:
            # If no background source catalogs are given, create a dummy point
            # source catalog and add it to the list. Without this, the
            # creation of the final SED file would fail. Put the source
            # down near the SEP and with a magnitude such that it won't
            # disturb anything
            filter_name = parameters['Readout']['filter'].lower()
            dummy_ptsrc = catalog_generator.PointSourceCatalog(ra=[0.], dec=[-89.])
            dummy_ptsrc.add_magnitude_column([40], instrument='nircam', filter_name=filter_name, magnitude_system='abmag')
            dummy_cat = 'dummy_ptsrc.cat'
            dummy_ptsrc.save(dummy_cat)
            self.catalog_files.append(dummy_cat)
            parameters['simSignals']['pointsource'] = dummy_cat

        self.instrument = parameters['Inst']['instrument'].lower()
        self.aperture = parameters['Readout']['array_name']
        try:
            self.namps = parameters['Readout']['namp']
        except KeyError:
            pass
        if self.instrument == 'niriss':
            self.module = None
            self.detetor = 'NIS'
        elif self.instrument == 'nircam':
            self.module = parameters['Readout']['array_name'][3]
            self.detector = parameters['Readout']['array_name'][0:5]
        else:
            raise ValueError("ERROR: Grism TSO mode not supported for {}".format(self.instrument))

        filter_name = parameters['Readout']['filter']
        pupil_name = parameters['Readout']['pupil']

        # In reality, the grism elements are in NIRCam's pupil wheel, and NIRISS's
        # filter wheel. But in the APT xml file, NIRISS grisms are in the pupil
        # wheel and the crossing filter is listed in the filter wheel. At that
        # point, NIRISS and NIRCam are consistent, so let's keep with this reversed
        # information
        if self.instrument == 'niriss':
            self.crossing_filter = pupil_name.upper()
            self.dispersion_direction = filter_name[-1].upper()
        elif self.instrument == 'nircam':
            if 'CLEAR' in pupil_name.upper():
                raise ValueError("ERROR: pupil cannot be CLEAR. It must be GRISMR or GRISMC.")
            self.crossing_filter = filter_name.upper()
            self.dispersion_direction = pupil_name[-1].upper()
        return parameters


    def make_lightcurves(self, catalog, frame_time, transmission_spec):
        """Given a transmission spectrum, create a series of lightcurves
        using ``batman``.

        Parameters
        ----------
        catalog : astropy.table.Table
            Table containing info from the TSO source catalog

        Returns
        -------
        lightcurves : numpy.ndarray
            2D array containing the light curve at each wavelength in the
            transmission spectrum
        """
        params = batman.TransitParams()

        # planet radius (in units of stellar radii)
        # DUMMY VALUE FOR MODEL INSTANTIATION
        params.rp = 0.1

        params.a = catalog['Semimajor_axis_in_stellar_radii']  # semi-major axis (in units of stellar radii)
        params.inc = catalog['Orbital_inclination_deg']        # orbital inclination (in degrees)
        params.ecc = catalog['Eccentricity']                   # eccentricity
        params.w = catalog['Longitude_of_periastron']          # longitude of periastron (in degrees)
        params.limb_dark = catalog['Limb_darkening_model']     # limb darkening model

        # Limb darkening coefficients [u1, u2, u3, u4]
        params.u = [np.float(e) for e in catalog['Limb_darkening_coeffs'][0].split(',')]

        # Get the time units from the catalog
        time_units = u.Unit(catalog['Time_units'][0])
        time = self.make_frame_times(catalog)

        params.t0 = (catalog['Time_of_inferior_conjunction'][0] * time_units).to(u.second).value  # time of inferior conjunction
        params.per = (catalog['Orbital_period'][0] * time_units).to(u.second).value       # orbital period

        model = batman.TransitModel(params, time)

        # Step along the transmission spectrum in wavelength space and
        # calculate a lightcurve at each step
        lightcurves = np.ones((len(time), len(transmission_spec['Transmission'])))
        for i, radius in enumerate(transmission_spec['Transmission']):
            params.rp = radius                          # updates planet radius
            new_flux = model.light_curve(params)        # recalculates light curve
            lightcurves[:, i] = new_flux

        # Save the 2D transmission data
        h0 = fits.PrimaryHDU(lightcurves)
        hdulist = fits.HDUList([h0])
        outfile = '{}{}'.format(self.basename, '_normalized_lightcurves_vs_time.fits')
        hdulist.writeto(outfile, overwrite=True)
        self.logger.info('2D array of lightcurves vs time saved to: {}'.format(outfile))

        return lightcurves, time


    def make_frame_times(self, catalog):
        """Create an array of times associated with all frames

        Parameters
        ----------
        catalog : astropy.table.Table
            Table containing info from the TSO source catalog

        Returns
        -------
        time : numpy.array
            Array of times corresponding to each frame
        """
        # Get the time units from the catalog
        time_units = u.Unit(catalog['Time_units'][0])
        start_time = catalog['Start_time'][0] * time_units
        end_time = catalog['End_time'][0] * time_units

        # Convert times to units of seconds to make working
        # with frametimes later easier
        start_time = start_time.to(u.second).value
        end_time = end_time.to(u.second).value

        # The time resolution must be one frametime since we will need one
        # lightcurve for each frame later
        time = np.arange(start_time, end_time, self.frametime)  # times at which to calculate light curve
        return time


    def param_checks(self):
        """Check validity of inputs
        """
        if self.orders not in [["+1"], ["+2"], ["+1", "+2"], None]:
            raise ValueError(("ERROR: Orders to be dispersed must be either None or some combination "
                              "of '+1', '+2'"))

    @staticmethod
    def crop_to_aperture(params, sub_bounds, array):
        """Create a mask showing the locations of the reference pixels

        Parameters
        ----------
        params : dict
            Nested dictionary contianing observation parameters, as read
            in from an input yaml file

        sub_bounds : list
            4 element list of subarray boundary locations in full frame
            coordinates. [xstart, ystart, xend, yend]. Generally created
            from siaf_interface.get_siaf_information

        array : numpy.ndarray
            2D full frame array

        Returns
        -------
        cropped : numpy.ndarray
            ```array``` cropped to the requested aperture
        """
        # Crop the mask to match the requested output array
        ap_suffix = params['Readout']['array_name'].split('_')[1]
        if ap_suffix not in ['FULL', 'CEN']:
            cropped = array[sub_bounds[1]:sub_bounds[3] + 1,
                            sub_bounds[0]:sub_bounds[2] + 1]
        return cropped


    def run_disperser(self, direct_file, orders=["+1", "+2"], add_background=True,
                      background_waves=None, background_fluxes=None, cache=False, finalize=False,
                      background_image_output='dispersed_background.fits'):
        """Run the disperser on the given direct seed image.

        Parameters
        ----------
        direct_file : str
            Name of file containing direct seed image

        orders : list
            Orders to include in the dispersed image.

        add_background : bool
            If True, a 2D dispersed background image is created and
            added

        background_waves : numpy.ndarray
            1D array of wavelengths (in mocrons) to be used when creating
            the background image

        background_fluxes : numpy.ndarray
            1D array of fluxes (in units of cgs f_lambda) to be used
            when creating the background image

        cache : bool
            Whether or not to cache the dispersed object. If it is cached,
            then it can be used later

        finalize : bool
            If True, call the finalize function of the disperser in order to
            create the final dispersed image of the object

        background_image_output : str
            Name of a fits file to which the final dispersed background image
            will be saved.

        Returns
        -------
        disp_seed : numpy.ndarray
            2D array containing the dispersed seed image
        """
        # Location of the configuration files needed for dispersion
        loc = os.path.join(self.datadir, "{}/GRISM_{}/current/".format(self.instrument,
                                                                       self.instrument.upper()))
        if not os.path.isdir(loc):
            raise ValueError(("{} directory is not present. GRISM_NIRCAM and/or GRISM_NIRISS portion of Mirage reference "
                              "file collection is out of date or not set up correctly. See "
                              "https://mirage-data-simulator.readthedocs.io/en/latest/reference_files.html foe details.".format(loc)))
        self.logger.info("Retrieving grism-related config files from: {}".format(loc))

        # Determine the name of the background file to use, as well as the
        # orders to disperse.
        dmode = 'mod{}_{}'.format(self.module, self.dispersion_direction)
        orders = self.orders

        # Create dispersed seed image from the direct images
        disp_seed = Grism_seed([direct_file], self.crossing_filter,
                               dmode, config_path=loc, instrument=self.instrument.upper(),
                               extrapolate_SED=self.extrapolate_SED, SED_file=self.final_SED_file,
                               SBE_save=self.source_stamps_file)
        disp_seed.observation(orders=orders)

        # Cache the dispersion if requested. This will allow you to
        #disperse the same object again later with a different lightcurve
        if cache:
            for order in orders:
                disp_seed.this_one[order].disperse_all(cache=True)
        else:
            disp_seed.disperse(orders=orders)

        # Only finalize and/or add the background if requested.
        if finalize:
            if add_background:
                # Create a 2D background image and find the maximum value. Then apply this
                # as the scaling when using the pre-computed 2D background image
                background_image = disp_seed.disperse_background_1D([background_waves, background_fluxes])
                scaling_factor = np.max(background_image)
                disp_seed.finalize(BackLevel=scaling_factor, tofits=background_image_output)
            else:
                disp_seed.finalize(Back=None, BackLevel=None)
        return disp_seed

    def save_seed(self, seed, segmentation_map, seed_header, params):
        """Save the 4D dispersed seed image to a fits file

        Parameters
        ----------
        seed : numpy.ndarray
            4D seed image

        segmentation_map : numpy.ndarray
            2D segmentation image

        seed_header : astropy.fits.PrimaryHDU
            Header object to use as the basis of the saved seed image

        params : dict
            Nested dictionary of instrument/observation parameters. From
            reading in a Mirage input yaml file.
        """
        arrayshape = seed.shape
        if len(arrayshape) == 4:
            units = 'ADU'
            integ, grps, yd, xd = arrayshape
            tgroup = self.frametime * (params['Readout']['nframe'] + params['Readout']['nskip'])
            self.logger.info('Seed image is 4D.')
        else:
            raise ValueError('Only 4D seed images supported. This seed image is {}D'.format(len(arrayshape)))

        xcent_fov = xd / 2
        ycent_fov = yd / 2

        seed_header['XCENTER'] = xcent_fov
        seed_header['YCENTER'] = ycent_fov
        seed_header['UNITS'] = units
        seed_header['TGROUP'] = tgroup

        if params['Readout']['pupil'][0].upper() == 'F':
            usefilt = 'pupil'
        else:
            usefilt = 'filter'

        if self.total_seed_segments_and_parts == 1:
            self.seed_file = os.path.join(self.basename + '_' + params['Readout'][usefilt] +
                                          '_seed_image.fits')
        else:
            seg_string = str(self.segment_number).zfill(3)
            part_string = str(self.segment_part_number).zfill(3)
            self.seed_file = os.path.join(self.basename + '_' + params['Readout'][usefilt] +
                                          '_seg{}_part{}_seed_image.fits'.format(seg_string, part_string))

        self.seed_files.append(self.seed_file)

        # These cannot come from the values in the direct seed image, since
        # they are for the padded array.
        seed_header['NOMXDIM'] = seed.shape[1]
        seed_header['NOMYDIM'] = seed.shape[0]
        seed_header['NOMXSTRT'] = 1
        seed_header['NOMXEND'] = seed.shape[1]
        seed_header['NOMYSTRT'] = 1
        seed_header['NOMYEND'] = seed.shape[0]

        # Observations with high data volumes (e.g. moving targets, TSO)
        # can be split into multiple "segments" in order to cap the
        # maximum file size
        seed_header['EXSEGTOT'] = self.total_seed_segments  # Total number of segments in exp
        seed_header['EXSEGNUM'] = self.segment_number       # Segment number of this file
        seed_header['PART_NUM'] = self.segment_part_number  # Segment part number of this file

        # Total number of integrations and groups in the current segment
        # (combining all parts of the segment)
        seed_header['SEGINT'] = self.segment_ints
        seed_header['SEGGROUP'] = self.segment_frames

        # Total number of integrations and groups in the exposure
        seed_header['EXPINT'] = params['Readout']['nint']
        seed_header['EXPGRP'] = params['Readout']['ngroup']

        # Indexes of the ints and groups where the data in this file go
        # Frame and integration indexes of the segment within the exposure
        seed_header['SEGFRMST'] = self.segment_frame_start_number
        seed_header['SEGFRMED'] = self.segment_frame_start_number + grps - 1
        seed_header['SEGINTST'] = self.segment_int_start_number
        seed_header['SEGINTED'] = self.segment_int_start_number + self.segment_ints - 1

        # Frame and integration indexes of the part within the segment
        seed_header['PTINTSRT'] = self.part_int_start_number
        seed_header['PTFRMSRT'] = self.part_frame_start_number

        # The 1b pipeline adds two header keywords to the data indicating the position
        # of the source on the detector. Currently the pipeline assumes that the
        # source is at the reference location of the aperture being used. Since Mirage
        # outputs are equivalent to what is produced by the 1b pipeline, we need to
        # add these keywords manually.
        siaf = pysiaf.Siaf('nircam')[params['Readout']['array_name']]
        seed_header['XREF_SCI'] = siaf.XSciRef
        seed_header['YREF_SCI'] = siaf.YSciRef
        seed_header['TSOVISIT'] = True

        # Save
        h0 = fits.PrimaryHDU()
        h1 = fits.ImageHDU(seed, name='DATA')
        h2 = fits.ImageHDU(segmentation_map)
        h2.header['EXTNAME'] = 'SEGMAP'

        # Put the keywords into the 0th and 1st extension headers
        for key in seed_header:
            h0.header[key] = seed_header[key]
            h1.header[key] = seed_header[key]

        hdulist = fits.HDUList([h0, h1, h2])
        hdulist.writeto(self.seed_file, overwrite=True)

        self.logger.info("Seed image and segmentation map saved as {}".format(self.seed_file))
        self.logger.info("Seed image, segmentation map, and metadata available as:")
        self.logger.info("self.seedimage, self.seed_segmap, self.seedinfo.\n\n")

    def split_param_file(self, params):
        """Create 2 copies of the input parameter file. One will contain
        all but the TSO source, while the other will contain only the TSO
        source.

        Parameters
        ----------
        params : dict
            Nested dictionary of instrument/observation parameters. From
            reading in a Mirage input yaml file.
        """
        file_dir, filename = os.path.split(self.paramfile)
        suffix = filename.split('.')[-1]

        # Check to see if there are any catalogs for background objects
        bkgd_cats = ['pointsource', 'galaxyListFile', 'extended']
        present = [True if params['simSignals'][cat] is not None else False for cat in bkgd_cats]

        # Copy #1 - contains all source catalogs other than TSO sources
        # and be set to wfss mode.
        if any(present):
            background_params = copy.deepcopy(params)
            background_params['simSignals']['tso_grism_catalog'] = 'None'
            background_params['Inst']['mode'] = 'wfss'
            background_params['Output']['grism_source_image'] = True
            background_params['Output']['file'] = params['Output']['file'].replace('.fits',
                                                                                   '_background_sources.fits')
            self.background_paramfile = self.paramfile.replace('.{}'.format(suffix),
                                                               '_background_sources.{}'.format(suffix))
            utils.write_yaml(background_params, self.background_paramfile)
        else:
            self.background_paramfile = None

        # Copy #2 - contaings only the TSO grism source catalog,
        # is set to wfss mode, and has no background
        if params['simSignals']['tso_grism_catalog'] is None:
            raise ValueError("tso_grism_catalog is not populated in the input yaml. No TSO source to simulate.")

        tso_params = copy.deepcopy(params)
        tso_params['Inst']['mode'] = 'wfss'
        tso_params['Output']['grism_source_image'] = True
        tso_params['simSignals']['bkgdrate'] = 0.
        tso_params['simSignals']['zodiacal'] = 'None'
        tso_params['simSignals']['scattered'] = 'None'
        other_catalogs = ['pointsource', 'galaxyListFile', 'extended', 'tso_imaging_catalog',
                          'movingTargetList', 'movingTargetSersic', 'movingTargetExtended',
                          'movingTargetToTrack']
        for catalog in other_catalogs:
            tso_params['simSignals'][catalog] = 'None'
        tso_params['simSignals']['pointsource'] = tso_params['simSignals']['tso_grism_catalog']
        tso_params['Output']['file'] = tso_params['Output']['file'].replace('.fits', '_tso_grism_sources.fits')

        self.tso_paramfile = self.paramfile.replace('.{}'.format(suffix),
                                                    '_tso_grism_sources.{}'.format(suffix))
        utils.write_yaml(tso_params, self.tso_paramfile)

    def tso_catalog_check(self, catalog, exp_time):
        """Check that the start and end times specified in the TSO catalog file (which are
        used to calculate the lightcurves) are long enough to encompass the entire exposure.
        If not, extend the end time to the required time.

        Parameters
        ----------
        catalog : astropy.table.Table
            Table contianing catalog of TSO source information

        exp_time : float
            Total time of the exposure, in seconds

        Returns
        -------
        catalog : astropy.table.Table
            Potentially modified catalog where the length of the lightcurve
            has been expanded to fill the entire ``exp_time``
        """
        time_unit_str = catalog['Time_units'][0]

        # Catch common unit errors
        if time_unit_str.lower() in ['seconds', 'minutes', 'hours', 'days']:
            time_unit_str = time_unit_str[0:-1]
            catalog['Time_units'][0] = time_unit_str
        catalog_time_units = u.Unit(time_unit_str)
        catalog_total_time = (catalog['End_time'][0] - catalog['Start_time'][0]) * catalog_time_units

        # Make sure the lightcurve time is at least as long as the exposure
        # time
        if exp_time > catalog_total_time:
            self.logger.info(('WARNING: Lightcurve duration specified in TSO catalog file is less than '
                              'the total duration of the exposure. Adding extra time to the end of the '
                              'lightcurve to match.'))
            catalog['End_time'][0] = catalog['Start_time'][0] + catalog_total_time.value

        # Make sure the time of inferior conjunction is betwen
        # the starting and ending times
        if ((catalog['Time_of_inferior_conjunction'][0] < catalog['Start_time'][0]) or
           (catalog['Time_of_inferior_conjunction'][0] > catalog['End_time'][0])):
            raise ValueError(("ERROR: the inferior conjucntion time in the TSO catalog is "
                              "outside the bounds of the starting and ending times."))
        return catalog
