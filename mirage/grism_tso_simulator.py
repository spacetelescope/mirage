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
import yaml

from astropy.io import ascii, fits
import astropy.units as u
import batman
import numpy as np
from NIRCAM_Gsim.grism_seed_disperser import Grism_seed
from scipy.interpolate import interp1d

from mirage import wfss_simulator
from mirage.catalogs import spectra_from_catalog
from mirage.seed_image import catalog_seed_image
from mirage.dark import dark_prep
from mirage.ramp_generator import obs_generator
from mirage.utils import read_fits
from mirage.utils.constants import CATALOG_YAML_ENTRIES
from mirage.utils.file_splitting import find_file_splits
from mirage.utils import utils
from mirage.yaml import yaml_update


class GrismTSO():
    def __init__(self, parameter_file, SED_file=None, SED_normalizing_catalog_column=None,
                 final_SED_file=None, SED_dict=None, save_dispersed_seed=True, source_stamps_file=None,
                 extrapolate_SED=True, override_dark=None, disp_seed_filename=None, orders=["+1", "+2"],
                 create_continuum_seds=True):

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

        # Set the user-input parameters
        self.paramfile = parameter_file
        self.SED_file = SED_file
        self.SED_normalizing_catalog_column = SED_normalizing_catalog_column
        self.final_SED_file = final_SED_file
        self.SED_dict = SED_dict
        self.save_dispersed_seed = save_dispersed_seed
        self.source_stamps_file = source_stamps_file
        self.extrapolate_SED = extrapolate_SED
        self.override_dark = override_dark
        self.disp_seed_filename = disp_seed_filename
        self.orders = orders
        self.create_continuum_seds = create_continuum_seds

        # Make sure the right combination of parameter files and SED file
        # are given
        self.param_checks()

        # Attempt to find the crossing filter and dispersion direction
        # from the input paramfiles. Adjust any imaging mode parameter
        # files to have the mode set to wfss. This will ensure the seed
        # images will be the proper (expanded) dimensions
        #self.paramfiles = self.find_param_info()

        # Make sure inputs are correct
        #self.check_inputs()

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
        # Make 2 copies of the input parameter file, separating the TSO
        # source from the other sources
        self.split_param_file(orig_parameters)

        # Run the catalog_seed_generator on the non-TSO (background) sources
        background_direct = catalog_seed_image.Catalog_seed()
        background_direct.paramfile = self.background_paramfile
        background_direct.make_seed()
        background_segmentation_map = background_direct.seed_segmap

        # Run the disperser on the background sources
        print('\n\nDispersing background sources\n\n')
        background_dispersed = self.run_disperser(background_direct.seed_file, orders=self.orders,
                                                  create_continuum_seds=True, add_background=True)  # finalize=True)

        # Run the catalog_seed_generator on the TSO source
        tso_direct = catalog_seed_image.Catalog_seed()
        tso_direct.paramfile = self.tso_paramfile
        tso_direct.make_seed()
        tso_segmentation_map = tso_direct.seed_segmap
        outside_tso_source = tso_segmentation_map == 0
        tso_segmentation_map[outside_tso_source] = background_segmentation_map[outside_tso_source]

        # Dimensions are (y, x)
        self.seed_dimensions = tso_direct.nominal_dims

        # Read in the transmission spectrum that goes with the TSO source
        tso_params = utils.read_yaml(self.tso_paramfile)
        tso_catalog_file = tso_params['simSignals']['tso_grism_catalog']

        print('TSO catalog file', tso_catalog_file)


        tso_catalog = ascii.read(tso_catalog_file)
        #self.check_tso_catalog_inputs(tso_catalog)

        transmission_file = tso_catalog['Transmission_spectrum'].data
        transmission_spectrum = ascii.read(transmission_file[0])

        # Calculate the total exposure time, including resets, to check
        # against the times provided in the catalog file.
        total_exposure_time = self.calculate_exposure_time() * u.second

        # Check to be sure the start and end times provided in the catalog
        # are enough to cover the length of the exposure.
        tso_catalog = self.tso_catalog_check(tso_catalog, total_exposure_time)

        # Use batman to create lightcurves from the transmission spectrum
        lightcurves, times = self.make_lightcurves(tso_catalog, self.frametime, transmission_spectrum)

        # Determine which frames of the exposure will take place with the unaltered stellar
        # spectrum. This will be all frames where the associated lightcurve is 1.0 everywhere.
        transit_frames, unaltered_frames = self.find_transit_frames(lightcurves)

        # Run the disperser using the original, unaltered stellar spectrum. Set 'cache=True'
        print('\n\nDispersing TSO source\n\n')
        grism_seed_object = self.run_disperser(tso_direct.seed_file, orders=self.orders, add_background=False)  # finalize=False)


        #for i, order in enumerate(self.orders):
        #    if i == 0:
        #        no_transit_signal = copy.deepcopy(grism_seed_object.this_one[order].simulated_image)
        #    else:
        #        no_transit_signal += grism_seed_object.this_one[order].simulated_image
        no_transit_signal = grism_seed_object.final


        print('Shape of no_transit_signal: {}'.format(no_transit_signal.shape))

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
        print('lightcurves shape:', lightcurves.shape)

        print('move below into seed_builder function?')

        ints_per_segment = self.int_segment_indexes[1:] - self.int_segment_indexes[:-1]
        groups_per_segment = self.grp_segment_indexes[1:] - self.grp_segment_indexes[:-1]
        print("Integration segment indexes: ", self.int_segment_indexes)
        print("ints_per_segment: ", ints_per_segment)
        print("Group segment indexes: ", self.grp_segment_indexes)
        print("groups_per_segment: ", groups_per_segment)

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

        # List of all output seed files
        self.seed_files = []

        for i, int_dim in enumerate(ints_per_segment):
            int_start = self.int_segment_indexes[i]
            int_end = self.int_segment_indexes[i+1]
            for j, grp_dim in enumerate(groups_per_segment):
                initial_frame = self.grp_segment_indexes[j]
                # int_dim and grp_dim are the number of integrations and
                # groups in the current segment PART
                print("Integrations: {}, Groups: {}".format(int_dim, grp_dim))
                segment_seed = np.zeros((int_dim, grp_dim, self.seed_dimensions[0], self.seed_dimensions[1]))

                #we need to deal with reset frames here. previous_frame signal should reset to zero,
                #and no dispersion nor segment_seed population is necessary


                for integ in np.arange(int_dim):
                    overall_integration_number = int_start + integ

                    #print('\n\nOverall integration number: {}'.format(overall_integration_number))
                    #print('from: {} + {}'.format(self.int_segment_indexes[i], integ))
                    #print('\n')


                    previous_frame = np.zeros(self.seed_dimensions)
                    for frame in np.arange(grp_dim):
                        print('TOTAL FRAME COUNTER: ', total_frame_counter)
                        print('integ and frame: ', integ, frame)
                        # If a frame is from the part of the lightcurve
                        # with no transit, then the signal in the frame
                        # comes from no_transit_signal
                        if total_frame_counter in unaltered_frames:
                            print("{} is unaltered.".format(total_frame_counter))
                            frame_only_signal = (background_dispersed.final + no_transit_signal) * self.frametime
                        # If the frame is from a part of the lightcurve
                        # where the transit is happening, then call the
                        # cached disperser with the appropriate lightcurve
                        elif total_frame_counter in transit_frames:
                            print("{} is during the transit".format(total_frame_counter))
                            frame_transmission = lightcurves[total_frame_counter, :]
                            trans_interp = interp1d(transmission_spectrum['Wavelength'], frame_transmission)
                            for order in self.orders:
                                grism_seed_object.this_one[order].disperse_all_from_cache(trans_interp)
                            frame_only_signal = (background_dispersed.final + grism_seed_object.final) * self.frametime

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
                    print('RESET FRAME! ', total_frame_counter)

                # At the end of the segment/part, save the segment_seed
                # to a fits file.
                segment_number = np.where(int_end <= self.file_segment_indexes)[0][0]
                if segment_number == previous_segment:
                    segment_part_number += 1
                    self.part_int_start_number = int_start - segment_starting_int_number
                    self.part_frame_start_number = initial_frame
                    self.segment_ints += int_dim
                    self.segment_frames += grp_dim
                else:
                    segment_part_number = 1
                    previous_segment = copy.deepcopy(segment_number)

                    self.segment_frame_start_number = initial_frame
                    self.segment_int_start_number = int_start
                    self.part_int_start_number = 0
                    self.part_frame_start_number = 0
                    segment_starting_int_number = copy.deepcopy(int_start)
                    self.segment_ints = int_dim
                    self.segment_frames = grp_dim




                #self.segment_int_start_number = self.int_segment_indexes[i]
                #self.segment_frame_start_number = self.grp_segment_indexes[j]
                #self.part_int_start_number = overall_integration_number - self.segment_int_start_number  # THIS ISNT RIGHT YET- segment_starting_int_number
                #self.part_frame_start_number = self.grp_segment_indexes[j]

                print('Segment and part numbers: ', segment_number, segment_part_number)
                print('Overall integration number: ', overall_integration_number)
                segment_file_name = '{}seg{}_part{}_seed_image.fits'.format(segment_file_base,
                                                                            str(segment_number).zfill(3),
                                                                            str(segment_part_number).zfill(3))


                print('\n\nSegment int and frame start numbers: {} {}'.format(self.segment_int_start_number, self.segment_frame_start_number))
                print('Part int and frame start numbers (ints and frames within the segment): {} {}'.format(self.part_int_start_number, self.part_frame_start_number))





                print('save the file here. Maybe need to move seed_catalog_images saveSeedImage into utils?')


                # when saving here: segmentation map doesn't matter much since the seed is already dispersed
                #
                tso_seed_header = fits.getheader(tso_direct.seed_file)
                self.save_seed(segment_seed, tso_segmentation_map, tso_seed_header, orig_parameters,
                               segment_number, segment_part_number)





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
        #self.total_seed_segments = (len(self.grp_segment_indexes) - 1) * (len(self.int_segment_indexes) - 1)

        # If the file needs to be split, check to see what the splitting
        # would be in the case of groups rather than frames. This will
        # help align the split files between the seed image and the dark
        # object later (which is split by groups).
        if self.split_seed:
            split_seed_g, group_segment_indexes_g, self.file_segment_indexes = find_file_splits(self.seed_dimensions[1],
                                                                                                self.seed_dimensions[0],
                                                                                                self.numgroups,
                                                                                                self.numints)
        else:
            self.file_segment_indexes = np.array([0, self.numints])
        self.total_seed_segments = len(self.file_segment_indexes) - 1

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

        cats = [parameters['simSignals'][cattype] for cattype in CATALOG_YAML_ENTRIES]
        cats = [e for e in cats if e.lower() != 'none']
        self.catalog_files.extend(cats)

        self.instrument = parameters['Inst']['instrument'].lower()
        self.aperture = parameters['Readout']['array_name']
        try:
            self.namps = parameters['Readout']['namp']
        except KeyError:
            pass
        if self.instrument == 'niriss':
            self.module = None
        elif self.instrument == 'nircam':
            self.module = parameters['Readout']['array_name'][3]
        else:
            raise ValueError("ERROR: Grism TSO mode not supported for {}".format(self.instrument))

        filter_name = parameters['Readout']['filter']
        pupil_name = parameters['Readout']['pupil']
        #dispname = ('{}_dispsersed_seed_image.fits'.format(parameters['Output']['file'].split('.fits')[0]))
        #self.default_dispersed_filename = os.path.join(parameters['Output']['directory'], dispname)

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
        return parameters


    @staticmethod
    def make_lightcurves(catalog, frame_time, transmission_spec):
        """Given a transmission spectrum, create a series of lightcurves
        using ``batman``.

        Parameters
        ----------
        catalog : astropy.table.Table
            Table containing info from the TSO source catalog

        Returns
        -------
        lightcurves : numpy.ndarray
            2D array containing the light curve at each wavelengthin the
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
        start_time = catalog['Start_time'][0] * time_units
        end_time = catalog['End_time'][0] * time_units

        # Convert times to units of seconds to make working
        # with frametimes later easier
        start_time = start_time.to(u.second).value
        end_time = end_time.to(u.second).value
        params.t0 = (catalog['Time_of_inferior_conjunction'][0] * time_units).to(u.second).value  # time of inferior conjunction
        params.per = (catalog['Orbital_period'][0] * time_units).to(u.second).value       # orbital period

        # The time resolution must be one frametime since we will need one
        # lightcurve for each frame later
        time = np.linspace(start_time, end_time, frame_time)  # times at which to calculate light curve
        model = batman.TransitModel(params, time)

        # Step along the transmission spectrum in wavelength space and
        # calculate a lightcurve at each step
        lightcurves = np.ones((len(time), len(transmission_spec['Transmission'])))
        for i, radius in enumerate(transmission_spec['Transmission']):
            params.rp = radius                          # updates planet radius
            new_flux = model.light_curve(params)        # recalculates light curve
            lightcurves[:, i] = new_flux
        return lightcurves, time

    def param_checks(self):
        """Check validity of inputs
        """
        if self.orders not in [["+1"], ["+2"], ["+1", "+2"], None]:
            raise ValueError(("ERROR: Orders to be dispersed must be either None or some combination "
                              "of '+1', '+2'"))

    def run_disperser(self, direct_file, orders=["+1", "+2"], create_continuum_seds=False, add_background=True):  # , finalize=False):
        """
        """
        # Stellar spectrum hdf5 file will be required, so no need to create one here.
        # Create hdf5 file with spectra of all sources if requested
        if create_continuum_seds:
            self.SED_file = spectra_from_catalog.make_all_spectra(self.catalog_files, input_spectra=self.SED_dict,
                                                                  input_spectra_file=self.SED_file,
                                                                  extrapolate_SED=self.extrapolate_SED,
                                                                  output_filename=self.final_SED_file,
                                                                  normalizing_mag_column=self.SED_normalizing_catalog_column)

        # Location of the configuration files needed for dispersion
        loc = os.path.join(self.datadir, "{}/GRISM_{}/".format(self.instrument,
                                                               self.instrument.upper()))

        # Determine the name of the background file to use, as well as the
        # orders to disperse.
        if self.instrument == 'nircam':
            dmode = 'mod{}_{}'.format(self.module, self.dispersion_direction)
            background_file = ("{}_{}_back.fits"
                               .format(self.crossing_filter, dmode))
        elif self.instrument == 'niriss':
            dmode = 'GR150{}'.format(self.dispersion_direction)
            background_file = "{}_{}_medium_background.fits".format(self.crossing_filter.lower(), dmode.lower())
            print('Background file is {}'.format(background_file))
        orders = self.orders


        print(direct_file, self.crossing_filter, dmode, loc, self.instrument)
        print(self.extrapolate_SED, self.SED_file, self.source_stamps_file)


        # Create dispersed seed image from the direct images
        disp_seed = Grism_seed([direct_file], self.crossing_filter,
                               dmode, config_path=loc, instrument=self.instrument.upper(),
                               extrapolate_SED=self.extrapolate_SED, SED_file=self.SED_file,
                               SBE_save=self.source_stamps_file)
        disp_seed.observation()#orders=orders)
        for order in orders:
            disp_seed.this_one[order].disperse_all(cache=True)

        print('DISPERSED SIZE: ', disp_seed.this_one[order].simulated_image.shape)


        #if finalize:
        if add_background:
            #print('FINALIZING')
            #print(background_file)
            #with fits.open(os.path.join('/ifs/jwst/wit/mirage_data/nircam/GRISM_NIRCAM', background_file)) as h:
            #    delme = h[0].data
            #print('background file size: ', delme.shape)
            disp_seed.finalize(Back=background_file)
        else:
            disp_seed.finalize(Back=None, BackLevel=None)
        return disp_seed

    def save_seed(self, seed, segmentation_map, seed_header, params, segment_num, part_num):
        """
        """
        arrayshape = seed.shape
        if len(arrayshape) == 4:
            units = 'ADU'
            integ, grps, yd, xd = arrayshape
            tgroup = self.frametime * (params['Readout']['nframe'] + params['Readout']['nskip'])
            print('Seed image is 4D.')
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

        if self.total_seed_segments == 1:
            self.seed_file = os.path.join(self.basename + '_' + params['Readout'][usefilt] +
                                          '_seed_image.fits')
        else:
            seg_string = str(segment_num).zfill(3)
            part_string = str(part_num).zfill(3)
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
        seed_header['EXSEGTOT'] = self.total_seed_segments
        seed_header['EXSEGNUM'] = segment_num
        seed_header['PART_NUM'] = part_num

        # Total number of integrations and groups in the current segment
        # (combining all parts of the segment)
        seed_header['SEGINT'] = self.segment_ints
        seed_header['SEGGROUP'] = self.segment_frames

        # Total number of integrations and groups in the exposure
        seed_header['EXPINT'] = params['Readout']['nint']
        seed_header['EXPGRP'] = params['Readout']['ngroup']

        # Indexes of the ints and groups where the data in this file go
        seed_header['SEGFRMST'] = self.segment_frame_start_number
        seed_header['SEGFRMED'] = self.segment_frame_start_number + grps - 1
        seed_header['SEGINTST'] = self.segment_int_start_number
        seed_header['SEGINTED'] = self.segment_int_start_number + integ - 1
        seed_header['PTINTSRT'] = self.part_int_start_number
        seed_header['PTFRMSRT'] = self.part_frame_start_number

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

        print("Seed image and segmentation map saved as {}".format(self.seed_file))
        print("Seed image, segmentation map, and metadata available as:")
        print("self.seedimage, self.seed_segmap, self.seedinfo.")

    """
    def seed_builder(self, not_in_transit, in_transit):
        not_in_transit is a list of frame numbers

        UGH, lots of stuff would have to be passed in here, or made
        into class variables

        background_dispersed
        no_transit_signal
        lightcurves
        times
        grism_seed_object


        ints_per_segment = self.int_segment_indexes[:-1] - self.int_segment_indexes[1:]
        groups_per_segment = self.grp_segment_indexes[:-1] - self.grp_segment_indexes[1:]
        print("Integration segment indexes: ", self.int_segment_indexes)
        print("ints_per_segment: ", ints_per_segment)
        print("Group segment indexes: ", self.grp_segment_indexes)
        print("groups_per_segment: ", groups_per_segment)
        total_frame_counter = 0
        previous_segment = 1
        segment_part_number = 0
        for i, int_dim in enumerate(ints_per_segment):
            int_end = self.int_segment_indexes[i+1]
            for j, grp_dim in enumerate(groups_per_segment):
                print("Integrations: {}, Groups: {}".format(int_dim, grp_dim))
                segment_seed = np.zeros((int_dim, grp_dim, self.seed_dimensions[0], self.seed_dimensions[1]))

                #we need to deal with reset frames here. previous_frame signal should reset to zero,
                #and no dispersion nor segment_seed population is necessary


                for integ in np.arange(int_dim):
                    previous_frame = np.zeros(self.seed_dimensions)
                    for frame in np.arange(grp_dim):
                        print('TOTAL FRAME COUNTER: ', total_frame_counter)
                        print('integ and frame: ', integ, frame)
                        # If a frame is from the part of the lightcurve
                        # with no transit, then the signal in the frame
                        # comes from no_transit_signal
                        if total_frame_counter in not_in_transit:
                            print("{} is unaltered.".format(total_frame_counter))
                            frame_only_signal = (background_dispersed.final + no_transit_signal) * self.frametime
                        # If the frame is from a part of the lightcurve
                        # where the transit is happening, then call the
                        # cached disperser with the appropriate lightcurve
                        elif total_frame_counter in in_transit:
                            print("{} is during the transit".format(total_frame_counter))
                            lightcurve = lightcurves[total_frame_counter, :]
                            lc_interp = interp1d(times, lightcurve)
                            for order in self.orders:
                                grism_seed_object.this_one[order].disperse_all_from_cache(lc_interp)
                            frame_only_signal = (background_dispersed.final + grism_seed_object.final) * self.frametime

                        # Now add the signal from this frame to that in the
                        # previous frame in order to arrive at the total
                        # cumulative signal
                        segment_seed[integ, frame, :, :] = previous_frame + frame_only_signal
                        previous_frame = copy.deepcopy(segment_seed[integ, frame, :, :])

                    # At the end of each integration, increment the
                    # total_frame_counter by the number of resets between
                    # integrations
                    total_frame_counter += self.numresets
                    print('RESET FRAME! ', total_frame_counter)

                # At the end of the segment/part, save the segment_seed
                # to a fits file.
                segment_number = np.where(int_end <= self.file_segment_indexes)[0][0]
                if segment_number == previous_segment:
                    segment_part_number += 1
                else:
                    segment_part_number = 1
                    previous_segment = copy.deepcopy(segment_number)

                print('Segment and part numbers: ', segment_number, segment_part_number)
                segment_file_name = '{}seg{}_part{}_seed_image.fits'.format(segment_file_base,
                                                                            str(segment_number).zfill(3),
                                                                            str(segment_part_number).zfill(3))
                #you have to save segment_seed here rather than returning it since you are in the loop
    """

    def split_param_file(self, params):
        """Create 2 copies of the input parameter file. One will contain
        all but the TSO source, while the other will contain only the TSO
        source.
        """
        # Read in the initial parameter file
        #params = read_yaml(self.paramfile)

        file_dir, filename = os.path.split(self.paramfile)
        suffix = filename.split('.')[-1]

        # Copy #1 - contains all source catalogs other than TSO sources
        # and be set to wfss mode.
        background_params = copy.deepcopy(params)
        background_params['simSignals']['tso_grism_catalog'] = 'None'
        background_params['Inst']['mode'] = 'wfss'
        background_params['Output']['grism_source_image'] = True
        background_params['Output']['file'] = params['Output']['file'].replace('.fits',
                                                                               '_background_sources.fits')
        self.background_paramfile = self.paramfile.replace('.{}'.format(suffix),
                                                           '_background_sources.{}'.format(suffix))
        utils.write_yaml(background_params, self.background_paramfile)

        # Copy #2 - contaings only the TSO grism source catalog,
        # is set to wfss mode, and has no background
        params['Inst']['mode'] = 'wfss'
        params['Output']['grism_source_image'] = True
        params['simSignals']['bkgdrate'] = 0.
        params['simSignals']['zodiacal'] = 'None'
        params['simSignals']['scattered'] = 'None'
        other_catalogs = ['pointsource', 'galaxyListFile', 'extended', 'tso_imaging_catalog',
                          'movingTargetList', 'movingTargetSersic', 'movingTargetExtended',
                          'movingTargetToTrack']
        for catalog in other_catalogs:
            params['simSignals'][catalog] = 'None'
        params['Output']['file'] = params['Output']['file'].replace('.fits', '_tso_grism_sources.fits')

        self.tso_paramfile = self.paramfile.replace('.{}'.format(suffix),
                                                    '_tso_grism_sources.{}'.format(suffix))
        utils.write_yaml(params, self.tso_paramfile)

    @staticmethod
    def tso_catalog_check(catalog, exp_time):
        """Check that the start and end times specified in the TSO catalog file (which are
        used to calculate the lightcurves) are long enough to encompass the entire exposure.
        If not, extend the end time to the required time.
        """

        print(catalog['Time_units'][0])
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
            print(('WARNING: Lightcurve duration specified in TSO catalog file is less than '
                   'the total duration of the exposure. Adding extra time to the end of the '
                   'lightcurve to match.'))
            catalog['End_time'][0] = catalog['Start_time'][0] + catalog_total_time.value#.to(catalog_time_units).value

        # Make sure the time of inferior conjunction is betwen
        # the starting and ending times
        if ((catalog['Time_of_inferior_conjunction'][0] < catalog['Start_time'][0]) or
           (catalog['Time_of_inferior_conjunction'][0] > catalog['End_time'][0])):
            raise ValueError(("ERROR: the inferior conjucntion time in the TSO catalog is "
                              "outside the bounds of the starting and ending times."))
        return catalog





