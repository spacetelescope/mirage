#! /usr/bin env python

'''
Reorganization of ramp simulator code. This class is used to construct
a "seed image" from a point source catalog. This seed image is a
noiseless countrate image containing only sources. No noise, no
cosmic rays.
'''

import argparse
import datetime
import sys
import glob
import os
import copy
import re
import shutil
from yaml.scanner import ScannerError

import math
import yaml
import time
import pkg_resources
import asdf
import scipy.signal as s1
from scipy.ndimage import rotate
import numpy as np
from photutils import detect_sources
from astropy.coordinates import SkyCoord
from astropy.io import fits, ascii
from astropy.table import Table, Column
from astropy.modeling.models import Shift, Sersic2D, Polynomial2D, Mapping
import astropy.units as u
import pysiaf

from . import moving_targets
from . import segmentation_map as segmap
import mirage
from mirage.catalogs.catalog_generator import TSO_GRISM_INDEX
from mirage.seed_image import tso
from ..reference_files import crds_tools
from ..utils import backgrounds
from ..utils import rotations, polynomial, read_siaf_table, utils
from ..utils import set_telescope_pointing_separated as set_telescope_pointing
from ..utils import siaf_interface, file_io
from ..utils.constants import CRDS_FILE_TYPES, MEAN_GAIN_VALUES
from ..utils.flux_cal import fluxcal_info
from ..psf.psf_selection import get_gridded_psf_library, get_psf_wings
from ..utils.constants import grism_factor, TSO_MODES
from mirage.utils.file_splitting import find_file_splits, SplitFileMetaData
from ..psf.segment_psfs import (get_gridded_segment_psf_library_list,
                                get_segment_offset, get_segment_library_list)


INST_LIST = ['nircam', 'niriss', 'fgs']
MODES = {'nircam': ["imaging", "ts_imaging", "wfss", "ts_grism"],
         'niriss': ["imaging", "ami", "pom", "wfss"],
         'fgs': ["imaging"]}
TRACKING_LIST = ['sidereal', 'non-sidereal']

inst_abbrev = {'nircam': 'NRC',
               'niriss': 'NIS',
               'fgs': 'FGS'}

ALLOWEDOUTPUTFORMATS = ['DMS']
WFE_OPTIONS = ['predicted', 'requirements']
WFEGROUP_OPTIONS = np.arange(5)


class Catalog_seed():
    def __init__(self, offline=False):
        """Instantiate the Catalog_seed class

        Parameters
        ----------
        offline : bool
            If True, the check for the existence of the MIRAGE_DATA
            directory is skipped. This is primarily for Travis testing
        """
        self.offline = offline

        # Locate the module files, so that we know where to look
        # for config subdirectory
        self.modpath = pkg_resources.resource_filename('mirage', '')

        # Get the location of the MIRAGE_DATA environment
        # variable, so we know where to look for darks, CR,
        # PSF files, etc later
        self.env_var = 'MIRAGE_DATA'
        datadir = utils.expand_environment_variable(self.env_var, offline=offline)

        # Check that CRDS-related environment variables are set correctly
        self.crds_datadir = crds_tools.env_variables()

        # self.coord_adjust contains the factor by which the
        # nominal output array size needs to be increased
        # (used for WFSS mode), as well as the coordinate
        # offset between the nominal output array coordinates,
        # and those of the expanded array. These are needed
        # mostly for WFSS observations, where the nominal output
        # array will not sit centered in the expanded output image.
        self.coord_adjust = {'x': 1., 'xoffset': 0, 'y': 1., 'yoffset': 0}

        # NIRCam rough noise values. Used to make educated guesses when
        # creating segmentation maps
        self.single_ron = 6.  # e-/read
        self.grism_background = 0.25  # e-/sec

        # keep track of the maximum index number of sources
        # in the various catalogs. We don't want to end up
        # with multiple sources having the same index numbers
        self.maxindex = 0

    def make_seed(self):
        """MAIN FUNCTION"""
        print('\n\nRunning catalog_seed_image..\n')

        # Read in input parameters and quality check
        self.seed_files = []
        self.readParameterFile()

        # Create dictionary to use when looking in CRDS for reference files
        self.crds_dict = crds_tools.dict_from_yaml(self.params)

        # Expand param entries to full paths where appropriate
        self.params = utils.full_paths(self.params, self.modpath, self.crds_dict, offline=self.offline)
        self.filecheck()
        self.basename = os.path.join(self.params['Output']['directory'],
                                     self.params['Output']['file'][0:-5].split('/')[-1])
        self.params['Output']['file'] = self.basename + self.params['Output']['file'][-5:]
        self.subdict = utils.read_subarray_definition_file(self.params['Reffiles']['subarray_defs'])
        self.check_params()
        self.params = utils.get_subarray_info(self.params, self.subdict)
        self.coord_transform = self.read_distortion_reffile()
        self.grism_direct_factor = grism_factor(self.params['Inst']['instrument'].lower())
        self.expand_catalog_for_segments = bool(self.params['simSignals']['expand_catalog_for_segments'])
        self.add_psf_wings = self.params['simSignals']['add_psf_wings']

        # If the output is a direct image to be dispersed, expand the size
        # of the nominal FOV so the disperser can account for sources just
        # outside whose traces will fall into the FOV
        if (self.params['Output']['grism_source_image']) or (self.params['Inst']['mode'] in ["pom"]):
            self.calcCoordAdjust()

        # Image dimensions
        self.nominal_dims = np.array([self.subarray_bounds[3] - self.subarray_bounds[1] + 1,
                                      self.subarray_bounds[2] - self.subarray_bounds[0] + 1])
        self.output_dims = (self.nominal_dims * np.array([self.coord_adjust['y'],
                                                          self.coord_adjust['x']])).astype(np.int)

        print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
        print('output dimensions are: {}'.format(self.output_dims))
        print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')

        # Read in the flat field file so it can be used later
        self.prepare_flat()

        # calculate the exposure time of a single frame, based on the size of the subarray
        self.frametime = utils.calc_frame_time(self.params['Inst']['instrument'],
                                               self.params['Readout']['array_name'],
                                               self.nominal_dims[1], self.nominal_dims[0],
                                               self.params['Readout']['namp'])
        print("Frametime is {}".format(self.frametime))

        # Get information on the number of frames
        self.get_frame_count_info()

        # If the total number of frames/pixels is larger than the file
        # size cap can handle, then we need to break the exposure up
        # into multiple chunks
        frames_per_group = self.frames_per_integration / self.params['Readout']['ngroup']
        if self.params['Inst']['mode'] in ['ts_imaging']:
            split_seed, group_segment_indexes, integration_segment_indexes = find_file_splits(self.output_dims[1],
                                                                                              self.output_dims[0],
                                                                                              self.frames_per_integration,
                                                                                              self.params['Readout']['nint'],
                                                                                              frames_per_group=frames_per_group)

            # If the file needs to be split, check to see what the splitting
            # would be in the case of groups rather than frames. This will
            # help align the split files between the seed image and the dark
            # object later (which is split by groups).
            if split_seed:
                split_seed_g, self.group_segment_indexes_g, self.file_segment_indexes = find_file_splits(self.output_dims[1],
                                                                                                    self.output_dims[0],
                                                                                                    self.params['Readout']['ngroup'],
                                                                                                    self.params['Readout']['nint'])


                #print('\nOutputs:')
                #print(split_seed, group_segment_indexes, integration_segment_indexes)
                #print(split_seed_g, self.group_segment_indexes_g, self.file_segment_indexes)
                #print(self.file_segment_indexes)
                #print('')

                # In order to avoid the case of having a single integration
                # in the final file, which leads to rate rather than rateints
                # files in the pipeline, check to be sure that the integration
                # splitting indexes indicate the last split isn't a single
                # integration
                if len(self.file_segment_indexes) > 2:
                    delta_int = self.file_segment_indexes[1:] - self.file_segment_indexes[0: -1]
                    if delta_int[-1] == 1 and delta_int[0] != 1:
                        self.file_segment_indexes[-2] -= 1
                        print('Adjusted to avoid single integration: ', self.file_segment_indexes)

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
                for intnum in integration_segment_indexes[1:]:
                    segnum_check.append(np.where(intnum <= self.file_segment_indexes)[0][0])
                maxseg = max(segnum_check)
                for i in range(1, maxseg + 1):
                    if i not in segnum_check:
                        integration_segment_indexes = copy.deepcopy(self.file_segment_indexes)

            else:
                self.file_segment_indexes = integration_segment_indexes

        else:
            split_seed = False
            group_segment_indexes = [0, self.params['Readout']['ngroup']]
            integration_segment_indexes = [0, self.params['Readout']['nint']]
            self.file_segment_indexes = copy.deepcopy(integration_segment_indexes)
            self.total_seed_segments = 1
            self.segment_number = 1
            self.segment_part_number = 1

        self.total_seed_segments = len(self.file_segment_indexes) - 1
        self.total_seed_segments_and_parts = (len(integration_segment_indexes) - 1) * (len(group_segment_indexes) - 1)

        # Read in the PSF library file corresponding to the detector and filter
        # For WFSS simulations, use the PSF libraries with the appropriate CLEAR element
        self.psf_pupil = self.params['Readout']['pupil']
        self.psf_filter = self.params['Readout']['filter']
        if self.params['Readout']['pupil'].lower() in ['grismr', 'grismc']:
            self.psf_pupil = 'CLEAR'
        if self.params['Readout']['filter'].lower() in ['gr150r', 'gr150c']:
            self.psf_filter = 'CLEAR'

        # If reading in a normal PSF, use get_gridded_psf_library to get a
        # single photutils.griddedPSFModel object
        if not self.expand_catalog_for_segments:
            self.psf_library = get_gridded_psf_library(
                self.params['Inst']['instrument'], self.detector, self.psf_filter, self.psf_pupil,
                self.params['simSignals']['psfwfe'],
                self.params['simSignals']['psfwfegroup'],
                self.params['simSignals']['psfpath'])
            self.psf_library_core_y_dim, self.psf_library_core_x_dim = self.psf_library.data.shape[-2:]
            self.psf_library_oversamp = self.psf_library.oversampling

        # If reading in segment PSFs, use get_gridded_segment_psf_library_list
        # to get a list of photutils.griddedPSFModel objects
        else:
            self.psf_library = get_gridded_segment_psf_library_list(
                self.params['Inst']['instrument'], self.detector, self.psf_filter,
                self.params['simSignals']['psfpath'], pupilname=self.psf_pupil)
            self.psf_library_core_y_dim, self.psf_library_core_x_dim = self.psf_library[0].data.shape[-2:]
            self.psf_library_oversamp = 1

        # Set the psf core dimensions to actually be 2 rows and columns
        # less than the dimensions in the library file. This is because
        # we will later evaluate the library using these core dimensions.
        # If we were to evaluate a library that is 51x51 pixels using a
        # 51x51 pixel grid, then if the source is centered close to the
        # edge of the central pixel, you can end up with an zeroed out
        # edge row or column in the evaluated array. So we do this to be
        # sure that we are evaluating the library with a slightly smaller
        # array than the array in the library.
        self.psf_library_core_x_dim = np.int(self.psf_library_core_x_dim / self.psf_library_oversamp) - \
            self.params['simSignals']['gridded_psf_library_row_padding']
        self.psf_library_core_y_dim = np.int(self.psf_library_core_y_dim / self.psf_library_oversamp) - \
            self.params['simSignals']['gridded_psf_library_row_padding']

        if self.add_psf_wings is True:
            self.psf_wings = get_psf_wings(self.params['Inst']['instrument'], self.detector,
                                           self.psf_filter, self.psf_pupil,
                                           self.params['simSignals']['psfwfe'],
                                           self.params['simSignals']['psfwfegroup'],
                                           os.path.join(self.params['simSignals']['psfpath'], 'psf_wings'))

            # Read in the file that defines PSF array sizes based on magnitude
            self.psf_wing_sizes = ascii.read(self.params['simSignals']['psf_wing_threshold_file'])
            max_wing_size = self.psf_wings.shape[0]
            too_large = np.where(np.array(self.psf_wing_sizes['number_of_pixels']) > max_wing_size)[0]
            if len(too_large) > 0:
                print(('Some PSF sizes in {} are larger than the PSF library file dimensions. '
                       'Resetting these values in the table to be equal to the PSF dimensions'
                       .format(os.path.basename(self.params['simSignals']['psf_wing_threshold_file']))))
                self.psf_wing_sizes['number_of_pixels'][too_large] = max_wing_size

        # Prepare for saving seed images. Set default values for the
        # keywords that must appear in the headers of any seed image
        # file

        # Attributes related to file splitting, which is done in cases
        # where the data volume is too large. Currently this is only
        # implemented for imaging TSO. The only other mode where it
        # will need to be implemented is moving tagets
        self.segment_number = 1
        self.segment_part_number = 1
        self.segment_frames = self.frames_per_integration  # self.seedimage.shape[1]
        self.segment_ints = self.params['Readout']['nint']  # self.seedimage.shape[0]
        self.segment_frame_start_number = 0
        self.segment_int_start_number = 0
        self.part_int_start_number = 0
        self.part_frame_start_number = 0

        # For imaging mode, generate the countrate image using the catalogs
        if self.params['Telescope']['tracking'].lower() != 'non-sidereal':
            print('Creating signal rate image of synthetic inputs.')
            self.seedimage, self.seed_segmap = self.create_sidereal_image()
            outapp = ''

        # If we are tracking a non-sidereal target, then
        # everything in the catalogs needs to be streaked across
        # the detector
        if self.params['Telescope']['tracking'].lower() == 'non-sidereal':
            print('Creating signal ramp of synthetic inputs')
            self.seedimage, self.seed_segmap = self.non_sidereal_seed()
            outapp = '_nonsidereal_target'

        # If non-sidereal targets are requested (KBOs, asteroids, etc,
        # create a RAPID integration which includes those targets
        mov_targs_ramps = []
        if (self.runStep['movingTargets'] | self.runStep['movingTargetsSersic']
                | self.runStep['movingTargetsExtended']):
            print(("Creating signal ramp of sources that are moving with "
                   "respect to telescope tracking."))
            trailed_ramp, trailed_segmap = self.make_trailed_ramp()
            outapp += '_trailed_sources'

            # Now we need to expand frameimage into a ramp
            # so we can add the trailed objects
            print('Combining trailed object ramp with that containing tracked targets')
            if self.params['Telescope']['tracking'].lower() != 'non-sidereal':
                self.seedimage = self.combineSimulatedDataSources('countrate', self.seedimage, trailed_ramp)
            else:
                self.seedimage = self.combineSimulatedDataSources('ramp', self.seedimage, trailed_ramp)
            self.seed_segmap += trailed_segmap

        # MASK IMAGE
        # Create a mask so that we don't add signal to masked pixels
        # Initially this includes only the reference pixels
        # Keep the mask image equal to the true subarray size, since this
        # won't be used to make a requested grism source image
        if self.params['Inst']['mode'] not in ['wfss', 'ts_grism']:
            self.maskimage = np.zeros((self.ffsize, self.ffsize), dtype=np.int)
            self.maskimage[4:self.ffsize-4, 4:self.ffsize-4] = 1.

            # crop the mask to match the requested output array
            ap_suffix = self.params['Readout']['array_name'].split('_')[1]
            if ap_suffix not in ['FULL', 'CEN']:
                self.maskimage = self.maskimage[self.subarray_bounds[1]:self.subarray_bounds[3] + 1,
                                                self.subarray_bounds[0]:self.subarray_bounds[2] + 1]

        # TSO imaging mode here
        if self.params['Inst']['mode'] in ['ts_imaging']:
            self.create_tso_seed(split_seed=split_seed, integration_splits=integration_segment_indexes,
                                 frame_splits=group_segment_indexes)
        else:

            # For seed images to be dispersed in WFSS mode,
            # embed the seed image in a full frame array. The disperser
            # tool does not work on subarrays
            aperture_suffix = self.params['Readout']['array_name'].split('_')[-1]
            if ((self.params['Inst']['mode'] in ['wfss', 'ts_grism']) & \
               (aperture_suffix not in ['FULL', 'CEN'])):
                self.seedimage, self.seed_segmap = self.pad_wfss_subarray(self.seedimage, self.seed_segmap)

            # For NIRISS POM data, extract the central 2048x2048
            if self.params['Inst']['mode'] in ["pom"]:
                self.seedimage, self.seed_segmap = self.extract_full_from_pom(self.seedimage, self.seed_segmap)

            if self.params['Inst']['mode'] not in ['wfss', 'ts_grism']:
                # Multiply the mask by the seed image and segmentation map in
                # order to reflect the fact that reference pixels have no signal
                # from external sources. Seed images to be dispersed do not have
                # this step applied.
                self.seedimage *= self.maskimage
                self.seed_segmap *= self.maskimage

            # Save the combined static + moving targets ramp
            if self.params['Inst']['instrument'].lower() != 'fgs':
                self.seed_file = os.path.join(self.basename + '_' + self.params['Readout'][self.usefilt] + '_final_seed_image.fits')
            else:
                self.seed_file = '{}_final_seed_image.fits'.format(self.basename)

            # Define self.seed_files for non-TSO observations
            if len(self.seed_files) == 0:
                self.seed_files = [self.seed_file]

            self.saveSeedImage(self.seedimage, self.seed_segmap, self.seed_file)
            print("Final seed image and segmentation map saved as {}".format(self.seed_file))
            print("Seed image, segmentation map, and metadata available as:")
            print("self.seedimage, self.seed_segmap, self.seedinfo.")

        # Return info in a tuple
        # return (self.seedimage, self.seed_segmap, self.seedinfo)

    def create_tso_seed(self, split_seed=False, integration_splits=-1, frame_splits=-1):
        """Create a seed image for time series sources. Just like for moving
        targets, the seed image will be 4-dimensional

        Parameters
        ----------
        split_seed : bool
            Whether or not the seed image will be split between multiple files

        integration_splits : list
            List containing the first integration number within each of the
            split files. In addition, the final number in the list should be
            the index of the final integration.

        frame_splits : list
            Similar to ``integration_splits`` but for frame numbers. This will
            be used in the case where the readout pattern contains many frames
            per group

        Returns
        -------
        seed : numpy.ndarray
            4D array seed image.

        segmap : numpy.ndarray
            2D segmentation map
        """
        # Read in the TSO catalog. This is only for TSO mode, not
        # grism TSO, which will need a different catalog format.
        if self.params['Inst']['mode'].lower() == 'ts_imaging':
            tso_cat = self.get_point_source_list(self.params['simSignals']['tso_imaging_catalog'], source_type='ts_imaging')

            # Create lists of seed images and segmentation maps for all
            # TSO objects
            tso_seeds = []
            tso_segs = []
            tso_lightcurves = []
            for source in tso_cat:
                # Let's set the TSO source to be a unique index number.
                # Otherwise the index number will be modified to be one
                # greater than the max value after working on the background
                # sources
                source['index'] = TSO_GRISM_INDEX

                # Place row in an empty table
                t = tso_cat[:0].copy()
                t.add_row(source)

                # Create the seed image contribution from each source
                ptsrc_seed, ptsrc_seg = self.make_point_source_image(t)

                # Multiply by the flat field
                ptsrc_seed *= self.flatfield

                # Add to the list of sources (even though the list will
                # always have only one item)
                tso_seeds.append(copy.deepcopy(ptsrc_seed))
                tso_segs.append(copy.deepcopy(ptsrc_seg))

                # Under the assumption that there will always be only one
                # TSO source, let's assume that the dataset number in the
                # lightcurve file is always 1. This seems easiest for the
                # users when creating the file.
                lightcurve = tso.read_lightcurve(source['lightcurve_file'], 1)
                tso_lightcurves.append(lightcurve)

            if not split_seed:
                # Add the TSO sources to the seed image and segmentation map
                seed, segmap = tso.add_tso_sources(self.seedimage, self.seed_segmap, tso_seeds, tso_segs,
                                                   tso_lightcurves, self.frametime, self.total_frames,
                                                   self.total_frames, self.frames_per_integration,
                                                   self.params['Readout']['nint'],
                                                   self.params['Readout']['resets_bet_ints'],
                                                   samples_per_frametime=5)

                self.segment_number = 1
                self.segment_part_number = 1
                self.segment_frames = seed.shape[1]
                self.segment_ints = seed.shape[0]
                self.segment_frame_start_number = 0
                self.segment_int_start_number = 0
                self.part_int_start_number = 0
                self.part_frame_start_number = 0

                # Zero out the reference pixels
                seed *= self.maskimage
                segmap *= self.maskimage

                # Save the seed image segment to a file
                if self.total_seed_segments_and_parts == 1:
                    self.seed_file = os.path.join(self.basename + '_' + self.params['Readout'][self.usefilt] +
                                                  '_seed_image.fits')
                else:
                    raise ValueError("ERROR: TSO seed file should not be split at this point.")
                    seg_string = str(self.segment_number).zfill(3)
                    part_string = str(self.segment_part_number).zfill(3)
                    self.seed_file = os.path.join(self.basename + '_' + self.params['Readout'][self.usefilt] +
                                                  '_seg{}_part{}_seed_image.fits'.format(seg_string, part_string))

                self.seed_files.append(self.seed_file)
                self.saveSeedImage(seed, segmap, self.seed_file)

            else:
                split_meta = SplitFileMetaData(integration_splits, frame_splits,
                                               self.file_segment_indexes, self.group_segment_indexes_g,
                                               self.frames_per_integration, self.frames_per_group, self.frametime)

                counter = 0
                i = 1
                for int_start in integration_splits[:-1]:
                    int_end = integration_splits[i]

                    j = 1
                    previous_frame = np.zeros(self.seedimage.shape)
                    for initial_frame in frame_splits[:-1]:
                        total_frames = split_meta.total_frames[counter]
                        total_ints = split_meta.total_ints[counter]
                        time_start = split_meta.time_start[counter]
                        frame_start = split_meta.frame_start[counter]
                        seed, segmap = tso.add_tso_sources(self.seedimage, self.seed_segmap,
                                                                           tso_seeds, tso_segs,
                                                                           tso_lightcurves, self.frametime,
                                                                           total_frames, self.total_frames,
                                                                           self.frames_per_integration,
                                                                           total_ints,
                                                                           self.params['Readout']['resets_bet_ints'],
                                                                           starting_time=time_start,
                                                                           starting_frame=frame_start,
                                                                           samples_per_frametime=5)
                        seed += previous_frame
                        previous_frame = seed[-1, -1, :, :]

                        # Zero out the reference pixels
                        seed *= self.maskimage
                        segmap *= self.maskimage

                        # Get metadata values to save in seed image header
                        self.segment_number = split_meta.segment_number[counter]
                        self.segment_ints = split_meta.segment_ints[counter]
                        self.segment_frames = split_meta.segment_frames[counter]
                        self.segment_part_number = split_meta.segment_part_number[counter]
                        self.segment_frame_start_number = split_meta.segment_frame_start_number[counter]
                        self.segment_int_start_number = split_meta.segment_int_start_number[counter]
                        self.part_int_start_number = split_meta.part_int_start_number[counter]
                        self.part_frame_start_number = split_meta.part_frame_start_number[counter]
                        counter += 1

                        # Save the seed image segment to a file
                        print('\n\n\ntotal_seed_segments-and_parts: ', self.total_seed_segments_and_parts)
                        seg_string = str(self.segment_number).zfill(3)
                        part_string = str(self.segment_part_number).zfill(3)
                        self.seed_file = os.path.join(self.basename + '_' + self.params['Readout'][self.usefilt] +
                                                      '_seg{}_part{}_seed_image.fits'.format(seg_string, part_string))

                        self.saveSeedImage(seed, segmap, self.seed_file)
                        self.seed_files.append(self.seed_file)
                        print('Adding file {} to self.seed_files\n'.format(self.seed_file))
                        j += 1

                    i += 1

            self.seedimage = seed
            self.seed_segmap = segmap
        else:
            raise ValueError('Grism TSO simulations should be done with grism_tso_simulator.py!')

    def get_frame_count_info(self):
        """Calculate information on the number of frames per group and
        per integration
        """
        numints = self.params['Readout']['nint']
        numgroups = self.params['Readout']['ngroup']
        numframes = self.params['Readout']['nframe']
        numskips = self.params['Readout']['nskip']
        numresets = self.params['Readout']['resets_bet_ints']

        self.frames_per_group = numframes + numskips
        self.frames_per_integration = numgroups * self.frames_per_group
        self.total_frames = numgroups * self.frames_per_group

        if numints > 1:
            # Frames for all integrations
            self.total_frames *= numints
            # Add the resets for all but the first integration
            self.total_frames += (numresets * (numints - 1))

    def extract_full_from_pom(self, seedimage, seed_segmap):
        """ Given the seed image and segmentation images for the NIRISS POM field of view,
        extract the central 2048x2048 pixel area where the detector sits.  The routine is only
        called when the mode is "pom".  The mode is set to "imaging" in these routine, as after
        the routine is called the subsequent results are the same as when the mode is set to
        "imaging" in the parameter file.

        Parameters:
        -----------

        seedimage : numpy.ndarray (float)   dimension 2322x2322 pixels
        seed_segmap : numpy.ndarray (int)    dimension 2322x2322 pixels

        Returns:
        ---------

        newseedimage : numpy.ndarray (float)  dimension 2048x2048
        newseed_segmap : numpy.ndarray (int)   dimension 2048x2048

        """
        # For the NIRISS POM mode, extact the central 2048x2048 pixels for the
        # ramp simulation.  Set the mode back to "imaging".
        newseedimage = np.copy(seedimage[self.coord_adjust['yoffset']:self.coord_adjust['yoffset']+2048,
                                         self.coord_adjust['xoffset']:self.coord_adjust['xoffset']+2048])
        newseed_segmap = np.copy(seed_segmap[self.coord_adjust['yoffset']:self.coord_adjust['yoffset']+2048,
                                             self.coord_adjust['xoffset']:self.coord_adjust['xoffset']+2048])
        self.params['Inst']['mode'] = "imaging"
        return newseedimage, newseed_segmap

    def add_detector_to_zeropoints(self, detector):
        """Manually add detector dependence to the zeropoint table for
        NIRCam and NIRISS simualtions. This is being done as a placeholder
        for the future, where we expect zeropoints to be detector-dependent.

        Parameters:
        -----------
        detector : str
            Name of detector to add to the table

        Returns:
        --------
        base_table : astropy.table.table
            Modified table of zeropoint info
        """
        # Add "Detector" to the list of column names
        base_table = copy.deepcopy(self.zps)
        num_entries = len(self.zps)
        det_column = Column(np.repeat(detector, num_entries), name="Detector")
        base_table.add_column(det_column, index=0)
        return base_table

    def basic_get_image(self, filename):
        """
        Read in image from a fits file

         Parameters:
        -----------
        filename : str
            Name of fits file to be read in

         Returns:
        --------
        data : obj
            numpy array of data within file

        header : obj
            Header from 0th extension of data file
        """
        data, header = fits.getdata(filename, header=True)
        if len(data.shape) != 2:
            data, header = fits.getdata(filename, 1)
        return data, header

    def prepare_flat(self):
        """
        Read in and prepare the flat field file, to be used
        during the seed creation process.

        Parameters:
        -----------
        None

        Returns:
        --------
        None
        """
        fname = self.params['Reffiles']['pixelflat']

        # Read in file
        try:
            flat, header = fits.getdata(fname, header=True)
        except:
            raise IOError('WARNING: unable to read in {}'.format(fname))

        if (self.params['Output']['grism_source_image']) or (self.params['Inst']['mode'] in ["pom", "wfss"]):
            # If the simulation is for WFSS or POM data, make sure the
            # reference pixels around the flat are also set to
            # non-zero, otherwise you will get a 4 pixel wide picture frame of 0's
            # in the expanded seed image
            bottom = flat[4, :]
            top = flat[2043, :]
            left = flat[:, 4]
            right = flat[:, 2043]
            for i in range(4):
                flat[i, :] = bottom
                flat[i+2044, :] = top
                flat[:, i] = left
                flat[:, i+2044] = right
            flat[0:4, 0:4] = flat[4, 4]
            flat[2044:, 0:4] = flat[2043, 4]
            flat[0:4, 2044:] = flat[4, 2043]
            flat[2044:, 2044:] = flat[2043, 2043]

        # Crop to expected subarray
        try:
            flat = flat[self.subarray_bounds[1]:self.subarray_bounds[3]+1,
                        self.subarray_bounds[0]:self.subarray_bounds[2]+1]
        except:
            raise ValueError("Unable to crop flat field to expected subarray.")

        # If we are making a NIRISS grism direct image, we need to embed the true flat field
        # in an array of the appropriate dimension, where any pixels outside the
        # actual aperture are set to 1.0
        if ((self.params['Output']['grism_source_image']) or (self.params['Inst']['mode'] in ["pom", "wfss"])) and self.instrument == 'niriss':
            mapshape = flat.shape
            # cannot use this: g, yd, xd = signalramp.shape
            # need to update dimensions: self.flat = np.ones((yd, xd))
            self.flatfield = np.ones(self.output_dims)
            ys = self.coord_adjust['yoffset']
            xs = self.coord_adjust['xoffset']
            self.flatfield[ys:ys+mapshape[0], xs:xs+mapshape[1]] = np.copy(flat)
        else:
            # In this case, we have either 1) imaging mode data, or
            # 2) nircam wfss data. In the former case, the nominal flat
            # is all we need. In the latter, the flat is not applied
            # until after dispersing, at which point the seed image is
            # back to the nominal detector/aperture shape.
            self.flatfield = flat

    def pad_wfss_subarray(self, seed, seg):
        """
        WFSS seed images that are to go into the disperser
        must be full frame (or larger). The disperser cannot
        work on subarray images. So embed the subarray seed image
        in a full-frame sized array.

        Parameters:
        -----------
        None

        Returns:
        --------
        Padded seed image and segmentation map
        """
        seeddim = seed.shape
        nx = np.int(2048 * self.grism_direct_factor)
        ffextra = np.int((nx - 2048) / 2)
        bounds = np.array(self.subarray_bounds)

        subextrax = np.int((seeddim[-1] - bounds[2]) / 2)
        subextray = np.int((seeddim[-2] - bounds[3]) / 2)

        extradiffy = ffextra - subextray
        extradiffx = ffextra - subextrax
        exbounds = [extradiffx, extradiffy, extradiffx+seeddim[-1]-1, extradiffy+seeddim[-2]-1]

        if len(seeddim) == 2:
            padded_seed = np.zeros((nx, nx))
            padded_seed[exbounds[1]:exbounds[3] + 1, exbounds[0]:exbounds[2] + 1] = seed
            padded_seg = np.zeros((nx, nx), dtype=np.int)
            padded_seg[exbounds[1]:exbounds[3] + 1, exbounds[0]:exbounds[2] + 1] = seg
        elif len(seeddim) == 4:
            padded_seed = np.zeros((seeddim[0], seeddim[1], nx, nx))
            padded_seed[:, :, exbounds[1]:exbounds[3] + 1, exbounds[0]:exbounds[2] + 1] = seed
            padded_seg = np.zeros((seeddim[0], seeddim[1], nx, nx), dtype=np.int)
            padded_seg[:, :, exbounds[1]:exbounds[3] + 1, exbounds[0]:exbounds[2] + 1] = seg
        else:
            raise ValueError("Seed image is not 2D or 4D. It should be.")
        return padded_seed, padded_seg

    def saveSeedImage(self, seed_image, segmentation_map, seed_file_name):
        """Save the seed image and accompanying segmentation map to a
        fits file.

        Parameters
        ----------
        seed_image : numpy.ndarray
            Array containing the seed image

        segmentation_map : numpy.ndimage
            Array containing the segmentation map

        seed_file_name : str
            Name of FITS file to save ``seed_image`` and `segmentation_map``
            into
        """
        arrayshape = seed_image.shape
        if len(arrayshape) == 2:
            units = 'ADU/sec'
            yd, xd = arrayshape
            grps = 0
            integ = 0
            tgroup = 0.
            arraygroup = 0
            arrayint = 0
            print('Seed image is 2D.')
        elif len(arrayshape) == 3:
            units = 'ADU'
            grps, yd, xd = arrayshape
            integ = 0
            tgroup = self.frametime * (self.params['Readout']['nframe'] + self.params['Readout']['nskip'])
            print('Seed image is 3D.')
        elif len(arrayshape) == 4:
            units = 'ADU'
            integ, grps, yd, xd = arrayshape
            tgroup = self.frametime * (self.params['Readout']['nframe'] + self.params['Readout']['nskip'])
            print('Seed image is 4D.')

        xcent_fov = xd / 2
        ycent_fov = yd / 2

        kw = {}
        kw['XCENTER'] = xcent_fov
        kw['YCENTER'] = ycent_fov
        kw['UNITS'] = units
        kw['TGROUP'] = tgroup

        # Set FGS filter to "N/A" in the output file
        # as this is the value DMS looks for.
        if self.params['Readout'][self.usefilt] == "NA":
            self.params['Readout'][self.usefilt] = "N/A"
        kw['FILTER'] = self.params['Readout'][self.usefilt]
        kw['PHOTFLAM'] = self.photflam
        kw['PHOTFNU'] = self.photfnu
        kw['PHOTPLAM'] = self.pivot * 1.e4  # put into angstroms
        kw['NOMXDIM'] = self.nominal_dims[1]
        kw['NOMYDIM'] = self.nominal_dims[0]
        kw['NOMXSTRT'] = np.int(self.coord_adjust['xoffset'] + 1)
        kw['NOMXEND'] = np.int(self.nominal_dims[1] + self.coord_adjust['xoffset'])
        kw['NOMYSTRT'] = np.int(self.coord_adjust['yoffset'] + 1)
        kw['NOMYEND'] = np.int(self.nominal_dims[0] + self.coord_adjust['yoffset'])

        # Files/inputs used during seed image production
        kw['YAMLFILE'] = self.paramfile
        kw['GAINFILE'] = self.params['Reffiles']['gain']
        kw['DISTORTN'] = self.params['Reffiles']['astrometric']
        kw['IPC'] = self.params['Reffiles']['ipc']
        kw['PIXARMAP'] = self.params['Reffiles']['pixelAreaMap']
        kw['CROSSTLK'] = self.params['Reffiles']['crosstalk']
        kw['FLUX_CAL'] = self.params['Reffiles']['flux_cal']
        kw['FTHRUPUT'] = self.params['Reffiles']['filter_throughput']
        kw['PTSRCCAT'] = self.params['simSignals']['pointsource']
        kw['GALAXCAT'] = self.params['simSignals']['galaxyListFile']
        kw['EXTNDCAT'] = self.params['simSignals']['extended']
        kw['MTPTSCAT'] = self.params['simSignals']['movingTargetList']
        kw['MTSERSIC'] = self.params['simSignals']['movingTargetSersic']
        kw['MTEXTEND'] = self.params['simSignals']['movingTargetExtended']
        kw['NONSDRAL'] = self.params['simSignals']['movingTargetToTrack']
        kw['BKGDRATE'] = self.params['simSignals']['bkgdrate']
        kw['TRACKING'] = self.params['Telescope']['tracking']
        kw['POISSON'] = self.params['simSignals']['poissonseed']
        kw['PSFWFE'] = self.params['simSignals']['psfwfe']
        kw['PSFWFGRP'] = self.params['simSignals']['psfwfegroup']
        kw['MRGEVRSN'] = mirage.__version__

        # Observations with high data volumes (e.g. moving targets, TSO)
        # can be split into multiple "segments" in order to cap the
        # maximum file size
        kw['EXSEGTOT'] = self.total_seed_segments  # Total number of segments in exp
        kw['EXSEGNUM'] = self.segment_number       # Segment number of this file
        kw['PART_NUM'] = self.segment_part_number  # Segment part number of this file

        # Total number of integrations and groups in the current segment
        # (combining all parts of the segment)
        kw['SEGINT'] = self.segment_ints
        kw['SEGGROUP'] = self.segment_frames

        # Total number of integrations and groups in the exposure
        kw['EXPINT'] = self.params['Readout']['nint']
        kw['EXPGRP'] = self.params['Readout']['ngroup']

        # Indexes of the ints and groups where the data in this file go
        # Frame and integration indexes of the segment within the exposure
        kw['SEGFRMST'] = self.segment_frame_start_number
        kw['SEGFRMED'] = self.segment_frame_start_number + grps - 1
        kw['SEGINTST'] = self.segment_int_start_number
        kw['SEGINTED'] = self.segment_int_start_number + integ - 1

        # Frame and integration indexes of the part within the segment
        kw['PTINTSRT'] = self.part_int_start_number
        kw['PTFRMSRT'] = self.part_frame_start_number

        # Seed images provided to disperser are always embedded in an array
        # with dimensions equal to full frame * self.grism_direct_factor
        if self.params['Inst']['mode'] in ['wfss', 'ts_grism']:
            kw['NOMXDIM'] = self.ffsize
            kw['NOMYDIM'] = self.ffsize
            kw['NOMXSTRT'] = np.int(self.ffsize * (self.grism_direct_factor - 1) / 2.)
            kw['NOMXEND'] = kw['NOMXSTRT'] + self.ffsize - 1
            kw['NOMYSTRT'] = np.int(self.ffsize * (self.grism_direct_factor - 1) / 2.)
            kw['NOMYEND'] = kw['NOMYSTRT'] + self.ffsize - 1

        # TSO-specific header keywords
        if self.params['Inst']['mode'] in ['ts_imaging', 'ts_grism']:
            kw['TSOVISIT'] = True
        else:
            kw['TSOVISIT'] = False

        kw['XREF_SCI'] = self.siaf.XSciRef
        kw['YREF_SCI'] = self.siaf.YSciRef

        kw['GRISMPAD'] = self.grism_direct_factor
        self.seedinfo = kw
        self.saveSingleFits(seed_image, seed_file_name, key_dict=kw, image2=segmentation_map, image2type='SEGMAP')

    def combineSimulatedDataSources(self, inputtype, input1, mov_tar_ramp):
        """Combine the exposure containing the trailed sources with the
        countrate image containing the static sources
        inputtype can be 'countrate' in which case input needs to be made
        into a ramp before combining with mov_tar_ramp, or 'ramp' in which
        case you can combine directly. Use 'ramp' with
        non-sidereal TRACKING data, and 'countrate' with sidereal TRACKING data
        """
        if inputtype == 'countrate':
            # First change the countrate image into a ramp
            yd, xd = input1.shape
            numints = self.params['Readout']['nint']
            num_frames = self.params['Readout']['ngroup'] * \
                (self.params['Readout']['nframe'] + self.params['Readout']['nskip'])
            print("Countrate image of synthetic signals being converted to "
                  "RAPID/NISRAPID integration with {} frames.".format(num_frames))
            input1_ramp = np.zeros((numints, num_frames, yd, xd))
            for i in range(num_frames):
                input1_ramp[0, i, :, :] = input1 * self.frametime * (i + 1)
            if numints > 1:
                for integ in range(1, numints):
                    input1_ramp[integ, :, :, :] = input1_ramp[0, :, :, :]

        else:
            # If input1 is a ramp rather than a countrate image
            input1_ramp = input1

        # Combine the input1 ramp and the moving target ramp, which are
        # now both RAPID mode
        totalinput = input1_ramp + mov_tar_ramp
        return totalinput

    def make_trailed_ramp(self):
        # Create a ramp for objects that are trailing through
        # the field of view during the integration
        mov_targs_ramps = []
        mov_targs_segmap = None

        if self.params['Telescope']['tracking'].lower() != 'non-sidereal':
            tracking = False
            ra_vel = None
            dec_vel = None
        else:
            tracking = True
            ra_vel = self.ra_vel
            dec_vel = self.dec_vel
            # print("Moving target mode, creating trailed object with")
            # print("RA velocity of {} and dec_val of {}".format(ra_vel, dec_vel))

        if self.runStep['movingTargets']:
            print("Adding moving point sources to seed image.")
            mov_targs_ptsrc, mt_ptsrc_segmap = self.movingTargetInputs(self.params['simSignals']['movingTargetList'],
                                                                       'pointSource',
                                                                       MT_tracking=tracking,
                                                                       tracking_ra_vel=ra_vel,
                                                                       tracking_dec_vel=dec_vel)
            # Multiply by flat field since these sources are trailed across detector
            if self.instrument == 'nircam' and self.params['Inst']['mode'] in ['wfss', 'ts_grism']:
                pass
            else:
                mov_targs_ptsrc *= self.flatfield
            mov_targs_ramps.append(mov_targs_ptsrc)
            mov_targs_segmap = np.copy(mt_ptsrc_segmap)

        # moving target using a sersic object
        if self.runStep['movingTargetsSersic']:
            print("Adding moving galaxies to seed image.")
            mov_targs_sersic, mt_galaxy_segmap = self.movingTargetInputs(self.params['simSignals']['movingTargetSersic'],
                                                                         'galaxies',
                                                                         MT_tracking=tracking,
                                                                         tracking_ra_vel=ra_vel,
                                                                         tracking_dec_vel=dec_vel)
            # Multiply by flat field
            if self.instrument == 'nircam' and self.params['Inst']['mode'] in ['wfss', 'ts_grism']:
                pass
            else:
                mov_targs_sersic *= self.flatfield
            mov_targs_ramps.append(mov_targs_sersic)
            if mov_targs_segmap is None:
                mov_targs_segmap = np.copy(mt_galaxy_segmap)
            else:
                mov_targs_segmap += mt_galaxy_segmap

        # moving target using an extended object
        if self.runStep['movingTargetsExtended']:
            print("Adding moving extended sources to seed image.")
            mov_targs_ext, mt_ext_segmap = self.movingTargetInputs(self.params['simSignals']['movingTargetExtended'],
                                                                   'extended',
                                                                   MT_tracking=tracking,
                                                                   tracking_ra_vel=ra_vel,
                                                                   tracking_dec_vel=dec_vel)
            # Multiply by flat field
            if self.instrument == 'nircam' and self.params['Inst']['mode'] in ['wfss', 'ts_grism']:
                pass
            else:
                mov_targs_ext *= self.flatfield
            mov_targs_ramps.append(mov_targs_ext)
            if mov_targs_segmap is None:
                mov_targs_segmap = np.copy(mt_ext_segmap)
            else:
                mov_targs_segmap += mt_ext_segmap

        mov_targs_integration = None
        if self.runStep['movingTargets'] or self.runStep['movingTargetsSersic'] or self.runStep['movingTargetsExtended']:
            # Combine the ramps of the moving targets if there is more than one type
            mov_targs_integration = mov_targs_ramps[0]
            if len(mov_targs_ramps) > 1:
                for i in range(1, len(mov_targs_ramps)):
                    mov_targs_integration += mov_targs_ramps[i]
        return mov_targs_integration, mov_targs_segmap

    def calcCoordAdjust(self):
        # Calculate the factors by which to expand the output array size, as well as the coordinate
        # offsets between the nominal output array and the input lists if the observation being
        # modeled is wfss

        dtor = math.radians(1.)
        instrument = self.params['Inst']['instrument']

        # Normal imaging with grism image requested
        if (instrument.lower() == 'nircam' and self.params['Output']['grism_source_image']) or \
           (instrument.lower() == 'niriss' and (self.params['Inst']['mode'] in ["pom", "wfss"] or self.params['Output']['grism_source_image'])):
            self.coord_adjust['x'] = self.grism_direct_factor
            self.coord_adjust['y'] = self.grism_direct_factor
            self.coord_adjust['xoffset'] = np.int((self.grism_direct_factor - 1.) *
                                                  (self.subarray_bounds[2] -
                                                   self.subarray_bounds[0] + 1) / 2.)
            self.coord_adjust['yoffset'] = np.int((self.grism_direct_factor - 1.) *
                                                  (self.subarray_bounds[3] -
                                                   self.subarray_bounds[1] + 1) / 2.)

    def non_sidereal_seed(self):
        """
        Create a seed EXPOSURE in the case where the instrument is tracking
        a non-sidereal target
        """

        # Create a count rate image containing only the non-sidereal target(s)
        # These will be stationary in the fov
        nonsidereal_countrate, nonsidereal_segmap, self.ra_vel, self.dec_vel, vel_flag \
            = self.nonsidereal_CRImage(self.params['simSignals']['movingTargetToTrack'])

        # Expand into a RAPID exposure and convert from signal rate to signals
        ns_yd, ns_xd = nonsidereal_countrate.shape
        ns_int = self.params['Readout']['nint']
        ns_group = self.params['Readout']['ngroup']
        ns_nframe = self.params['Readout']['nframe']
        ns_nskip = self.params['Readout']['nskip']
        totframes = ns_group * (ns_nframe + ns_nskip)
        tmptimes = self.frametime * np.arange(1, totframes + 1)

        non_sidereal_ramp = np.zeros((ns_int, totframes, ns_yd, ns_xd))
        for i in range(totframes):
            for integ in range(ns_int):
                non_sidereal_ramp[integ, i, :, :] = nonsidereal_countrate * tmptimes[i]

        # Now we need to collect all the other sources (point sources,
        # galaxies, extended) in the other input files, and treat them
        # as targets which will move across the field of view during
        # the exposure.
        mtt_data_list = []
        mtt_data_segmap = None

        if self.runStep['pointsource']:
            # Now ptsrc is a list, which we need to provide to
            # movingTargetInputs
            print("Adding moving background point sources to seed image.")
            mtt_ptsrc, mtt_ptsrc_segmap = self.movingTargetInputs(self.params['simSignals']['pointsource'],
                                                                  'pointSource',
                                                                  MT_tracking=True,
                                                                  tracking_ra_vel=self.ra_vel,
                                                                  tracking_dec_vel=self.dec_vel,
                                                                  trackingPixVelFlag=vel_flag)
            # Multiply by flat field since these sources are trailed
            # across detector
            if self.instrument == 'nircam' and self.params['Inst']['mode'] in ['wfss', 'ts_grism']:
                pass
            else:
                mtt_ptsrc *= self.flatfield
            mtt_data_list.append(mtt_ptsrc)
            if mtt_data_segmap is None:
                mtt_data_segmap = np.copy(mtt_ptsrc_segmap)
            else:
                mtt_data_segmap += mtt_ptsrc_segmap
            print(("Done with creating moving targets from {}"
                   .format(self.params['simSignals']['pointsource'])))

        if self.runStep['galaxies']:
            print("Adding moving background galaxies to seed image.")
            mtt_galaxies, mtt_galaxies_segmap = self.movingTargetInputs(self.params['simSignals']['galaxyListFile'],
                                                                        'galaxies',
                                                                        MT_tracking=True,
                                                                        tracking_ra_vel=self.ra_vel,
                                                                        tracking_dec_vel=self.dec_vel,
                                                                        trackingPixVelFlag=vel_flag)
            # Multiply by flat field
            if self.instrument == 'nircam' and self.params['Inst']['mode'] in ['wfss', 'ts_grism']:
                pass
            else:
                mtt_galaxies *= self.flatfield
            mtt_data_list.append(mtt_galaxies)
            if mtt_data_segmap is None:
                mtt_data_segmap = np.copy(mtt_galaxies_segmap)
            else:
                mtt_data_segmap += mtt_galaxies_segmap

            print(("Done with creating moving targets from {}".
                   format(self.params['simSignals']['galaxyListFile'])))

        if self.runStep['extendedsource']:
            print("Adding moving background extended sources to seed image.")
            mtt_ext, mtt_ext_segmap = self.movingTargetInputs(self.params['simSignals']['extended'],
                                                              'extended',
                                                              MT_tracking=True,
                                                              tracking_ra_vel=self.ra_vel,
                                                              tracking_dec_vel=self.dec_vel,
                                                              trackingPixVelFlag=vel_flag)
            # Multiply by flat field
            if self.instrument == 'nircam' and self.params['Inst']['mode'] in ['wfss', 'ts_grism']:
                pass
            else:
                mtt_ext *= self.flatfield
            mtt_data_list.append(mtt_ext)
            if mtt_data_segmap is None:
                mtt_data_segmap = np.copy(mtt_ext_segmap)
            else:
                mtt_data_segmap += mtt_ext_segmap

            print(("Done with creating moving targets from {}".
                   format(self.params['simSignals']['extended'])))

        # Add in the other objects which are not being tracked on
        # (i.e. the sidereal targets)
        if len(mtt_data_list) > 0:
            for i in range(len(mtt_data_list)):
                non_sidereal_ramp += mtt_data_list[i]
                # non_sidereal_zero += mtt_zero_list[i]
        if mtt_data_segmap is not None:
            nonsidereal_segmap.segmap += mtt_data_segmap
        return non_sidereal_ramp, nonsidereal_segmap.segmap

    def readMTFile(self, filename):
        """
        Read in moving target list file

        Arguments:
        ----------
        filename : str
            name of moving target catalog file

        Returns:
        --------
        returns : obj
            Table containing moving target entries
            pixelflag (boolean) -- If true, locations are in units of
                pixels. If false, locations are RA, Dec
            pixelvelflag (boolean) -- If true, moving target velocities
                are in units of pixels/hour. If false, arcsec/hour
            magsys -- magnitude system of the moving target magnitudes
        """
        mtlist = ascii.read(filename, comment='#')

        # Convert all relevant columns to floats
        for col in mtlist.colnames:
            if mtlist[col].dtype in ['int64', 'int']:
                mtlist[col] = mtlist[col].data * 1.

        # Check to see whether the position is in x,y or ra,dec
        pixelflag = False
        try:
            if 'position_pixels' in mtlist.meta['comments'][0:4]:
                pixelflag = True
        except:
            pass

        # If present, check whether the velocity entries are pix/sec
        # or arcsec/sec.
        pixelvelflag = False
        try:
            if 'velocity_pixels' in mtlist.meta['comments'][0:4]:
                pixelvelflag = True
        except:
            pass

        # If present, check whether the radius entries (for galaxies)
        # are in arcsec or pixels. If in arcsec, change to pixels
        if 'radius' in mtlist.colnames:
            if 'radius_pixels' not in mtlist.meta['comments'][0:4]:
                mtlist['radius'] /= self.siaf.XSciScale

        # If galaxies are present, change position angle from degrees
        # to radians
        if 'pos_angle' in mtlist.colnames:
            mtlist['pos_angle'] = mtlist['pos_angle'] * np.pi / 180.

        # Check to see if magnitude system is specified in comments
        # If not, assume AB mags
        msys = 'abmag'

        condition = ('stmag' in mtlist.meta['comments'][0:4]) | ('vegamag' in mtlist.meta['comments'][0:4])
        if condition:
            msys = [l for l in mtlist.meta['comments'][0:4] if 'mag' in l][0]
            msys = msys.lower()

        return mtlist, pixelflag, pixelvelflag, msys.lower()

    def get_index_numbers(self, catalog_table):
        """Get index numbers associated with the sources in a catalog

        Parameters
        ----------
        catalog_table : astropy.table.Table
            Source catalog read in from an ascii file

        Returns
        -------
        indexes : list
            List of index numbers corresponding to sources in the catalog
        """
        # If the input catalog has an index column
        # use that, otherwise add one
        if 'index' in catalog_table.colnames:
            indexes = catalog_table['index']
        else:
            indexes = np.arange(1, len(catalog_table['x_or_RA']) + 1)
        # Make sure there is no 0th object
        if np.min(indexes) == 0:
            indexes += 1
        # Make sure the index numbers don't overlap with any
        # sources already present. Increment the maxindex
        # value.
        if np.min(indexes) <= self.maxindex:
            indexes += self.maxindex
        self.maxindex = np.max(indexes)
        return indexes

    def movingTargetInputs(self, filename, input_type, MT_tracking=False,
                           tracking_ra_vel=None, tracking_dec_vel=None,
                           trackingPixVelFlag=False):
        """Read in listfile of moving targets and perform needed
        calculations to get inputs for moving_targets.py

        Parameters
        ----------

        filename : str
            Name of catalog contining moving target sources

        input_type : str
            Specifies type of sources. Can be 'pointSource','galaxies', or 'extended'

        MT_tracking : bool
            If True, observation is non-sidereal (i.e. telescope is tracking the moving target)

        tracking_ra_vel : float
            Velocity of the moving target in the detector x or right ascension direction.
            Units are pixels/hour or arcsec/hour depending on trackingPixVelFlag.

        tracking_dec_vel : float
            Velocity of moving target in the detector y or declination direction.
            Units are pixels/hour or arcsec/hour depending on trackingPixVelFlag.

        trackingPixVelFlag : bool
            If True, tracking_ra_vel and tracking_dec_vel are in units of pixels/hour, and
            velocities are in the detector x and y directions, respectively.
            If False, velocity untis are arcsec/hour and directions are RA, and Dec.

        Returns
        -------

        mt_integration : numpy.ndarray
            4D array containing moving target seed image

        moving_segmap.segmap : numpy.ndarray
            2D array containing segmentation map that goes with the seed image.
            Segmentation map is based on the final frame in the seed image.
        """
        # Read input file - should be able to use for all modes
        mtlist, pixelFlag, pixvelflag, magsys = self.readMTFile(filename)

        # If the input catalog has an index column
        # use that, otherwise add one
        indexes = self.get_index_numbers(mtlist)

        if MT_tracking is True:
            # Here, we are tracking a non-sidereal target.
            try:
                # If there are moving targets on top of the non-
                # sidereal tracking (e.g. tracking Io but Europa
                # comes into the fov), then the velocity vector
                # of the moving target needs to be adjusted.
                # If the input catalog already contains
                # 'x_or_RA_velocity' then we know we have a moving
                # target. If it doesn't, then we have sidereal
                # targets, and we can simply set their velocity
                # as the inverse of that being tracked.
                mtlist['x_or_RA_velocity'] -= tracking_ra_vel
                mtlist['y_or_Dec_velocity'] -= tracking_dec_vel
                pixvelflag = trackingPixVelFlag
            except:
                print('Setting velocity of targets equal to the non-sidereal tracking velocity')
                mtlist['x_or_RA_velocity'] = 0. - tracking_ra_vel
                mtlist['y_or_Dec_velocity'] = 0. - tracking_dec_vel
                pixvelflag = trackingPixVelFlag

        # Exposure times for all frames
        numints = self.params['Readout']['nint']
        numgroups = self.params['Readout']['ngroup']
        numframes = self.params['Readout']['nframe']
        numskips = self.params['Readout']['nskip']
        numresets = self.params['Readout']['resets_bet_ints']

        frames_per_group = numframes + numskips
        frames_per_integration = numgroups * frames_per_group
        total_frames = numgroups * frames_per_group
        # If only one integration per exposure, then total_frames
        # above is correct. For >1 integration, we need to add the reset
        # frame to each integration (except the first), and sum the number of
        # frames for all integrations

        if numints > 1:
            # Frames for all integrations
            total_frames *= numints
            # Add the resets for all but the first and last integrations
            total_frames += (numresets * (numints - 1))

        frameexptimes = self.frametime * np.arange(-1, total_frames)

        # output image dimensions
        dims = self.nominal_dims
        newdimsx = np.int(dims[1] * self.coord_adjust['x'])
        newdimsy = np.int(dims[0] * self.coord_adjust['y'])

        # Set up seed integration
        #mt_integration = np.zeros((numints, total_frames, newdimsy, newdimsx))
        mt_integration = np.zeros((numints, frames_per_integration, newdimsy, newdimsx))

        # Corresponding (2D) segmentation map
        moving_segmap = segmap.SegMap()
        moving_segmap.xdim = newdimsx
        moving_segmap.ydim = newdimsy
        moving_segmap.initialize_map()

        # Check the source list and remove any sources that are well outside the
        # field of view of the detector. These sources cause the coordinate
        # conversion to hang.
        indexes, mtlist = self.remove_outside_fov_sources(indexes, mtlist, pixelFlag, 4096)

        # Determine the name of the column to use for source magnitudes
        mag_column = self.select_magnitude_column(mtlist, filename)

        times = []
        obj_counter = 0
        time_reported = False
        for index, entry in zip(indexes, mtlist):
            start_time = time.time()
            # For each object, calculate x,y or RA,Dec of initial position
            pixelx, pixely, ra, dec, ra_str, dec_str = self.get_positions(
                entry['x_or_RA'], entry['y_or_Dec'], pixelFlag, 4096)

            # Now generate a list of x,y position in each frame
            if pixvelflag is False:
                # Calculate the RA,Dec in each frame
                # input velocities are arcsec/hour. ra/dec are in units of degrees,
                # so divide velocities by 3600^2.
                ra_frames = ra + (entry['x_or_RA_velocity'] / 3600. / 3600.) * frameexptimes
                dec_frames = dec + (entry['y_or_Dec_velocity'] / 3600. / 3600.) * frameexptimes

                x_frames = []
                y_frames = []
                for in_ra, in_dec in zip(ra_frames, dec_frames):
                    # Calculate the x,y position at each frame
                    px, py, pra, pdec, pra_str, pdec_str = self.get_positions(in_ra, in_dec, False, 4096)
                    x_frames.append(px)
                    y_frames.append(py)
                x_frames = np.array(x_frames)
                y_frames = np.array(y_frames)

            else:
                # If input velocities are pixels/hour, then generate the list of
                # x,y in each frame directly
                x_frames = pixelx + (entry['x_or_RA_velocity'] / 3600.) * frameexptimes
                y_frames = pixely + (entry['y_or_Dec_velocity'] / 3600.) * frameexptimes

            # Get countrate and PSF size info
            if entry[mag_column] is not None:
                rate = utils.magnitude_to_countrate(self.params['Inst']['mode'],
                                                    magsys, entry[mag_column],
                                                    photfnu=self.photfnu,
                                                    photflam=self.photflam)
            else:
                rate = 1.0

            psf_x_dim = self.find_psf_size(rate)
            psf_dimensions = (psf_x_dim, psf_x_dim)
            #psf_dimensions = (self.psf_library_x_dim, self.psf_library_y_dim)

            # If we have a point source, we can easily determine whether
            # it completely misses the detector, since we know the size
            # of the stamp already. For galaxies and extended sources,
            # we have to get the stamp image first to see if any part of
            # the stamp lands on the detector.
            status = 'on'
            if input_type == 'pointSource':

                status = self.on_detector(x_frames, y_frames, psf_dimensions,
                                          (newdimsx, newdimsy))
            if status == 'off':
                continue

            # Create the PSF
            eval_psf, minx, miny, wings_added = self.create_psf_stamp(pixelx, pixely, psf_x_dim, psf_x_dim,
                                                                      ignore_detector=True)

            # Skip sources that fall completely off the detector
            if eval_psf is None:
                continue

            if input_type == 'pointSource':
                stamp = eval_psf
                stamp *= rate

            elif input_type == 'extended':
                stamp, header = self.basic_get_image(entry['filename'])
                stamp = self.rotate_extended_image(stamp, entry['pos_angle'], ra, dec)

                # If no magnitude is given, use the extended image as-is
                if rate != 1.0:
                    stamp /= np.sum(stamp)
                    stamp *= rate

                # Convolve with instrument PSF if requested
                if self.params['simSignals']['PSFConvolveExtended']:
                    stamp_dims = stamp.shape
                    # If the stamp image is smaller than the PSF in either
                    # dimension, embed the stamp in an array that matches
                    # the psf size. This is so the upcoming convolution will
                    # produce an output that includes the wings of the PSF
                    psf_shape = eval_psf.shape
                    if ((stamp_dims[0] < psf_shape[0]) or (stamp_dims[1] < psf_shape[1])):
                        stamp = self.enlarge_stamp(stamp, psf_shape)
                        stamp_dims = stamp.shape

                    # Convolve stamp with PSF
                    stamp = s1.fftconvolve(stamp, eval_psf, mode='same')

            elif input_type == 'galaxies':
                pixelx, pixely, ra, dec, ra_str, dec_str = self.get_positions(entry['x_or_RA'],
                                                                              entry['y_or_Dec'],
                                                                              pixelFlag, 4096)

                pixelv2, pixelv3 = pysiaf.utils.rotations.getv2v3(self.attitude_matrix, ra, dec)

                xposang = self.calc_x_position_angle(pixelv2, pixelv3, entry['pos_angle'])

                # First create the galaxy
                stamp = self.create_galaxy(entry['radius'], entry['ellipticity'], entry['sersic_index'],
                                           xposang*np.pi/180., rate)

                # If the stamp image is smaller than the PSF in either
                # dimension, embed the stamp in an array that matches
                # the psf size. This is so the upcoming convolution will
                # produce an output that includes the wings of the PSF
                galdims = stamp.shape
                psf_shape = eval_psf.shape
                if ((galdims[0] < psf_shape[0]) or (galdims[1] < psf_shape[1])):
                    stamp = self.enlarge_stamp(stamp, psf_shape)
                    galdims = stamp.shape

                # Convolve the galaxy with the instrument PSF
                stamp = s1.fftconvolve(stamp, eval_psf, mode='same')

            # Now that we have stamp images for galaxies and extended
            # sources, check to see if they overlap the detector or not.
            # NOTE: this will only catch sources that never overlap the
            # detector for any of their positions.
            if input_type != 'pointSource':
                status = self.on_detector(x_frames, y_frames, stamp.shape,
                                          (newdimsx, newdimsy))
            if status == 'off':
                continue

            # Each entry will have stamp image as array, ra_init, dec_init,
            # ra_velocity, dec_velocity, frametime, numframes, subsample_factor,
            # outputarrayxsize, outputarrayysize
            # (maybe without the values that will be the same to each entry.

            # Need to feed info into moving_targets one integration at a time.
            # No need to feed in the reset frames, but they are necessary
            # before this point in order to get the timing and positions
            # correct.
            for integ in range(numints):
                framestart = integ * frames_per_integration + integ
                frameend = framestart + frames_per_integration + 1

                # Now check to see if the stamp image overlaps the output
                # aperture for this integration only. Above we removed sources
                # that never overlap the aperture. Here we want to get rid
                # of sources that overlap the detector in some integrations,
                # but not this particular integration
                status = 'on'
                status = self.on_detector(x_frames[framestart:frameend],
                                          y_frames[framestart:frameend],
                                          stamp.shape, (newdimsx, newdimsy))
                if status == 'off':
                    continue

                mt = moving_targets.MovingTarget()
                mt.subsampx = 3
                mt.subsampy = 3
                mt_source = mt.create(stamp, x_frames[framestart:frameend],
                                      y_frames[framestart:frameend],
                                      self.frametime, newdimsx, newdimsy)
                mt_integration[integ, :, :, :] += mt_source

                noiseval = self.single_ron / 100. + self.params['simSignals']['bkgdrate']
                if self.params['Inst']['mode'].lower() in ['wfss', 'ts_grism']:
                    noiseval += self.grism_background

                if input_type in ['pointSource', 'galaxies']:
                    moving_segmap.add_object_noise(mt_source[-1, :, :], 0, 0, index, noiseval)
                else:
                    indseg = self.seg_from_photutils(mt_source[-1, :, :], np.int(index), noiseval)
                    moving_segmap.segmap += indseg

            # Check the elapsed time for creating each object
            elapsed_time = time.time() - start_time
            times.append(elapsed_time)
            if obj_counter > 3 and not time_reported:
                avg_time = np.mean(times)
                total_time = len(indexes) * avg_time
                print(("Expected time to process {} sources: {:.2f} seconds "
                       "({:.2f} minutes)".format(len(indexes), total_time, total_time/60)))
                time_reported = True
            obj_counter += 1
        return mt_integration, moving_segmap.segmap

    def on_detector(self, xloc, yloc, stampdim, finaldim):
        """Given a set of x, y locations, stamp image dimensions,
        and final image dimensions, determine whether the stamp
        image will overlap at all with the final image, or
        completely miss it.

        Parameters
        ---------

        xloc : list
            X-coordinate locations of source

        yloc : list
            Y-coordinate locations of source

        stampdim : tuple
            (x,y) dimension lengths of stamp image

        finaldim : tuple
            (x,y) dimension lengths of final image

        Returns
        -------

        status : str
            state of stamp image:
            "on" -- stamp image fully or partially overlaps the final
                 image for at least one xloc, yloc pair
            "off" -- stamp image never overlaps with final image
        """
        status = 'on'
        stampx, stampy = stampdim
        stampminx = np.min(xloc - (stampx / 2))
        stampminy = np.min(yloc - (stampy / 2))
        stampmaxx = np.max(xloc + (stampx / 2))
        stampmaxy = np.max(yloc + (stampy / 2))
        finalx, finaly = finaldim
        if ((stampminx > finalx) or (stampmaxx < 0) or
           (stampminy > finaly) or (stampmaxy < 0)):
            status = 'off'
        return status

    def get_positions(self, input_x, input_y, pixel_flag, max_source_distance):
        """ Given an input position ( (x,y) or (RA,Dec) ) for a source, calculate
        the corresponding detector (x,y), RA, Dec, and provide RA and Dec strings.

        Parameters
        ----------

        input_x : str
            Detector x coordinate or RA of source. RA can be in decimal degrees or
            (e.g 10:23:34.2 or 10h23m34.2s)

        input_y : str
            Detector y coordinate or Dec of source. Dec can be in decimal degrees
            or (e.g. 10d:23m:34.2s)

        pixel_flag : bool
            True if input_x and input_y are in units of pixels. False if they are
            in the RA, Dec coordinate system.

        max_source_distance : float
            Maximum number of pixels from the aperture's reference location to keep
            a source. Sources very far off the detector will cause the calculation of
            RA, Dec to hang.

        Returns
        -------

        pixelx : float
            Detector x coordinate of source

        pixely : float
            Detector y coordinate of source

        ra : float
            RA of source (degrees)

        dec : float
            Dec of source (degrees)

        ra_string : str
            String representation of RA

        dec_string : str
            String representation of Dec
        """
        try:
            entry0 = float(input_x)
            entry1 = float(input_y)
            if not pixel_flag:
                ra_string, dec_string = self.makePos(entry0, entry1)
                ra_number = entry0
                dec_number = entry1
        except:
            # if inputs can't be converted to floats, then
            # assume we have RA/Dec strings. Convert to floats.
            ra_string = input_x
            dec_string = input_y
            ra_number, dec_number = utils.parse_RA_Dec(ra_string, dec_string)

        # Case where point source list entries are given with RA and Dec
        if not pixel_flag:

            # If distortion is to be included - either with or without the full set of coordinate
            # translation coefficients
            pixel_x, pixel_y = self.RADecToXY_astrometric(ra_number, dec_number)

        else:
            # Case where the point source list entry locations are given in units of pixels
            # In this case we have the source position, and RA/Dec are calculated only so
            # they can be written out into the output source list file.
            pixel_x = entry0
            pixel_y = entry1
            ra_number, dec_number, ra_string, dec_string = self.XYToRADec(pixel_x, pixel_y)

        return pixel_x, pixel_y, ra_number, dec_number, ra_string, dec_string

    def nonsidereal_CRImage(self, file):
        """
        Create countrate image of non-sidereal sources
        that are being tracked.

        Arguments:
        ----------
        file : str
            name of the file containing the tracked moving targets.

        Returns:
        --------
        returns : obj
            countrate image (2D) containing the tracked non-sidereal targets.
        """
        totalCRList = []
        totalSegList = []

        # Read in file containing targets
        targs, pixFlag, velFlag, magsys = self.readMTFile(file)

        # We need to keep track of the proper motion of the
        # target being tracked, because all other objects in
        # the field of view will be moving at the same rate
        # in the opposite direction
        track_ra_vel = targs[0]['x_or_RA_velocity']
        track_dec_vel = targs[0]['y_or_Dec_velocity']

        # Sort the targets by whether they are point sources,
        # galaxies, extended
        ptsrc_rows = []
        galaxy_rows = []
        extended_rows = []
        for i, line in enumerate(targs):
            if 'point' in line['object'].lower():
                ptsrc_rows.append(i)
            elif 'sersic' in line['object'].lower():
                galaxy_rows.append(i)
            else:
                extended_rows.append(i)

        # Re-use functions for the sidereal tracking situation
        if len(ptsrc_rows) > 0:
            ptsrc = targs[ptsrc_rows]
            if pixFlag:
                meta0 = 'position_pixels'
            else:
                meta0 = ''
            if velFlag:
                meta1 = 'velocity_pixels'
            else:
                meta1 = ''
            meta2 = magsys

            meta3 = ('Point sources with non-sidereal tracking. '
                     'File produced by catalog_seed_image.py')
            meta4 = ('from run using non-sidereal moving target '
                     'list {}.'.format(self.params['simSignals']['movingTargetToTrack']))
            ptsrc.meta['comments'] = [meta0, meta1, meta2, meta3, meta4]

            temp_ptsrc_filename = os.path.join(self.params['Output']['directory'],
                                               'temp_non_sidereal_point_sources.list')
            ptsrc.write(temp_ptsrc_filename, format='ascii', overwrite=True)

            ptsrc = self.get_point_source_list(temp_ptsrc_filename)
            ptsrcCRImage, ptsrcCRSegmap = self.make_point_source_image(ptsrc)

            if self.instrument == 'nircam' and self.params['Inst']['mode'] in ['wfss', 'ts_grism']:
                pass
            else:
                ptsrcCRImage *= self.flatfield

            totalCRList.append(ptsrcCRImage)
            totalSegList.append(ptsrcCRSegmap)

        if len(galaxy_rows) > 0:
            galaxies = targs[galaxy_rows]
            if pixFlag:
                meta0 = 'position_pixels'
            else:
                meta0 = ''
            if velFlag:
                meta1 = 'velocity_pixels'
            else:
                meta1 = ''
            meta2 = magsys
            meta3 = ('Galaxies (2d sersic profiles) with non-sidereal '
                     'tracking. File produced by ramp_simulator.py')
            meta4 = ('from run using non-sidereal moving target '
                     'list {}.'.format(self.params['simSignals']['movingTargetToTrack']))
            galaxies.meta['comments'] = [meta0, meta1, meta2, meta3, meta4]
            galaxies.write(os.path.join(self.params['Output']['directory'], 'temp_non_sidereal_sersic_sources.list'), format='ascii', overwrite=True)

            galaxyCRImage, galaxySegmap = self.make_galaxy_image('temp_non_sidereal_sersic_sources.list')

            if self.instrument == 'nircam' and self.params['Inst']['mode'] in ['wfss', 'ts_grism']:
                pass
            else:
                galaxyCRImage *= self.flatfield

            totalCRList.append(galaxyCRImage)
            totalSegList.append(galaxySegmap)

        if len(extended_rows) > 0:
            extended = targs[extended_rows]

            if pixFlag:
                meta0 = 'position_pixels'
            else:
                meta0 = ''
            if velFlag:
                meta1 = 'velocity_pixels'
            else:
                meta1 = ''
            meta2 = magsys
            meta3 = 'Extended sources with non-sidereal tracking. File produced by ramp_simulator.py'
            meta4 = 'from run using non-sidereal moving target list {}.'.format(self.params['simSignals']['movingTargetToTrack'])
            extended.meta['comments'] = [meta0, meta1, meta2, meta3, meta4]
            extended.write(os.path.join(self.params['Output']['directory'], 'temp_non_sidereal_extended_sources.list'), format='ascii', overwrite=True)

            extlist, extstamps = self.getExtendedSourceList('temp_non_sidereal_extended_sources.list')

            # translate the extended source list into an image
            extCRImage, extSegmap = self.make_extended_source_image(extlist, extstamps)

            # Multiply by the flat field
            if self.instrument == 'nircam' and self.params['Inst']['mode'] in ['wfss', 'ts_grism']:
                pass
            else:
                extCRImage *= self.flatfield

            totalCRList.append(extCRImage)
            totalSegList.append(extSegmap)

        # Now combine into a final countrate image of non-sidereal sources (that are being tracked)
        if len(totalCRList) > 0:
            totalCRImage = totalCRList[0]
            totalSegmap = totalSegList[0]
            if len(totalCRList) > 1:
                for i in range(1, len(totalCRList)):
                    totalCRImage += totalCRList[i]
                    totalSegmap += totalSegList[i] + i
        else:
            raise ValueError(("No non-sidereal countrate targets produced."
                              "You shouldn't be here."))
        return totalCRImage, totalSegmap, track_ra_vel, track_dec_vel, velFlag

    def create_sidereal_image(self):
        # Generate a signal rate image from input sources
        if (self.params['Output']['grism_source_image'] == False) and (not self.params['Inst']['mode'] in ["pom", "wfss"]):
            signalimage = np.zeros(self.nominal_dims)
            segmentation_map = np.zeros(self.nominal_dims)
        else:
            xd = np.int(self.nominal_dims[1] * self.coord_adjust['x'])
            yd = np.int(self.nominal_dims[0] * self.coord_adjust['y'])
            signalimage = np.zeros((yd, xd), dtype=np.float)
            segmentation_map = np.zeros((yd, xd))

        instrument_name = self.params['Inst']['instrument'].lower()
        # yd, xd = signalimage.shape
        arrayshape = signalimage.shape

        # POINT SOURCES
        # Read in the list of point sources to add
        # Adjust point source locations using astrometric distortion
        # Translate magnitudes to counts in a single frame
        if self.runStep['pointsource'] is True:
            if not self.expand_catalog_for_segments:

                # Translate the point source list into an image
                print('Calculating point source lists')
                pslist = self.get_point_source_list(self.params['simSignals']['pointsource'])
                psfimage, ptsrc_segmap = self.make_point_source_image(pslist)

            elif self.expand_catalog_for_segments:
                # Expand the point source list for each mirror segment, and add together
                # the 18 point source images that used different PSFs.
                print('Expanding the source catalog for 18 mirror segments')

                # Create empty image and segmentation map
                ptsrc_segmap = segmap.SegMap()
                ptsrc_segmap.ydim, ptsrc_segmap.xdim = self.output_dims
                ptsrc_segmap.initialize_map()
                psfimage = np.zeros(self.output_dims)

                library_list = get_segment_library_list(
                    self.params['Inst']['instrument'].lower(), self.detector, self.psf_filter,
                    self.params['simSignals']['psfpath'], pupil=self.psf_pupil
                )
                for i_segment in np.arange(1, 19):
                    print('\nCalculating point source lists for segment {}'.format(i_segment))
                    # Get the RA/Dec offset that matches the given segment
                    offset_vector = get_segment_offset(i_segment, self.detector, library_list)

                    # need to add a new offset option to get_point_source_list:
                    pslist = self.get_point_source_list(self.params['simSignals']['pointsource'],
                                                        segment_offset=offset_vector)

                    # Create a point source image, using the specific point
                    # source list and PSF for the given segment
                    seg_psfimage, ptsrc_segmap = self.make_point_source_image(pslist, segment_number=i_segment,
                                                                              ptsrc_segmap=ptsrc_segmap)

                    if self.params['Output']['save_intermediates'] is True:
                        seg_psfImageName = self.basename + '_pointSourceRateImage_seg{:02d}.fits'.format(i_segment)
                        h0 = fits.PrimaryHDU(seg_psfimage)
                        h0.writeto(seg_psfImageName, overwrite=True)
                        print("    Segment {} point source image and segmap saved as {}".format(i_segment,
                                                                                                seg_psfImageName))

                    psfimage += seg_psfimage

            # Multiply by the flat field unless it is a NIRCam observation
            # that is going to be dispersed
            if self.instrument == 'nircam' and self.params['Inst']['mode'] in ['wfss', 'ts_grism']:
                pass
            else:
                psfimage *= self.flatfield

            ptsrc_segmap = ptsrc_segmap.segmap

            # Add the point source image to the overall image
            signalimage = signalimage + psfimage
            segmentation_map += ptsrc_segmap

            # To avoid problems with overlapping sources between source
            # types in observations to be dispersed, make the point
            # source-only segmentation map available as a class variable
            self.point_source_seed = psfimage
            self.point_source_seg_map = ptsrc_segmap

            # For seed images to be dispersed in WFSS mode,
            # embed the seed image in a full frame array. The disperser
            # tool does not work on subarrays
            aperture_suffix = self.params['Readout']['array_name'].split('_')[-1]
            if ((self.params['Inst']['mode'] in ['wfss', 'ts_grism']) & \
                 (aperture_suffix not in ['FULL', 'CEN'])):
                self.point_source_seed, self.point_source_seg_map = self.pad_wfss_subarray(self.point_source_seed,
                                                                                           self.point_source_seg_map)

            # Save the point source seed image
            if instrument_name != 'fgs':
                self.ptsrc_seed_filename = os.path.join(self.basename + '_' + self.params['Readout'][self.usefilt] + '_ptsrc_seed_image.fits')
            else:
                self.ptsrc_seed_filename = '{}_ptsrc_seed_image.fits'.format(self.basename)
            self.saveSeedImage(self.point_source_seed, self.point_source_seg_map, self.ptsrc_seed_filename)
            print("Point source image and segmap saved as {}".format(self.ptsrc_seed_filename))

        else:
            self.point_source_seed = None
            self.point_source_seg_map = None
            self.ptsrc_seed_filename = None

        # Simulated galaxies
        # Read in the list of galaxy positions/magnitudes to simulate
        # and create a countrate image of those galaxies.
        if self.runStep['galaxies'] is True:
            galaxyCRImage, galaxy_segmap = self.make_galaxy_image(self.params['simSignals']['galaxyListFile'])

            # Multiply by the flat field unless it is a NIRCam observation
            # that is going to be dispersed
            if self.instrument == 'nircam' and self.params['Inst']['mode'] in ['wfss', 'ts_grism']:
                pass
            else:
                galaxyCRImage *= self.flatfield

            # To avoid problems with overlapping sources between source
            # types in observations to be dispersed, make the galaxy-
            # only segmentation map available as a class variable
            self.galaxy_source_seed = galaxyCRImage
            self.galaxy_source_seg_map = galaxy_segmap

            # For seed images to be dispersed in WFSS mode,
            # embed the seed image in a full frame array. The disperser
            # tool does not work on subarrays
            aperture_suffix = self.params['Readout']['array_name'].split('_')[-1]
            if ((self.params['Inst']['mode'] in ['wfss', 'ts_grism']) & \
                 (aperture_suffix not in ['FULL', 'CEN'])):
                self.galaxy_source_seed, self.galaxy_source_seg_map = self.pad_wfss_subarray(self.galaxy_source_seed,
                                                                                             self.galaxy_source_seg_map)

            # Save the galaxy source seed image
            if instrument_name != 'fgs':
                self.galaxy_seed_filename = os.path.join(self.basename + '_' + self.params['Readout'][self.usefilt] + '_galaxy_seed_image.fits')
            else:
                self.galaxy_seed_filename = '{}_galaxy_seed_image.fits'.format(self.basename)
            self.saveSeedImage(self.galaxy_source_seed, self.galaxy_source_seg_map, self.galaxy_seed_filename)
            print("Simulated galaxy image and segmap saved as {}".format(self.galaxy_seed_filename))

            # Add galaxy segmentation map to the master copy
            segmentation_map = self.add_segmentation_maps(segmentation_map, galaxy_segmap)

            # add the galaxy image to the signalimage
            signalimage = signalimage + galaxyCRImage

        else:
            self.galaxy_source_seed = None
            self.galaxy_source_seg_map = None
            self.galaxy_seed_filename = None

        # read in extended signal image and add the image to the overall image
        if self.runStep['extendedsource'] is True:
            extlist, extstamps = self.getExtendedSourceList(self.params['simSignals']['extended'])

            # translate the extended source list into an image
            extimage, ext_segmap = self.make_extended_source_image(extlist, extstamps)

            # Multiply by the flat field unless it is a NIRCam observation
            # that is going to be dispersed
            if self.instrument == 'nircam' and self.params['Inst']['mode'] in ['wfss', 'ts_grism']:
                pass
            else:
                extimage *= self.flatfield

            # To avoid problems with overlapping sources between source
            # types in observations to be dispersed, make the point
            # source-only segmentation map available as a class variable
            self.extended_source_seed = extimage
            self.extended_source_seg_map = ext_segmap

            # For seed images to be dispersed in WFSS mode,
            # embed the seed image in a full frame array. The disperser
            # tool does not work on subarrays
            aperture_suffix = self.params['Readout']['array_name'].split('_')[-1]
            if ((self.params['Inst']['mode'] in ['wfss', 'ts_grism']) & \
                 (aperture_suffix not in ['FULL', 'CEN'])):
                self.extended_source_seed, self.extended_source_seg_map = self.pad_wfss_subarray(self.extended_source_seed,
                                                                                                 self.extended_source_seg_map)

            # Save the extended source seed image
            if instrument_name != 'fgs':
                self.extended_seed_filename = os.path.join(self.basename + '_' + self.params['Readout'][self.usefilt] + '_extended_seed_image.fits')
            else:
                self.extended_seed_filename = '{}_extended_seed_image.fits'.format(self.basename)
            self.saveSeedImage(self.extended_source_seed, self.extended_source_seg_map, self.extended_seed_filename)
            print("Extended object image and segmap saved as {}".format(self.extended_seed_filename))

            # Add galaxy segmentation map to the master copy
            segmentation_map = add_segmentation_maps(segmentation_map, ext_segmap)

            # add the extended image to the synthetic signal rate image
            signalimage = signalimage + extimage

        else:
            self.extended_source_seed = None
            self.extended_source_seg_map = None
            self.extended_seed_filename = None

        # ZODIACAL LIGHT
        if self.runStep['zodiacal'] is True:
            print(("\n\nWARNING: A file has been provided for a zodiacal light contribution "
                   "but zodi is included in the background addition in imaging/imaging TSO modes "
                   "if bkgdrate is set to low/medium/high, or for any WFSS/Grism TSO observations.\n\n"))
            # zodiangle = self.eclipticangle() - self.params['Telescope']['rotation']
            zodiangle = self.params['Telescope']['rotation']
            zodiacalimage, zodiacalheader = self.getImage(self.params['simSignals']['zodiacal'], arrayshape, True, zodiangle, arrayshape/2)

            # Multiply by the flat field
            if self.instrument == 'nircam' and self.params['Inst']['mode'] in ['wfss', 'ts_grism']:
                pass
            else:
                zodiacalimage *= self.flatfield

            signalimage = signalimage + zodiacalimage*self.params['simSignals']['zodiscale']

        # SCATTERED LIGHT - no rotation here.
        if self.runStep['scattered']:
            scatteredimage, scatteredheader = self.getImage(self.params['simSignals']['scattered'],
                                                            arrayshape, False, 0.0, arrayshape/2)

            # Multiply by the flat field
            if self.instrument == 'nircam' and self.params['Inst']['mode'] in ['wfss', 'ts_grism']:
                pass
            else:
                scatteredimage *= self.flatfield

            signalimage = signalimage + scatteredimage*self.params['simSignals']['scatteredscale']

        # CONSTANT BACKGROUND - multiply by flat field unless it is a NIRCam
        # observation that is going to be dispersed
        if self.params['Inst']['mode'] not in ['wfss', 'ts_grism']:
            signalimage = signalimage + (self.params['simSignals']['bkgdrate'] * self.flatfield)

        # Save the final rate image of added signals
        if self.params['Output']['save_intermediates'] is True:
            rateImageName = self.basename + '_AddedSources_adu_per_sec.fits'
            self.saveSingleFits(signalimage, rateImageName)
            print("Signal rate image of all added sources saved as {}".format(rateImageName))

        return signalimage, segmentation_map

    @staticmethod
    def add_segmentation_maps(map1, map2):
        """Add two segmentation maps together. In the case of overlapping
        objects, the object in map1 is kept and the object in map2 is
        ignored.

        Parameters
        ----------
        map1 : numpy.ndarray
            2D segmentation map

        map2 : numpy.ndarray
            2D segmentation map

        Returns
        -------
        combined : numpy.ndarray
            Summed segmentation map
        """
        map1_zeros = map1 == 0
        combined = copy.deepcopy(map1)
        combined[map1_zeros] += map2[map1_zeros]
        return combined

    def get_point_source_list(self, filename, source_type='pointsources', segment_offset=None):
        """Read in the list of point sources to add, calculate positions
        on the detector, filter out sources outside the detector, and
        calculate countrates corresponding to the given magnitudes

        Parameters
        ----------
        filename : str
            Name of catalog file to examine

        source_type : str
            Flag specifying exactly what type of catalog is being read.
            This is because there are small differences between catalog
            types. Options are ``pointsources`` which is the default,
            ``ts_imaging`` for time series, and ``ts_grism`` for grism
            time series.

        Returns
        -------
        pointSourceList : astropy.table.Table
            Table containing source information
        """
        # Make sure that a valid PSF path has been provided
        if not os.path.isdir(self.params['simSignals']['psfpath']):
            raise ValueError('Invalid PSF path provided in YAML:',
                             self.params['simSignals']['psfpath'])

        pointSourceList = Table(names=('index', 'pixelx', 'pixely', 'RA', 'Dec', 'RA_degrees',
                                       'Dec_degrees', 'magnitude', 'countrate_e/s',
                                       'counts_per_frame_e', 'lightcurve_file'),
                                dtype=('i', 'f', 'f', 'S14', 'S14', 'f', 'f', 'f', 'f', 'f', 'S50'))

        try:
            lines, pixelflag, magsys = self.read_point_source_file(filename)
            if pixelflag:
                print("Point source list input positions assumed to be in units of pixels.")
            else:
                print("Point list input positions assumed to be in units of RA and Dec.")
        except:
            raise NameError("WARNING: Unable to open the point source list file {}".format(filename))

        # Create table of point source countrate versus psf size
        if self.add_psf_wings is True:
            self.translate_psf_table(magsys)

        # File to save adjusted point source locations
        psfile = self.params['Output']['file'][0:-5] + '_{}.list'.format(source_type)
        pslist = open(psfile, 'w')

        # If the input catalog has an index column
        # use that, otherwise add one
        indexes = self.get_index_numbers(lines)

        dtor = math.radians(1.)
        nx = (self.subarray_bounds[2] - self.subarray_bounds[0]) + 1
        ny = (self.subarray_bounds[3] - self.subarray_bounds[1]) + 1
        xc = (self.subarray_bounds[2] + self.subarray_bounds[0]) / 2.
        yc = (self.subarray_bounds[3] + self.subarray_bounds[1]) / 2.

        # Define the min and max source locations (in pixels) that fall onto the subarray
        # Include the effects of a requested grism_direct image, and also keep sources that
        # will only partially fall on the subarray
        # pixel coords here can still be negative and kept if the grism image is being made

        # First, coord limits for just the subarray
        miny = 0
        maxy = self.subarray_bounds[3] - self.subarray_bounds[1]
        minx = 0
        maxx = self.subarray_bounds[2] - self.subarray_bounds[0]

        # Expand the limits if a grism direct image is being made
        if (self.params['Output']['grism_source_image'] == True) or (self.params['Inst']['mode'] in ["pom", "wfss"]):
            extrapixy = np.int((maxy + 1)/2 * (self.coord_adjust['y'] - 1.))
            miny -= extrapixy
            maxy += extrapixy
            extrapixx = np.int((maxx + 1)/2 * (self.coord_adjust['x'] - 1.))
            minx -= extrapixx
            maxx += extrapixx

        # Write out the RA and Dec of the field center to the output file
        # Also write out column headers to prepare for source list
        pslist.write(("# Field center (degrees): %13.8f %14.8f y axis rotation angle "
                      "(degrees): %f  image size: %4.4d %4.4d\n" %
                      (self.ra, self.dec, self.params['Telescope']['rotation'], nx, ny)))
        pslist.write('#\n')
        pslist.write(("#    Index   RA_(hh:mm:ss)   DEC_(dd:mm:ss)   RA_degrees      "
                      "DEC_degrees     pixel_x   pixel_y    magnitude   counts/sec    counts/frame    TSO_lightcurve_catalog\n"))

        # If creating a segment-wise simulation, shift all of the RAs/Decs in
        # the list by the given offset
        if segment_offset is not None:
            lines = self.shift_sources_by_offset(lines, segment_offset, pixelflag)

        # Check the source list and remove any sources that are well outside the
        # field of view of the detector. These sources cause the coordinate
        # conversion to hang.
        indexes, lines = self.remove_outside_fov_sources(indexes, lines, pixelflag, 4096)

        # Determine the name of the column to use for source magnitudes
        mag_column = self.select_magnitude_column(lines, filename)

        print('Filtering point sources to keep only those on the detector')
        start_time = time.time()
        times = []
        time_reported = False
        # Loop over input lines in the source list
        for index, values in zip(indexes, lines):
            # Warn user of how long this calcuation might take...
            if len(times) < 100:
                elapsed_time = time.time() - start_time
                times.append(elapsed_time)
                start_time = time.time()
            elif len(times) == 100 and time_reported is False:
                avg_time = np.mean(times)
                total_time = len(indexes) * avg_time
                print(("Expected time to process {} sources: {:.2f} seconds "
                       "({:.2f} minutes)".format(len(indexes), total_time, total_time/60)))
                time_reported = True

            pixelx, pixely, ra, dec, ra_str, dec_str = self.get_positions(values['x_or_RA'],
                                                                          values['y_or_Dec'],
                                                                          pixelflag, 4096)

            # Get the input magnitude and countrate of the point source
            mag = float(values[mag_column])
            countrate = utils.magnitude_to_countrate(self.params['Inst']['mode'], magsys, mag,
                                                     photfnu=self.photfnu, photflam=self.photflam,
                                                     vegamag_zeropoint=self.vegazeropoint)

            psf_len = self.find_psf_size(countrate)
            edgex = np.int(psf_len // 2)
            edgey = np.int(psf_len // 2)

            if pixely > (miny-edgey) and pixely < (maxy+edgey) and pixelx > (minx-edgex) and pixelx < (maxx+edgex):
                # set up an entry for the output table
                entry = [index, pixelx, pixely, ra_str, dec_str, ra, dec, mag]

                # Calculate the countrate for the source
                framecounts = countrate * self.frametime

                # add the countrate and the counts per frame to pointSourceList
                # since they will be used in future calculations
                entry.append(countrate)
                entry.append(framecounts)

                # Add the TSO catalog name if present
                if source_type == 'ts_imaging':
                    tso_catalog = values['lightcurve_file']
                else:
                    tso_catalog = 'None'
                entry.append(tso_catalog)

                # add the good point source, including location and counts, to the pointSourceList
                pointSourceList.add_row(entry)

                # write out positions, distances, and counts to the output file
                pslist.write("%i %s %s %14.8f %14.8f %9.3f %9.3f  %9.3f  %13.6e   %13.6e  %s\n" %
                             (index, ra_str, dec_str, ra, dec, pixelx, pixely, mag, countrate, framecounts, tso_catalog))

        self.n_pointsources = len(pointSourceList)
        print("Number of point sources found within or close to the requested aperture: {}".format(self.n_pointsources))
        # close the output file
        pslist.close()

        # If no good point sources were found in the requested array, alert the user
        if len(pointSourceList) < 1:
            print("INFO: no point sources within the requested array.")
            # print("The point source image option is being turned off")
            # self.runStep['pointsource']=False
            # if self.runStep['extendedsource'] == False and self.runStep['cosmicray'] == False:
            #    print("Error: no input point sources, extended image, nor cosmic rays specified")
            #    print("Exiting...")
            #    sys.exit()

        return pointSourceList

    def translate_psf_table(self, magnitude_system):
        """Given a magnitude system, translate the table of PSF sizes
        versus magnitudes into PSF sizes versus countrates

        Parameters
        ----------
        magnitude_system : str
            Magnitude system of the sources: 'abmag', 'stmag', 'vegamag'
        """
        magnitudes = self.psf_wing_sizes[magnitude_system].data

        # Place table in order of ascending magnitudes
        sort = np.argsort(magnitudes)
        for colname in self.psf_wing_sizes.colnames:
            self.psf_wing_sizes[colname] = self.psf_wing_sizes[colname][sort]
        magnitudes = self.psf_wing_sizes[magnitude_system].data

        # Calculate corresponding countrates
        countrates = utils.magnitude_to_countrate(self.params['Inst']['mode'], magnitude_system,
                                                  magnitudes, photfnu=self.photfnu,
                                                  photflam=self.photflam,
                                                  vegamag_zeropoint=self.vegazeropoint)
        self.psf_wing_sizes['countrate'] = countrates

    def find_psf_size(self, countrate):
        """Determine the dimentions of the PSF to use based on an object's
        countrate.

        Parameters
        ----------
        countrate : float
            Source countrate

        Returns
        -------
        xdim : int
            Size of PSF in pixels in the x direction

        ydim : int
            Size of PSF in pixels in the y direction
        """
        if self.add_psf_wings is False:
            return self.psf_library_core_x_dim

        brighter = np.where(countrate >= self.psf_wing_sizes['countrate'])[0]
        if len(brighter) == 0:
            # Dimmest bin == size of pf library
            dimension = self.psf_library_core_x_dim
            #dim = np.max(self.psf_wing_sizes['number_of_pixels'])
        else:
            dimension = self.psf_wing_sizes['number_of_pixels'][brighter[0]]
        return dimension

    def shift_sources_by_offset(self, lines, segment_offset, pixelflag):
        print('    Shifting point source locations by arcsecond offset {}'.format(segment_offset))

        shifted_lines = lines.copy()
        shifted_lines.remove_rows(np.arange(0, len(shifted_lines)))

        V2ref_arcsec = self.siaf.V2Ref
        V3ref_arcsec = self.siaf.V3Ref
        position_angle = self.params['Telescope']['rotation']
        print('    Position angle = ', position_angle)
        attitude_ref = pysiaf.utils.rotations.attitude(V2ref_arcsec, V3ref_arcsec, self.ra, self.dec, position_angle)

        # Shift every source by the appropriate offset
        x_displacement_arcsec, y_displacement_arcsec = segment_offset
        for line in lines:
            x_or_RA, y_or_Dec = line['x_or_RA', 'y_or_Dec']
            # Convert the input source locations to V2/V3
            if not pixelflag:
                # Convert RA/Dec (sky frame) to V2/V3 (telescope frame)
                v2, v3 = pysiaf.utils.rotations.getv2v3(attitude_ref, x_or_RA, y_or_Dec)
            else:
                # Convert X/Y (detector frame) to V2/V3 (telescope frame)
                v2, v3 = self.siaf.det_to_tel(x_or_RA, y_or_Dec)

            # Add the arcsecond displacement to each V2/V3 source position
            v2 -= x_displacement_arcsec
            v3 += y_displacement_arcsec

            # Translate back to RA/Dec and save
            ra, dec = pysiaf.utils.rotations.pointing(attitude_ref, v2, v3)
            # TODO: need to be smarter about which magnitude to use here
            shifted_lines.add_row([ra, dec, line['magnitude']])

        return shifted_lines

    def remove_outside_fov_sources(self, index, source, pixflag, delta_pixels):
        """Filter out entries in the source catalog that are located well outside the field of
        view of the detector. This can be a fairly rough cut. We just need to remove sources
        that are very far from the detector.

        Parameters:
        -----------
        index : list
            List of index numbers corresponding to the sources

        source : Table
            astropy Table containing catalog information

        pixflag : bool
            Flag indicating whether catalog positions are given in units of
            pixels (True), or RA, Dec (False)

        delta_pixels : int
            Number of columns/rows outside of the nominal detector size (2048x2048)
            to keep sources in the source list. (e.g. delta_pixels=2048 will keep all
            sources located at -2048 to 4096.)

        Returns:
        --------
        index : list
            List of filtered index numbers corresponding to sources
            within or close to the field of view

        source : Table
            astropy Table containing filtered list of sources
        """
        catalog_x = source['x_or_RA']
        catalog_y = source['y_or_Dec']

        if pixflag:
            minx = 0 - delta_pixels
            maxx = self.output_dims[1] + delta_pixels
            miny = 0 - delta_pixels
            maxy = self.output_dims[0] + delta_pixels
            good = ((catalog_x > minx) & (catalog_x < maxx) & (catalog_y > miny) & (catalog_y < maxy))
        else:
            delta_degrees = (delta_pixels * self.siaf.XSciScale) / 3600. * u.deg
            reference = SkyCoord(ra=self.ra * u.deg, dec=self.dec * u.deg)

            # Need to determine the units of the RA values.
            # Dec units should always be degrees whether or not they are in decimal
            # or DD:MM:SS or DDd:MMm:SSs formats.
            dec_unit = u.deg
            try:
                # if it can be converted to a float, assume decimal degrees
                entry0 = float(catalog_x[0])
                ra_unit = u.deg
            except ValueError:
                # if it cannot be converted to a float, then the unit is 'hour'
                ra_unit = 'hour'

            # Assume that units are consisent within each column. (i.e. no mixing of
            # 12h:23m:34.5s and 189.87463 degrees within a column)
            catalog = SkyCoord(ra=catalog_x, dec=catalog_y, unit=(ra_unit, dec_unit))
            good = np.where(reference.separation(catalog) < delta_degrees)[0]

        filtered_sources = source[good]
        filtered_indexes = index[good]

        return filtered_indexes, filtered_sources

    def make_point_source_image(self, pointSources, segment_number=None, ptsrc_segmap=None):
        """Create a seed image containing all of the point sources
        provided by the source catalog

        Parameters
        ----------
        pointSources : astropy.table.Table
            Table of point sources
        segment_number : int, optional
            The number of the mirror segment to make an image for
        ptsrc_segmap : optional
            The point source segmentation map to keep adding to

        Returns
        -------
        psfimage : numpy.ndarray
            2D array containing the seed image with point sources

        seg.segmap : numpy.ndarray
            2D array containing the segmentation map that
            corresponds to ``psfimage``
        """
        dims = np.array(self.nominal_dims)

        # Create the empty image
        psfimage = np.zeros(self.output_dims)

        if ptsrc_segmap is None:
            # Create empty segmentation map
            ptsrc_segmap = segmap.SegMap()
            ptsrc_segmap.ydim, ptsrc_segmap.xdim = self.output_dims
            ptsrc_segmap.initialize_map()

        # Loop over the entries in the point source list
        for i, entry in enumerate(pointSources):

            # Find the PSF size to use based on the countrate
            psf_x_dim = self.find_psf_size(entry['countrate_e/s'])

            # Assume same PSF size in x and y
            psf_y_dim = psf_x_dim

            scaled_psf, min_x, min_y, wings_added = self.create_psf_stamp(
                entry['pixelx'], entry['pixely'], psf_x_dim, psf_y_dim,
                segment_number=segment_number
            )

            # Skip sources that fall completely off the detector
            if scaled_psf is None:
                continue

            scaled_psf *= entry['countrate_e/s']

            # PSF may not be centered in array now if part of the array falls
            # off of the aperture
            stamp_x_loc = psf_x_dim // 2 - min_x
            stamp_y_loc = psf_y_dim // 2 - min_y
            updated_psf_dimensions = scaled_psf.shape

            # If the source subpixel location is beyond 0.5 (i.e. the edge
            # of the pixel), then we shift the wing->core offset by 1.
            # We also need to shift the location of the wing array on the
            # detector by 1
            if wings_added:
                x_delta = int(np.modf(entry['pixelx'])[0] > 0.5)
                y_delta = int(np.modf(entry['pixely'])[0] > 0.5)
            else:
                x_delta = 0
                y_delta = 0

            # Get the coordinates that describe the overlap between the
            # PSF image and the output aperture
            xap, yap, xpts, ypts, (i1, i2), (j1, j2), (k1, k2), \
                (l1, l2) = self.create_psf_stamp_coords(entry['pixelx']+x_delta, entry['pixely']+y_delta,
                                                        updated_psf_dimensions,
                                                        stamp_x_loc, stamp_y_loc,
                                                        coord_sys='aperture')

            # Skip sources that fall completely off the detector
            if None in [i1, i2, j1, j2, k1, k2, l1, l2]:
                continue

            try:
                psfimage[j1:j2, i1:i2] += scaled_psf[l1:l2, k1:k2]

                # Divide readnoise by 100 sec, which is a 10 group RAPID ramp?
                noiseval = self.single_ron / 100. + self.params['simSignals']['bkgdrate']
                if self.params['Inst']['mode'].lower() in ['wfss', 'ts_grism']:
                    noiseval += self.grism_background
                ptsrc_segmap.add_object_noise(scaled_psf[l1:l2, k1:k2], j1, i1, entry['index'], noiseval)
            except IndexError:
                # In here we catch sources that are off the edge
                # of the detector. These may not necessarily be caught in
                # getpointsourcelist because if the PSF is not centered
                # in the webbpsf stamp, then the area to be pulled from
                # the stamp may shift off of the detector.
                pass

            if ((len(pointSources) > 100) and (np.mod(i, 100))) == 0:
                print('{}: Working on source {}'.format(str(datetime.datetime.now()), i))

        return psfimage, ptsrc_segmap

    def create_psf_stamp(self, x_location, y_location, psf_dim_x, psf_dim_y,
                         ignore_detector=False, segment_number=None):
        """From the gridded PSF model, location within the aperture, and
        dimensions of the stamp image (either the library PSF image, or
        the galaxy/extended stamp image with which the PSF will be
        convolved), evaluate the GriddedPSFModel at
        the appropriate location on the detector and return the PSF stamp

        Parameters
        ----------
        x_location : float
            X-coordinate of the PSF in the coordinate system of the
            aperture being simulated.

        y_location : float
            Y-coordinate of the PSF in the coordinate system of the
            aperture being simulated.

        psf_dim_x : int
            Number of columns of the array containing the PSF

        psf_dim_y : int
            Number of rows of the array containing the PSF

        ignore_detector : bool
            If True, the returned coordinates can have values outside the
            size of the subarray/detector (i.e. coords can be negative or
            larger than full frame). If False, coordinates are constrained
            to be on the detector.

        Returns
        -------
        full_psf : numpy.ndarray
            2D array containing the normalized PSF image. Total signal should
            be close to 1.0 (not exactly 1.0 due to asymmetries and distortion)
            Array will be copped based on how much falls on or off the detector

        k1 : int
            Row number on the PSF/stamp image corresponding to the bottom-most
            row that overlaps the detector/aperture

        l1 : int
            Column number on the PSF/stamp image corresponding to the left-most
            column that overlaps the detector/aperture

        add_wings : bool
            Whether or not PSF wings are to be added to the PSF core
        """
        # PSF will always be centered
        psf_x_loc = psf_dim_x // 2
        psf_y_loc = psf_dim_y // 2

        # Translation needed to go from PSF core (self.psf_library)
        # coordinate system to the PSF wing coordinate system (i.e.
        # center the PSF core in the wing image)
        psf_wing_half_width_x = np.int(psf_dim_x // 2)
        psf_wing_half_width_y = np.int(psf_dim_y // 2)
        psf_core_half_width_x = np.int(self.psf_library_core_x_dim // 2)
        psf_core_half_width_y = np.int(self.psf_library_core_y_dim // 2)
        delta_core_to_wing_x = psf_wing_half_width_x - psf_core_half_width_x
        delta_core_to_wing_y = psf_wing_half_width_y - psf_core_half_width_y

        # This assumes a square PSF shape!!!!
        # If no wings are to be added, then we can skip all the wing-
        # and pixel phase-related work below.
        if ((self.add_psf_wings is False) or (delta_core_to_wing_x <= 0)):
            add_wings = False

            if segment_number is not None:
                library = self.psf_library[segment_number - 1]
            else:
                library = self.psf_library

            # Get coordinates decribing overlap between the evaluated psf
            # core and the full frame of the detector. We really only need
            # the xpts_core and ypts_core from this in order to know how
            # to evaluate the library
            # Note that we don't care about the pixel phase here.
            psf_core_dims = (self.psf_library_core_y_dim, self.psf_library_core_x_dim)
            xc_core, yc_core, xpts_core, ypts_core, (i1c, i2c), (j1c, j2c), (k1c, k2c), \
                (l1c, l2c) = self.create_psf_stamp_coords(x_location, y_location, psf_core_dims,
                                                          psf_core_half_width_x, psf_core_half_width_y,
                                                          coord_sys='full_frame',
                                                          ignore_detector=ignore_detector)

            # Skip sources that fall completely off the detector
            if None in [i1c, i2c, j1c, j2c, k1c, k2c, l1c, l2c]:
                return None, None, None, False

            # Step 4
            full_psf = library.evaluate(x=xpts_core, y=ypts_core, flux=1.0,
                                        x_0=xc_core, y_0=yc_core)
            k1 = k1c
            l1 = l1c

        else:
            add_wings = True
            # If the source subpixel location is beyond 0.5 (i.e. the edge
            # of the pixel), then we shift the wing->core offset by 1.
            # We also need to shift the location of the wing array on the
            # detector by 1
            x_phase = np.modf(x_location)[0]
            y_phase = np.modf(y_location)[0]
            x_location_delta = int(x_phase > 0.5)
            y_location_delta = int(y_phase > 0.5)
            if x_phase > 0.5:
                delta_core_to_wing_x -= 1
            if y_phase > 0.5:
                delta_core_to_wing_y -= 1

            # offset_x, and y below will not change because that is
            # the offset between the full wing array and the user-specified
            # wing array size

            # Get the psf wings array - first the nominal size
            # Later we may crop if the source is only partially on the detector
            full_wing_y_dim, full_wing_x_dim = self.psf_wings.shape
            offset_x = np.int((full_wing_x_dim - psf_dim_x) / 2)
            offset_y = np.int((full_wing_y_dim - psf_dim_y) / 2)

            full_psf = copy.deepcopy(self.psf_wings[offset_y:offset_y+psf_dim_y, offset_x:offset_x+psf_dim_x])

            # Get coordinates describing overlap between PSF image and the
            # full frame of the detector
            # Step 1
            xcenter, ycenter, xpts, ypts, (i1, i2), (j1, j2), (k1, k2), \
                (l1, l2) = self.create_psf_stamp_coords(x_location+x_location_delta, y_location+y_location_delta,
                                                        (psf_dim_y, psf_dim_x), psf_x_loc, psf_y_loc,
                                                        coord_sys='full_frame', ignore_detector=ignore_detector)

            if None in [i1, i2, j1, j2, k1, k2, l1, l2]:
                return None, None, None, False

            # Step 2
            # If the core of the psf lands at least partially on the detector
            # then we need to evaluate the psf library
            if ((k1 < (psf_wing_half_width_x + psf_core_half_width_x)) and
               (k2 > (psf_wing_half_width_x - psf_core_half_width_x)) and
               (l1 < (psf_wing_half_width_y + psf_core_half_width_y)) and
               (l2 > (psf_wing_half_width_y - psf_core_half_width_y))):

                # Step 3
                # Get coordinates decribing overlap between the evaluated psf
                # core and the full frame of the detector. We really only need
                # the xpts_core and ypts_core from this in order to know how
                # to evaluate the library
                # Note that we don't care about the pixel phase here.
                psf_core_dims = (self.psf_library_core_y_dim, self.psf_library_core_x_dim)
                xc_core, yc_core, xpts_core, ypts_core, (i1c, i2c), (j1c, j2c), (k1c, k2c), \
                    (l1c, l2c) = self.create_psf_stamp_coords(x_location, y_location, psf_core_dims,
                                                              psf_core_half_width_x, psf_core_half_width_y,
                                                              coord_sys='full_frame', ignore_detector=ignore_detector)

                if None in [i1c, i2c, j1c, j2c, k1c, k2c, l1c, l2c]:
                    return None, None, None, False

                # Step 4
                psf = self.psf_library.evaluate(x=xpts_core, y=ypts_core, flux=1.,
                                                x_0=xc_core, y_0=yc_core)

                # Step 5
                wing_start_x = k1c + delta_core_to_wing_x
                wing_end_x = k2c + delta_core_to_wing_x
                wing_start_y = l1c + delta_core_to_wing_y
                wing_end_y = l2c + delta_core_to_wing_y

                full_psf[wing_start_y:wing_end_y, wing_start_x:wing_end_x] = psf

            # Whether or not the core is on the detector, crop the PSF
            # to the proper shape based on how much is on the detector
            full_psf = full_psf[l1:l2, k1:k2]

        return full_psf, k1, l1, add_wings

    def create_psf_stamp_coords(self, aperture_x, aperture_y, stamp_dims, stamp_x, stamp_y,
                                coord_sys='full_frame', ignore_detector=False):
        """Calculate the coordinates in the aperture coordinate system
        where the stamp image wil be placed based on the location of the
        stamp image in the aperture and the size of the stamp image.

        Parameters
        ----------
        aperture_x : float
            X-coordinate of the PSF in the coordinate system of the
            aperture being simulated.

        aperture_y : float
            Y-coordinate of the PSF in the coordinate system of the
            aperture being simulated.

        stamp_dims : tup
            (x, y) dimensions of the stamp image that will be placed
            into the final seed image. This stamp image can be either the
            PSF image itself, or the stamp image of the galaxy/extended
            source that the PSF is going to be convolved with.

        stamp_x : float
            Location in x of source within the stamp image

        stamp_y : float
            Location in y of source within the stamp image

        coord_sys : str
            Inidicates which coordinate system to return coordinates for.
            Options are 'full_frame' for full frame coordinates, or
            'aperture' for aperture coordinates (including any expansion
            for grism source image)

        ignore_detector : bool
            If True, the returned coordinates can have values outside the
            size of the subarray/detector (i.e. coords can be negative or
            larger than full frame). If False, coordinates are constrained
            to be on the detector.

        Returns
        -------
        x_points : numpy.ndarray
            2D array of x-coordinates in the aperture coordinate system
            where the stamp image will fall.

        y_points : numpy.ndarray
            2D array of y-coordinates in the aperture coordinate system
            where the stamp image will fall.

        (i1, i2) : tup
            Beginning and ending x coordinates (in the aperture coordinate
            system) where the stamp image will fall

        (j1, j2) : tup
            Beginning and ending y coordinates (in the aperture coordinate
            system) where the stamp image will fall

        (k1, k2) : tup
            Beginning and ending x coordinates (in the stamp's coordinate
            system) that overlap the aperture

        (l1, l2) : tup
            Beginning and ending y coordinates (in the stamp's coordinate
            system) that overlap the aperture
        """
        if coord_sys == 'full_frame':
            xpos = aperture_x + self.subarray_bounds[0]
            ypos = aperture_y + self.subarray_bounds[1]
            out_dims_x = self.ffsize
            out_dims_y = self.ffsize
        elif coord_sys == 'aperture':
            xpos = aperture_x + self.coord_adjust['xoffset']
            ypos = aperture_y + self.coord_adjust['yoffset']
            out_dims_x = self.output_dims[1]
            out_dims_y = self.output_dims[0]

        stamp_y_dim, stamp_x_dim = stamp_dims

        # Get coordinates that describe the overlap between the stamp
        # and the aperture
        (i1, i2, j1, j2, k1, k2, l1, l2) = self.cropped_coords(xpos, ypos, (out_dims_y, out_dims_x),
                                                               stamp_x, stamp_y, stamp_dims,
                                                               ignore_detector=ignore_detector)

        # If the stamp is completely off the detector, use dummy arrays
        # for x_points and y_points
        if j1 is None or j2 is None or i1 is None or i2 is None:
            x_points = np.zeros((2, 2))
            y_points = x_points
        else:
            y_points, x_points = np.mgrid[j1:j2, i1:i2]

        return xpos, ypos, x_points, y_points, (i1, i2), (j1, j2), (k1, k2), (l1, l2)

    def ensure_odd_lengths(self, x_dim, y_dim, x_center, y_center):
        """Given the dimensions and the coordinates of the center of an
        array, ensure the array has an odd number of rows and columns,
        calculate the updated half-width, and return the minimum and
        maximum row and column indexes.

        Parameters
        ----------
        x_dim : int
            Length of the array in the x-dimension

        y_dim : int
            Length of the array in the y-dimension

        x_center : float
            Coordinate of the center of the array, usually
            in some other coordinate system (e.g. full frame
            coords, while the array is a subarray)

        y_center : float
            Coordinate of the center of the array in the y
            direction, usually in some other coordinate
            system

        Returns
        -------
        x_min : int
            Minimum index in the x direction of the array
            in the coordinate system of ``x_center, y_center``.

        x_max : int
            Maximum index in the x direction of the array
            in the coordinate system of ``x_center, y_center``.

        y_min : int
            Minimum index in the y direction of the array
            in the coordinate system of ``x_center, y_center``.

        y_max : int
            Maximum index in the y direction of the array
            in the coordinate system of ``x_center, y_center``.
        """
        if x_dim % 2 == 0:
            x_dim -= 1
        if y_dim % 2 == 0:
            y_dim -= 1
        x_half_width = x_dim // 2
        y_half_width = y_dim // 2
        x_min = np.int(x_center) - x_half_width
        x_max = np.int(x_center) + x_half_width + 1
        y_min = np.int(y_center) - y_half_width
        y_max = np.int(y_center) + y_half_width + 1
        return x_min, x_max, y_min, y_max

    def cropped_coords(self, aperture_x, aperture_y, aperture_dims, stamp_x, stamp_y, stamp_dims,
                       ignore_detector=False):
        """Given the location of a source on the detector along with the size of
        the PSF/stamp image for that object, calcuate the limits of the detector
        coordinates onto which the object will fall.

        Parameters:
        -----------
        aperture_x : float
            Column location of source on detector (aperture coordinate system
            including any padding for WFSS seed image)

        aperture_y : float
            Row location of source on detector (aperture coordinate system
            including any padding for WFSS seed image)

        aperture_dims : tup
            (y, x) dimensions of the aperture on which the sources will be placed
            (e.g. full frame, full_frame+extra, subarray)

        stamp_x : float
            Location in x of source within the stamp image

        stamp_y : float
            Location in y of source within the stamp image

        stamp_dims : tup
            (y, x) dimensions of the source's stamp image. (e.g. the size of the PSF
            or galaxy image stamp)

        ignore_detector: bool
            If True, the returned coordinates can have values outside the
            size of the subarray/detector (i.e. coords can be negative or
            larger than full frame). If False, coordinates are constrained
            to be on the detector.

        Returns:
        --------
        i1 : int
            Column number on the detector/aperture corresponding to the left
            edge of the PSF/stamp image.

        i2 : int
            Column number on the detector/aperture corresponding to the right
            edge of the PSF/stamp image.

        j1 : int
            Row number on the detector/aperture corresponding to the bottom
            edge of the PSF/stamp image.

        j2 : int
            Row number on the detector/aperture corresponding to the top
            edge of the PSF/stamp image.

        l1 : int
            Column number on the PSF/stamp image corresponding to the left-most
            column that overlaps the detector/aperture

        l2 : int
            Column number on the PSF/stamp image corresponding to the right-most
            column that overlaps the detector/aperture

        k1 : int
            Row number on the PSF/stamp image corresponding to the bottom-most
            row that overlaps the detector/aperture

        k2 : int
            Row number on the PSF/stamp image corresponding to the top-most
            row that overlaps the detector/aperture
        """
        stamp_y_dim, stamp_x_dim = stamp_dims
        aperture_y_dim, aperture_x_dim = aperture_dims

        stamp_x = math.floor(stamp_x)
        stamp_y = math.floor(stamp_y)
        aperture_x = math.floor(aperture_x)
        aperture_y = math.floor(aperture_y)

        i1 = aperture_x - stamp_x
        j1 = aperture_y - stamp_y

        if ((i1 > (aperture_x_dim)) or (j1 > (aperture_y_dim))) and not ignore_detector:
            # In this case the stamp does not overlap the aperture at all
            return tuple([None]*8)
        delta_i1 = 0
        delta_j1 = 0

        if not ignore_detector:
            if i1 < 0:
                delta_i1 = copy.deepcopy(i1)
                i1 = 0
            if j1 < 0:
                delta_j1 = copy.deepcopy(j1)
                j1 = 0

        i2 = i1 + (stamp_x_dim + delta_i1)
        j2 = j1 + (stamp_y_dim + delta_j1)

        if ((i2 < 0) or (j2 < 0)) and not ignore_detector:
            # Stamp does not overlap the aperture at all
            return tuple([None]*8)

        delta_i2 = 0
        delta_j2 = 0
        if not ignore_detector:
            if i2 > aperture_x_dim:
                delta_i2 = i2 - aperture_x_dim
                i2 = aperture_x_dim
            if j2 > aperture_y_dim:
                delta_j2 = j2 - aperture_y_dim
                j2 = aperture_y_dim

        k1 = 0 - delta_i1
        k2 = stamp_x_dim - delta_i2
        l1 = 0 - delta_j1
        l2 = stamp_y_dim - delta_j2
        return (i1, i2, j1, j2, k1, k2, l1, l2)

    def cropPSF(self, psf):
        '''take an array containing a psf and crop it such that the brightest
        pixel is in the center of the array'''
        nyshift, nxshift = np.where(psf == np.max(psf))
        nyshift = nyshift[0]
        nxshift = nxshift[0]
        py, px = psf.shape

        xl = nxshift - 0
        xr = px - nxshift - 1
        if xl <= xr:
            xdist = xl
        if xr < xl:
            xdist = xr

        yl = nyshift - 0
        yr = py - nyshift - 1
        if yl <= yr:
            ydist = yl
        if yr < yl:
            ydist = yr

        return psf[nyshift - ydist:nyshift + ydist + 1, nxshift - xdist:nxshift + xdist + 1]

    def read_point_source_file(self, filename):
        """Read in the point source catalog file

         Parameters:
        -----------
        filename : str
            Filename of catalog file to be read in

         Returns:
        --------
        gtab : Table
            astropy Table containing catalog

         pflag : bool
            Flag indicating units of source locations. True for detector
            pixels, False for RA, Dec

         msys : str
            Magnitude system of the source brightnesses (e.g. 'abmag')
        """
        try:
            gtab = ascii.read(filename)
            # Look at the header lines to see if inputs
            # are in units of pixels or RA, Dec
            pflag = False
            try:
                if 'position_pixels' in gtab.meta['comments'][0:4]:
                    pflag = True
            except:
                pass
            # Check to see if magnitude system is specified
            # in the comments. If not default to AB mag
            msys = 'abmag'
            condition = ('stmag' in gtab.meta['comments'][0:4]) | ('vegamag' in gtab.meta['comments'][0:4])
            if condition:
                msys = [l for l in gtab.meta['comments'][0:4] if 'mag' in l][0]
                msys = msys.lower()

        except:
            raise IOError("WARNING: Unable to open the source list file {}".format(filename))

        return gtab, pflag, msys

    def select_magnitude_column(self, catalog, catalog_file_name):
        """Select the appropriate column to use for source magnitudes from the input source catalog. If there
        is a specific column name constructed as <instrument>_<filter>_magnitude to use for source magnitudes
        (e.g. nircam_f200w_magnitude) then use that. (NOTE: for FGS we search for a column name of
        'fgs_magnitude'). If not, check for a generic 'magnitude' column. If neither are present, raise an
        error.

        Parameters
        ----------

        catalog : astropy.table.Table
            Source catalog

        catalog_file_name : str
            Name of the catalog file. Used only when raising errors.

        Returns
        -------

        specific_mag_col : str
            The name of the catalog column to use for source magnitudes
        """
        # Determine the filter name to look for
        if self.params['Inst']['instrument'].lower() in ['nircam', 'niriss']:
        #    if self.params['Readout']['pupil'][0].upper() == 'F':
        #        usefilt = 'pupil'
        #    else:
        #        usefilt = 'filter'
            filter_name = self.params['Readout'][self.usefilt].lower()
            # Construct the column header to look for
            specific_mag_col = "{}_{}_magnitude".format(self.params['Inst']['instrument'].lower(),
                                                        filter_name)

        elif self.params['Inst']['instrument'].lower() == 'fgs':
            specific_mag_col = "{}_magnitude".format(self.params['Readout']['array_name'].split('_')[0].lower())
            filter_name = 'none'

        # Search catalog column names.
        if specific_mag_col in catalog.colnames:
            print("Using {} column in {} for magnitudes".format(specific_mag_col,
                                                                os.path.split(catalog_file_name)[1]))
            return specific_mag_col

        elif 'magnitude' in catalog.colnames:
            print(("WARNING: Catalog {} does not have a magnitude column called {}, "
                   "but does have a generic 'magnitude' column. Continuing simulation using that."
                   .format(os.path.split(catalog_file_name)[1], specific_mag_col)))
            return "magnitude"
        else:
            raise ValueError(("WARNING: Catalog {} has no magnitude column specifically for {} {}, "
                              "nor a generic 'magnitude' column. Unable to proceed."
                              .format(os.path.split(catalog_file_name)[1], self.params['Inst']['instrument'],
                                      filter_name.upper())))

    def makePos(self, alpha1, delta1):
        # given a numerical RA/Dec pair, convert to string
        # values hh:mm:ss
        if ((alpha1 < 0) or (alpha1 >= 360.)):
            alpha1 = alpha1 % 360
        if delta1 < 0.:
            sign = "-"
            d1 = abs(delta1)
        else:
            sign = "+"
            d1 = delta1
        decd = int(d1)
        value = 60. * (d1 - float(decd))
        decm = int(value)
        decs = 60. * (value - decm)
        a1 = alpha1 / 15.0
        radeg = int(a1)
        value = 60. * (a1 - radeg)
        ramin = int(value)
        rasec = 60. * (value - ramin)
        alpha2 = "%2.2d:%2.2d:%7.4f" % (radeg, ramin, rasec)
        delta2 = "%1s%2.2d:%2.2d:%7.4f" % (sign, decd, decm, decs)
        alpha2 = alpha2.replace(" ", "0")
        delta2 = delta2.replace(" ", "0")

        return alpha2, delta2

    def RADecToXY_astrometric(self, ra, dec):
        """Translate backwards, RA, Dec to V2, V3. If a distortion reference file is
        provided, use that. Otherwise fall back to pysiaf.

        Parameters:
        -----------
        ra : float
            Right ascention value, in degrees, to be translated.

        dec : float
            Declination value, in degrees, to be translated.

        Returns:
        --------
        pixelx : float
            X coordinate value in the aperture corresponding to the input location

        pixely : float
            Y coordinate value in the aperture corresponding to the input location
        """
        loc_v2, loc_v3 = pysiaf.utils.rotations.getv2v3(self.attitude_matrix, ra, dec)

        if self.coord_transform is not None:
            # Use the distortion reference file to translate from V2, V3 to RA, Dec
            pixelx, pixely = self.coord_transform.inverse(loc_v2, loc_v3)
            pixelx -= self.subarray_bounds[0]
            pixely -= self.subarray_bounds[1]
        else:
            # print('SIAF: using {} to transform from tel to sci'.format(self.siaf.AperName))
            pixelx, pixely = self.siaf.tel_to_sci(loc_v2, loc_v3)
            # Subtract 1 from SAIF-derived results since SIAF works in a 1-indexed coord system
            pixelx -= 1
            pixely -= 1
        return pixelx, pixely

    def object_separation(self, radec1, radec2, wcs):
        """
        Calculate the distance between two points on the sky given their
        RA, Dec values. Also calculate the angle (east of north?) between
        the two points.

        Parameters:
        -----------
        radec1 : list
            2-element list giving the RA, Dec (in decimal degrees) for
            the first object

        radec2 : list
            2-element list giving the RA, Dec (in decimal degrees) for
            the second object

        Returns:
        --------
        distance : float
            Angular separation (in degrees) between the two objects

        angle : float
            Angle (east of north?) separating the two sources
        """
        c1 = SkyCoord(radec1[0]*u.degree, radec1[1]*u.degree, frame='icrs')
        c2 = SkyCoord(radec2[0]*u.degree, radec2[1]*u.degree, frame='icrs')
        sepra, sepdec = c1.spherical_offsets_to(c2).to_pixel(wcs)
        return sepra, sepdec

    def XYToRADec(self, pixelx, pixely):
        """Translate a given x, y location on the detector to RA, Dec. If a
        distortion reference file is provided, use that. Otherwise fall back to
        using pysiaf.

        Parameters:
        -----------
        pixelx : float
            X coordinate value in the aperture

        pixely : float
            Y coordinate value in the aperture

        Returns:
        --------
        ra : float
            Right ascention value in degrees

        dec : float
            Declination value in degrees

        ra_str : str
            Right ascention value in HH:MM:SS

        dec_str : str
            Declination value in DD:MM:SS
        """
        if self.coord_transform is not None:
            pixelx += self.subarray_bounds[0]
            pixely += self.subarray_bounds[1]
            loc_v2, loc_v3 = self.coord_transform(pixelx, pixely)
        else:
            # Use SIAF to do the calculations if the distortion reffile is
            # not present. In this case, add 1 to the input pixel values
            # since SIAF works in a 1-indexed coordinate system.
            loc_v2, loc_v3 = self.siaf.sci_to_tel(pixelx + 1, pixely + 1)

        ra, dec = pysiaf.utils.rotations.pointing(self.attitude_matrix, loc_v2, loc_v3)

        # Translate the RA/Dec floats to strings
        ra_str, dec_str = self.makePos(ra, dec)

        return ra, dec, ra_str, dec_str

    def readGalaxyFile(self, filename):
        # Read in the galaxy source list
        try:
            # read table
            gtab = ascii.read(filename)

            # Look at the header lines to see if inputs
            # are in units of pixels or RA, Dec
            pflag = False
            rpflag = False
            try:
                if 'position_pixels' in gtab.meta['comments'][0:4]:
                    pflag = True
            except:
                pass
            try:
                if 'radius_pixels' in gtab.meta['comments'][0:4]:
                    rpflag = True
            except:
                pass
            # Check to see if magnitude system is specified in the comments
            # If not assume AB mags
            msys = 'abmag'
            condition = ('stmag' in gtab.meta['comments'][0:4]) | ('vegamag' in gtab.meta['comments'][0:4])
            if condition:
                msys = [l for l in gtab.meta['comments'][0:4] if 'mag' in l][0]
                msys = msys.lower()

        except:
            raise IOError("WARNING: Unable to open the galaxy source list file {}".format(filename))

        return gtab, pflag, rpflag, msys

    def filterGalaxyList(self, galaxylist, pixelflag, radiusflag, magsystem, catfile):
        # given a list of galaxies (location, size, orientation, magnitude)
        # keep only those which will fall fully or partially on the output array

        filteredList = Table(names=('index', 'pixelx', 'pixely', 'RA', 'Dec',
                                    'RA_degrees', 'Dec_degrees', 'V2', 'V3',
                                    'radius', 'ellipticity', 'pos_angle',
                                    'sersic_index', 'magnitude', 'countrate_e/s',
                                    'counts_per_frame_e'),
                             dtype=('i', 'f', 'f', 'S14', 'S14', 'f', 'f', 'f',
                                    'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f'))

        # Each entry in galaxylist is:
        # index x_or_RA  y_or_Dec  radius  ellipticity  pos_angle  sersic_index  magnitude
        # remember that x/y are interpreted as coordinates in the output subarray
        # NOT full frame coordinates. This is the same as the point source list coords

        # First, begin to define the pixel limits beyond which a galaxy will be completely
        # outside of the field of view
        # First, coord limits for just the subarray
        miny = 0
        maxy = self.subarray_bounds[3] - self.subarray_bounds[1]
        minx = 0
        maxx = self.subarray_bounds[2] - self.subarray_bounds[0]
        ny = self.subarray_bounds[3] - self.subarray_bounds[1] + 1
        nx = self.subarray_bounds[2] - self.subarray_bounds[0] + 1

        # Expand the limits if a grism direct image is being made
        if (self.params['Output']['grism_source_image'] == True) or (self.params['Inst']['mode'] in ["pom", "wfss"]):
            extrapixy = np.int((maxy + 1)/2 * (self.grism_direct_factor - 1.))
            miny -= extrapixy
            maxy += extrapixy
            extrapixx = np.int((maxx + 1)/2 * (self.grism_direct_factor - 1.))
            minx -= extrapixx
            maxx += extrapixx

            nx = np.int(nx * self.grism_direct_factor)
            ny = np.int(ny * self.grism_direct_factor)

        # If an index column is present use that, otherwise
        # create one
        indexes = self.get_index_numbers(galaxylist)

        # Check the source list and remove any sources that are well outside the
        # field of view of the detector. These sources cause the coordinate
        # conversion to hang.
        indexes, galaxylist = self.remove_outside_fov_sources(indexes, galaxylist, pixelflag, 4096)

        # Determine the name of the column to use for source magnitudes
        mag_column = self.select_magnitude_column(galaxylist, catfile)

        # Loop over galaxy sources
        for index, source in zip(indexes, galaxylist):

            # If galaxy radii are given in units of arcseconds, translate to pixels
            if radiusflag is False:
                source['radius'] /= self.siaf.XSciScale

            # how many pixels beyond the nominal subarray edges can a source be located and
            # still have it fall partially on the subarray? Galaxy stamps are nominally set to
            # have a length and width equal to 100 times the requested radius.
            edgex = source['radius'] * 100 / 2 - 1
            edgey = source['radius'] * 100 / 2 - 1

            # reset the field of view limits for the size of the current stamp image
            outminy = miny - edgey
            outmaxy = maxy + edgey
            outminx = minx - edgex
            outmaxx = maxx + edgex

            pixelx, pixely, ra, dec, ra_str, dec_str = self.get_positions(source['x_or_RA'],
                                                                          source['y_or_Dec'],
                                                                          pixelflag, 4096)

            # only keep the source if the peak will fall within the subarray
            if pixely > outminy and pixely < outmaxy and pixelx > outminx and pixelx < outmaxx:

                pixelv2, pixelv3 = pysiaf.utils.rotations.getv2v3(self.attitude_matrix, ra, dec)
                entry = [index, pixelx, pixely, ra_str, dec_str, ra, dec, pixelv2, pixelv3,
                         source['radius'], source['ellipticity'], source['pos_angle'], source['sersic_index']]

                # Now look at the input magnitude of the point source
                # append the mag and pixel position to the list of ra, dec
                mag = float(source[mag_column])
                entry.append(mag)

                # Convert magnitudes to countrate (ADU/sec) and counts per frame
                rate = utils.magnitude_to_countrate(self.params['Inst']['mode'], magsystem, mag,
                                                    photfnu=self.photfnu, photflam=self.photflam,
                                                    vegamag_zeropoint=self.vegazeropoint)
                framecounts = rate * self.frametime

                # add the countrate and the counts per frame to pointSourceList
                # since they will be used in future calculations
                entry.append(rate)
                entry.append(framecounts)

                # add the good point source, including location and counts, to the pointSourceList
                filteredList.add_row(entry)

        # Write the results to a file
        self.n_galaxies = len(filteredList)
        print(("Number of galaxies found within or close to the requested aperture: {}".format(self.n_galaxies)))

        if self.n_galaxies == 0:
            if self.n_pointsources == 0:
                raise ValueError(('No point sources or galaxies found in input catalog; empty seed image '
                                  'would be created.'))

        filteredList.meta['comments'] = [("Field center (degrees): %13.8f %14.8f y axis rotation angle "
                                          "(degrees): %f  image size: %4.4d %4.4d\n" %
                                          (self.ra, self.dec, self.params['Telescope']['rotation'], nx, ny))]
        filteredOut = self.basename + '_galaxySources.list'
        filteredList.write(filteredOut, format='ascii', overwrite=True)
        return filteredList

    def create_galaxy(self, radius, ellipticity, sersic, posang, totalcounts):
        """Create a model 2d sersic image with a given radius, eccentricity,
        position angle, and total counts.

        Parameters
        ----------
        radius : float
            Half light radius of the sersic profile, in units of pixels

        ellipticity : float
            Ellipticity of sersic profile

        sersic : float
            Sersic index

        posang : float
            Position angle in units of degrees

        totalcounts : float
            Total summed signal of the output image

        Returns
        -------
        img : numpy.ndarray
            2D array containing the 2D sersic profile
        """

        # create the grid of pixels
        meshmax = np.min([np.int(self.ffsize * self.coord_adjust['y']), np.int(radius * 100.)])

        # Make sure the grid has odd dimensions so that the galaxy will
        # be centered
        if meshmax % 2 == 0:
            meshmax += 1

        y, x = np.meshgrid(np.arange(meshmax), np.arange(meshmax))

        # Center the galaxy in the array
        xc = (meshmax // 2)
        yc = (meshmax // 2)

        # Create model
        mod = Sersic2D(amplitude=1, r_eff=radius, n=sersic, x_0=xc, y_0=yc,
                       ellip=ellipticity, theta=posang)

        # Create instance of model
        img = mod(x, y)

        # Scale such that the total number of counts in the galaxy matches the input
        summedcounts = np.sum(img)
        if summedcounts == 0:
            print('Zero counts in image in create_galaxy: ', radius, ellipticity, sersic, posang, totalcounts)
            factor = 0.
        else:
            factor = totalcounts / summedcounts
        img = img * factor

        # Crop image down such that it contains 99.95% of the total signal
        img = self.crop_galaxy_stamp(img, 0.9995)
        return img

    def crop_galaxy_stamp(self, stamp, threshold):
        """Crop an input stamp image containing a galaxy to a size that
        contains only threshold times the total signal. This is an
        attempt to speed up the simulator a bit, since the galaxy stamp
        images are often very large. Note that galaxy stamp images being
        fed into this function are currently always square.

        Parameters
        ----------
        stamp : numpy.ndarray
            2D stamp image of galaxy

        threshold : float
            Fraction of total flux to keep in the cropped image
            (e.g. 0.999 = 99.9%)

        Returns
        -------
        stamp : numpy.ndarray
            2D cropped image
        """
        totsignal = np.sum(stamp)
        # In the case of no signal, return the original stamp image.
        # This can happen for some Sersic profile parameters when the galaxy
        # is very small compared to the pixel size.
        if totsignal == 0.:
            return stamp
        yd, xd = stamp.shape
        mid = np.int(xd / 2)
        for rad in range(mid):
            signal = np.sum(stamp[mid - rad:mid + rad + 1, mid - rad:mid + rad + 1]) / totsignal
            if signal >= threshold:
                return stamp[mid - rad:mid + rad + 1, mid - rad:mid + rad + 1]
        # If we make it all the way through the stamp without
        # hitting the threshold, then return the full stamp image
        return stamp

    def make_galaxy_image(self, file):
        """Using the entries in the ``simSignals:galaxyList`` file, create a countrate image
        of model galaxies (2D sersic profiles)

        Parameters
        ----------
        catalog_file : str
            Name of ascii catalog file containing galaxy sources

        Returns
        -------
        galimage : numpy.ndarray
            2D array containing countrate image of galaxy sources

        segmentation.segmap : numpy.ndarray
            Segmentation map corresponding to ``galimage``
        """
        # Read in the list of galaxies (positions and magnitides)
        glist, pixflag, radflag, magsys = self.readGalaxyFile(file)
        if pixflag:
            print("Galaxy list input positions assumed to be in units of pixels.")
        else:
            print("Galaxy list input positions assumed to be in units of RA and Dec.")

        if radflag:
            print("Galaxy list input radii assumed to be in units of pixels.")
        else:
            print("Galaxy list input radii assumed to be in units of arcsec.")

        # Extract and save only the entries which will land (fully or partially) on the
        # aperture of the output
        galaxylist = self.filterGalaxyList(glist, pixflag, radflag, magsys, file)

        # galaxylist is a table with columns:
        # 'pixelx', 'pixely', 'RA', 'Dec', 'RA_degrees', 'Dec_degrees', 'radius', 'ellipticity',
        # 'pos_angle', 'sersic_index', 'magnitude', 'countrate_e/s', 'counts_per_frame_e'

        # final output image
        yd, xd = self.output_dims

        # create the final galaxy countrate image
        galimage = np.zeros((yd, xd))

        # Create corresponding segmentation map
        segmentation = segmap.SegMap()
        segmentation.xdim = xd
        segmentation.ydim = yd
        segmentation.initialize_map()

        # Create table of point source countrate versus psf size
        if self.add_psf_wings is True:
            self.translate_psf_table(magsys)

        # For each entry, create an image, and place it onto the final output image
        start_time = time.time()
        times = []
        time_reported = False
        for entry in galaxylist:
            # Warn user of how long this calcuation might take...
            if len(times) < 30:
                elapsed_time = time.time() - start_time
                times.append(elapsed_time)
                start_time = time.time()
            elif len(times) == 30 and time_reported is False:
                avg_time = np.mean(times)
                total_time = len(galaxylist) * avg_time
                print(("Expected time to process {} sources: {:.2f} seconds "
                       "({:.2f} minutes)".format(len(galaxylist), total_time, total_time / 60)))
                time_reported = True

            # Get position angle in the correct units. Inputs for each
            # source are degrees east of north. So we need to find the
            # angle between north and V3, and then the angle between
            # V3 and the y-axis on the detector. The former can be found
            # using rotations.posang(attitude_matrix, v2, v3). The latter
            # is just V3SciYAngle in the SIAF (I think???)
            # v3SciYAng is measured in degrees, from V3 towards the Y axis,
            # measured from V3 towards V2.
            xposang = self.calc_x_position_angle(entry['V2'], entry['V3'], entry['pos_angle'])

            # First create the galaxy
            stamp = self.create_galaxy(entry['radius'], entry['ellipticity'], entry['sersic_index'],
                                       xposang*np.pi/180., entry['countrate_e/s'])

            # If the stamp image is smaller than the PSF in either
            # dimension, embed the stamp in an array that matches
            # the psf size. This is so the upcoming convolution will
            # produce an output that includes the wings of the PSF
            galdims = stamp.shape

            # Using the PSF "core" normalized to 1 will keep more light near
            # the core of the galaxy, compared to the more rigorous
            # approach that uses the full convolution including the wings.
            # Whether this is a problem or not will depend on the relative
            # sizes of the photometry aperture versus the extended source.
            psf_dimensions = np.array(self.psf_library.data.shape[-2:])
            psf_shape = np.array((psf_dimensions / self.psf_library_oversamp) -
                                 self.params['simSignals']['gridded_psf_library_row_padding']).astype(np.int)

            if ((galdims[0] < psf_shape[0]) or (galdims[1] < psf_shape[1])):
                stamp = self.enlarge_stamp(stamp, psf_shape)
                galdims = stamp.shape

            # Get the PSF which will be convolved with the galaxy profile
            psf_image, min_x, min_y, wings_added = self.create_psf_stamp(entry['pixelx'], entry['pixely'],
                                                                         psf_shape[1], psf_shape[0], ignore_detector=True)

            # Skip sources that fall completely off the detector
            if psf_image is None:
                continue

            # Normalize the signal in the PSF stamp so that the final galaxy
            # signal will match the requested value
            psf_image = psf_image / np.sum(psf_image)

            # If the source subpixel location is beyond 0.5 (i.e. the edge
            # of the pixel), then we shift the wing->core offset by 1.
            # We also need to shift the location of the wing array on the
            # detector by 1
            if wings_added:
                x_delta = int(np.modf(entry['pixelx'])[0] > 0.5)
                y_delta = int(np.modf(entry['pixely'])[0] > 0.5)
            else:
                x_delta = 0
                y_delta = 0

            # Calculate the coordinates describing the overlap between
            # the PSF image and the galaxy image
            xap, yap, xpts, ypts, (i1, i2), (j1, j2), (k1, k2), \
                (l1, l2) = self.create_psf_stamp_coords(entry['pixelx']+x_delta, entry['pixely']+y_delta,
                                                        galdims, galdims[1] // 2, galdims[0] // 2,
                                                        coord_sys='aperture')

            # Make sure the stamp is at least partially on the detector
            if i1 is not None and i2 is not None and j1 is not None and j2 is not None:
                # Convolve the galaxy image with the PSF image
                stamp = s1.fftconvolve(stamp, psf_image, mode='same')

                # Now add the stamp to the main image
                if ((j2 > j1) and (i2 > i1) and (l2 > l1) and (k2 > k1) and (j1 < yd) and (i1 < xd)):
                    galimage[j1:j2, i1:i2] += stamp[l1:l2, k1:k2]
                    # Divide readnoise by 100 sec, which is a 10 group RAPID ramp
                    noiseval = self.single_ron / 100. + self.params['simSignals']['bkgdrate']
                    if self.params['Inst']['mode'].lower() in ['wfss', 'ts_grism']:
                        noiseval += self.grism_background
                    segmentation.add_object_noise(stamp[l1:l2, k1:k2], j1, i1, entry['index'], noiseval)

                else:
                    pass
                    # print("Source located entirely outside the field of view. Skipping.")
            else:
                pass

        return galimage, segmentation.segmap

    def calc_x_position_angle(self, v2_value, v3_value, position_angle):
        """Calcuate the position angle of the source relative to the x
        axis of the detector given the source's v2, v3 location and the
        user-input position angle (degrees east of north).

        Parameters
        ----------
        v2_value : float
            V2 location of source in units of arcseconds

        v3_value : float
            V3 location of source in units of arcseconds

        position_angle : float
            Position angle of source in degrees east of north

        Returns
        -------
        x_posang : float
            Position angle of source relative to detector x
            axis, in units of degrees
        """
        north_to_east_V3ang = rotations.posangle(self.attitude_matrix, v2_value, v3_value)
        x_posang = 0. - (self.siaf.V3SciXAngle - north_to_east_V3ang + self.local_roll - position_angle
                         + 90. + self.params['Telescope']['rotation'])
        return x_posang

    def getExtendedSourceList(self, filename):
        # read in the list of point sources to add, and adjust the
        # provided positions for astrometric distortion

        extSourceList = Table(names=('index', 'pixelx', 'pixely', 'RA', 'Dec',
                                     'RA_degrees', 'Dec_degrees', 'magnitude',
                                     'countrate_e/s', 'counts_per_frame_e'),
                              dtype=('i', 'f', 'f', 'S14', 'S14', 'f', 'f', 'f', 'f', 'f'))

        try:
            lines, pixelflag, magsys = self.read_point_source_file(filename)
            if pixelflag:
                print("Extended source list input positions assumed to be in units of pixels.")
            else:
                print("Extended list input positions assumed to be in units of RA and Dec.")
        except:
            raise FileNotFoundError("WARNING: Unable to open the extended source list file {}".format(filename))

        # File to save adjusted point source locations
        eoutcat = self.params['Output']['file'][0:-5] + '_extendedsources.list'
        eslist = open(eoutcat, 'w')

        dtor = math.radians(1.)
        nx = (self.subarray_bounds[2] - self.subarray_bounds[0]) + 1
        ny = (self.subarray_bounds[3] - self.subarray_bounds[1]) + 1
        xc = (self.subarray_bounds[2] + self.subarray_bounds[0]) / 2.
        yc = (self.subarray_bounds[3] + self.subarray_bounds[1]) / 2.

        # Write out the RA and Dec of the field center to the output file
        # Also write out column headers to prepare for source list
        eslist.write(("# Field center (degrees): %13.8f %14.8f y axis rotation angle "
                      "(degrees): %f  image size: %4.4d %4.4d\n" %
                      (self.ra, self.dec, self.params['Telescope']['rotation'], nx, ny)))
        eslist.write('# \n')
        eslist.write(("#    Index   RA_(hh:mm:ss)   DEC_(dd:mm:ss)   RA_degrees      "
                      "DEC_degrees     pixel_x   pixel_y    magnitude   counts/sec    counts/frame\n"))

        # Add an index column if not present
        indexes = self.get_index_numbers(lines)

        # Check the source list and remove any sources that are well outside the
        # field of view of the detector. These sources cause the coordinate
        # conversion to hang.
        indexes, lines = self.remove_outside_fov_sources(indexes, lines, pixelflag, 4096)

        print("After extended sources, max index is {}".format(self.maxindex))

        # Determine the name of the column to use for source magnitudes
        mag_column = self.select_magnitude_column(lines, filename)

        # Loop over input lines in the source list
        all_stamps = []
        for indexnum, values in zip(indexes, lines):
            if not os.path.isfile(values['filename']):
                raise FileNotFoundError('{} from extended source catalog does not exist.'.format(values['filename']))

            pixelx, pixely, ra, dec, ra_str, dec_str = self.get_positions(values['x_or_RA'],
                                                                          values['y_or_Dec'],
                                                                          pixelflag, 4096)
            # Get the input magnitude
            try:
                mag = float(values[mag_column])
            except ValueError:
                mag = None

            # Now find out how large the extended source image is, so we
            # know if all, part, or none of it will fall in the field of view
            ext_stamp = fits.getdata(values['filename'])
            if len(ext_stamp.shape) != 2:
                ext_stamp = fits.getdata(values['filename'], 1)

            # print('Extended source location, x, y, ra, dec:', pixelx, pixely, ra, dec)
            # print('extended source size:', ext_stamp.shape)

            # Rotate the stamp image if requested
            ext_stamp = self.rotate_extended_image(ext_stamp, values['pos_angle'], ra, dec)

            eshape = np.array(ext_stamp.shape)
            if len(eshape) == 2:
                edgey, edgex = eshape // 2
            else:
                raise ValueError(("WARNING, extended source image {} is not 2D! "
                                  "This is not supported.".format(values['filename'])))

            # Define the min and max source locations (in pixels) that fall onto the subarray
            # Inlude the effects of a requested grism_direct image, and also keep sources that
            # will only partially fall on the subarray
            # pixel coords here can still be negative and kept if the grism image is being made

            # First, coord limits for just the subarray
            miny = 0
            maxy = self.subarray_bounds[3] - self.subarray_bounds[1]
            minx = 0
            maxx = self.subarray_bounds[2] - self.subarray_bounds[0]

            # Expand the limits if a grism direct image is being made
            if (self.params['Output']['grism_source_image'] == True) or (self.params['Inst']['mode'] in ["pom", "wfss"]):
                extrapixy = np.int((maxy + 1)/2 * (self.coord_adjust['y'] - 1.))
                miny -= extrapixy
                maxy += extrapixy
                extrapixx = np.int((maxx + 1)/2 * (self.coord_adjust['x'] - 1.))
                minx -= extrapixx
                maxx += extrapixx

            # Now, expand the dimensions again to include point sources that fall only partially on the
            # subarray
            miny -= edgey
            maxy += edgey
            minx -= edgex
            maxx += edgex

            # Keep only sources within the appropriate bounds
            if pixely > miny and pixely < maxy and pixelx > minx and pixelx < maxx:

                # Set up an entry for the output table
                entry = [indexnum, pixelx, pixely, ra_str, dec_str, ra, dec, mag]

                # save the stamp image after normalizing to a total signal of 1.
                norm_factor = np.sum(ext_stamp)
                ext_stamp /= norm_factor
                all_stamps.append(ext_stamp)

                # If a magnitude is given then adjust the countrate to match it
                if mag is not None:
                    # Convert magnitudes to countrate (ADU/sec) and counts per frame
                    countrate = utils.magnitude_to_countrate(self.params['Inst']['mode'], magsys, mag,
                                                             photfnu=self.photfnu, photflam=self.photflam,
                                                             vegamag_zeropoint=self.vegazeropoint)
                    framecounts = countrate * self.frametime
                    magwrite = mag

                else:
                    # In this case, no magnitude is given in the extended input list
                    # Assume the input stamp image is in units of e/sec then.
                    print("No magnitude given for extended source in {}.".format(values['filename']))
                    print("Assuming the original file is in units of counts per sec.")
                    print("Multiplying original file values by 'extendedscale'.")
                    countrate = norm_factor * self.params['simSignals']['extendedscale']
                    framecounts = countrate * self.frametime
                    magwrite = 99.99999

                # add the countrate and the counts per frame to pointSourceList
                # since they will be used in future calculations
                # entry.append(scale)
                entry.append(countrate)
                entry.append(framecounts)

                # add the good point source, including location and counts, to the pointSourceList
                # self.pointSourceList.append(entry)
                extSourceList.add_row(entry)

                # Write out positions, distances, and counts to the output file
                eslist.write(("%i %s %s %14.8f %14.8f %9.3f %9.3f  %9.3f  %13.6e   %13.6e\n" %
                             (indexnum, ra_str, dec_str, ra, dec, pixelx, pixely, magwrite, countrate,
                              framecounts)))

        print("Number of extended sources found within or close to the requested aperture: {}".format(len(extSourceList)))
        # close the output file
        eslist.close()

        # If no good point sources were found in the requested array, alert the user
        if len(extSourceList) < 1:
            print("Warning: no non-sidereal extended sources within the requested array.")
            print("The extended source image option is being turned off")

        return extSourceList, all_stamps

    def rotate_extended_image(self, stamp_image, pos_angle, right_ascention, declination):
        """Given the user-input position angle for the extended source
        image, calculate the appropriate angle of the stamp image
        relative to the detector x axis, and rotate the stamp image.
        TO DO: if the stamp image contains a WCS, use that to
        determine rotation angle

        Parameters
        ----------
        stamp_image : numpy.ndarray
            2D stamp image of the extended source

        pos_angle : float
            Position angle of stamp image relative to north in degrees

        right_ascention : float
            RA of source, in decimal degrees

        declination : float
            Dec of source, in decimal degrees

        Returns
        -------
        rotated : numpy.ndarray
            Rotated stamp image
        """
        # Add later: check for WCS and use that
        # if no WCS:
        pixelv2, pixelv3 = pysiaf.utils.rotations.getv2v3(self.attitude_matrix, right_ascention, declination)
        x_pos_ang = self.calc_x_position_angle(pixelv2, pixelv3, pos_angle)

        # The definition of position angle is different between the galaxy
        # source inputs and the extended source inputs. calc_x_position_angle
        # was created for galaxies, so adjust for the different definition here.
        rotation_angle = x_pos_ang - 2. * self.params['Telescope']['rotation']
        rotated = rotate(stamp_image, rotation_angle, mode='constant', cval=0.)
        return rotated

    def make_extended_source_image(self, extSources, extStamps):
        # Create the empty image
        yd, xd = self.output_dims
        extimage = np.zeros(self.output_dims)

        # Create corresponding segmentation map
        segmentation = segmap.SegMap()
        segmentation.xdim = xd
        segmentation.ydim = yd
        segmentation.initialize_map()

        # Loop over the entries in the source list
        for entry, stamp in zip(extSources, extStamps):
            stamp_dims = stamp.shape

            stamp *= entry['countrate_e/s']

            # If the stamp needs to be convolved with the NIRCam PSF,
            # create the correct PSF  here and read it in
            if self.params['simSignals']['PSFConvolveExtended']:
                # If the stamp image is smaller than the PSF in either
                # dimension, embed the stamp in an array that matches
                # the psf size. This is so the upcoming convolution will
                # produce an output that includes the wings of the PSF
                psf_dimensions = np.array(self.psf_library.data.shape[-2:])
                psf_shape = np.array((psf_dimensions / self.psf_library_oversamp) -
                                     self.params['simSignals']['gridded_psf_library_row_padding']).astype(np.int)
                if ((stamp_dims[0] < psf_shape[0]) or (stamp_dims[1] < psf_shape[1])):
                    stamp = self.enlarge_stamp(stamp, psf_shape)
                    stamp_dims = stamp.shape

                # Create the PSF
                # Using the PSF "core" normalized to 1 will keep more light near
                # the core of the galaxy, compared to the more rigorous
                # approach that uses the full convolution including the wings.
                # Whether this is a problem or not will depend on the relative
                # sizes of the photometry aperture versus the extended source.
                psf_image, min_x, min_y, wings_added = self.create_psf_stamp(entry['pixelx'], entry['pixely'],
                                                                             psf_shape[1], psf_shape[0], ignore_detector=True)

                # Skip sources that fall completely off the detector
                if psf_image is None:
                    continue

                # Normalize the PSF so that the final signal in the extended
                # source mathces the requested signal
                psf_image = psf_image / np.sum(psf_image)

                # If the source subpixel location is beyond 0.5 (i.e. the edge
                # of the pixel), then we shift the wing->core offset by 1.
                # We also need to shift the location of the wing array on the
                # detector by 1
                if wings_added:
                    x_delta = int(np.modf(entry['pixelx'])[0] > 0.5)
                    y_delta = int(np.modf(entry['pixely'])[0] > 0.5)
                else:
                    x_delta = 0
                    y_delta = 0

                # Calculate the coordinates describing the overlap
                # between the extended image and the PSF image
                xap, yap, xpts, ypts, (i1, i2), (j1, j2), (k1, k2), \
                    (l1, l2) = self.create_psf_stamp_coords(entry['pixelx']+x_delta, entry['pixely']+y_delta,
                                                            stamp_dims, stamp_dims[1] // 2, stamp_dims[0] // 2,
                                                            coord_sys='aperture')

                if None in [i1, i2, j1, j2, k1, k2, l1, l2]:
                    continue

                # Convolve the extended image with the stamp image
                stamp = s1.fftconvolve(stamp, psf_image, mode='same')
            else:
                # If no PSF convolution is to be done, find the
                # coordinates describing the overlap between the
                # original stamp image and the aperture
                xap, yap, xpts, ypts, (i1, i2), (j1, j2), (k1, k2), \
                    (l1, l2) = self.create_psf_stamp_coords(entry['pixelx'], entry['pixely'],
                                                            stamp_dims, stamp_dims[1] // 2, stamp_dims[0] // 2,
                                                            coord_sys='aperture')

            # Make sure the stamp is at least partially on the detector
            if i1 is not None and i2 is not None and j1 is not None and j2 is not None:

                # Now add the stamp to the main image
                if ((j2 > j1) and (i2 > i1) and (l2 > l1) and (k2 > k1) and (j1 < yd) and (i1 < xd)):
                    extimage[j1:j2, i1:i2] += stamp[l1:l2, k1:k2]

                # Divide readnoise by 100 sec, which is a 10 group RAPID ramp?
                noiseval = self.single_ron / 100. + self.params['simSignals']['bkgdrate']
                if self.params['Inst']['mode'].lower() in ['wfss', 'ts_grism']:
                    noiseval += self.grism_background

                # Make segmentation map
                indseg = self.seg_from_photutils(stamp[l1:l2, k1:k2] * entry['countrate_e/s'],
                                                 entry['index'], noiseval)
                segmentation.segmap[j1:j2, i1:i2] += indseg
        return extimage, segmentation.segmap

    def enlarge_stamp(self, image, dims):
        """Place the given image within an enlarged array of zeros. If the
        requested dimension lengths are odd while ``image``'s dimension
        lengths are even, then the new array is expanded by one to also have
        even dimensions. This ensures that ``image`` will be centered
        within ``array``.

        Parameters
        ----------
        image : numpy.ndarray
            2D image

        dims : list
            2-element list of (y-dimension, x-dimension) to embed ``image``
            within

        Returns
        -------
        array : numpy.ndarray
            Expanded image
        """
        dim_y, dim_x = dims
        image_size_y, image_size_x = image.shape
        if dim_y % 2 != image_size_y % 2:
            dim_y += 1
        if dim_x % 2 != image_size_x % 2:
            dim_x += 1

        array = np.zeros((dim_y, dim_x))

        dx = dim_x - image_size_x
        dx = np.int(dx / 2)

        dy = dim_y - image_size_y
        dy = np.int(dy / 2)

        array[dy:dim_y-dy, dx:dim_x-dx] = image
        return array

    def seg_from_photutils(self, image, number, noise):
        # Create a segmentation map for the input image
        # using photutils. In this case, the input noise
        # represents the single frame noise for the appropriate
        # observing mode
        map = detect_sources(image, noise * 3., 8).data + number
        return map

    def makeFilterTable(self):
        # Create the table that contains the possible filter list, quantum yields, and countrates for a
        # star with vega magnitude of 15 in each filter. Do this by reading in phot_file
        # listed in the parameter file.

        # FUTURE WORK: If the countrates are left as 0, then pysynphot will
        # be used to calculate them later
        try:
            cvals_tab = ascii.read(self.params['Reffiles']['phot'])
            instrumentfilternames = cvals_tab['filter'].data
            stringcountrates = cvals_tab['countrate_for_vegamag15'].data
            instrumentmag15countrates = [float(s) for s in stringcountrates]
            strinstrumentqy = cvals_tab['quantum_yield'].data
            qy = [float(s) for s in strinstrumentqy]
            self.countvalues = dict(zip(instrumentfilternames, instrumentmag15countrates))
            self.qydict = dict(zip(instrumentfilternames, qy))
        except:
            raise IOError("WARNING: Unable to read in {}.".format(self.params['Reffiles']['phot']))

    def readParameterFile(self):
        """Read in the parameter file"""
        # List of fields in the yaml file to check for extra colons
        search_cats = ['title:', 'PI_Name:', 'Science_category:', 'observation_label:']

        # Open the yaml file and check for the presence of extra colons
        adjust_file = False
        with open(self.paramfile) as infile:
            read_data = infile.readlines()
            for i, line in enumerate(read_data):
                for search_term in search_cats:
                    if search_term in line:
                        idx = []
                        hashidx = [200]
                        for m in re.finditer(':', line):
                            idx.append(m.start())
                        for mm in re.finditer('#', line):
                            hashidx.append(mm.start())
                        num = np.sum(np.array(idx) < min(hashidx))
                        if num > 1:
                            adjust_file = True
                            later_string = line[idx[0]+1:]
                            later_string = later_string.replace(':', ',')
                            newline = line[0: idx[0]+1] + later_string
                            read_data[i] = newline

        if adjust_file:
            # Make a copy of the original file and then delete it
            yaml_copy = 'orig_{}'.format(self.paramfile)
            shutil.copy2(self.paramfile, yaml_copy)
            os.remove(self.paramfile)

            # Write the adjusted lines to a new copy of the input file
            with open(self.paramfile, 'w') as f:
                for item in read_data:
                    f.write("{}".format(item))

        # Load the yaml file
        try:
            with open(self.paramfile, 'r') as infile:
                self.params = yaml.safe_load(infile)
        except (ScannerError, FileNotFoundError, IOError) as e:
            print(e)

    def check_params(self):
        """Check input parameters for expected datatypes, values"""
        # Check instrument name
        if self.params['Inst']['instrument'].lower() not in INST_LIST:
            raise NotImplementedError("WARNING: {} instrument not implemented within ramp simulator")

        # Check entered mode:
        possibleModes = MODES[self.params['Inst']['instrument'].lower()]
        self.params['Inst']['mode'] = self.params['Inst']['mode'].lower()
        if self.params['Inst']['mode'] in possibleModes:
            pass
        else:
            raise ValueError(("WARNING: unrecognized mode {} for {}. Must be one of: {}"
                              .format(self.params['Inst']['mode'],
                                      self.params['Inst']['instrument'], possibleModes)))

        # Check telescope tracking entry
        self.params['Telescope']['tracking'] = self.params['Telescope']['tracking'].lower()
        if self.params['Telescope']['tracking'] not in TRACKING_LIST:
            raise ValueError(("WARNING: telescope tracking set to {}, but must be one "
                              "of {}.".format(self.params['Telescope']['tracking'],
                                              TRACKING_LIST)))

        # Non-sidereal WFSS observations are not yet supported
        if self.params['Telescope']['tracking'] == 'non-sidereal' and \
           self.params['Inst']['mode'] in ['wfss', 'ts_grism']:
            raise ValueError(("WARNING: wfss observations with non-sidereal "
                              "targets not yet supported."))

        # Set nframe and nskip according to the values in the
        # readout pattern definition file
        self.read_pattern_check()

        # Make sure that the requested number of groups is
        # less than or equal to the maximum allowed.
        # For full frame science operations, ngroup is going
        # to be limited to 10 for all readout patterns
        # except for the DEEP patterns, which can go to 20.
        # match = self.readpatterns['name'] == self.params['Readout']['readpatt'].upper()
        # if sum(match) == 1:
        #    maxgroups = self.readpatterns['maxgroups'].data[match][0]
        # if sum(match) == 0:
        #    print("Unrecognized readout pattern {}. Assuming a maximum allowed number of groups of 10.".format(self.params['Readout']['readpatt']))
        #    maxgroups = 10

        # if (self.params['Readout']['ngroup'] > maxgroups):
        #    print("WARNING: {} is limited to a maximum of {} groups. Proceeding with ngroup = {}.".format(self.params['Readout']['readpatt'], maxgroups, maxgroups))
        #    self.params['Readout']['readpatt'] = maxgroups

        # Check for entries in the parameter file that are None or blank,
        # indicating the step should be skipped. Create a dictionary of steps
        # and populate with True or False
        self.runStep = {}
        self.runStep['pixelflat'] = self.checkRunStep(self.params['Reffiles']['pixelflat'])
        self.runStep['illuminationflat'] = self.checkRunStep(self.params['Reffiles']['illumflat'])
        self.runStep['astrometric'] = self.checkRunStep(self.params['Reffiles']['astrometric'])
        # self.runStep['distortion_coeffs'] = self.checkRunStep(self.params['Reffiles']['distortion_coeffs'])
        self.runStep['ipc'] = self.checkRunStep(self.params['Reffiles']['ipc'])
        self.runStep['crosstalk'] = self.checkRunStep(self.params['Reffiles']['crosstalk'])
        self.runStep['occult'] = self.checkRunStep(self.params['Reffiles']['occult'])
        self.runStep['pointsource'] = self.checkRunStep(self.params['simSignals']['pointsource'])
        self.runStep['galaxies'] = self.checkRunStep(self.params['simSignals']['galaxyListFile'])
        self.runStep['extendedsource'] = self.checkRunStep(self.params['simSignals']['extended'])
        self.runStep['movingTargets'] = self.checkRunStep(self.params['simSignals']['movingTargetList'])
        self.runStep['movingTargetsSersic'] = self.checkRunStep(self.params['simSignals']['movingTargetSersic'])
        self.runStep['movingTargetsExtended'] = self.checkRunStep(self.params['simSignals']['movingTargetExtended'])
        self.runStep['MT_tracking'] = self.checkRunStep(self.params['simSignals']['movingTargetToTrack'])
        self.runStep['zodiacal'] = self.checkRunStep(self.params['simSignals']['zodiacal'])
        self.runStep['scattered'] = self.checkRunStep(self.params['simSignals']['scattered'])
        # self.runStep['fwpw'] = self.checkRunStep(self.params['Reffiles']['filtpupilcombo'])
        self.runStep['pixelAreaMap'] = self.checkRunStep(self.params['Reffiles']['pixelAreaMap'])

        # Notify user if no catalogs are provided
        if self.params['simSignals']['pointsource'] == 'None':
            print('No point source catalog provided in yaml file.')

        if self.params['simSignals']['galaxyListFile'] == 'None':
            print('No galaxy catalog provided in yaml file.')

        # Determine the instrument module and detector from the aperture name
        aper_name = self.params['Readout']['array_name']
        self.instrument = self.params["Inst"]["instrument"].lower()
        try:
            # previously detector was e.g. 'A1'. Let's make it NRCA1 to be more in
            # line with other instrument formats
            # detector = self.subdict[self.subdict['AperName'] == aper_name]['Detector'][0]
            # module = detector[0]
            detector = aper_name.split('_')[0]
            self.detector = detector
            if self.instrument == 'nircam':
                module = detector[3]
            elif self.instrument == 'niriss':
                module = detector[0]
            elif self.instrument == 'fgs':
                detector = detector.replace("FGS", "GUIDER")
                module = detector[0]
        except IndexError:
            raise ValueError('Unable to determine the detector/module in aperture {}'.format(aper_name))

        if self.instrument == 'niriss':
            newfilter, newpupil = utils.check_niriss_filter(self.params['Readout']['filter'],self.params['Readout']['pupil'])
            self.params['Readout']['filter'] = newfilter
            self.params['Readout']['pupil'] = newpupil

        # Make sure the requested filter is allowed. For imaging, all filters
        # are allowed. In the future, other modes will be more restrictive
        if self.params['Readout']['pupil'][0].upper() == 'F':
            self.usefilt = 'pupil'
        else:
            self.usefilt = 'filter'

        # If instrument is FGS, then force filter to be 'NA' for the purposes
        # of constructing the correct PSF input path name. Then change to be
        # the DMS-required "N/A" when outputs are saved
        if self.instrument == 'fgs':
            self.params['Readout']['filter'] = 'NA'
            self.params['Readout']['pupil'] = 'NA'

        # Get basic flux calibration information
        self.vegazeropoint, self.photflam, self.photfnu, self.pivot = \
            fluxcal_info(self.params, self.usefilt, detector, module)

        # Convert the input RA and Dec of the pointing position into floats
        # Check to see if the inputs are in decimal units or hh:mm:ss strings
        try:
            self.ra = float(self.params['Telescope']['ra'])

            self.dec = float(self.params['Telescope']['dec'])
        except:
            self.ra, self.dec = utils.parse_RA_Dec(self.params['Telescope']['ra'],
                                                   self.params['Telescope']['dec'])

        #if abs(self.dec) > 90. or self.ra < 0. or self.ra > 360. or \
        if abs(self.dec) > 90. or \
           self.ra is None or self.dec is None:
            raise ValueError("WARNING: bad requested RA and Dec {} {}".format(self.ra, self.dec))

        # make sure the rotation angle is a float
        try:
            self.params['Telescope']["rotation"] = float(self.params['Telescope']["rotation"])
        except ValueError:
            print(("ERROR: bad rotation value {}, setting to zero."
                   .format(self.params['Telescope']["rotation"])))
            self.params['Telescope']["rotation"] = 0.

        siaf_inst = self.params['Inst']['instrument']
        if siaf_inst.lower() == 'nircam':
            siaf_inst = 'NIRCam'
        instrument_siaf = siaf_interface.get_instance(siaf_inst)
        self.siaf = instrument_siaf[self.params['Readout']['array_name']]
        self.local_roll, self.attitude_matrix, self.ffsize, \
            self.subarray_bounds = siaf_interface.get_siaf_information(instrument_siaf,
                                                                       self.params['Readout']['array_name'],
                                                                       self.ra, self.dec,
                                                                       self.params['Telescope']['rotation'])

        print('SIAF: Requested {}   got {}'.format(self.params['Readout']['array_name'], self.siaf.AperName))
        # Set the background value if the high/medium/low settings
        # are used
        bkgdrate_options = ['high', 'medium', 'low']

        # For WFSS observations, we want the background in the direct
        # seed image to be zero. The dispersed background will be created
        # and added as part of the dispersion process
        if self.params['Inst']['mode'].lower() == 'wfss':
            self.params['simSignals']['bkgdrate'] = 0.

        if np.isreal(self.params['simSignals']['bkgdrate']):
            self.params['simSignals']['bkgdrate'] = float(self.params['simSignals']['bkgdrate'])
        else:
            if self.params['simSignals']['bkgdrate'].lower() in bkgdrate_options:
                print(("Calculating background rate using jwst_background "
                       "based on {} level".format(self.params['simSignals']['bkgdrate'])))

                # Find the appropriate filter throughput file
                if os.path.split(self.params['Reffiles']['filter_throughput'])[1] == 'placeholder.txt':
                    filter_file = utils.get_filter_throughput_file(self.params['Inst']['instrument'].lower(),
                                                                   self.params['Readout'][self.usefilt],
                                                                   fgs_detector=detector, nircam_module=module)
                else:
                    filter_file = self.params['Reffiles']['filter_throughput']

                print(("Using {} filter throughput file for background calculation."
                       .format(filter_file)))

                # To translate background signals from MJy/sr to e-/sec to
                # ADU/sec, we need a mean gain value
                if self.params['Inst']['instrument'].lower() == 'nircam':
                    if '5' in detector:
                        shorthand = 'lw{}'.format(module.lower())
                    else:
                        shorthand = 'sw{}'.format(module.lower())
                    gain_value = MEAN_GAIN_VALUES[self.params['Inst']['instrument'].lower()][shorthand]
                elif self.params['Inst']['instrument'].lower() == 'niriss':
                    gain_value = MEAN_GAIN_VALUES[self.params['Inst']['instrument'].lower()]
                elif self.params['Inst']['instrument'].lower() == 'fgs':
                    gain_value = MEAN_GAIN_VALUES[self.params['Inst']['instrument'].lower()][detector]

                if self.params['simSignals']['use_dateobs_for_background']:
                    bkgd_wave, bkgd_spec = backgrounds.day_of_year_background_spectrum(self.params['Telescope']['ra'],
                                                                                       self.params['Telescope']['dec'],
                                                                                       self.params['Output']['date_obs'])
                    self.params['simSignals']['bkgdrate'] = backgrounds.calculate_background(self.ra, self.dec,
                                                                                             filter_file, True,
                                                                                             gain_value, self.siaf,
                                                                                             back_wave=bkgd_wave,
                                                                                             back_sig=bkgd_spec)
                    print("Background rate determined using date_obs: {}".format(self.params['Output']['date_obs']))
                else:
                    # Here the background level is based on high/medium/low rather than date
                    orig_level = copy.deepcopy(self.params['simSignals']['bkgdrate'])
                    self.params['simSignals']['bkgdrate'] = backgrounds.calculate_background(self.ra, self.dec,
                                                                                             filter_file, False,
                                                                                             gain_value, self.siaf,
                                                                                             level=self.params['simSignals']['bkgdrate'].lower())
                    print("Background rate determined using {} level: {}".format(orig_level, self.params['simSignals']['bkgdrate']))
            else:
                raise ValueError(("WARNING: unrecognized background rate value. "
                                  "Must be either a number or one of: {}"
                                  .format(bkgdrate_options)))

        # Check that the various scaling factors are floats and within a reasonable range
        # self.params['cosmicRay']['scale'] = self.checkParamVal(self.params['cosmicRay']['scale'], 'cosmicRay', 0, 100, 1)
        self.params['simSignals']['extendedscale'] = self.checkParamVal(self.params['simSignals']['extendedscale'],
                                                                        'extendedEmission', 0, 10000, 1)
        self.params['simSignals']['zodiscale'] = self.checkParamVal(self.params['simSignals']['zodiscale'],
                                                                    'zodi', 0, 10000, 1)
        self.params['simSignals']['scatteredscale'] = self.checkParamVal(self.params['simSignals']['scatteredscale'],
                                                                         'scatteredLight', 0, 10000, 1)

        # make sure the requested output format is an allowed value
        if self.params['Output']['format'] not in ALLOWEDOUTPUTFORMATS:
            raise ValueError(("WARNING: unsupported output format {} requested. "
                              "Possible options are {}.".format(self.params['Output']['format'],
                                                                ALLOWEDOUTPUTFORMATS)))

        # Entries for creating the grism input image
        if not isinstance(self.params['Output']['grism_source_image'], bool):
            if self.params['Output']['grism_source_image'].lower() == 'none':
                self.params['Output']['grism_source_image'] = False
            else:
                raise ValueError("WARNING: grism_source_image needs to be True or False")

        # Location of extended image on output array, pixel x, y values.
        try:
            self.params['simSignals']['extendedCenter'] = np.fromstring(self.params['simSignals']['extendedCenter'],
                                                                        dtype=int, sep=", ")
        except:
            raise RuntimeError(("WARNING: not able to parse the extendedCenter list {}. "
                                "It should be a comma-separated list of x and y pixel positions."
                                .format(self.params['simSignals']['extendedCenter'])))

    def checkRunStep(self, filename):
        # check to see if a filename exists in the parameter file.
        if ((len(filename) == 0) or (filename.lower() == 'none')):
            return False
        else:
            return True

    def read_pattern_check(self):
        # Check the readout pattern that's entered and set nframe and nskip
        # accordingly
        self.params['Readout']['readpatt'] = self.params['Readout']['readpatt'].upper()

        # Read in readout pattern definition file
        # and make sure the possible readout patterns are in upper case
        self.readpatterns = ascii.read(self.params['Reffiles']['readpattdefs'])
        self.readpatterns['name'] = [s.upper() for s in self.readpatterns['name']]

        # If the requested readout pattern is in the table of options,
        # then adopt the appropriate nframe and nskip
        if self.params['Readout']['readpatt'] in self.readpatterns['name']:
            mtch = self.params['Readout']['readpatt'] == self.readpatterns['name']
            self.params['Readout']['nframe'] = self.readpatterns['nframe'][mtch].data[0]
            self.params['Readout']['nskip'] = self.readpatterns['nskip'][mtch].data[0]
            print(('Requested readout pattern {} is valid. '
                  'Using the nframe = {} and nskip = {}'
                   .format(self.params['Readout']['readpatt'],
                           self.params['Readout']['nframe'],
                           self.params['Readout']['nskip'])))
        else:
            # If the read pattern is not present in the definition file
            # then quit.
            raise ValueError(("WARNING: the {} readout pattern is not defined in {}."
                              .format(self.params['Readout']['readpatt'],
                                      self.params['Reffiles']['readpattdefs'])))

    def filecheck(self):
        # Make sure the requested input files exist
        # For reference files, assume first that they are located in
        # the directory tree under the directory specified by the MIRAGE_DATA
        # environment variable. If not, assume the input is a full path
        # and check there.
        rlist = [['Reffiles', 'astrometric']]
        plist = [['simSignals', 'psfpath']]
        ilist = [['simSignals', 'pointsource'],
                 ['simSignals', 'galaxyListFile'],
                 ['simSignals', 'extended'],
                 ['simSignals', 'movingTargetList'],
                 ['simSignals', 'movingTargetSersic'],
                 ['simSignals', 'movingTargetExtended'],
                 ['simSignals', 'movingTargetToTrack']]
        # for ref in rlist:
        #    self.ref_check(ref)
        for path in plist:
            self.path_check(path)
        for inp in ilist:
            self.input_check(inp)

    def ref_check(self, rele):
        """
        Check for the existence of the input reference file
        Assume first that the file is in the directory tree
        specified by the MIRAGE_DATA environment variable.

        Parameters:
        -----------
        rele -- Tuple containing the nested keys that point
                to the refrence file of interest. These come
                from the yaml input file

        Reutrns:
        --------
        Nothing
        """
        rfile = self.params[rele[0]][rele[1]]
        if rfile.lower() != 'none':
            c1 = os.path.isfile(rfile)
            if not c1:
                raise FileNotFoundError(("WARNING: Unable to locate the {}, {} "
                                         "input file! Not present in {}"
                                         .format(rele[0], rele[1], rfile)))

    def path_check(self, p):
        """
        Check for the existence of the input path.
        Assume first that the path is in relation to
        the directory tree specified by the MIRAGE_DATA
        environment variable

        Parameters:
        -----------
        p -- Tuple containing the nested keys that point
             to a directory in self.params

        Returns:
        --------
        Nothing
        """
        pth = self.params[p[0]][p[1]]
        c1 = os.path.exists(pth)
        if not c1:
            raise NotADirectoryError(("WARNING: Unable to find the requested path "
                                      "{}. Not present in directory tree specified by "
                                      "the {} environment variable."
                                      .format(pth, self.env_var)))

    def input_check(self, inparam):
        # Check for the existence of the input file. In
        # this case we do not check the directory tree
        # specified by the MIRAGE_DATA environment variable.
        # This is intended primarily for user-generated inputs like
        # source catalogs
        ifile = self.params[inparam[0]][inparam[1]]
        if ifile.lower() != 'none':
            c = os.path.isfile(ifile)
            if not c:
                raise FileNotFoundError(("WARNING: Unable to locate {} Specified "
                                         "by the {}:{} field in the input yaml file."
                                         .format(ifile, inparam[0], inparam[1])))

    def checkParamVal(self, value, typ, vmin, vmax, default):
        # make sure the input value is a float and between min and max
        try:
            value = float(value)
        except:
            raise ValueError("WARNING: {} for {} is not a float.".format(value, typ))

        if ((value >= vmin) & (value <= vmax)):
            return value
        else:
            print(("ERROR: {} for {} is not within reasonable bounds. "
                   "Setting to {}".format(value, typ, default)))
            return default

    def read_distortion_reffile(self):
        """Read in the CRDS-format distortion reference file and save
        the coordinate transformation model
        """
        coord_transform = None
        if self.runStep['astrometric']:
            with asdf.open(self.params['Reffiles']['astrometric']) as dist_file:
                coord_transform = dist_file.tree['model']
        # else:
        #    coord_transform = self.simple_coord_transform()
        return coord_transform

    """
    def calculate_background(self, ra, dec, ffile, back_wave=None, back_sig=None, level='medium'):
        '''Use the JWST background calculator to come up with
        an appropriate background level for the observation.
        Options for level include low, medium, high'''
        from jwst_backgrounds import jbt
        from astropy import units as u
        from astropy.units.equivalencies import si, cgs

        # Read in filter throughput file
        filt_wav, filt_thru = file_io.read_filter_throughput(ffile)

        # If the user wants a background signal from a particular day,
        # then extract that array here
        if self.params['simSignals']['use_dateobs_for_background']:
            # Interpolate background to match filter wavelength grid
            bkgd_interp = np.interp(filt_wav, back_wave, back_sig)

            # Combine
            filt_bkgd = bkgd_interp * filt_thru

            # Integrate over the filter bandpass
            filt_integ = np.trapz(filt_thru, x=filt_wav)

            # Integrate
            bval = np.trapz(filt_bkgd, x=filt_wav) / filt_integ * u.MJy / u.steradian

        else:
            # If the user has requested background in terms of low/medium/high,
            # then we need to examine all the background arrays.
            # Loop over each day (in the background)
            # info, convolve the background curve with the filter
            # throughput curve, and then integrate. THEN, we can
            # calculate the low/medium/high values.
            bval = backgrounds.low_medium_high_background_value(ra, dec, level, filt_wav, filt_thru)
            bval = bval * u.MJy / u.steradian

        # Convert from MJy/str to ADU/sec
        # then divide by area of pixel
        flambda = cgs.erg / si.angstrom / si.cm ** 2 / si.s
        # fnu = cgs.erg / si.Hz / si.cm ** 2 / si.s
        photflam = self.photflam * flambda
        # photnu = self.photfnu * fnu
        pivot = self.pivot * u.micron
        mjy = photflam.to(u.MJy, u.spectral_density(pivot))

        # Divide by pixel area in steradians to get
        # MJy/str per ADU/s
        pixel_area = self.siaf.XSciScale * u.arcsec * self.siaf.YSciScale * u.arcsec
        mjy_str = mjy / pixel_area.to(u.steradian)

        # Convert the background signal from MJy/str
        # to ADU/sec
        bval /= mjy_str
        return bval.value
        """

    def saveSingleFits(self, image, name, key_dict=None, image2=None, image2type=None):
        # Save an array into the first extension of a fits file
        h0 = fits.PrimaryHDU()
        h1 = fits.ImageHDU(image, name='DATA')
        if image2 is not None:
            h2 = fits.ImageHDU(image2)
            if image2type is not None:
                h2.header['EXTNAME'] = image2type

        # if a keyword dictionary is provided, put the
        # keywords into the 0th and 1st extension headers
        if key_dict is not None:
            for key in key_dict:
                h0.header[key] = key_dict[key]
                h1.header[key] = key_dict[key]

        if image2 is None:
            hdulist = fits.HDUList([h0, h1])
        else:
            hdulist = fits.HDUList([h0, h1, h2])
        hdulist.writeto(name, overwrite=True)

    def add_options(self, parser=None, usage=None):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage, description='Create seed image via catalogs')
        parser.add_argument("paramfile", help=('File describing the input parameters and instrument '
                                               'settings to use. (YAML format).'))
        parser.add_argument("--param_example", help='If used, an example parameter file is output.')
        return parser


if __name__ == '__main__':

    usagestring = 'USAGE: catalog_seed_image.py inputs.yaml'

    seed = Catalog_seed()
    parser = seed.add_options(usage=usagestring)
    args = parser.parse_args(namespace=seed)
    seed.make_seed()
