# ! /usr/bin/env python

'''
Class to produce yaml files that can be used as input for
the ramp simulator

Use
---
    To generate an observationlist.yaml file and generate all .yamls at once
    ::
        yam = yaml_generator.SimInput(input_xml=apt_file_xml, pointing_file=apt_file_pointing,
                                      catalogs=catalogs, verbose=True, output_dir=out_dir,
                                      simdata_output_dir=out_dir)
        yam.create_inputs()

    NOTE there is currently no way to generate yamls from an existing observationlist.yaml


Inputs
------
    xml file - Name of xml file exported from APT.
    pointing file - Name of associated pointing file exported from APT.

    Optional inputs:

    output_dir - Directory into which the output yaml files are written

    simdata_output_dir - Directory to place in the output_directory field of the yaml files.
                         This is the directory where the simulated ramps will be saved.

    table_file - Ascii table containing observation info. This is the
                 output from apt_inputs.py Use this if you are
                 not providing xml and pointing files from APT.

    datatype - Specifies the type of output data to save. Can be "raw", in which
               case the raw (uncalibrated) file is saved, "linear", where the
               linearized file is saved and ready to be run through the jump
               detection and ramp-fitting steps of the pipeline, or
               "linear, raw", where both versions are saved.

    use_nonstsci_names - set to True to override the use of the standard
                         STScI naming convention for output files

    subarray_def_file - Ascii file containing NIRCam subarray definitions

    readpatt_def_file - Ascii file containing NIRCam readout pattern definitions

    point_source - point source catalog file. Can be a single file, or a list of
                   catalogs. If it is a list, each filename is expected to contain
                   the filter name for which it is to be used. Catalogs and filters
                   will then be matched up in the output yaml files.

    galaxyListFile - galaxy (sersic) source catalog file. Can be a single name, or
                     a list of names. Behavior is identical to point_source above.

    extended - extended source catalog file. Behavior is identical to point_source
                above.

    convolveExtended - Set to True to convolve extended sources with NIRCam PSF

    movingTarg - Moving (point source) target catalog (sources moving through fov) names.
                 Behavior is the same as point_sources above.

    movingTargSersic - Moving galaxy (sersic) target catalog (sources moving through fov)
                       names. Behavior is the same as point_sources above.

    movingTargExtended - Moving extended source target catalog (sources moving through fov)
                         names. Behavior is the same as point_sources above.

    movingTargToTrack - Catalog of non-sidereal targets for non-sidereal tracking observations.
                        Behavior is the same as point_sources above.

    bkgdrate - Uniform background rate (e-/s) to add to observation.

    epoch_list - Ascii table file containing epoch start times and telescope roll angles
                 to use for each observation.


Dependencies
------------
    argparse, astropy, numpy, glob, copy
    apt_inputs.py - Functions for reading and parsing xml and pointing files from APT.

History
-------
    July 2017 - V0: Initial version. Bryan Hilbert
    Feb 2018 - V1: Updates to accomodate multiple filter pairs per
                   observation. Lauren Chambers
    August 2018 - V2: Replaced input SIAF CSV file with pysiaf dependence. Bryan Hilbert

'''

import sys
import os
import argparse
from glob import glob
from copy import deepcopy

from astropy.time import Time, TimeDelta
from astropy.table import Table
from astropy.io import ascii, fits
import numpy as np
import pkg_resources
import pysiaf

from ..apt import apt_inputs
from ..utils.utils import calc_frame_time, ensure_dir_exists, expand_environment_variable
from .generate_observationlist import get_observation_dict
from ..constants import NIRISS_PUPIL_WHEEL_ELEMENTS, NIRISS_FILTER_WHEEL_ELEMENTS

ENV_VAR = 'MIRAGE_DATA'


class SimInput:
    def __init__(self, input_xml=None, pointing_file=None, datatype='linear',
                 use_JWST_pipeline=True, catalogs=None, observation_list_file=None, verbose=False,
                 output_dir='./', simdata_output_dir='./', parameter_defaults=None, offline=False):
        """Initialize instance. Read APT xml and pointing files if provided.

        Also sets the reference files definitions for all instruments.

        Parameters
        ----------
        input_xml : str
            Filename of xml file exported by APT
        pointing_file : str
            Filename of pointing file exported by APT
        datatype : str
        use_JWST_pipeline : bool

        parameter_defaults : dict
            Default values of parameters like roll angle (PAV3) to pass on to observation list
            generator
        """
        self.info = {}
        self.input_xml = input_xml
        self.pointing_file = pointing_file
        self.datatype = datatype
        self.use_JWST_pipeline = use_JWST_pipeline
        self.catalogs = catalogs
        self.observation_list_file = observation_list_file
        self.verbose = verbose
        self.output_dir = output_dir
        self.simdata_output_dir = simdata_output_dir

        self.table_file = None
        self.use_nonstsci_names = False
        # self.subarray_def_file = 'config'
        # self.readpatt_def_file = 'config'
        # self.crosstalk = 'config'
        # self.filtpupil_pairs = 'config'
        # self.fluxcal = 'config'
        # self.dq_init_config = 'config'
        # self.saturation_config = 'config'
        # self.superbias_config = 'config'
        # self.refpix_config = 'config'
        # self.linearity_config = 'config'
        # self.filter_throughput = 'config'
        # self.observation_list_file = None
        self.use_linearized_darks = True
        self.psfwfe = 'predicted'
        self.psfwfegroup = 0
        self.resets_bet_ints = 1  # NIRCam should be 1
        self.tracking = 'sidereal'
        self.psf_paths = None
        self.expand_catalog_for_segments = False
        self.add_psf_wings = True

        # Expand the MIRAGE_DATA environment variable
        self.datadir = expand_environment_variable(ENV_VAR, offline=offline)

        # Get the path to the 'MIRAGE' package
        self.modpath = pkg_resources.resource_filename('mirage', '')

        self.set_global_definitions()
        self.path_defs()

        if (input_xml is not None) and (catalogs is not None):
            if self.observation_list_file is None:
                self.observation_list_file = os.path.join(self.output_dir, 'observation_list.yaml')
            self.apt_xml_dict = get_observation_dict(self.input_xml, self.observation_list_file, self.catalogs,
                                                     verbose=self.verbose, parameter_defaults=parameter_defaults)

        self.reffile_setup(offline=offline)

    def add_catalogs(self):
        """
        Add list(s) of source catalogs to the table containing the
        observation information
        """
        self.info['point_source'] = [None] * len(self.info['Module'])
        self.info['galaxyListFile'] = [None] * len(self.info['Module'])
        self.info['extended'] = [None] * len(self.info['Module'])
        self.info['convolveExtended'] = [False] * len(self.info['Module'])
        self.info['movingTarg'] = [None] * len(self.info['Module'])
        self.info['movingTargSersic'] = [None] * len(self.info['Module'])
        self.info['movingTargExtended'] = [None] * len(self.info['Module'])
        self.info['movingTargToTrack'] = [None] * len(self.info['Module'])

        for i in range(len(self.info['ShortFilter'])):
            if np.int(self.info['detector'][i][-1]) < 5:
                filtkey = 'ShortFilter'
                pupilkey = 'ShortPupil'
            else:
                filtkey = 'LongFilter'
                pupilkey = 'LongPupil'
            filt = self.info[filtkey][i]
            pup = self.info[pupilkey][i]

            if self.point_source[i] is not None:
                # In here, we assume the user provided a catalog to go with each filter
                # so now we need to find the filter for each entry and generate a list that makes sense
                self.info['point_source'][i] = os.path.abspath(os.path.expandvars(
                    self.catalog_match(filt, pup, self.point_source, 'point source')))
            else:
                self.info['point_source'][i] = None
            if self.galaxyListFile[i] is not None:
                self.info['galaxyListFile'][i] = os.path.abspath(os.path.expandvars(
                    self.catalog_match(filt, pup, self.galaxyListFile, 'galaxy')))
            else:
                self.info['galaxyListFile'][i] = None
            if self.extended[i] is not None:
                self.info['extended'][i] = os.path.abspath(os.path.expandvars(
                    self.catalog_match(filt, pup, self.extended, 'extended')))
            else:
                self.info['extended'][i] = None
            if self.movingTarg[i] is not None:
                self.info['movingTarg'][i] = os.path.abspath(os.path.expandvars(
                    self.catalog_match(filt, pup, self.movingTarg, 'moving point source target')))
            else:
                self.info['movingTarg'][i] = None
            if self.movingTargSersic[i] is not None:
                self.info['movingTargSersic'][i] = os.path.abspath(os.path.expandvars(
                    self.catalog_match(filt, pup, self.movingTargSersic, 'moving sersic target')))
            else:
                self.info['movingTargSersic'][i] = None
            if self.movingTargExtended[i] is not None:
                self.info['movingTargExtended'][i] = os.path.abspath(os.path.expandvars(
                    self.catalog_match(filt, pup, self.movingTargExtended, 'moving extended target')))
            else:
                self.info['movingTargExtended'][i] = None
            if self.movingTargToTrack[i] is not None:
                self.info['movingTargToTrack'][i] = os.path.abspath(os.path.expandvars(
                    self.catalog_match(filt, pup, self.movingTargToTrack, 'non-sidereal moving target')))
            else:
                self.info['movingTargToTrack'][i] = None
        if self.convolveExtended is True:
            self.info['convolveExtended'] = [True] * len(self.info['Module'])

    def catalog_match(self, filter, pupil, catalog_list, cattype):
        """
        Given a filter and pupil value, along with a list of input
        catalogs, find the catalog names that contain each filter/
        pupil name.

        Parameters
        ----------
        filter : str
            Name of a filter element
        pupil : str
            Name of a pupil element
        catalog_list : list
            List of catalog filenames
        cattype : str
            Type of catalog in the list.

        Returns
        -------
        match : str
            Name of catalog that contains the name of the
            input filter/pupil element
        """
        if pupil[0].upper() == 'F':
            match = [s for s in catalog_list if pupil.lower() in s.lower()]
            if len(match) == 0:
                self.no_catalog_match(pupil, cattype)
                return None
            elif len(match) > 1:
                self.multiple_catalog_match(pupil, cattype, match)
            return match[0]
        else:
            match = [s for s in catalog_list if filter.lower() in s.lower()]
            if len(match) == 0:
                self.no_catalog_match(filter, cattype)
                return None
            elif len(match) > 1:
                self.multiple_catalog_match(filter, cattype, match)
            return match[0]

    def create_inputs(self):
        """Create observation table """
        self.path_defs()

        if ((self.input_xml is not None) &
           (self.pointing_file is not None) &
           (self.observation_list_file is not None)):

            # Define directories and paths
            indir, infile = os.path.split(self.input_xml)
            final_file = os.path.join(self.output_dir,
                                      'Observation_table_for_' + infile +
                                      '_with_yaml_parameters.csv')

            # Read XML file and make observation table
            apt = apt_inputs.AptInput(input_xml=self.input_xml, pointing_file=self.pointing_file,
                                      output_dir=self.output_dir)
            # apt.input_xml = self.input_xml
            # apt.pointing_file = self.pointing_file
            apt.observation_list_file = self.observation_list_file
            apt.apt_xml_dict = self.apt_xml_dict

            apt.output_dir = self.output_dir
            apt.create_input_table()

            self.info = apt.exposure_tab

            # Add start time info to each element
            self.make_start_times()

            # Add a list of output yaml names to the dictionary
            self.make_output_names()

        elif self.table_file is not None:
            print('Reading table file: {}'.format(self.table_file))
            info = ascii.read(self.table_file)
            self.info = self.table_to_dict(info)
            final_file = self.table_file + '_with_yaml_parameters.csv'

        else:
            raise FileNotFoundError(("WARNING. You must include either an ascii table file of observations"
                                     " or xml and pointing files from APT plus the observation list file."
                                     "Aborting."))

        # For each element in the lists below, use the detector name to
        # find the appropriate reference files. Create lists, and add
        # these lists to the dictionary
        darks = []
        lindarks = []
        superbias = []
        linearity = []
        saturation = []
        gain = []
        astrometric = []
        ipc = []
        pam = []

        # detector_labels = self.info['detector']
        # for det in detector_labels:
        for instrument, det in zip([s.lower() for s in self.info['Instrument']], self.info['detector']):
            instrument = instrument.lower()
            darks.append(self.get_dark(instrument, det))
            lindarks.append(self.get_lindark(instrument, det))
            superbias.append(self.get_reffile(self.superbias_list[instrument], det))
            linearity.append(self.get_reffile(self.linearity_list[instrument], det))
            saturation.append(self.get_reffile(self.saturation_list[instrument], det))
            gain.append(self.get_reffile(self.gain_list[instrument], det))
            astrometric.append(self.get_reffile(self.astrometric_list[instrument], det))
            ipc.append(self.get_reffile(self.ipc_list[instrument], det))
            pam.append(self.get_reffile(self.pam_list[instrument], det))

        self.info['dark'] = darks
        # If linearized darks are to be used, set the darks to None
        if self.use_linearized_darks:
            self.info['dark'] = [None] * len(darks)
            self.info['lindark'] = lindarks
        else:
            self.info['dark'] = darks
            self.info['lindark'] = [None] * len(lindarks)
        self.info['superbias'] = superbias
        self.info['linearity'] = linearity
        self.info['saturation'] = saturation
        self.info['gain'] = gain
        self.info['astrometric'] = astrometric
        self.info['ipc'] = ipc
        self.info['pixelAreaMap'] = pam

        # Add setting describing whether JWST pipeline will be used
        self.info['use_JWST_pipeline'] = [self.use_JWST_pipeline] * len(darks)

        # add background rate to the table
        # self.info['bkgdrate'] = np.array([self.bkgdrate]*len(self.info['Mode']))

        # grism entries
        grism_source_image = ['False'] * len(self.info['Mode'])
        grism_input_only = ['False'] * len(self.info['Mode'])
        for i in range(len(self.info['Mode'])):
            if self.info['Mode'][i].lower() == 'wfss':
                if self.info['detector'][i] == 'NIS' or self.info['detector'][i][-1] == '5':
                    grism_source_image[i] = 'True'
                    grism_input_only[i] = 'True'
                # SW detectors shouldn't be wfss
                if self.info['Instrument'][i] == 'NIRCAM' and self.info['detector'][i][-1] != '5':
                    self.info['Mode'][i] = 'imaging'
        self.info['grism_source_image'] = grism_source_image
        self.info['grism_input_only'] = grism_input_only

        # level-3 associated keywords that are not present in APT file.
        # not quite sure how to populate these
        self.info['visit_group'] = ['01'] * len(self.info['Mode'])
        # self.info['sequence_id'] = ['1'] * len(self.info['Mode'])
        seq = []
        for par in self.info['CoordinatedParallel']:
            if par.lower() == 'true':
                seq.append('2')
            if par.lower() == 'false':
                seq.append('1')
        self.info['sequence_id'] = seq
        # self.info['obs_template'] = ['NIRCam Imaging'] * len(self.info['Mode'])

        # Deal with user-provided PSFs that differ across observations/visits/exposures
        self.info['psfpath'] = self.get_psf_path()

        # write out the updated table, including yaml filenames
        # start times, and reference files
        #if 0: #for debugging
        #    for key in self.info.keys():
        #        print('{:>40} has {:>3} entries'.format(key, len(self.info[key])))
        table = Table(self.info)
        table.write(final_file, format='csv', overwrite=True)
        # ascii.write(Table(self.info), final_file, format='csv', overwrite=True)
        print('Updated observation table file saved to {}'.format(final_file))

        # Now go through the lists one element at a time
        # and create a yaml file for each.
        yamls = []
        # for i in range(len(detector_labels)):
        for i, instrument in enumerate(self.info['Instrument']):
            instrument = instrument.lower()
            if instrument not in 'fgs nircam niriss'.split():
                # do not write files for MIRI and NIRSpec
                continue
            file_dict = {}
            for key in self.info:
                file_dict[key] = self.info[key][i]

            # break dither number into numbers for primary
            # and subpixel dithers
            tot_dith = np.int(file_dict['dither'])
            # primarytot = np.int(file_dict['PrimaryDithers'])
            primarytot = np.int(file_dict['number_of_dithers'])

            if file_dict['SubpixelPositions'].upper() == 'NONE':
                subpixtot = 1
            else:
                try:
                    subpixtot = np.int(file_dict['SubpixelPositions'])
                except:
                    subpixtot = np.int(file_dict['SubpixelPositions'][0])
            primary_dither = np.ceil(1. * tot_dith / subpixtot)
            file_dict['primary_dither_num'] = np.int(primary_dither)
            subpix_dither = (tot_dith-1) % subpixtot
            file_dict['subpix_dither_num'] = subpix_dither + 1
            file_dict['subarray_def_file'] = self.global_subarray_definition_files[instrument]
            file_dict['readpatt_def_file'] = self.global_readout_pattern_files[instrument]
            file_dict['crosstalk_file'] = self.global_crosstalk_files[instrument]
            file_dict['filtpupilcombo_file'] = self.global_filtpupilcombo_files[instrument]
            file_dict['flux_cal_file'] = self.global_flux_cal_files[instrument]
            file_dict['psf_wing_threshold_file'] = self.global_psf_wing_threshold_file[instrument]

            fname = self.write_yaml(file_dict)
            yamls.append(fname)

        # Write out summary of all written yaml files
        filenames = [y.split('/')[-1] for y in yamls]
        mosaic_numbers = sorted(list(set([f.split('_')[0] for f in filenames])))
        obs_ids = sorted(list(set([m[7:10] for m in mosaic_numbers])))

        print('\n')

        total_exposures = 0
        for obs in obs_ids:
            visit_list = list(set([m[10:] for m in mosaic_numbers
                                   if m[7:10] == obs]))
            n_visits = len(visit_list)
            activity_list = list(set([vf[17:29] for vf in self.info['yamlfile']
                                      if vf[7:10] == obs]))
            n_activities = len(activity_list)
            exposure_list = list(set([vf[20:25] for vf in self.info['yamlfile']
                                      if vf[7:10] == obs]))
            n_exposures = len(exposure_list) * n_activities
            total_exposures += n_exposures
            all_obs_files = [m for m in self.info['yamlfile'] if m[7:10] == obs]
            total_files = len(all_obs_files)

            obs_id_int = np.array([np.int(ele) for ele in self.info['ObservationID']])
            obs_indexes = np.where(obs_id_int == np.int(obs))[0]
            obs_entries = np.array(self.info['ParallelInstrument'])[obs_indexes]
            coord_par = self.info['CoordinatedParallel'][obs_indexes[0]]
            if coord_par:
                par_indexes = np.where(obs_entries)[0]
                if len(par_indexes) > 0:
                    parallel_instrument = self.info['Instrument'][obs_indexes[par_indexes[0]]]
                else:
                    parallel_instrument = 'NONE'
            else:
                parallel_instrument = 'NONE'

            # Note that changing the == below to 'is' results in an error
            pri_indexes = np.where(obs_entries == False)[0]
            prime_instrument = self.info['Instrument'][obs_indexes[pri_indexes[0]]]

            if prime_instrument.upper() == 'NIRCAM':
                module = self.info['Module'][obs_indexes[pri_indexes[0]]]
            elif parallel_instrument.upper() == 'NIRCAM':
                module = self.info['Module'][obs_indexes[par_indexes[0]]]

            if ((prime_instrument.upper() == 'NIRCAM') or (parallel_instrument.upper() == 'NIRCAM')):
                if module in ['A', 'B']:
                    n_det = 5
                if module == 'ALL':
                    n_det = 10
                    module = 'A and B'
            else:
                # number of detectors
                n_det = 1

            if coord_par == 'true':
                instrument_string = '    Prime: {}, Parallel: {}'.format(prime_instrument, parallel_instrument)
            else:
                instrument_string = '    Prime: {}'.format(prime_instrument)

            print('Observation {}:'.format(obs))
            print(instrument_string)
            print('    {} visit(s)'.format(n_visits))
            print('    {} activitiy(ies)'.format(n_activities))
            print('    {} exposure(s)'.format(n_exposures))
            if ((prime_instrument.upper() == 'NIRCAM') or (parallel_instrument.upper() == 'NIRCAM')):
                print('    {} NIRCam detector(s) in module {}'.format(n_det, module))
            print('    {} file(s)'.format(total_files))

        print('\n{} exposures total.'.format(total_exposures))
        print('{} output files written to: {}'.format(len(yamls), self.output_dir))

    def create_output_name(self, input_obj, index=0):
        """Put together the JWST formatted fits file name based on observation parameters

        Parameters
        ----------
        input_obj : dict
            Dictionary from apt_inputs giving information for each exposure

        Returns
        -------
        base : str
            JWST formatted filename base (excluding pipeline step suffix and ".fits")
        """
        proposal_id = '{0:05d}'.format(int(input_obj['ProposalID'][index]))
        observation = input_obj['obs_num'][index]
        visit_number = input_obj['visit_num'][index]
        visit_group = input_obj['visit_group'][index]
        parallel_sequence_id = input_obj['sequence_id'][index]
        activity_id = input_obj['act_id'][index]
        exposure = input_obj['exposure'][index]

        base = 'jw{}{}{}_{}{}{}_{}_'.format(proposal_id, observation, visit_number,
                                            visit_group, parallel_sequence_id, activity_id,
                                            exposure)
        return base

    def find_ipc_file(self, inputipc):
        """Given a list of potential IPC kernel files for a given
        detector, select the most appropriate one, and check to see
        whether the kernel needs to be inverted, in order to populate
        the invertIPC field. This is not intended to be terribly smart.
        The first inverted kernel found will be used. If none are found,
        the first kernel will be used and set to be inverted.

        Parameters
        ----------
        inputipc : list
           List of fits files containing IPC kernels for a single detector

        Returns
        -------
        (ipcfile, invstatus) : tup
           ipcfile is the name of the IPC kernel file to use, and invstatus
           lists whether the kernel needs to be inverted or not.
        """
        for ifile in inputipc:
            kernel = fits.getdata(ifile)
            kshape = kernel.shape

            # If kernel is 4 dimensional, extract the 3x3 kernel associated
            # with a single pixel
            if len(kernel.shape) == 4:
                kernel = kernel[:, :, np.int(kshape[2]/2), np.int(kshape[2]/2)]

            if kernel[1, 1] < 1.0:
                return (ifile, False)
        # If no inverted kernel was found, just return the first file
        return (inputipc[0], True)

    def get_dark(self, instrument, detector):
        """Return the name of a dark current file to use as input
        based on the detector being used

        Parameters
        ----------
        detector : str
            Name of detector being used
        Returns
        -------
        files : str
            Name of a dark current file to use for this detector
        """
        files = self.dark_list[instrument][detector]
        if len(files) == 1:
            return files[0]
        else:
            rand_index = np.random.randint(0, len(files) - 1)
            return files[rand_index]

    def get_lindark(self, instrument, detector):
        """
        Return the name of a linearized dark current file to
        use as input based on the detector being used

        Parameters
        ----------
        detector : str
            Name of detector being used

        Returns
        -------
        files : str
            Name of a linearized dark current file to use for this detector
        """
        files = self.lindark_list[instrument][detector]
        if len(files) == 1:
            return files[0]
        else:
            rand_index = np.random.randint(0, len(files) - 1)
            return files[rand_index]

    def get_readpattern_defs(self, filename=None):
        """Read in the readpattern definition file and return table.

        Returns
        -------
        tab : obj
            astropy.table.Table containing readpattern definitions
        filename : str
            Path to input file name
        """
        if filename is not None:
            return ascii.read(filename)

        tab = ascii.read(self.readpatt_def_file)
        return tab

    def get_reffile(self, refs, detector):
        """
        Return the appropriate reference file for detector
        and given reference file dictionary.

        Parameters
        ----------
        refs : dict
            dictionary in the form of:
             {'A1':'filenamea1.fits', 'A2':'filenamea2.fits'...}
             Containing reference file names
        detector : str
            Name of detector

        Returns
        -------
        refs: str
            Name of reference file appropriate for given detector
        """
        for key in refs:
            if detector in key:
                return refs[key]
        print("WARNING: no file found for detector {} in {}"
              .format(detector, refs))

    def get_subarray_defs(self, filename=None):
        """Read in subarray definition file and return table

        Returns
        -------
        sub : obj
            astropy.table.Table containing subarray definition information
        filename : str
            Path to input file name
        """
        if filename is not None:
            return ascii.read(filename)

        sub = ascii.read(self.subarray_def_file)
        return sub

    def make_output_names(self):
        """
        Create output yaml file names to go with all of the
        entries in the dictionary
        """
        yaml_names = []
        fits_names = []

        if self.use_nonstsci_names:
            for i in range(len(self.info['Module'])):
                act = str(self.info['act_id'][i]).zfill(2)
                if self.info['Instrument'][i].lower() == 'niriss':
                    det = 'NIS'
                elif self.info['Instrument'][i].lower() == 'fgs':
                    det = 'FGS'
                else:
                    det = self.info['detector'][i]
                mode = self.info['Mode'][i]
                dither = str(self.info['dither'][i]).zfill(2)

                yaml_names.append(os.path.abspath(os.path.join(self.output_dir, 'Act{}_{}_{}_Dither{}.yaml'
                                                                            .format(act, det, mode, dither))))
                fits_names.append('Act{}_{}_{}_Dither{}_uncal.fits'.format(act, det, mode, dither))

        else:
            for i in range(len(self.info['Module'])):
                if self.info['Instrument'][i].upper() == 'NIRCAM':
                    fulldetector = 'nrc{}'.format(self.info['detector'][i].lower())
                else:
                    fulldetector = self.info['detector'][i].lower()
                outfilebase = self.create_output_name(self.info, index=i)
                outfile = "{}{}{}".format(outfilebase, fulldetector, '_uncal.fits')
                yamlout = "{}{}{}".format(outfilebase, fulldetector, '.yaml')

                yaml_names.append(yamlout)
                fits_names.append(outfile)

        self.info['yamlfile'] = yaml_names
        self.info['outputfits'] = fits_names

    def set_global_definitions(self):
        """Store the subarray definitions of all supported instruments."""
        # TODO: Investigate how this could be combined with the creation of
        #  self.configfiles in reffile_setup()

        self.global_subarray_definitions = {}
        self.global_readout_patterns = {}
        self.global_subarray_definition_files = {}
        self.global_readout_pattern_files = {}

        self.global_crosstalk_files = {}
        self.global_filtpupilcombo_files = {}
        self.global_flux_cal_files = {}
        self.global_psf_wing_threshold_file = {}
        self.global_psfpath = {}
        # self.global_filter_throughput_files = {} ?

        for instrument in 'niriss fgs nircam miri nirspec'.split():
            if instrument.lower() == 'niriss':
                readout_pattern_file = 'niriss_readout_pattern.txt'
                subarray_def_file = 'niriss_subarrays.list'
                crosstalk_file = 'niriss_xtalk_zeros.txt'
                filtpupilcombo_file = 'niriss_dual_wheel_list.txt'
                flux_cal_file = 'niriss_zeropoints.list'
                psf_wing_threshold_file = 'niriss_psf_wing_rate_thresholds.txt'
                psfpath = os.path.join(self.datadir, 'niriss/gridded_psf_library')
            elif instrument.lower() == 'fgs':
                readout_pattern_file = 'guider_readout_pattern.txt'
                subarray_def_file = 'guider_subarrays.list'
                crosstalk_file = 'guider_xtalk_zeros.txt'
                filtpupilcombo_file = 'guider_filter_dummy.list'
                flux_cal_file = 'guider_zeropoints.list'
                psf_wing_threshold_file = 'fgs_psf_wing_rate_thresholds.txt'
                psfpath = os.path.join(self.datadir, 'fgs/gridded_psf_library')
            elif instrument.lower() == 'nircam':
                readout_pattern_file = 'nircam_read_pattern_definitions.list'
                subarray_def_file = 'NIRCam_subarray_definitions.list'
                crosstalk_file = 'xtalk20150303g0.errorcut.txt'
                filtpupilcombo_file = 'nircam_filter_pupil_pairings.list'
                flux_cal_file = 'NIRCam_zeropoints.list'
                psf_wing_threshold_file = 'nircam_psf_wing_rate_thresholds.txt'
                psfpath = os.path.join(self.datadir, 'nircam/gridded_psf_library')
            else:
                readout_pattern_file = 'N/A'
                subarray_def_file = 'N/A'
                crosstalk_file = 'N/A'
                filtpupilcombo_file = 'N/A'
                flux_cal_file = 'N/A'
                psf_wing_threshold_file = 'N/A'
                psfpath = 'N/A'
            if instrument in 'niriss fgs nircam'.split():
                self.global_subarray_definitions[instrument] = self.get_subarray_defs(filename=os.path.join(self.modpath, 'config', subarray_def_file))
                self.global_readout_patterns[instrument] = self.get_readpattern_defs(filename=os.path.join(self.modpath, 'config', readout_pattern_file))
            self.global_subarray_definition_files[instrument] = os.path.join(self.modpath, 'config', subarray_def_file)
            self.global_readout_pattern_files[instrument] = os.path.join(self.modpath, 'config', readout_pattern_file)
            self.global_crosstalk_files[instrument] = os.path.join(self.modpath, 'config', crosstalk_file)
            self.global_filtpupilcombo_files[instrument] = os.path.join(self.modpath, 'config', filtpupilcombo_file)
            self.global_flux_cal_files[instrument] = os.path.join(self.modpath, 'config', flux_cal_file)
            self.global_psf_wing_threshold_file[instrument] = os.path.join(self.modpath, 'config', psf_wing_threshold_file)
            self.global_psfpath[instrument] = psfpath

    def make_start_times(self):
        """Create exposure start times for each entry in the observation dictionary."""
        date_obs = []
        time_obs = []
        expstart = []
        nframe = []
        nskip = []
        namp = []

        # choose arbitrary start time for each epoch
        epoch_base_time = '16:44:12'
        epoch_base_time0 = deepcopy(epoch_base_time)

        if 'epoch_start_date' in self.info.keys():
            epoch_base_date = self.info['epoch_start_date'][0]
        else:
            epoch_base_date = self.info['Date'][0]
        base = Time(epoch_base_date + 'T' + epoch_base_time)
        base_date, base_time = base.iso.split()

        # Pick some arbirary overhead values
        act_overhead = 90  # seconds. (filter change)
        visit_overhead = 600  # seconds. (slew)

        # Get visit, activity_id info for first exposure
        actid = self.info['act_id'][0]
        visit = self.info['visit_num'][0]
        # obsname = self.info['obs_label'][0]

        # for i in range(len(self.info['Module'])):
        for i, instrument in enumerate(self.info['Instrument']):
            # Get dither/visit
            # Files with the same activity_id should have the same start time
            # Overhead after a visit break should be large, smaller between
            # exposures within a visit
            next_actid = self.info['act_id'][i]
            next_visit = self.info['visit_num'][i]
            next_obsname = self.info['obs_label'][i]

            # Find the readpattern of the file
            readpatt = self.info['ReadoutPattern'][i]
            groups = np.int(self.info['Groups'][i])
            integrations = np.int(self.info['Integrations'][i])

            if instrument.lower() in ['miri', 'nirspec']:
                nframe.append(0)
                nskip.append(0)
                namp.append(0)
                date_obs.append(base_date)
                time_obs.append(base_time)
                expstart.append(base.mjd)

            else:
                # Now read in readpattern definitions
                readpatt_def = self.global_readout_patterns[instrument.lower()]

                # Read in file containing subarray definitions
                subarray_def = self.global_subarray_definitions[instrument.lower()]

                match2 = readpatt == readpatt_def['name']
                if np.sum(match2) == 0:
                    raise RuntimeError(("WARNING!! Readout pattern {} not found in definition file."
                                        .format(readpatt)))

                # Now get nframe and nskip so we know how many frames in a group
                fpg = np.int(readpatt_def['nframe'][match2][0])
                spg = np.int(readpatt_def['nskip'][match2][0])
                nframe.append(fpg)
                nskip.append(spg)

                # Get the aperture name. For non-NIRCam instruments,
                # this is simply the self.info['aperture']. But for NIRCam,
                # we need to be careful of entries like NRCBS_FULL, which is used
                # for observations using all 4 shortwave B detectors. In that case,
                # we need to build the aperture name from the combination of detector
                # and subarray name.
                # if np.all(np.unique(self.info['Instrument']) == 'NIRISS'):
                #     # aperture = self.info['aperture']
                #
                # elif np.all(np.unique(self.info['Instrument']) == 'NIRCAM'):
                #     sub = self.info['Subarray'][i]
                #     det = 'NRC' + self.info['detector'][i]
                #     aperture = det + '_' + sub

                aperture = self.info['aperture'][i]
                if 'NRC' == aperture[0:3]:
                    sub = self.info['Subarray'][i]
                    det = 'NRC' + self.info['detector'][i]
                    aperture = det + '_' + sub



                # Get the number of amps from the subarray definition file
                match = aperture == subarray_def['AperName']

                # needed for NIRCam case
                if np.sum(match) == 0:
                    aperture = [apername for apername, name in
                                np.array(subarray_def['AperName', 'Name']) if
                                (sub in apername) or (sub in name)]

                    match = aperture == subarray_def['AperName']

                    if len(aperture) > 1 or len(aperture) == 0 or np.sum(match) == 0:
                        raise ValueError('Cannot combine detector {} and subarray {}\
                            into valid aperture name.'.format(det, sub))
                    # We don't want aperture as a list
                    aperture = aperture[0]

                amp = subarray_def['num_amps'][match][0]
                namp.append(amp)

                # same activity ID
                if next_actid == actid:
                    # in this case, the start time should remain the same
                    date_obs.append(base_date)
                    time_obs.append(base_time)
                    expstart.append(base.mjd)
                    continue

                epoch_date = self.info['epoch_start_date'][i]
                epoch_time = deepcopy(epoch_base_time0)
                # new epoch - update the base time
                if epoch_date != epoch_base_date:
                    epoch_base_date = deepcopy(epoch_date)
                    base = Time(epoch_base_date + 'T' + epoch_base_time)
                    base_date, base_time = base.iso.split()
                    basereset = True
                    date_obs.append(base_date)
                    time_obs.append(base_time)
                    expstart.append(base.mjd)
                    actid = deepcopy(next_actid)
                    visit = deepcopy(next_visit)
                    obsname = deepcopy(next_obsname)
                    continue

                # new visit
                if next_visit != visit:
                    # visit break. Larger overhead
                    overhead = visit_overhead
                elif ((next_actid > actid) & (next_visit == visit)):
                    # same visit, new activity. Smaller overhead
                    overhead = act_overhead
                else:
                    # should never get in here
                    raise NotImplementedError()

                # For cases where the base time needs to change
                # continue down here
                siaf_inst = self.info['Instrument'][i].upper()
                if siaf_inst == 'NIRCAM':
                    siaf_inst = "NIRCam"
                siaf_obj = pysiaf.Siaf(siaf_inst)[aperture]

                # Calculate the readout time for a single frame
                frametime = calc_frame_time(siaf_inst, aperture,
                                            siaf_obj.XSciSize, siaf_obj.YSciSize, amp)

                # Estimate total exposure time
                exptime = ((fpg + spg) * groups + fpg) * integrations * frametime

                # Delta should include the exposure time, plus overhead
                delta = TimeDelta(exptime + overhead, format='sec')
                base += delta
                base_date, base_time = base.iso.split()

                # Add updated dates and times to the list
                date_obs.append(base_date)
                time_obs.append(base_time)
                expstart.append(base.mjd)

                # increment the activity ID and visit
                actid = deepcopy(next_actid)
                visit = deepcopy(next_visit)
                obsname = deepcopy(next_obsname)

        self.info['date_obs'] = date_obs
        self.info['time_obs'] = time_obs
        # self.info['expstart'] = expstart
        self.info['nframe'] = nframe
        self.info['nskip'] = nskip
        self.info['namp'] = namp

    def multiple_catalog_match(self, filter, cattype, matchlist):
        """
        Alert the user if more than one catalog matches the filter/pupil

        Parameters
        ----------
        filter : str
          Name of filter element
        cattype : str
          Type of catalog (e.g. pointsource)
        matchlist : list
          Matching catalog names
        """
        print("WARNING: multiple {} catalogs matched! Using the first.".format(cattype))
        print("Observation filter: {}".format(filter))
        print("Matched point source catalogs: {}".format(matchlist))

    def no_catalog_match(self, filter, cattype):
        """
        Alert user if no catalog match was found.

        Parameters
        ----------
        filter : str
          Name of filter element
        cattype : str
          Type of catalog (e.g. pointsource)

        """
        print("WARNING: unable to find filter ({}) name".format(filter))
        print("in any of the given {} inputs".format(cattype))
        print("Using the first input for now. Make sure input catalog names have")
        print("the appropriate filter name in the filename to get matching to work.")

    def path_defs(self):
        """Expand input files to have full paths"""
        self.input_xml = os.path.abspath(os.path.expandvars(self.input_xml))
        self.pointing_file = os.path.abspath(os.path.expandvars(self.pointing_file))
        self.output_dir = os.path.abspath(os.path.expandvars(self.output_dir))
        self.simdata_output_dir = os.path.abspath(os.path.expandvars(self.simdata_output_dir))
        if self.table_file is not None:
            self.table_file = os.path.abspath(os.path.expandvars(self.table_file))

        ensure_dir_exists(self.output_dir)
        ensure_dir_exists(self.simdata_output_dir)

        # self.subarray_def_file = self.set_config(self.subarray_def_file, 'subarray_def_file')
        # self.readpatt_def_file = self.set_config(self.readpatt_def_file, 'readpatt_def_file')
        # self.filtpupil_pairs = self.set_config(self.filtpupil_pairs, 'filtpupil_pairs')
        # self.fluxcal = self.set_config(self.fluxcal, 'fluxcal')
        # self.filter_throughput = self.set_config(self.filter_throughput, 'filter_throughput')
        # self.dq_init_config = self.set_config(self.dq_init_config, 'dq_init_config')
        # self.refpix_config = self.set_config(self.refpix_config, 'refpix_config')
        # self.saturation_config = self.set_config(self.saturation_config, 'saturation_config')
        # self.superbias_config = self.set_config(self.superbias_config, 'superbias_config')
        # self.linearity_config = self.set_config(self.linearity_config, 'linearity_config')

        if self.observation_list_file is not None:
            self.observation_list_file = os.path.abspath(os.path.expandvars(self.observation_list_file))
        # if self.crosstalk not in [None, 'config']:
        #     self.crosstalk = os.path.abspath(os.path.expandvars(self.crosstalk))
        # elif self.crosstalk == 'config':
        #     self.crosstalk = os.path.join(self.modpath, 'config', self.configfiles['crosstalk'])

    def reffile_setup(self, offline=False):
        """Create lists of reference files associate with each detector.

        Parameters
        ----------
        instrument : str
            Name of instrument
        """
        # Prepare to find files listed as 'config'
        # and set up PSF path

        # set up as dictionary of dictionaries
        self.configfiles = {}
        self.psfpath = {}
        self.psfbasename = {}
        self.psfpixfrac = {}
        self.reference_file_dir = {}

        for instrument in 'nircam niriss fgs'.split():
            self.configfiles[instrument] = {}
            self.psfpath[instrument] = os.path.join(self.datadir, instrument, 'gridded_psf_library')
            self.psfbasename[instrument] = instrument
            self.reference_file_dir[instrument] = os.path.join(self.datadir, instrument, 'reference_files')

            # Set instrument-specific file paths
            if instrument == 'nircam':
                self.psfpixfrac[instrument] = 0.25
            elif instrument == 'niriss':
                self.psfpixfrac[instrument] = 0.1
            elif instrument == 'fgs':
                self.psfpixfrac[instrument] = 0.1

            # Set global file paths
            self.configfiles[instrument]['dq_init_config'] = os.path.join(self.modpath, 'config', 'dq_init.cfg')
            self.configfiles[instrument]['saturation_config'] = os.path.join(self.modpath, 'config', 'saturation.cfg')
            self.configfiles[instrument]['superbias_config'] = os.path.join(self.modpath, 'config', 'superbias.cfg')
            self.configfiles[instrument]['refpix_config'] = os.path.join(self.modpath, 'config', 'refpix.cfg')
            self.configfiles[instrument]['linearity_config'] = os.path.join(self.modpath, 'config', 'linearity.cfg')
            self.configfiles[instrument]['filter_throughput'] = os.path.join(self.modpath, 'config', 'placeholder.txt')

        for instrument in 'miri nirspec'.split():
            self.configfiles[instrument] = {}
            self.psfpixfrac[instrument] = 0
            self.psfbasename[instrument] = 'N/A'

        # create empty dictionaries
        list_names = 'superbias linearity gain saturation ipc astrometric pam dark lindark'.split()
        for list_name in list_names:
            setattr(self, '{}_list'.format(list_name), {})

        self.det_list = {}
        self.det_list['nircam'] = ['A1', 'A2', 'A3', 'A4', 'A5', 'B1', 'B2', 'B3', 'B4', 'B5']
        self.det_list['niriss'] = ['NIS']
        self.det_list['fgs'] = ['G1', 'G2']
        self.det_list['nirspec'] = ['NRS']
        self.det_list['miri'] = ['MIR']

        for instrument in 'nircam niriss fgs miri nirspec'.split():
            for list_name in list_names:
                getattr(self, '{}_list'.format(list_name))[instrument] = {}

            if offline:
                # no access to central store. Set all files to none.
                for list_name in list_names:
                    if list_name in 'dark lindark'.split():
                        default_value = ['None']
                    else:
                        default_value = 'None'
                    for det in self.det_list[instrument]:
                        getattr(self, '{}_list'.format(list_name))[instrument][det] = default_value

            elif instrument == 'nircam':
                sb_dir = os.path.join(self.datadir, 'nircam/reference_files/superbias')
                lin_dir = os.path.join(self.datadir, 'nircam/reference_files/linearity')
                gain_dir = os.path.join(self.datadir, 'nircam/reference_files/gain')
                sat_dir = os.path.join(self.datadir, 'nircam/reference_files/saturation')
                ipc_dir = os.path.join(self.datadir, 'nircam/reference_files/ipc')
                dist_dir = os.path.join(self.datadir, 'nircam/reference_files/distortion')
                pam_dir = os.path.join(self.datadir, 'nircam/reference_files/pam')
                rawdark_dir = os.path.join(self.datadir, 'nircam/darks/raw')
                lindark_dir = os.path.join(self.datadir, 'nircam/darks/linearized')
                for det in self.det_list[instrument]:
                    sbfiles = glob(os.path.join(sb_dir, '*fits'))
                    self.superbias_list[instrument][det] = [d for d in sbfiles if 'NRC' + det in d][0]
                    linfiles = glob(os.path.join(lin_dir, '*fits'))
                    longdet = deepcopy(det)
                    if '5' in det:
                        longdet = det.replace('5', 'LONG')
                    self.linearity_list[instrument][det] = [d for d in linfiles if 'NRC' + longdet in d][0]

                    gainfiles = glob(os.path.join(gain_dir, '*fits'))
                    self.gain_list[instrument][det] = [d for d in gainfiles if 'NRC' + det in d][0]

                    satfiles = glob(os.path.join(sat_dir, '*fits'))
                    self.saturation_list[instrument][det] = [d for d in satfiles if 'NRC' + det in d][0]

                    ipcfiles = glob(os.path.join(ipc_dir, 'Kernel_to_add_IPC*fits'))
                    self.ipc_list[instrument][det] = [d for d in ipcfiles if 'NRC' + det in d][0]

                    distfiles = glob(os.path.join(dist_dir, '*asdf'))
                    self.astrometric_list[instrument][det] = [d for d in distfiles if 'NRC' + det in d][0]

                    pamfiles = glob(os.path.join(pam_dir, '*fits'))
                    self.pam_list[instrument][det] = [d for d in pamfiles if det in d][0]

                    self.dark_list[instrument][det] = glob(os.path.join(rawdark_dir, det, '*.fits'))
                    self.lindark_list[instrument][det] = glob(os.path.join(lindark_dir, det, '*.fits'))

            elif instrument in ['nirspec', 'miri']:
                for key in 'subarray_def_file fluxcal filtpupil_pairs readpatt_def_file crosstalk ' \
                           'dq_init_config saturation_config superbias_config refpix_config ' \
                           'linearity_config filter_throughput'.split():
                    self.configfiles[instrument][key] = 'N/A'
                default_value = 'none'
                for list_name in list_names:
                    for det in self.det_list[instrument]:
                        getattr(self, '{}_list'.format(list_name))[instrument][det] = default_value

            else:  # niriss and fgs
                for det in self.det_list[instrument]:
                    if det == 'G1':
                        self.ipc_list[instrument][det] = glob(os.path.join(self.reference_file_dir[instrument], 'ipc/Kernel_to_add_IPC_effects_from_jwst_fgs_ipc_0003.fits'))[0]
                        self.dark_list[instrument][det] = glob(os.path.join(self.datadir, 'fgs/darks/raw',
                                               '*30632_1x88_FGSF03511-D-NR-G1-5346180117_1_497_SE_2015-12-12T19h00m12_dms_uncal*.fits'))
                        self.astrometric_list[instrument][det] = glob(
                            os.path.join(self.reference_file_dir[instrument],
                                         'distortion/*distortion_0004.asdf'))[0]
                        self.lindark_list[instrument][det] = glob(os.path.join(self.datadir, 'fgs/darks/linearized', '*_497_*fits'))

                    elif det == 'G2':
                        self.ipc_list[instrument][det] = glob(os.path.join(self.reference_file_dir[instrument], 'ipc/Kernel_to_add_IPC_effects_from_jwst_fgs_ipc_0003.fits'))[0]
                        self.dark_list[instrument][det] = glob(os.path.join(self.datadir, 'fgs/darks/raw',
                                               '*30670_1x88_FGSF03511-D-NR-G2-5346181816_1_498_SE_2015-12-12T21h31m01_dms_uncal*.fits'))
                        self.astrometric_list[instrument][det] = glob(
                            os.path.join(self.reference_file_dir[instrument],
                                         'distortion/*distortion_0003.asdf'))[0]
                        self.lindark_list[instrument][det] = glob(os.path.join(self.datadir, 'fgs/darks/linearized', '*_498_*fits'))

                    elif det == 'NIS':
                        self.ipc_list[instrument][det] = glob(os.path.join(self.reference_file_dir[instrument],
                                                                           'ipc/Kernel_to_add_IPC_effects_from_jwst_niriss_ipc_0007.fits'))[0]
                        self.dark_list[instrument][det] = glob(os.path.join(self.datadir, 'niriss/darks/raw',
                                                                            '*uncal.fits'))
                        self.lindark_list[instrument][det] = glob(os.path.join(self.datadir, 'niriss/darks/linearized',
                                                                               '*linear_dark_prep_object.fits'))
                        self.astrometric_list[instrument][det] = glob(
                            os.path.join(self.reference_file_dir[instrument],
                                         'distortion/*distortion*.asdf'))[0]
                    self.superbias_list[instrument][det] = glob(os.path.join(self.reference_file_dir[instrument], 'superbias/*superbias*.fits'))[0]
                    self.linearity_list[instrument][det] = glob(os.path.join(self.reference_file_dir[instrument], 'linearity/*linearity*.fits'))[0]
                    self.gain_list[instrument][det] = glob(os.path.join(self.reference_file_dir[instrument], 'gain/*gain*.fits'))[0]
                    self.saturation_list[instrument][det] = glob(os.path.join(self.reference_file_dir[instrument], 'saturation/*saturation*.fits'))[0]

                    # suspecting that the FGS wcs reference file has a problem
                    # self.astrometric_list[instrument][det] = 'none'

                    self.pam_list[instrument][det] = glob(os.path.join(self.reference_file_dir[instrument], 'pam/*area*.fits'))[0]

    def set_config(self, file, prop):
        """
        If a given file is listed as 'config'
        then set it in the yaml output as being in
        the config subdirectory.

        Parameters
        ----------
        file : str
            Name of the input file
        prop : str
            Type of file that file is.

        Returns:
        --------
        file : str
            Full path name to the input file
        """
        if file.lower() not in ['config']:
            file = os.path.abspath(file)
        elif file.lower() == 'config':
            file = os.path.join(self.modpath, 'config', self.configfiles[prop])
        return file


    def get_psf_path(self):
        """ Create a list of the path to the PSF library directory for
        each observation/visit/exposure in the APT program.

        Parameters:
        -----------
        psf_paths : list, str, or None
            Either a list of the paths to the PSF library(ies), with a
            length equal to the number of activities in the APT program,
            a string containing the path to one PSF library,
            or None. If a list, each path will be written
            chronologically into each yaml file. If a string, that path
            will be written into every yaml file. If None, the
            default PSF library path will be used for all yamls.

        Returns:
        --------
        paths_out : list
            The list of paths to the PSF library(ies), with a length
            equal to the number of activities in the APT program.
        """
        act_ids = sorted(list(set(self.info['act_id'])))
        act_id_indices = []
        for act_id in self.info['act_id']:
            act_id_indices.append(act_ids.index(act_id))
        n_activities = len(act_ids)

        # If no path explicitly provided, use the default path.
        if self.psf_paths is None:
            print('No PSF path provided. Using default path as PSF path for all yamls.')
            paths_out = []
            for instrument in self.info['Instrument']:
                default_path = self.global_psfpath[instrument.lower()]
                # default_path = os.path.join(self.datadir, instrument.lower(), 'webbpsf_library')
                paths_out.append(default_path)
            return paths_out

        elif isinstance(self.psf_paths, str):
            print('Using provided PSF path.')
            paths_out = [self.psf_paths] * len(self.info['act_id'])
            return paths_out

        elif isinstance(self.psf_paths, list) and len(self.psf_paths) != n_activities:
            raise ValueError('Invalid PSF paths parameter provided. Please '
                             'provide the psf_paths in the form of a list of '
                             'strings with a length equal to the number of '
                             'activities in the APT program ({}), not equal to {}.'
                             .format(n_activities, len(self.psf_paths)))

        elif isinstance(self.psf_paths, list):
            print('Using provided PSF paths.')
            paths_out = [sorted(self.psf_paths)[i] for i in act_id_indices]  # Why is this sorted..?? Seg number?
            return paths_out

        elif not isinstance(self.psf_paths, list) or not isinstance(self.psf_paths, str):
            raise TypeError('Invalid PSF paths parameter provided. Please '
                            'provide the psf_paths in the form of a list or string, not'
                            '{}'.format(type(self.psf_paths)))


    def table_to_dict(self, tab):
        """
        Convert the ascii table of observations to a dictionary

        Parameters
        ----------
        tab : obj
            astropy.table.Table containing observation information

        Returns
        -------
        dict : dict
            Dictionary of observation information
        """
        dict = {}
        for colname in tab.colnames:
            dict[colname] = tab[colname].data
        return dict

    def write_yaml(self, input):
        """
        Create yaml file for a single exposure/detector

        Parameters
        ----------
        input : dict
            dictionary containing all needed exposure
            information for one exposure
        """
        instrument = input['Instrument']
        # select the right filter
        if input['detector'] in ['NIS']:
            # if input['APTTemplate'] == 'NirissExternalCalibration': 'NirissImaging':
            filtkey = 'FilterWheel'
            pupilkey = 'PupilWheel'
            # set the FilterWheel and PupilWheel for NIRISS
            if input['APTTemplate'] in ['NirissAmi']:
                filter_name = input['Filter']
                input[filtkey] = filter_name
                input[pupilkey] = 'NRM'
            elif input['APTTemplate'] not in ['NirissExternalCalibration', 'NirissWfss']:
                filter_name = input['Filter']
                if filter_name in NIRISS_PUPIL_WHEEL_ELEMENTS:
                    input[pupilkey] = filter_name
                    input[filtkey] = 'CLEAR'
                elif filter_name in NIRISS_FILTER_WHEEL_ELEMENTS:
                    input[pupilkey] = 'CLEARP'
                    input[filtkey] = filter_name
                else:
                    raise RuntimeError('Filter {} not valid'.format(filter_name))
            catkey = ''
        elif input['detector'] in ['FGS']:
            filtkey = 'FilterWheel'
            pupilkey = 'PupilWheel'
            catkey = ''
        elif input['detector'] in ['NRS', 'MIR']:
            filtkey = 'FilterWheel'
            pupilkey = 'PupilWheel'
            catkey = ''
        else:
            if np.int(input['detector'][-1]) < 5:
                filtkey = 'ShortFilter'
                pupilkey = 'ShortPupil'
                catkey = 'sw'
            else:
                filtkey = 'LongFilter'
                pupilkey = 'LongPupil'
                catkey = 'lw'

        outfile = input['outputfits']
        yamlout = input['yamlfile']

        yamlout = os.path.join(self.output_dir, yamlout)
        with open(yamlout, 'w') as f:
            f.write('Inst:\n')
            f.write('  instrument: {}          # Instrument name\n'.format(instrument))
            f.write('  mode: {}                # Observation mode (e.g. imaging, WFSS)\n'.format(input['Mode']))
            f.write('  use_JWST_pipeline: {}   # Use pipeline in data transformations\n'.format(input['use_JWST_pipeline']))
            f.write('\n')
            f.write('Readout:\n')
            f.write('  readpatt: {}        # Readout pattern (RAPID, BRIGHT2, etc) overrides nframe, nskip unless it is not recognized\n'.format(input['ReadoutPattern']))
            f.write('  ngroup: {}              # Number of groups in integration\n'.format(input['Groups']))
            f.write('  nint: {}          # Number of integrations per exposure\n'.format(input['Integrations']))
            f.write('  resets_bet_ints: {} #Number of detector resets between integrations\n'.format(self.resets_bet_ints))

            if instrument.lower() == 'nircam':
                # if input['aperture'] in ['NRCA3_DHSPIL', 'NRCB4_DHSPIL']:
                if 'NRCA3_DHSPIL' in input['aperture'] or 'NRCB4_DHSPIL' in input['aperture']: # in ['NRCA3_DHSPIL', 'NRCB4_DHSPIL']:
                    full_ap = input['aperture']
                else:
                    apunder = input['aperture'].find('_')
                    full_ap = 'NRC' + input['detector'] + '_' + input['aperture'][apunder + 1:]
            if instrument.lower() in ['niriss', 'fgs']:
                full_ap = input['aperture']

            subarray_definitions = self.global_subarray_definitions[instrument.lower()]


            if full_ap not in subarray_definitions['AperName']:
                full_ap_new = [apername for apername, name in
                               np.array(subarray_definitions['AperName', 'Name']) if
                               (full_ap in apername) or (full_ap in name)]
                if len(full_ap_new) > 1 or len(full_ap_new) == 0:
                    raise ValueError('Cannot match {} with valid aperture name for observation {}.'
                                     .format(full_ap, input['obs_num']))
                else:
                    full_ap = full_ap_new[0]

            f.write('  array_name: {}    # Name of array (FULL, SUB160, SUB64P, etc) overrides subarray_bounds below\n'.format(full_ap))
            f.write('  filter: {}       # Filter of simulated data (F090W, F322W2, etc)\n'.format(input[filtkey]))
            f.write('  pupil: {}        # Pupil element for simulated data (CLEAR, GRISMC, etc)\n'.format(input[pupilkey]))
            f.write('\n')
            f.write('Reffiles:                                 # Set to None or leave blank if you wish to skip that step\n')
            f.write('  dark: {}   # Dark current integration used as the base\n'.format(input['dark']))
            f.write('  linearized_darkfile: {}   # Linearized dark ramp to use as input. Supercedes dark above\n'.format(input['lindark']))
            f.write('  badpixmask: None   # If linearized dark is used, populate output DQ extensions using this file\n')
            f.write('  superbias: {}     # Superbias file. Set to None or leave blank if not using\n'.format(input['superbias']))
            f.write('  linearity: {}    # linearity correction coefficients\n'.format(input['linearity']))
            f.write('  saturation: {}    # well depth reference files\n'.format(input['saturation']))
            f.write('  gain: {} # Gain map\n'.format(input['gain']))
            f.write('  pixelflat: None \n')
            f.write('  illumflat: None                               # Illumination flat field file\n')
            f.write('  astrometric: {}  # Astrometric distortion file (asdf)\n'.format(input['astrometric']))
            f.write('  ipc: {} # File containing IPC kernel to apply\n'.format(input['ipc']))
            f.write(('  invertIPC: {}      # Invert the IPC kernel before the convolution. True or False. Use True if the kernel is '
                     'designed for the removal of IPC effects, like the JWST reference files are.\n'.format(False)))
            f.write('  occult: None                                    # Occulting spots correction image\n')
            f.write(('  pixelAreaMap: {}      # Pixel area map for the detector. Used to introduce distortion into the output ramp.\n'
                     .format(input['pixelAreaMap'])))
            f.write(('  subarray_defs: {} # File that contains a list of all possible subarray names and coordinates\n'
                     .format(input['subarray_def_file'])))
            f.write(('  readpattdefs: {}  # File that contains a list of all possible readout pattern names and associated '
                     'NFRAME/NSKIP values\n'.format(input['readpatt_def_file'])))
            f.write('  crosstalk: {}   # File containing crosstalk coefficients\n'.format(input['crosstalk_file']))
            f.write(('  filtpupilcombo: {}   # File that lists the filter wheel element / pupil wheel element combinations. '
                     'Used only in writing output file\n'.format(input['filtpupilcombo_file'])))
            f.write(('  flux_cal: {} # File that lists flux conversion factor and pivot wavelength for each filter. Only '
                     'used when making direct image outputs to be fed into the grism disperser code.\n'.format(input['flux_cal_file'] )))
            f.write('  filter_throughput: {} #File containing filter throughput curve\n'.format(self.configfiles[instrument.lower()]['filter_throughput']))
            f.write('\n')
            f.write('nonlin:\n')
            f.write('  limit: 60000.0                           # Upper singal limit to which nonlinearity is applied (ADU)\n')
            f.write('  accuracy: 0.000001                        # Non-linearity accuracy threshold\n')
            f.write('  maxiter: 10                              # Maximum number of iterations to use when applying non-linearity\n')
            f.write('  robberto:  False                         # Use Massimo Robberto type non-linearity coefficients\n')
            f.write('\n')
            f.write('cosmicRay:\n')
            cosmic_ray_path = os.path.join(self.datadir, instrument.lower(), 'cosmic_ray_library')
            f.write('  path: {}               # Path to CR library\n'.format(cosmic_ray_path))
            f.write('  library: SUNMAX    # Type of cosmic rayenvironment (SUNMAX, SUNMIN, FLARE)\n')
            f.write('  scale: 1.0     # Cosmic ray scaling factor\n')
            # temporary tweak here to make it work with NIRISS
            detector_label = input['detector']

            if instrument.lower() in ['nircam', 'wfsc']:
                # detector_label = input['detector']
                f.write('  suffix: IPC_NIRCam_{}    # Suffix of library file names\n'.format(
                    detector_label))
            elif instrument.lower() == 'niriss':
                f.write('  suffix: IPC_NIRISS_{}    # Suffix of library file names\n'.format(
                    detector_label))
            elif instrument.lower() == 'fgs':
                if detector_label == 'G1':
                    detector_string = 'GUIDER1'
                elif detector_label == 'G2':
                    detector_string = 'GUIDER2'
                f.write('  suffix: IPC_FGS_{}    # Suffix of library file names\n'.format(
                    detector_string))
            f.write('  seed: {}                 # Seed for random number generator\n'.format(np.random.randint(1, 2**32-2)))
            f.write('\n')
            f.write('simSignals:\n')
            if instrument.lower() in ['nircam', 'wfsc']:
                PointSourceCatalog = input['{}_ptsrc'.format(catkey)]
                GalaxyCatalog = input['{}_galcat'.format(catkey)]
                ExtendedCatalog = input['{}_ext'.format(catkey)]
                ExtendedScale = input['{}_extscl'.format(catkey)]
                ExtendedCenter = input['{}_extcent'.format(catkey)]
                MovingTargetList = input['{}_movptsrc'.format(catkey)]
                MovingTargetSersic = input['{}_movgal'.format(catkey)]
                MovingTargetExtended = input['{}_movext'.format(catkey)]
                MovingTargetConvolveExtended = input['{}_movconv'.format(catkey)]
                MovingTargetToTrack = input['{}_solarsys'.format(catkey)]
                BackgroundRate = input['{}_bkgd'.format(catkey)]
            elif instrument.lower() in ['niriss', 'fgs']:
                PointSourceCatalog = input['PointSourceCatalog']
                GalaxyCatalog = input['GalaxyCatalog']
                ExtendedCatalog = input['ExtendedCatalog']
                ExtendedScale = input['ExtendedScale']
                ExtendedCenter = input['ExtendedCenter']
                MovingTargetList = input['MovingTargetList']
                MovingTargetSersic = input['MovingTargetSersic']
                MovingTargetExtended = input['MovingTargetExtended']
                MovingTargetConvolveExtended = input['MovingTargetConvolveExtended']
                MovingTargetToTrack = input['MovingTargetToTrack']
                BackgroundRate = input['BackgroundRate']

            f.write(('  pointsource: {}   #File containing a list of point sources to add (x, y locations and magnitudes)\n'
                     .format(PointSourceCatalog)))
            f.write('  gridded_psf_library_row_padding: 4  # Number of outer rows and columns to avoid when evaluating library. RECOMMEND 4.\n')
            f.write('  psf_wing_threshold_file: {}   # File defining PSF sizes versus magnitude\n'.format(input['psf_wing_threshold_file']))
            f.write('  add_psf_wings: {}  # Whether or not to place the core of the psf from the gridded library into an image of the wings before adding.\n'.format(self.add_psf_wings))
            f.write('  psfpath: {}   #Path to PSF library\n'.format(input['psfpath']))
            f.write('  psfwfe: {}   #PSF WFE value (predicted or requirements)\n'.format(self.psfwfe))
            f.write('  psfwfegroup: {}      #WFE realization group (0 to 4)\n'.format(self.psfwfegroup))
            f.write(('  galaxyListFile: {}    #File containing a list of positions/ellipticities/magnitudes of galaxies '
                     'to simulate\n'.format(GalaxyCatalog)))
            f.write('  extended: {}          #Extended emission count rate image file name\n'.format(ExtendedCatalog))
            f.write('  extendedscale: {}                          #Scaling factor for extended emission image\n'.format(ExtendedScale))
            f.write(('  extendedCenter: {}                   #x, y pixel location at which to place the extended image '
                     'if it is smaller than the output array size\n'.format(ExtendedCenter)))
            f.write(('  PSFConvolveExtended: True #Convolve the extended image with the PSF before adding to the output '
                     'image (True or False)\n'))
            f.write(('  movingTargetList: {}          #Name of file containing a list of point source moving targets (e.g. '
                     'KBOs, asteroids) to add.\n'.format(MovingTargetList)))
            f.write(('  movingTargetSersic: {}  #ascii file containing a list of 2D sersic profiles to have moving through '
                     'the field\n'.format(MovingTargetSersic)))
            f.write(('  movingTargetExtended: {}      #ascii file containing a list of stamp images to add as moving targets '
                     '(planets, moons, etc)\n'.format(MovingTargetExtended)))
            f.write(('  movingTargetConvolveExtended: {}       #convolve the extended moving targets with PSF before adding.\n'
                     .format(MovingTargetConvolveExtended)))
            f.write(('  movingTargetToTrack: {} #File containing a single moving target which JWST will track during '
                     'observation (e.g. a planet, moon, KBO, asteroid)	This file will only be used if mode is set to '
                     '"moving_target" \n'.format(MovingTargetToTrack)))
            f.write('  zodiacal:  None                          #Zodiacal light count rate image file \n')
            f.write('  zodiscale:  1.0                            #Zodi scaling factor\n')
            f.write('  scattered:  None                          #Scattered light count rate image file\n')
            f.write('  scatteredscale: 1.0                        #Scattered light scaling factor\n')
            f.write(('  bkgdrate: {}                         #Constant background count rate (ADU/sec/pixel) or '
                     '"high","medium","low" similar to what is used in the ETC\n'.format(BackgroundRate)))
            f.write(('  poissonseed: {}                  #Random number generator seed for Poisson simulation)\n'
                     .format(np.random.randint(1, 2**32-2))))
            f.write('  photonyield: True                         #Apply photon yield in simulation\n')
            f.write('  pymethod: True                            #Use double Poisson simulation for photon yield\n')
            f.write('  expand_catalog_for_segments: {}                     # Expand catalog for 18 segments and use distinct PSFs\n'
                .format(self.expand_catalog_for_segments))

            f.write('\n')
            f.write('Telescope:\n')
            f.write('  ra: {}                      # RA of simulated pointing\n'.format(input['ra_ref']))
            f.write('  dec: {}                    # Dec of simulated pointing\n'.format(input['dec_ref']))
            if 'pav3' in input.keys():
                pav3_value = input['pav3']
            else:
                pav3_value = input['PAV3']
            f.write('  rotation: {}                    # PA_V3 in degrees, i.e. the position angle of the V3 axis at V1 (V2=0, V3=0) measured from N to E.\n'.format(pav3_value))
            f.write('  tracking: {}   #Telescope tracking. Can be sidereal or non-sidereal\n'.format(self.tracking))
            f.write('\n')
            f.write('newRamp:\n')
            f.write('  dq_configfile: {}\n'.format(self.configfiles[instrument.lower()]['dq_init_config']))
            f.write('  sat_configfile: {}\n'.format(self.configfiles[instrument.lower()]['saturation_config']))
            f.write('  superbias_configfile: {}\n'.format(self.configfiles[instrument.lower()]['superbias_config']))
            f.write('  refpix_configfile: {}\n'.format(self.configfiles[instrument.lower()]['refpix_config']))
            f.write('  linear_configfile: {}\n'.format(self.configfiles[instrument.lower()]['linearity_config']))
            f.write('\n')
            f.write('Output:\n')
            # f.write('  use_stsci_output_name: {} # Output filename should follow STScI naming conventions (True/False)\n'.format(outtf))
            f.write('  directory: {}  # Output directory\n'.format(self.simdata_output_dir))
            f.write('  file: {}   # Output filename\n'.format(outfile))
            f.write(("  datatype: {} # Type of data to save. 'linear' for linearized ramp. 'raw' for raw ramp. 'linear, "
                     "raw' for both\n".format(self.datatype)))
            f.write('  format: DMS          # Output file format Options: DMS, SSR(not yet implemented)\n')
            f.write('  save_intermediates: False   # Save intermediate products separately (point source image, etc)\n')
            f.write('  grism_source_image: {}   # grism\n'.format(input['grism_source_image']))
            f.write('  unsigned: True   # Output unsigned integers? (0-65535 if true. -32768 to 32768 if false)\n')
            f.write('  dmsOrient: True    # Output in DMS orientation (vs. fitswriter orientation).\n')
            f.write('  program_number: {}    # Program Number\n'.format(input['ProposalID']))
            f.write('  title: {}   # Program title\n'.format(input['Title'].replace(':', ', ')))
            f.write('  PI_Name: {}  # Proposal PI Name\n'.format(input['PI_Name']))
            f.write('  Proposal_category: {}  # Proposal category\n'.format(input['Proposal_category']))
            f.write('  Science_category: {}  # Science category\n'.format(input['Science_category']))
            f.write("  observation_number: '{}'    # Observation Number\n".format(input['obs_num']))
            f.write('  observation_label: {}    # User-generated observation Label\n'.format(input['obs_label'].strip()))
            f.write("  visit_number: '{}'    # Visit Number\n".format(input['visit_num']))
            f.write("  visit_group: '{}'    # Visit Group\n".format(input['visit_group']))
            f.write("  visit_id: '{}'    # Visit ID\n".format(input['visit_id']))
            f.write("  sequence_id: '{}'    # Sequence ID\n".format(input['sequence_id']))
            f.write("  activity_id: '{}'    # Activity ID. Increment with each exposure.\n".format(input['act_id']))
            f.write("  exposure_number: '{}'    # Exposure Number\n".format(input['exposure']))
            f.write("  obs_id: '{}'   # Observation ID number\n".format(input['observation_id']))
            f.write("  date_obs: '{}'  # Date of observation\n".format(input['date_obs']))
            f.write("  time_obs: '{}'  # Time of observation\n".format(input['time_obs']))
            # f.write("  obs_template: '{}'  # Observation template\n".format(input['obs_template']))
            f.write("  primary_dither_type: {}  # Primary dither pattern name\n".format(input['PrimaryDitherType']))
            f.write("  total_primary_dither_positions: {}  # Total number of primary dither positions\n".format(input['PrimaryDithers']))
            f.write("  primary_dither_position: {}  # Primary dither position number\n".format(np.int(input['primary_dither_num'])))
            f.write("  subpix_dither_type: {}  # Subpixel dither pattern name\n".format(input['SubpixelDitherType']))
            # For WFSS we need to strip out the '-Points' from
            # the number of subpixel positions entry
            dash = input['SubpixelPositions'].find('-')
            if (dash == -1):
                val = input['SubpixelPositions']
            else:
                val = input['SubpixelPositions'][0:dash]
            if val == 'None':
                val = 1
            # try:
            #     dash = input['SubpixelPositions'].find('-')
            #     val = input['SubpixelPositions'][0:dash]
            # except:
            #     val = input['SubpixelPositions']
            f.write("  total_subpix_dither_positions: {}  # Total number of subpixel dither positions\n".format(val))
            f.write("  subpix_dither_position: {}  # Subpixel dither position number\n".format(np.int(input['subpix_dither_num'])))
            f.write("  xoffset: {}  # Dither pointing offset in x (arcsec)\n".format(input['idlx']))
            f.write("  yoffset: {}  # Dither pointing offset in y (arcsec)\n".format(input['idly']))

        return yamlout

    def add_options(self, parser=None, usage=None):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage, description='Simulate JWST ramp')
        parser.add_argument("--input_xml", help='XML file from APT describing the observations.')
        parser.add_argument("--pointing_file", help='Pointing file from APT describing observations.')
        parser.add_argument("--datatype", help='Type of data to save. Can be "linear", "raw" or "linear, raw"', default="linear")
        parser.add_argument("--output_dir", help='Directory into which the yaml files are output', default='./')
        parser.add_argument("--table_file", help='Ascii table containing observation info. Use this or xml + pointing files.', default=None)
        parser.add_argument("--use_nonstsci_names", help="Use STScI naming convention for output files", action='store_true')
        parser.add_argument("--subarray_def_file", help="Ascii file containing subarray definitions", default='config')
        parser.add_argument("--readpatt_def_file", help='Ascii file containing readout pattern definitions', default='config')
        parser.add_argument("--crosstalk", help="Crosstalk coefficient file", default='config')
        parser.add_argument("--filtpupil_pairs", help="List of paired filter/pupil elements", default='config')
        parser.add_argument("--fluxcal", help="File with zeropoints per filter", default='config')
        parser.add_argument("--dq_init_config", help="DQ Initialization config file", default='config')
        parser.add_argument("--saturation_config", help="Saturation config file", default='config')
        parser.add_argument("--superbias_config", help="Superbias subtraction config file", default='config')
        parser.add_argument("--refpix_config", help="Refpix subtraction config file", default='config')
        parser.add_argument("--linearity_config", help="Linearity config file", default='config')
        parser.add_argument("--observation_list_file", help="Table file containing epoch start times, telescope roll angles, catalogs for each observation", default=None)
        parser.add_argument("--use_JWST_pipeline", help='True/False', action='store_true')
        parser.add_argument("--use_linearized_darks", help='True/False', action='store_true')
        parser.add_argument("--simdata_output_dir", help='Output directory for simulated exposure files', default='./')
        parser.add_argument("--psfpath", help='Directory containing PSF library',
                            default=os.path.join(self.datadir, 'gridded_psf_library'))
        parser.add_argument("--psfbasename", help="Basename of the files in the PSF library", default='nircam')
        parser.add_argument("--psfpixfrac", help="Subpixel centering resolution of files in PSF library", default=0.25)
        parser.add_argument("--psfwfe", help="Wavefront error value to use for PSFs", default='predicted')
        parser.add_argument("--psfwfegroup", help="Realization index number for PSF files", default=0)
        parser.add_argument("--resets_bet_ints", help="Number of detector resets between integrations", default=1)
        parser.add_argument("--tracking", help="Type of telescope tracking: 'sidereal' or 'non-sidereal'", default='sidereal')

        return parser


if __name__ == '__main__':

    usagestring = 'USAGE: yaml_generator.py NIRCam_obs.xml NIRCam_obs.pointing'

    input = SimInput()
    parser = input.add_options(usage=usagestring)
    args = parser.parse_args(namespace=input)
    input.reffile_setup()
    input.create_inputs()
