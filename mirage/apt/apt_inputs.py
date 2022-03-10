# ! /usr/bin/env python

'''
Given APT output files, read in data relevant to the data simulator,
organize, and create input files for the simulator.

Inputs:

xml file - Name of xml file exported from APT.
pointing file - Name of associated pointing file exported from APT.

Optional inputs:

epoch_list - Name of ascii file which lists observation labels and
             associated starting observation time, as well as telescope
             roll angle (PAV3). If you wish to have observations back-
             to-back, give them all the same starting time.

Outputs:

output_csv - Ascii table containing all relevant information needed
             to construct a ramp simulator input file for each
             observation/dither/detector combination.

Dependencies:

argparse, lxml, astropy, numpy, collections

JWST Calibration pipeline (only for the use of set_telescope_pointing.py)
in order to translate PAV3 values to local roll angles for each detector.


HISTORY:

July 2017 - V0: Initial version. Bryan Hilbert
Feb 2018  - V1: Updated to work for multiple filter pairs per observation
            Lauren Chambers
August 2018 - V2: Replaced manual Ra, Dec calculations with pysiaf functionality
October 2018 - Major modifications to read programs of all science instruments and parallels
               Johannes Sahlmann
'''
import copy
import os
import logging
import re
import argparse
import pkg_resources
import warnings

from astropy.table import Table, vstack
from astropy.time import Time, TimeDelta
from astropy.io import ascii
import numpy as np
from pysiaf import JWST_PRD_VERSION, rotations, Siaf
import yaml

from . import read_apt_xml
from ..logging import logging_functions
from ..utils import siaf_interface, constants, utils
from mirage.utils.constants import NIRCAM_UNSUPPORTED_PUPIL_VALUES, LOG_CONFIG_FILENAME, STANDARD_LOGFILE_NAME, \
                                   NIRCAM_LW_GRISMTS_APERTURES, NIRCAM_SW_GRISMTS_APERTURES


classpath = os.path.abspath(os.path.join(os.path.dirname(__file__), '../'))
log_config_file = os.path.join(classpath, 'logging', LOG_CONFIG_FILENAME)
logging_functions.create_logger(log_config_file, STANDARD_LOGFILE_NAME)


class AptInput:
    """Summary

    Attributes:
        exposure_tab (TYPE): Description
        input_xml (str): Description
        observation_list_file (str): Description
        obstab (TYPE): Description
        output_csv (TYPE): Description
        pointing_file (str): Description
    """

    def __init__(self, input_xml=None, pointing_file=None, output_dir=None, output_csv=None,
                 observation_list_file=None, offline=False):
        self.logger = logging.getLogger('mirage.apt.apt_inputs')

        self.input_xml = input_xml
        self.pointing_file = pointing_file
        self.output_dir = output_dir
        self.output_csv = output_csv
        self.observation_list_file = observation_list_file
        self.offline = offline

        # Locate the module files, so that we know where to look
        # for config subdirectory
        self.config_path = os.path.join(pkg_resources.resource_filename('mirage', ''), 'config')

    def add_epochs(self, intab):
        """NOT CURRENTLY USED"""
        # add information on the epoch of each observation
        # if the user entered a list of epochs, read that in
        default_date = '2020-10-14'

        if self.epoch_list is not None:
            epochs = ascii.read(self.epoch_list, header_start=0, data_start=1)
        else:
            epochs = Table()
            epochs['observation'] = intab['obs_label']
            epochs['date'] = ['2018-10-14'] * len(intab['obs_label'])
            epochs['pav3'] = [0.] * len(intab['obs_label'])

        # insert epoch info for each observation into dictionary
        epoch_start = []
        epoch_pav3 = []
        for obs in intab['obs_label']:
            match = obs == epochs['observation'].data
            if np.sum(match) == 0:
                self.logger.error("No valid epoch line found for observation {}".format(obs))
                self.logger.error('{}'.format(epochs['observation'].data))
                epoch_start.append(default_date)
                epoch_pav3.append(0.)
            else:
                epoch_start.append(epochs['date'][match].data[0])
                epoch_pav3.append(epochs['pav3'][match].data[0])
        intab['epoch_start_date'] = epoch_start
        intab['pav3'] = epoch_pav3
        return intab

    def add_observation_info(self, intab):
        """Add information about each observation.

        Catalog names, dates, PAV3 values, etc., which are retrieved from the observation list
        yaml file.

        Parameters
        ----------
        intab : obj
            astropy.table.Table containing exposure information

        Returns
        -------
        intab : obj
            Updated table with information from the observation list
            yaml file added.

        """
        with open(self.observation_list_file, 'r') as infile:
            self.obstab = yaml.safe_load(infile)

        OBSERVATION_LIST_FIELDS = 'Date PAV3 Filter PointSourceCatalog GalaxyCatalog ' \
                                  'ExtendedCatalog ExtendedScale ExtendedCenter MovingTargetList ' \
                                  'MovingTargetSersic MovingTargetExtended ' \
                                  'MovingTargetConvolveExtended MovingTargetToTrack ' \
                                  'ImagingTSOCatalog GrismTSOCatalog ' \
                                  'BackgroundRate DitherIndex CosmicRayLibrary CosmicRayScale'.split()

        nircam_mapping = {'ptsrc': 'PointSourceCatalog',
                          'galcat': 'GalaxyCatalog',
                          'ext': 'ExtendedCatalog',
                          'extscl': 'ExtendedScale',
                          'extcent': 'ExtendedCenter',
                          'movptsrc': 'MovingTargetList',
                          'movgal': 'MovingTargetSersic',
                          'movext': 'MovingTargetExtended',
                          'movconv': 'MovingTargetConvolveExtended',
                          'solarsys': 'MovingTargetToTrack',
                          'img_tso': 'ImagingTSOCatalog',
                          'grism_tso': 'GrismTSOCatalog',
                          'bkgd': 'BackgroundRate',
                          }

        unique_instrument_names = [name.lower() for name in np.unique(intab['Instrument'])]

        # initialize dictionary keys
        for key in OBSERVATION_LIST_FIELDS:
            intab[key] = []

        if 'nircam' in unique_instrument_names:
            for channel in ['SW', 'LW']:
                for name, item in nircam_mapping.items():
                    key = '{}_{}'.format(channel.lower(), name)
                    intab[key] = []

        # loop over entries in input dictionary
        for index, instrument in enumerate(intab['Instrument']):
            instrument = instrument.lower()

            # retrieve corresponding entry from observation list
            entry = get_entry(self.obstab, intab['entry_number'][index])

            if instrument == 'nircam':
                # keep the number of entries in the dictionary consistent
                for key in OBSERVATION_LIST_FIELDS:
                    if key in ['Date', 'PAV3', 'Instrument', 'CosmicRayLibrary', 'CosmicRayScale']:
                        value = str(entry[key])
                    else:
                        value = str(None)

                    intab[key].append(value)

                for channel in ['SW', 'LW']:
                    for name, item in nircam_mapping.items():
                        key = '{}_{}'.format(channel.lower(), name)
                        if item in 'ExtendedScale ExtendedCenter MovingTargetConvolveExtended BackgroundRate'.split():
                            intab[key].append(entry['FilterConfig'][channel][item])
                        else:
                            intab[key].append(self.full_path(entry['FilterConfig'][channel][item]))

            else:
                for key in OBSERVATION_LIST_FIELDS:
                    value = str(entry[key])

                    # Expand catalog names to contain full paths
                    catalog_names = 'PointSourceCatalog GalaxyCatalog ' \
                                    'ExtendedCatalog MovingTargetList ' \
                                    'MovingTargetSersic MovingTargetExtended ' \
                                    'MovingTargetToTrack ImagingTSOCatalog ' \
                                    'GrismTSOCatalog'.split()
                    if key in catalog_names:
                        value = self.full_path(value)
                    intab[key].append(value)

                # keep the number of entries in the dictionary consistent
                if 'nircam' in unique_instrument_names:
                    for channel in ['SW', 'LW']:
                        for name, item in nircam_mapping.items():
                            key = '{}_{}'.format(channel.lower(), name)
                            intab[key].append(str(None))

        intab['epoch_start_date'] = intab['Date']
        return intab

    def base36encode(self, integer):
        """
        Translate a base 10 integer to base 36

        Parameters
        ----------
        integer : int
            a base 10 integer

        Returns
        -------
        integer : int
            The integer translated to base 36
        """
        chars, encoded = '0123456789abcdefghijklmnopqrstuvwxyz', ''

        while integer > 0:
            integer, remainder = divmod(integer, 36)
            encoded = chars[remainder] + encoded

        return encoded.zfill(2)

    def combine_dicts(self, dict1, dict2):
        """Combine two dictionaries into a single dictionary.

        Parameters
        ----------
        dict1 : dict
            dictionary
        dict2 : dict
            dictionary

        Returns
        -------
        combined : dict
            Combined dictionary
        """
        combined = dict1.copy()
        combined.update(dict2)
        return combined

    def create_input_table(self, skip_observations=None, verbose=False):
        """
        Main function for creating a table of parameters for each
        exposure

        Parameters
        ----------
        skip_observations : list
            List of observation numbers to be skipped when reading in pointing file

        verbose : bool
            If True, extra information is printed to the log
        """
        self.skip_observations = skip_observations

        # Expand paths to full paths
        # self.input_xml = os.path.abspath(self.input_xml)
        # self.pointing_file = os.path.abspath(self.pointing_file)
        if self.output_csv is not None:
            self.output_csv = os.path.abspath(self.output_csv)
        if self.observation_list_file is not None:
            self.observation_list_file = os.path.abspath(self.observation_list_file)

        # if APT.xml content has already been generated during observation list creation
        # (generate_observationlist.py) load it here
        if self.apt_xml_dict is None:
            raise RuntimeError('self.apt_xml_dict is not defined')

        # Read in the pointing file and produce dictionary
        pointing_dictionary = self.get_pointing_info(self.pointing_file, propid=self.apt_xml_dict['ProposalID'][0],
                                                     skipped_obs_from_xml=self.skip_observations)

        # Check that the .xml and .pointing files agree
        assert len(self.apt_xml_dict['ProposalID']) == len(pointing_dictionary['obs_num']),\
            ('Inconsistent table size from XML file ({}) and pointing file ({}). Something was not '
             'processed correctly in apt_inputs.'.format(len(self.apt_xml_dict['ProposalID']),
                                                         len(pointing_dictionary['obs_num'])))

        # Combine the dictionaries
        observation_dictionary = self.combine_dicts(self.apt_xml_dict, pointing_dictionary)

        # Add epoch and catalog information
        observation_dictionary = self.add_observation_info(observation_dictionary)

        if verbose:
            self.logger.info('Summary of observation dictionary:')
            for key in observation_dictionary.keys():
                self.logger.info('{:<25}: number of elements is {:>5}'.format(key, len(observation_dictionary[key])))

        # Global Alignment observations need to have the pointing information for the
        # FGS exposures updated
        if 'WfscGlobalAlignment' in observation_dictionary['APTTemplate']:
            observation_dictionary = self.global_alignment_pointing(observation_dictionary)

        # Expand the dictionary to have one entry for each detector in each exposure
        self.exposure_tab = self.expand_for_detectors(observation_dictionary)

        # For fiducial point overrides, save the pointing aperture and actual aperture separately
        self.check_aperture_override()

        # Add start times for each exposure
        # Ignore warnings as astropy.time.Time will give a warning
        # related to unknown leap seconds if the date is too far in
        # the future.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.exposure_tab = make_start_times(self.exposure_tab, offline=self.offline)

        # Fix data for filename generation
        # Set parallel seq id
        for j, isparallel in enumerate(self.exposure_tab['ParallelInstrument']):
            if isparallel:
                self.exposure_tab['sequence_id'][j] = '2'

        # set exposure number (new sequence for every combination of seq id and act id and observation number and detector)
        temp_table = Table([self.exposure_tab['sequence_id'], self.exposure_tab['exposure'], self.exposure_tab['act_id'], self.exposure_tab['obs_num'], self.exposure_tab['detector']], names=('sequence_id', 'exposure', 'act_id', 'obs_num', 'detector'))

        for obs_num in np.unique(self.exposure_tab['obs_num']):
            for act_id in np.unique(temp_table['act_id']):
                # prime_index = np.where((temp_table['sequence_id'] == '1') & (temp_table['act_id'] == act_id) & (temp_table['obs_num'] == obs_num))[0]
                # parallel_index = np.where((temp_table['sequence_id'] == '2') & (temp_table['act_id'] == act_id) & (temp_table['obs_num'] == obs_num))[0]
                for detector in np.unique(temp_table['detector']):
                    prime_index = np.where((temp_table['sequence_id'] == '1') & (temp_table['act_id'] == act_id) & (temp_table['obs_num'] == obs_num) & (temp_table['detector']==detector))[0]
                    parallel_index = np.where((temp_table['sequence_id'] == '2') & (temp_table['act_id'] == act_id) & (temp_table['obs_num'] == obs_num) & (temp_table['detector']==detector))[0]

                    temp_table['exposure'][prime_index] = ['{:05d}'.format(n+1) for n in np.arange(len(prime_index))]
                    temp_table['exposure'][parallel_index] = ['{:05d}'.format(n+1) for n in np.arange(len(parallel_index))]
        self.exposure_tab['exposure'] = list(temp_table['exposure'])

        if verbose:
            for key in self.exposure_tab.keys():
                self.logger.info('{:>20} has {:>10} items'.format(key, len(self.exposure_tab[key])))

        # Create a pysiaf.Siaf instance for each instrument in the proposal
        self.siaf = {}
        for instrument_name in np.unique(observation_dictionary['Instrument']):
            self.siaf[instrument_name] = siaf_interface.get_instance(instrument_name)

        # Output to a csv file.
        if self.output_csv is None:
            indir, infile = os.path.split(self.input_xml)
            self.output_csv = os.path.join(self.output_dir, 'Observation_table_for_' + infile.split('.')[0] + '.csv')
        ascii.write(Table(self.exposure_tab), self.output_csv, format='csv', overwrite=True)
        self.logger.info('csv exposure list written to {}'.format(self.output_csv))

    def check_aperture_override(self):
        if bool(self.exposure_tab['FiducialPointOverride']) is True:
            instruments = self.exposure_tab['Instrument']
            apertures = self.exposure_tab['aperture']

            aperture_key = constants.instrument_abbreviations

            fixed_apertures = []
            for i, (instrument, aperture) in enumerate(zip(instruments, apertures)):
                inst_match_ap = aperture.startswith(aperture_key[instrument.lower()])
                if not inst_match_ap:
                    # Handle the one case we understand, for now
                    if instrument.lower() == 'fgs' and aperture[:3] == 'NRC':
                        obs_num = self.exposure_tab['obs_num'][i]

                        if self.exposure_tab['APTTemplate'][i] == 'WfscGlobalAlignment':
                            guider_number = self.exposure_tab['aperture'][i][3]
                        elif self.exposure_tab['APTTemplate'][i] == 'FgsExternalCalibration':
                            guider_number = read_apt_xml.get_guider_number(self.input_xml, obs_num)
                        else:
                            raise ValueError("WARNING: unsupported APT template with Fiducial Override.")
                        guider_aperture = 'FGS{}_FULL'.format(guider_number)
                        fixed_apertures.append(guider_aperture)
                    else:
                        self.logger.error('{} {} {} {}'.format(instrument, aperture, inst_match_ap, aperture_key[instrument.lower()]))
                        raise ValueError('Unknown FiducialPointOverride in program. Instrument = {} but aperture = {}.'.format(instrument, aperture))
                else:
                    fixed_apertures.append(aperture)

            # Add new dictionary entry to document the FiducialPointOverride (pointing aperture)
            self.exposure_tab['pointing_aperture'] = self.exposure_tab['aperture']
            # Rewrite the existing imaging aperture entry to match the primary instrument
            self.exposure_tab['aperture'] = fixed_apertures
        else:
            # Add new dictionary entry to document that the imaging aperture is the
            # same as the pointing aperture
            self.exposure_tab['pointing_aperture'] = self.exposure_tab['aperture']

    def expand_for_detectors(self, input_dictionary):
        """Expand dictionary to have one entry per detector, rather than the
        one line per module that is in the input

        Parameters
        ----------
        input_dictionary : dict
            dictionary containing one entry per module

        Returns
        -------
        observation_dictionary : dict
            dictionary expanded to have one entry per detector
        """
        observation_dictionary = {}
        for key in input_dictionary:
            observation_dictionary[key] = []
        observation_dictionary['detector'] = []

        # Read in tables of aperture information
        nircam_subarray_file = os.path.join(self.config_path, 'NIRCam_subarray_definitions.list')
        nircam_apertures = utils.read_subarray_definition_file(nircam_subarray_file)

        for index, instrument in enumerate(input_dictionary['Instrument']):
            instrument = instrument.lower()
            if instrument == 'nircam' and input_dictionary['Mode'][index] != 'coron':
                # NIRCam case: Expand for detectors. Create one entry in each list for each
                # detector, rather than a single entry for 'ALL' or 'BSALL'

                # Determine module and subarray of the observation
                sub = input_dictionary['Subarray'][index]
                module = input_dictionary['Module'][index]
                if module == 'ALL':
                    module = ['A', 'B']
                else:
                    module = list(module)

                # Match up `sub` with aperture names in the aperture table
                # FULL matches up with the standard full frame imaging
                # apertures as well as full frame Grism apertures
                matches = np.where(nircam_apertures['Name'] == sub)[0]

                if len(matches) == 0:
                    raise ValueError("ERROR: aperture {} in not present in the subarray definition file {}"
                                     .format(sub, nircam_subarray_file))

                # Keep only apertures in the correct module
                matched_allmod_apertures = nircam_apertures['AperName'][matches].data
                matched_apertures = []
                for mod in module:
                    good = [ap for ap in matched_allmod_apertures if 'NRC{}'.format(mod) in ap]
                    matched_apertures.extend(good)
                if sub in ['FULL', 'SUB160', 'SUB320', 'SUB640', 'SUB64P', 'SUB160P', 'SUB400P', 'FULLP']:
                    mode = input_dictionary['Mode'][index]
                    template = input_dictionary['APTTemplate'][index]
                    if (sub == 'FULL'):
                        if mode in ['imaging', 'ts_imaging', 'wfss']:
                            # This block should catch full-frame observations
                            # in either imaging (including TS imaging) or
                            # wfss mode
                            matched_aps = np.array([ap for ap in matched_apertures if 'GRISM' not in ap if 'MASK' not in ap])
                            matched_apertures = []
                            detectors = []
                            for ap in matched_aps:
                                detectors.append(ap[3:5])
                                split = ap.split('_')
                                if len(split) == 3:
                                    ap_string = '{}_{}'.format(split[1], split[2])
                                elif len(split) == 2:
                                    ap_string = split[1]
                                matched_apertures.append(ap_string)

                        elif mode == 'ts_grism':
                            # This block should get Grism Time Series
                            # observations that use the full frame
                            matched_apertures = np.array([ap for ap in matched_apertures if 'GRISM' in ap])
                            filtered_aperture = input_dictionary['aperture'][index]
                            filtered_splits = filtered_aperture.split('_')
                            filtered_ap_no_det = '{}_{}'.format(filtered_splits[1], filtered_splits[2])
                            detectors = [filtered_splits[0][3:5]]

                            # Get correct apertures
                            apertures_to_add = []

                            final_matched_apertures = []
                            final_detectors = []
                            filtered_ap_det = filtered_aperture[3:5]
                            for ap in matched_apertures:
                                det = ap[3:5]
                                if ap == filtered_aperture or det != filtered_ap_det:
                                    final_matched_apertures.append(ap)
                                    final_detectors.append(det)
                            matched_apertures = final_matched_apertures
                            detectors = final_detectors
                    else:
                        # 'Standard' imaging subarrays: SUB320, SUB400P, etc
                        matched_apertures = [ap for ap in matched_apertures if sub in ap]
                        detectors = [ap.split('_')[0][3:5] for ap in matched_apertures]
                        matched_apertures = [sub] * len(detectors)

                elif 'SUBGRISM' in sub:
                    # This should catch only Grism Time Series observations
                    # and engineering imaging observations, which are the
                    # only 2 templates that can use SUBGRISM apertures
                    long_filter = input_dictionary['LongFilter'][index]
                    filter_dependent_apertures = [ap for ap in matched_apertures if len(ap.split('_')) == 3]

                    filtered_aperture = input_dictionary['aperture'][index]
                    filtered_splits = filtered_aperture.split('_')
                    filtered_ap_no_det = '{}_{}'.format(filtered_splits[1], filtered_splits[2])
                    detectors = [filtered_splits[0][3:5]]

                    # Get correct apertures
                    apertures_to_add = []

                    final_matched_apertures = []
                    final_detectors = []
                    filtered_ap_det = filtered_aperture[3:5]
                    for ap in matched_apertures:
                        det = ap[3:5]
                        if ap == filtered_aperture or det != filtered_ap_det:
                            final_matched_apertures.append(ap)
                            final_detectors.append(det)
                    matched_apertures = []
                    for ap in final_matched_apertures:
                        split = ap.split('_')
                        if len(split) == 3:
                            ap_string = '{}_{}'.format(split[1], split[2])
                        elif len(split) == 2:
                            ap_string = split[1]
                        matched_apertures.append(ap_string)
                    detectors = final_detectors
                else:
                    # TA, WFSC apertures
                    stripped_apertures = []
                    detectors = []
                    for ap in matched_apertures:
                        detectors.append(ap[3:5])
                        split = ap.split('_')
                        if len(split) == 3:
                            ap_string = '{}_{}'.format(split[1], split[2])
                        elif len(split) == 2:
                            ap_string = split[1]
                        stripped_apertures.append(ap_string)
                    matched_apertures = stripped_apertures

                full_apertures = ['NRC{}_{}'.format(det, sub) for det, sub in zip(detectors, matched_apertures)]

                # Add entries to observation dictionary
                num_entries = len(detectors)
                #observation_dictionary['Subarray'].extend(matched_apertures) extend? or replace?
                for key in input_dictionary:
                    #if key not in ['Subarray']:
                    if key not in ['aperture', 'detector', 'Subarray']:
                        observation_dictionary[key].extend(([input_dictionary[key][index]] * num_entries))
                observation_dictionary['detector'].extend(detectors)
                observation_dictionary['aperture'].extend(full_apertures)
                observation_dictionary['Subarray'].extend(matched_apertures)

            else:
                if instrument == 'niriss':
                    detectors = ['NIS']

                elif instrument == 'nirspec':
                    detectors = ['NRS']

                elif instrument == 'fgs':
                    if input_dictionary['APTTemplate'][index] == 'WfscGlobalAlignment':
                        guider_number = input_dictionary['aperture'][index][3]
                    elif input_dictionary['APTTemplate'][index] == 'FgsExternalCalibration':
                        guider_number = read_apt_xml.get_guider_number(self.input_xml, input_dictionary['obs_num'][index])
                    detectors = ['G{}'.format(guider_number)]

                elif instrument == 'miri':
                    detectors = ['MIR']

                elif instrument == 'nircam' and input_dictionary['Mode'][index] == 'coron':
                    # For coronagraphic observations, there is no need to expand for detectors.
                    # Coronographic observations will always use only a single detector, which
                    # we already know from the 'aperture' key in the input dictionary
                    detectors = [input_dictionary['aperture'][index][3:5]]

                    # Reset the mode of the coronagraphic observations to be imaging, since
                    # 'coron' is not a supported mode, and those running Mirage with these
                    # files currently have to manually switch the mode over to 'imaging'
                    input_dictionary['Mode'][index] = 'imaging'

                n_detectors = len(detectors)
                for key in input_dictionary:
                    observation_dictionary[key].extend(([input_dictionary[key][index]] * n_detectors))
                observation_dictionary['detector'].extend(detectors)

        """
        # Correct NIRCam aperture names for commissioning subarrays
        for index, instrument in enumerate(observation_dictionary['Instrument']):
            instrument = instrument.lower()
            if instrument == 'nircam':
                detector = observation_dictionary['detector'][index]
                sub = observation_dictionary['Subarray'][index]

                # this should probably better be handled by using the subarray_definitions file upstream
                if 'DHSPIL' in sub:
                    subarray, module = sub.split('DHSPIL')
                    subarray_size = subarray[3:]
                    detector = 'A3' if module == 'A' else 'B4'
                    aperture_name = 'NRC{}_DHSPIL_SUB{}'.format(detector, subarray_size)
                elif 'FP1' in sub:
                    subarray, module = sub.split('FP1')
                    subarray_size = subarray[3:]
                    detector = 'A3' if module == 'A' else 'B4'
                    aperture_name = 'NRC{}_FP1_SUB{}'.format(detector, subarray_size)
                else:
                    aperture_name = 'NRC' + detector + '_' + sub

                observation_dictionary['aperture'][index] = aperture_name
                observation_dictionary['detector'][index] = detector
        """

        return observation_dictionary

    def extract_grism_aperture(self, apertures, filter_name):
        """In the case of a Grism observation (WFSS or GRISM TSO), where a
        given crossing filter is used, find the appropriate aperture to
        use.

        Paramters
        ---------
        apertures : list
            List of possible grism apertures

        filter_name : str
            Name of crossing filter

        Returns
        -------
        apertures : list
            Modified list containig the correct aperture
        """
        filter_match = [True if filter_name in mtch else False for mtch in apertures]
        if any(filter_match):
            self.logger.debug('EXACT FILTER MATCH')
            self.logger.debud('{}'.format(filter_match))
            apertures = list(np.array(apertures)[filter_match])
        else:
            self.logger.debug('NO EXACT FILTER MATCH')
            filter_int = int(filter_name[1:4])
            aperture_int = np.array([int(ap.split('_')[-1][1:4]) for ap in apertures])
            wave_diffs = np.abs(aperture_int - filter_int)
            min_diff_index = np.where(wave_diffs == np.min(wave_diffs))[0]
            apertures = list(apertures[min_diff_index])

            self.logger.debug('{} {} {} {}'.format(filter_int, aperture_int, min_diff_index, apertures))

        return apertures

    def extract_value(self, line):
        """Extract text from xml line

        Parameters
        ----------
        line : str
            Line from xml file

        Returns
        -------
        line : str
            Text between > and < in the input line
        """
        gt = line.find('>')
        lt = line.find('<', gt)
        return line[gt + 1:lt]

    def filter_unsuppoted_obs_numbers(self, pointing_dict, skipped_obs):
        """Remove observations from the pointing dictionary that are for
        unsupported observing modes

        Parameters
        ----------
        pointing_dict : dict
            Dictionary of pointing information

        skipped_obs : list
            List of observation numbers that should be removed

        Returns
        -------
        pointing_dict : dict
            Dictionary with the appropriate entries removed
        """
        if skipped_obs is None:
            return pointing_dict

        for obnum in skipped_obs:
            #matches = pointing_dict['obs_num'] != obnum
            matches = [True if obnum != o else False for o in pointing_dict['obs_num']]

            for key in pointing_dict:
                values = np.array(pointing_dict[key])
                values = values[matches]
                pointing_dict[key] = list(values)

        return pointing_dict


    def full_path(self, in_path):
        """
        If the input path is not None, expand
        any environment variables and make an
        absolute path. Return the updated path.

        Parameters
        ----------
        in_path : str
            Path to be expanded

        Returns
        -------
        in_path : str
            Expanded, absolute path
        """
        if in_path.lower() == 'none':
            return in_path
        else:
            return os.path.abspath(os.path.expandvars(in_path))

    def get_pointing_info(self, file, propid=0, skipped_obs_from_xml=None, verbose=False):
        """Read in information from APT's pointing file.

        Parameters
        ----------
        file : str
            Name of APT-exported pointing file to be read

        propid : int
            Proposal ID number (integer). This is used to
            create various ID fields

        skipped_obs_from_xml : list
            List of observation numbers to be removed

        Returns
        -------
        pointing : dict
            Dictionary of pointing-related information

        TODO
        ----
            extract useful information from header?
            check visit numbers
            set parallel proposal number correctly

        """
        tar = []
        tile = []
        exp = []
        dith = []
        aperture = []
        pps_aperture = []
        targ1 = []
        targ2 = []
        ra = []
        dec = []
        basex = []
        basey = []
        dithx = []
        dithy = []
        v2 = []
        v3 = []
        idlx = []
        idly = []
        level = []
        type_str = []
        expar = []
        dkpar = []
        ddist = []
        observation_number = []
        visit_number = []
        visit_id = []
        visit_grp = []
        activity_id = []
        observation_label = []
        observation_id = []
        seq_id = []

        act_counter = 1
        with open(file) as f:
            for line in f:

                # Skip comments and new lines except for the line with the version of the PRD
                if (line[0] == '#') or (line in ['\n']) or ('=====' in line):

                    # Compare the version of the PRD from APT and pysiaf
                    if 'PRDOPSSOC' in line:
                        apt_prd_version = line.split(' ')[-2]
                        if apt_prd_version != JWST_PRD_VERSION:
                            self.logger.warning(('\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'
                                                 'The pointing file from APT was created using PRD version: {},\n'
                                                 'while the current installation of pysiaf uses PRD version: {}.\n'
                                                 'This inconsistency may lead to errors in source locations or\n'
                                                 'the WCS of simulated data if the apertures being simulated are\n'
                                                 'shifted between the two versions. We highly recommend using a\n'
                                                 'consistent version of the PRD between APT and pysiaf.\n'
                                                 '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'
                                                 .format(apt_prd_version, JWST_PRD_VERSION)))
                    else:
                        continue
                # Extract proposal ID
                elif line.split()[0] == 'JWST':
                    propid_header = line.split()[7]
                    try:
                        propid = int(propid_header)
                    except ValueError:
                        # adopt value passed to function
                        pass
                    if verbose:
                        self.logger.info('Extracted proposal ID {}'.format(propid))
                    continue

                elif (len(line) > 1):
                    elements = line.split()

                    # Look for lines that give visit/observation numbers
                    if line[0:2] == '* ':
                        paren = line.rfind('(')
                        if paren == -1:
                            obslabel = line[2:]
                            obslabel = obslabel.strip()
                        else:
                            obslabel = line[2:paren-1]
                            obslabel = obslabel.strip()
                        if (' (' in obslabel) and (')' in obslabel):
                            obslabel = re.split(r' \(|\)', obslabel)[0]

                    skip = False

                    if line[0:2] == '**':
                        v = elements[2]
                        obsnum, visitnum = v.split(':')
                        obsnum = str(obsnum).zfill(3)
                        visitnum = str(visitnum).zfill(3)
                        if (skip is True) and (verbose):
                            self.logger.info('Skipping observation {} ({})'.format(obsnum, obslabel))

                    try:
                        # Skip the line at the beginning of each
                        # visit that gives info on the target,
                        # but is not actually an observation
                        # These lines have 'Exp' values of 0,
                        # while observations have a value of 1
                        # (that I've seen so far)

                        # The exception to this rule is TA images.
                        # These have Exp values of 0. Sigh.


                        if ((int(elements[1]) > 0) & ('NRC' in elements[4]
                                                         or 'NIS' in elements[4]
                                                         or 'FGS' in elements[4]
                                                         or 'NRS' in elements[4]
                                                         or 'MIR' in elements[4])
                            ) or (('TA' in elements[4]) & ('NRC' in elements[4]
                                                           or 'NIS' in elements[4])):
                            if (elements[18] == 'PARALLEL') and ('MIRI' in elements[4]):
                                skip = True

                            if skip:
                                # act_counter += 1
                                continue
                            act = self.base36encode(act_counter)
                            activity_id.append(act)
                            observation_label.append(obslabel)
                            observation_number.append(obsnum)
                            visit_number.append(visitnum)
                            prop_5digit = "{0:05d}".format(int(propid))
                            vid = "{}{}{}".format(prop_5digit, obsnum, visitnum)
                            visit_id.append(vid)
                            # Visit group hard coded to 1. It's not clear how APT divides visits up into visit
                            # groups. For now just keep everything in a single visit group.
                            vgrp = '01'
                            visit_grp.append(vgrp)
                            # Parallel sequence is hard coded to 1 (Simulated instrument as prime rather than
                            # parallel) at the moment. Future improvements may allow the proper sequence
                            # number to be constructed.
                            seq = '1'
                            seq_id.append(seq)
                            tar.append(int(elements[0]))
                            tile.append(int(elements[1]))
                            exnum = str(elements[2]).zfill(5)
                            exp.append(exnum)
                            dith.append(int(elements[3]))
                            pps_aperture.append(elements[4])

                            ap = elements[4]
                            if ('GRISMR_WFSS' in elements[4]):
                                ap = ap.replace('GRISMR_WFSS', 'FULL')
                            elif ('GRISMC_WFSS' in elements[4]):
                                ap = ap.replace('GRISMC_WFSS', 'FULL')

                            aperture.append(ap)
                            targ1.append(int(elements[5]))
                            targ2.append(elements[6])
                            ra.append(elements[7])
                            dec.append(elements[8])
                            basex.append(elements[9])
                            basey.append(elements[10])
                            dithx.append(float(elements[11]))
                            dithy.append(float(elements[12]))
                            v2.append(float(elements[13]))
                            v3.append(float(elements[14]))
                            idlx.append(float(elements[15]))
                            idly.append(float(elements[16]))
                            level.append(elements[17])
                            type_str.append(elements[18])
                            expar.append(int(elements[19]))
                            dkpar.append(int(elements[20]))

                            # For the moment we assume that the instrument being simulated is not being
                            # run in parallel, so the parallel proposal number will be all zeros,
                            # as seen in the line below.
                            observation_id.append("V{}P{}{}{}{}".format(vid, '00000000', vgrp, seq, act))
                            # act_counter += 1

                    except ValueError as e:
                        if verbose:
                            self.logger.info('Skipping line:\n{}\nproducing error:\n{}'.format(line, e))
                        pass

        pointing = {'exposure': exp, 'dither': dith, 'aperture': aperture, 'pps_aperture': pps_aperture,
                    'targ1': targ1, 'targ2': targ2, 'ra': ra, 'dec': dec,
                    'basex': basex, 'basey': basey, 'dithx': dithx,
                    'dithy': dithy, 'v2': v2, 'v3': v3, 'idlx': idlx,
                    'idly': idly, 'obs_label': observation_label,
                    'obs_num': observation_number, 'visit_num': visit_number,
                    'act_id': activity_id, 'visit_id': visit_id, 'visit_group': visit_grp,
                    'sequence_id': seq_id, 'observation_id': observation_id}

        pointing = self.filter_unsuppoted_obs_numbers(pointing, skipped_obs_from_xml)
        return pointing


    def global_alignment_pointing(self, obs_dict):
        """Adjust the pointing dictionary information for global alignment
        observations. Some of the entries need to be changed from NIRCam to
        FGS. Remember that not all observations in the dictionary will
        necessarily be WfscGlobalAlignment template. Be sure the leave all
        other templates unchanged.

        Parameters
        ----------
        obs_dict : dict
            Dictionary of observation parameters, as returned from add_observation_info()

        Returns
        -------
        obs_dict : dict
            Dictionary with modified values for FGS pointing in Global Alignment templates
        """

        # We'll always be changing NIRCam to FGS, so set up the NIRCam siaf
        # instance outside of loop
        nrc_siaf = Siaf('nircam')['NRCA3_FULL']

        ga_index = np.array(obs_dict['APTTemplate']) == 'WfscGlobalAlignment'

        observation_numbers = np.unique(np.array(obs_dict['obs_num'])[ga_index])

        for obs_num in observation_numbers:
            obs_indexes = np.where(np.array(obs_dict['obs_num']) == obs_num)[0]

            # Get the subarray and aperture entries for the observation
            aperture_values = np.array(obs_dict['aperture'])[obs_indexes]
            subarr_values = np.array(obs_dict['Subarray'])[obs_indexes]

            # Subarray values, which come from the xml file, are correct. The aperture
            # values, which come from the pointing file, are not correct. We need to
            # copy over the FGS values from the Subarray column to the aperture column
            to_fgs = [True if 'FGS' in subarr else False for subarr in subarr_values]
            aperture_values[to_fgs] = subarr_values[to_fgs]
            fgs_aperture = aperture_values[to_fgs][0]

            all_aperture_values = np.array(obs_dict['aperture'])
            all_aperture_values[obs_indexes] = aperture_values
            obs_dict['aperture'] = all_aperture_values

            # Update the pointing info for the FGS exposures
            fgs = Siaf('fgs')[fgs_aperture]
            basex, basey = fgs.tel_to_idl(nrc_siaf.V2Ref, nrc_siaf.V3Ref)
            dithx = np.array(obs_dict['dithx'])[obs_indexes[to_fgs]]
            dithy = np.array(obs_dict['dithy'])[obs_indexes[to_fgs]]
            idlx = basex + dithx
            idly = basey + dithy

            basex_col = np.array(obs_dict['basex'])
            basey_col = np.array(obs_dict['basey'])
            idlx_col = np.array(obs_dict['idlx'])
            idly_col = np.array(obs_dict['idly'])

            basex_col[obs_indexes[to_fgs]] = basex
            basey_col[obs_indexes[to_fgs]] = basey
            idlx_col[obs_indexes[to_fgs]] = idlx
            idly_col[obs_indexes[to_fgs]] = idly

            obs_dict['basex'] = basex_col
            obs_dict['basey'] = basey_col
            obs_dict['idlx'] = idlx_col
            obs_dict['idly'] = idly_col

        return obs_dict

    def tight_dithers(self, input_dict):
        """
        In NIRCam, when the 'FULL' dither pattern is
        used, it is possible to set the number of primary
        dithers to '3TIGHT' rather than just a number
        (e.g. '3'). If the number of dithers is set to '3TIGHT'
        remove 'TIGHT' from the entries and leave
        only the number behind.

        Parameters
        ----------

        input_dict : dict
           Dictionary where each key points to a list containing
           observation details. For example, input_dict['PrimarDither']
           is a list of the number of primary dithers for all observations

        Returns
        -------

        input_dict : dict
            Updated dictionary where 'TIGHT' has been removed from
            PrimaryDither list
        """
        inlist = input_dict['PrimaryDithers']
        modlist = [v if 'TIGHT' not in v else v.strip('TIGHT') for v in inlist]
        input_dict['PrimaryDithers'] = modlist
        return input_dict

    def add_options(self, parser=None, usage=None):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage, description='Simulate JWST ramp')
        parser.add_argument("input_xml", help='XML file from APT describing the observations.')
        parser.add_argument("pointing_file", help='Pointing file from APT describing observations.')
        parser.add_argument("--output_csv", help="Name of output file containing list of observations.",
                            default=None)
        parser.add_argument("--observation_list_file", help=('Ascii file containing a list of '
                                                             'observations, times, and roll angles, '
                                                             'catalogs'), default=None)
        return parser


def get_entry(dict, entry_number):
    """Return a numbered entry from a dictionary that corresponds to the observataion_list.yaml.

    Parameters
    ----------
    dict
    entry_number

    Returns
    -------

    """
    entry_key = 'EntryNumber{}'.format(entry_number)
    for key, observation in dict.items():
        if entry_key in observation.keys():
            return observation[entry_key]


def get_filters(pointing_info):
    """Return a dictionary of instruments and filters contained within a
    pointing dictionary from an APT file.

    This function is aware that sometimes filters are installed in pupil wheels,
    and also that some wheels contain non-filter items such as the NIRCam
    weak lenses. It will return a list of those items that are indeed spectral bandpass filters.

    Parameters
    ----------
    pointing_info : dict
        This is the dictionary held in self.exposure_tab from the apt_inputs.APTInput
        class.

    Returns
    -------
    filters : dict
        Dictionary with instrument names as keys. The value for each is a
        list containing the names of the filters for that instrument present
        in the pointing dictionary.
    """
    instrument_list = set(pointing_info['Instrument'])
    filters = {}
    for inst in instrument_list:
        good = np.where(np.array(pointing_info['Instrument']) == inst.upper())[0]
        if inst.upper() == 'NIRCAM':
            filter_list = []
            short_filters = np.array(pointing_info['ShortFilter'])[good]
            long_filters = np.array(pointing_info['LongFilter'])[good]
            short_pupils = np.array(pointing_info['ShortPupil'])[good]
            long_pupils = np.array(pointing_info['LongPupil'])[good]

            for s_filt, s_pup, l_filt, l_pup in zip(short_filters, short_pupils, long_filters, long_pupils):
                if s_pup not in NIRCAM_UNSUPPORTED_PUPIL_VALUES:
                    filter_list.append('{}/{}'.format(s_filt, s_pup))
                if l_pup not in NIRCAM_UNSUPPORTED_PUPIL_VALUES:
                    if 'GRISM' in l_pup:
                        l_pup = 'CLEAR'
                    filter_list.append('{}/{}'.format(l_filt, l_pup))

        elif inst.upper() == 'FGS':
            filter_list = ['guider1', 'guider2']

        else:
            short_filters = np.array(pointing_info['FilterWheel'])[good]
            short_pupils = np.array(pointing_info['PupilWheel'])[good]

            short_filter_only = np.where(((short_pupils == 'CLEAR') | (short_pupils == 'CLEARP')))[0]
            filter_list = list(set(short_pupils))

            # Remove CLEAR and non-imaging elements if present
            niriss_clears = ['CLEAR', 'CLEARP', 'None', 'NRM', 'GR700XD', 'GR150R', 'GR150C']
            for clear in niriss_clears:
                if clear in filter_list:
                    filter_list.remove(clear)

            filter_list.extend(list(set(short_filters[short_filter_only])))

        filters[inst.upper()] = filter_list
    return filters


def make_start_times(obs_info, offline=False):
    """Create exposure start times for each entry in the observation dictionary.

    Parameters
    ----------
    obs_info : dict
        Dictionary of exposures. Development was around a dictionary containing
        APT xml-derived properties as well as pointing file properties. Should
        be before expanding to have one entry for each detector in each exposure.

    offline : bool
        Whether the class is being called with or without access to
        Mirage reference data. Used primarily for testing.

    Returns
    -------
    obs_info : dict
        Modified dictionary with observation dates and times added
    """
    logger = logging.getLogger('mirage.apt.apt_inputs')

    date_obs = []
    time_obs = []
    expstart = []
    nframe = []
    nskip = []
    namp = []

    # Read in file containing subarray definitions
    config_information = utils.organize_config_files(offline=offline)

    if 'epoch_start_date' in obs_info.keys():
        epoch_base_date = obs_info['epoch_start_date'][0]
    else:
        epoch_base_date = obs_info['Date'][0]

    base = Time(obs_info['epoch_start_date'][0])
    #base = Time(epoch_base_date + 'T' + epoch_base_time)
    base_date, base_time = base.iso.split()

    # Pick some arbirary overhead values
    act_overhead = 90  # seconds. (filter change)
    visit_overhead = 600  # seconds. (slew)

    # Get visit, activity_id, dither_id info for first exposure
    ditherid = obs_info['dither'][0]
    actid = obs_info['act_id'][0]
    visit = obs_info['visit_num'][0]
    obsid = obs_info['ObservationID'][0]
    exp = obs_info['exposure'][0]
    entry_num = obs_info['entry_number'][0]

    for i, instrument in enumerate(obs_info['Instrument']):
        # Get dither/visit
        # Files with the same activity_id should have the same start time
        # Overhead after a visit break should be large, smaller between
        # exposures within a visit
        next_actid = obs_info['act_id'][i]
        next_visit = obs_info['visit_num'][i]
        next_obsname = obs_info['obs_label'][i]
        next_ditherid = obs_info['dither'][i]
        next_obsid = obs_info['ObservationID'][i]
        next_exp = obs_info['exposure'][i]
        next_entry_num = obs_info['entry_number'][i]

        # Find the readpattern of the file
        readpatt = obs_info['ReadoutPattern'][i]
        groups = int(obs_info['Groups'][i])
        integrations = int(obs_info['Integrations'][i])

        if instrument.lower() in ['miri', 'nirspec']:
            nframe.append(0)
            nskip.append(0)
            namp.append(0)
            date_obs.append(base_date)
            time_obs.append(base_time)
            expstart.append(base.mjd)

        else:
            readpatt_def = config_information['global_readout_patterns'][instrument.lower()]
            subarray_def = config_information['global_subarray_definitions'][instrument.lower()]

            match2 = readpatt == readpatt_def['name']
            if np.sum(match2) == 0:
                raise RuntimeError(("WARNING!! Readout pattern {} not found in definition file."
                                    .format(readpatt)))

            # Now get nframe and nskip so we know how many frames in a group
            fpg = int(readpatt_def['nframe'][match2][0])
            spg = int(readpatt_def['nskip'][match2][0])
            nframe.append(fpg)
            nskip.append(spg)

            # Get the aperture name. For non-NIRCam instruments,
            # this is simply the obs_info['aperture']. But for NIRCam,
            # we need to be careful of entries like NRCBS_FULL, which is used
            # for observations using all 4 shortwave B detectors. In that case,
            # we need to build the aperture name from the combination of detector
            # and subarray name.
            aperture = obs_info['aperture'][i]

            # Get the number of amps from the subarray definition file
            match = aperture == subarray_def['AperName']

            if np.sum(match) == 0:
                if '_MASKLWB' in aperture or '_MASKSWB' in aperture:
                    apsplit = aperture.split('_')
                    no_filter = '{}_{}'.format(apsplit[0], apsplit[1])
                    match = no_filter == subarray_def['AperName']

            # needed for NIRCam case
            if np.sum(match) == 0:
                logger.info(('Aperture: {} does not match any entries in the subarray definition file. Guessing at the '
                             'aperture for the purpose of calculating the exposure time and number of amps.'.format(aperture)))
                sub = aperture.split('_')[1]
                aperture = [apername for apername, name in
                            np.array(subarray_def['AperName', 'Name']) if
                            (sub in apername) or (sub in name)]

                match = aperture == subarray_def['AperName']

                if len(aperture) > 1 or len(aperture) == 0 or np.sum(match) == 0:
                    raise ValueError('Cannot combine detector {} and subarray {}\
                                     into valid aperture name.'.format(det, sub))
                # We don't want aperture as a list
                aperture = aperture[0]

            # For grism tso observations, get the number of
            # amplifiers to use from the APT file.
            # For other modes, check the subarray def table.
            try:
                amp = int(obs_info['NumOutputs'][i])
            except ValueError:
                amp = subarray_def['num_amps'][match][0]

            # Default to amps=4 for subarrays that can have 1 or 4
            # if the number of amps is not defined. Hopefully we
            # should never enter this code block given the lines above.
            if amp == 0:
                amp = 4
                logger.info(('Aperture {} can be used with 1 or 4 readout amplifiers. Defaulting to use 4.'
                             'In the future this information should be made a user input.'.format(aperture)))
            namp.append(amp)

            # same activity ID
            # Remove this for now, since Mirage was not correctly
            # specifying activities. At the moment all exposures have
            # the same activity ID, which means we must allow the
            # the epoch_start_date to change even if the activity ID
            # does not. This will change back in the future when we
            # figure out more realistic activity ID values.
            #if next_actid == actid:
            #    # in this case, the start time should remain the same
            #    date_obs.append(base_date)
            #    time_obs.append(base_time)
            #    expstart.append(base.mjd)
            #    continue

            epoch_date = obs_info['epoch_start_date'][i]
            #epoch_time = copy.deepcopy(epoch_base_time0)

            # new epoch - update the base time
            if epoch_date != epoch_base_date:
                epoch_base_date = copy.deepcopy(epoch_date)
                #base = Time(epoch_base_date + 'T' + epoch_base_time)
                base = Time(obs_info['epoch_start_date'][i])
                base_date, base_time = base.iso.split()
                basereset = True
                date_obs.append(base_date)
                time_obs.append(base_time)
                expstart.append(base.mjd)
                actid = copy.deepcopy(next_actid)
                visit = copy.deepcopy(next_visit)
                obsid = copy.deepcopy(next_obsid)
                obsname = copy.deepcopy(next_obsname)
                ditherid = copy.deepcopy(next_ditherid)
                exp = copy.deepcopy(next_exp)
                entry_num = copy.deepcopy(next_entry_num)
                continue

            # new observation or visit (if a different epoch time has
            # not been provided)
            if ((next_obsid != obsid) | (next_visit != visit)):
                # visit break. Larger overhead
                overhead = visit_overhead
            elif ((next_actid > actid) & (next_visit == visit)):
                # This block should be updated when we have more realistic
                # activity IDs
                # same visit, new activity. Smaller overhead
                overhead = act_overhead
            elif ((next_ditherid != ditherid) & (next_visit == visit)):
                # same visit, new dither position. Smaller overhead
                overhead = act_overhead
            else:
                # same observation, activity, dither. Filter changes
                # will still fall in here, which is not accurate
                overhead = 0.  # Reset frame captured in exptime below

            # For cases where the base time needs to change
            # continue down here
            siaf_inst = obs_info['Instrument'][i].upper()
            siaf_obj = Siaf(siaf_inst)[aperture]

            # Calculate the readout time for a single frame
            frametime = utils.calc_frame_time(siaf_inst, aperture,
                                              siaf_obj.XSciSize, siaf_obj.YSciSize, amp)

            # Estimate total exposure time
            exptime = ((fpg + spg) * groups + fpg) * integrations * frametime

            if ((next_obsid == obsid) & (next_visit == visit) & (next_actid == actid) & (next_ditherid == ditherid) & (next_entry_num == entry_num)):
                # If we are in the same exposure (but with a different detector),
                # then we should keep the start time the same
                delta = TimeDelta(0., format='sec')
            else:
                # If we are moving on to the next exposure, activity, or visit
                # then move the start time by the expoure time of the current
                # exposure, plus the overhead
                delta = TimeDelta(exptime + overhead, format='sec')

            base += delta
            base_date, base_time = base.iso.split()

            # Add updated dates and times to the list
            date_obs.append(base_date)
            time_obs.append(base_time)
            expstart.append(base.mjd)

            # increment the activity ID and visit
            actid = copy.deepcopy(next_actid)
            visit = copy.deepcopy(next_visit)
            obsname = copy.deepcopy(next_obsname)
            ditherid = copy.deepcopy(next_ditherid)
            obsid = copy.deepcopy(next_obsid)
            exp = copy.deepcopy(next_exp)
            entry_num = copy.deepcopy(next_entry_num)

    obs_info['date_obs'] = date_obs
    obs_info['time_obs'] = time_obs
    obs_info['nframe'] = nframe
    obs_info['nskip'] = nskip
    obs_info['namp'] = namp
    return obs_info


def ra_dec_update(exposure_dict, siaf_instances, verbose=False):
    """Given the V2, V3 values for the reference locations associated
    with detector apertures, calculate corresponding RA, Dec.

    Parameters
    ----------
    exposure_dict : dict
        Dictionary of exposure parameters, like self.exposure_tab, after expanding for
        detectors

    siaf_instances : dict
        Dictionary of instrument level SIAF instances. Instrument names are the
        keys and the SIAF instances are the values

    Returns
    -------
    exposure_dict : dict
        Modified exposure dictionary with updated RA, Dec values for the pointing
    """
    #sw_grismts_apertures = ['NRCA1_GRISMTS256', 'NRCA1_GRISMTS128', 'NRCA1_GRISMTS64', 'NRCA1_GRISMTS',
    #                        'NRCA3_GRISMTS256', 'NRCA3_GRISMTS128', 'NRCA3_GRISMTS64', 'NRCA3_GRISMTS']

    #lw_grismts_apertures = ['NRCA5_GRISM256_F277W', 'NRCA5_GRISM128_F277W', 'NRCA5_GRISM64_F277W', 'NRCA5_GRISM_F277W',
    #                        'NRCA5_GRISM256_F322W2', 'NRCA5_GRISM128_F322W2', 'NRCA5_GRISM64_F322W2', 'NRCA5_GRISM_F322W2',
    #                        'NRCA5_GRISM256_F356W', 'NRCA5_GRISM128_F356W', 'NRCA5_GRISM64_F356W', 'NRCA5_GRISM_F356W',
    #                        'NRCA5_GRISM256_F444W', 'NRCA5_GRISM128_F444W', 'NRCA5_GRISM64_F444W', 'NRCA5_GRISM_F444W']

    #intermediate_lw_grismts_apertures = ['NRCA5_TAGRISMTS_SCI_F444W', 'NRCA5_TAGRISMTS_SCI_F322W2']

    aperture_ra = []
    aperture_dec = []
    intermediate_apertures = []

    lw_grismts_aperture = None
    for i in range(len(exposure_dict['Module'])):
        siaf_instrument = exposure_dict["Instrument"][i]
        aperture_name = exposure_dict['aperture'][i]
        pointing_ra = float(exposure_dict['ra'][i])
        pointing_dec = float(exposure_dict['dec'][i])
        pointing_v2 = float(exposure_dict['v2'][i])
        pointing_v3 = float(exposure_dict['v3'][i])

        # When we run across a LW grism TS aperture, save
        # the aperture name, because we'll need it when looking
        # at the accompanying SW apertuers to follow. THIS
        # RELIES ON THE LW ENTRY COMING BEFORE THE SW ENTRIES.
        # Also note that F277W, F322W2, and F356W all use the same
        # intermeidate aperture, since their nominal apertures
        # (i.e. NRCA5_GRISM64_F322W2, NRCA5_GRISM64_F277W) are all
        # exactly the same as well.
        if aperture_name in NIRCAM_LW_GRISMTS_APERTURES:
            lw_grismts_aperture = copy.deepcopy(aperture_name)
            lw_intermediate_aperture = utils.get_lw_grism_tso_intermeidate_aperture(aperture_name)
            intermediate_apertures.append(lw_intermediate_aperture)
        elif aperture_name in NIRCAM_SW_GRISMTS_APERTURES:
            intermediate_apertures.append(lw_intermediate_aperture)
        else:
            intermediate_apertures.append(None)

        if 'pav3' in exposure_dict.keys():
            pav3 = float(exposure_dict['pav3'][i])
        else:
            pav3 = float(exposure_dict['PAV3'][i])

        telescope_roll = pav3

        aperture = siaf_instances[siaf_instrument][aperture_name]

        if 'NRCA5_GRISM' in aperture_name and 'WFSS' not in aperture_name:
            ra = pointing_ra
            dec = pointing_dec
            #exposure_dict['aperture'][i] = lw_intermediate_aperture  This does not work, because the intermediate aperture is a small box, and mirage happily
            #        creates a small square subarray when the seed image generator is run with this. we need something in catalog seed generator that knows
            #        about the link between the grismts apertures and the intermediates and uses the intermediate in the Siaf instance. Also, there are
            #        currently no intermediate apertures in pysiaf for the F277W, F356W filters. What happens (in APT/reality) in those cases?
        else:
            if aperture_name in NIRCAM_SW_GRISMTS_APERTURES:
                # Special case. When looking at grism time series observation
                # we force the pointing to be at the reference location of the
                # LW *intermediate* aperture, rather than paying attention to
                # the V2, V3 in the pointing file. V2, V3 from the intermediate
                # aperture is where the source would land on the detector if
                # the grism were not in the beam. This is exactly what we want
                # for the SW detectors, where this is no grism.

                # Generate an attitude matrix from this and
                # use to get the RA, Dec in the SW apertures
                lw_gts = siaf_instances[siaf_instrument][lw_intermediate_aperture]
                pointing_v2 = lw_gts.V2Ref
                pointing_v3 = lw_gts.V3Ref
                aperture.V2Ref = pointing_v2
                aperture.V3Ref = pointing_v3

            local_roll, attitude_matrix, fullframesize, subarray_boundaries = \
                siaf_interface.get_siaf_information(siaf_instances[siaf_instrument], aperture_name,
                                                    pointing_ra, pointing_dec, telescope_roll,
                                                    v2_arcsec=pointing_v2, v3_arcsec=pointing_v3)

            # Calculate RA, Dec of reference location for the detector
            # Add in any offsets from the pointing file in the BaseX, BaseY columns
            ra, dec = rotations.pointing(attitude_matrix, aperture.V2Ref, aperture.V3Ref)

        aperture_ra.append(ra)
        aperture_dec.append(dec)

    exposure_dict['grismts_intermediate_aperture'] = intermediate_apertures
    exposure_dict['ra_ref'] = aperture_ra
    exposure_dict['dec_ref'] = aperture_dec
    return exposure_dict



# if __name__ == '__main__':
#
#     usagestring = 'USAGE: apt_inputs.py NIRCam_obs.xml NIRCam_obs.pointing'
#
#     input = AptInput()
#     parser = input.add_options(usage=usagestring)
#     args = parser.parse_args(namespace=input)
#     input.create_input_table()
