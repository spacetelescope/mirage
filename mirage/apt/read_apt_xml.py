# ! /usr/bin/env python


import copy
import os
import re
from collections import OrderedDict

from lxml import etree
from astropy.io import ascii
import numpy as np

from ..utils.utils import append_dictionary

APT_DIR = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
PACKAGE_DIR = os.path.dirname(APT_DIR)


class ReadAPTXML():
    """Class to open and parse XML files from APT. Can read templates for
    NircamImaging, NircamEngineeringImaging, WfscCommissioning,
    WfscGlobaLAlignment, WfscCoarsePhasing, WfscFinePhasing (incomplete),
    and NircamWfss modes.

    Attributes
    ----------
    apt: str
        APT namespace for XML files
    APTObservationParams: dict
        Dictionary of APT parameters that accumulates all parameters in all
        tiles and all observation. Passed out to be further parsed in the
        apt_inputs script
    obs_tuple_list: list
        Compiled list of all the parameters for all tiles in a single observation
    """

    def __init__(self):
        # Define the APT namespace
        self.apt = '{http://www.stsci.edu/JWST/APT}'

        # Set up dictionary of observation parameters to be populated
        ProposalParams_keys = ['PI_Name', 'Proposal_category', 'ProposalID',
                               'Science_category', 'Title']
        ObsParams_keys = ['Module', 'Subarray', 'Instrument',
                          'PrimaryDitherType', 'PrimaryDithers', 'SubpixelPositions',
                          'SubpixelDitherType', 'CoordinatedParallel', 'ParallelInstrument',
                          'ObservationID', 'TileNumber', 'APTTemplate',
                          'ApertureOverride', 'ObservationName',
                          'DitherPatternType', 'ImageDithers', # NIRISS
                          'number_of_dithers', # uniform name across instruments
                          'FiducialPointOverride',
                          ]
        FilterParams_keys = ['ShortFilter', 'LongFilter', 'ShortPupil', 'LongPupil',
                             'ReadoutPattern', 'Groups', 'Integrations',
                             'FilterWheel', 'PupilWheel' # for NIRISS
                             ]
        OtherParams_keys = ['Mode', 'Grism',
                            'IntegrationsShort', 'GroupsShort', 'Dither', # MIRI
                            'GroupsLong', 'ReadoutPatternShort', 'IntegrationsLong',
                            'Exposures', 'Wavelength', 'ReadoutPatternLong', 'Filter',
                            'EtcIdLong', 'EtcIdShort', 'EtcId',
                            ]

        self.APTObservationParams_keys = ProposalParams_keys + ObsParams_keys + \
            FilterParams_keys + OtherParams_keys
        self.APTObservationParams = {}
        for key in self.APTObservationParams_keys:
            self.APTObservationParams[key] = []
        self.empty_exposures_dictionary = copy.deepcopy(self.APTObservationParams)
        self.observation_info = OrderedDict()

    def read_xml(self, infile, verbose=False):
        """Read in the .xml file from APT, and output dictionary of parameters.

        Arguments
        ---------
        infile (str):
            Path to input .xml file

        Returns
        -------
        dict:
            Dictionary with extracted observation parameters

        Raises
        ------
        ValueError:
            If an .xml file is provided that includes an APT template that is not
            supported.
            If the .xml file includes a fiducial pointing override with an
            unknown subarray specification
        """

        # Open XML file, get element tree of the APT proposal
        with open(infile) as f:
            tree = etree.parse(f)

        # Get high-level information: proposal info - - - - - - - - - - - - - -

        # Set default values
        propid_default = 42
        proptitle_default = 'Looking for my towel'
        scicat_default = 'Planets and Planet Formation'
        piname_default = 'D.N. Adams'
        propcat_default = 'GO'

        # Get just the element with the proposal information
        proposal_info = tree.find(self.apt + 'ProposalInformation')

        # Title
        try:
            prop_title = proposal_info.find(self.apt + 'Title').text
        except:
            prop_title = proptitle_default

        # Proposal ID
        try:
            prop_id = '{:05d}'.format(np.int(proposal_info.find(self.apt + 'ProposalID').text))
        except:
            prop_id = '{:05d}'.format(propid_default)

        # Proposal Category
        try:
            prop_category = proposal_info.find(self.apt + 'ProposalCategory')[0]
            prop_category = etree.QName(prop_category).localname
        except:
            prop_category = propcat_default

        # Science Category
        try:
            science_category = proposal_info.find(self.apt + 'ScientificCategory').text
        except:
            science_category = scicat_default

        # Principal Investigator Name
        try:
            pi_firstname = proposal_info.find('.//' + self.apt + 'FirstName').text
            pi_lastname = proposal_info.find('.//' + self.apt + 'LastName').text
            pi_name = ' '.join([pi_firstname, pi_lastname])
        except:
            pi_name = piname_default

        # Get parameters for each observation  - - - - - - - - - - - - - - - -

        # Find all observations (but use only those that use NIRCam or are WFSC)
        observation_data = tree.find(self.apt + 'DataRequests')
        observation_list = observation_data.findall('.//' + self.apt + 'Observation')

        # Loop through observations, get parameters
        for i_obs, obs in enumerate(observation_list):
            observation_number = np.int(obs.find(self.apt + 'Number').text)

            # Create empty list that will be populated with a tuple of parameters
            # for every observation
            self.obs_tuple_list = []

            # Determine what template is used for the observation
            template = obs.find(self.apt + 'Template')[0]
            template_name = etree.QName(template).localname

            # Are all the templates in the XML file something that we can handle?
            known_APT_templates = ['NircamImaging', 'NircamWfss', 'WfscCommissioning',
                                   'NircamEngineeringImaging', 'WfscGlobalAlignment',
                                   'WfscCoarsePhasing', 'WfscFinePhasing',
                                   'NirissExternalCalibration', 'NirissWfss',  # NIRISS
                                   'NirspecImaging', 'NirspecInternalLamp',  # NIRSpec
                                   'MiriMRS',  # MIRI
                                   'FgsExternalCalibration',  # FGS
                                   'NircamDark', 'NirissDark',  # Darks
                                   'NircamTimeSeries',
                                   ]
            if template_name not in known_APT_templates:
                # If not, turn back now.
                raise ValueError('No protocol written to read {} template.'.format(template_name))

            # Get observation label
            label_ele = obs.find(self.apt + 'Label')
            if label_ele is not None:
                label = label_ele.text
                if (' (' in label) and (')' in label):
                    label = re.split(r' \(|\)', label)[0]
            else:
                label = 'None'

            if verbose:
                print('+'*100)
                print('Observation `{}` labelled `{}` uses template `{}`'.format(observation_number, label, template_name))
                number_of_entries = len(self.APTObservationParams['Instrument'])
                print('APTObservationParams Dictionary holds {} entries before reading template'.format(number_of_entries))

            # Get coordinated parallel
            coordparallel = obs.find(self.apt + 'CoordinatedParallel').text
            CoordinatedParallelSet = None
            if coordparallel == 'true':
                try:
                    CoordinatedParallelSet = obs.find(self.apt + 'CoordinatedParallelSet').text
                except AttributeError:
                    raise RuntimeError('Program does not specify parallels correctly.')
                    # if verbose:
                    #     print('No CoordinatedParallelSet found. Set to default.')

            try:
                obs_label = obs.find(self.apt + 'Label').text
            except AttributeError:
                # label tag not present
                obs_label = 'Observation 1'

            # extract visit numbers
            visit_numbers = [np.int(element.items()[0][1]) for element in obs if
                             element.tag.split(self.apt)[1] == 'Visit']

            prop_params = [pi_name, prop_id, prop_title, prop_category,
                           science_category, coordparallel, observation_number, obs_label]

            proposal_parameter_dictionary = {'PI_Name': pi_name, 'ProposalID': prop_id,
                                             'Title': prop_title,
                                             'Proposal_category': prop_category,
                                             'Science_category': science_category,
                                             'CoordinatedParallel': coordparallel,
                                             'ObservationID': observation_number,
                                             'ObservationName': obs_label,
                                             }

            if template_name in ['NircamImaging', 'NircamEngineeringImaging', 'NirissExternalCalibration',
                                 'NirspecImaging', 'MiriMRS', 'FgsExternalCalibration']:
                exposures_dictionary = self.read_generic_imaging_template(template, template_name, obs,
                                                                          proposal_parameter_dictionary,
                                                                          verbose=verbose)
                if coordparallel == 'true':
                    parallel_template_name = etree.QName(obs.find(self.apt + 'FirstCoordinatedTemplate')[0]).localname
                    if parallel_template_name in ['MiriImaging']:
                        pass
                    else:
                        parallel_exposures_dictionary = self.read_parallel_exposures(obs, exposures_dictionary,
                                                                                     proposal_parameter_dictionary,
                                                                                     verbose=verbose)
                        exposures_dictionary = append_dictionary(exposures_dictionary, parallel_exposures_dictionary, braid=True)

            # If template is WFSC Commissioning
            elif template_name in ['WfscCommissioning']:
                exposures_dictionary, num_WFCgroups = self.read_commissioning_template(template, template_name, obs, prop_params)

            # If template is WFSC Global Alignment
            elif template_name in ['WfscGlobalAlignment']:
                exposures_dictionary, n_exp = self.read_globalalignment_template(template, template_name, obs, prop_params)

            # If template is WFSC Coarse Phasing
            elif template_name in ['WfscCoarsePhasing']:
                exposures_dictionary, n_tiles_phasing = self.read_coarsephasing_template(template, template_name, obs, prop_params)

            # If template is WFSC Fine Phasing
            elif template_name in ['WfscFinePhasing']:
                exposures_dictionary, n_tiles_phasing = self.read_finephasing_template(template, template_name, obs, prop_params)

            # If template is WFSS
            elif template_name == 'NircamWfss':
                exposures_dictionary = self.read_wfss_template(template, template_name, obs, prop_params)
                print('Looking in APT, it seems that NIRCam WFSS cannot have parallel observations')
                print('but this check was for a science porposal, and not COM or CAL')

            elif template_name == 'NirissWfss':
                exposures_dictionary, exp_dictionary_length = self.read_niriss_wfss_template(template, template_name, obs, proposal_parameter_dictionary)

                if coordparallel == 'true':
                    parallel_template_name = etree.QName(obs.find(self.apt + 'FirstCoordinatedTemplate')[0]).localname
                    if parallel_template_name in ['MiriImaging']:
                        pass
                    elif parallel_template_name in ['NircamImaging']:
                        parallel_exposures_dictionary = self.read_parallel_exposures(obs, exposures_dictionary,
                                                                                     proposal_parameter_dictionary,
                                                                                     verbose=verbose,
                                                                                     exposures_dictionary_length=exp_dictionary_length)
                        exposures_dictionary = append_dictionary(exposures_dictionary, parallel_exposures_dictionary, braid=True)
                    else:
                        raise ValueError('Parallel template {} (with primary template {}) not supported.'
                                         .format(parallel_template_name, template_name))
            else:
                print('SKIPPED: Observation `{}` labelled `{}` uses template `{}`'.format(observation_number,
                                                                                          label,
                                                                                          template_name))
                continue

            if verbose:
                print('Dictionary read from template has {} entries.'.format(len(exposures_dictionary['Instrument'])))
                print(exposures_dictionary['number_of_dithers'])
                # print(exposures_dictionary['Instrument'])

            # set default number of dithers to one, for downstream processing
            for i, n_dither in enumerate(exposures_dictionary['number_of_dithers']):
                if (template_name == 'NircamEngineeringImaging') and (n_dither == '2PLUS'):
                    exposures_dictionary['number_of_dithers'][i] = '2'
                elif int(n_dither) == 0:
                    exposures_dictionary['number_of_dithers'][i] = '1'

            # add the exposure dictionary to the main dictionary
            self.APTObservationParams = append_dictionary(self.APTObservationParams,
                                                          exposures_dictionary)

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            # Now we need to look for mosaic details, if any
            mosaic_tiles = obs.findall('.//' + self.apt + 'MosaicTiles')

            # count only tiles that are included
            tile_state = np.array([mosaic_tiles[i].find('.//' + self.apt + 'TileState').text for i in range(len(mosaic_tiles))])
            n_tiles = np.sum(tile_state=='Tile Included')

            label = obs_label

            if verbose:
                print("Found {} tile(s) for observation {} {}".format(n_tiles, observation_number, label))
                print('Found {} visits with numbers: {}'.format(len(visit_numbers), visit_numbers))

            if n_tiles > 1:
                for i in range(n_tiles - 1):
                    self.APTObservationParams = append_dictionary(self.APTObservationParams, exposures_dictionary)

            self.observation_info[observation_number] = {}
            self.observation_info[observation_number]['visit_numbers'] = visit_numbers

            if verbose:
                number_of_entries_after = len(self.APTObservationParams['Instrument'])
                print('APTObservationParams Dictionary holds {} entries after reading template ({:+d} entries)'
                      .format(number_of_entries_after, number_of_entries_after-number_of_entries))

        if verbose:
            print('Finished reading APT xml file.')
            print('+'*100)
        return self.APTObservationParams

    def add_exposure(self, dictionary, tup):
        # add an exposure to the dictionary
        dictionary['PI_Name'].append(tup[0])
        dictionary['ProposalID'].append(tup[1])
        dictionary['Title'].append(tup[2])
        dictionary['Proposal_category'].append(tup[3])
        dictionary['Science_category'].append(tup[4])
        dictionary['Mode'].append(tup[5])
        dictionary['Module'].append(tup[6])
        dictionary['Subarray'].append(tup[7])
        dictionary['PrimaryDitherType'].append(tup[8])
        dictionary['PrimaryDithers'].append(tup[9])
        dictionary['SubpixelDitherType'].append(tup[10])
        dictionary['SubpixelPositions'].append(tup[11])
        dictionary['ShortFilter'].append(tup[12])
        dictionary['LongFilter'].append(tup[13])
        dictionary['ReadoutPattern'].append(tup[14])
        dictionary['Groups'].append(tup[15])
        dictionary['Integrations'].append(tup[16])
        dictionary['ShortPupil'].append(tup[17])
        dictionary['LongPupil'].append(tup[18])
        dictionary['Grism'].append(tup[19])
        dictionary['CoordinatedParallel'].append(tup[20])
        dictionary['ObservationID'].append(tup[21])
        dictionary['TileNumber'].append(tup[22])
        dictionary['APTTemplate'].append(tup[23])
        dictionary['Instrument'].append(tup[24])
        dictionary['ObservationName'].append(tup[25])
        return dictionary

    def read_generic_imaging_template(self, template, template_name, obs, proposal_parameter_dictionary, verbose=False, parallel=False):
        """Read imaging template content regardless of instrument.

        Save content to object attributes. Support for coordinated parallels is included.

        Parameters
        ----------
        template : etree xml element
            xml content of template
        template_name : str
            name of the template
        obs : etree xml element
            xml content of observation
        proposal_parameter_dictionary : dict
            Dictionary of proposal parameters to extract from template

        Returns
        -------
        exposures_dictionary : OrderedDict
            Dictionary containing relevant exposure parameters

        """
        if parallel:
            # boolean indicating which instrument is not prime but parallel
            parallel_instrument = True
            if template_name == 'FgsExternalCalibration':
                instrument = 'FGS'
            elif template_name == 'MiriImaging':
                instrument = 'MIRI'
            elif template_name == 'NirissImaging':
                instrument = 'NIRISS'
            elif template_name == 'NircamImaging':
                instrument = 'NIRCAM'
            prime_instrument = obs.find(self.apt + 'Instrument').text
            print('IN READ_GENERIC_IMAGING_TEMPLATE')
            print('Prime: {}   Parallel: {}'.format(prime_instrument, instrument))
            prime_template = obs.find(self.apt + 'Template')[0]
            prime_template_name = etree.QName(prime_template).localname
            prime_ns = "{{{}/Template/{}}}".format(self.apt.replace('{','').replace('}',''), prime_template_name)
            print('PRIME TEMPLATE NAME IS: {}'.format(prime_template_name))
        else:
            instrument = obs.find(self.apt + 'Instrument').text
            parallel_instrument = False
            prime_instrument = instrument

        # if instrument.lower() in 'miri nirspec':
        #     return {}

        exposures_dictionary = copy.deepcopy(self.empty_exposures_dictionary)
        ns = "{{{}/Template/{}}}".format(self.apt.replace('{','').replace('}',''), template_name)

        DitherPatternType = None

        if ((prime_instrument in ['NIRCAM', 'FGS']) or
           (prime_instrument == 'NIRISS' and prime_template_name == 'NirissWfss')):
            dither_key_name = 'PrimaryDithers'
        elif prime_instrument in ['NIRISS', 'MIRI', 'NIRSPEC']:
            dither_key_name = 'ImageDithers'

        # number of dithers defaults to 1
        number_of_dithers = 1

        if instrument.lower() == 'nircam':
            # NIRCam uses FilterConfig structure to specifiy exposure parameters

            # store Module, Subarray, ... fields
            observation_dict = {}
            for field in template:
            #     print('{} {}'.format(field.tag, field.text))
                # if 'Filters ' in key:
                #     continue
                key = field.tag.split(ns)[1]
                value = field.text
                observation_dict[key] = value

            # Determine if there is an aperture override
            override = obs.find('.//' + self.apt + 'FiducialPointOverride')
            FiducialPointOverride = True if override is not None else False

            if dither_key_name in observation_dict.keys():
                number_of_dithers = observation_dict[dither_key_name]
                if observation_dict['SubpixelDitherType'] in ['3-POINT-WITH-MIRI-F770W']:
                    number_of_dithers = str(np.int(number_of_dithers) * 3)
                elif observation_dict['SubpixelDitherType'] in ['STANDARD']:
                    number_of_dithers = str(np.int(number_of_dithers) * np.int(observation_dict['SubpixelPositions']))
            elif parallel and prime_instrument == 'NIRISS' and prime_template_name == 'NirissWfss':
                observation_dict['PrimaryDithers'] = prime_template.find(prime_ns + dither_key_name).text
                observation_dict['DitherSize'] = prime_template.find(prime_ns + 'DitherSize').text
                number_of_dithers = observation_dict['PrimaryDithers'][0]
            else:
                print('Element {} not found, use default value.'.format(dither_key_name))

            # Find filter parameters for all filter configurations within obs
            filter_configs = template.findall('.//' + ns + 'FilterConfig')
            # loop over filter configurations
            for filter_config_index, filter_config in enumerate(filter_configs):
                filter_config_dict = {}
                # print('Filter config index {}'.format(filter_config_index))
                for element in filter_config:
                    key = element.tag.split(ns)[1]
                    value = element.text
                    # if verbose:
                    #     print('{} {}'.format(key, value))
                    if key == 'ShortFilter':
                        if ' + ' in value:
                            split_ind = value.find(' + ')
                            ShortPupil = value[0:split_ind]
                            value = value[split_ind + 1:]
                        else:
                            ShortPupil = 'CLEAR'
                        filter_config_dict['ShortPupil'] = ShortPupil
                    elif key == 'LongFilter':
                        if ' + ' in value:
                            split_ind = value.find(' + ')
                            LongPupil = value[0:split_ind]
                            value = value[split_ind + 1:]
                        else:
                            LongPupil = 'CLEAR'
                        filter_config_dict['LongPupil'] = LongPupil

                    filter_config_dict[key] = value

                for key in self.APTObservationParams_keys:
                    if key in filter_config_dict.keys():
                        value = filter_config_dict[key]
                    elif key in observation_dict.keys():
                        value = observation_dict[key]
                    elif key in proposal_parameter_dictionary.keys():
                        value = proposal_parameter_dictionary[key]
                    elif key == 'Instrument':
                        value = instrument
                    elif key == 'ParallelInstrument':
                        value = parallel_instrument
                    elif key == 'number_of_dithers':
                        value = str(number_of_dithers)
                    elif key == 'FiducialPointOverride':
                        value = str(FiducialPointOverride)
                    elif key == 'APTTemplate':
                        value = template_name
                    else:
                        value = str(None)

                    if (key == 'Mode'):
                        value = 'imaging'

                    exposures_dictionary[key].append(value)

        else:
            for element in template:
                element_tag_stripped = element.tag.split(ns)[1]
                # if verbose:
                #     print('{} {}'.format(element_tag_stripped, element.text))
                # loop through exposures and collect exposure parameters
                if element_tag_stripped == 'DitherPatternType':
                    DitherPatternType = element.text
                elif element_tag_stripped == 'ImageDithers':
                    # ImageDithers = element.text
                    number_of_dithers = element.text
                elif element_tag_stripped == 'PrimaryDithers':
                    number_of_dithers = element.text
                elif element_tag_stripped == 'Dithers':
                    DitherPatternType = element.find(ns + 'MrsDitherSpecification').find(ns + 'DitherType').text
                    number_of_dithers = int(DitherPatternType[0])

                # Determine if there is an aperture override
                override = obs.find('.//' + self.apt + 'FiducialPointOverride')
                FiducialPointOverride = True if override is not None else False
                #
                # observation_dict['FiducialPointOverride'] = str(FiducialPointOverride)

                # Different SI conventions of how to list exposure parameters
                if ((instrument.lower() == 'niriss') and (element_tag_stripped == 'ExposureList')) | \
                        ((instrument.lower() == 'fgs') and (element_tag_stripped == 'Exposures'))| \
                        ((instrument.lower() == 'miri') and (element_tag_stripped == 'ExposureList'))| \
                        ((instrument.lower() == 'nirspec') and (element_tag_stripped == 'Exposures')):
                    for exposure in element.findall(ns + 'Exposure'):
                        exposure_dict = {}
                        exposure_dict['DitherPatternType'] = DitherPatternType
                        if number_of_dithers is None:
                            number_of_dithers = 1
                        exposure_dict[dither_key_name] = np.int(number_of_dithers)
                        for exposure_parameter in exposure:
                            parameter_tag_stripped = exposure_parameter.tag.split(ns)[1]
                            # if verbose:
                            #     print('{} {}'.format(parameter_tag_stripped, exposure_parameter.text))
                            exposure_dict[parameter_tag_stripped] = exposure_parameter.text

                        exposure_dict['number_of_dithers'] = exposure_dict[dither_key_name]
                        # fill dictionary to return
                        for key in self.APTObservationParams_keys:
                            if key in exposure_dict.keys():
                                value = exposure_dict[key]
                                # print(key)
                            elif key in proposal_parameter_dictionary.keys():
                                value = proposal_parameter_dictionary[key]
                                # print(key)
                            elif key == 'Instrument':
                                value = instrument
                            elif key == 'ParallelInstrument':
                                value = parallel_instrument
                            elif key == 'number_of_dithers':
                                value = str(number_of_dithers)
                            elif key == 'FiducialPointOverride':
                                value = str(FiducialPointOverride)
                            elif key == 'APTTemplate':
                                value = template_name
                            else:
                                value = str(None)

                            if (key in ['PrimaryDithers', 'ImageDithers']) and ((value is None) or (value == 'None')):
                                value = '1'

                            if (key == 'Mode'):# and (template_name in ['NirissExternalCalibration', 'FgsExternalCalibration']):
                                value = 'imaging'

                            exposures_dictionary[key].append(value)

                        # add keys that were not defined in self.APTObservationParams_keys
                        # (to be fixed in Class.__init__ later )
                        for key in exposure_dict.keys():
                            if key not in self.APTObservationParams_keys:
                                # if key not yet present, create entry
                                if key not in exposures_dictionary.keys():
                                    exposures_dictionary[key] = [str(exposure_dict[key])]
                                else:
                                    exposures_dictionary[key].append(str(exposure_dict[key]))

        for key in exposures_dictionary.keys():
            if type(exposures_dictionary[key]) is not list:
                exposures_dictionary[key] = list(exposures_dictionary[key])

        # make sure all list items in the returned dictionary have the same length
        for key, item in exposures_dictionary.items():
            if len(item) == 0:
                exposures_dictionary[key] = [0] * len(exposures_dictionary['Instrument'])

        return exposures_dictionary

    def read_imaging_template(self, template, template_name, obs, prop_params):
        """Read NIRCam imaging template.

        Parameters
        ----------
        template
        template_name
        obs
        prop_params

        Returns
        -------

        """
        # Get proposal parameters
        pi_name, prop_id, prop_title, prop_category, science_category, coordparallel, i_obs, obs_label = prop_params

        # dictionary that holds the content of this observation only
        exposures_dictionary = copy.deepcopy(self.empty_exposures_dictionary)

        # Set namespace
        ns = "{{http://www.stsci.edu/JWST/APT/Template/{}}}".format(template_name)
        instrument = obs.find(self.apt + 'Instrument').text

        # Set parameters that are constant for all imaging obs
        #typeflag = template_name
        typeflag = 'imaging'
        grismval = 'N/A'
        short_pupil = 'CLEAR'

        # Find observation-specific parameters
        mod = template.find(ns + 'Module').text
        subarr = template.find(ns + 'Subarray').text
        pdithtype = template.find(ns + 'PrimaryDitherType').text

        # Determine if there is an aperture override
        mod, subarr = self.check_for_aperture_override(obs, mod, subarr, i_obs)

        try:
            pdither = template.find(ns + 'PrimaryDithers').text
        except:
            pdither = '1'

        sdithtype = template.find(ns + 'SubpixelDitherType').text

        try:
            sdither = template.find(ns + 'SubpixelPositions').text
        except:
            try:
                stemp = template.find(ns + 'CoordinatedParallelSubpixelPositions').text
                sdither = np.int(stemp[0])
            except:
                sdither = '1'

        # Find filter parameters for all filter configurations within obs
        filter_configs = template.findall('.//' + ns + 'FilterConfig')

        for filt in filter_configs:
            sfilt = filt.find(ns + 'ShortFilter').text
            lfilt = filt.find(ns + 'LongFilter').text
            rpatt = filt.find(ns + 'ReadoutPattern').text
            grps = filt.find(ns + 'Groups').text
            ints = filt.find(ns + 'Integrations').text

            # Separate pupil and filter in case of filter that is
            # mounted in the pupil wheel
            if ' + ' in sfilt:
                split_ind = sfilt.find(' + ')
                short_pupil = sfilt[0:split_ind]
                sfilt = sfilt[split_ind + 1:]
            else:
                short_pupil = 'CLEAR'

            if ' + ' in lfilt:
                p = lfilt.find(' + ')
                long_pupil = lfilt[0:p]
                lfilt = lfilt[p + 1:]
            else:
                long_pupil = 'CLEAR'

            # Add all parameters to dictionary
            tup_to_add = (pi_name, prop_id, prop_title, prop_category,
                          science_category, typeflag, mod, subarr, pdithtype,
                          pdither, sdithtype, sdither, sfilt, lfilt,
                          rpatt, grps, ints, short_pupil,
                          long_pupil, grismval, coordparallel,
                          i_obs, 1, template_name, instrument, obs_label)
            exposures_dictionary = self.add_exposure(exposures_dictionary, tup_to_add)
            # self.APTObservationParams = self.add_exposure(self.APTObservationParams, tup_to_add)
            self.obs_tuple_list.append(tup_to_add)

        # make sure all list items in the returned dictionary have the same length
        for key, item in exposures_dictionary.items():
            if len(item) == 0:
                exposures_dictionary[key] = [0] * len(exposures_dictionary['Instrument'])
        return exposures_dictionary

    def read_commissioning_template(self, template, template_name, obs, prop_params):
        # Get proposal parameters
        pi_name, prop_id, prop_title, prop_category, science_category, coordparallel, i_obs, obs_label = prop_params

        # dictionary that holds the content of this observation only
        exposures_dictionary = copy.deepcopy(self.empty_exposures_dictionary)

        # Set namespace
        ns = "{http://www.stsci.edu/JWST/APT/Template/WfscCommissioning}"

        # Set parameters that are constant for all WFSC obs
        #typeflag = template_name
        typeflag = 'imaging'
        grismval = 'N/A'
        subarr = 'FULL'
        pdithtype = 'NONE'
        pdither = '1'
        sdithtype = 'STANDARD'
        sdither = '1'

        # Find observation-specific parameters
        mod = template.find(ns + 'Module').text
        num_WFCgroups = int(template.find(ns + 'ExpectedWfcGroups').text)

        # Determine if there is an aperture override
        mod, subarr = self.check_for_aperture_override(obs, mod, subarr, i_obs)

        # Find filter parameters for all filter configurations within obs
        filter_configs = template.findall('.//' + ns + 'FilterConfig')

        for filt in filter_configs:
            sfilt = filt.find(ns + 'ShortFilter').text
            try:
                lfilt = filt.find(ns + 'LongFilter').text
            except AttributeError:
                lfilt = 'F480M'
            rpatt = filt.find(ns + 'ReadoutPattern').text
            grps = filt.find(ns + 'Groups').text
            ints = filt.find(ns + 'Integrations').text

            # Separate pupil and filter in case of filter that is
            # mounted in the pupil wheel
            if ' + ' in sfilt:
                split_ind = sfilt.find(' + ')
                short_pupil = sfilt[0:split_ind]
                sfilt = sfilt[split_ind + 1:]
            else:
                short_pupil = 'CLEAR'

            if ' + ' in lfilt:
                p = lfilt.find(' + ')
                long_pupil = lfilt[0:p]
                lfilt = lfilt[p + 1:]
            else:
                long_pupil = 'CLEAR'

            # Repeat for the number of expected WFSC groups + 1
            for j in range(num_WFCgroups + 1):
                # Add all parameters to dictionary
                tup_to_add = (pi_name, prop_id, prop_title, prop_category,
                              science_category, typeflag, mod, subarr, pdithtype,
                              pdither, sdithtype, sdither, sfilt, lfilt,
                              rpatt, grps, ints, short_pupil,
                              long_pupil, grismval, coordparallel,
                              i_obs , j + 1, template_name, 'NIRCAM', obs_label)

                exposures_dictionary = self.add_exposure(exposures_dictionary, tup_to_add)
                self.obs_tuple_list.append(tup_to_add)

        # make sure all list items in the returned dictionary have the same length
        for key, item in exposures_dictionary.items():
            if len(item) == 0:
                exposures_dictionary[key] = [0] * len(exposures_dictionary['Instrument'])
                # self.APTObservationParams = self.add_exposure(self.APTObservationParams, tup_to_add)
                # self.obs_tuple_list.append(tup_to_add)

        return exposures_dictionary, num_WFCgroups

    def read_globalalignment_template(self, template, template_name, obs, prop_params):
        # Get proposal parameters
        pi_name, prop_id, prop_title, prop_category, science_category, coordparallel, i_obs, obs_label = prop_params

        # dictionary that holds the content of this observation only
        exposures_dictionary = copy.deepcopy(self.empty_exposures_dictionary)

        ns = "{http://www.stsci.edu/JWST/APT/Template/WfscGlobalAlignment}"

        # Set parameters that are constant for all WFSC obs
        #typeflag = template_name
        typeflag = 'imaging'
        grismval = 'N/A'
        short_pupil = 'CLEAR'
        subarr = 'FULL'
        pdither = '1'
        pdithtype = 'NONE'
        sdithtype = 'STANDARD'
        sdither = '1'

        # Determine the Global Alignment Iteration Type
        GA_iteration = obs.find('.//' + ns + 'GaIteration').text

        if GA_iteration == 'ADJUST1':
            n_exp = 3
        elif GA_iteration == 'ADJUST2':
            n_exp = 6  # technically 5, but 3 is repeated?
        elif GA_iteration == 'BSCORRECT':
            # Technically has 2 dithers, but that doesn't seem to be incorporated...
            n_exp = 2
        elif GA_iteration == 'CORRECT+ADJUST':
            n_exp = 6  # technically 5, but 3 is repeated?
        elif GA_iteration == 'CORRECT':
            n_exp = 3

        # Find observation-specific parameters
        mod = template.find(ns + 'Module').text
        # num_WFCgroups = int(template.find(ns + 'ExpectedWfcGroups').text)

        # Determine if there is an aperture override
        mod, subarr = self.check_for_aperture_override(obs, mod, subarr, i_obs)

        # Find filter parameters for all filter configurations within obs
        ga_nircam_configs = template.findall('.//' + ns + 'NircamParameters')

        for conf in ga_nircam_configs:
            sfilt = conf.find(ns + 'ShortFilter').text
            try:
                lfilt = conf.find(ns + 'LongFilter').text
            except AttributeError:
                lfilt = 'F480M'
            rpatt = conf.find(ns + 'ReadoutPattern').text
            grps = conf.find(ns + 'Groups').text
            ints = conf.find(ns + 'Integrations').text

            # Separate pupil and filter in case of filter that is
            # mounted in the pupil wheel
            if ' + ' in sfilt:
                split_ind = sfilt.find(' + ')
                short_pupil = sfilt[0:split_ind]
                sfilt = sfilt[split_ind + 1:]
            else:
                short_pupil = 'CLEAR'

            if ' + ' in lfilt:
                p = lfilt.find(' + ')
                long_pupil = lfilt[0:p]
                lfilt = lfilt[p + 1:]
            else:
                long_pupil = 'CLEAR'

        # Repeat for the number of exposures + 1
        for j in range(n_exp + 1):
            # Add all parameters to dictionary
            tup_to_add = (pi_name, prop_id, prop_title, prop_category,
                          science_category, typeflag, mod, subarr, pdithtype,
                          pdither, sdithtype, sdither, sfilt, lfilt,
                          rpatt, grps, ints, short_pupil,
                          long_pupil, grismval, coordparallel,
                          i_obs, j + 1, template_name, 'NIRCAM', obs_label)

            exposures_dictionary = self.add_exposure(exposures_dictionary, tup_to_add)
            self.obs_tuple_list.append(tup_to_add)

        # make sure all list items in the returned dictionary have the same length
        for key, item in exposures_dictionary.items():
            if len(item) == 0:
                exposures_dictionary[key] = [0] * len(exposures_dictionary['Instrument'])

            # self.APTObservationParams = self.add_exposure(self.APTObservationParams, tup_to_add)
            # self.obs_tuple_list.append(tup_to_add)

        return exposures_dictionary, n_exp

    def read_coarsephasing_template(self, template, template_name, obs, prop_params):
        # Get proposal parameters
        pi_name, prop_id, prop_title, prop_category, science_category, coordparallel, i_obs, obs_label = prop_params

        # dictionary that holds the content of this observation only
        exposures_dictionary = copy.deepcopy(self.empty_exposures_dictionary)

        ns = "{http://www.stsci.edu/JWST/APT/Template/WfscCoarsePhasing}"

        # Set parameters that are constant for all WFSC obs
        #typeflag = template_name
        typeflag = 'imaging'
        grismval = 'N/A'
        pdither = '1'
        pdithtype = 'NONE'
        sdithtype = 'STANDARD'
        sdither = '1'

        # Find the module and derive the subarrays
        mod = template.find(ns + 'Module').text
        if mod == 'A':
            mods = ['SUB96DHSPILA'] + ['DHSPILA'] * 6
            subarrs = ['NRCA3_DHSPIL_SUB96'] + ['NRCA3_DHSPIL'] * 6
        if mod == 'B':
            mods = ['SUB96DHSPILB'] + ['DHSPILB'] * 6
            subarrs = ['NRCB4_DHSPIL_SUB96'] + ['NRCB4_DHSPIL'] * 6

        # Find the exposure parameters for the In Focus, DHS, and Defocus modes
        readouts = [r.text for r in obs.findall('.//' + ns + 'ReadoutPattern')]
        groups = [g.text for g in obs.findall('.//' + ns + 'Groups')]
        integrations = [i.text for i in obs.findall('.//' + ns + 'Integrations')]
        inds = [0, 1, 1, 1, 1, 2, 2]
        readouts = np.array(readouts)[inds]
        groups = np.array(groups)[inds]
        integrations = np.array(integrations)[inds]

        # List the pupils and filters in the appropriate order
        sw_pupils = ['CLEAR', 'GDHS0', 'GDHS0', 'GDHS60', 'GDHS60', 'WLP8', 'WLM8']
        sw_filts = ['F212N', 'F150W2', 'F150W2', 'F150W2', 'F150W2', 'F212N', 'F212N']
        lw_pupils = ['F405N'] * 7
        lw_filts = ['F444W'] * 7

        for i in range(7):
            mod = mods[i]
            subarr = subarrs[i]

            sfilt = sw_filts[i]
            lfilt = lw_filts[i]
            short_pupil = sw_pupils[i]
            long_pupil = lw_pupils[i]

            rpatt = readouts[i]
            grps = groups[i]
            ints = integrations[i]

            # Repeat for two dithers
            for j in range(2):
                # Add all parameters to dictionary
                tup_to_add = (pi_name, prop_id, prop_title, prop_category,
                              science_category, typeflag, mod, subarr, pdithtype,
                              pdither, sdithtype, sdither, sfilt, lfilt,
                              rpatt, grps, ints, short_pupil,
                              long_pupil, grismval, coordparallel,
                              i_obs, j + 1, template_name, 'NIRCAM', obs_label)

                exposures_dictionary = self.add_exposure(exposures_dictionary, tup_to_add)
                self.obs_tuple_list.append(tup_to_add)

        # make sure all list items in the returned dictionary have the same length
        for key, item in exposures_dictionary.items():
            if len(item) == 0:
                exposures_dictionary[key] = [0] * len(exposures_dictionary['Instrument'])

                    # self.APTObservationParams = self.add_exposure(self.APTObservationParams, tup_to_add)
                    # self.obs_tuple_list.append(tup_to_add)
        n_tiles_phasing = 14

        return exposures_dictionary, n_tiles_phasing

    def read_finephasing_template(self, template, template_name, obs, prop_params):
        # Get proposal parameters
        pi_name, prop_id, prop_title, prop_category, science_category, coordparallel, i_obs, obs_label = prop_params

        # dictionary that holds the content of this observation only
        exposures_dictionary = copy.deepcopy(self.empty_exposures_dictionary)

        ns = "{http://www.stsci.edu/JWST/APT/Template/WfscFinePhasing}"

        # Set parameters that are constant for all WFSC obs
        #typeflag = template_name
        typeflag = 'imaging'
        grismval = 'N/A'
        pdither = '1'
        pdithtype = 'NONE'
        sdithtype = 'STANDARD'
        sdither = '1'

        # Find the module and derive the subarrays
        mod = template.find(ns + 'Module').text
        if mod == 'A':
            mod = 'A3'
        elif mod == 'B':
            mod = 'B4'

        # Determine the sensing type, and list the pupils and filters
        # in the appropriate order
        sensing_type = obs.find('.//' + ns + 'SensingType').text

        n_configs = 0
        n_dithers = []
        subarrs = []
        mods = []
        sw_pupils = []
        sw_filts = []
        lw_pupils = []
        lw_filts = []
        readouts = []
        groups = []
        integrations = []

        if sensing_type in ['LOS Jitter', 'Both']:
            ########## IF WE WANT TO MODEL TARGET ACQ:
            # n_configs += 2
            # n_dithers += [1] * n_configs

            # subarrs += ['SUB64FP1' + mod, 'SUB8FP1' + mod]
            # mods += subarrs

            # sw_pupils += ['CLEAR', 'CLEAR']
            # sw_filts += ['F212N', 'F200W']
            # lw_pupils += ['F405N'] * n_configs # default?
            # lw_filts += ['F444W'] * n_configs # default?

            # # Find/define the exposure parameters for the target
            # # acquisition and LOS imaging modes
            # readouts += ['RAPID', 'RAPID']
            # acq_groups = obs.find('.//' + ns + 'AcqNumGroups').text
            # LOSimg_groups = obs.find('.//' + ns + 'LosImgNumGroups').text
            # groups += [acq_groups, LOSimg_groups]
            # LOSimg_ints = obs.find('.//' + ns + 'LosImgNumInts').text
            # integrations += [1, LOSimg_ints]

            n_configs += 1
            n_dithers += [1] * n_configs

            subarrs += ['SUB8FP1{}'.format(mod[0])]
            mods += [mod]

            sw_pupils += ['CLEAR']
            sw_filts += ['F200W']
            lw_pupils += ['F405N'] * n_configs # default?
            lw_filts += ['F444W'] * n_configs # default?

            # Find/define the exposure parameters for the target
            # acquisition and LOS imaging modes
            readouts += ['RAPID']
            groups += [obs.find('.//' + ns + 'LosImgNumGroups').text]
            integrations += [obs.find('.//' + ns + 'LosImgNumInts').text]

        if sensing_type in ['Fine Phasing', 'Both']:
            # Deterimine what diversity of sensing
            diversity = obs.find('.//' + ns + 'Diversity').text
            if diversity == 'ALL':
                n_configs_fp = 5
                n_configs += n_configs_fp
            elif diversity == 'ALL+187N':
                n_configs_fp = 7
                n_configs += n_configs_fp
            elif diversity == 'PM8':
                n_configs_fp = 2
                n_configs += n_configs_fp

            n_dithers += [2] * n_configs_fp

            subarrs += ['FP1'] * n_configs_fp
            mods += [mod] * n_configs_fp

            sw_pupils += ['WLM8', 'WLP8', 'WLP8', 'WLM8', 'CLEAR', 'WLM8', 'WLP8'][:n_configs_fp]
            sw_filts += ['F212N', 'F212N', 'WLP4', 'WLP4', 'WLP4', 'F187N', 'F187N'][:n_configs_fp]
            lw_pupils += ['F405N'] * n_configs_fp
            lw_filts += ['F444W'] * n_configs_fp

            # Find the exposure parameters for the +/- 8, + 12, and +/-4 modes
            readouts_fp = [r.text for r in obs.findall('.//' + ns + 'ReadoutPattern')]
            groups_fp = [g.text for g in obs.findall('.//' + ns + 'Groups')]
            integrations_fp = [i.text for i in obs.findall('.//' + ns + 'Integrations')]
            inds = [0, 0, 1, 2, 2, 3, 3][:n_configs_fp]
            readouts += list(np.array(readouts_fp)[inds])
            groups += list(np.array(groups_fp)[inds])
            integrations += list(np.array(integrations_fp)[inds])


        sensing = obs.find('.//' + self.apt + 'WavefrontSensing').text
        if sensing == 'SENSING_ONLY':
            n_repeats = 1
        else:
            n_repeats = 2

        for z in range(n_repeats):
            for i in range(n_configs):
                subarr = subarrs[i]
                mod = mods[i]

                sfilt = sw_filts[i]
                lfilt = lw_filts[i]
                short_pupil = sw_pupils[i]
                long_pupil = lw_pupils[i]

                rpatt = readouts[i]
                grps = groups[i]
                ints = integrations[i]

                n_dith = n_dithers[i]

                # Repeat for designated number of dithers
                for j in range(n_dith):
                    # Add all parameters to dictionary
                    tup_to_add = (pi_name, prop_id, prop_title, prop_category,
                                  science_category, typeflag, mod, subarr, pdithtype,
                                  pdither, sdithtype, sdither, sfilt, lfilt,
                                  rpatt, grps, ints, short_pupil,
                                  long_pupil, grismval, coordparallel,
                                  i_obs, j + 1, template_name, 'NIRCAM', obs_label)

                    exposures_dictionary = self.add_exposure(exposures_dictionary, tup_to_add)
                    self.obs_tuple_list.append(tup_to_add)

                # make sure all list items in the returned dictionary have the same length
        for key, item in exposures_dictionary.items():
            if len(item) == 0:
                exposures_dictionary[key] = [0] * len(exposures_dictionary['Instrument'])
                    # self.APTObservationParams = self.add_exposure(self.APTObservationParams, tup_to_add)
                    # self.obs_tuple_list.append(tup_to_add)

        n_tiles_phasing = sum(n_dithers) * n_repeats

        return exposures_dictionary, n_tiles_phasing

    def read_wfss_template(self, template, template_name, obs, prop_params):
        # Get proposal parameters
        pi_name, prop_id, prop_title, prop_category, science_category, coordparallel, i_obs, obs_label = prop_params

        # dictionary that holds the content of this observation only
        exposures_dictionary = copy.deepcopy(self.empty_exposures_dictionary)

        # Set namespace
        ns = "{http://www.stsci.edu/JWST/APT/Template/NircamWfss}"

        instrument = obs.find(self.apt + 'Instrument').text

        mod = template.find(ns + 'Module').text
        subarr = template.find(ns + 'Subarray').text
        grismval = template.find(ns + 'Grism').text
        if grismval == 'BOTH':
            grismval = ['GRISMR', 'GRISMC']
        else:
            grismval = [grismval]

        explist = template.find(ns + 'ExposureList')
        expseqs = explist.findall(ns + 'ExposureSequences')

        # Determine if there is an aperture override
        mod, subarr = self.check_for_aperture_override(obs, mod, subarr, i_obs)

        # if BOTH was specified for the grism,
        # then we need to repeat the sequence of
        # grism/direct/grism/direct/outoffield for each grism
        print("Is this in the wrong order? should we loop over exposuresequence and then grism?")
        for gnum in range(len(grismval)):
            for expseq in expseqs:
                # sequence = grism, direct, grism, direct, outoffield
                # if grism == both, sequence is done for grismr,
                # then repeated for grismc
                grism = grismval[gnum]
                # need to switch the order of the grism and direct
                # exposures in order for them to be chronological
                grismexp = expseq.find(ns + 'GrismExposure')
                typeflag = 'WFSS'
                sfilt = grismexp.find(ns + 'ShortFilter').text
                lfilt = grismexp.find(ns + 'LongFilter').text
                rpatt = grismexp.find(ns + 'ReadoutPattern').text
                groups = grismexp.find(ns + 'Groups').text
                integrations = grismexp.find(ns + 'Integrations').text

                pdithtype = template.find(ns + 'PrimaryDitherType').text
                try:
                    pdither = template.find(ns + 'PrimaryDithers').text
                except AttributeError:
                    pdither = 1
                sdither = template.find(ns + 'SubpixelPositions').text
                sdithtype = template.find(ns + 'SubpixelPositions').text

                # separate pupil and filter in case of filter
                # that is mounted in the pupil wheel
                if ' + ' in sfilt:
                    p = sfilt.find(' + ')
                    short_pupil = sfilt[0:p]
                    sfilt = sfilt[p + 1:]
                else:
                    short_pupil = 'CLEAR'

                long_pupil = grism
                tup_to_add = (pi_name, prop_id, prop_title, prop_category,
                              science_category, typeflag, mod, subarr,
                              pdithtype, pdither, sdithtype,
                              sdither, sfilt, lfilt, rpatt, groups,
                              integrations, short_pupil, long_pupil,
                              grism, coordparallel,
                              i_obs, 1, template_name, instrument, obs_label)

                exposures_dictionary = self.add_exposure(exposures_dictionary, tup_to_add)
                # self.APTObservationParams = self.add_exposure(self.APTObservationParams, tup_to_add)
                self.obs_tuple_list.append(tup_to_add)


                # TEST
                print('TEST add number_of_dithers for nircam wfss GRISM OBS')
                if 'Point' in sdither:
                    sdither_value = int(sdither[0])
                else:
                    sdither_value = int(sdither)
                print(pdither, sdither)
                exposures_dictionary['number_of_dithers'].append(sdither_value*int(pdither))

                # Check to see if the direct image will be collected
                # This will be either true or false
                direct_done = expseq.find(ns + 'DirectImage').text

                # We need to collect the information about the direct image even if it is
                # not being collected, because the same info will be used for the out of field
                # observations.
                directexp = expseq.find(ns + 'DiExposure')
                typeflag = 'imaging'
                pdither = '1'  # direct image has no dithers
                sdither = '1'  # direct image has no dithers
                sdithtype = '1'  # direct image has no dithers
                grism = 'N/A'
                sfilt = directexp.find(ns + 'ShortFilter').text
                lfilt = directexp.find(ns + 'LongFilter').text
                rpatt = directexp.find(ns + 'ReadoutPattern').text
                grps = directexp.find(ns + 'Groups').text
                ints = directexp.find(ns + 'Integrations').text

                # separate pupil and filter in case of filter
                # that is mounted in the pupil wheel
                if ' + ' in sfilt:
                    p = sfilt.find(' + ')
                    short_pupil = sfilt[0:p]
                    sfilt = sfilt[p + 1:]
                else:
                    short_pupil = 'CLEAR'

                if ' + ' in lfilt:
                    p = lfilt.find(' + ')
                    long_pupil = lfilt[0:p]
                    lfilt = lfilt[p + 1:]
                else:
                    long_pupil = 'CLEAR'

                direct_tup_to_add = (pi_name, prop_id, prop_title, prop_category,
                                     science_category, typeflag, mod, subarr, pdithtype,
                                     pdither, sdithtype, sdither, sfilt, lfilt,
                                     rpatt, grps, ints, short_pupil, long_pupil,
                                     grism, coordparallel,
                                     i_obs, 1, template_name, instrument, obs_label)

                # Only add the direct image if APT says it will be observed
                if direct_done == 'true':
                    self.obs_tuple_list.append(tup_to_add)
                    exposures_dictionary = self.add_exposure(exposures_dictionary, direct_tup_to_add)

                    # TEST
                    print('TEST add number_of_dithers for nircam wfss DIRECT OBS')
                    exposures_dictionary['number_of_dithers'].append(1)


            # Now we need to add the two out-of-field exposures, which are
            # not present in the APT file (but are in the associated pointing
            # file from APT. We can just
            # duplicate the entries for the direct images taken immediately
            # prior.
            exposures_dictionary = self.add_exposure(exposures_dictionary, direct_tup_to_add)
            exposures_dictionary = self.add_exposure(exposures_dictionary, direct_tup_to_add)
            self.obs_tuple_list.append(tup_to_add)
            self.obs_tuple_list.append(tup_to_add)

            # TEST
            print('TEST add number_of_dithers for nircam wfss OUT OF FIELDx2')
            exposures_dictionary['number_of_dithers'].extend([1, 1])

        # Make sure all entries are lists
        for key in exposures_dictionary.keys():
            if type(exposures_dictionary[key]) is not list:
                exposures_dictionary[key] = list(exposures_dictionary[key])

        # make sure all list items in the returned dictionary have the same length
        for key, item in exposures_dictionary.items():
            if len(item) == 0:
                exposures_dictionary[key] = [0] * len(exposures_dictionary['Instrument'])
        return exposures_dictionary

    def read_niriss_wfss_template(self, template, template_name, obs, proposal_param_dict, parallel=False):
        '''prime_template is needed if WFSS is parallel to a nircam imaging observation. In that case we need to look
        into the nircam observation to see if the niriss direct images are to be dithered '''
        # Get proposal parameters
        #pi_name, prop_id, prop_title, prop_category, science_category, coordparallel, i_obs, obs_label = proposal_param_dict
        #i_obs = proposal_param_dict['ObservationID']

        instrument = 'NIRISS'

        # Dummy module name for NIRISS. Needed for consistency in dictionary entry
        mod = 'N'
        subarr = 'FULL'
        long_filter = 'N/A'
        long_pupil = 'N/A'

        # dictionary that holds the content of this observation only
        exposures_dictionary = copy.deepcopy(self.empty_exposures_dictionary)

        # Set namespace
        ns = "{http://www.stsci.edu/JWST/APT/Template/NirissWfss}"

        if parallel:
            prime_template = obs.find(self.apt + 'Template')[0]
            prime_template_name = etree.QName(prime_template).localname
            prime_ns = "{{{}/Template/{}}}".format(self.apt.replace('{','').replace('}',''), prime_template_name)

            # boolean indicating which instrument is not prime but parallel
            parallel_instrument = True
            prime_instrument = obs.find(self.apt + 'Instrument').text
            print('IN READ_NIRISS_WFSS_TEMPLATE')
            print('Prime: {}   Parallel: {}'.format(prime_instrument, instrument))
            dither_direct = prime_template.find(prime_ns + 'DitherNirissWfssDirectImages').text
            sdither_type_grism = prime_template.find(prime_ns + 'CoordinatedParallelSubpixelPositions').text
            sdither_grism = sdither_type_grism[0]
        else:
            parallel_instrument = False
            prime_instrument = instrument
            dither_direct = 'NO_DITHERING'
            sdither_type_grism = '1'
            sdither_grism = '1'

        # Can be various types
        # WFSS stand-alone observation or WFSS as parallel, to NIRCam imaging: Integer: number
        # of primary dithers WFSS as prime, with NIRCam imaging as parallel: str
        # (e.g. '2-POINT-LARGE-NIRCam')
        dvalue = template.find(ns + 'PrimaryDithers').text
        try:
            dvalue_int = np.int(dvalue)
            pdither_grism = dvalue
        except ValueError:
            # When NIRISS is prime with NIRCam parallel, the PrimaryDithers field can be
            # (e.g. '2-POINT-LARGE-NIRCAM'), where the first character is always the number
            # of dither positions. Not sure how to save both this name as well as the DitherSize
            # value. I don't think there are header keywords for both, with PATTTYPE being the
            # only keyword for dither pattern names.
            pdither_grism = dvalue[0]
            # pdither_type = dvalue

        # Dither size can be SMALL, MEDIUM, LARGE
        pdither_type_grism = template.find(ns + 'DitherSize').text


        print('Primary dither info:', pdither_grism, pdither_type_grism)


        explist = template.find(ns + 'ExposureList')
        expseqs = explist.findall(ns + 'ExposureSequences')

        # Determine if there is an aperture override
        mod, subarr = self.check_for_aperture_override(obs, mod, subarr, proposal_param_dict['ObservationID'])

        delta_exp_dict_length = 0
        for expseq in expseqs:
            # Grism values are listed for each ExposureSequence
            grismval = expseq.find(ns + 'Sequence').text
            if grismval == 'BOTH':
                grismval = ['GRISM150R', 'GRISM150C']
                # Due to the way APT groups the direct images between the grism exposures
                # in cases where grism is set to 'BOTH' we need to adjust the exposure_dictionary
                # length to compare to when there are parallel observations
                delta_exp_dict_length -= 1
            else:
                grismval = [grismval]
            filter_name = expseq.find(ns + 'Filter').text

            # Loop over grism selections
            for grism_number, grism in enumerate(grismval):
                # sequence = direct pre, dithered grism, direct post
                # but if grism == both: direct pre1, dithered grism, direct pre2,
                # dithered grism2, direct post

                # Mini dictionary just for exposure sequence
                exp_seq_dict = {}

                # Collect info on the direct exposure
                directexp = expseq.find(ns + 'DiExposure')
                typeflag = 'imaging'
                if dither_direct == 'NO_DITHERING':
                    pdither = '1'  # direct image has no dithers
                    pdither_type = '1'  # direct image has no dithers
                    sdither = '1'
                    sdither_type = '1'
                else:
                    pdither = pdither_grism
                    pdither_type = pdither_type_grism
                    sdither = sdither_grism
                    sdither_type = sdither_type_grism

                print('for direct image, primary dither info: ', pdither, pdither_type)

                tile = '1'
                direct_grismvalue = 'N/A'
                pupil = 'CLEARP'  # NIRISS filter MUST be in filter wheel, not PUPIL wheel
                rpatt = directexp.find(ns + 'ReadoutPattern').text
                grps = directexp.find(ns + 'Groups').text
                ints = directexp.find(ns + 'Integrations').text

                #direct_tup_to_add = (pi_name, prop_id, prop_title, prop_category,
                #                     science_category, typeflag, mod, subarr, pdithtype,
                #                     pdither, sdithtype, sdither, filter_name, long_filter,
                #                     rpatt, grps, ints, pupil, long_pupil,
                #                     direct_grismvalue, coordparallel,
                #                     i_obs, 1, template_name, instrument, obs_label)

                # Add the 'pre' direct image
                print("Adding pre-direct image for grismval: {}".format(grism))
                #self.obs_tuple_list.append(direct_tup_to_add)
                #exposures_dictionary = self.add_exposure(exposures_dictionary, direct_tup_to_add)

                print('NIRISS predirect image, filter and pupil are: {} {} {} {}'.format(filter_name, long_filter, pupil, long_pupil))


                # Collect info on grism exposure
                grismexp = expseq.find(ns + 'GrismExposure')
                grism_typeflag = 'wfss'
                grism_pupil = grism
                grism_rpatt = grismexp.find(ns + 'ReadoutPattern').text
                grism_grps = grismexp.find(ns + 'Groups').text
                grism_ints = grismexp.find(ns + 'Integrations').text

                #grism_tup_to_add = (pi_name, prop_id, prop_title, prop_category,
                #                    science_category, grism_typeflag, mod, subarr,
                #                    pdithtype_grism, pdither_grism, sdithtype,
                #                    sdither, filter_name, long_filter, grism_rpatt, grism_grps,
                #                    grism_ints, grism_pupil, long_pupil, grism, coordparallel,
                #                    i_obs, 1, template_name, instrument, obs_label)

                # Add the grism exposure
                print("Adding grism image for grismval: {}".format(grism))
                #exposures_dictionary = self.add_exposure(exposures_dictionary, grism_tup_to_add)
                #self.obs_tuple_list.append(grism_tup_to_add)

                print('NIRISS grism image, filter and pupil are: {} {} {} {}'.format(filter_name, long_filter, grism_pupil, long_pupil))


                # Now add the entry for the 'post' direct image
                print("Adding post-direct image for grismval: {}".format(grism))
                #self.obs_tuple_list.append(direct_tup_to_add)
                #exposures_dictionary = self.add_exposure(exposures_dictionary, direct_tup_to_add)

                print('NIRISS postdirect image, filter and pupil are: {} {} {} {}'.format(filter_name, long_filter, pupil, long_pupil))

                # Update to add values via dictionary
                # Adding exposures this way (i.e. including separate entries for direct images before and after each grism)
                # makes a table with a length that agrees with the length of the pointing file. BUT it does not work in terms
                # of matching up the number of prime and parallel exposures. In that case, the 2 direct images that are after
                # the first grism and before the second are grouped into a single exposure with 2 dithers. It is the length
                # of this 'condensed' table that needs to match up with the length of the table from the other instrument.
                # Note sure how to handle this since all three tables are eventually combined. This is only a  problem with
                # NIRISS WFSS as prime because when it is parallel you are only allowed to choose one grism at at time.
                exp_seq_dict['Mode'] = [typeflag, grism_typeflag, typeflag]
                exp_seq_dict['Module'] = ['N'] * 3
                exp_seq_dict['Subarray'] = ['FULL'] * 3  # Niriss WFSS is always full frame
                exp_seq_dict['PrimaryDitherType'] = [pdither_type, pdither_type_grism, pdither_type]
                exp_seq_dict['PrimaryDithers'] = [pdither, pdither_grism, pdither]
                exp_seq_dict['SubpixelPositions'] = [sdither, sdither_grism, sdither]
                exp_seq_dict['SubpixelDitherType'] = [sdither_type, sdither_type_grism, sdither_type]
                exp_seq_dict['CoordinatedParallel'] = [parallel] * 3
                exp_seq_dict['Instrument'] = [instrument] * 3
                exp_seq_dict['ParallelInstrument'] = [parallel_instrument] * 3
                exp_seq_dict['ShortFilter'] = [filter_name] * 3
                exp_seq_dict['LongFilter'] = [long_filter] * 3
                exp_seq_dict['ReadoutPattern'] = [rpatt, grism_rpatt, rpatt]
                exp_seq_dict['Groups'] = [grps, grism_grps, grps]
                exp_seq_dict['Integrations'] = [ints, grism_ints, ints]
                exp_seq_dict['ShortPupil'] = [pupil, grism_pupil, pupil]
                exp_seq_dict['LongPupil'] = [long_pupil] * 3
                exp_seq_dict['Grism'] = [direct_grismvalue, grism, direct_grismvalue]
                exp_seq_dict['ObservationID'] = [proposal_param_dict['ObservationID']] * 3
                exp_seq_dict['TileNumber'] = [tile] * 3
                exp_seq_dict['APTTemplate'] = [template_name] * 3
                exp_seq_dict['ObservationName'] = [proposal_param_dict['ObservationName']] * 3
                exp_seq_dict['number_of_dithers'] = [1, int(pdither_grism)*int(sdither), 1]
                exp_seq_dict['FilterWheel'] = [filter_name] * 3
                exp_seq_dict['PupilWheel'] = [pupil, grism_pupil, pupil]

                #######################################################################
                # Update exposure dictionary to return
                print('below we need to make sure we are not double counting any keywords')
                for key in self.APTObservationParams_keys:
                    if key in exp_seq_dict.keys():
                        value = exp_seq_dict[key]
                    elif key in proposal_param_dict.keys():
                        value = [proposal_param_dict[key]] * 3
                    else:
                        value = [str(None)] * 3

                    if (key in ['PrimaryDithers', 'ImageDithers']) and ((value is None) or (value == 'None')):
                        value = ['1'] * 3

                    exposures_dictionary[key].extend(value)

                # add keys that were not defined in self.APTObservationParams_keys
                # (to be fixed in Class.__init__ later )
                for key in exp_seq_dict.keys():
                    if key not in self.APTObservationParams_keys:
                        # if key not yet present, create entry
                        if key not in exposures_dictionary.keys():
                            print('key {} not present in APTObservationParams nor exposures_dictionary'.format(key))
                            stop
                            exposures_dictionary[key] = [str(exposure_dict[key])]
                        else:
                            print('key {} not present in APTObservationParams'.format(key))
                            stop
                            exposures_dictionary[key].append(str(exposure_dict[key]))

        # Make sure all entries are lists
        for key in exposures_dictionary.keys():
            if type(exposures_dictionary[key]) is not list:
                exposures_dictionary[key] = list(exposures_dictionary[key])

        # Make sure all list items in the returned dictionary have the same length
        for key, item in exposures_dictionary.items():
            if len(item) == 0:
                exposures_dictionary[key] = [0] * len(exposures_dictionary['Instrument'])

        # Dictionary length to use when comparing to that from a parallel observation
        exp_len = len(exposures_dictionary['number_of_dithers']) + delta_exp_dict_length
        return exposures_dictionary, exp_len

    def check_for_aperture_override(self, obs, mod, subarr, i_obs):
        '''Determine if there is an aperture override'''

        override = obs.find('.//' + self.apt + 'FiducialPointOverride')
        if override is not None:
            mod = override.text
            if 'FULL' not in mod:
                subarray_list_file = os.path.join(PACKAGE_DIR, 'config',
                                                  'NIRCam_subarray_definitions.list')
                config = ascii.read(subarray_list_file)
                try:
                    i_sub = list(config['AperName']).index(mod)
                except ValueError:
                    i_sub = [mod in name for name in np.array(config['AperName'])]
                    i_sub = np.where(i_sub)[0]
                    if len(i_sub) > 1 or len(i_sub) == 0:
                        raise ValueError('Unable to match FiducialPointOverride {} to valid aperture in observation {}.'.format(mod, i_obs))

                subarr = config[i_sub]['Name']
                if type(subarr) != np.str_: # Don't know why, but astropy tables aren't behaving
                    subarr = subarr[0]

                print('Aperture override: subarray {}'.format(subarr))

            return mod, subarr

        else:
            return mod, subarr

    def read_parallel_exposures(self, obs, exposures_dictionary, proposal_parameter_dictionary, verbose=False,
                                exposures_dictionary_length=None):
        """Read the exposures of the parallel instrument.

        Parameters
        ----------
        obs : APT xml element
            Observation section of xml file
        exposures_dictionary : dict
            Exposures of the prime instrument
        proposal_parameter_dictionary : dict
            Parameters to extract
        verbose : bool
            Verbosity

        Returns
        -------
        parallel_exposures_dictionary : dict
            Parallel exposures.

        """
        # Determine what template is used for the parallel observation


        print('IN READ_PARALLEL_EXPOSURES')


        template = obs.find(self.apt + 'FirstCoordinatedTemplate')[0]
        template_name = etree.QName(template).localname
        print(template_name)
        if template_name in ['NircamImaging', 'NircamEngineeringImaging', 'NirissExternalCalibration',
                             'NirspecImaging', 'MiriMRS', 'FgsExternalCalibration']:
            parallel_exposures_dictionary = self.read_generic_imaging_template(template,
                                                                               template_name, obs,
                                                                               proposal_parameter_dictionary,
                                                                               parallel=True,
                                                                               verbose=verbose)
            parallel_length = len(parallel_exposures_dictionary['number_of_dithers'])

        elif template_name == 'NirissWfss':
            parallel_exposures_dictionary, parallel_length = self.read_niriss_wfss_template(template, template_name, obs,
                                                                                            proposal_parameter_dictionary,
                                                                                            parallel=True)
        else:
            raise ValueError('Parallel observation template {} not supported.'.format(template_name))

        # If an adjusted length of the exposures dictionary (from NIRISS WFSS) is not supplied,
        # then calculate it
        if exposures_dictionary_length is None:
            exposures_dictionary_length = len(exposures_dictionary['number_of_dithers'])

        if parallel_length != exposures_dictionary_length:
            print('Primary exposures: {}, Parallel exposures: {}'.format(len(exposures_dictionary['number_of_dithers']), parallel_length))
            raise RuntimeError('Mismatch in the number of parallel observations.')

        return parallel_exposures_dictionary
