# ! /usr/bin/env python


import copy
import os
import re
from collections import OrderedDict

from lxml import etree
from astropy.io import ascii
import numpy as np

APT_DIR = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
PACKAGE_DIR = os.path.dirname(APT_DIR)

flatten_list = lambda l: [item for sublist in l for item in sublist]
# def flatten_list(l):
#     if l == []:
#         return []
#     else:
#         fl = [item for sublist in l for item in sublist]
#         return fl[0]


def append_dictionary(base_dictionary, added_dictionary, braid=False):
    """Append the content of added_dictionary key-by-key to the base_dictionary.

    This assumes that the keys refer to lists.

    Parameters
    ----------
    base_dictionary : dict
    added_dictionary : dict
    braid : bool
        If true, the elements of added_dictionary are added in alternating sequence.
        This is used to synchronize parallel observations with the pointing file.

    Returns
    -------
    new_dictionary : dict
        Dictionary where every key holds a list of lists

    """
    new_dictionary = copy.deepcopy(base_dictionary)

    # extract an arbitrary key name
    first_key = [key for i, key in enumerate(base_dictionary.keys()) if i == 0][0]

    # Insert keys from added_dictionary that are not yet present in base_dictionary
    for key in added_dictionary.keys():
        if key not in base_dictionary.keys():
            new_dictionary[key] = ['None'] * len(base_dictionary[first_key])

    # # Insert keys from added_dictionary that are not yet present in base_dictionary
    # for key in added_dictionary.keys():
    #     if key not in base_dictionary.keys():
    #         new_dictionary[key] = ['None'] * len(base_dictionary[first_key])

    # Append the items
    for key in new_dictionary.keys():
        if key not in added_dictionary.keys():
            continue
        # print('{} {}'.format(key, new_dictionary[key]))
        if len(new_dictionary[key]) == 0:
            new_dictionary[key] = added_dictionary[key]
        else:
            if braid:
                # solution from https://stackoverflow.com/questions/3678869/pythonic-way-to-combine-two-lists-in-an-alternating-fashion
                new_dictionary[key] = [sub[i] for i in range(len(added_dictionary[key])) for sub in
                                       [new_dictionary[key], added_dictionary[key]]]
            else:
                new_dictionary[key] = new_dictionary[key] + added_dictionary[key]

    return new_dictionary


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
                          # 'Number', # the observation number  = ObservationID
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
        # obs_results = observation_data.findall('.//' + self.apt + 'Observation')
        observation_list = observation_data.findall('.//' + self.apt + 'Observation')
        # observation_numbers = [np.int(obs.find(self.apt + 'Number').text) for obs in observation_list]


        # Loop through observatiobs, get parameters
        # for i_obs, obs in zip(i_observations, observations):
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
                                   'NirissExternalCalibration',  # NIRISS
                                   'NirspecImaging', 'NirspecInternalLamp',  # NIRSpec
                                   'MiriMRS',  # MIRI
                                   'FgsExternalCalibration',  # FGS
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
                # number_of_exposures = len(flatten_list(self.APTObservationParams['Instrument']))
                number_of_entries = len(self.APTObservationParams['Instrument'])
                print('APTObservationParams Dictionary holds {} entries before reading template'.format(number_of_entries))
                # if i_obs==2:
                #     1/0

            # Get coordinated parallel
            coordparallel = obs.find(self.apt + 'CoordinatedParallel').text
            CoordinatedParallelSet = None
            if coordparallel == 'true':
                CoordinatedParallelSet = obs.find(self.apt + 'CoordinatedParallelSet').text

            try:
                obs_label = obs.find(self.apt + 'Label').text
            except AttributeError:
                # label tag not present
                obs_label = 'Observation 1'

            # DitherPatternType = obs.find(self.apt + 'DitherPatternType').text
            # ImageDithers = obs.find(self.apt + 'ImageDithers').text

            # extract visit numbers
            visit_numbers = [np.int(element.items()[0][1]) for element in obs if
                             element.tag.split(self.apt)[1] == 'Visit']



            # Determine pointing offset?
            offset = obs.find('.//' + self.apt + 'Offset')
            try:
                offset_x = offset.get('Xvalue')
                offset_y = offset.get('Yvalue')
            except AttributeError:
                offset_x, offset_y = 0, 0
            if (offset_x != 0) or (offset_y != 0):
                print('* * * OFFSET OF ({}, {}) IN OBS {} NOT APPLIED ***'.format(offset_x,
                                                                                  offset_y,
                                                                                  observation_number))

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

            # 1/0
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # If template is NircamImaging or NircamEngineeringImaging
            # if template_name in ['ircamEngineeringImaging']:
            # # if template_name in ['NircamImaging', 'NircamEngineeringImaging']:
            #     exposures_dictionary = self.read_imaging_template(template, template_name, obs, prop_params)
            #     # exposures_dictionary = copy.deepcopy(self.APTObservationParams)
            #     if coordparallel == 'true':
            #
            #         parallel_exposures_dictionary = self.read_parallel_exposures(obs,
            #                                                                 exposures_dictionary,
            #                                                                 proposal_parameter_dictionary,
            #                                                                 verbose=verbose)
            #         exposures_dictionary = self.append_dictionary(exposures_dictionary,
            #                                                       parallel_exposures_dictionary,
            #                                                       braid=True)
            #     self.APTObservationParams = self.append_dictionary(self.APTObservationParams,
            #                                                            exposures_dictionary)

            if template_name in ['NircamImaging', 'NircamEngineeringImaging', 'NirissExternalCalibration', 'NirspecImaging', 'MiriMRS', 'FgsExternalCalibration']:
                exposures_dictionary = self.read_generic_imaging_template(template, template_name, obs, proposal_parameter_dictionary, verbose=verbose)
                if coordparallel == 'true':
                    parallel_exposures_dictionary = self.read_parallel_exposures(obs, exposures_dictionary, proposal_parameter_dictionary, verbose=verbose)
                    exposures_dictionary = append_dictionary(exposures_dictionary, parallel_exposures_dictionary, braid=True)

                # self.APTObservationParams = append_dictionary(self.APTObservationParams, exposures_dictionary)

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


            if verbose:
                print('Dictionary read from template has {} entries.'.format(len(exposures_dictionary['Instrument'])))
                print(exposures_dictionary['number_of_dithers'])
                print(exposures_dictionary['Instrument'])

            # set default number of dithers to one, for downstream processing
            for i, n_dither in enumerate(exposures_dictionary['number_of_dithers']):
                if int(n_dither) == 0:
                    exposures_dictionary['number_of_dithers'][i] = '1'


            # add the exposure dictionary to the main dictionary
            # if template_name in ['WfscCommissioning', 'WfscGlobalAlignment', 'WfscCoarsePhasing', 'WfscFinePhasing', 'NircamWfss']:
            self.APTObservationParams = append_dictionary(self.APTObservationParams,
                                                          exposures_dictionary)

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


            # Now we need to look for mosaic details, if any
            mosaic_tiles = obs.findall('.//' + self.apt + 'MosaicTiles')

            # count only tiles that are included
            tile_state = np.array([mosaic_tiles[i].find('.//' + self.apt + 'TileState').text for i in range(len(mosaic_tiles))])
            n_tiles = np.sum(tile_state=='Tile Included')

            # If WFSC, look at expected groups rather than mosaic tiles:
            # if n_tiles == 0 and template_name == 'WfscCommissioning':
            #     if num_WFCgroups:
            #         n_tiles = num_WFCgroups + 1
            # if n_tiles == 0 and template_name == 'WfscGlobalAlignment':
            #     n_tiles = n_exp + 1
            # if n_tiles == 0 and 'Phasing' in template_name:
            #     n_tiles = n_tiles_phasing

            label = obs_label

            print("Found {} tile(s) for observation {} {}".format(n_tiles, observation_number, label))
            print('Found {} visits with numbers: {}'.format(len(visit_numbers), visit_numbers))


            if n_tiles > 1:
                for i in range(n_tiles - 1):
                    self.APTObservationParams = append_dictionary(self.APTObservationParams, exposures_dictionary)

            self.observation_info[observation_number] = {}
            self.observation_info[observation_number]['visit_numbers'] = visit_numbers

            # all_visit_numbers = [self.observation_info[key]['visit_numbers'] for key in self.observation_info.keys()]
            # flatten list of lists
            # all_visit_numbers = [item for sublist in all_visit_numbers for item in sublist]

            # print(self.APTObservationParams['number_of_dithers'])
            # 1/0
            # number_of_exposures = np.sum([np.int(s) for s in self.APTObservationParams['number_of_dithers']])
            # number_of_exposures = np.sum([np.int(s) for s in self.APTObservationParams['PrimaryDithers']])
            # assert number_of_exposures == len(all_visit_numbers)
            # 1/0

            if verbose:
                number_of_entries_after = len(self.APTObservationParams['Instrument'])
                print('APTObservationParams Dictionary holds {} entries after reading template ({:+d} entries)'.format(number_of_entries_after, number_of_entries_after-number_of_entries))
                # assert number_of_entries_after == len(self.APTObservationParams['ReadoutPatternLong'])
                # print(self.APTObservationParams['Instrument'])

        if verbose:
            print('Finished reading APT xml file.')
            print('+'*100)

        # for key, item in self.APTObservationParams.items():
        #     if len(item) == 0:
        #         self.APTObservationParams[key] = [0] * len(self.APTObservationParams['Instrument'])

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
            prime_instrument = obs.find(self.apt + 'Instrument').text
            print('Prime: {}   Parallel: {}'.format(prime_instrument, instrument))
        else:
            instrument = obs.find(self.apt + 'Instrument').text
            parallel_instrument = False
            prime_instrument = instrument

        # verbose = True

        # Get proposal parameters
        # pi_name, prop_id, prop_title, prop_category, science_category, coordparallel, i_obs = prop_params
        exposures_dictionary = copy.deepcopy(self.empty_exposures_dictionary)
        # exposures_dictionary = OrderedDict()
        # for key in self.APTObservationParams_keys:
        #     exposures_dictionary[key] = []

        # Set namespace
        # if parallel:
        #     ns = "{{http://www.stsci.edu/JWST/APT/FirstCoordinatedTemplate/{}}}".format(template_name)
        # else:
        # ns = "{{http://www.stsci.edu/JWST/APT/Template/{}}}".format(template_name)
        ns = "{{{}/Template/{}}}".format(self.apt.replace('{','').replace('}',''), template_name)

        DitherPatternType = None

        if prime_instrument in ['NIRCAM', 'FGS']:
            dither_key_name = 'PrimaryDithers'
        elif prime_instrument in ['NIRISS', 'MIRI', 'NIRSPEC']:
            dither_key_name = 'ImageDithers'
        # number of dithers defaults to 1
        number_of_dithers = 1
        # ImageDithers = 1


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

            # 1/0

            # Set parameters that are constant for all imaging obs
            # typeflag = template_name
            # typeflag = 'imaging'
            # grismval = 'N/A'
            # # short_pupil = 'CLEAR'
            #
            # # Find observation-specific parameters
            # mod = template.find(ns + 'Module').text
            # subarr = template.find(ns + 'Subarray').text
            # pdithtype = template.find(ns + 'PrimaryDitherType').text

            # Determine if there is an aperture override
            # mod, subarr = self.check_for_aperture_override(obs, mod, subarr, proposal_parameter_dictionary['ObservationID'])
            override = obs.find('.//' + self.apt + 'FiducialPointOverride')
            FiducialPointOverride = False
            if override is not None:
                FiducialPointOverride = True

            observation_dict['FiducialPointOverride'] = str(FiducialPointOverride)

            if dither_key_name in observation_dict.keys():
                number_of_dithers = observation_dict[dither_key_name]
            else:
                # try:
                #     number_of_dithers = template.find(ns + dither_key_name).text
                # except:
                print('Element {} not found, use default value.'.format(dither_key_name))

            # sdithtype = template.find(ns + 'SubpixelDitherType').text
            #
            # try:
            #     sdither = template.find(ns + 'SubpixelPositions').text
            # except:
            #     try:
            #         stemp = template.find(ns + 'CoordinatedParallelSubpixelPositions').text
            #         sdither = np.int(stemp[0])
            #     except:
            #         sdither = '1'

            # Find filter parameters for all filter configurations within obs
            filter_configs = template.findall('.//' + ns + 'FilterConfig')
            # 1/0
            verbose = True
            # loop over filter configurations
            for filter_config_index, filter_config in enumerate(filter_configs):
                filter_config_dict = {}
                print('Filter config index {}'.format(filter_config_index))
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
                    else:
                        value = str(None)

                    if (key == 'Mode'):
                        value = 'imaging'

                    exposures_dictionary[key].append(value)


            # print(exposures_dictionary)
            # 1/0

        else:
            # verbose = True
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
                    # number of dithers
                    # ImageDithers = int(DitherPatternType[0])
                    number_of_dithers = int(DitherPatternType[0])


                # Different SI conventions of how to list exposure parameters
                if ((instrument.lower()=='niriss') and (element_tag_stripped == 'ExposureList')) | \
                        ((instrument.lower() == 'fgs') and (element_tag_stripped == 'Exposures'))| \
                        ((instrument.lower() == 'miri') and (element_tag_stripped == 'ExposureList'))| \
                        ((instrument.lower() == 'nirspec') and (element_tag_stripped == 'Exposures')):
                    for exposure in element.findall(ns + 'Exposure'):
                        exposure_dict = {}
                        exposure_dict['DitherPatternType'] = DitherPatternType
                        # exposure_dict['ImageDithers'] = np.int(ImageDithers)
                        if number_of_dithers is None:
                            number_of_dithers = 1
                        exposure_dict[dither_key_name] = np.int(number_of_dithers)
                        for exposure_parameter in exposure:
                            parameter_tag_stripped = exposure_parameter.tag.split(ns)[1]
                            if verbose:
                                print('{} {}'.format(parameter_tag_stripped, exposure_parameter.text))
                            exposure_dict[parameter_tag_stripped] = exposure_parameter.text

                        exposure_dict['number_of_dithers'] = exposure_dict[dither_key_name]
                        # 1/0
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
                            # elif key == 'number_of_dithers':
                            #     value = number_of_dithers
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

            # self.APTObservationParams = self.add_exposure(self.APTObservationParams, tup_to_add)
            # self.APTObservationParams = exposures_dictionary

        for key in exposures_dictionary.keys():
            if type(exposures_dictionary[key]) is not list:
                exposures_dictionary[key] = list(exposures_dictionary[key])

        # make sure all list items in the returned dictionary have the same length
        for key, item in exposures_dictionary.items():
            if len(item) == 0:
                exposures_dictionary[key] = [0] * len(exposures_dictionary['Instrument'])


        return exposures_dictionary
        # print(exposures_dictionary)
        # return



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
        short_pupil = 'CLEAR'
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
        # pdithtype = template.find(ns + 'PrimaryDitherType').text
        # pdither = template.find(ns + 'PrimaryDithers').text
        # sdither = template.find(ns + 'SubpixelPositions').text
        # sdithtype = template.find(ns + 'SubpixelPositions').text
        explist = template.find(ns + 'ExposureList')
        expseqs = explist.findall(ns + 'ExposureSequences')

        # Determine if there is an aperture override
        mod, subarr = self.check_for_aperture_override(obs, mod, subarr, i_obs)

        # if BOTH was specified for the grism,
        # then we need to repeat the sequence of
        # grism/direct/grism/direct/outoffield for each grism
        for gnum in range(len(grismval)):
            for expseq in expseqs:
                # sequence = grism, direct, grism, direct, outoffield
                # if grism == both, sequence is done for grismr,
                # then repeated for grismc
                grismvalue = grismval[gnum]
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

                long_pupil = grismvalue
                tup_to_add = (pi_name, prop_id, prop_title, prop_category,
                              science_category, typeflag, mod, subarr,
                              pdithtype, pdither, sdithtype,
                              sdither, sfilt, lfilt, rpatt, groups,
                              integrations, short_pupil, long_pupil,
                              grismvalue, coordparallel,
                              i_obs, 1, template_name, instrument, obs_label)

                exposures_dictionary = self.add_exposure(exposures_dictionary, tup_to_add)
                # self.APTObservationParams = self.add_exposure(self.APTObservationParams, tup_to_add)
                self.obs_tuple_list.append(tup_to_add)

                directexp = expseq.find(ns + 'DiExposure')
                #typeflag = template_name
                typeflag = 'imaging'
                pdither = '1'  # direct image has no dithers
                sdither = '1'  # direct image has no dithers
                sdithtype = '1'  # direct image has no dithers
                grismvalue = 'N/A'
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
                                     grismvalue, coordparallel,
                                     i_obs, 1, template_name, instrument, obs_label)
                # self.APTObservationParams = self.add_exposure(self.APTObservationParams, direct_tup_to_add)
                self.obs_tuple_list.append(tup_to_add)
                exposures_dictionary = self.add_exposure(exposures_dictionary, direct_tup_to_add)


            # Now we need to add the two out-of-field exposures, which are
            # not present in the APT file (but are in the associated pointing
            # file from APT. We can just
            # duplicate the entries for the direct images taken immediately
            # prior. BUT, will there ever be a case where there is no preceding
            # direct image?
            exposures_dictionary = self.add_exposure(exposures_dictionary, direct_tup_to_add)
            exposures_dictionary = self.add_exposure(exposures_dictionary, direct_tup_to_add)
            # self.APTObservationParams = self.add_exposure(self.APTObservationParams, direct_tup_to_add)
            # self.APTObservationParams = self.add_exposure(self.APTObservationParams, direct_tup_to_add)
            self.obs_tuple_list.append(tup_to_add)
            self.obs_tuple_list.append(tup_to_add)

        # make sure all list items in the returned dictionary have the same length
        for key, item in exposures_dictionary.items():
            if len(item) == 0:
                exposures_dictionary[key] = [0] * len(exposures_dictionary['Instrument'])
        return exposures_dictionary

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

    def read_parallel_exposures(self, obs, exposures_dictionary, proposal_parameter_dictionary, verbose=False):
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
        template = obs.find(self.apt + 'FirstCoordinatedTemplate')[0]
        template_name = etree.QName(template).localname
        parallel_exposures_dictionary = self.read_generic_imaging_template(template,
                                                                           template_name, obs,
                                                                           proposal_parameter_dictionary,
                                                                           parallel=True,
                                                                           verbose=verbose)

        # if exposures_dictionary['Instrument'][0] in ['NIRISS']:
        #     dither_key_name = 'ImageDithers'
        # elif exposures_dictionary['Instrument'][0] in ['NIRSPEC', 'MIRI']:
        #     dither_key_name = 'ImageDithers'
        # elif exposures_dictionary['Instrument'][0] in ['NIRCAM']:
        #     dither_key_name = 'PrimaryDithers'

        # print('{} == {}'.format(exposures_dictionary['number_of_dithers'], parallel_exposures_dictionary['number_of_dithers']))

        # if exposures_dictionary['Instrument'][0] not in ['NIRCAM']:
            # check that replicating the dithers for the parallel instrument worked
            # if len(parallel_exposures_dictionary[dither_key_name]) != len(
            #         exposures_dictionary[dither_key_name]):
        if len(parallel_exposures_dictionary['number_of_dithers']) != len(
                exposures_dictionary['number_of_dithers']):
            raise RuntimeError('Mismatch in the number of parallel observations.')

        return parallel_exposures_dictionary