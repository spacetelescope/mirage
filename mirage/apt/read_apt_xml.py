# ! /usr/bin/env python


import copy
import os
import re
from collections import OrderedDict

from lxml import etree
from astropy.io import ascii
from astropy.table import Table, Column
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
                          'DitherPatternType', 'ImageDithers',  # NIRISS
                          'number_of_dithers',  # uniform name across instruments
                          'FiducialPointOverride', 'TargetID', 'TargetRA', 'TargetDec'
                          ]
        FilterParams_keys = ['ShortFilter', 'LongFilter', 'ShortPupil', 'LongPupil',
                             'ReadoutPattern', 'Groups', 'Integrations',
                             'FilterWheel', 'PupilWheel',  # for NIRISS
                             'NumOutputs'
                             ]
        OtherParams_keys = ['Mode', 'Grism',
                            'IntegrationsShort', 'GroupsShort', 'Dither',  # MIRI
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
        """Main function. Read in the .xml file from APT, and output dictionary of parameters.

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

        # Get target names - - - - - - - - - - - - - - - - - - - - - - - - - -
        targs = tree.find(self.apt + 'Targets')
        target_elements = targs.findall(self.apt + 'Target')
        self.target_info = {}
        for target in target_elements:
            t_name = target.find(self.apt + 'TargetName').text
            t_coords = target.find(self.apt + 'EquatorialCoordinates').items()[0][1]
            ra_hour, ra_min, ra_sec, dec_deg, dec_arcmin, dec_arcsec = t_coords.split(' ')
            ra = '{}:{}:{}'.format(ra_hour, ra_min, ra_sec)
            dec = '{}:{}:{}'.format(dec_deg, dec_arcmin, dec_arcsec)
            self.target_info[t_name] = (ra, dec)
        print('target_info:')
        print(self.target_info)

        # Get parameters for each observation  - - - - - - - - - - - - - - - -

        # Find all observations (but use only those that use NIRCam or are WFSC)
        observation_data = tree.find(self.apt + 'DataRequests')
        observation_list = observation_data.findall('.//' + self.apt + 'Observation')

        # Loop through observations, get parameters
        for i_obs, obs in enumerate(observation_list):
            observation_number = obs.find(self.apt + 'Number').text.zfill(3)

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
                                   'NircamGrismTimeSeries', 'NircamTimeSeries',
                                   'NirissExternalCalibration', 'NirissWfss', 'NirissAmi',  # NIRISS
                                   'NirspecImaging', 'NirspecInternalLamp',  # NIRSpec
                                   'MiriMRS',  # MIRI
                                   'FgsExternalCalibration',  # FGS
                                   'NircamDark', 'NirissDark'  # Darks
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

            # Get coordinated parallel
            coordparallel = obs.find(self.apt + 'CoordinatedParallel').text

            if verbose:
                print('+'*100)
                print('Observation `{}` labelled `{}` uses template `{}`'.format(observation_number, label,
                                                                                 template_name))
                number_of_entries = len(self.APTObservationParams['Instrument'])
                print('APTObservationParams Dictionary holds {} entries before reading template'
                      .format(number_of_entries))
                if coordparallel == 'true':
                    print('Coordinated parallel observation')

            CoordinatedParallelSet = None
            if coordparallel == 'true':
                try:
                    CoordinatedParallelSet = obs.find(self.apt + 'CoordinatedParallelSet').text
                except AttributeError:
                    raise RuntimeError('Program does not specify parallels correctly.')

            try:
                obs_label = obs.find(self.apt + 'Label').text
            except AttributeError:
                # label tag not present
                obs_label = 'Observation 1'

            # Get target name
            try:
                targ_name = obs.find(self.apt + 'TargetID').text.split(' ')[1]
            except IndexError as e:
                print("No target ID for observation: {}".format(obs))
                targ_name = obs.find(self.apt + 'TargetID').text.split(' ')[0]

            # For NIRSpec Internal Lamp
            if targ_name == 'NONE':
                self.target_info[targ_name] = ('0', '0')

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
                                             'TargetID': targ_name,
                                             'TargetRA': self.target_info[targ_name][0],
                                             'TargetDec': self.target_info[targ_name][1]
                                             }

            if template_name in ['NircamImaging', 'NircamEngineeringImaging', 'NirissExternalCalibration',
                                 'NirspecImaging', 'MiriMRS', 'FgsExternalCalibration', 'NirissAmi']:
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
                exposures_dictionary, num_WFCgroups = self.read_commissioning_template(template,
                                                                                       template_name, obs,
                                                                                       prop_params)

            # If template is WFSC Global Alignment
            elif template_name in ['WfscGlobalAlignment']:
                exposures_dictionary, n_exp = self.read_globalalignment_template(template, template_name, obs,
                                                                                 prop_params)

            # If template is WFSC Coarse Phasing
            elif template_name in ['WfscCoarsePhasing']:
                exposures_dictionary, n_tiles_phasing = self.read_coarsephasing_template(template,
                                                                                         template_name, obs,
                                                                                         prop_params)

            # If template is WFSC Fine Phasing
            elif template_name in ['WfscFinePhasing']:
                exposures_dictionary, n_tiles_phasing = self.read_finephasing_template(template,
                                                                                       template_name, obs,
                                                                                       prop_params)

            # If template is NIRCam Grism Time Series
            elif template_name == 'NircamGrismTimeSeries':
                exposures_dictionary = self.read_nircam_grism_time_series(template, template_name, obs,
                                                                          proposal_parameter_dictionary)

            # If template is NIRCam Imaging Time Series
            elif template_name == 'NircamTimeSeries':
                exposures_dictionary = self.read_nircam_imaging_time_series(template, template_name, obs,
                                                                            proposal_parameter_dictionary)

            # If template is WFSS
            elif template_name == 'NircamWfss':
                exposures_dictionary = self.read_nircam_wfss_template(template, template_name, obs,
                                                                      proposal_parameter_dictionary)

            elif template_name == 'NirissWfss':
                exposures_dictionary = self.read_niriss_wfss_template(template, template_name, obs,
                                                                      proposal_parameter_dictionary,
                                                                      verbose=verbose)
                if coordparallel == 'true':
                    parallel_template_name = etree.QName(obs.find(self.apt + 'FirstCoordinatedTemplate')[0]).localname
                    if parallel_template_name in ['MiriImaging']:
                        pass
                    elif parallel_template_name in ['NircamImaging']:
                        parallel_exposures_dictionary = self.read_parallel_exposures(obs, exposures_dictionary,
                                                                                     proposal_parameter_dictionary,
                                                                                     verbose=verbose)
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

            # # set default number of dithers, for downstream processing
            # for i, n_dither in enumerate(exposures_dictionary['number_of_dithers']):
            #     if (template_name == 'NircamEngineeringImaging') and (n_dither == '2PLUS'):
            #         exposures_dictionary['number_of_dithers'][i] = '2'
            #     elif int(n_dither) == 0:
            #         exposures_dictionary['number_of_dithers'][i] = '1'

            # add the exposure dictionary to the main dictionary
            self.APTObservationParams = append_dictionary(self.APTObservationParams,
                                                          exposures_dictionary)
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            # Now we need to look for mosaic details, if any
            mosaic_tiles = obs.findall('.//' + self.apt + 'MosaicTiles')

            # count only tiles that are included
            tile_state = np.array([mosaic_tiles[i].find('.//' + self.apt + 'TileState').text for i in range(len(mosaic_tiles))])
            n_tiles = np.sum(np.array([True if state == 'Tile Included' else False for state in tile_state]))

            label = obs_label

            if verbose:
                print("Found {} tile(s) for observation {} {}".format(n_tiles, observation_number, label))
                if len(visit_numbers) > 0:
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
        # Temporary for creating truth tables to use in tests
        #bool_cols = ['ParallelInstrument']
        #int_cols = ['Groups', 'Integrations']
        #final_table = Table()
        #for key in self.APTObservationParams.keys():
        #    if key in bool_cols:
        #        data_type = bool
        #    elif key in int_cols:
        #        data_type = str
        #    else:
        #        data_type = str
        #    new_col = Column(data=self.APTObservationParams[key], name=key, dtype=data_type)
        #    final_table.add_column(new_col)
        #tmpoutdir = '/Users/hilbert/python_repos/mirage/tests/test_data'
        #filebase = os.path.split(infile)[1]
        #tmpoutfile = '{}{}'.format(filebase.split('.xml')[0], '.txt')
        #ascii.write(final_table, os.path.join(tmpoutdir, tmpoutfile), overwrite=True)

        return self.APTObservationParams

    def add_exposure(self, dictionary, tup):
        """Add an exposure to the exposure dictionary

        Parameters
        ----------
        dictionary : dict
            Information on individual exposures

        tup : tuple
            A tuple contianing information to add to dictionary as
            the next exposure

        Returns
        -------
        dictionary : dict
            With new exposure added
        """
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

    def read_generic_imaging_template(self, template, template_name, obs, proposal_parameter_dictionary,
                                      verbose=False, parallel=False):
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
            if verbose:
                print('Prime: {}   Parallel: {}'.format(prime_instrument, instrument))
            prime_template = obs.find(self.apt + 'Template')[0]
            prime_template_name = etree.QName(prime_template).localname
            prime_ns = "{{{}/Template/{}}}".format(self.apt.replace('{', '').replace('}', ''), prime_template_name)
            if verbose:
                print('PRIME TEMPLATE NAME IS: {}'.format(prime_template_name))
        else:
            instrument = obs.find(self.apt + 'Instrument').text
            parallel_instrument = False
            prime_instrument = instrument
            prime_template = template
            prime_template_name = template_name

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
        number_of_subpixel_positions = 1

        number_of_primary_dithers = 1
        number_of_subpixel_dithers = 1

        if instrument.lower() == 'nircam':
            # NIRCam uses FilterConfig structure to specifiy exposure parameters

            # store Module, Subarray, ... fields
            observation_dict = {}
            for field in template:
                key = field.tag.split(ns)[1]
                value = field.text
                observation_dict[key] = value

            # Determine if there is an aperture override
            override = obs.find('.//' + self.apt + 'FiducialPointOverride')
            FiducialPointOverride = True if override is not None else False

            # Get the number of primary and subpixel dithers
            primary_dithers_present = dither_key_name in observation_dict.keys()
            if primary_dithers_present:
                number_of_primary_dithers = observation_dict[dither_key_name]

                if (template_name == 'NircamEngineeringImaging') and (number_of_primary_dithers == '2PLUS'):
                    # Handle the special case for 2PLUS
                    number_of_primary_dithers = 2

                if observation_dict[dither_key_name] in ['2TIGHTGAPS']:
                    number_of_primary_dithers = observation_dict[dither_key_name][0]

                # Special case for 8NIRSPEC dither pattern
                if number_of_primary_dithers == '8NIRSPEC':
                    number_of_primary_dithers = '8'

            else:
                print('Primary dither element {} not found, use default primary dithers value (1).'.format(dither_key_name))

            # Find the number of subpixel dithers
            if not parallel:
                if observation_dict['SubpixelDitherType'] in ['3-POINT-WITH-MIRI-F770W']:
                    # Handle the special case for MIRI
                    number_of_subpixel_dithers = 3
                elif "-WITH-NIRISS" in observation_dict['SubpixelDitherType']:
                    number_of_subpixel_dithers = np.int(observation_dict['SubpixelDitherType'][0])
                elif observation_dict['SubpixelDitherType'] in ['STANDARD', 'IMAGING']:
                    number_of_subpixel_dithers = np.int(observation_dict['SubpixelPositions'])
            else:
                # For parallel instrument we ignore any dither info and set values to 0
                number_of_primary_dithers = 0
                number_of_subpixel_dithers = 0

                # The exception is if NIRISS WFSS is prime, in which case the dither pattern
                # (from NIRISS) can be imposed on only the second of each trio of exposures
                # (the exposure that is matched up to the NIRISS grism exposure), OR to all
                # three of each trio of exposures (meaning those matched up to the direct,
                # grism, and direct NIRISS exposures).

                # This doesn't actually matter at the moment since parallel instrument dither values
                # in the table are ignored.
                if prime_instrument == 'NIRISS' and prime_template_name == 'NirissWfss':
                    observation_dict['PrimaryDithers'] = prime_template.find(prime_ns + dither_key_name).text

                    # In the case where PrimaryDithers is e.g. 2-POINT-WITH-NIRCam,
                    # extract the '2' and place it in the PrimaryDithers field
                    try:
                        int_dithers = np.int(observation_dict['PrimaryDithers'])
                    except ValueError:
                        observation_dict['PrimaryDitherType'] = copy.deepcopy(observation_dict['PrimaryDithers'])
                        observation_dict['PrimaryDithers'] = observation_dict['PrimaryDithers'][0]
                    observation_dict['DitherSize'] = prime_template.find(prime_ns + 'DitherSize').text
                    number_of_primary_dithers = np.int(observation_dict['PrimaryDithers'][0])
                    number_of_subpixel_dithers = 1

            # Combine primary and subpixel dithers
            number_of_dithers = str(np.int(number_of_primary_dithers) * number_of_subpixel_dithers)
            print('Number of dithers: {} primary * {} subpixel = {}'.format(number_of_primary_dithers,
                                                                            number_of_subpixel_dithers,
                                                                            number_of_dithers))

            # Find filter parameters for all filter configurations within obs
            filter_configs = template.findall('.//' + ns + 'FilterConfig')
            # loop over filter configurations
            for filter_config_index, filter_config in enumerate(filter_configs):
                filter_config_dict = {}
                # print('Filter config index {}'.format(filter_config_index))
                for element in filter_config:
                    key = element.tag.split(ns)[1]
                    value = element.text
                    if key == 'ShortFilter':
                        ShortPupil, ShortFilter = self.separate_pupil_and_filter(value)
                        filter_config_dict['ShortPupil'] = ShortPupil
                    elif key == 'LongFilter':
                        LongPupil, LongFilter = self.separate_pupil_and_filter(value)
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

            ##########################################################
            # If NIRCam is prime with NIRISS WFSS parallel, then a DITHER_DIRECT
            # field will be added to the xml, describing whether the direct images
            # on either side of the grism exposure should be ditered. In a rare case
            # of the prime instrument having to conform to what the parallel instrument
            # is doing, this means that the NIRCam exposures taken at the same time as
            # the NIRISS direct images will also be dithered or not to match NIRISS. In
            # this special mode, NIRISS direct images are taken before and after each
            # grism exposure. Therefore for the NIRCam prime exposures, you can collect
            # them into groups of 3, and the dither_direct value will affect the first
            # and third exposures in each group. If dither_direct is false then all dithering,
            # primary and subpixel, are skipped. If dither_direct is true, then primary
            # and subpixel dithers are both done. It is guaranteed (by APT) in this case that
            # then number of exposures is a multiple of 3.
            try:
                dither_direct = observation_dict['DitherNirissWfssDirectImages']
                if dither_direct == 'NO_DITHERING':
                    if verbose:
                        print(('NIRISS WFSS parallel and NO_DITHERING set for direct imgages. Adjusting '
                               'number_of_dithers to 1 for the matching NIRCam exposures.'))
                    num_dithers = exposures_dictionary['number_of_dithers']
                    for counter in range(0, len(num_dithers), 3):
                        num_dithers[counter: counter+3] = ['1', num_dithers[counter+1], '1']
            except:
                pass
            ############################################################

        else:

            # Determine if there is an aperture override
            override = obs.find('.//' + self.apt + 'FiducialPointOverride')
            FiducialPointOverride = True if override is not None else False

            for element in template:
                element_tag_stripped = element.tag.split(ns)[1]

                # loop through exposures and collect dither parameters
                if element_tag_stripped == 'DitherPatternType':
                    DitherPatternType = element.text
                elif element_tag_stripped == 'ImageDithers':
                    number_of_primary_dithers = int(element.text)
                elif element_tag_stripped == 'SubpixelPositions':
                    if element.text != 'NONE':
                        number_of_subpixel_positions = np.int(element.text)
                elif element_tag_stripped == 'PrimaryDithers':
                    if (element.text is not None) & (element.text != 'NONE'):
                        number_of_primary_dithers = int(element.text)
                elif element_tag_stripped == 'Dithers':
                    DitherPatternType = element.find(ns + 'MrsDitherSpecification').find(ns + 'DitherType').text
                    number_of_primary_dithers = int(DitherPatternType[0])
                elif element_tag_stripped == 'SubpixelDithers':
                    if element.text is not None:
                        number_of_subpixel_dithers = int(element.text)


                # handle the NIRISS AMI case
                if number_of_subpixel_positions > number_of_subpixel_dithers:
                    number_of_subpixel_dithers = np.copy(number_of_subpixel_positions)

                # Determine if there is an aperture override
                override = obs.find('.//' + self.apt + 'FiducialPointOverride')
                FiducialPointOverride = True if override is not None else False

                # To reduce confusion, if this is the parallel instrument,
                # set the number of dithers to zero, since the prime
                # instrument controls the number of dithers
                if parallel:
                    number_of_primary_dithers = 0
                    number_of_subpixel_dithers = 0

                # Combine primary and subpixel dithers
                number_of_dithers = str(number_of_primary_dithers * number_of_subpixel_dithers)


                # Different SI conventions of how to list exposure parameters
                if ((instrument.lower() == 'niriss') and (element_tag_stripped == 'ExposureList')) | \
                        ((instrument.lower() == 'fgs') and (element_tag_stripped == 'Exposures'))| \
                        ((instrument.lower() == 'miri') and (element_tag_stripped == 'ExposureList'))| \
                        ((instrument.lower() == 'nirspec') and (element_tag_stripped == 'Exposures')):
                    for exposure in element.findall(ns + 'Exposure'):
                        exposure_dict = {}

                        # Load dither information into dictionary
                        exposure_dict['DitherPatternType'] = DitherPatternType

                        if (number_of_dithers is None) | (number_of_dithers == 'NONE'):
                            number_of_dithers = 1 * number_of_subpixel_positions

                        exposure_dict[dither_key_name] = np.int(number_of_dithers)
                        exposure_dict['number_of_dithers'] = exposure_dict[dither_key_name]

                        for exposure_parameter in exposure:
                            parameter_tag_stripped = exposure_parameter.tag.split(ns)[1]
                            # if verbose:
                            #     print('{} {}'.format(parameter_tag_stripped, exposure_parameter.text))
                            exposure_dict[parameter_tag_stripped] = exposure_parameter.text

                        # fill dictionary to return
                        for key in self.APTObservationParams_keys:
                            if key in exposure_dict.keys():
                                value = exposure_dict[key]
                            elif key in proposal_parameter_dictionary.keys():
                                value = proposal_parameter_dictionary[key]
                            elif key == 'Instrument':
                                value = instrument
                            elif key == 'ParallelInstrument':
                                value = parallel_instrument
                            elif key == 'FiducialPointOverride':
                                value = str(FiducialPointOverride)
                            elif key == 'APTTemplate':
                                value = template_name
                            else:
                                value = str(None)

                            if (key in ['PrimaryDithers', 'ImageDithers']) and (str(value) == 'None'):
                                value = '1'

                            if (key == 'Mode'):
                                if template_name not in ['NirissAmi']:
                                    value = 'imaging'
                                else:
                                    value = 'ami'

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

            if not parallel:
                print('Number of dithers: {} primary * {} subpixel = {}'.format(number_of_primary_dithers,
                                                                                number_of_subpixel_dithers,
                                                                                number_of_dithers))

        for key in exposures_dictionary.keys():
            if type(exposures_dictionary[key]) is not list:
                exposures_dictionary[key] = list(exposures_dictionary[key])

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

            # Add the number of dithers
            number_of_dithers = int(pdither) * int(sdither)
            exposures_dictionary['number_of_dithers'] = [str(number_of_dithers)] * len(exposures_dictionary['Instrument'])

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

        # Add the number of dithers
        number_of_dithers = int(pdither) * int(sdither)
        exposures_dictionary['number_of_dithers'] = [str(number_of_dithers)] * len(
            exposures_dictionary['Instrument'])

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
        mods = [mod] * 12
        if mod == 'A':
            subarrs = ['SUB96DHSPILA'] + ['FULL'] * 6
        if mod == 'B':
            subarrs = ['SUB96DHSPILB'] + ['FULL'] * 6

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

        # Add the number of dithers
        number_of_dithers = int(pdither) * int(sdither)
        exposures_dictionary['number_of_dithers'] = [str(number_of_dithers)] * len(
            exposures_dictionary['Instrument'])

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

            subarrs += ['SUB8FP1{}'.format(mod)]
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

            subarrs += ['FULL'] * n_configs_fp
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

                # Add all parameters to dictionary
                tup_to_add = (pi_name, prop_id, prop_title, prop_category,
                              science_category, typeflag, mod, subarr, pdithtype,
                              pdither, sdithtype, sdither, sfilt, lfilt,
                              rpatt, grps, ints, short_pupil,
                              long_pupil, grismval, coordparallel,
                              i_obs, 1, template_name, 'NIRCAM', obs_label)

                exposures_dictionary = self.add_exposure(exposures_dictionary, tup_to_add)
                exposures_dictionary['number_of_dithers'] += str(n_dith)

                self.obs_tuple_list.append(tup_to_add)

        # make sure all list items in the returned dictionary have the same length
        for key, item in exposures_dictionary.items():
            if len(item) == 0:
                exposures_dictionary[key] = [0] * len(exposures_dictionary['Instrument'])
                    # self.APTObservationParams = self.add_exposure(self.APTObservationParams, tup_to_add)
                    # self.obs_tuple_list.append(tup_to_add)

        n_tiles_phasing = sum(n_dithers) * n_repeats

        return exposures_dictionary, n_tiles_phasing

    def read_nircam_grism_time_series(self, template, template_name, obs, proposal_parameter_dictionary):
        """Parse a NIRCam Grism Time Series observation template from an APT xml file.
        Produce an exposure dictionary that lists all exposures (excluding dithers)
        from the template

        Parameters
        ----------
        template : lxml.etree._Element
            Template section from APT xml

        template_name : str
            The type of template (e.g. 'NirissWfss')

        obs : lxml.etree._Element
            Observation section from APT xml

        proposal_param_dict : dict
            Dictionary of proposal level information from the xml file
            (e.g. PI, Science Category, etc)

        Returns
        -------
        exposures_dictionary : dict
            Dictionary containing details on all exposures contained within
            the template. These details include things like filter, pupil,
            readout pattern, subarray, etc. Specifically for Grism Time Series,
            there will be entries for the TA exposure and the Time Series
            exposure.
        """
        # Dictionary that holds the content of this observation only
        exposures_dictionary = copy.deepcopy(self.empty_exposures_dictionary)

        # Set namespace
        ns = "{http://www.stsci.edu/JWST/APT/Template/NircamGrismTimeSeries}"

        # Observation-wide info
        instrument = obs.find(self.apt + 'Instrument').text

        # Get target name
        try:
            targ_name = obs.find(self.apt + 'TargetID').text.split(' ')[1]
        except IndexError as e:
            print("No target ID for observation: {}".format(obs))
            targ_name = obs.find(self.apt + 'TargetID').text.split(' ')[0]

        # Mode specific info, including target acq
        acq_target = template.find(ns + 'AcqTargetID').text
        if acq_target == 'Same Target as Observation':
            acq_target = targ_name

        acq_readout_pattern = template.find(ns + 'AcqReadoutPattern').text
        acq_groups = template.find(ns + 'AcqGroups').text
        acq_filter = 'F335M'
        acq_subarray = 'SUB32TATSGRISM'
        acq_integrations = '1'
        module = template.find(ns + 'Module').text
        subarray = template.find(ns + 'Subarray').text
        readout_pattern = template.find(ns + 'ReadoutPattern').text
        groups = template.find(ns + 'Groups').text
        integrations = template.find(ns + 'Integrations').text
        num_exps = template.find(ns + 'NumExps').text
        num_outputs = template.find(ns + 'NumOutputs').text  # Number of amplifiers used

        # Neither TA exposures nor Grism Time Series exposures allow dithering
        number_of_dithers = '1'
        number_of_subpixel_positions = '1'
        number_of_primary_dithers = '1'
        number_of_subpixel_dithers = '1'

        short_pupil_filter = template.find(ns + 'ShortPupilFilter').text
        long_pupil_filter = template.find(ns + 'LongPupilFilter').text
        short_pupil, short_filter = self.separate_pupil_and_filter(short_pupil_filter)
        long_pupil, long_filter = self.separate_pupil_and_filter(long_pupil_filter)

        # Populate observation dictionary with TA and Time Series exposures
        exposures_dictionary['Instrument'] = [instrument] * 2
        exposures_dictionary['Module'] = [module] * 2
        exposures_dictionary['TargetID'] = [acq_target, targ_name]
        exposures_dictionary['APTTemplate'] = [template_name] * 2
        exposures_dictionary['Mode'] = ['imaging', 'ts_grism']
        exposures_dictionary['Subarray'] = [acq_subarray, subarray]
        exposures_dictionary['ReadoutPattern'] = [acq_readout_pattern, readout_pattern]
        exposures_dictionary['Groups'] = [acq_groups, groups]
        exposures_dictionary['Integrations'] = [acq_integrations, integrations]
        exposures_dictionary['Exposures'] = ['1', num_exps]
        exposures_dictionary['NumOutputs'] = ['1', num_outputs]
        exposures_dictionary['PrimaryDithers'] = [number_of_primary_dithers, number_of_primary_dithers]
        exposures_dictionary['SubpixelPositions'] = [number_of_subpixel_dithers, number_of_subpixel_dithers]
        exposures_dictionary['ImageDithers'] = [number_of_dithers, number_of_dithers]
        exposures_dictionary['number_of_dithers'] = [number_of_dithers, number_of_dithers]
        exposures_dictionary['ShortFilter'] = [str(None), short_filter]
        exposures_dictionary['ShortPupil'] = [str(None), short_pupil]
        exposures_dictionary['LongFilter'] = [acq_filter, long_filter]
        exposures_dictionary['LongPupil'] = ['CLEAR', long_pupil]
        exposures_dictionary['FiducialPointOverride'] = [str(False)] * 2
        exposures_dictionary['ParallelInstrument'] = [False] * 2

        # Populate other keywords with None
        for key in self.APTObservationParams_keys:
            value = 'reset_value'
            if key in proposal_parameter_dictionary.keys() and key != 'TargetID':
                value = [proposal_parameter_dictionary[key]] * 2
            elif exposures_dictionary[key] == []:
                value = [str(None)] * 2
            if value != 'reset_value':
                exposures_dictionary[key].extend(value)

        return exposures_dictionary

    def read_nircam_imaging_time_series(self, template, template_name, obs, proposal_parameter_dictionary):
        """Parse a NIRCam Imaging Time Series observation template from an APT xml file.
        Produce an exposure dictionary that lists all exposures (excluding dithers)
        from the template

        Parameters
        ----------
        template : lxml.etree._Element
            Template section from APT xml

        template_name : str
            The type of template (e.g. 'NirissWfss')

        obs : lxml.etree._Element
            Observation section from APT xml

        proposal_param_dict : dict
            Dictionary of proposal level information from the xml file
            (e.g. PI, Science Category, etc)

        Returns
        -------
        exposures_dictionary : dict
            Dictionary containing details on all exposures contained within
            the template. These details include things like filter, pupil,
            readout pattern, subarray, etc. Specifically for Grism Time Series,
            there will be entries for the TA exposure and the Time Series
            exposure.
        """
        # Dictionary that holds the content of this observation only
        exposures_dictionary = copy.deepcopy(self.empty_exposures_dictionary)

        # Set namespace
        ns = "{http://www.stsci.edu/JWST/APT/Template/NircamTimeSeries}"

        # Observation-wide info
        instrument = obs.find(self.apt + 'Instrument').text

        # Get target name
        try:
            targ_name = obs.find(self.apt + 'TargetID').text.split(' ')[1]
        except IndexError as e:
            print("No target ID for observation: {}".format(obs))
            targ_name = obs.find(self.apt + 'TargetID').text.split(' ')[0]

        # Mode specific info, including target acq
        acq_target = template.find(ns + 'AcqTargetID').text
        if acq_target == 'Same Target as Observation':
            acq_target = targ_name

        acq_readout_pattern = template.find(ns + 'AcqReadoutPattern').text
        acq_groups = template.find(ns + 'AcqGroups').text
        acq_filter = 'F335M'
        acq_subarray = 'SUB32TATS'
        acq_integrations = '1'
        module = template.find(ns + 'Module').text
        subarray = template.find(ns + 'Subarray').text
        readout_pattern = template.find(ns + 'ReadoutPattern').text
        groups = template.find(ns + 'Groups').text
        integrations = template.find(ns + 'Integrations').text
        num_exps = template.find(ns + 'NumExps').text

        if subarray == 'FULL':
            num_outputs = 4
        else:
            num_outputs = 1

        # Neither TA exposures nor Grism Time Series exposures allow dithering
        number_of_dithers = '1'
        number_of_subpixel_positions = '1'
        number_of_primary_dithers = '1'
        number_of_subpixel_dithers = '1'

        short_pupil = template.find(ns + 'ShortPupil').text
        short_filter = template.find(ns + 'ShortFilter').text
        long_pupil = template.find(ns + 'LongPupil').text
        long_filter = template.find(ns + 'LongFilter').text

        # Populate observation dictionary with TA and Time Series exposures
        exposures_dictionary['Instrument'] = [instrument] * 2
        exposures_dictionary['Module'] = [module] * 2
        exposures_dictionary['TargetID'] = [acq_target, targ_name]
        exposures_dictionary['APTTemplate'] = [template_name] * 2
        exposures_dictionary['Mode'] = ['imaging', 'ts_imaging']
        exposures_dictionary['Subarray'] = [acq_subarray, subarray]
        exposures_dictionary['ReadoutPattern'] = [acq_readout_pattern, readout_pattern]
        exposures_dictionary['Groups'] = [acq_groups, groups]
        exposures_dictionary['Integrations'] = [acq_integrations, integrations]
        exposures_dictionary['Exposures'] = ['1', num_exps]
        exposures_dictionary['NumOutputs'] = ['1', num_outputs]
        exposures_dictionary['PrimaryDithers'] = [number_of_primary_dithers, number_of_primary_dithers]
        exposures_dictionary['SubpixelPositions'] = [number_of_subpixel_dithers, number_of_subpixel_dithers]
        exposures_dictionary['ImageDithers'] = [number_of_dithers, number_of_dithers]
        exposures_dictionary['number_of_dithers'] = [number_of_dithers, number_of_dithers]
        exposures_dictionary['ShortFilter'] = [str(None), short_filter]
        exposures_dictionary['ShortPupil'] = [str(None), short_pupil]
        exposures_dictionary['LongFilter'] = [acq_filter, long_filter]
        exposures_dictionary['LongPupil'] = ['CLEAR', long_pupil]
        exposures_dictionary['FiducialPointOverride'] = [str(False)] * 2
        exposures_dictionary['ParallelInstrument'] = [False] * 2

        # Populate other keywords with None
        for key in self.APTObservationParams_keys:
            value = 'reset_value'
            if key in proposal_parameter_dictionary.keys() and key != 'TargetID':
                value = [proposal_parameter_dictionary[key]] * 2
            elif exposures_dictionary[key] == []:
                value = [str(None)] * 2
            if value != 'reset_value':
                exposures_dictionary[key].extend(value)

        return exposures_dictionary

    def read_nircam_wfss_template(self, template, template_name, obs, proposal_param_dict):
        """Parse a NIRCam WFSS observation template from an APT xml file. Produce an exposure dictionary
        that lists all exposures (excluding dithers) from the template.

        Parameters
        ----------
        template : lxml.etree._Element
            Template section from APT xml

        template_name : str
            The type of template (e.g. 'NirissWfss')

        obs : lxml.etree._Element
            Observation section from APT xml

        proposal_param_dict : dict
            Dictionary of proposal level information from the xml file
            (e.g. PI, Science Category, etc)

        Returns
        -------
        exposures_dictionary : dict
            Dictionary containing details on all exposures contained within the template. These details
            include things like filter, pupil, readout pattern, subarray, etc
        """
        # Dictionary that holds the content of this observation only
        exposures_dictionary = copy.deepcopy(self.empty_exposures_dictionary)

        # Set namespace
        ns = "{http://www.stsci.edu/JWST/APT/Template/NircamWfss}"

        # Observation-wide info
        instrument = obs.find(self.apt + 'Instrument').text
        module = template.find(ns + 'Module').text
        subarr = template.find(ns + 'Subarray').text

        # Determine if there is an aperture override
        override = obs.find('.//' + self.apt + 'FiducialPointOverride')
        FiducialPointOverride = True if override is not None else False

        # Get primary and subpixel dither values for the grism exposures
        primary_dither_type_grism = template.find(ns + 'PrimaryDitherType').text
        if primary_dither_type_grism.lower() != 'none':
            primary_dither_grism = template.find(ns + 'PrimaryDithers').text
        else:
            primary_dither_grism = '1'
        subpix_dither_type_grism = template.find(ns + 'SubpixelPositions').text
        if subpix_dither_type_grism.lower() != 'none':
            subpix_dither_grism = subpix_dither_type_grism[0]
        else:
            subpix_dither_grism = '1'
        grism_number_of_dithers = str(int(primary_dither_grism) * int(subpix_dither_grism))

        # Direct and out of field images are never dithered
        primary_dither_direct = '1'
        primary_dither_type_direct = 'None'
        subpix_dither_direct = '1'
        subpix_dither_type_direct = 'None'
        direct_number_of_dithers = '1'

        # Which grism(s) are to be used
        grismval = template.find(ns + 'Grism').text
        if grismval == 'BOTH':
            grismval = ['GRISMR', 'GRISMC']
        else:
            grismval = [grismval]

        explist = template.find(ns + 'ExposureList')
        expseqs = explist.findall(ns + 'ExposureSequences')

        # if BOTH was specified for the grism,
        # then we need to repeat the sequence of
        # grism/direct/grism/direct/outoffield for each grism
        for grism_number, grism in enumerate(grismval):
            out_of_field_dict = {}
            for expseq in expseqs:
                # sequence = grism, direct, grism, direct, outoffield, outoffield
                # if grism == both, entire sequence is done for grismr,
                # then repeated for grismc

                # Mini dictionary just for exposure sequence
                exp_seq_dict = {}

                # Switch the order of the grism and direct
                # exposures from that within xml file in order for them to be chronological
                grismexp = expseq.find(ns + 'GrismExposure')
                grism_typeflag = 'wfss'
                grism_short_filter = grismexp.find(ns + 'ShortFilter').text
                grism_long_filter = grismexp.find(ns + 'LongFilter').text
                grism_readpatt = grismexp.find(ns + 'ReadoutPattern').text
                grism_groups = grismexp.find(ns + 'Groups').text
                grism_integrations = grismexp.find(ns + 'Integrations').text

                # Separate SW pupil and filter in case of filter
                # that is mounted in the pupil wheel
                grism_short_pupil, grism_short_filter = self.separate_pupil_and_filter(grism_short_filter)
                grism_long_pupil = grism

                # Check to see if the direct image will be collected
                # This will be either 'true' or 'false'
                direct_done = expseq.find(ns + 'DirectImage').text

                # We need to collect the information about the direct image even if it is
                # not being observed, because the same info will be used for the out of field
                # observations later.
                directexp = expseq.find(ns + 'DiExposure')
                direct_typeflag = 'imaging'
                direct_grism = 'N/A'
                direct_short_filter = directexp.find(ns + 'ShortFilter').text
                direct_long_filter = directexp.find(ns + 'LongFilter').text
                direct_readpatt = directexp.find(ns + 'ReadoutPattern').text
                direct_groups = directexp.find(ns + 'Groups').text
                direct_integrations = directexp.find(ns + 'Integrations').text

                # Separate pupil and filter in case of filter
                # that is mounted in the pupil wheel
                direct_short_pupil, direct_short_filter = self.separate_pupil_and_filter(direct_short_filter)
                direct_long_pupil, direct_long_filter = self.separate_pupil_and_filter(direct_long_filter)

                # Only add the direct image if APT says it will be observed
                if direct_done == 'true':
                    exp_seq_dict['Mode'] = [grism_typeflag, direct_typeflag]
                    exp_seq_dict['Module'] = [module] * 2
                    exp_seq_dict['Subarray'] = [subarr] * 2
                    exp_seq_dict['PrimaryDitherType'] = [primary_dither_type_grism, primary_dither_type_direct]
                    exp_seq_dict['PrimaryDithers'] = [primary_dither_grism, primary_dither_direct]
                    exp_seq_dict['SubpixelPositions'] = [subpix_dither_grism, subpix_dither_direct]
                    exp_seq_dict['SubpixelDitherType'] = [subpix_dither_type_grism, subpix_dither_type_direct]
                    exp_seq_dict['CoordinatedParallel'] = ['false'] * 2
                    exp_seq_dict['Instrument'] = [instrument] * 2
                    exp_seq_dict['ParallelInstrument'] = [False] * 2
                    exp_seq_dict['ShortFilter'] = [grism_short_filter, direct_short_filter]
                    exp_seq_dict['LongFilter'] = [grism_long_filter, direct_long_filter]
                    exp_seq_dict['ReadoutPattern'] = [grism_readpatt, direct_readpatt]
                    exp_seq_dict['Groups'] = [grism_groups, direct_groups]
                    exp_seq_dict['Integrations'] = [grism_integrations, direct_integrations]
                    exp_seq_dict['ShortPupil'] = [grism_short_pupil, direct_short_pupil]
                    exp_seq_dict['LongPupil'] = [grism_long_pupil, direct_long_pupil]
                    exp_seq_dict['Grism'] = [grism, direct_grism]
                    exp_seq_dict['ObservationID'] = [proposal_param_dict['ObservationID']] * 2
                    exp_seq_dict['TileNumber'] = ['1'] * 2
                    exp_seq_dict['APTTemplate'] = [template_name] * 2
                    exp_seq_dict['ObservationName'] = [proposal_param_dict['ObservationName']] * 2
                    exp_seq_dict['number_of_dithers'] = [grism_number_of_dithers, direct_number_of_dithers]
                    exp_seq_dict['FilterWheel'] = ['none'] * 2  # used for NIRISS
                    exp_seq_dict['PupilWheel'] = ['none'] * 2  # used for NIRISS
                    exp_seq_dict['FiducialPointOverride'] = [FiducialPointOverride] * 2
                else:
                    exp_seq_dict['Mode'] = [grism_typeflag]
                    exp_seq_dict['Module'] = [module]
                    exp_seq_dict['Subarray'] = [subarr]
                    exp_seq_dict['PrimaryDitherType'] = [primary_dither_type_grism]
                    exp_seq_dict['PrimaryDithers'] = [primary_dither_grism]
                    exp_seq_dict['SubpixelPositions'] = [subpix_dither_grism]
                    exp_seq_dict['SubpixelDitherType'] = [subpix_dither_type_grism]
                    exp_seq_dict['CoordinatedParallel'] = ['false']
                    exp_seq_dict['Instrument'] = [instrument]
                    exp_seq_dict['ParallelInstrument'] = [False]
                    exp_seq_dict['ShortFilter'] = [grism_short_filter]
                    exp_seq_dict['LongFilter'] = [grism_long_filter]
                    exp_seq_dict['ReadoutPattern'] = [grism_readpatt]
                    exp_seq_dict['Groups'] = [grism_groups]
                    exp_seq_dict['Integrations'] = [grism_integrations]
                    exp_seq_dict['ShortPupil'] = [grism_short_pupil]
                    exp_seq_dict['LongPupil'] = [grism_long_pupil]
                    exp_seq_dict['Grism'] = [grism]
                    exp_seq_dict['ObservationID'] = [proposal_param_dict['ObservationID']]
                    exp_seq_dict['TileNumber'] = ['1']
                    exp_seq_dict['APTTemplate'] = [template_name]
                    exp_seq_dict['ObservationName'] = [proposal_param_dict['ObservationName']]
                    exp_seq_dict['number_of_dithers'] = [grism_number_of_dithers]
                    exp_seq_dict['FilterWheel'] = ['none']  # used for NIRISS
                    exp_seq_dict['PupilWheel'] = ['none']  # used for NIRISS
                    exp_seq_dict['FiducialPointOverride'] = [FiducialPointOverride]

                # Add exp_seq_dict to the exposures_dictionary
                exposures_dictionary = self.append_to_exposures_dictionary(exposures_dictionary,
                                                                           exp_seq_dict,
                                                                           proposal_param_dict)

            # Now we need to add the two out-of-field exposures, which are
            # not present in the APT file (but are in the associated pointing
            # file from APT.) We can just duplicate the entries for the direct
            # images taken immediately prior. Out of field exposures are collected
            # for each grism in the observation
            out_of_field_dict['Mode'] = [direct_typeflag] * 2
            out_of_field_dict['Module'] = [module] * 2
            out_of_field_dict['Subarray'] = [subarr] * 2
            out_of_field_dict['PrimaryDitherType'] = [primary_dither_type_direct] * 2
            out_of_field_dict['PrimaryDithers'] = [primary_dither_direct] * 2
            out_of_field_dict['SubpixelPositions'] = [subpix_dither_direct] * 2
            out_of_field_dict['SubpixelDitherType'] = [subpix_dither_type_direct] * 2
            out_of_field_dict['CoordinatedParallel'] = ['false'] * 2
            out_of_field_dict['Instrument'] = [instrument] * 2
            out_of_field_dict['ParallelInstrument'] = [False] * 2
            out_of_field_dict['ShortFilter'] = [direct_short_filter] * 2
            out_of_field_dict['LongFilter'] = [direct_long_filter] * 2
            out_of_field_dict['ReadoutPattern'] = [direct_readpatt] * 2
            out_of_field_dict['Groups'] = [direct_groups] * 2
            out_of_field_dict['Integrations'] = [direct_integrations] * 2
            out_of_field_dict['ShortPupil'] = [direct_short_pupil] * 2
            out_of_field_dict['LongPupil'] = [direct_long_pupil] * 2
            out_of_field_dict['Grism'] = [direct_grism] * 2
            out_of_field_dict['ObservationID'] = [proposal_param_dict['ObservationID']] * 2
            out_of_field_dict['TileNumber'] = ['1'] * 2
            out_of_field_dict['APTTemplate'] = [template_name] * 2
            out_of_field_dict['ObservationName'] = [proposal_param_dict['ObservationName']] * 2
            out_of_field_dict['number_of_dithers'] = [direct_number_of_dithers] * 2
            out_of_field_dict['FilterWheel'] = ['none'] * 2  # used for NIRISS
            out_of_field_dict['PupilWheel'] = ['none'] * 2  # used for NIRISS
            out_of_field_dict['FiducialPointOverride'] = [FiducialPointOverride] * 2

            # Add out_of_field_dict to the exposures_dictionary
            exposures_dictionary = self.append_to_exposures_dictionary(exposures_dictionary,
                                                                       out_of_field_dict,
                                                                       proposal_param_dict)

        # Make sure all entries are lists
        for key in exposures_dictionary.keys():
            if type(exposures_dictionary[key]) is not list:
                exposures_dictionary[key] = list(exposures_dictionary[key])

        # Make sure all list items in the returned dictionary have the same length
        for key, item in exposures_dictionary.items():
            if len(item) == 0:
                exposures_dictionary[key] = [0] * len(exposures_dictionary['Instrument'])

        return exposures_dictionary

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
        if template_name in ['NircamImaging', 'NircamEngineeringImaging', 'NirissExternalCalibration',
                             'NirspecImaging', 'MiriMRS', 'FgsExternalCalibration']:
            parallel_exposures_dictionary = self.read_generic_imaging_template(template,
                                                                               template_name, obs,
                                                                               proposal_parameter_dictionary,
                                                                               parallel=True,
                                                                               verbose=verbose)
        elif template_name == 'NirissWfss':
            parallel_exposures_dictionary = self.read_niriss_wfss_template(template, template_name, obs,
                                                                           proposal_parameter_dictionary,
                                                                           parallel=True)
        else:
            raise ValueError('Parallel observation template {} not supported.'.format(template_name))

        # Find length of the exposures dictionary to compare
        parallel_length = len(parallel_exposures_dictionary['number_of_dithers'])
        exposures_dictionary_length = len(exposures_dictionary['number_of_dithers'])

        if parallel_length != exposures_dictionary_length:
            raise RuntimeError('Mismatch in the number of parallel observations.')

        return parallel_exposures_dictionary

    def append_to_exposures_dictionary(self, exp_dictionary, exposure_seq_dict, prop_param_dict):
        """Append exposure(s) information from a dictionary to an existing exposures dictionary

        Parameters
        ----------
        exp_dictionary : dict
            Dictionary containing information on multiple exposures

        exposure_seq_dict : dict
            Dictionary containing information on a single exposure. This dictionary should have
            the same keys as exp_dictionary. The contents of this dictionary will be added to
            exp_dictionary

        prop_param_dict : dict
            A dictionary containing proposal-wide information, such as title and PI name

        Reutrns
        -------
        exp_dictionary : dict
            With the new exposure(s) added
        """
        keys = list(exposure_seq_dict.keys())
        number_of_exposures = len(exposure_seq_dict[keys[0]])

        for key in self.APTObservationParams_keys:
            if key in exposure_seq_dict.keys():
                value = exposure_seq_dict[key]
            elif key in prop_param_dict.keys():
                value = [prop_param_dict[key]] * number_of_exposures
            else:
                value = [str(None)] * number_of_exposures

            if (key in ['PrimaryDithers', 'ImageDithers']) and ((value is None) or (value == 'None')):
                value = ['1'] * number_of_exposures
            exp_dictionary[key].extend(value)

        # add keys that were not defined in self.APTObservationParams_keys
        for key in exposure_seq_dict.keys():
            if key not in self.APTObservationParams_keys:
                # if key not yet present, create entry
                if key not in exp_dictionary.keys():
                    print('Key {} not present in APTObservationParams nor exposures_dictionary'.format(key))
                    exp_dictionary[key] = [str(exposure_seq_dict[key])]
                else:
                    print('Key {} not present in APTObservationParams'.format(key))
                    exp_dictionary[key].append(str(exposure_seq_dict[key]))
        return exp_dictionary

    def read_niriss_wfss_template(self, template, template_name, obs, proposal_param_dict, parallel=False,
                                  verbose=False):
        """Parse a NIRISS WFSS observation template from an APT xml file. Produce an exposure dictionary
        that lists all exposures (excluding dithers) from the template.

        Parameters
        ----------
        template : lxml.etree._Element
            Template section from APT xml

        template_name : str
            The type of template (e.g. 'NirissWfss')

        obs : lxml.etree._Element
            Observation section from APT xml

        proposal_param_dict : dict
            Dictionary of proposal level information from the xml file
            (e.g. PI, Science Category, etc)

        parallel : bool
            If True, template should be for parallel observations. If False, NIRISS WFSS
            observation is assumed to be prime

        Returns
        -------
        exposures_dictionary : dict
            Dictionary containing details on all exposures contained within the template. These details
            include things like filter, pupil, readout pattern, subarray, etc

        exp_len : int
            Dictionary length to use when comparing to that from a parallel observation. This is not
            necessarily the same as the true length of the dictionary due to the way in which APT
            groups overvations
        """
        instrument = 'NIRISS'

        # Dummy module name for NIRISS. Needed for consistency in dictionary entry
        mod = 'N'
        subarr = 'FULL'
        long_filter = 'N/A'
        long_pupil = 'N/A'

        # Dictionary that holds the content of this observation only
        exposures_dictionary = copy.deepcopy(self.empty_exposures_dictionary)

        # Set namespace
        ns = "{http://www.stsci.edu/JWST/APT/Template/NirissWfss}"

        # Template from the prime instrument is needed if WFSS is parallel to a nircam imaging observation.
        # In that case we need to look into the nircam observation to see if the niriss direct images are
        # to be dithered
        if verbose:
            print("Reading NIRISS WFSS template")
        if parallel:
            prime_template = obs.find(self.apt + 'Template')[0]
            prime_template_name = etree.QName(prime_template).localname
            prime_ns = "{{{}/Template/{}}}".format(self.apt.replace('{', '').replace('}', ''), prime_template_name)

            # Boolean indicating which instrument is not prime but parallel
            parallel_instrument = True
            prime_instrument = obs.find(self.apt + 'Instrument').text
            if verbose:
                print('Prime: {}   Parallel: {}'.format(prime_instrument, instrument))
            pdither_grism = prime_template.find(prime_ns + 'PrimaryDithers').text
            pdither_type_grism = prime_template.find(prime_ns + 'PrimaryDitherType').text
            dither_direct = prime_template.find(prime_ns + 'DitherNirissWfssDirectImages').text
            sdither_type_grism = prime_template.find(prime_ns + 'CoordinatedParallelSubpixelPositions').text
            try:
                sdither_grism = str(np.int(sdither_type_grism[0]))
            except ValueError:
                sdither_grism = prime_template.find(prime_ns + 'SubpixelPositions').text
        else:
            parallel_instrument = False
            prime_instrument = instrument
            dither_direct = 'NO_DITHERING'
            sdither_type_grism = 'None'
            sdither_grism = '1'
            # Dither size can be SMALL, MEDIUM, LARGE. Only look for this if NIRISS is prime
            #pdither_type_grism = template.find(ns + 'DitherSize').text
            try:
                pdither_type_grism = template.find(ns + 'PrimaryDitherType').text
            except AttributeError:
                pdither_type_grism = 'None'

            # Can be various types
            # WFSS stand-alone observation or WFSS as parallel to NIRCam prime imaging:
            # this will be the number of primary dithers.
            # WFSS as prime, with NIRCam imaging as parallel, this will be a string.
            # (e.g. '2-POINT-LARGE-NIRCam')
            dvalue = template.find(ns + 'PrimaryDithers').text
            try:
                pdither_grism = str(np.int(dvalue))
            except ValueError:
                # When NIRISS is prime with NIRCam parallel, the PrimaryDithers field can be
                # (e.g. '2-POINT-LARGE-NIRCAM'), where the first character is always the number
                # of dither positions. Not sure how to save both this name as well as the DitherSize
                # value. I don't think there are header keywords for both, with PATTTYPE being the
                # only keyword for dither pattern names.
                pdither_grism = str(np.int(dvalue[0]))

        # Check if this observation has parallels
        coordinated_parallel = obs.find(self.apt + 'CoordinatedParallel').text

        explist = template.find(ns + 'ExposureList')
        expseqs = explist.findall(ns + 'ExposureSequences')

        # Determine if there is an aperture override
        override = obs.find('.//' + self.apt + 'FiducialPointOverride')
        FiducialPointOverride = True if override is not None else False

        delta_exp_dict_length = 0
        for expseq in expseqs:
            # Grism values are listed for each ExposureSequence
            grismval = expseq.find(ns + 'Sequence').text
            if grismval == 'BOTH':
                grismval = ['GRISM150R', 'GRISM150C']
                both_grisms = True
                entry_repeats = [2, 3]
            else:
                grismval = [grismval]
                both_grisms = False
                entry_repeats = [3]
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
                    pdither_type = 'None'  # direct image has no dithers
                    sdither = '1'
                    sdither_type = 'None'
                else:
                    pdither = pdither_grism
                    pdither_type = pdither_type_grism
                    sdither = sdither_grism
                    sdither_type = sdither_type_grism

                tile = '1'
                direct_grismvalue = 'N/A'
                #pupil = 'CLEARP'  # NIRISS filter MUST be in filter wheel, not PUPIL wheel
                pupil = 'CLEAR'
                rpatt = directexp.find(ns + 'ReadoutPattern').text
                grps = directexp.find(ns + 'Groups').text
                ints = directexp.find(ns + 'Integrations').text

                # Collect info on grism exposure
                grismexp = expseq.find(ns + 'GrismExposure')
                grism_typeflag = 'wfss'
                grism_pupil = grism
                grism_rpatt = grismexp.find(ns + 'ReadoutPattern').text
                grism_grps = grismexp.find(ns + 'Groups').text
                grism_ints = grismexp.find(ns + 'Integrations').text

                # Update values in dictionary
                repeats = entry_repeats[grism_number]

                exp_seq_dict['Module'] = ['N'] * repeats
                exp_seq_dict['Subarray'] = ['FULL'] * repeats  # Niriss WFSS is always full frame
                exp_seq_dict['CoordinatedParallel'] = [coordinated_parallel] * repeats
                exp_seq_dict['Instrument'] = [instrument] * repeats
                exp_seq_dict['ParallelInstrument'] = [parallel_instrument] * repeats
                exp_seq_dict['ShortFilter'] = [filter_name] * repeats
                exp_seq_dict['LongFilter'] = [long_filter] * repeats
                exp_seq_dict['LongPupil'] = [long_pupil] * repeats
                exp_seq_dict['ObservationID'] = [proposal_param_dict['ObservationID']] * repeats
                exp_seq_dict['TileNumber'] = [tile] * repeats
                exp_seq_dict['APTTemplate'] = [template_name] * repeats
                exp_seq_dict['ObservationName'] = [proposal_param_dict['ObservationName']] * repeats
                #exp_seq_dict['FilterWheel'] = [filter_name] * repeats
                exp_seq_dict['PupilWheel'] = [filter_name] * repeats
                exp_seq_dict['FiducialPointOverride'] = [FiducialPointOverride] * repeats

                if not both_grisms:
                    exp_seq_dict['Mode'] = [typeflag, grism_typeflag, typeflag]
                    exp_seq_dict['PrimaryDitherType'] = [pdither_type, pdither_type_grism, pdither_type]
                    exp_seq_dict['PrimaryDithers'] = [pdither, pdither_grism, pdither]
                    exp_seq_dict['SubpixelPositions'] = [sdither, sdither_grism, sdither]
                    exp_seq_dict['SubpixelDitherType'] = [sdither_type, sdither_type_grism, sdither_type]
                    exp_seq_dict['ReadoutPattern'] = [rpatt, grism_rpatt, rpatt]
                    exp_seq_dict['Groups'] = [grps, grism_grps, grps]
                    exp_seq_dict['Integrations'] = [ints, grism_ints, ints]
                    exp_seq_dict['ShortPupil'] = [pupil, grism_pupil, pupil]
                    exp_seq_dict['Grism'] = [direct_grismvalue, grism, direct_grismvalue]
                    exp_seq_dict['number_of_dithers'] = [str(int(pdither)*int(sdither)),
                                                         str(int(pdither_grism)*int(sdither_grism)),
                                                         str(int(pdither)*int(sdither))]
                    #exp_seq_dict['PupilWheel'] = [pupil, grism_pupil, pupil]
                    exp_seq_dict['FilterWheel'] = [pupil, grism_pupil, pupil]
                else:
                    if grism_number == 0:
                        exp_seq_dict['Mode'] = [typeflag, grism_typeflag]
                        exp_seq_dict['PrimaryDitherType'] = [pdither_type, pdither_type_grism]
                        exp_seq_dict['PrimaryDithers'] = [pdither, pdither_grism]
                        exp_seq_dict['SubpixelPositions'] = [sdither, sdither_grism]
                        exp_seq_dict['SubpixelDitherType'] = [sdither_type, sdither_type_grism]
                        exp_seq_dict['ReadoutPattern'] = [rpatt, grism_rpatt]
                        exp_seq_dict['Groups'] = [grps, grism_grps]
                        exp_seq_dict['Integrations'] = [ints, grism_ints]
                        exp_seq_dict['ShortPupil'] = [pupil, grism_pupil]
                        exp_seq_dict['Grism'] = [direct_grismvalue, grism]
                        exp_seq_dict['number_of_dithers'] = [str(int(pdither)*int(sdither)),
                                                             str(int(pdither_grism)*int(sdither_grism))]
                        #exp_seq_dict['PupilWheel'] = [pupil, grism_pupil]
                        exp_seq_dict['FilterWheel'] = [pupil, grism_pupil]
                    elif grism_number == 1:
                        exp_seq_dict['Mode'] = [typeflag, grism_typeflag, typeflag]
                        exp_seq_dict['PrimaryDitherType'] = [pdither_type, pdither_type_grism, pdither_type]
                        exp_seq_dict['PrimaryDithers'] = [str(int(pdither)+1), pdither_grism, pdither]
                        exp_seq_dict['SubpixelPositions'] = [sdither, sdither_grism, sdither]
                        exp_seq_dict['SubpixelDitherType'] = [sdither_type, sdither_type_grism, sdither_type]
                        exp_seq_dict['ReadoutPattern'] = [rpatt, grism_rpatt, rpatt]
                        exp_seq_dict['Groups'] = [grps, grism_grps, grps]
                        exp_seq_dict['Integrations'] = [ints, grism_ints, ints]
                        exp_seq_dict['ShortPupil'] = [pupil, grism_pupil, pupil]
                        exp_seq_dict['Grism'] = [direct_grismvalue, grism, direct_grismvalue]
                        exp_seq_dict['number_of_dithers'] = [str((int(pdither)*int(sdither))*2),
                                                             str(int(pdither_grism)*int(sdither_grism)),
                                                             str(int(pdither)*int(sdither))]
                        #exp_seq_dict['PupilWheel'] = [pupil, grism_pupil, pupil]
                        exp_seq_dict['FilterWheel'] = [pupil, grism_pupil, pupil]
                #######################################################################
                # Update exposure dictionary to return
                # Add out_of_field_dict to the exposures_dictionary
                exposures_dictionary = self.append_to_exposures_dictionary(exposures_dictionary,
                                                                           exp_seq_dict,
                                                                           proposal_param_dict)
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
        return exposures_dictionary

    def separate_pupil_and_filter(self, filter_string):
        """Filters listed for NIRCam observations can take the form 'F164N+F444W' in cases
        where filters in the filter wheel and the pupil wheel are used in combination. This
        function separates the two values.

        Parameters
        ----------
        filter_string : str
            Filter name as given in xml file from APT

        Returns
        -------
        filter_name : str
            Name of the filter in the filter wheel

        pupil_name : str
            Name of the filter in the pupil wheel
        """
        if '+' in filter_string:
            pupil_name, filter_name = filter_string.split('+')
            pupil_name = pupil_name.strip()
            filter_name = filter_name.strip()
        else:
            pupil_name = 'CLEAR'
            filter_name = filter_string
        return pupil_name, filter_name


def get_guider_number(xml_file, observation_number):
    """"Parse the guider number for a particular FGSExternalCalibration observation.
    """
    observation_number = int(observation_number)
    apt_namespace = '{http://www.stsci.edu/JWST/APT}'
    fgs_namespace = '{http://www.stsci.edu/JWST/APT/Template/FgsExternalCalibration}'

    with open(xml_file) as f:
        tree = etree.parse(f)

    observation_data = tree.find(apt_namespace + 'DataRequests')
    observation_list = observation_data.findall('.//' + apt_namespace + 'Observation')
    for obs in observation_list:
        if int(obs.findtext(apt_namespace + 'Number')) == observation_number:
            detector = obs.findtext('.//' + fgs_namespace + 'Detector')
            number = detector[-1]
            return number

    raise RuntimeError('Could not find guider number in observation {} in {}'.format(observation_number, xml_file))
