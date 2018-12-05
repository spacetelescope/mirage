"""Generate observation list files based on default values and APT output files.

Authors
-------
    - Lauren Chambers
    - Johannes Sahlmann

Use
---
    ::
        from mirage.yaml import generate_observationlist
        generate_observationlist.get_observation_dict(xml_file, yaml_file, catalogs,
            parameter_defaults=None, verbose=False):

TODO
----
    - Determine solution to set default parameters explicitly in configuration file or similar.
    - Clarify use and role of FilterConfig

"""
import collections
import copy
import os

from astropy.table import Table, vstack
import numpy as np

from ..apt import read_apt_xml


def dictionary_slice(dictionary, index):
    """Return a dictionary with only the i'th element from every list stored in a key.

    Parameters
    ----------
    dictionary
    index

    Returns
    -------

    """
    new_dict = {}
    for key in dictionary.keys():
        new_dict[key] = [dictionary[key][index]]
    return new_dict


def expand_for_dithers(indict, verbose=True):
    """Expand a given dictionary to create one entry for each dither.

    Supports parallel observations.

    Moved here and modified from apt_inputs.py

    Parameters
    ----------
    indict : dict
        dictionary of observations

    Returns
    -------
    expanded : dict
        Dictionary, expanded to include a separate entry for
        each dither
    """
    expanded = {}
    for key in indict:
        expanded[key] = []

    # Loop over entries in dict and duplicate by the
    # number of dither positions
    # keys = indict.keys()
    # if np.all(np.unique(indict['Instrument']) == 'NIRCAM'):
    # # if not np.all(np.array(indict['PrimaryDithers']).astype(int) == 1):
    #     # NIRCam exposures
    #     for i in range(len(indict['PrimaryDithers'])):
    #         arr = np.array([item[i] for item in indict.values()])
    #         entry = dict(zip(keys, arr))
    #
    #         # In WFSS, SubpixelPositions will be either '4-Point' or '9-Point'
    #         subpix = entry['SubpixelPositions']
    #         if subpix in ['0', 'NONE']:
    #             subpix = [[1]]
    #         if subpix == '4-Point':
    #             subpix = [[4]]
    #         if subpix == '9-Point':
    #             subpix = [[9]]
    #
    #         primary = entry['PrimaryDithers']
    #         if primary == '0':
    #             primary = [1]
    #         reps = np.int(subpix[0][0]) * np.int(primary[0])
    #         for key in keys:
    #             for j in range(reps):
    #                 expanded[key].append(indict[key][i])


    # use astropy table operations to expand dithers while maintaining parallels in sync
    # implementation assumes that only one instrument is used in parallel
    table = Table(indict)
    table['row'] = np.arange(len(table))
    expanded_table = None #copy.deepcopy(table)

    # if verbose:
    #     table['ObservationID', 'Instrument', 'CoordinatedParallel', 'ParallelInstrument', 'number_of_dithers'].pprint()


    # complication here is to handle cases with unsupported instruments (MIRI, NIRSpec) in parallel
    for i, row in enumerate(table['row']):
        number_of_dithers = np.int(table['number_of_dithers'][i])
        expand_prime_dithers_only = False
        expand_parallel_dithers = False
        # 1/0

        # skip over parallel observations because they are already accounted for
        if table['ParallelInstrument'][i]:
            print('skipping row {}'.format(i))
            print(table[i]['Instrument'])
            continue

        try:
            if (table['CoordinatedParallel'][i] == 'true') and (not table['ParallelInstrument'][i]) \
                    and (table['ParallelInstrument'][i + 1]) and (table['Instrument'][i] != table['Instrument'][i + 1]):
                expand_parallel_dithers = True
            else:
                expand_prime_dithers_only = True

        except IndexError:  # last row in table is not a parallel
            expand_prime_dithers_only = True

        if (table['CoordinatedParallel'][i] == 'false'):
            expand_prime_dithers_only = True

        if expand_prime_dithers_only and expand_parallel_dithers:
            raise RuntimeError('Possible conflict found when expanding for dithers.')

        print('table row: {}.'.format(i))
        print(table[i]['Instrument'])
        print('expand_prime_dithers_only: {}, expand_parallel_dithers: {}'.format(expand_prime_dithers_only, expand_parallel_dithers))


        if expand_parallel_dithers:
            dither_table = table[i:i + 2]

            print('expand_parallel_dithers:')
            print(dither_table, number_of_dithers)

            if (number_of_dithers > 1):
                #replicate parallel observation n times
                dither_table = vstack([dither_table]*number_of_dithers)

            if expanded_table is None:
                expanded_table = dither_table
            else:
                # if verbose:
                #     print('Parallel: Adding {:>3d} rows to table with {:>3d} rows'.format(len(dither_table), len(expanded_table)))
                expanded_table = vstack((expanded_table, dither_table))

        elif expand_prime_dithers_only:

            print('expand_prime_dithers_only:')
            print(dither_table, number_of_dithers)

            # add row multiplied by number of dithers
            dither_table = vstack([table[i]]*number_of_dithers)
            if expanded_table is None:
                expanded_table = dither_table
            else:
                # print('Prime:    Adding {:>3d} rows to table with {:>3d} rows'.format(len(dither_table),
                #                                                             len(expanded_table)))
                expanded_table = vstack((expanded_table, dither_table))

    # set number of dithers to 1 after expansion
    expanded_table['number_of_dithers'] = np.ones(len(expanded_table)).astype(np.int)

    # NIRCam cannot handle when PrimaryDithers=None
    for index, value in enumerate(expanded_table['PrimaryDithers']):
        if value == 'None':
            expanded_table['PrimaryDithers'][index] = expanded_table['number_of_dithers'][index]

    expanded = {}
    for key in expanded_table.colnames:
        expanded[key] = np.array(expanded_table[key]).tolist()

    if verbose:
        print('Number of entries before expanding dithers: {}'.format(len(table)))
        print('Number of entries after expanding dithers:  {}'.format(len(expanded_table)))

    if verbose:
        for obs_id in np.unique(expanded_table['ObservationID']):
            print('Expanded table for Observation {} has {} entries'.format(obs_id, len(np.where(expanded_table['ObservationID']==obs_id)[0])))

    return expanded


def get_observation_dict(xml_file, yaml_file, catalogs, parameter_defaults=None, verbose=False):
    """Write observation list file (required mirage input) on the basis of APT files.

    Parameters
    ----------
    xml_file : str
        path to APT .xml file
    yaml_file : str
        output_file
    catalogs : dict
        Dictionary of catalog files, one entry per instrument. For NIRCam the entry has to be a
        dictionary itself, e.g. catalogs['nircam']['lw'] = somefile
        If the user prvides a list of catalogs, that list has to have one entry per observation in
        the program, accounting for any instrument used.
    parameter_defaults : dict
        Dictionary of default parameter value, e.g. date, roll angle, ...

    Returns
    -------
    xml_dict : dict
        Expanded dictionary that holds exposure information

    TODO
    ----
        Read default values from configuration file

    """
    # avoid side effects
    catalogs = copy.deepcopy(catalogs)

    # Read in filters from APT .xml file
    readxml_obj = read_apt_xml.ReadAPTXML()

    xml_dict = readxml_obj.read_xml(xml_file, verbose=verbose)

    print('Primary dithers output from read_xml: {}'.format(xml_dict['PrimaryDithers']))
    print('Number of dithers output from read_xml: {}'.format(xml_dict['number_of_dithers']))

    print('IN generate_observationlist: ', xml_dict['Instrument'])
    print(np.unique(xml_dict['Instrument']))

    # if verbose:
    print('Summary of observation dictionary:')
    for key in xml_dict.keys():
        print('{:<25}: number of elements is {:>5}'.format(key, len(xml_dict[key])))

    # create an expanded dictionary that contains lists of parameters expanded for dithers
    xml_dict = expand_for_dithers(xml_dict)
    print('Summary of observation dictionary after expanding for dithers:')
    for key in xml_dict.keys():
        print('{:<25}: number of elements is {:>5}'.format(key, len(xml_dict[key])))

    return_dict = None

    # array of unique instrument names
    used_instruments = np.unique(xml_dict['Instrument'])
    unique_observation_ids = np.unique(xml_dict['ObservationID']).tolist()

    print('LIST OF USED INSTRUMENTS:', used_instruments)

    # Only require the number of catalogs equal to the number of observations
    # for each instrument. Keep in mind that multiple instruments can be involved in
    # a given observation due to parallels (right?) But in the case of serial observations
    # with different instruments, we don't want to over-count observations and
    # require more catalogs than are really necessary
    print('TEST TEST: make sure required number of catalogs is correct')
    number_of_obs = {}
    for instrument_name in used_instruments:
        print('CHECKING FOR ', instrument_name.lower())
        print(type(xml_dict['Instrument']))
        inst_observations = np.array(np.array(xml_dict['Instrument']) == instrument_name)
        print(inst_observations)
        print(np.array(xml_dict['ObservationID'])[inst_observations])
        unique_inst_obs = np.unique(np.array(xml_dict['ObservationID'])[inst_observations])
        print(unique_inst_obs, len(unique_inst_obs))
        number_of_obs[instrument_name.lower()] = len(unique_inst_obs)
    print('TEST TEST')


    number_of_observations = len(unique_observation_ids)
    number_of_exposures = len(xml_dict['ObservationID'])

    # Check that the appropriate catalogs have been included
    for inst in used_instruments:
        if inst.lower() not in catalogs.keys():
            raise KeyError('Missing a catalog entry for {} in the catalog dictionary.'.format(inst))

    entry_numbers = []

    # ensure that catalog files are lists with number of elements matching the number of observations
    if not isinstance(catalogs, collections.Mapping):
        raise ValueError('Please provide a catalog dictionary.')

    for key in catalogs.keys():
        catalog_file_list = None
        catalog_files = None
        if key.lower() == 'nircam':
            # check that a dictionary is provided for nircam
            if not isinstance(catalogs[key], collections.Mapping):
                raise ValueError('Please provide a lw/sw dictionary for nircam.')
            else:
                for module_key in catalogs[key].keys():
                    catalog_files = catalogs[key][module_key]
                    if isinstance(catalog_files, str):
                        catalog_file_list = [catalog_files] * number_of_observations
                        catalogs[key][module_key] = catalog_file_list

                    #if len(catalogs[key][module_key]) != number_of_observations:
                    print('TEST HERE TOO, original line commented out above')
                    print(key, number_of_obs)
                    print(len(catalogs[key][module_key]))
                    print(number_of_obs[key])
                    if len(catalogs[key][module_key]) != number_of_obs[key]:
                        raise RuntimeError('Please specify one catalog per observation for {}'.format(key.lower()))

        else:
            catalog_files = catalogs[key]
            if isinstance(catalog_files, str):
                catalog_file_list = [catalog_files] * number_of_observations
                catalogs[key] = catalog_file_list

            #if len(catalogs[key]) != number_of_observations:
            print('TEST FOR OTHER INST CATALOG LENGTHS (orig line is commented out above')
            if len(catalogs[key]) != number_of_obs[key]:
                raise RuntimeError(
                    'Please specify one catalog per observation for {}'.format(key.lower()))

    # if verbose:
    #     print('Summary of dictionary extracted from {}'.format(xml_file))
    #     for key in xml_dict.keys():
    #         print('{:<25}: number of elements is {:>5}'.format(key, len(xml_dict[key])))


    # set default values. These are overwritten if defaults argument is present
    default_values = {}
    default_values['Date'] = '2019-07-04'
    default_values['PAV3'] = '111.'
    default_values['GalaxyCatalog'] = 'None'
    default_values['ExtendedCatalog'] = 'None'
    default_values['ExtendedScale'] = '1.0'
    default_values['ExtendedCenter'] = '1024,1024'
    default_values['MovingTargetList'] = 'None'
    default_values['MovingTargetSersic'] = 'None'
    default_values['MovingTargetExtended'] = 'None'
    default_values['MovingTargetConvolveExtended'] = 'True'
    default_values['MovingTargetToTrack'] = 'None'
    default_values['BackgroundRate_sw'] = '0.5'
    default_values['BackgroundRate_lw'] = '1.2'
    default_values['BackgroundRate'] = '0.5'


    default_parameter_name_list = [key for key, item in default_values.items() if key not in 'Date PAV3 BackgroundRate BackgroundRate_sw BackgroundRate_lw'.split()]

    # set default parameters if given as argument
    if parameter_defaults is not None:
        for key in parameter_defaults.keys():
            if key in default_values.keys():
                default_values[key] = parameter_defaults[key]


    # assemble string that will constitute the yaml content
    text_out = ["# Observation list created by generate_observationlist.py\n\n"]

    text = ['']
    entry_number = 0  # running number for every entry in the observation list


    print(xml_dict['ObservationID'])
    print(xml_dict['Instrument'])


    # Create an instrument-specific counter to be used with input catalogs
    all_instruments = np.unique(xml_dict['Instrument'])
    counter = {}
    for inst in all_instruments:
        counter[inst] = 0

    observation_numbers = np.unique(xml_dict['ObservationID'])
    for observation_index, observation_number in enumerate(observation_numbers):
        first_index = xml_dict['ObservationID'].index(observation_number)
        text += [
            "Observation{}:\n".format(observation_number),
            "  Name: '{}'\n".format(xml_dict['ObservationName'][first_index])
            ]
        observation_rows = np.where(np.array(xml_dict['ObservationID']) == observation_number)[0]
        print(observation_rows, type(observation_rows))
        instruments_in_observation = np.unique(np.array(xml_dict['Instrument'])[observation_rows])
        print('Observation {}, instruments used: {}'.format(observation_index, instruments_in_observation))
        for index in observation_rows:
            number_of_dithers = np.int(xml_dict['number_of_dithers'][index])
            instrument = xml_dict['Instrument'][index]
            for dither_index in range(number_of_dithers):
                text += [
                    "  EntryNumber{}:\n".format(entry_number),
                    "    Instrument: {}\n".format(instrument),
                    "    Date: {}\n".format(default_values['Date']),
                    "    PAV3: {}\n".format(default_values['PAV3']),
                    "    DitherIndex: {}\n".format(dither_index),
                ]
                if return_dict is None:
                    return_dict = dictionary_slice(xml_dict, index)
                else:
                    return_dict = read_apt_xml.append_dictionary(return_dict, dictionary_slice(xml_dict, index))
                if instrument.lower() in ['nircam', 'wfsc']:

                    sw_filt = xml_dict['ShortFilter'][index]
                    lw_filt = xml_dict['LongFilter'][index]

                    #text += [
                    #    "    FilterConfig:\n",
                    #    "      SW:\n",
                    #    "        Filter: {}\n".format(sw_filt),
                    #    "        PointSourceCatalog: {}\n".format(catalogs[instrument.lower()]['sw'][observation_index]),
                    #    "        BackgroundRate: {}\n".format(default_values['BackgroundRate_sw']),
                    #    ]

                    text += [
                        "    FilterConfig:\n",
                        "      SW:\n",
                        "        Filter: {}\n".format(sw_filt),
                        "        PointSourceCatalog: {}\n".format(catalogs[instrument.lower()]['sw'][counter[instrument]]),
                        "        BackgroundRate: {}\n".format(default_values['BackgroundRate_sw']),
                        ]

                    for key in default_parameter_name_list:
                        text += ["        {}: {}\n".format(key, default_values[key])]
                    #text += [
                    #    "      LW:\n",
                    #    "        Filter: {}\n".format(lw_filt),
                    #    "        PointSourceCatalog: {}\n".format(
                    #        catalogs[instrument.lower()]['lw'][observation_index]),
                    #    "        BackgroundRate: {}\n".format(
                    #        default_values['BackgroundRate_lw']),
                    #    ]

                    text += [
                        "      LW:\n",
                        "        Filter: {}\n".format(lw_filt),
                        "        PointSourceCatalog: {}\n".format(
                            catalogs[instrument.lower()]['lw'][counter[instrument]]),
                        "        BackgroundRate: {}\n".format(
                            default_values['BackgroundRate_lw']),
                        ]

                    for key in default_parameter_name_list:
                        text += ["        {}: {}\n".format(key, default_values[key])]

                elif instrument.lower() in ['niriss', 'fgs', 'nirspec', 'miri']:
                    if (instrument.lower() == 'niriss') and (xml_dict['APTTemplate'][index] in ['NirissExternalCalibration']):
                        filter_wheel_value = xml_dict['FilterWheel'][index]
                        pupil_wheel_value = xml_dict['PupilWheel'][index]
                        text += [
                            "    FilterWheel: {}\n".format(filter_wheel_value),
                            "    PupilWheel: {}\n".format(pupil_wheel_value),
                            ]
                        if 'CLEAR' in filter_wheel_value:
                            filter_value = pupil_wheel_value
                        elif ('CLEAR' in pupil_wheel_value) or ('NRM' in pupil_wheel_value):
                            filter_value = filter_wheel_value
                    else:
                        filter_value = xml_dict['Filter'][index]
                    #text += [
                    #    "    Filter: {}\n".format(filter_value),
                    #    "    PointSourceCatalog: {}\n".format(catalogs[instrument.lower()][observation_index]),
                    #    "    BackgroundRate: {}\n".format(default_values['BackgroundRate']),
                    #    ]

                    text += [
                        "    Filter: {}\n".format(filter_value),
                        "    PointSourceCatalog: {}\n".format(catalogs[instrument.lower()][counter[instrument]]),
                        "    BackgroundRate: {}\n".format(default_values['BackgroundRate']),
                        ]

                    for key in default_parameter_name_list:
                        text += ["    {}: {}\n".format(key, default_values[key])]

                entry_numbers.append(entry_number)
                entry_number += 1

        # Update the catalog counters for the instruments used in this observation
        print("TESTING CATALOG COUNTER UPDATES")
        for inst_name in instruments_in_observation:
            counter[inst_name] += 1
        print('updated counters are: {}'.format(counter))

    print('BUT WHAT ABOUT THE CASE WHERE THERE ARE MULTIPLE FILTER CONFIGS *WITHIN* AN OBSERVATION??')
    print('HOW DO WE GET DIFFERENT CATALOGS IN THERE? FORCE THE USER TO MANUALLY EDIT THE OBSERVATION LIST FILE?')
    print('STRONLY SUGGEST THAT USERS MAKE CATALOGS WITH MAGNITUDE COLUMNS TO COVER ALL THE FILTERS IN THE OBSERVATION')
    print('WHAT ABOUT CASE WHERE SAME CATALOG IS TO BE USED FOR ALL OBS? USER SHOULD ONLY HAVE TO PROVIDE IT ONCE,')
    print('RATHER THAN A LONG LIST OF THE SAME CATALOG NAME OVER AND OVER.')

    text_out += text

    return_dict['entry_number'] = entry_numbers

    # If the directory to hold the observation file does not yet exist, create it
    obs_dir = os.path.dirname(yaml_file)
    if obs_dir is not '' and os.path.isdir(obs_dir) is False:
        try:
            os.mkdir(obs_dir)
        except OSError:
            print("Creation of the directory {} failed".format(obs_dir))
        else:
            print("Successfully created the directory {} to hold the observation list file.".format(obs_dir))

    f = open(yaml_file, 'w')
    for line in text_out:
        f.write(line)
    f.close()
    print('\nWrote {} observations and {} entries to {}'.format(len(observation_numbers), entry_number, yaml_file))

    return return_dict
