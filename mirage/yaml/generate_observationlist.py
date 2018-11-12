"""Generate observation list files based on default values and APT output files.

Authors
-------
    - Lauren Chambers
    - Johannes Sahlmann

Use
---
    ::
        from mirage.yaml import write_observationlist
        write_observationlist.write_yaml(xml_file, pointing_file,
            yaml_file, ps_cat_sw, ps_cat_lw)

TODO
----
    - Determine solution to set default parameters explicitly in configuration file or similar.
    - Clarify use and role of FilterConfig

"""
import collections
import copy

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
    indict :dict
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

        if expand_parallel_dithers:
            dither_table = table[i:i + 2]

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
    # if verbose:
    #     print('Summary of observation dictionary:')
    #     for key in xml_dict.keys():
    #         print('{:<25}: number of elements is {:>5}'.format(key, len(xml_dict[key])))

    # create an expanded dictionary that contains lists of parameters expanded for dithers
    xml_dict = expand_for_dithers(xml_dict)

    return_dict = None


    # array of unique instrument names
    used_instruments = np.unique(xml_dict['Instrument'])
    unique_observation_ids = np.unique(xml_dict['ObservationID']).tolist()
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

                    if len(catalogs[key][module_key]) != number_of_observations:
                        raise RuntimeError('Please specify one catalog per observation for {}'.format(key.lower()))

        else:
            catalog_files = catalogs[key]
            if isinstance(catalog_files, str):
                catalog_file_list = [catalog_files] * number_of_observations
                catalogs[key] = catalog_file_list

            if len(catalogs[key]) != number_of_observations:
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


    observation_numbers = np.unique(xml_dict['ObservationID'])
    for observation_index, observation_number in enumerate(observation_numbers):
        first_index = xml_dict['ObservationID'].index(observation_number)
        text += [
            "Observation{}:\n".format(observation_number),
            "  Name: '{}'\n".format(xml_dict['ObservationName'][first_index])
            ]
        for index in np.where(np.array(xml_dict['ObservationID']) == observation_number)[0]:
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

                    text += [
                        "    FilterConfig:\n",
                        "      SW:\n",
                        "        Filter: {}\n".format(sw_filt),
                        "        PointSourceCatalog: {}\n".format(catalogs[instrument.lower()]['sw'][observation_index]),
                        "        BackgroundRate: {}\n".format(default_values['BackgroundRate_sw']),
                        ]
                    for key in default_parameter_name_list:
                        text += ["        {}: {}\n".format(key, default_values[key])]
                    text += [
                        "      LW:\n",
                        "        Filter: {}\n".format(lw_filt),
                        "        PointSourceCatalog: {}\n".format(
                            catalogs[instrument.lower()]['lw'][observation_index]),
                        "        BackgroundRate: {}\n".format(
                            default_values['BackgroundRate_lw']),
                        ]
                    for key in default_parameter_name_list:
                        text += ["        {}: {}\n".format(key, default_values[key])]

                elif instrument.lower() in ['niriss', 'fgs', 'nirspec', 'miri']:
                    text += [
                        "    Filter: {}\n".format(xml_dict['FilterWheel'][index]),
                        "    Pupil: {}\n".format(xml_dict['PupilWheel'][index]),
                        "    PointSourceCatalog: {}\n".format(catalogs[instrument.lower()][observation_index]),
                        "    BackgroundRate: {}\n".format(default_values['BackgroundRate']),
                        ]
                    for key in default_parameter_name_list:
                        text += ["    {}: {}\n".format(key, default_values[key])]

                entry_numbers.append(entry_number)
                entry_number += 1

    text_out += text

    return_dict['entry_number'] = entry_numbers

    f = open(yaml_file, 'w')
    for line in text_out:
        f.write(line)
    f.close()
    print('\nWrote {} observations and {} entries to {}'.format(len(observation_numbers), entry_number, yaml_file))

    return return_dict
