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


# Get the mapping between the user-input catalog dictionary keys and
# the keys used in the full table of observation parameters
CAT_TYPE_MAPPING = {'point_souce': 'PointsourceCatalog', 'galaxy': 'GalaxyCatalog',
                    'extended': 'ExtendedCatalog', 'moving_pointsource': 'MovingTargetList',
                    'moving_sersic': 'MovingTargetSersic', 'moving_extended': 'MovingTargetExtended',
                    'moving_target_to_track': 'MovingTargetToTrack'}

POSSIBLE_CATS = list(CAT_TYPE_MAPPING.keys())


def catalog_dictionary_per_observation(cats, obs_nums, targets, defaults):
    """Translate a dictionary of catalogs from a case of either:

    1. Separate catalogs for each target name
    2. Separate catalogs for each target name and instrument

    into a dictionary of catalogs for each instrument and observation

    Parameters
    ----------
    cats : dict
        Dictionary of catalogs. Can be:

        Same catalogs for all instruments within each observation
        catalogs = {'my_targ_1': {'point_source': 'ptsrc1.cat',
                                  'galaxy': 'galaxy1.cat',
                                  'extended': 'ex1.cat'},
                    'my_targ_2': {'point_source': 'ptsrc2.cat',
                                  'galaxy': 'galaxy2.cat',
                                  'extended': 'ex2.cat'}}

        Different catalogs for each instrument in each observation
        catalogs = {'my_targ_1': {'nircam': {'point_source': 'ptsrc1.cat',
                                             'galaxy': 'galaxy1.cat',
                                             'extended': 'ex1.cat'},
                                  'niriss': {'pointsource': 'ptsrc_nis.cat',
                                             'galaxy': 'galaxy_nis.cat'}},
                    'my_targ_2': {'nircam': {'point_source': 'ptsrc2.cat',
                                             'galaxy': 'galaxy2.cat',
                                             'extended': 'ex2.cat'}}}

    obs_nums : numpy.ndarray
        1D array of observation ID numbers

    targets : numpy.ndarray
        1d array of target names, with a 1:1 correspondence to obs_nums

    defaults : dict
        Dictionary of default catalog values

    Returns
    -------
    obs_cats : dict
        Dictionary of catalogs per observation, with keys that match
        those in the defaults

        obs_cats = {'001': {'nircam': {'PointsourceCatalog': 'ptsrc1.cat',
                                       'GalaxyCatalog': 'galaxy1.cat',
                                       'ExtendedCatalog': 'ex1.cat'},
                            'niriss': {'PointsourceCatalog': 'ptsrc_nis.cat',
                                       'GalaxyCatalog': 'galaxy_nis.cat'},
                            }
                    '002': {'nircam': {'PointsourceCatalog': 'ptsrc2.cat',
                                       'GalaxyCatalog': 'galaxy2.cat',
                                       'ExtendedCatalog': 'ex2.cat'},
                            'niriss': {'PointsourceCatalog': 'ptsrc_nis2.cat',
                                       'GalaxyCatalog': 'galaxy_nis2.cat'}
                            }
                   }

    """
    # Set up the output dictionary. Populate with keys for all observations
    # and default catalog values to cover any entries in obs_cats that are
    # note present
    obs_cats = {}
    for number in obs_nums:
        obs_cats[number] = {'nircam': {}, 'niriss': {}, 'fgs': {}}
        for cat_type in POSSIBLE_CATS:
            obs_cats[number]['nircam'][CAT_TYPE_MAPPING[cat_type]] = defaults[CAT_TYPE_MAPPING[cat_type]]

    # Loop over the keys in the top level of the input dictionary
    for key1 in cats:

        # Find the observation numbers that use this target
        match = targets == key1

        # Check to see if the second level of the input dictionary is
        # a dictionary of catalogs, or a dictionary of instruments
        keys2 = cats[key1].keys()
        keys_present = [True if poss in cat_keys else False for poss in POSSIBLE_CATS]
        if any(keys_present):
            # Dictionary contains catalog names, so we use the same catalogs
            # for all instruments

            # Loop over the observation numbers that use this target and
            # populate the entries for each with the catalog names. In
            # this case the catalog names are the same for all instruments
            for obs_number in obs_nums[match]:
                for key2 in keys2:
                    obs_cats[obs_number]['nircam'][CAT_TYPE_MAPPING[key2]] = cats[key1][key2]
                    obs_cats[obs_number]['niriss'][CAT_TYPE_MAPPING[key2]] = cats[key1][key2]
                    obs_cats[obs_number]['fgs'][CAT_TYPE_MAPPING[key2]] = cats[key1][key2]
        else:
            # Dictionary contains instrument names
            # Loop over observation numbers that use this target and
            # populate the different catalogs for each instrument
            for obs_number in obs_nums[match]:
                for instrument in keys2:
                    ctypes = cats[key1][key2].keys()
                    for ctype in ctypes:
                        obs_cats[obs_number][instrument][CAT_TYPE_MAPPING[ctype]] = cats[key1][instrument][ctype]
    return obs_cats


def convert_background_dict(bkgd):
    """Given a dictionary of background rates where the keys are observation
    numbers and the values are strings or numbers, expand the dictionary to
    contain entries for each instrument, and channels in the case of nircam.
    This function primarily makes sense as a way to enable users to keep a
    single background level 'high', 'medium', 'low' that applies to all
    insturuments in a given observation.

    Parameters
    ----------
    bkgd : dict
        For example:
        background = {'001': 'high', '002': 'medium', '003': 2.3}

    Return
    ------
    new_bkgd : dict
        For example:
        background = {'001': {'nircam': {'sw': high, 'lw': high}, 'niriss': high, 'fgs': high},
                      '002': {'nircam': {'sw': medium, 'lw': medium}, 'niriss': medium, 'fgs': medium},
                      '003': {'nircam': {'sw': 2.3, 'lw': 2.3}, 'niriss': 2.3, 'fgs': 2.3}
    """
    new_bkgd = {}
    for key, value in bkgd.items():
        new_bkgd[key] = {'nircam': {'sw': value, 'lw': value}, 'niriss': value, 'fgs': value}
    return new_bkgd


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

    # use astropy table operations to expand dithers while maintaining parallels in sync
    # implementation assumes that only one instrument is used in parallel
    table = Table(indict)
    table['row'] = np.arange(len(table))
    expanded_table = None #copy.deepcopy(table)

    # complication here is to handle cases with unsupported instruments (MIRI, NIRSpec) in parallel
    for i, row in enumerate(table['row']):
        number_of_dithers = np.int(table['number_of_dithers'][i])
        expand_prime_dithers_only = False
        expand_parallel_dithers = False

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


def get_observation_dict(xml_file, yaml_file, catalogs, parameter_overrides=None, verbose=False):
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
    # Read in filters from APT .xml file
    readxml_obj = read_apt_xml.ReadAPTXML()

    xml_dict = readxml_obj.read_xml(xml_file, verbose=verbose)

    # if verbose:
    #print('Summary of observation dictionary:')
    #for key in xml_dict.keys():
    #    print('{:<25}: number of elements is {:>5}'.format(key, len(xml_dict[key])))

    # create an expanded dictionary that contains lists of parameters expanded for dithers
    xml_dict = expand_for_dithers(xml_dict, verbose=verbose)
    #print('Summary of observation dictionary after expanding for dithers:')
    #for key in xml_dict.keys():
    #    print('{:<25}: number of elements is {:>5}'.format(key, len(xml_dict[key])))
    return_dict = None

    # array of unique instrument names
    used_instruments = np.unique(xml_dict['Instrument'])
    all_observation_ids = xml_dict['ObservationID']
    unique_observation_ids = np.unique(all_observation_ids).tolist()

    # List of target names from the proposal
    all_targets = readxml_obj.target_names

    # Only require the number of catalogs equal to the number of observations
    # for each instrument. Keep in mind that multiple instruments can be involved in
    # a given observation due to parallels. But in the case of serial observations
    # with different instruments, we don't want to over-count observations and
    # require more catalogs than are really necessary
    #number_of_obs = {}
    #for instrument_name in used_instruments:
    #    inst_observations = np.array(np.array(xml_dict['Instrument']) == instrument_name)
    #    unique_inst_obs = np.unique(np.array(xml_dict['ObservationID'])[inst_observations])
    #    number_of_obs[instrument_name.lower()] = len(unique_inst_obs)


    # Translate the dictionary of catalogs into a per-obervation, per-instrument dictionary

    # Same catalogs for all observations and instruments
    # catalogs = {'pointsource': 'ptsrc.cat', 'galaxy': 'gal.cat'}

    # Same catalogs for all instruments within each observation
    # catalogs = {'my_targ_1': {'pointsource': 'ptsrc1.cat', 'galaxy': 'galaxy1.cat', 'extended': 'ex1.cat'},
    #             'my_targ_2': {'pointsource': 'ptsrc2.cat', 'galaxy': 'galaxy2.cat', 'extended': 'ex2.cat'}}

    # Different catalogs for each instrument in each observation
    # catalogs = {'my_targ_1': {'nircam': {'pointsource': 'ptsrc1.cat', 'galaxy': 'galaxy1.cat', 'extended': 'ex1.cat'},
    #                           'niriss': {'pointsource': 'ptsrc_nis.cat', 'galaxy': 'galaxy_nis.cat'}},
    #             'my_targ_2': {'nircam': {'pointsource': 'ptsrc2.cat', 'galaxy': 'galaxy2.cat', 'extended': 'ex2.cat'}}}









    """
    ######################PREVIOUS CODE#########################
    # ensure that catalog files are lists with number of elements matching the number of observations
    if not isinstance(catalogs, collections.Mapping):
        raise TypeError('Please provide a catalog dictionary.')

    # Determine what kind of catalog has been provided
    # input_catalogs = copy.deepcopy(catalogs)
    if any(x in ['nircam', 'fgs', 'niriss'] for x in catalogs.keys()):
        cat_type = 'inst'
    elif set(catalogs.keys()) == set(['gal', 'ptsrc']):
        cat_type = 'object'
    else:
        raise ValueError('Invalid catalog dictionary provided.')

    # Process instrument-specific catalogs (assuming only point source)
    input_catalogs = copy.deepcopy(catalogs)
    ptsrc_catalog_lists = {}
    gal_catalog_lists = {}
    if cat_type == 'inst':
        # Strip out catalogs for any instruments that aren't used
        for key in catalogs:
            if key.upper() not in used_instruments:
                del input_catalogs[key]

        # Check that the appropriate catalogs have been included
        for inst in used_instruments:
            if inst.lower() not in input_catalogs.keys():
                raise KeyError('Missing a catalog entry for {} in the catalog dictionary.'.format(inst))

        for inst in input_catalogs.keys():
            if inst.lower() == 'nircam':
                # check that a dictionary is provided for nircam
                if not isinstance(input_catalogs[inst], collections.Mapping):
                    raise ValueError('Please provide a lw/sw dictionary for nircam.')
                else:
                    ptsrc_catalog_lists[inst] = {}
                    gal_catalog_lists[inst] = {}
                    for module_key in input_catalogs[inst].keys():
                        catalog_files = input_catalogs[inst][module_key]

                        # If the provided catalog is a string, multiply it by
                        # n_exposures
                        if isinstance(catalog_files, str):
                            catalog_file_list = [catalog_files] * number_of_obs[inst]
                            ptsrc_catalog_lists[inst][module_key] = catalog_file_list
                        # If the provided catalog is a list, use that
                        elif isinstance(catalog_files, list):
                            ptsrc_catalog_lists[inst][module_key] = catalog_files
                        else:
                            raise TypeError('Must provide catalogs as file path '
                                            'strings, or a list of file path strings.')
                        gal_catalog_lists[inst][module_key] = ['None'] * number_of_obs[inst]

                        if len(ptsrc_catalog_lists[inst][module_key]) != number_of_obs[inst]:
                            raise RuntimeError(('Please specify one catalog per observation for {}. '
                                                'Current catalog is {}').format(key.lower(),
                                                                                catalog_file_list))

            else:
                catalog_files = input_catalogs[inst]
                if isinstance(catalog_files, str):
                    catalog_file_list = [catalog_files] * number_of_obs[inst]
                    ptsrc_catalog_lists[inst] = catalog_file_list
                    gal_catalog_lists[inst] = ['None'] * number_of_obs[inst]
                else:
                    raise TypeError('Must provide catalogs as file path strings.')

                if len(catalog_file_list) != number_of_obs[inst]:
                    raise RuntimeError(('Please specify one catalog per observation for {}. '
                                        'Current catalog is {}').format(inst.lower(), catalog_file_list))

    # Process object-specific catalogs (expecting one point source and one galactic catalog,
    # each containing columns for all instruments)
    elif cat_type == 'object':
        ptsrc_catalog = input_catalogs['ptsrc']
        gal_catalog = input_catalogs['gal']
        if not (isinstance(ptsrc_catalog, str) and isinstance(gal_catalog, str)):
            raise TypeError('Must provide catalogs as file path strings.')

        for inst in used_instruments:
            inst = inst.lower()
            ptsrc_catalog_file_list = [ptsrc_catalog] * number_of_obs[inst]
            gal_catalog_file_list = [gal_catalog] * number_of_obs[inst]

            if inst == 'nircam':
                ptsrc_catalog_lists[inst] = {}
                gal_catalog_lists[inst] = {}
                ptsrc_catalog_lists[inst]['sw'] = ptsrc_catalog_file_list
                ptsrc_catalog_lists[inst]['lw'] = ptsrc_catalog_file_list
                gal_catalog_lists[inst]['sw'] = gal_catalog_file_list
                gal_catalog_lists[inst]['lw'] = gal_catalog_file_list
            else:
                ptsrc_catalog_lists[inst] = ptsrc_catalog_file_list
                gal_catalog_lists[inst] = gal_catalog_file_list
    ######################PREVIOUS CODE#########################
    """







    # if verbose:
    #     print('Summary of dictionary extracted from {}'.format(xml_file))
    #     for key in xml_dict.keys():
    #         print('{:<25}: number of elements is {:>5}'.format(key, len(xml_dict[key])))


    # Set default values. These are overwritten if there is an appropriate
    # entry in parameter_defaults
    default_values = {}
    default_values['Date'] = '2020-02-20'
    default_values['PAV3'] = '0.'
    default_values['PointsourceCatalog'] = 'None'
    default_values['GalaxyCatalog'] = 'None'
    default_values['ExtendedCatalog'] = 'None'
    default_values['ExtendedScale'] = '1.0'
    default_values['ExtendedCenter'] = '1024,1024'
    default_values['MovingTargetList'] = 'None'
    default_values['MovingTargetSersic'] = 'None'
    default_values['MovingTargetExtended'] = 'None'
    default_values['MovingTargetConvolveExtended'] = 'True'
    default_values['MovingTargetToTrack'] = 'None'
    default_values['BackgroundRate_sw'] = 'low'
    default_values['BackgroundRate_lw'] = 'low'
    default_values['BackgroundRate'] = '0.5'
    default_values['CosmicRayLibrary'] = 'SUNMAX'
    default_values['CosmicRayScale'] = 1.0
    default_parameter_name_list = ['MovingTargetConvolveExtended', 'ExtendedScale', 'ExtendedCenter']

    # Cosmic rays
    # Can be:
    # cr = {'library': 'SUNMAX', 'scale': 1.0}
    # cr = {'001': {'library': 'SUNMAX', 'scale': 1.2}}
    cosmic_rays = parameter_overrides['cosmic_rays']
    # Case where one value of library and scale are to be used for
    # all observations
    if cosmic_rays is not None:
        if 'library' in cosmic_rays.keys():
            default_values['CosmicRayLibrary'] = cosmic_rays['library']
            default_values['CosmicRayScale'] = cosmic_rays['scale']
            # Now set cosmic_rays to None so that it won't be used when looping
            # over observations below
            cosmic_rays = None
        else:
            # Case where different values are given for different observations
            # Just use cosmic_rays below when looping over observations
            pass

    # Background levels
    # background = 'high'
    # background = 22.2
    # background = {'001': 'high', '002': 'medium', '003': 22.3}
    # background = {'001': {'nircam': {'sw': 0.2, 'lw':0.3}, 'niriss': 0.4, 'fgs': 0.2}}
    background = parameter_overrides['background']
    if background is not None:
        if isinstance(background, str) or isinstance(background, float) or isinstance(background, int):
            default_values['BackgroundRate_sw'] = background
            default_values['BackgroundRate_lw'] = background
            default_values['BackgroundRate'] = background
            # Now set background to None so that it won't be used when looping
            # over observations below
            background = None
        else:
            bkeys = list(background.keys())
            if (isinstance(background[bkeys[0]], str) or isinstance(background[bkeys[0]], float) or isinstance(background[bkeys[0]], int)):
                background = convert_background_dict(background)
            else:
                # Just use background as-is below when looping over observations
                pass

    # Dates
    # dates = '2019-5-25'
    # dates = {'001': '2019-05-25', '002': '2019-11-15'}
    dates = parameter_overrides['dates']
    if dates is not None:
        if isinstance(dates, str):
            default_values['Date'] = dates
            # Now set dates to None so that it won't be used when looping
            # over observations below
            dates = None
        else:
            # Just use dates below when looping over observations
            pass

    #for key in parameter_defaults.keys():
    #    if key in default_values.keys():
    #        default_values[key] = parameter_defaults[key]

    # Roll angle, aka PAV3
    # pav3 = 34.5
    # pav3 = {'001': 34.5, '002': 154.5}
    pav3 = parameter_overrides['roll_angle']
    if pav3 is not None:
        if isinstance(pav3, float) or isinstance(pav3, int):
            default_values['PAV3'] = pav3
            # Now set pav3 to None so that it won't be used when looping
            # over observations below
            pav3 = None
        else:
            # Just use pav3 below when looping over observations
            pass

    # Catalogs
    # In the easy case where the same catalogs are to be used,
    # just populate the default values
    if catalogs is not None:
        cat_keys = catalogs.keys()
        keys_present = [True if poss in cat_keys else False for poss in POSSIBLE_CATS]
        if any(keys_present):
            for cat_key in cat_keys[keys_present]:
                if cat_key == 'point_source':
                    default_values['PointsourceCatalog'] = catalogs[cat_key]
                if cat_key == 'galaxy':
                    default_values['GalaxyCatalog'] = catalogs[cat_key]
                if cat_key == 'extended':
                    default_values['ExtendedCatalog'] = catalogs[cat_key]
                if cat_key == 'moving_pointsource':
                    default_values['MovingTargetList'] = catalogs[cat_key]
                if cat_key == 'moving_sersic':
                    default_values['MovingTargetSersic'] = catalogs[cat_key]
                if cat_key == 'moving_extended':
                    default_values['MovingTargetExtended'] = catalogs[cat_key]
                if cat_key == 'moving_target_to_track':
                    default_values['MovingTargetToTrack'] = catalogs[cat_key]
            # Now that we have modified the default values, set catalogs to
            # None so that it is not accessed later
            catalogs_per_observation = None
        else:
            # If the catalog dictionary is more complex, specifying different
            # catalogs for each target, or different catalogs for each instrument
            # and target, then translate this dictionary into a dictionary of
            # catalogs for each observation and instrument
            catalogs_per_observation = catalog_dictionary_per_observation(catalogs,
                                                                          np.array(all_observation_ids),
                                                                          np.array(all_targets),
                                                                          default_values)
    else:
        catalogs_per_observation = None

    # assemble string that will constitute the yaml content
    text_out = ["# Observation list created by generate_observationlist.py\n\n"]

    text = ['']
    entry_number = 0  # running number for every entry in the observation list

    # Create an instrument-specific counter to be used with input catalogs
    all_instruments = np.unique(xml_dict['Instrument'])
    counter = {}
    for inst in all_instruments:
        counter[inst] = 0

    entry_numbers = []
    observation_numbers = np.unique(xml_dict['ObservationID'])
    for observation_index, observation_number in enumerate(observation_numbers):
        first_index = xml_dict['ObservationID'].index(observation_number)
        text += [
            "Observation{}:\n".format(observation_number),
            "  Name: '{}'\n".format(xml_dict['ObservationName'][first_index])
            ]
        observation_rows = np.where(np.array(xml_dict['ObservationID']) == observation_number)[0]
        instruments_in_observation = np.unique(np.array(xml_dict['Instrument'])[observation_rows])
        for index in observation_rows:
            number_of_dithers = np.int(xml_dict['number_of_dithers'][index])
            instrument = xml_dict['Instrument'][index]
            for dither_index in range(number_of_dithers):

                # Get the proper date value
                if dates is None:
                    date_value = default_values['Date']
                else:
                    try:
                        date_value = dates[observation_number]
                    except KeyError:
                        print(("ERROR: No date value specified for Observation {} in date dictionary. "
                               "Quitting.".format(observation_number)))
                        raise KeyError

                # Get the proper PAV3 value
                if pav3 is None:
                    pav3_value = default_values['PAV3']
                else:
                    try:
                        pav3_value = pav3[observation_number]
                    except KeyError:
                        print(("ERROR: No roll angle value specified for Observation {} in roll_angle "
                               "dictionary. Quitting.".format(observation_number)))
                        raise KeyError


                # Get the proper catalog values
                if catalogs_per_observation is None:
                    ptsrc_catalog_value = default_values['PointsourceCatalog']
                    galaxy_catalog_value = default_values['GalaxyCatalog']
                    extended_catalog_value = default_values['ExtendedCatalog']
                    mov_ptsrc_catalog_value = default_values['MovingTargetList']
                    mov_sersic_catalog_value = default_values['MovingTargetSersic']
                    mov_extended_catalog_value = default_values['MovingTargetExtended']
                    mov_tracked_catalog_value = default_values['MovingTargetToTrack']
                else:
                    try:
                        catalogs_to_use = catalogs_per_observation[observation_number][instrument]
                    except KeyError:
                        print(("ERROR: Missing observation number or instrument entry in catalog "
                               "dictionary. Failed to find catalogs[{}][{}]".format(observation_number, instrument)))
                        raise KeyError
                    ptsrc_catalog_value = catalogs_to_use['PointsourceCatalog']
                    galaxy_catalog_value = catalogs_to_use['GalaxyCatalog']
                    extended_catalog_value = catalogs_to_use['ExtendedCatalog']
                    mov_ptsrc_catalog_value = catalogs_to_use['MovingTargetList']
                    mov_sersic_catalog_value = catalogs_to_use['MovingTargetSersic']
                    mov_extended_catalog_value = catalogs_to_use['MovingTargetExtended']
                    mov_tracked_catalog_value = catalogs_to_use['MovingTargetToTrack']

                # Get the proper cosmic ray values
                if cosmic_rays is None:
                    cr_library_value = default_values['CosmicRayLibrary']
                    cr_scale_value = default_values['CosmicRayScale']
                else:
                    try:
                        cr_library_value = cosmic_rays[observation_number]['library']
                        cr_scale_value = cosmic_rays[observation_number]['scale']
                    except KeyError:
                        print(("ERROR: No cosmic ray library and/or scale value specified for "
                               "Observation {} in cosmic_ray dictionary. Quitting."
                               .format(observation_number)))
                        raise KeyError

                text += [
                    "  EntryNumber{}:\n".format(entry_number),
                    "    Instrument: {}\n".format(instrument),
                    "    Date: {}\n".format(date_value),
                    "    PAV3: {}\n".format(pav3_value),
                    "    DitherIndex: {}\n".format(dither_index),
                    "    CosmicRayLibrary: {}\n".format(cr_library_value),
                    "    CosmicRayScale: {}\n".format(cr_scale_value),
                ]
                if return_dict is None:
                    return_dict = dictionary_slice(xml_dict, index)
                else:
                    return_dict = read_apt_xml.append_dictionary(return_dict, dictionary_slice(xml_dict, index))
                if instrument.lower() in ['nircam', 'wfsc']:

                    sw_filt = xml_dict['ShortFilter'][index]
                    lw_filt = xml_dict['LongFilter'][index]

                    # Get the proper background rate
                    if background is None:
                        background_sw_value = default_values['BackgroundRate_sw']
                        background_lw_value = default_values['BackgroundRate_lw']
                    else:
                        try:
                            background_sw_value = background[observation_number]['nircam']['sw']
                            background_lw_value = background[observation_number]['nircam']['lw']
                        except KeyError:
                            print(("ERROR: Missing entry in the background dictionary for NIRCam SW and/or "
                                   "LW channels, observation number: {}".format(observation_number)))
                            raise KeyError

                    text += [
                        "    FilterConfig:\n",
                        "      SW:\n",
                        "        Filter: {}\n".format(sw_filt),
                        "        PointSourceCatalog: {}\n".format(ptsrc_catalog_value),
                        "        GalaxyCatalog: {}\n".format(galaxy_catalog_value),
                        "        ExtendedCatalog: {}\n".format(extended_catalog_value),
                        "        MovingTargetList: {}\n".format(mov_ptsrc_catalog_value),
                        "        MovingTargetSersic: {}\n".format(mov_sersic_catalog_value),
                        "        MovingTargetExtended: {}\n".format(mov_extended_catalog_value),
                        "        MovingTargetToTrack: {}\n".format(mov_tracked_catalog_value),
                        "        BackgroundRate: {}\n".format(background_sw_value),
                        ]

                    for key in default_parameter_name_list:
                        text += ["        {}: {}\n".format(key, default_values[key])]

                    text += [
                        "      LW:\n",
                        "        Filter: {}\n".format(lw_filt),
                        "        PointSourceCatalog: {}\n".format(ptsrc_catalog_value),
                        "        GalaxyCatalog: {}\n".format(galaxy_catalog_value),
                        "        ExtendedCatalog: {}\n".format(extended_catalog_value),
                        "        MovingTargetList: {}\n".format(mov_ptsrc_catalog_value),
                        "        MovingTargetSersic: {}\n".format(mov_sersic_catalog_value),
                        "        MovingTargetExtended: {}\n".format(mov_extended_catalog_value),
                        "        MovingTargetToTrack: {}\n".format(mov_tracked_catalog_value),
                        "        BackgroundRate: {}\n".format(background_lw_value),
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

                    # Get the proper background rate
                    if background is None:
                        background_value = default_values['BackgroundRate']
                    else:
                        try:
                            background_value = background[observation_number][instrument]
                        except KeyError:
                            print(("ERROR: Missing entry in the background dictionary for observation "
                                   "number: {}, instrument: {}".format(observation_number, instrument)))
                            raise KeyError

                    text += [
                        "    Filter: {}\n".format(filter_value),
                        "    PointSourceCatalog: {}\n".format(ptsrc_catalog_value),
                        "    GalaxyCatalog: {}\n".format(galaxy_catalog_value),
                        "    ExtendedCatalog: {}\n".format(extended_catalog_value),
                        "    MovingTargetList: {}\n".format(mov_ptsrc_catalog_value),
                        "    MovingTargetSersic: {}\n".format(mov_sersic_catalog_value),
                        "    MovingTargetExtended: {}\n".format(mov_extended_catalog_value),
                        "    MovingTargetToTrack: {}\n".format(mov_tracked_catalog_value),
                        "    BackgroundRate: {}\n".format(background_value),
                        ]

                    for key in default_parameter_name_list:
                        text += ["    {}: {}\n".format(key, default_values[key])]

                entry_numbers.append(entry_number)
                entry_number += 1

        # Update the catalog counters for the instruments used in this observation
        #for inst_name in instruments_in_observation:
        #    counter[inst_name] += 1

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
