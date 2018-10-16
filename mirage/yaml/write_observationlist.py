"""Generate observation list files based on default values and APT
output files

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
"""
import re

from lxml import etree
import numpy as np

from ..apt import read_apt_xml


def write_yaml(xml_file, yaml_file, catalog_files, ps_cat_sw=None, ps_cat_lw=None, verbose=False):
    """Write observation list file (required mirage input) on the basis of APT files.

    Parameters
    ----------
    xml_file : str
        path to APT .xml file
    pointing_file
    yaml_file : str
        output_file
    catalog_files : str or list(str)
        path to file(s) that contain the catalog of sources.
    ps_cat_sw
    ps_cat_lw

    Returns
    -------
    xml_dict : dict
        Expanded dictionary that holds exposure information

    """
    # Read in observations from XML file
    # with open(xml_file) as f:
    #     tree = etree.parse(f)
    #
    if type(catalog_files) is str:
        catalog_files = [catalog_files]
    #
    # apt = '{http://www.stsci.edu/JWST/APT}'
    # observation_data = tree.find(apt + 'DataRequests')
    # obs_results = observation_data.findall('.//' + apt + 'Observation')
    #
    # # Collect observation names and numbers for NIRCam/WFSC observations
    # observations = []
    # i_observations = []
    # obs_names = []
    # instruments = []
    # for i_obs, o in enumerate(obs_results):
    #     instrument = o.find(apt + 'Instrument').text
    #     if instrument in ['NIRCAM', 'WFSC', 'NIRISS', 'FGS']:
    #         observations.append(o)
    #         i_observations.append(i_obs)
    #         instruments.append(instrument)
    #
    #         # Find observation name, and modify to exclude parantheticals
    #         label = o.find(apt + 'Label').text
    #         if (' (' in label) and (')' in label):
    #             label = re.split(r' \(|\)', label)[0]
    #         obs_names.append(label)
    #
    # num_obs = len(observations)

    # Read in filters from APT .xml file
    readxml_obj = read_apt_xml.ReadAPTXML()

    xml_dict = readxml_obj.read_xml(xml_file, verbose=verbose)

    # array of unique instrument names
    used_instruments = np.unique(xml_dict['Instrument'])


    if verbose:
        print('Summary of dictionary extracted from {}'.format(xml_file))
        for key in xml_dict.keys():
            print('{:<25}: number of elements is {:>5}'.format(key, len(xml_dict[key])))

    # set defaults
    date = '2019-07-04'
    PAV3 = '0.'
    GalaxyCatalog = 'None'
    ExtendedCatalog = 'None'
    ExtendedScale = '1.0'
    ExtendedCenter = '1024,1024'
    MovingTargetList = 'None'
    MovingTargetSersic = 'None'
    MovingTargetExtended = 'None'
    MovingTargetConvolveExtended = 'True'
    MovingTargetToTrack = 'None'
    BackgroundRate = '0.5'


    number_of_exposures = len(xml_dict['ObservationID'])

    # if only one catalog is provided, that catalog will be used for all observations
    if len(catalog_files) == 1:
        catalog_files = catalog_files * number_of_exposures


    # assemble string that will constitute the yaml content
    text_out = ["# Observation list created by write_observationlist.py\n\n"]

    from astropy.table import Table

    # loop over exposures
    for instrument in used_instruments:
        exposure_indices = np.where(np.array(xml_dict['Instrument']) == instrument)[0]
        # 1/0
        for exposure_index in exposure_indices:
            text = [
                "Observation{}:\n".format(xml_dict['ObservationID'][exposure_index]),
                "  Name: '{}'\n".format(xml_dict['ObservationName'][exposure_index]),
                "  Date: {}\n".format(date),
                "  PAV3: {}\n".format(PAV3),
            ]
            if instrument in ['NIRCAM', 'WFSC']:
                # LW and SW stuff
                # Define "default" values (probably should be changed eventually)
                date = '2019-07-04'
                PAV3 = '0.'
                GalaxyCatalog = 'None'
                ExtendedCatalog = 'None'
                ExtendedScale = '1.0'
                ExtendedCenter = '1024,1024'
                MovingTargetList = 'None'
                MovingTargetSersic = 'None'
                MovingTargetExtended = 'None'
                MovingTargetConvolveExtended = 'True'
                MovingTargetToTrack = 'None'
                BackgroundRate_sw = '0.5'
                BackgroundRate_lw = '1.2'

                sw_filters = {}
                lw_filters = {}
                sw_filters_all = np.array(xml_dict['ShortFilter'])
                lw_filters_all = np.array(xml_dict['LongFilter'])
                tile_nums = xml_dict['TileNumber']
                observation_ids = xml_dict['ObservationID']

                for i_obs_all in set(observation_ids):
                    # i_obs_all = int(i_obs_all)
                    current_obs_indices = [i == i_obs_all for i in observation_ids]
                    if len(set(np.array(sw_filters_all)[current_obs_indices])) > 1:
                        print('Note: Multiple filters in observation {}'.format(
        i_obs_all))
                        # At some point could use the tile_nums to fix this
                    sw_filters[i_obs_all] = sw_filters_all[current_obs_indices]
                    lw_filters[i_obs_all] = lw_filters_all[current_obs_indices]

                # # Check that all parameters have the right length
                # all_param_lengths = [len(ps_cat_sw), len(ps_cat_lw), len(sw_filters),
                #                      len(lw_filters), len(observations),len(i_observations),
                #                      len(obs_names)]
                #
                # if len(set(all_param_lengths)) > 1:
                #     print(all_param_lengths)
                #     raise ValueError('Not all provided parameters have compatible '
                #                      'dimensions. Will not write {}'.format(yaml_file))

                for i_filt, (sw_filt, lw_filt) in enumerate(zip(sw_filters[obs_number], lw_filters[obs_number])):
                        text += [
                            "  FilterConfig{}:\n".format(i_filt + 1),
                            "    SW:\n",
                            "      Filter: {}\n".format(sw_filt),
                            "      PointSourceCatalog: {}\n".format(ps_cat_sw[i_obs]),
                            "      GalaxyCatalog: {}\n".format(GalaxyCatalog),
                            "      ExtendedCatalog: {}\n".format(ExtendedCatalog),
                            "      ExtendedScale: {}\n".format(ExtendedScale),
                            "      ExtendedCenter: {}\n".format(ExtendedCenter),
                            "      MovingTargetList: {}\n".format(MovingTargetList),
                            "      MovingTargetSersic: {}\n".format(MovingTargetSersic),
                            "      MovingTargetExtended: {}\n".format(MovingTargetExtended),
                            "      MovingTargetConvolveExtended: {}\n".format(MovingTargetConvolveExtended),
                            "      MovingTargetToTrack: {}\n".format(MovingTargetToTrack),
                            "      BackgroundRate: {}\n".format(BackgroundRate_sw),
                            "    LW:\n",
                            "      Filter: {}\n".format(lw_filt),
                            "      PointSourceCatalog: {}\n".format(ps_cat_lw[i_obs]),
                            "      GalaxyCatalog: {}\n".format(GalaxyCatalog),
                            "      ExtendedCatalog: {}\n".format(ExtendedCatalog),
                            "      ExtendedScale: {}\n".format(ExtendedScale),
                            "      ExtendedCenter: {}\n".format(ExtendedCenter),
                            "      MovingTargetList: {}\n".format(MovingTargetList),
                            "      MovingTargetSersic: {}\n".format(MovingTargetSersic),
                            "      MovingTargetExtended: {}\n".format(MovingTargetExtended),
                            "      MovingTargetConvolveExtended: {}\n".format(MovingTargetConvolveExtended),
                            "      MovingTargetToTrack: {}\n".format(MovingTargetToTrack),
                            "      BackgroundRate: {}\n\n".format(BackgroundRate_lw)
                        ]



                text_out += text * xml_dict['PrimaryDithers'][exposure_index]

            elif instrument in ['NIRISS', 'FGS', 'NIRSPEC', 'MIRI']:
                text += [
                    "  Filter: {}\n".format(xml_dict['FilterWheel'][exposure_index]),
                    "  PointSourceCatalog: {}\n".format(catalog_files[exposure_index]),
                    "  GalaxyCatalog: {}\n".format(GalaxyCatalog),
                    "  ExtendedCatalog: {}\n".format(ExtendedCatalog),
                    "  ExtendedScale: {}\n".format(ExtendedScale),
                    "  ExtendedCenter: {}\n".format(ExtendedCenter),
                    "  MovingTargetList: {}\n".format(MovingTargetList),
                    "  MovingTargetSersic: {}\n".format(MovingTargetSersic),
                    "  MovingTargetExtended: {}\n".format(MovingTargetExtended),
                    "  MovingTargetConvolveExtended: {}\n".format(MovingTargetConvolveExtended),
                    "  MovingTargetToTrack: {}\n".format(MovingTargetToTrack),
                    "  BackgroundRate: {}\n".format(BackgroundRate),
                ]
            # multiply by the number of dithers
                text_out += text * xml_dict['ImageDithers'][exposure_index]
    #         1/0
    # Table(xml_table).pprint()
    # 1/0

    # if len(used_instruments) > 1:
    #     # Several instruments used within one APT file, e.g. parallels
    #     raise NotImplementedError
    # else:
    #     if used_instruments[0] in ['NIRCAM', 'WFSC']:
    #         # Define "default" values (probably should be changed eventually)
    #         date = '2019-07-04'
    #         PAV3 = '0.'
    #         GalaxyCatalog = 'None'
    #         ExtendedCatalog = 'None'
    #         ExtendedScale = '1.0'
    #         ExtendedCenter = '1024,1024'
    #         MovingTargetList = 'None'
    #         MovingTargetSersic = 'None'
    #         MovingTargetExtended = 'None'
    #         MovingTargetConvolveExtended = 'True'
    #         MovingTargetToTrack = 'None'
    #         BackgroundRate_sw = '0.5'
    #         BackgroundRate_lw = '1.2'
    #
    #         sw_filters = {}
    #         lw_filters = {}
    #         sw_filters_all = np.array(xml_table['ShortFilter'])
    #         lw_filters_all = np.array(xml_table['LongFilter'])
    #         tile_nums = xml_table['TileNumber']
    #         observation_ids = xml_table['ObservationID']
    #
    #         for i_obs_all in set(observation_ids):
    #             # i_obs_all = int(i_obs_all)
    #             current_obs_indices = [i == i_obs_all for i in observation_ids]
    #             if len(set(np.array(sw_filters_all)[current_obs_indices])) > 1:
    #                 print('Note: Multiple filters in observation {}'.format(i_obs_all))
    #                 # At some point could use the tile_nums to fix this
    #             sw_filters[i_obs_all] = sw_filters_all[current_obs_indices]
    #             lw_filters[i_obs_all] = lw_filters_all[current_obs_indices]
    #
    #         # Check that all parameters have the right length
    #         all_param_lengths = [len(ps_cat_sw), len(ps_cat_lw), len(sw_filters),
    #                              len(lw_filters), len(observations), len(i_observations),
    #                              len(obs_names)]
    #
    #     elif used_instruments[0] in ['NIRISS']:
    #         # Define "default" values (probably should be changed eventually)
    #         date = '2019-07-04'
    #         PAV3 = '0.'
    #         GalaxyCatalog = 'None'
    #         ExtendedCatalog = 'None'
    #         ExtendedScale = '1.0'
    #         ExtendedCenter = '1024,1024'
    #         MovingTargetList = 'None'
    #         MovingTargetSersic = 'None'
    #         MovingTargetExtended = 'None'
    #         MovingTargetConvolveExtended = 'True'
    #         MovingTargetToTrack = 'None'
    #         BackgroundRate = '0.5'
    #
    #         filters = {}
    #         filters_all = np.array(xml_table['PupilWheel'])
    #         # lw_filters_all = np.array(xml_table['LongFilter'])
    #         tile_nums = xml_table['TileNumber']
    #         observation_ids = xml_table['ObservationID']
    #
    #         for i_obs_all in set(observation_ids):
    #             # i_obs_all = int(i_obs_all)
    #             current_obs_indices = [i == i_obs_all for i in observation_ids]
    #             if len(set(np.array(filters_all)[current_obs_indices])) > 1:
    #                 print('Note: Multiple filters in observation {}'.format(i_obs_all))
    #                 # At some point could use the tile_nums to fix this
    #             filters[i_obs_all] = filters_all[current_obs_indices]
    #             # lw_filters[i_obs_all] = lw_filters_all[current_obs_indices]
    #
    #         # if only one catalog is provided, that catalog will be used for all observations
    #         if len(catalog_files) == 1:
    #             catalog_files = catalog_files*len(observations)
    #
    #
    #         # Check that all parameters have the right length
    #         all_param_lengths = [len(catalog_files), len(filters), len(observations),
    #                              len(i_observations), len(obs_names)]
    #
    # if len(set(all_param_lengths)) > 1:
    #     print(all_param_lengths)
    #     raise ValueError('Not all provided parameters have compatible '
    #                      'dimensions. Will not write {}'.format(yaml_file))
    #
    # write = ["# Observation list created by write_observationlist.py  in\n",
    #          "# mirage/yaml. Note: all values except filters and\n",
    #          "# observation names are default.\n\n"]
    #
    # if used_instruments[0] in ['NIRCAM', 'WFSC']:
    #     for i_obs in range(num_obs):
    #         obs_number = i_observations[i_obs] + 1
    #         write += [
    #             "Observation{}:\n".format(obs_number),
    #             "  Name: '{}'\n".format(obs_names[i_obs]),
    #             "  Date: {}\n".format(date),
    #             "  PAV3: {}\n".format(PAV3),
    #         ]
    #
    #         for i_filt, (sw_filt, lw_filt) in enumerate(zip(sw_filters[obs_number],
    #                                                         lw_filters[obs_number])):
    #             write += [
    #                 "  FilterConfig{}:\n".format(i_filt + 1),
    #                 "    SW:\n",
    #                 "      Filter: {}\n".format(sw_filt),
    #                 "      PointSourceCatalog: {}\n".format(ps_cat_sw[i_obs]),
    #                 "      GalaxyCatalog: {}\n".format(GalaxyCatalog),
    #                 "      ExtendedCatalog: {}\n".format(ExtendedCatalog),
    #                 "      ExtendedScale: {}\n".format(ExtendedScale),
    #                 "      ExtendedCenter: {}\n".format(ExtendedCenter),
    #                 "      MovingTargetList: {}\n".format(MovingTargetList),
    #                 "      MovingTargetSersic: {}\n".format(MovingTargetSersic),
    #                 "      MovingTargetExtended: {}\n".format(MovingTargetExtended),
    #                 "      MovingTargetConvolveExtended: {}\n".format(MovingTargetConvolveExtended),
    #                 "      MovingTargetToTrack: {}\n".format(MovingTargetToTrack),
    #                 "      BackgroundRate: {}\n".format(BackgroundRate_sw),
    #                 "    LW:\n",
    #                 "      Filter: {}\n".format(lw_filt),
    #                 "      PointSourceCatalog: {}\n".format(ps_cat_lw[i_obs]),
    #                 "      GalaxyCatalog: {}\n".format(GalaxyCatalog),
    #                 "      ExtendedCatalog: {}\n".format(ExtendedCatalog),
    #                 "      ExtendedScale: {}\n".format(ExtendedScale),
    #                 "      ExtendedCenter: {}\n".format(ExtendedCenter),
    #                 "      MovingTargetList: {}\n".format(MovingTargetList),
    #                 "      MovingTargetSersic: {}\n".format(MovingTargetSersic),
    #                 "      MovingTargetExtended: {}\n".format(MovingTargetExtended),
    #                 "      MovingTargetConvolveExtended: {}\n".format(MovingTargetConvolveExtended),
    #                 "      MovingTargetToTrack: {}\n".format(MovingTargetToTrack),
    #                 "      BackgroundRate: {}\n\n".format(BackgroundRate_lw)
    #             ]
    #
    # elif used_instruments[0] in ['NIRISS', 'FGS']:
    #     for i_obs in range(num_obs):
    #         obs_number = i_observations[i_obs] + 1
    #         write += [
    #             "Observation{}:\n".format(obs_number),
    #             "  Name: '{}'\n".format(obs_names[i_obs]),
    #             "  Date: {}\n".format(date),
    #             "  PAV3: {}\n".format(PAV3),
    #         ]
    #
    #         # for i_filt, filt in enumerate(filters[obs_number]):
    #         for i_filt, filt in enumerate(filters[i_obs]):
    #             write += [
    #                 # "  FilterConfig{}:\n".format(i_filt + 1),
    #                 # "    SW:\n",
    #                 "  Filter: {}\n".format(filt),
    #                 "  PointSourceCatalog: {}\n".format(catalog_files[i_obs]),
    #                 "  GalaxyCatalog: {}\n".format(GalaxyCatalog),
    #                 "  ExtendedCatalog: {}\n".format(ExtendedCatalog),
    #                 "  ExtendedScale: {}\n".format(ExtendedScale),
    #                 "  ExtendedCenter: {}\n".format(ExtendedCenter),
    #                 "  MovingTargetList: {}\n".format(MovingTargetList),
    #                 "  MovingTargetSersic: {}\n".format(MovingTargetSersic),
    #                 "  MovingTargetExtended: {}\n".format(MovingTargetExtended),
    #                 "  MovingTargetConvolveExtended: {}\n".format(MovingTargetConvolveExtended),
    #                 "  MovingTargetToTrack: {}\n".format(MovingTargetToTrack),
    #                 "  BackgroundRate: {}\n".format(BackgroundRate),
    #             ]

    f = open(yaml_file, 'w')
    for line in text_out:
        f.write(line)
    f.close()
    print('\nSuccessfully wrote {} exposures to {}'.format(number_of_exposures, yaml_file))

    return xml_dict

# if __name__ == '__main__':
#     xml_file = '../OTECommissioning/OTE01/OTE01-1134.xml'
#     pointing_file = '../OTECommissioning/OTE01/OTE01-1134.pointing'
#     yaml_file = 'test.yaml'
#
#     write_yaml(xml_file, pointing_file, yaml_file)
