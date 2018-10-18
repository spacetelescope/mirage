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
import numpy as np

from ..apt import read_apt_xml


def write_yaml(xml_file, yaml_file, catalog_files=None, ps_cat_sw=None, ps_cat_lw=None,
               verbose=False):
    """Write observation list file (required mirage input) on the basis of APT files.

    Parameters
    ----------
    xml_file : str
        path to APT .xml file
    yaml_file : str
        output_file
    catalog_files : str or list(str)
        path to file(s) that contain the catalog of sources.
    ps_cat_sw : str or list(str)
        NIRCam SW catalog, one per observations
    ps_cat_lw : str or list(str)
        NIRCam LW catalog, one per observations

    Returns
    -------
    xml_dict : dict
        Expanded dictionary that holds exposure information

    """
    if (catalog_files is not None) and (type(catalog_files) is str):
        catalog_files = [catalog_files]

    # Read in filters from APT .xml file
    readxml_obj = read_apt_xml.ReadAPTXML()

    xml_dict = readxml_obj.read_xml(xml_file, verbose=verbose)
    if verbose:
        print('Summary of observation dictionary:')
        for key in xml_dict.keys():
            print('{:<25}: number of elements is {:>5}'.format(key, len(xml_dict[key])))

    # array of unique instrument names
    used_instruments = np.unique(xml_dict['Instrument'])
    unique_observation_ids = np.unique(xml_dict['ObservationID']).tolist()
    number_of_observations = len(unique_observation_ids)
    number_of_exposures = len(xml_dict['ObservationID'])

    # temporary fix for NIRCam
    if ('NIRCAM' in used_instruments) and (ps_cat_sw is None):
        ps_cat_sw = catalog_files * number_of_observations
    if ('NIRCAM' in used_instruments) and (ps_cat_lw is None):
        ps_cat_lw = catalog_files * number_of_observations

    # non-NIRCam: if only one catalog is provided, that catalog will be used for all observations
    if (catalog_files is not None) and (len(catalog_files) == 1):
        catalog_files = catalog_files * number_of_exposures

    if verbose:
        print('Summary of dictionary extracted from {}'.format(xml_file))
        for key in xml_dict.keys():
            print('{:<25}: number of elements is {:>5}'.format(key, len(xml_dict[key])))

    # set default values (should be changed eventually)
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

    # assemble string that will constitute the yaml content
    text_out = ["# Observation list created by write_observationlist.py\n\n"]

    # loop over instruments and exposures (allows for multiple instruments and parallels)
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

            # 'ImageDithers' and 'PrimaryDithers' are either 0 (when not set), or 1 (default), or N (when set)
            number_of_dithers = np.max([np.int(xml_dict['ImageDithers'][exposure_index]), np.int(xml_dict['PrimaryDithers'][exposure_index])])
            # if np.int(xml_dict['PrimaryDithers'][exposure_index]) == 0:
            #     number_of_dithers = np.int(xml_dict['ImageDithers'][exposure_index])
            # elif np.int(xml_dict['ImageDithers'][exposure_index]) == 0:
            #     number_of_dithers = np.int(xml_dict['PrimaryDithers'][exposure_index])
            # else:
            #     raise ValueError('Both PrimaryDithers and ImageDithers are not zero.')

            if instrument in ['NIRCAM', 'WFSC']:

                # set default values
                BackgroundRate_sw = '0.5'
                BackgroundRate_lw = '1.2'

                sw_filt = xml_dict['ShortFilter'][exposure_index]
                lw_filt = xml_dict['LongFilter'][exposure_index]
                # i_obs = xml_dict['ObservationID'][exposure_index]-1
                unique_observation_index = unique_observation_ids.index(xml_dict['ObservationID'][exposure_index])

                text += [
                    # "  FilterConfig{}:\n".format(i_obs+1),
                    "  FilterConfig{}:\n".format(exposure_index),
                    "    SW:\n",
                    "      Filter: {}\n".format(sw_filt),
                    "      PointSourceCatalog: {}\n".format(ps_cat_sw[unique_observation_index]),
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
                    "      PointSourceCatalog: {}\n".format(ps_cat_lw[unique_observation_index]),
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
                # multiply by the number of dithers
                text_out += text * number_of_dithers

            elif instrument in ['NIRISS', 'FGS', 'NIRSPEC', 'MIRI']:
                # set default values
                BackgroundRate = '0.5'

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
                text_out += text * number_of_dithers

    f = open(yaml_file, 'w')
    for line in text_out:
        f.write(line)
    f.close()
    print('\nSuccessfully wrote {} exposures to {}'.format(number_of_exposures, yaml_file))

    return xml_dict
