#! /usr/bin/env python
"""Test the interface between mirage and APT files for automatic observation list generation.

Authors
-------
    - Johannes Sahlmann

Use
---
    >>> pytest -s test_apt_interface.py


"""

import glob
import os
import pytest
import shutil

from astropy.io import ascii
import numpy as np

from mirage.yaml import generate_observationlist, yaml_generator
from mirage.apt.read_apt_xml import ReadAPTXML
# for debugging
# from mirage.apt import read_apt_xml, apt_inputs
# from mirage.utils import siaf_interface
# import importlib
# importlib.reload( yaml_generator )
# importlib.reload( generate_observationlist )
# importlib.reload( read_apt_xml )
# importlib.reload( apt_inputs )
# importlib.reload( siaf_interface )
# from mirage.yaml import generate_observationlist, yaml_generator

TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), 'test_data')
TEMPORARY_DIR = os.path.join(os.path.dirname(__file__), 'temp_data')

# @pytest.fixture(scope="module")
def temporary_directory(test_dir=TEMPORARY_DIR):
    """Create a test directory for permission management.

    Parameters
    ----------
    test_dir : str
        Path to directory used for testing

    Yields
    -------
    test_dir : str
        Path to directory used for testing
    """
    if os.path.isdir(test_dir):
        shutil.rmtree(test_dir)
        os.mkdir(test_dir)  # creates directory with default mode=511
    else:
        os.mkdir(test_dir)


def test_observation_list_generation_minimal():

    # generate output directory
    temporary_directory()

    instrument = 'NIRISS'

    catalogs = {}
    if instrument == 'NIRISS':
        apt_dir = os.path.join(TEST_DATA_DIR, instrument)
        apt_file_seed = '1087_minimal'
        source_list_file_name = os.path.join(apt_dir, 'niriss_point_sources.list')
        catalogs[instrument.lower()] = source_list_file_name

    # Write observationlist.yaml
    observation_list_file = os.path.join(TEMPORARY_DIR, '{}_observation_list.yaml'.format(instrument.lower()))
    apt_file_xml = os.path.join(apt_dir, '{}.xml'.format(apt_file_seed))
    outputs = generate_observationlist.get_observation_dict(apt_file_xml, observation_list_file, catalogs)

    assert os.path.isfile(observation_list_file)


def test_complete_input_generation():
    """Exercise mirage input generation from APT files (.xml and .pointing)."""

    # generate output directory
    temporary_directory()


    for instrument in ['NIRCam', 'NIRISS', 'NIRSpec', 'MIRI', 'misc', 'FGS']:
    # for instrument in ['NIRISS', 'NIRSpec', 'MIRI', 'FGS', 'NIRCam']:
    # for instrument in ['NIRISS', 'NIRSpec', 'MIRI', 'FGS']:
    # for instrument in ['NIRISS']:
    # for instrument in ['misc']:
    # for instrument in ['NIRSpec']:
    # for instrument in ['MIRI']:
    # for instrument in ['FGS']:
    # for instrument in ['NIRCam']:


        apt_dir = os.path.join(TEST_DATA_DIR, instrument)
        if instrument == 'NIRISS':
            apt_file_seeds = ['1087_minimal', '1088', '1087', 'm31_field_test_observation']
            # apt_file_seeds = ['1087']
            # apt_file_seeds = ['1087_minimal']
            source_list_file_name = os.path.join(apt_dir, 'niriss_point_sources.list')
        elif instrument == 'NIRCam':
            apt_file_seeds = ['1069', '1144-OTE-10', 'NIRCamTest']
            # apt_file_seeds = ['NIRCamTest']
            # apt_file_seeds = ['1069']
            # apt_file_seeds = ['1144-OTE-10']
            source_list_file_name = os.path.join(apt_dir, 'seed_im_from_catalog_test_ptsrc_catalog.list')
        elif instrument == 'NIRSpec':
            apt_file_seeds = ['1164']
            source_list_file_name = os.path.join(apt_dir, 'niriss_point_sources.list')
        elif instrument == 'MIRI':
            apt_file_seeds = ['1171']
            source_list_file_name = os.path.join(apt_dir, 'niriss_point_sources.list')
        elif instrument == 'FGS':
            apt_file_seeds = ['MIRAGE-test']
            source_list_file_name = os.path.join(apt_dir, 'niriss_point_sources.list')

        elif instrument == 'misc':
            apt_xml_files = glob.glob(os.path.join(apt_dir, '*/*.xml'))
            apt_file_seeds = [f.split(apt_dir)[1] for f in apt_xml_files]
            source_list_file_name = os.path.join(apt_dir, 'niriss_point_sources.list')

        catalogs = {}
        for instrument_name in 'fgs nircam niriss miri nirspec'.split():
            if instrument_name.lower() == 'nircam':
                catalogs[instrument_name.lower()] = {}
                catalogs[instrument_name.lower()]['sw'] = source_list_file_name
                catalogs[instrument_name.lower()]['lw'] = source_list_file_name
            else:
                catalogs[instrument_name.lower()] = source_list_file_name

        for i, apt_file_seed in enumerate(apt_file_seeds):
            print('\n\n' + '=' * 100 + '\n')

            # For the moment, skip tests that contain NirissImaging observations
            # or that take too long for Travis
            skip_for_now = ['DeepField', 'NCam010', '1071']
            skip_bool = any([True if prop in apt_file_seed else False for prop in skip_for_now])
            if skip_bool:
                continue
                
            obs_yaml_files = glob.glob(os.path.join(TEMPORARY_DIR, 'jw*.yaml'))
            for file in obs_yaml_files:
                os.remove(file)

            if '.xml' in apt_file_seed:
                apt_file_xml = os.path.join(apt_dir, apt_file_seed[1:])
                apt_file_pointing = os.path.join(apt_dir, apt_file_seed[1:].replace('.xml', '.pointing'))
                observation_list_file = os.path.join(TEMPORARY_DIR, '{}_observation_list.yaml'.format(apt_file_seed.replace('/', '_').split('.')[0]))

            else:
                observation_list_file = os.path.join(TEMPORARY_DIR, '{}_observation_list.yaml'.format(apt_file_seed))
                apt_file_xml = os.path.join(apt_dir, '{}.xml'.format(apt_file_seed))
                apt_file_pointing = os.path.join(apt_dir, '{}.pointing'.format(apt_file_seed))

            print('Processing program {}'.format(apt_file_xml))

            yam = yaml_generator.SimInput(input_xml=apt_file_xml, pointing_file=apt_file_pointing,
                                          catalogs=catalogs, observation_list_file=observation_list_file,
                                          verbose=True, output_dir=TEMPORARY_DIR, simdata_output_dir=TEMPORARY_DIR,
                                          offline=True)

            try:
                yam.create_inputs()
            except RuntimeError as e:
                print('\nERROR detected. Skipping {} because of error\n{}\n\n'.format(apt_file_xml, e))
                continue

            yfiles = glob.glob(os.path.join(yam.output_dir,'jw*{}*.yaml'.format(yam.info['ProposalID'][0])))
            valid_instrument_list = [s for s in yam.info['Instrument'] if s.lower() in 'fgs nircam niriss'.split()]
            assert len(valid_instrument_list) == len(yfiles)


def read_xml(xml_file):
    """Read in and parse xml file from APT

    Parameters
    ----------
    xml_file : str
        XML file from APT

    Returns
    -------
    exposure_dictionary : dict
        Dictionary containing exposure information for observations in xml_file
    """
    xml = ReadAPTXML()
    exposure_dict = xml.read_xml(xml_file)
    return exposure_dict


def test_xml_reader():
    """Tests for the xml reader in read_apt_xml.py"""

    programs_to_test = ['08888', '12345', '54321']

    for program in programs_to_test:
        apt_dir = os.path.join(TEST_DATA_DIR, 'misc/{}'.format(program))
        xml_filename = glob.glob(os.path.join(apt_dir, '*.xml'))[0]
        comp_file = xml_filename.replace('.xml', '.txt')

        # Create a new exposure dictionary. Compare to the truth version
        exposure_dict = read_xml(xml_filename)
        comparison_dict = ascii.read(comp_file)

        # Columns to convert to strings since astropy defaults to reading them in as integers
        int_to_string_cols = ['ProposalID', 'PrimaryDithers', 'SubpixelPositions', 'ObservationID',
                              'number_of_dithers', 'Groups', 'Integrations', 'TileNumber']
        # Boolean columns to convert to strings in order for the comparison to work
        bool_to_string_cols = ['ParallelInstrument']

        # Convert integer columns to string
        for col in int_to_string_cols:
            data = comparison_dict[col].data
            if col == 'ProposalID':
                test = str(data[0])
                leading_zeros = 0
                if len(test) < 5:
                    leading_zeros = 5 - len(test)
                data = ['{}{}'.format('0'*leading_zeros, entry) for entry in data]
            else:
                data = [str(entry) for entry in data]
            comparison_dict[col] = data
        # Convert boolean columns to string
        for col in bool_to_string_cols:
            data = exposure_dict[col]
            data = [str(entry) for entry in data]
            exposure_dict[col] = data

        # Check table lengths
        assert len(exposure_dict['Instrument']) == len(comparison_dict['Instrument'])

        # Compare values in each column (skip ETC-related columns)
        for col in comparison_dict.colnames:
            if col[0:3] != 'Etc':
                data = comparison_dict[col].data
                if isinstance(data[0], np.int64):
                    data = [str(d) for d in data]
                    comparison_dict[col] = data
                if isinstance(data[0], int):
                    assert all(exposure_dict[col] == comparison_dict[col].data), print(col,
                                                                                       exposure_dict[col],
                                                                                       comparison_dict[col].data)
                else:
                    assert all(exposure_dict[col] == comparison_dict[col].data), print(col,
                                                                                       exposure_dict[col],
                                                                                       comparison_dict[col].data)


# for debugging
if __name__ == '__main__':
    test_complete_input_generation()
    #test_xml_reader()
