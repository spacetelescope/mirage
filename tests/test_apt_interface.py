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
from mirage.utils.utils import ensure_dir_exists

TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), 'test_data')
TEMPORARY_DIR = os.path.join(os.path.dirname(__file__), 'temp_data')

# Determine if tests are being run on Travis
ON_TRAVIS = 'travis' in os.path.expanduser('~')

if not ON_TRAVIS:
    orig_mirage_data = os.environ['MIRAGE_DATA']
os.environ['MIRAGE_DATA'] = '/test/'


@pytest.fixture(scope="module")
def temporary_directory(test_dir=TEMPORARY_DIR):
    """Create a temporary directory for testing.

    Parameters
    ----------
    test_dir : str
        Path to directory used for testing

    Yields
    -------
    test_dir : str
        Path to directory used for testing
    """
    # Create directory and yield its name
    ensure_dir_exists(test_dir)  # creates directory with default mode=511
    yield test_dir

    # Remove directory
    if os.path.isdir(test_dir):
        shutil.rmtree(test_dir)


def test_observation_list_generation_minimal(temporary_directory):

    instrument = 'NIRISS'

    catalogs = {}
    if instrument == 'NIRISS':
        apt_dir = os.path.join(TEST_DATA_DIR, instrument)
        apt_file_seed = '1087_minimal'
        source_list_file_name = os.path.join(apt_dir, 'niriss_point_sources.list')
        catalogs = {}
        catalogs['point_source'] = source_list_file_name

    # Write observationlist.yaml
    observation_list_file = os.path.join(temporary_directory, '{}_observation_list.yaml'.format(instrument.lower()))
    apt_file_xml = os.path.join(apt_dir, '{}.xml'.format(apt_file_seed))
    outputs = generate_observationlist.get_observation_dict(apt_file_xml, observation_list_file, catalogs)

    assert os.path.isfile(observation_list_file)


@pytest.mark.skipif(ON_TRAVIS,
                   reason="Cannot access mirage data in the central storage directory from Travis CI.")
def test_complete_input_generation(temporary_directory):
    """Exercise mirage input generation from APT files (.xml and .pointing)."""

    for instrument in ['NIRCam', 'NIRISS', 'NIRSpec', 'MIRI', 'misc', 'FGS']:
        apt_dir = os.path.join(TEST_DATA_DIR, instrument)
        if instrument == 'NIRISS':
            apt_file_seeds = ['com1093', '1087_minimal', '1088', '1087', 'm31_field_test_observation']
            source_list_file_name = os.path.join(apt_dir, 'niriss_point_sources.list')
        elif instrument == 'NIRCam':
            apt_file_seeds = ['1069', '1144-OTE-10', 'NIRCamTest']
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

        for i, apt_file_seed in enumerate(apt_file_seeds):
            print('\n\n' + '=' * 100 + '\n')

            catalogs = None

            # For the moment, skip tests that contain NirissImaging observations
            # or that take too long for Travis
            skip_for_now = ['DeepField', 'NCam010', '1071']
            skip_bool = any([True if prop in apt_file_seed else False for prop in skip_for_now])
            if skip_bool:
                continue

            obs_yaml_files = glob.glob(os.path.join(temporary_directory, 'jw*.yaml'))
            for file in obs_yaml_files:
                os.remove(file)

            if '.xml' in apt_file_seed:
                apt_file_xml = os.path.join(apt_dir, apt_file_seed[1:])
                apt_file_pointing = os.path.join(apt_dir, apt_file_seed[1:].replace('.xml', '.pointing'))
            else:
                apt_file_xml = os.path.join(apt_dir, '{}.xml'.format(apt_file_seed))
                apt_file_pointing = os.path.join(apt_dir, '{}.pointing'.format(apt_file_seed))

            print('Processing program {}'.format(apt_file_xml))

            yam = yaml_generator.SimInput(input_xml=apt_file_xml, pointing_file=apt_file_pointing,
                                          catalogs=catalogs, verbose=True, output_dir=TEMPORARY_DIR,
                                          simdata_output_dir=TEMPORARY_DIR,
                                          offline=True)
            try:
                yam.create_inputs()
            except RuntimeError as e:
                print('\nERROR detected. Skipping {} because of error\n{}\n\n'.format(apt_file_xml, e))
                continue

            yfiles = glob.glob(os.path.join(yam.output_dir, 'jw*{}*.yaml'.format(yam.info['ProposalID'][0])))
            valid_instrument_list = [s for s in yam.info['Instrument'] if s.lower() in 'fgs nircam niriss'.split()]
            assert len(valid_instrument_list) == len(yfiles)

            if os.path.basename(apt_file_xml) in ['54321_niriss_wfss_prime_nircam_imaging_parallel.xml',
                                                  '12345_nircam_imaging_prime_niriss_wfss_parallel.xml']:
                prog = os.path.basename(apt_file_xml).split('_')[0]
                expected_nircam_files = {'54321': [60, 75], '12345': [40, 90, 40, 20]}
                expected_niriss_files = {'54321': [12, 15], '12345': [8, 18, 8, 0]}
                nircam_files_from_list = np.array([True if s.lower() == 'nircam' else False for s in valid_instrument_list])
                num_nircam_files = np.sum(nircam_files_from_list)
                num_niriss_files = len(yfiles) - num_nircam_files
                assert num_nircam_files == np.sum(np.array(expected_nircam_files[prog]))
                assert num_niriss_files == np.sum(np.array(expected_niriss_files[prog]))

                num_obs = len(expected_nircam_files[prog])
                for i in range(1, num_obs+1):
                    obs_num = str(i).zfill(3)
                    obs_yfiles_nircam = glob.glob(os.path.join(yam.output_dir, 'jw{}{}*nrc*.yaml'.format(prog, obs_num)))
                    obs_yfiles_niriss = glob.glob(os.path.join(yam.output_dir, 'jw{}{}*nis*.yaml'.format(prog, obs_num)))
                    assert len(obs_yfiles_nircam) == expected_nircam_files[prog][i-1]
                    assert len(obs_yfiles_niriss) == expected_niriss_files[prog][i-1]


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
        if program == '54321':
            xml_filename = '54321_niriss_wfss_prime_nircam_imaging_parallel.xml'
        elif program == '12345':
            xml_filename = '12345_nircam_imaging_prime_niriss_wfss_parallel.xml'
        elif program == '08888':
            xml_filename = '08888_niriss_nircam_wfss.xml'
        xml_filename = os.path.join(apt_dir, xml_filename)
        comp_file = xml_filename.replace('.xml', '.txt')

        # Create a new exposure dictionary. Compare to the truth version
        exposure_dict = read_xml(xml_filename)
        comparison_dict = ascii.read(comp_file)

        # Columns to convert to strings since astropy defaults to reading them in as integers
        int_to_string_cols = ['ProposalID', 'PrimaryDithers', 'SubpixelPositions', 'ObservationID',
                              'number_of_dithers', 'Groups', 'Integrations', 'TileNumber']
        # Boolean columns to convert to strings in order for the comparison to work
        bool_to_string_cols = ['ParallelInstrument', 'FiducialPointOverride']

        # Convert integer columns to string
        for col in int_to_string_cols:
            data = comparison_dict[col].data
            if col == 'ProposalID':
                data = [str(entry).zfill(5) for entry in data]
            elif col == 'ObservationID':
                data = [str(entry).zfill(3) for entry in data]
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
                assert all(exposure_dict[col] == comparison_dict[col].data), print(program, col,
                                                                                   exposure_dict[col],
                                                                                   comparison_dict[col].data)


# Return environment variable to original value. This is helpful when
# calling many tests at once, some of which need the real value.
if not ON_TRAVIS:
    os.environ['MIRAGE_DATA'] = orig_mirage_data
