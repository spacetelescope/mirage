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

from mirage.yaml import generate_observationlist, yaml_generator
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

@pytest.fixture()
def temporary_directory(test_dir=TEMPORARY_DIR):
    """Create a test directory.

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

    yield test_dir

    print("teardown tests/temp_data directory")
    if os.path.isdir(test_dir):
        shutil.rmtree(test_dir)


def test_observation_list_generation_minimal(temporary_directory):
    instrument = 'NIRISS'

    catalogs = {}
    if instrument == 'NIRISS':
        apt_dir = os.path.join(TEST_DATA_DIR, instrument)
        apt_file_seed = '1087_minimal'
        source_list_file_name = os.path.join(apt_dir, 'niriss_point_sources.list')
        catalogs[instrument.lower()] = source_list_file_name

    # Write observationlist.yaml
    observation_list_file = os.path.join(temporary_directory, '{}_observation_list.yaml'.format(instrument.lower()))
    apt_file_xml = os.path.join(apt_dir, '{}.xml'.format(apt_file_seed))
    generate_observationlist.get_observation_dict(apt_file_xml, observation_list_file, catalogs=catalogs)

    assert os.path.isfile(observation_list_file)

# Assemble list of programs to test
complete_input_generation_parameters = [('NIRISS', '1087_minimal', 'niriss_point_sources.list'),
                                        ('NIRISS', '1088', 'niriss_point_sources.list'),
                                        ('NIRISS', '1087', 'niriss_point_sources.list'),
                                        ('NIRISS', 'm31_field_test_observation', 'niriss_point_sources.list'),
                                        ('NIRCam', '1069', 'seed_im_from_catalog_test_ptsrc_catalog.list'),
                                        ('NIRCam', '1144-OTE-10', 'seed_im_from_catalog_test_ptsrc_catalog.list'),
                                        ('NIRCam', 'NIRCamTest', 'seed_im_from_catalog_test_ptsrc_catalog.list'),
                                        ('NIRSpec', '1164', 'niriss_point_sources.list'),
                                        ('MIRI', '1171', 'niriss_point_sources.list'),
                                        ('FGS', 'MIRAGE-test', 'niriss_point_sources.list')]
# Add miscellaneous programs
apt_xml_files = glob.glob(os.path.join(TEST_DATA_DIR, 'misc', '*/*.xml'))
for file in apt_xml_files:
    apt_dir = os.path.join(TEST_DATA_DIR, 'misc')
    seed = file.split(apt_dir)[1]
    complete_input_generation_parameters.append(('misc', seed, 'niriss_point_sources.list'))

@pytest.mark.parametrize('instrument, apt_file_seed, source_list_file_name', complete_input_generation_parameters)
def test_complete_input_generation(temporary_directory, instrument, apt_file_seed, source_list_file_name):
    """Exercise mirage input generation from APT files (.xml and .pointing)."""

    # Determine source list filename
    apt_dir = os.path.join(TEST_DATA_DIR, instrument)
    source_list_file_name = os.path.join(apt_dir, source_list_file_name)

    # Populate the catalogs
    catalogs = {}
    for instrument_name in 'fgs nircam niriss miri nirspec'.split():
        if instrument_name.lower() == 'nircam':
            catalogs[instrument_name.lower()] = {}
            catalogs[instrument_name.lower()]['sw'] = source_list_file_name
            catalogs[instrument_name.lower()]['lw'] = source_list_file_name
        else:
            catalogs[instrument_name.lower()] = source_list_file_name

    print('\n\n' + '=' * 100 + '\n')
    if 'DeepField' in apt_file_seed: # exlude until support for NircamTimeSeries is included
        return

    # # Remove any existing yaml files
    # obs_yaml_files = glob.glob(os.path.join(temporary_directory, 'jw*.yaml'))
    # for file in obs_yaml_files:
    #     os.remove(file)

    if '.xml' in apt_file_seed:
        apt_file_xml = os.path.join(apt_dir, apt_file_seed[1:])
        apt_file_pointing = os.path.join(apt_dir, apt_file_seed[1:].replace('.xml', '.pointing'))
        observation_list_file = os.path.join(temporary_directory, '{}_observation_list.yaml'.format(apt_file_seed.replace('/', '_').split('.')[0]))

    else:
        observation_list_file = os.path.join(temporary_directory, '{}_observation_list.yaml'.format(apt_file_seed))
        apt_file_xml = os.path.join(apt_dir, '{}.xml'.format(apt_file_seed))
        apt_file_pointing = os.path.join(apt_dir, '{}.pointing'.format(apt_file_seed))

    print('Processing program {}'.format(apt_file_xml))
    yam = yaml_generator.SimInput(input_xml=apt_file_xml, pointing_file=apt_file_pointing,
                                  catalogs=catalogs, observation_list_file=observation_list_file,
                                  verbose=True, output_dir=temporary_directory, simdata_output_dir=temporary_directory,
                                  offline=True)

    yam.create_inputs()

    yfiles = glob.glob(os.path.join(yam.output_dir,'jw*{}*.yaml'.format(yam.info['ProposalID'][0])))
    valid_instrument_list = [s for s in yam.info['Instrument'] if s.lower() in 'fgs nircam niriss'.split()]
    assert len(valid_instrument_list) == len(yfiles)

# for debugging
if __name__ == '__main__':
    test_complete_input_generation()
