"""Test the interface between mirage and APT files for automatic observation list generation.

Authors
-------
    - Johannes Sahlmann

Use
---
    >>> pytest -s test_niriss_imaging.py


"""

import glob
import os
import pytest

# from mirage import imaging_simulator as im
from mirage.yaml import write_observationlist, yaml_generator


# os.environ['MIRAGE_DATA'] = ''
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
    if not os.path.isdir(test_dir):
        os.mkdir(test_dir)  # creates directory with default mode=511
    # os.chmod(test_dir, 777)
    # yield test_dir
    # print("teardown test directory")
    # if os.path.isdir(test_dir):
    #     os.remove(test_dir)



def test_observation_list_generation_minimal():

    # generate output directory
    temporary_directory()

    instrument = 'NIRISS'

    if instrument == 'NIRISS':
        apt_dir = os.path.join(TEST_DATA_DIR, instrument)
        apt_file_seed = '1087_minimal'
        source_list_file_name = os.path.join(apt_dir, 'niriss_point_sources.list')

    # Write observationlist.yaml
    observation_list_file = os.path.join(TEMPORARY_DIR, '{}_observation_list.yaml'.format(instrument.lower()))

    apt_file_xml = os.path.join(apt_dir, '{}.xml'.format(apt_file_seed))
    # apt_file_pointing = os.path.join(apt_dir, '{}.pointing'.format(apt_file_seed))
    # apt_file_pointing,
    write_observationlist.write_yaml(apt_file_xml, observation_list_file, source_list_file_name)

    assert os.path.isfile(source_list_file_name)


def test_complete_input_generation():

    # generate output directory
    temporary_directory()

    instrument = 'NIRISS'

    if instrument == 'NIRISS':
        apt_dir = os.path.join(TEST_DATA_DIR, instrument)
        apt_file_seed = '1087'
        source_list_file_name = os.path.join(apt_dir, 'niriss_point_sources.list')

    # Write observationlist.yaml
    observation_list_file = os.path.join(TEMPORARY_DIR, '{}_observation_list.yaml'.format(instrument.lower()))

    apt_file_xml = os.path.join(apt_dir, '{}.xml'.format(apt_file_seed))
    apt_file_pointing = os.path.join(apt_dir, '{}.pointing'.format(apt_file_seed))

    apt_xml_dict = write_observationlist.write_yaml(apt_file_xml, observation_list_file, source_list_file_name)

    yam = yaml_generator.SimInput()
    yam.input_xml = os.path.join(apt_dir, '{}.xml'.format(apt_file_seed))
    yam.pointing_file = os.path.join(apt_dir, '{}.pointing'.format(apt_file_seed))
    yam.output_dir = TEMPORARY_DIR
    yam.simdata_output_dir = TEMPORARY_DIR
    yam.observation_table = observation_list_file
    yam.use_JWST_pipeline = True
    yam.use_linearized_darks = False
    yam.datatype = 'linear'
    yam.reffile_setup(instrument=instrument, offline=True) # only PRIME, ignores FGS
    yam.set_global_definitions()
    yam.create_inputs(apt_xml_dict=apt_xml_dict)

    yfiles = glob.glob(os.path.join(yam.output_dir,'jw*.yaml'))
    # print(len(yfiles))
    # print(len(apt_xml_dict['Instrument']))

    # assert os.path.isfile(source_list_file_name)
