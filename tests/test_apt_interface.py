#! /usr/bin/env python
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
import shutil

from mirage.yaml import write_observationlist, yaml_generator
from mirage.apt import read_apt_xml, apt_inputs

import importlib
importlib.reload( yaml_generator )
importlib.reload( write_observationlist )
importlib.reload( read_apt_xml )
importlib.reload( apt_inputs )


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
    if os.path.isdir(test_dir):
        shutil.rmtree(test_dir)
        os.mkdir(test_dir)  # creates directory with default mode=511
    else:
        os.mkdir(test_dir)
    # os.chmod(test_dir, 777)
    # yield test_dir
    # print("teardown test directory")
    # if os.path.isdir(test_dir):
    #     os.remove(test_dir)



def est_observation_list_generation_minimal():

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
    write_observationlist.write_yaml(apt_file_xml, observation_list_file, source_list_file_name)

    assert os.path.isfile(source_list_file_name)


def test_complete_input_generation():
    """Exercise mirage input generation from APT files (.xml and .pointing)."""

    # generate output directory
    temporary_directory()

    # for instrument in ['NIRCam', 'NIRISS']:
    # for instrument in ['NIRISS']:
    for instrument in ['NIRSpec']:

        apt_dir = os.path.join(TEST_DATA_DIR, instrument)
        if instrument == 'NIRISS':
            apt_file_seeds = ['1088', '1087', 'm31_field_test_observation']
            source_list_file_name = os.path.join(apt_dir, 'niriss_point_sources.list')
        elif instrument == 'NIRCam':
            apt_file_seeds = ['NIRCamTest']
            source_list_file_name = os.path.join(apt_dir, 'seed_im_from_catalog_test_ptsrc_catalog.list')
        elif instrument == 'NIRSpec':
            apt_file_seeds = ['1164']
            source_list_file_name = os.path.join(apt_dir, 'niriss_point_sources.list')

        for apt_file_seed in apt_file_seeds:
            # Write observationlist.yaml
            observation_list_file = os.path.join(TEMPORARY_DIR, '{}_{}_observation_list.yaml'.format(instrument.lower(), apt_file_seed))

            apt_file_xml = os.path.join(apt_dir, '{}.xml'.format(apt_file_seed))
            apt_file_pointing = os.path.join(apt_dir, '{}.pointing'.format(apt_file_seed))

            apt_xml_dict = write_observationlist.write_yaml(apt_file_xml, observation_list_file, source_list_file_name, verbose=True)
            print(apt_xml_dict)
            yam = yaml_generator.SimInput()
            yam.input_xml = apt_file_xml
            yam.pointing_file = apt_file_pointing
            yam.output_dir = TEMPORARY_DIR
            yam.simdata_output_dir = TEMPORARY_DIR
            yam.observation_table = observation_list_file
            yam.use_JWST_pipeline = True
            yam.use_linearized_darks = False
            yam.datatype = 'linear'
            yam.reffile_setup(instrument=instrument, offline=True)
            yam.set_global_definitions()
            yam.create_inputs(apt_xml_dict=apt_xml_dict)

            yfiles = glob.glob(os.path.join(yam.output_dir,'jw*{}*.yaml'.format(yam.info['ProposalID'][0])))
            assert len(yam.info['Instrument']) == len(yfiles)

if __name__ == '__main__':
    test_complete_input_generation()
