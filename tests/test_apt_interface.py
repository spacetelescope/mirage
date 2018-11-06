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
# importlib.reload(generate_observationlist)
# importlib.reload( read_apt_xml )
# importlib.reload( apt_inputs )
# importlib.reload( siaf_interface )


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
    generate_observationlist.get_observation_dict(apt_file_xml, observation_list_file, catalogs=catalogs)

    assert os.path.isfile(observation_list_file)


def test_complete_input_generation():
    """Exercise mirage input generation from APT files (.xml and .pointing)."""

    # generate output directory
    temporary_directory()


    # for instrument in ['NIRCam', 'NIRISS', 'NIRSpec', 'MIRI', 'misc', 'FGS']:
    for instrument in ['NIRISS', 'NIRSpec', 'MIRI', 'FGS', 'NIRCam']:
    # for instrument in ['NIRISS', 'NIRSpec', 'MIRI', 'FGS']:
    # for instrument in ['NIRISS']:
    # for instrument in ['misc']:
    # for instrument in ['NIRSpec']:
    # for instrument in ['MIRI']:
    # for instrument in ['FGS']:
    # for instrument in ['NIRCam']:

        print('='*100)
        apt_dir = os.path.join(TEST_DATA_DIR, instrument)
        if instrument == 'NIRISS':
            apt_file_seeds = ['1088', '1087', 'm31_field_test_observation']
            # apt_file_seeds = ['1088']
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
                                          verbose=True, output_dir=TEMPORARY_DIR, simdata_output_dir=TEMPORARY_DIR)

            yam.reffile_setup(offline=True)
            yam.create_inputs()

            yfiles = glob.glob(os.path.join(yam.output_dir,'jw*{}*.yaml'.format(yam.info['ProposalID'][0])))
            valid_instrument_list = [s for s in yam.info['Instrument'] if s.lower() in 'fgs nircam niriss'.split()]
            assert len(valid_instrument_list) == len(yfiles)

# for debugging
# if __name__ == '__main__':
#     test_complete_input_generation()
