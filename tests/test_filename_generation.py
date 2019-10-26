#! /usr/bin/env python
"""Test the file name generation in terms of formatting.

Authors
-------
    - Johannes Sahlmann

Use
---
    >>> pytest -s test_filename_generation.py


"""

import glob
import os

from astropy.table import Table
import numpy as np
import pytest

from mirage.yaml import generate_observationlist, yaml_generator
from .test_apt_interface import ON_TRAVIS, temporary_directory, TEST_DATA_DIR

@pytest.mark.skipif(ON_TRAVIS,
                   reason="Cannot access mirage data in the central storage directory from Travis CI.")
def test_filenames(temporary_directory):
    """Exercise mirage input generation from APT files (.xml and .pointing)."""

    for instrument in ['misc']:
        apt_dir = os.path.join(TEST_DATA_DIR, instrument)

        if instrument == 'misc':
            apt_file_seeds = ['12345/12345_nircam_imaging_prime_niriss_wfss_parallel.xml']
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

            obs_yaml_files = glob.glob(os.path.join(temporary_directory, 'jw*.yaml'))
            for file in obs_yaml_files:
                os.remove(file)

            if '.xml' in apt_file_seed:
                apt_file_xml = os.path.join(apt_dir, apt_file_seed)
                apt_file_pointing = os.path.join(apt_dir, apt_file_seed.replace('.xml', '.pointing'))
                observation_list_file = os.path.join(temporary_directory,
                                                     '{}_observation_list.yaml'.format(apt_file_seed.replace('/', '_').split('.')[0]))


            print('Processing program {}'.format(apt_file_xml))

            yam = yaml_generator.SimInput(input_xml=apt_file_xml, pointing_file=apt_file_pointing,
                                          observation_list_file=observation_list_file,
                                          verbose=True, output_dir=temporary_directory, simdata_output_dir=temporary_directory,
                                          offline=True)
            try:
                yam.create_inputs()
            except RuntimeError as e:
                print('\nERROR detected. Skipping {} because of error\n{}\n\n'.format(apt_file_xml, e))
                continue


            Table([yam.info['yamlfile']]).pprint(max_lines=-1)

            exposure_numbers = np.array([np.int(s.split('_')[2]) for s in yam.info['yamlfile']])
            assert np.max(exposure_numbers) <= 18
