'''Define unit tests for parsing APT templates with pytest.

Authors
-------
    - Lauren Chambers

Use
---
    Ensure you have pytest installed. Then, simply run pytest in any
    parent directory of nircam_simulator/tests/:
    >>> pytest
'''

import os

from nircam_simulator.scripts import write_observationlist, yaml_generator, utils
from nircam_simulator.scripts.get_catalog import get_all_catalogs

INSTRUMENTS = ['NIRCam']
PROPOSAL_ID = '1111'

TESTS_DIR = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


def test_APT():
    '''Parse the given APT files and create a set of .yamls for each
    instrument in the INSTRUMENTS list.
    '''
    for instrument in INSTRUMENTS:
        # Define .pointing and .xml file locations
        pointing_file = os.path.join(TESTS_DIR, 'test_data', instrument + 'Test.pointing')
        xml_file = os.path.join(TESTS_DIR, 'test_data', instrument + 'Test.xml')

        # Point to appropriate output directory
        out_dir = os.path.join(TESTS_DIR, 'test_data', 'APT_{}_out'.format(instrument))

        # Create/locate catalogs for target(s)
        sw_cats, lw_cats = get_all_catalogs(pointing_file, PROPOSAL_ID)

        # Write observationlist.yaml
        observationlist_file = os.path.join(out_dir, instrument + '_observationlist.yaml')
        write_observationlist.write_yaml(xml_file, pointing_file, observationlist_file,
                                         ps_cat_sw=sw_cats, ps_cat_lw=lw_cats)

        # Create a series of data simulator input yaml files
        yam = yaml_generator.SimInput()
        yam.input_xml = xml_file
        yam.pointing_file = pointing_file
        yam.siaf = utils.get_siaf()['NIRCam']
        yam.output_dir = out_dir
        yam.simdata_output_dir = out_dir
        yam.observation_table = observationlist_file
        yam.use_JWST_pipeline = False
        yam.use_linearized_darks = True
        yam.datatype = 'linear'
        yam.reffile_setup()
        yam.create_inputs()
