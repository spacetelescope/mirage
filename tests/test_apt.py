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
import glob
import shutil

from lxml import etree
import pytest

from mirage.yaml import generate_observationlist, yaml_generator

# for debugging
# from mirage.apt import read_apt_xml, apt_inputs
# from mirage.utils import siaf_interface
# import importlib
# importlib.reload( yaml_generator )
# importlib.reload( write_observationlist )
# importlib.reload( read_apt_xml )
# importlib.reload( apt_inputs )
# importlib.reload( siaf_interface )



# to be set before running test
# os.environ['MIRAGE_DATA'] = ''

PROPOSAL_ID = '1111'
APT_NAMESPACE = '{http://www.stsci.edu/JWST/APT}'

TESTS_DIR = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

# Determine if tests are being run on Github Actions CI
ON_GITHUB = '/home/runner' in os.path.expanduser('~')
if not ON_GITHUB:
    orig_mirage_data = os.environ['MIRAGE_DATA']
os.environ['MIRAGE_DATA'] = '/test/'

# @pytest.mark.xfail
def RunAllAPTTemplates(instrument):
    '''Parse the given APT files and create a set of .yamls for a given
    instrument
    '''
    # Define .pointing and .xml file locations
    pointing_file = os.path.join(TESTS_DIR, 'test_data', instrument, instrument + 'Test.pointing')
    xml_file = os.path.join(TESTS_DIR, 'test_data',  instrument, instrument + 'Test.xml')

    # Open XML file, get element tree of the APT proposal to determine how
    # many observations there are
    with open(xml_file) as f:
        tree = etree.parse(f)
    observation_data = tree.find(APT_NAMESPACE + 'DataRequests')
    obs_results = observation_data.findall('.//' + APT_NAMESPACE + 'Observation')
    n_obs = len(obs_results)

    # Locate catalogs for target(s) (one catalog per observation and channel)
    #sw_cats = [os.path.join(TESTS_DIR, 'test_data', '2MASS_RA273.09deg_Dec65.60deg.list')] * n_obs
    #lw_cats = [os.path.join(TESTS_DIR, 'test_data', 'WISE_RA273.09deg_Dec65.60deg.list')] * n_obs
    #cat_dict = {'nircam': {'lw': lw_cats,
    #                       'sw': sw_cats}}
    cat_dict = None

    # Point to appropriate output directory
    out_dir = os.path.join(TESTS_DIR, 'test_data',  instrument, 'APT_{}_out'.format(instrument))

    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
        os.makedirs(out_dir)
    else:
        os.makedirs(out_dir)

    # Create a series of data simulator input yaml files
    yam = yaml_generator.SimInput(input_xml=xml_file, pointing_file=pointing_file,
                                  catalogs=cat_dict, use_JWST_pipeline=False,
                                  verbose=True, output_dir=out_dir, simdata_output_dir=out_dir,
                                  offline=True)
    yam.create_inputs()

    # Ensure that some of the expected files have been created
    assert os.path.exists(os.path.join(out_dir, 'Observation_table_for_' +
                                                instrument +
                                                'Test.xml_with_yaml_parameters.csv')), \
        'Observation table not created.'

    number_of_yaml_files  = len(glob.glob(os.path.join(out_dir, 'jw{:05d}*.yaml'.format(int(PROPOSAL_ID)))))
    print('PROPOSAL_ID: {}'.format(PROPOSAL_ID))
    print('number of observations: {}'.format(n_obs))
    print('number of files written: {}'.format(number_of_yaml_files))
    assert n_obs == 17
    # assert number_of_yaml_files == 150
    assert number_of_yaml_files >= n_obs, 'Fewer yaml files created than observations'

    # If a reference observationlist.yaml file exists, ensure that the
    # file that was just created matches it
    # NOT USED because out_dir is being recreated at runtime
    # reference_yaml = os.path.join(out_dir, 'REFERENCE_' + instrument + '_observationlist.yaml')
    # print(reference_yaml)
    # if os.path.exists(reference_yaml):
    #     assert yaml.load(reference_yaml) == yaml.load(observationlist_file),\
    #         'The created observationlist.yaml file does not match the reference' +\
    #         'observationlist.yaml file. Either the APT parser is malfunctioning,' +\
    #         'or the reference yaml is out of date.'


@pytest.mark.skipif(ON_GITHUB,
                    reason="Cannot access mirage data in the central storage directory from Github CI.")
def test_environment_variable():
    '''Ensure the MIRAGE_DATA environment variable has been set
    '''
    failure_msg = \
        "MIRAGE_DATA environment variable is not set. This must be set to the " \
        "base directory containing the darks, cosmic ray, PSF, etc input files " \
        "needed for the simulation. These files must be downloaded separately " \
        "from the Mirage package."

    try:
        env_var = os.environ['MIRAGE_DATA']
        assert env_var != '', failure_msg
    except KeyError:
        assert False, failure_msg


def test_RunNIRCamAPTTemplates():
    '''Parse the given APT files and create a set of .yamls for NIRCam
    '''
    RunAllAPTTemplates('NIRCam')


# Return environment variable to original value. This is helpful when
# calling many tests at once, some of which need the real value.
if not ON_GITHUB:
    os.environ['MIRAGE_DATA'] = orig_mirage_data

# for debugging
# if __name__ == '__main__':
#     RunAllAPTTemplates('NIRCam')
