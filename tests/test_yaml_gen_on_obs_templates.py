"""Test the yaml generator on a collection of APT files representing various supported
observation templates/combinations of templates

Authors
-------
    - Bryan Hilbert

Use
---
    >>> pytest -s test_yaml_gen_on_obs_templates.py
"""

from glob import glob
import numpy as np
import os
import pytest

from mirage.yaml.yaml_generator import SimInput

# Define directory and file locations
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


# Determine if tests are being run on Github Actions CI
ON_GITHUB = '/home/runner' in os.path.expanduser('~')

if not ON_GITHUB:
    orig_mirage_data = os.environ['MIRAGE_DATA']

# Reset the MIRAGE_DATA env variable to be a real location so yaml_generator
# doesn't croak
os.environ["MIRAGE_DATA"] = __location__
os.environ["CRDS_PATH"] = os.path.join(__location__, "temp")

# Define input data location
test_input_dir = os.path.join(os.path.dirname(__file__), 'test_data/obs_mode_tests')

def call_yaml_generator(xml_file):
    """Call the yaml generator. Assume the pointing file is in the same directory
    and has the same name as the xml file
    """
    temp_output_dir = os.path.join(__location__, "temp")
    pointing_file = xml_file.replace('.xml', '.pointing')

    yam = SimInput(xml_file, pointing_file, verbose=True,
                   offline=True, output_dir=temp_output_dir,
                   simdata_output_dir=temp_output_dir,
                   datatype='raw', reffile_defaults='crds'
                    )
    yam.use_linearized_darks = True
    yam.create_inputs()
    return yam


def test_fgs_external():
    xml = os.path.join(test_input_dir, 'fgs_external.xml')
    yamls = call_yaml_generator(xml)
    assert len(yamls.yaml_files) == 2

def test_miri_img_niriss_wfss():
    xml = os.path.join(test_input_dir, 'miri_img_niriss_wfss.xml')
    yamls = call_yaml_generator(xml)
    assert len(yamls.yaml_files) == 30

def test_miri_img_nrc_img():
    xml = os.path.join(test_input_dir, 'miri_img_nrc_img.xml')
    yamls = call_yaml_generator(xml)
    assert len(yamls.yaml_files) == 75

def test_nircam_img_niriss_wfss():
    xml = os.path.join(test_input_dir, 'nircam_img_niriss_wfss.xml')
    yamls = call_yaml_generator(xml)
    assert len(yamls.yaml_files) == 96
    nrc_yamls = np.array([True if 'nrc' in ele else False for ele in yamls.yaml_files])
    nis_yamls = np.array([True if 'nis' in ele else False for ele in yamls.yaml_files])
    assert np.sum(nrc_yamls) == 80
    assert np.sum(nis_yamls) == 16

def test_nircam_wfss_miri_img():
    xml = os.path.join(test_input_dir, 'nircam_wfss_miri_img.xml')
    yamls = call_yaml_generator(xml)
    assert len(yamls.yaml_files) == 140

def test_nircam_wfss_niriss_img():
    xml = os.path.join(test_input_dir, 'nircam_wfss_niriss_img.xml')
    yamls = call_yaml_generator(xml)
    assert len(yamls.yaml_files) == 336
    nrc_yamls = np.array([True if 'nrc' in ele else False for ele in yamls.yaml_files])
    nis_yamls = np.array([True if 'nis' in ele else False for ele in yamls.yaml_files])
    assert np.sum(nrc_yamls) == 280
    assert np.sum(nis_yamls) == 56

def test_niriss_img():
    xml = os.path.join(test_input_dir, 'niriss_img.xml')
    yamls = call_yaml_generator(xml)
    assert len(yamls.yaml_files) == 1

def test_niriss_wfss_miri_img():
    xml = os.path.join(test_input_dir, 'niriss_wfss_miri_img.xml')
    yamls = call_yaml_generator(xml)
    assert len(yamls.yaml_files) == 10

def test_niriss_wfss_nrc_img():
    xml = os.path.join(test_input_dir, 'niriss_wfss_nrc_img.xml')
    yamls = call_yaml_generator(xml)
    assert len(yamls.yaml_files) == 48
    nrc_yamls = np.array([True if 'nrc' in ele else False for ele in yamls.yaml_files])
    nis_yamls = np.array([True if 'nis' in ele else False for ele in yamls.yaml_files])
    assert np.sum(nrc_yamls) == 40
    assert np.sum(nis_yamls) == 8

def test_nrc_img_and_miri_img():
    xml = os.path.join(test_input_dir, 'nrc_img_and_miri_img.xml')
    yamls = call_yaml_generator(xml)
    assert len(yamls.yaml_files) == 60

def test_nrc_img_niriss_img():
    xml = os.path.join(test_input_dir, 'nrc_img_niriss_img.xml')
    yamls = call_yaml_generator(xml)
    assert len(yamls.yaml_files) == 36
    nrc_yamls = np.array([True if 'nrc' in ele else False for ele in yamls.yaml_files])
    nis_yamls = np.array([True if 'nis' in ele else False for ele in yamls.yaml_files])
    assert np.sum(nrc_yamls) == 30
    assert np.sum(nis_yamls) == 6

def test_nrs_mos_and_nrc_img():
    xml = os.path.join(test_input_dir, 'nrs_mos_and_nrc_img.xml')
    yamls = call_yaml_generator(xml)
    assert len(yamls.yaml_files) == 20

def test_skipped_obs():
    xml = os.path.join(test_input_dir, 'skipped_obs.xml')
    yamls = call_yaml_generator(xml)
    assert len(yamls.yaml_files) == 21
    nrc_yamls = np.array([True if 'nrc' in ele else False for ele in yamls.yaml_files])
    nis_yamls = np.array([True if 'nis' in ele else False for ele in yamls.yaml_files])
    assert np.sum(nrc_yamls) == 20
    assert np.sum(nis_yamls) == 1













