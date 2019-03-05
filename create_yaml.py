#! /usr/bin/env python

from mirage.yaml import generate_observationlist, yaml_generator
from mirage.apt import read_apt_xml, apt_inputs
from mirage.utils import siaf_interface


xml_name = '793.xml'
pointing_name = '793.pointing'
output_directory = './yaml_files/'
simdata_output_directory = './simulated_mirage_data/'
data_type = 'linear, raw'
#catalogues = {'niriss' : 'stars_field19_combined_allfilters.list' }
catalogues = {'niriss' : 'field19.list' }
observation_file = '793_observation_list.yaml'
generate_observationlist.get_observation_dict(xml_name,observation_file,catalogs=catalogues)
parameter_defaults = {'PAV3': 275.,'DATE':'2020-09-20'}
yam = yaml_generator.SimInput(input_xml=xml_name, pointing_file=pointing_name,
                              catalogs=catalogues, observation_list_file=observation_file,
                              verbose=True, output_dir=output_directory, simdata_output_dir=simdata_output_directory,
                              use_JWST_pipeline=True,offline=True, parameter_defaults=parameter_defaults)

yam.reffile_setup()
yam.create_inputs()

"""
from glob import glob
from scipy.stats import sigmaclip
import numpy as np
from astropy.io import fits
from astropy.visualization import simple_norm
import matplotlib.pyplot as plt


from mirage.seed_image import catalog_seed_image
from mirage.dark import dark_prep
from mirage.ramp_generator import obs_generator

from mirage.apt import apt_inputs
from mirage.yaml import yaml_generator
# Create a series of data simulator input yaml files
# from APT files
yam = yaml_generator.SimInput()
yam.input_xml = '793.xml'
yam.pointing_file = '793.pointing'
yam.output_dir = './yaml_files/'
yam.simdata_output_dir = './simulated_mirage_data/'
yam.observation_table = '793_observation_list.yaml'
yam.use_JWST_pipeline = True
yam.use_linearized_darks = False
yam.datatype = 'linear,raw'
yam.parameter_defaults = {'PAV3': 272.,'DATE':'2020-10-15'}
#yam.catalogs = {'niriss' : 'stars_field20_combined_allfilters.list' }
yam.reffile_setup()
yam.create_inputs()

"""
"""
if MAKE_YAML:
    # set me up with yaml files first...
    xml_name = '/Users/anand/Documents/NIRISS/COM/com1093smalltest.xml'
    pointing_name = '/Users/anand/Documents/NIRISS/COM/com1093smalltest.pointing'
    output_directory = '/Users/anand/Documents/NIRISS/COM/com1093smalltest_mirout/'
    data_type = 'linear, raw'
    # AB Dor is field20 from Kevin's excel file... /ifs/jwst/wit/niriss/kevin/cars_simulations/targets.xlsx
    catalogues = {'niriss' : 'stars_field20_combined_allfilters.list' }
    observation_file = 'com1093smalltest_observation_list.yaml'
    generate_observationlist.get_observation_dict(xml_name,observation_file,catalogs=catalogues)
    parameter_defaults = {'PAV3': 272.,'DATE':'2020-10-15'}
    yam = yaml_generator.SimInput(input_xml=xml_name, pointing_file=pointing_name,
                                  catalogs=catalogues, observation_list_file=observation_file,
                                  verbose=True, output_dir=output_directory, simdata_output_dir=output_directory,
                                  offline=True, parameter_defaults=parameter_defaults)
    yam.create_inputs()
"""
