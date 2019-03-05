#! /usr/bin/env python

from mirage.yaml import generate_observationlist, yaml_generator
from mirage.apt import read_apt_xml, apt_inputs
from mirage.utils import siaf_interface


xml_name = 'input_files/793.xml'
pointing_name = 'input_files/793.pointing'
output_directory = './yaml_files/'
simdata_output_directory = './simulated_mirage_data/'
data_type = 'linear, raw'
#AB Dor is in field 19 in Kevin's targets.xlsx file
catalogues = {'niriss' : 'stars_field19_combined_allfilters.list'}
observation_file = '793_observation_list.yaml'
generate_observationlist.get_observation_dict(xml_name,observation_file,catalogs=catalogues)
parameter_defaults = {'PAV3': 275.,'DATE':'2020-09-20'}
yam = yaml_generator.SimInput(input_xml=xml_name, pointing_file=pointing_name,
                              catalogs=catalogues, observation_list_file=observation_file,
                              verbose=True, output_dir=output_directory, simdata_output_dir=simdata_output_directory,
                              use_JWST_pipeline=True,offline=True, parameter_defaults=parameter_defaults)

yam.reffile_setup()
yam.create_inputs()
