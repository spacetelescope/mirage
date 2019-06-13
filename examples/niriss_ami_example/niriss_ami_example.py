#! /usr/bin/env python
import glob
import io
import os
import yaml

from mirage.yaml import generate_observationlist, yaml_generator
from mirage import imaging_simulator

home_dir = os.environ['HOME']

ami_example_dir = os.path.dirname(__file__)
print(ami_example_dir,"===================")

xml_name = os.path.join(ami_example_dir, 'input_files/793_mirage_example.xml')
pointing_name = os.path.join(ami_example_dir, 'input_files/793_mirage_example.pointing')

output_directory = os.path.join(home_dir, 'ami_mirage_simulation_example')
simdata_output_directory = output_directory

# AB Dor is in field 19 and HD37093 is in field 20 in Kevin's targets.xlsx file
catalogues = {'niriss': os.path.join(ami_example_dir, 'stars_field19_20_combined_allfilters.list')}
parameter_defaults = {'PAV3': 275., 'DATE': '2020-09-20'}

observation_file = os.path.join(output_directory, '793_observation_list.yaml')

generate_observationlist.get_observation_dict(xml_name, observation_file, catalogs=catalogues)
yam = yaml_generator.SimInput(input_xml=xml_name, pointing_file=pointing_name,
                              catalogs=catalogues, observation_list_file=observation_file,
                              verbose=True, output_dir=output_directory, simdata_output_dir=simdata_output_directory,
                              use_JWST_pipeline=True, offline=False, parameter_defaults=parameter_defaults)

datatype = 'linear, raw'
yam.datatype = datatype
yam.create_inputs()


yaml_files = glob.glob(os.path.join(output_directory, 'jw*.yaml'))

for file in yaml_files:

    # set astrometric reference file to None to use pysiaf
    with open(file, 'r') as infile:
        yaml_content = yaml.load(infile)
    yaml_content['Reffiles']['astrometric'] = 'None'
    modified_file = file.replace('.yaml', '_mod.yaml')
    with io.open(modified_file, 'w') as outfile:
        yaml.dump(yaml_content, outfile, default_flow_style=False)

    t1 = imaging_simulator.ImgSim()
    t1.paramfile = str(modified_file)
    t1.create()

