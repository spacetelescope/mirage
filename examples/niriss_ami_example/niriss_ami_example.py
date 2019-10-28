#! /usr/bin/env python
import glob
import io
import os
import yaml
import sys

from mirage import imaging_simulator
from mirage.yaml import yaml_generator


home_dir = os.environ['HOME']

ami_example_dir = os.path.dirname(__file__)
print(ami_example_dir)

xml_name = os.path.join(ami_example_dir, 'input_files/793_mirage_example.xml')
pointing_name = os.path.join(ami_example_dir, 'input_files/793_mirage_example.pointing')

output_directory = os.path.join(home_dir, 'ami_mirage_simulation_example')
simdata_output_directory = output_directory

# AB Dor is in field 19 and HD37093 is in field 20 in Kevin's targets.xlsx file
catalogues = {'AB-DOR': {'point_source': os.path.join(ami_example_dir,'stars_field19_20_combined_allfilters.list')
                        },
              'HD-37093': {'point_source': os.path.join(ami_example_dir,'stars_field19_20_combined_allfilters.list')
                        }
             }
pav3 = 275
dates = '2020-09-20'
reffile_defaults = 'crds'
datatype = 'linear, raw'


yam = yaml_generator.SimInput(input_xml=xml_name, pointing_file=pointing_name,
                              catalogs=catalogues, roll_angle=pav3,
                              dates=dates, reffile_defaults=reffile_defaults,
                              verbose=True, output_dir=output_directory,
                              simdata_output_dir=simdata_output_directory,
                              datatype=datatype)

yam.create_inputs()


#Create all files
yaml_files = glob.glob(os.path.join(output_directory, 'jw*.yaml'))
#create one file for testing
#yaml_files = glob.glob(os.path.join(output_directory, 'jw00793001001_01101_00001_nis.yaml'))
print(yaml_files)

for file in yaml_files:

    # set astrometric reference file to None to use pysiaf
    with open(file, 'r') as infile:
        yaml_content = yaml.load(infile)
    yaml_content['Reffiles']['astrometric'] = 'None'
    yaml_content['psf_wing_threshold_file'] = 'config'
    modified_file = file.replace('.yaml', '_mod.yaml')
    with io.open(modified_file, 'w') as outfile:
        yaml.dump(yaml_content, outfile, default_flow_style=False)

    t1 = imaging_simulator.ImgSim()
    t1.paramfile = str(modified_file)
    t1.create()

