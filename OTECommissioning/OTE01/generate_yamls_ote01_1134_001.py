import glob
import os

from mirage import image_simulator
from mirage.catalogs import get_catalog
from mirage.yaml import yaml_generator

# Import APT files
prop_id = '1134'
root = 'OTE01-1134-reduced_mosaic'
ote_dir = '/user/lchambers/OTECommSims/OTE01_reducedmosaic/'
pointing_file = os.path.join(ote_dir, root + '.pointing')
xml_file = os.path.join(ote_dir, root + '.xml')

# Get SW and LW catalogs
cats = get_catalog.get_all_catalogs(pointing_file, prop_id)
target_coords, catalog_filenames_sw, catalog_filenames_lw = cats

# Create a series of data simulator input yaml files from APT files
print('Beginning to generate YAML files')
yam = yaml_generator.SimInput()
yam.input_xml = xml_file
yam.pointing_file = pointing_file
siaf_file = os.path.expandvars('$MIRAGE_DATA/nircam/reference_files/SIAF/NIRCam_SIAF_2018-01-08.csv')
yam.siaf = siaf_file
yam.output_dir = ote_dir #os.path.join(os.getcwd(), ote_dir)
yam.simdata_output_dir = ote_dir #os.path.join(os.getcwd(), ote_dir)
yam.observation_table = os.path.join(ote_dir, 'OTE01-1134-reduced_mosaic_observationlist.yaml')

yam.use_JWST_pipeline = False # changed to False
yam.use_linearized_darks = True # changed to True
yam.datatype = 'linear'

yam.reffile_setup()
yam.create_inputs()

# Print out how many YAMLs were generated
yaml_list = glob.glob(ote_dir + 'jw01134001*.yaml')
print('Found {} yaml files to generate images.'.format(len(yaml_list)))
