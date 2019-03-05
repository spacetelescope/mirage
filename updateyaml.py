#! /usr/bin/env python

import glob
import pathlib
from ruamel.yaml import YAML

files = glob.glob('yaml_files/jw*yaml')
def update(value1,value2,value3,value4):
    yaml = YAML()
    yaml.preserve_quotes = True
    for i in files:
        mf = pathlib.Path(i)
        doc = yaml.load(mf)
        doc['Reffiles']['astrometric'] = value1
        doc['Reffiles']['distortion_coeffs'] = value2
        doc['Output']['datatype'] = value3
        doc['Output']['save_intermediates'] = value4
        yaml.dump(doc, mf)

update('None',\
       '/ifs/jwst/wit/mirage_data/niriss/reference_files/SIAF/NIRISS_SIAF_09-28-2017.csv',\
       'linear,raw',\
       'true')






