#! /usr/bin/env python

'''
Given an imaging mode yaml file, convert relevent
parameters to those needed for WFSS simulations,
and output the updated yaml file.

This is meant primarily for WFSS observations, where
WFSS yaml files are needed for the initial seed
image generation, but then later, when running the
observation generator,

quantites that need to be updated:
mode - NO, CAN BE WFSS FOR ALL
filter
pupil?
output:file
grism_source_image
'''
import sys
import yaml


class YamlUpdate():
    def __init__(self):
        self.file = None
        self.filter = None
        self.pupil = None
        self.raw_outfile = None
        self.outname = None

    def read_yaml(self,file):
        try:
            with open(file,'r') as f:
                data = yaml.safe_load(f)
        except:
            print("WARNING: unable to open {}".format(self.paramfile))
            sys.exit()
        return data


    def wfss_seed_to_dispersed(self,indata):
        indata['Readout']['filter'] = self.filter
        indata['Readout']['pupil'] = self.pupil
        if self.raw_outfile is None:
            ofile = indata['Output']['file']
            suffix = ('_dispersed_{}_crossing_{}_uncal.fits'
                            .format(self.pupil,self.filter))
            self.raw_outfile = ofile.replace('.fits',suffix)
        indata['Output']['file'] = self.raw_outfile
        indata['Output']['grism_source_image'] = False
        return indata


    def run(self):
        # Read in yaml contents
        y = self.read_yaml(self.file)

        # Adjust the fields necessary to create a yaml
        # file that can be used in the obs_generator
        # for a WFSS observation
        y2 = self.wfss_seed_to_dispersed(y)

        # Save adjusted data into a new yaml file
        if self.outname is None:
            suffix = ('_dispersed_{}_crossing_{}.yaml'
                      .format(self.pupil,self.filter))
            self.outname = self.file.replace('.yaml',suffix)
        with open(self.outname,'w') as output:
            yaml.dump(y2,output,default_flow_style=False)
