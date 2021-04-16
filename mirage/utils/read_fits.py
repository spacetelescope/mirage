#! /usr/bin/env python

'''
Module for reading in a given fits file using either RampModel
or astropy, packaging information into a common format,
and returning

metadata needed within the simulator (these need to be in
a common format regardless of read-in method). Other meta
data can stay in whatever format produced when it is read
in, under the assmption that it will be written using the
same format

read pattern - exposure.readpatt READPATT
nints - exposure.nints NINTS
ngroups - exposure.ngroups NGROUPS
nframes - exposure.nframes NFRAMES
nskip - exposure.nskip NSKIP
groupgap - exposure.groupgap GROUPGAP
exptype = exposure.type EXP_TYPE
detector - instrument.detector DETECTOR
instrument - instrument.name INSTRUME
'''
import logging
import numpy as np
import os

from astropy.io import fits
from jwst.datamodels import RampModel

from mirage.logging import logging_functions
from mirage.utils.constants import LOG_CONFIG_FILENAME, STANDARD_LOGFILE_NAME


classdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../'))
log_config_file = os.path.join(classdir, 'logging', LOG_CONFIG_FILENAME)
logging_functions.create_logger(log_config_file, STANDARD_LOGFILE_NAME)


class Read_fits():
    def __init__(self):
        self.translate = {}
        self.translate['READPATT'] = 'exposure.readpatt'
        self.translate['NINTS'] = 'exposure.nints'
        self.translate['NGROUPS'] = 'expsoure.ngroups'
        self.translate['NFRAMES'] = 'exposure.nframes'
        self.translate['NSKIP'] = 'exposure.nskip'
        self.translate['GROUPGAP'] = 'exposure.groupgap'
        self.translate['TFRAME'] = 'exposure.frame_time'
        self.translate['TGROUP'] = 'exposure.group_time'
        self.translate['EXP_TYPE'] = 'exposure.type'
        self.translate['DETECTOR'] = 'instrument.detector'
        self.translate['INSTRUME'] = 'instrument.name'
        self.translate['FASTAXIS'] = 'subarray.fastaxis'
        self.translate['SLOWAXIS'] = 'subarray.slowaxis'

    def rampmodel_to_obj(self):
        # convert a RampModel instance to a read_fits object
        self.data = self.model.data
        self.zeroframe = self.model.zeroframe
        self.sbAndRefpix = None
        self.zero_sbAndRefpix = None

        self.header = {}
        for key in self.translate:
            try:
                self.header[key] = self.model.meta[self.translate[key]]
            except:
                self.header[key] = None

    def read_astropy(self):

        h = fits.open(self.file)

        self.data = None
        self.zeroframe = None
        self.sbAndRefpix = None
        self.zero_sbAndRefpix = None
        for i in range(len(h)):
            name = h[i].name
            if name == 'SCI':
                self.data = h[i].data
            if name == 'ZEROFRAME':
                self.zeroframe = h[i].data
            if name == 'SBANDREFPIX':
                self.sbAndRefpix = h[i].data
            if name == 'ZEROSBANDREFPIX':
                self.zero_sbAndRefpix = h[i].data

        #to match what happens with the RampModel version,
        #populate any of the remaining None extensions with
        #arrays of zeros
        #if self.zeroframe is None:
        #    self.zeroframe = np.zeros((ngroup,ny,nx),dtype=np.float)

        self.header = {}
        for key in self.translate:
            try:
                self.header[key] = h[0].header[key]
            except:
                self.header[key] = None

    def read_datamodel(self):
        logger = logging.getLogger('mirage.utils.read_fits.read_datamodel')

        h = RampModel(self.file)

        #remove any non-pipeline related keywords (e.g. CV3 temps/voltages)
        if 'extra_fits' in dir(h):
            h.__delattr__('extra_fits')

        self.data = h.data

        #Currently a bug in level1bmodel when zeroframe is
        #not present and a cube of zeros is returned
        #If the datamodel returns a default zeroframe of all
        #zeros, then set it to None here
        #self.zeroframe = h.zeroframe
        if np.all(h.zeroframe == 0):
            logger.info("Zeroframe in {} is all zeros. Returning None.".format(self.file))
            self.zeroframe = None
        else:
            self.zeroframe = h.zeroframe

        self.sbAndRefpix = None

        self.header = {}
        for key in self.translate:
            try:
                self.header[key] = h.meta[self.translate[key]]
            except:
                self.header[key] = None

    def insert_into_datamodel(self,subfile):
        #read in a dummy/substitute file as a datamodel,
        #and insert the data and self.header metadata
        #into it

        h = RampModel(subfile)
        h.data = self.data
        try:
            h.zeroframe = self.zeroframe
        except:
            pass

        h.err = np.zeros_like(self.data)
        h.groupdq = np.zeros_like(self.data)
        nint,ng,ny,nx = self.data.shape
        h.pixeldq = np.zeros((ny,nx))

        h.meta.exposure.readpatt = self.header['READPATT']
        h.meta.exposure.nints = self.header['NINTS']
        h.meta.exposure.ngroups = self.header['NGROUPS']
        h.meta.exposure.nframes = self.header['NFRAMES']
        h.meta.exposure.nskip = self.header['NSKIP']
        h.meta.exposure.groupgap = self.header['GROUPGAP']
        h.meta.exposure.type = self.header['EXP_TYPE']
        h.meta.instrument.detector = self.header['DETECTOR']
        h.meta.instrument.name = self.header['INSTRUME']
        h.meta.subarray.fastaxis = self.header['FASTAXIS']
        h.meta.subarray.slowaxis = self.header['SLOWAXIS']

        return h
