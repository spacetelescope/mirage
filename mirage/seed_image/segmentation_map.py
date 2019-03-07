#! /usr/bin/env python

'''
Segmentation map creation. Developed in conjunction with
the seed image generator code for catalogs
catalog_seed_image.py
'''

import numpy as np


class SegMap():
    def __init__(self):
        self.xdim = 2048
        self.ydim = 2048
        self.zdim = None

    def initialize_map(self):
        if self.zdim is None:
            self.segmap = np.zeros((self.ydim, self.xdim), dtype=np.int64)
        else:
            self.segmap = np.zeros((self.zdim, self.ydim, self.xdim), dtype=np.int64)

    def add_object_basic(self, ystart, yend, xstart, xend, number):
        # Add an object to the segmentation map
        # in simplest way possible. All pixels
        # in the box are set to the index number
        # regardless of signal
        self.segmap[ystart:yend, xstart:xend] = number

    def add_object_perccut(self, image, ystart, xstart, number, perc):
        # Add an object to the segmentation map
        # In this case, only flag pixels whose
        # signal is higher than the given cutoff
        yd, xd = image.shape
        stamp = self.segmap[ystart:ystart+yd, xstart:xstart+xd]
        maxsig = np.max(image)
        cutoff = maxsig * perc
        flag = image >= cutoff
        stamp[flag] = number

    def add_object_noise(self, image, ystart, xstart, number, noise):
        # Add an object to the segmentation map
        # In this case, only flag pixels whose
        # signal is higher than a calculated cutoff
        # that is based on average noise/background
        # values for NIRCam
        yd, xd = image.shape
        stamp = self.segmap[ystart:ystart+yd, xstart:xstart+xd]
        flag = image >= noise
        stamp[flag] = number
