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
        self.intdim = None

    def initialize_map(self):
        if self.zdim is None and self.intdim is None:
            self.segmap = np.zeros((self.ydim, self.xdim), dtype=np.int64)
        elif self.zdim is not None and self.intdim is None:
            self.segmap = np.zeros((self.zdim, self.ydim, self.xdim), dtype=np.int64)
        elif self.zdim is not None and self.intdim is not None:
            self.segmap = np.zeros((self.intdim, self.zdim, self.ydim, self.xdim), dtype=np.int64)

    def add_object_basic(self, ystart, yend, xstart, xend, number):
        # Add an object to the segmentation map
        # in simplest way possible. All pixels
        # in the box are set to the index number
        # regardless of signal
        ndim = len(self.segmap.shape)
        if ndim == 2:
            self.segmap[ystart:yend, xstart:xend] = number
        elif ndim == 3:
            self.segmap[:, ystart:yend, xstart:xend] = number
        elif ndim == 4:
            self.segmap[:, :, ystart:yend, xstart:xend] = number

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

    def add_object_threshold(self, image, ystart, xstart, number, threshold):
        """Add an object to the segmentation map
        In this case, only flag pixels whose
        signal is higher than the given threshold value

        Paramters
        ---------
        image : numpy.ndarray
            Image containing the source to be added to the segmentation map

        ystart : int
            Y coordinate, in the coordinate system of the full seed image
            of the lower left corner of '''image'''

        xstart : imt
            X coordinate, in the coordinate system of the full seed image
            of the lower left corner of '''image'''

        number : int
            Value to be placed into the segmentation map for this object

        threshold : float
            Pixels with signal values higher than this will be added to
            the segmentation map as part of this object
        """
        yd, xd = image.shape
        stamp = self.segmap[ystart:ystart+yd, xstart:xstart+xd]
        flag = image >= threshold
        stamp[flag] = number
