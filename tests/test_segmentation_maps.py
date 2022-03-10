#! /usr/bin/env python

"""Test aspects of the segmentation map creation
"""
import numpy as np
from mirage.seed_image.segmentation_map import SegMap


def test_add_object_basic():
    """Set given pixels to contain an object
    """
    object_num = 8
    map = SegMap()
    map.initialize_map()
    map.add_object_basic(400, 410, 200, 210, object_num)
    assert np.unique(map.segmap[400:410, 200:210]) == object_num


def test_add_object_percut():
    """Add pixels above a given percentage cutoff
    """
    object_num = 8
    image = np.zeros((10, 10))
    image[:, 4] = 10.
    image[:, 3] = 6.

    map = SegMap()
    map.initialize_map()
    map.add_object_perccut(image, 300, 100, object_num, 0.65)

    # Only the pixels in column 4 should be added to the map
    assert np.unique(map.segmap[300:310, 104]) == object_num
    assert np.unique(map.segmap[300:310, 103]) == 0


def test_add_object_threshold():
    """Add pixels above a signal threshold to the segmap
    """
    object_num = 8
    image = np.zeros((10, 10))
    image[:, 4] = 10.
    image[:, 3] = 6.

    map = SegMap()
    map.initialize_map()
    map.add_object_threshold(image, 300, 100, object_num, 7.5)

    # Only the pixels in column 4 should be added to the map
    assert np.unique(map.segmap[300:310, 104]) == object_num
    assert np.unique(map.segmap[300:310, 103]) == 0


def test_dimensionality():
    """Create segmaps with various numbers of dimensions
    """
    map = SegMap()
    map.initialize_map()
    assert map.segmap.shape == ((2048, 2048))

    map = SegMap()
    map.zdim = 3
    map.initialize_map()
    assert map.segmap.shape == ((3, 2048, 2048))

    map = SegMap()
    map.zdim = 3
    map.intdim = 4
    map.initialize_map()
    assert map.segmap.shape == ((4, 3, 2048, 2048))
