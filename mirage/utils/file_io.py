#! /usr/bin/env python

"""
This module contains utilities used for reading and writing various files
used or created by Mirage.

Author
------

    - Bryan Hilbert

Use
---

    This module can be imported and used as such:

    ::

        from mirage.utils import file_io
        gain, gain_head = file_io.read_gain_file(filename)
"""
import os

from astropy.io import ascii, fits
import numpy as np


def read_filter_throughput(filename):
    """Read in the ascii file containing a filter throughput curve

    Parameters
    ----------
    filename : str
        Name of ascii file containing throughput info

    Returns
    -------
    (wavelengths, throughput) : tup
        Tuple of 1D numpy arrays containing the wavelength and throughput
        values from the file
    """
    tab = ascii.read(filename)
    return tab['Wavelength_microns'].data, tab['Throughput'].data


def read_gain_file(filename):
    """
    Read in CRDS-formatted gain reference file

    Paramters
    ---------
    filename : str
        Name of (fits) gain reference file

    Returns
    -------
    image : numpy.ndarray
        2D array gain map

    header : dict
        Information contained in the header of the file
    """
    try:
        with fits.open(filename) as h:
            image = h[1].data
            header = h[0].header
    except (FileNotFoundError, IOError) as e:
        print(e)

    mngain = np.nanmedian(image)

    # Set pixels with a gain value of 0 equal to mean
    image[image == 0] = mngain

    # Set any pixels with non-finite values equal to mean
    image[~np.isfinite(image)] = mngain
    return image, header
