"""Various utility functions.

Authors
-------
    - Lauren Chambers

Use
---
    This module can be imported as such:

    >>> import utils
    siaf_files = get_siaf()

Dependencies
------------
    The user must have a configuration file named ``siaf_config.json``
    placed in mirage/config/ directory.
"""

import json
import os
import re

from astropy.io import ascii as asc


def calc_frame_time(instrument, aperture, xdim, ydim, amps):
    """Calculate the readout time for a single frame
    of a given size and number of amplifiers. Note that for
    NIRISS and FGS, the fast readout direction is opposite to 
    that in NIRCam, so we switch xdim and ydim so that we can
    keep a single equation.

    Parameters:
    -----------
    instrument : str
        Name of the instrument being simulated

    aperture : str
        Name of aperture being simulated (e.g "NRCA1_FULL")
        Currently this is only used to check for the FGS
        ACQ1 aperture, which uses a unique value of colpad
        below.

    xdim : int
        Number of columns in the frame

    ydim : int
        Number of rows in the frame
    
    amps : int
        Number of amplifiers used to read out the frame

    Returns:
    --------
    frametime : float
        Readout time in seconds for the frame
    """
    instrument = instrument.lower()
    if instrument == "nircam":
        xs = xdim
        ys = ydim
        colpad = 12

        # Fullframe
        if amps == 4:
            rowpad = 1
            fullpad = 1
        else:
            # All subarrays
            rowpad = 2
            fullpad = 0
            if ((xdim <= 8) & (ydim <= 8)):
                # The smallest subarray
                rowpad = 3

    elif instrument == "niriss":
        xs = ydim
        ys = xdim
        colpad = 12

        # Fullframe
        if amps == 4:
            rowpad = 1
            fullpad = 1
        else:
            rowpad = 2
            fullpad = 0

    elif instrument == 'fgs':
        xs = ydim
        ys = xdim
        colpad = 6
        if 'acq1' in aperture.lower():
            colpad = 12
        rowpad = 1
        if amps == 4:
            fullpad = 1
        else:
            fullpad = 0

    return ((1.0 * xs / amps + colpad) * (ys + rowpad) + fullpad) * 1.e-5
            
def get_siaf():
    '''Return a dictionary that holds the contents of the SIAF config
    file.

    Returns
    -------
    siaf_files : dict
        A dictionary that holds the contents of the config file.
    '''
    scripts_dir = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    package_dir = os.path.dirname(scripts_dir)
    config_file = os.path.join(package_dir, 'config', 'siaf_config.json')

    with open(config_file, 'r') as config_file:
        siaf_files = json.load(config_file)

    return siaf_files


def get_aperture_definition(aperture_name, instrument):
    '''Parses the SIAF to get the definition of a given aperture
    '''

    siaf_file = get_siaf()[instrument]
    siaf_table = asc.read(siaf_file, header_start=1)
    row_ind = list(siaf_table['AperName']).index(aperture_name)
    siaf_row = siaf_table[row_ind]

    detector = aperture_name[3:5]  # Won't work for ALL, A, or B case

    if re.search(r"_F\d\d\d[MWN]", aperture_name):
        filt = aperture_name.split('_')[-1]
    else:
        filt = 'ANY'

    refpix_x = siaf_row['XSciRef']
    refpix_y = siaf_row['YSciRef']
    refpix_v2 = siaf_row['V2Ref']
    refpix_v3 = siaf_row['V3Ref']

    x_size = siaf_row['XSciSize']
    y_size = siaf_row['YSciSize']

    xstart = refpix_x - x_size / 2
    xend = refpix_x + x_size / 2
    ystart = refpix_y - y_size / 2
    yend = refpix_y + y_size / 2

    # Need to add num_amps

    aperture_definition = [aperture_name, detector, filt, xstart, ystart, xend,
                           yend, refpix_x, refpix_y, refpix_v2, refpix_v3]
    return aperture_definition
