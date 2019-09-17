#! /usr/bin/env python

"""Create test RA, Dec data to be used in test_ra_dec_to_x_y.py. Travis does
not have access to the reference files, so this must be done manually, and
the results placed into test_ra_dec_to_x_y.py.
"""
import asdf
import glob
import os
import numpy as np
import pysiaf

from mirage.seed_image import catalog_seed_image
from mirage.utils import siaf_interface
from mirage.utils.utils import expand_environment_variable


data_dir = expand_environment_variable('MIRAGE_DATA', offline=False)
instruments = ['nircam', 'niriss', 'fgs']
xyvals = {}

# Create list of points at which to calculate x, y
pointing_ra = 12.0
pointing_dec = 12.0
rotation = 0.0
ra_list = []
dec_list = []
delta_pix = np.array([0, 10, 20, 40, 75, 180, 500, 1000])
delta_deg = delta_pix * 0.031 / 3600.
for delta in delta_deg:
    ra_entry = pointing_ra + delta
    dec_entry = pointing_dec + delta
    ra_list.append(ra_entry)
    dec_list.append(dec_entry)

c = catalog_seed_image.Catalog_seed()
for instrument in instruments:
    siaf = siaf_interface.get_instance(instrument)

    if instrument.lower() == 'nircam':
        a_detectors = ['NRCA{}'.format(i+1) for i in range(5)]
        b_detectors = ['NRCB{}'.format(i+1) for i in range(5)]
        detectors = a_detectors + b_detectors
    elif instrument.lower() == 'niriss':
        detectors = ['NIS']
    elif instrument.lower() == 'fgs':
        detectors = ['FGS1', 'FGS2']

    print('\n\nDetectors: {}\n\n'.format(detectors))

    # Find reference files
    reffile_dir = os.path.join(data_dir, '{}/reference_files/distortion'.format(instrument))
    reffiles = glob.glob(os.path.join(reffile_dir, '*.asdf'))

    for detector in detectors:

        # Look through reference files, find file with matching detector
        for reffile in reffiles:
            with asdf.open(reffile) as f:
                file_det = f.tree['meta']['instrument']['detector'].upper()
                c.coord_transform = f.tree['model']

            if '5' in detector:
                match_detector = detector.replace('5', 'LONG')
            else:
                match_detector = detector

            if 'FGS' in detector:
                match_detector = detector.replace('FGS', 'GUIDER')

            if file_det == match_detector:
                apertures = siaf.apernames
                for aperture in apertures:

                    # Skip the raw orientation apertures
                    #if '_SUB' in aperture:
                    if 'OSS' not in aperture and detector in aperture:
                        c.local_roll, c.attitude_matrix, c.ffsize, \
                         c.subarray_bounds = siaf_interface.get_siaf_information(siaf,
                                                                                 aperture,
                                                                                 pointing_ra, pointing_dec,
                                                                                 rotation)
                        xvals_from_reffile = []
                        yvals_from_reffile = []
                        for ra, dec in zip(ra_list, dec_list):
                            x, y = c.RADecToXY_astrometric(ra, dec)
                            xvals_from_reffile.append(x)
                            yvals_from_reffile.append(y)
                        xyvals['{}_{}'.format(instrument.lower(), aperture)] = (xvals_from_reffile, yvals_from_reffile)
for key in xyvals:
    print("'{}': {},".format(key, xyvals[key]))
