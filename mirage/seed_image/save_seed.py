#! /usr/bin/env python

"""
Utility for saving seed images
"""
import logging
import os

from astropy.io import fits
import numpy as np

import mirage
from mirage.logging import logging_functions
from mirage.utils.constants import LOG_CONFIG_FILENAME, STANDARD_LOGFILE_NAME


classdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../'))
log_config_file = os.path.join(classdir, 'logging', LOG_CONFIG_FILENAME)
logging_functions.create_logger(log_config_file, STANDARD_LOGFILE_NAME)



def save(seed_image, param_file, parameters, photflam, photfnu, pivot_wavelength,
         fullframe_size, nominal_dimensions, coord_adjust, grism_direct_factor,
         filename=None, segmentation_map=None, frametime=None, base_unit='ADU'):
    """Save a seed image
    """
    logger = logging.getLogger('mirage.seed_image.save_seed.save')

    arrayshape = seed_image.shape
    if len(arrayshape) == 2:
        units = '{}/sec'.format(base_unit)
        yd, xd = arrayshape
        tgroup = 0.
        logger.info('Seed image is 2D.')
    elif len(arrayshape) == 3:
        units = base_unit
        g, yd, xd = arrayshape
        tgroup = frametime * (parameters['Readout']['nframe'] + parameters['Readout']['nskip'])
        logger.info('Seed image is 3D.')
    elif len(arrayshape) == 4:
        units = base_unit
        integ, g, yd, xd = arrayshape
        tgroup = frametime * (parameters['Readout']['nframe'] + parameters['Readout']['nskip'])
        logger.info('Seed image is 4D.')

    xcent_fov = xd / 2
    ycent_fov = yd / 2

    kw = {}
    kw['xcenter'] = xcent_fov
    kw['ycenter'] = ycent_fov
    kw['units'] = kw['UNITS'] = units
    kw['TGROUP'] = tgroup
    if parameters['Readout']['pupil'][0].upper() == 'F':
        usefilt = 'pupil'
    else:
        usefilt = 'filter'

    if filename is None:
        basename = os.path.join(parameters['Output']['directory'],
                                parameters['Output']['file'][0:-5].split('/')[-1])
        filename = '{}_{}_seed_image.fits'.format(basename, parameters['Readout'][usefilt])

    # Set FGS filter to "N/A" in the output file
    # as this is the value DMS looks for.
    if parameters['Readout'][usefilt] == "NA":
        parameters['Readout'][usefilt] = "N/A"
    kw['filter'] = parameters['Readout'][usefilt]
    kw['PHOTFLAM'] = photflam
    kw['PHOTFNU'] = photfnu
    kw['PHOTPLAM'] = pivot_wavelength * 1.e4  # put into angstroms
    kw['NOMXDIM'] = nominal_dimensions[1]
    kw['NOMYDIM'] = nominal_dimensions[0]
    kw['NOMXSTRT'] = int(coord_adjust['xoffset'] + 1)
    kw['NOMXEND'] = int(nominal_dimensions[1] + coord_adjust['xoffset'])
    kw['NOMYSTRT'] = int(coord_adjust['yoffset'] + 1)
    kw['NOMYEND'] = int(nominal_dimensions[0] + coord_adjust['yoffset'])

    # Files/inputs used during seed image production
    kw['YAMLFILE'] = param_file
    kw['GAINFILE'] = parameters['Reffiles']['gain']
    kw['DISTORTN'] = parameters['Reffiles']['astrometric']
    kw['IPC'] = parameters['Reffiles']['ipc']
    kw['PIXARMAP'] = parameters['Reffiles']['pixelAreaMap']
    kw['CROSSTLK'] = parameters['Reffiles']['crosstalk']
    kw['FLUX_CAL'] = parameters['Reffiles']['flux_cal']
    kw['FTHRUPUT'] = parameters['Reffiles']['filter_throughput']
    kw['PTSRCCAT'] = parameters['simSignals']['pointsource']
    kw['GALAXCAT'] = parameters['simSignals']['galaxyListFile']
    kw['EXTNDCAT'] = parameters['simSignals']['extended']
    kw['MTPTSCAT'] = parameters['simSignals']['movingTargetList']
    kw['MTSERSIC'] = parameters['simSignals']['movingTargetSersic']
    kw['MTEXTEND'] = parameters['simSignals']['movingTargetExtended']
    kw['NONSDRAL'] = parameters['simSignals']['movingTargetToTrack']
    kw['BKGDRATE'] = parameters['simSignals']['bkgdrate']
    kw['TRACKING'] = parameters['Telescope']['tracking']
    kw['POISSON'] = parameters['simSignals']['poissonseed']
    kw['PSFWFE'] = parameters['simSignals']['psfwfe']
    kw['PSFWFGRP'] = parameters['simSignals']['psfwfegroup']
    kw['MRGEVRSN'] = mirage.__version__

    # Seed images provided to disperser are always embedded in an array
    # with dimensions equal to full frame * self.grism_direct_factor
    if parameters['Inst']['mode'] in ['wfss', 'ts_wfss']:
        kw['NOMXDIM'] = fullframe_size
        kw['NOMYDIM'] = fullframe_size
        kw['NOMXSTRT'] = int(fullframe_size * (grism_direct_factor - 1) / 2.)
        kw['NOMXEND'] = kw['NOMXSTRT'] + fullframe_size - 1
        kw['NOMYSTRT'] = int(fullframe_size * (grism_direct_factor - 1) / 2.)
        kw['NOMYEND'] = kw['NOMYSTRT'] + fullframe_size - 1

    kw['GRISMPAD'] = grism_direct_factor
    seedinfo = kw
    save_single_fits(seed_image, filename, key_dict=kw, image2=segmentation_map, image2type='SEGMAP')

    # Keep this print statement in the code that calls this function
    #print("Seed image and segmentation map saved as {}".format(self.seed_file))
    #print("Seed image, segmentation map, and metadata available as:")
    #print("self.seedimage, self.seed_segmap, self.seedinfo.")
    return filename, seedinfo


def save_single_fits(image, name, key_dict=None, image2=None, image2type=None):
        # Save an array into the first extension of a fits file
        h0 = fits.PrimaryHDU()
        h1 = fits.ImageHDU(image, name='DATA')
        if image2 is not None:
            h2 = fits.ImageHDU(image2)
            if image2type is not None:
                h2.header['EXTNAME'] = image2type

        # if a keyword dictionary is provided, put the
        # keywords into the 0th and 1st extension headers
        if key_dict is not None:
            for key in key_dict:
                h0.header[key] = key_dict[key]
                h1.header[key] = key_dict[key]

        if image2 is None:
            hdulist = fits.HDUList([h0, h1])
        else:
            hdulist = fits.HDUList([h0, h1, h2])
        hdulist.writeto(name, overwrite=True)

