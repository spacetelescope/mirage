#! /usr/bin/env python

'''
Create a mirage-format seed image from
an input fits file. This fits file must contain a
valid WCS, and be distortion-free.

This code will extract a sub-image from the input
fits file. The subimage is centered at the input
RA and Dec (crop_center_ra, crop_center_dec) and
its size is determined by the aperture input.

The extracted sub-image is then blotted onto the
requested NIRCam detector's FOV, including the
distortion model.

The resulting image is then saved in such a format
that it can be used as input in the obs_generator.py,
image_simulator.py, or disperser code in the
mirage package. Essentially this script
is designed to take the place of catalog_seed_image.py
in the case where the user has a distortion-free
mosaic image that they wish to break up into
NIRCam images.

Another way to use this code would be for an input
fits image of some extended object that the user
wishes to distort according to the NIRCam distortion
model. This distorted file can then be placed in
an extended object catalog which will be used by
the catalog_seed_image.py step. In this method,
multiple extended objets could be distorted and
added to the same frame.

Inputs:
mosaicfile: The name of the fits file containing the
            distortion-free image to be blotted. The code
            currently assumes that the mosaicfile is
            oriented north-up.

From here, there are two ways to call fits_seed_image.py

1. parameter file: a yaml input file matching the format
   (optional)      required by other steps of the nircam
                   data simulator

2. Manual inputs:
   aperture -
              aperture name matching an entry in the
              subarray definition file (e.g. NRCB5_FULL)

   crop_center_ra, crop_center_dec -
              The RA and Dec at the center of the sub-image
              to crop from mosaicfile. This, in combination
              with the array size (determined from aperture)
              define the sub-image to use.

   blot_center_ra, blot_center_dec -
              The RA and Dec at the center of the blotted
              sub-image. If these are equal to crop_center_ra
              and dec, then the center of the cropped image
              will remain at the center of the blotted image.

   blot_pav3 -
              Position angle of the blotted image.

   flux_cal_file -
              Ascii file listing the zeropoints for all
              NIRCam filters. By setting this value equal
              to 'config', the code will use the flux
              calibration file that is packaged with the
              NIRCam data simulator. This information is
              used only to populate the PHOTFLAM, PHOTFNU,
              and PHOTPLAM header keywords.

   filter -
              The NIRCam filter to use. This is used to
              select the appropriate row from the flux
              calibration file.

   pupil -
              The NIRCam pupil value to use. This is used
              to select the appropriate row from the flux
              calibration file in the case where a filter
              that is located in the pupil wheel is used.
              If the pupil value does not start with 'F'
              (e.g. F405N), then it is ignored.

   grism_source_image -
              True/False. If you intend to send the seed
              image output to the disperser software, then
              set this value to True. The output image will
              then be made larger than the nominal aperture
              size by a factor of sqrt(2), in order for the
              disperser to keep track of sources just outside
              the nominal FOV, whose dispersed spectra may
              fall onto the detector.

   outfile -
              Name of the file to contain the output seed
              image.

   outdir -
              Directory in which to place outfile

Outputs:
self.seedimage, self.seed_segmap, self.seedinfo
contain the seed (countrate) image, the associated
segmentation map, and header information required
by subsequent steps of the nircam data simulator.
The segmentation map is only used in the case where
the seed image is dispersed by a call to the
disperser software.

Example calls:

1. Using parameter file
s = fits_seed_image.ImgSeed('fits_seed_image_test.yaml')
s.mosaicfile = 'large_distortionFree_mosaic_drz.fits'
s.crop_and_blot()

2. Manual inputs
s = fits_seed_image.ImgSeed()
s.mosaicfile = 'large_distortionFree_mosaic_drz.fits'
s.aperture = 'NRCB5_FULL'
s.crop_center_ra = 53.1
s.crop_center_dec = -27.8
s.blot_center_ra = 53.1
s.blot_center_de = -27.8
s.blot_pav3 = 0.
s.subarray_defs = 'config'
s.flux_cal_file = 'config'
s.filter = 'F250M'
s.pupil = 'CLEAR'
s.grism_source_image = False
s.outfile = 'test_mosaic_seed.fits'
s.outdir = '/my/simulation/outputs/'
s.crop_and_blot()
'''
import argparse
import os
import pkg_resources
import sys

from astropy.io import fits, ascii
from math import radians
import numpy as np
from photutils import detect_sources
import pysiaf
import yaml

from . import crop_mosaic, blot_image
from mirage.reference_files import crds_tools
from mirage.utils.siaf_interface import get_siaf_information


config_files = {'nircam': {'flux_cal':'NIRCam_zeropoints.list'
                           },
                'niriss': {'flux_cal': 'niriss_zeropoints.list'
                           }
                }

class ImgSeed:
    def __init__(self, paramfile=None, mosaicfile=None):
        self.mosaicfile = mosaicfile
        self.data_extension_number = 0.
        self.wcs_extension_number = 0.
        self.crop_center_ra = 0.
        self.crop_center_dec = 0.
        self.channel = None
        self.detector = None
        self.blot_center_ra = 0.
        self.blot_center_dec = 0.
        self.blot_pav3 = 0.
        self.aperture = ''
        self.flux_cal_file = ''
        self.filter = ''
        self.pupil = ''
        self.grism_source_image = False
        self.cropped_file = 'cropped_image.fits'
        self.outfile = 'seed_image_from_mosaic.fits'
        self.outdir = './'
        self.grism_direct_factor = np.sqrt(2.)

        # Locate the module files, so that we know where to look
        # for config subdirectory
        self.modpath = pkg_resources.resource_filename('mirage','')

        # self.coords contains the factor by which the
        # nominal output array size needs to be increased
        # (used for WFSS mode), as well as the coordinate
        # offset between the nominal output array coordinates,
        # and those of the expanded array. These are needed
        # mostly for WFSS observations, where the nominal output
        # array will not sit centered in the expanded output image.
        self.coords = {'x':1., 'xoffset':0., 'y':1., 'yoffset':0.}

        # If a paramfile is provided, read in and
        # set params
        if paramfile is not None:
            self.read_param_file(paramfile)

    def read_param_file(self, file):
        """Read the input yaml file into a nested dictionary

        Parameters
        ----------
        file : str
            Name of the input yaml file
        """
        if os.path.isfile(file):
            with open(file, 'r') as f:
                params = yaml.safe_load(f)
        else:
            raise FileNotFoundError(("ERROR: {} does not exist."
                   .format(file)))

        self.instrument = params['Inst']['instrument']
        self.aperture = params['Readout']['array_name']
        self.detector = self.aperture[3:5]
        if '5' in self.detector:
            self.channel = 'long'
        else:
            self.channel = 'short'

        self.crop_center_ra = params['Telescope']['ra']
        self.crop_center_dec = params['Telescope']['dec']
        self.blot_center_ra = params['Telescope']['ra']
        self.blot_center_dec = params['Telescope']['dec']
        self.blot_pav3 = params['Telescope']['rotation']

        self.flux_cal_file = params['Reffiles']['flux_cal']
        self.distortion_file = params['Reffiles']['astrometric']
        self.filter = params['Readout']['filter']
        self.pupil = params['Readout']['pupil']

        # If input is to be dispersed later, we need to
        # expand the region to be blotted
        if params['Output']['grism_source_image']:
            self.grism_coord_adjust()
        self.outfile = params['Output']['file']
        self.outdir = params['Output']['directory']
        self.siaf = pysiaf.Siaf(self.instrument)

        # If we have a yaml file where the distortion reference
        # file is specified only as 'crds', then we need to
        # retrieve this file from CRDS
        if self.distortion_file == 'crds':
            reference_dictionary = crds_tools.dict_from_yaml(params)
            mapping = crds_tools.get_reffiles(reference_dictionary, ['distortion'], download=True)
            self.distortion_file = mapping['distortion']

    def expand_config_paths(self, file, filetype):
        if file.lower() == 'config':
            sfile = config_files[self.instrument.lower()][filetype]
            return os.path.join(self.modpath, 'config', sfile)

    def fluxcal_info(self, fluxfile, filter, module):
        # Read in list of zeropoints/photflam/photfnu
        self.zps = ascii.read(fluxfile)
        mtch = ((self.zps['Filter'] == filter) \
                & (self.zps['Module'] == module))
        self.photflam = self.zps['PHOTFLAM'][mtch][0]
        self.photfnu = self.zps['PHOTFNU'][mtch][0]
        self.pivot = self.zps['Pivot_wave'][mtch][0]

    def aperture_info(self):
        """
        Return basic information on the reuqested aperture from SIAF
        """
        local_roll, attitude_matrix, framesize, subarray_bounds = get_siaf_information(self.siaf,
                                                                                       self.aperture,
                                                                                       self.crop_center_ra,
                                                                                       self.crop_center_dec,
                                                                                       self.blot_pav3)
        xstart, ystart, xend, yend = subarray_bounds

        self.subarr_bounds = {'xstart':xstart, 'ystart':ystart, 'xend':xend, 'yend':yend}

    def crop_and_blot(self):
        """Extract an area of the input fits file that is slightly larger
        than the FOV, and use the JWST pipeline version of blot (called
        resample) to resample to the correct pixel scale and introduce the
        proper distortion
        """

        # Locations of conifg files
        #self.flux_cal_file = self.expand_config_paths(self.flux_cal_file, 'flux_cal')

        if self.detector is None:
            self.detector = self.aperture[3:5]
        if self.channel is None:
            if '5' in self.detector:
                self.channel = 'long'
            else:
                self.channel = 'short'

        # Get information on the aperture
        self.aperture_info()
        xstart = self.subarr_bounds['xstart']
        ystart = self.subarr_bounds['ystart']
        xdim = self.subarr_bounds['xend'] - xstart + 1
        ydim = self.subarr_bounds['yend'] - ystart + 1
        module = self.detector[0].upper()

        # Flux calibration information
        if self.pupil[0].upper() == 'F':
            usefilt = self.pupil.upper()
        else:
            usefilt = self.filter.upper()
        self.fluxcal_info(self.flux_cal_file, usefilt, module)

        # Extract
        pixel_scale = self.siaf[self.aperture].XSciScale
        crop = crop_mosaic.Extraction(mosaicfile=self.mosaicfile,
                                      data_extension_number=self.data_extension_number,
                                      wcs_extension_number=self.wcs_extension_number,
                                      center_ra=self.crop_center_ra,
                                      center_dec=self.crop_center_dec,
                                      dimensions=(xdim, ydim),
                                      outfile=self.cropped_file,
                                      jwst_pixel_scale=pixel_scale)
        crop.extract()

        # Blot
        blot = blot_image.Blot(instrument=self.instrument, aperture=self.aperture,
                               ra=[self.blot_center_ra], dec=[self.blot_center_dec],
                               pav3=[self.blot_pav3], blotfile=crop.cropped,
                               distortion_file=self.distortion_file)
        blot.blot()

        for dm in blot.blotted_datamodels:
            # Create segmentation map
            # (used only when dispsersing image for WFSS)
            self.seed_segmap = self.make_segmap(dm)

            # Save
            self.save_seed_image(dm)

    def make_segmap(self, model):
        """Create a segmentation map of the input image

        Parameters
        ----------
        model : jwst.datamodels.ImageModel

        Returns
        -------
        map : numpy.ndarray
            2D segmentation map
        """
        noise = np.median(model.data)
        seg_map = detect_sources(model.data, noise*5., 8)
        if seg_map is None:
            print('No segmentation map created. Returning empty segmap')
            map_image = np.zeros_like(model.data)
        else:
            map_image = seg_map.data
        return map_image

    def save_seed_image(self, model):
        """Save the seed image in the proper format that it can be used
        in subsequent steps of mirage

        Parameters
        ----------
        model : jwst.datamodels.ImageModel
        """
        self.seed_image = model.data
        array_shape = self.seed_image.shape
        dim = len(array_shape)
        if dim == 2:
            units = 'e-/sec'
            yd, xd = array_shape
            tgroup = 0.
        else:
            raise ValueError(("WARNING: {}D seed images from "
                   "input fits file not yet supported."
                   .format(dim)))

        dot = self.outfile.rfind('.fits')
        if dot == -1:
            dot = len(self.outfile)

        outfile = self.outfile.replace('.fits', '_{}_seed_image.fits'.format(self.filter))
        if outfile == self.outfile:
          outfile = '{}_{}_seed_image.fits'.format(self.outfile, self.filter)
        self.seed_file = os.path.join(self.outdir, outfile)

        xcent_fov = xd / 2
        ycent_fov = yd / 2
        kw = {}
        kw['xcenter'] = xcent_fov
        kw['ycenter'] = ycent_fov
        kw['units'] = units
        kw['TGROUP'] = tgroup
        if self.pupil.upper() == 'F':
            kw['filter'] = self.pupil
        else:
            kw['filter'] = self.filter
        kw['PHOTFLAM'] = self.photflam
        kw['PHOTFNU'] = self.photfnu
        kw['PHOTPLAM'] = self.pivot * 1.e4 #put into angstroms
        kw['NOMXDIM'] = self.subarr_bounds['xend'] - self.subarr_bounds['xstart'] + 1
        kw['NOMYDIM'] = self.subarr_bounds['yend'] - self.subarr_bounds['ystart'] + 1
        kw['NOMXSTRT'] = self.coords['xoffset'] + 1
        kw['NOMXEND'] = self.coords['xoffset'] + kw['NOMXDIM']
        kw['NOMYSTRT'] = self.coords['yoffset'] + 1
        kw['NOMYEND'] = self.coords['yoffset'] + kw['NOMYDIM']
        kw['GRISMPAD'] = self.grism_direct_factor
        self.seedinfo = kw
        self.save_fits(self.seed_image, self.seed_file, key_dict=kw,
                       image2=self.seed_segmap, image2type='SEGMAP')
        print('Blotted image saved to: {}'.format(self.seed_file))

    def save_fits(self, image, output_file, key_dict=None, image2=None,
                       image2type=None):
        """Save the given arrays into a fits file

        Parameters
        ----------
        image : numpy.ndarray
            2D image array

        output_file : str
            Name of the fits file to save the images into

        key_dict : dict
            Dictionary of header keywords to add to the primary and
            first extension headers

        image2 : numpy.ndarray
            Image to save in a second extension of the fits file
            (e.g. segmentation map)

        image2type : str
            Name to use for the extension into which ``image2`` is written
        """
        h0 = fits.PrimaryHDU()
        h1 = fits.ImageHDU(image, name='DATA')
        if image2 is not None:
            h2 = fits.ImageHDU(image2)
            if image2type is not None:
                h2.header['EXTNAME'] = image2type

        # If a keyword dictionary is provided, put the
        # keywords into the 0th and 1st extension headers
        if key_dict is not None:
            for key in key_dict:
                h0.header[key] = key_dict[key]
                h1.header[key] = key_dict[key]

        if image2 is None:
            hdulist = fits.HDUList([h0,h1])
        else:
            hdulist = fits.HDUList([h0,h1,h2])
        hdulist.writeto(output_file, overwrite=True)

    def grism_coord_adjust(self):
        """Calculate the factors by which to expand the
        output array size, as well as the coordinate
        offsets between the nominal output array and
        the input lists if the observation being
        modeled is to be dispersed with the grism
        """
        dtor = radians(1.)

        # Normal imaging with grism image requested
        self.coords['x'] = self.grism_direct_factor
        self.coords['y'] = self.grism_direct_factor
        self.coords['xoffset'] = np.int((self.grism_direct_factor - 1.)
                                        * (self.subarray_bounds[2] -
                                           self.subarray_bounds[0]+1) / 2.)
        self.coords['yoffset'] = np.int((self.grism_direct_factor - 1.)
                                        * (self.subarray_bounds[3] -
                                           self.subarray_bounds[1]+1) / 2.)

    def add_options(self, parser = None, usage = None):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage,\
                        description='Create seed image via catalogs')
        parser.add_argument("paramfile",help='File describing the input parameters and instrument settings to use. (YAML format).')
        parser.add_argument("mosaicfile",help="Fits file containing the image to create the seed image from")
        parser.add_argument("crop_center_ra",help="RA of center of the image to crop from mosaicimage")
        parser.add_argument("crop_center_dec",help="Dec of center of the image to crop from mosaicimage")
        parser.add_argument("channel",help="NIRCam channel: 'short' or 'long'")
        parser.add_argument("detector",help="NIRCam detector: e.g. 'A1','B5'")
        parser.add_argument("blot_center_ra",help="RA of center of the output blotted image")
        parser.add_argument("blot_center_dec",help="Dec of center of the output blotted image")
        parser.add_argument("blot_pav3",help="Position angle of the output blotted image")
        parser.add_argument("aperture",help="Aperture to use. (e.g. NRCA1_FULL)")
        parser.add_argument("subarray_defs",help="Ascii file containing subarray definitions")
        return parser


if __name__ == '__main__':

    usagestring = 'USAGE: fits_seed_image.py inputs.yaml'

    seed = Imgseed()
    parser = seed.add_options(usage = usagestring)
    args = parser.parse_args(namespace=seed)
    seed.read_param_file(args.paramfile)
    seed.crop_and_blot()
