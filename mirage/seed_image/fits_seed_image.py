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

   subarray_defs -
              Name of the ascii file containing the definitions
              of the possible subarrays. By setting this value
              equal to 'config', the code will use the
              subarray definition file packaged with the
              NIRCam data simulator.

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

import os
import sys
import argparse

import yaml
import pkg_resources
import numpy as np
from math import radians
from astropy.io import fits, ascii
from photutils import detect_sources


from . import crop_mosaic, blot_image

config_files = {'subarray_defs':'NIRCam_subarray_definitions.list',
                'flux_cal':'NIRCam_zeropoints.list'}

class ImgSeed:
    def __init__(self,paramfile = None):
        self.mosaicfile = None
        self.crop_center_ra = 0.
        self.crop_center_dec = 0.
        self.channel = None
        self.detector = None
        self.blot_center_ra = 0.
        self.blot_center_dec = 0.
        self.blot_pav3 = 0.
        self.aperture = ''
        self.subarray_defs = ''
        self.flux_cal_file = ''
        self.filter = ''
        self.pupil = ''
        self.grism_source_image = False
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
        self.coords = {'x':1.,'xoffset':0.,'y':1.,'yoffset':0.}

        # If a paramfile is provided, read in and
        # set params
        if paramfile is not None:
            self.read_param_file(paramfile)



    def read_param_file(self, file):
        '''If an input yaml file is given
        which matches the format of the input
        yaml files for the other mirage
        steps, read in and set needed parameters'''

        if os.path.isfile(file):
            with open(file,'r') as f:
                params = yaml.safe_load(f)
        else:
            print(("WARNING: {} does not exist."
                   .format(file)))
            sys.exit()

        self.aperture = params['Readout']['array_name']
        self.detector = self.aperture[3:5]
        if '5' in self.detector:
            self.channel = 'long'
        else:
            self.channel = 'short'
        #module = self.detector[0].upper()
        self.crop_center_ra = params['Telescope']['ra']
        self.crop_center_dec = params['Telescope']['dec']
        self.blot_center_ra = params['Telescope']['ra']
        self.blot_center_dec = params['Telescope']['dec']
        self.blot_pav3 = params['Telescope']['rotation']
        self.subarray_defs = params['Reffiles']['subarray_defs']
        self.flux_cal_file = params['Reffiles']['flux_cal']
        self.filter = params['Readout']['filter']
        self.pupil = params['Readout']['pupil']
        # If input is to be dispersed later, we need to
        # expand the region to be blotted
        if params['Output']['grism_source_image']:
            self.grism_coord_adjust()
        self.outfile = params['Output']['file']
        self.outdir = params['Output']['directory']


    def expand_config_paths(self,file,filetype):
        if file.lower() == 'config':
            sfile = config_files[filetype]
            return os.path.join(self.modpath, 'config', sfile)


    def fluxcal_info(self, fluxfile, filter, module):
        # Read in list of zeropoints/photflam/photfnu
        self.zps = ascii.read(fluxfile)
        mtch = ((self.zps['Filter'] == filter) \
                & (self.zps['Module'] == module))
        self.photflam = self.zps['PHOTFLAM'][mtch][0]
        self.photfnu = self.zps['PHOTFNU'][mtch][0]
        self.pivot = self.zps['Pivot_wave'][mtch][0]


    def aperture_info(self, apname):
        '''
        Return basic information on the reuqested
        aperture. This information comes from the
        subarray defintion file in the config
        directory
        '''
        # Read in the subarray definition file
        print(self.subarray_defs)
        self.read_subarray_definition_file(self.subarray_defs)

        # Find the matching entry in the definition file
        if self.aperture in self.subdict['AperName']:
            mtch = self.aperture == self.subdict['AperName']
            xstart = self.subdict['xstart'].data[mtch][0]
            ystart = self.subdict['ystart'].data[mtch][0]
            xend = self.subdict['xend'].data[mtch][0]
            yend = self.subdict['yend'].data[mtch][0]
            self.array = {'xstart':xstart,
                               'ystart':ystart,
                               'xend':xend,
                               'yend':yend}
        else:
            print(("WARNING: Array name {} does not "
                   .format(apname)))
            print("to be present in subarray definition file.")
            sys.exit()


    def read_subarray_definition_file(self, file):
        '''Read in the file that contains a list
        of subarray names and positions on the detector
        '''
        try:
            self.subdict = ascii.read(file,data_start=1,header_start=0)
        except:
            print("Error: could not read in subarray definitions file:")
            print(file)
            sys.exit()


    def crop_and_blot(self):
        '''Extract an area of the input fits
        file that is slightly larger than the
        FOV, and use the JWST pipeline version
        of blot to resample to the correct
        pixel scale and introduce the proper
        distortion'''

        # Locations of conifg files
        self.subarray_defs = self.expand_config_paths(self.subarray_defs, 'subarray_defs')
        self.flux_cal_file = self.expand_config_paths(self.flux_cal_file, 'flux_cal')

        if self.detector is None:
            self.detector = self.aperture[3:5]
        if self.channel is None:
            if '5' in self.detector:
                self.channel = 'long'
            else:
                self.channel = 'short'

        # Get information on the aperture
        self.aperture_info(self.aperture)
        xstart = self.array['xstart']
        ystart = self.array['ystart']
        xdim = self.array['xend'] - xstart + 1
        ydim = self.array['yend'] - ystart + 1
        module = self.detector[0].upper()

        # Flux calibration information
        if self.pupil[0].upper() == 'F':
            usefilt = self.pupil.upper()
        else:
            usefilt = self.filter.upper()
        self.fluxcal_info(self.flux_cal_file, usefilt, module)

        # Extract
        crop = crop_mosaic.Extraction()
        crop.mosaicfile = self.mosaicfile
        crop.center_ra = self.crop_center_ra
        crop.center_dec = self.crop_center_dec
        crop.channel = self.channel
        crop.corner = (xstart, ystart)
        crop.dimensions = (xdim, ydim)
        crop.extract()

        # Blot
        blot = blot_image.Blot()
        blot.blotfile = crop.cropped
        blot.detector = [self.detector]
        blot.center_ra = [self.blot_center_ra]
        blot.center_dec = [self.blot_center_dec]
        blot.pav3 = [self.blot_pav3]
        blot.blot()

        for dm in blot.blotted_datamodels:
            # Create segmentation map
            # (used only when dispsersing image for WFSS)
            self.seed_segmap = self.make_segmap(dm)

            # Save
            self.save_seed_image(dm)


    def make_segmap(self, model):
        '''Create a segmentation map of the
        input image'''
        noise = np.median(model.data)
        map = detect_sources(model.data,noise*3.,8).data
        return map


    def save_seed_image(self, model):
        '''Save the seed image in the proper
        format that it can be used in subsequent
        steps of mirage'''
        self.seed_image = model.data
        dim = len(self.seed_image.shape)
        if dim == 2:
            units = 'e-/sec'
            yd,xd = arrayshape
            tgroup = 0.
        else:
            print(("WARNING: {}D seed images from "
                   "input fits file not yet supported."
                   .format(dim)))
            sys.exit()

        dot = self.outfile.rfind('.fits')
        if dot == -1:
            dot = len(self.outfile)
        self.basename = os.path.join(self.outdir,
                                     self.outfile[0:dot])
        self.seed_file = self.basename + '_' \
                      + self.filter + '_seed_image.fits'

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
        kw['NOMXDIM'] = self.array['xend'] - self.array['xstart'] + 1
        kw['NOMYDIM'] = self.array['yend'] - self.array['ystart'] + 1
        kw['NOMXSTRT'] = self.coords['xoffset'] + 1
        kw['NOMXEND'] = self.coords['xoffset'] + kw['NOMXDIM']
        kw['NOMYSTRT'] = self.coords['yoffset'] + 1
        kw['NOMYEND'] = self.coords['yoffset'] + kw['NOMYDIM']
        kw['GRISMPAD'] = self.grism_direct_factor
        self.seedinfo = kw
        self.save_fits(self.seedimage, self.seed_file, key_dict = kw,
                       image2 = self.seed_segmap, image2type = 'SEGMAP')


    def saveSingleFits(self, image, name, key_dict = None,
                       image2 = None, image2type = None):
        # Save an array into the first extension of a fits file
        h0 = fits.PrimaryHDU()
        h1 = fits.ImageHDU(image,name='DATA')
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
        hdulist.writeto(name,overwrite=True)


    def grism_coord_adjust(self):
        # Calculate the factors by which to expand the
        # output array size, as well as the coordinate
        # offsets between the nominal output array and
        # the input lists if the observation being
        # modeled is to be dispersed with the grism
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

    seed = ImgSeed()
    parser = seed.add_options(usage = usagestring)
    args = parser.parse_args(namespace=seed)
    seed.read_param_file(args.paramfile)
    seed.crop_and_blot()
