#! /usr/bin/env python

"""
Extract a roughly SCA-sized area from a distortion-free mosaic.
This extracted area can then be blotted to match NIRCam distorion
and then used as an input to the ramp simulator.
"""

import argparse

import numpy as np
from astropy import wcs
from astropy.io import fits

class Extraction():
    """
    Parameters
    ----------
    data_extension_number : int
        Extension number within the fits file where the data are stored

    wcs_extension_number : int
        Extension number within the fits file where the WCS data are
        stored in the header

    center_ra : float
        RA value in decimal degrees at the cropped image's reference location

    center_dec : float
        Dec value in decimal degrees at the cropped image's reference location
    """



    def __init__(self, mosaicfile=None, data_extension_number=0, wcs_extension_number=0,
                 center_ra=None, center_dec=None, outfile=None, dimensions=(0, 0),
                 jwst_pixel_scale=None):
        self.mosaicfile = mosaicfile
        self.data_extension_number = data_extension_number
        self.wcs_extension_number = wcs_extension_number
        self.center_ra = center_ra
        self.center_dec = center_dec
        #self.channel = 'short'
        self.outfile = outfile
        #self.nrc_scale = {'short':0.031,
        #                  'long':0.063}
        self.dimensions = dimensions
        self.jwst_pixel_scale = jwst_pixel_scale

        # Buffer around the extraction dimensions to include.
        # This is necessary in case the cropped image will be
        # rotated during later resampling, in order to avoid
        # having the rotated image corners fall outside the
        # extracted image. This is a multiplicative value.
        self.buffer = np.sqrt(2.)


    def extract(self):
        # get WCS info from mosaic file
        mosaic = fits.open(self.mosaicfile)
        mosaic_wcs = wcs.WCS(mosaic[0].header)
        mosaic_shape = mosaic[0].data.shape

        # Assume the input mosaic has been drizzled, so remove
        # any SIP coefficients that may be present
        print('\n\n******************************************************************')
        print('Assuming input files have been drizzled. Removing SIP coefficients')
        print('******************************************************************\n\n')
        mosaic_wcs.sip = None

        radec = np.array([[self.center_ra, self.center_dec]])
        mosaic_xy_at_center_radec = mosaic_wcs.wcs_world2pix(radec, 1)
        mosaic_center_x = mosaic_xy_at_center_radec[0][0] - 1
        mosaic_center_y = mosaic_xy_at_center_radec[0][1] - 1

        # These are only for populating the header of the cropped file
        #mosaic_ra_ref = mosaic[0].header['CRVAL1']
        #mosaic_dec_ref = mosaic[0].header['CRVAL2']

        mosaic_ra_ref, mosaic_dec_ref = mosaic_wcs.wcs.crval
        mosaic_x_ref, mosaic_y_ref = mosaic_wcs.wcs.crpix
        self.mosaic_x_type, self.mosaic_y_type = mosaic_wcs.wcs.ctype
        self.cd11, self.cd12 = mosaic_wcs.wcs.cd[0]
        self.cd21, self.cd22 = mosaic_wcs.wcs.cd[1]
        self.mosaic_scale_x = np.abs(self.cd11) * 3600.  # arcsec/pix
        self.mosaic_scale_y = np.abs(self.cd22) * 3600.  # arcsec/pix

        #mosaic_x_ref = mosaic[0].header['CRPIX1']
        #mosaic_y_ref = mosaic[0].header['CRPIX2']
        #self.mosaic_x_type = mosaic[0].header['CTYPE1']
        #self.mosaic_y_type = mosaic[0].header['CTYPE2']
        #self.mosaic_scale_x = mosaic[0].header['CDELT1']*3600. #arcsec/pix
        #self.mosaic_scale_y = mosaic[0].header['CDELT2']*3600. #arcsec/pix

        #self.cd11 = mosaic[0].header['CD1_1']
        #self.cd12 = mosaic[0].header['CD1_2']
        #self.cd21 = mosaic[0].header['CD2_1']
        #self.cd22 = mosaic[0].header['CD2_2']
        self.mosaic_roll = np.arccos(self.cd11 / (self.mosaic_scale_x / 3600.))

        # Define size of region to extract
        # self.dimensions are in units of jwst pixwls
        #outscale_x = self.nrc_scale[self.channel.lower()]
        xlen, ylen = self.dimensions

        # Expand by the buffer
        xlen *= self.buffer
        ylen *= self.buffer

        # Determine the dimensions of the mosaic aperture to be cropped in
        # units of mosaic pixels
        self.nx = np.absolute(np.int(xlen * self.jwst_pixel_scale / self.mosaic_scale_x))
        self.ny = np.absolute(np.int(ylen * self.jwst_pixel_scale / self.mosaic_scale_y))

        # Need to calculate what pixel in the mosaic corresponds to the
        # RA, Dec of the center of the cutout (center_ra, center_dec)
        # For the moment, assume mosaic is north up, which Anton's are...
        #delta_ra = mosaic_ra_ref - self.center_ra
        #delta_dec = mosaic_dec_ref - self.center_dec

        #delta_x = delta_ra * self.mosaic_scale_x
        #delta_y = delta_dec * self.mosaic_scale_y

        #centerx = mosaic_x_ref + delta_x
        #centery = mosaic_y_ref + delta_y

        #intcenterx = np.int(centerx)
        #intcentery = np.int(centery)

        # Extract the SCA-sized area from the moasic
        half_height = self.ny // 2 + 1
        half_width = self.nx // 2 + 1

        miny = np.int(mosaic_center_y - half_height)
        maxy = np.int(mosaic_center_y + half_height + 1)
        minx = np.int(mosaic_center_x - half_width)
        maxx = np.int(mosaic_center_x + half_width + 1)
        if miny < 0:
            miny = 0
        if minx < 0:
            minx = 0
        if maxy > mosaic_shape[0]:
            maxy = mosaic_shape[0]
        if maxx > mosaic_shape[1]:
            maxx = mosaic_shape[1]

        print("Coords of center of cropped area", mosaic_center_x, mosaic_center_y)
        print("X-min, X-max coords: ", minx, maxx)
        print("Y-min, Y-max coords: ", miny, maxy)

        crop = mosaic[0].data[miny: maxy, minx : maxx]

        # Place into a data model to prepare for blotting
        self.cropped = self.populate_datamodel(crop)

        # Save only if outfile is not None
        if self.outfile is not None:
            self.savefits(crop, self.outfile)
            print("Extracted image saved to {}".format(self.outfile))

    def savefits(self,array,ofile):
        h0 = fits.PrimaryHDU()
        h1 = fits.ImageHDU()
        h1.data = array

        h1.header['CDELT1'] = self.mosaic_scale_x / 3600.
        h1.header['CDELT2'] = self.mosaic_scale_y / 3600.
        h1.header['CRPIX1'] = self.nx / 2 + 0.5
        h1.header['CRPIX2'] = self.ny / 2 + 0.5
        h1.header['CRVAL1'] = self.center_ra
        h1.header['CRVAL2'] = self.center_dec
        h1.header['CTYPE1'] = self.mosaic_x_type
        h1.header['CTYPE2'] = self.mosaic_y_type
        h1.header['CUNIT1'] = 'deg'
        h1.header['CUNIT2'] = 'deg'
        h1.header['RA_REF'] = self.center_ra
        h1.header['DEC_REF'] = self.center_dec
        h1.header['ROLL_REF'] = self.mosaic_roll * 180./np.pi
        h1.header['WCSAXES'] = 2
        h1.header['EXTNAME'] = 'SCI    '
        h1.header['EXTVER'] = 1
        h1.header['PC1_1'] = self.cd11 / (self.mosaic_scale_x / 3600.)
        h1.header['PC1_2'] = self.cd12 / (self.mosaic_scale_x / 3600.)
        h1.header['PC2_1'] = self.cd21 / (self.mosaic_scale_y / 3600.)
        h1.header['PC2_2'] = self.cd22 / (self.mosaic_scale_y / 3600.)

        hlist = fits.HDUList([h0,h1])
        hlist.writeto(ofile,overwrite=True)


    def populate_datamodel(self, array):
        from jwst import datamodels

        # now we need to update the WCS of the mosaic piece
        cropped = datamodels.ImageModel()
        cropped.meta.wcsinfo.cdelt1 = self.mosaic_scale_x / 3600.
        cropped.meta.wcsinfo.cdelt2 = self.mosaic_scale_y / 3600.
        cropped.meta.wcsinfo.crpix1 = self.nx / 2 + 0.5
        cropped.meta.wcsinfo.crpix2 = self.ny / 2 + 0.5
        cropped.meta.wcsinfo.crval1 = self.center_ra
        cropped.meta.wcsinfo.crval2 = self.center_dec
        cropped.meta.wcsinfo.ctype1 = self.mosaic_x_type
        cropped.meta.wcsinfo.ctype2 = self.mosaic_y_type
        cropped.meta.wcsinfo.cunit1 = 'deg'
        cropped.meta.wcsinfo.cunit2 = 'deg'
        cropped.meta.wcsinfo.ra_ref = self.center_ra
        cropped.meta.wcsinfo.dec_ref = self.center_dec
        cropped.meta.wcsinfo.roll_ref = self.mosaic_roll * 180./np.pi
        cropped.meta.wcsinfo.wcsaxes = 2
        cropped.meta.wcsinfo.pc1_1 = self.cd11 / (self.mosaic_scale_x / 3600.)
        cropped.meta.wcsinfo.pc1_2 = self.cd12 / (self.mosaic_scale_x / 3600.)
        cropped.meta.wcsinfo.pc2_1 = self.cd21 / (self.mosaic_scale_y / 3600.)
        cropped.meta.wcsinfo.pc2_2 = self.cd22 / (self.mosaic_scale_y / 3600.)

        cropped.data = array
        return cropped


    def add_options(self,parser=None,usage=None):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage,description='Extract SCA-sized area from moasic')
        parser.add_argument("mosaicfile",help="Filename of fits file containing mosaic.")
        parser.add_argument("center_ra",help="RA at the center of the extracted area",type=np.float)
        parser.add_argument("center_dec",help="Dec at the center of the extracted area",type=np.float)
        parser.add_argument("--channel",help="NIRCam channel the extracted data are meant to simulate (short, long)",default='short')
        parser.add_argument("--outfile",help="Name of output fits file containing extracted image",default=None)
        parser.add_argument("--dimensions",help="Tuple of (x size, y size) of cropped image")
        parser.add_argument("--buffer",help="Width of buffer (in NIRCam pixels) to add around area to be cropped.",default=300)
        return parser

if __name__ == '__main__':
    usagestring = 'USAGE: crop_mosaic.py filename.fits 23.43 -21.2'

    ex = Extraction()
    parser = ex.add_options(usage = usagestring)
    args = parser.parse_args(namespace=ex)
    ex.extract()

