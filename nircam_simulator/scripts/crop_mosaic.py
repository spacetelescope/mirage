#! /usr/bin/env python

'''
Extract a roughly SCA-sized area from a distortion-free mosaic.
This extracted area can then be blotted to match NIRCam distorion
and then used as an input to the ramp simulator.
'''

from astropy.io import fits
import argparse
import numpy as np

class Extraction():

    def __init__(self):
        self.mosaicfile = ''
        self.center_ra = 0.0
        self.center_dec = 0.0
        self.channel = 'short'
        self.outfile = None
        self.nrc_scale = {'short':0.031,
                          'long':0.063}


    def extract(self):
        # get WCS info from mosaic file
        mosaic = fits.open(self.mosaicfile)
        mosaic_ra_ref = mosaic[0].header['CRVAL1']
        mosaic_dec_ref = mosaic[0].header['CRVAL2']
        mosaic_x_ref = mosaic[0].header['CRPIX1']
        mosaic_y_ref = mosaic[0].header['CRPIX2']
        self.mosaic_x_type = mosaic[0].header['CTYPE1']
        self.mosaic_y_type = mosaic[0].header['CTYPE2']
        self.mosaic_scale_x = mosaic[0].header['CDELT1']*3600. #arcsec/pix
        self.mosaic_scale_y = mosaic[0].header['CDELT2']*3600. #arcsec/pix
        #self.mosaic_roll = 0. #North up
        self.cd11 = mosaic[0].header['CD1_1']
        self.cd12 = mosaic[0].header['CD1_2']
        self.cd21 = mosaic[0].header['CD2_1']
        self.cd22 = mosaic[0].header['CD2_2']
        self.mosaic_roll = np.arccos(self.cd11/(self.mosaic_scale_x/3600.))

        # define size of region to extract
        outscale_x = self.nrc_scale[self.channel.lower()]
        
        # This assumes full frame!!! Need some buffer to account for distortion
        nominal_size = 2897. + 300. #300 pixel buffer
        self.nx = np.absolute(np.int(nominal_size * outscale_x / self.mosaic_scale_x))
        self.ny = np.absolute(np.int(nominal_size * outscale_x / self.mosaic_scale_y))

        # Need to calculate what pixel in the mosaic corresponds to the
        # RA, Dec of the center of the cutout (center_ra,center_dec)
        # For the moment, assume mosaic is north up, which Anton's are...
        delta_ra = mosaic_ra_ref - self.center_ra
        delta_dec = mosaic_dec_ref - self.center_dec

        delta_x = delta_ra * self.mosaic_scale_x
        delta_y = delta_dec * self.mosaic_scale_y
        
        centerx = mosaic_x_ref + delta_x
        centery = mosaic_y_ref + delta_y

        intcenterx = np.int(centerx)
        intcentery = np.int(centery)

        # extract the SCA-sized area from the moasic
        ny_over_2 = np.int(self.ny / 2)
        nx_over_2 = np.int(self.nx / 2)
        print("Coords of center of cropped area",intcenterx,intcentery)
        print("X-min, Y-max coords:",intcenterx-nx_over_2,intcenterx+nx_over_2)
        print("Y-min, Y-max coords:",intcentery-ny_over_2,intcentery+ny_over_2)
        
        crop = mosaic[0].data[intcentery-ny_over_2:intcentery+ny_over_2,intcenterx-nx_over_2:intcenterx+nx_over_2]

        # save extracted area - tests show saving with fits is much much
        # faster than saving with a datamodel. So save with fits, and then
        # place into a datamodel to prepare for later blotting
        #self.saveDMS(crop,self.outfile)

        # Place into a data model to prepare for blotting
        self.cropped = self.populate_datamodel(crop)
        
        # Save only if outfile is not None
        if self.outfile is not None:
            #self.outfile = 'cropped_from_'+self.mosaicfile+'_ra_'+str(self.center_ra)+'_dec_'+str(self.center_dec)+'.fits'
            self.savefits(crop,self.outfile)
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
        # These assume north up orientation
        #h0.header['PC1_1'] = self.mosaic_scale_x / 3600.
        #h0.header['PC1_2'] = 0.
        #h0.header['PC2_1'] = 0.
        #h0.header['PC2_2'] = self.mosaic_scale_y / 3600.

        h1.header['PC1_1'] = self.cd11 / (self.mosaic_scale_x/3600.)
        h1.header['PC1_2'] = self.cd12 / (self.mosaic_scale_x/3600.)
        h1.header['PC2_1'] = self.cd21 / (self.mosaic_scale_y/3600.)
        h1.header['PC2_2'] = self.cd22 / (self.mosaic_scale_y/3600.)

        #for key in h1.header:
        #    print(key,h1.header[key])

        hlist = fits.HDUList([h0,h1])
        hlist.writeto(ofile,overwrite=True)


    def populate_datamodel(self,array):
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

        # These assume north up orientation
        #cropped.meta.wcsinfo.pc1_1 = self.mosaic_scale_x / 3600.
        #cropped.meta.wcsinfo.pc1_2 = 0.
        #cropped.meta.wcsinfo.pc2_1 = 0.
        #cropped.meta.wcsinfo.pc2_2 = self.mosaic_scale_y / 3600.

        cropped.meta.wcsinfo.pc1_1 = self.cd11 / (self.mosaic_scale_x/3600.)
        cropped.meta.wcsinfo.pc1_2 = self.cd12 / (self.mosaic_scale_x/3600.)
        cropped.meta.wcsinfo.pc2_1 = self.cd21 / (self.mosaic_scale_y/3600.)
        cropped.meta.wcsinfo.pc2_2 = self.cd22 / (self.mosaic_scale_y/3600.)

        cropped.data = array

        #save the cropped view
        #cropped.save(ofile)
        return cropped


    def add_options(self,parser=None,usage=None):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage,description='Extract SCA-sized area from moasic')
        parser.add_argument("mosaicfile",help="Filename of fits file containing mosaic.")
        parser.add_argument("center_ra",help="RA at the center of the extracted area",type=np.float)
        parser.add_argument("center_dec",help="Dec at the center of the extracted area",type=np.float)
        parser.add_argument("--channel",help="NIRCam channel the extracted data are meant to simulate (short, long)",default='short')
        parser.add_argument("--outfile",help="Name of output fits file containing extracted image",default=None)
        return parser

if __name__ == '__main__':
    usagestring = 'USAGE: crop_mosaic.py filename.fits 23.43 -21.2'

    ex = Extraction()
    parser = ex.add_options(usage = usagestring)
    args = parser.parse_args(namespace=ex)
    ex.extract()

