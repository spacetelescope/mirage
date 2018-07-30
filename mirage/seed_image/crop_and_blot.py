#! /usr/bin/env python

'''
test to combine crop_mosaic.py and blot_image.py
'''

from . import crop_mosaic, blot_image

crop = crop_mosaic.Extraction()
crop.mosaicfile = 'goodss_candels_udf_wfc3_f105w_010mas_test_drz.fits'
crop.center_ra = 53.122751
crop.center_dec = -27.805089
crop.channel = 'short'
crop.extract()

#crop.cropped now has the datamodel with the mosaic piece...

blot = blot_image.Blot()
blot.blotfile = crop.cropped
blot.detector = ['A1','A1']
blot.center_ra = [53.122751,53.122751]
blot.center_dec = [-27.805089,-27.805089]
blot.pav3 = [0.,12.]
blot.blot()

