#! /usr/bin/env python

'''
Quick test of the imaging_simulator.py module
'''

from nircam_simulator.scripts import imaging_simulator

m = imaging_simulator.ImgSim()
m.paramfile = 'V42424001001P0000000001104_B3_F115W.yaml'
#m.override_dark = 'V42424001001P0000000001104_B3_F115W_uncal_linear_dark_prep_object.fits'
m.create()
