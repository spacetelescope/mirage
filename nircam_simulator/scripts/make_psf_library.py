#! /usr/bin/env

'''
create PSF library files for use by the ramp simulator
'''

import make_nircam_psfs as m
import numpy as np

# Set parameters
wfe_list = ['predicted', 'requirements']
filter_list = ['F070W', 'F090W', 'F115W', 'F140M', 'F150W2', 'F150W', 'F162M',
               'F164N', 'F182M', 'F187N', 'F200W', 'F210M', 'F212N', 'F250M',
               'F277W', 'F300M', 'F322W2', 'F323N', 'F335M', 'F356W', 'F360M',
               'F405N', 'F410M', 'F430M', 'F444W', 'F460M', 'F466N', 'F470N',
               'F480M']
detector_list = ['A3', 'A4', 'B1', 'B2', 'B3', 'B4', 'B5']
xpos_list = [1024]
ypos_list = [1024]
realizations = np.arange(5)
num_subpix_pos = 25
subpix_rows = np.int(np.sqrt(num_subpix_pos))
subpix_offsets = np.linspace(-0.5, 0.5, subpix_rows)

# Create PSFs
lib = m.PSFs()
lib.offsets = subpix_offsets
lib.wfe_list = wfe_list
lib.filter_list = filter_list
lib.detector_list = detector_list
lib.xpos_list = xpos_list
lib.ypos_list = ypos_list
lib.realization_list = realizations
lib.run()
