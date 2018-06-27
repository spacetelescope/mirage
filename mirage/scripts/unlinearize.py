#! /usr/bin/env python

'''
Run the linearity correction step backwards on an input
image or ramp. Introduce non-linearity into a linear
input ramp.
'''

import sys
import numpy as np

def unlinearize(image,coeffs,sat,lin_satmap,maxiter=10,accuracy=0.000001,robberto=False
                ,save_accuracy_map=False,accuracy_file='unlinearize_no_convergence.fits'):
    # Insert non-linearity into the linear synthetic sources

    # If the sizes of the satmap or coeffs are different than 
    # the data, return an error
    if sat.shape != image.shape[-2:]:
        print("WARNING: image y,x shape is {},".format(image.shape[-2:]))
        print("but input saturation map shape is {}".format(sat.shape))
        sys.exit()

    # REWORK SO THAT THE LINARIZED SATURATION MAP IS AN INPUT
    # RATHER THAN BEING CREATED HERE. THIS IS BECAUSE THE SUPERBIAS
    # AND REFPIX SIGNALS MUST BE SUBTRACTED BEFORE LINEARIZING
    # THE ORIGINAL RAW SATURATION MAP
    # Sat is the original saturation map data for non-linear ramps. 
    # Translate to saturation maps for the linear ramps here so 
    # that we can pay attention only to non-saturated pixels in
    # the input linear image. Do this by applying the standard
    # linearity correction to the saturation map
    #lin_satmap = nonLinFunc(sat,coeffs,sat)

    # Find pixels with "good" signals, to have the nonlin applied. 
    # Negative pix or pix with signals above the requested max 
    # value will not be changed.
    x = np.copy(image)
    i1 = np.where((image > 0.) & (image < lin_satmap))
    dev = np.zeros_like(image,dtype=float)
    dev[i1] = 1.
    i2 = np.where((image <= 0.) | (image >= lin_satmap))
    numhigh = np.where(image >= lin_satmap)
    i = 0

    # Initial run of the nonlin function - when calling the 
    # non-lin function, give the original satmap for the 
    # non-linear signal values
    val = nonLinFunc(image,coeffs,sat)
    val[i2]=1.

    if robberto:            
        x = image * val
    else:
        x[i1] = (image[i1]+image[i1]/val[i1]) / 2. 
        while i < maxiter:
            i=i+1
            val = nonLinFunc(x,coeffs,sat)
            val[i2]=1.
            dev[i1] = abs(image[i1]/val[i1]-1.)
            inds = np.where(dev[i1] > accuracy)
            if inds[0].size < 1:
                break
            val1 = nonLinDeriv(x,coeffs,sat)
            val1[i2] = 1.
            x[i1] = x[i1] + (image[i1]-val[i1])/val1[i1]

    # If we max out the number of iterations, 
    # save the array of accuracy values. Spot 
    # checks reveal the pix that don't meet the 
    # accuracy reqs are randomly located on the detector, 
    # and don't seem to be correlated with point source
    # locations.
    if i == maxiter and save_accuracy_map:
        from astropy.io import fits
        print(("WARNING: some pixels failed to unlinearize correctly within "
               "the maximum number of iterations. Map of accuracy of the "
               "unlinearized values saved to {}.".format(accuracy_file)))
        devcheck = np.copy(dev)
        devcheck[i2] = -1.
        h0 = fits.PrimaryHDU()
        h1 = fits.ImageHDU(devcheck)
        hl = fits.HDUList([h0,h1])
        hl.writeto(accuracy_file,overwrite=True)           
    return x


def nonLinFunc(image,coeffs,limits):
    # Apply linearity correction coefficients
    # to image.
    values = np.copy(image)
    bady = 0
    badx = 1
    if len(image.shape) == 3:
        bady = 1
        badx = 2
    elif len(image.shape) == 4:
        bady = 2
        badx = 3

    bad = np.where(values > limits)
    values[bad] = limits[bad[bady],bad[badx]]
    ncoeff = coeffs.shape[0]
    t = np.copy(coeffs[-1,:,:])
    for i in range(ncoeff-2,-1,-1):
        t = coeffs[i,:,:] + values*t
    return t


def nonLinDeriv(image,coeffs,limits):
    # First derivative of non-lin correction
    # function
    values = np.copy(image)
    bady = 0
    badx = 1
    if len(image.shape) == 3:
        bady = 1
        badx = 2
        
    bad = np.where(values > limits)
    values[bad] = limits[bad[bady],bad[badx]]
    ncoeff = coeffs.shape[0]
    t = (ncoeff-1) * np.copy(coeffs[-1,:,:])
    for i in range(ncoeff-3,-1,-1):
        t = (i+1) * coeffs[i+1,:,:] + values*t
    return t
