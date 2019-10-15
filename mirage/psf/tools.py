#! /usr/bin/env python

"""
Tools for working with PSFs
"""
import os

from astropy.io import fits
from astropy.modeling.models import Gaussian2D
import numpy as np
from photutils.datasets import load_irac_psf
import webbpsf

from mirage.psf.psf_selection import get_psf_wings
from mirage.utils.siaf_interface import get_instance


def gaussian_psf(fwhm, xdim, ydim):
    """Return a 2D Gaussian kernel

    Parameters
    ----------
    fwhm : float
        FWHM of the Gaussian kernel to return

    xdim : int
        Length of the created PSF in the x-dimension

    ydim : int
        Length of the created PSF in the y-dimension

    Returns
    -------
    kernel : numpy.ndarray
        2D array containing 2D Gaussian kernel
    """
    # Center the PSF in the array
    xmean = xdim / 2.
    ymean = ydim / 2.

    # Translate FWHM to standard deviation
    dev = fwhm / 2.35

    # Create model
    model = Gaussian2D(amplitude=1., x_mean=xmean, y_mean=ymean, x_stddev=dev, y_stddev=dev)

    # Evalulate model
    x_grid, y_grid = np.mgrid[0: xdim, 0: ydim]
    kernel = model.evaluate(x=x_grid, y=y_grid, amplitude=1., x_mean=xmean,
                            y_mean=ymean, x_stddev=dev, y_stddev=dev, theta=0.)
    return kernel


def get_HST_PSF(meta):
    """Get a representative PSF for the HST instrument/detector in
    question
    """
    if meta['instrument'] == 'WFC3':
        if meta['filter'] == 'F160W':
            grid_size = something
            y, x = np.mgrid[0:51, 0:51]
            gm1 = Gaussian2D(100, 25, 25, 3, 3)


def get_IRAC_PSF(meta):
    """Use photutils to get a representative PSF for the IRAC data
    specified in ``meta``.

    Parameters
    ----------
    meta : dict
        Dictionary of instrument information. Most important is
        meta['channel'] which contains the integer IRAC channel
        number. PSFs are loaded from photutils

    Returns
    -------
    hdu.data : numpy.ndarray
        2D array containing the IRAC PSF
    """
    if isinstance(meta['channel'], int):
        if meta['channel'] >= 1 and meta['channel'] <= 4:
            hdu = load_irac_psf(channel=meta['channel'])
        else:
            raise ValueError("ERROR: IRAC Channel number must be an integer from 1 to 4.")
    else:
        raise ValueError("ERROR: IRAC Channel number must be an integer from 1 to 4.")

    pixscale = hdu.header['SECPIX']
    mosaic_pix_scale1 = meta['pix_scale1']
    if not np.isclose(mosaic_pix_scale1, pixscale, rtol=0.05):
        raise ValueError(("Reported pixel scale in the mosaic does not match that in "
                          "the IRAC PSF from photutils. If the pixel scale of the mosaic "
                          "is really not {}, please provide your own fits file containing "
                          "a representative PSF at the desired pixel scale.".format(pixscale)))

    return hdu.data


def get_JWST_pixel_scale(meta, aperture=None):
    """Use pysiaf to find the nominal pixel scale for a given JWST
    instrument/detector

    Parameters
    ----------
    meta : dict
        Dictionary of instrument information. Most important is
        meta['channel'] which contains the integer IRAC channel
        number. PSFs are loaded from photutils

    aperture : str
        Name of aperture to look up pixel scale for. If None,
        then fall back to using FULL or CEN as appropriate.

    Returns
    -------
    siaf_aper.XSciScale : float
        Pixel scale (arcsec/pixel) in x-direction on detector

    siaf_aper.YSciScale : float
        Pixel scale (arcsec/pixel) in y-direction on detector
    """
    # Create SIAF instance
    siaf = get_instance(meta['instrument'])

    # If needed, get the name to append to the detector to create
    # the full aperture
    # name
    if aperture is None:
        if meta['instrument'].lower() == 'nircam':
            aperture = 'FULL'
        elif meta['instrument'].lower() == 'niriss':
            aperture = 'CEN'
        elif meta['instrument'].lower() == 'fgs':
            aperture = 'FULL'

        # In SIAF, the FGS detectors are FGS1,2 rather than GUIDER1,2
        detector = meta['detector'].upper()
        if 'GUIDER' in detector:
            detector = detector.replace('GUIDER', 'FGS')

        aperture = '{}_{}'.format(detector, aperture)

    # Get the aperture-specific SIAF info
    siaf_aper = siaf[aperture]

    # Return the x and y pixel scale
    return siaf_aper.XSciScale, siaf_aper.YSciScale


def get_JWST_PSF(meta):
    """Use webbpsf to get a representative PSF for the JWST data specified
    in ``meta``.

    Parameters
    ----------
    meta : dict
        Dictionary of instrument information. Most important is
        meta['channel'] which contains the integer IRAC channel
        number. PSFs are loaded from photutils

    Returns
    -------
    psf_model : numpy.ndarray
        2D array containing the PSF
    """
    # Find the nominal pixel scale, which is what the psf_wings
    # library is saved at
    nominal_pix_scale_x, nominal_pix_scale_y = get_JWST_pixel_scale(meta)

    # If the pixel scale reported in the mosaic file is the same as the
    # nominal pixel scale for that instrument/detector, then just read in
    # the PSF from the psf_wings file
    if meta['pix_scale1'] is not None:
        if np.isclose(meta['pix_scale1'], nominal_pix_scale_x, atol=0., rtol=0.05):
            instrument = meta['instrument'].lower()
            psf_model = get_psf_wings(meta['instrument'], meta['detector'], meta['filter'], meta['pupil'],
                                      'predicted', 0, os.path.join(os.path.expandvars('$MIRAGE_DATA'),
                                                                   instrument, 'gridded_psf_library/psf_wings'))
        else:
            # If the mosaic's reported pixel scale does not agree with the
            # nominal pixel scale (such as if the data have been drizzled)
            # then call webbpsf and create a PSF with the appropriate scale
            inst = meta['instrument'].upper()
            if inst == 'NIRCAM':
                psf = webbpsf.NIRCam()
            elif inst == 'NIRISS':
                psf = webbpsf.NIRISS()
            elif inst == 'FGS':
                psf = webbpsf.FGS()

            psf.detector = meta['detector'].upper()
            try:
                psf.filter = meta['filter'].upper()
            except ValueError:
                try:
                    psf.filter = meta['pupil'].upper()
                except ValueError:
                    raise ValueError("ERROR: Neither {} nor {} are valid filter values for {} in WebbPSF."
                                     .format(meta['filter'], meta['pupil'], inst))

            opd_list = psf.opd_list
            predicted_opd = [opd for opd in opd_list if 'predicted' in opd][0]
            psf.pupilopd = (predicted_opd, 0)
            oversample = meta['pix_scale1'] / nominal_pix_scale_x
            something = psf.calc_psf(oversample=oversample, fov_arcsec=999)
            print('this is not done yet')

    return psf_model


def get_psf_metadata(filename):
    """Retrieve basic PSF-related metadata

    Parameters
    ----------
    filename : str
        Name of fits file containing mosaic

    Returns
    -------
    metadata : dict
        Dictionary of instrument-specific data
    """
    metadata = {}
    with fits.open(filename) as hdulist:

        # We want basic metadata on the instrument. Assume it's in the
        # primary header
        header = hdulist[0].header

        try:
            telescope = header['TELESCOP'].upper()
        except KeyError:
            telescope = None
        metadata['telescope'] = telescope

        try:
            instrument = header['INSTRUME'].upper()
        except KeyError:
            instrument = None
        metadata['instrument'] = instrument

        if telescope == 'JWST':
            try:
                metadata['detector'] = header['DETECTOR']
                metadata['filter'] = header['FILTER']
                metadata['pupil'] = header['PUPIL']
                metadata['pix_scale1'] = np.abs(header['CD1_1']) * 3600.
                metadata['pix_scale2'] = np.abs(header['CD2_2']) * 3600.
            except KeyError:
                metadata['detector'] = None
                metadata['filter'] = None
                metadata['pupil'] = None

        if telescope == 'HST':
            try:
                metadata['detector'] = header['DETECTOR']
                metadata['filter'] = header['FILTER']
                metadata['pa'] = header['PA_APER']  # PA of reference aperture center
                metadata['pix_scale1'] = np.abs(header['CD1_1']) * 3600.
                metadata['pix_scale2'] = np.abs(header['CD2_2']) * 3600.
            except KeyError:
                metadata['detector'] = None
                metadata['filter'] = None
                metadata['pa'] = None
                metadata['pix_scale1'] = None
                metadata['pix_scale2'] = None

        if instrument == 'IRAC':
            try:
                metadata['channel'] = int(header['CHNLNUM'])
                metadata['pix_scale1'] = np.abs(header['PXSCAL1'])
                metadata['pix_scale2'] = np.abs(header['PXSCAL2'])
                metadata['pa'] = header['PA']  #  [deg] Position angle of axis 2 (E of N)
            except KeyError:
                metadata['channel'] = None
                metadata['pix_scale1'] = None
                metadata['pix_scale2'] = None
                metadata['pa'] = None

        if telescope not in 'JWST HST SPITZER'.split():
            try:
                metadata['pix_scale1'] = np.abs(header['CD1_1']) * 3600.
                metadata['pix_scale2'] = np.abs(header['CD2_2']) * 3600.
            except KeyError:
                metadata['pix_scale1'] = None
                metadata['pix_scale2'] = None
    return metadata


def same_array_size(array1, array2):
    """Crop the input arrays such that they are the same size in both
    dimensions

    Parameters
    ----------
    array1 : numpy.ndarray
        2D array

    array2 : numpy.ndarray
        2D array

    Returns
    -------
    array1 : numpy.ndarray
        Potentially cropped 2D array

    array2 : numpy.ndarray
        Potentially cropped 2D array

    """
    shape1 = array1.shape
    shape2 = array2.shape
    miny = min([shape1[0], shape2[0]])
    minx = min([shape1[1], shape2[1]])

    # Crop y dimension
    if shape1[0] == miny:
        smaller = array1
        larger = array2
        smaller_shape = shape1
        larger_shape = shape2
    else:
        smaller = array2
        larger = array1
        smaller_shape = shape2
        larger_shape = shape1

    diff = abs(larger_shape[0] - miny)
    low_delta = diff // 2
    if diff % 2 == 1:
        high_delta = low_delta + 1
    else:
        high_delta = low_delta
    larger = larger[low_delta: larger_shape[0]-high_delta, :]

    if shape1[0] == miny:
        array1 = smaller
        array2 = larger
    else:
        array1 = larger
        array2 = smaller
    shape1 = array1.shape
    shape2 = array2.shape

    # Crop x dimension
    if shape1[1] == minx:
        smaller = array1
        larger = array2
        smaller_shape = shape1
        larger_shape = shape2
    else:
        smaller = array2
        larger = array1
        smaller_shape = shape2
        larger_shape = shape1

    diff = abs(larger_shape[1] - minx)
    low_delta = diff // 2
    if diff % 2 == 1:
        high_delta = low_delta + 1
    else:
        high_delta = low_delta
    larger = larger[:, low_delta: larger_shape[1]-high_delta]

    if shape1[1] == minx:
        array1 = smaller
        array2 = larger
    else:
        array1 = larger
        array2 = smaller
    shape1 = array1.shape
    shape2 = array2.shape

    return array1, array2
