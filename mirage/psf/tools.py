#! /usr/bin/env python

"""
Tools for working with PSFs
"""
import os

from astropy.io import fits
from astropy.modeling import models, fitting
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
    xmean = xdim // 2
    ymean = ydim // 2

    # Translate FWHM to standard deviation
    dev = fwhm / 2.35

    # Create model
    model = Gaussian2D(amplitude=1., x_mean=xmean, y_mean=ymean, x_stddev=dev, y_stddev=dev)

    # Evalulate model
    x_grid, y_grid = np.mgrid[0: xdim, 0: ydim]
    kernel = model.evaluate(x=x_grid, y=y_grid, amplitude=1., x_mean=xmean,
                            y_mean=ymean, x_stddev=dev, y_stddev=dev, theta=0.)
    return kernel


def get_HST_PSF(fwhm_pixels):
    """Get a representative PSF for the HST instrument/detector in
    question. Currently this is just a 2D Gaussian of the appropriate
    FWHM.

    Parameters
    ----------
    fwhm_pixels : float
        FWHM of the Gaussian in units of pixels

    Returns
    -------
    psf : numpy.ndarray
        2D array containing PSF

    """
    # Set the dimensions of the PSF array to be +/- 50 * FWHM
    dim = np.int(np.round(100 * fwhm_pixels))

    # Make sure the array has an odd number of rows and columns
    if dim % 2 == 0:
        dim += 1

    # Create 2D Gaussian
    return gaussian_psf(fwhm_pixels, dim, dim)


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
        Dictionary of instrument information.

    Returns
    -------
    psf_model : numpy.ndarray
        2D array containing the PSF
    """
    # Find the nominal pixel scale, which is what the psf_wings
    # library is saved at
    nominal_pix_scale_x, nominal_pix_scale_y = get_JWST_pixel_scale(meta)

    mosaic_pix_scale1 = meta['pix_scale1']
    if not np.isclose(mosaic_pix_scale1, nominal_pix_scale_x, rtol=0.05):
        raise ValueError(("Reported pixel scale in the mosaic does not match the nominal "
                          "pixel scale for {} {}/{}. If the pixel scale of the mosaic "
                          "is really not {}, please provide your own fits file containing "
                          "a representative PSF at the desired pixel scale.".format(meta['instrument'],
                                                                                    meta['filter'],
                                                                                    meta['pupil'],
                                                                                    nominal_pix_scale_x)))
    # If the pixel scale reported in the mosaic file is the same as the
    # nominal pixel scale for that instrument/detector, then just read in
    # the PSF from the psf_wings file
    instrument = meta['instrument'].lower()
    psf_model = get_psf_wings(meta['instrument'], meta['detector'], meta['filter'], meta['pupil'],
                              'predicted', 0, os.path.join(os.path.expandvars('$MIRAGE_DATA'),
                                                           instrument, 'gridded_psf_library/psf_wings'))

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

        metadata['telescope'] = metadata_entry('TELESCOP', header).upper()
        metadata['instrument'] = metadata_entry('INSTRUME', header).upper()

        if metadata['telescope'] == 'JWST':
            metadata['detector'] = metadata_entry('DETECTOR', header)
            metadata['filter'] = metadata_entry('FILTER', header)
            metadata['pupil'] = metadata_entry('PUPIL', header)
            metadata['pix_scale1'] = np.abs(metadata_entry('CD1_1', header)) * 3600.
            metadata['pix_scale2'] = np.abs(metadata_entry('CD2_2', header)) * 3600.

        if metadata['telescope'] == 'HST':
            metadata['detector'] = metadata_entry('DETECTOR', header)
            metadata['filter'] = metadata_entry('FILTER', header)
            metadata['pa'] = metadata_entry('PA_APER', header)  # PA of reference aperture center
            metadata['pix_scale1'] = np.abs(metadata_entry('CD1_1', header)) * 3600.
            metadata['pix_scale2'] = np.abs(metadata_entry('CD2_2', header)) * 3600.

        if metadata['instrument'] == 'IRAC':
            metadata['channel'] = int(metadata_entry('CHNLNUM', header))
            metadata['pix_scale1'] = np.abs(metadata_entry('PXSCAL1', header))
            metadata['pix_scale2'] = np.abs(metadata_entry('PXSCAL2', header))
            metadata['pa'] = metadata_entry('PA', header)  # deg] Position angle of axis 2 (E of N)

        if metadata['telescope'] not in 'JWST HST SPITZER'.split():
            metadata['pix_scale1'] = np.abs(metadata_entry('CD1_1', header)) * 3600.
            metadata['pix_scale2'] = np.abs(metadata_entry('CD2_2', header)) * 3600.

    return metadata


def measure_fwhm(array):
    """Fit a Gaussian2D model to a PSF and return the FWHM

    Parameters
    ----------
    array : numpy.ndarray
        Array containing PSF

    Returns
    -------
    x_fwhm : float
        FWHM in x direction in units of pixels

    y_fwhm : float
        FWHM in y direction in units of pixels
    """
    yp, xp = array.shape
    y, x, = np.mgrid[:yp, :xp]
    p_init = models.Gaussian2D()
    fit_p = fitting.LevMarLSQFitter()
    fitted_psf = fit_p(p_init, x, y, array)
    return fitted_psf.x_fwhm, fitted_psf.y_fwhm


def metadata_entry(keyword, header):
    """Get the ``header`` keyword value from the ``header`` object. If
    it is not present, return None

    Paramters
    ---------
    keyword : str
        Name of header keyword to examine

    header : astropy.io.fits.header
        Header object from fits file

    Returns
    -------
    value : str, int, float
        Value of the header keyword
    """
    try:
        value = header[keyword]
    except KeyError:
        value = None
    return value


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
