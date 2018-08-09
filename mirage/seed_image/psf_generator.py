#! /suyr/bin/env python

"""
Class to generate PSFs to be used with Mirage
"""
import os

import numpy as np
from astropy.io import fits
from photutils.psf import FittableImageModel


class PSF():
    def __init__(self, x_position, y_position, psf_base,
        interval=0.25, oversampling=1):
        """
        Populate FittableImageModel instance with data from the appropriate PSF file

        Parameters:
        -----------
        x_position : float
            Column number on the detector where the PSF will be located

        y_position : float
            Row number on the detector where the PSF will be located

        psf_base : str
            Base of the PSF filename (e.g. "nircam"

        interval : float
            Sub-pixel centering resolution of the PSF library (e.g. 0.25)

        oversampling : int
            Oversampling factor of the data in the PSF file

        Returns:
        --------
        None
        """
        self.interval = interval
        self.oversampling = oversampling

        # Determine the filename of the appropriate PSF file
        psf_filename = self.find_subpix_psf_filename(x_position, y_position, psf_base)

        # Check for the existence of the file
        if os.path.isfile(psf_filename):
            try:
                # Place PSF in FittableImageModel
                self.model = self.populate_epsfmodel(psf_filename, oversample=oversampling)
            except:
                raise RuntimeError("ERROR: Could not load PSF file {} from library"
                                   .format(psf_filename))
        else:
            raise FileNotFoundError("PSF file {} not found.".format(psf_filename))

    def find_subpix_psf_filename(self, xloc, yloc, basename):
        """Given an x, y location on the
        detector, determine the filename for the most appropriate
        PSF file to use. This function only looks for the sub-pixel
        position, and doesn't know about PSF variation across the
        detector. Therefore only the fractional part of (xloc, yloc)
        is really important.

        Parameters:
        -----------
        xloc : int
            Column number on full detector of source location

        yloc : int
            Row number on full detector of source location

        Returns:
        --------
        psf_filename : str
            Name of fits file containing PSF to use
        """
        # Find sub-pixel offsets in position from the center of the pixel
        xfract, xoff = np.modf(xloc)
        yfract, yoff = np.modf(yloc)

        # Resolution of PSF sub pixel positions
        #interval = self.params['simSignals']['psfpixfrac']
        numperpix = int(1. / self.interval)

        # Now we need to determine the proper PSF
        # file to read in from the library
        # This depends on the sub-pixel offsets above
        a_in = self.interval * int(numperpix * xfract + 0.5) - 0.5
        b_in = self.interval * int(numperpix * yfract + 0.5) - 0.5
        astr = "{0:.{1}f}".format(a_in, 2)
        bstr = "{0:.{1}f}".format(b_in, 2)

        # Generate the psf file name based on the center of the point source
        # in units of fraction of a pixel
        frag = astr + '_' + bstr
        frag = frag.replace('-', 'm')
        frag = frag.replace('.', 'p')

        # Generate fits filename
        psf_filename = basename + '_' + frag + '.fits'
        return psf_filename

    def minimal_psf_evaluation(self):
        """
        Create a PSF by evaluating a FittableImageModel instance. Return
        an array only just big enough to contain the PSF data.

        Parameters:
        -----------
        model : obj
            FittableImageModel instance containing the PSF model

        Returns:
        --------
        eval_psf : ndarray
            2D numpy array containing the evaluated PSF, with normalized signal
        """
        eval_xshape = np.int(np.ceil(self.model.shape[1] / self.model.oversampling))
        eval_yshape = np.int(np.ceil(self.model.shape[0] / self.model.oversampling))
        y, x = np.mgrid[0:eval_yshape, 0:eval_xshape]
        eval_psf = self.model.evaluate(x=x, y=y, flux=1.,
                                            x_0=eval_xshape / 2,
                                            y_0=eval_yshape / 2)

        # For testing
        #h0=fits.PrimaryHDU(eval_psf)
        #hlist = fits.HDUList([h0])
        #hlist.writeto("temp_minimal_eval_psf.fits", overwrite=True)

        return eval_psf

    def populate_epsfmodel(self, infile, oversample=1):
        """Create an instance of EPSFModel and populate the data
        from the given fits file. Also populate information
        about the oversampling rate.

        Parameters:
        -----------
        infile : str
            Name of fits file containing PSF data

        oversample : int
            Factor by which the PSF data are oversampled

        Returns:
        --------
        psf : obj
            FittableImageModel instance
        """
        psf_data = fits.getdata(infile)

        # Normalize
        psf_data /= np.sum(psf_data)

        # Create instance. Assume the PSF is centered
        # in the array
        psf = FittableImageModel(psf_data, oversampling=oversample)
        #psf = EPSFModel(psf_data, oversampling=oversample)
        return psf
