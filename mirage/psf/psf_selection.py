#! /usr/bin/env python

import numpy as np
from astropy.io import fits


class PSFCollection:
    """Class to contain a PSF library across a single detector for a
    single filter. Through interpolation, the PSF at any location on
    the detector can be created."""

    def __init__(self, instrument, filtername, library_path):
        """Upon instantiation of the class, read in the PSF library
        contained in the given fits file. Also pull out relevant
        information such as the oversampling factor and the locations
        on the detector of the fiducial PSFs. This assumes the PSF
        library file is in the format created by ```Mirage's```
        ```CreatePSFLibrary``` class in psf_library.py.

        Parameters:
        -----------
        library_file : str
            Name of a fits file containing the PSF library.

        Returns:
        --------
        None
        """
        library_file = self.get_library_file(instrument, filtername, library_path)

        try:
            with fits.open(library_file) as hdu:
                self.library = hdu[0].data
                self.library_info = hdu[0].header
        except:
            raise OSError("Unable to open {}.".format(psffile))

    def create_interpolator(self, xloc, yloc, data):
        """Create an interpolator function for the detector-dependent
        PSF position. Need a separate interpolator function for each pixel?

        Parameters:
        -----------


        Returns:
        --------
        """
        # Get the locations on the detector for each PSF
        # in the library
        num_psfs = self.library_info['NUM_PSFS']
        x_det = []
        y_det = []
        x_index = []
        y_index = []
        for num in range(len(num_psfs)):
            xval, yval = self.library_info['DET_XY'+num]
            x_det.append(xval)
            y_det.append(yval)
            xi, yi = self.library_info['DET_IJ'+num]
            x_index.append(xi)
            y_index.append(yi)

            
        # Have data of dimension xpsf, ypsf at locations x_det, y_det
        psfydim = header['NAXIS2']
        psfxdim = header['NAXIS1'] 
        x_psf = np.arange(psfxdim)
        y_psf = np.arange(psfydim)

        xs = np.unique(x_det)
        ys = np.unique(y_det)
        interpolator = np.zeros((psfydim, psfxdim))
        for ypix in range(psfydim):
            for xpix in range(psfxdim):
                pixeldata = data[:, :, ypix, xpix]
                # Interpolate this pixel to find value at location
                # finalx, finaly on the detector
                interp2d_fncn = interp2d(xs, ys, pixeldata, kind='linear')
                interpolator[ypix, xpix] = interp2d_fncn

                #and then later to evaluate
                intdata3 = interp2d_fncn(finalx, finaly)
        return interpolator

                
    def get_library_file(self, instrument, filt, library_path):
        """Given an instrument and filter name along with the path of
        the PSF library, find the appropriate library file to load.

        Parameters:
        -----------
        instrument : str
            Name of instrument the PSFs are from

        filt : str
            Name of filter used for PSF library creation

        library_path : str
            Path pointing to the location of the PSF library

        Returns:
        --------
        lib_file : str
            Name of the PSF library file for the instrument and filtername
        """
        lib_file = glob(os.path.join(library_path), instrument + '_' + filt.lower() + '*.fits')

        # If no matching files are found, or more than 1 matching file is found, raise
        # an error.
        if len(lib_file) == 0:
            raise FileNotFoundError("No files matching {}, {} found in {}."
                                    .format(instrument, filt, library_path))
        elif len(lib_file) > 1:
            raise ValueError("Multiple files matching {}, {} found in {}."
                             .format(instrument, filt, library_path))

        return lib_file[0]

    def populate_epsfmodel(self):
        """Create an instance of EPSFModel and populate the data
        from the given fits file. Also populate information
        about the oversampling rate.

        Parameters:
        -----------
        psf_data : numpy.ndarray
            Numpy array containing the 2d image of a PSF

        Returns:
        --------
        psf : obj
            FittableImageModel instance
        """
        # Normalize
        self.at_location /= np.sum(self.at_location)

        # Create instance. Assume the PSF is centered
        # in the array
        oversample = self.library_info['OVERSAMP']
        psf = FittableImageModel(self.at_location, oversampling=oversample)
        #psf = EPSFModel(psf_data, oversampling=oversample)
        return psf




    
    def position_interpolation(self, x, y, method=spline):
        """Interpolate the PSF library to construct the PSF at a
        a given location on the detector.

        Parameters:
        -----------
        x :  int (float?)
            X-coordinate value for the new PSF
        
        y : int (float?)
            Y-coordinate value for the new PSF

        method : str
            Type of interpolation to use

        
        Returns:
        --------
        None
        """
        # Get the locations on the detector for each PSF
        # in the library
        num_psfs = self.library_info['NUM_PSFS']
        x_det = []
        y_det = []
        x_index = []
        y_index = []
        for num in range(len(num_psfs)):
            xval, yval = self.library_info['DET_XY'+num]
            x_det.append(xval)
            y_det.append(yval)
            xi, yi = self.library_info['DET_IJ'+num]
            x_index.append(xi)
            y_index.append(yi)
            
        yd_lib, xd_lib = self.library.shape
        
        # Do we really need to loop over and interpolate each pixel
        # individually???
        for psfx in range(xd_lib):
            for psfy in range(yd_lib):
                data = []
                self.library[detector....]
        self.at_location = interpolated_psf
