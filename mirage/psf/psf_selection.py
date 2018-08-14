#! /usr/bin/env python

from ast import literal_eval
from glob import glob
import os

from astropy.io import fits
import numpy as np
#from photutils.utils import ShepardIDWInterpolator as idw
from photutils import FittableImageModel
from scipy.interpolate import interp2d


class PSFCollection:
    """Class to contain a PSF library across a single detector for a
    single filter. Through interpolation, the PSF at any location on
    the detector can be created."""

    def __init__(self, instrument, detector, filtername, library_path):
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
                self.library_info = hdu[0].header
                det_index = self.select_detector(detector, library_file)
                self.library = hdu[0].data[det_index, :, :, :, :]
        except OSError:
            print("OSError: Unable to open {}.".format(library_file))
        except IndexError:
            print(("IndexError: File {} has no data in extension 0."
                   .format(library_file)))

        # Get some basic information about the library
        self.num_psfs = self.library_info['NUM_PSFS']
        self.x_det = []
        self.y_det = []
        self.x_index = []
        self.y_index = []
        for num in range(self.num_psfs):
            xval, yval = literal_eval(self.library_info['DET_XY' + str(num)])
            self.x_det.append(xval)
            self.y_det.append(yval)
            xi, yi = literal_eval(self.library_info['DET_IJ' + str(num)])
            self.x_index.append(xi)
            self.y_index.append(yi)
        self.interpolator = self.create_interpolator()
        print("check_over_x_y_order_in_all_functions")

    def create_interpolator(self, interp_type='interp2d'):
        """Create an interpolator function for the detector-dependent
        PSF position. Need a separate interpolator function for each pixel?

        Parameters:
        -----------
        interp_type: str
            Can be 'interp2d' or 'idw'

        pass self.x_det as xloc, or just use the class variable??

        Returns:
        --------
        """
        # Get the locations on the detector for each PSF
        # in the library

        # Have data of dimension xpsf, ypsf at locations x_det, y_det
        # psfydim = self.library_info['NAXIS2']
        # psfxdim = self.library_info['NAXIS1']
        # x_psf = np.arange(psfxdim)
        # y_psf = np.arange(psfydim)

        # Create interpolator for each pixel
        if interp_type.lower() == 'interp2d':
            interpolator = self.make_interp2d_functions(self.x_det,
                                                        self.y_det,
                                                        self.library)
        elif interp_type.lower() == 'idw':
            interpolator = None
        else:
            raise ValueError(("interp_type of {} not supported."
                              .format(interp_type)))
        return interpolator

    def evaluate_interp2d(self, splines, xout, yout):
        """Evaluate the splines produced in make_interp2d_functions for each
        pixel in the PSF given a location on the detector

        Parameters:
        -----------
        splines : list
            List of lists of splines. Output from make_interp2d_functions

        xout : int
            x coordinate on detector at which to evaluate

        yout : int
            y coordinate on detector at which to evaluate

        Returns:
        --------
        psf : numpy.ndarray
            2D array containing the interpolated PSF
        """
        xlen = len(splines[0])
        ylen = len(splines)
        psf = np.zeros((ylen, xlen))
        for y, rowlist in enumerate(splines):
            for x, spline in enumerate(rowlist):
                pix = spline(xout, yout)
                psf[y, x] = pix
        return psf

    def find_nearest_neighbors(self, x_out, y_out, num_nearest):
        """Determine the location of the num_nearest nearest reference PSFs
        to a given location on the detector.

        Parameters:
        -----------
        x_out : int
            x-coordinate on the detector to find distances from

        y_out : int
            y-coordinate on the detector to find distances from

        num_nearest : int
            Number of nearest reference PSFs to find

        Returns:
        --------
        match_indexes : list
            List containing 4 tuples; one for each of the 4 nearest
            reference PSFs. Each tuple contains (index of x_psfs,y_psfs
            list, radial distance to x_out,yout)
        """
        x_psfs = np.array(self.x_det)
        y_psfs = np.array(self.y_det)
        distances = np.sqrt(np.abs(x_psfs - x_out)**2 +
                            np.abs(y_psfs - y_out)**2)
        ranking = np.argsort(distances)

        # If the requested PSF position falls exactly on the location
        # of one of the reference PSFs, then use just that single PSF
        # with a weight of 1.0
        if min(distances) == 0:
            num_nearest = 1

        match_index = ranking[0:num_nearest]
        x_nearest = x_psfs[match_index]
        y_nearest = y_psfs[match_index]
        result = (x_nearest, y_nearest, distances[match_index])
        return result

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
        lib_file = glob(os.path.join(library_path, instrument.lower() + '_' +
                                     filt.lower() + '*.fits'))

        # If no matching files are found, or more than 1 matching file is
        # found, raise an error.
        if len(lib_file) == 0:
            raise FileNotFoundError("No files matching {}, {} found in {}."
                                    .format(instrument, filt, library_path))
        elif len(lib_file) > 1:
            raise ValueError("Multiple files matching {}, {} found in {}."
                             .format(instrument, filt, library_path))

        return lib_file[0]

    def make_interp2d_functions(self, xpos, ypos, psfdata):
        """Create a list of lists containing splines for each pixel in the PSF
        image. Note that this assumes a regular grid of locations for the
        reference PSFs. (e.g. Assume np.unique creates an input x list of
        [0, 680, 1365, 2047], and a y list with the same values, then it is
        expected that there are 16 PSFs in the libaray at locations
        corresponding to all combinations of the values in the x and y lists.

        Parameters:
        -----------
        xpos : list
            X coordinate on the detector corresponding to each PSF

        ypos : list
            Y coordinate on the detector corresponding to each PSF

        psfdata : numpy.ndarray
            PSF library. 4D array with dimensions (x loc on detector, y loc on
            detector, x PSF size, y PSF size)

        Returns:
        --------
        interp_functions : list
            List of lists containing the spline for each pixel in the PSF
        """
        xs = np.unique(xpos)
        ys = np.unique(ypos)
        interp_functions = []
        loc_ydim, loc_xdim, psf_ydim, psf_xdim = psfdata.shape

        # Need to loop over the pixels and produce an interpolating function
        # for each.
        for ypix in range(psf_ydim):
            for xpix in range(psf_xdim):
                xlist = []
                pixeldata = self.library[:, :, ypix, xpix]
                interp2d_fncn = interp2d(xs, ys, pixeldata, kind='linear')
                xlist.append(interp2d_fncn)
            interp_functions.append(xlist)
        return interp_functions

    def populate_epsfmodel(self, psf_data):
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
        psf_data /= np.sum(psf_data)

        # Create instance. Assume the PSF is centered
        # in the array
        oversample = self.library_info['OVERSAMP']

        print(type(psf_data))
        print(psf_data.shape)
        print(oversample)
        stop

        psf = FittableImageModel(psf_data, oversampling=oversample)
        # psf = EPSFModel(psf_data, oversampling=oversample)
        return psf

    def position_interpolation(self, x, y, method="spline", idw_number_nearest=4, idw_alpha=-2):
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

        idw_number_nearest : int
            For weighted averaging (method='idw'), the number of nearest
            reference PSFs to the location of the interpolated PSF to use
            when calculating the weighted average. Default = 4.

        idw_alpha : float
            Exponent used for turning distances into weights. Default is -2.

        Returns:
        --------
        None
        """
        # Get the locations on the detector for each PSF
        # in the library
        if method == "spline":
            interp_psf = self.evaluate_interp2d(self.interpolator, x, y)
        elif method == "idw":
            nearest_x, nearest_y, nearest_dist = self.find_nearest_neighbors(x, y, idw_number_nearest)
            if len(nearest_dist) > 1:
                psfs_to_avg = self.library[nearest_x, nearest_y, :, :]
                interp_psf = self.weighted_avg(psfs_to_avg, nearest_dist,
                                               idw_alpha)
            elif len(nearest_dist) == 1:
                interp_psf = self.library[nearest_x, nearest_y, :, :]
        else:
            raise ValueError(("{} interpolation method not supported."
                              .format(method)))

        # Place resulting PSF in instance of FittableImageModel (or EPSFModel)
        self.psf = self.populate_epsfmodel(interp_psf)

    def select_detector(self, det_name, input_file):
        """Given a PSF library, select only the PSFs associated with a
        given detector.

        Parameters:
        -----------
        det_name : str
            Name of the detector whose PSFs are to be returned.

        input_file : str
            Name of the file containing the PSF library. Only used for
            printing error messages.

        Returns:
        --------
        match : int
            Index number indicating which slice of self.library corresponds
            to the given detector
        """
        det_list = []
        for i in range(self.library_info['NAXIS5']):
            det_list.append(self.library_info['DETNAME' + str(i)])

        if det_name in det_list:
            match = np.where(np.array(det_list) == det_name)[0][0]
            return match
        else:
            raise ValueError(("No PSFs for detector {} in {}."
                              .format(det_name, input_file)))

    def weighted_avg(self, psfs, distances, alpha):
        """Calculate the weighted average of the input PSFs given their distances
        to the location of interest on the detector

        Parameters:
        -----------
        psfs : numpy.ndarray
            Array containing the PSFs to be used to calculate the weighted
            average. 3D array

        distances : list
            List of radial distances from the location of interest to the PSFs
            in psfs

        alpha : float
            Exponent to use in determining weights. weights = distance**alpha
            Should be between -2 and 0.

        Returns:
        --------
        avg_psf : numpy.ndarray
            Array containing the weighted average PSF
        """
        alpha = np.float(alpha)
        weights = distances**alpha
        output = np.average(psfs, axis=0, weights=weights)
        return output
