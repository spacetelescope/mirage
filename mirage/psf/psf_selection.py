#! /usr/bin/env python
"""Open a PSF library file and interpolate to a field-dependent location.

Authors
-------
 - Bryan Hilbert
 - Lauren Chambers

Use
---\

"""

from ast import literal_eval
from glob import glob
import os

from astropy.io import fits
import numpy as np
#from photutils.utils import ShepardIDWInterpolator as idw
from photutils import FittableImageModel
import pysiaf
from scipy.interpolate import interp2d, RectBivariateSpline
import webbpsf


class PSFCollection:
    """Class to contain a PSF library across a single detector for a
    single filter. Through interpolation, the PSF at any location on
    the detector can be created."""

    def __init__(self, instrument, detector, filtername, library_path, expand_for_segments):
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

        # Open the library file (or files)
        library_list = self.get_library_list(instrument, filtername, library_path, expand_for_segments)

        # Load the library into a list
        self.library_list = []
        for file in library_list:
            try:
                with fits.open(file) as hdu:
                    header = hdu[0].header
                    det_index = self.select_detector(header, detector, file)
                    library_dict = {'header': header,
                                    'data': hdu[0].data[det_index, :, :, :, :]}
                self.library_list.append(library_dict)
            except OSError:
                print("OSError: Unable to open {}.".format(file))
            except IndexError:
                print(("IndexError: File {} has no data in extension 0."
                       .format(file)))

        # Get some basic information about the library
        library = self.library_list[0]
        self.psf_y_dim, self.psf_x_dim = library['data'].shape[-2:]
        self.psf_x_dim /= library['header']['OVERSAMP']
        self.psf_y_dim /= library['header']['OVERSAMP']
        self.num_psfs = library['header']['NUM_PSFS']
        self.x_det = []
        self.y_det = []
        self.x_index = []
        self.y_index = []
        for num in range(self.num_psfs):
            yval, xval = literal_eval(library['header']['DET_YX' + str(num)])
            self.x_det.append(xval)
            self.y_det.append(yval)
            yi, xi = literal_eval(library['header']['DET_JI' + str(num)])
            self.x_index.append(xi)
            self.y_index.append(yi)

        # Get the locations on the detector for each PSF in the library
        if self.num_psfs > 1:
            self.interpolator_list = self.create_interpolator()
        else:
            self.interpolator_list = [None] * len(self.library_list)

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

        interpolator_list = []
        for library in self.library_list:
            # Create interpolator for each pixel
            if interp_type.lower() == 'interp2d':
                interpolator = self.make_interp2d_functions(self.x_det,
                                                            self.y_det,
                                                            library['data'])
            elif interp_type.lower() == 'idw':
                interpolator = None
            else:
                raise ValueError(("interp_type of {} not supported."
                                  .format(interp_type)))
            interpolator_list.append(interpolator)

        return interpolator_list

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

        x_loc_index = np.array(self.x_index)
        y_loc_index = np.array(self.y_index)
        match_index = ranking[0:num_nearest]
        x_nearest = x_loc_index[match_index]
        y_nearest = y_loc_index[match_index]
        result = (x_nearest, y_nearest, distances[match_index])
        return result

    def get_library_list(self, instrument, filt, library_path, expand_for_segments):
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
        lib_file : list
            Name of the PSF library file(s) for the instrument and filtername
        """
        # Find the matching file (or files, if expanding for 18 segments)
        if not expand_for_segments:
            lib_list = glob(os.path.join(library_path, instrument.lower() + '_' +
                                         filt.lower() + '*.fits'))

        elif expand_for_segments:
            lib_list = glob(os.path.join(library_path, instrument.lower() + '_' +
                                         filt.lower() + '*seg*.fits'))
            lib_list = sorted(lib_list)

        # If no matching files are found, or more than 1/18 matching file/s is/are
        # found, raise an error.
        if len(lib_list) == 0:
            raise FileNotFoundError("No files matching {}, {} found in {}."
                                    .format(instrument, filt, library_path))
        elif len(lib_list) > 1 and not expand_for_segments:
            raise ValueError("Multiple files matching {}, {} found in {}."
                             .format(instrument, filt, library_path))
        elif len(lib_list) != 18 and  expand_for_segments:
            raise ValueError("Did not find 18 segment files matching {}, {} found in {}. Instead found {}."
                             .format(instrument, filt, library_path, len(lib_list)))

        return lib_list

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
            xlist = []
            for xpix in range(psf_xdim):
                pixeldata = psfdata[:, :, ypix, xpix]
                # interp2d_fncn = interp2d(xs, ys, pixeldata, kind='linear')
                interp2d_fncn = RectBivariateSpline(xs, ys, pixeldata, kx=1, ky=1)
                xlist.append(interp2d_fncn)
            interp_functions.append(xlist)
        return interp_functions

    def minimal_psf_evaluation(self, model, deltax=0., deltay=0.):
        """
        Create a PSF by evaluating a FittableImageModel instance. Return
        an array only just big enough to contain the PSF data.

         Parameters:
        -----------
        model : obj
            FittableImageModel instance containing the PSF model

        deltax : float
            Offset in the x-direction, in units of nominal pixels, to
            shift the interpolated PSF within the output (interpolated)
            frame.

        deltay : float
            Offset in the y-direction, in units of nominal pixels, to
            shift the interpolated PSF within the output (interpolated)
            frame.

         Returns:
        --------
        eval_psf : ndarray
            2D numpy array containing the evaluated PSF, with normalized signal
        """
        eval_xshape = np.int(np.ceil(model.shape[1] / model.oversampling))
        eval_yshape = np.int(np.ceil(model.shape[0] / model.oversampling))
        y, x = np.mgrid[0:eval_yshape, 0:eval_xshape]
        eval_psf = model.evaluate(x=x, y=y, flux=1.,
                                  x_0=(eval_xshape - 1) / 2 + deltax,
                                  y_0=(eval_yshape - 1) / 2 + deltay)
        return eval_psf

    def populate_epsfmodel(self, psf_data, header):
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
        oversample = header['OVERSAMP']
        psf = FittableImageModel(psf_data, oversampling=oversample)
        # psf = EPSFModel(psf_data, oversampling=oversample)
        return psf

    def position_interpolation(self, x, y, method="spline", idw_number_nearest=4,
                               idw_alpha=-2, segment_number=None):
        """Interpolate the PSF library to construct the PSF at a
        a given location on the detector. Note that the coordinate system used
        in this case has (0.0, 0.0) centered in the lower left pixel. (0.5, 0.5)
        corresponds to the upper right corner of that pixel, and (-0.5, -0.5) the
        lower left corner. This is important more for when the FittableImageModel
        that is returned here is evaluated.

        Parameters
        ----------
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

        Returns
        -------
        out : FittableImageModel
            Instance of FittableImageModel containing the interpolated PSF
        """

        # If needed, restrict the library to the data for the desired segment.
        if segment_number is None:
            library = self.library_list[0]['data']
            header = self.library_list[0]['header']
            interpolator = self.interpolator_list[0]
        else:
            library = self.library_list[segment_number - 1]['data']
            header = self.library_list[segment_number - 1]['header']
            assert int(header['SEGID']) == int(segment_number), \
                "Uh-oh. The segment ID of the library does not match the requested " \
                "segment. The library_list was not assembled correctly."
            interpolator = self.interpolator_list[segment_number - 1] # Really, this is always none....

        if self.num_psfs == 1:
            interp_psf = library[0, 0, :, :]
        else:
            if method == "spline":
                interp_psf = self.evaluate_interp2d(interpolator, x, y)
            elif method == "idw":
                nearest_x, nearest_y, nearest_dist = self.find_nearest_neighbors(x, y, idw_number_nearest)
                if len(nearest_dist) > 1:
                    psfs_to_avg = library[nearest_x, nearest_y, :, :]
                    interp_psf = self.weighted_avg(psfs_to_avg, nearest_dist,
                                                   idw_alpha)
                elif len(nearest_dist) == 1:
                    interp_psf = library[nearest_x, nearest_y, :, :]
            else:
                raise ValueError(("{} interpolation method not supported."
                                .format(method)))

        # Return resulting PSF in instance of FittableImageModel (or EPSFModel)
        return self.populate_epsfmodel(interp_psf, header)

    def select_detector(self, header, det_name, input_file):
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
            Index number indicating which slice of the library file corresponds
            to the given detector
        """
        det_list = []
        for i in range(header['NAXIS5']):
            det_list.append(header['DETNAME' + str(i)])

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

    def get_segment_offset(self, segment_number, detector):
        """Convert vectors coordinates in the local segment control
        coordinate system to NIRCam detector X and Y coordinates,
        at least proportionally, in order to calculate the location of

        Parameters
        ----------
        segment : int
            Segment ID, i.e 3

        Returns
        -------
        x_displacement
            The shift of the segment PSF in NIRCam SW x pixels
        y_displacement
            The shift of the segment PSF in NIRCam SW y pixels
        """
        # Verify that the segment number in the header matches the index
        assert int(self.library_list[segment_number - 1]['header']['SEGID']) == int(segment_number),\
            "Uh-oh. The segment ID of the library does not match the requested " \
            "segment. The library_list was not assembled correctly."

        xtilt = self.library_list[segment_number - 1]['header']['XTILT']
        ytilt = self.library_list[segment_number - 1]['header']['YTILT']
        segment = self.library_list[segment_number - 1]['header']['SEGNAME'][:2]

        # These conversion factors were empirically calculated by measuring the
        # relation between tilt and the pixel displacement
        tilt_to_pixel_slope = 13.4
        tilt_to_pixel_intercept = 0

        # segment = webbpsf.constants.SEGNAMES_WSS_ORDER[segment_number - 1][:2]

        control_xaxis_rotations = {
            'A1': 180, 'A2': 120, 'A3': 60, 'A4': 0, 'A5': -60,
            'A6': -120, 'B1': 0, 'C1': 60, 'B2': -60, 'C2': 0,
            'B3': -120, 'C3': -60, 'B4': -180, 'C4': -120,
            'B5': -240, 'C5': -180, 'B6': -300, 'C6': -240
        }

        x_rot = control_xaxis_rotations[segment]  # degrees
        x_rot_rad = x_rot * np.pi / 180  # radians

        # Note that y is defined as the x component and x is defined as the y component.
        # This is because "xtilt" moves the PSF in the y direction, and vice versa.
        tilt_onto_y = (xtilt * np.cos(x_rot_rad)) - (ytilt * np.sin(x_rot_rad))
        tilt_onto_x = (xtilt * np.sin(x_rot_rad)) + (ytilt * np.cos(x_rot_rad))

        # TODO: IS THE SLOPE DIFFERENT FOR LW DETECTORS????
        x_displacement = -(tilt_onto_x * tilt_to_pixel_slope)  # pixels
        y_displacement = -(tilt_onto_y * tilt_to_pixel_slope)  # pixels

        # Get the appropriate pixel scale from pysiaf
        siaf = pysiaf.Siaf('nircam')
        aperture = siaf['NRC{}_FULL'.format(detector[-2:])]
        nircam_x_pixel_scale = aperture.XSciScale  # arcsec/pixel
        nircam_y_pixel_scale = aperture.YSciScale  # arcsec/pixel

        # Convert the pixel displacement into angle
        x_arcsec = x_displacement * nircam_x_pixel_scale  # arcsec
        y_arcsec = y_displacement * nircam_y_pixel_scale  # arcsec

        return x_arcsec, y_arcsec

        # return x_displacement, y_displacement