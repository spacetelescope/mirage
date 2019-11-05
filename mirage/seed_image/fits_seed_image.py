#! /usr/bin/env python

"""
Create a mirage-format seed image from
an input fits file. This fits file must contain a
valid WCS, and be distortion-free.

This code will extract a sub-image from the input
fits file. The subimage is centered at the input
RA and Dec (crop_center_ra, crop_center_dec) and
its size is determined by the aperture input.

The extracted sub-image is then blotted (resampled)
onto the requested NIRCam detector's FOV, including
the distortion model.

The resulting image is then saved in such a format
that it can be used as input in the obs_generator.py,
image_simulator.py, or disperser code in the
mirage package. Essentially this script
is designed to take the place of catalog_seed_image.py
in the case where the user has a distortion-free
mosaic image that they wish to break up into
NIRCam images.

Another way to use this code would be for an input
fits image of some extended object that the user
wishes to distort according to the NIRCam distortion
model. This distorted file can then be placed in
an extended object catalog which will be used by
the catalog_seed_image.py step. In this method,
multiple extended objets could be distorted and
added to the same frame.

Inputs:
mosaic_file: The name of the fits file containing the
            distortion-free image to be blotted. The code
            currently assumes that the mosaic_file is
            oriented north-up.

From here, there are two ways to call fits_seed_image.py

1. parameter file: a yaml input file matching the format
   (optional)      required by other steps of the nircam
                   data simulator

2. Manual inputs:
   aperture -
              aperture name matching an entry in the
              subarray definition file (e.g. NRCB5_FULL)

   crop_center_ra, crop_center_dec -
              The RA and Dec at the center of the sub-image
              to crop from mosaic_file. This, in combination
              with the array size (determined from aperture)
              define the sub-image to use.

   blot_center_ra, blot_center_dec -
              The RA and Dec at the center of the blotted
              sub-image. If these are equal to crop_center_ra
              and dec, then the center of the cropped image
              will remain at the center of the blotted image.

   blot_pav3 -
              Position angle of the blotted image.

   flux_cal_file -
              Ascii file listing the zeropoints for all
              NIRCam filters. By setting this value equal
              to 'config', the code will use the flux
              calibration file that is packaged with the
              NIRCam data simulator. This information is
              used only to populate the PHOTFLAM, PHOTFNU,
              and PHOTPLAM header keywords.

   filter -
              The NIRCam filter to use. This is used to
              select the appropriate row from the flux
              calibration file.

   pupil -
              The NIRCam pupil value to use. This is used
              to select the appropriate row from the flux
              calibration file in the case where a filter
              that is located in the pupil wheel is used.
              If the pupil value does not start with 'F'
              (e.g. F405N), then it is ignored.

   grism_source_image -
              True/False. If you intend to send the seed
              image output to the disperser software, then
              set this value to True. The output image will
              then be made larger than the nominal aperture
              size by a factor of sqrt(2), in order for the
              disperser to keep track of sources just outside
              the nominal FOV, whose dispersed spectra may
              fall onto the detector.

   outfile -
              Name of the file to contain the output seed
              image.

   outdir -
              Directory in which to place outfile

Outputs:
self.seedimage, self.seed_segmap, self.seedinfo
contain the seed (countrate) image, the associated
segmentation map, and header information required
by subsequent steps of the nircam data simulator.
The segmentation map is only used in the case where
the seed image is dispersed by a call to the
disperser software.

Example calls:
See the Simulated_data_from_mosaic_image.ipynb notebook
in the examples directory
"""
import argparse
import copy
import os
import pkg_resources
import sys

from astropy.io import fits, ascii
import datetime
from math import radians
import numpy as np
from photutils import detect_sources
from photutils import TopHatWindow, TukeyWindow, CosineBellWindow, SplitCosineBellWindow, HanningWindow
from photutils.centroids import centroid_2dg
from photutils.psf import resize_psf
from photutils.psf.matching import create_matching_kernel
import pysiaf
from scipy.ndimage import shift
from scipy.signal import fftconvolve
import yaml

from . import crop_mosaic, blot_image
from mirage.psf.psf_selection import get_psf_wings
from mirage.psf import tools
from mirage.seed_image.save_seed import save
from mirage.reference_files import crds_tools
from mirage.utils.constants import EXPTYPES
from mirage.utils.siaf_interface import get_siaf_information


config_files = {'nircam': {'flux_cal': 'NIRCam_zeropoints.list'},
                'niriss': {'flux_cal': 'niriss_zeropoints.list'}
                }

KNOWN_PSF_TELESCOPES = {"JWST", "HST", "SPITZER"}


class ImgSeed:
    def __init__(self, paramfile=None, mosaic_file=None, data_extension_number=0, wcs_extension_number=0,
                 cropped_file='cropped_image.fits', outdir=None, blotted_file=None, psf_file=None,
                 mosaic_fwhm=None, mosaic_fwhm_units='arcsec', gaussian_psf=True, save_intermediates=False):
        """Create a seed image from a distortionless input mosaic. The user can supply
        a PSF associated with the mosaic in the ``psf_file`` parameter. If ``psf_file``
        is None and ``gaussian_psf`` is True, then astropy's Gaussian2D model will be
        used to construct a PSF for the mosaic. If ``psf_file`1` is None and ``gaussian_psf``
        is False, and the mosaic comes from a "known" telescope, then a custom function
        will be used to construct the PSF. Currently the "known" telescopes include
        JWST, HST, and Spitzer IRAC.

        Note that this code cannot be used in cases where the mosaic's PSF FWHM is
        larger than the JWST PSF, as translating the mosaic to the JWST PSF would then
        be sharpening the image.

        Parameters
        ----------
            paramfile : str
                Name of Mirage yaml input file used to supply instrument info

            mosaic_file : str
                Name of the fits file containing the mosaic file to use

            data_extension_number : int
                Extension number within the mosaic FITS file where the data
                array is located

            wcs_extension_number : int
                Extension number within the mosaic FITS file where the
                world coordinate system information is located

            cropped_file : str
                Name of the file to save the cropped image into. If None,
                the cropped image is not saved

            outdir : str
                Name of the output directory to save the output files into.
                If None, the output directory listed in ``paramfile`` will
                be used.

            blotted_file : str
                Name of FITS file to save resampled image into

            psf_file : str
                Name of FITS file containing a PSF corresponding to the mosaic
                data. This will be used to create a matching kernel to allow
                convolution of the mosaic data to match the JWST PSF. If None,
                a PSF will be created using an astropy Gaussian2D model or,
                in the case of a mosaic from JWST, a webbpsf-generated PSF
                will be read in.

            mosaic_fwhm : float
                FWHM of the PSF in the mosaic. Can be in pixels or arcseconds,
                as specified by the ``mosaic_fwhm_units`` keyword. If None,
                an astropy Gaussian2D model will be fit to the PSF to estimate
                the FWHM

            mosaic_fwhm_units : str
                Units of ``mosaic_fwhm``. Can be 'pixels' or 'arcsec'

            gaussian_psf : bool
                If ``psf_file`` is None and this is True, the mosaic PSF will
                be created using an astropy Gaussian2D model.

            save_intermediates : bool
                If True, intermediate outputs, including the PSFs and the
                matching PSF kernel, as well as the cropped image after
                convolution with the matching PSF kernel, will be saved.
        """
        allowed_mosaic_fwhm_units = ['pixels', 'arcsec']

        self.mosaic_file = mosaic_file
        self.data_extension_number = data_extension_number
        self.wcs_extension_number = wcs_extension_number
        self.crop_center_ra = 0.
        self.crop_center_dec = 0.
        self.channel = None
        self.detector = None
        self.blot_center_ra = 0.
        self.blot_center_dec = 0.
        self.blot_pav3 = 0.
        self.aperture = ''
        self.flux_cal_file = 'config'
        self.distortion_file = 'crds'
        self.filter = ''
        self.pupil = ''
        self.grism_source_image = False
        self.cropped_file = cropped_file
        self.blotted_file = blotted_file
        self.outdir = outdir
        self.grism_direct_factor = np.sqrt(2.)
        self.psf_file = psf_file
        self.gaussian_psf = gaussian_psf
        self.mosaic_fwhm = mosaic_fwhm
        if mosaic_fwhm_units not in allowed_mosaic_fwhm_units:
            raise ValueError(("ERROR: mosaic_fwhm_units must be one of: {}"
                              .format(allowed_mosaic_fwhm_units)))
        self.mosaic_fwhm_units = mosaic_fwhm_units
        self.save_intermediates = save_intermediates

        # Locate the module files, so that we know where to look
        # for config subdirectory
        self.modpath = pkg_resources.resource_filename('mirage', '')

        # self.coords contains the factor by which the
        # nominal output array size needs to be increased
        # (used for WFSS mode), as well as the coordinate
        # offset between the nominal output array coordinates,
        # and those of the expanded array. These are needed
        # mostly for WFSS observations, where the nominal output
        # array will not sit centered in the expanded output image.
        self.coords = {'x': 1., 'xoffset': 0., 'y': 1., 'yoffset': 0.}

        # If a paramfile is provided, read in and
        # set params
        self.paramfile = paramfile
        if paramfile is not None:
            self.read_param_file(paramfile)

    def aperture_info(self):
        """
        Get basic information on the reuqested aperture from SIAF
        """
        local_roll, attitude_matrix, self.framesize, subarray_bounds = get_siaf_information(self.siaf,
                                                                                            self.aperture,
                                                                                            self.crop_center_ra,
                                                                                            self.crop_center_dec,
                                                                                            self.blot_pav3)
        xstart, ystart, xend, yend = subarray_bounds

        self.subarr_bounds = {'xstart': xstart, 'ystart': ystart, 'xend': xend, 'yend': yend}

        self.nominal_dims = np.array([yend - ystart + 1, xend - xstart + 1])

    def crop_and_blot(self):
        """MAIN FUNCTION

        Extract an area of the input mosaic fits file that is slightly larger
        than the desired aperture, and use the JWST pipeline version of blot (called
        resample) to resample this to the correct pixel scale and introduce the
        proper distortion for the requested detector
        """
        self.detector_channel_value()
        self.distortion_file_value()
        module = self.module_value()
        self.siaf = pysiaf.Siaf(self.instrument)

        # Get information on the aperture
        self.aperture_info()
        xstart = self.subarr_bounds['xstart']
        ystart = self.subarr_bounds['ystart']
        xdim = self.subarr_bounds['xend'] - xstart + 1
        ydim = self.subarr_bounds['yend'] - ystart + 1

        # Flux calibration information
        self.flux_cal_file = self.expand_config_paths(self.flux_cal_file, 'flux_cal')
        if self.pupil[0].upper() == 'F':
            usefilt = self.pupil.upper()
        else:
            usefilt = self.filter.upper()
        self.fluxcal_info(usefilt, module)

        # Get information on the instrument used to create the mosaic
        self.mosaic_metadata = tools.get_psf_metadata(self.mosaic_file)

        # Check that the FWHM of the mosaic PSF is smaller than the PSF
        # of the JWST PSF. If the opposite is true, we can't create a
        # matching kernel
        if self.mosaic_fwhm_units == 'arcsec':
            mosaic_fwhm_arcsec = copy.deepcopy(self.mosaic_fwhm)
            self.mosaic_fwhm /= self.mosaic_metadata['pix_scale1']

        # Identify the correct PSF for the JWST instrument/detector
        # Let's just read in the psf_wings file, which contains the PSF
        # at the center of the detector. It doesn't make much sense to
        # worry about spatially-dependent PSFs in this case
        self.jwst_psf = get_psf_wings(self.instrument, self.detector, self.filter, self.pupil,
                                      self.params['simSignals']['psfwfe'],
                                      self.params['simSignals']['psfwfegroup'],
                                      os.path.join(self.params['simSignals']['psfpath'], 'psf_wings'))
        self.outscale1, self.outscale2 = self.find_jwst_pixel_scale()

        # Find the FWHM of the JWST PSF
        j_yd, j_xd = self.jwst_psf.shape
        mid_y = j_yd // 2
        mid_x = j_xd // 2
        box = self.jwst_psf[mid_y-10:mid_y+11, mid_x-10:mid_x+11]
        jwst_x_fwhm, jwst_y_fwhm = tools.measure_fwhm(box)
        jwst_y_fwhm_arcsec = jwst_y_fwhm * self.outscale2
        print('JWST FWHM in pix: ', jwst_y_fwhm)
        print('JWST FWHM in arcsec: ', jwst_y_fwhm_arcsec)

        # If the FWHM of the mosaic image is larger than that of the JWST
        # PSF, then we cannot continue because we cannot create a matching
        # PSF kernel to use for convolution
        if mosaic_fwhm_arcsec > jwst_y_fwhm_arcsec:
            raise ValueError(("ERROR: FWHM of the mosaic image is larger than that of "
                              "the JWST PSF. Unable to create a matching PSF kernel "
                              "using photutils. \nMosaic FWHM: {}\nJWST FWHM: {}"
                              .format(self.mosaic_fwhm, jwst_y_fwhm_arcsec)))

        # The JWST PSF as read in is large (399 x 399 pixels). Crop to
        # something reasonable so that the convolution doesn't take
        # forever
        half_width = 50
        self.jwst_psf = self.jwst_psf[mid_y-half_width: mid_y+half_width+1, mid_x-half_width: mid_x+half_width+1]

        # Collect the metadata relating to the mosaic and (optionally) PSF
        # and read in/create the PSF
        self.prepare_psf()

        # Crop
        print("Cropping subarray from mosaic")
        pixel_scale = self.siaf[self.aperture].XSciScale
        crop = crop_mosaic.Extraction(mosaicfile=self.mosaic_file,
                                      data_extension_number=self.data_extension_number,
                                      wcs_extension_number=self.wcs_extension_number,
                                      center_ra=self.crop_center_ra,
                                      center_dec=self.crop_center_dec,
                                      dimensions=(xdim, ydim),
                                      outfile=self.cropped_file,
                                      jwst_pixel_scale=pixel_scale)
        crop.extract()

        # Convolve the cropped image with the appropriate PSF
        print("Convolving with PSF kernel")
        crop.cropped = self.psf_convolution(crop.cropped)

        # Blot
        print("Resampling subarray onto JWST pixel grid")
        blot = blot_image.Blot(instrument=self.instrument, aperture=self.aperture,
                               ra=[self.blot_center_ra], dec=[self.blot_center_dec],
                               pav3=[self.blot_pav3], blotfile=crop.cropped,
                               distortion_file=self.distortion_file)
        blot.blot()

        # Blot the PSF associated with the mosaic
        #blot_psf = blot_image.Blot(instrument=self.instrument, aperture=self.aperture,
        #                           ra=[self.blot_center_ra], dec=[self.blot_center_dec],
        #                           pav3=[self.blot_pav3], blotfile=mosaic_psf_model,
        #                           distortion_file=self.distortion_file)
        #blot_psf.blot()

        for dm in blot.blotted_datamodels:
            # Create segmentation map
            # (used only when dispsersing image for WFSS)
            self.seed_image = dm.data
            self.seed_segmap = self.make_segmap(dm)

            # Save
            if self.blotted_file is not None:
                self.blotted_file = os.path.join(self.outdir, self.blotted_file)
            self.seed_file, self.seedinfo = save(dm.data, self.paramfile, self.params,
                                                 self.photflam, self.photfnu, self.pivot,
                                                 self.framesize, self.nominal_dims,
                                                 self.coords, self.grism_direct_factor,
                                                 segmentation_map=self.seed_segmap,
                                                 filename=self.blotted_file, base_unit='e-')
            print('Blotted image saved to: {}'.format(self.seed_file))

    def detector_channel_value(self):
        """Get the detector and optional channel value based on the aperture
        name
        """
        if self.detector is None:
            self.detector = self.aperture[0:5]
        self.detector = self.detector.upper()
        self.instrument = self.instrument.upper()

        # Make sure the FGS detector is GUIDER[12]
        if 'FGS' in self.detector:
            self.detector = 'GUIDER{}'.format(self.detctor[-1])

        if 'NRC' in self.detector:
            if self.channel is None:
                if '5' in self.detector:
                    self.channel = 'LONG'
                else:
                    self.channel = 'SHORT'

    def distortion_file_value(self):
        """Make sure the distortion file name is a valid file. If it has a
        value of 'crds', then query CRDS and retrieve the appropriate file
        """
        if self.distortion_file == 'crds':
            if self.paramfile is None:
                self.reference_dictionary = {}
                self.reference_dictionary['INSTRUME'] = self.instrument
                self.reference_dictionary['DETECTOR'] = self.detector
                if '5' in self.detector:
                    self.reference_dictionary['DETECTOR'] = self.detector.replace('5', 'LONG')

                if self.instrument.upper() == 'NIRCAM':
                    self.reference_dictionary['CHANNEL'] = self.channel

                if 'FGS' in self.reference_dictionary['DETECTOR']:
                    self.reference_dictionary['DETECTOR'] = 'GUIDER{}'.format(self.reference_dictionary['DETECTOR'][-1])

                # Use the current date and time in order to get the most recent
                # reference file
                self.reference_dictionary['DATE-OBS'] = datetime.date.today().isoformat()
                current_date = datetime.datetime.now()
                self.reference_dictionary['TIME-OBS'] = current_date.time().isoformat()

                # For the purposes of choosing reference files, the exposure type should
                # always be set to imaging, since it is used to locate sources in the
                # seed image, prior to any dispersion.
                self.reference_dictionary['EXP_TYPE'] = EXPTYPES[self.instrument.lower()]["imaging"]

                # This assumes that filter and pupil names match up with reality,
                # as opposed to the more user-friendly scheme of allowing any
                # filter to be in the filter field.
                self.reference_dictionary['FILTER'] = self.filter
                self.reference_dictionary['PUPIL'] = self.pupil

            mapping = crds_tools.get_reffiles(self.reference_dictionary, ['distortion'], download=True)
            self.distortion_file = mapping['distortion']

    def expand_config_paths(self, filename, filetype):
        """If a reference file has the name "config" in the input yaml
        file, then convert that to the appropriate file name. This is
        done using a global dictionary.

        Parameters
        ----------
        filename : str
            The filename as listed in the yaml file

        filetype : str
            The type of reference file in question. Currently
            the dictionary only has an entry for 'flux_cal'

        Returns
        -------
        filename : str
            The updated filename for the reference file type in question
        """
        if filename.lower() == 'config':
            sfile = config_files[self.instrument.lower()][filetype]
            return os.path.join(self.modpath, 'config', sfile)
        else:
            return filename

    def find_jwst_pixel_scale(self):
        """Find the pixel scale for the JWST instrument and aperture to
        be used for the seed image
        """
        meta_output = {}
        meta_output['instrument'] = self.instrument
        return tools.get_JWST_pixel_scale(meta_output, aperture=self.aperture)

    def fluxcal_info(self, filter_name, module):
        """Get the flux calibration-related information from the flux
        cal reference file.

        Parameters
        ----------
        filter_name : str
            Name of the filter to look up the flux cal info for

        module : str
            Module name to use for the look-up 'A' or 'B' for
            NIRCam, and 'N' for NIRISS
        """
        # Read in list of zeropoints/photflam/photfnu
        self.zps = ascii.read(self.flux_cal_file)
        mtch = ((self.zps['Filter'] == filter_name) \
                & (self.zps['Module'] == module))
        self.photflam = self.zps['PHOTFLAM'][mtch][0]
        self.photfnu = self.zps['PHOTFNU'][mtch][0]
        self.pivot = self.zps['Pivot_wave'][mtch][0]

    def grism_coord_adjust(self):
        """Calculate the factors by which to expand the
        output array size, as well as the coordinate
        offsets between the nominal output array and
        the input lists if the observation being
        modeled is to be dispersed with the grism
        """
        dtor = radians(1.)

        # Normal imaging with grism image requested
        self.coords['x'] = self.grism_direct_factor
        self.coords['y'] = self.grism_direct_factor
        self.coords['xoffset'] = np.int((self.grism_direct_factor - 1.)
                                        * (self.subarr_bounds['xend'] -
                                           self.subarr_bounds['xstart'] + 1) / 2.)
        self.coords['yoffset'] = np.int((self.grism_direct_factor - 1.)
                                        * (self.subarr_bounds['yend'] -
                                           self.subarr_bounds['ystart']+1) / 2.)

    def make_segmap(self, model):
        """Create a segmentation map of the input image

        Parameters
        ----------
        model : jwst.datamodels.ImageModel

        Returns
        -------
        map_image : numpy.ndarray
            2D segmentation map
        """
        noise = np.median(model.data)
        seg_map = detect_sources(model.data, noise*5., 8)
        if seg_map is None:
            print('No segmentation map created. Returning empty segmap')
            map_image = np.zeros_like(model.data)
        else:
            map_image = seg_map.data
        return map_image

    @staticmethod
    def matching_kernel(psf1, psf2, window_type='TukeyWindow', alpha=None, beta=None):
        """Use photutils to create a matching kernel given two PSFs and
        the window type and parameters

        Parameters
        ----------
        psf1 : numpy.ndarray
            2D array containing the first PSF

        psf2 : numpy.ndarray
            2D array containing the second PSF

        window_type : str
            Name of the window function to use when filtering the matching kernel

        alpha : float
            Optional input for some of the window functions

        beta : float
            Optional input for some of the window functions

        Returns
        -------
        matched_kernel : numpy.ndarray
            2D array containing the matching PSF kernel
        """
        # Create the filtering window
        orig_window_type = copy.deepcopy(window_type)
        window_type = window_type.lower()
        if window_type == 'tophatwindow':
            window = TopHatWindow(beta=beta)
        elif window_type == 'cosinebellwindow':
            window = CosineBellWindow(alpha=alpha)
        elif window_type == 'splitcosinebellwindow':
            window = SplitCosineBellWindow(alpha=alpha, beta=beta)
        elif window_type == 'tukeywindow':
            window = TukeyWindow(alpha=alpha)
        elif window_type == 'hanningwindow':
            window = HanningWindow()
        else:
            raise ValueError("ERROR: Unrecognized window_type: {}".format(orig_window_type))

        # Create the matching kernel
        matched_kernel = create_matching_kernel(psf1, psf2, window=window)

        return matched_kernel

    def module_value(self):
        """Determine the module value based on the instrument and detector
        """
        if self.instrument == 'NIRCAM':
            module = self.detector[3].upper()
        elif self.instrument == 'NIRISS':
            module = 'N'
        elif self.instrument == 'FGS':
            module = 'G'
        return module

    def prepare_psf(self):
        """Get metadata relevant to the instrument/detector as well as the
        pixel scale from either the mosaic file, the PSF file, or both.
        Read in or create the PSF that is representative of the mosaic
        data.
        """
        # Get the PSF associated with the mosaic
        if self.psf_file is None:

            # If no PSF file is given and there is no pixel scale listed in
            # the mosaic file, then we cannot continue
            if self.mosaic_metadata['pix_scale1'] is None:
                raise ValueError(("ERROR: pixel scale value not present in mosaic file "
                                  "(in CD1_1 header keyword). This information is needed "
                                  "to be able to convolve the mosaic with the proper PSF "
                                  "kernel."))

            # Get the dimensions of the JWST PSF kernel representing the
            # final PSF to use in the simulation
            jwst_ydim, jwst_xdim = self.jwst_psf.shape

            # If no psf file is given, and the user requests a 2D Gaussian,
            # then make a 2D Gaussian
            if self.gaussian_psf:
                print("Creating 2D Gaussian for mosiac PSF")

                # If a Gaussian FWHM value is given, then construct the
                # PSF using astropy's Gaussian2D kernel

                # Rough guess on dimensions to use
                scale_ratio1 = self.outscale1 / self.mosaic_metadata['pix_scale1']
                scale_ratio2 = self.outscale2 / self.mosaic_metadata['pix_scale2']
                gauss_xdim = int(np.round(jwst_xdim * scale_ratio1))
                gauss_ydim = int(np.round(jwst_ydim * scale_ratio2))

                # Make sure the array has an odd number of rows and columns
                if gauss_xdim % 2 == 0:
                    gauss_xdim += 1
                if gauss_ydim % 2 == 0:
                    gauss_ydim += 1

                print('Mosaic pixel scale:', self.mosaic_metadata['pix_scale1'], self.mosaic_metadata['pix_scale2'])
                print('JWST pixel scale:', self.outscale1, self.outscale2)
                print('scale ratios:', scale_ratio1, scale_ratio2)
                print('JWST PSF dims:', jwst_xdim, jwst_ydim)
                print('Gaussian PSF dimensions: ', gauss_xdim, gauss_ydim)


                # Create 2D Gaussian
                self.mosaic_psf = tools.gaussian_psf(self.mosaic_fwhm, gauss_xdim, gauss_ydim)

                print('Temporarily saving psf for development')
                h0 = fits.PrimaryHDU(self.mosaic_psf)
                hlist = fits.HDUList([h0])
                hlist.writeto('gaussian_2d_psf.fits', overwrite=True)
            elif self.mosaic_metadata['telescope'] in KNOWN_PSF_TELESCOPES and not self.gaussian_psf:
                # If no PSF file is given and the user does not want a
                # generic 2D Gaussian (self.gaussian_psf is False),
                # check to see if the mosaic is from one of the telescopes
                # with a known PSF. If so, use the appropriate function to
                # construct a PSF
                if self.mosaic_metadata['telescope'] == 'HST':
                    print('Creating HST PSF, using 2D Gaussian')
                    print('HST FWHM in arcsec: {}'.format(self.mosaic_fwhm * self.mosaic_metadata['pix_scale2']))
                    self.mosaic_psf = tools.get_HST_PSF(self.mosaic_metadata, self.mosaic_fwhm)
                elif self.mosaic_psf_metadata['telescope'] == 'JWST':
                    print("Retrieving JWST PSF")
                    self.mosaic_psf = tools.get_JWST_PSF(self.mosaic_metadata)
                elif self.mosaic_psf_metadata['instrument'] == 'IRAC':
                    print("Retrieving IRAC PSF")
                    self.mosaic_psf = tools.get_IRAC_PSF(self.mosaic_metadata)
            else:
                raise ValueError(("For telescopes other than {}, you must either provide a "
                                  "representative PSF using the psf_file keyword, or set "
                                  "gaussian_psf=True in order to construct a 2D Gaussian PSF."
                                  "This is neeeded to create the proper PSF kernal to transform "
                                  "to the appropriate JWST PSF."))
        else:
            # If a PSF file is provided, check for any metadata. Metadata
            # from the mosaic takes precidence over metadata in the PSF file.
            psf_metadata = tools.get_psf_metadata(psf_filename)

            # If the mosaic has no pixel scale info but the PSF file does,
            # use the value from the PSF file.
            if self.mosaic_metadata['pix_scale1'] is None:
                if psf_metadata['pix_scale1'] is not None:
                    self.mosaic_metadata['pix_scale1'] = psf_metadata['pix_scale1']
                else:
                    raise ValueError(("ERROR: no pixel scale value present in mosaic file nor PSF "
                                      "file metadata (in CD1_1 header keyword). This information is "
                                      "needed to be able to convolve the mosaic with the proper PSF "
                                      "kernel."))
            self.mosaic_psf = fits.getdata(psf_filename)

    def psf_convolution(self, model):
        """Convolve the cropped image with the appropriate PSF for the
        JWST detector being simulated.

        Parameters
        ----------
        model : jwst.datamodels.ImageModel
            Data model instance containing the cropped image

        Returns
        -------
        model : jwst.datamodels.ImageModel
            Data model with image convolved by the PSF
        """
        # The jwst_psf and the mosaic_psf must have the same array size
        # and the same pixel scale. First deal with the pixel scale.

        # Rescale one of the PSFs if necessary, in order to get matching pixel scales.
        # Since the matching kernel is going to be convolved with the mosaic image,
        # then it seems like we should match the PSFs at the mosaic pixel scale.
        if not np.isclose(self.outscale1, self.mosaic_metadata['pix_scale1'], atol=0., rtol=0.01):
            orig_jwst = copy.deepcopy(self.jwst_psf)
            self.jwst_psf = resize_psf(self.jwst_psf, self.outscale1, self.mosaic_metadata['pix_scale1'], order=3)

            resized_y_dim, resized_x_dim = self.jwst_psf.shape
            if ((resized_y_dim % 2 == 0) or (resized_x_dim % 2 == 0)):
                if resized_y_dim % 2 == 0:
                    new_y_dim = resized_y_dim + 1
                else:
                    new_y_dim = resized_y_dim
                if resized_x_dim % 2 == 0:
                    new_x_dim = resized_x_dim + 1
                else:
                    new_x_dim = resized_x_dim

                # Add a column/row to the resized array,
                jwst_psf_padded = np.zeros((new_y_dim, new_x_dim))
                jwst_psf_padded[0: resized_y_dim, 0: resized_x_dim] = self.jwst_psf
                # Rather than zeros, make the top row/leftmost column a
                # copy of the row/column next to it
                if new_y_dim > resized_y_dim:
                    jwst_psf_padded[-1, 0: resized_x_dim] = self.jwst_psf[-1, :]
                if new_x_dim > resized_x_dim:
                    jwst_psf_padded[0: resized_y_dim, -1] = self.jwst_psf[:, -1]
                if ((new_y_dim > resized_y_dim) and (new_x_dim > resized_x_dim)):
                    jwst_psf_padded[-1, -1] = self.jwst_psf[-1, 1]

                # Shift the source to be centered in the center pixel
                centerx, centery = centroid_2dg(jwst_psf_padded)
                jwst_shifted = shift(jwst_psf_padded, [0.5, 0.5], order=1)
                centerx, centery = centroid_2dg(jwst_shifted)
                self.jwst_psf = jwst_shifted
            else:
                jwst_psf_padded = self.jwst_psf

        jwst_shape = self.jwst_psf.shape
        mosaic_shape = self.mosaic_psf.shape
        if ((jwst_shape[0] % 2 != mosaic_shape[0] % 2) or (jwst_shape[1] % 2 != mosaic_shape[1] % 2)):
            raise ValueError(("ERROR: Mosaic PSF and JWST PSF have different shapes in terms "
                              "of odd/even numbers of rows and/or columns. Try adding or subtracting "
                              "rows/columns to the mosaic PSF. Mosaic PSF shape: {}, JWST PSF shape: {}"
                              .format(mosaic_shape, jwst_shape)))

        # Now crop either the resized JWST PSF or the mosaic PSF in
        # order to get them both to the same array size
        print("Crop PSFs to have the same array size")
        self.jwst_psf, self.mosaic_psf = tools.same_array_size(self.jwst_psf, self.mosaic_psf)

        # Now we make a matching kernel. The mosaic can then be
        # convolved with this kernel in order to adjust the PSFs to match
        # those from JWST.
        print("Create matching kernel")
        kernel = self.matching_kernel(self.mosaic_psf, self.jwst_psf, window_type='TukeyWindow',
                                      alpha=1.5, beta=1.5)

        if self.save_intermediates:
            print('Save JWST psf and matching psf in outgoing_and_matching_kernel.fits')
            ha = fits.PrimaryHDU(orig_jwst)
            h0 = fits.ImageHDU(self.jwst_psf)
            h1 = fits.ImageHDU(self.mosaic_psf)
            h2 = fits.ImageHDU(kernel)
            hlist = fits.HDUList([ha, h0, h1, h2])
            outfile = os.path.join(self.outdir,
                                   '{}_outgoing_and_matching_kernel.fits'.format(self.output_base))
            hlist.writeto(outfile, overwrite=True)

        print('Convolve image cropped from mosaic with the matching PSF kernel')
        start_time = datetime.datetime.now()
        convolved_mosaic = fftconvolve(model.data, kernel, mode='same')
        end_time = datetime.datetime.now()
        delta_time = end_time - start_time
        print("Convolution took {} seconds".format(delta_time.seconds))
        model.data = convolved_mosaic

        if self.save_intermediates:
            print('Saving convolved mosaic as convolved_mosaic.fits')
            h0 = fits.PrimaryHDU(convolved_mosaic)
            hlist = fits.HDUList([h0])
            outfile = os.path.join(self.outdir,
                                   '{}_convolved_mosaic.fits'.format(self.output_base))
            hlist.writeto(outfile, overwrite=True)

        return model

    def read_param_file(self, file):
        """Read the input yaml file into a nested dictionary

        Parameters
        ----------
        file : str
            Name of the input yaml file
        """
        if os.path.isfile(file):
            with open(file, 'r') as f:
                self.params = yaml.safe_load(f)
        else:
            raise FileNotFoundError(("ERROR: {} does not exist."
                                     .format(file)))

        self.instrument = self.params['Inst']['instrument'].upper()
        self.aperture = self.params['Readout']['array_name'].upper()
        self.detector_channel_value()

        self.crop_center_ra = self.params['Telescope']['ra']
        self.crop_center_dec = self.params['Telescope']['dec']
        self.blot_center_ra = self.params['Telescope']['ra']
        self.blot_center_dec = self.params['Telescope']['dec']
        self.blot_pav3 = self.params['Telescope']['rotation']

        self.flux_cal_file = self.params['Reffiles']['flux_cal']
        self.distortion_file = self.params['Reffiles']['astrometric']
        self.filter = self.params['Readout']['filter']
        self.pupil = self.params['Readout']['pupil']

        self.output_base = self.params['Output']['file']

        # If input is to be dispersed later, we need to
        # expand the region to be blotted
        if self.params['Output']['grism_source_image']:
            self.grism_coord_adjust()

        if self.outdir is None:
            self.outdir = self.params['Output']['directory']

        # Create a dictionary for future CRDS query
        self.reference_dictionary = crds_tools.dict_from_yaml(self.params)

    def add_options(self, parser=None, usage=None):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage, \
                                             description='Create seed image via catalogs')
        parser.add_argument("paramfile", help='File describing the input parameters and instrument settings to use. (YAML format).')
        parser.add_argument("mosaic_file", help="Fits file containing the image to create the seed image from")
        parser.add_argument("data_extension_number", help="Extension in fits file that contains mosaic image", default=0)
        parser.add_argument("wcs_extension_number", help="Extension in mosaic fits file that contains WCS information", default=0)
        parser.add_argument("cropped_file", help="Filename used to save image cropped from mosaic", default=None)
        parser.add_argument("blotted_file", help="Filename used to save resampled image", default=None)
        parser.add_argument("outdir", help="Output directory. If None, the output dir from paramfile is used.", default=None)
        parser.add_argument("psf_file", help="Name of fits file containing PSF corresponding to mosaic image", default=None)
        parser.add_argument("mosaic_fwhm", help="FWHM of PSF in mosaic image, in pixels or arcsec", default=None)
        parser.add_argument("mosaic_fwhm_units", help="Units of PSF FWHM in mosaic image. 'pixels' or 'arcsec'", default='arcsec')
        parser.add_argument("gaussian_psf", help="If True, create mosaic PSF from Gaussian2D model", default=False, action='store_true')
        parser.add_argument("save_intermediates", help="If True, save intermediate outputs", default=False, action='store_true')
        return parser


if __name__ == '__main__':

    usagestring = 'USAGE: fits_seed_image.py inputs.yaml'
    seed = Imgseed()
    parser = seed.add_options(usage=usagestring)
    args = parser.parse_args(namespace=seed)
    seed.crop_and_blot()
