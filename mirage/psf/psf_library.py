""" Generate "PSF Library" files to run with MIRaGe.

Authors
-------
    - Shannon Osborne
    - Lauren Chambers
"""

import glob
import itertools
import os
import time

import astropy.convolution
from astropy.io import fits
import numpy as np
import webbpsf


class CreatePSFLibrary:

    # Class variables for NIRCam short vs long wave information:
    nrca_short_filters = ['F070W', 'F090W', 'F115W', 'F140M', 'F150W2', 'F150W', 'F162M',
                          'F164N', 'F182M', 'F187N', 'F200W', 'F210M', 'F212N']
    nrca_long_filters = ['F250M', 'F277W', 'F300M', 'F322W2', 'F323N', 'F335M', 'F356W', 'F360M',
                         'F405N', 'F410M', 'F430M', 'F444W', 'F460M', 'F466N', 'F470N', 'F480M']

    nrca_short_detectors = ['NRCA1', 'NRCA2', 'NRCA3', 'NRCA4', 'NRCB1', 'NRCB2', 'NRCB3', 'NRCB4']
    nrca_long_detectors = ['NRCA5', 'NRCB5']

    def __init__(self, instrument, filters="all", detectors="all", num_psfs=16, psf_location=(1023, 1023),
                 add_distortion=True, fov_pixels=101, oversample=5, opd_type="requirements", opd_number=0,
                 save=True, fileloc=None, filename=None, overwrite=True, ote=None, header_addons=None,
                 **kwargs):
        """
        Description
        -----------
        Class to create a PSF library in the following format:
            For a given instrument, the output file (1 per filter) will contain a 5D array
            with axes [SCA, j, i, y, x] where SCA is the detector, (j,i) is the PSF position
            on the detector grid (integer positions) and (y,x) is the 2D PSF. This class must
            be run separately for each instrument.

        Parameters
        ----------
        instrument : str
            The name of the instrument you want to run. Can be any capitalization. Can
            only run 1 instrument at a time. This class is only set up for NIRCam,
            NIRISS, and FGS.

        filters : str or list, optional
            Which filter(s) you want to create a library for.
            Can be a string of 1 filter name, a list of filter names (as strings), or
            the default "all" will run through all the filters for that instrument.
            Spelling/capitalization must match what WebbPSF expects. See also special
            case for NIRCam. Default is "all".

        detectors : str or list, optional
            Which detector(s) you want to create a library for.
            Can be a string of 1 detector name, a list of detector names (as strings), or
            the default "all" will run through all the detectors for that instrument.
            Spelling/capitalization must match what WebbPSF expects. See also special
            case for NIRCam. Default is "all".

        num_psfs : int, optional
            The total number of fiducial PSFs to be created and saved in the files. This
            must be a square number. Default is 16. E.g. num_psfs = 16 will create a 4x4
            grid of fiducial PSFs.

        psf_location : tuple, optional
            If num_psfs = 1, then this is used to set the (y,x) location of that PSF.
            Default is (1023,1023).

        add_distortion : bool, optional
            If True, the PSF will have distortions applied: the geometric distortion from
            the detectors (using data from SIAF) and the rotation of the detectors with
            respect to the focal plane. Default is True.

        fov_pixels : int, optional
            The field of view in undersampled detector pixels used by WebbPSF when
            creating the PSFs. Default is 101 pixels.

        oversample : int, optional
            The oversampling factor used by WebbPSF when creating the PSFs. Default is 5.

        opd_type : str, optional
            The type of OPD map you would like to use to create the PSFs. Options are
            "predicted" or "requirements" where the predicted map is of the expected
            WFE and the requirements map is slightly more conservative (has slightly
            larger WFE). Default is "requirements".

        opd_number : int, optional
            The realization of the OPD map pulled from the OPD file. Options are an
            integer from 0 to 9, one for each of the 10 Monte Carlo realizations of
            the telescope included in the OPD map file. Default is 0.

        save : bool, optional
            True/False boolean if you want to save your file

        fileloc : str, optional
            Where to save your file if "save" keyword is set to True. Default of None
            will save in the current directory

        filename : str, optional
            The name to save your current file under if "save" keyword is set to True.
            Default of None will save it in the form: "{instr}_{filt}_fovp{}_samp{}_npsf{}.fits"

        ote: webbpsf.opds.OTE_Linear_Model_WSS, optional
            The OPD map to use to generate all PSFs, in the form of a
            OTE_Linear_Model_WSS object. If provided, overrides the ``opd_type`` and
            ``opd_number`` options.

        header_addons: astropy.io.fits.header.Header, optional
            Additional header entries to include in PSF library header

        overwrite : bool, optional
            True/False boolean to overwrite the output file if it already exists. Default
            is True.

        **kwargs
            This can be used to add any extra arguments to the webbpsf calc_psf() method
            call.

        Special Case for NIRCam:
        For NIRCam, you can set detectors and filters with multiple options.
        You may set both filters and detectors = "all" just like the other instruments,
        and the short and long wave filter/detectors will be separated and run in the
        correct pairings.
        If you choose only certain filters (either by name or with "shortwave" or
        "longwave"), you may set detectors to "shortwave" or "longwave" or you can
        set it to be "all" and it will pull the all appropriate SW/LW detectors.

        Use
        ---
        c = CreatePSFLibrary(instrument, filters, detectors, num_psfs, add_distortion,
                             fov_pixels, oversample, save, fileloc, filename, overwrite)
        c.create_files()

        nis = CreatePSFLibrary("NIRISS") # will run all filters/detectors
        nis.create_files()

        """

        # Pull correct capitalization of instrument name
        webbpsf_name_dict = {"NIRCAM": "NIRCam", "NIRSPEC": "NIRSpec",
                             "NIRISS": "NIRISS", "MIRI": "MIRI", "FGS": "FGS"}
        self.instr = webbpsf_name_dict[instrument.upper()]

        # Create instance of instrument in WebbPSF
        self.webb = getattr(webbpsf, self.instr)()

        # Set the filters and detectors based on the inputs
        self.filter_input = filters
        self.detector_input = detectors

        # A list of filters and a list of list of detectors (1 sublist per filter)
        self.filter_list = self._set_filters()
        self.detector_list = [self._set_detectors(filter) for filter in self.filter_list]

        # Set the locations on the detector of the fiducial PSFs
        self.ij_list, self.loc_list, self.location_list = self._set_psf_locations(num_psfs, psf_location)

        # For NIRCam: Check if filters/detectors match in terms of if they are long/shortwave
        if self.instr == "NIRCam":
            for filt, det_list in zip(self.filter_list, self.detector_list):
                for det in det_list:
                    if "5" in det and filt in CreatePSFLibrary.nrca_short_filters:
                        self._raise_value_error("short filter", det, filt)
                    elif "5" not in det and filt in CreatePSFLibrary.nrca_long_filters:
                        self._raise_value_error("long filter", det, filt)

        # Set PSF attributes
        self.add_distortion = add_distortion
        self.fov_pixels = fov_pixels
        self.oversample = oversample
        self.opd_type = opd_type
        self.opd_number = opd_number
        self.ote = ote
        self._kwargs = kwargs

        # Set saving attributes
        self.save = save
        self.overwrite = overwrite
        self.fileloc = fileloc
        self.filename = filename
        self.header_addons = header_addons

    @staticmethod
    def _raise_value_error(msg_type, det, filt):
        """Raise an error based on mis-matched short/long wave detectors/filters"""

        if "short filter" in msg_type.lower():
            message = "You are trying to apply a shortwave filter ({}) to a " \
                      "longwave detector ({}). ".format(filt, det)
        if "long filter" in msg_type.lower():
            message = "You are trying to apply a longwave filter ({}) to a " \
                      "shortwave detector ({}). ".format(filt, det)

        raise ValueError(message + "Please change these entries so the filter "
                                   "falls within the detector band.")

    def _set_detectors(self, filter):
        """Set the list of detectors to include in the PSF library files"""

        if self.detector_input == "all":
            if self.instr != "NIRCam":
                detector_list = self.webb.detector_list
            elif self.instr == "NIRCam" and filter in CreatePSFLibrary.nrca_short_filters:
                detector_list = CreatePSFLibrary.nrca_short_detectors
            elif self.instr == "NIRCam" and filter in CreatePSFLibrary.nrca_long_filters:
                detector_list = CreatePSFLibrary.nrca_long_detectors
        elif self.detector_input == "shortwave":
            detector_list = CreatePSFLibrary.nrca_short_detectors
        elif self.detector_input == "longwave":
            detector_list = CreatePSFLibrary.nrca_long_detectors
        elif type(self.detector_input) is str:
            detector_list = self.detector_input.split()
        elif type(self.detector_input) is list:
            detector_list = self.detector_input
        else:
            raise TypeError("Method of setting detectors is not valid.")

        # If the user hand chose a detector list, check it's valid for the chosen instrument
        if self.detector_input not in ["all", "shortwave", "longwave"]:
            det = set(detector_list).difference(set(self.webb.detector_list))
            if det != set():
                raise ValueError("Instrument {} doesn't have the detector(s) "
                                 "{}.".format(self.instr, det))

        return detector_list

    def _set_filters(self):
        """Set the list of filters to create PSF library files for"""

        if self.filter_input == "all":
            filter_list = self.webb.filter_list
        elif self.filter_input == "shortwave":
            filter_list = CreatePSFLibrary.nrca_short_filters
        elif self.filter_input == "longwave":
            filter_list = CreatePSFLibrary.nrca_long_filters
        elif type(self.filter_input) is str:
            filter_list = self.filter_input.split()
        elif type(self.filter_input) is list:
            filter_list = self.filter_input
        else:
            raise TypeError("Method of setting filters is not valid.")

        # If the user hand chose a filter list, check it's valid for the chosen instrument
        if self.filter_input not in ["all", "shortwave", "longwave"]:
            filt = set(filter_list).difference(set(self.webb.filter_list))
            if filt != set():
                raise ValueError("Instrument {} doesn't have the filter(s) "
                                 "{}.".format(self.instr, filt))

        return filter_list

    def _set_psf_locations(self, num_psfs, psf_location):
        """Set the locations on the detector of the fiducial PSFs. Assumes a 2048x2048 detector"""

        # The locations these PSFs should be centered on for a 2048x2048 detector
        self.num_psfs = num_psfs

        if np.sqrt(self.num_psfs).is_integer():
            self.length = int(np.sqrt(self.num_psfs))
        else:
            raise ValueError("You must choose a square number of fiducial PSFs to create (E.g. 9, 16, etc.)")

        # Set the center values
        if num_psfs == 1:
            # (1023.5, 1023.5) is the center, but we want an integer location- so default is (1023,1023)
            ij_list = [(0, 0)]
            loc_list = list(psf_location[::-1])  # list of x,y location
            location_list = [psf_location[::-1]]  # tuple of (x,y)
        else:
            ij_list = list(itertools.product(range(self.length), range(self.length)))
            loc_list = [int(round(num * 2047)) for num in np.linspace(0, 1, self.length, endpoint=True)]
            location_list = list(itertools.product(loc_list, loc_list))  # list of tuples (x,y) (for WebbPSF)

        return ij_list, loc_list, location_list

    def _set_opd(self):
        """Define the telescope's OPD, either grabbing the requested
        FITS file (requirements or predicted) or loading an OPD as a
        OTE_Linear_Model_WSS instance.

        For now, sets one OPD for all files being created.
        """
        if self.ote is None:
            # Set OPD Map (pull most recent version with self.webb.opd_list call) - always predicted then requirements
            if self.opd_type.lower() == "requirements":
                opd = self.webb.opd_list[1]
            elif self.opd_type.lower() == "predicted":
                opd = self.webb.opd_list[0]

            self.opd_name = opd
            self.opd_realization = self.opd_number

            self.webb.pupilopd = (opd, self.opd_number)

        elif not isinstance(self.ote, webbpsf.opds.OTE_Linear_Model_WSS):
            raise TypeError('Must provide pupil OPD as an OTE_Linear_Model_WSS object,'
                            'not {}'.format(type(self.ote)))

        else:
            self.webb.pupil = self.ote  # Also defines self.webb.pupilopd

            self.opd_name = "Modified OPD"
            self.opd_realization = "N/A"

    def create_files(self):
        """
        For a given instrument, the output file (1 per filter) will contain a 5D array
        with axes [SCA, j, i, y, x] where SCA is the detector, (j,i) is the PSF position
        on the detector grid (integer positions) and (y,x) is the 2D PSF.

        Returns
        -------
        Returns a list of all the hdulist objects created (a 5D array of [SCA, j, i, y, x],
        1 per filter) and saves out the library files if requested.

        """

        # Set the pupil OPD based on the input
        self._set_opd()

        # Set extension to read based on distortion choice
        if self.add_distortion:
            ext = "OVERDIST"
        else:
            ext = "OVERSAMP"

        # Create kernel to smooth pixel based on oversample
        kernel = astropy.convolution.Box2DKernel(width=self.oversample)

        # Set output mode
        self.webb.options['output_mode'] = 'Oversampled Image'

        # For every filter
        final_list = []
        for filt, det_list in zip(self.filter_list, self.detector_list):
            print("\nStarting filter: {}".format(filt))
            self.webb.filter = filt

            # Create an array to fill ([SCA, j, i, y, x])
            psf_size = self.fov_pixels * self.oversample
            psf_arr = np.empty((len(det_list), self.length, self.length, psf_size, psf_size))

            # For every detector
            for k, det in enumerate(det_list):
                print("  Running detector: {}".format(det))
                self.webb.detector = det

                # For each of the 9 locations on the detector (loc = tuple = (x,y))
                for (i, j), loc in zip(self.ij_list, self.location_list):
                    self.webb.detector_position = loc  # (X,Y) - line 286 in webbpsf_core.py

                    # Create PSF
                    psf = self.webb.calc_psf(add_distortion=self.add_distortion,
                                             fov_pixels=self.fov_pixels,
                                             oversample=self.oversample,
                                             **self._kwargs)

                    # Convolve PSF with a square kernel
                    psf_conv = astropy.convolution.convolve(psf[ext].data, kernel)

                    # Add PSF to 5D array
                    psf_arr[k, j, i, :, :] = psf_conv

            # Write header
            header = fits.Header()

            header["INSTRUME"] = (self.instr, "Instrument")
            header["FILTER"] = (filt, "Filter name")
            header["PUPILOPD"] = (self.opd_name, "Pupil OPD source name")
            header["OPD_REAL"] = (self.opd_realization, "Pupil OPD source realization from file")

            for i, det in enumerate(det_list):
                header["DETNAME{}".format(i)] = (det, "The #{} detector included in this file".format(i))

            header["FOVPIXEL"] = (self.fov_pixels, "Field of view in pixels (full array)")
            header["FOV"] = (psf[ext].header["FOV"], "Field of view in arcsec (full array) ")
            header["OVERSAMP"] = (self.oversample, "Oversampling factor for FFTs in computation")
            header["NWAVES"] = (psf[ext].header["NWAVES"], "Number of wavelengths used in calculation")

            for k, ij in enumerate(self.ij_list):  # these were originally written out in (i,j) and (x,y)

                # Even arrays are shifted by 0.5 so they are centered correctly during calc_psf computation
                # But this needs to be expressed correctly in the header
                loc = np.asarray(self.location_list[k], dtype=float)
                if self.fov_pixels % 2 == 0:
                    loc += 0.5

                header["DET_JI{}".format(k)] = (str(ij[::-1]), "The #{} PSF's (j,i) detector position".format(k))
                header["DET_YX{}".format(k)] = (str(tuple(loc[::-1])),
                                                "The #{} PSF's (y,x) detector pixel position".format(k))

            header["NUM_PSFS"] = (self.num_psfs, "The total number of fiducial PSFs")

            # The range of location values
            if self.num_psfs == 1:
                # In this case, loc_list is the single x and y value
                header["I0_X"] = (self.loc_list[0], "The x pixel value for i=0 (AXIS4)")
                header["J0_Y"] = (self.loc_list[1], "The y pixel value for j=0 (AXIS3)")
            else:
                last = len(self.loc_list) - 1
                header["I0_X"] = (self.loc_list[0], "The x pixel value for i=0 (AXIS4)")
                header["I{}_X".format(last)] = (self.loc_list[-1],
                                                "The x pixel value for i={} (final value; AXIS4)".format(last))
                header["J0_Y"] = (self.loc_list[0], "The y pixel value for j=0 (AXIS3)")
                header["J{}_Y".format(last)] = (self.loc_list[-1],
                                                "The y pixel value for j={} (final value; AXIS3)".format(last))

            # Distortion information
            if self.add_distortion:
                header["ROTATION"] = (psf[ext].header["ROTATION"], "PSF rotated to match detector rotation")
                header["DISTORT"] = (psf[ext].header["DISTORT"], "SIAF distortion coefficients applied")
                header["SIAF_VER"] = (psf[ext].header["SIAF_VER"], "SIAF PRD version used")

                for key in list(psf[ext].header.keys()):
                    if "COEF_" in key:
                        header[key] = (psf[ext].header[key], "SIAF distortion coefficient for {}".format(key))

            # Pull values from the last made psf
            header["WAVELEN"] = (psf[ext].header["WAVELEN"], "Weighted mean wavelength in meters")
            header["DIFFLMT"] = (psf[ext].header["DIFFLMT"], "Diffraction limit lambda/D in arcsec")
            header["FFTTYPE"] = (psf[ext].header["FFTTYPE"], "Algorithm for FFTs: numpy or fftw")
            header["NORMALIZ"] = (psf[ext].header["NORMALIZ"], "PSF normalization method")
            header["JITRTYPE"] = (psf[ext].header["JITRTYPE"], "Type of jitter applied")
            header["JITRSIGM"] = (psf[ext].header["JITRSIGM"], "Gaussian sigma for jitter [arcsec]")
            header["TEL_WFE"] = (psf[ext].header["TEL_WFE"], "[nm] Telescope pupil RMS wavefront error")

            header["DATE"] = (psf[ext].header["DATE"], "Date of calculation")
            header["AUTHOR"] = (psf[ext].header["AUTHOR"], "username@host for calculation")
            header["VERSION"] = (psf[ext].header["VERSION"], "WebbPSF software version")
            header["DATAVERS"] = (psf[ext].header["DATAVERS"], "WebbPSF reference data files version ")

            # Include any header add-ons from the user, if applicable
            if self.header_addons is not None:
                # Append new keywords on to the existing header
                header += self.header_addons

            # Add descriptor for how the file was made
            # the output file (1 per filter) will contain a 5D array with axes[SCA, j, i, y, x]
            header["COMMENT"] = "For a given instrument, the output file (1 per filter) will contain "
            header["COMMENT"] = "a 5D array with axes[SCA, j, i, y, x] where (j,i) is the PSF position "
            header["COMMENT"] = "on the detector grid (integer positions) and (y,x) is the 2D PSF. The "
            header["COMMENT"] = "order of the detectors can be found under the  header DETNAME* "
            header["COMMENT"] = "keywords and the order of the fiducial PSFs ((j,i) and (y,x)) under "
            header["COMMENT"] = "the header DET_JI*/DET_YX* keywords"

            # Add header labels
            header.insert("INSTRUME", ('', ''))
            header.insert("INSTRUME", ('COMMENT', '/ PSF Library Information'))

            header.insert("NORMALIZ", ('', ''))
            header.insert("NORMALIZ", ('COMMENT', '/ WebbPSF Creation Information'))

            if self.header_addons is None:
                header.insert("DATAVERS", ('COMMENT', '/ File Description'), after=True)
                header.insert("DATAVERS", ('', ''), after=True)
            else:
                header.insert("DATAVERS", ('COMMENT', '/ User-Supplied Information'), after=True)
                header.insert("DATAVERS", ('', ''), after=True)

                last_added_keyword = list(self.header_addons.keys())[-1]
                header.insert(last_added_keyword, ('COMMENT', '/ File Description'), after=True)
                header.insert(last_added_keyword, ('', ''), after=True)

            # Combine the header and data
            hdu = fits.HDUList([fits.PrimaryHDU(psf_arr, header=header)])

            # Write file out
            if self.save:

                if self.fileloc is None:
                    self.fileloc = os.path.expandvars('$MIRAGE_DATA/{}/'
                                                      'test_webbpsf_library'.format(self.instr.lower()))
                if self.filename is None:
                    name = "{}_{}_fovp{}_samp{}_npsf{}.fits".format(self.instr.lower(), filt.lower(),
                                                                    self.fov_pixels, self.oversample,
                                                                    self.num_psfs)
                    filepath = os.path.join(self.fileloc, name)
                else:
                    filepath = os.path.join(self.fileloc, self.filename)

                print("  Saving file: {}".format(filepath))

                hdu.writeto(filepath, overwrite=self.overwrite)

            final_list.append(hdu)

        return final_list


def generate_segment_psfs(ote, segment_tilts, out_dir, filters=['F212N', 'F480M'],
                          detectors='all', fov_pixels=1024, overwrite=False):
    """Generate NIRCam PSF libraries for all 18 mirror segments given a perturbed OTE
    mirror state. Saves each PSF library as a FITS file named in the following format:
        nircam_{filter}_fovp{fov size}_samp1_npsf1_seg{segment number}.fits

    Parameters
    ----------
    ote : webbpsf.opds.OTE_Linear_Model_WSS object
        WebbPSF OTE object describing perturbed OTE state with tip and tilt removed

    segment_tilts : numpy.ndarray
        List of X and Y tilts for each mirror segment, in microradians

    out_dir : str
        Directory in which to save FITS files

    filters : str or list, optional
        Which filters to generate PSF libraries for. Default is ['F212N', 'F480M']
        (the two filters used for most commissioning activities).

    detectors : str or list, optional
        Which detectors to generate PSF libraries for. Default is 'all'.

    fov_pixels : int, optional
        Size of the PSF to generate, in pixels. Default is 1024.

    overwrite : bool, optional
            True/False boolean to overwrite the output file if it already
            exists. Default is True.

    """
    for i in range(18):
        start_time = time.time()

        # Restrict the pupil to the current segment
        i_segment = i + 1
        segname = webbpsf.webbpsf_core.segname(i_segment)
        print('GENERATING SEGMENT {} DATA'.format(segname))
        print('------------------------------')

        pupil = webbpsf.webbpsf_core.one_segment_pupil(i_segment)
        ote.amplitude = pupil[0].data

        # Add header keywords about segment
        hdr = fits.Header()
        hdr['SEGID'] = (i_segment, 'ID of the mirror segment')
        hdr['SEGNAME'] = (segname, 'Name of the mirror segment')
        hdr['XTILT'] = (round(segment_tilts[i, 0], 2), 'X tilt of the segment in microns')
        hdr['YTILT'] = (round(segment_tilts[i, 1], 2), 'Y tilt of the segment in microns')

        # Generate the library file(s)
        c = CreatePSFLibrary('NIRCam', filters=filters, detectors=detectors,
                             fov_pixels=fov_pixels, oversample=1, num_psfs=1,
                             fileloc=out_dir, ote=ote, overwrite=overwrite,
                             header_addons=hdr)
        c.create_files()

        # Update file names to include segment number
        created_files = glob.glob(os.path.join(out_dir, 'nircam_*_samp1_npsf1.fits'))
        for f in created_files:
            old_name = os.path.basename(f)
            old_root = old_name.split('.fits')[0]
            new_path = os.path.join(out_dir, old_root + '_seg{:02d}.fits'.format(i_segment))
            os.rename(f, new_path)
            print('  Renamed file:', new_path)

        print('Elapsed time:', time.time() - start_time)
        print()
