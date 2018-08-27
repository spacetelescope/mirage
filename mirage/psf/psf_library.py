import os
import itertools

import webbpsf
import numpy as np
import astropy.convolution
from astropy.io import fits


class CreatePSFLibrary:
    """
    Class Description:
    -----------------
    Class to create a PSF library in the following format:
        For a given instrument, 1 file per filter in the form [SCA, j, i, y, x] where
        (j,i) is the PSF position on the detector grid (integer positions) and (y,x)
        is the 2D PSF.

    This class must be run separately for each instrument.

    Special Case for NIRCam:
    For NIRCam, you can set detectors and filters both to "shortwave" or both to
    "longwave" to run all the shortwave/longwave filters/detectors.
    Or if you want to run all of the filters and detectors, setting both to "all"
    will separate the short and long wave filter/detectors and run them in the
    correct pairings.
    If you choose to specify specific filters and detectors, they must all be either
    shortwave or longwave. Mismatched lists of short and long wave filters and
    detectors will result in an error.

    Parameters:
    -----------
    instrument: str
        The name of the instrument you want to run. Can be any capitalization. Can
        only run 1 instrument at a time. Right now this is only set up for NIRCam,
        NIRISS, and FGS (they are 2048x2048)

    filters: str or list
        Which filter(s) you want to create a library for.

        Can be a string of 1 filter name, a list of filter names (as strings), or
        the default "all" will run through all the filters in the filter_list
        attribute of webbpsf.INSTR(). Spelling/capitalization must match what
        webbpsf expects. See also special case for NIRCam. Default is "all"

    detectors: str or list
        Which detector(s) you want to create a library for.

        Can be a string of 1 detector name, a list of detector names (as strings), or
        the default "all" will run through all the detectors in the detector_list
        attribute of webbpsf.INSTR(). Spelling/capitalization must match what
        webbpsf expects. See also special case for NIRCam. Default is "all"

    fov_pixels: int
        The field of view in undersampled detector pixels used by WebbPSF when
        creating the PSFs. Default is 101 pixels.

    oversample: int
        The oversampling factor used by WebbPSF when creating the PSFs. Default is 5.

    num_psfs: int
        The total number of fiducial PSFs to be created and saved in the files. This
        number must be a square number. Default is 16.
        E.g. num_psfs = 16 will have the class create a 4x4 grid of fiducial PSFs.

    opd_type: str
        The type of OPD map you would like to use to create the PSFs. Options are
        "predicted" or "requirements" where the predicted map is of the expected
        WFE and the requirements map is slightly more conservative (has slightly
        larger WFE). Default is "requirements"

    opd_number: int
        The realization of the OPD map pulled from the OPD file. Options are an
        integer from 0 to 9, one for each of the 10 Monte Carlo realizations of
        the telescope included in the OPD map file. Default is 0.

    save: bool
        True/False boolean if you want to save your file

    fileloc: str
        Where to save your file if "save" keyword is set to True. Default of None
        will save in the current directory

    filename: str
        The name to save your current file under if "save" keyword is set to True.
        Default of None will save it as "INSTRNAME_FILTERNAME.fits"

    overwrite: bool
        True/False boolean to overwrite the output file if it already exists.


    Use:
    ----
    c = CreatePSFLibrary(instrument, filters, detectors, fov_pixels, oversample,
                         num_psfs, save, fileloc, filename, overwrite)
    c.create_files()

    nis = CreatePSFLibrary("NIRISS") # will run all filters/detectors
    nis.create_files()

    """

    # Class variables for NIRCam short vs long wave information:
    nrca_short_filters = ['F070W', 'F090W', 'F115W', 'F140M', 'F150W2', 'F150W', 'F162M', 'F164N', 'F182M', 'F187N',
                          'F200W', 'F210M', 'F212N']
    nrca_long_filters = ['F250M', 'F277W', 'F300M', 'F322W2', 'F323N', 'F335M', 'F356W', 'F360M', 'F405N', 'F410M',
                         'F430M', 'F444W', 'F460M', 'F466N', 'F470N', 'F480M']

    nrca_short_detectors = ['NRCA1', 'NRCA2', 'NRCA3', 'NRCA4', 'NRCB1', 'NRCB2', 'NRCB3', 'NRCB4']
    nrca_long_detectors = ['NRCA5', 'NRCB5']

    def _set_filters(self):
        """ Get the list of filters to create PSF library files for """

        # Set filter list to loop over
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

        # If the user hand chose a filter list, check it's a valid list for the chosen instrument
        if self.filter_input not in ["all", "shortwave", "longwave"]:
            for filt in filter_list:
                if filt not in self.webb.filter_list:
                    raise ValueError("Instrument {} doesn't have a filter called {}.".format(self.instr, filt))

        return filter_list

    def _set_detectors(self):
        """ Get the list of detectors to include in the PSF library files """

        # Set detector list to loop over
        if self.detector_input == "all":
            detector_list = self.webb.detector_list
        elif self.detector_input == "shortwave":
            detector_list = CreatePSFLibrary.nrca_short_detectors
        elif self.detector_input == "longwave":
            detector_list = CreatePSFLibrary.nrca_long_detectors
        elif type(self.detector_input) is str:
            detector_list = self.detector_input.split()
        elif type(self.detector_input) is list:
            detector_list = self.detector_input

        # If the user hand chose a detector list, check it's a valid list for the chosen instrument
        if self.detector_input not in ["all", "shortwave", "longwave"]:
            for det in detector_list:
                if det not in self.webb.detector_list:
                    raise ValueError("Instrument {} doesn't have a detector called {}.".format(self.instr, det))

        return detector_list

    @staticmethod
    def _raise_value_error(msg_type, det, filt):
        """Raise a specific ValueError based on mis-matched short/long wave detectors/filters"""

        if "short filter" in msg_type.lower():
            message = "You are trying to apply a shortwave filter ({}) to a longwave detector ({}). ".format(filt, det)
        if "long filter" in msg_type.lower():
            message = "You are trying to apply a longwave filter ({}) to a shortwave detector ({}). ".format(filt, det)

        raise ValueError(message + "Please change these entries so the filter falls within the detector band.")

    def _set_psf_locations(self, num_psfs):
        """ Set the locations on the detector of the fiducial PSFs. Assumes a 2048x2048 detector"""

        # The locations these PSFs should be centered on for a 2048x2048 detector
        self.num_psfs = num_psfs

        if np.sqrt(self.num_psfs).is_integer():
            self.length = int(np.sqrt(self.num_psfs))
        else:
            raise ValueError("You must choose a square number of fiducial PSFs to create (E.g. 9, 16, etc.)")

        # Set the values
        ij_list = list(itertools.product(range(self.length), range(self.length)))
        self.loc_list = [int(round(num * 2047)) for num in np.linspace(0, 1, self.length, endpoint=True)]
        location_list = list(itertools.product(self.loc_list, self.loc_list))  # list of tuples

        return ij_list, location_list

    def __init__(self, instrument, filters="all", detectors="all", fov_pixels=101, oversample=5, num_psfs=16,
                 opd_type="requirements", opd_number=0, save=True, fileloc=None, filename=None, overwrite=True):

        # Pull correct capitalization of instrument name
        webbpsf_name_dict = {"NIRCAM": "NIRCam", "NIRSPEC": "NIRSpec", "NIRISS": "NIRISS",
                             "MIRI": "MIRI", "FGS": "FGS"}

        self.instr = webbpsf_name_dict[instrument.upper()]

        # Create instance of instrument in WebbPSF (same as webbpsf.instr)
        self.webb = getattr(webbpsf, self.instr)()

        # Set the filters and detectors based on the inputs
        self.filter_input = filters
        self.detector_input = detectors

        self.filter_list = self._set_filters()
        self.detector_list = self._set_detectors()

        # Set the locations on the detector of the fiducial PSFs
        self.ij_list, self.location_list = self._set_psf_locations(num_psfs)

        # For NIRCam: Check if filters/detectors match in terms of if they are longwave/shortwave
        # The "all" case will be updated later
        if self.instr == "NIRCam" and (self.filter_input, self.detector_input) != ("all", "all"):
            for det in self.detector_list:
                if "5" in det:
                    [self._raise_value_error("short filter", det, filt) for filt in
                     self.filter_list if filt in CreatePSFLibrary.nrca_short_filters]

                else:
                    [self._raise_value_error("long filter", det, filt) for filt in
                     self.filter_list if filt in CreatePSFLibrary.nrca_long_filters]

        # Set PSF attributes
        self.fov_pixels = fov_pixels
        self.oversample = oversample
        self.opd_type = opd_type
        self.opd_number = opd_number
        self.save = save
        self.overwrite = overwrite
        self.fileloc = fileloc
        self.filename = filename

    def create_files(self):
        """
        This method creates the following:

        For a given instrument, 1 file per filter in the form [SCA, j, i, y, x] where
        (j,i) is the PSF position on the detector grid (integer positions) and (y,x)
        is the 2D PSF.

        All variables needed to run this method are defined during the creation for the
        instance

        Returns:
        -------
        This saves out the library files if requested and then returns a list of all the
        hdulist objects created (each in the form of [SCA, j, i, y, x], 1 per filter
        requested).

        """

        # If someone wants to run all of NIRCam, run the function 2x: shortwave and longwave
        if self.instr == "NIRCam" and (self.filter_input, self.detector_input) == ("all", "all"):
            short_list = self._run_files(CreatePSFLibrary.nrca_short_filters, CreatePSFLibrary.nrca_short_detectors)
            long_list = self._run_files(CreatePSFLibrary.nrca_long_filters, CreatePSFLibrary.nrca_long_detectors)

            final_list = short_list + long_list

        else:
            final_list = self._run_files(self.filter_list, self.detector_list)

        return final_list

    def _run_files(self, filter_list, detector_list):
        """
        This method is called in the create_files() method

        For a given instrument, 1 file per filter in the form [SCA, j, i, y, x] where
        (j,i) is the PSF position on the detector grid (integer positions) and (y,x)
        is the 2D PSF.

        Parameters filter_list and detector_list are defined in the create_files method

        Returns:
        -------
        This saves out the library files if requested and then returns a list of all the
        hdulist objects created (each in the form of [SCA, j, i, y, x], 1 per filter
        requested).

        """

        # Create kernel to smooth pixel based on oversample
        kernel = astropy.convolution.Box2DKernel(width=self.oversample)

        # Set output mode
        self.webb.options['output_mode'] = 'Oversampled Image'

        # Set OPD Map (pull most recent version with self.webb.opd_list call) - always predicted then requirements
        if self.opd_type.lower() == "requirements":
            opd = self.webb.opd_list[1]
        elif self.opd_type.lower() == "predicted":
            opd = self.webb.opd_list[0]
        self.webb.pupilopd = (opd, self.opd_number)

        # For every filter
        final_list = []
        for filt in filter_list:
            print("\nStarting filter: {}".format(filt))

            # Set filter
            self.webb.filter = filt

            # Create an array to fill ([SCA, j, i, y, x])
            psf_size = self.fov_pixels * self.oversample
            psf_arr = np.empty((len(detector_list), self.length, self.length, psf_size, psf_size))

            # For every detector
            for k, det in enumerate(detector_list):
                print("  Running detector {}".format(det))

                self.webb.detector = det

                # For each of the 9 locations on the detector (loc = tuple = (x,y))
                for (i, j), loc in zip(self.ij_list, self.location_list):
                    self.webb.detector_position = loc  # (X,Y) - line 286 in webbpsf_core

                    # Create PSF
                    psf = self.webb.calc_psf(fov_pixels=self.fov_pixels, oversample=self.oversample)

                    # Convolve PSF with a square kernel
                    psf_conv = astropy.convolution.convolve(psf["OVERDIST"].data, kernel)

                    # Add PSF to 5D array
                    psf_arr[k, j, i, :, :] = psf_conv

            # Write header
            header = fits.Header()

            header["INSTRUME"] = (self.instr, "Instrument")
            header["FILTER"] = (filt, "Filter name")
            header["PUPILOPD"] = (self.webb.pupilopd[0], "Pupil OPD source name")
            header["OPD_REAL"] = (self.webb.pupilopd[1], "Pupil OPD source realization from file")

            for i, det in enumerate(detector_list):
                header["DETNAME{}".format(i)] = (det, "The #{} detector included in this file".format(i))

            header["FOVPIXEL"] = (self.fov_pixels, "Field of view in pixels (full array)")
            header["OVERSAMP"] = (self.oversample, "Oversampling factor for FFTs in computation")

            for k, ij in enumerate(self.ij_list):  # these were originally written out in (i,j) and (x,y)
                header["DET_JI{}".format(k)] = (str((ij[1], ij[0])), "The #{} PSF's (j,i) detector position".format(k))
                header["DET_YX{}".format(k)] = (str((self.location_list[k][1], self.location_list[k][0])),
                                                "The #{} PSF's (y,x) detector pixel position".format(k))

            header["NUM_PSFS"] = (self.num_psfs, "The total number of fiducial PSFs")

            last = len(self.loc_list) - 1
            header["I0_X"] = (self.loc_list[0], "The x pixel value for i=0 (AXIS4)")
            header["I{}_X".format(last)] = (self.loc_list[-1],
                                            "The x pixel value for i={} (final value; AXIS4)".format(last))
            header["J0_Y"] = (self.loc_list[0], "The y pixel value for j=0 (AXIS3)")
            header["J{}_Y".format(last)] = (self.loc_list[-1],
                                            "The y pixel value for j={} (final value; AXIS3)".format(last))

            # Pull values from the last made psf
            header["NORMALIZ"] = (psf[0].header["NORMALIZ"], "PSF normalization method")
            header["DATE"] = (psf[0].header["DATE"], "Date of calculation")
            header["AUTHOR"] = (psf[0].header["AUTHOR"], "username@host for calculation")
            header["VERSION"] = (psf[0].header["VERSION"], "WebbPSF software version")
            header["DATAVERS"] = (psf[0].header["DATAVERS"], "WebbPSF reference data files version ")

            # Add descriptor for how the file was made
            header["COMMENT"] = "For a given instrument, 1 file per filter in the form [SCA, j, i, y, x]"
            header["COMMENT"] = "where (j,i) is the PSF position on the detector grid (integer "
            header["COMMENT"] = "positions) and (y,x) is the 2D PSF. The order of the detectors can be "
            header["COMMENT"] = "found under the  header DETNAME* keywords and the order of the fiducial "
            header["COMMENT"] = "PSFs ((j,i) and (y,x)) under the header DET_JI*/DET_YX* keywords"

            # Add header labels
            header.insert("INSTRUME", ('', ''))
            header.insert("INSTRUME", ('COMMENT', '/ PSF Library Information'))

            header.insert("NORMALIZ", ('', ''))
            header.insert("NORMALIZ", ('COMMENT', '/ WebbPSF Creation Information'))

            header.insert("DATAVERS", ('COMMENT', '/ File Description'), after=True)
            header.insert("DATAVERS", ('', ''), after=True)

            # Combine the header and data
            hdu = fits.HDUList([fits.PrimaryHDU(psf_arr, header=header)])

            # Write file out
            if self.save:

                # Set file information
                if self.fileloc is None:
                    self.fileloc = os.path.expandvars('$MIRAGE_DATA/{}/test_webbpsf_library'.format(self.instr.lower()))

                if self.filename is None:
                    # E.g. filename: nircam_f090w_fovp1000_samp5_npsf16.fits
                    name = "{}_{}_fovp{}_samp{}_npsf{}.fits".format(self.instr.lower(), filt.lower(), self.fov_pixels,
                                                                    self.oversample, self.num_psfs)
                    self.filepath = os.path.join(self.fileloc, name)
                else:
                    self.filepath = os.path.join(self.fileloc, self.filename)

                print("  Saving file: {}".format(self.filepath))

                hdu.writeto(self.filepath, overwrite=self.overwrite)

            # Create something to return
            final_list.append(hdu)

        return final_list
