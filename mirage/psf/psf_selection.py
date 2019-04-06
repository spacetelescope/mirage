#! /usr/bin/env python

from ast import literal_eval
from glob import glob
import os

from astropy.io import fits
import numpy as np
#from photutils.utils import ShepardIDWInterpolator as idw
from photutils import FittableImageModel
from scipy.interpolate import interp2d, RectBivariateSpline
from webbpsf.utils import to_griddedpsfmodel


def get_gridded_psf_library(instrument, detector, filtername, pupilname, width, oversamp,
                            number_of_psfs, wavefront_error, wavefront_error_group, library_path):
    """
    Find the filename for the appropriate gridded PSF library and read it in
    """
    library_file = get_library_file(instrument, detector, filtername, pupilname,
                                    wavefront_error, wavefront_error_group, library_path)
    print("PSFs will be generated using: {}".format(os.path.basename(library_file)))

    try:
        library = to_griddedpsfmodel(library_file)
    except OSError:
        print("OSError: Unable to open {}.".format(library_file))
    return library


def get_library_file_old(instrument, detector, filt, pupil, fov_pix, oversample, num_psf, wfe,
                         wfe_group, library_path):
        """Given an instrument and filter name along with the path of
        the PSF library, find the appropriate library file to load.

        Parameters:
        -----------
        instrument : str
            Name of instrument the PSFs are from

        detector : str
            Name of the detector within ```instrument```

        filt : str
            Name of filter used for PSF library creation

        pupil : str
            Name of pupil wheel element used for PSF library creation

        fov_pix : int
            With of the PSF stamps in units of nominal pixels

        oversample : int
            Oversampling factor of the PSF stamps. (e.g oversamp=2 means
            NIRCam SW PSFs of 0.031 / 2 arcsec per pixel)

        num_psf : int
            Number of PSFs across the detector in the library. (e.g. for
            a 3x3 library of PSFs, num_psf=9)

        wfe : str
            Wavefront error. Can be 'predicted' or 'requirements'

        wfe_group : int
            Wavefront error realization group. Must be an integer from 0 - 9.

        library_path : str
            Path pointing to the location of the PSF library

        Returns:
        --------
        lib_file : str
            Name of the PSF library file for the instrument and filtername
        """
        filename = '{}_{}_{}_{}_fovp{}_samp{}_npsf{}_wfe_{}_wfegroup{}.fits'.format(instrument.lower(),
                                                                                   detector.lower(),
                                                                                   filt.lower(),
                                                                                   pupil.lower(),
                                                                                   fov_pix, oversample,
                                                                                   num_psf, wfe, wfe_group)
        print('PSF file to use: {}'.format(filename))
        lib_file = os.path.join(library_path, filename)

        # If no matching files are found, or more than 1 matching file is
        # found, raise an error.
        if not os.path.isfile:
            raise FileNotFoundError("PSF library file {} does not exist."
                                    .format(lib_file))
        return lib_file


def get_library_file(instrument, detector, filt, pupil, wfe, wfe_group, library_path):
    """Given an instrument and filter name along with the path of
        the PSF library, find the appropriate library file to load.

        Parameters:
        -----------
        instrument : str
            Name of instrument the PSFs are from

        detector : str
            Name of the detector within ```instrument```

        filt : str
            Name of filter used for PSF library creation

        pupil : str
            Name of pupil wheel element used for PSF library creation

        wfe : str
            Wavefront error. Can be 'predicted' or 'requirements'

        wfe_group : int
            Wavefront error realization group. Must be an integer from 0 - 9.

        library_path : str
            Path pointing to the location of the PSF library

        Returns:
        --------
        lib_file : str
            Name of the PSF library file for the instrument and filtername
    """
    psf_files = glob(os.path.join(library_path, '*.fits'))

    # Create a dictionary of header information for all PSF library files
    psf_table = {}
    matches = []

    instrument = instrument.upper()
    detector = detector.upper()
    filt = filt.upper()
    pupil = pupil.upper()
    wfe = wfe.lower()

    for filename in psf_files:
        header = fits.getheader(filename)
        file_inst = header['INSTRUME'].upper()
        file_det = header['DETECTOR'].upper()
        file_filt = header['FILTER'].upper()
        #file_pupil = header['PUPIL'].upper()
        file_pupil = 'CLEAR'
        print('PUPIL VALUE SET TO CLEAR WHILE AWAITING KEYWORD')
        opd = header['PUPILOPD']
        if 'requirements' in opd:
            file_wfe = 'requirements'
        elif 'predicted' in opd:
            file_wfe = 'predicted'

        if 'slice' in opd:
            file_wfe_grp = np.int(opd.split(' ')[-1])
        else:
            file_wfe_grp = 0

        match = (file_inst == instrument and file_det == detector and file_filt == filt and
                 file_pupil == pupil and file_wfe == wfe and file_wfe_grp == wfe_group)

        if match:
            matches.append(filename)
        psf_table[filename] = [file_inst, file_det, file_filt, file_pupil, file_wfe, file_wfe_grp, match]

    # Find files matching the requested inputs
    if len(matches) == 1:
        return matches[0]
    elif len(matches) == 0:
        raise ValueError("No PSF library file found matching requested parameters.")
    elif len(matches) > 1:
        raise ValueError("More than one PSF library file matches requested parameters: {}".format(matches))
