#! /usr/bin/env python

from copy import copy
from glob import glob
import os

from astropy.io import fits
import numpy as np
from webbpsf.utils import to_griddedpsfmodel

from mirage.utils.constants import NIRISS_PUPIL_WHEEL_FILTERS


def get_gridded_psf_library(instrument, detector, filtername, pupilname, wavefront_error,
                            wavefront_error_group, library_path):
    """Find the filename for the appropriate gridded PSF library and
    read it in to a griddedPSFModel

    Parameters
    ----------
    instrument : str
        Name of instrument the PSFs are from

    detector : str
        Name of the detector within ```instrument```

    filtername : str
        Name of filter used for PSF library creation

    pupilname : str
        Name of pupil wheel element used for PSF library creation

    wavefront_errpr : str
        Wavefront error. Can be 'predicted' or 'requirements'

    wavefront_error__group : int
        Wavefront error realization group. Must be an integer from 0 - 9.

    library_path : str
        Path pointing to the location of the PSF library

    Returns:
    --------
    library : photutils.griddedPSFModel
        Object containing PSF library

    """
    library_file = get_library_file(instrument, detector, filtername, pupilname,
                                    wavefront_error, wavefront_error_group, library_path)
    print("PSFs will be generated using: {}".format(os.path.basename(library_file)))

    try:
        library = to_griddedpsfmodel(library_file)
    except OSError:
        print("OSError: Unable to open {}.".format(library_file))
    return library


def get_library_file(instrument, detector, filt, pupil, wfe, wfe_group, library_path):
    """Given an instrument and filter name along with the path of
    the PSF library, find the appropriate library file to load.

    Parameters
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

    Returns
    --------
    lib_file : str
        Name of the PSF library file for the instrument and filtername
    """
    psf_files = glob(os.path.join(library_path, '*.fits'))

    # Create a dictionary of header information for all PSF library files
    # psf_table = {}
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
        if file_inst.upper() == 'NIRCAM':
            file_pupil = 'CLEAR'
        elif file_inst.upper() == 'NIRISS':
            file_pupil = 'CLEARP'
        print('PUPIL VALUE SET TO CLEAR WHILE AWAITING KEYWORD')

        # NIRISS has many filters in the pupil wheel. Webbpsf does
        # not make a distinction, but Mirage does. Adjust the info
        # to match Mirage's expectations
        if file_inst.upper() == 'NIRISS' and file_filt in NIRISS_PUPIL_WHEEL_FILTERS:
            save_filt = copy(file_filt)
            if file_pupil == 'CLEARP':
                file_filt = 'CLEAR'
            else:
                raise ValueError(('Pupil value is something other than '
                                  'CLEARP, but the filter being used is '
                                  'in the pupil wheel.'))
            file_pupil = save_filt

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
        # psf_table[filename] = [file_inst, file_det, file_filt, file_pupil, file_wfe, file_wfe_grp, match]

    # Find files matching the requested inputs
    if len(matches) == 1:
        return matches[0]
    elif len(matches) == 0:
        raise ValueError("No PSF library file found matching requested parameters.")
    elif len(matches) > 1:
        raise ValueError("More than one PSF library file matches requested parameters: {}".format(matches))
