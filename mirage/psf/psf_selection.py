#! /usr/bin/env python

"""This module contains code to locate the appropriate PSF library files
to use for a given simulation. It supports the selection of one PSF
"core" file, which is assumed to contain a 3D array of PSFs that is read
into a ``griddedPSFmodel`` instance, as well as a single PSF "wings" file,
which contains a single PSF instance.

For both of these files, once they are identified, they are read in via
the appropriate mechanisms for the data they contain, and the resulting
object is returned.

Author
------

    - Bryan Hilbert

Use
---

    This module can be imported and called as such:
    ::
        from mirage.psf import psf_selection
        library = psf_selection.get_gridded_psf_library('nircam', 'nrcb1',
                                                        'f200w', 'clear',
                                                        'predicted', 0,
                                                        '/path/to/library/')
"""


from copy import copy
from glob import glob
import logging
import os
import warnings
import pdb

from astropy.io import fits
import numpy as np
from webbpsf.utils import to_griddedpsfmodel

from mirage.logging import logging_functions
from mirage.utils.constants import NIRISS_PUPIL_WHEEL_FILTERS, NIRCAM_PUPIL_WHEEL_FILTERS, \
                                   LOG_CONFIG_FILENAME, STANDARD_LOGFILE_NAME, PSF_NORM_MIN, \
                                   PSF_NORM_MAX, NIRISS_NRM_PSF_THROUGHPUT_REDUCTION, \
                                   NIRISS_CLEARP_PSF_THROUGHPUT_REDUCTION, NIRCAM_WLP12_PSF_THROUGHPUT_REDUCTION
from mirage.utils.utils import expand_environment_variable, standardize_filters

classdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../'))
log_config_file = os.path.join(classdir, 'logging', LOG_CONFIG_FILENAME)
logging_functions.create_logger(log_config_file, STANDARD_LOGFILE_NAME)

def check_normalization(lib, lower_limit=0.80, upper_limit=1.0):
    """Check that the gridded PSF library is properly normalized. We expect
    the total signal of the PSF to be roughly 1.0 (minus up to several percent
    since it should be normalized to 1.0 at the pupil).

    Parameters
    ----------
    lib : photutils.psf.models.GriddedPSFModel
        Gridded PSF model instance

    lower_limit : float
        Lower limit for the total signal in the PSF

    upper_limit : float
        Upper limit for the total signal in the PSF

    Returns
    -------
    result : tup
        2-element tuple. The first is a boolean, which is
        True if the normalization is as expected. False otherwise.
        The second element is a short string describing the
        result.
    """
    ndims = len(lib.data.shape)
    if ndims == 3:
        total_signal = np.sum(lib.data[0, :, :])
    elif ndims == 2:
        total_signal = np.sum(lib.data)
    total_signal /= lib.meta['oversamp'][0]**2

    if total_signal > upper_limit:
        result = False, 'too high'
    elif total_signal < lower_limit:
        result = False, 'too low'
    else:
        result = True, 'correct'
    return result


def confirm_gridded_properties(filename, instrument, detector, filtername, pupilname,
                               wavefront_error_type, wavefront_error_group, file_path,
                               extname='PRIMARY'):
    """Examine the header of the gridded PSF model file to confirm that
    the properties of the data match those expected.

    Parameters
    ----------
    filename : str
        Base name of the PSF library file to be checked

    instrument : str
        Name of instrument the PSFs are from

    detector : str
        Name of the detector within ```instrument```

    filtername : str
        Name of filter used for PSF library creation

    pupilname : str
        Name of pupil wheel element used for PSF library creation

    wavefront_error_type : str
        Wavefront error. Can be 'predicted' or 'requirements'

    wavefront_error_group : int
        Wavefront error realization group. Must be an integer from 0 - 9.

    file_path : str
        Path pointing to the location of the PSF library

    extname : str
        Name of the extension within ``filename`` to check

    Returns
    -------
    full_filename : str
        Full path and filename if the file properties are as expected.
        None if the properties do not match.
    """

    # Determine if the PSF path is default or not
    mirage_dir = expand_environment_variable('MIRAGE_DATA')
    default_psf = file_path == os.path.join(mirage_dir,
                                               '{}/gridded_psf_library'.format(instrument.lower()))

    full_filename = os.path.join(file_path, filename)
    with fits.open(full_filename) as hdulist:
        header = hdulist[extname.upper()].header

    inst = header['INSTRUME']
    try:
        det = header['DETECTOR']
    except KeyError:
        det = header['DET_NAME']
    filt = header['FILTER']
    try:
        pupil = header['PUPIL']
    except KeyError:
        # If no pupil mask value is present, then assume the CLEAR is
        # being used
        if instrument.upper() == 'NIRCAM':
            pupil = 'CLEAR'
        elif instrument.upper() == 'NIRISS':
            pupil = 'CLEARP'

    # NIRISS has many filters in the pupil wheel. WebbPSF does
    # not make a distinction, but Mirage does. Adjust the info
    # to match Mirage's expectations
    if inst.upper() == 'NIRISS' and filt in NIRISS_PUPIL_WHEEL_FILTERS:
        save_filt = copy(filt)
        if pupil == 'CLEARP':
            filt = 'CLEAR'
        else:
            raise ValueError(('Pupil value is something other than '
                              'CLEARP, but the filter being used is '
                              'in the pupil wheel.'))
        pupil = save_filt

    # Same for NIRCam
    if inst.upper() == 'NIRCAM' and filt in NIRCAM_PUPIL_WHEEL_FILTERS:
        save_filt = copy(filt)
        if pupil == 'CLEAR':
            if save_filt[0:2] == 'F4':
                filt = 'F444W'
            elif save_filt[0:2] == 'F3':
                filt = 'F322W2'
            elif save_filt[0:2] == 'F1':
               filt = 'F150W2'
        else:
            raise ValueError(('Pupil value is something other than '
                              'CLEAR, but the filter being used is '
                              'in the pupil wheel.'))
        pupil = save_filt

    elif (inst.upper() == 'NIRCAM') and (pupil.upper() in ['WLP4', 'WLM4', 'WLP12']):
        # If there is a weak lens listed in the library file's pupil value,
        # disentangle it into the correct filter and pupil values
        filt, pupil = fix_weak_lens_filter_pupil_values(filt, pupil)

    opd_file = header['OPD_FILE']
    if default_psf:
        if 'predicted' in opd_file:
            wfe_type = 'predicted'
        elif 'requirements' in opd_file:
            wfe_type = 'requirements'
        realization = header['OPDSLICE']

    # make the check below pass for FGS
    if instrument.lower() == 'fgs':
        pupil = 'N/A'
        pupilname = 'N/A'
        filt = 'N/A'
        filtername = 'N/A'

    match = inst.lower() == instrument.lower() and \
            det.lower() == detector.lower() and \
            filt.lower() == filtername.lower() and \
            pupil.lower() == pupilname.lower()
    if match and not default_psf:
        return full_filename
    elif match and wfe_type == wavefront_error_type.lower() and \
        realization == wavefront_error_group:
            return full_filename
    else:
        return None


def fix_weak_lens_filter_pupil_values(filter_value, pupil_value):
    """By default, WebbPSF modifies weak lens-related filter
    and pupil values that are in the headers of the gridded PSF
    library files. In cases where WLP4 (which is in the filter wheel)
    is used in combination with one of the +/-8 weak lenses (in the
    pupil wheel), WebbPSF will 'add' them together, and report a
    pupil value of WLM4 or WLP12. It will also set the filter value
    to F212N (becuase WLP4 has a F212N-like coating on it). For our
    purposes, we need to know exactly what is in the filter and pupil
    wheels. This function disentangles those values.

    Parameters
    ----------
    filter_value : str
        Value from the FILTER header keyword in the PSF library file

    pupil_value : str
        Value from the PUPIL header keyword in the PSF library file
    """
    save_filt = copy(filter_value)
    file_filt = 'WLP4' ## it's built in filter wheel
    save_pupil = copy(pupil_value)

    # Reverse engineer if it was paired with CLEAR, WLP8 or WLM8
    if save_pupil == 'WLP4':
        file_pupil = 'CLEAR'
    elif save_pupil == 'WLP12':
        file_pupil = 'WLP8'
    elif save_pupil == 'WLM4':
        file_pupil = 'WLM8'
    else:
        raise ValueError(('The pupil paired with WLP4, {}, does not '
                          'match expectations'.format(save_pupil)))
    return file_filt, file_pupil


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

    wavefront_error : str
        Wavefront error. Can be 'predicted' or 'requirements'

    wavefront_error_group : int
        Wavefront error realization group. Must be an integer from 0 - 9.

    library_path : str
        Path pointing to the location of the PSF library

    Returns:
    --------
    library : photutils.griddedPSFModel
        Object containing PSF library

    """
    logger = logging.getLogger('mirage.psf.psf_selection.get_gridded_psf_library')

    # In the default case, we expect the (total PSF signal)/ (oversample factor**2)
    # to be close to 1.0. In certain cases for NIRISS, this expectation is lower by
    # some factor. Here we set the defaul factor to lower expectations to 1.0
    # (i.e. don't lower expectations).
    grid_min_factor = 1.0

    # First, as a way to save time, let's assume a file naming convention
    # and search for the appropriate file that way. If we find a match,
    # confirm the properties of the file via the header. This way we don't
    # need to open and examine every file in the gridded library, which
    # saves at least a handful of seconds.
    if instrument.lower() == 'fgs':
        default_file_pattern = '{}_{}_fovp*_samp*_npsf*_{}_realization{}.fits'.format(instrument.lower(),
                                                                                      detector.lower(),
                                                                                      wavefront_error.lower(),
                                                                                      wavefront_error_group)
    else:
        # NIRISS gridded library names don't follow standard filter/pupil rules.
        # The filenames are all <filter>_<clear>, where <clear> is clear if it
        # is in the filter wheel and clearp if in the pupil wheel.
        if instrument.lower() == 'niriss':
            if filtername.lower() == 'clear':
                filename_filter = pupilname
                filename_pupil = filtername
            elif pupilname.lower() == 'clearp':
                # If CLEARP is used, we need to alert the normalization check
                # below to expect PSFs with a total signal/oversamp**2 lower
                # than otherwise expected by a factor of ~0.84
                filename_filter = filtername
                filename_pupil = pupilname
                grid_min_factor = NIRISS_CLEARP_PSF_THROUGHPUT_REDUCTION
            # filter=clear, pupil=nrm is currently not allowed
            if pupilname.lower() == 'nrm':
                # In NRM mode, the ~0.15 throughput factor associated with the
                # NRM mask is baked into the PSFs from Webbpsf, so we need to
                # adjust the expectations of the normalization check below.
                filename_filter = filtername
                filename_pupil = 'mask_nrm'
                grid_min_factor = NIRISS_NRM_PSF_THROUGHPUT_REDUCTION
        elif instrument.lower() == 'nircam':
            # In the case of a WLP4+WLP8 observation, the gridded PSF library field
            # of view is small enough (and the PSF is large enough) that a significant
            # fraction of the flux is outside the field of view (of the PSF "core" library.
            # In this case, we need to loosen then requirement for the minimum normalized
            # flux in the array.
            if ((filtername.lower() == 'wlp4') & (pupilname.lower() == 'wlp8')):
                grid_min_factor = NIRCAM_WLP12_PSF_THROUGHPUT_REDUCTION
            filename_filter = filtername
            filename_pupil = pupilname if 'GDHS' not in pupilname else 'clear'  # for WFSC team practice purposes we don't produce DHS "PSFs"

        default_file_pattern = '{}_{}_{}_{}_fovp*_samp*_npsf*_{}_realization{}.fits'.format(instrument.lower(),
                                                                                            detector.lower(),
                                                                                            filename_filter.lower(),
                                                                                            filename_pupil.lower(),
                                                                                            wavefront_error.lower(),
                                                                                            wavefront_error_group)
    default_matches = glob(os.path.join(library_path, default_file_pattern))
    library_file = None
    if len(default_matches) == 1:
        library_file = confirm_gridded_properties(default_matches[0], instrument, detector, filtername,
                                                  pupilname, wavefront_error, wavefront_error_group,
                                                  library_path)

    # If the above search found no matching files, or multiple matching
    # files (based only on filename), or if the matching file's gridded
    # PSF model properties don't match what's expected, then resort to
    # opening and examining all files in the library.
    if library_file is None:
        logger.info(("No matching gridded PSF library file found based on filename pattern. Checking "
                 "library files' metadata. This will be slower."))
        library_file = get_library_file(instrument, detector, filtername, pupilname,
                                        wavefront_error, wavefront_error_group, library_path)

    logger.info("PSFs will be generated using: {}".format(os.path.abspath(library_file)))

    lib_head = fits.getheader(library_file)
    itm_sim = lib_head.get('ORIGIN', '') == 'ITM'

    if not itm_sim:
        try:
            library = to_griddedpsfmodel(library_file)
        except OSError:
            logger.error("OSError: Unable to open {}.".format(library_file))
    else:
        # Handle input ITM images
        library = _load_itm_library(library_file)

    # Check that the gridded PSF library is normalized as expected
    check_max = PSF_NORM_MAX * grid_min_factor
    check_min = PSF_NORM_MIN * grid_min_factor
    correct_norm, reason = check_normalization(library, lower_limit=check_min, upper_limit=check_max)
    if correct_norm:
        return library
    else:
        raise ValueError(("Gridded PSF library in {} appears to be improperly normalized."
                          "The total signal in a PSF is {}".format(library_file, reason)))

def get_library_file(instrument, detector, filt, pupil, wfe, wfe_group,
                     library_path, wings=False, segment_id=None):
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

    wings : bool, optional
        Must the library file contain PSF wings or PSF cores? Default is False.

    segment_id : int or None, optional
        If specified, returns a segment PSF library file and denotes the ID
        of the mirror segment

    Returns
    --------
    matches : str
        Name of the PSF library file for the instrument and filter name
    """
    logger = logging.getLogger('mirage.psf.psf_selection.get_library_file')

    psf_files = glob(os.path.join(library_path, '*.fits'))

    # Determine if the PSF path is default or not
    mirage_dir = expand_environment_variable('MIRAGE_DATA')
    gridded_dir = os.path.join(mirage_dir, '{}/gridded_psf_library'.format(instrument.lower()))
    if wings:
        gridded_dir = os.path.join(gridded_dir, 'psf_wings')
    default_psf = library_path == gridded_dir

    # Create a dictionary of header information for all PSF library files
    matches = []

    instrument = instrument.upper()
    detector = detector.upper()
    filt = filt.upper()
    pupil = pupil.upper()
    wfe = wfe.lower()

    # set default
    file_wfe = ''

    # handle the NIRISS NRM case
    if pupil == 'NRM':
        pupil = 'MASK_NRM'

    # Handle the DHS for Coarse Phasing - this is a workaround for webbpsf not
    # implementing this. We're going to load an ITM image in any case in this mode
    # so the PSF is entirely unused, but we need to load something or else MIRAGE errors.
    if pupil == 'GDHS0' or pupil == 'GDHS60':
        pupil = 'CLEAR'

    for filename in psf_files:
        try:
            header = fits.getheader(filename)

            # Determine if it is an ITM file
            itm_sim = header.get('ORIGIN', '') == 'ITM'

            # Compare the header entries to the user input
            file_inst = header['INSTRUME'].upper()
            try:
                file_det = header['DETECTOR'].upper()
            except KeyError:
                file_det = header['DET_NAME'].upper()
            file_filt = header['FILTER'].upper()

            try:
                file_pupil = header['PUPIL'].upper()
            except KeyError:
                # If no pupil mask value is present, then assume the CLEAR is
                # being used
                if file_inst.upper() == 'NIRCAM':
                    file_pupil = 'CLEAR'
                elif file_inst.upper() == 'NIRISS':
                    try:
                        file_pupil = header['PUPIL'].upper()  # can be 'MASK_NRM'
                    except KeyError:
                        file_pupil = 'CLEARP'

            # NIRISS has many filters in the pupil wheel. WebbPSF does
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

            # Same for NIRCam
            if file_inst.upper() == 'NIRCAM' and file_filt in NIRCAM_PUPIL_WHEEL_FILTERS:
                save_filt = copy(file_filt)
                if file_pupil == 'CLEAR':
                    if save_filt[0:2] == 'F4':
                        file_filt = 'F444W'
                    elif save_filt[0:2] == 'F3':
                        file_filt = 'F322W2'
                    elif save_filt[0:2] == 'F1':
                        file_filt = 'F150W2'
                else:
                    raise ValueError(('Pupil value is something other than '
                                      'CLEAR, but the filter being used is '
                                      'in the pupil wheel.'))
                file_pupil = save_filt

            elif (file_inst.upper() == 'NIRCAM') and ((file_pupil.upper() == 'WLP4') |
                                                      (file_pupil.upper() == 'WLM4') |
                                                      (file_pupil.upper() == 'WLP12')):
                # If there is a weak lens listed in the library file's pupil value,
                # disentangle it into the correct filter and pupil values
                file_filt, file_pupil = fix_weak_lens_filter_pupil_values(file_filt, file_pupil)

            if segment_id is None and not itm_sim:
                opd = header['OPD_FILE']
                if 'requirements' in opd:
                    file_wfe = 'requirements'
                elif 'predicted' in opd:
                    file_wfe = 'predicted'

                file_wfe_grp = header['OPDSLICE']

            if segment_id is not None:
                segment_id = int(segment_id)
                file_segment_id = int(header['SEGID'])

            if segment_id is None and itm_sim:
                # If we have an ITM library, then wfe is
                # meaningless, so force it to match
                file_wfe = 'predicted'
                wfe = 'predicted'

            # allow check below to pass for FGS
            if instrument.lower() == 'fgs':
                file_filt = 'N/A'
                filt = 'N/A'
                file_pupil = 'N/A'
                pupil = 'N/A'

            # Evaluate if the file matches the given parameters
            match = (file_inst == instrument
                     and file_det == detector
                     and file_filt == filt
                     and file_pupil == pupil
                     and file_wfe == wfe)

            if not wings and segment_id is None and not itm_sim and default_psf:
                match = match and file_wfe_grp == wfe_group
            if segment_id is not None:
                match = match and file_segment_id == segment_id
            elif not itm_sim and default_psf:
                match = match and file_wfe == wfe

            # If so, add to the list of all matches
            if match:
                matches.append(filename)
        except KeyError as e:
            warnings.warn('While searching for PSF file, error raised when examining {}:\n{}\nContinuing.'.format(os.path.basename(filename), e))
            continue

    # Find files matching the requested inputs
    if len(matches) == 1:
        return matches[0]
    elif len(matches) == 0:
        logger.info('Requested parameters:\ninstrument {}\ndetector {}\nfilt {}\npupil {}\nwfe {}\n'
              'wfe_group {}\nlibrary_path {}\n'.format(instrument, detector, filt, pupil, wfe,
                                                       wfe_group, library_path))
        raise ValueError("No PSF library file found matching requested parameters.")
    elif len(matches) > 1:
        raise ValueError("More than one PSF library file matches requested parameters: {}".format(matches))


def get_psf_wings(instrument, detector, filtername, pupilname, wavefront_error, wavefront_error_group,
                  library_path):
    """Locate the file containing PSF wing image and read them in. The
    idea is that there will only be one file for a given detector/filter/
    pupil/WFE/realization combination. This file will contain a PSF
    sampled at detector resolution and covering some large area in pixels.
    Later, when making the seed image, the appropriate subarray will be
    pulled out of this array for each input source depending on its
    magnitude.

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

    wavefront_error : str
        Wavefront error. Can be 'predicted' or 'requirements'

    wavefront_error_group : int
        Wavefront error realization group. Must be an integer from 0 - 9.

    library_path : str
        Path pointing to the location of the PSF library

    Returns
    -------
    psf_wings : numpy.ndarray
        Array containing the PSF wing data. Note that the outermost row
        and column are not returned, in order to avoid edge effects

    """
    logger = logging.getLogger('mirage.psf.psf_selection.get_psf_wings')

    # First, as a way to save time, let's assume a file naming convention
    # and search for the appropriate file that way. If we find a match,
    # confirm the properties of the file via the header. This way we don't
    # need to open and examine every file in the gridded library, which
    # saves at least a handful of seconds.
    default_file_pattern = '{}_{}_{}_{}_fovp*_samp*_{}_realization{}.fits'.format(instrument.lower(),
                                                                                  detector.lower(),
                                                                                  filtername.lower(),
                                                                                  pupilname.lower(),
                                                                                  wavefront_error.lower(),
                                                                                  wavefront_error_group)
    default_matches = glob(os.path.join(library_path, default_file_pattern))

    wings_file = None
    if len(default_matches) == 1:
        wings_file = confirm_gridded_properties(default_matches[0], instrument, detector, filtername,
                                                pupilname, wavefront_error, wavefront_error_group,
                                                library_path, extname='DET_DIST')

    # If the above search found no matching files, or multiple matching
    # files (based only on filename), or if the matching file's gridded
    # PSF model properties don't match what's expected, then resort to
    # opening and examining all files in the library.
    if wings_file is None:
        logger.info(("No matching gridded PSF wing library file found based on filename pattern. Checking "
                     "library files' metadata. This will be slower."))
        # Find the file containing the PSF wings
        wings_file = get_library_file(instrument, detector, filtername, pupilname,
                                      wavefront_error, wavefront_error_group, library_path, wings=True)

    logger.info("PSF wings will be from: {}".format(os.path.basename(wings_file)))
    with fits.open(wings_file) as hdulist:
        psf_wing = hdulist['DET_DIST'].data
    # Crop the outer row and column in order to remove any potential edge
    # effects leftover from creation
    psf_wing = psf_wing[1:-1, 1:-1]

    for shape in psf_wing.shape:
        if shape % 2 == 0:
            logger.error(("WARNING: PSF wing file contains an even number of rows or columns. "
                   "These must be even."))
            raise ValueError
    return psf_wing


def _load_itm_library(library_file):
    """Load ITM FITS file

    Parameters
    ----------
    library_path : str
        Path pointing to the location of the PSF library

    Returns
    -------
    library : photutils.griddedPSFModel
        Object containing PSF library
    """
    data = fits.getdata(library_file)
    hdr = fits.getheader(library_file)
    if data.shape == (2048, 2048):
        # Normalize the data
        data /= np.sum(data)

        # Add PSF location and oversampling keywords
        hdr['DET_YX0'] = ('(1023, 1023)', "The #0 PSF's (y,x) detector pixel position")
        hdr['OVERSAMP'] = (1, 'Oversampling factor for FFTs in computation')

        # Convert to HDUList and create library
        phdu = fits.PrimaryHDU(data, hdr)
        hdulist = fits.HDUList(phdu)
        library = to_griddedpsfmodel(hdulist)

        return library
    else:
        raise ValueError('Expecting ITM data of size (2048, 2048), not {}'.format(data.shape))
