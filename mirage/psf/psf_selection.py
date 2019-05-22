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
        from mirage.psf import psf_selction
        library = psf_selection.get_gridded_psf_library('nircam', 'nrcb1',
                                                        'f200w', 'clear',
                                                        'predicted', 0,
                                                        '/path/to/library/')
"""


from copy import copy
from glob import glob
import os

from astropy.io import fits
import numpy as np
import pysiaf
from webbpsf.utils import to_griddedpsfmodel

from mirage.utils.constants import NIRISS_PUPIL_WHEEL_FILTERS


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

    opd_file = header['OPD_FILE']
    if 'predicted' in opd_file:
        wfe_type = 'predicted'
    elif 'requirements' in opd_file:
        wfe_type = 'requirements'
    realization = header['OPDSLICE']

    if inst.lower() == instrument.lower() and det.lower() == detector.lower() and \
       filt.lower() == filtername.lower() and pupil.lower() == pupilname.lower() and \
       wfe_type == wavefront_error_type.lower() and realization == wavefront_error_group:
        return full_filename
    else:
        return None


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
    # First, as a way to save time, let's assume a file naming convention
    # and search for the appropriate file that way. If we find a match,
    # confirm the properties of the file via the header. This way we don't
    # need to open and examine every file in the gridded library, which
    # saves at least a handful of seconds.
    default_file_pattern = '{}_{}_{}_{}_fovp*_samp*_npsf*_{}_realization{}.fits'.format(instrument.lower(),
                                                                                        detector.lower(),
                                                                                        filtername.lower(),
                                                                                        pupilname.lower(),
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
        library_file = get_library_file(instrument, detector, filtername, pupilname,
                                        wavefront_error, wavefront_error_group, library_path)

    print("PSFs will be generated using: {}".format(os.path.basename(library_file)))

    try:
        library = to_griddedpsfmodel(library_file)
    except OSError:
        print("OSError: Unable to open {}.".format(library_file))
    return library


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
        ?

    segment_id : int or None, optional
        If specified, returns a segment PSF library file and denotes the ID
        of the mirror segment

    Returns
    --------
    matches : str
        Name of the PSF library file for the instrument and filter name
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
        # Open the file header
        header = fits.getheader(filename)

        # Compare the header entries to the user input
        file_inst = header['INSTRUME'].upper()

        if segment_id is None:
            try:
                file_det = header['DETECTOR'].upper()
            except KeyError:
                file_det = header['DET_NAME'].upper()
            det_match = file_det == detector
        else:
            det_keys = [k for k in header.keys() if 'DETNAME' in k]
            dets = [header[k] for k in det_keys]
            det_match = detector in dets
        file_filt = header['FILTER'].upper()

        try:
            file_pupil = header['PUPIL_MASK'].upper()
        except KeyError:
            # If no pupil mask value is present, then assume the CLEAR is
            # being used
            if file_inst.upper() == 'NIRCAM':
                file_pupil = 'CLEAR'
            elif file_inst.upper() == 'NIRISS':
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

        if segment_id is None:
            opd = header['OPD_FILE']
            if 'requirements' in opd:
                file_wfe = 'requirements'
            elif 'predicted' in opd:
                file_wfe = 'predicted'

            file_wfe_grp = header['OPDSLICE']

        if segment_id is not None:
            segment_id = int(segment_id)
            file_segment_id = int(header['SEGID'])

        # Evaluate if the file matches the given parameters
        match = (file_inst == instrument
                 and det_match
                 and file_filt == filt
                 and file_pupil == pupil)
        if not wings and segment_id is None:
            match = match and file_wfe_grp == wfe_group
        if segment_id is not None:
            match = match and file_segment_id == segment_id
        else:
            match = match and file_wfe == wfe

        # If so, add to the list of all matches
        if match:
            matches.append(filename)

    # Find files matching the requested inputs
    if len(matches) == 1:
        return matches[0]
    elif len(matches) == 0:
        raise ValueError("No PSF library file found matching requested parameters.")
    elif len(matches) > 1:
        raise ValueError("More than one PSF library file matches requested parameters: {}".format(matches))


def get_segment_library_list(instrument, detector, filt,
                             library_path, pupil='CLEAR'):
    """Given an instrument and filter name along with the path of
    the PSF library, find the appropriate 18 segment PSF library files.

    Parameters
    -----------
    instrument : str
        Name of instrument the PSFs are from

    detector : str
        Name of the detector within ```instrument```

    filt : str
        Name of filter used for PSF library creation

    library_path : str
        Path pointing to the location of the PSF library

    pupil : str, optional
        Name of pupil wheel element used for PSF library creation. Default is
        'CLEAR'.

    segment_id : int or None, optional
        If specified, returns a segment PSF library file and denotes the ID
        of the mirror segment

    Returns
    --------
    library_list : list
        List of the names of the segment PSF library files for the instrument
        and filter name
    """
    library_list = []
    for seg_id in np.arange(1, 19):
         segment_file = get_library_file(
             instrument, detector, filt, pupil, '', 0, library_path,
             segment_id=seg_id
         )
         library_list.append(segment_file)

    return library_list


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
        # Find the file containing the PSF wings
        wings_file = get_library_file(instrument, detector, filtername, pupilname,
                                      wavefront_error, wavefront_error_group, library_path, wings=True)

    print("PSF wings will be from: {}".format(os.path.basename(wings_file)))
    with fits.open(wings_file) as hdulist:
        psf_wing = hdulist['DET_DIST'].data
    # Crop the outer row and column in order to remove any potential edge
    # effects leftover from creation
    psf_wing = psf_wing[1:-1, 1:-1]

    for shape in psf_wing.shape:
        if shape % 2 == 0:
            print(("WARNING: PSF wing file contains an even number of rows or columns. "
                   "These must be even."))
            raise ValueError
    return psf_wing


def get_segment_offset(segment_number, detector, library_list):
    """Convert vectors coordinates in the local segment control
    coordinate system to NIRCam detector X and Y coordinates,
    at least proportionally, in order to calculate the location
    of the segment PSFs on the given detector.

    Parameters
    ----------
    segment : int
        Segment ID, i.e 3
    detector : str
        Name of NIRCam detector
    library_list : list
        List of the names of the segment PSF library files

    Returns
    -------
    x_displacement
        The shift of the segment PSF in NIRCam SW x pixels
    y_displacement
        The shift of the segment PSF in NIRCam SW y pixels
    """

    # Verify that the segment number in the header matches the index
    seg_index = int(segment_number) - 1
    header = fits.getheader(library_list[seg_index])

    assert int(header['SEGID']) == int(segment_number), \
        "Uh-oh. The segment ID of the library does not match the requested " \
        "segment. The library_list was not assembled correctly."
    xtilt = header['XTILT']
    ytilt = header['YTILT']
    segment = header['SEGNAME'][:2]

    # These conversion factors were empirically calculated by measuring the
    # relation between tilt and the pixel displacement
    tilt_to_pixel_slope = 13.4
    tilt_to_pixel_intercept = 0

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
    x_displacement = -(tilt_onto_x * tilt_to_pixel_slope) + tilt_to_pixel_intercept  # pixels
    y_displacement = -(tilt_onto_y * tilt_to_pixel_slope) + tilt_to_pixel_intercept  # pixels

    # Get the appropriate pixel scale from pysiaf
    siaf = pysiaf.Siaf('nircam')
    aperture = siaf['NRC{}_FULL'.format(detector[-2:].upper())]
    nircam_x_pixel_scale = aperture.XSciScale  # arcsec/pixel
    nircam_y_pixel_scale = aperture.YSciScale  # arcsec/pixel

    # Convert the pixel displacement into angle
    x_arcsec = x_displacement * nircam_x_pixel_scale  # arcsec
    y_arcsec = y_displacement * nircam_y_pixel_scale  # arcsec

    return x_arcsec, y_arcsec
