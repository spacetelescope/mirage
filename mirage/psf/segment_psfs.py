#! /usr/bin/env python

"""This module contains code to make, locate, and parse the appropriate segment
PSF library files to use for a given simulation. These files are made using
the WebbPSF psf_grid() method, and are turned into photutils GriddedPSFModel
objects using webbpsf.utils.to_griddedpsfmodel.

Author
------

    - Lauren Chambers

Use
---

    This module can be imported and called as such:
    ::
        from mirage.psf.segment_psfs import get_gridded_segment_psf_library_list
        lib = get_gridded_segment_psf_library_list(instrument, detector, filter,
                out_dir, pupilname="CLEAR")
"""
import os
import time

from astropy.io import fits
import numpy as np
import pysiaf
import webbpsf
from webbpsf.gridded_library import CreatePSFLibrary
from webbpsf.utils import to_griddedpsfmodel

from mirage.psf.psf_selection import get_library_file


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
    # Create webbpsf NIRCam instance
    nc = webbpsf.NIRCam()

    # Create dummy CreatePSFLibrary instance to get lists of filter and detectors
    lib = CreatePSFLibrary

    # Define the filter list to loop through
    if isinstance(filters, str):
        filters = [filters]
    elif not isinstance(filters, list):
        raise TypeError('Please define filters as a string or list, not {}'.format(type(filters)))

    # Define the detector list to loop through
    if detectors == 'all':
        detectors = ['NRCA1', 'NRCA2', 'NRCA3', 'NRCA4', 'NRCA5',
                     'NRCB1', 'NRCB2', 'NRCB3', 'NRCB4', 'NRCB5',]
    elif isinstance(detectors, str):
        detectors = [detectors]
    elif not isinstance(detectors, list):
        raise TypeError('Please define detectors as a string or list, not {}'.format(type(detectors)))


    # Create PSF grids for all segments, detectors, and filters
    for i in range(18):
        start_time = time.time()
        i_segment = i + 1
        segname = webbpsf.webbpsf_core.segname(i_segment)
        print('GENERATING SEGMENT {} DATA'.format(segname))
        print('------------------------------')

        det_filt_match = False
        for det in sorted(detectors):
            for filt in list(filters):
                # Make sure the detectors and filters match
                if (det in lib.nrca_short_detectors and filt not in lib.nrca_short_filters) \
                        or (det in lib.nrca_long_detectors and filt not in lib.nrca_long_filters):
                    continue

                det_filt_match = True

                # Define the filter and detector
                nc.filter = filt
                nc.detector = det

                # Restrict the pupil to the current segment
                pupil = webbpsf.webbpsf_core.one_segment_pupil(i_segment)
                ote.amplitude = pupil[0].data
                nc.pupil = ote

                # Generate the PSF grid
                # NOTE: we are choosing a polychromatic simulation here to better represent the
                # complexity of simulating unstacked PSFs. See the WebbPSF website for more details.
                grid = nc.psf_grid(num_psfs=1, save=False, all_detectors=False,
                                   use_detsampled_psf=True, fov_pixels=fov_pixels,
                                   oversample=1, overwrite=overwrite, add_distortion=False,
                                   nlambda=10)

                # Remove and add header keywords about segment
                del grid.meta["grid_xypos"]
                del grid.meta["oversampling"]
                grid.meta['SEGID'] = (i_segment, 'ID of the mirror segment')
                grid.meta['SEGNAME'] = (segname, 'Name of the mirror segment')
                grid.meta['XTILT'] = (round(segment_tilts[i, 0], 2), 'X tilt of the segment in microns')
                grid.meta['YTILT'] = (round(segment_tilts[i, 1], 2), 'Y tilt of the segment in microns')

                # Write out file
                filename = 'nircam_{}_{}_fovp{}_samp1_npsf1_seg{:02d}.fits'.format(det.lower(), filt.lower(),
                                                                                   fov_pixels, i_segment)
                filepath = os.path.join(out_dir, filename)
                primaryhdu = fits.PrimaryHDU(grid.data)
                tuples = [(a, b, c) for (a, (b, c)) in sorted(grid.meta.items())]
                primaryhdu.header.extend(tuples)
                hdu = fits.HDUList(primaryhdu)
                hdu.writeto(filepath, overwrite=overwrite)
                print('Saved gridded library file to {}'.format(filepath))

        if det_filt_match == False:
            raise ValueError('No matching filters and detectors given - all '
                             'filters are longwave but detectors are shortwave, '
                             'or vice versa.')

        print('\nElapsed time:', time.time() - start_time, '\n')


def get_gridded_segment_psf_library_list(instrument, detector, filtername,
                                         library_path, pupilname="CLEAR"):
    """Find the filenames for the appropriate gridded segment PSF libraries and
    read them into griddedPSFModel objects

    Parameters
    ----------
    instrument : str
        Name of instrument the PSFs are from

    detector : str
        Name of the detector within ```instrument```

    filtername : str
        Name of filter used for PSF library creation

    library_path : str
        Path pointing to the location of the PSF library

    pupilname : str, optional
        Name of pupil wheel element used for PSF library creation. Default is "CLEAR".

    Returns:
    --------
    libraries : list of photutils.griddedPSFModel
        List of object containing segment PSF libraries

    """
    library_list = get_segment_library_list(instrument, detector, filtername, library_path, pupil=pupilname)

    print("Segment PSFs will be generated using:")
    for filename in library_list:
        print(os.path.basename(filename))

    libraries = []
    for filename in library_list:
        with fits.open(filename) as hdulist:
        #     hdr = hdulist[0].header
        #     d = hdulist[0].data
        #
        # data = d[0][0]
        # phdu = fits.PrimaryHDU(data, header=hdr)
        # hdulist = fits.HDUList(phdu)

            lib_model = to_griddedpsfmodel(hdulist)
            libraries.append(lib_model)

    return libraries


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
