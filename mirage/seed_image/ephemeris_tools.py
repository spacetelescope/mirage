#! /usr/bin/env python

"""
This module contains functions that work with JPL Horizons ephemeris files
"""


from astropy.coordinates import SkyCoord
from astropy.table import Table
import calendar
from datetime import datetime, timezone
import numpy as np
from scipy.interpolate import interp1d


def calculate_nested_positions(x_or_ra_frames, y_or_dec_frames, times_list, spatial_frequency,
                               ra_ephemeris=None, dec_ephemeris=None, position_units='angular'):
    """Given the (x, y) or (RA, Dec) locations of a source at each frame within an integration,
    use an ephemeris or interpolation to calculate the source location at sub-frametime scales,
    in order to support moving target additions.

    Parameters
    ----------
    x_or_ra_frames : list
        List of x-pixel or RA values corresponding to the source location at each frame of the integration

    y_or_dec_frames : list
        List of y-pixel or Dec values corresponding to the source location at each frame of the integration

    times_list : list
        List of times corresponding to each frame of the integration. When using an ephemeris, the units of
        the time values are assumed to be calendar.timegm, so that they can be used by the ephemeris functions.
        When working without an ephemeris, any unit that can be used by np.interp should be ok.

    spatial_frequency : float
        Number of pixels or arcseconds between sub-frame images. Historically, this has been 0.3 pixels. The
        output subframe positions and times will be calculated at this spatial frequency.

    ra_ephemeris : scipy.interpolate.interp1d
        Interpolation function describing the RA portion of the ephemeris of the source.i.e. output from
        create_interpol_function().

    dec_ephemeris : scipy.interpolate.interp1d
        Interpolation function describing the Dec portion of the ephemeris of the source. i.e. output from
        create_interpol_function().

    position_units : str
        Describes the units of ``x_or_ra_frames``, ``y_or_dec_frames``, and ``spatial_frequency``. Can be
        'angular' if positions are in units of RA, Dec (degrees) and ``spatial_frequency`` is in arcsec, or
        'pixels' if positions and ``spatital_frequency`` are in units of detector pixels.

    Returns
    -------
    ra_frames_nested : list
        Nested list of floats giving the pixel x, or RA, positions of the source. Each element of the list
        is a list describing the RA position of the source betwen the current frame and the next frame, with
        a spatial frequency of ``spatial_frequency``.

    dec_frames_nested : list
        Nested list of floats giving the pixel y, or Dec, positions of the source. Each element of the list
        is a list describing the Dec position of the source betwen the current frame and the next frame, with
        a spatial frequency of ``spatial_frequency``.

    subframe_times_nested : list
        Nested list of floats giving the times associated with all elements of ``ra_frames_nested`` and
        ``dec_frames_nested``
    """
    subframe_times_nested = []
    ra_frames_nested = []
    dec_frames_nested = []

    # Create an interpolation function using the frame locations for the case where ephemeris functions
    # are not provided
    if ra_ephemeris is None:
        ra_interpol = interp1d(times_list, x_or_ra_frames)
        dec_interpol = interp1d(times_list, y_or_dec_frames)

    # Here we are looping over all frames in the entire exposure.
    for i, (ra_frame, dec_frame, frame_time) in enumerate(zip(x_or_ra_frames[1:], y_or_dec_frames[1:], times_list[1:])):
        # Note that within this loop, i starts at 0, but is pointing to the first (not zeroth)
        # elements of ra_frames, dec_frames, all_times.

        if position_units == 'angular':
            # If inputs are in RA, Dec (which are in degrees), translate to arcseconds
            delta_ra = (ra_frame - x_or_ra_frames[i]) * 3600.
            delta_dec = (dec_frame - y_or_dec_frames[i]) * 3600.
        elif position_units == 'pixels':
            # If input positions are in units of pixels, there's no need to translate
            pass

        # Calculate distance between the location in the current frame and the preceding frame
        delta_pos = np.sqrt(delta_ra**2 + delta_dec**2)

        # How many points do we need to follow the given spatial scale?
        num_sub_frame_points = int(np.ceil(delta_pos / spatial_frequency))
        if num_sub_frame_points == 1:
            # If the source moves less than the spatial frequency limit, then we'll only need to evaluate
            # the PSF once, using the end time of the frame
            sub_frame_times = [frame_time]
        elif num_sub_frame_points == 2:
            # If the source moves just over the spatial frequency limit, then we'll need to evaluate the
            # PSF twice. Use the start and end time of the frame
            sub_frame_times = [times_list[i], frame_time]
        elif num_sub_frame_points > 2:
            # Here we need to evaluate the PSF three or more times. Use the frame start and end time,
            # and then spread the remaining times evenly throughout the frame time
            frame_start_dt = datetime.fromtimestamp(times_list[i], tz=timezone.utc)
            frame_end_dt = datetime.fromtimestamp(frame_time, tz=timezone.utc)
            subframe_delta_time = (frame_end_dt - frame_start_dt) / (num_sub_frame_points - 1)
            sub_frame_times = [to_timestamp(frame_start_dt + subframe_delta_time * n) for n in range(num_sub_frame_points)]
        elif num_sub_frame_points == 0:
            # In this case the source didn't move at all. Not even a fraction of a pixel.
            sub_frame_times = [frame_time]

        if ra_ephemeris is not None and position_units == 'angular':
            # If ephemeris functions are provided, and positions are in RA, Dec, then calculate the sub
            # frame locations using those
            subframe_ra = ra_ephemeris(sub_frame_times)
            subframe_dec = dec_ephemeris(sub_frame_times)
        else:
            # If no ephemeris functions are provided, then use interpolation to get the sub frame locations
            subframe_ra = ra_interpol(sub_frame_times)
            subframe_dec = dec_interpol(sub_frame_times)

        ra_frames_nested.append(subframe_ra)
        dec_frames_nested.append(subframe_dec)
        subframe_times_nested.append(sub_frame_times)
        return ra_frames_nested, dec_frames_nested, subframe_times_nested


def create_interpol_function(ephemeris):
    """Given an ephemeris, create an interpolation function that can be
    used to get the RA Dec at an arbitrary time

    Parameters
    ----------
    ephemeris : astropy.table.Table

    Returns
    -------
    eph_interp : scipy.interpolate.interp1d
    """
    # In order to create an interpolation function, we need to translate
    # the datetime objects into calendar timestamps
    time = [to_timestamp(entry) for entry in ephemeris['Time']]
    ra_interp = interp1d(time, ephemeris['RA'].data,
                         bounds_error=True, kind='quadratic')
    dec_interp = interp1d(time, ephemeris['Dec'].data,
                          bounds_error=True, kind='quadratic')
    return ra_interp, dec_interp


def get_ephemeris(method):
    """Wrapper function to simplify the creation of an ephemeris

    Parameters
    ----------
    method : str
        Method to use to create the ephemeris. Can be one of two
        options:
        1) Name of an ascii file containing an ephemeris.
        2) 'create' - Horizons is queried in order to produce an ephemeris

    Returns
    -------
    ephemeris : tup
        Tuple of interpolation functions for (RA, Dec). Interpolation
        functions are for RA (or Dec) in degrees as a function of
        calendar timestamp
    """
    if method.lower() != 'create':
        ephem = read_ephemeris_file(method)
    else:
        raise NotImplementedError('Horizons query not yet working')
        start_date = datetime.datetime.strptime(starting_date, '%Y-%m-%d')
        earlier = start_date - datetime.timedelta(days=1)
        later = start_date + datetime.timedelta(days=1)
        step_size = 0.1  # days
        ephem = query_horizons(target_name, earlier, later, step_size)

    ephemeris = create_interpol_function(ephem)
    return ephemeris


def to_timestamp(date):
    """Convert a datetime object into a calendar timestamp object

    Parameters
    ----------
    date : datetime.datetime
        Datetime object e.g. datetime.datetime(2020, 10, 31, 0, 0)

    Returns
    -------
    cal : calendar.timegm
        Calendar timestamp corresponding to the input datetime
    """
    return  date.timestamp()


def read_ephemeris_file(filename):
    """Read in an ephemeris file from Horizons

    Parameters
    ----------
    filename : str
        Name of ascii file to be read in

    Returns
    -------
    ephemeris : astropy.table.Table
        Table of ephemeris data containing SkyCoord positions and
        datetime times.
    """
    with open(filename) as fobj:
        lines = fobj.readlines()

    use_line = False
    ra = []
    dec = []
    time = []
    for i, line in enumerate(lines):
        newline = " ".join(line.split())
        if 'Date__(UT)__HR:MN' in line:
            use_line = True
            start_line = i
        if use_line:
            try:
                date_val, time_val, ra_h, ra_m, ra_s, dec_d, dec_m, dec_s, *others = newline.split(' ')
                ra_str = '{}h{}m{}s'.format(ra_h, ra_m, ra_s)
                dec_str = '{}d{}m{}s'.format(dec_d, dec_m, dec_s)
                location = SkyCoord(ra_str, dec_str, frame='icrs')
                ra.append(location.ra.value)
                dec.append(location.dec.value)
                dt = datetime.strptime("{} {}".format(date_val, time_val), "%Y-%b-%d %H:%M")
                time.append(dt)
            except:
                pass

            if (('*****' in line) and (i > (start_line+2))):
                use_line = False

    ephemeris = Table()
    ephemeris['Time'] = time
    ephemeris['RA'] = ra
    ephemeris['Dec'] = dec
    return ephemeris


def query_horizons(object_name, start_date, stop_date, step_size):
    """Use astroquery to query Horizons and produce an ephemeris

    Paramteres
    ----------
    start_date : str
        e.g. '2020-10-10'

    stop_date : str
        e.g. '2020-11-27'

    step_size : str
        e.g. '10d'

    Returns
    -------
    eph :
        XXXXXX
    """
    from astroquery.jplhorizons import Horizons

    obj = Horizons(id=object_name, location='568',
                   epochs={'start':start_date, 'stop':stop_date,
                           'step':step_size})
    #eph = obj.ephemerides()
    raise NotImplementedError('Horizons query from within Mirage not yet implemented.')
    return eph


