#! /usr/bin/env python

"""
This module contains functions that work with JPL Horizons ephemeris files
"""


from astropy.coordinates import SkyCoord
from astropy.table import Table
import calendar
from datetime import datetime
from scipy.interpolate import interp1d


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
    return calendar.timegm(date.timetuple())


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


