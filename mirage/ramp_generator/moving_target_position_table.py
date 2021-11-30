#! /usr/bin/env python

"""Functions relating to the addition of the moving_target_positions table to the uncal file
in cases where we are dealing with a moving target.

Here is an example of a moving target position table:

h = fits.open('example_OTB_file_uncal.fits')
h['MOVING_TARGET_POSITION'].data
Out[5]:
FITS_rec([(59003.01650053, 1024, 1023, 1024., 1023., 0., 0., -4.62085139e+08, 3.3316837e+08, -6.86525173e+08, -1152263.83524321, -1152263.83524321, -1152263.83524321, -1152263.83524321, -1041098.21566225, -578924.63900379, 8.91446768e+08, 8.92099015e+08, -99.),
          (59003.017184  , 1024, 1023, 1024., 1023., 0., 0., -4.62085139e+08, 3.3316837e+08, -6.86525173e+08, -1152263.37236178, -1152263.37236178, -1152263.37236178, -1152263.37236178, -1041098.54461875, -578924.54186454, 8.91446769e+08, 8.92099015e+08, -99.)],
         dtype=(numpy.record, [('time', '>f8'), ('moving_target_x', '>i4'), ('moving_target_y', '>i4'), ('ref_pixel_RA', '>f8'), ('ref_pixel_Dec', '>f8'), ('moving_target_RA', '>f8'), ('moving_target_Dec', '>f8'), ('mt_x_helio', '>f8'), ('mt_y_helio', '>f8'), ('mt_z_helio', '>f8'), ('jwst_x_helio', '>f8'), ('jwst_y_helio', '>f8'), ('jwst_z_helio', '>f8'), ('mt_x_jwst', '>f8'), ('mt_y_jwst', '>f8'), ('mt_z_jwst', '>f8'), ('mt_jwst_distance', '>f8'), ('mt_sun_distance', '>f8'), ('phase_angle', '>f8')]))

mov = h['MOVING_TARGET_POSITION'].data

In [7]: mov.columns
Out[7]:
ColDefs(
    name = 'time'; format = 'D'; unit = 'seconds'
    name = 'moving_target_x'; format = 'J'; unit = 'pixel'
    name = 'moving_target_y'; format = 'J'; unit = 'pixel'
    name = 'ref_pixel_RA'; format = 'D'; unit = 'deg'
    name = 'ref_pixel_Dec'; format = 'D'; unit = 'deg'
    name = 'moving_target_RA'; format = 'D'; unit = 'deg'
    name = 'moving_target_Dec'; format = 'D'; unit = 'deg'
    name = 'mt_x_helio'; format = 'D'; unit = 'km'
    name = 'mt_y_helio'; format = 'D'; unit = 'km'
    name = 'mt_z_helio'; format = 'D'; unit = 'km'
    name = 'jwst_x_helio'; format = 'D'; unit = 'km'
    name = 'jwst_y_helio'; format = 'D'; unit = 'km'
    name = 'jwst_z_helio'; format = 'D'; unit = 'km'
    name = 'mt_x_jwst'; format = 'D'; unit = 'km'
    name = 'mt_y_jwst'; format = 'D'; unit = 'km'
    name = 'mt_z_jwst'; format = 'D'; unit = 'km'
    name = 'mt_jwst_distance'; format = 'D'; unit = 'km'
    name = 'mt_sun_distance'; format = 'D'; unit = 'km'
    name = 'phase_angle'; format = 'D'; unit = 'deg'
)

"""
from astropy.time import Time, TimeDelta
from datetime import datetime
import numpy as np

from mirage.seed_image.ephemeris_tools import to_timestamp


def create_mt_pos_entry(time_string, movtarg_x, movtarg_y, refpix_ra, refpix_dec,
                        movtarg_ra, movtarg_dec, mt_x_helio, mt_y_helio, mt_z_helio, jwst_x_helio,
                        jwst_y_helio, jwst_z_helio, mt_x_jwst, mt_y_jwst, mt_z_jwst, mt_jwst_distance,
                        mt_sun_distance, phase_angle):
        """Create a single entry for the MOVING_TARGET_POSITION table. The assumption here
        is that inputs will be coming from the GROUP table, which is why the time-based
        columns are so screwy.

        Parameters
        ----------
        time_string : str
            Time of the entry in isot format
            (e.g. '2022-12-18T00:11:41.996')

        movtarg_x : float
            The X location of the moving target in the aperture

        movtarg_y : float
            The Y location of the moving target in the aperture

        refpix_ra : float
            RA of the reference location of the aperture

        refpix_dec : float
            Dec of the reference location of the aperture

        movtarg_ra : float
            RA of the moving target

        movtarg_dec : float
            Dec of the moving target

        mt_x_helio : float
            X-coord of moving target heliocentric distance

        mt_y_helio : float
            Y-coord of moving target heliocentric distance

        mt_z_helio : float
            Z-coord of moving target heliocentric distance

        jwst_x_helio : float
            X-coord of JWST heliocentric distance

        jwst_y_helio : float
            Y-coord of JWST heliocentric distance

        jwst_z_helio : float
            Z-coord of JWST heliocentric distance

        mt_x_jwst : float
            X-coord of moving target-JWST distance

        mt_y_jwst : float
            Y-coord of moving target-JWST distance

        mt_z_jwst : float
            Z-coord of moving target-JWST distance

        mt_jwst_distance : float
            Moving target-JWST distance

        mt_sun_distance : float
            Moving target-Sun distance

        phase_angle : float
            Phase angle

        Returns
        -------
        group : nump.ndarray
            Input values organized into format needed for group entry in
            JWST formatted file
        """
        mjd_day = Time(time_string).mjd

        position = np.ndarray(
            (1, ),
            dtype=[
                ('time', '>f8'),
                ('mt_detector_x', '>i4'),
                ('mt_detector_y', '>i4'),
                ('ref_pixel_RA', '>f8'),
                ('ref_pixel_Dec', '>f8'),
                ('mt_apparent_RA', '>f8'),
                ('mt_apparent_Dec', '>f8'),
                ('mt_apparent_x_helio', '>f8'),
                ('mt_apparent_y_helio', '>f8'),
                ('mt_apparent_z_helio', '>f8'),
                ('jwst_x_helio', '>f8'),
                ('jwst_y_helio', '>f8'),
                ('jwst_z_helio', '>f8'),
                ('mt_apparent_x_jwst', '>f8'),
                ('mt_apparent_y_jwst', '>f8'),
                ('mt_apparent_z_jwst', '>f8'),
                ('mt_apparent_jwst_distance', '>f8'),
                ('mt_apparent_sun_distance', '>f8'),
                ('apparent_phase_angle', '>f8'),
                ('apparent_light_travel_time', '>f8'),
                ('mt_true_x_helio', '>f8'),
                ('mt_true_y_helio', '>f8'),
                ('mt_true_z_helio', '>f8'),
                ('mt_true_x_jwst', '>f8'),
                ('mt_true_y_jwst', '>f8'),
                ('mt_true_z_jwst', '>f8')
            ]
        )
        position[0]['time'] = mjd_day
        position[0]['mt_detector_x'] = movtarg_x
        position[0]['mt_detector_y'] = movtarg_y
        position[0]['ref_pixel_RA'] = refpix_ra
        position[0]['ref_pixel_Dec'] = refpix_dec
        position[0]['mt_apparent_RA'] = movtarg_ra
        position[0]['mt_apparent_Dec'] = movtarg_dec
        position[0]['mt_apparent_x_helio'] = mt_x_helio
        position[0]['mt_apparent_y_helio'] = mt_y_helio
        position[0]['mt_apparent_z_helio'] = mt_z_helio
        position[0]['jwst_x_helio'] = jwst_x_helio
        position[0]['jwst_y_helio'] = jwst_y_helio
        position[0]['jwst_z_helio'] = jwst_z_helio
        position[0]['mt_apparent_x_jwst'] = mt_x_jwst
        position[0]['mt_apparent_y_jwst'] = mt_y_jwst
        position[0]['mt_apparent_z_jwst'] = mt_z_jwst
        position[0]['mt_apparent_jwst_distance'] = mt_jwst_distance
        position[0]['mt_apparent_sun_distance'] = mt_sun_distance
        position[0]['apparent_phase_angle'] = phase_angle
        position[0]['apparent_light_travel_time'] = np.sqrt(mt_x_jwst**2 + mt_y_jwst**2 + mt_z_jwst**2) / 3.e8
        position[0]['mt_true_x_helio'] = mt_x_helio
        position[0]['mt_true_y_helio'] = mt_y_helio
        position[0]['mt_true_z_helio'] = mt_z_helio
        position[0]['mt_true_x_jwst'] = mt_x_jwst
        position[0]['mt_true_y_jwst'] = mt_y_jwst
        position[0]['mt_true_z_jwst'] = mt_z_jwst

        return position


def populate_moving_target_table(grouptable, ephem_interp_func, movtarg_x, movtarg_y, refpix_ra, refpix_dec):
    """Given an instance of the Group table from a datamodel, along with
    an interpolation function for the target's RA and Dec, construct a
    basic moving_target_position table that can be added to the
    moving_target_position attribute of the datamodel. In this initial
    case, we use dummy data to populate the columns relating to distances
    between the target and JWST and the target and the sun.

    Parameters
    ----------
    grouptable : numpy.array
        Group table array, in the format used as input for the Group
        attribute. Each element of grouptable is a tuple containing
        one entry for each column of the table.

    ephem_interp_func : tup
        2-tuple of scipy.interpolate.interp1d objects. These are the
        interpolation functions that will give the RA, Dec of the target
        at an input time

    movtarg_x : float
        The x-coordinate of the location of the moving target (pixels).

    movtarg_y : float
        The y-coordinate of the location of the moving target (pixels).

    refpix_ra : float
        RA value corresponding to the aperture's reference location (deg)

    refpix_dec : float
        Dec value corresponding to the aperture's reference location (deg)

    Returns
    -------
    mt_position : numpy.array
        Array containing the moving_target_position table, formatted such
        that the moving_target_position datamodel attribute can be directly
        set as this array.
    """
    ra_func, dec_func = ephem_interp_func

    # Create the table with a first row populated by garbage
    mt_position = create_mt_pos_entry('2000-01-01', 1024, 1024, 0., 0., 0., 0.,
                                      999., 999., 999., 888., 888., 888.,
                                      777., 777., 777., 1234., 1234., 90.)

    # Dummy data to use for the time being
    mt_y_helio = 5.30e8
    mt_x_helio = 5.30e8
    mt_z_helio = 1.23e5
    jwst_x_helio = 1.56e8
    jwst_y_helio = 1.6e4
    jwst_z_helio = 1.23e5
    mt_x_jwst = 3.74e8
    mt_y_jwst = 5.3e8
    mt_z_jwst = 1.4e2
    mt_jwst_distance = 6.49e8
    mt_sun_distance = 7.5e8
    phase_angle = 9.77

    for line in grouptable:
        line_day = Time(line[0][5])
        line_day_datetime = obstime_to_datetime(line_day.mjd)
        line_day_calstamp = to_timestamp(line_day_datetime)
        interp_ra = ra_func(line_day_calstamp)
        interp_dec = dec_func(line_day_calstamp)

        entry = create_mt_pos_entry(line[0][5], movtarg_x, movtarg_y, refpix_ra, refpix_dec,
                        interp_ra, interp_dec, mt_x_helio, mt_y_helio, mt_z_helio, jwst_x_helio,
                        jwst_y_helio, jwst_z_helio, mt_x_jwst, mt_y_jwst, mt_z_jwst, mt_jwst_distance,
                        mt_sun_distance, phase_angle)

        mt_position = np.vstack([mt_position, entry])

    # Now remove the top garbage row from the table
    mt_position = mt_position[1:]
    return mt_position


def obstime_to_datetime(obstime):
    """Convert the observation time (from the datamodel instance) into a
    datetime (which can be used by the ephemeris interpolation function)

    Parameters
    ----------
    obstime : float
        Observation time in MJD. If populating the MT_RA, MT_DEC keywords,
        this will be the mid-time of the exposure.

    Returns
    -------
    time_datetime : datetime.datetime
        Datetime associated with ```obstime```
    """
    mid_time_astropy = Time(obstime, format='mjd')
    isot = mid_time_astropy.isot
    date_val, time_val = isot.split('T')
    time_val_parts = time_val.split(':')
    fullsec = float(time_val.split(':')[2])
    microsec = '{}'.format(int((fullsec - int(fullsec)) * 1e6))
    new_time_val = '{}:{}:{}:{}'.format(time_val_parts[0], time_val_parts[1], int(fullsec), microsec)
    time_datetime = datetime.strptime("{} {}".format(date_val, new_time_val), "%Y-%m-%d %H:%M:%S:%f")
    return time_datetime


