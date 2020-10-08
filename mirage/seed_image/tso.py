#! /usr/bin/env python

"""
This module contains tools for adding a source with time-varying signal
to a 4-dimensional seed image
"""
import copy
import logging
import numpy as np
import h5py
import os
from scipy.integrate import romb

from mirage.catalogs.hdf5_catalog import open_tso
from mirage.logging import logging_functions
from mirage.utils.constants import LOG_CONFIG_FILENAME, STANDARD_LOGFILE_NAME


classdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../'))
log_config_file = os.path.join(classdir, 'logging', LOG_CONFIG_FILENAME)
logging_functions.create_logger(log_config_file, STANDARD_LOGFILE_NAME)



def add_tso_sources(seed_image, seed_segmentation_map, psf_seeds, segmentation_maps, lightcurves, frametime,
                    total_frames, exposure_total_frames, frames_per_integration, number_of_ints, resets_bet_ints,
                    starting_time=0, starting_frame=0, samples_per_frametime=5):
    """inputs are lists, so that there can be multiple sources, although this
    situation is extremely unlikely

    Parameters
    ----------
    seed_image : numpy.ndarray
        2D seed image in units of counts per second

    seed_segmentation_map : numpy.ndarray
        2D segmentation map associated with ``seed_image``

    psf_seeds : list
        List of 2D seed images, each containing only the signal
        (counts/sec) for a single TSO target.

    segmentation_maps : list
        List of 2D segmentation maps associated with ``psf_seeds``

    lightcurves : list
        List of dictionaries containing lightcurve data. (i.e. output of
        ``read_lightcurve``

    frametime : float
        Exposure time associated with a single frame

    total_frames : int
        Total number of frames in the observation (this includes
        resets between integrations)

    exposure_total_frames : int
        Number of frames in the exposure (excluding resets between
        integrations)

    frames_per_integration : int
        Number of frames per integration of the exposure

    number_of_ints : int
        Number of integrations in the exposure to be simulated

    resets_bet_ints : int
        Number of reset frames between integrations

    samples_per_frametime : int
        Number of points per ``frametime`` to interpolate the lightcurve
        to. For each frame, the lightcurve will be integrated across
        these points using Romberg Integration. This means that
        samples_per_frametime must be 2^k + 1 for some integer k.

    Returns
    -------
    final_seed : numpy.ndarray
        4D array containing the 2D seed image for each frame of each
        integration in the exposure

    final_seed_segmentation_map : numpy.ndarray
        2D array containing the segmentation map
    """
    logger = logging.getLogger('mirage.seed_image.tso.add_tso_sources')

    yd, xd = seed_image.shape
    total_exposure_time = exposure_total_frames * frametime

    # Make sure that samples_per_frametime has a value that will allow
    # Romberg integration
    k_remainder, k = np.modf(np.log2(samples_per_frametime - 1))
    if k_remainder > 0.01:
        raise ValueError("samples_per_frametime must be 2^k + 1 for some integer k.")

    # Put seed image in units of counts per frame
    # seed_image *= frametime
    seed_image_per_frame = seed_image * frametime

    # Frame seed contains only the signal for that particular frame,
    # rather than the cumulative signal since the beginning of the
    # integration
    frame_seed = np.zeros((total_frames, yd, xd))
    final_seed_segmentation_map = np.zeros((yd, xd))
    # Loop over TSO objects
    for source_number, (psf, seg_map, lightcurve) in enumerate(zip(psf_seeds, segmentation_maps, lightcurves)):
        # Check that the provided lightcurve is long enough to cover the
        # exposure time of the observation. If not...extend the lightcurve
        # with values of 1.0 until is is long enough. If the lightcurve
        # starts at some non-zero time, then extend the lightcurve with
        # 1.0's back to time=0.
        lightcurve = check_lightcurve_time(lightcurve, total_exposure_time, frametime)

        # Scale the TSO source's seed image contribution to be for one
        # frametime rather than 1 second
        ft_psf = psf * frametime

        # Interpolate the lightcurve to prepare for integration
        interp_lightcurve = interpolate_lightcurve(copy.deepcopy(lightcurve), samples_per_frametime, frametime)
        dx = frametime / (samples_per_frametime - 1)

        # Integrate the lightcurve for each frame
        logger.info('\nIntegrating lightcurve signal for each frame ')
        for frame_number in np.arange(total_frames) + starting_frame:
            frame_index = frame_number - starting_frame
            #print("Loop 1, frame number and index: ", frame_number, frame_index)
            min_index = frame_number * (samples_per_frametime - 1)
            indexes = np.arange(min_index, min_index + samples_per_frametime)

            # Normalize the integrated signal by the frametime, as that
            # is the integral of a flat line at 1.0 over one frametime
            relative_signal = romb(interp_lightcurve['fluxes'].value[indexes], dx) / frametime
            frame_psf = ft_psf * relative_signal
            tmpy, tmpx = psf.shape
            frame_seed[frame_index, :, :] += frame_psf

        # Add the TSO target to the segmentation map
        final_seed_segmentation_map = update_segmentation_map(seed_segmentation_map, seg_map.segmap)

    # Translate the frame-by-frame seed into the final, cumulative seed
    # image. Rearrange into integrations, resetting the signal for each
    # new integration.
    logger.info('Translate the frame-by-frame transit seed into the final, cumulative seed image.')
    integration_starts = np.arange(number_of_ints) * (frames_per_integration + resets_bet_ints) + starting_frame
    reset_frames = integration_starts[1:] - 1

    if total_frames-len(reset_frames) > frames_per_integration:
        dimension = frames_per_integration
    else:
        dimension = total_frames-len(reset_frames)
    final_seed = np.zeros((number_of_ints, dimension, yd, xd))

    for frame in np.arange(total_frames) + starting_frame:
        int_number = np.where(frame >= integration_starts)[0][-1]
        rel_frame = frame - integration_starts[int_number]

        if frame in integration_starts:
            final_seed[int_number, 0, :, :] = copy.deepcopy(frame_seed[frame-starting_frame, :, :]) + seed_image_per_frame
        elif frame not in reset_frames:
            final_seed[int_number, rel_frame, :, :] = final_seed[int_number, rel_frame-1, :, :] + \
                frame_seed[frame-starting_frame, :, :] + seed_image_per_frame

    return final_seed, final_seed_segmentation_map


def check_lightcurve_time(light_curve, exposure_time, frame_time):
    """Check to be sure the provided lightcurve is long enough to cover
    the supplied total exposure time. If not, lengthen at the beginning
    or end such that it does. Times will only be added to the beginning
    if the first time entry in the lightcurve is > 0. Lightcurves where
    the initial time entry is < 0 will have all times < 0 chopped. This
    will allow the user to simulate lightcurves where the exposure starts
    somewhere in the middle of the lightcurve.

    Parameters
    ----------
    light_curve : dict
        Dictionary of lightcurve. "fluxes" and "times" keys contain
        arrays of those values

    exposure_time : float
        Total exposure time for the full exposure being simulated
        (in seconds)

    frame_time : float
        Exposure time of a single frame of the observation

    Returns
    -------
    light_curve : dict
        Potentially modified with added or removed elements
    """
    logger = logging.getLogger('mirage.seed_image.tso.check_lightcurve_time')

    times = copy.deepcopy(light_curve["times"].value)
    fluxes = copy.deepcopy(light_curve["fluxes"].value)
    time_units = light_curve["times"].unit
    flux_units = light_curve["fluxes"].unit
    adjusted = False

    # Remove elements where time < 0.
    if np.min(times) < 0.:
        positive_times = times >= 0.
        times = times[positive_times]
        fluxes = fluxes[positive_times]
        adjusted = True

    # If the times begin at values significantly > 0,
    # then add entries to bring the start back to time = 0
    if np.min(times) > 0.:
        logger.info(("Lightcurve time values do not start at zero. Prepending an entry with time=0 "
                     "and flux = 1."))
        times = np.insert(times, 0, 0.)
        fluxes = np.insert(fluxes, 0, 1.)
        adjusted = True

    # If the ending time is less than the exposure's total
    # observation time, then add entries with flux=1
    if np.max(times) < exposure_time:
        logger.info(("Lightcurve time values extend only to {} seconds. This is not long enough "
                     "to cover the entire exposure time of {} seconds. Extending to cover the full "
                     "exposure time with flux = 1.".format(np.max(times), exposure_time)))
        times = np.append(times, exposure_time + 5 * frame_time)
        fluxes = np.append(fluxes, 1.)
        adjusted = True

    if adjusted:
        light_curve["times"] = times * time_units
        light_curve["fluxes"] = fluxes * flux_units

    return light_curve


def interpolate_lightcurve(light_curve, samples_per_frame_time, frame_time):
    """Given a lightcurve with arbitrary sampling times, interpolate
    such that it has ``samples_per_frame_time`` samples per frametime.
    In order to use Romberg, we need to integrate over 2^k + 1 samples each time.

    Parameters
    ----------
    light_curve : dict
        Dictionary containing light curve. Times are in 'times' keyword,
        and fluxes are in 'fluxes' keyword. Each should be an array with
        associated unit.

    samples_per_frame_time : int
        Number of times to sample the lightcurve per frame exposure time

    frame_time : float
        Exposure time of one frame, in seconds.

    Returns
    -------
    light_curve : dict
        Modified dictionary containing the resampled lightcurve
    """
    time_units = light_curve['times'].unit
    flux_units = light_curve['fluxes'].unit
    divisor = samples_per_frame_time - 1.
    points = np.arange(light_curve['times'][0].value, light_curve['times'][-1].value, frame_time/divisor)
    light_curve["fluxes"] = np.interp(points, light_curve['times'].value, light_curve['fluxes'].value) * flux_units
    light_curve["times"] = points * time_units
    return light_curve


def read_lightcurve(filename, index_value):
    """Read in hdf5 file containing lightcurve data and return the
    lightcurve in the given index.

    Parameters
    ----------
    filename : str
        Name of hdf5 file containing lightcurves

    index_value : int
        Dataset value within ``filename`` continaing the lightcurve of
        interest

    Returns
    -------
    dataset : dict
        Dictionary containing lightcurve data for the object in dataset
        number ``index_value``.
    """
    logger = logging.getLogger('mirage.seed_image.tso.read_lightcurve')
    file_contents = open_tso(filename)

    try:
        dataset = file_contents[index_value]
    except KeyError as e:
        logger.info(e)
    return dataset


def update_segmentation_map(segmap, object_map):
    """Add the segmentation map for an individual object to that for the
    seed image. Note that multiple objects cannot be recorded in a single
    pixel. In this case, the TSO object of interest will be added to the
    segmentation map after the background objects, so it will overwrite
    any pixels that are affected by background objects.

    Parameters
    ----------
    segmap : numpy.ndarray
        2D segmenatation map of the scene

    object_map : numpy.ndarray
        2D segmentation map of object to add to the scene

    Returns
    -------
    segmap : numpy.ndarray
        Updated 2D segmentation map of the scene
    """
    obj_pix = object_map != 0
    segmap[obj_pix] = object_map[obj_pix]
    return segmap
