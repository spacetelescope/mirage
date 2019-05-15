#! /usr/bin/env python

"""
This module contains code for adding a source with time-varying signal
to a seed image


Inputs:
* seed image (counts/second)
* PSF scaled to counts/second for given magnitude
* coordinates within aperture to place the PSF, coordinates in the PSF that correspond to the area on the PSF to use
* hdf5 file of the source

method:
* step along the light curve in time, at a cadence equal to the frame time.
* for each step, integrate under the light curve to find out the relative brightness of the source
* scale the input PSF based on integration results
* scale the seed image to be in counts/frame
* add PSF to the seed image
* add seed image to the signal in the seed image corresponding to the previous frame
* add Poisson noise/cosmic rays
* move on to the next frame

outputs:
* 4d seed image with time varying source (similar to moving targets’ 4d seed image


"""


def tso_seed(seed_image, tso_catalog_file, psf_library, samples_per_frametime=5):
    """

    NOT USED!!


    Will there ever be a case when more than one object is varying in brightness?
    It seems unlikely.
    We can make lists of varying sources I suppose and then pass the lists to the
    function for adding sources. Since adding a source involves stepping along in
    time and creating a 4D seed image, we don't want to do that multiple times
    """
    # Read in the TSO catalog file - the point source catalog reader should work
    #tso_cat = read_tso_catalog(tso_catalog_file)
    #tso_cat = readPointSourceFile(tso_catalog_file)

    # We may be able simply to call getPointSourceList
    # This will take care of sources off the detector, and will check the
    # source index values compared to those in the other catalogs

    # add warning in catalog_Seed_image that moving targets and tso cannot be combined

    # Check source index numbers. We don't want to repeat
    # numbers from the other catalogs
    #check_index_numbers()

    # Read in the lightcurve file
    #lightcurve = hdf5.open_tso(lightcurve_file)
    #lightcurve_list = lightcurve[tso_cat['index']]


    # ?? might be slightly slower, but would avoid repeating code
    tso_seeds = []
    tso_segs = []
    tso_lightcurves = []
    for source in tso_cat:
        seed, seg = make_point_source_image(tso_cat[source])
        tso_seeds.append(seed)
        tso_segs.append(seg)
        lightcurve = read_lightcurve(tso_cat['lightcurve_file'])
        tso_lightcurves.append(lightcurve)
    seed_image = add_tso_sources(seed_image, tso_seeds, tso_segs, lightcurve_list, samples_per_frametime)
    return seed_image

    OR USE THE CASE BELOW:
    # Lists to hold source info
    psf_list = []
    image_coord_list = []
    psf_coord_list = []
    lightcurve_list = []
    for source in tso_cat:
        # Get x, y locations for sources
        getPos(source['ra'], source['dec'])

        # Create the PSF to use
        psf = psf_library.evaluate()
        psf_list.append(psf)

        # Scale PSF to appropraite magnitude
        psf *= scaling_factor

        # Get the coordinates of the aperture/PSF overlap
        i,j,k,l,etc = mirage.stamp_coords()
        image_coord_list.append(coordinates)
        psf_coord_list.append(psf_coordinates)

    # Add the sources to the seed image
    new_seed_image = add_tso_sources(seed_image, psf_list, image_coord_list, psf_coord_list, lightcurve_list)
    return new_seed_image


def add_tso_sources(seed_image, seed_segmentation_map, psf_seeds, sementation_maps, lightcurves, frametime,
                    samples_per_frametime=5):
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

    samples_per_frametime : int
        Number of points per ``frametime`` to interpolate the lightcurve
        to. For each frame, the lightcurve will be integrated across
        these points using Romberg Integration. This means that
        samples_per_frametime must be 2^k + 1 for some integer k.
    """
    # Make sure that samples_per_frametime has a value that will allow
    # Romberg integration
    k_remainder, k = np.modf(np.log2(samples_per_frametime - 1))
    if k_remainder > 0.01:
        raise ValueError("samples_per_frametime must be 2^k + 1 for some integer k.")

    # Put seed image in units of counts per frame
    seed_image *= frametime

    # Loop over TSO objects
    for psf, seg_map, lightcurve in zip(psf_seeds, segmentation_maps, lightcurves):
        psf *= frametime

        # Interpolate the lightcurve to prepare for integration
        frame_timestamps = np.arange(lightcurve['time'].data[0], lightcurve['time'].data[-1], frametime)
        interp_lightcurve = interpolate_lightcurve(lightcurve, samples_per_frametime)

        # Integrate the lightcurve for each frame
        for frame_number in range(total_frames):
            min_index = frame_number * (samples_per_frametime - 1)
            indexes = np.arange(min_index, min_index + samples_per_frametime)
            # Normalize the integrated signal by the frametime, as that
            # is the integral of a flat line at 1.0 over one frametime
            relative_signal = np.romb(interp_lightcurve[indexes]) / frametime
            frame_psf = psf * relative_signal
            if frame_number == 0:
                frame_seed[frame_number, :, :] = seed_image + frame_psf
            else:
                frame_seed[frame_number, :, :] = frame_seed[frame_number-1, :, :] + seed_image + frame_psf

        # Add the TSO target to the segmentation map
        seed_segmentation_map = update_segmentation_map(seed_segmentation_map, seg_map)
    return frame_seed, seed_segmentation_map


def interpolate_lightcurve(lightcurve, samples_per_frame_time):
    """Given a lightcurve with arbitrary sampling times, interpolate
    such that it has 3 (or 5?) samples per frametime. Start/mid/end, where
    start is the same as end from the previous frametime, and end is
    the same as start from the subsequent frametime. In order to use
    Romberg, we need to integrate over 2^k + 1 samples each time.
    """
    divisor = samples_per_frame_time - 1.
    points = np.arange(0, end_time, frametime/divisor)
    lightcurve["fluxes"] = np.interpolate(points)
    lightcurve["times"] = points
    return lightcurve


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
    file_contents = hdf5.open_tso(filename)
    try:
        dataset = file_contents[index_value]
    except KeyError as e:
        print(e)
    return dataset


def update_segmentation_map(segmap, object_map):
    """Add the segmentation map for an individual object to that for the
    seed image. Note that multiple objects cannot be recorded in a single
    pixel. In this case, the TSO object of interest will be added to the
    segmentation map after the background objects, so it will overwrite
    any pixels that are affected by background objects.

    Parameters
    ----------
    segmap :

    object_map :

    Returns
    -------
    segmap :
    """
    obj_pix = object_map != 0
    segmap[obj_pix] = object_map[obj_pix]
    return segmap
