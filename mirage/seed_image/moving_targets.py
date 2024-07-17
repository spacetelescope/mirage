#! /usr/bin/env python

'''
Class that creates an integration containing multiple frames and
shows a source that is moving relative to the detector.

Arguments:
----------
stamp -- 2D stamp image containing target
xframes -- list of x-coordinate pixel position of target
           in each frame
yframes -- list of y-coordinate pixel position of target
           in each frame
timeframes -- nested list of the time associated with each subframe
frametime -- exposure time in seconds corresponding to one
             detector readout (varies with subarray size)
outx -- x-dimension size of the output aperture (2048 for
        full-frame)
outy -- y-dimension size of the output aperture (2048 for
        full-frame)

Returns:
--------
3D array containing the signal of the source in each frame
of the integration

Author:
-------
Bryan Hilbert
'''

import logging
import os
import sys

import numpy as np
from astropy.io import fits

from ..logging import logging_functions
from ..utils.constants import LOG_CONFIG_FILENAME, STANDARD_LOGFILE_NAME
from ..utils.utils import flatten_nested_list


classdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../'))
log_config_file = os.path.join(classdir, 'logging', LOG_CONFIG_FILENAME)
logging_functions.create_logger(log_config_file, STANDARD_LOGFILE_NAME)


class MovingTarget():

    def __init__(self):
        self.logger = logging.getLogger(__name__)
        self.verbose = False
        self.subsampx = 3
        self.subsampy = 3

    def create(self, stamps_nested, xframes_nested, yframes_nested, frametime_nested, outx, outy):
        """
        MAIN FUNCTION

        Arguments:
        ----------
        stamp -- 2D stamp image containing target. Set to None for point sources
        xframes -- list of x-coordinate pixel position of target
                   in each frame
        yframes -- list of y-coordinate pixel position of target
                   in each frame
        frametime -- exposure time in seconds corresponding to one
                     detector readout (varies with subarray size)
        outx -- x-dimension size of the output aperture (2048 for
                full-frame)
        outy -- y-dimension size of the output aperture (2048 for
                full-frame)

        Returns:
        --------
        3D array containing the signal of the source in each frame
        of the integration
        """

        #in this update, stamps_nested, xframes, and yframes are now all nested lists.
        #sublists are for positions within a given frame.

        # Retrieve the list of nominal source locations in each frame,
        # defined as the postiion in the finel element of each nested list.
        xframes = [e[-1] for e in xframes_nested]
        yframes = [e[-1] for e in yframes_nested]


        #there should be no subsampling needed at this point. The nested list of psfs
        #have been evaluated at the proper subpixel locations

        print(f'xframes is {xframes}')
        print(f'yframes is {yframes}')


        # Make sure subsampling factor is an integer
        self.subsampx = int(self.subsampx)
        self.subsampy = int(self.subsampy)

        # Quick fix for the case where xinit,yinit are integers
        xinit = float(xframes_nested[1][0])
        yinit = float(yframes_nested[1][0])

        # List of times for all frames
        numframes = len(xframes_nested)-1
        #times = frametime * np.arange(-1, numframes)

        # Generate a list of locations at dist-pixel increments
        # between the beginning and ending locations
        #xs, ys = self.equidistantXY(xframes[0], yframes[0], xframes[-1],
        #                            yframes[-1], 1./self.subsampx)

        # List of sub-frame locations spanning the entire integration.
        # Replaces the call to self.equidistantXY() above.
        xs = flatten_nested_list(xframes_nested)
        ys = flatten_nested_list(yframes_nested)
        stamps = flatten_nested_list(stamps_nested)
        frametimes = flatten_nested_list(frametime_nested)

        # Subsample the stamp image
        #if stamp is not None:
        #    substamp = self.subsample(stamp, self.subsampx, self.subsampy)
        #else:
        #    substamp = None
        ##substamplen = substamp.shape


        #print(f'In moving_targets, initial stamp total signal is: {np.sum(stamp)}')
        #print(f'Total signal in subsampled stamp is: {np.sum(substamp)}')



        # Create the initial output frame
        #ystamplen, xstamplen = stamp.shape
        #minx = int(np.min([np.floor(xframes[0]) - np.ceil(xstamplen/2.),\
        #                   np.floor(xframes[-1])-np.ceil(xstamplen/2.)]))
        #maxx = int(np.max([np.floor(xframes[-1] + xstamplen/2.),\
        #                   np.floor(xframes[0]+xstamplen/2.)]))
        #miny = int(np.min([np.floor(yframes[0]) - np.ceil(ystamplen/2.),\
        #                   np.floor(yframes[-1])-np.ceil(ystamplen/2.)]))
        #maxy = int(np.max([np.floor(yframes[-1] + ystamplen/2.),\
        #                   np.floor(yframes[0]+ystamplen/2.)]))

        ystamplen, xstamplen = stamps_nested[0][0].shape
        minx = int(np.min([np.floor(xframes_nested[0][0]) - np.ceil(xstamplen/2.),\
                           np.floor(xframes_nested[-1][-1])-np.ceil(xstamplen/2.)]))
        maxx = int(np.max([np.floor(xframes_nested[-1][-1] + xstamplen/2.),\
                           np.floor(xframes_nested[0][0]+xstamplen/2.)]))
        miny = int(np.min([np.floor(yframes_nested[0][0]) - np.ceil(ystamplen/2.),\
                           np.floor(yframes_nested[-1][-1])-np.ceil(ystamplen/2.)]))
        maxy = int(np.max([np.floor(yframes_nested[-1][-1] + ystamplen/2.),\
                           np.floor(yframes_nested[0][0]+ystamplen/2.)]))


        # Don't let stamps fall off the edges of the output array
        mnx = minx
        mxx = maxx
        mny = miny
        mxy = maxy
        if minx < 0:
            mnx = 0
        if maxx > outx:
            mxx = int(outx-1)
        if miny < 0:
            mny = 0
        if maxy > outy:
            mxy = int(outy-1)

        # Subsample the output frame
        #totxpoints = np.min([outx,mxx-mnx+1])
        #totypoints = np.min([outy,mxy-mny+1])
        #outputframe0 = np.zeros((int(totypoints*self.subsampy),\
        #                         int(totxpoints*self.subsampx)))
        #outputframe1 = np.zeros((int(totypoints*self.subsampy),\
        #                         int(totxpoints*self.subsampx)))
        #outfull = np.zeros((numframes, outy, outx))
        #outsubshape = outputframe0.shape

        # Determine what size array is needed to hold all locations of the source
        # across all frames
        totxpoints = np.min([outx, mxx-mnx+1])
        totypoints = np.min([outy, mxy-mny+1])

        # Create outputframes that are just large enough to contain the source
        outputframe0 = np.zeros((int(totypoints), int(totxpoints)))
        outputframe1 = np.zeros((int(totypoints), int(totxpoints)))

        # Create the output integration at the full aperture size
        outfull = np.zeros((numframes, outy, outx))
        outsubshape = outputframe0.shape



        """
        minx - is the minimum x coordinate of the stamp image across the entire integration in the aperture coordinates
        maxx - is the maximum x coordinate of the stamp image across the entire integration in the aperture coordinates
        minx cannot be less than 0 or more than 2048

        totxpoints is the lesser of 2048 or maxx-minx

        xframessub is the (integer) distance between the x location of the source in aperture coords and minx (minimum x coord of stamp image across integration)
        xssub is the (integer) distance between the source location (in each subframe) and minx (minimum x coord of stamp images across integration)
        """


        # Translate the source location x and y values to the coordinates
        # of the output frame
        #deltacenterx = np.round(self.subsampx / 2. - 1 + 0.000001)
        #deltacentery = np.round(self.subsampy / 2. - 1 + 0.000001)
        #xframessub = np.round((xframes-mnx) * self.subsampx) + deltacenterx
        #yframessub = np.round((yframes-mny) * self.subsampy) + deltacentery
        #xssub = np.round((xs-mnx) * self.subsampx) + deltacenterx
        #yssub = np.round((ys-mny) * self.subsampy) + deltacentery


        # Same as above, but no need to deal with subsampled coordinates.

        # xframesub, yframesub are the coordinates of the source in each subframe, shifted to be
        # in the coordinate system of outputframe0/outputframe1
        xframessub = np.round((xframes_nested - mnx))
        yframessub = np.round((yframes_nested - mny))

        # xs is just the flattend version of xframes_nested. So xssub, yssub is the flattened
        # version of xframessub, yframessub
        xssub = np.round((xs - mnx))
        yssub = np.round((ys - mny))


        for i in range(1, numframes+1):
            # Find the velocity of the source during this frame
            outputframe1 = np.copy(outputframe0)




            ##################################################
            ### This block I think is just finding the subframe entries that
            ### correspond to the current frame. With the nested lists, we
            ### already know this info
            ### So we should be able to skip this
            #if xframessub[i-1] < xframessub[i]:
            #    goodxs = ((xssub > (xframessub[i-1]+1e-7)) & (xssub < (xframessub[i]-1e-7)))
            #else:
            #    goodxs = ((xssub > (xframessub[i]+1e-7)) & (xssub < (xframessub[i-1]-1e-7)))
            #if yframessub[i-1] < yframessub[i]:
            #    goodys = ((yssub > (yframessub[i-1]+1e-7)) & (yssub < (yframessub[i]-1e-7)))
            #else:
            #    goodys = ((yssub > (yframessub[i]+1e-7)) & (yssub < (yframessub[i-1]-1e-7)))

            #xsum = np.sum(goodxs)
            #ysum = np.sum(goodys)
            #if xsum >= ysum:
            #    good = goodxs
            #else:
            #    good = goodys
            ###################################################


            if (np.all((xframes[i-1:i+1]-xstamplen) > outx) or \
                (np.all((yframes[i-1:i+1]-ystamplen) > outy))):
                # If the stamp is completely above or to the right of the aperture,
                # outputframe1 is the same as the previous frame. (No stamp added.)
                outputframe1 = np.copy(outputframe0)
            elif (np.all((xframes[i-1:i+1]+xstamplen) < 0) or \
                  (np.all((yframes[i-1:i+1]+ystamplen) < 0))):
                # If the stamp is completely below or to the leftt of the aperture,
                # outputframe1 is the same as the previous frame. (No stamp added.)
                outputframe1 = np.copy(outputframe0)
            else:
                # If the stamp is at least partially within the aperture, add the moving stamp to the frame

                print(f'Frame {i}: xssub[good] is {xssub[good]} and yssub[good] is {yssub[good]}')


                #outputframe1 = self.inputMotion(outputframe1, substamp, xframessub[i-1:i+1],
                #                                yframessub[i-1:i+1],xssub[good],yssub[good],
                #                                frametime)
                outputframe1 = self.inputMotion(outputframe1, substamp, xframessub[i-1:i+1],
                                                yframessub[i-1:i+1],xssub[good],yssub[good],
                                                frametimes[i-1:i+1])

            outputframe0 = np.copy(outputframe1)


            print(f'Total signal in outputframe1 divided by frametime is: {np.sum(outputframe1)/frametime}')



            # Put the output frames back to the original resolution
            resampled = self.resample(outputframe1, self.subsampx, self.subsampy)


            print(f'Total signal in resampled outputframe divided by frametime is: {np.sum(resampled)/frametime}')



            resampylen, resampxlen = resampled.shape
            #outfull[i-1,mny:mxy+1,mnx:mxx+1] = resampled
            maxfully = mny + resampylen
            maxfullx = mnx + resampxlen
            maxrey = resampylen
            maxrex = resampxlen
            if (mny + resampylen) > outy:
                diffind = (mny + resampylen) - outy
                maxfully = outy
                maxrey -= diffind
            if (mnx + resampxlen) > outx:
                diffind = (mnx + resampxlen) - outx
                maxfullx = outx
                maxrex -= diffind
            outfull[i-1, mny:maxfully, mnx:maxfullx] = resampled[0:maxrey, 0:maxrex]
        return outfull

    def create_psf_stamp(self, x_location, y_location, psf_dim_x, psf_dim_y,
                         ignore_detector=False, segment_number=None):
        """From the gridded PSF model, location within the aperture, and
        dimensions of the stamp image (either the library PSF image, or
        the galaxy/extended stamp image with which the PSF will be
        convolved), evaluate the GriddedPSFModel at
        the appropriate location on the detector and return the PSF stamp

        Parameters
        ----------
        x_location : float
            X-coordinate of the PSF in the coordinate system of the
            aperture being simulated.

        y_location : float
            Y-coordinate of the PSF in the coordinate system of the
            aperture being simulated.

        psf_dim_x : int
            Number of columns of the array containing the PSF

        psf_dim_y : int
            Number of rows of the array containing the PSF

        ignore_detector : bool
            If True, the returned coordinates can have values outside the
            size of the subarray/detector (i.e. coords can be negative or
            larger than full frame). If False, coordinates are constrained
            to be on the detector.

        Returns
        -------
        full_psf : numpy.ndarray
            2D array containing the normalized PSF image. Total signal should
            be close to 1.0 (not exactly 1.0 due to asymmetries and distortion)
            Array will be cropped based on how much falls on or off the detector

        k1 : int
            Row number on the PSF/stamp image corresponding to the bottom-most
            row that overlaps the detector/aperture

        l1 : int
            Column number on the PSF/stamp image corresponding to the left-most
            column that overlaps the detector/aperture

        add_wings : bool
            Whether or not PSF wings are to be added to the PSF core
        """
        # PSF will always be centered
        psf_x_loc = psf_dim_x // 2
        psf_y_loc = psf_dim_y // 2

        # Translation needed to go from PSF core (self.psf_library)
        # coordinate system to the PSF wing coordinate system (i.e.
        # center the PSF core in the wing image)
        psf_wing_half_width_x = int(psf_dim_x // 2)
        psf_wing_half_width_y = int(psf_dim_y // 2)
        psf_core_half_width_x = int(self.psf_library_core_x_dim // 2)
        psf_core_half_width_y = int(self.psf_library_core_y_dim // 2)
        delta_core_to_wing_x = psf_wing_half_width_x - psf_core_half_width_x
        delta_core_to_wing_y = psf_wing_half_width_y - psf_core_half_width_y

        # This assumes a square PSF shape!!!!
        # If no wings are to be added, then we can skip all the wing-
        # and pixel phase-related work below.
        if ((self.add_psf_wings is False) or (delta_core_to_wing_x <= 0)):
            add_wings = False

            if segment_number is not None:
                library = self.psf_library[segment_number - 1]
            else:
                library = self.psf_library

            # Get coordinates decribing overlap between the evaluated psf
            # core and the full frame of the detector. We really only need
            # the xpts_core and ypts_core from this in order to know how
            # to evaluate the library
            # Note that we don't care about the pixel phase here.
            psf_core_dims = (self.psf_library_core_y_dim, self.psf_library_core_x_dim)
            xc_core, yc_core, xpts_core, ypts_core, (i1c, i2c), (j1c, j2c), (k1c, k2c), \
                (l1c, l2c) = self.create_psf_stamp_coords(x_location, y_location, psf_core_dims,
                                                          psf_core_half_width_x, psf_core_half_width_y,
                                                          coord_sys='full_frame',
                                                          ignore_detector=ignore_detector)

            # Skip sources that fall completely off the detector
            if None in [i1c, i2c, j1c, j2c, k1c, k2c, l1c, l2c]:
                return None, None, None, False

            # Step 4
            full_psf = library.evaluate(x=xpts_core, y=ypts_core, flux=1.0,
                                        x_0=xc_core, y_0=yc_core)
            k1 = k1c
            l1 = l1c

        else:
            add_wings = True
            # If the source subpixel location is beyond 0.5 (i.e. the edge
            # of the pixel), then we shift the wing->core offset by 1.
            # We also need to shift the location of the wing array on the
            # detector by 1
            x_phase = np.modf(x_location)[0]
            y_phase = np.modf(y_location)[0]
            x_location_delta = int(x_phase > 0.5)
            y_location_delta = int(y_phase > 0.5)
            if x_phase > 0.5:
                delta_core_to_wing_x -= 1
            if y_phase > 0.5:
                delta_core_to_wing_y -= 1

            # offset_x, and y below will not change because that is
            # the offset between the full wing array and the user-specified
            # wing array size

            # Get the psf wings array - first the nominal size
            # Later we may crop if the source is only partially on the detector
            full_wing_y_dim, full_wing_x_dim = self.psf_wings.shape
            offset_x = int((full_wing_x_dim - psf_dim_x) / 2)
            offset_y = int((full_wing_y_dim - psf_dim_y) / 2)

            full_psf = copy.deepcopy(self.psf_wings[offset_y:offset_y+psf_dim_y, offset_x:offset_x+psf_dim_x])

            # Get coordinates describing overlap between PSF image and the
            # full frame of the detector
            # Step 1
            xcenter, ycenter, xpts, ypts, (i1, i2), (j1, j2), (k1, k2), \
                (l1, l2) = self.create_psf_stamp_coords(x_location+x_location_delta, y_location+y_location_delta,
                                                        (psf_dim_y, psf_dim_x), psf_x_loc, psf_y_loc,
                                                        coord_sys='full_frame', ignore_detector=ignore_detector)

            if None in [i1, i2, j1, j2, k1, k2, l1, l2]:
                return None, None, None, False

            # Step 2
            # If the core of the psf lands at least partially on the detector
            # then we need to evaluate the psf library
            if ((k1 < (psf_wing_half_width_x + psf_core_half_width_x)) and
               (k2 > (psf_wing_half_width_x - psf_core_half_width_x)) and
               (l1 < (psf_wing_half_width_y + psf_core_half_width_y)) and
               (l2 > (psf_wing_half_width_y - psf_core_half_width_y))):

                # Step 3
                # Get coordinates decribing overlap between the evaluated psf
                # core and the full frame of the detector. We really only need
                # the xpts_core and ypts_core from this in order to know how
                # to evaluate the library
                # Note that we don't care about the pixel phase here.
                psf_core_dims = (self.psf_library_core_y_dim, self.psf_library_core_x_dim)
                xc_core, yc_core, xpts_core, ypts_core, (i1c, i2c), (j1c, j2c), (k1c, k2c), \
                    (l1c, l2c) = self.create_psf_stamp_coords(x_location, y_location, psf_core_dims,
                                                              psf_core_half_width_x, psf_core_half_width_y,
                                                              coord_sys='full_frame', ignore_detector=ignore_detector)

                if None in [i1c, i2c, j1c, j2c, k1c, k2c, l1c, l2c]:
                    return None, None, None, False

                # Step 4
                psf = self.psf_library.evaluate(x=xpts_core, y=ypts_core, flux=1.,
                                                x_0=xc_core, y_0=yc_core)

                # Step 5
                wing_start_x = k1c + delta_core_to_wing_x
                wing_end_x = k2c + delta_core_to_wing_x
                wing_start_y = l1c + delta_core_to_wing_y
                wing_end_y = l2c + delta_core_to_wing_y

                full_psf[wing_start_y:wing_end_y, wing_start_x:wing_end_x] = psf

            # Whether or not the core is on the detector, crop the PSF
            # to the proper shape based on how much is on the detector
            full_psf = full_psf[l1:l2, k1:k2]

        return full_psf, k1, l1, add_wings

    def create_psf_stamp_coords(self, aperture_x, aperture_y, stamp_dims, stamp_x, stamp_y,
                                coord_sys='full_frame', ignore_detector=False):
        """Calculate the coordinates in the aperture coordinate system
        where the stamp image wil be placed based on the location of the
        stamp image in the aperture and the size of the stamp image.

        Parameters
        ----------
        aperture_x : float
            X-coordinate of the PSF in the coordinate system of the
            aperture being simulated.

        aperture_y : float
            Y-coordinate of the PSF in the coordinate system of the
            aperture being simulated.

        stamp_dims : tup
            (x, y) dimensions of the stamp image that will be placed
            into the final seed image. This stamp image can be either the
            PSF image itself, or the stamp image of the galaxy/extended
            source that the PSF is going to be convolved with.

        stamp_x : float
            Location in x of source within the stamp image

        stamp_y : float
            Location in y of source within the stamp image

        coord_sys : str
            Inidicates which coordinate system to return coordinates for.
            Options are 'full_frame' for full frame coordinates, or
            'aperture' for aperture coordinates (including any expansion
            for grism source image)

        ignore_detector : bool
            If True, the returned coordinates can have values outside the
            size of the subarray/detector (i.e. coords can be negative or
            larger than full frame). If False, coordinates are constrained
            to be on the detector.

        Returns
        -------
        x_points : numpy.ndarray
            2D array of x-coordinates in the aperture coordinate system
            where the stamp image will fall.

        y_points : numpy.ndarray
            2D array of y-coordinates in the aperture coordinate system
            where the stamp image will fall.

        (i1, i2) : tup
            Beginning and ending x coordinates (in the aperture coordinate
            system) where the stamp image will fall

        (j1, j2) : tup
            Beginning and ending y coordinates (in the aperture coordinate
            system) where the stamp image will fall

        (k1, k2) : tup
            Beginning and ending x coordinates (in the stamp's coordinate
            system) that overlap the aperture

        (l1, l2) : tup
            Beginning and ending y coordinates (in the stamp's coordinate
            system) that overlap the aperture
        """
        if coord_sys == 'full_frame':
            xpos = aperture_x + self.subarray_bounds[0]
            ypos = aperture_y + self.subarray_bounds[1]
            out_dims_x = self.ffsize
            out_dims_y = self.ffsize
        elif coord_sys == 'aperture':
            xpos = aperture_x + self.coord_adjust['xoffset']
            ypos = aperture_y + self.coord_adjust['yoffset']
            out_dims_x = self.output_dims[1]
            out_dims_y = self.output_dims[0]

        stamp_y_dim, stamp_x_dim = stamp_dims

        # Get coordinates that describe the overlap between the stamp
        # and the aperture
        (i1, i2, j1, j2, k1, k2, l1, l2) = self.cropped_coords(xpos, ypos, (out_dims_y, out_dims_x),
                                                               stamp_x, stamp_y, stamp_dims,
                                                               ignore_detector=ignore_detector)

        # If the stamp is completely off the detector, use dummy arrays
        # for x_points and y_points
        if j1 is None or j2 is None or i1 is None or i2 is None:
            x_points = np.zeros((2, 2))
            y_points = x_points
        else:
            y_points, x_points = np.mgrid[j1:j2, i1:i2]

        return xpos, ypos, x_points, y_points, (i1, i2), (j1, j2), (k1, k2), (l1, l2)


    def resample(self, frame, sampx, sampy):
        """
        Return subsampled image back to original resolution

        Arguments:
        ----------
        frame -- subsampled image
        sampx -- x-dimension subsampling factor (e.g. 3 means 3x oversampled
                 compared to original image
        sampy -- y-dimension subsampling factor

        Returns:
        --------
        resampled image
        """
        framey, framex = frame.shape
        newframe = np.zeros((int(framey/sampy), int(framex/sampx)))
        newframey, newframex = newframe.shape

        for j in range(newframey):
            for i in range(newframex):
                newframe[j, i] = np.sum(frame[sampy*j:sampy*(j+1), sampx*i:sampx*(i+1)])
        return newframe

    def coordCheck(self, center, len_stamp, len_out):
        """
        Find indexes of stamp and frame to use
        given that the stamp may fall off the edge
        of the frame. Works on only one coordinate dimension

        Arguments:
        ----------
        center -- coordinate (in full aperture coords of the center of the
                  stamp image.
        len_stamp -- Size of the stamp image
        len_out -- Size of the full aperture image

        Returns:
        --------
        x and y coordinates corresponding to the beginning and ending
        (i.e. top and bottom for y-dimension, or left and right for x-
        dimension) of the stamp image on the full frame aperture, as
        well as the beginning and ending coordinates within the stamp
        image that fall onto the full frame aperture.
        """
        outxmin = center - len_stamp/2
        outxmax = outxmin + len_stamp
        stampxmin = 0
        stampxmax = len_stamp

        # Left edge of stamp is off the edge of output
        if outxmin < 0:
            if outxmin >= (0.-len_stamp):
                stampxmin = 0 - outxmin
                outxmin = 0
            else:
                # Here the image is completely off the output frame
                stampxmin = np.nan
                outxmin = np.nan
                stampxmax = np.nan
                outxmax = np.nan

        # Right edge of stamp is off the edge of the output
        if outxmax > len_out:
            if outxmax <= (len_out + len_stamp):
                delta = outxmax - len_out
                stampxmax = len_stamp - delta
                outxmax = len_out
            else:
                # Here the image is completely off the left side of output
                outxmax = np.nan
                stampxmax = np.nan
                outxmin = np.nan
                stampxmin = np.nan

        indexes = [outxmin, outxmax, stampxmin, stampxmax]
        if np.all(np.isfinite(indexes)):
            ioutxmin = int(outxmin)
            ioutxmax = int(outxmax)
            istampxmin = int(stampxmin)
            istampxmax = int(stampxmax)
            dout = ioutxmax - ioutxmin
            dstamp = istampxmax - istampxmin

            if dout == dstamp:
                pass
            elif dout == (dstamp+1):
                if istampxmin > 0:
                    istampxmin -= 1
                else:
                    istampxmax += 1
            elif dstamp == (dout+1):
                if ioutxmin > 0:
                    ioutxmin -= 1
                else:
                    ioutxmax += 1
            else:
                self.logger.error("WARNING: bad stamp/output match. Quitting.")
                raise ValueError('Bad stamp/output match.')
            return ioutxmin, ioutxmax, istampxmin, istampxmax
        else:
            # If values are NaN then we can't change them to integers
            return outxmin, outxmax, stampxmin, stampxmax

    def inputMotion(self, inframe, source, xbounds, ybounds, xs, ys, frame_times):
        """
        Smear out the source to create an output frame image
        given the necessary info about the source location and velocity

        Parameters
        ----------
        inframe : numpy.ndarray
            2D array representing the image

        source : numpy.ndarray
            2D stamp image containing the source

        xbounds : list
            2-element list containing the starting and ending x-dimension
            coordinates of the source (i.e. location corresponding to the
            beginning and ending of the frame)

        ybounds : list
            2-element list containing the starting and ending y-dimension
            coordinates of the source

        xs : list
            x-coordinate positions of the source

        ys : list
            y-coordinate positions of the source

        frame_time : float
            exposure time of a single frame

        Returns
        -------
        inframe : numpy.ndarray
            With streaked source added
        """
        stop
        frameylen,framexlen = inframe.shape
        srcylen,srcxlen = source.shape
        xlist = np.append(xs, xbounds[1])
        ylist = np.append(ys, ybounds[1])
        xlist = np.insert(xlist, 0, xbounds[0])
        ylist = np.insert(ylist, 0, ybounds[0])
        xlist = np.round(xlist)
        ylist = np.round(ylist)

        # exposure time between two of these equidistant points
        pt_exptime = frame_time / len(xlist)

        print(f'In inputMotion, list of xlist is {len(xlist)} elements')
        print(f'xs: {xs}')
        print(f'xbounds: {xbounds}')
        print(f'xlist: {xlist}')

        print(f'Partial frametime is {pt_exptime}. Total frametime is {frame_time}')

        for i in range(0,len(xlist)):
            outxmin, outxmax, stampxmin, stampxmax = self.coordCheck(xlist[i], srcxlen, framexlen)
            outymin, outymax, stampymin, stampymax = self.coordCheck(ylist[i], srcylen, frameylen)
            outcoords = np.array([outxmin,outxmax,outymin,outymax])

            # If any of the coordinates are set to NaN, then the stamp image is completely off
            # the output frame and it shouldn't be added
            if np.all(np.isfinite(outcoords)):
                #dist = np.sqrt((xlist[i]-xlist[i-1])**2 + (ylist[i]-ylist[i-1])**2)

                # Evaluate gridded psf model to create PSF for this subpixel location
                #create_psf_stamp is in catalog_seed_image...
                eval_psf, minx, miny, wings_added = self.create_psf_stamp(pixelx, pixely, psf_x_dim, psf_x_dim,
                                                                      ignore_detector=True)

                inframe[outymin:outymax, outxmin:outxmax] += (source[stampymin:stampymax, stampxmin:stampxmax] * pt_exptime)


                #if i==len(xlist)-1:
                print(f'     In inputMotion, total of smeared signal divided by frametime is: {np.sum(inframe[outymin:outymax, outxmin:outxmax])}')


        if len(xs) > 0:
            stop
        return inframe

    def subsample(self, image, factorx, factory):
        """
        Subsample the input image

        Arguments:
        ----------
        image -- 2D image
        factorx -- factor in the x-dimension to subsample the image
                 (e.g. factorx=2 will break each pixel into 2 pixels
                  in the x dimension)
        factory -- factor in the y-dimension to subsample the image

        Setting factorx = 2, factory = 2 will break each pixel in the
        original image into a 2x2 grid of pixels

        Returns:
        --------
        Subsampled image
        """
        ydim, xdim = image.shape
        substamp = np.zeros((ydim*factory, xdim*factorx))

        for i in range(xdim):
            for j in range(ydim):
                substamp[factory*j:factory*(j+1), factorx*i:factorx*(i+1)] = image[j, i] / (factorx * factory)
        return substamp

    def equidistantXY(self,xstart, ystart, xend, yend, dist):
        """
        Return a list of x,y positions that are equidistant
        between the beginning and ending positions, with
        a distance of dist pixels between them

        Arguments:
        ----------
        xstart -- beginning x coordinate
        ystart -- beginning y coordinate
        xend -- ending x coordinate
        yend -- ending y coordinate
        dist -- distance in pixels between adjacent positions
        """
        xlen = 0
        ylen = 0

        deltax = xend - xstart
        deltay = yend - ystart
        ang = np.arctan2(deltay, deltax)

        dx = np.cos(ang) * dist
        dy = np.sin(ang) * dist

        if dx != 0.:
            xs = np.arange(xstart, xend+dx/2, dx)
            xlen = len(xs)
        else:
            # Motion parallel to y axis
            xs = np.array([xstart])

        if dy != 0:
            ys = np.arange(ystart, yend+dy/2, dy)
            ylen = len(ys)
        else:
            # Motion parallel to x asis
            ys = np.array([ystart])

        # Make sure lengths agree
        if xlen == 0:
            xs = np.zeros(ylen) + xstart
        if ylen == 0:
            ys = np.zeros(xlen) + ystart

        return xs, ys

    def radecPerFrame(self, ra0, dec0, ravel, decvel, time):
        """
        Generate a list of RA,Dec locations for a source in
        a series of frames

        Arguments:
        ----------
        ra0 -- Initial RA value of source
        dec0 -- Initial Dec value of source
        ravel -- Velocity of source in the RA direction
        decvel -- Velocity of source in the Dec direction
        time -- List of times corresponding to all frames

        Returns:
        --------
        List of RA, Dec positions corresponding to all input times
        """
        ra = ra0 + (ravel * time)
        dec = dec0 + (decvel * time)
        return ra, dec

    def xyPerFrame(self, velocity, time, ang, x0, y0):
        """
        Generate list of x,y positions for a source given an
        initial position, velocity, and velocity angle

        Arguments:
        ----------
        velocity -- Velocity of source
        time -- List of times at which we want to find positions
        ang -- Angle (in radians) at which source is traveling.
               An ang of zero corresponds to moving along the +x axis.
               An ang of np.pi/2 corresponds to moveing along the +y axis.
        x0 -- Initial x coordinate of source
        y0 -- Initial y coordinate of source

        Returns:
        --------
        Tuple of x-list, y-list souce locations
        """
        ratex = velocity * np.cos(ang)
        ratey = velocity * np.sin(ang)

        # x,y in each frame
        xs = x0 + ratex*time
        ys = y0 + ratey*time

        return xs,ys
