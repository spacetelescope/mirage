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

import sys

import numpy as np
from astropy.io import fits

class MovingTarget():

    def __init__(self):
        self.verbose = False
        self.subsampx = 3
        self.subsampy = 3

    def create(self, stamp, xframes, yframes, frametime, outx, outy):
        """
        MAIN FUNCTION

        Arguments:
        ----------
        stamp -- 2D stamp image containing target
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

        # Make sure subsampling factor is an integer
        self.subsampx = np.int(self.subsampx)
        self.subsampy = np.int(self.subsampy)

        # Quick fix for the case where xinit,yinit are integers
        xinit = np.float(xframes[1])
        yinit = np.float(yframes[1])

        # Change position angle to radians
        #posang = posang * np.pi / 180.

        # List of times for all frames
        numframes = len(xframes)-1
        times = frametime * np.arange(-1, numframes)

        # Generate a list of locations at dist-pixel increments
        # between the beginning and ending locations
        xs, ys = self.equidistantXY(xframes[0], yframes[0], xframes[-1],
                                    yframes[-1], 1./self.subsampx)

        # Subsample the PSF image
        substamp = self.subsample(stamp, self.subsampx, self.subsampy)
        substamplen = substamp.shape

        # Create the initial output frame
        ystamplen, xstamplen = stamp.shape
        #minx = np.min([np.floor(xinit) - np.ceil(xstamplen/2.),np.floor(xframes[-1])-np.ceil(xstamplen/2.)])
        minx = np.int(np.min([np.floor(xframes[0]) - np.ceil(xstamplen/2.),\
                              np.floor(xframes[-1])-np.ceil(xstamplen/2.)]))
        maxx = np.int(np.max([np.floor(xframes[-1] + xstamplen/2.),\
                              np.floor(xframes[0]+xstamplen/2.)]))
        #miny = np.min([np.floor(yinit) - np.ceil(ystamplen/2.),np.floor(yframes[-1])-np.ceil(ystamplen/2.)])
        miny = np.int(np.min([np.floor(yframes[0]) - np.ceil(ystamplen/2.),\
                              np.floor(yframes[-1])-np.ceil(ystamplen/2.)]))
        maxy = np.int(np.max([np.floor(yframes[-1] + ystamplen/2.),\
                              np.floor(yframes[0]+ystamplen/2.)]))

        # Don't let stamps fall off the edges of the output array
        mnx = minx
        mxx = maxx
        mny = miny
        mxy = maxy
        if minx < 0:
            mnx = 0
        if maxx > outx:
            mxx = np.int(outx-1)
        if miny < 0:
            mny = 0
        if maxy > outy:
            mxy = np.int(outy-1)

        # Subsample the output frame
        totxpoints = np.min([outx,mxx-mnx+1])
        totypoints = np.min([outy,mxy-mny+1])
        outputframe0 = np.zeros((np.int(totypoints*self.subsampy),\
                                 np.int(totxpoints*self.subsampx)))
        outputframe1 = np.zeros((np.int(totypoints*self.subsampy),\
                                 np.int(totxpoints*self.subsampx)))
        outfull = np.zeros((numframes, outy, outx))
        outsubshape = outputframe0.shape

        # Translate the source location x and y values to the coordinates
        # of the output frame
        deltacenterx = np.round(self.subsampx / 2. - 1 + 0.000001)
        deltacentery = np.round(self.subsampy / 2. - 1 + 0.000001)
        xframessub = np.round((xframes-mnx) * self.subsampx) + deltacenterx
        yframessub = np.round((yframes-mny) * self.subsampy) + deltacentery
        xssub = np.round((xs-mnx) * self.subsampx) + deltacenterx
        yssub = np.round((ys-mny) * self.subsampy) + deltacentery

        for i in range(1,numframes+1):
            # Find the velocity of the source during this frame
            xvelocity = (xframes[i] - xframes[i-1]) / frametime
            yvelocity = (yframes[i] - yframes[i-1]) / frametime
            secPerPix = 1. / np.sqrt(xvelocity*xvelocity + yvelocity*yvelocity)
            outputframe1 = np.copy(outputframe0)

            if xframessub[i-1] < xframessub[i]:
                goodxs = ((xssub > (xframessub[i-1]+1e-7)) & (xssub < (xframessub[i]-1e-7)))
            else:
                goodxs = ((xssub > (xframessub[i]+1e-7)) & (xssub < (xframessub[i-1]-1e-7)))
            if yframessub[i-1] < yframessub[i]:
                goodys = ((yssub > (yframessub[i-1]+1e-7)) & (yssub < (yframessub[i]-1e-7)))
            else:
                goodys = ((yssub > (yframessub[i]+1e-7)) & (yssub < (yframessub[i-1]-1e-7)))

            xsum = np.sum(goodxs)
            ysum = np.sum(goodys)
            if xsum >= ysum:
                good = goodxs
            else:
                good = goodys

            if (np.all((xframes[i-1:i+1]-xstamplen) > outx) or \
                (np.all((yframes[i-1:i+1]-ystamplen) > outy))):
                outputframe1 = np.copy(outputframe0)
            elif (np.all((xframes[i-1:i+1]+xstamplen) < 0) or \
                  (np.all((yframes[i-1:i+1]+ystamplen) < 0))):
                outputframe1 = np.copy(outputframe0)
            else:
                outputframe1 = self.inputMotion(outputframe1, substamp, xframessub[i-1:i+1],
                                                yframessub[i-1:i+1],xssub[good],yssub[good],
                                                secPerPix)

            outputframe0 = np.copy(outputframe1)

            # Put the output frames back to the original resolution
            resampled = self.resample(outputframe1, self.subsampx, self.subsampy)
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
        newframe = np.zeros((np.int(framey/sampy), np.int(framex/sampx)))
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
            ioutxmin = np.int(outxmin)
            ioutxmax = np.int(outxmax)
            istampxmin = np.int(stampxmin)
            istampxmax = np.int(stampxmax)
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
                print("WARNING: bad stamp/output match. Quitting.")
                sys.exit()
            return ioutxmin, ioutxmax, istampxmin, istampxmax
        else:
            # If values are NaN then we can't change them to integers
            return outxmin, outxmax, stampxmin, stampxmax

    def inputMotion(self, inframe, source, xbounds, ybounds, xs, ys, secperpix):
        """
        Smear out the source to create an output frame image
        given the necessary info about the source location and velocity

        Arguments:
        ----------
        inframe -- 2D array representing the image
        source -- 2D stamp image containing the source
        xbounds -- 2-element list containing the starting and ending x-dimension
                   coordinates of the source (i.e. location corresponding to the
                   beginning and ending of the frame)
        ybounds -- 2-element list containing the starting and ending y-dimension
                   coordinates of the source
        xs -- list of x-coordinate positions of the source
        ys -- list of y-coordinate positions of the source
        secperpix -- Inverse velocity of the source, in seconds per pixel
        """
        frameylen,framexlen = inframe.shape
        srcylen,srcxlen = source.shape
        xlist = np.append(xs, xbounds[1])
        ylist = np.append(ys, ybounds[1])
        xlist = np.insert(xlist, 0, xbounds[0])
        ylist = np.insert(ylist, 0, ybounds[0])
        xlist = np.round(xlist)
        ylist = np.round(ylist)

        for i in range(1,len(xlist)):
            outxmin, outxmax, stampxmin, stampxmax = self.coordCheck(xlist[i], srcxlen, framexlen)
            outymin, outymax, stampymin, stampymax = self.coordCheck(ylist[i], srcylen, frameylen)
            outcoords = np.array([outxmin,outxmax,outymin,outymax])

            # If any of the coordinates are set to NaN, then the stamp image is completely off
            # the output frame and it shouldn't be added
            if np.all(np.isfinite(outcoords)):
                dist = np.sqrt((xlist[i]-xlist[i-1])**2 + (ylist[i]-ylist[i-1])**2)
                inframe[outymin:outymax, outxmin:outxmax] += (source[stampymin:stampymax, stampxmin:stampxmax]*secperpix*dist)
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
                substamp[factory*j:factory*(j+1), factorx*i:factorx*(i+1)] = image[j, i]
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
