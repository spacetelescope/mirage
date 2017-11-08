#! /usr/bin/env python

'''
Class that simulates moving targets as viewed on a JWST instrument. Written
to accompany ramp_simulator.py

'''

from astropy.io import fits
import numpy as np
import sys


class MovingTarget():

    def __init__(self):
        self.verbose = False
        self.subsampx = 3
        self.subsampy = 3

    #def create(self,stamp,xinit,yinit,rate,posang,frametime,numframes,subsample_factor,outx,outy):
    #def create(self,stamp,rainit,decinit,rarate,decrate,frametime,numframes,subsample_factor,outx,outy):
    def create(self,stamp,xframes,yframes,frametime,outx,outy):
        #inputs include the stamp image of the object, the initial
        #x and y location, the rate at which the object is moving
        #(in pixels/sec) and the position angle (degrees), to describe the
        #direction in which the object is moving.
        
        #assume xinit is ra,dec. time calc stays same.
        #velocity is arcsec/sec
        #assume an ra_velocity and dec_velocity input
        #xframes,yframes then in ra,dec
        #translate xframes,yframes to x,y

        #make sure subsampling factor is an integer
        self.subsampx = np.int(self.subsampx)
        self.subsampy = np.int(self.subsampy)

        #quick fix for the case where xinit,yinit are integers
        xinit = np.float(xframes[1])
        yinit = np.float(yframes[1])

        #change position angle to radians
        #posang = posang * np.pi / 180.

        #inverse rate
        #rate = np.sqrt(rarate*rarate + decrate*decrate)
        #secPerPix = 1./ rate
        #CALCULATE THIS WITHIN EACH FRAME
        
        #list of times for all frames
        numframes = len(xframes)-1
        times = frametime * np.arange(-1,numframes)

        #list of RA,Dec positions for all frames
        #raframes,decframes = self.radecPerFrame(xinit,yinit,rarate,decrate,time)
        #raframes = rainit + rarate*times
        #decrames = decinit + decrate*times
        
        #calculate the x,y location of the object in each frame
        #Include the x,y location at the time of detector reset,
        #so that we can generate a realistic frame 0.
        #xframes,yframes = self.xyPerFrame(rate,times,posang,xinit,yinit)

        #generate a list of locations at dist-pixel increments
        #between the beginning and ending locations
        xs, ys = self.equidistantXY(xframes[0],yframes[0],xframes[-1],yframes[-1],1./self.subsampx) #,posang)
        #print('xs,ys',xs,ys)

        #subsample the PSF image
        substamp = self.subsample(stamp,self.subsampx,self.subsampy)
        substamplen = substamp.shape

        #create the initial output frame
        ystamplen,xstamplen = stamp.shape
        #minx = np.min([np.floor(xinit) - np.ceil(xstamplen/2.),np.floor(xframes[-1])-np.ceil(xstamplen/2.)])
        minx = np.int(np.min([np.floor(xframes[0]) - np.ceil(xstamplen/2.),np.floor(xframes[-1])-np.ceil(xstamplen/2.)]))        

        #print(xinit,xstamplen/2,xframes[-1])
        #print("intially, mnx is: {}".format(minx))
        maxx = np.int(np.max([np.floor(xframes[-1] + xstamplen/2.),np.floor(xframes[0]+xstamplen/2.)]))
        #print('maxx is {}. xframes[0] and [-1] are {},{}, stamplen is {}'.format(maxx,xframes[0],xframes[-1],xstamplen))
        #miny = np.min([np.floor(yinit) - np.ceil(ystamplen/2.),np.floor(yframes[-1])-np.ceil(ystamplen/2.)])
        miny = np.int(np.min([np.floor(yframes[0]) - np.ceil(ystamplen/2.),np.floor(yframes[-1])-np.ceil(ystamplen/2.)]))

        maxy = np.int(np.max([np.floor(yframes[-1] + ystamplen/2.),np.floor(yframes[0]+ystamplen/2.)]))

        #Don't let stamps fall off the edges of the output array
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

        #print("Now, minx is {}".format(mnx))
        #sys.exit()

        #subsample the output frame
        totxpoints = np.min([outx,mxx-mnx+1])
        totypoints = np.min([outy,mxy-mny+1])
        #print('totxpoints,outx,mxx,mnx',totxpoints,outx,mxx,mnx)
        #print('totypoints,outy,mxy,mny',totypoints,outy,mxy,mny)
        #print(totypoints,self.subsampy,totxpoints,self.subsampx,outx,mnx,mxx,mny,mxy)
        #print(minx,maxx,miny,maxy)
        outputframe0 = np.zeros((np.int(totypoints*self.subsampy),np.int(totxpoints*self.subsampx)))
        #sys.exit()
        outputframe1 = np.zeros((np.int(totypoints*self.subsampy),np.int(totxpoints*self.subsampx)))
        outfull = np.zeros((numframes,outy,outx))
        outsubshape = outputframe0.shape

        #print('x,y stamp size',xstamplen,ystamplen)
        #print("miny,maxy in full frame coords: ",miny,maxy)
        #print("minx,maxx in full frame coords: ",minx,maxx)
        #print('xframes in full frame coords:',xframes)
        #print('yframes in full frame coords:',yframes)

        #translate the source location x and y values to the coordinates
        #of the output frame
        deltacenterx = np.round(self.subsampx / 2. - 1 + 0.000001)
        deltacentery = np.round(self.subsampy / 2. - 1 + 0.000001)
        xframessub = np.round((xframes-mnx) * self.subsampx) + deltacenterx
        yframessub = np.round((yframes-mny) * self.subsampy) + deltacentery
        xssub = np.round((xs-mnx) * self.subsampx) + deltacenterx
        yssub = np.round((ys-mny) * self.subsampy) + deltacentery

        #print('xssub',xssub)
        #print('yssub',yssub)

        #now create the multiple-frame integration using subsampled frames

        #check to see if the source is off, or partially off the detector. 
        #Adjust coordinate indexes as necessary, and keep track of coordinates
        #in the stamp image
        #outxmin,outxmax,stampxmin,stampxmax = self.coordCheck(xframessub[1],substamplen[1],outsubshape[1])
        #outymin,outymax,stampymin,stampymax = self.coordCheck(xframessub[0],substamplen[0],outsubshape[0])

        #print(outxmin,outxmax,stampxmin,stampxmax)
        #print(outymin,outymax,stampymin,stampymax)
        #outcoords = np.array([outxmin,outxmax,outymin,outymax])
#       # outputframe[0,yframessub[0]-substamplen[0]/2:yframessub[0]+substamplen[0]/2,xframessub[0]-substamplen[1]/2:xframessub[0]+substamplen[1]/2] += substamp

        ##if any of the coordinates are set to NaN, then the stamp image is completely off
        ##the output frame and it shouldn't be added
        #if np.all(np.isfinite(outcoords)):
        #    outputframe[0,outymin:outymax,outxmin,outxmax] += substamp[stampymin:stampymax,stampxmin,stampxmax]

        
        for i in range(1,numframes+1):
            #for frames after the 0th, start with the previous frame
            #if i != 0:
                #outputframe[i-1,:,:] = np.copy(outputframe[i-2,:,:])

            #print("Working on frame number: {}".format(i-1))

            #find the velocity of the source during this frame
            xvelocity = (xframes[i]-xframes[i-1]) / frametime
            yvelocity = (yframes[i]-yframes[i-1]) / frametime
            secPerPix = 1. / np.sqrt(xvelocity*xvelocity + yvelocity*yvelocity)
            #print('velocities:',xvelocity,yvelocity,secPerPix)
            
            #print("Trying to calculate xs, xssub separate for each frame! Check for correctness!!!!")
            #xs, ys = self.equidistantXY(xframes[i-1],yframes[i-1],xframes[1],yframes[1],1./self.subsampx)
            #xssub = np.round((xs-mnx) * self.subsampx) + deltacenterx
            #yssub = np.round((ys-mnx) * self.subsampy) + deltacentery

            outputframe1 = np.copy(outputframe0)
            
            #print('LIMITS TO FIND STEPS: {}, {}'.format(xframessub[i-1],xframessub[i]))
            #print(yssub,yframessub[i-1],yframessub[i])
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

            #print('goodxs:',xssub[good])
            #print('goodys:',yssub[good])

            #print("making frame:")
            #print(xframessub[i-1:i+1],yframessub[i-1:i+1],xssub[good],yssub[good],secPerPix)

            #outputframe[i-1,:,:] = self.inputMotion(outputframe[i-1,:,:],substamp,xframessub[i-1:i+1],yframessub[i-1:i+1],xssub[good],yssub[good],secPerPix)

            #print("  ")
            #print('xstart xend, ystart, yend, xmidpts, ymidpts')
            #print(xframessub[i-1:i+1],yframessub[i-1:i+1],xssub[good],yssub[good])
            #print("  ")

            #goodones = xssub[good]
            #print(xframessub[i-1],goodones[0])
            #if xframessub[i-1] > goodones[0]:
            #    print("frame limit is greater.")
            #elif xframessub[i-1] < goodones[0]:
            #    print("xssub is greater")
            #elif xframessub[i-1] == goodones[0]:
            #    print('they are equal')

            #print('Going into inputMotion:')
            #print('outputframe1 shape',outputframe1.shape)
            #print('subsampled stamp shape',substamp.shape)
            #print('source location: x: {} to {}, y: {} to {}'.format(xframessub[i-1],xframessub[i],yframessub[i-1],yframessub[i]))
            #print('locations to use:',xssub[good],yssub[good])

            #if (np.all((xframessub[i-1:i+1]-substamplen[1]) > outsubshape[1]) or (np.all((yframessub[i-1:i+1]-substamplen[0]) > outsubshape[0]))):
            #    print(xframessub[i-1:i+1])
            #    print(substamplen)
            #    print(outsubshape[1])
            #    print("For frame {}, entire stamp is left or above the fov".format(i-1))
            #    outputframe1 = np.copy(outputframe0)
            #elif (np.all((xframessub[i-1:i+1]+substamplen[1]) < 0) or (np.all((yframessub[i-1:i+1]+substamplen[0]) < 0))):
            #    print(xframessub[i-1:i+1])
            #    print(substamplen)
            #    print("0.")
            #    print("For frame {}, entire stamp is right or below the fov".format(i-1))
            #    outputframe1 = np.copy(outputframe0)

            if (np.all((xframes[i-1:i+1]-xstamplen) > outx) or (np.all((yframes[i-1:i+1]-ystamplen) > outy))):
                #print(xframes[i-1:i+1])
                #print(xstamplen)
                #print(outx)
                #print("For frame {}, entire stamp is left or above the fov".format(i-1))
                outputframe1 = np.copy(outputframe0)
            elif (np.all((xframes[i-1:i+1]+xstamplen) < 0) or (np.all((yframes[i-1:i+1]+ystamplen) < 0))):
                #print(xframes[i-1:i+1])
                #print(xstamplen)
                #print("0.")
                #print("For frame {}, entire stamp is right or below the fov".format(i-1))
                outputframe1 = np.copy(outputframe0)
            else:
                #print(xframes[i-1:i+1])
                #print(xstamplen)
                #print(xframes[i-1:i+1]-xstamplen,xframes[i-1:i+1]+xstamplen)
                #print("For frame {}, stamp is at least partially within the fov".format(i-1))
                #print("going into inputmotion")
                #print(xframessub[i-1:i+1],yframessub[i-1:i+1],xssub[good],yssub[good],secPerPix)
                outputframe1 = self.inputMotion(outputframe1,substamp,xframessub[i-1:i+1],yframessub[i-1:i+1],xssub[good],yssub[good],secPerPix)

                #h0 = fits.PrimaryHDU()
                #h1 = fits.ImageHDU(outputframe1)
                #hl = fits.HDUList([h0,h1])
                #oname = 'test.fits'
                #hl.writeto(oname,clobber=True)
                #print(mnx,mxx,mny,mxy)
                #print(xframes)
                #sys.exit()


            outputframe0 = np.copy(outputframe1)
            
            #oname = 'test_movingtarg_out_frame_'+str(i)+'.fits'
            #h0 = fits.PrimaryHDU(outputframe1)
            #hl = fits.HDUList([h0])
            #hl.writeto(oname,clobber=True)


            #put the output frames back to the original resolution
            #print("miny,maxy: ",miny,maxy)
            #print(outfull[i-1,mny:mxy+1,mnx:mxx+1].shape,test.shape)
            #print(outputframe1.shape,test.shape)
            #print(mnx,mxx,mny,mxy)
            #check for sources that fall off the edges of the output array
            #print(mny,mxy,mnx,mxx)

            #print(outputframe1.shape)
            #print(self.subsampx,self.subsampy)
            outfull[i-1,mny:mxy+1,mnx:mxx+1] = self.resample(outputframe1,self.subsampx,self.subsampy)
        
        #h0 = fits.PrimaryHDU()
        #h1 = fits.ImageHDU(outfull)
        #hl = fits.HDUList([h0,h1])
        #oname = 'test_movingtarg_out_frame_origres'+str(i)+'.fits'
        #hl.writeto(oname,clobber=True)
        #sys.exit()

        return outfull


    def resample(self,frame,sampx,sampy):
        #return subsampled image back to original resolution
        framey,framex = frame.shape
        newframe = np.zeros((np.int(framey/sampy),np.int(framex/sampx)))
        newframey,newframex = newframe.shape

        for j in range(newframey):
            for i in range(newframex):
                newframe[j,i] = np.sum(frame[sampy*j:sampy*(j+1),sampx*i:sampx*(i+1)])
        return newframe
        

    def coordCheck(self,center,len_stamp,len_out):
        # Find indexes of stamp and frame to use
        # given that the stamp may fall off the edge
        # of the frame
        outxmin = center - len_stamp/2
        outxmax = outxmin + len_stamp
        stampxmin = 0
        stampxmax = len_stamp

        #print('center,len',center,len_stamp)
        #print('before checks!!!',outxmin,outxmax,stampxmin,stampxmax)

        #left edge of stamp is off the edge of output
        if outxmin < 0:
            #print('left edge of stamp is off')
            if outxmin >= (0.-len_stamp):
                stampxmin = 0 - outxmin
                outxmin = 0
            else:
                #here the image is completely off the output frame
                stampxmin = np.nan
                outxmin = np.nan
                stampxmax = np.nan
                outxmax = np.nan

        #right edge of stamp is off the edge of the output
        if outxmax > len_out:
            #print('right edge of stamp is off')
            if outxmax <= (len_out+len_stamp):
                delta = outxmax - len_out
                stampxmax = len_stamp - delta
                outxmax = len_out
            else:
                #here the image is completely off the left side of output
                outxmax = np.nan
                stampxmax = np.nan
                outxmin = np.nan
                stampxmin = np.nan

        #print('center,len',center,len_stamp)
        #print('after edge checks!!!',outxmin,outxmax,stampxmin,stampxmax)

                
        indexes = [outxmin,outxmax,stampxmin,stampxmax]
        if np.all(np.isfinite(indexes)):
            ioutxmin = np.int(outxmin)
            ioutxmax = np.int(outxmax)
            istampxmin = np.int(stampxmin)
            istampxmax = np.int(stampxmax)
            dout = ioutxmax - ioutxmin
            dstamp = istampxmax - istampxmin
            #print('before fix!!!',ioutxmin,ioutxmax,istampxmin,istampxmax)
            #print(dout,dstamp)

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

            #print('center,len',center,len_stamp)
            #print('after fix!!!',ioutxmin,ioutxmax,istampxmin,istampxmax)
            #print(ioutxmax-ioutxmin,istampxmax-istampxmin)   
            return ioutxmin,ioutxmax,istampxmin,istampxmax
        else:
            # if values are NaN then we can't change them to integers
            return outxmin,outxmax,stampxmin,stampxmax


    def inputMotion(self,inframe,source,xbounds,ybounds,xs,ys,secperpix):
        #smear out the source to create an output frame image
        #given the necessary info about the source location and velocity
        frameylen,framexlen = inframe.shape
        srcylen,srcxlen = source.shape
        xlist = np.append(xs,xbounds[1])
        ylist = np.append(ys,ybounds[1])
        xlist = np.insert(xlist,0,xbounds[0])
        ylist = np.insert(ylist,0,ybounds[0])


        #print('xlist',xlist)
        #print('ylist',ylist)
        xlist = np.round(xlist)
        ylist = np.round(ylist)
        #print('xlist',xlist)
        #print('ylist',ylist)

        for i in range(1,len(xlist)):
            #print('Working on location {},{}'.format(xlist[i],ylist[i]))
            outxmin,outxmax,stampxmin,stampxmax = self.coordCheck(xlist[i],srcxlen,framexlen)
            outymin,outymax,stampymin,stampymax = self.coordCheck(ylist[i],srcylen,frameylen)
            #print('x',outxmin,outxmax,stampxmin,stampxmax,srcxlen,frameylen)
            #print('y',outymin,outymax,stampymin,stampymax,srcylen,frameylen)

            outcoords = np.array([outxmin,outxmax,outymin,outymax])

            #if any of the coordinates are set to NaN, then the stamp image is completely off
            #the output frame and it shouldn't be added
            if np.all(np.isfinite(outcoords)):       

                dist = np.sqrt((xlist[i]-xlist[i-1])**2 + (ylist[i]-ylist[i-1])**2)
                #print('distance is {}'.format(dist))
                #print('x-coords',outxmin,outxmax,stampxmin,stampxmax)
                #print('y-coords',outymin,outymax,stampymin,stampymax)
                #inframe[ylist[i]-np.ceil(srcylen/2.):ylist[i]+np.ceil(srcylen/2.),xlist[i]-np.ceil(srcxlen/2.):xlist[i]+np.ceil(srcxlen/2.)] += (source*secperpix*dist)

                #print('inframe shape',inframe[outymin:outymax,outxmin:outxmax].shape)
                #print('source shape',source[stampymin:stampymax,stampxmin:stampxmax].shape)
                
                inframe[outymin:outymax,outxmin:outxmax] += (source[stampymin:stampymax,stampxmin:stampxmax]*secperpix*dist)
        return inframe


    def subsample(self,image,factorx,factory):
        #subsample the input image
        ydim,xdim = image.shape
        substamp = np.zeros((ydim*factory,xdim*factorx))
        
        for i in range(xdim):
            for j in range(ydim):
                #substamp[j*factory:j*(factory+1),i*factorx:i*(factorx+1)] = image[j,i]
                substamp[factory*j:factory*(j+1),factorx*i:factorx*(i+1)] = image[j,i]

        return substamp


    def equidistantXY(self,xstart,ystart,xend,yend,dist): #,ang):
        #return a list of x,y positions that are equidistant
        #between the beginning and ending positions, with 
        #a distance of dist pixels between them
        xlen = 0
        ylen = 0

        deltax = xend - xstart
        deltay = yend - ystart
        ang = np.arctan2(deltay,deltax)
        
        dx = np.cos(ang) * dist
        dy = np.sin(ang) * dist

        if dx != 0.:
            xs = np.arange(xstart,xend+dx/2,dx)
            xlen = len(xs)
        else:
            #motion parallel to y axis
            xs = np.array([xstart])

        if dy != 0:
            ys = np.arange(ystart,yend+dy/2,dy)
            ylen = len(ys)
        else:
            #motion parallel to x asis
            ys = np.array([ystart])

        #make sure lengths agree
        if xlen == 0:
            xs = np.zeros(ylen) + xstart
        if ylen == 0:
            ys = np.zeros(xlen) + ystart
            
        return xs,ys
        

    def radecPerFrame(self,ra0,dec0,ravel,decvel,time):
        #generate a list of RA,Dec locations for all
        #frames
        ra = ra0 + (ravel * time)
        dec = dec0 + (decvel * time)
        return ra,dec

    
    def xyPerFrame(self,velocity,time,ang,x0,y0):

        #rate of movement in x and y (pix/sec)
        ratex = velocity * np.cos(ang)
        ratey = velocity * np.sin(ang)

        #x,y in each frame
        xs = x0 + ratex*time
        ys = y0 + ratey*time
        
        return xs,ys
    



        #return: ramp containing moving target. Also return the x,y position
        #of the target in each frame.
