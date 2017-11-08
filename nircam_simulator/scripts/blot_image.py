#! /usr/bin/env python

'''
Blot an image back to some input WCSs.

Use in conjunction with crop_mosaic.py
'''
import sys, os
from copy import copy
import glob
from astropy.io import fits
from jwst import datamodels
import gwcs
from jwst.outlier_detection import outlier_detection
from jwst.assign_wcs import AssignWcsStep
from jwst.datamodels import container
from . import set_telescope_pointing_separated as stp


class Blot():

    def __init__(self):
        self.detector = ['']
        self.center_ra = [0.]
        self.center_dec = [0.]
        self.pav3 = [0.]
        self.outfile = None

        distdir = '/ifs/jwst/wit/witserv/data4/nrc/hilbert/distortion_reference_file/jwreftools/nircam/'
        self.distfiles = glob.glob(os.path.join(distdir,'NRC*FULL_distortion.asdf'))

        # Pixel scales for the NIRCam detectors
        self.nrc_scale = {'A1':(0.0311,0.0313),
                     'A2':(0.0308,0.0309),
                     'A3':(0.0313,0.0315),
                     'A4':(0.0309,0.0309),
                     'A5':(0.0628,0.0631),
                     'B1':(0.0307,0.0308),
                     'B2':(0.0311,0.0313),
                     'B3':(0.0308,0.0309),
                     'B4':(0.0313,0.0314),
                     'B5':(0.0629,0.0632)}

        self.v2refs = {'A1':120.671376,
                       'A2':120.112090,
                       'A3':51.934455,
                       'A4':52.276773,
                       'A5':86.103458,
                       'B1':-120.968203,
                       'B2':-121.144349,
                       'B3':-53.123823,
                       'B4':-52.818167,
                       'B5':-89.389195}

        self.v3refs = {'A1':-527.387665,
                       'A2':-459.680558,
                       'A3':-527.803414,
                       'A4':-459.809697,
                       'A5':-493.227512,
                       'B1':-457.752680,
                       'B2':-525.458174,
                       'B3':-457.780381,
                       'B4':-525.727255,
                       'B5':-491.443963}

        self.v3yangles = {'A1':-0.573843,
                          'A2':-0.212057,
                          'A3':0.185477,
                          'A4':0.057459,
                          'A5':-0.089504,
                          'B1':0.375008,
                          'B2':0.830079,
                          'B3':-0.485114,
                          'B4':-0.344918,
                          'B5':-0.008844}

    def blot(self):
        
        # Make sure detector, ra, dec, and roll have same number
        # of elements
        if ((len(self.detector) != len(self.center_ra)) | \
                (len(self.detector) != len(self.center_dec)) | \
                (len(self.detector) != len(self.pav3))):
            print('WARNING: detector, center_ra, center_dec')
            print('and pav3 all must have the same number')
            print('of elements.')
            sys.exit()

        if type(self.blotfile) == str:
            input_mod = datamodels.ImageModel(self.blotfile)
            outbase = self.blotfile
            
            # Create a GWCS object from the input file's header
            #input_header = fits.getheader(self.blotfile)
            #transform = gwcs.utils.make_fitswcs_transform(header)
            #input_mod.meta.wcs = gwcs.WCS(transform)

        elif type(self.blotfile) == datamodels.image.ImageModel:
            # Not a great solution at the moment. Save
            # imagemodel instance so that you have a fits
            # file from which the header can be retrieved
            input_mod = copy(self.blotfile)
            self.blotfile.save('temp.fits')
            self.blotfile = 'temp.fits'
            outbase = 'mosaic'
            
        else:
            print('WARNING: unrecognized type for blotfile')
            sys.exit()

        # Create a GWCS object from the input file's header
        input_header = fits.getheader(self.blotfile,ext=1)
        transform = gwcs.utils.make_fitswcs_transform(input_header)
        input_mod.meta.wcs = gwcs.WCS(forward_transform=transform,output_frame='world')
        # Filter and pupil information
        filter = input_mod.meta.instrument.filter
        try:
            pupil = input_mod.meta.instrument.pupil
        except:
            pupil = 'CLEAR'

        # get position angle of input data
        input_pav3 = input_mod.meta.wcsinfo.roll_ref

        # parity is always -1 for nircam
        parity = -1

        # Name of temporary file output for set_telescope_pointing
        # to work on
        shellname = 'temp_wcs_container.fits'
        
        blist = []
        for (det,ra,dec,roll) in \
            zip(self.detector,self.center_ra,\
                self.center_dec,self.pav3):
            # get detector-specific info 
            v2ref,v3ref,v3ang = self.get_siaf_info(det)

            # create datamodel with appropriate metadata
            bmodel = self.make_model(det,ra,dec,v2ref,v3ref,v3ang,parity,input_pav3,filter,pupil)
            #shellname = 'wcs_model_to_blot_to_{}_{}_{}_{}.fits'.format(det,ra,dec,roll)
            bmodel.save(shellname,overwrite=True)
            
            #tmpname = 'wcs_model_to_blot_to_BASEMODEL_{}_{}_{}_{}.fits'.format(det,ra,dec,roll)
            #bmodel.save(tmpname,overwrite=True)
            
            # use set_telescope_pointing to compute local roll
            # angle and PC matrix
            stp.add_wcs(shellname,roll=roll)

            bmodel = datamodels.open(shellname)

            # Now we need to run assign_wcs step so that these
            # files get a gwcs object attached to them.
            dist_reffile = [s for s in self.distfiles if det in s][0]
            bmodel = AssignWcsStep.call(bmodel,override_distortion=dist_reffile)

            #tmpname = 'wcs_model_to_blot_to_ASSIGNWCS_{}_{}_{}_{}.fits'.format(det,ra,dec,roll)
            #bmodel.save(tmpname,overwrite=True)

            # Add to the list of data model instances to blot to
            blist.append(bmodel)
            
        #place the model instances to blot to in a ModelContainer
        blot_list = container.ModelContainer(blist)
        
        #blot the image to each of the WCSs in the blot_list
        pars = {'sinscl':1.0, 'interp':'poly5'}
        reffiles = {}
        mm = outlier_detection.OutlierDetection(blot_list,reffiles=reffiles,**pars)
        blotted_datamodels = outlier_detection.blot_median(input_mod,blot_list,**mm.outlierpars)

        for (bltted,det,ra,dec,roll) in \
            zip(blotted_datamodels,self.detector,self.center_ra,\
                self.center_dec,self.pav3):
            if self.outfile is None:
                self.outfile = 'blotted_from_{}_to_{}_{}_{}_{}.fits'.format(outbase,det,ra,dec,roll)
            bltted.save(self.outfile)


    def get_siaf_info(self,detname):
        # get v2,v3 reference values and y3yangle for a given detector
        vv = self.v2refs[detname]
        vvv = self.v3refs[detname]
        ang = self.v3yangles[detname]
        return vv,vvv,ang


    def make_model(self,detname,raval,decval,v2,v3,v3angle,vpar,rollval,filt,pup):
        #define the WCS of the blotted output
        blot_to = datamodels.ImageModel((2048,2048))
        blot_to.meta.wcsinfo.cdelt1 = self.nrc_scale[detname][0] / 3600.
        blot_to.meta.wcsinfo.cdelt2 = self.nrc_scale[detname][1] / 3600.
        blot_to.meta.wcsinfo.crpix1 = 1024.5
        blot_to.meta.wcsinfo.crpix2 = 1024.5
        blot_to.meta.wcsinfo.crval1 = raval
        blot_to.meta.wcsinfo.crval2 = decval
        blot_to.meta.wcsinfo.ctype1 = 'RA---TAN'
        blot_to.meta.wcsinfo.ctype2 = 'DEC--TAN'
        blot_to.meta.wcsinfo.cunit1 = 'deg'
        blot_to.meta.wcsinfo.cunit2 = 'deg'
        blot_to.meta.wcsinfo.ra_ref = raval
        blot_to.meta.wcsinfo.dec_ref = decval
        blot_to.meta.wcsinfo.roll_ref = rollval 
        blot_to.meta.wcsinfo.wcsaxes = 2
        blot_to.meta.wcsinfo.v2_ref = v2
        blot_to.meta.wcsinfo.v3_ref = v3
        blot_to.meta.wcsinfo.v3yangle = v3angle
        blot_to.meta.wcsinfo.vparity = vpar
        blot_to.meta.instrument.channel = 'SHORT'
        blot_to.meta.instrument.detector = 'NRC'+detname
        blot_to.meta.instrument.filter = filt
        blot_to.meta.instrument.module = detname[0]
        blot_to.meta.instrument.name = 'NIRCAM'
        blot_to.meta.instrument.pupil = pup
        blot_to.meta.telescope = 'JWST'
        blot_to.meta.exposure.start_time = 57410.24546415885
        blot_to.meta.exposure.end_time = 57410.2477009838
        blot_to.meta.exposure.type = 'NRC_IMAGE'
        blot_to.meta.target.ra = raval
        blot_to.meta.target.dec = decval
        return blot_to


    def add_options(self,parser=None,usage=None):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage,description='Extract SCA-sized area from moasic')
        parser.add_argument("blotfile",help="Filename or model instance name of fits file containing mosaic.")
        #parser.add_argument("--blotmodel",help="Datamodel containing array to blot",default=None)
        parser.add_argument("--detector",help="NIRCam detectors to blot to. Multiple inputs ok. (e.g. A1 A5)",nargs='*') 
        parser.add_argument("center_ra",help="RA at the center of the extracted area",type=np.float,nargs='*')
        parser.add_argument("center_dec",help="Dec at the center of the extracted area",type=np.float,nargs='*')
        parser.add_argument("pav3",help="Position angle for outputs to use when blotting",nargs='*')
        parser.add_argument("--outfile",help="Name of output fits file containing extracted image",default=None,nargs='*')
        return parser


if __name__ == '__main__':
    usagestring = 'USAGE: blot_image.py filename.fits --detector A1 A1 --center_ra 10.2 10.2001 --center_dec 12.9 12.91 --local_roll 0 0'

    b = Blot()
    parser = b.add_options(usage = usagestring)
    args = parser.parse_args(namespace=b)
    b.blot()

