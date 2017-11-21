#! /usr/bin/env python

'''Create NIRCam PSF files to be used for data simulation (via
Kevin Volk's rampsim.py). For each filter, create PSFs where the
source is centered on a pixel, as well as in steps of 0.25 pixels in
the x and y directions. In addition (??) create sets of PSFs for the 3
wavefront senseing error options that Kevin uses in his script.
'''

import webbpsf
import numpy as np
import os

#offsets = np.arange(-0.5, 0.6, .1)
#wfe (wavefront error) options:
#Units are wavefront error in microns
#NIRCam has two options for each measurement:
#nominal case, and optimistic case, in terms of how well nircam is aligned
#115 - OTE and ISIM intrinsic wfe, optimistic case
#123 - OTE and ISIM intrinsic wfe, nominal case
#132 - above, plus a slight defocus to blur image, to appox image motion. optimistic case
#136 - above, plus a slight defocus to blur image, to appox image motion. nomincal case
#150 - above, plus additional WFE due to NIRCam internal optics. optimistic case
#155 - above, plus additional WFE due to NIRCam internal optics. nominal case

#116 - REVISED requirements value for NIRCam+OTE
#90 - REVISED expected value for NIRCam+OTE



#wfe_list = [0, 123, 136, 155] #for now just do the nominal cases
#filter_list = ['F070W', 'F090W', 'F115W', 'F140M', 'F150W2', 'F150W', 'F162M', 'F164N', 'F182M', 'F187N', 'F200W', 'F210M', 'F212N', 'F250M', 'F277W', 'F300M', 'F322W2', 'F323N', 'F335M', 'F356W', 'F360M', 'F405N', 'F410M', 'F430M', 'F444W', 'F460M', 'F466N', 'F470N', 'F480M']
filter_list = ['F070W', 'F090W', 'F115W', 'F140M', 'F150W2', 'F150W', 'F162M', 'F164N', 'F182M',
               'F187N', 'F200W', 'F210M', 'F212N', 'F250M', 'F277W', 'F300M', 'F322W2', 'F323N',
               'F335M', 'F356W', 'F360M', 'F405N', 'F410M', 'F430M', 'F444W', 'F460M', 'F466N',
               'F470N', 'F480M']
detector_list = ['A1', 'A2', 'A3', 'A4', 'A5', 'B1', 'B2', 'B3', 'B4', 'B5']

#lists for script testing
#offsets = np.arange(-0.1, 0, 0.1)
#wfe_list = [0]
#filter_list = ['F070W']

wfefiles = {}
wfefiles['predicted'] = 'OPD_RevW_ote_for_NIRCam_predicted.fits.gz'
wfefiles['requirements'] = 'OPD_RevW_ote_for_NIRCam_requirements.fits.gz'

pixscales = {}
pixscales['A1'] = (0.0311, 0.0313)
pixscales['A2'] = (0.0308, 0.0309)
pixscales['A3'] = (0.0313, 0.0315)
pixscales['A4'] = (0.0309, 0.0309)
pixscales['A5'] = (0.0628, 0.0631)
pixscales['B1'] = (0.0307, 0.0308)
pixscales['B2'] = (0.0311, 0.0313)
pixscales['B3'] = (0.0308, 0.0309)
pixscales['B4'] = (0.0313, 0.0314)
pixscales['B5'] = (0.0629, 0.0632)

class PSFs:
    def __init__(self):
        self.infile = None

    def make_offset_string(self, offset):
        '''create the offset string that will go into the PSF output file name'''
        #round offset to zero if that's what it is supposed to be but
        #python has kept it not quite zero
        if np.absolute(offset) < 1e-6:
            offset = 0

        if offset < 0:
            first = 'm0p'
        if offset >= 0:
            first = '0p'
        ostr = str(offset)
        dot = ostr.find('.')
        second = ostr[dot+1:]
        return first+second


    def run(self):
        if len(self.xpos_list) != len(self.ypos_list):
            print("WARNING: xpos_list and ypos_list have different")
            print("lengths. They must be the same length.")
            sys.exit()

        if ((max(self.realization_list) > 9) | (min(self.realization_list) < 0)):
            print("WARNING: realization numbers must be between 0 and 9, inclusive.")
            sys.exit()

        offsets_r = []
        offsets_t = []
        offsets_x = []
        offsets_y = []
        #outdir = '/ifs/jwst/wit/witserv/data4/nrc/hilbert/simulated_data/psf_files/'
        #outdir = '/ifs/jwst/wit/nircam/nircam_simulator_data/new_psf_data'
        outdir = './'

        for wfe in self.wfe_list:
            #wfestr = str(wfe)
            #if wfe == 0:
            #    wfestr = 'zero'
            #wfefile = wfebase + wfestr + '.fits'

            wfefile = wfefiles[wfe]

            for filter in self.filter_list:
                fwave = np.int(filter[1:4])
                if fwave < 250:
                    good_detectors = [d for d in self.detector_list if '5' not in d]
                else:
                    good_detectors = [d for d in self.detector_list if '5' in d]

                nc = webbpsf.NIRCam()
                nc.pupilopd = wfefile
                if wfe == 0:
                    nc.pupilopd = None

                #loop over all 10 different statistical
                #realizations of OPD files
                opdname = nc.pupilopd
                if wfe != 0:
                    for i in self.realization_list:

                        # Loop over detector. For the moment, just create PSFs in
                        # the center of the detector. No field dependence yet,
                        # because we need an interpolation scheme.
                        for det in good_detectors:
                            subdir = '/'+det+'/'+filter.lower()+'/'+wfe+'/'
                            if not os.path.exists(outdir+subdir):
                                print("Creating subdirectory {}".format(outdir+subdir))
                                os.makedirs(outdir+subdir)
                            scalex, scaley = pixscales[det]

                            for ox in self.offsets:
                                for oy in self.offsets:
                                    offsetr = np.sqrt((ox*scalex)**2+(oy*scaley)**2)
                                    offsets_r.append(offsetr)
                                    if oy != 0:
                                        theta = np.arctan(np.absolute(ox)/np.absolute(oy))*180/np.pi
                                    if ox == 0:
                                        if oy >= 0:
                                            theta=0
                                        if oy < 0:
                                            theta=180.
                                    if ox < 0:
                                        if oy < 0:
                                            theta = 180.-theta
                                        if oy == 0:
                                            theta = 90.
                                    if ox > 0:
                                        if oy > 0:
                                            theta = 360.-theta
                                        if oy < 0:
                                            theta = 180.+theta
                                        if oy == 0:
                                            theta = 270.
                                    offsets_t.append(theta)
                                    offsets_x.append(ox)
                                    offsets_y.append(oy)
                                    offsett = theta
                                    offsetx = ox
                                    offsety = oy
                                    ox_str = self.make_offset_string(offsetx)
                                    oy_str = self.make_offset_string(offsety)

                                    # Loop over position on detector
                                    for ypos, xpos in zip(self.ypos_list, self.xpos_list):
                                        nc.pupilopd = (opdname, i)
                                        ofilename = "{}{}nircam_{}_x{}_y{}_{}_{}_{}_{}_{}.fits".format(outdir, subdir, det, str(xpos), str(ypos), filter.lower(), wfe, str(i), ox_str, oy_str)
                                        nc.detector = det
                                        nc.detector_position = (ypos, xpos)
                                        nc.options['output_mode'] = 'detector_sampled'
                                        nc.options['source_offset_r'] = offsetr #arcseconds
                                        nc.options['source_offset_theta'] = offsett #degrees ccw from +Y
                                        nc.filter = filter
                                        nc.calcPSF(ofilename, fov_pixels=301)
                                        print(ofilename)
                #else:
                #    ofilename = "{}{}nircam_{}_{}_{}_{}.fits".format(outdir, subdir, filter.lower(), wfestr, ox_str, oy_str)
                #    nc.options['output_mode'] = 'detector_sampled'
                #    nc.options['source_offset_r'] = offsetr #arcseconds
                #    nc.options['source_offset_theta'] = offsett #degrees ccw from +Y
                #    nc.filter = filter
                #    nc.calcPSF(ofilename, fov_pixels=301)


    def add_options(self, parser=None, usage=None):
        if parser == None:
            parser = argparse.ArgumentParser(usage=usage, description="generate PSF library for NIRCam")
        parser.add_argument("offsets", help="list of position offsets, applied to x and y", default= np.arange(-0.5, 0.6, .1))
        #parser.add_argument("wfe_list", help="list of wavefront error values (in microns) to apply", default=[0, 123, 136, 155])
        parser.add_argument("wfe_list", help="list of wavefront error values to apply", default=['predicted', 'requirements'])
        parser.add_argument("filter_list", help="list of filters to use", default=filter_list)
        parser.add_argument("detector_list", help="list of detectors to use", default=detector_list)
        parser.add_argument("xpos_list", help="List of x pixel coordinates of PSF locations", default=[1024])
        parser.add_argument("ypos_list", help="List of y pixel coordinates of PSF locations", default=[1024])
        parser.add_argument("realization_list", help="index numbers of realizations to use", default=np.arange(10))
        return parser


if __name__ == '__main__':
    usagestring = 'USAGE: make_nircam_psfs.py offsets wfelist filterlist'

    psflib = PSFs()
    parser = add_options(usage=usagestring)
    args = parser.parse_args(namespace=psflib)

    psflib.run()
