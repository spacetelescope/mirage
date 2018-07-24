#! /usr/bin env python

'''
Reorganization of ramp simulator code. This class is used to construct
a "seed image" from a point source catalog. This seed image is a
noiseless countrate image containing only sources. No noise, no
cosmic rays.
'''


import argparse, sys, glob, os
import copy
import pkg_resources
import scipy.signal as s1
import numpy as np
import math
from photutils import detect_threshold, detect_sources
from astropy.io import fits, ascii
from astropy.table import Table, Column
from astropy.modeling.models import Shift, Sersic2D, Polynomial2D, Mapping
from astropy.convolution import convolve
from asdf import AsdfFile
import yaml
import time

from mirage.scripts import rotations  # this is Colin's set of rotation functions
from mirage.scripts import polynomial # more of Colin's functions
from mirage.scripts import read_siaf_table
from mirage.scripts import set_telescope_pointing_separated as set_telescope_pointing
from mirage.scripts import moving_targets
from mirage.scripts import segmentation_map as segmap
from mirage.scripts.psf_generator import PSF

INST_LIST = ['nircam', 'niriss', 'fgs']
MODES = {'nircam': ["imaging", "ts_imaging", "wfss", "ts_wfss"],
         'niriss': ["imaging"],
         'fgs': ["imaging"]}  
TRACKING_LIST = ['sidereal','non-sidereal']
inst_abbrev = {'nircam': 'NRC',
               'niriss': 'NIS',
               'fgs': 'FGS'}
PIXELSCALE = {'nircam': {'sw':0.031, 'lw':0.063},
              'niriss': 0.065,
              'fgs': 0.065}
FULL_ARRAY_SIZE = {'nircam': 2048,
                   'niriss': 2048,
                   'fgs':2048}
ALLOWEDOUTPUTFORMATS = ['DMS']
WFE_OPTIONS = ['predicted', 'requirements']
WFEGROUP_OPTIONS = np.arange(5)

class Catalog_seed():
    def __init__(self):
        # Locate the module files, so that we know where to look
        # for config subdirectory
        self.modpath = pkg_resources.resource_filename('mirage', '')

        # Get the location of the MIRAGE_DATA environment
        # variable, so we know where to look for darks, CR,
        # PSF files, etc later
        self.env_var = 'MIRAGE_DATA'
        datadir = os.environ.get(self.env_var)

        if datadir is None:
            raise ValueError(("WARNING: {} environment variable is not set."
                              "This must be set to the base directory"
                              "containing the darks, cosmic ray, PSF, etc"
                              "input files needed for the simulation."
                              "These files must be downloaded separately"
                              "from the Mirage package.".format(self.env_var)))

        # if a grism signal rate image is requested, expand
        # the width and height of the signal rate image by this
        # factor, so that the grism simulation software can
        # track sources that are outside the requested subarray
        # in order to calculate contamination.
        self.grism_direct_factor = np.sqrt(2.)

        # self.coord_adjust contains the factor by which the
        # nominal output array size needs to be increased
        # (used for WFSS mode), as well as the coordinate
        # offset between the nominal output array coordinates,
        # and those of the expanded array. These are needed
        # mostly for WFSS observations, where the nominal output
        # array will not sit centered in the expanded output image.
        self.coord_adjust = {'x':1., 'xoffset':0., 'y':1., 'yoffset':0.}

        # NIRCam rough noise values. Used to make educated guesses when
        # creating segmentation maps
        self.single_ron = 6. # e-/read
        self.grism_background = 0.25 # e-/sec

        # keep track of the maximum index number of sources
        # in the various catalogs. We don't want to end up
        # with multiple sources having the same index numbers
        self.maxindex = 0

    def make_seed(self):
        """MAIN FUNCTION"""
        # Read in input parameters and quality check
        self.readParameterFile()
        self.fullPaths()
        self.filecheck()
        self.basename = os.path.join(self.params['Output']['directory'],
                                     self.params['Output']['file'][0:-5].split('/')[-1])
        self.params['Output']['file'] = self.basename + self.params['Output']['file'][-5:]

        self.readSubarrayDefinitionFile()
        self.checkParams()
        self.getSubarrayBounds()
        self.instrument_specific_dicts(self.params['Inst']['instrument'].lower())

        # If the output is a direct image to be dispersed, expand the size
        # of the nominal FOV so the disperser can account for sources just
        # outside whose traces will fall into the FOV
        if self.params['Output']['grism_source_image']:
            self.calcCoordAdjust()

        #image dimensions
        self.nominal_dims = np.array([self.subarray_bounds[3] - self.subarray_bounds[1] + 1,
                                      self.subarray_bounds[2] - self.subarray_bounds[0] + 1])
        self.output_dims = (self.nominal_dims * np.array([self.coord_adjust['y'],
                                                          self.coord_adjust['x']])).astype(np.int)

        # calculate the exposure time of a single frame, based on the size of the subarray
        self.calcFrameTime()

        # Read in the pixel area map, which will be needed for certain
        # sources in the seed image
        self.prepare_PAM()
        
        # For imaging mode, generate the countrate image using the catalogs
        if self.params['Telescope']['tracking'].lower() != 'non-sidereal':
            print('Creating signal rate image of synthetic inputs.')
            self.seedimage, self.seed_segmap = self.addedSignals()
            outapp = ''

        # If we are tracking a non-sidereal target, then
        # everything in the catalogs needs to be streaked across
        # the detector
        if self.params['Telescope']['tracking'].lower() == 'non-sidereal':
            print('Creating signal ramp of synthetic inputs')
            self.seedimage, self.seed_segmap = self.non_sidereal_seed()
            outapp = '_nonsidereal_target'

        # If non-sidereal targets are requested (KBOs, asteroids, etc,
        # create a RAPID integration which includes those targets
        mov_targs_ramps = []
        if (self.runStep['movingTargets'] | self.runStep['movingTargetsSersic']
            | self.runStep['movingTargetsExtended']):
            print(("Creating signal ramp of sources that are moving with "
                   "respect to telescope tracking."))
            trailed_ramp, trailed_segmap = self.make_trailed_ramp()
            outapp += '_trailed_sources'

            # Now we need to expand frameimage into a ramp
            # so we can add the trailed objects
            print('Combining trailed object ramp with that containing tracked targets')
            if self.params['Telescope']['tracking'].lower() != 'non-sidereal':
                self.seedimage = self.combineSimulatedDataSources('countrate', self.seedimage, trailed_ramp)
            else:
                self.seedimage = self.combineSimulatedDataSources('ramp', self.seedimage, trailed_ramp)
            self.seed_segmap += trailed_segmap

        # For seed images to be dispersed in WFSS mode,
        # embed the seed image in a full frame array. The disperser
        # tool does not work on subarrays
        if ((self.params['Inst']['mode'] in ['wfss','ts_wfss']) & \
            ('FULL' not in self.params['Readout']['array_name'])):
            self.seedimage, self.seed_segmap = self.pad_wfss_subarray(self.seedimage, self.seed_segmap)

        # Save the combined static + moving targets ramp
        self.saveSeedImage()
        # Return info in a tuple
        # return (self.seedimage, self.seed_segmap, self.seedinfo)

    def add_detector_to_zeropoints(self, detector):
        """Manually add detector dependence to the zeropoint table for
        NIRCam and NIRISS simualtions. This is being done as a placeholder
        for the future, where we expect zeropoints to be detector-dependent.

        Parameters:
        -----------
        detector : str
            Name of detector to add to the table

        Returns:
        --------
        Nothing
        """
        # Add "Detector" to the list of column names
        base_table = copy.deepcopy(self.zps)
        num_entries = len(self.zps)
        det_column = Column(np.repeat(detector, num_entries), name="Detector")
        base_table.add_column(det_column, index = 0)
        return base_table

    def prepare_PAM(self):
        """
        Read in and prepare the pixel area map (PAM), to be used
        during the seed creation process.

        Parameters:
        -----------
        None

        Returns:
        --------
        None
        """
        fname = self.params['Reffiles']['pixelAreaMap']

        # Read in PAM
        try:
            pam, header = fits.getdata(fname, header=True)
        except:
            raise IOError('WARNING: unable to read in {}'.format(fname))

        # Crop to expected subarray
        try:
            pam = pam[self.subarray_bounds[1]:self.subarray_bounds[3]+1,
                      self.subarray_bounds[0]:self.subarray_bounds[2]+1]
        except:
            raise ValueError("Unable to crop pixel area map to expected subarray.")

        # If we are making a grism direct image, we need to embed the true pixel area
        # map in an array of the appropriate dimension, where any pixels outside the
        # actual aperture are set to 1.0
        if self.params['Output']['grism_source_image']:
            mapshape = pam.shape
            #cannot use this: g, yd, xd = signalramp.shape
            #need to update dimensions: self.pam = np.ones((yd, xd))
            self.pam = np.ones(self.output_dims)
            ys = self.coord_adjust['yoffset']
            xs = self.coord_adjust['xoffset']
            self.pam[ys:ys+mapshape[0], xs:xs+mapshape[1]] = np.copy(pam)
        else:
            self.pam = pam

    def pad_wfss_subarray(self, seed, seg):
        """
        WFSS seed images that are to go into the disperser
        must be full frame (or larger). The disperser cannot 
        work on subarray images. So embed the subarray seed image
        in a full-frame sized array.

        Parameters:
        -----------
        None

        Returns:
        --------
        Padded seed image and segmentation map
        """
        seeddim = seed.shape
        nx = np.int(2048 * self.grism_direct_factor)
        ffextra = np.int((nx - 2048) / 2)
        bounds = np.array(self.subarray_bounds)

        subextrax = np.int((seeddim[-1] - bounds[2]) / 2)
        subextray = np.int((seeddim[-2] - bounds[3]) / 2)

        extradiffy = ffextra - subextray
        extradiffx = ffextra - subextrax

        #print(self.subarray_bounds)
        #print(ffextra)
        #print(subextrax, subextray)
        #print(extradiffx, extradiffy)
        #print(seeddim)
        #print(nx)
        #print(extradiffx, extradiffy, extradiffx+seeddim[-1]-1, extradiffy+seeddim[-2]-1)

        exbounds = [extradiffx, extradiffy, extradiffx+seeddim[-1]-1, extradiffy+seeddim[-2]-1]

        if len(seeddim) == 2:
            padded_seed = np.zeros((nx, nx))
            padded_seed[exbounds[1]:exbounds[3] + 1, exbounds[0]:exbounds[2] + 1] = seed
            padded_seg = np.zeros((nx, nx), dtype=np.int)
            padded_seg[exbounds[1]:exbounds[3] + 1, exbounds[0]:exbounds[2] + 1] = seg
        elif len(seeddim) == 4:
            padded_seed = np.zeros((seeddim[0], seeddim[1], nx, nx))
            padded_seed[:, :, exbounds[1]:exbounds[3] + 1, exbounds[0]:exbounds[2] + 1] = seed
            padded_seg = np.zeros((seeddim[0], seeddim[1], nx, nx), dtype=np.int)
            padded_seg[:, :, exbounds[1]:exbounds[3] + 1, exbounds[0]:exbounds[2] + 1] = seg
        else:
            raise ValueError("Seed image is not 2D or 4D. It should be.")
        return padded_seed, padded_seg
        
    def saveSeedImage(self):
        # Create the grism direct image or ramp to be saved
        arrayshape = self.seedimage.shape
        if len(arrayshape) == 2:
            units = 'ADU/sec'
            yd, xd = arrayshape
            tgroup = 0.
            print('Seed image is 2D.')
        elif len(arrayshape) == 3:
            units = 'ADU'
            g, yd, xd = arrayshape
            tgroup = self.frametime * (self.params['Readout']['nframe'] + self.params['Readout']['nskip'])
            print('Seed image is 3D.')
        elif len(arrayshape) == 4:
            units = 'ADU'
            integ, g, yd, xd = arrayshape
            tgroup = self.frametime * (self.params['Readout']['nframe'] + self.params['Readout']['nskip'])
            print('Seed image is 4D.')

        xcent_fov = xd / 2
        ycent_fov = yd / 2
        
        kw = {}
        kw['xcenter'] = xcent_fov
        kw['ycenter'] = ycent_fov
        kw['units'] = units
        kw['TGROUP'] = tgroup
        if self.params['Readout']['pupil'][0].upper() == 'F':
            usefilt = 'pupil'
        else:
            usefilt = 'filter'

        self.seed_file = os.path.join(self.basename + '_' + self.params['Readout'][usefilt] + '_seed_image.fits')
            
        # Set FGS filter to "N/A" in the output file
        # as this is the value DMS looks for.
        if self.params['Readout'][usefilt] == "NA":
            self.params['Readout'][usefilt] = "N/A"   
        kw['filter'] = self.params['Readout'][usefilt]
        kw['PHOTFLAM'] = self.photflam
        kw['PHOTFNU'] = self.photfnu
        kw['PHOTPLAM'] = self.pivot * 1.e4 # put into angstroms
        kw['NOMXDIM'] = self.nominal_dims[1]
        kw['NOMYDIM'] = self.nominal_dims[0]
        kw['NOMXSTRT'] = self.coord_adjust['xoffset'] + 1
        kw['NOMXEND'] = self.nominal_dims[1] + self.coord_adjust['xoffset']
        kw['NOMYSTRT'] = self.coord_adjust['yoffset'] + 1
        kw['NOMYEND'] = self.nominal_dims[0] + self.coord_adjust['yoffset']

        # Seed images provided to disperser are always embedded in an array
        # with dimensions equal to full frame * self.grism_direct_factor
        if self.params['Inst']['mode'] in ['wfss', 'ts_wfss']:
            kw['NOMXDIM'] = self.ffsize
            kw['NOMYDIM'] = self.ffsize
            kw['NOMXSTRT'] = np.int(self.ffsize * (self.grism_direct_factor - 1) / 2.)
            kw['NOMXEND'] = kw['NOMXSTRT'] + self.ffsize - 1
            kw['NOMYSTRT'] = np.int(self.ffsize * (self.grism_direct_factor - 1) / 2.)
            kw['NOMYEND'] = kw['NOMYSTRT'] + self.ffsize - 1
        
        kw['GRISMPAD'] = self.grism_direct_factor
        self.seedinfo = kw
        self.saveSingleFits(self.seedimage, self.seed_file, key_dict=kw, image2=self.seed_segmap, image2type='SEGMAP')
        print("Seed image and segmentation map saved as {}".format(self.seed_file))
        print("Seed image, segmentation map, and metadata available as:")
        print("self.seedimage, self.seed_segmap, self.seedinfo.")


    def fullPaths(self):
        # Expand all input paths to be full paths
        # This is to allow for easier Condor-ization of
        # many runs
        pathdict = {'Reffiles':['dark', 'linearized_darkfile', 'superbias',
                                'subarray_defs', 'linearity',
                                'saturation', 'gain', 'pixelflat',
                                'illumflat', 'astrometric', 'distortion_coeffs', 'ipc',
                                'crosstalk', 'occult', 'pixelAreaMap',
                                'flux_cal', 'readpattdefs', 'filter_throughput'],
                    'simSignals':['pointsource', 'psfpath', 'galaxyListFile', 'extended',
                                  'movingTargetList', 'movingTargetSersic',
                                  'movingTargetExtended', 'movingTargetToTrack'],
                    'Output':['file', 'directory']}

        all_config_files = {'nircam': {'Reffiles-subarray_defs': 'NIRCam_subarray_definitions.list',
                                       'Reffiles-flux_cal': 'NIRCam_zeropoints.list',
                                       'Reffiles-crosstalk': 'xtalk20150303g0.errorcut.txt',
                                       'Reffiles-readpattdefs': 'nircam_read_pattern_definitions.list',
                                       'Reffiles-filter_throughput': 'placeholder.txt'},
                            'niriss': {'Reffiles-subarray_defs': 'niriss_subarrays.list',
                                       'Reffiles-flux_cal': 'niriss_zeropoints.list',
                                       'Reffiles-crosstalk': 'niriss_xtalk_zeros.txt',
                                       'Reffiles-readpattdefs': 'niriss_readout_pattern.txt',
                                       'Reffiles-filter_throughput': 'placeholder.txt'},
                            'fgs': {'Reffiles-subarray_defs': 'guider_subarrays.list',
                                    'Reffiles-flux_cal': 'guider_zeropoints.list',
                                    'Reffiles-crosstalk': 'guider_xtalk_zeros.txt',
                                    'Reffiles-readpattdefs': 'guider_readout_pattern.txt',
                                    'Reffiles-filter_throughput': 'placeholder.txt'}}
        config_files = all_config_files[self.params['Inst']['instrument'].lower()]

        for key1 in pathdict:
            for key2 in pathdict[key1]:
                if self.params[key1][key2].lower() not in ['none', 'config']:
                    self.params[key1][key2] = os.path.abspath(os.path.expandvars(self.params[key1][key2]))
                elif self.params[key1][key2].lower() == 'config':
                    cfile = config_files['{}-{}'.format(key1, key2)]
                    fpath = os.path.join(self.modpath, 'config', cfile)
                    self.params[key1][key2] = fpath
                    print("'config' specified: Using {} for {}:{} input file".format(fpath, key1, key2))
                    
    def mag_to_countrate(self, magsys, mag, photfnu=None, photflam=None):
        # Convert object magnitude to counts/sec
        if magsys.lower() == 'abmag':
            try:
                return 10**((mag + 48.6) / -2.5) / photfnu
            except:
                raise ValueError(("AB mag to countrate conversion failed."
                                  "magnitude = {}, photfnu = {}".format(mag, photfnu)))
        if magsys.lower() == 'vegamag':
            try:
                return 10**((self.vegazeropoint - mag) / 2.5)
            except:
                raise ValueError(("Vega mag to countrate conversion failed."
                                  "magnitude = {}".format(mag)))
        if magsys.lower() == 'stmag':
            try:
                return 10**((mag + 21.1) / -2.5) / photflam
            except:
                raise ValueError(("ST mag to countrate conversion failed."
                                  "magnitude = {}, photflam = {}".format(mag, photflam)))

    def combineSimulatedDataSources(self, inputtype, input1, mov_tar_ramp):
        """Combine the exposure containing the trailed sources with the
        countrate image containing the static sources
        inputtype can be 'countrate' in which case input needs to be made
        into a ramp before combining with mov_tar_ramp, or 'ramp' in which
        case you can combine directly. Use 'ramp' with
        non-sidereal TRACKING data, and 'countrate' with sidereal TRACKING data
        """
        if inputtype == 'countrate':
            # First change the countrate image into a ramp
            yd, xd = input1.shape
            numints = self.params['Readout']['nint']
            num_frames = self.params['Readout']['ngroup'] * \
                         (self.params['Readout']['nframe'] + self.params['Readout']['nskip'])
            print("Countrate image of synthetic signals being converted to "
                  "RAPID/NISRAPID integration with {} frames.".format(num_frames))
            input1_ramp = np.zeros((numints, num_frames, yd, xd))
            for i in range(num_frames):
                input1_ramp[0, i, :, :] = input1 * self.frametime * (i + 1)
            if numints > 1:
                for integ in range(1, numints):
                    input1_ramp[integ, :, :, :] = input1_ramp[0, :, :, :]

        else:
            # If input1 is a ramp rather than a countrate image
            input1_ramp = input1

        # Combine the input1 ramp and the moving target ramp, which are
        # now both RAPID mode
        totalinput = input1_ramp + mov_tar_ramp
        return totalinput

    def make_trailed_ramp(self):
        # Create a ramp for objects that are trailing through
        # the field of view during the integration
        mov_targs_ramps = []
        mov_targs_segmap = None

        if self.params['Telescope']['tracking'].lower() != 'non-sidereal':
            tracking = False
            ra_vel = None
            dec_vel = None
        else:
            tracking = True
            ra_vel = self.ra_vel
            dec_vel = self.dec_vel
            # print("Moving target mode, creating trailed object with")
            # print("RA velocity of {} and dec_val of {}".format(ra_vel, dec_vel))

        if self.runStep['movingTargets']:
            #print('Starting moving targets for point sources!')
            mov_targs_ptsrc, mt_ptsrc_segmap = self.movingTargetInputs(self.params['simSignals']['movingTargetList'],
                                                                       'pointSource',
                                                                       MT_tracking=tracking,
                                                                       tracking_ra_vel=ra_vel,
                                                                       tracking_dec_vel=dec_vel)
            # Multiply by pixel area map since these sources are trailed across detector
            mov_targs_ptsrc *= self.pam
            mov_targs_ramps.append(mov_targs_ptsrc)
            #print("Moving target segmap, min, max {}, {}".format(np.min(mt_ptsrc_segmap), np.max(mt_ptsrc_segmap)))
            mov_targs_segmap = np.copy(mt_ptsrc_segmap)

        # moving target using a sersic object
        if self.runStep['movingTargetsSersic']:
            #print("Moving targets, sersic!")
            mov_targs_sersic, mt_galaxy_segmap = self.movingTargetInputs(self.params['simSignals']['movingTargetSersic'],
                                                                         'galaxies',
                                                                         MT_tracking=tracking,
                                                                         tracking_ra_vel=ra_vel,
                                                                         tracking_dec_vel=dec_vel)
            # Multiply by pixel area map 
            mov_targs_sersic *= self.pam
            mov_targs_ramps.append(mov_targs_sersic)
            if mov_targs_segmap is None:
                mov_targs_segmap = np.copy(mt_galaxy_segmap)
            else:
                mov_targs_segmap += mt_galaxy_segmap

        # moving target using an extended object
        if self.runStep['movingTargetsExtended']:
            #print("Extended moving targets!!!")
            mov_targs_ext, mt_ext_segmap = self.movingTargetInputs(self.params['simSignals']['movingTargetExtended'],
                                                                   'extended',
                                                                   MT_tracking=tracking,
                                                                   tracking_ra_vel=ra_vel,
                                                                   tracking_dec_vel=dec_val)
            # Multiply by pixel area map 
            mov_targs_ext *= self.pam
            mov_targs_ramps.append(mov_targs_ext)
            if mov_targs_segmap is None:
                mov_targs_segmap = np.copy(mt_ext_segmap)
            else:
                mov_targs_segmap += mt_ext_segmap

        mov_targs_integration = None
        if self.runStep['movingTargets'] or self.runStep['movingTargetsSersic'] or self.runStep['movingTargetsExtended']:
            # Combine the ramps of the moving targets if there is more than one type
            mov_targs_integration = mov_targs_ramps[0]
            if len(mov_targs_ramps) > 1:
                for i in range(1, len(mov_targs_ramps)):
                    mov_targs_integration += mov_targs_ramps[0]
        return mov_targs_integration, mov_targs_segmap

    def calcFrameTime(self):
        # calculate the exposure time of a single frame of the proposed output ramp
        # based on the size of the cropped dark current integration
        # numint, numgrp, yd, xd = self.dark.data.shape
        yd, xd = self.nominal_dims
        # self.frametime = (xd/self.params['Readout']['namp'] + 12.) * (yd + 1) * 10.00 * 1.e-6
        # UPDATED VERSION, 16 Sept 2017
        if 'nircam' in self.params['Inst']['instrument'].lower():
            colpad = 12
            rowpad = 2
            if ((xd <= 8) & (yd <= 8)):
                rowpad = 3
            self.frametime = ((1.0 * xd / self.params['Readout']['namp'] + colpad) * (yd + rowpad)) * 1.e-5
        elif self.params['Inst']['instrument'].lower() in ['niriss', 'fgs']:
            # the following applies to NIRISS and Guider full frame imaging and
            # NIRISS sub-arrays.
            #
            # According JDox the NIRCam full frame time is 10.73677 seconds the
            # same as for NIRISS, but right now the change does not apply to NIRCam.
            #
            #
            # note that the Guider frame time may be different for small sub-arrays
            # less than 64 pixels square, but that needs to be confirmed.
            colpad = 12
            if self.params['Readout']['namp'] == 4:
                pad1 = 1
                pad2 = 1
            else:
                pad1 = 2
                pad2 = 0
            self.frametime = (pad2 + (yd / self.params['Readout']['namp'] + colpad) 
                              * (xd + pad1)) * 0.00001

    def calcCoordAdjust(self):
        # Calculate the factors by which to expand the output array size, as well as the coordinate
        # offsets between the nominal output array and the input lists if the observation being
        # modeled is wfss

        dtor = math.radians(1.)

        # Normal imaging with grism image requested
        if self.params['Output']['grism_source_image']:
            self.coord_adjust['x'] = self.grism_direct_factor
            self.coord_adjust['y'] = self.grism_direct_factor
            self.coord_adjust['xoffset'] = np.int((self.grism_direct_factor - 1.) * (self.subarray_bounds[2] - self.subarray_bounds[0] + 1) / 2.)
            self.coord_adjust['yoffset'] = np.int((self.grism_direct_factor - 1.) * (self.subarray_bounds[3] - self.subarray_bounds[1] + 1) / 2.)

    def non_sidereal_seed(self):
        """
        Create a seed EXPOSURE in the case where the instrument is tracking
        a non-sidereal target
        """

        # Create a count rate image containing only the non-sidereal target(s)
        # These will be stationary in the fov
        nonsidereal_countrate, nonsidereal_segmap, self.ra_vel, self.dec_vel, vel_flag \
            = self.nonsidereal_CRImage(self.params['simSignals']['movingTargetToTrack'])

        # Expand into a RAPID exposure and convert from signal rate to signals
        ns_yd, ns_xd = nonsidereal_countrate.shape
        ns_int = self.params['Readout']['nint']
        ns_group = self.params['Readout']['ngroup']
        ns_nframe = self.params['Readout']['nframe']
        ns_nskip = self.params['Readout']['nskip']
        totframes = ns_group * (ns_nframe + ns_nskip)
        tmptimes = self.frametime * np.arange(1, totframes + 1)

        #non_sidereal_ramp = np.zeros((totframes, ns_yd, ns_xd))
        non_sidereal_ramp = np.zeros((ns_int, ns_group, ns_yd, ns_xd))
        for i in range(totframes):
            for integ in range(ns_int):
                non_sidereal_ramp[integ, i, :, :] = nonsidereal_countrate * tmptimes[i]

        # Now we need to collect all the other sources (point sources,
        # galaxies, extended) in the other input files, and treat them
        # as targets which will move across the field of view during
        # the exposure.
        mtt_data_list = []
        mtt_data_segmap = None
        #mtt_zero_list = []

        if self.runStep['pointsource']:
            # Now ptsrc is a list, which we need to provide to
            # movingTargetInputs
            mtt_ptsrc, mtt_ptsrc_segmap = self.movingTargetInputs(self.params['simSignals']['pointsource'],
                                                                  'pointSource',
                                                                  MT_tracking=True,
                                                                  tracking_ra_vel=self.ra_vel,
                                                                  tracking_dec_vel=self.dec_vel,
                                                                  trackingPixVelFlag=vel_flag)
            # Multiply by pixel area map since these sources are trailed
            # across detector
            mtt_ptsrc *= self.pam
            mtt_data_list.append(mtt_ptsrc)
            if mtt_data_segmap is None:
                mtt_data_segmap = np.copy(mtt_ptsrc_segmap)
            else:
                mtt_data_segmap += mtt_ptsrc_segmap
            print(("Done with creating moving targets from {}"
                   .format(self.params['simSignals']['pointsource'])))

        if self.runStep['galaxies']:
            mtt_galaxies, mtt_galaxies_segmap = self.movingTargetInputs(self.params['simSignals']['galaxyListFile'],
                                                                        'galaxies',
                                                                        MT_tracking=True,
                                                                        tracking_ra_vel=self.ra_vel,
                                                                        tracking_dec_vel=self.dec_vel,
                                                                        trackingPixVelFlag=vel_flag)
            # Multiply by pixel area map 
            mtt_galaxies *= self.pam
            mtt_data_list.append(mtt_galaxies)
            if mtt_data_segmap is None:
                mtt_data_segmap = np.copy(mtt_galaxies_segmap)
            else:
                mtt_data_segmap += mtt_galaxies_segmap

            print(("Done with creating moving targets from {}".
                   format(self.params['simSignals']['galaxyListFile'])))

        if self.runStep['extendedsource']:
            mtt_ext, mtt_ext_segmap = self.movingTargetInputs(self.params['simSignals']['extended'],
                                                              'extended',
                                                              MT_tracking=True,
                                                              tracking_ra_vel=self.ra_vel,
                                                              tracking_dec_vel=self.dec_vel,
                                                              trackingPixVelFlag=vel_flag)
            # Multiply by pixel area map 
            mtt_ext *= self.pam
            mtt_data_list.append(mtt_ext)
            if mtt_data_segmap is None:
                mtt_data_segmap = np.copy(mtt_ext_segmap)
            else:
                mtt_data_segmap += mtt_ext_segmap

            print(("Done with creating moving targets from {}".
                   format(self.params['simSignals']['extended'])))

        # Add in the other objects which are not being tracked on
        # (i.e. the sidereal targets)
        if len(mtt_data_list) > 0:
            for i in range(len(mtt_data_list)):
                non_sidereal_ramp += mtt_data_list[i]
                # non_sidereal_zero += mtt_zero_list[i]
        if mtt_data_segmap is not None:
            nonsidereal_segmap += mtt_data_segmap
        return non_sidereal_ramp, nonsidereal_segmap


    def readMTFile(self, file):
        """
        Read in moving target list file

        Arguments:
        ----------
        file : str
            name of moving target catalog file

        Returns:
        --------
        returns : obj
            Table containing moving target entries
            pixelflag (boolean) -- If true, locations are in units of
                pixels. If false, locations are RA, Dec
            pixelvelflag (boolean) -- If true, moving target velocities
                are in units of pixels/hour. If false, arcsec/hour
            magsys -- magnitude system of the moving target magnitudes
        """
        mtlist = ascii.read(file, comment='#')

        # Convert all relevant columns to floats
        for col in mtlist.colnames:
            if mtlist[col].dtype in ['int64', 'int']:
                mtlist[col] = mtlist[col].data * 1.

        # Check to see whether the position is in x,y or ra,dec
        pixelflag = False
        try:
            if 'position_pixels' in mtlist.meta['comments'][0:4]:
                pixelflag = True
        except:
            pass

        # If present, check whether the velocity entries are pix/sec
        # or arcsec/sec.
        pixelvelflag = False
        try:
            if 'velocity_pixels' in mtlist.meta['comments'][0:4]:
                pixelvelflag = True
        except:
            pass

        # If present, check whether the radius entries (for galaxies)
        # are in arcsec or pixels. If in arcsec, change to pixels
        if 'radius' in mtlist.colnames:
            if 'radius_pixels' not in mtlist.meta['comments'][0:4]:
                mtlist['radius'] /= self.pixscale[0]

        # If galaxies are present, change position angle from degrees
        # to radians
        if 'pos_angle' in mtlist.colnames:
            mtlist['pos_angle'] = mtlist['pos_angle'] * np.pi / 180.

        # Check to see if magnitude system is specified in comments
        # If not, assume AB mags
        msys = 'abmag'

        condition=('stmag' in mtlist.meta['comments'][0:4]) | ('vegamag' in mtlist.meta['comments'][0:4])
        if condition:
            msys = [l for l in mtlist.meta['comments'][0:4] if 'mag' in l][0]
            msys = msys.lower()

        return mtlist, pixelflag, pixelvelflag, msys.lower()

    def movingTargetInputs(self, file, input_type, MT_tracking=False,
                           tracking_ra_vel=None, tracking_dec_vel=None,
                           trackingPixVelFlag=False):
        """Read in listfile of moving targets and perform needed
        calculations to get inputs for moving_targets.py

        input_type can be 'pointSource','galaxies', or 'extended'
        """
        # Read input file - should be able to use for all modes
        mtlist, pixelFlag, pixvelflag, magsys = self.readMTFile(file)

        # If the input catalog has an index column
        # use that, otherwise add one
        if 'index' in mtlist.colnames:
            indexes = mtlist['index']
        else:
            indexes = np.arange(1, len(mtlist['x_or_RA']) + 1)
        # Make sure there is no 0th object
        if np.min(indexes) == 0:
            indexes += 1
        # Make sure the index numbers don't overlap with any
        # sources already present. Increment the maxindex
        # value.
        if np.min(indexes) <= self.maxindex:
            indexes += self.maxindex
        self.maxindex = np.max(indexes)

        if MT_tracking == True:
            # Here, we are tracking a non-sidereal target.
            try:
                # If there are moving targets on top of the non-
                # sidereal tracking (e.g. tracking Io but Europa
                # comes into the fov), then the velocity vector
                # of the moving target needs to be adjusted.
                # If the input catalog already contains
                # 'x_or_RA_velocity' then we know we have a moving
                # target. If it doesn't, then we have sidereal
                # targets, and we can simply set their velocity
                # as the inverse of that being tracked.
                mtlist['x_or_RA_velocity'] -= tracking_ra_vel #* (1./365.25/24.)
                mtlist['y_or_Dec_velocity'] -= tracking_dec_vel #* (1./365.25/24.)
                pixvelflag = trackingPixVelFlag
            except:
                print('Setting velocity of targets equal to the non-sidereal tracking velocity')
                mtlist['x_or_RA_velocity'] = 0. - tracking_ra_vel #* (1./365.25/24.)
                mtlist['y_or_Dec_velocity'] = 0. - tracking_dec_vel #* (1./365.25/24.)
                pixvelflag = trackingPixVelFlag

        # Get necessary information for coordinate transformations
        coord_transform = None
        if self.runStep['astrometric']:
            # Read in the CRDS-format distortion reference file
            with AsdfFile.open(self.params['Reffiles']['astrometric']) as dist_file:
                coord_transform = dist_file.tree['model']
        else:
            coord_transform = self.simple_coord_transform()
                
        # Using the requested RA,Dec of the reference pixel, along with the
        # V2,V3 of the reference pixel, and the requested roll angle of the telescope,
        # create a matrix that can be used to translate between V2,V3 and RA,Dec
        # for any pixel.
        # v2,v3 need to be in arcsec, and RA, Dec, and roll all need to be in degrees
        attitude_matrix = self.getAttitudeMatrix()

        # Exposure times for all frames
        numints = self.params['Readout']['nint']
        numgroups = self.params['Readout']['ngroup']
        numframes = self.params['Readout']['nframe']
        numskips = self.params['Readout']['nskip']
        numresets = self.params['Readout']['resets_bet_ints']

        frames_per_group = numframes + numskips
        total_frames = numgroups * frames_per_group
        # If only one integration per exposure, then total_frames
        # above is correct. For >1 integration, we need to add the reset
        # frame to each integration (except the first), and sum the number of
        # frames for all integrations

        if numints > 1:
            # Frames for all integrations
            total_frames *= numints
            # Add the resets for all but the first and last integrations
            total_frames += (numresets * (numints - 1))

        frameexptimes = self.frametime * np.arange(-1,total_frames)

        #output image dimensions
        #dims = np.array(self.dark.data[0,0,:,:].shape)
        dims = self.nominal_dims
        newdimsx = np.int(dims[1] * self.coord_adjust['x'])
        newdimsy = np.int(dims[0] * self.coord_adjust['y'])

        # Set up seed integration
        #mt_integration = np.zeros((len(frameexptimes)-1, newdimsy, newdimsx))
        mt_integration = np.zeros((numints, numgroups * frames_per_group, newdimsy, newdimsx))

        # Corresponding (2D) segmentation map
        moving_segmap = segmap.SegMap()
        moving_segmap.xdim = newdimsx
        moving_segmap.ydim = newdimsy
        moving_segmap.initialize_map()

        for index, entry in zip(indexes, mtlist):
            # For each object, calculate x,y or RA,Dec of initial position
            pixelx, pixely, ra, dec, ra_str, dec_str = self.getPositions(
                entry['x_or_RA'], entry['y_or_Dec'], attitude_matrix,
                coord_transform, pixelFlag)

            # Now generate a list of x,y position in each frame
            if pixvelflag == False:
                # Calculate the RA,Dec in each frame
                # input velocities are arcsec/hour. ra/dec are in units of degrees,
                # so divide velocities by 3600^2.
                ra_frames = ra + (entry['x_or_RA_velocity'] / 3600. / 3600.) * frameexptimes
                dec_frames = dec + (entry['y_or_Dec_velocity'] / 3600. / 3600.) * frameexptimes

                x_frames = []
                y_frames = []
                for in_ra, in_dec in zip(ra_frames, dec_frames):
                    # Calculate the x,y position at each frame
                    px, py, pra, pdec, pra_str, pdec_str = self.getPositions(
                        in_ra, in_dec, attitude_matrix, coord_transform, False)
                    x_frames.append(px)
                    y_frames.append(py)
                x_frames = np.array(x_frames)
                y_frames = np.array(y_frames)

            else:
                # If input velocities are pixels/hour, then generate the list of
                # x,y in each frame directly
                x_frames = pixelx + (entry['x_or_RA_velocity'] / 3600.) * frameexptimes
                y_frames = pixely + (entry['y_or_Dec_velocity'] / 3600.) * frameexptimes


            # Load PSF from the appropriate file given the source location
            # Remember that this source location should be in the coordinate
            # system of the aperture being simulated, and not include any padding
            # for WFSS seed images.
            psf_obj = PSF(entry['pixelx'], entry['pixely'], self.psfname,
                            interval=self.params['simSignals']['psfpixfrac'], oversampling=1)

            # If we have a point source, we can easily determine whether
            # it completely misses the detector, since we know the size
            # of the stamp already. For galaxies and extended sources,
            # we have to get the stamp image first to see if any part of
            # the stamp lands on the detector.
            status = 'on'
            if input_type == 'pointSource':
                status = self.on_detector(x_frames, y_frames, psf_obj.model.shape,
                                          (newdimsx, newdimsy))
            if status == 'off':
                continue

            # So now we have xinit,yinit, a list of x,y positions for
            # each frame, and the frametime.
            # Subsample factor can be hardwired for now. outx and outy
            # are also known. So all we need is the stamp image, then
            # we can call moving_targets.py and feed it these things,
            # which contain all the info needed

            # Evaluate the PSF model. Place in stamp image just large enough to contain it
            eval_psf = psf_obj.minimal_psf_evaluation()
            
            if input_type == 'pointSource':
                stamp = eval_psf

            elif input_type == 'extended':
                stamp, header = self.basicGetImage(entry['filename'])
                if entry['pos_angle'] != 0.:
                    stamp = self.basicRotateImage(stamp, entry['pos_angle'])

                # Convolve with instrument PSF if requested
                if self.params['simSignals']['PSFConvolveExtended']:
                    stamp = s1.fftconvolve(stamp, eval_psf, mode='same')

            elif input_type == 'galaxies':
                stamp = self.create_galaxy(entry['radius'], entry['ellipticity'], entry['sersic_index'], entry['pos_angle'], 1.)
                # Convolve the galaxy with the instrument PSF
                stamp = s1.fftconvolve(stamp, eval_psf, mode='same')

            # Normalize the PSF to a total signal of 1.0
            totalsignal = np.sum(stamp)
            stamp /= totalsignal

            # Scale the stamp image to the requested magnitude
            rate = self.mag_to_countrate(magsys, entry['magnitude'],
                                         photfnu=self.photfnu,
                                         photflam=self.photflam)
            stamp *= rate

            # Now that we have stamp images for galaxies and extended
            # sources, check to see if they overlap the detector or not.
            # NOTE: this will only catch sources that never overlap the
            # detector for any of their positions.
            if input_type != 'pointSource':
                status = self.on_detector(x_frames, y_frames, stamp.shape,
                                          (newdimsx, newdimsy))
            if status == 'off':
                continue

            # Each entry will have stamp image as array, ra_init, dec_init,
            # ra_velocity, dec_velocity, frametime, numframes, subsample_factor,
            # outputarrayxsize, outputarrayysize
            # (maybe without the values that will be the same to each entry.

            #entryList = (stamp,ra,dec,entry[3]/3600.,entry[4]/3600.,self.frametime,numframes,subsample_factor,outx,outy)
            #entryList = (stamp,x_frames,y_frames,self.frametime,subsample_factor,outx,outy)

            # Need to feed info into moving_targets one integration at a time.
            # No need to feed in the reset frames, but they are necessary
            # before this point in order to get the timing and positions
            # correct.
            for integ in range(numints):
                framestart = integ * (frames_per_group * numgroups) + integ
                frameend = framestart + (frames_per_group * numgroups) + 1

                # Now check to see if the stamp image overlaps the output
                # aperture for this integration only. Above we removed sources
                # that never overlap the aperture. Here we want to get rid
                # of sources that overlap the detector in some integrations,
                # but not this particular integration
                status = 'on'
                status = self.on_detector(x_frames[framestart:frameend],
                                          y_frames[framestart:frameend],
                                          stamp.shape, (newdimsx, newdimsy))
                if status == 'off':
                    continue

                mt = moving_targets.MovingTarget()
                mt.subsampx = 3
                mt.subsampy = 3
                #mt_source = mt.create(stamp, x_frames, y_frames, self.frametime, newdimsx, newdimsy)
                #print("Prepping to input to moving_targets: integ: {}, framestart: {}, frameend: {}, xframelen: {}, yframelen: {}, {}, {}, {}".format(integ,framestart,frameend,len(x_frames),len(y_frames),newdimsx,newdimsy,stamp.shape))

                mt_source = mt.create(stamp, x_frames[framestart:frameend],
                                      y_frames[framestart:frameend],
                                      self.frametime, newdimsx, newdimsy)
                #mt_integration += mt_source
                mt_integration[integ, :, :, :] += mt_source

                noiseval = self.single_ron / 100. + self.params['simSignals']['bkgdrate']
                if self.params['Inst']['mode'].lower() in ['wfss','ts_wfss']:
                    noiseval += self.grism_background

                if input_type in ['pointSource', 'galaxies']:
                    moving_segmap.add_object_noise(mt_source[-1, :, :], 0, 0, index, noiseval)
                else:
                    indseg = self.seg_from_photutils(mt_source[-1, :, :], index, noiseval)
                    moving_segmap.segmap += indseg
        return mt_integration, moving_segmap.segmap

    def simple_coord_transform(self):
        """
        For the case where the distortion reference file is not provided, 
        generate a simple coordinate transform function that does not include
        distortion

        "Forward" transform = science -> ideal -> V2,V3

        Parameters:
        -----------
        None

        Returns:
        --------
        returns : astropy.modeling.Polynomial2D
            Compound model containing the necessary functions
            to transform coordinates assuming no distortion
        """
        xshift = Shift(0. - self.refpix_pos["x"])
        yshift = Shift(0. - self.refpix_pos["y"])
        pixelscalex = self.xsciscale
        pixelscaley = self.ysciscale

        # Science to Ideal transform
        cx = {}
        cx["c0_0"] = 0.
        cx["c1_0"] = pixelscalex
        cx["c1_1"] = 0.
        cx["c0_1"] = 0.
        x_sci_to_ideal_func = Polynomial2D(2, **cx)

        cy = {}
        cy["c0_0"] = 0.
        cy["c1_0"] = 0.
        cy["c1_1"] = 0.
        cy["c0_1"] = pixelscaley
        y_sci_to_ideal_func = Polynomial2D(2, **cy)

        # Ideal to V2, V3 transform
        parity = self.parity
        v3_ideal_y_angle = self.v3yang

        xc = {}
        yc = {}
        xc['c0_0'] = 0.
        yc['c0_0'] = 0.
        xc['c1_0'] = parity * np.cos(v3_ideal_y_angle)
        xc['c0_1'] = np.sin(v3_ideal_y_angle)
        yc['c1_0'] = (0.-parity) * np.sin(v3_ideal_y_angle)
        yc['c0_1'] = np.cos(v3_ideal_y_angle)
        x_ideal_to_v2v3_func = Polynomial2D(1, **xc)
        y_ideal_to_v2v3_func = Polynomial2D(1, **yc)

        # Ideal to Science transform
        cx = {}
        cx["c0_0"] = 0.
        cx["c1_0"] = 1./ pixelscalex
        cx["c1_1"] = 0.
        cx["c0_1"] = 0.
        x_ideal_to_sci_func = Polynomial2D(2, **cx)

        cy = {}
        cy["c0_0"] = 0.
        cy["c1_0"] = 0.
        cy["c1_1"] = 0.
        cy["c0_1"] = 1./ pixelscaley
        y_ideal_to_sci_func = Polynomial2D(2, **cy)

        # V2, V3 to Ideal transform
        xc = {}
        yc = {}
        xc['c0_0'] = 0.
        yc['c0_0'] = 0.
        xc['c1_0'] = parity * np.cos(v3_ideal_y_angle)
        xc['c0_1'] = parity * (0. - np.sin(v3_ideal_y_angle))
        yc['c1_0'] = np.sin(v3_ideal_y_angle)
        yc['c0_1'] = np.cos(v3_ideal_y_angle)
        x_v2v3_to_ideal_func = Polynomial2D(1, **xc)
        y_v2v3_to_ideal_func = Polynomial2D(1, **yc)

        # Shift by V2ref, V3ref
        v2shift = Shift(self.v2_ref)
        v3shift= Shift(self.v2_ref)
        
        #Now create a compound model for each with the appropriate
        #inverse
        sci2idl = Mapping([0,1,0,1]) | x_sci_to_ideal_func & y_sci_to_ideal_func
        sci2idl.inverse = Mapping([0,1,0,1]) | x_ideal_to_sci_func & y_ideal_to_sci_func

        idl2v2v3 = Mapping([0,1,0,1]) | x_ideal_to_v2v3_func & y_ideal_to_v2v3_func
        idl2v2v3.inverse = Mapping([0,1,0,1]) | x_v2v3_to_ideal_func & y_v2v3_to_ideal_func

        core_model = sci2idl | idl2v2v3
        model = xshift & yshift | core_model | v2shift & v3shift
        return model

    def on_detector(self,xloc, yloc, stampdim, finaldim):
        """Given a set of x, y locations, stamp image dimensions,
        and final image dimensions, determine whether the stamp
        image will overlap at all with the final image, or
        completely miss it.

        Parameters:
        ----------
        xloc : list
            x-coordinate locations of source
        yloc : list
            y-coordinate locations of source
        stampdim : tuple
            x,y dimension lengths of stamp image
        finaldim : tuple
            x,y dimension lengths of final image

        Returns:
        --------
        state of stamp image:
        "on" -- stamp image fully or partially overlaps the final
             image for at least one xloc, yloc pair
        "off" -- stamp image never overlaps with final image
        """
        status = 'on'
        stampx, stampy = stampdim
        stampminx = np.min(xloc - (stampx / 2))
        stampminy = np.min(yloc - (stampy / 2))
        stampmaxx = np.max(xloc + (stampx / 2))
        stampmaxy = np.max(yloc + (stampy / 2))
        finalx, finaly = finaldim
        if ((stampminx > finalx) or (stampmaxx < 0) or
            (stampminy > finaly) or (stampmaxy < 0)):
            status = 'off'
        return status

    def getPositions(self, inx, iny, matrix, transform, pixelflag):
        #input a row containing x,y or ra,dec values, and figure out
        #x,y, RA, Dec, and RA string and Dec string
        try:
            entry0 = float(inx)
            entry1 = float(iny)
            if not pixelflag:
                ra_str, dec_str = self.makePos(entry0, entry1)
                ra = entry0
                dec = entry1
        except:
            # if inputs can't be converted to floats, then
            # assume we have RA/Dec strings. Convert to floats.
            ra_str = inx
            dec_str = iny
            ra, dec = self.parseRADec(ra_str, dec_str)

        # Case where point source list entries are given with RA and Dec
        if not pixelflag:

            # If distortion is to be included - either with or without the full set of coordinate
            # translation coefficients
            if self.runStep['astrometric']:
                pixelx, pixely = self.RADecToXY_astrometric(ra, dec, matrix, transform)
            else:
                # No distortion at all - "manual mode"
                pixelx, pixely = self.RADecToXY_manual(ra, dec)

        else:
            # Case where the point source list entry locations are given in units of pixels
            # In this case we have the source position, and RA/Dec are calculated only so
            # they can be written out into the output source list file.

            pixelx = entry0
            pixely = entry1

            ra, dec, ra_str, dec_str = self.XYToRADec(pixelx, pixely, matrix, transform)
        return pixelx, pixely, ra, dec, ra_str, dec_str

    def nonsidereal_CRImage(self, file):
        """
        Create countrate image of non-sidereal sources
        that are being tracked.

        Arguments:
        ----------
        file : str
            name of the file containing the tracked moving targets.

        Returns:
        --------
        returns : obj
            countrate image (2D) containing the tracked non-sidereal targets.
        """
        totalCRList = []
        totalSegList = []

        # Read in file containing targets
        targs, pixFlag, velFlag, magsys = self.readMTFile(file)

        # We need to keep track of the proper motion of the
        # target being tracked, because all other objects in
        # the field of view will be moving at the same rate
        # in the opposite direction
        track_ra_vel = targs[0]['x_or_RA_velocity']
        track_dec_vel = targs[0]['y_or_Dec_velocity']

        # Sort the targets by whether they are point sources,
        # galaxies, extended
        ptsrc_rows = []
        galaxy_rows = []
        extended_rows = []
        for i, line in enumerate(targs):
            if 'point' in line['object'].lower():
                ptsrc_rows.append(i)
            elif 'sersic' in line['object'].lower():
                galaxy_rows.append(i)
            else:
                extended_rows.append(i)

        # Re-use functions for the sidereal tracking situation
        if len(ptsrc_rows) > 0:
            ptsrc = targs[ptsrc_rows]
            if pixFlag:
                meta0 = 'position_pixel'
            else:
                meta0 = ''
            if velFlag:
                meta1 = 'velocity_pixel'
            else:
                meta1 = ''
            meta2 = magsys

            meta3 = ('Point sources with non-sidereal tracking. '
                     'File produced by catalog_seed_image.py')
            meta4 = ('from run using non-sidereal moving target '
                     'list {}.'.format(self.params['simSignals']['movingTargetToTrack']))
            ptsrc.meta['comments'] = [meta0, meta1, meta2, meta3, meta4]
            ptsrc.write(os.path.join(self.params['Output']['directory'], 'temp_non_sidereal_point_sources.list'), format='ascii', overwrite=True)

            ptsrc = self.getPointSourceList('temp_non_sidereal_point_sources.list')
            ptsrcCRImage, ptsrcCRSegmap = self.makePointSourceImage(ptsrc)
            totalCRList.append(ptsrcCRImage)
            totalSegList.append(ptsrcCRSegmap)

        if len(galaxy_rows) > 0:
            galaxies = targs[galaxy_rows]
            if pixFlag:
                meta0 = 'position_pixel'
            else:
                meta0 = ''
            if velFlag:
                meta1 = 'velocity_pixel'
            else:
                meta1 = ''
            meta2 = magsys
            meta3 = ('Galaxies (2d sersic profiles) with non-sidereal '
                     'tracking. File produced by ramp_simulator.py')
            meta4 = ('from run using non-sidereal moving target '
                     'list {}.'.format(self.params['simSignals']['movingTargetToTrack']))
            galaxies.meta['comments'] = [meta0, meta1, meta2, meta3, meta4]
            galaxies.write(os.path.join(self.params['Output']['directory'], 'temp_non_sidereal_sersic_sources.list'), format='ascii', overwrite=True)

            galaxyCRImage, galaxySegmap = self.makeGalaxyImage('temp_non_sidereal_sersic_sources.list')
            galaxyCRImage *= self.pam
            totalCRList.append(galaxyCRImage)
            totalSegList.append(galaxySegmap)

        if len(extended_rows) > 0:
            extended = targs[extended_rows]

            if pixFlag:
                meta0 = 'position_pixel'
            else:
                meta0 = ''
            if velFlag:
                meta1 = 'velocity_pixel'
            else:
                meta1 = ''
            meta2 = magsys
            meta3 = 'Extended sources with non-sidereal tracking. File produced by ramp_simulator.py'
            meta4 = 'from run using non-sidereal moving target list {}.'.format(self.params['simSignals']['movingTargetToTrack'])
            extended.meta['comments'] = [meta0, meta1, meta2, meta3, meta4]
            extended.write(os.path.join(self.params['Output']['directory'], 'temp_non_sidereal_extended_sources.list'), format='ascii', overwrite=True)

            extlist, extstamps = self.getExtendedSourceList('temp_non_sidereal_extended_sources.list')

            # translate the extended source list into an image
            extCRImage, extSegmap = self.makeExtendedSourceImage(extlist, extstamps)

            # Multiply by the pixel area map
            extCRImage *= self.pam

            totalCRList.append(extCRImage)
            totalSegList.append(extSegmap)

        # Now combine into a final countrate image of non-sidereal sources (that are being tracked)
        if len(totalCRList) > 0:
            totalCRImage = totalCRList[0]
            totalSegmap = totalSegList[0]
            if len(totalCRList) > 1:
                for i in range(1, len(totalCRList)):
                    totalCRImage += totalCRList[i]
                    totalSegmap += totalSegList[i] + i
        else:
            raise ValueError(("No non-sidereal countrate targets produced."
                              "You shouldn't be here."))
        return totalCRImage, totalSegmap, track_ra_vel, track_dec_vel, velFlag

    def addedSignals(self):
        # Generate a signal rate image from input sources
        if self.params['Output']['grism_source_image'] == False:
            signalimage = np.zeros(self.nominal_dims)
            segmentation_map = np.zeros(self.nominal_dims)
        else:
            xd = np.int(self.nominal_dims[1] * self.coord_adjust['x'])
            yd = np.int(self.nominal_dims[0] * self.coord_adjust['y'])
            signalimage = np.zeros((yd, xd), dtype=np.float)
            segmentation_map = np.zeros((yd, xd))

        # yd, xd = signalimage.shape
        arrayshape = signalimage.shape

        #MASK IMAGE
        #Create a mask so that we don't add signal to masked pixels
        #Initially this includes only the reference pixels
        #Keep the mask image equal to the true subarray size, since this
        #won't be used to make a requested grism source image
        maskimage = np.zeros((self.ffsize, self.ffsize), dtype=np.int)
        maskimage[4:self.ffsize-4, 4:self.ffsize-4] = 1.

        #crop the mask to match the requested output array
        if "FULL" not in self.params['Readout']['array_name']:
            maskimage = maskimage[self.subarray_bounds[1]:self.subarray_bounds[3] + 1, self.subarray_bounds[0]:self.subarray_bounds[2] + 1]


        # POINT SOURCES
        # Read in the list of point sources to add
        # Adjust point source locations using astrometric distortion
        # Translate magnitudes to counts in a single frame
        if self.runStep['pointsource'] == True:
            pslist = self.getPointSourceList(self.params['simSignals']['pointsource'])

            # translate the point source list into an image
            psfimage, ptsrc_segmap = self.makePointSourceImage(pslist)

            # save the point source image for examination by user
            if self.params['Output']['save_intermediates'] == True:
                psfImageName = self.basename + '_pointSourceRateImage_adu_per_sec.fits'
                h0 = fits.PrimaryHDU(psfimage)
                h1 = fits.ImageHDU(ptsrc_segmap)
                hlist = fits.HDUList([h0, h1])
                hlist.writeto(psfImageName, overwrite=True)
                # ptsrc_seg_file = psfImageName[0:-5] + '_segmap.fits'
                # self.saveSingleFits(psfimage, psfImageName)
                print("Point source image and segmap saved as {}".format(psfImageName))
                # self.saveSingleFits(ptsrc_segmap, ptsrc_seg_file)
                # print("Point source segmentation map saved as {}".format(ptsrc_seg_file))

            # Add the point source image to the overall image
            signalimage = signalimage + psfimage
            segmentation_map += ptsrc_segmap

        # Simulated galaxies
        # Read in the list of galaxy positions/magnitudes to simulate
        # and create a countrate image of those galaxies.
        if self.runStep['galaxies'] == True:
            galaxyCRImage, galaxy_segmap = self.makeGalaxyImage(self.params['simSignals']['galaxyListFile'])

            # Multiply by the pixel area map
            galaxyCRImage *= self.pam
            
            # Add galaxy segmentation map to the master copy
            segmentation_map += galaxy_segmap

            # save the galaxy image for examination by the user
            if self.params['Output']['save_intermediates'] == True:
                galImageName = self.basename + '_galaxyRateImage_adu_per_sec.fits'
                h0 = fits.PrimaryHDU(galaxyCRImage)
                h1 = fits.ImageHDU(galaxy_segmap)
                hlist = fits.HDUList([h0, h1])
                hlist.writeto(galImageName, overwrite=True)
                # self.saveSingleFits(galaxyCRImage, galImageName)
                print("Simulated galaxy image and segmap saved as {}".format(galImageName))

            # add the galaxy image to the signalimage
            signalimage = signalimage + galaxyCRImage

        # read in extended signal image and add the image to the overall image
        if self.runStep['extendedsource'] == True:
            extlist, extstamps = self.getExtendedSourceList(self.params['simSignals']['extended'])

            # translate the extended source list into an image
            extimage, ext_segmap = self.makeExtendedSourceImage(extlist, extstamps)

            # Multiply by the pixel area map
            extimage *= self.pam
            
            # Add galaxy segmentation map to the master copy
            segmentation_map += ext_segmap

            # Save extended source image and segmap
            if self.params['Output']['save_intermediates'] == True:
                extImageName = self.basename + '_extendedObject_adu_per_sec.fits'
                h0 = fits.PrimaryHDU(extimage)
                h1 = fits.ImageHDU(ext_segmap)
                hlist = fits.HDUList([h0, h1])
                hlist.writeto(extImageName, overwrite=True)
                print("Extended object image and segmap saved as {}".format(extImageName))

            # add the extended image to the synthetic signal rate image
            signalimage = signalimage + extimage

        # ZODIACAL LIGHT
        if self.runStep['zodiacal'] == True:
            #zodiangle = self.eclipticangle() - self.params['Telescope']['rotation']
            zodiangle = self.params['Telescope']['rotation']
            zodiacalimage, zodiacalheader = self.getImage(self.params['simSignals']['zodiacal'], arrayshape, True, zodiangle, arrayshape/2)

            # Multiply by the pixel area map
            zodiacalimage *= self.pam
            
            signalimage = signalimage + zodiacalimage*self.params['simSignals']['zodiscale']

        # SCATTERED LIGHT - no rotation here.
        if self.runStep['scattered']:
            scatteredimage, scatteredheader = self.getImage(self.params['simSignals']['scattered'],
                                                            arrayshape, False, 0.0, arrayshape/2)

            # Multiply by the pixel area map
            scatteredimage *= self.pam
            
            signalimage = signalimage + scatteredimage*self.params['simSignals']['scatteredscale']

        # CONSTANT BACKGROUND - multiply by pixel area map first
        signalimage = signalimage + self.params['simSignals']['bkgdrate'] * self.pam

        # Save the final rate image of added signals
        if self.params['Output']['save_intermediates'] == True:
            rateImageName = self.basename + '_AddedSources_adu_per_sec.fits'
            self.saveSingleFits(signalimage, rateImageName)
            print("Signal rate image of all added sources saved as {}".format(rateImageName))

        return signalimage, segmentation_map

    def getDistortionCoefficients(self, table, from_sys, to_sys, aperture):
        '''from the table of distortion coefficients, get the coeffs that correspond
        to the requested transformation and return as a list for x and another for y
        '''
        match = table['AperName'] == aperture
        if np.any(match) == False:
            raise ValueError("Aperture name {} not found in input CSV file.".format(aperture))

        row = table[match]

        if ((from_sys == 'science') & (to_sys == 'ideal')):
            label = 'Sci2Idl'
        elif ((from_sys == 'ideal') & (to_sys == 'science')):
            label = 'Idl2Sci'
        else:
            raise ValueError(("WARNING: from_sys of {} and to_sys of {} not "
                              "a valid transformation.".format(from_sys, to_sys)))
            
        # get the coefficients, return as list
        X_cols = [c for c in row.colnames if label + 'X' in c]
        Y_cols = [c for c in row.colnames if label + 'Y' in c]
        x_coeffs = [row[c].data[0] for c in X_cols]
        y_coeffs = [row[c].data[0] for c in Y_cols]

        # Strip off masked coefficients, where the column exists but there
        # is no corresponding coefficient in the SIAF
        x_coeffs = [c for c in x_coeffs if np.isfinite(c)]
        y_coeffs = [c for c in y_coeffs if np.isfinite(c)]
        #x_coeffs = [c if np.isfinite(c) else 0. for c in x_coeffs]
        #y_coeffs = [c if np.isfinite(c) else 0. for c in y_coeffs]

        # Also get the V2, V3 values of the reference pixel
        v2ref = row['V2Ref'].data[0]
        v3ref = row['V3Ref'].data[0]

        # Get parity and V3 Y angle info
        parity = row['VIdlParity'].data[0]
        yang = row['V3IdlYAngle'].data[0]
        v3scixang = row['V3SciXAngle'].data[0]

        # Get pixel scale info - not used but needs to be in output
        xsciscale = row['XSciScale'].data[0]
        ysciscale = row['YSciScale'].data[0]

        return x_coeffs, y_coeffs, v2ref, v3ref, parity, yang, xsciscale, ysciscale, v3scixang

    def getPointSourceList(self, filename):
        # read in the list of point sources to add, and adjust the
        # provided positions for astrometric distortion

        # find the array sizes of the PSF files in the library. Assume they are all the same.
        # We want the distance from the PSF peak to the edge, assuming the peak is centered
        if self.params['simSignals']['psfwfe'] != 0:
            numstr = str(self.params['simSignals']['psfwfe'])
        else:
            numstr = 'zero'
        #psflibfiles = glob.glob(self.params['simSignals']['psfpath'] + '*')
        psflibfiles = glob.glob(os.path.join(self.params['simSignals']['psfpath'], '*.fits'))

        # If a PSF library is specified, then just get the dimensions from one of the files
        if self.params['simSignals']['psfpath'] != None:
            h = fits.open(psflibfiles[0])
            edgex = h[0].header['NAXIS1'] / 2 - 1
            edgey = h[0].header['NAXIS2'] / 2 - 1
            self.psfhalfwidth = np.array([edgex, edgey])
            h.close()
        else:
            # if no PSF library is specified, then webbpsf will be creating the PSF on the
            # fly. In this case, we assume webbpsf's default output size of 301x301 pixels?
            edgex = int(301 / 2)
            edgey = int(301 / 2)
            print("INFO: no PSF library specified, but point sources are to be added to")
            print("the output. PSFs will be generated by WebbPSF on the fly")
            raise NotImplementedError("Not yet implemented.")

        pointSourceList = Table(names=('index', 'pixelx', 'pixely', 'RA', 'Dec', 'RA_degrees',
                                       'Dec_degrees', 'magnitude', 'countrate_e/s',
                                       'counts_per_frame_e'),
                                dtype=('i', 'f', 'f', 'S14', 'S14', 'f', 'f', 'f', 'f', 'f'))

        try:
            lines, pixelflag, magsys = self.readPointSourceFile(filename)
            if pixelflag:
                print("Point source list input positions assumed to be in units of pixels.")
            else:
                print("Point list input positions assumed to be in units of RA and Dec.")
        except:
            raise NameError("WARNING: Unable to open the point source list file {}".format(filename))

        # File to save adjusted point source locations
        psfile = self.params['Output']['file'][0:-5] + '_pointsources.list'
        pslist = open(psfile, 'w')

        # If the input catalog has an index column
        # use that, otherwise add one
        if 'index' in lines.colnames:
            print('Using point source catalog index numbers')
            indexes = lines['index']
        else:
            print('No point source catalog index numbers. Adding to output: {}.'.format(psfile))
            indexes = np.arange(1, len(lines['x_or_RA']) + 1)
        # Make sure there is no 0th object
        if np.min(indexes) == 0:
            indexes += 1
        # Make sure the index numbers don't overlap with any
        # sources already present. Increment the maxindex
        # value.
        if np.min(indexes) <= self.maxindex:
            indexes += self.maxindex
        self.maxindex = np.max(indexes)
        print("After point sources, max index is {}".format(self.maxindex))

        dtor = math.radians(1.)
        nx = (self.subarray_bounds[2] - self.subarray_bounds[0]) + 1
        ny = (self.subarray_bounds[3] - self.subarray_bounds[1]) + 1
        xc = (self.subarray_bounds[2] + self.subarray_bounds[0]) / 2.
        yc = (self.subarray_bounds[3] + self.subarray_bounds[1]) / 2.

        # Location of the subarray's reference pixel.
        xrefpix = self.refpix_pos['x']
        yrefpix = self.refpix_pos['y']

        # center positions, sub-array sizes in pixels
        # now offset the field center to array center for astrometric distortion corrections
        coord_transform = None
        if self.runStep['astrometric']:
            # Read in the CRDS-format distortion reference file
            with AsdfFile.open(self.params['Reffiles']['astrometric']) as dist_file:
                coord_transform = dist_file.tree['model']
        else:
            coord_transform = self.simple_coord_transform()
            
        # Using the requested RA, Dec of the reference pixel, along with the
        # V2, V3 of the reference pixel, and the requested roll angle of the telescope
        # create a matrix that can be used to translate between V2, V3 and RA, Dec
        # for any pixel.
        # v2, v3 need to be in arcsec, and RA, Dec, and roll all need to be in degrees
        attitude_matrix = self.getAttitudeMatrix()

        #Define the min and max source locations (in pixels) that fall onto the subarray
        #Include the effects of a requested grism_direct image, and also keep sources that
        #will only partially fall on the subarray
        #pixel coords here can still be negative and kept if the grism image is being made
        
        # First, coord limits for just the subarray
        miny = 0
        maxy = self.subarray_bounds[3] - self.subarray_bounds[1]
        minx = 0
        maxx = self.subarray_bounds[2] - self.subarray_bounds[0]

        # Expand the limits if a grism direct image is being made
        if self.params['Output']['grism_source_image'] == True:
            extrapixy = np.int((maxy + 1)/2 * (self.coord_adjust['y'] - 1.))
            miny -= extrapixy
            maxy += extrapixy
            extrapixx = np.int((maxx + 1)/2 * (self.coord_adjust['x'] - 1.))
            minx -= extrapixx
            maxx += extrapixx

        # Now, expand the dimensions again to include point
        # sources that fall only partially on the subarray
        miny -= edgey
        maxy += edgey
        minx -= edgex
        maxx += edgex

        # Write out the RA and Dec of the field center to the output file
        # Also write out column headers to prepare for source list
        pslist.write(("# Field center (degrees): %13.8f %14.8f y axis rotation angle "
                      "(degrees): %f  image size: %4.4d %4.4d\n" %
                      (self.ra, self.dec, self.params['Telescope']['rotation'], nx, ny)))
        pslist.write('#\n')
        pslist.write(("#    Index   RA_(hh:mm:ss)   DEC_(dd:mm:ss)   RA_degrees      "
                      "DEC_degrees     pixel_x   pixel_y    magnitude   counts/sec    counts/frame\n"))

        start_time = time.time()
        times = []
        # Loop over input lines in the source list
        for index, values in zip(indexes, lines):
            try:
            #line below (if 1>0) used to keep the block
            # of code below at correct indent for the try: above
            # the try: is commented out for code testing.
            #if 1 > 0:
                # Warn user of how long this calcuation might take...
                if index < 100:
                    elapsed_time = time.time() - start_time
                    times.append(elapsed_time)
                    start_time = time.time()
                elif index == 100:
                    avg_time = np.mean(times)
                    total_time = len(indexes) * avg_time
                    print(("Expected time to process {} sources: {:.2f} seconds "
                           "({:.2f} minutes)".format(len(indexes), total_time, total_time/60)))

                try:
                    entry0 = float(values['x_or_RA'])
                    entry1 = float(values['y_or_Dec'])
                    if not pixelflag:
                        ra_str, dec_str = self.makePos(entry0, entry1)
                        ra = entry0
                        dec = entry1
                except:
                    # if inputs can't be converted to floats, then
                    # assume we have RA/Dec strings. Convert to floats.
                    ra_str = values['x_or_RA']
                    dec_str = values['y_or_Dec']
                    ra, dec = self.parseRADec(ra_str, dec_str)

                # Case where point source list entries are given with RA and Dec
                if not pixelflag:
                     
                    # If distortion is to be included - either with or without the full set of coordinate
                    # translation coefficients
                    #if self.runStep['astrometric']:

                    # Same function call regardless of whether distortion file is provided or not
                    pixelx, pixely = self.RADecToXY_astrometric(ra, dec, attitude_matrix, coord_transform)
                    #else:
                    #    # No distortion at all - "manual mode"
                    #    pixelx, pixely = self.RADecToXY_manual(ra, dec)

                else:
                    # Case where the point source list entry locations are given in units of pixels
                    # In this case we have the source position, and RA/Dec are calculated only so
                    # they can be written out into the output source list file.

                    # Assume that the input x and y values are coordinate values
                    # WITHIN THE SPECIFIED SUBARRAY. So for example, a source in the file
                    # at 0, 0 when you are making a SUB160 ramp will fall on the lower left
                    # corner of the SUB160 subarray, NOT the lower left corner of the full
                    # frame.

                    pixelx = entry0
                    pixely = entry1

                    ra, dec, ra_str, dec_str = self.XYToRADec(pixelx, pixely, attitude_matrix,
                                                              coord_transform)

                # Get the input magnitude of the point source
                mag = float(values['magnitude'])

                if pixely > miny and pixely < maxy and pixelx > minx and pixelx < maxx:
                    # set up an entry for the output table
                    entry = [index, pixelx, pixely, ra_str, dec_str, ra, dec, mag]

                    # translate magnitudes to countrate
                    # scale = 10.**(0.4*(15.0-mag))

                    # get the countrate that corresponds to a 15th magnitude star for this filter
                    # if self.params['Readout']['pupil'][0].upper() == 'F':
                    #    usefilt = 'pupil'
                    # else:
                    #    usefilt = 'filter'
                    # cval = self.countvalues[self.params['Readout'][usefilt]]
                    #
                    # DEAL WITH THIS LATER, ONCE PYSYNPHOT IS INCLUDED WITH PIPELINE DIST?
                    # if cval == 0:
                    #    print("Countrate value for {} is zero in {}.".format(self.params['Readout'][usefilt], self.parameters['phot_file']))
                    #    print("Eventually attempting to calculate value using pysynphot.")
                    #    print("but pysynphot is not present in jwst build 6, so pushing off to later...")
                    #    sys.exit()
                    #    cval = self.findCountrate(self.params['Readout'][usefilt])
                    #
                    # translate to counts in single frame at requested array size
                    # framecounts = scale*cval*self.frametime
                    # countrate = scale*cval

                    # Calculate the countrate for the source
                    countrate = self.mag_to_countrate(magsys, mag, photfnu=self.photfnu, photflam=self.photflam)
                    framecounts = countrate * self.frametime

                    # add the countrate and the counts per frame to pointSourceList
                    # since they will be used in future calculations
                    entry.append(countrate)
                    entry.append(framecounts)

                    # add the good point source, including location and counts, to the pointSourceList
                    pointSourceList.add_row(entry)

                    # write out positions, distances, and counts to the output file
                    pslist.write("%i %s %s %14.8f %14.8f %9.3f %9.3f  %9.3f  %13.6e   %13.6e\n" % (index, ra_str, dec_str, ra, dec, pixelx, pixely, mag, countrate, framecounts))
            except:
                pass
        self.n_pointsources = len(pointSourceList)
        print("Number of point sources found within the requested aperture: {}".format(self.n_pointsources))
        # close the output file
        pslist.close()

        # If no good point sources were found in the requested array, alert the user
        if len(pointSourceList) < 1:
            print("INFO: no point sources within the requested array.")
            # print("The point source image option is being turned off")
            # self.runStep['pointsource']=False
            # if self.runStep['extendedsource'] == False and self.runStep['cosmicray'] == False:
            #    print("Error: no input point sources, extended image, nor cosmic rays specified")
            #    print("Exiting...")
            #    sys.exit()

        return pointSourceList

    def makePointSourceImage(self, pointSources):
        dims = np.array(self.nominal_dims)

        # Offset that needs to be applied to the x, y positions of the
        # source list to account for case where we make a
        # seed image that is extra-large, to be used as a grism
        # direct image
        deltax = 0
        deltay = 0

        newdimsx = np.int(dims[1] * self.coord_adjust['x'])
        newdimsy = np.int(dims[0] * self.coord_adjust['y'])
        deltax = self.coord_adjust['xoffset']
        deltay = self.coord_adjust['yoffset']
        dims = np.array([newdimsy, newdimsx])

        # Create the empty image
        psfimage = np.zeros((dims[0], dims[1]))

        # Create empty segmentation map
        seg = segmap.SegMap()
        seg.xdim = newdimsx
        seg.ydim = newdimsy
        seg.initialize_map()

        #Loop over the entries in the point source list
        for entry in pointSources:
            # Adjust x, y position if the grism output image is requested
            xpos = entry['pixelx'] + deltax
            ypos = entry['pixely'] + deltay
            nx = math.floor(xpos)
            ny = math.floor(ypos)

            # Desired counts per second in the point source
            counts = entry['countrate_e/s'] # / self.frametime

            # Load PSF from the appropriate file given the source location
            # Remember that this source location should be in the coordinate
            # system of the aperture being simulated, and not include any padding
            # for WFSS seed images.
            psf_obj = PSF(entry['pixelx'], entry['pixely'], self.psfname,
                            interval=self.params['simSignals']['psfpixfrac'], oversampling=1)

            # Calculate the coordinate limits of the aperture/PSF stamp overlap
            psf_ydim, psf_xdim = psf_obj.model.shape
            (i1, i2, j1, j2, k1, k2, l1, l2) = self.cropped_coords(xpos, ypos, psf_xdim, psf_ydim,
                                                                   newdimsx, newdimsy)

            try:
                ypts, xpts = np.mgrid[j1:j2, i1:i2]
                scaled_psf = psf_obj.model.evaluate(x=xpts, y=ypts, flux=counts, x_0=nx, y_0=ny)
                psfimage[ypts, xpts] += scaled_psf
                # Divide readnoise by 100 sec, which is a 10 group RAPID ramp?
                noiseval = self.single_ron / 100. + self.params['simSignals']['bkgdrate']
                if self.params['Inst']['mode'].lower() in ['wfss','ts_wfss']:
                    noiseval += self.grism_background
                seg.add_object_noise(scaled_psf[l1:l2, k1:k2], j1, i1, entry['index'], noiseval)
            except:
                # In here we catch sources that are off the edge
                # of the detector. These may not necessarily be caught in
                # getpointsourcelist because if the PSF is not centered
                # in the webbpsf stamp, then the area to be pulled from
                # the stamp may shift off of the detector.
                pass

        return psfimage, seg.segmap

    def cropped_coords(self, x_det, y_det, psfxdim, psfydim, aperturexdim, apertureydim):
        """
        Given the location of a source on the detector along with the size of
        the PSF/stamp image for that object, calcuate the limits of the detector
        coordinates onto which the object will fall.

        Parameters:
        -----------
        x_det : flt
            Column location of source on detector (aperture coordinate system
            including any padding for WFSS seed image)
        
        y_det : flt
            Row location of source on detector (aperture coordinate system
            including any padding for WFSS seed image)

        psfxdim : int
            Number of columns in the PSF/stamp image that are being placed
            onto the detector/aperture

        psfydim : int
            Number of rows in the PSF/stamp image that are being placed
            onto the detector/aperture

        aperturexdim : int
            Number of columns in the final seed image

        apertureydim : int
            Number of rows in the final seed image

        Returns:
        --------
        i1 : int
            Column number on the detector/aperture corresponding to the left
            edge of the PSF/stamp image.

        i2 : int
            Column number on the detector/aperture corresponding to the right
            edge of the PSF/stamp image.

        j1 : int
            Row number on the detector/aperture corresponding to the bottom
            edge of the PSF/stamp image.

        j2 : int
            Row number on the detector/aperture corresponding to the top
            edge of the PSF/stamp image.

        l1 : int
            Column number on the PSF/stamp image corresponding to the left-most
            column that overlaps the detector/aperture

        l2 : int
            Column number on the PSF/stamp image corresponding to the right-most
            column that overlaps the detector/aperture

        k1 : int
            Row number on the PSF/stamp image corresponding to the bottom-most
            row that overlaps the detector/aperture

        k2 : int
            Row number on the PSF/stamp image corresponding to the top-most
            row that overlaps the detector/aperture
        """
        # Find sub-pixel offsets in position from the center of the pixel
        nx = math.floor(x_det)
        ny = math.floor(y_det)
    
        # Assume the center of the PSF is the center of the array
        nxshift = psfxdim // 2
        nyshift = psfydim // 2

        #nx = int(xoff)
        #ny = int(yoff)
        i1 = max(nx - nxshift, 0)
        i2 = min(nx + 1 + nxshift, aperturexdim)
        j1 = max(ny - nyshift, 0)
        j2 = min(ny + 1 + nyshift, apertureydim)
        k1 = nxshift - (nx - i1)
        k2 = nxshift + (i2 - nx)
        l1 = nyshift - (ny - j1)
        l2 = nyshift + (j2 - ny)
        
        # If the cutout for the psf is larger than
        # the psf array, truncate it, along with the array
        # in the source image where it will be placed
        # This is used only for segmap now that epsf is
        # implemented.
        if l2 > psfydim:
            l2 = psfydim
            j2 = j1 + (l2 - l1)
        
        if k2 > psfxdim:
            k2 = psfxdim
            i2 = i1 + (k2 - k1)

        # At this point coordinates are in the final output array coordinate system, so there
        # should be no negative values, nor values larger than the output array size
        if j1 < 0 or i1 < 0 or l1 < 0 or k1 < 0:
            raise ValueError(("WARNING: Coordinates of overlap between detector aperture"
                              "and PSF/stamp image are < 0."))
        if j2 > (apertureydim + 1) or i2 > (aperturexdim + 1) or l2 > (psfxdim + 1) or k2 > (psfydim + 1):
            raise ValueError(("WARNING: Coordinates of overlap between detector aperture:"
                              "and PSF/stamp are larger than the dimensions of the image."))
        return (i1, i2, j1, j2, k1, k2, l1, l2)

    def cropPSF(self, psf):
        '''take an array containing a psf and crop it such that the brightest
        pixel is in the center of the array'''
        nyshift, nxshift = np.where(psf == np.max(psf))
        nyshift = nyshift[0]
        nxshift = nxshift[0]
        py, px = psf.shape

        xl = nxshift - 0
        xr = px - nxshift - 1
        if xl <= xr:
            xdist = xl
        if xr < xl:
            xdist = xr

        yl = nyshift - 0
        yr = py - nyshift - 1
        if yl <= yr:
            ydist = yl
        if yr < yl:
            ydist = yr

        return psf[nyshift - ydist:nyshift + ydist + 1, nxshift - xdist:nxshift + xdist + 1]

    def readPointSourceFile(self, filename):
        # Read in the point source list
        try:
            gtab = ascii.read(filename)
            # Look at the header lines to see if inputs
            # are in units of pixels or RA, Dec
            pflag = False
            try:
                if 'position_pixels' in gtab.meta['comments'][0:4]:
                    pflag = True
            except:
                pass
            # Check to see if magnitude system is specified
            # in the comments. If not default to AB mag
            msys = 'abmag'
            condition=('stmag' in gtab.meta['comments'][0:4]) | ('vegamag' in gtab.meta['comments'][0:4])
            if condition:
                msys = [l for l in gtab.meta['comments'][0:4] if 'mag' in l][0]
                msys = msys.lower()

        except:
            raise IOError("WARNING: Unable to open the source list file {}".format(filename))

        return gtab, pflag, msys

    def getAttitudeMatrix(self):
        # create an attitude_matrix from the distortion reference file model and other info
        # calculate a local roll angle for the aperture
        self.local_roll = set_telescope_pointing.compute_local_roll(self.params['Telescope']['rotation'], self.ra, self.dec, self.v2_ref, self.v3_ref)
        # create attitude_matrix
        attitude_matrix = rotations.attitude(self.refpix_pos['v2'], self.refpix_pos['v3'], self.ra, self.dec, self.local_roll)
        return attitude_matrix

    def makePos(self, alpha1, delta1):
        # given a numerical RA/Dec pair, convert to string
        # values hh:mm:ss
        if alpha1 < 0.:
            alpha1 = alpha1 + 360.
        if delta1 < 0.:
            sign = "-"
            d1 = abs(delta1)
        else:
            sign = " + "
            d1 = delta1
        decd = int(d1)
        value = 60. * (d1 - float(decd))
        decm = int(value)
        decs = 60. * (value - decm)
        a1 = alpha1/15.0
        radeg = int(a1)
        value = 60. * (a1 - radeg)
        ramin = int(value)
        rasec = 60. * (value - ramin)
        alpha2 = "%2.2d:%2.2d:%7.4f" % (radeg, ramin, rasec)
        delta2 = "%1s%2.2d:%2.2d:%7.4f" % (sign, decd, decm, decs)
        alpha2 = alpha2.replace(" ", "0")
        delta2 = delta2.replace(" ", "0")
        return alpha2, delta2

    def parseRADec(self, rastr, decstr):
        # convert the input RA and Dec strings to floats
        try:
            rastr = rastr.lower()
            rastr = rastr.replace("h", ":")
            rastr = rastr.replace("m", ":")
            rastr = rastr.replace("s", "")
            decstr = decstr.lower()
            decstr = decstr.replace("d", ":")
            decstr = decstr.replace("m", ":")
            decstr = decstr.replace("s", "")

            values = rastr.split(":")
            ra0 = 15.*(int(values[0]) + int(values[1])/60. + float(values[2])/3600.)

            values = decstr.split(":")
            if "-" in values[0]:
                sign = -1
                values[0] = values[0].replace("-", " ")
            else:
                sign =  +1
            dec0 = sign*(int(values[0]) + int(values[1])/60. + float(values[2])/3600.)
            return ra0, dec0
        except:
            raise ValueError("Error parsing RA, Dec strings: {} {}".format(rastr, decstr))

    def RADecToXY_astrometric(self, ra, dec, attitude_matrix, coord_transform):
        # Translate backwards, RA, Dec to V2, V3
        pixelv2, pixelv3 = rotations.getv2v3(attitude_matrix, ra, dec)

        if self.runStep['distortion_coeffs']:
            # If the full set of distortion coefficients are provided, then
            # use those to make the exact transformation from the 'ideal'
            # to 'science' coordinate systems

            # Now V2, V3 to undistorted angular distance from the reference pixel
            xidl = self.v2v32idlx(pixelv2 - self.v2_ref, pixelv3 - self.v3_ref)
            yidl = self.v2v32idly(pixelv2 - self.v2_ref, pixelv3 - self.v3_ref)

            # Finally, undistorted distances to distorted pixel values
            ncoeff = len(self.x_sci2idl)
            if ncoeff == 21:
                polynomial_order = 5
            elif ncoeff == 15:
                polynomial_order = 4
            elif ncoeff == 10:
                polynomial_order = 3
            else:
                raise ValueError(("WARNING: {} science to ideal coefficients read in from "
                                  "SIAF. Not sure what polynomial order this corresponds to."
                                  .format(ncoeff)))

            deltapixelx, deltapixely, err, iter = polynomial.invert(self.x_sci2idl, self.y_sci2idl,
                                                                    xidl, yidl, polynomial_order)
            pixelx = deltapixelx + self.refpix_pos['x']
            pixely = deltapixely + self.refpix_pos['y']

        else:
            # If the full set of distortion coefficients are not provided,
            # then we fall back to the coordinate transform provided by the
            # distortion reference file. These results are not exact, and
            # become less accurate the farther the source is from the center
            # of the detector. Results can be incorrect by ~20 pixels in the
            # corners of the detector.

            # Now go backwards from V2, V3 to distorted pixels
            # deltapixelx, deltapixely = coord_transform.inverse(pixelv2-self.refpix_pos['v2'], pixelv3-self.refpix_pos['v3'])
            pixelx, pixely = coord_transform.inverse(pixelv2, pixelv3)

        return pixelx, pixely

    def RADecToXY_manual(self, ra, dec):
        # In this case, the sources are provided as an RA, Dec list,
        # but no astrometry information is provided. So assume an average
        # pixel scale and calculate the pixel position of the source from that.
        # This obviously does not include distortion, and is kind of a last
        # resort.
        ra_source = ra * 3600.
        dec_source = dec * 3600.

        #dist_between, deltaang = self.object_separation([self.ra, self.dec],
        #                                               [ra_source, dec_source])

        #simple_wcs_obj = {}
        #simple_wcs_obj["CRPIX1"] = self.refpix_pos["x"]
        #simple_wcs_obj["CRPIX2"] = self.refpix_pos["y"]
        #simple_wcs_obj["CRVAL1"] =
        #simple_wcs_obj["CRVAL2"] =
        #simple_wcs_obj["CTYPE1"] =
        #simple_wcs_obj["CTYPE2"] =
        #simple_wcs_obj["CD1_1"] =
        #simple_wcs_obj["CD1_2"] =
        #simple_wcs_obj["CD2_1"] =
        #simple_wcs_obj["CD2_2"] = 

        #wcs_obj = astropy.wcs.WCS(header = simple_wcs_obj)
        deltara, deltadec = self.object_separation([self.ra, self.dec],
                                                   [ra_source, dec_source], wcs_obj)
        # Now translate to deltax and deltay if the
        # position angle is non-zero
        #tot_ang = deltaang + (0. - self.params['Telescope']['rotation'] * np.pi / 180.)

        #deltax = dist_between * np.sin(tot_ang) / self.pixscale[0]
        #deltay = dist_between * np.cos(tot_ang) / self.pixscale[0]

        #pixelx = self.refpix_pos['x'] + deltax
        #pixely = self.refpix_pos['y'] + deltay

        
        
        
        return pixelx, pixely

    def object_separation(self, radec1, radec2, wcs):
        """
        Calculate the distance between two points on the sky given their
        RA, Dec values. Also calculate the angle (east of north?) between
        the two points.

        Parameters:
        -----------
        radec1 : list
            2-element list giving the RA, Dec (in decimal degrees) for 
            the first object

        radec2 : list
            2-element list giving the RA, Dec (in decimal degrees) for
            the second object

        Returns:
        --------
        distance : float
            Angular separation (in degrees) between the two objects

        angle : float
            Angle (east of north?) separating the two sources
        """
        c1 = SkyCoord(radec1[0]*u.degree, radec1[1]*u.degree, frame='icrs')
        c2 = SkyCoord(radec2[0]*u.degree, radec2[1]*u.degree, frame='icrs')
        sepra, sepdec = c1.spherical_offsets_to(c2).to_pixel(wcs)
        return sepra, sepdec
        
    def XYToRADec(self, pixelx, pixely, attitude_matrix, coord_transform):
        # Translate a given x, y location on the detector
        # to RA, Dec

        # If distortion is to be included
        # if self.runStep['astrometric']:
        if coord_transform is not None:
            # Transform distorted pixels to V2, V3
            # deltav2, deltav3 = coord_transform(pixelx-self.refpix_pos['x'], pixely-self.refpix_pos['y'])
            # pixelv2 = deltav2 + self.refpix_pos['v2']
            # pixelv3 = deltav3 + self.refpix_pos['v3']
            pixelv2, pixelv3 = coord_transform(pixelx, pixely)

            # Now translate V2, V3 to RA, Dec
            ra, dec = rotations.pointing(attitude_matrix, pixelv2, pixelv3)

        else:
            # Without including distortion.
            # Fall back to "manual" calculations
            dist_between = np.sqrt((pixelx - self.refpix_pos['x'])**2 + (pixely - self.refpix_pos['y'])**2)
            deltaang = np.arctan2(pixely, pixelx)

            tot_ang = deltaang + (self.parms['Telescope']['rotation'] * np.pi / 180.)

            deltara = dist_between * np.sin(tot_ang) / self.pixoscale[0]
            deltadec = dist_between * np.cos(tot_ang) / self.pixscale[0]

            ra = self.ra + deltara
            dec = self.dec + deltadec

        # Translate the RA/Dec floats to strings
        ra_str, dec_str = self.makePos(ra, dec)

        return ra, dec, ra_str, dec_str

    def readGalaxyFile(self, filename):
        # Read in the galaxy source list
        try:
            # read table
            gtab = ascii.read(filename)

            # Look at the header lines to see if inputs
            # are in units of pixels or RA, Dec
            pflag = False
            rpflag = False
            try:
                if 'position_pixel' in gtab.meta['comments'][0:4]:
                    pflag = True
            except:
                pass
            try:
                if 'radius_pixel' in gtab.meta['comments'][0:4]:
                    rpflag = True
            except:
                pass
            # Check to see if magnitude system is specified in the comments
            # If not assume AB mags
            msys = 'abmag'
            condition=('stmag' in gtab.meta['comments'][0:4]) | ('vegamag' in gtab.meta['comments'][0:4])
            if condition:
                msys = [l for l in gtab.meta['comments'][0:4] if 'mag' in l][0]
                msys = msys.lower()

        except:
            raise FileNotFoundError("WARNING: Unable to open the galaxy source list file {}".format(filename))

        return gtab, pflag, rpflag, msys


    def filterGalaxyList(self, galaxylist, pixelflag, radiusflag, magsystem):
        # given a list of galaxies (location, size, orientation, magnitude)
        # keep only those which will fall fully or partially on the output array

        filteredList = Table(names=('index', 'pixelx', 'pixely', 'RA', 'Dec',
                                    'RA_degrees', 'Dec_degrees', 'V2', 'V3',
                                    'radius', 'ellipticity', 'pos_angle',
                                    'sersic_index', 'magnitude', 'countrate_e/s',
                                    'counts_per_frame_e'),
                             dtype=('i', 'f', 'f', 'S14', 'S14', 'f', 'f', 'f',
                                    'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f'))

        #each entry in galaxylist is:
        #index x_or_RA  y_or_Dec  radius  ellipticity  pos_angle  sersic_index  magnitude
        #remember that x/y are interpreted as coordinates in the output subarray
        #NOT full frame coordinates. This is the same as the point source list coords

        # each entry in galaxylist is:
        # x_or_RA  y_or_Dec  radius  ellipticity  pos_angle  sersic_index  magnitude
        # remember that x/y are interpreted as coordinates in the output subarray
        # NOT full frame coordinates. This is the same as the point source list coords

        # First, begin to define the pixel limits beyond which a galaxy will be completely
        # outside of the field of view
        # First, coord limits for just the subarray
        miny = 0
        maxy = self.subarray_bounds[3] - self.subarray_bounds[1]
        minx = 0
        maxx = self.subarray_bounds[2] - self.subarray_bounds[0]
        ny = self.subarray_bounds[3] - self.subarray_bounds[1]
        nx = self.subarray_bounds[2] - self.subarray_bounds[0]

        #Expand the limits if a grism direct image is being made
        if self.params['Output']['grism_source_image'] == True:
            extrapixy = np.int((maxy + 1)/2 * (self.grism_direct_factor - 1.))
            miny -= extrapixy
            maxy += extrapixy
            extrapixx = np.int((maxx + 1)/2 * (self.grism_direct_factor - 1.))
            minx -= extrapixx
            maxx += extrapixx

            nx = np.int(nx * self.grism_direct_factor)
            ny = np.int(ny * self.grism_direct_factor)

        # Create transform matrix for galaxy sources
        # Read in the CRDS-format distortion reference file
        coord_transform = None
        if self.runStep['astrometric']:
            with AsdfFile.open(self.params['Reffiles']['astrometric']) as dist_file:
                coord_transform = dist_file.tree['model']
        else:
            coord_transform = self.simple_coord_transform()
                
        # Using the requested RA, Dec of the reference pixel, along with the
        # V2, V3 of the reference pixel, and the requested roll angle of the telescope
        # create a matrix that can be used to translate between V2, V3 and RA, Dec
        # for any pixel
        # v2, v3 need to be in arcsec, and RA, Dec, and roll all need to be in degrees
        # attitude_matrix = rotations.attitude(self.refpix_pos['v2'], self.refpix_pos['v3'], self.ra, self.dec, self.params['Telescope']["rotation"])
        attitude_matrix = self.getAttitudeMatrix()

        # If an index column is present use that, otherwise
        # create one
        if 'index' in galaxylist.colnames:
            indexes = galaxylist['indexes']
        else:
            indexes = np.arange(1, len(galaxylist['radius']) + 1)
        # Make sure there is no 0th object
        if np.min(indexes) == 0:
            indexes += 1
        # Increment the index numbers so that these
        # sources don't overlap with any others already
        # in place
        if np.min(indexes) <= self.maxindex:
            indexes += self.maxindex
        self.maxindex = np.max(indexes)
        print("After galaxies, max index is {}".format(self.maxindex))

        # Loop over galaxy sources
        for index, source in zip(indexes, galaxylist):

            # If galaxy radii are given in units of arcseconds, translate to pixels
            if radiusflag == False:
                source['radius'] /= self.pixscale[0]

            # how many pixels beyond the nominal subarray edges can a source be located and
            # still have it fall partially on the subarray? Galaxy stamps are nominally set to
            # have a length and width equal to 100 times the requested radius.
            edgex = source['radius'] * 100 / 2 - 1
            edgey = source['radius'] * 100 / 2 - 1

            # reset the field of view limits for the size of the current stamp image
            outminy = miny - edgey
            outmaxy = maxy + edgey
            outminx = minx - edgex
            outmaxx = maxx + edgex

            try:
                entry0 = float(source['x_or_RA'])
                entry1 = float(source['y_or_Dec'])
                if not pixelflag:
                    ra_str, dec_str = self.makePos(entry0, entry1)
                    ra = entry0
                    dec = entry1
            except:
                # if inputs can't be converted to floats, then
                # assume we have RA/Dec strings. Convert to floats.
                ra_str = source['x_or_RA']
                dec_str = source['y_or_Dec']
                ra, dec = self.parseRADec(ra_str, dec_str)

            # case where point source list entries are given with RA and Dec
            if not pixelflag:

                # if distortion is to be included
                #if self.runStep['astrometric']:

                # Call is the same regardless of whether distortion reference file
                # is given or not
                pixelx, pixely = self.RADecToXY_astrometric(ra, dec, attitude_matrix, coord_transform)

                #else:
                #    # No distortion. Fall back to "manual" calculations
                #    pixelx, pixely = self.RADecToXY_manual(ra, dec)

            else:
                # case where the point source list entry locations are given in units of pixels
                # In this case we have the source position, and RA/Dec are calculated only so
                # they can be written out into the output source list file.

                # Assume that the input x and y values are coordinate values
                # WITHIN THE SPECIFIED SUBARRAY. So for example, a source in the file
                # at 0, 0 when you are making a SUB160 ramp will fall on the lower left
                # corner of the SUB160 subarray, NOT the lower left corner of the full
                # frame.

                pixelx = entry0
                pixely = entry1

                ra, dec, ra_str, dec_str = self.XYToRADec(pixelx, pixely, attitude_matrix, coord_transform)

            # only keep the source if the peak will fall within the subarray
            if pixely > outminy and pixely < outmaxy and pixelx > outminx and pixelx < outmaxx:

                pixelv2, pixelv3 = rotations.getv2v3(attitude_matrix, ra, dec)
                entry = [index, pixelx, pixely, ra_str, dec_str, ra, dec, pixelv2, pixelv3, source['radius'], source['ellipticity'], source['pos_angle'], source['sersic_index']]

                # Now look at the input magnitude of the point source
                # append the mag and pixel position to the list of ra, dec
                mag = float(source['magnitude'])
                entry.append(mag)

                # translate magnitudes to countrate
                # scale = 10.**(0.4*(15.0-mag))

                # get the countrate that corresponds to a 15th magnitude star for this filter
                # if self.params['Readout']['pupil'][0].upper() == 'F':
                #    usefilt = 'pupil'
                # else:
                #    usefilt = 'filter'
                # cval = self.countvalues[self.params['Readout'][usefilt]]

                # DEAL WITH THIS LATER, ONCE PYSYNPHOT IS INCLUDED WITH PIPELINE DIST?
                # if cval == 0:
                #    print("Countrate value for {} is zero in {}.".format(self.params['Readout'][usefilt], self.parameters['phot_file']))
                #    print("Eventually attempting to calculate value using pysynphot.")
                #    print("but pysynphot is not present in jwst build 6, so pushing off to later...")
                #    sys.exit()
                #    cval = self.findCountrate(self.params['Readout'][usefilt])

                # translate to counts in single frame at requested array size
                # framecounts = scale*cval*self.frametime
                # rate = scale*cval

                # Convert magnitudes to countrate (ADU/sec) and counts per frame
                rate = self.mag_to_countrate(magsystem, mag, photfnu=self.photfnu, photflam=self.photflam)
                framecounts = rate * self.frametime

                # add the countrate and the counts per frame to pointSourceList
                # since they will be used in future calculations
                entry.append(rate)
                entry.append(framecounts)

                # add the good point source, including location and counts, to the pointSourceList
                filteredList.add_row(entry)

        # Write the results to a file
        self.n_galaxies = len(filteredList)
        print(("Number of galaxies found within the requested aperture: {}"
               .format(self.n_galaxies)))

        if self.n_galaxies == 0:
            if self.n_pointsources == 0:
                raise ValueError('No point sources or galaxies found in input catalog; empty seed image would be created.')

        filteredList.meta['comments'] = ["Field center (degrees): %13.8f %14.8f y axis rotation angle (degrees): %f  image size: %4.4d %4.4d\n" % (self.ra,self.dec,self.params['Telescope']['rotation'],nx,ny)]
        filteredOut = self.basename + '_galaxySources.list'
        filteredList.write(filteredOut, format='ascii', overwrite=True)
        return filteredList

    def create_galaxy(self, radius, ellipticity, sersic, posang, totalcounts):
        # given relevent parameters, create a model sersic image with a given radius, eccentricity,
        # position angle, and total counts.

        # create the grid of pixels
        meshmax = np.min([np.int(self.ffsize * self.coord_adjust['y']), radius * 100.])
        x, y = np.meshgrid(np.arange(meshmax), np.arange(meshmax))

        # Center the galaxy in the array
        xc = meshmax / 2
        yc = meshmax / 2

        # Create model
        mod = Sersic2D(amplitude=1, r_eff=radius, n=sersic, x_0=xc, y_0=yc,
                       ellip=ellipticity, theta=posang)

        # Create instance of model
        img = mod(x, y)

        # Check to see if you've cropped too small and there is still significant signal
        # at the edges
        mxedge = np.max(np.array([np.max(img[:,-1]),np.max(img[:,0]),np.max(img[0,:]),np.max(img[-1,:])]))
        if mxedge > 0.001:
            print('Too small!')

        # Scale such that the total number of counts in the galaxy matches the input
        summedcounts = np.sum(img)
        if summedcounts == 0:
            print('in create_galaxy: ', radius, ellipticity, sersic, posang, totalcounts)
        factor = totalcounts / summedcounts
        img = img * factor

        # Crop image down such that it contains 99.95% of the total signal
        img = self.crop_galaxy_stamp(img,0.9995)
        return img

    def crop_galaxy_stamp(self, stamp, threshold):
        """Crop an input stamp image containing a galaxy to a size that
        contains only threshold times the total signal. This is an
        attempt to speed up the simulator a bit, since the galaxy stamp
        images are often very large. Note that galaxy stamp images being
        fed into this function are currently always square.

        Arguments:
        ----------
        stamp -- 2D stamp image of galaxy
        threshold -- fraction of total flux to keep in the cropped image
                      (e.g. 0.999 = 99.9%)

        Returns:
        --------
        cropped image
        """
        totsignal = np.sum(stamp)
        yd, xd = stamp.shape
        mid = np.int(xd / 2)
        for rad in range(mid):
            signal = np.sum(stamp[mid - rad:mid + rad + 1, mid - rad:mid + rad + 1]) / totsignal
            if signal >= threshold:
                return stamp[mid - rad:mid + rad + 1, mid - rad:mid + rad + 1]
        # If we make it all the way through the stamp without
        # hitting the threshold, then return the full stamp image
        return stamp

    def makeGalaxyImage(self, file):
        # Using the entries in the 'simSignals' 'galaxyList' file, create a countrate image
        # of model galaxies (sersic profile)

        # Read in the list of galaxies (positions and magnitides)
        glist, pixflag, radflag, magsys = self.readGalaxyFile(file)
        if pixflag:
            print("Galaxy list input positions assumed to be in units of pixels.")
        else:
            print("Galaxy list input positions assumed to be in units of RA and Dec.")

        if radflag:
            print("Galaxy list input radii assumed to be in units of pixels.")
        else:
            print("Galaxy list input radii assumed to be in units of arcsec.")

        # Extract and save only the entries which will land (fully or partially) on the
        # aperture of the output
        galaxylist = self.filterGalaxyList(glist, pixflag, radflag, magsys)

        # galaxylist is a table with columns:
        # 'pixelx', 'pixely', 'RA', 'Dec', 'RA_degrees', 'Dec_degrees', 'radius', 'ellipticity', 'pos_angle', 'sersic_index', 'magnitude', 'countrate_e/s', 'counts_per_frame_e'

        # final output image
        origyd, origxd = self.nominal_dims
        # origyd, origxd = self.dark.data[0, 0, :, :].shape
        yd = origyd
        xd = origxd

        # expand if a grism source image is being made
        xfact = 1
        yfact = 1
        if self.params['Output']['grism_source_image']:
            # xfact = self.grism_direct_factor
            # yfact = self.grism_direct_factor
            # elif
            yd = np.int(origyd * self.coord_adjust['y'])
            xd = np.int(origxd * self.coord_adjust['x'])

        # create the final galaxy countrate image
        galimage = np.zeros((yd, xd))
        dims = galimage.shape

        # Create corresponding segmentation map
        segmentation = segmap.SegMap()
        segmentation.xdim = xd
        segmentation.ydim = yd
        segmentation.initialize_map()

        # Adjust the coordinate system of the galaxy list if working with a grism direct image output
        deltax = 0
        deltay = 0
        if self.params['Output']['grism_source_image']:
            deltax = np.int((dims[1] - origxd) / 2)
            deltay = np.int((dims[0] - origyd) / 2)

        # create attitude matrix so we can calculate the North->V3 angle for
        # each galaxy
        attitude_matrix = self.getAttitudeMatrix()

        # For each entry, create an image, and place it onto the final output image
        for entry in galaxylist:

            # Get position angle in the correct units. Inputs for each
            # source are degrees east of north. So we need to find the
            # angle between north and V3, and then the angle between
            # V3 and the y-axis on the detector. The former can be found
            # using rotations.posang(attitude_matrix, v2, v3). The latter
            # is just V3SciYAngle in the SIAF (I think???)
            # v3SciYAng is measured in degrees, from V3 towards the Y axis,
            # measured from V3 towards V2.
            north_to_east_V3ang = rotations.posangle(attitude_matrix, entry['V2'], entry['V3'])
            # xposang = (0-self.v3scixang) - (north_to_east_V3ang - entry['pos_angle'])
            xposang = 0. - (self.v3scixang - north_to_east_V3ang + self.local_roll - entry['pos_angle'] + 90. + self.params['Telescope']['rotation'])

            # First create the galaxy image
            stamp = self.create_galaxy(entry['radius'], entry['ellipticity'], entry['sersic_index'], xposang*np.pi/180., entry['counts_per_frame_e'])

            # Source location in the padded seed image
            xpos = entry['pixelx'] + deltax
            ypos = entry['pixely'] + deltay
            nx = math.floor(xpos)
            ny = math.floor(ypos)

            #-------FUTURE IMPROVEMENT-------------
            # Identify the PSF that needs to be convolved with the stamp
            # First translate the given pixel coordinates, which are for the subarray
            # being simulated, into full frame coordinates
            #fullframex_equiv, fullframey_equiv = utils.subarray_to_full_coords(entry['pixelx'], entry['pixely'])

            # Using the full frame coordinates, choose/interpolate to find the most
            # appropriate PSF for the galaxy location.
            #psf = self.make_psf_model(fullframex_equiv, fullframey_equiv)
            #-------FUTURE IMPROVEMENT-------------

            # Load PSF model of the appropriate PSF given the source location
            psf_obj = PSF(entry['pixelx'], entry['pixely'], self.psfname,
                          interval=self.params['simSignals']['psfpixfrac'], oversampling=1)

            # Evaluate the EPSF model. Place into an array just large enough for it
            # as the next step will be convolution with the galaxy stamp image
            psf_image = psf_obj.minimal_psf_evaluation()

            # convolve the galaxy with the instrument PSF
            stamp = s1.fftconvolve(stamp, psf_image, mode='same')

            # Now add the stamp to the main image
            # Extract the appropriate subarray from the galaxy image if necessary
            galdims = stamp.shape

            # Calculate the coordinate limits of the aperture/PSF stamp overlap
            galydim, galxdim = stamp.shape
            (i1, i2, j1, j2, k1, k2, l1, l2) = self.cropped_coords(xpos, ypos, galxdim, galydim,
                                                                   dims[1], dims[0])

            if ((j2 > j1) and (i2 > i1) and (l2 > l1) and (k2 > k1) and (j1 < dims[0]) and (i1 < dims[0])):
                galimage[j1:j2, i1:i2] = galimage[j1:j2, i1:i2] + stamp[l1:l2, k1:k2]
                # Divide readnoise by 100 sec, which is a 10 group RAPID ramp?
                noiseval = self.single_ron / 100. + self.params['simSignals']['bkgdrate']
                if self.params['Inst']['mode'].lower() in ['wfss','ts_wfss']:
                    noiseval += self.grism_background
                segmentation.add_object_noise(stamp[l1:l2, k1:k2], j1, i1, entry['index'], noiseval)

            else:
                pass
                #print("Source located entirely outside the field of view. Skipping.")

        return galimage, segmentation.segmap

    def getExtendedSourceList(self, filename):
        # read in the list of point sources to add, and adjust the
        # provided positions for astrometric distortion

        extSourceList = Table(names=('index', 'pixelx', 'pixely', 'RA', 'Dec',
                                     'RA_degrees', 'Dec_degrees', 'magnitude',
                                     'countrate_e/s', 'counts_per_frame_e'),
                              dtype=('i', 'f', 'f', 'S14', 'S14', 'f', 'f', 'f', 'f', 'f'))

        try:
            lines, pixelflag, magsys = self.readPointSourceFile(filename)
            if pixelflag:
                print("Extended source list input positions assumed to be in units of pixels.")
            else:
                print("Extended list input positions assumed to be in units of RA and Dec.")
        except:
            raise FileNotFoundError("WARNING: Unable to open the extended source list file {}".format(filename))

        # File to save adjusted point source locations
        eoutcat = self.params['Output']['file'][0:-5] + '_extendedsources.list'
        eslist = open(eoutcat, 'w')

        dtor = math.radians(1.)
        nx = (self.subarray_bounds[2] - self.subarray_bounds[0]) + 1
        ny = (self.subarray_bounds[3] - self.subarray_bounds[1]) + 1
        xc = (self.subarray_bounds[2] + self.subarray_bounds[0]) / 2.
        yc = (self.subarray_bounds[3] + self.subarray_bounds[1]) / 2.

        # Location of the subarray's reference pixel.
        xrefpix = self.refpix_pos['x']
        yrefpix = self.refpix_pos['y']

        # center positions, sub-array sizes in pixels
        # now offset the field center to array center for astrometric distortion corrections
        coord_transform = None
        if self.runStep['astrometric']:
            # Read in the CRDS-format distortion reference file
            with AsdfFile.open(self.params['Reffiles']['astrometric']) as dist_file:
                coord_transform = dist_file.tree['model']
        else:
            coord_transform = self.simple_coord_transform()
                
        # Using the requested RA, Dec of the reference pixel, along with the
        # V2, V3 of the reference pixel, and the requested roll angle of the telescope
        # create a matrix that can be used to translate between V2, V3 and RA, Dec
        # for any pixel.
        # v2, v3 need to be in arcsec, and RA, Dec, and roll all need to be in degrees
        attitude_matrix = self.getAttitudeMatrix()

        # Write out the RA and Dec of the field center to the output file
        # Also write out column headers to prepare for source list
        eslist.write(("# Field center (degrees): %13.8f %14.8f y axis rotation angle "
                      "(degrees): %f  image size: %4.4d %4.4d\n" %
                      (self.ra, self.dec, self.params['Telescope']['rotation'], nx, ny)))
        eslist.write('# \n')
        eslist.write(("#    Index   RA_(hh:mm:ss)   DEC_(dd:mm:ss)   RA_degrees      "
                      "DEC_degrees     pixel_x   pixel_y    magnitude   counts/sec    counts/frame\n"))

        # Add an index column if not present
        if 'index' in lines.colnames:
            pass
        else:
            print('No extended object catalog index numbers. Adding to output: {}.'.format(eoutcat))
            indexes = np.arange(1, len(lines['filename']) + 1)
            lines['index'] = indexes

        # Make sure there is no 0th source
        if np.min(indexes) == 0:
            indexes += 1
        # Make sure the index numbers don't overlap with
        # sources that are already added
        if np.min(indexes) <= self.maxindex:
            indexes += self.maxindex
        self.maxindex = np.max(indexes)
        print("After extended sources, max index is {}".format(self.maxindex))

        #Loop over input lines in the source list
        all_stamps = []
        for indexnum, values in zip(indexes, lines):
            try:
            #line below (if 1>0) used to keep the block of code below at correct indent for the try: above
            #the try: is commented out for code testing.
            #if 1>0:
                try:
                    entry0 = float(values['x_or_RA'])
                    entry1 = float(values['y_or_Dec'])

                    if not pixelflag:
                        ra_str, dec_str = self.makePos(entry0, entry1)
                        ra = entry0
                        dec = entry1
                except:
                    # if inputs can't be converted to floats, then
                    # assume we have RA/Dec strings. Convert to floats.
                    ra_str = values['x_or_RA']
                    dec_str = values['y_or_Dec']
                    ra, dec = self.parseRADec(ra_str, dec_str)

                # Case where point source list entries are given with RA and Dec
                if not pixelflag:

                    # If distortion is to be included - either with or without the full set of coordinate
                    # translation coefficients
                    #if self.runStep['astrometric']:
                    # Same function call regardless of whether distortion file is provided or not
                    pixelx, pixely = self.RADecToXY_astrometric(ra, dec, attitude_matrix, coord_transform)
                    #else:
                    #    # No distortion at all - "manual mode"
                    #    pixelx, pixely = self.RADecToXY_manual(ra, dec)

                else:
                    # Case where the point source list entry locations are given in units of pixels
                    # In this case we have the source position, and RA/Dec are calculated only so
                    # they can be written out into the output source list file.

                    # Assume that the input x and y values are coordinate values
                    # WITHIN THE SPECIFIED SUBARRAY. So for example, a source in the file
                    # at 0, 0 when you are making a SUB160 ramp will fall on the lower left
                    # corner of the SUB160 subarray, NOT the lower left corner of the full
                    # frame.

                    pixelx = entry0
                    pixely = entry1

                    ra, dec, ra_str, dec_str = self.XYToRADec(pixelx, pixely, attitude_matrix,
                                                              coord_transform)

                # Get the input magnitude
                try:
                    mag = float(values['magnitude'])
                except:
                    mag = None

                # Now find out how large the extended source image is, so we
                # know if all, part, or none of it will fall in the field of view
                ext_stamp = fits.getdata(values['filename'])
                if len(ext_stamp.shape) != 2:
                    ext_stamp = fits.getdata(values['filename'], 1)

                eshape = np.array(ext_stamp.shape)
                if len(eshape) == 2:
                    edgey, edgex = eshape / 2
                else:
                    print(("WARNING, extended source image {} is not 2D! "
                           "Not sure how to proceed. Quitting.".format(values['filename'])))
                    sys.exit()

                #Define the min and max source locations (in pixels) that fall onto the subarray
                #Inlude the effects of a requested grism_direct image, and also keep sources that
                #will only partially fall on the subarray
                #pixel coords here can still be negative and kept if the grism image is being made

                #First, coord limits for just the subarray

                miny = 0
                maxy = self.subarray_bounds[3] - self.subarray_bounds[1]
                minx = 0
                maxx = self.subarray_bounds[2] - self.subarray_bounds[0]

                # Expand the limits if a grism direct image is being made
                if self.params['Output']['grism_source_image'] == True:
                    extrapixy = np.int((maxy + 1)/2 * (self.coord_adjust['y'] - 1.))
                    miny -= extrapixy
                    maxy += extrapixy
                    extrapixx = np.int((maxx + 1)/2 * (self.coord_adjust['x'] - 1.))
                    minx -= extrapixx
                    maxx += extrapixx

                # Now, expand the dimensions again to include point sources that fall only partially on the
                # subarray
                miny -= edgey
                maxy += edgey
                minx -= edgex
                maxx += edgex

                # Keep only sources within the appropriate bounds
                if pixely > miny and pixely < maxy and pixelx > minx and pixelx < maxx:

                    #set up an entry for the output table
                    entry = [indexnum, pixelx, pixely, ra_str, dec_str, ra, dec, mag]

                    # save the stamp image after normalizing to a total signal of 1.
                    norm_factor = np.sum(ext_stamp)
                    ext_stamp /= norm_factor
                    all_stamps.append(ext_stamp)

                    # If a magnitude is given then adjust the countrate to match it
                    if mag is not None:
                        # translate magnitudes to countrate
                        # scale = 10.**(0.4*(15.0-mag))
                        #
                        # get the countrate that corresponds to a 15th magnitude star for this filter
                        # if self.params['Readout']['pupil'][0].upper() == 'F':
                        #   usefilt = 'pupil'
                        # else:
                        #    usefilt = 'filter'
                        # cval = self.countvalues[self.params['Readout'][usefilt]]
                        #
                        # DEAL WITH THIS LATER, ONCE PYSYNPHOT IS INCLUDED WITH PIPELINE DIST?
                        # if cval == 0:
                        #    print("Countrate value for {} is zero in {}.".format(self.params['Readout'][usefilt], self.parameters['phot_file']))
                        #    print("Eventually attempting to calculate value using pysynphot.")
                        #    print("but pysynphot is not present in jwst build 6, so pushing off to later...")
                        #    sys.exit()
                        #    cval = self.findCountrate(self.params['Readout'][usefilt])

                        # translate to counts in single frame at requested array size
                        # framecounts = scale*cval*self.frametime
                        # countrate = scale*cval
                        # magwrite = mag

                        # Convert magnitudes to countrate (ADU/sec) and counts per frame
                        countrate = self.mag_to_countrate(magsys, mag, photfnu=self.photfnu, photflam=self.photflam)
                        framecounts = countrate * self.frametime
                        magwrite = mag

                    else:
                        # In this case, no magnitude is given in the extended input list
                        # Assume the input stamp image is in units of e/sec then.
                        print("No magnitude given for extended source in {}.".format(values['filename']))
                        print("Assuming the original file is in units of counts per sec.")
                        print("Multiplying original file values by 'extendedscale'.")
                        countrate = norm_factor * self.params['simSignals']['extendedscale']
                        framecounts = countrate*self.frametime
                        magwrite = 99.99999

                    # add the countrate and the counts per frame to pointSourceList
                    # since they will be used in future calculations
                    # entry.append(scale)
                    entry.append(countrate)
                    entry.append(framecounts)

                    # add the good point source, including location and counts, to the pointSourceList
                    # self.pointSourceList.append(entry)
                    extSourceList.add_row(entry)

                    #write out positions, distances, and counts to the output file
                    eslist.write("%i %s %s %14.8f %14.8f %9.3f %9.3f  %9.3f  %13.6e   %13.6e\n" % (indexnum, ra_str, dec_str, ra, dec, pixelx, pixely, magwrite, countrate, framecounts))
            except:
                #print("ERROR: bad point source line %s. Skipping." % (line))
                pass
        print("Number of extended sources found within the requested aperture: {}".format(len(extSourceList)))
        # close the output file
        eslist.close()

        # If no good point sources were found in the requested array, alert the user
        if len(extSourceList) < 1:
            print("Warning: no non-sidereal extended sources within the requested array.")
            print("The extended source image option is being turned off")

        return extSourceList, all_stamps

    def makeExtendedSourceImage(self, extSources, extStamps):
        dims = np.array(self.nominal_dims)

        # offset that needs to be applied to the x, y positions of the
        # source list to account for case where we make a point
        # source image that is extra-large, to be used as a grism
        # direct image
        newdimsx = np.int(dims[1] * self.coord_adjust['x'])
        newdimsy = np.int(dims[0] * self.coord_adjust['y'])
        deltax = self.coord_adjust['xoffset']
        deltay = self.coord_adjust['yoffset']
        dims = np.array([newdimsy, newdimsx])

        # create the empty image
        extimage = np.zeros((dims[0], dims[1]))

        # Create corresponding segmentation map
        segmentation = segmap.SegMap()
        segmentation.xdim = newdimsx
        segmentation.ydim = newdimsy
        segmentation.initialize_map()

        # Loop over the entries in the source list
        for entry, stamp in zip(extSources, extStamps):
            # adjust x, y position if the grism output image is requested
            xpos = entry['pixelx'] + deltax
            ypos = entry['pixely'] + deltay
            xoff = math.floor(xpos)
            yoff = math.floor(ypos)
            
            # desired counts per second in the source
            counts = entry['countrate_e/s'] # / self.frametime

            # If the stamp needs to be convolved with the NIRCam PSF,
            # locate the correct PSF file here and read it in
            if self.params['simSignals']['PSFConvolveExtended']:
            
                # Load the appropriate PSF file given the source location
                psf_obj = PSF(entry['pixelx'], entry['pixely'], self.psfname,
                              interval=self.params['simSignals']['psfpixfrac'], oversampling=1)
                
                # Evaluate the EPSF model. Place into an array just large enough for it
                # as the next step will be convolution with the galaxy stamp image
                psf_image = psf_obj.minimal_psf_evaluation()
                
                # Convolve the galaxy with the instrument PSF
                stamp = s1.fftconvolve(stamp, psf_image, mode='same')
  
            # Extract the appropriate subarray from the image if necessary
            # Assume that the brightest pixel corresponds to the peak of the source
            #psfdims = stamp.shape
            

            # Calculate the coordinate limits of the aperture/PSF stamp overlap
            stampydim, stampxdim = stamp.shape
            (i1, i2, j1, j2, k1, k2, l1, l2) = self.cropped_coords(xpos, ypos, stampxdim, stampydim,
                                                                   dims[1], dims[0])
                
            # Add stamp image to the extended source countrate image
            extimage[j1:j2, i1:i2] = extimage[j1:j2, i1:i2] + stamp[l1:l2, k1:k2] * counts
            # Divide readnoise by 100 sec, which is a 10 group RAPID ramp?
            noiseval = self.single_ron / 100. + self.params['simSignals']['bkgdrate']
            if self.params['Inst']['mode'].lower() in ['wfss','ts_wfss']:
                noiseval += self.grism_background

            #segmentation.add_object_noise(stamp[l1:l2, k1:k2]*counts, j1, i1, entry['index'], noiseval)
            indseg = self.seg_from_photutils(stamp[l1:l2, k1:k2] * counts, entry['index'], noiseval)
            segmentation.segmap[j1:j2, i1:i2] += indseg
        return extimage, segmentation.segmap

    def seg_from_photutils(self, image, number, noise):
        # Create a segmentation map for the input image
        # using photutils. In this case, the input noise
        # represents the single frame noise for the appropriate
        # observing mode
        map = detect_sources(image, noise * 3., 8).data + number
        return map

    def makeFilterTable(self):
        # Create the table that contains the possible filter list, quantum yields, and countrates for a
        # star with vega magnitude of 15 in each filter. Do this by reading in phot_file
        # listed in the parameter file.

        # FUTURE WORK: If the countrates are left as 0, then pysynphot will
        # be used to calculate them later
        try:
            cvals_tab = ascii.read(self.params['Reffiles']['phot'])
            instrumentfilternames = cvals_tab['filter'].data
            stringcountrates = cvals_tab['countrate_for_vegamag15'].data
            instrumentmag15countrates = [float(s) for s in stringcountrates]
            strinstrumentqy = cvals_tab['quantum_yield'].data
            qy = [float(s) for s in strinstrumentqy]
            self.countvalues = dict(zip(instrumentfilternames, instrumentmag15countrates))
            self.qydict = dict(zip(instrumentfilternames, qy))

        except:
            raise IOError("WARNING: Unable to read in {}.".format(self.params['Reffiles']['phot']))

    def readParameterFile(self):
        # read in the parameter file
        try:
            with open(self.paramfile, 'r') as infile:
                self.params = yaml.load(infile)
        except:
            raise IOError("WARNING: unable to open {}".format(self.paramfile))

    def checkParams(self):
        """Check input parameters for expected datatypes, values"""
        # Check instrument name
        if self.params['Inst']['instrument'].lower() not in INST_LIST:
            raise NotImplementedError("WARNING: {} instrument not implemented within ramp simulator")

        # Check entered mode:
        possibleModes = MODES[self.params['Inst']['instrument'].lower()]
        self.params['Inst']['mode'] = self.params['Inst']['mode'].lower()
        if self.params['Inst']['mode'] in possibleModes:
            pass
        else:
            raise ValueError(("WARNING: unrecognized mode {} for {}. Must be one of: {}"
                   .format(self.params['Inst']['mode'],
                           self.params['Inst']['instrument'],possibleModes)))

        # Check telescope tracking entry
        self.params['Telescope']['tracking'] = self.params['Telescope']['tracking'].lower()
        if self.params['Telescope']['tracking'] not in TRACKING_LIST:
            raise ValueError(("WARNING: telescope tracking set to {}, but must be one "
                              "of {}.".format(self.params['Telescope']['tracking'],
                                              TRACKING_LIST)))

        # Non-sidereal WFSS observations are not yet supported
        if self.params['Telescope']['tracking'] == 'non-sidereal' and \
           self.params['Inst']['mode'] in ['wfss','ts_wfss']:
            raise ValueError(("WARNING: wfss observations with non-sidereal "
                              "targets not yet supported."))
        
        # Set nframe and nskip according to the values in the
        # readout pattern definition file
        self.read_pattern_check()

        #Make sure that the requested number of groups is
        #less than or equal to the maximum allowed.
        #For full frame science operations, ngroup is going
        #to be limited to 10 for all readout patterns
        #except for the DEEP patterns, which can go to 20.
        #match = self.readpatterns['name'] == self.params['Readout']['readpatt'].upper()
        #if sum(match) == 1:
        #    maxgroups = self.readpatterns['maxgroups'].data[match][0]
        # if sum(match) == 0:
        #    print("Unrecognized readout pattern {}. Assuming a maximum allowed number of groups of 10.".format(self.params['Readout']['readpatt']))
        #    maxgroups = 10

        # if (self.params['Readout']['ngroup'] > maxgroups):
        #    print("WARNING: {} is limited to a maximum of {} groups. Proceeding with ngroup = {}.".format(self.params['Readout']['readpatt'], maxgroups, maxgroups))
        #    self.params['Readout']['readpatt'] = maxgroups


        # check for entries in the parameter file that are None or blank,
        # indicating the step should be skipped. Create a dictionary of steps
        # and populate with True or False
        self.runStep = {}
        self.runStep['pixelflat'] = self.checkRunStep(self.params['Reffiles']['pixelflat'])
        self.runStep['illuminationflat'] = self.checkRunStep(self.params['Reffiles']['illumflat'])
        self.runStep['astrometric'] = self.checkRunStep(self.params['Reffiles']['astrometric'])
        self.runStep['distortion_coeffs'] = self.checkRunStep(self.params['Reffiles']['distortion_coeffs'])
        self.runStep['ipc'] = self.checkRunStep(self.params['Reffiles']['ipc'])
        self.runStep['crosstalk'] = self.checkRunStep(self.params['Reffiles']['crosstalk'])
        self.runStep['occult'] = self.checkRunStep(self.params['Reffiles']['occult'])
        self.runStep['pointsource'] = self.checkRunStep(self.params['simSignals']['pointsource'])
        self.runStep['galaxies'] = self.checkRunStep(self.params['simSignals']['galaxyListFile'])
        self.runStep['extendedsource'] = self.checkRunStep(self.params['simSignals']['extended'])
        self.runStep['movingTargets'] = self.checkRunStep(self.params['simSignals']['movingTargetList'])
        self.runStep['movingTargetsSersic'] = self.checkRunStep(self.params['simSignals']['movingTargetSersic'])
        self.runStep['movingTargetsExtended'] = self.checkRunStep(self.params['simSignals']['movingTargetExtended'])
        self.runStep['MT_tracking'] = self.checkRunStep(self.params['simSignals']['movingTargetToTrack'])
        self.runStep['zodiacal'] = self.checkRunStep(self.params['simSignals']['zodiacal'])
        self.runStep['scattered'] = self.checkRunStep(self.params['simSignals']['scattered'])
        # self.runStep['fwpw'] = self.checkRunStep(self.params['Reffiles']['filtpupilcombo'])
        self.runStep['pixelAreaMap'] = self.checkRunStep(self.params['Reffiles']['pixelAreaMap'])

        # Notify user if no catalogs are provided
        if self.params['simSignals']['pointsource'] == 'None':
            print('No point source catalog provided in yaml file.')

        if self.params['simSignals']['galaxyListFile'] == 'None':
            print('No galaxy catalog provided in yaml file.')

        # create table that will contain filters/quantum yield/and vegamag=15 countrates
        # self.makeFilterTable()

        # Read in list of zeropoints/photflam/photfnu
        self.zps = ascii.read(self.params['Reffiles']['flux_cal'])

        # Determine the instrument module and detector from the aperture name
        aper_name = self.params['Readout']['array_name']
        try:
            detector = self.subdict[self.subdict['AperName'] == aper_name]['Detector'][0]
            module = detector[0]
        except IndexError:
            raise ValueError('Unable to determine the detector/module in aperture {}'.format(aper_name))

        # In the future we expect zeropoints to be detector dependent, as they currently
        # are for FGS. So if we are working with NIRCAM or NIRISS, manually add a Detector key
        # to the dictionary as a placeholder.
        if self.params["Inst"]["instrument"].lower() in ["nircam", "niriss"]:
            self.zps = self.add_detector_to_zeropoints(detector)
        
        # make sure the requested filter is allowed. For imaging, all filters are allowed.
        # In the future, other modes will be more restrictive
        if self.params['Readout']['pupil'][0].upper() == 'F':
            usefilt = 'pupil'
        else:
            usefilt = 'filter'

        # If instrument is FGS, then force filter to be 'NA' for the purposes
        # of constructing the correct PSF input path name. Then change to be
        # the DMS-required "N/A" when outputs are saved
        if self.params['Inst']['instrument'].lower() == 'fgs':
            self.params['Readout']['filter'] = 'NA'
            self.params['Readout']['pupil'] = 'NA'
            
        if self.params['Readout'][usefilt] not in self.zps['Filter']:
            raise ValueError(("WARNING: requested filter {} is not in the list of "
                   "possible filters.".format(self.params['Readout'][usefilt])))

        # Get the photflambda and photfnu values that go with
        # the filter
        mtch = ((self.zps['Detector'] == detector) &
                (self.zps['Filter'] == self.params['Readout'][usefilt]) &
                (self.zps['Module'] == module))
        self.vegazeropoint=self.zps['VEGAMAG'][mtch][0]
        self.photflam = self.zps['PHOTFLAM'][mtch][0]
        self.photfnu = self.zps['PHOTFNU'][mtch][0]
        self.pivot = self.zps['Pivot_wave'][mtch][0]

        #PSF: generate the name of the PSF file to use
        #if the psf path has been left blank or set to 'None'
        #then assume the user does not want to add point sources
        if self.params['simSignals']['psfpath'] is not None:
            if self.params['simSignals']['psfpath'][-1] != '/':
                self.params['simSignals']['psfpath']=self.params['simSignals']['psfpath'] + '/'

            wfe = self.params['simSignals']['psfwfe']
            if wfe not in WFE_OPTIONS:
                raise ValueError(("WARNING: invalid wavefront error (psfwfe) input: {}"
                                  "psfwfe must be one of: {}".format(wfe, WFE_OPTIONS)))
            wfegroup = self.params['simSignals']['psfwfegroup']
            if wfegroup not in WFEGROUP_OPTIONS:
                raise ValueError(("WARNING: invalid wavefront group (psfwfegroup) "
                                  "value: {}. psfwfegroup must be one of: {}"
                                  .format(wfegroup, WFEGROUP_OPTIONS)))
            basename = self.params['simSignals']['psfbasename'] + '_'
            if wfe == 0:
                psfname = basename + self.params['simSignals'][usefilt].lower() + '_zero'
                self.params['simSignals']['psfpath'] = self.params['simSignals']['psfpath'] + \
                                                       self.params['simSignals'][usefilt].lower() + \
                                                       '/zero/'
            else:
                if self.params['Inst']['instrument'].lower() != 'fgs':
                    psfname = '{}{}_x{}_y{}_{}_{}_{}'.format(basename, detector,
                                                             'psfxpos', 'psfypos',
                                                             self.params['Readout'][usefilt].lower(),
                                                             str(wfe), str(wfegroup))
                else:
                    psfname = '{}{}_x{}_y{}_{}_{}'.format(basename, detector,
                                                             'psfxpos', 'psfypos',
                                                             str(wfe), str(wfegroup))

                psfname = psfname.replace('psfxpos', '1024')
                psfname = psfname.replace('psfypos', '1024')
                if self.params['Inst']['instrument'].lower() != 'fgs':
                    pathaddition = "{}/{}/{}".format(detector,
                                                     self.params['Readout'][usefilt].lower(),
                                                     str(wfe))
                else:
                    pathaddition = "{}/{}".format(detector, str(wfe))
                    
                self.params['simSignals']['psfpath'] = os.path.join(self.params['simSignals']['psfpath'], pathaddition)
                self.psfname = os.path.join(self.params['simSignals']['psfpath'], psfname)
        else:
            # case where psfPath is None. In this case, create a PSF on the fly to use
            # for adding sources
            raise IOError("WARNING: no PSF library path given.")
            #print("update this to include a call to WebbPSF????????")
            #self.psfimage = np.zeros((5, 5), dtype=np.float32)
            #sum1 = 0
            #for i in range(5):
            #    for j in range(5):
            #        self.psfimage[i, j] = (0.02**(abs((i - 2))) * (0.02**abs(j - 2)))
            #        sum1 = sum1 + self.psfimage[i, j]
            #self.psfimage = self.psfimage / sum1
            #self.psfname = None

        # ASTROMETRY
        # Read in the distortion coefficients file if present. These will provide a more exact
        # transform from RA, Dec to x, y than the astrometric distortion reference file above.
        # The file above can be off by ~20 pixels in the corners of the array. This file will give
        # exact answers
        if self.runStep['distortion_coeffs'] == True:
            if os.path.isfile(self.params['Reffiles']['distortion_coeffs']):
                distortionTable = ascii.read(self.params['Reffiles']['distortion_coeffs'], header_start=1, format='csv')
            else:
                raise FileNotFoundError(("WARNING: Input distortion coefficients file {} "
                                         "does not exist."
                                         .format(self.params['Reffiles']['distortion_coeffs'])))

            # read in coefficients for the forward 'science' to 'ideal' coordinate transformation.
            # 'science' is in units of distorted pixels, while 'ideal' is the undistorted
            # angular distance from the reference pixel
            ap_name = self.params['Readout']['array_name']

            self.x_sci2idl, self.y_sci2idl, self.v2_ref, self.v3_ref, \
                self.parity, self.v3yang, self.xsciscale, self.ysciscale, \
                self.v3scixang = self.getDistortionCoefficients(distortionTable,
                                                                'science', 'ideal', ap_name)

            #Generate the coordinate transform for V2, V3 to 'ideal'
            siaf = ascii.read(self.params['Reffiles']['distortion_coeffs'], header_start=1, format='csv')

            match = siaf['AperName'] == ap_name
            if not np.any(match):
                raise ValueError("Aperture name {} not found in input CSV file.".format(ap_name))

            siaf_row = siaf[match]

            self.v2v32idlx, self.v2v32idly = read_siaf_table.\
                                             get_siaf_v2v3_transform(siaf_row,
                                                                     ap_name,
                                                                     to_system='ideal')

        #convert the input RA and Dec of the pointing position into floats
        #check to see if the inputs are in decimal units or hh:mm:ss strings
        try:
            self.ra = float(self.params['Telescope']['ra'])

            self.dec = float(self.params['Telescope']['dec'])
        except:
            self.ra, self.dec = self.parseRADec(self.params['Telescope']['ra'],
                                                self.params['Telescope']['dec'])

        if abs(self.dec) > 90. or self.ra < 0. or self.ra > 360. or \
           self.ra is None or self.dec is None:
            raise ValueError("WARNING: bad requested RA and Dec {} {}".format(self.ra, self.dec))

        # make sure the rotation angle is a float
        try:
            self.params['Telescope']["rotation"] = float(self.params['Telescope']["rotation"])
        except:
            print(("ERROR: bad rotation value {}, setting to zero."
                   .format(self.params['Telescope']["rotation"])))
            self.params['Telescope']["rotation"] = 0.

        # Set the background value if the high/medium/low settings
        # are used
        bkgdrate_options = ['high', 'medium', 'low']

        try:
            self.params['simSignals']['bkgdrate'] = float(self.params['simSignals']['bkgdrate'])
        except:
            if self.params['simSignals']['bkgdrate'].lower() in bkgdrate_options:
                print(("Calculating background rate using jwst_background "
                       "based on {} level".format(self.params['simSignals']['bkgdrate'])))

                # Find the appropriate filter throughput file
                if os.path.split(self.params['Reffiles']['filter_throughput'])[1] == 'placeholder.txt':
                    instrm = self.params['Inst']['instrument'].lower()
                    if instrm == 'nircam':
                        filter_file = ("{}_nircam_plus_ote_throughput_mod{}_sorted.txt"
                                       .format(self.params['Readout'][usefilt].upper(), module.lower()))
                    elif instrm == 'niriss':
                        filter_file = ("{}_niriss_throughput1.txt"
                                       .format(self.params['Readout'][usefilt].lower()))
                    elif instrm == 'fgs':
                        #det = self.params['Readout']['array_name'].split('_')[0]
                        filter_file = "{}_throughput_py.txt".format(detector.lower())
                    filt_dir = os.path.split(self.params['Reffiles']['filter_throughput'])[0]
                    filter_file = os.path.join(filt_dir, filter_file)

                else:
                    filter_file = self.params['Reffiles']['filter_throughput']

                print(("Using {} filter throughput file for background calculation."
                       .format(filter_file)))

                self.params['simSignals']['bkgdrate'] = \
                                self.calculate_background(self.ra,
                                                          self.dec,
                                                          filter_file,
                                                          level=self.params['simSignals']['bkgdrate'].lower())
                print('Background level set to: {}'.format(self.params['simSignals']['bkgdrate']))
            else:
                raise ValueError(("WARNING: unrecognized background rate value. "
                                  "Must be either a number or one of: {}"
                                  .format(bkgdrate_options)))

        #check that the various scaling factors are floats and within a reasonable range
        #self.params['cosmicRay']['scale'] = self.checkParamVal(self.params['cosmicRay']['scale'], 'cosmicRay', 0, 100, 1)
        self.params['simSignals']['extendedscale'] = self.checkParamVal(self.params['simSignals']['extendedscale'], 'extendedEmission', 0, 10000, 1)
        self.params['simSignals']['zodiscale'] = self.checkParamVal(self.params['simSignals']['zodiscale'], 'zodi', 0, 10000, 1)
        self.params['simSignals']['scatteredscale'] = self.checkParamVal(self.params['simSignals']['scatteredscale'], 'scatteredLight', 0, 10000, 1)

        # make sure the requested output format is an allowed value
        if self.params['Output']['format'] not in ALLOWEDOUTPUTFORMATS:
            raise ValueError(("WARNING: unsupported output format {} requested. "
                              "Possible options are {}.".format(self.params['Output']['format'],
                                                                ALLOWEDOUTPUTFORMATS)))

        # Entries for creating the grims input image
        if not isinstance(self.params['Output']['grism_source_image'], bool):
            if self.params['Output']['grism_source_image'].lower() == 'none':
                self.params['Output']['grism_source_image'] = False
            else:
                raise ValueError("WARNING: grism_source_image needs to be True or False")

        # Location of extended image on output array, pixel x, y values.
        try:
            self.params['simSignals']['extendedCenter'] = np.fromstring(self.params['simSignals']['extendedCenter'], dtype=int, sep=", ")
        except:
            raise RuntimeError(("WARNING: not able to parse the extendedCenter list {}. "
                                "It should be a comma-separated list of x and y pixel positions."
                                .format(self.params['simSignals']['extendedCenter'])))

        # check the output metadata, including visit and observation numbers, obs_id, etc
        #
        # kwchecks = ['program_number', 'visit_number', 'visit_group',
        #            'sequence_id', 'activity_id', 'exposure_number', 'observation_number', 'obs_id', 'visit_id']
        # for quality in kwchecks:
        #    try:
        #        self.params['Output'][quality] = str(self.params['Output'][quality])
        #    except:
        #        print("WARNING: unable to convert {} to string. This is required.".format(self.params['Output'][quality]))
        #        sys.exit()

    def checkRunStep(self, filename):
        # check to see if a filename exists in the parameter file.
        if ((len(filename) == 0) or (filename.lower() == 'none')):
            return False
        else:
            return True

    def read_pattern_check(self):
        # Check the readout pattern that's entered and set nframe and nskip
        # accordingly
        self.params['Readout']['readpatt'] = self.params['Readout']['readpatt'].upper()

        # Read in readout pattern definition file
        # and make sure the possible readout patterns are in upper case
        self.readpatterns = ascii.read(self.params['Reffiles']['readpattdefs'])
        self.readpatterns['name'] = [s.upper() for s in self.readpatterns['name']]

        # If the requested readout pattern is in the table of options,
        # then adopt the appropriate nframe and nskip
        if self.params['Readout']['readpatt'] in self.readpatterns['name']:
            mtch = self.params['Readout']['readpatt'] == self.readpatterns['name']
            self.params['Readout']['nframe'] = self.readpatterns['nframe'][mtch].data[0]
            self.params['Readout']['nskip'] = self.readpatterns['nskip'][mtch].data[0]
            print(('Requested readout pattern {} is valid. '
                  'Using the nframe = {} and nskip = {}'
                   .format(self.params['Readout']['readpatt'],
                           self.params['Readout']['nframe'],
                           self.params['Readout']['nskip'])))
        else:
            # If the read pattern is not present in the definition file
            # then quit.
            raise ValueError(("WARNING: the {} readout pattern is not defined in {}."
                              .format(self.params['Readout']['readpatt'],
                                      self.params['Reffiles']['readpattdefs'])))

    def filecheck(self):
        # Make sure the requested input files exist
        # For reference files, assume first that they are located in
        # the directory tree under the directory specified by the MIRAGE_DATA
        # environment variable. If not, assume the input is a full path
        # and check there.
        rlist = [['Reffiles', 'astrometric'],
                 ['Reffiles', 'distortion_coeffs']]
        plist = [['simSignals', 'psfpath']]
        ilist = [['simSignals', 'pointsource'],
                 ['simSignals', 'galaxyListFile'],
                 ['simSignals', 'extended'],
                 ['simSignals', 'movingTargetList'],
                 ['simSignals', 'movingTargetSersic'],
                 ['simSignals', 'movingTargetExtended'],
                 ['simSignals', 'movingTargetToTrack']]
        for ref in rlist:
            self.ref_check(ref)
        for path in plist:
            self.path_check(path)
        for inp in ilist:
            self.input_check(inp)

    def ref_check(self, rele):
        """
        Check for the existence of the input reference file
        Assume first that the file is in the directory tree
        specified by the MIRAGE_DATA environment variable.

        Parameters:
        -----------
        rele -- Tuple containing the nested keys that point
                to the refrence file of interest. These come
                from the yaml input file

        Reutrns:
        --------
        Nothing
        """
        rfile = self.params[rele[0]][rele[1]]
        if rfile.lower() != 'none':
            c1 = os.path.isfile(rfile)
            if not c1:
                raise FileNotFoundError(("WARNING: Unable to locate the {}, {} "
                                         "input file! Not present in {}"
                                         .format(rele[0], rele[1], rfile)))

    def path_check(self, p):
        """
        Check for the existence of the input path.
        Assume first that the path is in relation to
        the directory tree specified by the MIRAGE_DATA
        environment variable
        
        Parameters:
        -----------
        p -- Tuple containing the nested keys that point
             to a directory in self.params

        Returns:
        --------
        Nothing
        """
        pth = self.params[p[0]][p[1]]
        c1 = os.path.exists(pth)
        if not c1:
            raise NotADirectoryError(("WARNING: Unable to find the requested path "
                                      "{}. Not present in directory tree specified by "
                                      "the {} environment variable."
                                      .format(pth, self.env_var)))

    def input_check(self, inparam):
        # Check for the existence of the input file. In
        # this case we do not check the directory tree
        # specified by the MIRAGE_DATA environment variable.
        # This is intended primarily for user-generated inputs like
        # source catalogs
        ifile = self.params[inparam[0]][inparam[1]]
        if ifile.lower() != 'none':
            c = os.path.isfile(ifile)
            if not c:
                raise FileNotFoundError(("WARNING: Unable to locate {} Specified "
                                         "by the {}:{} field in the input yaml file."
                                         .format(ifile, inparam[0], inparam[1])))

    def checkParamVal(self, value, typ, vmin, vmax, default):
        # make sure the input value is a float and between min and max
        try:
            value = float(value)
        except:
            raise ValueError("WARNING: {} for {} is not a float.".format(value, typ))

        if ((value >= vmin) & (value <= vmax)):
            return value
        else:
            print(("ERROR: {} for {} is not within reasonable bounds. "
                   "Setting to {}".format(value, typ, default)))
            return default

    def readSubarrayDefinitionFile(self):
        # read in the file that contains a list of subarray names and positions on the detector

        try:
            self.subdict = ascii.read(self.params['Reffiles']['subarray_defs'], data_start=1, header_start=0)
        except:
            raise RuntimeError("Error: could not read in subarray definitions file: {}".format(self.params['Reffiles']['subarray_defs']))

    def getSubarrayBounds(self):
        # find the bounds of the requested subarray
        if self.params['Readout']['array_name'] in self.subdict['AperName']:
            mtch = self.params['Readout']['array_name'] == self.subdict['AperName']
            self.subarray_bounds = [self.subdict['xstart'].data[mtch][0], self.subdict['ystart'].data[mtch][0], self.subdict['xend'].data[mtch][0], self.subdict['yend'].data[mtch][0]]
            self.refpix_pos = {'x':self.subdict['refpix_x'].data[mtch][0], 'y':self.subdict['refpix_y'][mtch][0], 'v2':self.subdict['refpix_v2'].data[mtch][0], 'v3':self.subdict['refpix_v3'].data[mtch][0]}

            namps = self.subdict['num_amps'].data[mtch][0]
            if namps != 0:
                self.params['Readout']['namp'] = namps
            else:
                if ((self.params['Readout']['namp'] == 1) or
                    (self.params['Readout']['namp'] == 4)):
                    print(("CAUTION: Aperture {} can be used with either "
                           "a 1-amp".format(self.subdict['AperName'].data[mtch][0])))
                    print("or a 4-amp readout. The difference is a factor of 4 in")
                    print(("readout time. You have requested {} amps."
                           .format(self.params['Readout']['namp'])))
                else:
                    raise ValueError(("WARNING: {} requires the number of amps "
                                      "to be 1 or 4. You have requested {}."
                                      .format(self.params['Readout']['array_name'],
                                              self.params['Readout']['namp'])))
        else:
            raise ValueError(("WARNING: subarray name {} not found in the "
                              "subarray dictionary {}."
                              .format(self.params['Readout']['array_name'],
                                      self.params['Reffiles']['subarray_defs'])))

    def instrument_specific_dicts(self, instrument):
        # get instrument-specific values for things that
        # don't need to be in the parameter file

        # array size of a full frame image
        self.ffsize = FULL_ARRAY_SIZE[instrument]

        # pixel scale - return as a 2-element list, with pixscale for x and y.
        if instrument.lower() == 'nircam':
            filt = self.params['Readout']['filter']
            fnum = int(filt[1:4])
            if fnum < 230:
                channel = 'sw'
            else:
                channel = 'lw'
            self.pixscale = [PIXELSCALE[instrument][channel], PIXELSCALE[instrument][channel]]
        else:
            self.pixscale = [PIXELSCALE[instrument], PIXELSCALE[instrument]]

    def read_filter_throughput(self, file):
        '''Read in the ascii file containing the filter
        throughput curve'''
        tab = ascii.read(file)
        return tab['Wavelength_microns'].data, tab['Throughput'].data

    def calculate_background(self, ra, dec, ffile, level='medium'):
        '''Use the JWST background calculator to come up with
        an appropriate background level for the observation.
        Options for level include low, medium, high'''
        from jwst_backgrounds import jbt
        from astropy import units as u
        from astropy.units.equivalencies import si, cgs

        # Read in filter throughput file
        filt_wav, filt_thru = self.read_filter_throughput(ffile)

        # Get background information
        # Any wavelength will return the 2D array that includes
        # all wavelengths, so just use a dummy value of 2.5 microns
        bg = jbt.background(ra, dec, 2.5)

        # Now we need to loop over each day (in the background)
        # info, convolve the background curve with the filter
        # throughput curve, and then integrate. THEN, we can
        # calculate the low/medium/high values.
        bsigs = np.zeros(len(bg.bkg_data['total_bg'][:, 0]))
        for i in range(len(bg.bkg_data['total_bg'][:, 0])):
            back_wave = bg.bkg_data['wave_array']
            back_sig = bg.bkg_data['total_bg'][i, :]

            # Interpolate background to match filter wavelength grid
            bkgd_interp = np.interp(filt_wav, back_wave, back_sig)

            # Combine
            filt_bkgd = bkgd_interp * filt_thru

            # Integrate
            bsigs[i] = np.trapz(filt_bkgd, x=filt_wav)

        # Now sort and determine the low/medium/high levels
        x = np.sort(bsigs)
        y = np.arange(1, len(x) + 1) / len(x)

        if level.lower() == 'low':
            perc = 0.1
        elif level.lower() == 'medium':
            perc = 0.5
        elif level.lower() == 'high':
            perc = 0.9
        else:
            raise ValueError("Unrecognized background level string")

        # Interpolate to the requested level
        bval = np.interp(perc, y, x) * u.MJy / u.steradian

        # Convert from MJy/str to ADU/sec
        # then divide by area of pixel
        flambda = cgs.erg / si.angstrom / si.cm ** 2 / si.s
        #fnu = cgs.erg / si.Hz / si.cm ** 2 / si.s
        photflam = self.photflam * flambda
        #photnu = self.photfnu * fnu
        pivot = self.pivot * u.micron
        mjy = photflam.to(u.MJy, u.spectral_density(pivot))

        # Divide by pixel area in steradians to get
        # MJy/str per ADU/s
        pixel_area = self.xsciscale * u.arcsec * self.ysciscale * u.arcsec
        mjy_str = mjy / pixel_area.to(u.steradian)

        # Convert the background signal from MJy/str
        # to ADU/sec
        bval /= mjy_str
        return bval.value

    def saveSingleFits(self, image, name, key_dict=None, image2=None, image2type=None):
        #save an array into the first extension of a fits file
        h0 = fits.PrimaryHDU()
        h1 = fits.ImageHDU(image, name='DATA')
        if image2 is not None:
            h2 = fits.ImageHDU(image2)
            if image2type is not None:
                h2.header['EXTNAME'] = image2type

        # if a keyword dictionary is provided, put the
        # keywords into the 0th and 1st extension headers
        if key_dict is not None:
            for key in key_dict:
                h0.header[key] = key_dict[key]
                h1.header[key] = key_dict[key]

        if image2 is None:
            hdulist = fits.HDUList([h0, h1])
        else:
            hdulist = fits.HDUList([h0, h1, h2])
        hdulist.writeto(name, overwrite=True)

    def add_options(self, parser=None, usage=None):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage, description='Create seed image via catalogs')
        parser.add_argument("paramfile", help='File describing the input parameters and instrument settings to use. (YAML format).')
        parser.add_argument("--param_example", help='If used, an example parameter file is output.')
        return parser


if __name__ == '__main__':

    usagestring = 'USAGE: catalog_seed_image.py inputs.yaml'

    seed = Catalog_seed()
    parser = seed.add_options(usage=usagestring)
    args = parser.parse_args(namespace=seed)
    seed.make_seed()
