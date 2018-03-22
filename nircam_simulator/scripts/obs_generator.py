#! /usr/bin/env python

'''
Convert a signal rate seed image into a signal ramp.
Add cosmic rays, poisson noise, etc.
'''

import sys, os
import pkg_resources
import random
import copy
from math import radians
import datetime
import numpy as np
from astropy.io import fits, ascii
from astropy.table import Table
from astropy.time import Time, TimeDelta
import scipy.signal as s1
import yaml
from . import set_telescope_pointing_separated as stp
from . import unlinearize
from . import read_fits


INST_LIST = ['nircam', 'niriss', 'fgs']
MODES = {"nircam": ["imaging", "ts_imaging", "wfss", "ts_wfss"],
         "niriss": ["imaging"],
         "fgs": ["imaging"]}

class Observation():
    def __init__(self):
        self.linDark = None
        self.seed = None
        self.segmap = None
        self.seedheader = None
        self.seedunits = 'ADU/sec'
        
        # self.coord_adjust contains the factor by which the
        # nominal output array size needs to be increased
        # (used for WFSS mode), as well as the coordinate
        # offset between the nominal output array coordinates.
        self.coord_adjust = {'x': 1., 'xoffset': 0., 'y': 1., 'yoffset': 0.}

        # Array size of a full frame image
        self.ffsize = 2048 #pixels on a side

        # Locate the module files, so that we know where to look
        # for config subdirectory
        self.modpath = pkg_resources.resource_filename('nircam_simulator', '')

        # Get the location of the NIRCAM_SIM_DATA environment
        # variable, so we know where to look for darks, CR,
        # PSF files, etc later
        self.env_var = 'NIRCAM_SIM_DATA'
        self.datadir = os.environ.get(self.env_var)
        if self.datadir is None:
            localpath = '/ifs/jwst/wit/nircam/nircam_simulator_data'
            local = os.path.exists(localpath)
            if local:
                self.datadir = localpath
                os.environ['NIRCAM_SIM_DATA'] = localpath
            else:
                print(("WARNING: {} environment variable is not set."
                       .format(self.env_var)))
                print("This must be set to the base directory")
                print("containing the darks, cosmic ray, PSF, etc")
                print("input files needed for the simulation.")
                print("This should be set correctly if you installed")
                print("the nircam_sim_data conda package.")
                sys.exit()

    def create(self):
        """MAIN FUNCTION"""
        print('')
        print("Running observation generator....")
        print('')

        # Read in the parameter file
        self.readParameterFile()

        # Expand all paths in order to be more condor-friendly
        self.expand_env_var()
        self.filecheck()
        self.fullPaths()

        # Get the input dark if a filename is supplied
        if self.linDark is None:
            self.linDark = self.params['Reffiles']['linearized_darkfile']
            print('Reading in dark file: {}'.format(self.linDark))
        if type(self.linDark) == type('string'):
            print('Reading in dark file: {}'.format(self.linDark))
            self.linDark = self.readDarkFile(self.linDark)

        # Finally, collect information about the detector,
        # which will be needed for astrometry later
        self.detector = self.linDark.header['DETECTOR']
        self.instrument = self.linDark.header['INSTRUME']
        self.fastaxis = self.linDark.header['FASTAXIS']
        self.slowaxis = self.linDark.header['SLOWAXIS']

        # Get detector/channel specific values
        self.channel_specific_dicts()

        # Get the input seed image if a filename is supplied
        if type(self.seed) == type('string'):
            self.seed, self.segmap, self.seedheader = self.readSeed(self.seed)

        # Some basic checks on the inputs to make sure
        # the script won't have to abort due to bad inputs
        # self.checkParams()
        self.readSubarrayDefinitionFile()
        self.getSubarrayBounds()
        self.checkParams()

        # Read in cosmic ray library files if
        # CRs are to be added to the data later
        if self.runStep['cosmicray']:
            self.readCRFiles()

        # Read in gain map to be used for adding Poisson noise
        # and to scale CRs to be in ADU
        self.readGainMap()

        # If seed image is in units of electrons/sec then divide
        # by the gain to put in ADU/sec
        if 'units' in self.seedheader.keys():
            if self.seedheader['units'] in ["e-/sec","e-"]:
                print(("Seed image is in units of {}. Dividing by gain."
                       .format(self.seedheader['units'])))
                self.seed /= self.gainim
        else:
            raise ValueError(("'units' keyword not present in header of "
                             "seed image. Unable to determine whether the "
                             "seed image is in units of ADU or electrons."))
                
        # Calculate the exposure time of a single frame, based on
        # the size of the subarray
        seeddim = len(self.seed.shape)
        if seeddim == 4:
            temp_frame = self.seed[0, 0, :, :]
        elif seeddim == 3:
            temp_frame = self.seed[0, :, :]
        elif seeddim == 2:
            temp_frame = self.seed
        self.calcFrameTime(temp_frame)

        # Calculate the rate of cosmic ray hits expected per frame
        self.getCRrate()

        # Read in saturation file
        if self.params['Reffiles']['saturation'] is not None:
            self.readSaturationFile()
        else:
            print('CAUTION: no saturation map provided. Using')
            print('{} for all pixels.'.format(self.params['nonlin']['limit']))
            dy, dx = self.dark.data.shape[2:]
            self.satmap = np.zeros((dy, dx)) + self.params['nonlin']['limit']

        # Translate to ramp if necessary,
        # Add poisson noise and cosmic rays
        # Rearrange into requested read pattern
        # All done in one function to save memory
        simexp, simzero = self.add_crs_and_noise(self.seed)

        # Multiply flat fields
        simexp = self.add_flatfield_effects(simexp)
        simzero = self.add_flatfield_effects(np.expand_dims(simzero, axis=1))[:, 0, :, :]

        # Mask any reference pixels
        if self.params['Output']['grism_source_image'] == False:
            simexp, simzero = self.maskRefPix(simexp, simzero)

        # Add IPC effects
        # (Dark current ramp already has IPC in it)
        if self.runStep['ipc']:
            simexp = self.addIPC(simexp)
            simzero = self.addIPC(np.expand_dims(simzero, axis=1))[:, 0, :, :]

        # Add the simulated source ramp to the dark ramp
        lin_outramp, lin_zeroframe, lin_sbAndRefpix = self.addSyntheticToDark(simexp,
                                                                              self.linDark,
                                                                              syn_zeroframe=simzero)

        # Add other detector effects (Crosstalk/PAM)
        lin_outramp = self.add_detector_effects(lin_outramp)
        lin_zeroframe = self.add_detector_effects(np.expand_dims(lin_zeroframe, axis=1))[:, 0, :, :]

        # Read in non-linearity correction coefficients. We need these
        # regardless of whether we are saving the linearized data or going
        # on to make raw data
        nonlincoeffs = self.get_nonlinearity_coeffs()

        # We need to first subtract superbias and refpix signals from the
        # original saturation limits, and then linearize them
        # Refpix signals will vary from group to group, but only by a few
        # ADU. So let's cheat and just use the refpix signals from group 0

        # Create a linearized saturation map
        limits = np.zeros_like(self.satmap) + 1.e6

        if self.linDark.sbAndRefpix is not None:
            lin_satmap = unlinearize.nonLinFunc(self.satmap - self.linDark.sbAndRefpix[0, 0, :, :],
                                                nonlincoeffs, limits)
        elif ((self.linDark.sbAndRefpix is None) & (self.runStep['superbias'])):
            # If the superbias and reference pixel signal is not available
            # but the superbias reference file is, then just use that.
            self.readSuperbiasFile()
            lin_satmap = unlinearize.nonLinFunc(self.satmap - self.superbias,
                                                nonlincoeffs, limits)

        elif ((self.linDark.sbAndRefpix is None) & (self.runStep['superbias'] == False)):
            # If superbias and refpix signal is not available and
            # the superbias reffile is also not available, fall back to
            # a superbias value that is roughly correct. Error in this value
            # will cause errors in saturation flagging for the highest signal
            # pixels.
            manual_sb = np.zeros_like(self.satmap) + 12000.
            lin_satmap = unlinearize.nonLinFunc(self.satmap - manual_sb,
                                                nonlincoeffs, limits)

        # Save the ramp if requested. This is the linear ramp,
        # ready to go into the Jump step of the pipeline
        if 'linear' in self.params['Output']['datatype'].lower():
            # Output filename: append 'linear'
            try:
                linearrampfile = self.params['Output']['file'].replace('uncal', 'linear')
            except:
                linearrampfile = self.params['Output']['file'].replace('.fits', '_linear.fits')

            # Full path of output file
            linearrampfile = linearrampfile.split('/')[-1]
            linearrampfile = os.path.join(self.params['Output']['directory'], linearrampfile)

            # Saturation flagging - to create the pixeldq extension
            # and make data ready for ramp fitting
            # Since we subtracted the superbias and refpix signal from the
            # saturation map prior to linearizing, we can now compare that map
            # to lin_outramp, which also does not include superbias nor refpix
            # signal, and is linear.
            groupdq = self.flag_saturation(lin_outramp, lin_satmap)

            # Create the error and groupdq extensions
            err, pixeldq = self.create_other_extensions(copy.deepcopy(lin_outramp))

            if self.params['Inst']['use_JWST_pipeline']:
                self.saveDMS(lin_outramp, lin_zeroframe, linearrampfile, mod='ramp',
                             err_ext=err, group_dq=groupdq, pixel_dq=pixeldq)
            else:
                self.savefits(lin_outramp, lin_zeroframe, linearrampfile, mod='ramp',
                              err_ext=err, group_dq=groupdq, pixel_dq=pixeldq)

            stp.add_wcs(linearrampfile, roll=self.params['Telescope']['rotation'])
            print("Final linearized exposure saved to:")
            print("{}".format(linearrampfile))

        # If the raw version is requested, we need to unlinearize
        # the ramp
        if 'raw' in self.params['Output']['datatype'].lower():
            if self.linDark.sbAndRefpix is not None:
                if self.params['Output']['save_intermediates']:
                    base_name = self.params['Output']['file'].split('/')[-1]
                    ofile = os.path.join(self.params['Output']['directory'],
                                         base_name[0:-5] + '_doNonLin_accuracy.fits')
                    savefile = True
                else:
                    ofile = None
                    savefile = False

                raw_outramp = unlinearize.unlinearize(lin_outramp, nonlincoeffs, self.satmap,
                                                      lin_satmap,
                                                      maxiter=self.params['nonlin']['maxiter'],
                                                      accuracy=self.params['nonlin']['accuracy'],
                                                      save_accuracy_map=savefile,
                                                      accuracy_file=ofile)
                raw_zeroframe = unlinearize.unlinearize(lin_zeroframe, nonlincoeffs, self.satmap,
                                                        lin_satmap,
                                                        maxiter=self.params['nonlin']['maxiter'],
                                                        accuracy=self.params['nonlin']['accuracy'],
                                                        save_accuracy_map=False)

                # Add the superbias and reference pixel signal back in
                #raw_outramp = self.add_sbAndRefPix(raw_outramp,self.linDark.sbAndRefpix)
                raw_outramp = self.add_sbAndRefPix(raw_outramp, lin_sbAndRefpix)
                raw_zeroframe = self.add_sbAndRefPix(raw_zeroframe, self.linDark.zero_sbAndRefpix)

                # Make sure all signals are < 65535
                raw_outramp[raw_outramp > 65535] = 65535
                raw_zeroframe[raw_zeroframe > 65535] = 65535

                # Save the raw ramp
                base_name = self.params['Output']['file'].split('/')[-1]
                rawrampfile = os.path.join(self.params['Output']['directory'], base_name)
                if self.params['Inst']['use_JWST_pipeline']:
                    self.saveDMS(raw_outramp, raw_zeroframe, rawrampfile, mod='1b')
                else:
                    self.savefits(raw_outramp, raw_zeroframe, rawrampfile, mod='1b')
                stp.add_wcs(rawrampfile,roll=self.params['Telescope']['rotation'])
                print("Final raw exposure saved to")
                print("{}".format(rawrampfile))
            else:
                print("WARNING: raw output ramp requested, but the signal associated")
                print("with the superbias and reference pixels is not present in")
                print("the dark current data object. Quitting.")
                sys.exit()

        print("Observation generation complete.")

    def expand_env_var(self):
        """
        Replace the environment variable name in any inputs
        where it is used.
        """
        for key1 in self.params:
            for key2 in self.params[key1]:
                if self.env_var in str(self.params[key1][key2]):
                    self.params[key1][key2] = self.params[key1][key2].replace('$'+self.env_var+'/', self.datadir)

    def filecheck(self):
        """
        Make sure the requested input files exist
        For reference files, assume first that they are located in
        the directory tree under the datadir (from the NIRCAM_SIM_DATA
        environment variable). If not, assume the input is a full path
        and check there.
        """
        rlist = [['Reffiles', 'badpixmask'],
                 ['Reffiles', 'linearity'],
                 ['Reffiles', 'saturation'],
                 ['Reffiles', 'ipc'],
                 ['Reffiles', 'pixelAreaMap'],
                 ['Reffiles', 'gain']]
        plist = [['cosmicRay', 'path']]
        for ref in rlist:
            self.ref_check(ref)
        for path in plist:
            self.path_check(path)

    def ref_check(self, rele):
        """
        Check for the existence of the input reference file
        Assume first that the file is in the directory tree
        specified by the NIRCAM_SIM_DATA environment variable.

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
            rfile = os.path.abspath(rfile)
            c1 = os.path.isfile(rfile)
            if c1:
                self.params[rele[0]][rele[1]] = rfile
            else:
                raise FileNotFoundError(("WARNING: Unable to locate the "
                                         "{}, {} input file! Not present in {}."
                                         .format(rele[0], rele[1], rfile)))

    def path_check(self, p):
        """
        Check for the existence of the input path.
        Assume first that the path is in relation to
        the directory tree specified by the NIRCAM_DATA_SIM
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
        pth = os.path.abspath(pth)
        c1 = os.path.exists(pth)
        if c1:
            self.params[p[0]][p[1]] = pth
        else:
            raise NotADirectoryError(("WARNING: Unable to find the requested path "
                                      "{}. Not present in directory tree specified "
                                      "by the {} environment variable."
                                      .format(pth, self.env_var)))

    def input_check(self, inparam):
        # Check for the existence of the input file. In
        # this case we do not check the directory tree
        # specified by the NIRCAM_SIM_DATA environment variable.
        # This is intended primarily for user-generated inputs like
        # source catalogs
        ifile = self.params[inparam[0]][inparam[1]]
        if ifile.lower() != 'none':
            ifile = os.path.abspath(ifile)
            c = os.path.isfile(ifile)
            if c:
                self.params[inparam[0]][inparam[1]] = ifile
            else:
                print("WARNING: Unable to locate {}".format(ifile))
                print("Specified by the {}:{} field in".format(inparam[0], inparam[1]))
                print("the input yaml file.")
                sys.exit()

    def fullPaths(self):
        # Expand all input paths to be full paths
        # This is to allow for easier Condor-ization of
        # many runs
        pathdict = {'Reffiles':['dark', 'linearized_darkfile',
                                'superbias',
                                'subarray_defs', 'readpattdefs',
                                'linearity', 'saturation', 'gain',
                                'pixelflat', 'illumflat',
                                'astrometric', 'distortion_coeffs',
                                'ipc', 'crosstalk', 'occult',
                                'filtpupilcombo', 'pixelAreaMap',
                                'flux_cal'],
                    'cosmicRay':['path'],
                    'simSignals':['pointsource', 'psfpath',
                                  'galaxyListFile', 'extended',
                                  'movingTargetList',
                                  'movingTargetSersic',
                                  'movingTargetExtended',
                                  'movingTargetToTrack'],
                    'Output':['file', 'directory']}

        config_files = {'Reffiles-readpattdefs':'nircam_read_pattern_definitions.list'
                        , 'Reffiles-subarray_defs':'NIRCam_subarray_definitions.list'
                        , 'Reffiles-flux_cal':'NIRCam_zeropoints.list'
                        , 'Reffiles-crosstalk':'xtalk20150303g0.errorcut.txt'
                        , 'Reffiles-filtpupilcombo':'nircam_filter_pupil_pairings.list'}

        all_config_files = {'nircam': {'Reffiles-readpattdefs':'nircam_read_pattern_definitions.list'
                                       , 'Reffiles-subarray_defs':'NIRCam_subarray_definitions.list'
                                       , 'Reffiles-flux_cal':'NIRCam_zeropoints.list'
                                       , 'Reffiles-crosstalk':'xtalk20150303g0.errorcut.txt'
                                       , 'Reffiles-filtpupilcombo':'nircam_filter_pupil_pairings.list'},
                            'niriss': {'Reffiles-readpattdefs':'niriss_readout_pattern.txt'
                                       , 'Reffiles-subarray_defs':'niriss_subarrays.list'
                                       , 'Reffiles-flux_cal':'niriss_zeropoint_values.out'
                                       , 'Reffiles-crosstalk':'niriss_xtalk_zeros.txt'
                                       , 'Reffiles-filtpupilcombo':'niriss_dual_wheel_list.txt'},
                            'fgs': {'Reffiles-readpattdefs':'nircam_read_pattern_definitions.list'
                                    , 'Reffiles-subarray_defs':'NIRCam_subarray_definitions.list'
                                    , 'Reffiles-flux_cal':'NIRCam_zeropoints.list'
                                    , 'Reffiles-crosstalk':'xtalk20150303g0.errorcut.txt'
                                    , 'Reffiles-filtpupilcombo':'nircam_filter_pupil_pairings.list'}}
        config_files = all_config_files[self.params['Inst']['instrument'].lower()]
                            
        for key1 in pathdict:
            for key2 in pathdict[key1]:
                if self.params[key1][key2].lower() not in ['none', 'config']:
                    self.params[key1][key2] = os.path.abspath(self.params[key1][key2])
                elif self.params[key1][key2].lower() == 'config':
                    cfile = config_files['{}-{}'.format(key1, key2)]
                    fpath = os.path.join(self.modpath, 'config', cfile)
                    self.params[key1][key2] = fpath
                    print("'config' specified: Using {} for {}:{} input file".format(fpath, key1, key2))

    def readParameterFile(self):
        # Read in the parameter file
        try:
            with open(self.paramfile, 'r') as infile:
                self.params = yaml.load(infile)
        except:
            raise IOError("WARNING: unable to open {}".format(self.paramfile))

    def readSubarrayDefinitionFile(self):
        # Read in the file that contains a list of subarray
        # names and positions on the detector
        try:
            self.subdict = ascii.read(self.params['Reffiles']['subarray_defs'], data_start=1, header_start=0)
        except:
            raise IOError("Error: could not read in subarray definitions file.")

    def readSaturationFile(self):
        # Read in saturation map
        if self.runStep['saturation_lin_limit']:
            if 1>0:
            #try:
                self.satmap, self.satheader = self.readCalFile(self.params['Reffiles']['saturation'])
                bad = ~np.isfinite(self.satmap)
                self.satmap[bad] = 1.e6
            #except:
            else:
                print(('WARNING: unable to open saturation file {}.'
                       .format(self.params['Reffiles']['saturation'])))
                print(("Please provide a valid file, or place 'none' "
                       "in the saturation entry in the parameter file, "))
                print(("in which case the nonlin limit value in the "
                       "parameter file ({}) will be used for all pixels."
                       .format(self.params['nonlin']['limit'])))
                sys.exit()
        else:
            print(('CAUTION: no saturation map provided. Using '
                   '{} for all pixels.'.format(self.params['nonlin']['limit'])))
            dy, dx = self.dark.data.shape[2:]
            self.satmap = np.zeros((dy, dx)) + self.params['nonlin']['limit']

    def readSuperbiasFile(self):
        # Read in superbias
        if self.runStep['superbias']:
            try:
                self.superbias, self.superbiasheader = self.readCalFile(self.params['Reffiles']['superbias'])
            except:
                raise IOError(("WARNING: unable to open superbias file {}. "
                               "Please provide a valid file in the superbias "
                               "entry in the parameter file."
                               .format(self.params['Reffiles']['superbias'])))
        else:
            raise ValueError('CAUTION: no superbias provided. Quitting.')

    def readpattern_compatible(self):
        # Make sure the input dark has a readout pattern
        # that is compatible with the requested output
        # readout pattern
        rapids = ["RAPID", "NISRAPID"]
        darkpatt = self.linDark.header['READPATT']
        if ((darkpatt != self.params['Readout']['readpatt']) &
            (darkpatt not in rapids)):
            raise ValueError(("WARNING: Unable to convert input dark with a "
                              "readout pattern of {}, to the requested readout "
                              "pattern of {}. The readout pattern of the dark "
                              "must be RAPID, NISRAPID or match the requested output "
                              "readout pattern.".format(darkpatt,
                                           self.params['Readout']['readpatt'])))

    def channel_specific_dicts(self):
        # Get detector/channel-specific values for things that
        # don't need to be in the parameter file

        # Pixel scale - return as a 2-element list, with pixscale for x and y.
        filt = self.params['Readout']['filter']
        fnum = int(filt[1:4])
        if fnum < 230:
            channel = 'sw'
            self.pixscale = [0.031, 0.031]
        else:
            channel = 'lw'
            self.pixscale = [0.063, 0.063]

    def flag_saturation(self, data, sat):
        # flag saturated data. return a dq map
        # with the appropriate dq value (2)
        # for saturation
        satmap = (data > sat).astype(np.uint32)
        satmap[satmap > 0] = 2
        return satmap

    def add_sbAndRefPix(self, ramp, sbref):
        # Add superbias and reference pixel-associated
        # signal to the input ramp
        rampdim = ramp.shape
        sbrefdim = sbref.shape
        if len(rampdim) != len(sbrefdim):
            if len(rampdim) == (len(sbrefdim) + 1):
                newramp = ramp + sbref
            else:
                raise ValueError(("WARNING: input ramp and superbias+refpix "
                                  "arrays have different dimensions. Cannot combine. "
                                  "Ramp dim: {}, SBRef dim: {}"
                                  .format(len(rampdim), len(sbrefdim))))
        else:
            # Inputs arrays have the same number of dimensions
            newramp = ramp + sbref
        return newramp

    def apply_lincoeff(self, data, cof):
        # Linearize the input data
        # data will be 2d
        # cof will be 3d, or 1d
        # cof[0] + num*cof[1] + cof[2]*num**2 + cof[3]*num**3 +...
        apply = 0.
        if len(cof.shape) == 1:
            for i in range(len(cof)):
                apply += (cof[i] * data**i)
        elif len(cof.shape) == 3:
            for i in range(len(cof[:, 0, 0])):
                apply += (cof[i, :, :] * data**i)
        return apply

    def getDistortionCoefficients(self, table, from_sys, to_sys, aperture):
        # From the table of distortion coefficients,
        # get the coeffs that correspond to the requested
        # transformation and return as a list for x and another for y
        match = table['AperName'] == aperture
        if np.any(match) == False:
            raise ValueError(("Aperture name {} not found in input CSV "
                              "file.".format(aperture)))

        row = table[match]

        if ((from_sys == 'science') & (to_sys == 'ideal')):
            label = 'Sci2Idl'
        elif ((from_sys == 'ideal') & (to_sys == 'science')):
            label = 'Idl2Sci'
        else:
            raise ValueError(("WARNING: from_sys of {} and to_sys of {} not a "
                              "valid transformation.".format(from_sys, to_sys)))

        # Get the coefficients, return as list
        X_cols = [c for c in row.colnames if label+'X' in c]
        Y_cols = [c for c in row.colnames if label+'Y' in c]
        x_coeffs = [row[c].data[0] for c in X_cols]
        y_coeffs = [row[c].data[0] for c in Y_cols]

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

    def get_nonlinearity_coeffs(self):
        if self.params['Reffiles']['linearity'] is not None:
            try:
                nonlin = self.get_nonlin_coeffs(self.params['Reffiles']['linearity'])
            except:
                print("Unable to read in non-linearity correction coefficients")
                print("from {}.".format(self.params['Reffiles']['linearity']))
                print("Using a set of mean coefficients.")
                nonlin = np.array([0., 1.0, 9.69903112e-07, 3.85263835e-11,
                                   1.09267058e-16, -5.30613939e-20, 9.27963411e-25])
        else:
            print("No linearity coefficient file provided. Proceeding using a")
            print("set of mean coefficients derived from CV3 data.")
            nonlin = np.array([0., 1.0, 9.69903112e-07, 3.85263835e-11,
                               1.09267058e-16, -5.30613939e-20, 9.27963411e-25])
        # print('Nonlinearity coefficients: ', nonlin)
        return nonlin


    def create_other_extensions(self, data):
        # If the linearized version of the file is to be saved, we
        # need to create the error, pixel dq, and group dq extensions

        # error extension - keep it simple. sqrt of signal
        toolow = data < 0.
        data[toolow] = 0.
        err = np.sqrt(data)

        # pixel dq extension - populate using the mask reference file
        if self.runStep['badpixfile']:
            mask_hdu = fits.open(self.params['Reffiles']['badpixmask'])
            mask = mask_hdu[1].data
            dqdef = mask_hdu[2].data
            mask_hdu.close()

            # Crop to match output subarray size
            if "FULL" not in self.params['Readout']['array_name']:
                mask = self.crop_to_subarray(mask)

            # If the JWST pipeline is available,
            # convert mask values to those used by the pipeline
            # based on the names in dq_def. This function is basically
            # a copy of the dynamic_mask function in dynamicdq.py
            # in the JWST pipeline
            if self.params['Inst']['use_JWST_pipeline']:
                pixeldq = self.convert_mask(mask, dqdef)
            else:
                # If the pipeline is not to be used, then the
                # best we can do is assume that the input bad
                # pixel value definitions match what the pipeline
                # expects, and keep the mask as read in.
                pixeldq = mask
        else:
            print("No bad pixel mask provided. Setting all pixels in")
            print("pixel data quality extension to 0, indicating they")
            print("are good.")
            pixeldq = np.zeros(data.shape[2:]).astype(np.uint32)

        # group dq extension - populate by flagging saturated pix
        #do this earlier (at the unlinearize step)
        return err, pixeldq


    def create_group_entry(self, integration, groupnum, endday, endmilli, endsubmilli, endgroup, xd, yd, gap, comp_code, comp_text, barycentric, heliocentric):
        #add the GROUP extension to the output file

        #from an example Mark Kyprianou sent:
        #integration_number, group_number - integers
        #end_day - integer (e.g. 5861) days since Jan 1 2000
        #end_milliseconds - integer (e.g. 9806061) milliseconds of the day
        #end_submilliseconds - integer (e.g 177) since last millisecond?
        #group_end_time e.g. '2016-01-18T02:43:26.061'
        #number_of_columns 2048
        #number_of_rows 2048
        #number_of_gaps 0 gaps in telemetry
        #completion_code_number 0 (nominal?)
        #Completion_code_text 'COMPLETE'-from howard
        #                     'Normal Completion' - from mark
        #bary_end_time (mjd) 57405.11165225
        #helio_end_time (mjd) 57405.1163058

        group = np.ndarray(
            (1, ),
            dtype=[
                ('integration_number', '<i2'),
                ('group_number', '<i2'),
                ('end_day', '<i2'),
                ('end_milliseconds', '<i4'),
                ('end_submilliseconds', '<i2'),
                ('group_end_time', 'S26'),
                ('number_of_columns', '<i2'),
                ('number_of_rows', '<i2'),
                ('number_of_gaps', '<i2'),
                ('completion_code_number', '<i2'),
                ('completion_code_text', 'S36'),
                ('bary_end_time', '<f8'),
                ('helio_end_time', '<f8')
            ]
        )
        group[0]['integration_number'] = integration
        group[0]['group_number'] = groupnum
        group[0]['end_day'] = endday
        group[0]['end_milliseconds'] = endmilli
        group[0]['end_submilliseconds'] = endsubmilli
        group[0]['group_end_time'] = endgroup
        group[0]['number_of_columns'] = xd
        group[0]['number_of_rows'] = yd
        group[0]['number_of_gaps'] = gap
        group[0]['completion_code_number'] = comp_code
        group[0]['completion_code_text'] = comp_text
        group[0]['bary_end_time'] = barycentric
        group[0]['helio_end_time'] = heliocentric
        return group


    def populate_group_table(self, starttime, grouptime, ramptime, numint, numgroup, ny, nx):
        #create some reasonable values to fill the GROUP extension table.
        #These will not be completely correct because access to other ssb
        #scripts and more importantly, databses, is necessary. But they should be
        #close.

        #print("Going in to populate_group_table:")
        #print(starttime, grouptime, ramptime, numint, numgroup, ny, nx)

        #Create the table with a first row populated by garbage
        grouptable = self.create_group_entry(999, 999, 0, 0, 0, 'void', 0, 0, 0, 0, 'void', 1., 1.)

        #Fixed for all exposures
        compcode = 0
        comptext = 'Normal Completion'
        numgap = 0
        baseday = Time('2000-01-01T00:00:00')

        #integration start times
        rampdelta = TimeDelta(ramptime, format='sec')
        groupdelta = TimeDelta(grouptime, format='sec')
        intstarts = starttime + (np.arange(numint)*rampdelta)

        for integ in range(numint):
            groups = np.arange(1, numgroup+1)
            groupends = intstarts[integ] + (np.arange(1, numgroup+1)*groupdelta)
            endday = (groupends - baseday).jd
            enddayint = [np.int(s) for s in endday]

            #now to get end_milliseconds, we need milliseconds from the beginning
            #of the day
            inday = TimeDelta(endday - enddayint, format='jd')
            endmilli = inday.sec * 1000.

            #submilliseconds - just use a random number
            endsubmilli = np.random.randint(0, 1000, len(endmilli))

            #group end time. need to remove : and - and make lowercase t
            #gstr = groupends.isot
            groupending = groupends.isot
            #groupending = [s.replace(':', '').replace('-', '').lower() for s in gstr]

            #approximate these as just the group end time in mjd
            barycentric = groupends.mjd
            heliocentric = groupends.mjd

            for grp, day, milli, submilli, grpstr, bary, helio in zip(groups, endday, endmilli, endsubmilli, groupending, barycentric, heliocentric):
                entry = self.create_group_entry(integ+1, grp, day, milli, submilli, grpstr, nx, ny, numgap, compcode, comptext, bary, helio)
                grouptable = np.vstack([grouptable, entry])

        # Now remove the top garbage row from the table
        grouptable = grouptable[1:]
        return grouptable


    def convert_mask(self, inmask, dq_table):
        #
        # Return a mask with values converted those those used
        # by the JWST pipeline
        from jwst.datamodels import dqflags

        # Get the DQ array and the flag definitions
        if (dq_table is not None and
            not np.isscalar(dq_table) and
            len(dq_table.shape) and
            len(dq_table)):
            #
            # Make an empty mask
            dqmask = np.zeros(inmask.shape, dtype=np.uint32)
            for record in dq_table:
                bitplane = record['VALUE']
                dqname = record['NAME'].strip()
                try:
                    standard_bitvalue = dqflags.pixel[dqname]
                except KeyError:
                    print(('Keyword {} does not correspond to an existing DQ '
                           'mnemonic, so will be ignored'.format(dqname)))
                    continue
                just_this_bit = np.bitwise_and(inmask, bitplane)
                pixels = np.where(just_this_bit != 0)
                dqmask[pixels] = np.bitwise_or(dqmask[pixels], standard_bitvalue)
        else:
            dqmask = inmask

        return dqmask


    def saveDMS(self, ramp, zeroframe, filename, mod='1b', err_ext = None
                , group_dq = None, pixel_dq = None):
        # Save the new, simulated integration in DMS format (i.e. DMS orientation
        # rather than raw fitswriter orientation, and using SSB's datamodels)
        if mod == '1b':
            from jwst.datamodels import Level1bModel as DataModel
        elif mod == 'ramp':
            from jwst.datamodels import RampModel as DataModel
        else:
            raise ValueError(("Model type to use for saving output is "
                              "not recognized. Must be either '1b' or 'ramp'."))
        outModel = DataModel()

        #make sure the ramp to be saved has the right number of dimensions
        imshape = ramp.shape
        if len(imshape) == 3:
            ramp = np.expand_dims(ramp, axis=0)

        #insert data into model
        outModel.data = ramp

        if mod == 'ramp':
            outModel.err = err_ext
            outModel.groupdq = group_dq
            outModel.pixeldq = pixel_dq

        #if saving the zeroth frame is requested, insert into the model instance
        if zeroframe is not None:
            #if the zeroframe is a 2D image, then add a dimension,
            #as the model expects 3D
            if len(zeroframe.shape) == 2:
                zeroframe = np.expand_dims(zeroframe, 0)
            outModel.zeroframe = zeroframe#.astype(np.uint16)
        else:
            print("Zeroframe not present. Setting to all zeros")
            numint, numgroup, ys, xs = ramp.shape
            outModel.zeroframe = np.zeros((numint, ys, xs))

        #EXPTYPE OPTIONS
        #exptypes = ['NRC_IMAGE', 'NRC_GRISM', 'NRC_TACQ', 'NRC_CORON',
        #            'NRC_DARK', 'NRC_TSIMAGE', 'NRC_TSGRISM']
        #nrc_tacq and nrc_coron are not currently implemented.

        exptype = {"nircam": {"imaging": "NRC_IMAGE", "ts_imaging": "NRC_TSIMAGE",
                              "wfss": "NRC_GRISM", "ts_wfss": "NRC_TSGRISM"},
                   "niriss": {"imaging": "NIS_IMAGE"},
                   "fgs": {"imaging": "FGS_IMAGE"}}

        try:
            outModel.meta.exposure.type = exptype[self.params['Inst']['instrument'].lower()]\
                                          [self.params['Inst']['mode'].lower()]
        except:
            raise ValueError('EXPTYPE mapping not complete for this!!! FIX ME!')

        #update various header keywords
        dims = outModel.data.shape
        dtor = radians(1.)
        pixelsize = self.pixscale[0] / 3600.0

        current_time = datetime.datetime.utcnow()
        start_time_string = self.params['Output']['date_obs'] + 'T' + self.params['Output']['time_obs']
        ct = Time(start_time_string)

        outModel.meta.date = start_time_string
        outModel.meta.telescope = 'JWST'
        outModel.meta.instrument.name = self.params['Inst']['instrument'].upper()
        if self.instrument.upper() == 'NIRCAM':
            outModel.meta.instrument.module = self.detector[3]
            channel = 'SHORT'
            if 'LONG' in self.detector:
                channel = 'LONG'
            outModel.meta.instrument.channel = channel
                
        outModel.meta.instrument.detector = self.detector
        outModel.meta.coordinates.reference_frame = 'ICRS'

        outModel.meta.subarray.fastaxis = self.fastaxis
        outModel.meta.subarray.slowaxis = self.slowaxis

        outModel.meta.origin = 'STScI'
        outModel.meta.filename = filename
        outModel.meta.filetype = 'raw'
        outModel.meta.observation.obs_id = self.params['Output']['obs_id']
        outModel.meta.observation.visit_id = self.params['Output']['visit_id']
        outModel.meta.observation.visit_number = self.params['Output']['visit_number']
        outModel.meta.observation.program_number = self.params['Output']['program_number']
        outModel.meta.observation.observation_number = self.params['Output']['observation_number']
        outModel.meta.observation.observation_label = self.params['Output']['observation_label']
        outModel.meta.observation.visit_group = self.params['Output']['visit_group']
        outModel.meta.observation.sequence_id = self.params['Output']['sequence_id']
        outModel.meta.observation.activity_id =self.params['Output']['activity_id']
        outModel.meta.observation.exposure_number = self.params['Output']['exposure_number']

        outModel.meta.program.pi_name = self.params['Output']['PI_Name']
        outModel.meta.program.title = self.params['Output']['title']
        outModel.meta.program.category = self.params['Output']['Proposal_category']
        outModel.meta.program.sub_category = 'UNKNOWN'
        outModel.meta.program.science_category = self.params['Output']['Science_category']
        outModel.meta.program.continuation_id = 0

        outModel.meta.target.catalog_name = 'UNKNOWN'
        outModel.meta.coordinates.reference_frame = 'ICRS'
        
        outModel.meta.wcsinfo.wcsaxes = 2
        outModel.meta.wcsinfo.crval1 = self.ra
        outModel.meta.wcsinfo.crval2 = self.dec
        outModel.meta.wcsinfo.crpix1 = self.refpix_pos['x']+1.
        outModel.meta.wcsinfo.crpix2 = self.refpix_pos['y']+1.
        outModel.meta.wcsinfo.ctype1 = 'RA---TAN'
        outModel.meta.wcsinfo.ctype2 = 'DEC--TAN'
        outModel.meta.wcsinfo.cunit1 = 'deg'
        outModel.meta.wcsinfo.cunit2 = 'deg'
        outModel.meta.wcsinfo.v2_ref = self.v2_ref
        outModel.meta.wcsinfo.v3_ref = self.v3_ref
        outModel.meta.wcsinfo.vparity = self.parity
        outModel.meta.wcsinfo.v3yangle = self.v3yang
        outModel.meta.wcsinfo.cdelt1 = self.xsciscale / 3600.
        outModel.meta.wcsinfo.cdelt2 = self.ysciscale / 3600.
        outModel.meta.wcsinfo.roll_ref = self.local_roll
        outModel.meta.target.ra = self.ra
        outModel.meta.target.dec = self.dec

        #ra_v1, dec_v1, and pa_v3 are not used by the level 2 pipelines
        outModel.meta.pointing.ra_v1 = self.ra
        outModel.meta.pointing.dec_v1 = self.dec
        outModel.meta.pointing.pa_v3 = self.params['Telescope']['rotation']

        ramptime = self.frametime*(1+self.params['Readout']['ngroup']*(self.params['Readout']['nframe']+self.params['Readout']['nskip']))
        #Add time for the reset frame....
        rampexptime = self.frametime * (self.params['Readout']['ngroup']*(self.params['Readout']['nframe']+self.params['Readout']['nskip']))

        outModel.meta.observation.date = self.params['Output']['date_obs']
        outModel.meta.observation.time = self.params['Output']['time_obs']

        if self.runStep['fwpw']:
            fwpw = ascii.read(self.params['Reffiles']['filtpupilcombo'])
        else:
            print("WARNING: Filter wheel element/pupil wheel element combo reffile not specified. Proceeding by")
            print("saving {} in FILTER keyword, and {} in PUPIL keyword".format(self.params['Readout']['filter'], self.params['Readout']['pupil']))
            fwpw = Table()
            fwpw['filter_wheel'] = self.params['Readout']['filter']
            fwpw['pupil_wheel'] = self.params['Readout']['pupil']

        #get the proper filter wheel and pupil wheel values for the header
        if self.params['Inst']['mode'].lower() not in ['wfss','ts_wfss']:
            mtch = fwpw['filter'] == self.params['Readout']['filter'].upper()
            fw = str(fwpw['filter_wheel'].data[mtch][0])
            pw = str(fwpw['pupil_wheel'].data[mtch][0])
            #grism='N/A'
        else:
            pw = self.params['Readout']['pupil']
            fw = self.params['Readout']['filter']
            #grism = pw

        outModel.meta.instrument.filter = fw
        outModel.meta.instrument.pupil = pw

        outModel.meta.dither.primary_type = self.params['Output']['primary_dither_type']
        outModel.meta.dither.position_number = self.params['Output']['primary_dither_position']
        outModel.meta.dither.total_points = self.params['Output']['total_primary_dither_positions']
        outModel.meta.dither.pattern_size = 0.0
        outModel.meta.dither.subpixel_type = self.params['Output']['subpix_dither_type']
        outModel.meta.dither.subpixel_number = self.params['Output']['subpix_dither_position']
        outModel.meta.dither.subpixel_total_points = self.params['Output']['total_subpix_dither_positions']
        outModel.meta.dither.xoffset = self.params['Output']['xoffset']
        outModel.meta.dither.yoffset = self.params['Output']['yoffset']

        # pixel coordinates in FITS header start from 1 not from 0
        xc = (self.subarray_bounds[2] + self.subarray_bounds[0])/2.+1.
        yc = (self.subarray_bounds[3] + self.subarray_bounds[1])/2.+1.

        outModel.meta.exposure.readpatt = self.params['Readout']['readpatt']

        #The subarray name needs to come from the "Name" column in the
        #subarray definitions dictionary
        mtch = self.subdict["AperName"] == self.params["Readout"]['array_name']
        outModel.meta.subarray.name = str(self.subdict["Name"].data[mtch][0])

        #subarray_bounds indexed to zero, but values in header should be
        #indexed to 1.
        outModel.meta.subarray.xstart = self.subarray_bounds[0]+1
        outModel.meta.subarray.ystart = self.subarray_bounds[1]+1
        outModel.meta.subarray.xsize = self.subarray_bounds[2]-self.subarray_bounds[0]+1
        outModel.meta.subarray.ysize = self.subarray_bounds[3]-self.subarray_bounds[1]+1

        nlrefpix = max(4-self.subarray_bounds[0], 0)
        nbrefpix = max(4-self.subarray_bounds[1], 0)
        nrrefpix=max(self.subarray_bounds[2]-(self.ffsize-4), 0)
        ntrefpix=max(self.subarray_bounds[3]-(self.ffsize-4), 0)

        outModel.meta.exposure.nframes = self.params['Readout']['nframe']
        outModel.meta.exposure.ngroups = self.params['Readout']['ngroup']
        outModel.meta.exposure.nints = self.params['Readout']['nint']

        outModel.meta.exposure.sample_time = 10
        outModel.meta.exposure.frame_time = self.frametime
        outModel.meta.exposure.group_time = self.frametime*(self.params['Readout']['nframe']+self.params['Readout']['nskip'])
        outModel.meta.exposure.groupgap = self.params['Readout']['nskip']

        outModel.meta.exposure.nresets_at_start = 1
        outModel.meta.exposure.nresets_between_ints = 1
        outModel.meta.exposure.integration_time = rampexptime
        outModel.meta.exposure.exposure_time = rampexptime * self.params['Readout']['nint']
        outModel.meta.model_type = 'RampModel'

        #set the exposure start time
        outModel.meta.exposure.start_time = ct.mjd
        endingTime = ct.mjd + outModel.meta.exposure.exposure_time/3600./24.
        outModel.meta.exposure.end_time = endingTime
        outModel.meta.exposure.mid_time = ct.mjd + outModel.meta.exposure.exposure_time/3600./24./2.
        outModel.meta.exposure.duration = ramptime

        #populate the GROUP extension table
        n_int, n_group, n_y, n_x = outModel.data.shape
        outModel.group = self.populate_group_table(ct, outModel.meta.exposure.group_time, rampexptime, n_int, n_group, n_y, n_x)

        outModel.save(filename)

        # Now we need to adjust the datamodl header keyword
        # If we leave it as Level1bModel, the pipeline doesn't
        # work properly
        if mod == '1b':
            temp = fits.open(filename, mode='update')
            temp[0].header['DATAMODL'] = 'RampModel'
            temp.flush()

        #print("Final output integration saved to {}".format(filename))
        return


    def savefits(self, ramp, zeroframe, filename, mod='1b', err_ext = None
                , group_dq = None, pixel_dq = None):
        #save the new, simulated integration in DMS format (i.e. DMS orientation
        #rather than raw fitswriter orientation, and using SSBs RampModel)

        #make sure the ramp to be saved has the right number of dimensions
        imshape = ramp.shape
        if len(imshape) == 3:
            ramp = np.expand_dims(ramp, axis=0)

        if mod == '1b':
            toohigh = ramp > 65535
            ramp[toohigh] = 65535

        #if saving the zeroth frame is requested, insert into the model instance
        if zeroframe is not None:
            #if the zeroframe is a 2D image, then add a dimension
            if len(zeroframe.shape) == 2:
                zeroframe = np.expand_dims(zeroframe, 0)
        else:
            print("Zeroframe not present. Setting to all zeros")
            numint, numgroup, ys, xs = ramp.shape

        # Place the arrays in the correct extensions of the HDUList
        #using int16 below causes problems! anything set to 65535
        #gets reset to -1, which screws up saturation flagging
        #I think the answer is to save as uint16...

        if mod == 'ramp':
            ex0 = fits.PrimaryHDU()
            ex1 = fits.ImageHDU(ramp.astype(np.float32), name='SCI')
            ex2 = fits.ImageHDU(pixel_dq.astype(np.uint32), name='PIXELDQ')
            ex3 = fits.ImageHDU(group_dq.astype(np.uint8), name='GROUPDQ')
            ex4 = fits.ImageHDU(err_ext.astype(np.float32), name='ERR')
            ex5 = fits.ImageHDU(zeroframe.astype(np.float32), name='ZEROFRAME')
            ex6 = fits.BinTableHDU(name='GROUP')
            outModel = fits.HDUList([ex0, ex1, ex2, ex3, ex4, ex5, ex6])
            groupextnum = 6

        elif mod == '1b':
            ex0 = fits.PrimaryHDU()
            ex1 = fits.ImageHDU(ramp.astype(np.uint16), name='SCI')
            ex2 = fits.ImageHDU(zeroframe.astype(np.uint16), name='ZEROFRAME')
            ex3 = fits.BinTableHDU(name='GROUP')
            outModel = fits.HDUList([ex0, ex1, ex2, ex3])
            groupextnum = 3

        #self.dark[1].data = ramp.astype(np.uint16)
        #self.dark[2].data = zeroframe.astype(np.uint16)

        #if there are other extensions, remove them
        #if len(self.dark) > 3:
        #    self.dark = self.dark[0:3]

        #EXPTYPE OPTIONS
        #exptype = {'nircam': {'imaging': 'NRC_IMAGE', 'ts_imaging': 'NRC_TSIMAGE',
        #           'wfss': 'NRC_GRISM', 'ts_wfss': 'NRC_TSGRISM'}}
        exptype = {"nircam": {"imaging": "NRC_IMAGE", "ts_imaging": "NRC_TSIMAGE",
                              "wfss": "NRC_GRISM", "ts_wfss": "NRC_TSGRISM"},
                   "niriss": {"imaging": "NIS_IMAGE"},
                   "fgs": {"imaging": "FGS_IMAGE"}}

        try:
            outModel[0].header['EXP_TYPE'] = exptype[self.params['Inst']['instrument'].lower()]\
                                             [self.params['Inst']['mode'].lower()]
        except:
            raise ValueError('EXPTYPE mapping not complete for this!!! FIX ME!')

        #update various header keywords
        dims = outModel[1].data.shape
        dtor = radians(1.)
        pixelsize = self.pixscale[0] / 3600.0

        current_time = datetime.datetime.utcnow()
        start_time_string = self.params['Output']['date_obs'] + 'T' + self.params['Output']['time_obs']
        ct = Time(start_time_string)

        outModel[0].header['DATE'] = start_time_string
        outModel[0].header['TELESCOP'] = 'JWST'

        outModel[0].header['INSTRUME'] = self.params['Inst']['instrument'].upper()
        outModel[0].header['DETECTOR'] = self.detector

        if self.instrument.upper() == 'NIRCAM':
            outModel[0].header['MODULE'] = self.detector[3]
            channel = 'SHORT'
            if 'LONG' in self.detector:
                channel = 'LONG'
            outModel[0].header['CHANNEL'] = channel
            
        outModel[0].header['FASTAXIS'] = self.fastaxis
        outModel[0].header['SLOWAXIS'] = self.slowaxis

        outModel[1].header['RADESYS'] = 'ICRS'

        outModel[0].header['ORIGIN'] = 'STScI'
        outModel[0].header['FILENAME'] = os.path.split(filename)[1]
        outModel[0].header['FILETYPE'] = 'raw'
        outModel[0].header['OBS_ID'] = self.params['Output']['obs_id']
        outModel[0].header['VISIT_ID'] = self.params['Output']['visit_id']
        outModel[0].header['VISIT'] = self.params['Output']['visit_number']
        outModel[0].header['PROGRAM'] = self.params['Output']['program_number']
        outModel[0].header['OBSERVTN'] = self.params['Output']['observation_number']
        outModel[0].header['OBSLABEL'] = self.params['Output']['observation_label']
        outModel[0].header['VISITGRP'] = self.params['Output']['visit_group']
        outModel[0].header['SEQ_ID'] = self.params['Output']['sequence_id']
        outModel[0].header['ACT_ID'] = self.params['Output']['activity_id']
        outModel[0].header['EXPOSURE'] = self.params['Output']['exposure_number']

        outModel[0].header['PI_NAME'] = self.params['Output']['PI_Name']
        outModel[0].header['TITLE'] = self.params['Output']['title']
        outModel[0].header['CATEGORY'] = self.params['Output']['Proposal_category']
        outModel[0].header['SUBCAT'] = 'UNKNOWN'
        outModel[0].header['SCICAT'] = self.params['Output']['Science_category']
        outModel[0].header['CONT_ID'] = 0

        outModel[0].header['TARGNAME'] = 'UNKNOWN'

        outModel[1].header['WCSAXES'] = 2
        outModel[1].header['CRVAL1'] = self.ra
        outModel[1].header['CRVAL2'] = self.dec
        outModel[1].header['CRPIX1'] = self.refpix_pos['x']+1.
        outModel[1].header['CRPIX2'] = self.refpix_pos['y']+1.
        outModel[1].header['CTYPE1'] = 'RA---TAN'
        outModel[1].header['CTYPE2'] = 'DEC--TAN'
        outModel[1].header['CUNIT1'] = 'deg'
        outModel[1].header['CUNIT2'] = 'deg'
        outModel[1].header['V2_REF'] = self.v2_ref
        outModel[1].header['V3_REF'] = self.v3_ref
        outModel[1].header['VPARITY'] = self.parity
        outModel[1].header['V3I_YANG'] = self.v3yang
        outModel[1].header['CDELT1'] = self.xsciscale / 3600.
        outModel[1].header['CDELT2'] = self.ysciscale / 3600.
        outModel[1].header['ROLL_REF'] = self.local_roll

        outModel[0].header['TARG_RA'] = self.ra #not correct
        outModel[0].header['TARG_DEC'] = self.dec #not correct

        #ra_v1, dec_v1, and pa_v3 are not used by the level 2 pipelines
        outModel[0].header['RA_V1'] = self.ra
        outModel[0].header['DEC_V1'] = self.dec
        outModel[0].header['PA_V3'] = self.params['Telescope']['rotation']

        ramptime = self.frametime*(1+self.params['Readout']['ngroup']*(self.params['Readout']['nframe']+self.params['Readout']['nskip']))
        #Add time for the reset frame....
        rampexptime = self.frametime * (self.params['Readout']['ngroup']*(self.params['Readout']['nframe']+self.params['Readout']['nskip']))

        # elapsed time from the end and from the start of the supposid ramp, in seconds
        # put the end of the ramp 1 second before the time the file is written
        # these only go in the fake ramp, not in the signal images....
        outModel[0].header['DATE-OBS'] = self.params['Output']['date_obs']
        outModel[0].header['TIME-OBS'] = self.params['Output']['time_obs']

        #Hmmm. Do we need a file with filter/pupil combos? Or should we rely
        #on the user to enter the correct values in the parameter file?
        if self.runStep['fwpw']:
            fwpw = ascii.read(self.params['Reffiles']['filtpupilcombo'])
        else:
            print("WARNING: Filter wheel element/pupil wheel element combo reffile not specified. Proceeding by")
            print("saving {} in FILTER keyword, and {} in PUPIL keyword".format(self.params['Readout']['filter'], self.params['Readout']['pupil']))
            fwpw = Table()
            fwpw['filter_wheel'] = self.params['Readout']['filter']
            fwpw['pupil_wheel'] = self.params['Readout']['pupil']

        #get the proper filter wheel and pupil wheel values for the header
        if self.params['Inst']['mode'].lower() not in ['wfss','ts_wfss']:
            mtch = fwpw['filter'] == self.params['Readout']['filter'].upper()
            fw = str(fwpw['filter_wheel'].data[mtch][0])
            pw = str(fwpw['pupil_wheel'].data[mtch][0])
        else:
            pw = self.params['Readout']['pupil']
            fw = self.params['Readout']['filter']

        outModel[0].header['FILTER'] = fw
        outModel[0].header['PUPIL'] = pw

        outModel[0].header['PATTTYPE'] = self.params['Output']['primary_dither_type']
        outModel[0].header['PATT_NUM'] = self.params['Output']['primary_dither_position']
        outModel[0].header['NUMDTHPT'] = self.params['Output']['total_primary_dither_positions']
        outModel[0].header['PATTSIZE'] = 0.0
        outModel[0].header['SUBPXTYP'] = self.params['Output']['subpix_dither_type']
        outModel[0].header['SUBPXNUM'] = self.params['Output']['subpix_dither_position']
        outModel[0].header['SUBPXPNS'] = self.params['Output']['total_subpix_dither_positions']
        outModel[0].header['XOFFSET'] = self.params['Output']['xoffset']
        outModel[0].header['YOFFSET'] = self.params['Output']['yoffset']

        # pixel coordinates in FITS header start from 1 not from 0
        xc=(self.subarray_bounds[2]+self.subarray_bounds[0])/2.+1.
        yc=(self.subarray_bounds[3]+self.subarray_bounds[1])/2.+1.

        outModel[0].header['READPATT'] = self.params['Readout']['readpatt']

        #The subarray name needs to come from the "Name" column in the
        #subarray definitions dictionary
        mtch = self.subdict["AperName"] == self.params["Readout"]['array_name']
        outModel[0].header['SUBARRAY'] = str(self.subdict["Name"].data[mtch][0])

        #subarray_bounds indexed to zero, but values in header should be
        #indexed to 1.
        outModel[0].header['SUBSTRT1'] = self.subarray_bounds[0]+1
        outModel[0].header['SUBSTRT2'] = self.subarray_bounds[1]+1
        outModel[0].header['SUBSIZE1'] = self.subarray_bounds[2]-self.subarray_bounds[0]+1
        outModel[0].header['SUBSIZE2'] = self.subarray_bounds[3]-self.subarray_bounds[1]+1

        nlrefpix = max(4-self.subarray_bounds[0], 0)
        nbrefpix = max(4-self.subarray_bounds[1], 0)
        nrrefpix=max(self.subarray_bounds[2]-(self.ffsize-4), 0)
        ntrefpix=max(self.subarray_bounds[3]-(self.ffsize-4), 0)

        outModel[0].header['NFRAMES'] = self.params['Readout']['nframe']
        outModel[0].header['NGROUPS'] = self.params['Readout']['ngroup']
        outModel[0].header['NINTS'] = self.params['Readout']['nint']

        outModel[0].header['TSAMPLE'] = 10
        outModel[0].header['TFRAME'] = self.frametime
        outModel[0].header['TGROUP'] = self.frametime*(self.params['Readout']['nframe']+self.params['Readout']['nskip'])
        outModel[0].header['GROUPGAP'] = self.params['Readout']['nskip']

        outModel[0].header['NRSTSTRT'] = 1
        outModel[0].header['NRESETS'] = 1
        outModel[0].header['EFFINTTM'] = rampexptime
        outModel[0].header['EFFEXPTM'] = rampexptime * self.params['Readout']['nint']

        #set the exposure start time as the current time
        outModel[0].header['EXPSTART'] = ct.mjd
        outModel[0].header['EXPEND'] = ct.mjd + outModel[0].header['EFFEXPTM']/3600./24.
        outModel[0].header['EXPMID'] = ct.mjd + outModel[0].header['EFFEXPTM']/3600./24./2.

        outModel[0].header['DURATION'] = ramptime

        #populate the GROUP extension table
        n_int, n_group, n_y, n_x = outModel[1].data.shape
        outModel[groupextnum].data = self.populate_group_table(ct, outModel[0].header['TGROUP'], rampexptime, n_int, n_group, n_y, n_x)


        print('wcsaxes is {}'.format(outModel[1].header['WCSAXES']))
        outModel.writeto(filename, overwrite=True)

        #print("Final output integration saved to {}".format(filename))
        return filename


    def populate_group_table(self, starttime, grouptime, ramptime, numint, numgroup, ny, nx):
        #create some reasonable values to fill the GROUP extension table.
        #These will not be completely correct because access to other ssb
        #scripts and more importantly, databses, is necessary. But they should be
        #close.

        #print("Going in to populate_group_table:")
        #print(starttime, grouptime, ramptime, numint, numgroup, ny, nx)

        #Create the table with a first row populated by garbage
        grouptable = self.create_group_entry(999, 999, 0, 0, 0, 'void', 0, 0, 0, 0, 'void', 1., 1.)

        #Fixed for all exposures
        compcode = 0
        comptext = 'Normal Completion'
        numgap = 0
        baseday = Time('2000-01-01T00:00:00')

        #integration start times
        rampdelta = TimeDelta(ramptime, format='sec')
        groupdelta = TimeDelta(grouptime, format='sec')
        intstarts = starttime + (np.arange(numint)*rampdelta)

        for integ in range(numint):
            groups = np.arange(1, numgroup+1)
            groupends = intstarts[integ] + (np.arange(1, numgroup+1)*groupdelta)
            endday = (groupends - baseday).jd
            enddayint = [np.int(s) for s in endday]

            #now to get end_milliseconds, we need milliseconds from the beginning
            #of the day
            inday = TimeDelta(endday - enddayint, format='jd')
            endmilli = inday.sec * 1000.

            #submilliseconds - just use a random number
            endsubmilli = np.random.randint(0, 1000, len(endmilli))

            #group end time. need to remove : and - and make lowercase t
            #gstr = groupends.isot
            groupending = groupends.isot
            #groupending = [s.replace(':', '').replace('-', '').lower() for s in gstr]

            #approximate these as just the group end time in mjd
            barycentric = groupends.mjd
            heliocentric = groupends.mjd

            for grp, day, milli, submilli, grpstr, bary, helio in zip(groups, endday, endmilli, endsubmilli, groupending, barycentric, heliocentric):
                entry = self.create_group_entry(integ+1, grp, day, milli, submilli, grpstr, nx, ny, numgap, compcode, comptext, bary, helio)
                grouptable = np.vstack([grouptable, entry])

        # Now remove the top garbage row from the table
        grouptable = grouptable[1:]
        return grouptable

    def crop_to_subarray(self, data):
        return data[self.subarray_bounds[1]:self.subarray_bounds[3] + 1,
                    self.subarray_bounds[0]:self.subarray_bounds[2] + 1]

    def get_nonlin_coeffs(self, linfile):
        # Read in non-linearity coefficients
        nonlin, nonlinheader = self.readCalFile(linfile)
        # Set NaN coefficients such that no correction will be made
        nans = np.isnan(nonlin[1, :, :])
        numnan = np.sum(nans)
        print('The linearity coefficients of {} pixels are NaNs. Setting these'.format(numnan))
        print('coefficients such that no linearity correction is made.')
        for i, cof in enumerate(range(nonlin.shape[0])):
            tmp = nonlin[cof, :, :]
            if i == 1:
                tmp[nans] = 1.
            else:
                tmp[nans] = 0.
            nonlin[cof, :, :] = tmp

        # # Crop to appropriate subarray - ALREADY DONE IN readCalFile
        # if "FULL" not in self.params['Readout']['array_name']:
        #     nonlin = self.crop_to_subarray(nonlin)
        return nonlin

    def readDarkFile(self, filename):
        # Read in a prepared dark current exposure
        obj = read_fits.Read_fits()
        obj.file = filename
        obj.read_astropy()
        # obj now contains: obj.data, obj.sbAndRefpix,
        #                   obj.zeroframe, obj.zero_sbAndRefpix,
        #                   obj.header
        # Values are None for objects that don't exist
        return obj

    def readSeed(self, filename):
        # Read in the file containing the seed image/ramp
        with fits.open(filename) as h:
            seed = h[1].data
            seedheader = h[0].header
            try:
                segmap = h[2].data
            except:
                segmap = None
        return seed, segmap, seedheader

    def calcFrameTime(self, frame):
        #calculate the exposure time of a single frame of the proposed output ramp
        #based on the size of the croped dark current integration
        yd, xd = frame.shape
        #self.frametime = (xd/self.params['Readout']['namp'] + 12.) * (yd+1) * 10.00 * 1.e-6
        #UPDATED VERSION, 16 Sept 2017
        colpad = 12
        rowpad = 2
        if ((xd <= 8) & (yd <= 8)):
            rowpad = 3
        self.frametime = ((1.0 * xd/self.params['Readout']['namp'] + colpad) * (yd+rowpad)) * 1.e-5
        print('Exposure time of a single frame: ', self.frametime)


    def getSubarrayBounds(self):
        # Find the bounds of the requested subarray
        if self.params['Readout']['array_name'] in self.subdict['AperName']:
            mtch = self.params['Readout']['array_name'] == self.subdict['AperName']
            self.subarray_bounds = [self.subdict['xstart'].data[mtch][0], self.subdict['ystart'].data[mtch][0], self.subdict['xend'].data[mtch][0], self.subdict['yend'].data[mtch][0]]
            self.refpix_pos = {'x':self.subdict['refpix_x'].data[mtch][0], 'y':self.subdict['refpix_y'][mtch][0], 'v2':self.subdict['refpix_v2'].data[mtch][0], 'v3':self.subdict['refpix_v3'].data[mtch][0]}

            namps = self.subdict['num_amps'].data[mtch][0]
            if namps != 0:
                self.params['Readout']['namp'] = namps
            else:
                if ((self.params['Readout']['namp'] == 1) or (self.params['Readout']['namp'] == 4)):
                    print(("CAUTION: Aperture {} can be used with either "
                           "a 1-amp"
                           .format(self.subdict['AperName'].data[mtch][0])))
                    print("or a 4-amp readout. The difference is a factor of 4 in")
                    print(("readout time. You have requested {} amps."
                           .format(self.params['Readout']['namp'])))
                else:
                    print(("WARNING: {} requires the number of amps to "
                           "be 1 or 4. You have requested {}."
                           .format(self.params['Readout']['array_name']
                                   , self.params['Readout']['namp'])))
                    sys.exit()

        else:
            print(("WARNING: subarray name {} not found in the subarray "
                   "dictionary {}."
                   .format(self.params['Readout']['array_name']
                           , self.params['Reffiles']['subarray_defs'])))
            sys.exit()


    def getCRrate(self):
        #get the base cosmic ray impact probability
        #these numbers are from Kevin's original script. Not sure where
        #they come from or their units
        self.crrate=0.
        if "SUNMAX" in self.params["cosmicRay"]["library"]:
            self.crrate=1.6955e-04
        if "SUNMIN" in self.params["cosmicRay"]["library"]:
            self.crrate=6.153e-05
        if "FLARES" in self.params["cosmicRay"]["library"]:
            self.crrate=0.10546

        self.crrate=self.crrate/self.frametime
        #print('self.crrate is {}'.format(self.crrate))
        if self.crrate > 0.:
            print("Base cosmic ray probability per pixel per second: {}".format(self.crrate))


    def addSyntheticToDark(self, synthetic, dark, syn_zeroframe=None):
        """Add the synthetic data (now an exposure) to the dark current
        exposure.

        If zeroframe is provided, the function uses that to create the
        dark+synthetic zeroframe that is returned. If not provided, the
        function attempts to use the 0th frame of the input synthetic ramp

        Combine the cube of synthetic signals to the real dark current ramp.
        Be sure to adjust the dark current ramp if nframe/nskip is different
        than the nframe/nskip values that the dark was taken with.

        Only RAPID, NISRAPID darks will be re-averaged into different 
        readout patterns. But a BRIGHT2 dark can be used to create a 
        BRIGHT2 simulated ramp

        Arguments:
        ----------
        synthetic -- simulated signals, 4D array
        dark -- dark current exposure, 4D array
        syn_zeroframe -- zeroframe data associated with simulated data

        Returns:
        --------
        4D exposure containing combined simulated + dark data
        """
        # Get the info for the dark integration
        darkpatt = dark.header['READPATT']
        dark_nframe = dark.header['NFRAMES']
        mtch = self.readpatterns['name'].data == darkpatt
        dark_nskip = self.readpatterns['nskip'].data[mtch][0]

        # If the zeroframes for the dark and the synthetic data
        # are present, combine. Otherwise the zeroframe will be
        # None.
        zeroframe = None
        if ((syn_zeroframe is not None) & (dark.zeroframe is not None)):
            zeroframe = dark.zeroframe + syn_zeroframe

        # To hold reordered superbias + refpix signals from the dark
        reorder_sbandref = np.zeros_like(synthetic)

        # We have already guaranteed that either the readpatterns match
        # or the dark is RAPID, so no need to worry about checking for
        # other cases here.
        rapids = ["RAPID", "NISRAPID"]
        if ((darkpatt in rapids) and (self.params['Readout']['readpatt'] not in rapids)):
            deltaframe = self.params['Readout']['nskip'] + \
                         self.params['Readout']['nframe']
            frames = np.arange(0, self.params['Readout']['nframe'])
            accumimage = np.zeros_like(synthetic[0, :, :], dtype=np.int32)
            sbaccumimage = np.zeros_like(synthetic[0, :, :], dtype=np.int32)

            # Loop over integrations
            for integ in range(self.params['Readout']['nint']):

                # Loop over groups
                for i in range(self.params['Readout']['ngroup']):
                    # average together the appropriate frames,
                    # skip the appropriate frames
                    print(('Averaging dark current ramp in addSyntheticToDark.'
                           'Frames {}, to become group {}'.format(frames, i)))

                    # If averaging needs to be done
                    if self.params['Readout']['nframe'] > 1:
                        accumimage = np.mean(dark.data[integ, frames, :, :], axis=0)
                        sbaccumimage = np.mean(dark.sbAndRefpix[integ, frames, :, :],
                                               axis=0)

                        # If no averaging needs to be done
                    else:
                        accumimage = dark.data[integ, frames[0], :, :]
                        sbaccumimage = dark.sbAndRefpix[integ, frames[0], :, :]

                    # Now add the averaged dark frame to the synthetic data,
                    # which has already been placed into the correct readout pattern
                    synthetic[integ, i, :, :] += accumimage
                    reorder_sbandref[integ, i, :, :] = sbaccumimage

                    # Increment the frame indexes
                    frames = frames + deltaframe

        elif (darkpatt == self.params['Readout']['readpatt']):
            # If the input dark is not RAPID, or if the readout pattern
            # of the input dark and the output ramp match, then no averaging
            # needs to be done and we can simply add the synthetic groups to
            # the dark current groups.
            synthetic = synthetic + dark.data[:, 0:self.params['Readout']['ngroup'], :, :]
            print("Number of pixels with exactly 0 signal in synthetic: {}".format(np.sum(synthetic==0)))
            reorder_sbandref = dark.sbAndRefpix
        return synthetic, zeroframe, reorder_sbandref

    def maskRefPix(self, ramp, zero):
        # Make sure that reference pixels have no signal
        # in the simulated source ramp
        maskimage = np.zeros((self.ffsize, self.ffsize), dtype=np.int)
        maskimage[4:self.ffsize - 4, 4:self.ffsize - 4] = 1.

        #crop the mask to match the requested output array
        if "FULL" not in self.params['Readout']['array_name']:
            maskimage = self.crop_to_subarray(maskimage)

        ramp *= maskimage
        zero *= maskimage
        return ramp, zero

    def add_flatfield_effects(self,ramp):
        """Add flat field effects to the exposure"""
        # ILLUMINATION FLAT
        if self.runStep['illuminationflat']:
            illuminationflat, illuminationflatheader = self.readCalFile(self.params['Reffiles']['illumflat'])
            ramp *= illuminationflat

        #PIXEL FLAT
        if self.runStep['pixelflat']:
            pixelflat, pixelflatheader = self.readCalFile(self.params['Reffiles']['pixelflat'])
            ramp *= pixelflat
        return ramp

    def add_detector_effects(self, ramp):
        if self.runStep['crosstalk']:
            ramp = self.addCrosstalk(ramp)

        if self.runStep['pixelAreaMap']:
            ramp = self.addPAM(ramp)

        return ramp

    #def addIPC(self, exposure):
    #    # Input is always 4D
    #    ints, groups, yd, xd = exposure.shape

    #    # Prepare IPC effects
    #    ipcimage = fits.getdata(self.params['Reffiles']['ipc'])

    #    # If the IPC kernel is designed for the
    #    # removal of IPC, we need to invert it
    #    if self.params['Reffiles']['invertIPC']:
    #        print("Inverting IPC kernel prior to convolving with image")
    #        yk, xk = ipcimage.shape
    #        newkernel = 0. - ipcimage
    #        newkernel[int((yk - 1) / 2), int((xk - 1) / 2)] = 1. - (ipcimage[1, 1] - np.sum(ipcimage))
    #        ipcimage = newkernel

    #    # Apply the kernel to each frame of each integration
    #    for integ in range(ints):
    #        for group in range(groups):
    #            exposure[integ, group, :, :] = s1.fftconvolve(exposure[integ, group, :, :], ipcimage, mode="same")
    #    return exposure

    def addIPC(self, data):
        """
        Add interpixel capacitance effects to the data. This is done by 
        convolving the data with a kernel. The kernel is read in from the
        file specified by self.params['Reffiles']['ipc']. The core of this
        function was copied from the IPC correction step in the JWST
        calibration pipeline.

        Parameters:
        -----------
        data : obj
            4d numpy ndarray containing the data to which the 
            IPC effects will be added

        Returns:
        --------
        returns : obj
            4d numpy ndarray of the modified data
        """
        output_data = np.copy(data)
        # Shape of the data, which may include reference pix
        shape = output_data.shape

        # Find the number of reference pixel rows and columns
        # in output_data
        if self.subarray_bounds[0] < 4:
            left_columns = 4 - self.subarray_bounds[0]
        else:
            left_columns = 0
        if self.subarray_bounds[2] > 2043:
            right_columns = 4 - (2047 - self.subarray_bounds[2])
        else:
            right_columns = 0
        if self.subarray_bounds[1] < 4:
            bottom_rows = 4 - self.subarray_bounds[1]
        else:
            bottom_rows = 0
        if self.subarray_bounds[3] > 2043:
            top_rows = 4 - (2047 - self.subarray_bounds[3])
        else:
            top_rows = 0

        # Get IPC kernel data
        try:
            # If addIPC has already been called, then the correct
            # IPC kernel already exists, in self.kernel
            kernel = np.copy(self.kernel)
        except:
            # If addIPC has not been called yet, then read in the
            # kernel from the specified file.
            kernel = fits.getdata(self.params['Reffiles']['ipc'])
            # Invert the kernel if requested, to go from a kernel
            # designed to remove IPC effects to one designed to
            # add IPC effects
            if self.params['Reffiles']['invertIPC']:
                print("Iverting IPC kernel prior to convolving with image")
                kernel = self.invert_ipc_kernel(kernel)
            self.kernel = np.copy(kernel)
        kshape = kernel.shape

        # These axes lengths exclude reference pixels, if there are any.
        ny = shape[-2] - (bottom_rows + top_rows)
        nx = shape[-1] - (left_columns + right_columns)

        # The temporary array temp is larger than the science part of
        # output_data by a border (set to zero) that's about half of the
        # kernel size, so the convolution can be done without checking for
        # out of bounds.
        # b_b, t_b, l_b, and r_b are the widths of the borders on the
        # bottom, top, left, and right, respectively.
        b_b = kshape[0] // 2
        t_b = kshape[0] - b_b - 1
        l_b = kshape[1] // 2
        r_b = kshape[1] - l_b - 1
        tny = ny + b_b + t_b
        yoff = bottom_rows           # offset in output_data
        tnx = nx + l_b + r_b
        xoff = left_columns          # offset in output_data

        # Loop over integrations and groups
        for integration in range(shape[0]):
            for group in range(shape[1]):
        
                # Copy the science portion (not the reference pixels) of
                # output_data to this temporary array, then make
                # subsequent changes in-place to output_data.
                temp = np.zeros((tny, tnx), dtype=output_data.dtype)
                temp[b_b:b_b + ny, l_b:l_b + nx] = \
                    output_data[integration, group, yoff:yoff + ny, xoff:xoff + nx].copy()

                # After setting this slice to zero, we'll incrementally add
                # to it.
                output_data[integration, group, yoff:yoff + ny, xoff:xoff + nx] = 0.

                if len(kshape) == 2:
                    # 2-D IPC kernel.  Loop over pixels of the deconvolution
                    # kernel. In this section, `part` has the same shape
                    # as `temp`.
                    middle_j = kshape[0] // 2
                    middle_i = kshape[1] // 2
                    for j in range(kshape[0]):
                        jstart = kshape[0] - j - 1
                        for i in range(kshape[1]):
                            if i == middle_i and j == middle_j:
                                continue  # the middle pixel is done last
                            part = kernel[j, i] * temp
                            istart = kshape[1] - i - 1
                            output_data[integration, group, yoff:yoff + ny, xoff:xoff + nx] += \
                                        part[jstart:jstart + ny, istart:istart + nx]
                    # The middle pixel of the IPC kernel is expected to be
                    # the largest, so add that last.
                    part = kernel[middle_j, middle_i] * temp
                    output_data[integration, group, yoff:yoff + ny, xoff:xoff + nx] += \
                                part[middle_j:middle_j + ny, middle_i:middle_i + nx]

                else:
                    # 4-D IPC kernel.  Extract a subset of the kernel:
                    # all of the first two axes, but only the portion
                    # of the last two axes corresponding to the science
                    # data (i.e. possibly a subarray,
                    # and certainly excluding reference pixels).
                    k_temp = np.zeros((kshape[0], kshape[1], tny, tnx),
                                      dtype=kernel.dtype)
                    k_temp[:, :, b_b:b_b + ny, l_b:l_b + nx] = \
                            kernel[:, :, yoff:yoff + ny, xoff:xoff + nx]

                    # In this section, `part` has shape (ny, nx), which is
                    # smaller than `temp`.
                    middle_j = kshape[0] // 2
                    middle_i = kshape[1] // 2
                    for j in range(kshape[0]):
                        jstart = kshape[0] - j - 1
                        for i in range(kshape[1]):
                            if i == middle_i and j == middle_j:
                                continue   # the middle pixel is done last
                            istart = kshape[1] - i - 1
                            # The slice of k_temp includes different pixels
                            # for the first or second axes within each loop,
                            # but the same slice for the last two axes.
                            # The slice of temp (a copy of the science data)
                            # includes a different offset for each loop.
                            part = k_temp[j, i, b_b:b_b + ny, l_b:l_b + nx] * \
                                   temp[jstart:jstart + ny, istart:istart + nx]
                            output_data[integration, group, yoff:yoff + ny, xoff:xoff + nx] += part
                    # Add the product for the middle pixel last.
                    part = k_temp[middle_j, middle_i, b_b:b_b + ny, l_b:l_b + nx] * \
                           temp[middle_j:middle_j + ny, middle_i:middle_i + nx]
                    output_data[integration, group, yoff:yoff + ny, xoff:xoff + nx] += part
        return output_data

    def invert_ipc_kernel(self, kern):
        """
        Invert the IPC kernel such that it goes from being used to remove
        IPC effects from data, to being used to add IPC effects to data,
        or vice versa.

        Parameters:
        -----------
        kern : obj
            numpy ndarray, either 2D or 4D, containing the kernel

        Returns:
        --------
        returns : obj
            numpy ndarray containing iInverted" kernel
        """
        shape = kern.shape
        ys = 0
        ye = shape[-2]
        xs = 0
        xe = shape[-1]
        if shape[-1] == 2048:
            xs = 4
            xe = 2044
        if shape[-2] == 2048:
            ys = 4
            ye = 2044
        if len(shape) == 2:
            subkernel = kern[ys:ye, xs:xe]
        elif len(shape) == 4:
            subkernel = kern[:, : , ys:ye, xs:xe]    

        dims = subkernel.shape
        # Force subkernel to be 4D to make the function cleaner
        # Dimensions are (kernely, kernelx, detectory, detectorx)
        if len(dims) == 2:
            subkernel = np.expand_dims(subkernel, axis=3)
            subkernel = np.expand_dims(subkernel, axis=4)
            dims = subkernel.shape
            
        # Make sure the total signal in the kernel = 1
        #ave = np.average(subkernel, axis=(0, 1))
        #npix = dims[0] * dims[1]
        #tflux = ave * npix
        #renorm = 1. / tflux
        #subkernel *= renorm

        delta = subkernel * 0.
        delta[nyc, nxc] = 1.
        a1 = np.fft.fft2(subkernel, axes=(0, 1)) 
        a2 = np.fft.fft2(delta, axes=(0, 1))
        aout = a2 / a1
        imout = np.fft.ifft2(aout, axes=(0, 1))
        realout = np.real(imout)
        imout1 = np.fft.fftshift(imout, axes=(0, 1))
        realout1 = np.real(imout1)

        # If the input kernel was 2D, make the output 2D
        # If the input was 4D and had reference pixels, then
        # surround the inverted kernel with reference pixels
        if len(shape) == 2:
            newkernel = realout1[:, :, 0, 0]
        elif len(shape) == 4:
            newkernel = np.copy(kern)
            newkernel[:, :, ys:ye, xs:xe] = realout1

        # Save the inverted kernel for future simulator runs
        h0 = fits.PrimaryHDU()
        h1 = fits.ImageHDU(newkernel)
        h1.header["DETECTOR"] = self.detector
        h1.header["INSTRUME"] = self.params["Inst"]["instrument"]
        hlist = fits.HDUList([h0, h1])
        indir, infile = os.path.split(self.params["Reffiles"]["ipc"])
        outname = os.path.join(indir, "Kernel_to_add_IPC_effects_from_" + infile)
        hlist.writeto(outname, overwrite=True)
        print(("Inverted IPC kernel saved to {} for future simulator "
               "runs.".format(outname)))
        return newkernel

    def addCrosstalk(self, exposure):
        # Input is always 4D
        ints, groups, yd, xd = exposure.shape
        if self.params['Readout']['namp'] == 4:
            if self.instrument.upper() == 'NIRCAM':
                xdet = self.detector[3:5].upper()
                if xdet[1] == 'L':
                    xdet = xdet[0] + '5'
            else:
                xdet = self.detector
            xtcoeffs = self.readCrossTalkFile(self.params['Reffiles']['crosstalk'], xdet)
            # Only sources on the detector will create crosstalk.
            # If signalimage is larger than full frame
            # because we are creating a grism image, then extract the
            # pixels corresponding to the actual
            # detector, and only create crosstalk values for those.
            xs = 0
            xe = xd
            ys = 0
            ye = yd
            #if self.params['Output']['grism_source_image']:
            #    xs, xe, ys, ye = self.extractFromGrismImage(exposure[0, 0, :, :])
            for integ in range(ints):
                for group in range(groups):
                    xtinput = exposure[integ, group, ys:ye, xs:xe]
                    xtimage = self.crossTalkImage(xtinput, xtcoeffs)

                    #Now add the crosstalk image to the signalimage
                    exposure[integ, group, ys:ye, xs:xe] += xtimage
        else:
            print("Crosstalk calculation requested, but the chosen subarray")
            print("is read out using only 1 amplifier.")
            print("Therefore there will be no crosstalk. Skipping this step.")
        return exposure


    def readCrossTalkFile(self, file, detector):
        # Read in appropriate line from the xtalk coefficients
        # file and return the coeffs
        xtcoeffs = ascii.read(file, header_start=0)

        coeffs = []
        mtch = xtcoeffs['Det'] == detector.upper()
        if np.any(mtch) == False:
            print('Detector {} not found in xtalk file {}'.format(detector, file))
            sys.exit()

        return xtcoeffs[mtch]


    def extractFromGrismImage(self, array):
        # Return the indexes that will allow you to extract the nominal output
        # image from the extra-large grism source image
        print("we shouldn't need to use this function")
        sys.exit()
        arrayshape = array.shape
        nominaly, nominalx = self.dark.data[0, 0, :, :].shape
        diffx = (arrayshape[1] - nominalx) / 2
        diffy = (arrayshape[0] - nominaly) / 2
        x1 = diffx
        y1 = diffy
        x2 = arrayshape[1] - diffx
        y2 = arrayshape[0] - diffy
        return x1, x2, y1, y2


    def crossTalkImage(self, orig, coeffs):
        # Using Xtalk coefficients, generate an image of the crosstalk signal
        xtalk_corr_im = np.zeros_like(orig)
        subamp_shift = {"0":1, "1":-1, "2":1, "3":-1}

        #List of starting columns for all quadrants.
        xtqstart = [0, 512, 1024, 1536, 2048]

        for amp in range(4):
            to_mult = orig[:, xtqstart[amp]:xtqstart[amp+1]]
            receivers = []
            for i in range(4):
                if i != amp:
                    receivers.append(i)
            # Reverse the values to multply if the amps being used
            # are adjacent or 3 amps apart
            for subamp in receivers:
                index = 'xt'+str(amp+1)+str(subamp+1)
                if ((np.absolute(amp-subamp) == 1) | (np.absolute(amp-subamp) == 3)):
                    corr_amp = np.fliplr(to_mult) * coeffs[index]
                if (np.absolute(amp-subamp) == 2):
                    corr_amp = to_mult * coeffs[index]

                xtalk_corr_im[:, xtqstart[subamp]:xtqstart[subamp+1]] += corr_amp

            # Per Armin's instructions, now repeat the process
            # using his xt??post coefficients, but shift the arrays
            # by one pixel according to readout direction.
            for subamp in receivers:
                index = 'xt'+str(amp+1)+str(subamp+1)+'post'
                if ((np.absolute(amp-subamp) == 1) | (np.absolute(amp-subamp) == 3)):
                    corr_amp = np.fliplr(to_mult) * coeffs[index]
                    corr_amp = np.roll(corr_amp, subamp_shift[str(subamp)], axis=1)
                if (np.absolute(amp-subamp) == 2):
                    corr_amp = to_mult * coeffs[index]
                    corr_amp = np.roll(corr_amp, subamp_shift[str(subamp)])

                xtalk_corr_im[:, xtqstart[subamp]:xtqstart[subamp+1]] += corr_amp

        # Save the crosstalk correction image
        if self.params['Output']['save_intermediates'] == True:
            phdu = fits.PrimaryHDU(xtalk_corr_im)
            base_name = self.params['Output']['file'].split('/')[-1]
            xtalkout = os.path.join(self.params['Output']['directory'], base_name[0:-5] + '_xtalk_correction_image.fits')
            phdu.writeto(xtalkout, overwrite=True)

        return xtalk_corr_im


    def addPAM(self, signalramp):
        #Pixel Area Map
        pixAreaMap = self.simpleGetImage(self.params['Reffiles']['pixelAreaMap'])

        # If we are making a grism direct image, we need to embed the true pixel area
        # map in an array of the appropriate dimension, where any pixels outside the
        # actual aperture are set to 1.0
        if self.params['Output']['grism_source_image']:
            mapshape = pixAreaMap.shape
            g, yd, xd = signalramp.shape
            pam = np.ones((yd, xd))
            ys = self.coord_adjust['yoffset']
            xs = self.coord_adjust['xoffset']
            pam[ys:ys+mapshape[0], xs:xs+mapshape[1]] = np.copy(pixAreaMap)
            pixAreaMap = pam

        signalramp *= pixAreaMap
        return signalramp


    def simpleGetImage(self, name):
        # Read in an array from a fits file and crop using subarray_bounds
        try:
           image, header = fits.getdata(name, header=True)
        except:
            print('WARNING: unable to read in {}'.format(name))
            sys.exit()

        #assume that the input is 2D, since we are using it to build a signal rate frame
        imageshape = image.shape
        if len(imageshape) != 2:
            self.printfunc("Error: image %s is not two-dimensional" % (name))
            return None, None

        imageshape = image.shape

        try:
            image = image[self.subarray_bounds[1]:self.subarray_bounds[3]+1, self.subarray_bounds[0]:self.subarray_bounds[2]+1]
        except:
            print("Unable to crop image from {}".format(name))
            sys.exit

        return image

    def add_crs_and_noise(self,seed):
        """Given a noiseless seed ramp, add cosmic
        rays and poisson noise"""
        yd, xd = seed.shape[-2:]
        seeddim = len(seed.shape)

        # Run one integration at a time
        # because each needs its own collection
        # of cosmic rays and poisson noise realization
        nint = self.params['Readout']['nint']
        ngroups = self.params['Readout']['ngroup']
        sim_exposure = np.zeros((nint, ngroups, yd, xd))
        sim_zero = np.zeros((nint, yd, xd))

        for integ in range(nint):
            print("Integration {}:".format(integ))
            if seeddim == 2:
                inseed = seed
            elif seeddim == 4:
                inseed = seed[integ, :, :, :]
            if self.runStep['cosmicray']:
                ramp, rampzero = self.frameToRamp(inseed)
            else:
                ramp, rampzero = self.frameToRamp_noCR(inseed)
            sim_exposure[integ, :, :, :] = ramp
            sim_zero[integ, :, :] = rampzero
        return sim_exposure, sim_zero

    def frameToRamp(self, data):
        """Convert rate image to ramp, add poisson noise
        and cosmic rays

        Arguments:
        ----------
        data -- seed image. Should be a 2d frame or 3d integration.
        If the original seed image is a 4d exposure, call frameToRamp
        with one integration at a time.

        Returns:
        --------
        outramp -- 3d integration with cosmic rays and poisson noise
        zeroframe -- 2d zeroframe
        """

        # Output ramp will be in requested readout pattern!
        ndim = len(data.shape)

        #hopefully we don't need this and can find deltaimage on the fly...
        if ndim == 4:
            print("Shouldn't be here! No 4D seed images!")
            sys.exit()
            nintin, ngroupin, yd, xd = data.shape
        elif ndim == 3:
            ngroupin, yd, xd = data.shape
        #    deltaimage = np.zeros((ngroupin-1,yd,xd))
        #    for frame in xrange(1,ngroupin):
        #        deltaimage[frame-1,:,:] = data[frame,:,:] - data[frame-1,:,:]
        elif ndim == 2:
            yd, xd = data.shape
        #    ngroupin = None
        #    deltaimage = np.copy(data)
        #-----------------------

        # If a ramp is given, create a -1st frame that is all zeros
        # so that we can create deltaframes for all frames later
        # This should be the case only for data containing
        # moving targets.
        if ndim == 3:
            print('Moving target data shape', data.shape, yd, xd)
            #newdata = np.zeros((nintin, ngroupin+1, yd, xd))
            #newdata[:, 1:, :, :] = data
            #data = newdata
            data = np.vstack((np.zeros((1, yd, xd)), data))

        outramp = np.zeros((self.params['Readout']['ngroup'], yd, xd), dtype=np.float)
        #totalsignalimage = np.zeros((yd,xd),dtype=np.float)

        # Set up functions to apply cosmic rays later
        # Need the total number of active pixels in the
        # output array to multiply the CR rate by
        if self.runStep['cosmicray']:
            npix = int(yd * xd + 0.02)

            # Reinitialize the cosmic ray functions for each integration
            crhits, crs_perframe = self.CRfuncs(npix, seed=self.params['cosmicRay']['seed'])

            #open output file to contain the list of cosmic rays
            base_name = self.params['Output']['file'].split('/')[-1]
            crlistout = os.path.join(self.params['Output']['directory'], base_name[0:-5] + '_cosmicrays.list')
            self.openCRListFile(crlistout, crhits)

            #counter for use in cosmic ray addition while looping over frames
            #framenum = 0

        # Difference between the latest outimage frame and the
        # latest newsignalimage frame. This is important when nframe>1
        #delta = 0.
        if ndim == 2:
            totalsignalimage = data
        elif ndim == 3:
            totalsignalimage = data[1, :, :]

        # Define signal in the previous frame
        # Needed in loop below
        previoussignal = np.zeros((yd, xd))

        # Container for zeroth frame
        zeroframe = None

        # Total frames per group (including skipped frames)
        framesPerGroup = self.params['Readout']['nframe']+self.params['Readout']['nskip']

        # Loop over each group
        for i in range(self.params['Readout']['ngroup']):

            # Hold the averaged group signal
            accumimage = np.zeros((yd, xd))#, dtype=np.int32)

            # Group 0: the initial nskip frames don't exist,
            # so adjust indexes accordingly
            rstart = 0
            if i == 0:
                rstart = self.params['Readout']['nskip']

            # Loop over frames within each group if necessary
            for j in range(rstart,framesPerGroup):
                # Frame index number in input data
                frameindex = (i * framesPerGroup) + j - self.params['Readout']['nskip']

                # Signal only since previous frame
                if ndim == 3:
                    deltaframe = data[frameindex+1] - data[frameindex]
                elif ndim == 2:
                    deltaframe = data * self.frametime

                # Add poisson noise
                poissonsignal = self.doPoisson(deltaframe, self.params['simSignals']['poissonseed'])

                # Increment poisson seed value so that the next frame doesn't have identical
                # noise
                self.params['simSignals']['poissonseed'] += 1
                
                # Create the frame by adding the delta signal
                # and poisson noise associated with the delta signal
                # to the previous frame
                framesignal = previoussignal + poissonsignal

                # Add cosmic rays
                if self.runStep['cosmicray']:
                    framesignal = self.doCosmicRays(framesignal, i, j, self.params['Readout']['nframe'], crs_perframe[frameindex], self.params['cosmicRay']['seed'])
                    # Increment the seed, so that every frame doesn't have identical
                    # cosmic rays
                    self.params['cosmicRay']['seed'] += 1

                # Keep track of the total signal in the ramp,
                # so that we don't neglect signal which comes
                # in during the frames that are skipped.
                #totalsignalimage += framesignal
                previoussignal = copy.deepcopy(framesignal)

                if ((i==0) & (j==0)):
                    zeroframe = copy.deepcopy(framesignal)

                # Now round off and truncate to integers,
                # simulating the A/D conversion
                # workimage is the signal accumulated in the current frame
                # NOTE: any NaN values will be translated into -2147483648
                #framesignal = np.around(framesignal)
                #framesignal = framesignal.astype(np.int32)

                # Add the frame to the group signal image
                #if ((self.params['Readout']['nskip'] > 0) & (j >= self.params['Readout']['nskip'])):
                if j >= self.params['Readout']['nskip']:
                    print('    Averaging frame {} into group {}'.format(frameindex, i))
                    accumimage += framesignal
                elif j < self.params['Readout']['nskip']:
                #elif ((self.params['Readout']['nskip'] > 0) & (j < self.params['Readout']['nskip'])):
                    print('    Skipping frame {}'.format(frameindex))

            # divide by nframes if > 1
            if self.params['Readout']['nframe'] > 1:
                accumimage /= self.params['Readout']['nframe']
            outramp[i, :, :] = accumimage

        if self.runStep['cosmicray']:
            # Close the cosmic ray list file
            self.cosmicraylist.close()

        return outramp, zeroframe

    def frameToRamp_noCR(self, data):
        # Convert input seed image/ramp to a
        # ramp that includes poisson noise. No
        # cosmic rays are added

        # Output ramp will be in requested readout pattern!
        ndim = len(data.shape)


        if ndim == 3:
            ngroupin, yd, xd = data.shape
        elif ndim == 2:
            yd, xd = data.shape

        # Define output ramp
        outramp = np.zeros((self.params['Readout']['ngroup'], yd, xd))

        # If a ramp is given, create a -1st frame that is all zeros
        # so that we can create deltaframes for all frames later
        if ndim == 3:
            data = np.vstack(np.zeros((1, yd, xd)), data)

        # Container for zeroth frame
        zeroframe = None

        if ndim == 2:
            totalsignal = np.zeros((yd, xd))

        # Total frames per group (including skipped frames)
        framesPerGroup = self.params['Readout']['nframe']+self.params['Readout']['nskip']
        # Loop over each group
        for i in range(self.params['Readout']['ngroup']):
            accumimage = np.zeros((yd, xd))#, dtype=np.int32)

            # Loop over frames within each group if necessary
            # create each frame
            for j in range(framesPerGroup):

                # Frame index number in input data
                frameindex = (i * framesPerGroup) + j

                # Add poisson noise
                if ndim == 3:
                    framesignal = self.doPoisson(data[frameindex+1],
                                                 self.params['simSignals']['poissonseed'])
                elif ndim == 2:
                    framesignal = self.doPoisson(data*frameindex,
                                                 self.params['simSignals']['poissonseed'])

                # Increment poisson seed value so that the next frame doesn't have identical
                # noise
                self.params['simSignals']['poissonseed'] += 1

                # Keep track of the total signal in the ramp,
                # so that we don't neglect signal which comes
                # in during the frames that are skipped.
                #totalsignalimage += framesignal
                #previoussignal = copy.deepcopy(framesignal)

                if ((i==0) & (j==0)):
                    zeroframe = copy.deepcopy(framesignal)

                # Now round off and truncate to integers,
                # simulating the A/D conversion
                # workimage is the signal accumulated in the current frame
                # NOTE: any NaN values will be translated into -2147483648
                #framesignal = np.around(framesignal)
                #framesignal = framesignal.astype(np.int32)

                # Add the frame to the group signal image
                if ((self.params['Readout']['nskip'] > 0) & (j >= self.params['Readout']['nskip'])):
                    print('    Averaging frame {} into group {}'.format(frameindex, i))
                    accumimage += framesignal
                elif ((self.params['Readout']['nskip'] > 0) & (j < self.params['Readout']['nskip'])):
                    print('    Skipping frame {}'.format(frameindex))

            # divide by nframes if > 1
            if self.params['Readout']['nframe'] > 1:
                accumimage /= self.params['Readout']['nframe']
            outramp[i, :, :] = accumimage
        return outramp, zeroframe

    def doPoisson(self, signalimage, seedval):
        """Add poisson noise to an input image. Input is assumed
        to be in units of ADU, meaning it must be multiplied by
        the gain when calcuating Poisson noise. Then divide by the
        gain in order for the returned image to also be in ADU

        Arguments:
        ----------
        signalimage -- 2D array of signals in ADU
        seedval -- Integer seed value for the random number generator

        Returns:
        --------
        signalimage with Poisson noise added
        """

        #newimage = np.zeros_like(signalimage,dtype=np.float)
        #ndim = signalimage.shape

        # Set the seed
        np.random.seed(seedval)

        # Find the appropriate quantum yield value for the filter
        #if self.params['simSignals']['photonyield']:
        #    try:
        #        if self.params['Readout']['pupil'][0].upper() == 'F':
        #            usefilt = 'pupil'
        #        else:
        #            usefilt = 'filter'
        #        pym1=self.qydict[self.params['Readout'][usefilt]] - 1.
        #    except:
        #        pym1=0.

        # Quantum yield is 1.0 for all NIRCam filters
        pym1 = 0.

        # Can't add Poisson noise to pixels with negative values
        # Set those to zero when adding noise, then replace with
        # original value
        signalgain = signalimage * self.gainim
        highpix = np.where(signalgain == np.nanmax(signalgain))
        if np.nanmin(signalgain) < 0.:
            neg = signalgain < 0.
            negatives = copy.deepcopy(signalgain)
            negatives[neg] = signalgain[neg]
            signalgain[neg] = 0.
        
        # Add poisson noise    
        newimage = np.random.poisson(signalgain,signalgain.shape).astype(np.float)

        if np.nanmin(signalgain) < 0.:
            newimage[neg] = negatives[neg]

        newimage /= self.gainim

                # Quantum yield for NIRCam is always 1.0 (so psym1=0)
                #if self.params['simSignals']['photonyield'] and pym1 > 0.000001 and newimage[i, j] > 0:
                #    if self.params['simSignals']['pymethod']:
                #        # Calculate the values to make the poisson
                #        # results the same with/without photon
                #        # Yield (but not for pymethod true and false)
                #        # ...use yield -1 because the value
                #        # cannot be less than 1
                #        values = np.random.poisson(pym1, newimage[i, j])
                #        newimage[i, j] = newimage[i, j] + values.sum()
                #    else:
                #        newimage[i, j] = newimage[i, j] * self.qydict[self.params['Readout'][usefilt]]
                #        fract = newimage[i, j] - int(newimage[i, j])
                #        if self.generator2.random() < fract:
                #            newimage[i, j] = newimage[i, j] + 1
        return newimage

    def doCosmicRays(self, image, ngroup, iframe, nframe, ncr, seedval):
        # Change the seed each time this is run, or else simulated
        # exposures that have more than 1 integration will have the
        # same cosmic rays in each integration
        self.generator1 = random.Random()
        self.generator1.seed(seedval)

        # Add cosmic rays to a frame
        nray = int(ncr)

        i = 0
        dims = image.shape
        while i < nray:
            i = i+1
            j = int(self.generator1.random()*dims[0])
            k = int(self.generator1.random()*dims[1])
            n = int(self.generator1.random()*10.0)
            m = int(self.generator1.random()*1000.0)
            crimage = np.copy(self.cosmicrays[n][m, :, :])
            i1 = max(j-10, 0)
            i2 = min(j+11, dims[0])
            j1 = max(k-10, 0)
            j2 = min(k+11, dims[1])
            k1 = 10-(j-i1)
            k2 = 10+(i2-j)
            l1 = 10-(k-j1)
            l2 = 10+(j2-k)

            # Insert cosmic ray (divided by gain to put into ADU)
            image[i1:i2, j1:j2] = image[i1:i2, j1:j2] + crimage[k1:k2, l1:l2] / self.gainim[k1:k2, l1:l2]

            self.cosmicraylist.write("{} {} {} {} {} {} {}\n".format((j2-j1)/2+j1, (i2-i1)/2+i1, ngroup, iframe, n, m, np.max(crimage[k1:k2, l1:l2])))
        return image


    def openCRListFile(self, filename, hits):
        # Open a file and print header info for the file
        # that will contain the list and positions of inserted
        # cosmic rays
        self.cosmicraylist = open(filename, "w")
        self.cosmicraylist.write("# Cosmic ray list (file set %s random seed %d)\n" % (self.crfile, self.params['cosmicRay']['seed']))
        self.cosmicraylist.write('# Cosmic ray rate per frame: %13.6e (scale factor %f)\n' % (hits, self.params['cosmicRay']['scale']))
        self.cosmicraylist.write('Image_x    Image_y    Group   Frame   CR_File_Index   CR_file_frame   Max_CR_Signal\n')


    def CRfuncs(self, npix, seed=4242):
        # Set up functions that will be used to generate
        # cosmic ray hits
        crhits = npix * self.crrate * self.params['cosmicRay']['scale'] * self.frametime
        np.random.seed(seed)
        #self.generator1 = random.Random()
        #self.generator1.seed(seed)
        # Need a set of CRs for all frames, including those
        # that are skipped, in order for the rate of CRs to
        # be consistent.
        crs_perframe = np.random.poisson(crhits, self.params['Readout']['nint'] * self.params['Readout']['ngroup'] * (self.params['Readout']['nframe']+self.params['Readout']['nskip']))
        return crhits, crs_perframe


    def poissonFuncs(self, seed=4243):
        """Not used"""
        np.random.seed(seed)
        self.generator2 = random.Random()
        self.generator2.seed(seed)


    def readCRFiles(self):
        # Read in the 10 files that comprise the cosmic ray library
        self.cosmicrays=[]
        self.cosmicraysheader=[]
        for i in range(10):
            idx = '_%2.2d_' % (i)
            str1 = idx + self.params['cosmicRay']['suffix'] + '.fits'
            #str1='_%2.2d_IPC.fits' % (i)
            name = self.crfile + str1
            #name = self.params['cosmicRay']['path'] + self.fileList['cosmicrays'] + str1
            with fits.open(name) as h:
                im = h[1].data
                head = h[0].header
            #im, head = self.readCalFile(name)
            self.cosmicrays.append(im)
            self.cosmicraysheader.append(head)

    def readSaturationMap(self):
        # Read in saturation map
        try:
            self.satmap, self.satheader = self.readCalFile(self.params['Reffiles']['saturation'])
            bad = ~np.isfinite(self.satmap)
            self.satmap[bad] = 1.e6
        except:
            print(('WARNING: unable to open saturation file {}.'
                   .format(self.params['Reffiles']['saturation'])))
            print("Please provide a valid file, or place 'none' in the")
            print("saturation entry in the parameter file, ")
            print("in which case the nonlin limit value in the parameter file")
            print("({}) will be used for all pixels.".format(self.params['nonlin']['limit']))
            sys.exit()


    def readGainMap(self):
        # Read in the gain map. This will be used to
        # translate signals from e/s to ADU/sec
        if self.runStep['gain']:
            self.gainim, self.gainhead = self.readCalFile(self.params['Reffiles']['gain'])
            #set any NaN's to 1.0
            bad = ((~np.isfinite(self.gainim)) | (self.gainim == 0))
            self.gainim[bad] = 1.0

            #Pixels that have a gain value of 0
            # will be reset to have values of 1.0
            #zs = gainim == 0
            #gainim[zs] = 1.0


    def readCalFile(self, filename):
        #read in the specified calibration file
        try:
            with fits.open(filename) as h:
                image = h[1].data
                header = h[0].header
        except:
            print("WARNING: Unable to open {}".format(filename))
            sys.exit()

        #extract the appropriate subarray if necessary
        if ((self.subarray_bounds[0] != 0) or
            (self.subarray_bounds[2] != (self.ffsize - 1)) or
            (self.subarray_bounds[1] != 0) or
            (self.subarray_bounds[3] != (self.ffsize - 1))):

            if len(image.shape) == 2:
                image = image[self.subarray_bounds[1]:self.subarray_bounds[3] + 1,
                              self.subarray_bounds[0]:self.subarray_bounds[2] + 1]

            if len(image.shape) == 3:
                image = image[:, self.subarray_bounds[1]:self.subarray_bounds[3] + 1,
                              self.subarray_bounds[0]:self.subarray_bounds[2] + 1]

        return image, header


    def checkParams(self):
        # check instrument name
        if self.params['Inst']['instrument'].lower() not in INST_LIST:
            raise ValueError(("WARNING: instrument {} not in the list of "
                              "available instruments: {}"
                              .format(self.params['Inst']['instrument'].lower(),INST_LIST)))

        # check output filename - make sure it's fits
        if self.params['Output']['file'][-5:].lower() != '.fits':
            self.params['Output']['file'] += '.fits'

        # check mode:
        possibleModes = MODES[self.params['Inst']['instrument'].lower()]
        self.params['Inst']['mode'] = self.params['Inst']['mode'].lower()
        if self.params['Inst']['mode'] in possibleModes:
            pass
        else:
            raise ValueError(("WARNING: unrecognized mode {}. Must be one of: {}"
                              .format(self.params['Inst']['mode'], possibleModes)))

        # Make sure input readout pattern, nframe/nkip combination
        # is valid
        self.readpatternCheck()

        # Check that readout patterns of input dark and requested output
        # are compatible
        self.readpattern_compatible()
        
        # Make sure ngroup and nint are integers
        try:
            self.params['Readout']['ngroup'] = int(self.params['Readout']['ngroup'])
        except:
            print("WARNING: Input value of ngroup is not an integer.")
            sys.exit

        try:
            self.params['Readout']['nint'] = int(self.params['Readout']['nint'])
        except:
            print("WARNING: Input value of nint is not an integer.")
            sys.exit

        # Make sure that the requested number of groups is less than or
        # equal to the maximum allowed.
        # For full frame science operations, ngroup is going to be limited
        # to 10 for all readout patterns
        # except for the DEEP patterns, which can go to 20.
        match = self.readpatterns['name'] == self.params['Readout']['readpatt'].upper()
        if sum(match) == 1:
            if 'FULL' in self.params['Readout']['array_name']:
                maxgroups = self.readpatterns['maxgroups'].data[match][0]
            else:
                # I'm not sure what the limit is for subarrays, if any
                maxgroups = 999

        if (self.params['Readout']['ngroup'] > maxgroups):
            print(("WARNING: {} is limited to a maximum of {} groups. "
                   "Proceeding with ngroup = {}."
                   .format(self.params['Readout']['readpatt'], maxgroups, maxgroups)))
            self.params['Readout']['readpatt'] = maxgroups

        # Check for entries in the parameter file that are None or blank,
        # indicating the step should be skipped. Create a dictionary of steps
        # and populate with True or False
        self.runStep = {}
        self.runStep['superbias'] = self.checkRunStep(self.params['Reffiles']['superbias'])
        self.runStep['nonlin'] = self.checkRunStep(self.params['Reffiles']['linearity'])
        self.runStep['gain'] = self.checkRunStep(self.params['Reffiles']['gain'])
        #self.runStep['phot'] = self.checkRunStep(self.params['Reffiles']['phot'])
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
        self.runStep['linearity'] = self.checkRunStep(self.params['Reffiles']['linearity'])
        self.runStep['cosmicray'] = self.checkRunStep(self.params['cosmicRay']['path'])
        self.runStep['saturation_lin_limit'] = self.checkRunStep(self.params['Reffiles']['saturation'])
        self.runStep['fwpw'] = self.checkRunStep(self.params['Reffiles']['filtpupilcombo'])
        self.runStep['linearized_darkfile'] = self.checkRunStep(self.params['Reffiles']['linearized_darkfile'])
        self.runStep['badpixfile'] = self.checkRunStep(self.params['Reffiles']['badpixmask'])
        self.runStep['pixelAreaMap'] = self.checkRunStep(self.params['Reffiles']['pixelAreaMap'])

        # NON-LINEARITY
        # Make sure the input accuracy is a float with reasonable bounds
        self.params['nonlin']['accuracy'] = self.checkParamVal(self.params['nonlin']['accuracy'], 'nlin accuracy', 1e-12, 1e-6, 1e-6)
        self.params['nonlin']['maxiter'] = self.checkParamVal(self.params['nonlin']['maxiter'], 'nonlin max iterations', 5, 40, 10)
        self.params['nonlin']['limit'] = self.checkParamVal(self.params['nonlin']['limit'], 'nonlin max value', 30000., 1.e6, 66000.)

        # Make sure the CR random number seed is an integer
        try:
            self.params['cosmicRay']['seed'] = int(self.params['cosmicRay']['seed'])
        except:
            self.params['cosmicRay']['seed'] = 66231289
            print(("ERROR: cosmic ray random number generator seed is bad. "
                   "Using the default value of {}."
                   .format(self.params['cosmicRay']['seed'])))

        # Also make sure the poisson random number seed is an integer
        try:
            self.params['simSignals']['poissonseed'] = int(self.params['simSignals']['poissonseed'])
        except:
            self.params['simSignals']['poissonseed'] = 815813492
            print(("ERROR: cosmic ray random number generator seed is bad. "
                   "Using the default value of {}."
                   .format(self.params['simSignals']['poissonseed'])))

        # COSMIC RAYS:
        # Generate the name of the actual CR file to use
        if self.params['cosmicRay']['path'] is None:
            self.crfile = None
        else:
            if self.params['cosmicRay']['path'][-1] != '/':
                self.params['cosmicRay']['path'] += '/'
            if self.params['cosmicRay']["library"].upper() in ["SUNMAX", "SUNMIN", "FLARES"]:
                self.crfile = self.params['cosmicRay']['path'] + "CRs_MCD1.7_"+self.params['cosmicRay']["library"].upper()
            else:
                self.crfile = None
                print(("Warning: unrecognised cosmic ray library {}"
                       .format(self.params['cosmicRay']["library"])))
                sys.exit()


        # Read in distortion and WCS-related data. These will be placed
        # in the header of the output file.
        ap_name = self.params['Readout']['array_name']

        #Read in the distortion coefficients file if present. These will provide a more exact
        #transform from RA, Dec to x, y than the astrometric distortion reference file above.
        #The file above can be off by ~20 pixels in the corners of the array. This file will give
        #exact answers
        if os.path.isfile(self.params['Reffiles']['distortion_coeffs']):
            distortionTable = ascii.read(self.params['Reffiles']['distortion_coeffs'],
                                         header_start=1, format='csv')
        else:
            print(("WARNING: Input distortion coefficients file {} "
                   "does not exist."
                   .format(self.params['Reffiles']['distortion_coeffs'])))
            sys.exit()

        self.x_sci2idl, self.y_sci2idl, self.v2_ref, self.v3_ref, self.parity, self.v3yang, self.xsciscale, self.ysciscale, self.v3scixang = self.getDistortionCoefficients(distortionTable, 'science', 'ideal', ap_name)

        # Convert the input RA and Dec of the pointing position into floats
        # check to see if the inputs are in decimal units or hh:mm:ss strings
        try:
            self.ra = float(self.params['Telescope']['ra'])
            self.dec = float(self.params['Telescope']['dec'])
        except:
            self.ra, self.dec = self.parseRADec(self.params['Telescope']['ra']
                                               , self.params['Telescope']['dec'])

        if abs(self.dec) > 90. or self.ra < 0. or self.ra > 360. or self.ra is None or self.dec is None:
            print("WARNING: bad requested RA and Dec {} {}".format(self.ra, self.dec))
            sys.exit()

        # Make sure the rotation angle is a float
        try:
            self.params['Telescope']["rotation"] = float(self.params['Telescope']["rotation"])
        except:
            print(("ERROR: bad rotation value {}, setting to zero."
                   .format(self.params['Telescope']["rotation"])))
            self.params['Telescope']["rotation"] = 0.

        # From the pointing info above, calculate the local roll angle
        # This is used only when saving the output
        self.local_roll = stp.compute_local_roll(
            self.params['Telescope']['rotation'], self.ra, self.dec, self.v2_ref, self.v3_ref)

        # Check that the various scaling factors are floats and
        # within a reasonable range
        self.params['cosmicRay']['scale'] = self.checkParamVal(self.params['cosmicRay']['scale'], 'cosmicRay', 0, 100, 1)
        self.params['simSignals']['extendedscale'] = self.checkParamVal(self.params['simSignals']['extendedscale'], 'extendedEmission', 0, 10000, 1)
        self.params['simSignals']['zodiscale'] = self.checkParamVal(self.params['simSignals']['zodiscale'], 'zodi', 0, 10000, 1)
        self.params['simSignals']['scatteredscale'] = self.checkParamVal(self.params['simSignals']['scatteredscale'], 'scatteredLight', 0, 10000, 1)

        # Make sure the requested output format is an allowed value
        if self.params['Output']['format'] not in ['DMS']:
            print(("WARNING: unsupported output format {} requested. "
                   "Possible options are {}."
                   .format(self.params['Output']['format'], ['DMS'])))
            sys.exit()

        # Check the output metadata, including visit and observation
        # numbers, obs_id, etc
        kwchecks = ['program_number', 'visit_number', 'visit_group',
                    'sequence_id', 'activity_id', 'exposure_number',
                    'observation_number', 'obs_id', 'visit_id']
        for quality in kwchecks:
            try:
                self.params['Output'][quality] = str(self.params['Output'][quality])
            except:
                print(("WARNING: unable to convert {} to string. "
                       "This is required.".format(self.params['Output'][quality])))
                sys.exit()


    def checkRunStep(self, filename):
        # Check to see if a filename exists in the parameter file.
        if ((len(filename) == 0) or (filename.lower() == 'none')):
            return False
        else:
            return True


    def readpatternCheck(self):
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
                  'Using nframe = {} and nskip = {}'
                   .format(self.params['Readout']['readpatt'],
                           self.params['Readout']['nframe'],
                           self.params['Readout']['nskip'])))
        else:
            # If the read pattern is not present in the definition file
            # then quit.
            print(("WARNING: the {} readout pattern is not defined in {}."
                   .format(self.params['Readout']['readpatt'],
                           self.params['Reffiles']['readpattdefs'])))
            print("Quitting.")
            sys.exit()


            # If read pattern is not present in the definition file but
            # the nframe/nskip combo is, then reset
            # readpatt to the appropriate value from the definition file
            #readpatt_nframe = self.readpatterns['nframe'].data
            #readpatt_nskip = self.readpatterns['nskip'].data
            #readpatt_name = self.readpatterns['name'].data
            #if self.params['Readout']['nframe'] in readpatt_nframe:
            #    nfmtch = self.params['Readout']['nframe'] == readpatt_nframe
            #    nskip_subset = readpatt_nskip[nfmtch]
            #    name_subset = readpatt_name[nfmtch]
            #    if self.params['Readout']['nskip'] in nskip_subset:
            #        finalmtch = self.params['Readout']['nskip'] == nskip_subset
            #        finalname = name_subset[finalmtch][0]
            #        print(("CAUTION: requested readout pattern {} not recognized."
            #               .format(self.params['Readout']['readpatt'])))
            #        print(("but the requested nframe/nskip combination ({},{}), "
            #               "matches those values for".
            #               format(self.params['Readout']['nframe'],
            #                      self.params['Readout']['nskip'])))
            #        print(("the {} readout pattern, as listed in {}."
            #              .format(finalname,self.params['Reffiles']['readpattdefs'])))
            #        print('Continuing on using that as the readout pattern.')
            #        self.params['Readout']['readpatt'] = finalname
            #    else:
            #        # Case where readpatt is not recognized, nframe is present
            #        # in the definition file, but nskip is not
            #        print(('Unrecognized readout pattern {}, and the input '
            #               'nframe/nskip combination {},{} does not'
            #               .format(self.params['Readout']['readpatt'],
            #                       self.params['Readout']['nframe'],
            #                       self.params['Readout']['nskip'])))
            #        print(('match any present in {}. This is not a valid NIRCam '
            #              'readout pattern. Quitting.'
            #              .format(self.params['Reffiles']['readpattdefs'])))
            #        sys.exit()
            #else:
            #    # Case where readpatt and nframe are not recognized
            #    print(('Unrecognized readout pattern {}, and the input '
            #           'nframe/nskip combination {},{} does not'
            #           .format(self.params['Readout']['readpatt'],
            #                   self.params['Readout']['nframe'],
            #                   self.params['Readout']['nskip'])))
            #    print(('match any present in {}. This is not a valid NIRCam '
            #           'readout pattern. Quitting.'
            #           .format(self.params['Reffiles']['readpattdefs'])))
            #    sys.exit()


    def checkParamVal(self, value, typ, vmin, vmax, default):
        # Make sure the input value is a float and between min and max
        try:
            value = float(value)
        except:
            print("WARNING: {} for {} is not a float.".format(value, typ))
            sys.exit()

        if ((value >= vmin) & (value <= vmax)):
            return value
        else:
            print(("ERROR: {} for {} is not within reasonable bounds. "
                   "Setting to {}".format(value, typ, default)))
            return default

    def add_options(self, parser=None, usage=None):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage,
                                             description='Simulate JWST ramp')
        parser.add_argument("paramfile", help = 'File describing the input parameters and instrument settings to use. (YAML format).')
        parser.add_argument("linDark", help = 'File containing linearized dark ramp.')
        parser.add_argument("seed", help = 'File containing seed image and segmentation map')
        return parser


if __name__ == '__main__':

    usagestring = ('USAGE: obs_generator.py inputs.yaml '
                   'lindark.fits seedimg.fits')

    obs = Observation()
    parser = obs.add_options(usage = usagestring)
    args = parser.parse_args(namespace = obs)
    obs.create()
