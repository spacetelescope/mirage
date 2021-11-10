#! /usr/bin/env python

'''
Module for preparing a given dark current exposure for
integration with a ramp of simulated sources created
using (e.g.) catalog_seed_image.py.

Authors:
--------

    - Bryan Hilbert, Kevin Volk

Use:
----

    This module can be imported as such:

    ::

        from mirage.dark.dark_prep import DarkPrep
        dark = DarkPrep()
        dark.paramfile = 'my_parameters.yaml'
        dark.prepare()
'''

import sys
import os
import argparse
import datetime
import logging
from math import floor
from glob import glob
import shutil

import yaml
import pkg_resources
import numpy as np
from astropy.io import fits, ascii
import astropy.units as u

import mirage
from mirage.logging import logging_functions
from mirage.utils import read_fits, utils, siaf_interface
from mirage.utils.constants import FGS1_DARK_SEARCH_STRING, FGS2_DARK_SEARCH_STRING, \
                                   LOG_CONFIG_FILENAME, STANDARD_LOGFILE_NAME
from mirage.utils.file_splitting import find_file_splits
from mirage.utils.timer import Timer
from mirage.reference_files import crds_tools


# Allowed instruments
INST_LIST = ['nircam', 'niriss', 'fgs']

classdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../'))
log_config_file = os.path.join(classdir, 'logging', LOG_CONFIG_FILENAME)
logging_functions.create_logger(log_config_file, STANDARD_LOGFILE_NAME)


class DarkPrep():
    def __init__(self, file_splitting=True, offline=False):
        """Instantiate DarkPrep object.

        Parameters
        ----------
        file_splitting : bool
            If True, the output dark current objects are split into
            segment files, following the logic of the JWST calibration
            pipeline. If False, no file-splitting is done. False should
            only be used in very limited cases, such as when creating new
            linearized dark current files to add to the collection of
            Mirage reference files

        offline : bool
            If True, the check for the existence of the MIRAGE_DATA
            directory is skipped. This is primarily for Travis testing
        """
        # Initialize the log using dictionary from the yaml file
        self.logger = logging.getLogger(__name__)

        self.offline = offline
        self.file_splitting = file_splitting

        # Locate the module files, so that we know where to look
        # for config subdirectory
        self.modpath = pkg_resources.resource_filename('mirage', '')

        # Get the location of the MIRAGE_DATA environment
        # variable, so we know where to look for darks, CR,
        # PSF files, etc later
        self.env_var = 'MIRAGE_DATA'
        datadir = utils.expand_environment_variable(self.env_var, offline=offline)

        # Check that CRDS-related environment variables are set correctly
        self.crds_datadir = crds_tools.env_variables()

        # Initialize timer
        self.timer = Timer()

    def check_params(self):
        """Check for acceptible values for the input parameters in the
        yaml file.
        """
        # Check instrument name
        if self.params['Inst']['instrument'].lower() not in INST_LIST:
            raise ValueError(("WARNING: {} instrument not implemented within "
                              "simulator".format(self.params['Inst']['instrument'])))

        # If user requests not to use the pipeline,
        # make sure the input is a linearized dark, and not
        # a raw dark
        if self.params['Inst']['use_JWST_pipeline'] is False:
            if self.params['Reffiles']['linearized_darkfile'].lower() == 'none':
                raise ValueError(("WARNING: You specified no use of the JWST pipeline, but "
                                  "have not provided a linearized dark file to use. Without the "
                                  "pipeline, a raw dark cannot be used."))

        # If no raw dark nor linearized dark is given, then there is no
        # input to work with.
        if (self.params['Reffiles']['linearized_darkfile'].lower() == 'none' and
            self.params['Reffiles']['dark'].lower() == 'none'):
            raise ValueError(("WARNING: No raw nor linearized dark file given. Unable to proceed."))

        # Make sure nframe, nskip, ngroup are all integers
        try:
            self.params['Readout']['nframe'] = int(self.params['Readout']['nframe'])
        except ValueError:
            self.logger.error("WARNING: Input value of nframe is not an integer.")

        try:
            self.params['Readout']['nskip'] = int(self.params['Readout']['nskip'])
        except ValueError:
            self.logger.error("WARNING: Input value of nskip is not an integer.")

        try:
            self.params['Readout']['ngroup'] = int(self.params['Readout']['ngroup'])
        except ValueError:
            self.logger.error("WARNING: Input value of ngroup is not an integer.")

        try:
            self.params['Readout']['nint'] = int(self.params['Readout']['nint'])
        except ValueError:
            self.logger.error("WARNING: Input value of nint is not an integer.")

        # Check for entries in the parameter file that are None or blank,
        # indicating the step should be skipped. Create a dictionary of steps
        # and populate with True or False
        self.runStep = {}
        # self.runStep['linearity'] = self.check_run_step(self.params['Reffiles']['linearity'])
        self.runStep['badpixmask'] = self.check_run_step(self.params['Reffiles']['badpixmask'])
        self.runStep['linearized_darkfile'] = self.check_run_step(self.params['Reffiles']['linearized_darkfile'])
        self.runStep['saturation_lin_limit'] = self.check_run_step(self.params['Reffiles']['saturation'])
        self.runStep['superbias'] = self.check_run_step(self.params['Reffiles']['superbias'])
        self.runStep['linearity'] = self.check_run_step(self.params['Reffiles']['linearity'])

        # If the pipeline is going to be used to create
        # the linearized dark current ramp, make sure the
        # specified configuration files are present.
        if not self.runStep['linearized_darkfile']:
            dqcheck = self.check_run_step(self.params['newRamp']['dq_configfile'])
            if not dqcheck:
                raise FileNotFoundError(("WARNING: DQ pipeline step configuration file not provided. "
                                         "This file is needed to run the pipeline."))
            satcheck = self.check_run_step(self.params['newRamp']['sat_configfile'])
            if not satcheck:
                raise FileNotFoundError(("WARNING: Saturation pipeline step configuration file not provided. "
                                         "This file is needed to run the pipeline."))
            sbcheck = self.check_run_step(self.params['newRamp']['superbias_configfile'])
            if not sbcheck:
                raise FileNotFoundError(("WARNING: Superbias pipeline step configuration file not provided. "
                                        "This file is needed to run the pipeline."))
            refpixcheck = self.check_run_step(self.params['newRamp']['refpix_configfile'])
            if not refpixcheck:
                raise FileNotFoundError(("WARNING: Refpix pipeline step configuration file not provided. "
                                         "This file is needed to run the pipeline."))
            lincheck = self.check_run_step(self.params['newRamp']['linear_configfile'])
            if not lincheck:
                raise FileNotFoundError(("WARNING: Linearity pipeline step configuration file not provided. "
                                         "This file is needed to run the pipeline."))

    def check_run_step(self, filename):
        """Check to see if a filename exists in the parameter file
        or if it is set to none.

        Parameters
        ----------
        filename : str
            Name of the entry in the yaml file

        Returns
        -------
        state : bool
            True if entry is a filename, False if it is non or empty
        """
        if ((len(filename) == 0) or (filename.lower() == 'none')):
            return False
        else:
            return True

    def crop_dark(self, model):
        """Cut the dark current array down to the size dictated
        by the subarray bounds

        Parameters
        ----------
        model : obj
            read_fits object holding dark current as well as superbias and refpix signal

        Returns
        -------
        model : obj
            read_fits object with cropped data arrays
        """
        modshape = model.data.shape
        yd = modshape[-2]
        xd = modshape[-1]

        if ((self.subarray_bounds[0] != 0) or (self.subarray_bounds[2] != (xd - 1))
           or (self.subarray_bounds[1] != 0) or (self.subarray_bounds[3] != (yd - 1))):

            if len(modshape) == 4:
                model.data = model.data[:, :, self.subarray_bounds[1]:self.subarray_bounds[3] + 1,
                                        self.subarray_bounds[0]:self.subarray_bounds[2] + 1]
                try:
                    model.sbAndRefpix = model.sbAndRefpix[:, :, self.subarray_bounds[1]:
                                                          self.subarray_bounds[3] + 1,
                                                          self.subarray_bounds[0]:self.subarray_bounds[2] + 1]
                except:
                    pass
                try:
                    model.zeroframe = model.zeroframe[:, self.subarray_bounds[1]:
                                                      self.subarray_bounds[3] + 1,
                                                      self.subarray_bounds[0]:self.subarray_bounds[2] + 1]
                except:
                    pass
            if len(modshape) == 3:
                model.data = model.data[:, self.subarray_bounds[1]:self.subarray_bounds[3] + 1,
                                        self.subarray_bounds[0]:self.subarray_bounds[2] + 1]
                try:
                    model.sbAndRefpix = model.sbAndRefpix[:, self.subarray_bounds[1]:
                                                          self.subarray_bounds[3] + 1,
                                                          self.subarray_bounds[0]:self.subarray_bounds[2] + 1]
                except:
                    pass
                try:
                    model.zeroframe = model.zeroframe[:, self.subarray_bounds[1]:
                                                      self.subarray_bounds[3] + 1,
                                                      self.subarray_bounds[0]:self.subarray_bounds[2] + 1]
                except:
                    pass
            elif len(modshape) == 2:
                model.data = model.data[self.subarray_bounds[1]:self.subarray_bounds[3] + 1,
                                        self.subarray_bounds[0]:self.subarray_bounds[2] + 1]
                try:
                    model.sbAndRefpix = model.sbAndRefpix[self.subarray_bounds[1]:self.subarray_bounds[3] + 1,
                                                          self.subarray_bounds[0]:self.subarray_bounds[2] + 1]
                except:
                    pass
                try:
                    model.zeroframe = model.zeroframe[self.subarray_bounds[1]:
                                                      self.subarray_bounds[3] + 1,
                                                      self.subarray_bounds[0]:self.subarray_bounds[2] + 1]
                except:
                    pass

        # Make sure that if the output is supposedly a
        # 4-amplifier file, that the number of
        # pixels in the x direction is a multiple of 4.
        nfast = self.subarray_bounds[2] - self.subarray_bounds[0] + 1
        nslow = self.subarray_bounds[3] - self.subarray_bounds[1] + 1
        nramp = (self.params['Readout']['nframe'] + self.params['Readout']['nskip']) * \
            self.params['Readout']['ngroup']

        if self.params['Readout']['namp'] != 4 and self.params['Readout']['namp'] != 1:
            raise ValueError(('ERROR: amplifier mode specified ({}) is not allowed'
                              .format(self.params['Readout']['namp'])))

        if self.params['Readout']['namp'] == 4:
            n = int(nfast / 4) * 4
            if n != nfast:
                raise ValueError(("ERROR: 4 amplifier mode specified but the number of pixels "
                                  "in the fast read direction ({}) is not a multiple of 4."
                                  .format(nfast)))
        return model

    def darkints(self):
        """Check the number of integrations in the dark
        current file and compare with the requested
        number of integrations. Add/remove integrations
        if necessary

        Parameters
        ----------

        None


        Returns
        -------

        None
        """
        ndarkints, ngroups, ny, nx = self.dark.data.shape
        reqints = self.params['Readout']['nint']

        if reqints <= ndarkints:
            # Fewer requested integrations than are in the
            # input dark.
            self.logger.info("{} output integrations requested.".format(reqints))
            self.dark.data = self.dark.data[0:reqints, :, :, :]
            if self.dark.sbAndRefpix is not None:
                self.dark.sbAndRefpix = self.dark.sbAndRefpix[0:reqints, :, :, :]
            if self.dark.zeroframe is not None:
                self.dark.zeroframe = self.dark.zeroframe[0:reqints, :, :]

        elif reqints > ndarkints:
            # More requested integrations than are in input dark.
            logger.info(('Requested output has {} integrations, while input dark has only {}. '
                         'Adding integrations to input dark by making copies of the input.'
                         .format(reqints, ndarkints)))
            self.integration_copy(reqints, ndarkints)

    def data_volume_check(self, obj):
        """Make sure that the input integration has
        enough frames/groups to create the requested
        number of frames/groups of the output

        Parameters
        ----------
        obj : obj
            Instance of read_fits class containing dark current data and other info

        Returns
        -------
        None
        """
        ngroup = int(self.params['Readout']['ngroup'])
        nframe = int(self.params['Readout']['nframe'])
        nskip = int(self.params['Readout']['nskip'])

        inputframes = obj.data.shape[1]
        if ngroup * (nskip + nframe) > inputframes:
            self.logger.warning(("Not enough frames in the input integration to "
                   "create the requested number of output groups. Input has "
                   "{} frames. Requested output is {} groups each created from "
                   "{} frames plus skipping {} frames between groups for a total "
                   "of {} frames."
                   .format(inputframes, ngroup, nframe, nskip, ngroup * (nframe + nskip))))
            self.logger.warning(("Making copies of {} dark current frames and adding them to "
                   "the end of the dark current integration."
                   .format(ngroup * (nskip + nframe) - inputframes)))

            # Figure out how many more frames we need,
            # in terms of how many copies of the original dark
            div = floor((ngroup * (nskip + nframe)) / inputframes)
            mod = (ngroup * (nskip + nframe)) % inputframes

            # If more frames are needed than there are frames
            # in the original dark, then make copies of the
            # entire thing as many times as necessary, adding
            # the signal from the previous final frame to each.
            for ints in range(div - 1):
                extra_frames = np.copy(obj.data)
                obj.data = np.hstack((obj.data, extra_frames + obj.data[:, -1, :, :]))
                if obj.sbAndRefpix is not None:
                    extra_sb = np.copy(obj.sbAndRefpix)
                    obj.sbAndRefpix = np.hstack((obj.sbAndRefpix, extra_sb + obj.sbAndRefpix[:, -1, :, :]))

            # At this point, if more frames are needed, but fewer
            # than an entire copy of self.dark.data, then add the
            # appropriate number of frames here.
            extra_frames = np.copy(obj.data[:, 1:mod + 1, :, :]) - obj.data[:, 0, :, :]
            obj.data = np.hstack((obj.data, extra_frames + obj.data[:, -1, :, :]))
            if obj.sbAndRefpix is not None:
                extra_sb = np.copy(obj.sbAndRefpix[:, 1:mod + 1, :, :]) - obj.sbAndRefpix[:, 0, :, :]
                obj.sbAndRefpix = np.hstack((obj.sbAndRefpix, extra_sb + obj.sbAndRefpix[:, -1, :, :]))

        elif ngroup * (nskip + nframe) < inputframes:
            # If there are more frames in the dark than we'll need,
            # crop the extras in order to reduce memory use
            obj.data = obj.data[:, 0:ngroup * (nskip + nframe), :, :]
            if obj.sbAndRefpix is not None:
                obj.sbAndRefpix = obj.sbAndRefpix[:, 0:ngroup * (nskip + nframe), :, :]
        obj.header['NGROUPS'] = ngroup * (nskip + nframe)

    def filecheck(self):
        """Make sure the requested input files exist
        or reference files, assume first that they are located in
        the directory tree under the datadir (from the MIRAGE_DATA
        environment variable). If not, assume the input is a full path
        and check there.
        """
        rlist = [['Reffiles', 'dark'],
                 ['Reffiles', 'linearized_darkfile'],
                 ['Reffiles', 'superbias'],
                 ['Reffiles', 'linearity'],
                 ['Reffiles', 'saturation']]
        for ref in rlist:
            self.ref_check(ref)

    def get_base_dark(self, input_file):
        """Read in the dark current ramp that will serve as the
        base for the simulated ramp"""

        # First make sure that the file exists
        local = os.path.isfile(input_file)
        # Future improvement: download a missing dark file
        if not local:
            try:
                self.logger.info("Local copy of dark current file {} not found. Attempting to download".format(input_file))
                darkfile = os.path.split(input_file)[1]
                hh = fits.open(os.path.join(self.dark_url, darkfile))
                hh.writeto(input_file)
            except Exception:
                self.logger.error("Local copy of dark current file {} not found, and unable to download. Quitting.".format(input_file))

        self.dark = read_fits.Read_fits()
        self.dark.file = input_file

        # Depending on the method indicated in the input file
        # read in the dark using astropy or RampModel
        if self.params['Inst']['use_JWST_pipeline']:
            self.dark.read_datamodel()

            try:
                # Remove non-pipeline related keywords
                # (e.g. CV3 temps/voltages)
                self.dark.__delattr__('extra_fits')
            except:
                pass

        else:
            self.dark.read_astropy()

        # We assume that the input dark current integration is raw, which means
        # the data are in the original ADU measured by the detector. So the
        #  data range should be no larger than 65536
        darkrange = 1.*self.dark.data.max() - 1.*self.dark.data.min()
        if darkrange > 65535.:
            raise ValueError(("WARNING: Range of data values in the input dark is too large. "
                              "We assume the input dark is raw ADU values, with a range no more than 65536."))

        # If the inputs are signed integers, change to unsigned.
        if self.dark.data.min() < 0.:
            self.dark.data += 32768

        # If the input is any readout pattern other than RAPID, then
        # make sure that the output readout patten matches. Only RAPID
        # can be averaged and transformed into another readout pattern
        if self.dark.header['READPATT'] not in ['RAPID', 'NISRAPID', 'FGSRAPID']:
            if self.params['Readout']['readpatt'].upper() != self.dark.header['READPATT']:
                self.logger.error(("WARNING: cannot transform input {} integration into "
                                   "output {} integration.".format(self.dark.header['READPATT'],
                                                                   self.params['Readout']['readpatt'])))
                raise ValueError(("Only RAPID, NISRAPID, or FGSRAPID inputs can be translated to a "
                                  "different readout pattern"))
        else:
            pass

        # Finally, collect information about the detector, which will be needed for astrometry later
        self.detector = self.dark.header['DETECTOR']
        self.instrument = self.dark.header['INSTRUME']
        self.fastaxis = self.dark.header['FASTAXIS']
        self.slowaxis = self.dark.header['SLOWAXIS']

    def get_frame_count_info(self):
        """Calculate information on the number of frames per group and
        per integration
        """
        self.numints = self.params['Readout']['nint']
        self.numgroups = self.params['Readout']['ngroup']
        numframes = self.params['Readout']['nframe']
        numskips = self.params['Readout']['nskip']
        numresets = self.params['Readout']['resets_bet_ints']

        self.frames_per_group = numframes + numskips
        self.frames_per_integration = self.numgroups * self.frames_per_group
        self.total_frames = self.numgroups * self.frames_per_group

        if self.numints > 1:
            # Frames for all integrations
            self.total_frames *= self.numints
            # Don't worry about counting reset frames between integrations
            # We're not concerned with timing here.

    def get_reffile_metadata(self, keyword):
        """Reetrive the name of the reference file used to calibate
        self.linDark

        Parameters
        ----------
        keyword : str
            Header keyword to be retrived (e.g. 'R_LINEAR')

        Returns
        -------
        value : str
            Header keyword value
        """
        try:
            value = self.linDark.header[keyword]
        except KeyError:
            value = 'N/A'
        return value

    def integration_copy(self, req, darkint):
        """Use copies of integrations in the dark current input to
        make up integrations in the case where the output has
        more integrations than the input

        Parameters
        ----------
        req : int
            Requested number of dark current integrations for the output

        darkint : int
            Number of integrations in the input dark current exposure

        Returns
        -------
        None
        """

        ncopies = np.int((req - darkint) / darkint)
        extras = (req - darkint) % darkint

        # Full copies of the dark exposure (all integrations)
        if ncopies > 0:
            copydata = self.dark.data
            if self.dark.sbAndRefpix is not None:
                copysb = self.dark.sbAndRefpix
            if self.dark.zeroframe is not None:
                copyzero = self.dark.zeroframe
            for i in range(ncopies):
                self.dark.data = np.vstack((self.dark.data, copydata))
                if self.dark.sbAndRefpix is not None:
                    self.dark.sbAndRefpix = np.vstack((self.dark.sbAndRefpix, copysb))
                if self.dark.zeroframe is not None:
                    self.dark.zeroframe = np.vstack((self.dark.zeroframe, copyzero))

        # Partial copy of dark (some integrations)
        if extras > 0:
            self.dark.data = np.vstack((self.dark.data, self.dark.data[0:extras, :, :, :]))
            if self.dark.sbAndRefpix is not None:
                self.dark.sbAndRefpix = np.vstack((self.dark.sbAndRefpix,
                                                   self.dark.sbAndRefpix[0:extras, :, :, :]))
            if self.dark.zeroframe is not None:
                self.dark.zeroframe = np.vstack((self.dark.zeroframe, self.dark.zeroframe[0:extras, :, :]))

        self.dark.header['NINTS'] = req

    def linearize_dark(self, darkobj, save_refs=True):
        """Beginning with the input dark current ramp, run the dq_init, saturation, superbias
        subtraction, refpix and nonlin pipeline steps in order to produce a linearized
        version of the ramp. This will be used when combining the dark ramp with the
        simulated signal ramp.

        Parameters
        -----------
        darkobj : obj
            Instance of read_fits class containing dark current data and info

        save_refs : bool
            Whether or not to save the names of the reference files used in the
            pipeline run.

        Returns
        -------
        linDarkObj : obj
            Modified read_fits instance with linearized dark current data
        """
        from jwst.dq_init import DQInitStep
        from jwst.saturation import SaturationStep
        from jwst.superbias import SuperBiasStep
        from jwst.refpix import RefPixStep
        from jwst.linearity import LinearityStep

        # First we need to place the read_fits object into a RampModel instance
        if self.runStep['linearized_darkfile']:
            subfile = self.params['Reffiles']['linearized_darkfile']
        else:
            subfile = self.params['Reffiles']['dark']
        dark = darkobj.insert_into_datamodel(subfile)

        self.logger.info('Creating a linearized version of the dark current input ramp using JWST calibration pipeline.')

        # Run the DQ_Init step
        if self.runStep['badpixmask']:
            linDark = DQInitStep.call(dark,
                                      config_file=self.params['newRamp']['dq_configfile'],
                                      override_mask=self.params['Reffiles']['badpixmask'])
        else:
            linDark = DQInitStep.call(dark, config_file=self.params['newRamp']['dq_configfile'])

        # If the saturation map is provided, use it. If not, default to whatever is in CRDS
        if self.runStep['saturation_lin_limit']:
            linDark = SaturationStep.call(linDark,
                                          config_file=self.params['newRamp']['sat_configfile'],
                                          override_saturation=self.params['Reffiles']['saturation'])
        else:
            linDark = SaturationStep.call(linDark,
                                          config_file=self.params['newRamp']['sat_configfile'])

        # If the superbias file is provided, use it. If not, default to whatever is in CRDS
        if self.runStep['superbias']:
            linDark = SuperBiasStep.call(linDark,
                                         config_file=self.params['newRamp']['superbias_configfile'],
                                         override_superbias=self.params['Reffiles']['superbias'])
        else:
            linDark = SuperBiasStep.call(linDark,
                                         config_file=self.params['newRamp']['superbias_configfile'])

        # Reference pixel correction
        linDark = RefPixStep.call(linDark,
                                  config_file=self.params['newRamp']['refpix_configfile'])

        # Save a copy of the superbias- and reference pixel-subtracted
        # dark. This will be used later to add these effects back in
        # after the synthetic signals have been added and the non-linearity
        # effects are added back in when using the PROPER combine method.
        sbAndRefpixEffects = dark.data - linDark.data

        # Linearity correction - save the output so that you won't need to
        # re-run the pipeline when using the same dark current file in the
        # future. Use the linearity coefficient file if provided
        base_name = self.params['Output']['file'].split('/')[-1]
        linearoutfile = base_name[0:-5] + '_linearized_dark_current_ramp.fits'
        if self.runStep['linearity']:
            linDark = LinearityStep.call(linDark,
                                         config_file=self.params['newRamp']['linear_configfile'],
                                         override_linearity=self.params['Reffiles']['linearity'])
        else:
            linDark = LinearityStep.call(linDark,
                                         config_file=self.params['newRamp']['linear_configfile'])

        # Now we need to put the data back into a read_fits object
        linDarkobj = read_fits.Read_fits()

        # Save the reference files used during the pipeline run
        if save_refs:
            self.linearity_reffile = linDark.meta.ref_file.linearity.name
            self.mask_reffile = linDark.meta.ref_file.mask.name
            self.saturation_reffile = linDark.meta.ref_file.saturation.name
            self.superbias_reffile = linDark.meta.ref_file.superbias.name

        # Populate the data
        linDarkobj.model = linDark
        linDarkobj.rampmodel_to_obj()
        linDarkobj.sbAndRefpix = sbAndRefpixEffects

        return linDarkobj

    @logging_functions.log_fail
    def prepare(self):
        """MAIN FUNCTION"""
        # Read in the yaml parameter file
        self.read_parameter_file()

        # Get the log caught up on what's already happened
        self.logger.info('\n\nRunning dark_prep..\n')
        self.logger.info('Reading parameter file: {}\n'.format(self.paramfile))
        self.logger.info('Original log file name: ./{}'.format(STANDARD_LOGFILE_NAME))

        # Make filter/pupil values respect the filter/pupil wheel they are in
        self.params['Readout']['filter'], self.params['Readout']['pupil'] = \
            utils.normalize_filters(self.params['Inst']['instrument'], self.params['Readout']['filter'], self.params['Readout']['pupil'])

        # Create dictionary to use when looking in CRDS for reference files
        self.crds_dict = crds_tools.dict_from_yaml(self.params)

        # Expand param entries to full paths where appropriate
        self.params = utils.full_paths(self.params, self.modpath, self.crds_dict, offline=self.offline)
        self.filecheck()

        # Base name for output files
        base_name = self.params['Output']['file'].split('/')[-1]
        self.basename = os.path.join(self.params['Output']['directory'],
                                     base_name[0:-5])

        # Check the entered read pattern info
        self.readpattern_check()

        # Check input parameters for any bad values
        self.check_params()

        # Find out how many groups/integrations we need
        self.get_frame_count_info()

        # If the simulation will have more than one integration, generate
        # a list of dark current files that can be read in to provide
        # more data
        if self.numints > 1:
            if self.runStep['linearized_darkfile']:
                dark_dir = os.path.split(self.params['Reffiles']['linearized_darkfile'])[0]
            else:
                dark_dir = os.path.split(self.params['Reffiles']['dark'])[0]
            search_string = '*.fits'
            if 'FGS1' in self.params['Readout']['array_name']:
                search_string = FGS1_DARK_SEARCH_STRING
            elif 'FGS2' in self.params['Readout']['array_name']:
                search_string = FGS2_DARK_SEARCH_STRING
            dark_list = glob(os.path.join(dark_dir, search_string))
        else:
            if self.runStep['linearized_darkfile']:
                dark_list = [self.params['Reffiles']['linearized_darkfile']]
            else:
                dark_list = [self.params['Reffiles']['dark']]
        dark_list = np.array(dark_list)

        # Read in the subarray definition file
        self.subdict = utils.read_subarray_definition_file(self.params['Reffiles']['subarray_defs'])
        self.params = utils.get_subarray_info(self.params, self.subdict)

        # Get the subarray boundaries from pysiaf
        siaf_inst = self.params['Inst']['instrument']
        instrument_siaf = siaf_interface.get_instance(siaf_inst)
        self.siaf = instrument_siaf[self.params['Readout']['array_name']]
        junk0, junk1, self.ffsize, \
            self.subarray_bounds = siaf_interface.get_siaf_information(instrument_siaf,
                                                                       self.params['Readout']['array_name'],
                                                                       0.0, 0.0,
                                                                       self.params['Telescope']['rotation'])

        # If there will be too many frames, then the file will need
        # to be split into pieces to save memory
        xdim = self.subarray_bounds[2] - self.subarray_bounds[0] + 1
        ydim = self.subarray_bounds[3] - self.subarray_bounds[1] + 1
        if self.file_splitting:
            # Mimic what is done in catalog_seed_image, so that we get the same answer. First,
            # split the file based on frames (i.e. RAPID). After that, split based on groups.
            frames_per_group = self.frames_per_integration / self.numgroups
            split_seed_frames, seg_ind, int_ind = find_file_splits(xdim, ydim, self.frames_per_integration,
                                                                   self.numints, frames_per_group=frames_per_group)



            forced_ints_per_file = int(self.frames_per_integration / self.numgroups) * (int_ind[1] - int_ind[0])
            split_seed, group_segment_indexes, integration_segment_indexes = find_file_splits(xdim, ydim,
                                                                                              self.params['Readout']['ngroup'],
                                                                                              self.params['Readout']['nint'],
                                                                                              force_delta_int=forced_ints_per_file)
            self.logger.info('File splitting enabled. Starting integration numbers for all file segments: {}'.format(integration_segment_indexes))
        else:
            self.logger.info('File-splitting disabled. Output will be forced into a single file.')
            split_seed = False
            group_segment_indexes = np.array([0, self.params['Readout']['ngroup']])
            integration_segment_indexes = np.array([0, self.params['Readout']['nint']])

        # In order to avoid the case of having a single integration
        # in the final file, which leads to rate rather than rateints
        # files in the pipeline, check to be sure that the integration
        # splitting indexes indicate the last split isn't a single
        # integration
        if len(integration_segment_indexes) > 2:
            delta_int = integration_segment_indexes[1:] - integration_segment_indexes[0: -1]
            if delta_int[-1] == 1 and delta_int[0] != 1:
                integration_segment_indexes[-2] -= 1
                self.logger.info('Adjusted final segment file to avoid single integration: {}'.format(integration_segment_indexes))

        # This currently pays attention only to splitting the file on integrations, and not
        # within an integration. It seems like we can keep it this way as the dark data is
        # reordered to the requested readout pattern before saving. It seems unlikely
        # that observations will be taken that require splitting the file within integrations
        # in that case. BUT at one point we do still have all the frames in memory
        # (within data_volume_check). So at that point, a DEEP exposure with 20 groups
        # will have 400 frames open. This is also an unlikely situation, but it could
        # be a problem.
        if len(integration_segment_indexes[:-1]) > 1:
            self.logger.info(('An estimate of processing time remaining will be provided after the first segment '
                               'has been completed.\n\n'))

        self.dark_files = []
        for seg_index, segment in enumerate(integration_segment_indexes[:-1]):
            # Start timer
            self.timer.start()

            # Get the number of integrations being simulated
            if split_seed:
                number_of_ints = integration_segment_indexes[seg_index+1] - segment
            else:
                number_of_ints = self.numints

            self.logger.info('Segment number: {}'.format(seg_index))
            self.logger.info('Number of integrations: {}'.format(number_of_ints))

            # Generate the list of dark current files to use
            if number_of_ints == 1:
                use_all_files = False
                files_to_use = [dark_list[0]]
            elif number_of_ints > 1 and number_of_ints <= len(dark_list):
                use_all_files = False
                try:
                    random_indices = np.random.choice(len(dark_list), number_of_ints, replace=False)
                    files_to_use = dark_list[random_indices]
                except TypeError:
                    self.logger.debug("Error in random choice")
                    self.logger.debug(dark_list)
                    self.logger.debug(number_of_ints)
                    self.logger.debug(random_indices)
                    self.logger.debug(type(random_indices))
                    self.logger.debug(type(random_indices[0]))
                    raise ValueError("Error in random choice.")
            else:
                use_all_files = True
                files_to_use = dark_list

            self.logger.info('Dark files to use:')
            for dummy in files_to_use:
                self.logger.info(dummy)

            # Create a list describing which dark current file goes with
            # which integration. Force the 0th element to be the first
            # file in the list, to be assured that the final_dark variable
            # is defined.
            mapping = np.random.choice(len(files_to_use), size=number_of_ints)
            mapping[0] = 0

            # Loop over dark current files
            for file_index, filename in enumerate(files_to_use):

                self.logger.info('Working on dark file: {}'.format(filename))
                if not use_all_files:
                    frames = np.arange(len(files_to_use))
                else:
                    # In this case, we need to repeat at least some of
                    # the dark ramps because the number of integrations
                    # is larger than the number of dark current files
                    frames = np.where(mapping == file_index)[0]

                # If the random number picker didn't pick this file for
                # any integrations, then skip reading it in.
                if len(frames) == 0:
                    continue

                # Read in the input dark current frame
                if not self.runStep['linearized_darkfile']:
                    self.get_base_dark(filename)
                    self.linDark = None
                else:
                    self.read_linear_dark(filename)
                    self.dark = self.linDark

                # Make sure there is enough data (frames/groups)
                # in the input integration to produce
                # the proposed output integration
                self.data_volume_check(self.dark)

                # Put the input dark (or linearized dark) into the
                # requested readout pattern
                self.dark, sbzeroframe = self.reorder_dark(self.dark)
                self.logger.info(('DARK has been reordered to {} to match the input readpattern of {}'
                                  .format(self.dark.data.shape, self.dark.header['READPATT'])))

                # If a raw dark was read in, create linearized version
                # here using the SSB pipeline. Better to do this
                # on the full dark before cropping, so that reference
                # pixels can be used in the processing.
                if ((self.params['Inst']['use_JWST_pipeline']) & (self.runStep['linearized_darkfile'] is False)):

                    # Linearize the dark ramp via the SSB pipeline.
                    # Also save a diff image of the original dark minus
                    # the superbias and refpix subtracted dark, to use later.

                    # In order to linearize the dark, the JWST pipeline must
                    # be present, and self.dark will have to be translated back
                    # into a RampModel instance
                    # print('Working on {}'.format(self.dark))
                    self.linDark = self.linearize_dark(self.dark)
                    self.logger.info("Linearized dark shape: {}".format(self.linDark.data.shape))

                    if self.params['Readout']['readpatt'].upper() in ['RAPID', 'NISRAPID', 'FGSRAPID']:
                        self.logger.info(("Output is {}, grabbing zero frame from linearized dark"
                                          .format(self.params['Readout']['readpatt'].upper())))
                        self.zeroModel = read_fits.Read_fits()
                        self.zeroModel.data = self.linDark.data[:, 0, :, :]
                        self.zeroModel.sbAndRefpix = self.linDark.sbAndRefpix[:, 0, :, :]
                    elif ((self.params['Readout']['readpatt'].upper() not in ['RAPID', 'NISRAPID', 'FGSRAPID']) &
                          (self.dark.zeroframe is not None)):
                        self.logger.info(("Now we need to linearize the zeroframe because the "
                                          "output readpattern is not RAPID, NISRAPID, or FGSRAPID"))
                        # Now we need to linearize the zeroframe. Place it
                        # into a RampModel instance before running the
                        # pipeline steps
                        self.zeroModel = read_fits.Read_fits()
                        self.zeroModel.data = np.expand_dims(self.dark.zeroframe, axis=1)
                        self.zeroModel.header = self.linDark.header
                        self.zeroModel.header['NGROUPS'] = 1
                        self.zeroModel = self.linearize_dark(self.zeroModel, save_refs=False)
                        # Return the zeroModel data to 3 dimensions
                        # integrations, y, x
                        self.zeroModel.data = self.zeroModel.data[:, 0, :, :]
                        self.zeroModel.sbAndRefpix = self.zeroModel.sbAndRefpix[:, 0, :, :]
                        # In this case the zeroframe has changed from what
                        # was read in. So let's remove the original zeroframe
                        # to avoid confusion
                        self.linDark.zeroframe = np.zeros(self.linDark.zeroframe.shape)
                    else:
                        # In this case, the input dark has no zero frame data.
                        # Let's add an (admittedly crude) approximation. Mirage
                        # will use this to create the zeroframe in the simulated
                        # data. The calibration pipeline does not currently do
                        # anything with the zeroframe. In the future it will be
                        # be used to determine signal rate for pixels that are
                        # saturated in all groups. Our approximation here will be
                        # the signal in the first group of the dark data normalized
                        # to the exposure time of a single frame. This block of the
                        # code should only be run if the input dark is a non-RAPID
                        # dark from ground testing (i.e. a converted FITSWriter file)
                        self.logger.info(('Non-RAPID input dark with no zeroframe extension.'
                                          ' Creating an approximation.'))
                        self.zeroModel = read_fits.Read_fits()
                        nint, ng, ny, nx = self.dark.data.shape
                        frametime = self.linDark.header['TFRAME']
                        grouptime = self.linDark.header['TGROUP']
                        zframe = self.linDark.data[0, 0, :, :] / grouptime * frametime

                        # Make 3D, replicate the zeroframe to cover the needed number
                        # of integrations
                        self.zeroModel.data = np.expand_dims(zframe, axis=0)
                        self.zeroModel.sbAndRefpix = np.expand_dims(self.linDark.sbAndRefpix[0, 0, :, :], axis=0)
                        for i in range(2, nint+1):
                            self.zeroModel.data = np.vstack((self.zeroModel.data, np.expand_dims(zframe, axis=0)))
                            self.zeroModel.sbAndRefpix = np.vstack((self.zeroModel.sbAndRefpix, np.expand_dims(self.linDark.sbAndRefpix[0, 0, :, :], axis=0)))

                        self.zeroModel.header = self.linDark.header
                        self.zeroModel.header['NGROUPS'] = 1

                    # Now crop self.linDark, self.dark, and zeroModel
                    # to requested subarray
                    self.dark = self.crop_dark(self.dark)
                    self.linDark = self.crop_dark(self.linDark)

                    if self.zeroModel is not None:
                        self.zeroModel = self.crop_dark(self.zeroModel)

                elif self.runStep['linearized_darkfile']:
                    # If no pipeline is run
                    self.zeroModel = read_fits.Read_fits()
                    self.zeroModel.data = self.dark.zeroframe
                    self.zeroModel.sbAndRefpix = sbzeroframe

                    # Crop the linearized dark to the requested
                    # subarray size
                    # THIS WILL CROP self.dark AS WELL SINCE
                    # self.linDark IS JUST A REFERENCE IN THE NON
                    # PIPELINE CASE!!
                    self.linDark = self.crop_dark(self.linDark)
                    if self.zeroModel.data is not None:
                        self.zeroModel = self.crop_dark(self.zeroModel)
                else:
                    raise NotImplementedError(("Mode not yet supported! Must use either: use_JWST_pipeline "
                                               "= True and a raw or linearized dark or supply a linearized dark. "
                                               "Cannot yet skip the pipeline and provide a raw dark."))

                if file_index == 0:
                    #print('self.linDark shape is {}. Expecting it to be 4D'.format(self.linDark.data.shape))
                    junk, num_grps, ydim, xdim = self.linDark.data.shape
                    final_dark = np.zeros((number_of_ints, num_grps, ydim, xdim), dtype=np.float32)
                    final_sbandrefpix = np.zeros((number_of_ints, num_grps, ydim, xdim), dtype=np.float32)
                    final_zerodata = np.zeros((number_of_ints, ydim, xdim), dtype=np.float32)
                    final_zero_sbandrefpix = np.zeros((number_of_ints, ydim, xdim), dtype=np.float32)

                if not use_all_files:
                    self.logger.info('Setting integration {} to use file {}\n'.format(file_index, os.path.split(filename)[1]))
                    final_dark[file_index, :, :, :] = self.linDark.data
                    final_sbandrefpix[file_index, :, :, :] = self.linDark.sbAndRefpix
                    final_zerodata[file_index, :, :] = self.zeroModel.data
                    final_zero_sbandrefpix[file_index, :, :] = self.zeroModel.sbAndRefpix
                else:
                    # In this case, we need to repeat at least some of
                    # the dark ramps because the number of integrations
                    # is larger than the number of dark current files
                    #frames = np.where(mapping == file_index)[0]
                    self.logger.info('File number {} will be used for integrations {}'.format(file_index, frames))
                    final_dark[frames, :, :, :] = self.linDark.data

                    #final_sbandrefpix[frames, :, :, :] = self.linDark.sbAndRefpix
                    for frame in frames:
                        final_sbandrefpix[frame, :, :, :] = self.linDark.sbAndRefpix


                    final_zerodata[frames, :, :] = self.zeroModel.data
                    final_zero_sbandrefpix[frames, :, :] = self.zeroModel.sbAndRefpix

            # Save the linearized dark
            # if self.params['Output']['save_intermediates']:
            h0 = fits.PrimaryHDU()
            h1 = fits.ImageHDU(final_dark, name='SCI')
            h2 = fits.ImageHDU(final_sbandrefpix, name='SBANDREFPIX')
            h3 = fits.ImageHDU(final_zerodata, name='ZEROFRAME')
            h4 = fits.ImageHDU(final_zero_sbandrefpix, name='ZEROSBANDREFPIX')

            # Populate basic info in the 0th extension header
            nints, ngroups, yd, xd = final_dark.shape
            h0.header['READPATT'] = self.params['Readout']['readpatt'].upper()
            h0.header['NINTS'] = nints
            h0.header['NGROUPS'] = ngroups
            h0.header['NFRAMES'] = self.params['Readout']['nframe']
            h0.header['NSKIP'] = self.params['Readout']['nskip']
            h0.header['DETECTOR'] = self.detector
            h0.header['INSTRUME'] = self.instrument
            h0.header['SLOWAXIS'] = self.slowaxis
            h0.header['FASTAXIS'] = self.fastaxis
            h0.header['R_LINEAR'] = self.linearity_reffile
            h0.header['R_MASK'] = self.mask_reffile
            h0.header['R_SATURA'] = self.saturation_reffile
            h0.header['R_SUPERB'] = self.superbias_reffile

            # Add some basic Mirage-centric info
            h0.header['MRGEVRSN'] = (mirage.__version__, 'Mirage version used')
            h0.header['YAMLFILE'] = (self.paramfile, 'Mirage input yaml file')

            hl = fits.HDUList([h0, h1, h2, h3, h4])
            if split_seed:
                objname = self.basename + '_seg{}_linear_dark_prep_object.fits'.format(str(seg_index+1).zfill(3))
            else:
                objname = self.basename + '_linear_dark_prep_object.fits'
            objname = os.path.join(self.params['Output']['directory'], objname)
            hl.writeto(objname, overwrite=True)
            self.logger.info(("Linearized dark frame plus superbias and reference"
                              "pixel signals, as well as zeroframe, saved to {}. "
                              "This can be used as input to the observation "
                              "generator.".format(objname)))
            self.dark_files.append(objname)

            # Timing information
            self.timer.stop(name='seg_{}'.format(str(seg_index+1).zfill(4)))

            # If there is more than one segment, provide an estimate of processing time
            self.logger.info('\n\nSegment {} out of {} complete.'.format(seg_index+1, len(integration_segment_indexes[:-1])))
            if len(integration_segment_indexes[:-1]) > 1:
                time_per_segment = self.timer.sum(key_str='seg_') / (seg_index+1)
                estimated_remaining_time = time_per_segment * (len(integration_segment_indexes[:-1]) - (seg_index+1)) * u.second
                time_remaining = np.around(estimated_remaining_time.to(u.minute).value, decimals=2)
                finish_time = datetime.datetime.now() + datetime.timedelta(minutes=time_remaining)
                self.logger.info(('Estimated time remaining in dark_prep: {} minutes. '
                                  'Projected finish time: {}'.format(time_remaining, finish_time)))

        # If only one dark current file is needed, return just the file
        # name rather than a list.
        if len(self.dark_files) == 1:
            self.dark_files = self.dark_files[0]

        # important variables
        # self.linDark
        # self.zeroModel
        # Make a read_fits instance that contains all the important
        # information, for ease in connecting with next step of the
        # simulator
        self.prepDark = read_fits.Read_fits()
        self.prepDark.data = final_dark
        self.prepDark.zeroframe = final_zerodata
        self.prepDark.sbAndRefpix = final_sbandrefpix
        self.prepDark.zero_sbAndRefpix = final_zero_sbandrefpix
        self.prepDark.header = self.linDark.header

        logging_functions.move_logfile_to_standard_location(self.paramfile, STANDARD_LOGFILE_NAME,
                                                            yaml_outdir=self.params['Output']['directory'])

    def read_linear_dark(self, input_file):
        """Read in the linearized version of the dark current ramp
        using the read_fits class"""
        try:
            self.logger.info(('Reading in linearized dark current ramp from {}'
                              .format(input_file)))
            self.linDark = read_fits.Read_fits()
            #self.linDark.file = self.params['Reffiles']['linearized_darkfile']
            self.linDark.file = input_file
            self.linDark.read_astropy()
        except:
            raise IOError('WARNING: Unable to read in linearized dark ramp.')

        # Finally, collect information about the detector, which will be needed for astrometry later
        self.detector = self.linDark.header['DETECTOR']
        self.instrument = self.linDark.header['INSTRUME']
        self.fastaxis = self.linDark.header['FASTAXIS']
        self.slowaxis = self.linDark.header['SLOWAXIS']

        self.linearity_reffile = self.get_reffile_metadata('R_LINEAR')
        self.mask_reffile = self.get_reffile_metadata('R_MASK')
        self.saturation_reffile = self.get_reffile_metadata('R_SATURA')
        self.superbias_reffile = self.get_reffile_metadata('R_SUPERB')

    def read_parameter_file(self):
        """Read in the yaml parameter file"""
        try:
            with open(self.paramfile, 'r') as infile:
                self.params = yaml.safe_load(infile)
        except FileNotFoundError:
            self.logger.error("Unable to open {}".format(self.paramfile))

    def readpattern_check(self):
        """Check the readout pattern that's entered and set the number of averaged frames
        as well as the number of skipped frames per group"""
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
            self.logger.info(('Requested readout pattern {} is valid. '
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

    def ref_check(self, rele):
        """Check for the existence of the input reference file
        Assume first that the file is in the directory tree
        specified by the MIRAGE_DATA environment variable.

        Parameters
        ----------
        rele : tup
            Keywords needed to access a given entry in self.params from the
            input yaml file

        Returns
        -------
        None
        """
        rfile = self.params[rele[0]][rele[1]]
        if rfile.lower() != 'none':
            c1 = os.path.isfile(rfile)
            if not c1:
                raise FileNotFoundError(("WARNING: Unable to locate the {}, {} "
                                         "input file! Not present in {}"
                                         .format(rele[0], rele[1], rfile)))

    def reorder_dark(self, dark):
        """Reorder the input dark ramp using the requested
        readout pattern (nframe, nskip). If the initial
        dark ramp is RAPID, NISRAPID, or FGSRAPID, then save and return
        the 0th frame.

        Parameters
        ----------
        dark : obj
            Instance of read_fits class containing dark current data and info

        Returns
        -------
        dark : obj
            Instance of read_fits class where dark data have been reordered to match
            the requested readout pattern

        sbzero : numpy.ndarray
            Zeroth frame data from the dark current data. This is saved separately
            because averaging for non-RAPID readout patterns will destroy the frame
        """
        if self.params['Reffiles']['linearized_darkfile']:
            datatype = np.float
        else:
            datatype = np.int32

        # Get the info for the dark integration
        darkpatt = dark.header['READPATT']
        dark_nframe = dark.header['NFRAMES']
        mtch = self.readpatterns['name'].data == darkpatt
        dark_nskip = self.readpatterns['nskip'].data[mtch][0]

        nint, ngroup, yd, xd = dark.data.shape
        #outdark = np.zeros((self.params['Readout']['nint'], self.params['Readout']['ngroup'], yd, xd))
        outdark = np.zeros((1, self.params['Readout']['ngroup'], yd, xd))

        if dark.sbAndRefpix is not None:
            #outsb = np.zeros((self.params['Readout']['nint'], self.params['Readout']['ngroup'], yd, xd))
            outsb = np.zeros((1, self.params['Readout']['ngroup'], yd, xd))

        # We can only keep a zero frame around if the input dark
        # is RAPID, NISRAPID, or FGSRAPID. Otherwise that information is lost.
        darkzero = None
        sbzero = None
        rapids = ['RAPID', 'NISRAPID', 'FGSRAPID']
        if ((darkpatt in rapids) & (dark.zeroframe is None)):
            dark.zeroframe = dark.data[:, 0, :, :]
            self.logger.info("Saving 0th frame from data to the zeroframe extension")
            if dark.sbAndRefpix is not None:
                sbzero = dark.sbAndRefpix[:, 0, :, :]
        elif ((darkpatt not in rapids) & (dark.zeroframe is None)):
            self.logger.info("Unable to save the zeroth frame because the input dark current ramp is not RAPID.")
            sbzero = None
        elif ((darkpatt in rapids) & (dark.zeroframe is not None)):
            # In this case we already have the zeroframe
            if dark.sbAndRefpix is not None:
                sbzero = dark.sbAndRefpix[:, 0, :, :]
        elif ((darkpatt not in rapids) & (dark.zeroframe is not None)):
            # Non RAPID dark, zeroframe is present, but
            # we can't get sbAndRefpix for the zeroth frame
            # because the pattern is not RAPID.
            sbzero = None

        # We have already guaranteed that either the readpatterns match
        # or the dark is RAPID, so no need to worry about checking for
        # other cases here.

        if ((darkpatt in rapids) and (self.params['Readout']['readpatt'] not in rapids)):

            framesPerGroup = self.params['Readout']['nframe'] + self.params['Readout']['nskip']
            accumimage = np.zeros_like(outdark[0, 0, :, :], dtype=datatype)

            if dark.sbAndRefpix is not None:
                zeroaccumimage = np.zeros_like(outdark[0, 0, :, :], dtype=np.float)

            # Loop over integrations
            #for integ in range(self.params['Readout']['nint']):
            for integ in range(dark.data.shape[0]):
                frames = np.arange(0, self.params['Readout']['nframe'])

                # Loop over groups
                for i in range(self.params['Readout']['ngroup']):
                    # Average together the appropriate frames,
                    # skip the appropriate frames
                    self.logger.info(("Averaging dark current ramp. Frames {}, to become "
                                      "group {}".format(frames, i)))

                    # If averaging needs to be done
                    if self.params['Readout']['nframe'] > 1:
                        accumimage = np.mean(dark.data[integ, frames, :, :], axis=0)
                        if dark.sbAndRefpix is not None:
                            zeroaccumimage = np.mean(dark.sbAndRefpix[integ, frames, :, :], axis=0)

                        # If no averaging needs to be done
                    else:
                        accumimage = dark.data[integ, frames[0], :, :]
                        if dark.sbAndRefpix is not None:
                            zeroaccumimage = dark.sbAndRefpix[integ, frames[0], :, :]

                    outdark[integ, i, :, :] += accumimage
                    if dark.sbAndRefpix is not None:
                        outsb[integ, i, :, :] += zeroaccumimage

                    # Increment the frame indexes
                    frames = frames + framesPerGroup

        elif (self.params['Readout']['readpatt'] == darkpatt):
            # If the input dark is not RAPID, or if the readout
            # pattern of the input dark and
            # the output ramp match, then no averaging needs to
            # be done
            outdark = dark.data[:, 0:self.params['Readout']['ngroup'], :, :]
            if dark.sbAndRefpix is not None:
                outsb = dark.sbAndRefpix[:, 0:self.params['Readout']['ngroup'], :, :]
        else:
            # This check should already have been done,
            # but just to be sure...
            raise ValueError(("WARNING: dark current readout pattern is {} and requested "
                              "output is {}. Cannot convert between the two."
                              .format(darkpatt, self.params['Readout']['readpatt'])))

        # Now place the reorganized dark into the data object and
        # update the appropriate metadata
        dark.data = outdark
        if dark.sbAndRefpix is not None:
            dark.sbAndRefpix = outsb
        dark.header['READPATT'] = self.params['Readout']['readpatt']
        dark.header['NFRAMES'] = self.params['Readout']['nframe']
        dark.header['NSKIP'] = self.params['Readout']['nskip']
        dark.header['NGROUPS'] = self.params['Readout']['ngroup']
        dark.header['NINTS'] = self.params['Readout']['nint']
        return dark, sbzero

    def add_options(self, parser=None, usage=None):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage, description='Simulate JWST ramp')
        parser.add_argument("paramfile", help=("File describing the input parameters and instrument "
                                               "settings to use. (YAML format)."))
        return parser


if __name__ == '__main__':

    usagestring = 'USAGE: dark_prep.py inputs.yaml'

    indark = DarkPrep()
    parser = indark.add_options(usage=usagestring)
    args = parser.parse_args(namespace=indark)
    indark.prepare()
