#! /usr/bin/env python

"""
Create a final simulated exposure.

This module contains code that will combine a seed image and a
dark current exposure into a final simulated exposure. Cosmic rays,
Poisson noise, and other detector effects are addded. This is the
final step when creating simulated data with Mirage. It can be run
after catalog_Seed_image.py and dark_prep.py

Authors:
--------

    - Bryan Hilbert, Kevin Volk

Use:
----

    This module can be imported as such:

    ::

        from mirage.ramp_generator.obs_generator import Observation
        ob = Observation()
        ob.paramfile = 'my_parameters.yaml'
        ob.create()
"""

import sys
import os
import random
import copy
from math import radians
import datetime
import warnings
import argparse

import yaml
import pkg_resources
import numpy as np
from astropy.io import fits, ascii
from astropy.table import Table
from astropy.time import Time, TimeDelta
import astropy.units as u
import pysiaf

import mirage
from mirage.ramp_generator import unlinearize
from mirage.reference_files import crds_tools
from mirage.utils import read_fits, utils, siaf_interface
from mirage.utils import set_telescope_pointing_separated as stp
from mirage.utils.constants import EXPTYPES, MEAN_GAIN_VALUES


INST_LIST = ['nircam', 'niriss', 'fgs']
MODES = {"nircam": ["imaging", "ts_imaging", "wfss", "ts_grism"],
         "niriss": ["imaging", "ami", "pom", "wfss"],
         "fgs": ["imaging"]}


class Observation():
    def __init__(self, offline=False):
        """Instantiate the Observation class

        Parameters
        ----------
        offline : bool
            If True, the check for the existence of the MIRAGE_DATA
            directory is skipped. This is primarily for Travis testing
        """
        self.linDark = None
        self.seed = None
        self.segmap = None
        self.seedheader = None
        self.seedunits = 'ADU/sec'
        self.offline = offline

        # self.coord_adjust contains the factor by which the
        # nominal output array size needs to be increased
        # (used for WFSS mode), as well as the coordinate
        # offset between the nominal output array coordinates.
        self.coord_adjust = {'x': 1., 'xoffset': 0., 'y': 1., 'yoffset': 0.}

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

    def add_crosstalk(self, exposure):
        """Add crosstalk effects to the input exposure

        Parameters
        ----------
        exposure : numpy.ndarray
            Always 4D

        Returns
        -------
        exposure : numpy.ndarray
            Exposure with crosstalk effects added
        """
        ints, groups, yd, xd = exposure.shape
        if self.params['Readout']['namp'] == 4:
            if self.instrument.upper() == 'NIRCAM':
                xdet = self.detector[3:5].upper()
                if xdet[1] == 'L':
                    xdet = xdet[0] + '5'
            else:
                xdet = self.detector
            xtcoeffs = self.read_crosstalk_file(self.params['Reffiles']['crosstalk'], xdet)
            # Only sources on the detector will create crosstalk.
            # If signalimage is larger than full frame
            # because we are creating a grism image, then extract the
            # pixels corresponding to the actual
            # detector, and only create crosstalk values for those.
            xs = 0
            xe = xd
            ys = 0
            ye = yd

            for integ in range(ints):
                for group in range(groups):
                    xtinput = exposure[integ, group, ys:ye, xs:xe]
                    xtimage = self.crosstalk_image(xtinput, xtcoeffs)

                    # Now add the crosstalk image to the signalimage
                    exposure[integ, group, ys:ye, xs:xe] += xtimage
        else:
            print("Crosstalk calculation requested, but the chosen subarray")
            print("is read out using only 1 amplifier.")
            print("Therefore there will be no crosstalk. Skipping this step.")
        return exposure

    def add_crs_and_noise(self, seed):
        """Given a noiseless seed ramp, add cosmic
        rays and poisson noise

        Paramters:
        ----------
        seed : numpy.ndarray
            Exposure to add CRs and noise to

        Returns
        --------
        sim_exposure : numpy.ndarray
            Exposure with CRs and noise added

        sim_zero : numpy.ndarray
            Zeroth read(s) of exposure
        """
        yd, xd = seed.shape[-2:]
        seeddim = len(seed.shape)

        # Run one integration at a time
        # because each needs its own collection
        # of cosmic rays and poisson noise realization
        if seeddim == 2:
            nint = self.params['Readout']['nint']
            ngroups = self.params['Readout']['ngroup']
        elif seeddim == 4:
            nint = seed.shape[0]
            ngroups = int(seed.shape[1] / (self.params['Readout']['nframe'] + self.params['Readout']['nskip']))

        sim_exposure = np.zeros((nint, ngroups, yd, xd))
        sim_zero = np.zeros((nint, yd, xd))

        for integ in range(nint):
            print("Integration {}:".format(integ))
            if seeddim == 2:
                inseed = seed
            elif seeddim == 4:
                inseed = seed[integ, :, :, :]
            if self.runStep['cosmicray']:
                ramp, rampzero = self.frame_to_ramp(inseed)
            else:
                ramp, rampzero = self.frame_to_ramp_no_cr(inseed)

            sim_exposure[integ, :, :, :] = ramp
            sim_zero[integ, :, :] = rampzero
        return sim_exposure, sim_zero

    def add_detector_effects(self, ramp):
        """Add detector-based effects to input data.
        Currently only crosstalk effects are added.

        Parameters
        ----------
        ramp : numpy.ndarray
            Array containing the exposure

        Returns
        -------
        ramp : numpy.ndarray
            Exposure with effects added
        """
        if self.runStep['crosstalk']:
            ramp = self.add_crosstalk(ramp)
        return ramp

    def add_flatfield_effects(self, ramp):
        """Add flat field effects to the exposure

        Paramters:
        ----------
        ramp : numpy.ndarray
            Array containing exposure

        Returns
        --------
        ramp : numpy.ndarray
            Exposure with flat field applied
        """
        # ILLUMINATION FLAT
        if self.runStep['illuminationflat']:
            illuminationflat, illuminationflatheader = self.read_cal_file(self.params['Reffiles']['illumflat'])
            ramp *= illuminationflat

        # PIXEL FLAT
        if self.runStep['pixelflat']:
            pixelflat, pixelflatheader = self.read_cal_file(self.params['Reffiles']['pixelflat'])
            ramp *= pixelflat
        return ramp

    def add_ipc(self, data):
        """
        Add interpixel capacitance effects to the data. This is done by
        convolving the data with a kernel. The kernel is read in from the
        file specified by self.params['Reffiles']['ipc']. The core of this
        function was copied from the IPC correction step in the JWST
        calibration pipeline.

        Parameters
        ----------
        data : obj
            4d numpy ndarray containing the data to which the
            IPC effects will be added

        Returns
        -------
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
            # If add_ipc has already been called, then the correct
            # IPC kernel already exists, in self.kernel
            kernel = np.copy(self.kernel)
        except:
            # If add_ipc has not been called yet, then read in the
            # kernel from the specified file.
            kernel = fits.getdata(self.params['Reffiles']['ipc'])
            # Invert the kernel if requested, to go from a kernel
            # designed to remove IPC effects to one designed to
            # add IPC effects
            if self.params['Reffiles']['invertIPC']:
                print("Inverting IPC kernel prior to convolving with image")
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

    def add_mirage_info(self):
        """Place Mirage-related information in a FITS hdulist so that it can
        be saved with the output data

        Returns
        -------
        hdulist : astroy.io.fits.HDUList
            HDU List containing Mirage-related info in the primary header
        """
        hdulist = fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU()])
        hdulist[0].header['MRGEVRSN'] = (mirage.__version__, 'Mirage version used')
        hdulist[0].header['YAMLFILE'] = (self.paramfile, 'Mirage input yaml file')
        #hdulist[0].header['GAINFILE'] = (self.params['Reffiles']['gain'], 'Gain file used by Mirage')
        hdulist[0].header['GAIN'] = (self.gain, 'Gain value used by Mirage')
        hdulist[0].header['DISTORTN'] = (self.params['Reffiles']['astrometric'],
                                         'Distortion reffile used by Mirage')
        hdulist[0].header['IPC'] = (self.params['Reffiles']['ipc'], 'IPC kernel used by Mirage')
        hdulist[0].header['PIXARMAP'] = (self.params['Reffiles']['pixelAreaMap'],
                                         'Pixel area map used by Mirage')
        hdulist[0].header['CROSSTLK'] = (self.params['Reffiles']['crosstalk'],
                                         'Crosstalk file used by Mirage')
        hdulist[0].header['FLUX_CAL'] = (self.params['Reffiles']['flux_cal'],
                                         'Flux calibration file used by Mirage')
        hdulist[0].header['FTHRUPUT'] = (self.params['Reffiles']['filter_throughput'],
                                         'Filter throughput file used by Mirage')
        hdulist[0].header['PTSRCCAT'] = (self.params['simSignals']['pointsource'],
                                         'Point source catalog used by Mirage')
        hdulist[0].header['GALAXCAT'] = (self.params['simSignals']['galaxyListFile'],
                                         'Galaxy source catalog used by Mirage')
        hdulist[0].header['EXTNDCAT'] = (self.params['simSignals']['extended'],
                                         'Extended source catalog used by Mirage')
        hdulist[0].header['MTPTSCAT'] = (self.params['simSignals']['movingTargetList'],
                                         'Moving point source catalog used by Mirage')
        hdulist[0].header['MTSERSIC'] = (self.params['simSignals']['movingTargetSersic'],
                                         'Moving Sersic catalog used by Mirage')
        hdulist[0].header['MTEXTEND'] = (self.params['simSignals']['movingTargetExtended'],
                                         'Moving extended target catalog used by Mirage')
        hdulist[0].header['NONSDRAL'] = (self.params['simSignals']['movingTargetToTrack'],
                                         'Non-Sidereal catalog used by Mirage')
        hdulist[0].header['BKGDRATE'] = (self.params['simSignals']['bkgdrate'],
                                         'Background rate used by Mirage')
        hdulist[0].header['TRACKING'] = (self.params['Telescope']['tracking'],
                                         'Telescope tracking type for Mirage')
        hdulist[0].header['POISSON'] = (self.params['simSignals']['poissonseed'],
                                        'Random num generator seed for Poisson noise in Mirage')
        hdulist[0].header['PSFWFE'] = (self.params['simSignals']['psfwfe'],
                                       'WebbPSF Wavefront error used by Mirage')
        hdulist[0].header['PSFWFGRP'] = (self.params['simSignals']['psfwfegroup'],
                                         'WebbPSF wavefront error group used by Mirage')
        hdulist[0].header['CRLIB'] = (self.params['cosmicRay']['library'],
                                      'Cosmic ray library used by Mirage')
        hdulist[0].header['CRSCALE'] = (self.params['cosmicRay']['scale'],
                                        'Cosmic ray rate scaling factor used by Mirage')
        hdulist[0].header['CRSEED'] = (self.params['cosmicRay']['seed'],
                                       'Random number generator seed for cosmic rays in Mirage')
        return hdulist

    def add_pam(self, signalramp):
        """ Apply Pixel Area Map to exposure

        Paramters:
        ----------
        signalramp : numpy.ndarray
            Array containing exposure

        Returns
        --------
        signalramp : numpy.ndarary
            Array after multiplying by the pixel area map
        """
        pixAreaMap = self.simple_get_image(self.params['Reffiles']['pixelAreaMap'])

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

    def add_superbias_and_refpix(self, ramp, sbref):
        """Add superbias and reference pixel-associated
        signal to the input ramp

        Parameters
        ----------
        ramp : numpy.ndarray
            Array containing exposure data

        sbref : numpy.ndarray
            Array containing superbias and reference pixel associated signal

        Returns
        -------
        newramp : numpy.ndarray
            Ramp with superbias and refpix signal added
        """
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

    def add_synthetic_to_dark(self, synthetic, dark, syn_zeroframe=None):
        """Add the synthetic data (now an exposure) to the dark current
        exposure.

        If zeroframe is provided, the function uses that to create the
        dark+synthetic zeroframe that is returned. If not provided, the
        function attempts to use the 0th frame of the input synthetic ramp

        Combine the cube of synthetic signals to the real dark current ramp.
        Be sure to adjust the dark current ramp if nframe/nskip is different
        than the nframe/nskip values that the dark was taken with.

        Only RAPID, NISRAPID, FGSRAPID darks will be re-averaged into different
        readout patterns. But a BRIGHT2 dark can be used to create a
        BRIGHT2 simulated ramp

        Arguments:
        ----------
        synthetic : numpy.ndarray
            simulated signals, 4D array

        dark : numpy.ndarray
            dark current exposure, 4D array

        syn_zeroframe : numpy.ndarray
            zeroframe data associated with simulated data

        Returns
        --------
        synthetic : numpy.ndarray
            4D exposure containing combined simulated + dark data
        zeroframe : numpy.ndarray
            Zeroth read(s) of simulated + dark data
        reorder_sbandref : numpy.ndarray
            superbias and refpix signal from the dark
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
        rapids = ["RAPID", "NISRAPID", "FGSRAPID"]
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
                    print(('Averaging dark current ramp in add_synthetic_to_dark.'
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
            reorder_sbandref = dark.sbAndRefpix
        return synthetic, zeroframe, reorder_sbandref

    def apply_lincoeff(self, data, cof):
        """Linearize the input data
        cof[0] + num*cof[1] + cof[2]*num**2 + cof[3]*num**3 +...

        Parameters
        ----------
        data : numpy.ndarray
            data will be 2d

        cof : numpy.ndarray
            Non-linearity coefficients. cof will be 3d, or 1d

        Returns
        -------
        apply : numpy.ndarray
            Linearized data
        """
        apply = 0.
        if len(cof.shape) == 1:
            for i in range(len(cof)):
                apply += (cof[i] * data**i)
        elif len(cof.shape) == 3:
            for i in range(len(cof[:, 0, 0])):
                apply += (cof[i, :, :] * data**i)
        return apply

    def check_param_val(self, value, typ, vmin, vmax, default):
        """Make sure the input value is a float and between given min and max

        Parameters
        ----------
        value : float
            Value to be checked

        typ : str
            Description of variable contents

        vmin : float
            Minimum acceptible value

        vmax : float
            Maximum acceptible value

        default : float
            Value to set value if it is outside accepible bounds

        Returns
        -------
        value : float
            Acceptible value of value
        """
        try:
            value = float(value)
        except ValueError:
            print("WARNING: {} for {} is not a float.".format(value, typ))

        if ((value >= vmin) & (value <= vmax)):
            return value
        else:
            print(("ERROR: {} for {} is not within reasonable bounds. "
                   "Setting to {}".format(value, typ, default)))
            return default

    def check_params(self):
        """Check that the values of various input parameters are acceptible"""
        # check instrument name
        if self.params['Inst']['instrument'].lower() not in INST_LIST:
            raise ValueError(("WARNING: instrument {} not in the list of "
                              "available instruments: {}"
                              .format(self.params['Inst']['instrument'].lower(), INST_LIST)))

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
        self.readpattern_check()

        # Check that readout patterns of input dark and requested output
        # are compatible
        self.readpattern_compatible()

        # Make sure ngroup and nint are integers
        try:
            self.params['Readout']['ngroup'] = int(self.params['Readout']['ngroup'])
        except:
            raise ValueError("WARNING: Input value of ngroup is not an integer.")

        try:
            self.params['Readout']['nint'] = int(self.params['Readout']['nint'])
        except:
            raise ValueError("WARNING: Input value of nint is not an integer.")

        # If instrument is FGS, then force filter to be 'N/A'
        if self.params['Inst']['instrument'].lower() == 'fgs':
            self.params['Readout']['filter'] = 'NA'
            self.params['Readout']['pupil'] = 'NA'

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
        self.runStep['superbias'] = self.check_run_step(self.params['Reffiles']['superbias'])
        self.runStep['nonlin'] = self.check_run_step(self.params['Reffiles']['linearity'])
        self.runStep['gain'] = self.check_run_step(self.params['Reffiles']['gain'])
        # self.runStep['phot'] = self.check_run_step(self.params['Reffiles']['phot'])
        self.runStep['pixelflat'] = self.check_run_step(self.params['Reffiles']['pixelflat'])
        self.runStep['illuminationflat'] = self.check_run_step(self.params['Reffiles']['illumflat'])
        self.runStep['astrometric'] = self.check_run_step(self.params['Reffiles']['astrometric'])
        self.runStep['ipc'] = self.check_run_step(self.params['Reffiles']['ipc'])
        self.runStep['crosstalk'] = self.check_run_step(self.params['Reffiles']['crosstalk'])
        self.runStep['occult'] = self.check_run_step(self.params['Reffiles']['occult'])
        self.runStep['pointsource'] = self.check_run_step(self.params['simSignals']['pointsource'])
        self.runStep['galaxies'] = self.check_run_step(self.params['simSignals']['galaxyListFile'])
        self.runStep['extendedsource'] = self.check_run_step(self.params['simSignals']['extended'])
        self.runStep['movingTargets'] = self.check_run_step(self.params['simSignals']['movingTargetList'])
        self.runStep['movingTargetsSersic'] = self.check_run_step(self.params['simSignals']['movingTargetSersic'])
        self.runStep['movingTargetsExtended'] = self.check_run_step(self.params['simSignals']['movingTargetExtended'])
        self.runStep['MT_tracking'] = self.check_run_step(self.params['simSignals']['movingTargetToTrack'])
        self.runStep['zodiacal'] = self.check_run_step(self.params['simSignals']['zodiacal'])
        self.runStep['scattered'] = self.check_run_step(self.params['simSignals']['scattered'])
        self.runStep['linearity'] = self.check_run_step(self.params['Reffiles']['linearity'])
        self.runStep['cosmicray'] = self.check_run_step(self.params['cosmicRay']['path'])
        self.runStep['saturation_lin_limit'] = self.check_run_step(self.params['Reffiles']['saturation'])
        self.runStep['fwpw'] = self.check_run_step(self.params['Reffiles']['filtpupilcombo'])
        self.runStep['linearized_darkfile'] = self.check_run_step(self.params['Reffiles']['linearized_darkfile'])
        self.runStep['badpixfile'] = self.check_run_step(self.params['Reffiles']['badpixmask'])
        self.runStep['pixelAreaMap'] = self.check_run_step(self.params['Reffiles']['pixelAreaMap'])

        # NON-LINEARITY
        # Make sure the input accuracy is a float with reasonable bounds
        self.params['nonlin']['accuracy'] = self.check_param_val(self.params['nonlin']['accuracy'],
                                                               'nlin accuracy', 1e-12, 1e-6, 1e-6)
        self.params['nonlin']['maxiter'] = self.check_param_val(self.params['nonlin']['maxiter'],
                                                              'nonlin max iterations', 5, 40, 10)
        self.params['nonlin']['limit'] = self.check_param_val(self.params['nonlin']['limit'],
                                                            'nonlin max value', 30000., 1.e6, 66000.)

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
                self.crfile = os.path.join(self.params['cosmicRay']['path'],
                                           "CRs_MCD1.7_" + self.params['cosmicRay']["library"].upper())
            else:
                self.crfile = None
                raise FileNotFoundError(("Warning: unrecognised cosmic ray library {}"
                                         .format(self.params['cosmicRay']["library"])))

        # Read in distortion and WCS-related data. These will be placed
        # in the header of the output file.
        ap_name = self.params['Readout']['array_name']

        # Convert the input RA and Dec of the pointing position into floats
        # check to see if the inputs are in decimal units or hh:mm:ss strings
        try:
            self.ra = float(self.params['Telescope']['ra'])
            self.dec = float(self.params['Telescope']['dec'])
        except:
            self.ra, self.dec = utils.parse_RA_Dec(self.params['Telescope']['ra'],
                                                   self.params['Telescope']['dec'])

        #if abs(self.dec) > 90. or self.ra < 0. or self.ra > 360. or self.ra is None or self.dec is None:
        if abs(self.dec) > 90. or self.ra is None or self.dec is None:
            raise ValueError("WARNING: bad requested RA and Dec {} {}".format(self.ra, self.dec))

        # Make sure the rotation angle is a float
        try:
            self.params['Telescope']["rotation"] = float(self.params['Telescope']["rotation"])
        except:
            print(("ERROR: bad rotation value {}, setting to zero."
                   .format(self.params['Telescope']["rotation"])))
            self.params['Telescope']["rotation"] = 0.

        # Get SIAF-related information and subarray bounds
        siaf_inst = self.params['Inst']['instrument']
        if siaf_inst.lower() == 'nircam':
            siaf_inst = 'NIRCam'
        instrument_siaf = siaf_interface.get_instance(siaf_inst)
        self.siaf = instrument_siaf[self.params['Readout']['array_name']]
        self.local_roll, self.attitude_matrix, self.ffsize, \
            self.subarray_bounds = siaf_interface.get_siaf_information(instrument_siaf,
                                                                       self.params['Readout']['array_name'],
                                                                       self.ra, self.dec,
                                                                       self.params['Telescope']['rotation'])

        # Check that the various scaling factors are floats and
        # within a reasonable range
        self.params['cosmicRay']['scale'] = self.check_param_val(self.params['cosmicRay']['scale'],
                                                                 'cosmicRay', 0, 100, 1)
        self.params['simSignals']['extendedscale'] = self.check_param_val(self.params['simSignals']
                                                                                     ['extendedscale'],
                                                                          'extendedEmission', 0, 10000, 1)
        self.params['simSignals']['zodiscale'] = self.check_param_val(self.params['simSignals']['zodiscale'],
                                                                      'zodi', 0, 10000, 1)
        self.params['simSignals']['scatteredscale'] = self.check_param_val(self.params['simSignals']
                                                                                      ['scatteredscale'],
                                                                           'scatteredLight', 0, 10000, 1)

        # Make sure the requested output format is an allowed value
        if self.params['Output']['format'] not in ['DMS']:
            raise NotImplementedError(("WARNING: unsupported output format {} requested. "
                                        "Possible options are {}."
                                        .format(self.params['Output']['format'], ['DMS'])))

        # Check the output metadata, including visit and observation
        # numbers, obs_id, etc
        kwchecks = ['program_number', 'visit_number', 'visit_group',
                    'sequence_id', 'activity_id', 'exposure_number',
                    'observation_number', 'obs_id', 'visit_id']
        for quality in kwchecks:
            try:
                self.params['Output'][quality] = str(self.params['Output'][quality])
            except ValueError:
                print(("WARNING: unable to convert {} to string. "
                       "This is required.".format(self.params['Output'][quality])))

        # Get the filter wheel and pupil wheel resolver positions for the
        # filter and pupil to use. This information will be placed in the
        # header of the output file
        if self.instrument.upper() in ['NIRCAM', 'NIRISS']:
            fw_positions = ascii.read(self.params['Reffiles']['filter_wheel_positions'])
            if self.instrument.upper() == 'NIRISS':
                f_match = self.params['Readout']['filter'] == fw_positions['Name']
                p_match = self.params['Readout']['pupil'] == fw_positions['Name']
            elif self.instrument.upper() == 'NIRCAM':
                if '5' in self.detector or 'LONG' in self.detector.upper():
                    channel = 'LW'
                else:
                    channel = 'SW'
                f_match = ((self.params['Readout']['filter'] == fw_positions['Name']) & (channel == fw_positions['Channel']))
                p_match = ((self.params['Readout']['pupil'] == fw_positions['Name']) & (channel == fw_positions['Channel']))
            self.filter_wheel_position = fw_positions['Filter_Resolver_Reading_Wheel_Degrees'][f_match].data[0]
            self.pupil_wheel_position = fw_positions['Pupil_Resolver_Reading_Wheel_Degrees'][p_match].data[0]
        elif self.instrument.upper() == 'FGS':
            self.filter_wheel_position = 999.
            self.pupil_wheel_position = 999.

    def check_run_step(self, filename):
        """Check to see if a filename exists in the parameter file.

        Parameters
        ----------
        filename : str
            Name of file to be checked

        Returns
        -------
        state : bool
            Indicates whether or not filename is 'none'
        """
        if ((len(filename) == 0) or (filename.lower() == 'none')):
            return False
        else:
            return True

    def combine_seeds(self, filenames):
        """If the seed image is broken amongst a list of files, read
        in those files and combine the data into a single seed image.

        Parameters
        ----------
        filenames : list
            List of files to be read in

        Returns
        -------
        seed : numpy.ndarray

        segmap : numpy.ndarray

        header : dict
        """
        print("Reconstructing seed image from multiple files")
        for i, filename in enumerate(filenames):

            print('File: ', filename)

            # Read in the data from one file
            seed_data, seg_data, header_data = self.read_seed(filename)
            ints, groups, ydim, xdim = seed_data.shape

            if i == 0:
                nints = header_data['SEGINT']
                ngroups = header_data['SEGGROUP']

                print('Final seed image shape: ({}, {}, {}, {})'.format(nints, ngroups, ydim, xdim))
                seed = np.zeros((nints, ngroups, ydim, xdim))
                segmap = seg_data

            # Place the data in the reconstructed seed image array based
            # on the integration and group counters
            int_start = header_data['PTINTSRT']
            grp_start = header_data['PTFRMSRT']
            seed[int_start: int_start+ints, grp_start: grp_start+groups, :, :] = seed_data

            #print('data placed into integration: ', int_start, int_start+ints)
            #print('                  groups: ', grp_start, grp_start+groups)
            #print(seed_data.shape)
            #print(seed.shape)

        return seed, segmap, header_data

    def convert_mask(self, inmask, dq_table):
        """Convert a bad pixel mask to contain values used
        by the JWST pipeline

        Parameters
        ----------
        inmask : numpy.ndarray
            Input bad pixel mask

        dq_table : astropy.fits.record
            Table containing data quality definitions

        Returns
        -------
        dqmask : numpy.ndarray
            Data quality array modified to use values in dq_table
        """
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

    def cr_funcs(self, npix, seed=4242):
        """Set up functions that will be used to generate
        cosmic ray hits

        Parameters
        ----------
        npix : int
            Number of pixels across which we are generating a collection of CR hits

        seed : int
            Seed value for random number generator

        Returns
        -------
        crhits : int
            Number of cosmic ray hits in the given number of pixels and time

        crs_perframe L numpy.ndarray
            Array of random values from Poisson distribution for the cosmic ray hits
            per readout frame
        """
        crhits = npix * self.crrate * self.params['cosmicRay']['scale'] * self.frametime
        np.random.seed(seed)

        # Need a set of CRs for all frames, including those
        # that are skipped, in order for the rate of CRs to
        # be consistent.
        crs_perframe = np.random.poisson(crhits, self.params['Readout']['nint'] *
                                         self.params['Readout']['ngroup'] *
                                         (self.params['Readout']['nframe']+self.params['Readout']['nskip']))
        return crhits, crs_perframe

    def create(self):
        """MAIN FUNCTION"""
        print("\nRunning observation generator....\n")

        # Read in the parameter file
        self.read_parameter_file()

        # Create dictionary to use when looking in CRDS for reference files
        self.crds_dict = crds_tools.dict_from_yaml(self.params)

        # Expand param entries to full paths where appropriate
        self.params = utils.full_paths(self.params, self.modpath, self.crds_dict, offline=self.offline)

        self.file_check()

        #print('self.linDark:', self.linDark)
        #print('self.seed:', self.seed)

        # Get the input dark if a filename is supplied
        if self.linDark is None:
            # If no linearized dark is provided, assume the entry in the
            # yaml file is the proper format
            self.linDark = [self.params['Reffiles']['linearized_darkfile']]
            seed_dict = {self.linDark[0]: self.seed}
            print('Reading in dark file from param file: {}'.format(self.linDark))
            self.linear_dark = self.read_dark_file(self.linDark[0])
        elif isinstance(self.linDark, list):
            # Case where dark is split amongst multiple files due to high
            # data volume
            print('Dark file list: {}'.format(self.linDark))
            seed_dict = self.map_seeds_to_dark()
            self.linear_dark = self.read_dark_file(self.linDark[0])
        elif isinstance(self.linDark, str):
            # If a single filename is given, read in the file
            print('Reading in dark file: {}'.format(self.linDark))
            self.linDark = [self.linDark]
            seed_dict = {self.linDark[0]: self.seed}
            self.linear_dark = self.read_dark_file(self.linDark[0])
        else:
            # Case where user has provided a catalogSeed object
            print('Dark object provided')
            self.linear_dark = copy.deepcopy(self.linDark)
            self.linDark = ['none']
            seed_dict = {self.linDark[0]: self.seed}

        # Finally, collect information about the detector,
        # which will be needed for astrometry later
        self.detector = self.linear_dark.header['DETECTOR']
        self.instrument = self.linear_dark.header['INSTRUME']
        self.fastaxis = self.linear_dark.header['FASTAXIS']
        self.slowaxis = self.linear_dark.header['SLOWAXIS']

        # Some basic checks on the inputs to make sure
        # the script won't have to abort due to bad inputs
        # self.check_params()
        self.subdict = utils.read_subarray_definition_file(self.params['Reffiles']['subarray_defs'])
        self.params = utils.get_subarray_info(self.params, self.subdict)

        self.check_params()

        # Read in cosmic ray library files if
        # CRs are to be added to the data later
        if self.runStep['cosmicray']:
            self.read_cr_files()

        # Read in gain map to be used for adding Poisson noise
        # and to scale CRs to be in ADU
        #self.read_gain_map()
        # For the time being, use the mean gain value in constants.py in
        # order to avoid discontinuities that can arise when gain reference
        # files are made using binned areas of pixels, as they are now.
        if self.instrument.lower() == 'niriss':
            self.gain = MEAN_GAIN_VALUES['niriss']
        elif self.instrument.lower() == 'nircam':
            det = copy.deepcopy(self.detector.lower())
            if 'long' in det:
                det = det.replace('long', '5')
            self.gain = MEAN_GAIN_VALUES['nircam'][det]
        elif self.instrument.lower() == 'fgs':
            self.gain = MEAN_GAIN_VALUES['fgs'][self.detector.lower()]

        # Calculate the exposure time of a single frame, based on
        # the size of the subarray
        #tmpy, tmpx = self.seed.shape[-2:]
        tmpy, tmpx = self.linear_dark.data.shape[-2:]
        self.frametime = utils.calc_frame_time(self.instrument, self.params['Readout']['array_name'],
                                               tmpx, tmpy, self.params['Readout']['namp'])
        print("Frametime is {}".format(self.frametime))

        # Calculate the rate of cosmic ray hits expected per frame
        self.get_cr_rate()

        # Read in saturation file
        if self.params['Reffiles']['saturation'] is not None:
            self.read_saturation_file()
        else:
            print('CAUTION: no saturation map provided. Using')
            print('{} for all pixels.'.format(self.params['nonlin']['limit']))
            dy, dx = self.linear_dark.data.shape[2:]
            self.satmap = np.zeros((dy, dx)) + self.params['nonlin']['limit']


        # Read in non-linearity correction coefficients. We need these
        # regardless of whether we are saving the linearized data or going
        # on to make raw data
        nonlincoeffs = self.get_nonlinearity_coeffs()

        # Read in superbias file if present
        self.read_superbias_file()

        for i, linDark in enumerate(self.linDark):
            temp_outdir, basename = os.path.split(self.params['Output']['file'])

            # Get the segment number of the file if present
            seg_location = linDark.find('_seg')
            if seg_location != -1:
                seg_str = linDark[seg_location+4:seg_location+7]
                #print('first segment string: ', seg_str)
            else:
                try:
                    seg_location = seed_dict[linDark][0].find('_seg')
                except AttributeError:
                    seg_location = -1
                if seg_location != -1:
                    seg_str = seed_dict[linDark][0][seg_location+4:seg_location+7]
                else:
                    seg_str = ''
                #print('second segment string: ', seg_str)

            if seg_str != '':
                # Assume standard JWST filename format
                try:
                    print("Creating output file name with segment number.")
                    parts = basename.split('_')
                    basename = '{}_{}_{}-seg{}_{}_{}'.format(parts[0], parts[1], parts[2], seg_str,
                                                             parts[3], parts[4])
                except IndexError:
                    # Non-standard filename format
                    basename = basename.replace('.fits', '-seg{}.fits'.format(seg_str))
            #basename = os.path.join(temp_outdir, basename)

            if i > 0:
                self.linear_dark = self.read_dark_file(self.linDark[i])

            seed_files = seed_dict[linDark]
            if isinstance(seed_files[0], str):
                print('\nSeed files:')
                print(seed_files)
            # Get the corresponding input seed image(s)
            if isinstance(seed_files, str):
                # If a single filename is suppliedself.read_seed(seed_files)
                self.seed_image, self.segmap, self.seedheader = self.read_seed(seed_files)
            elif isinstance(seed_files, list):
                # If the seed image is a list of files (due to high data
                # volume)
                self.seed_image, self.segmap, self.seedheader = self.combine_seeds(seed_files)
            else:
                # self.seed is a catalogSeed object.
                # In this case we assume that self.segmap and
                # self.seedheader have also been provided as the
                # appropriate objects, since they are saved in
                # the same file as the seed image
                self.seed_image = copy.deepcopy(seed_files)

            # If seed image is in units of electrons/sec then divide
            # by the gain to put in ADU/sec
            if 'UNITS' in self.seedheader.keys():
                if self.seedheader['UNITS'] in ["e-/sec", "e-"]:
                    print(("Seed image is in units of {}. Dividing by gain."
                           .format(self.seedheader['units'])))
                    self.seed_image /= self.gain
            else:
                raise ValueError(("'UNITS' keyword not present in header of "
                                  "seed image. Unable to determine whether the "
                                  "seed image is in units of ADU or electrons."))

            # Translate to ramp if necessary,
            # Add poisson noise and cosmic rays
            # Rearrange into requested read pattern
            # All done in one function to save memory
            simexp, simzero = self.add_crs_and_noise(self.seed_image)

            # Multiply flat fields
            simexp = self.add_flatfield_effects(simexp)
            simzero = self.add_flatfield_effects(np.expand_dims(simzero, axis=1))[:, 0, :, :]

            # Mask any reference pixels
            if self.params['Output']['grism_source_image'] is False:
                simexp, simzero = self.mask_refpix(simexp, simzero)

            # Add IPC effects
            # (Dark current ramp already has IPC in it)
            if self.runStep['ipc']:
                simexp = self.add_ipc(simexp)
                simzero = self.add_ipc(np.expand_dims(simzero, axis=1))[:, 0, :, :]

            # Add the simulated source ramp to the dark ramp
            lin_outramp, lin_zeroframe, lin_sbAndRefpix = self.add_synthetic_to_dark(simexp,
                                                                                     self.linear_dark,
                                                                                     syn_zeroframe=simzero)

            # Add other detector effects (Crosstalk/PAM)
            print('Adding crosstalk')
            lin_outramp = self.add_detector_effects(lin_outramp)
            lin_zeroframe = self.add_detector_effects(np.expand_dims(lin_zeroframe, axis=1))[:, 0, :, :]

            # We need to first subtract superbias and refpix signals from the
            # original saturation limits, and then linearize them
            # Refpix signals will vary from group to group, but only by a few
            # ADU. So let's cheat and just use the refpix signals from group 0

            # Create a linearized saturation map
            limits = np.zeros_like(self.satmap) + 1.e6

            if self.linear_dark.sbAndRefpix is not None:
                lin_satmap = unlinearize.nonLinFunc(self.satmap - self.linear_dark.sbAndRefpix[0, 0, :, :],
                                                    nonlincoeffs, limits)
            elif ((self.linear_dark.sbAndRefpix is None) & (self.runStep['superbias'])):
                # If the superbias and reference pixel signal is not available
                # but the superbias reference file is, then just use that.
                lin_satmap = unlinearize.nonLinFunc(self.satmap - self.superbias,
                                                    nonlincoeffs, limits)

            elif ((self.linear_dark.sbAndRefpix is None) & (self.runStep['superbias'] is False)):
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
            self.linear_output = None
            if 'linear' in self.params['Output']['datatype'].lower():
                # Output filename: append 'linear'
                if 'uncal' in basename:
                    linearrampfile = basename.replace('uncal', 'linear')
                else:
                    linearrampfile = basename.replace('.fits', '_linear.fits')

                # Full path of output file
                #linearrampfile = linearrampfile.split('/')[-1]
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
                    self.save_DMS(lin_outramp, lin_zeroframe, linearrampfile, mod='ramp',
                                  err_ext=err, group_dq=groupdq, pixel_dq=pixeldq)
                else:
                    self.save_fits(lin_outramp, lin_zeroframe, linearrampfile, mod='ramp',
                                   err_ext=err, group_dq=groupdq, pixel_dq=pixeldq)

                stp.add_wcs(linearrampfile, roll=self.params['Telescope']['rotation'])
                print("Final linearized exposure saved to:")
                print("{}".format(linearrampfile))
                self.linear_output = linearrampfile

            # If the raw version is requested, we need to unlinearize
            # the ramp
            self.raw_output = None
            if 'raw' in self.params['Output']['datatype'].lower():
                if self.linear_dark.sbAndRefpix is not None:
                    if self.params['Output']['save_intermediates']:
                        #base_name = self.params['Output']['file'].split('/')[-1]
                        ofile = os.path.join(self.params['Output']['directory'],
                                             basename[0:-5] + '_doNonLin_accuracy.fits')
                        savefile = True
                    else:
                        ofile = None
                        savefile = False

                    print('Unlinearizing exposure.')
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
                    print('Adding superbias and reference pixel signals.')
                    raw_outramp = self.add_superbias_and_refpix(raw_outramp, lin_sbAndRefpix)
                    raw_zeroframe = self.add_superbias_and_refpix(raw_zeroframe, self.linear_dark.zero_sbAndRefpix)

                    # Make sure all signals are < 65535
                    raw_outramp[raw_outramp > 65535] = 65535
                    raw_zeroframe[raw_zeroframe > 65535] = 65535

                    # Save the raw ramp
                    #base_name = self.params['Output']['file'].split('/')[-1]
                    rawrampfile = os.path.join(self.params['Output']['directory'], basename)
                    if self.params['Inst']['use_JWST_pipeline']:
                        self.save_DMS(raw_outramp, raw_zeroframe, rawrampfile, mod='1b')
                    else:
                        self.save_fits(raw_outramp, raw_zeroframe, rawrampfile, mod='1b')
                    stp.add_wcs(rawrampfile, roll=self.params['Telescope']['rotation'])
                    print("Final raw exposure saved to: ")
                    print("{}".format(rawrampfile))
                    self.raw_output = rawrampfile
                else:
                    raise ValueError(("WARNING: raw output ramp requested, but the signal associated "
                                      "with the superbias and reference pixels is not present in "
                                      "the dark current data object. Quitting."))

        print("Observation generation complete.")

    def create_group_entry(self, integration, groupnum, endday, endmilli, endsubmilli, endgroup,
                           xd, yd, gap, comp_code, comp_text, barycentric, heliocentric):
        """Add the GROUP extension to the output file
        From an example Mark Kyprianou sent:

        Parameters
        ----------
        integration : int
            Integration number

        group_number : int
            Group number

        endday : int
            Days since Jan 1 2000

        endmilli : integer
            Milliseconds of the day for given time

        endsubmilli : int
            Time since last millisecond?

        endgroup : str
            End group time, e.g. '2016-01-18T02:43:26.061'

        xd : int
            Number_of_columns e.g. 2048

        yd : int
            Number_of_rows e.g. 2048

        gap : int
            Number of gaps in telemetry

        comp_code : int
            Completion code number e.g. 0 (nominal?)

        comp_text : str
            Completion code text e.g. 'COMPLETE'-from howard
                                      'Normal Completion' - from mark

        barycentric : float
            Barycentric end time (mjd) 57405.11165225

        heliocentric : float
            Heliocentric end time (mjd) 57405.1163058

        Returns
        -------
        group : nump.ndarray
            Input values organized into format needed for group entry in
            JWST formatted file
        """
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

    def create_other_extensions(self, data):
        """If the linearized version of the file is to be saved, we
        need to create the error, pixel dq, and group dq extensions

        Parameters
        ----------
        data : numpy.ndarray
            Array containing the exposure data

        Returns
        -------
        err : numpy.ndarray
            Array containing error values

        pixeldq : numpy.ndarray
            Array containing data quality flags
        """
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

        return err, pixeldq

    def crop_to_subarray(self, data):
        """Crop the given (full frame) array to specified subarray

        Parameters
        ----------
        data : numpy.ndarray
            Array containing exposure data

        Returns
        -------
        data : numpy.ndarray
            Array cropped to requested size, using self.subarray_bounds
        """
        return data[self.subarray_bounds[1]:self.subarray_bounds[3] + 1,
                    self.subarray_bounds[0]:self.subarray_bounds[2] + 1]

    def crosstalk_image(self, orig, coeffs):
        """Using Xtalk coefficients, generate an image of the crosstalk signal

        Parameters
        ----------
        orig : numpy.ndarray
            Array to add crosstalk to

        coeffs : numpy.ndarray
            Crosstalk coefficients from the input coefficint file

        Returns
        -------
        xtalk_corr_im : numpy.ndarray
            Input data modified to have crosstalk
        """
        xtalk_corr_im = np.zeros_like(orig)
        subamp_shift = {"0": 1, "1": -1, "2": 1, "3": -1}

        # List of starting columns for all quadrants.
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
        if self.params['Output']['save_intermediates'] is True:
            phdu = fits.PrimaryHDU(xtalk_corr_im)
            base_name = self.params['Output']['file'].split('/')[-1]
            xtalkout = os.path.join(self.params['Output']['directory'], base_name[0:-5] +
                                    '_xtalk_correction_image.fits')
            phdu.writeto(xtalkout, overwrite=True)

        return xtalk_corr_im

    def do_cosmic_rays(self, image, ngroup, iframe, ncr, seedval):
        """Add cosmic rays to input data

        Parameters
        ----------
        image : numpy.ndarray
            2D array containing exposure data to add cosmic rays to

        ngroup : int
            Group number of the image. Only used when writing out summary info

        iframe : int
            Frame number of the image. Only used when writing out summary info

        ncr : int
            Number of cosmic rays to add to the frame

        seedval : int
            Seed to use for random number generator

        Returns
        -------
        image : numpy.ndarray
            Input image with cosmic rays added
        """
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
            image[i1:i2, j1:j2] = image[i1:i2, j1:j2] + crimage[k1:k2, l1:l2] / self.gain  # self.gainim[k1:k2, l1:l2]

            self.cosmicraylist.write("{} {} {} {} {} {} {}\n".format((j2-j1)/2+j1, (i2-i1)/2+i1, ngroup,
                                     iframe, n, m, np.max(crimage[k1:k2, l1:l2])))
        return image

    def do_poisson(self, signalimage, seedval):
        """Add poisson noise to an input image. Input is assumed
        to be in units of ADU, meaning it must be multiplied by
        the gain when calcuating Poisson noise. Then divide by the
        gain in order for the returned image to also be in ADU

        Parameters
        ----------
        signalimage : numpy.ndarray
            2D array of signals in ADU

        seedval : int
            Seed value for the random number generator

        Returns
        -------
        newimage : numpy.ndarray
            signalimage with Poisson noise added
        """
        # Set the seed
        np.random.seed(seedval)

        # Find the appropriate quantum yield value for the filter
        # if self.params['simSignals']['photonyield']:
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
        signalgain = signalimage * self.gain
        highpix = np.where(signalgain == np.nanmax(signalgain))
        if np.nanmin(signalgain) < 0.:
            neg = signalgain < 0.
            negatives = copy.deepcopy(signalgain)
            negatives[neg] = signalgain[neg]
            signalgain[neg] = 0.

        # Add poisson noise
        newimage = np.random.poisson(signalgain, signalgain.shape).astype(np.float)

        if np.nanmin(signalgain) < 0.:
            newimage[neg] = negatives[neg]

        newimage /= self.gain

        # Quantum yield for NIRCam is always 1.0 (so psym1=0)
        # if self.params['simSignals']['photonyield'] and pym1 > 0.000001 and newimage[i, j] > 0:
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

    def file_check(self):
        """
        Make sure the requested input files exist
        For reference files, assume first that they are located in
        the directory tree under the datadir (from the MIRAGE_DATA
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

    def flag_saturation(self, data, sat):
        """Flag saturated pixels in intput data. Return a dq map
        with the appropriate dq value (2) for saturation

        Parameters
        ----------
        data : numpy.ndarray
            Exposure data

        sat : numpy.ndarray
            2D Saturation map

        Returns
        -------
        satmap : numpy.ndarray
            Map containing flagged saturated pixels
        """
        satmap = (data > sat).astype(np.uint32)
        satmap[satmap > 0] = 2
        return satmap

    def frame_to_ramp(self, data):
        """Convert rate image to ramp, add poisson noise
        and cosmic rays

        Parameters
        ----------
        data : numpy.ndarray
            Seed image. Should be a 2d frame or 3d integration.
            If the original seed image is a 4d exposure, call frame_to_ramp
            with one integration at a time.

        Returns
        -------
        outramp : numpy.ndarray
            3d integration with cosmic rays and poisson noise

        zeroframe : numpy.ndarray
            2d zeroth frame
        """
        # Output ramp will be in requested readout pattern!
        ndim = len(data.shape)

        if ndim == 4:
            raise ValueError("Seed image shouldn't be 4D!")
            nintin, ngroupin, yd, xd = data.shape
        elif ndim == 3:
            ngroupin, yd, xd = data.shape
        elif ndim == 2:
            yd, xd = data.shape

        # If a ramp is given, create a -1st frame that is all zeros
        # so that we can create deltaframes for all frames later
        # This should be the case only for TSO observations or
        # moving targets.
        if ndim == 3:
            data = np.vstack((np.zeros((1, yd, xd)), data))

        outramp = np.zeros((self.params['Readout']['ngroup'], yd, xd), dtype=np.float)

        # Set up functions to apply cosmic rays later
        # Need the total number of active pixels in the
        # output array to multiply the CR rate by
        if self.runStep['cosmicray']:
            npix = int(yd * xd + 0.02)

            # Reinitialize the cosmic ray functions for each integration
            crhits, crs_perframe = self.cr_funcs(npix, seed=self.params['cosmicRay']['seed'])

            # open output file to contain the list of cosmic rays
            base_name = self.params['Output']['file'].split('/')[-1]
            crlistout = os.path.join(self.params['Output']['directory'], base_name[0:-5] + '_cosmicrays.list')
            self.open_cr_list_file(crlistout, crhits)

        # Difference between the latest outimage frame and the
        # latest newsignalimage frame. This is important when nframe>1
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
            accumimage = np.zeros((yd, xd))

            # Group 0: the initial nskip frames don't exist,
            # so adjust indexes accordingly
            rstart = 0
            if i == 0:
                rstart = self.params['Readout']['nskip']

            # Loop over frames within each group if necessary
            for j in range(rstart, framesPerGroup):
                # Frame index number in input data
                frameindex = (i * framesPerGroup) + j - self.params['Readout']['nskip']

                # Signal only since previous frame
                if ndim == 3:
                    deltaframe = data[frameindex+1] - data[frameindex]
                elif ndim == 2:
                    deltaframe = data * self.frametime

                # Add poisson noise
                poissonsignal = self.do_poisson(deltaframe, self.params['simSignals']['poissonseed'])

                # Increment poisson seed value so that the next frame doesn't have identical
                # noise
                self.params['simSignals']['poissonseed'] += 1

                # Create the frame by adding the delta signal
                # and poisson noise associated with the delta signal
                # to the previous frame
                framesignal = previoussignal + poissonsignal

                # Add cosmic rays
                if self.runStep['cosmicray']:
                    framesignal = self.do_cosmic_rays(framesignal, i, j,
                                                    crs_perframe[frameindex],
                                                    self.params['cosmicRay']['seed'])
                    # Increment the seed, so that every frame doesn't have identical
                    # cosmic rays
                    self.params['cosmicRay']['seed'] += 1

                # Keep track of the total signal in the ramp,
                # so that we don't neglect signal which comes
                # in during the frames that are skipped.
                previoussignal = copy.deepcopy(framesignal)

                if ((i == 0) & (j == 0)):
                    zeroframe = copy.deepcopy(framesignal)

                # Add the frame to the group signal image
                if j >= self.params['Readout']['nskip']:
                    print('    Averaging frame {} into group {}'.format(frameindex, i))
                    accumimage += framesignal
                elif j < self.params['Readout']['nskip']:
                    print('    Skipping frame {}'.format(frameindex))

            # divide by nframes if > 1
            if self.params['Readout']['nframe'] > 1:
                accumimage /= self.params['Readout']['nframe']
            outramp[i, :, :] = accumimage

        if self.runStep['cosmicray']:
            # Close the cosmic ray list file
            self.cosmicraylist.close()

        return outramp, zeroframe

    def frame_to_ramp_no_cr(self, data):
        """Convert input seed image/ramp to a
        ramp that includes poisson noise. No
        cosmic rays are added

        Parameters
        ----------
        data : numpy.ndarray
            Seed image. Should be a 2d frame or 3d integration.
            If the original seed image is a 4d exposure, call frame_to_ramp
            with one integration at a time.

        Returns
        -------
        outramp : numpy.ndarray
            3d integration with cosmic rays and poisson noise

        zeroframe : numpy.ndarray
            2d zeroth frame
        """
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
            accumimage = np.zeros((yd, xd))

            # Loop over frames within each group if necessary
            # create each frame
            for j in range(framesPerGroup):

                # Frame index number in input data
                frameindex = (i * framesPerGroup) + j

                # Add poisson noise
                if ndim == 3:
                    framesignal = self.do_poisson(data[frameindex+1],
                                                 self.params['simSignals']['poissonseed'])
                elif ndim == 2:
                    framesignal = self.do_poisson(data*frameindex,
                                                 self.params['simSignals']['poissonseed'])

                # Increment poisson seed value so that the next frame doesn't have identical
                # noise
                self.params['simSignals']['poissonseed'] += 1

                if ((i == 0) & (j == 0)):
                    zeroframe = copy.deepcopy(framesignal)

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

    def get_cr_rate(self):
        """Get the base cosmic ray impact probability.

        The following values are based on JWST-STScI-001928, "A library of simulated cosmic ray events impacting
        JWST HgCdTe detectors by Massimo Robberto", Table 1, times the pixel area of 18 microns square = 3.24e-06
        square cm.  Values are in nucleon events per pixel per second.  Corresponding values from the report are
        4.8983 nucleons/cm^2/second, 1.7783 nucleons/cm^2/second, and 3046.83 nucleons/cm^2/second.  The expected
        rates per full frame read (10.73677 seconds) over the whole set of 2048x2048 pixels are 715, 259, and
        444609 events respectively.

        Note that the SUNMIN rate is lower than the SUNMAX rate.  The MIN and MAX labels refer to the solar activity,
        and the galactic cosmic ray contribution at L2 is reduced at solar maximum compared to solar minimum.  The
        FLARE case is for the largest solar flare event on record (see the Robberto report) and corresponds to conditions
        under which JWST would presumably not be operating.
        """
        self.crrate = 0.
        # The previous values were per full frame read and there was a transcription issue in Volk's code.  These
        # have been corrected.  Values are cosmic ray "hit" rates per pixel per second.
        if "SUNMIN" in self.params["cosmicRay"]["library"]:
            self.crrate = 1.587e-05
        if "SUNMAX" in self.params["cosmicRay"]["library"]:
            self.crrate = 5.762e-06
        if "FLARES" in self.params["cosmicRay"]["library"]:
            self.crrate = 0.0098729

        if self.crrate > 0.:
            print("Base cosmic ray probability per pixel per second: {}".format(self.crrate))

    def get_nonlin_coeffs(self, linfile):
        """Read in non-linearity coefficients from given file

        Parameters
        ----------
        linfile : str
            Name of fits file containing linearity coefficients

        Returns
        -------
        nonlin : numpy.ndarray
            Collection of nonlinearity correction coefficients
        """
        nonlin, nonlinheader = self.read_cal_file(linfile)
        # Set NaN coefficients such that no correction will be made
        nans = np.isnan(nonlin[1, :, :])
        numnan = np.sum(nans)
        if numnan > 0:
            print(("The linearity coefficients of {} pixels are NaNs. "
                   "Setting these coefficients such that no linearity "
                   "correction is made.".format(numnan)))

        for i, cof in enumerate(range(nonlin.shape[0])):
            tmp = nonlin[cof, :, :]
            if i == 1:
                tmp[nans] = 1.
            else:
                tmp[nans] = 0.
            nonlin[cof, :, :] = tmp

        # # Crop to appropriate subarray - ALREADY DONE IN read_cal_file
        # if "FULL" not in self.params['Readout']['array_name']:
        #     nonlin = self.crop_to_subarray(nonlin)
        return nonlin

    def get_nonlinearity_coeffs(self):
        """Wrapper around get_nonlin_coeffs. If the file can't
        be opened, or no file is given, the code falls back to some
        average non-linearity coefficients. This would probably be
        bad to use...
        """
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

    def invert_ipc_kernel(self, kern):
        """
        Invert the IPC kernel such that it goes from being used to remove
        IPC effects from data, to being used to add IPC effects to data,
        or vice versa.

        Parameters
        ----------
        kern : obj
            numpy ndarray, either 2D or 4D, containing the kernel

        Returns
        -------
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
            subkernel = kern[:, :, ys:ye, xs:xe]

        dims = subkernel.shape
        # Force subkernel to be 4D to make the function cleaner
        # Dimensions are (kernely, kernelx, detectory, detectorx)
        if len(dims) == 2:
            subkernel = np.expand_dims(subkernel, axis=2)
            subkernel = np.expand_dims(subkernel, axis=3)
        dims = subkernel.shape

        delta = subkernel * 0.
        nyc = dims[0] // 2
        nxc = dims[1] // 2
        delta[nyc, nxc, :, :] = 1.

        a1 = np.fft.fft2(subkernel, axes=(0, 1))
        a2 = np.fft.fft2(delta, axes=(0, 1))
        aout = a2 / a1
        imout = np.fft.ifft2(aout, axes=(0, 1))
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

    def map_seeds_to_dark(self):
        """
        Create a mapping of which seed image filenames belong to each dark
        file. This is needed primarily for cases where the dark and/or
        seeds are split into multiple files due to data volume limitataions

        Returns
        -------
        mapping : dict
            Dictionary that gives the list of seed images associated with
            each dark current file
        """
        self.seed.sort()
        mapping = {}
        for dark in self.linDark:
            seg = dark.find('_seg')
            seg_str = dark[seg+1:seg+7]
            seeds = [name for name in self.seed if seg_str in name]
            mapping[dark] = seeds
        return mapping

    def mask_refpix(self, ramp, zero):
        """Make sure that reference pixels have no signal
        in the simulated source ramp

        Parameters
        ----------
        ramp : numpy.ndarray
            Array containing exposure data

        zero : numpy.ndarray
            Zeroth frame data

        Returns
        -------
        ramp : numpy.ndarray
            Exposure data with reference pixels zeroed out

        zero : numpy.ndarray
            Zeroth frame data with reference pixels zeroed out
        """
        maskimage = np.zeros((self.ffsize, self.ffsize), dtype=np.int)
        maskimage[4:self.ffsize - 4, 4:self.ffsize - 4] = 1.

        # Crop the mask to match the requested output array
        if "FULL" not in self.params['Readout']['array_name']:
            maskimage = self.crop_to_subarray(maskimage)

        ramp *= maskimage
        zero *= maskimage
        return ramp, zero

    def open_cr_list_file(self, filename, hits):
        """Open a new file and print header info. This file
        that will contain the list and positions of inserted
        cosmic rays.

        Parameters
        ----------
        filename : str
            Name of ascii file to contain the summary of added cosmic rays

        hits : int
            Number of cosmic ray hits per frame

        Returns
        -------
        None
        """
        self.cosmicraylist = open(filename, "w")
        self.cosmicraylist.write("# Cosmic ray list (file set %s random seed %d)\n" %
                                 (self.crfile, self.params['cosmicRay']['seed']))
        self.cosmicraylist.write('# Cosmic ray rate per frame: %13.6e (scale factor %f)\n' %
                                 (hits, self.params['cosmicRay']['scale']))
        self.cosmicraylist.write(("Image_x    Image_y    Group   Frame   CR_File_Index   CR_file_frame   "
                                  "Max_CR_Signal\n"))

    def path_check(self, p):
        """
        Check for the existence of the input path.
        Assume first that the path is in relation to
        the directory tree specified by the MIRAGE_DATA
        environment variable

        Parameters
        ----------
        p : tup
            Nested keys that point to a directory in self.params

        Returns
        -------
        Nothing
        """
        pth = self.params[p[0]][p[1]]
        c1 = os.path.exists(pth)
        if not c1:
            raise NotADirectoryError(("WARNING: Unable to find the requested path "
                                      "{}. Not present in directory tree specified by "
                                      "the {} environment variable."
                                      .format(pth, self.env_var)))

    def populate_group_table(self, starttime, grouptime, ramptime, numint, numgroup, ny, nx):
        """Create some reasonable values to fill the GROUP extension table.
        These will not be completely correct because access to other ssb
        scripts and more importantly, databases, is necessary. But they should be
        close.

        Parameters
        ----------
        starttime : astropy.time.Time
            Starting time of exposure

        grouptime : float
            Exposure time of a single group (seconds)

        ramptime : float
            Exposure time of the entire exposure (seconds)

        numint : int
            Number of integrations in data

        numgroup : int
            Number of groups per integration

        ny : int
            Number of pixels in the y dimension

        nx : int
            Number of pixels in the x dimension

        Returns
        -------
        grouptable : numpy.ndarray
            Group extension data for all groups in the exposure
        """
        # Create the table with a first row populated by garbage
        grouptable = self.create_group_entry(999, 999, 0, 0, 0, 'void', 0, 0, 0, 0, 'void', 1., 1.)

        # Quantities that are fixed for all exposures
        compcode = 0
        comptext = 'Normal Completion'
        numgap = 0

        # Ignore warnings as astropy.time.Time will give a warning
        # related to unknown leap seconds if the date is too far in
        # the future.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            baseday = Time('2020-01-01T00:00:00')

        # Integration start times
        rampdelta = TimeDelta(ramptime, format='sec')
        groupdelta = TimeDelta(grouptime, format='sec')
        intstarts = starttime + (np.arange(numint)*rampdelta)

        for integ in range(numint):
            groups = np.arange(1, numgroup+1)
            groupends = intstarts[integ] + (np.arange(1, numgroup+1)*groupdelta)
            endday = (groupends - baseday).jd

            # If the integration has a single group, force endday to be an array
            if isinstance(endday, float):
                endday = np.array([endday])
            enddayint = [np.int(s) for s in endday]

            # Now to get end_milliseconds, we need milliseconds from the beginning
            # of the day
            inday = TimeDelta(endday - enddayint, format='jd')
            endmilli = inday.sec * 1000.

            # Submilliseconds - just use a random number
            endsubmilli = np.random.randint(0, 1000, len(endmilli))

            # Group end time. need to remove : and - and make lowercase t
            groupending = groupends.isot

            # Approximate these as just the group end time in mjd
            barycentric = groupends.mjd
            heliocentric = groupends.mjd

            # For the case of an integration with a single group, force quantities to be
            # arrays so that everything is consistent
            if isinstance(groupending, str):
                groupending = np.array([groupending])
                barycentric = np.array([barycentric])
                heliocentric = np.array([heliocentric])

            for grp, day, milli, submilli, grpstr, bary, helio in zip(groups, endday, endmilli,
                                                                      endsubmilli, groupending,
                                                                      barycentric, heliocentric):
                entry = self.create_group_entry(integ+1, grp, day, milli, submilli, grpstr, nx, ny,
                                                numgap, compcode, comptext, bary, helio)
                grouptable = np.vstack([grouptable, entry])

        # Now remove the top garbage row from the table
        grouptable = grouptable[1:]
        return grouptable

    def read_cal_file(self, filename):
        """Read in the specified calibration fits file. This is for files that contain
        images (e.g. flats, superbias, etc)

        Parameters
        ----------
        filename : str
            Name of file to be opened

        Returns
        -------
        image : numpy.ndarray
            Array data from input file

        header : list
            Information from file header
        """
        try:
            with fits.open(filename) as h:
                image = h[1].data
                header = h[0].header
        except FileNotFoundError:
            print("WARNING: Unable to open {}".format(filename))

        # extract the appropriate subarray if necessary
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

    def read_cr_files(self):
        """Read in the 10 files that comprise the cosmic ray library"""
        self.cosmicrays = []
        self.cosmicraysheader = []
        for i in range(10):
            idx = '_%2.2d_' % (i)
            str1 = idx + self.params['cosmicRay']['suffix'] + '.fits'
            name = self.crfile + str1
            with fits.open(name) as h:
                im = h[1].data
                head = h[0].header
            self.cosmicrays.append(im)
            self.cosmicraysheader.append(head)

    def read_crosstalk_file(self, file, detector):
        """Read in appropriate line from the xtalk coefficients
        file for the given detector and return the coeffs

        Parameters
        ----------
        file : str
            Name of ascii file containing the crosstalk coefficients

        detector : str
            Name of the detector being simulated

        Returns
        -------
        xtcoeffs : list
            Collection of crosstalk coefficients associated with detector
        """
        xtcoeffs = ascii.read(file, header_start=0)

        coeffs = []
        mtch = xtcoeffs['Det'] == detector.upper()
        if np.any(mtch) is False:
            raise ValueError('Detector {} not found in xtalk file {}'.format(detector, file))

        return xtcoeffs[mtch]

    def read_dark_file(self, filename):
        """Read in a prepared dark current exposure

        Parameters
        ----------
        filename : str
            Name of fits file containing dark current data

        Returns
        -------
        obj : read_fits object
            obj contains: obj.data, obj.sbAndRefpix,
                          obj.zeroframe, obj.zero_sbAndRefpix,
                          obj.header
            Values are None for objects that don't exist
        """
        obj = read_fits.Read_fits()
        obj.file = filename
        obj.read_astropy()
        return obj

    def read_gain_map(self):
        """Read in the gain map. This will be used to
        translate signals from e/s to ADU/sec
        """
        if self.runStep['gain']:
            self.gainim, self.gainhead = self.read_cal_file(self.params['Reffiles']['gain'])
            # set any NaN's to 1.0
            bad = ((~np.isfinite(self.gainim)) | (self.gainim == 0))
            self.gainim[bad] = 1.0

            # Pixels that have a gain value of 0
            # will be reset to have values of 1.0
            # zs = gainim == 0
            # gainim[zs] = 1.0

    def read_parameter_file(self):
        """Read in the yaml parameter file (main input to Mirage)."""
        try:
            with open(self.paramfile, 'r') as infile:
                self.params = yaml.safe_load(infile)
        except FileNotFoundError as e:
            print("WARNING: unable to open {}".format(self.paramfile))
        if self.params['Inst']['instrument'].lower() == 'niriss':
            newfilter,newpupil = utils.check_niriss_filter(self.params['Readout']['filter'],self.params['Readout']['pupil'])
            self.params['Readout']['filter'] = newfilter
            self.params['Readout']['pupil'] = newpupil

    def read_saturation_file(self):
        """Read in saturation map from fits file"""
        if self.runStep['saturation_lin_limit']:
            try:
                self.satmap, self.satheader = self.read_cal_file(self.params['Reffiles']['saturation'])
                bad = ~np.isfinite(self.satmap)
                self.satmap[bad] = 1.e6
            except Exception:
                print(('WARNING: unable to open saturation file {}.'
                       .format(self.params['Reffiles']['saturation'])))
                print(("Please provide a valid file, or place 'none' "
                       "in the saturation entry in the parameter file, "))
                print(("in which case the nonlin limit value in the "
                       "parameter file ({}) will be used for all pixels."
                       .format(self.params['nonlin']['limit'])))
        else:
            print(('CAUTION: no saturation map provided. Using '
                   '{} for all pixels.'.format(self.params['nonlin']['limit'])))
            dy, dx = self.linear_dark.data.shape[2:]
            self.satmap = np.zeros((dy, dx)) + self.params['nonlin']['limit']

    def read_seed(self, filename):
        """Read in the fits file containing the seed image/ramp

        Parameters
        ----------
        filename : str
            Fits file containing the seed image

        Returns
        -------
        seed : numpy.ndarray
            Seed image

        segmap : numnpy.ndarray
            Segmentation map

        seedheader : list
            Information from the header of the seed image file
        """
        with fits.open(filename) as h:
            seed = h[1].data
            seedheader = h[0].header
            try:
                segmap = h[2].data
            except:
                segmap = None
        return seed, segmap, seedheader

    def read_superbias_file(self):
        """Read in superbias from fits file"""
        if self.runStep['superbias']:
            try:
                self.superbias, self.superbiasheader = self.read_cal_file(self.params['Reffiles']['superbias'])
            except:
                raise IOError(("WARNING: unable to open superbias file {}. "
                               "Please provide a valid file in the superbias "
                               "entry in the parameter file."
                               .format(self.params['Reffiles']['superbias'])))
        else:
            raise ValueError('CAUTION: no superbias provided. Quitting.')

    def readpattern_check(self):
        """Check the readout pattern that's entered and set number of used frames and
        number of skipped frames per group, from the readout pattern definition file
        """
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
            raise ValueError(("WARNING: the {} readout pattern is not defined in {}."
                              .format(self.params['Readout']['readpatt'],
                                      self.params['Reffiles']['readpattdefs'])))

    def readpattern_compatible(self):
        """Make sure the input dark has a readout pattern
        that is compatible with the requested output
        readout pattern. The dark must be a flavor of RAPID,
        or have a readout pattern that matches the output.
        """
        rapids = ["RAPID", "NISRAPID", "FGSRAPID"]
        darkpatt = self.linear_dark.header['READPATT']
        if ((darkpatt != self.params['Readout']['readpatt']) &
           (darkpatt not in rapids)):
            raise ValueError(("WARNING: Unable to convert input dark with a "
                              "readout pattern of {}, to the requested readout "
                              "pattern of {}. The readout pattern of the dark "
                              "must be RAPID, NISRAPID, FGSRAPID, or match the requested output "
                              "readout pattern.".format(darkpatt, self.params['Readout']['readpatt'])))

    def ref_check(self, rele):
        """
        Check for the existence of the input reference file
        Assume first that the file is in the directory tree
        specified by the MIRAGE_DATA environment variable.

        Parameters
        ----------
        rele : tup
            Nested keys that point to the refrence file of
            interest. These come from the yaml input file

        Reutrns:
        -------
        Nothing
        """
        rfile = self.params[rele[0]][rele[1]]
        if rfile.lower() != 'none':
            rfile = os.path.abspath(rfile)
            c1 = os.path.isfile(rfile)
            if not c1:
                raise FileNotFoundError(("WARNING: Unable to locate the {}, {} "
                                         "input file! Not present in {}"
                                         .format(rele[0], rele[1], rfile)))

    def int_times_table(self, integration_time, date_obs, time_obs, num_ints):
        """Create and populate the INT_TIMES table, which is saved as a
        separate extension in the output data file

        Parameters
        ----------

        integration_time : float
            Exposure time for a single integration, including the reset
            frame, in seconds

        date_obs : str
            Date string of observation ('2020-02-28')

        time_obs : str
            Time string of observation ('12:24:56')

        num_ints : int
            Number of integrations to put in the table

        Returns
        -------
        int_times_tab : astropy.table.Table
            Table of starting, mid, and end times for each integration
        """
        integration_numbers = np.arange(self.params['Readout']['nint'])

        start_time_string = date_obs + 'T' + time_obs
        start_time = Time(start_time_string)

        integ_time_delta = TimeDelta(integration_time * u.second)
        start_times = start_time + (integ_time_delta * integration_numbers)

        integration_time_exclude_reset = TimeDelta((integration_time - self.frametime) * u.second)
        end_times = start_times + integration_time_exclude_reset

        mid_times = start_times + integration_time_exclude_reset / 2.

        # For now, let's keep the BJD (Barycentric?) times identical
        # to the MJD times.
        start_times_bjd = start_times
        mid_times_bjd = mid_times
        end_times_bjd = end_times

        # Create table
        nrows = len(integration_numbers)
        data_list = [(integration_numbers[i] + 1, start_times.mjd[i], mid_times.mjd[i], end_times.mjd[i],
                      start_times_bjd.mjd[i], mid_times_bjd.mjd[i], end_times_bjd.mjd[i]) for i in range(nrows)]

        int_times_tab = np.array(data_list,
                                 dtype=[('integration_number','<i2'),
                                        ('int_start_MJD_UTC','<f8'),
                                        ('int_mid_MJD_UTC', '<f8'),
                                        ('int_end_MJD_UTC','<f8'),
                                        ('int_start_BJD_TDB','<f8'),
                                        ('int_mid_BJD_TDB','<f8'),
                                        ('int_end_BJD_TDB','<f8')])
        return int_times_tab

    def save_DMS(self, ramp, zeroframe, filename, mod='1b', err_ext=None,
                group_dq=None, pixel_dq=None):
        """Save the new, simulated integration in DMS format (i.e. DMS orientation
        rather than raw fitswriter orientation) using JWST data models

        Parameters
        ----------
        ramp : numpy.ndarray
            Array containing the exposure to be saved

        zeroframe : numpy.ndarray
            The zeroth frame(s) for the exposure

        filename : str
            Name of output file

        mod : str
            Format in which to save the data. Can be '1b' or 'ramp'
            '1b' will save the file in JWST Level 1B format. 'ramp'
            will save the data as if it has gone beyond level 1b, and
            contains the error and dq extensions

        err_ext : numpy.ndarray
            Array containing error values to save in error extension, if
            using mod='ramp'

        group_dq : numpy.ndarray
            Array containing group data quality values. Used only if mod='ramp'

        pixel_dq : numpy.ndarray
            Array containing pixel data quality values. Used only if mod='ramp'

        Returns
        -------
        None
        """
        extra_fits_hdulist = self.add_mirage_info()

        if mod == '1b':
            from jwst.datamodels import Level1bModel as DataModel
        elif mod == 'ramp':
            from jwst.datamodels import RampModel as DataModel
        else:
            raise ValueError(("Model type to use for saving output is "
                              "not recognized. Must be either '1b' or 'ramp'."))
        outModel = DataModel(extra_fits_hdulist)

        # make sure the ramp to be saved has the right number of dimensions
        imshape = ramp.shape
        if len(imshape) == 3:
            ramp = np.expand_dims(ramp, axis=0)

        # insert data into model
        outModel.data = ramp

        if mod == 'ramp':
            outModel.err = err_ext
            outModel.groupdq = group_dq
            outModel.pixeldq = pixel_dq

        # if saving the zeroth frame is requested, insert into the model instance
        if zeroframe is not None:
            # if the zeroframe is a 2D image, then add a dimension,
            # as the model expects 3D
            if len(zeroframe.shape) == 2:
                zeroframe = np.expand_dims(zeroframe, 0)
            outModel.zeroframe = zeroframe
        else:
            print("Zeroframe not present. Setting to all zeros")
            numint, numgroup, ys, xs = ramp.shape
            outModel.zeroframe = np.zeros((numint, ys, xs))

        try:
            outModel.meta.exposure.type = EXPTYPES[self.params['Inst']['instrument'].lower()]\
                [self.params['Inst']['mode'].lower()]
        except:
            raise ValueError('EXPTYPE mapping not complete for this!!! FIX ME!')

        # update various header keywords
        dims = outModel.data.shape
        dtor = radians(1.)

        # Ignore warnings as astropy.time.Time will give a warning
        # related to unknown leap seconds if the date is too far in
        # the future.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
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
        outModel.meta.observation.activity_id = self.params['Output']['activity_id']
        outModel.meta.observation.exposure_number = self.params['Output']['exposure_number']

        outModel.meta.program.pi_name = self.params['Output']['PI_Name']
        outModel.meta.program.title = self.params['Output']['title']
        outModel.meta.program.category = self.params['Output']['Proposal_category']
        outModel.meta.program.sub_category = 'UNKNOWN'
        outModel.meta.program.science_category = self.params['Output']['Science_category']
        outModel.meta.program.continuation_id = 0

        outModel.meta.target.catalog_name = 'UNKNOWN'
        outModel.meta.target.ra = self.params['Output']['target_ra']
        outModel.meta.target.dec = self.params['Output']['target_dec']
        outModel.meta.target.proposer_name = self.params['Output']['target_name']
        outModel.meta.coordinates.reference_frame = 'ICRS'

        outModel.meta.wcsinfo.wcsaxes = 2
        outModel.meta.wcsinfo.crval1 = self.ra
        outModel.meta.wcsinfo.crval2 = self.dec
        outModel.meta.wcsinfo.crpix1 = self.siaf.XSciRef
        outModel.meta.wcsinfo.crpix2 = self.siaf.YSciRef
        outModel.meta.wcsinfo.ctype1 = 'RA---TAN'
        outModel.meta.wcsinfo.ctype2 = 'DEC--TAN'
        outModel.meta.wcsinfo.cunit1 = 'deg'
        outModel.meta.wcsinfo.cunit2 = 'deg'
        outModel.meta.wcsinfo.v2_ref = self.siaf.V2Ref
        outModel.meta.wcsinfo.v3_ref = self.siaf.V3Ref
        outModel.meta.wcsinfo.vparity = self.siaf.VIdlParity
        outModel.meta.wcsinfo.v3yangle = self.siaf.V3IdlYAngle
        outModel.meta.wcsinfo.cdelt1 = self.siaf.XSciScale / 3600.
        outModel.meta.wcsinfo.cdelt2 = self.siaf.YSciScale / 3600.
        outModel.meta.wcsinfo.roll_ref = self.local_roll

        # Grism TSO data have the XREF_SCI and YREF_SCI keywords populated.
        # These are used to describe the location of the source on the detector.
        print('\n\nPopulating xref_sci in output file:')
        print(self.seedheader['XREF_SCI'])

        try:
            outModel.meta.wcsinfo.siaf_xref_sci = self.seedheader['XREF_SCI']
            outModel.meta.wcsinfo.siaf_yref_sci = self.seedheader['YREF_SCI']
        except KeyError:
            print('Unable to propagate XREF_SCI, YREF_SCI from seed image to simualted data file.')

        # ra_v1, dec_v1, and pa_v3 are not used by the level 2 pipelines
        # compute pointing of V1 axis
        pointing_ra_v1, pointing_dec_v1 = pysiaf.rotations.pointing(self.attitude_matrix, 0., 0.)
        outModel.meta.pointing.ra_v1 = pointing_ra_v1
        outModel.meta.pointing.dec_v1 = pointing_dec_v1
        outModel.meta.pointing.pa_v3 = self.params['Telescope']['rotation']

        ramptime = self.frametime * (1 + self.params['Readout']['ngroup'] *
                                     (self.params['Readout']['nframe'] + self.params['Readout']['nskip']))
        # Add time for the reset frame....
        rampexptime = self.frametime * (self.params['Readout']['ngroup'] *
                                        (self.params['Readout']['nframe']+self.params['Readout']['nskip']))

        outModel.meta.observation.date = self.params['Output']['date_obs']
        outModel.meta.observation.time = self.params['Output']['time_obs']

        # Create INT_TIMES table, to be saved in INT_TIMES extension
        int_times = self.int_times_table(ramptime, self.params['Output']['date_obs'], self.params['Output']['time_obs'],
                                         outModel.data.shape[0])
        outModel.int_times = int_times

        if self.runStep['fwpw']:
            fwpw = ascii.read(self.params['Reffiles']['filtpupilcombo'])
        else:
            print(("WARNING: Filter wheel element/pupil wheel element combo reffile not specified. "
                   "Proceeding by saving {} in FILTER keyword, and {} in PUPIL keyword".
                   format(self.params['Readout']['filter'], self.params['Readout']['pupil'])))
            fwpw = Table()
            fwpw['filter_wheel'] = self.params['Readout']['filter']
            fwpw['pupil_wheel'] = self.params['Readout']['pupil']

        # get the proper filter wheel and pupil wheel values for the header
        if ((self.params['Inst']['instrument'].lower() == 'nircam') and
           (self.params['Inst']['mode'].lower() not in ['wfss', 'ts_grism'])):
            mtch = fwpw['filter'] == self.params['Readout']['filter'].upper()
            fw = str(fwpw['filter_wheel'].data[mtch][0])
            pw = str(fwpw['pupil_wheel'].data[mtch][0])
        else:
            pw = self.params['Readout']['pupil']
            fw = self.params['Readout']['filter']

        # Get FGS filter/pupil in proper format
        if fw == 'NA':
            fw = 'N/A'
        if pw == 'NA':
            pw = 'N/A'

        # Filter and pupil info
        outModel.meta.instrument.filter = fw
        outModel.meta.instrument.pupil = pw
        if self.instrument.upper() == 'NIRISS':
            outModel.meta.instrument.filter_position = self.filter_wheel_position
            outModel.meta.instrument.pupil_position = self.pupil_wheel_position

        # Specify whether the exposure is part of a TSO observation
        if self.params['Inst']['mode'].lower() not in ['ts_imaging', 'ts_grism']:
            outModel.meta.visit.tsovisit = False
        else:
            outModel.meta.visit.tsovisit = True

        outModel.meta.dither.primary_type = self.params['Output']['primary_dither_type'].upper()
        outModel.meta.dither.position_number = self.params['Output']['primary_dither_position']
        outModel.meta.dither.total_points = self.params['Output']['total_primary_dither_positions']
        outModel.meta.dither.pattern_size = 'DEFAULT'
        outModel.meta.dither.subpixel_type = self.params['Output']['subpix_dither_type']
        outModel.meta.dither.subpixel_number = self.params['Output']['subpix_dither_position']
        outModel.meta.dither.subpixel_total_points = self.params['Output']['total_subpix_dither_positions']
        outModel.meta.dither.xoffset = self.params['Output']['xoffset']
        outModel.meta.dither.yoffset = self.params['Output']['yoffset']

        # pixel coordinates in FITS header start from 1 not from 0
        xc = (self.subarray_bounds[2] + self.subarray_bounds[0])/2.+1.
        yc = (self.subarray_bounds[3] + self.subarray_bounds[1])/2.+1.

        outModel.meta.exposure.readpatt = self.params['Readout']['readpatt']

        # The subarray name needs to come from the "Name" column in the
        # subarray definitions dictionary
        mtch = self.subdict["AperName"] == self.params["Readout"]['array_name']
        outModel.meta.subarray.name = str(self.subdict["Name"].data[mtch][0])

        # subarray_bounds indexed to zero, but values in header should be
        # indexed to 1.
        outModel.meta.subarray.xstart = self.subarray_bounds[0]+1
        outModel.meta.subarray.ystart = self.subarray_bounds[1]+1
        outModel.meta.subarray.xsize = self.subarray_bounds[2]-self.subarray_bounds[0]+1
        outModel.meta.subarray.ysize = self.subarray_bounds[3]-self.subarray_bounds[1]+1

        nlrefpix = max(4-self.subarray_bounds[0], 0)
        nbrefpix = max(4-self.subarray_bounds[1], 0)
        nrrefpix = max(self.subarray_bounds[2]-(self.ffsize-4), 0)
        ntrefpix = max(self.subarray_bounds[3]-(self.ffsize-4), 0)

        outModel.meta.exposure.nframes = self.params['Readout']['nframe']
        outModel.meta.exposure.ngroups = self.params['Readout']['ngroup']
        outModel.meta.exposure.nints = self.params['Readout']['nint']
        outModel.meta.exposure.integration_start = self.seedheader['SEGINTST'] + 1
        outModel.meta.exposure.integration_end = self.seedheader['SEGINTED'] + 1

        outModel.meta.exposure.sample_time = 10
        outModel.meta.exposure.frame_time = self.frametime
        outModel.meta.exposure.group_time = self.frametime * (self.params['Readout']['nframe'] +
                                                              self.params['Readout']['nskip'])
        outModel.meta.exposure.groupgap = self.params['Readout']['nskip']

        outModel.meta.exposure.nresets_at_start = 1
        outModel.meta.exposure.nresets_between_ints = 1
        outModel.meta.exposure.integration_time = rampexptime
        outModel.meta.exposure.exposure_time = rampexptime * self.params['Readout']['nint']
        outModel.meta.model_type = 'RampModel'

        # set the exposure start time
        outModel.meta.exposure.start_time = ct.mjd
        endingTime = ct.mjd + outModel.meta.exposure.exposure_time/3600./24.
        outModel.meta.exposure.end_time = endingTime
        outModel.meta.exposure.mid_time = ct.mjd + outModel.meta.exposure.exposure_time/3600./24./2.
        outModel.meta.exposure.duration = ramptime

        # populate the GROUP extension table
        n_int, n_group, n_y, n_x = outModel.data.shape
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            outModel.group = self.populate_group_table(ct, outModel.meta.exposure.group_time, rampexptime,
                                                       n_int, n_group, n_y, n_x)

        outModel.save(filename)

        # Now we need to adjust the datamodl header keyword
        # If we leave it as Level1bModel, the pipeline doesn't
        # work properly
        if mod == '1b':
            temp = fits.open(filename, mode='update')
            temp[0].header['DATAMODL'] = 'RampModel'
            temp.flush()
        return

    def save_fits(self, ramp, zeroframe, filename, mod='1b', err_ext=None,
                 group_dq=None, pixel_dq=None):
        """Save the new, simulated integration in DMS format (i.e. DMS orientation
        rather than raw fitswriter orientation) using astropy rather than
        JWST data models

        Parameters
        ----------
        ramp : numpy.ndarray
            Array containing the exposure to be saved

        zeroframe : numpy.ndarray
            The zeroth frame(s) for the exposure

        filename : str
            Name of output file

        mod : str
            Format in which to save the data. Can be '1b' or 'ramp'
            '1b' will save the file in JWST Level 1B format. 'ramp'
            will save the data as if it has gone beyond level 1b, and
            contains the error and dq extensions

        err_ext : numpy.ndarray
            Array containing error values to save in error extension, if
            using mod='ramp'

        group_dq : numpy.ndarray
            Array containing group data quality values. Used only if mod='ramp'

        pixel_dq : numpy.ndarray
            Array containing pixel data quality values. Used only if mod='ramp'

        Returns
        -------
        filename : str
            Same as input filename
        """
        # Make sure the ramp to be saved has the right number of dimensions
        imshape = ramp.shape
        if len(imshape) == 3:
            ramp = np.expand_dims(ramp, axis=0)

        if mod == '1b':
            toohigh = ramp > 65535
            ramp[toohigh] = 65535

        # If saving the zeroth frame is requested, insert into the model instance
        if zeroframe is not None:
            # If the zeroframe is a 2D image, then add a dimension
            if len(zeroframe.shape) == 2:
                zeroframe = np.expand_dims(zeroframe, 0)
        else:
            print("Zeroframe not present. Setting to all zeros")
            numint, numgroup, ys, xs = ramp.shape

        # Place the arrays in the correct extensions of the HDUList
        # using int16 below causes problems! anything set to 65535
        # gets reset to -1, which screws up saturation flagging
        # I think the answer is to save as uint16...

        # Create HDU List of Mirage-centric info
        extra_fits_hdulist = self.add_mirage_info()
        extra_header0 = extra_fits_hdulist[0].header

        if mod == 'ramp':
            ex0 = fits.PrimaryHDU(header=extra_header0)
            ex1 = fits.ImageHDU(ramp.astype(np.float32), name='SCI')
            ex2 = fits.ImageHDU(pixel_dq.astype(np.uint32), name='PIXELDQ')
            ex3 = fits.ImageHDU(group_dq.astype(np.uint8), name='GROUPDQ')
            ex4 = fits.ImageHDU(err_ext.astype(np.float32), name='ERR')
            ex5 = fits.ImageHDU(zeroframe.astype(np.float32), name='ZEROFRAME')
            ex6 = fits.BinTableHDU(name='GROUP')
            ex7 = fits.BinTableHDU(name='INT_TIMES')
            outModel = fits.HDUList([ex0, ex1, ex2, ex3, ex4, ex5, ex6, ex7])
            groupextnum = 6

        elif mod == '1b':
            ex0 = fits.PrimaryHDU(header=extra_header0)
            ex1 = fits.ImageHDU(ramp.astype(np.uint16), name='SCI')
            ex2 = fits.ImageHDU(zeroframe.astype(np.uint16), name='ZEROFRAME')
            ex3 = fits.BinTableHDU(name='GROUP')
            ex4 = fits.BinTableHDU(name='INT_TIMES')
            outModel = fits.HDUList([ex0, ex1, ex2, ex3, ex4])
            groupextnum = 3

        try:
            outModel[0].header['EXP_TYPE'] = EXPTYPES[self.params['Inst']['instrument'].lower()]\
                                             [self.params['Inst']['mode'].lower()]
        except:
            raise ValueError('EXPTYPE mapping not complete for this!!! FIX ME!')

        # update various header keywords
        dims = outModel[1].data.shape
        dtor = radians(1.)

        # Ignore warnings as astropy.time.Time will give a warning
        # related to unknown leap seconds if the date is too far in
        # the future.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
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
        outModel[0].header['TARGPROP'] = self.params['Output']['target_name']
        outModel[0].header['TARG_RA'] = self.params['Output']['target_ra']
        outModel[0].header['TARG_DEC'] = self.params['Output']['target_dec']

        outModel[1].header['WCSAXES'] = 2
        outModel[1].header['CRVAL1'] = self.ra
        outModel[1].header['CRVAL2'] = self.dec
        outModel[1].header['CRPIX1'] = self.siaf.XSciRef
        outModel[1].header['CRPIX2'] = self.siaf.YSciRef
        outModel[1].header['CTYPE1'] = 'RA---TAN'
        outModel[1].header['CTYPE2'] = 'DEC--TAN'
        outModel[1].header['CUNIT1'] = 'deg'
        outModel[1].header['CUNIT2'] = 'deg'
        outModel[1].header['V2_REF'] = self.siaf.V2Ref
        outModel[1].header['V3_REF'] = self.siaf.V3Ref
        outModel[1].header['VPARITY'] = self.siaf.VIdlParity
        outModel[1].header['V3I_YANG'] = self.siaf.V3IdlYAngle
        outModel[1].header['CDELT1'] = self.siaf.XSciScale / 3600.
        outModel[1].header['CDELT2'] = self.siaf.YSciScale / 3600.
        outModel[1].header['ROLL_REF'] = self.local_roll

        # ra_v1, dec_v1, and pa_v3 are not used by the level 2 pipelines
        # compute pointing of V1 axis
        pointing_ra_v1, pointing_dec_v1 = pysiaf.rotations.pointing(self.attitude_matrix, 0., 0.)
        outModel[0].header['RA_V1'] = pointing_ra_v1
        outModel[0].header['DEC_V1'] = pointing_dec_v1
        outModel[0].header['PA_V3'] = self.params['Telescope']['rotation']

        ramptime = self.frametime * (1 + self.params['Readout']['ngroup'] *
                                     (self.params['Readout']['nframe'] + self.params['Readout']['nskip']))
        # Add time for the reset frame....
        rampexptime = self.frametime * (self.params['Readout']['ngroup'] * (self.params['Readout']['nframe'] +
                                                                            self.params['Readout']['nskip']))

        # elapsed time from the end and from the start of the supposid ramp, in seconds
        # put the end of the ramp 1 second before the time the file is written
        # these only go in the fake ramp, not in the signal images....
        outModel[0].header['DATE-OBS'] = self.params['Output']['date_obs']
        outModel[0].header['TIME-OBS'] = self.params['Output']['time_obs']

        # Create INT_TIMES table, to be saved in INT_TIMES extension
        int_times = self.int_times_table(ramptime, self.params['Output']['date_obs'], self.params['Output']['time_obs'],
                                         outModel['SCI'].data.shape[0])
        outModel['INT_TIMES'].data = int_times

        if self.runStep['fwpw']:
            fwpw = ascii.read(self.params['Reffiles']['filtpupilcombo'])
        else:
            print(("WARNING: Filter wheel element/pupil wheel element combo reffile not specified. "
                   "Proceeding by saving {} in FILTER keyword, and {} in PUPIL keyword".
                   format(self.params['Readout']['filter'], self.params['Readout']['pupil'])))
            fwpw = Table()
            fwpw['filter_wheel'] = self.params['Readout']['filter']
            fwpw['pupil_wheel'] = self.params['Readout']['pupil']

        # get the proper filter wheel and pupil wheel values for the header
        if self.params['Inst']['mode'].lower() not in ['wfss', 'ts_grism']:
            mtch = fwpw['filter'] == self.params['Readout']['filter'].upper()
            fw = str(fwpw['filter_wheel'].data[mtch][0])
            pw = str(fwpw['pupil_wheel'].data[mtch][0])
        else:
            pw = self.params['Readout']['pupil']
            fw = self.params['Readout']['filter']

        # Get FGS filter/pupil in proper format
        if fw == 'NA':
            fw = 'N/A'
        if pw == 'NA':
            pw = 'N/A'

        # Filter and pupil info
        outModel[0].header['FILTER'] = fw
        outModel[0].header['PUPIL'] = pw
        if self.instrument.upper() == 'NIRISS':
            outModel[0].header['FWCPOS'] = self.filter_wheel_position
            outModel[0].header['PWCPOS'] = self.pupil_wheel_position

        # Specify whether the exposure is part of a TSO observation
        if self.params['Inst']['mode'].lower() not in ['ts_imaging', 'ts_grism']:
            outModel[0].header['TSOVISIT'] = False
        else:
            outModel[0].header['TSOVISIT'] = True

        outModel[0].header['PATTTYPE'] = self.params['Output']['primary_dither_type']
        outModel[0].header['PATT_NUM'] = self.params['Output']['primary_dither_position']
        outModel[0].header['NUMDTHPT'] = self.params['Output']['total_primary_dither_positions']
        outModel[0].header['PATTSIZE'] = 'DEFAULT'
        outModel[0].header['SUBPXTYP'] = self.params['Output']['subpix_dither_type']
        outModel[0].header['SUBPXNUM'] = self.params['Output']['subpix_dither_position']
        outModel[0].header['SUBPXPNS'] = self.params['Output']['total_subpix_dither_positions']
        outModel[0].header['XOFFSET'] = self.params['Output']['xoffset']
        outModel[0].header['YOFFSET'] = self.params['Output']['yoffset']

        # pixel coordinates in FITS header start from 1 not from 0
        xc = (self.subarray_bounds[2]+self.subarray_bounds[0])/2.+1.
        yc = (self.subarray_bounds[3]+self.subarray_bounds[1])/2.+1.

        outModel[0].header['READPATT'] = self.params['Readout']['readpatt']

        # The subarray name needs to come from the "Name" column in the
        # subarray definitions dictionary
        mtch = self.subdict["AperName"] == self.params["Readout"]['array_name']
        outModel[0].header['SUBARRAY'] = str(self.subdict["Name"].data[mtch][0])

        # subarray_bounds indexed to zero, but values in header should be
        # indexed to 1.
        outModel[0].header['SUBSTRT1'] = self.subarray_bounds[0]+1
        outModel[0].header['SUBSTRT2'] = self.subarray_bounds[1]+1
        outModel[0].header['SUBSIZE1'] = self.subarray_bounds[2]-self.subarray_bounds[0]+1
        outModel[0].header['SUBSIZE2'] = self.subarray_bounds[3]-self.subarray_bounds[1]+1

        nlrefpix = max(4-self.subarray_bounds[0], 0)
        nbrefpix = max(4-self.subarray_bounds[1], 0)
        nrrefpix = max(self.subarray_bounds[2]-(self.ffsize-4), 0)
        ntrefpix = max(self.subarray_bounds[3]-(self.ffsize-4), 0)

        outModel[0].header['NFRAMES'] = self.params['Readout']['nframe']
        outModel[0].header['NGROUPS'] = self.params['Readout']['ngroup']
        outModel[0].header['NINTS'] = self.params['Readout']['nint']

        outModel[0].header['TSAMPLE'] = 10
        outModel[0].header['TFRAME'] = self.frametime
        outModel[0].header['TGROUP'] = self.frametime * (self.params['Readout']['nframe'] +
                                                         self.params['Readout']['nskip'])
        outModel[0].header['GROUPGAP'] = self.params['Readout']['nskip']

        outModel[0].header['NRSTSTRT'] = 1
        outModel[0].header['NRESETS'] = 1
        outModel[0].header['EFFINTTM'] = rampexptime
        outModel[0].header['EFFEXPTM'] = rampexptime * self.params['Readout']['nint']

        # set the exposure start time as the current time
        outModel[0].header['EXPSTART'] = ct.mjd
        outModel[0].header['EXPEND'] = ct.mjd + outModel[0].header['EFFEXPTM']/3600./24.
        outModel[0].header['EXPMID'] = ct.mjd + outModel[0].header['EFFEXPTM']/3600./24./2.

        outModel[0].header['DURATION'] = ramptime

        # populate the GROUP extension table
        n_int, n_group, n_y, n_x = outModel[1].data.shape
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            outModel[groupextnum].data = self.populate_group_table(ct, outModel[0].header['TGROUP'], rampexptime,
                                                                   n_int, n_group, n_y, n_x)
        outModel.writeto(filename, overwrite=True)
        return filename

    def simple_get_image(self, name):
        """Read in an array from a fits file and crop using subarray_bounds

        Parameters
        ----------
        name : str
            Name of fits file to be read in

        Returns
        -------
        image : numpy.ndarray
            Array populated with the file contents
        """
        try:
            image, header = fits.getdata(name, header=True)
        except:
            raise FileNotFoundError('WARNING: unable to read in {}'.format(name))

        # assume that the input is 2D, since we are using it to build a signal rate frame
        imageshape = image.shape
        if len(imageshape) != 2:
            self.printfunc("Error: image %s is not two-dimensional" % (name))
            return None, None

        imageshape = image.shape

        try:
            image = image[self.subarray_bounds[1]:self.subarray_bounds[3]+1,
                          self.subarray_bounds[0]:self.subarray_bounds[2]+1]
        except:
            raise ValueError("Unable to crop image from {}".format(name))

        return image

    def add_options(self, parser=None, usage=None):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage,
                                             description='Simulate JWST ramp')
        parser.add_argument("paramfile", help=('File describing the input parameters and instrument '
                                               'settings to use. (YAML format).'))
        parser.add_argument("linDark", help='File containing linearized dark ramp.')
        parser.add_argument("seed", help='File containing seed image and segmentation map')
        return parser


if __name__ == '__main__':

    usagestring = ('USAGE: obs_generator.py inputs.yaml '
                   'lindark.fits seedimg.fits')

    obs = Observation()
    parser = obs.add_options(usage=usagestring)
    args = parser.parse_args(namespace=obs)
    obs.create()
