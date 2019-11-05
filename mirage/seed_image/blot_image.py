#! /usr/bin/env python

'''
Blot an image back to some input WCSs.

Use in conjunction with crop_mosaic.py
'''
import sys
import os
import glob
from copy import copy
import argparse

import numpy as np
import gwcs
from astropy.io import fits
from jwst import datamodels
from jwst.outlier_detection import outlier_detection
from jwst.assign_wcs import AssignWcsStep
from jwst.datamodels import container
import pysiaf

from ..utils import set_telescope_pointing_separated as stp
from ..utils.siaf_interface import get_instance


class Blot():
    def __init__(self, instrument=None, aperture=None, ra=None, dec=None,
                 pav3=None, output_file=None, blotfile=None, distortion_file=None):
        """Blot (i.e. resample) the given image to be centered on the input
        list of RA, Dec, pav3 values using the appropriate instrument/aperture
        distortion.

        Parameters
        ----------
        instrument : str
            JWST instrument name

        aperture : str
            Aperture describing the subarray to use. (e.g. 'NRCA1_FULL')

        ra : float or list
            List of RA values in decimal degrees, describing the locations
            of the centers of the blotted images

        dec : float or list
            List of Declination values in decimal degress, describing the
            locations of the centers of the blotted images

        pav3 : float or list
            List of PAV3 values in decimal degrees, describing the locations
            of the centers of the blotted images

        output_file : str
            Not currently used

        distortion_file : str
            Name of a CRDS distortion reference file to be used when running
            the assign_wcs pipeline step on the blotted data
        """
        if isinstance(ra, float):
            self.center_ra = [ra]
        else:
            self.center_ra = ra
        if isinstance(dec, float):
            self.center_dec = [dec]
        else:
            self.center_dec = dec
        if isinstance(pav3, float):
            self.pav3 = [pav3]
        else:
            self.pav3 = pav3

        self.instrument = instrument
        self.aperture = aperture
        self.detector = ['']
        self.outfile = output_file
        self.blotfile = blotfile
        self.distortion_file = distortion_file

    def blot(self):
        """MAIN FUNCTION
        Resamples the input image onto the WCS associated with the given
        instrument/apertures
        """
        # Get detector names from the aperture list
        if isinstance(self.aperture, str):
            self.aperture = [self.aperture]

        self.detector = [element.split('_')[0] for element in self.aperture]

        # Make sure detector, ra, dec, and roll have same number
        # of elements
        if ((len(self.center_dec) != len(self.center_ra)) | \
                (len(self.center_dec) != len(self.pav3))):
            raise ValueError(('WARNING: center_ra, center_dec '
                              'and pav3 all must have the same number '
                              'of elements.'))

        if type(self.blotfile) == str:
            input_mod = datamodels.ImageModel(self.blotfile)
            outbase = self.blotfile

        elif type(self.blotfile) == datamodels.image.ImageModel:
            # Not a great solution at the moment. Save
            # imagemodel instance so that you have a fits
            # file from which the header can be retrieved
            input_mod = copy(self.blotfile)
            self.blotfile.save('temp.fits')
            self.blotfile = 'temp.fits'
            outbase = 'mosaic'

        else:
            raise ValueError('WARNING: unrecognized type for blotfile')

        # Create a GWCS object from the input file's header
        input_header = fits.getheader(self.blotfile, ext=1)
        transform = gwcs.utils.make_fitswcs_transform(input_header)
        input_mod.meta.wcs = gwcs.WCS(forward_transform=transform, output_frame='world')

        # Filter and pupil information
        filtername = input_mod.meta.instrument.filter
        try:
            pupil = input_mod.meta.instrument.pupil
        except:
            pupil = 'CLEAR'

        # Get position angle of input data
        input_pav3 = input_mod.meta.wcsinfo.roll_ref

        # Name of temporary file output for set_telescope_pointing
        # to work on
        shellname = 'temp_wcs_container.fits'

        blist = []
        for (aper, det, ra, dec, roll) in \
            zip(self.aperture, self.detector, self.center_ra, self.center_dec, self.pav3):
            # Get aperture-specific info
            self.siaf = get_instance(self.instrument)[aper]

            # Create datamodel with appropriate metadata
            bmodel = self.make_model(det, ra, dec, input_pav3, filtername, pupil)
            bmodel.save(shellname, overwrite=True)

            # Use set_telescope_pointing to compute local roll
            # angle and PC matrix
            stp.add_wcs(shellname, roll=roll)

            bmodel = datamodels.open(shellname)

            # Now we need to run assign_wcs step so that these
            # files get a gwcs object attached to them.
            if self.distortion_file is not None:
                bmodel = AssignWcsStep.call(bmodel, override_distortion=self.distortion_file)
            else:
                bmodel = AssignWcsStep.call(bmodel)

            # Add to the list of data model instances to blot to
            blist.append(bmodel)

        # Place the model instances to blot to in a ModelContainer
        blot_list = container.ModelContainer(blist)

        # Blot the image to each of the WCSs in the blot_list
        pars = {'sinscl': 1.0, 'interp': 'poly5'}
        reffiles = {}
        blotter = outlier_detection.OutlierDetection(blot_list, reffiles=reffiles, **pars)
        blotter.input_models = blot_list
        self.blotted_datamodels = blotter.blot_median(input_mod)

    def get_siaf_info(self, instrument_name, aperture_name):
        """Get v2,v3 reference values and y3yangle for a given aperture

        Parameters
        ----------
        instrument_name : str
            JWST instrument name

        aperture_name : str
            Aperture name to use in pysiaf (e.g. 'NRCA1_FULL')
        """
        self.siaf = pysiaf.Siaf(instrument_name)[aperture_name]

    def make_model(self, detector_name, ra_val, dec_val, roll_val,
                   filter_element, pupil_element):
        """Define the WCS of the blotted output

        Parameters
        ----------
        detector_name : str
            Detector name

        ra_val : float
            RA value in decimal degrees

        dec_val : float
            Declination value in decimal degrees

        roll_val : float
            Local roll angle for the aperture

        filter_element : str
            Filter name

        pupil_element : str
            Pupil name

        Returns
        -------
        blot_to : jwst.datamodels.ImageModel
            Empty ImageModel instance with WCS values populated
        """
        blot_to = datamodels.ImageModel((2048,2048))
        blot_to.meta.wcsinfo.cdelt1 = self.siaf.XSciScale * 3600.
        blot_to.meta.wcsinfo.cdelt2 = self.siaf.YSciScale * 3600.
        blot_to.meta.wcsinfo.crpix1 = 1024.5
        blot_to.meta.wcsinfo.crpix2 = 1024.5
        blot_to.meta.wcsinfo.crval1 = ra_val
        blot_to.meta.wcsinfo.crval2 = dec_val
        blot_to.meta.wcsinfo.ctype1 = 'RA---TAN'
        blot_to.meta.wcsinfo.ctype2 = 'DEC--TAN'
        blot_to.meta.wcsinfo.cunit1 = 'deg'
        blot_to.meta.wcsinfo.cunit2 = 'deg'
        blot_to.meta.wcsinfo.ra_ref = ra_val
        blot_to.meta.wcsinfo.dec_ref = dec_val
        blot_to.meta.wcsinfo.roll_ref = roll_val
        blot_to.meta.wcsinfo.wcsaxes = 2
        blot_to.meta.wcsinfo.v2_ref = self.siaf.V2Ref
        blot_to.meta.wcsinfo.v3_ref = self.siaf.V3Ref
        blot_to.meta.wcsinfo.v3yangle = self.siaf.V3SciYAngle
        blot_to.meta.wcsinfo.vparity = self.siaf.VIdlParity
        blot_to.meta.instrument.channel = 'SHORT'
        if '5' not in detector_name:
            blot_to.meta.instrument.detector = detector_name
        else:
            blot_to.meta.instrument.detector = detector_name.replace('5', 'LONG')

        blot_to.meta.instrument.filter = filter_element
        blot_to.meta.instrument.module = detector_name[0]
        blot_to.meta.instrument.name = 'NIRCAM'
        blot_to.meta.instrument.pupil = pupil_element
        blot_to.meta.telescope = 'JWST'
        blot_to.meta.exposure.start_time = 57410.24546415885
        blot_to.meta.exposure.end_time = 57410.2477009838
        blot_to.meta.exposure.type = 'NRC_IMAGE'
        blot_to.meta.target.ra = ra_val
        blot_to.meta.target.dec = dec_val
        return blot_to

    def add_options(self,parser=None,usage=None):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage,description='Extract SCA-sized area from moasic')
        parser.add_argument("--instrument", help="JWST instrument name to which image will be resampled")
        parser.add_argument("--aperture", help="Instrument aperture to which image will be resampled")
        parser.add_argument("--blotfile", help="Filename or model instance name of fits file containing mosaic.")
        parser.add_argument("--ra", help="RA at the center of the resampled area", type=np.float,nargs='*')
        parser.add_argument("--dec", help="Dec at the center of the resampled area", type=np.float,nargs='*')
        parser.add_argument("--pav3", help="Position angle for outputs to use when blotting", nargs='*')
        parser.add_argument("--output_file", help="Name of output fits file containing extracted image", default=None,nargs='*')
        return parser


if __name__ == '__main__':
    usagestring = 'USAGE: blot_image.py filename.fits --detector A1 A1 --center_ra 10.2 10.2001 --center_dec 12.9 12.91 --local_roll 0 0'

    b = Blot()
    parser = b.add_options(usage = usagestring)
    args = parser.parse_args(namespace=b)
    b.blot()
