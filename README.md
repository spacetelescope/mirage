This repository contains code that can be used to generate
simulated NIRCam data.

The code needed to create simulated data is split into
four general parts:

1. Generate a "seed image"
2 (OPTIONAL) to create WFSS data, disperse the seed image
3. Prepare an existing dark current ramp
4. Combine the seed image with the dark 

More details on these steps are given below:

1. Generate a "seed image"

This portion of the code generates a "seed image" from
input source catalogs. This seed image will
then be used as input to later steps in the NIRCam Data
Simulator (not included in this repo yet).

The seed image is a noiseless countrate image containing
all of the sources specified in the input catalogs. Sources
are scaled to the requested magnitudes and placed in the
appropriate location given the NIRCam distortion model.

The main input to the code is a yaml file (examples are
in the 'inputs' subdirectory). This file contains more
information than is actually needed to run the code,
but the information is all needed to run the entire
NIRCam Data Simulator, so for the moment we've kept all
of the inputs.

Input fields used here include:

Inst:

  mode: examples - 'imaging', 'WFSS', 'moving_target'

Readout:

  array_name e.g. NRCB5_FULL
  
  filter

  pupil

Reffiles:

  subarray_defs - file located in 'config' subdirectory

  astrometric - CRDS-formatted distortion reference file

  distortion_coeffs - SIAF

  flux_cal - file with zeropoints located in 'config' subdirectory

simSignals:

  pointsource - ptsrc catalog, example in 'catalogs' subdirectory

  psfpath - path to PSF library

  psfbasename - base name of files in library

  psfpixfrac - subpixel resolution of PSF library files

  psfwfe - wavefront error to use

  psfwfegroup - realization number for a given WFE

  galaxyListFile - galaxy catalog, example in 'catalogs' subdirectory

  extended - catalog of extended sources

  extendedscale - multiplicative factor to scale extended target brightness (only if the magnitude of the extended target is set to None in the extended catalog)

  extendedCenter - obsolete. No longer used.

  PSFConvolveExtended - True/False, convolve extended objects with NIRCam PSF

  movingTargetList - catalog of targets to trail across FOV (e.g KBOs while observing a sidereal target)

  movingTargetSersic - catalog of "galaxies" (2D Sersic profiles) to trail across FOV

  movingTargetExtended - catalog of extended sources to trail across FOV

  movingTargetToTrack - catalog of a non-sidereal target to track. In this case, all targets in the pointsource, galaxyListFile, and extended catalogs will be trailed across the FOV. To use this, set the 'mode' keyword at the top to 'moving_target'
  
  bkgdrate - Uniform background signal to add, in units of ADU/sec

Telescope:

  ra: RA at the reference location of the detector/aperture being simulated

  dec: Dec at the reference location of the detector/aperture being simulated

  rotation: PAV3 of the telescope (degrees E of N)

Output:

  file: filename used as a base for the seed image output

  directory: directory in which to place the output

  save_intermediates: True/False, save intermediate data products

  grism_source_image: True/False. If true, the length and width of the seed image are increased to larger than the detector FOV



To use the code:

1) From the command line:

python catalog_seed_image.py myfile.yaml

2) Within python:

from nircam_simulator.nircam_simulator.scripts import catalog_seed_image as csi

cat = csi.Catalog_seed()

cat.paramfile = 'myfile.yaml'

cat.run()


Outputs:

Multi-extension fits file with name ending in 'seed_image.fits', containing:

Extension 0: empty

Extension 1: seed image

Extension 2: segmentation map

Also, the seed image, segmentation map, and exposure info dictionary are available as:

self.seedimage, self.seed_segmap, and self.seedinfo



3. Prepare an existing dark current ramp

The input dark current exposure will be reorganized into the
requested readout pattern (if possible). If the input is not
a linearized exposure, then it will be run through the
initial stages of the JWST calibration pipeline in order to
linearize the data. This includes superbias subtraction and
reference pixel subtraction, followed by the linearization
step.

The signal associated with the superbias and reference pixels
is saved along side the linearized dark ramp such that it
can be added back in later, if the user requests a raw output
ramp from the NIRCam Data Simulator.

Dependencies:

If the:

Inst:
  use_JWST_pipeline

input is set to true, then the JWST calibration pipeline is needed.


Output:

The linearized dark current and zeroth frame as saved to a fits file
that uses the name from the Output:file entry in the input yaml file
and ending with '_linearizedDark.fits'.

These are also available as self.linDark and self.zeroModel


To use:

python dark_prep.py myinputs.yaml

or:

from nircam_simulator.nircam_simluator.scripts import dark_prep
dark = dark_prep.DarkPrep()
dark.paramfile = 'myinputs.yaml'
dark.run()


4. Combine the seed image with the dark

This step takes as input a seed image (with segmentation
map), which is a countrate image generated by (for example)
catalog_seed_image. It also takes a linearized dark current
exposure created by dark_prep.

This step converts the seed image into a 4d signal ramp,
adds poisson noise, cosmic rays, and other detector effects.

The ramp is then reorganized into the requested readout
pattern and added to the dark current ramp.

To use:

python obs_generator.py myinputs.yaml

or:

from nircam_simulator.nircam_simulator.scripts import obs_generator
obs = obs_generator.Observation()
obs.linDark = 'V42424024002P000000000112o_B5_F250M_uncal_linear_dark_prep_object.fits'
obs.seed = 'V42424024002P000000000112o_B5_F250M_uncal_F250M_seed_image.fits'
obs.paramfile = 'myinputs.yaml'
obs.create()


To create a simulated exposure, string together all of the steps:

from nircam_simulator.nircam_simulator.scripts import catalog_seed_image
from nircam_simulator.nircam_simulator.scripts import dark_prep
from nircam_simulator.nircam_simulator.scripts import obs_generator

yamlfile = 'seed_catalog_test.yaml'

cat = catalog_seed_image.Catalog_seed()
cat.paramfile = yamlfile
seedimage, segmap, seedinfo = cat.make_seed()

d = dark_prep.DarkPrep()
d.paramfile = yamlfile
d.prepare()

obs = obs_generator.Observation()
obs.linDark = d.prepDark
obs.seed = cat.seedimage 
obs.segmap = cat.seed_segmap 
obs.seedheader = cat.seedinfo 
#obs.seed = 'V42424024002P000000000112o_B5_F250M_uncal_F250M_seed_image.fits'
obs.paramfile = yamlfile
obs.create()

