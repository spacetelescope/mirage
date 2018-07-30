[![Build Status](https://travis-ci.org/spacetelescope/mirage.svg?branch=master)](https://travis-ci.org/spacetelescope/mirage)

# MIRaGe = Multi Instrument Ramp Generator

This repository contains code that can be used to generate
simulated NIRCam, NIRISS, or FGS data. These data can be in one
of two formats:

`raw` - No calirbations applied. Detector level effects such as non-linearity,
superbias, etc are still present.

`linearized` - Detector level effects have been removed, and data have been
linearized, but are still in ramp format, where multiple non-destructive
reads of the detector are present.


## Installation

To install `mirage`, first clone the repository:
`git clone https://github.com/spacetelescope/mirage.git`

Then, install the package:
```
cd mirage
pip install .
```

## Dependencies

To simulate wide field slitless spectroscopy (WFSS) data:

* [NIRCam_Gsim][d1]: to disperse imaging data
* [GRISM_NIRCAM][d2]: NIRCam-specific grism configuration and sensitivity files
* [GRISMCONF][d3]: grism dispersion polynomials

Background signals:

* [JWST backgrounds][d4]: Generate JWST Exposure Time Calculator-type backgrounds (zodiacal+thermal)

Calibration pipeline:

* [JWST calibration pipeline][d5]. Necessary if using raw dark current exposures as input. Optional otherwise.

[d1]: https://github.com/npirzkal/NIRCAM_Gsim
[d2]: https://github.com/npirzkal/GRISM_NIRCAM
[d3]: https://github.com/npirzkal/GRISMCONF
[d4]: https://github.com/spacetelescope/jwst_backgrounds
[d5]: https://github.com/STScI-JWST/jwst


## Examples

See the notebooks in the `examples` subdirectory. There are notebooks
for imaging simulations, WFSS simulations, moving target
(non-sidereal) simulations, and simulations of OTE commissioning.


## Functionality

The code needed to create simulated data is split into
four general parts:

1. Generate a "seed image"
2. (OPTIONAL) to create WFSS data, disperse the seed image
3. Prepare an existing dark current ramp
4. Combine the seed image with the dark

More details on these steps are given below:

### Generate a "seed image"

This portion of the code generates a "seed image" from
either input source catalogs or an input large field-of-view
observation (fits file).

The seed image is a noiseless countrate image containing
all of the sources specified in the input catalogs. Sources
are scaled to the requested magnitudes and placed in the
appropriate location according to the instrument distortion model.

The main input to the code is a yaml file (examples are
in the 'inputs' subdirectory). This file contains more
information than is actually needed to run the code,
but the information is all needed to run the entire
data simulation process, so for the moment we've kept all
of the inputs.


**To use the code:**

1) From the command line:
```
python catalog_seed_image.py myfile.yaml
```

2) Within python:
```
from mirage.seed_image import catalog_seed_image as csi
cat = csi.Catalog_seed()
cat.paramfile = 'myfile.yaml'
cat.make_seed()
```

**Outputs:**

Multi-extension fits file with name ending in 'seed_image.fits', containing:

Extension 0: empty

Extension 1: seed image

Extension 2: segmentation map

Also, the seed image, segmentation map, and exposure info dictionary are available as:
`self.seedimage`, `self.seed_segmap`, and `self.seedinfo`

### Disperse the seed image

Requires multiple imaging seed images as input. Output is a single, dispersed
seed image that can be passed to later simulator steps just as imaging seed
images. See WFSS notebook for an example.

**To use:**
```
from NIRCAM_Gsim.grism_seed_disperser import Grism_seed
crossing_filter = 'F444W'
module = 'A'    # 'A' or 'B'
direction = 'R' # 'R' for row or 'C' for column
image_seeds = [seed1.seed_file, seed2.seed_file, seed3.seed_file, seed4.seed_file]
dmode = 'mod{}_{}'.format(module.upper(),direction.upper())
loc = os.path.join(datadir,"GRISM_NIRCAM/")
background_file = "{}_{}_back.fits".format(crossing_filter,dmode)
t = Grism_seed(image_seeds,crossing_filter,dmode,config_path=loc)
t.observation()
t.finalize(Back = background_file)
```

### Prepare an existing dark current ramp

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
ramp from Mirage.

**Output:**

The linearized dark current and zeroth frame as saved to a fits file
that uses the name from the Output:file entry in the input yaml file
and ending with '_linearizedDark.fits'.

These are also available as self.linDark and self.zeroModel

**To use:**
python dark_prep.py myinputs.yaml

**or:**
```
from mirage.dark import dark_prep
dark = dark_prep.DarkPrep()
dark.paramfile = 'myinputs.yaml'
dark.prepare()
```

### Combine the seed image with the dark

This step takes as input a seed image (with associated segmentation
map) and a linearized dark current exposure.

The seed image is converted into a 4d signal ramp,
and poisson noise, cosmic rays, and other detector effects
are added.

The ramp is then reorganized into the requested readout
pattern and added to the dark current ramp.

**To use:**
`python obs_generator.py myinputs.yaml`

**or:**
```
from mirage.ramp_generator import obs_generator
obs = obs_generator.Observation()
obs.linDark = 'V42424024002P000000000112o_B5_F250M_uncal_linear_dark_prep_object.fits'
obs.seed = 'V42424024002P000000000112o_B5_F250M_uncal_F250M_seed_image.fits'
obs.paramfile = 'myinputs.yaml'
obs.create()
```

To create a simulated exposure, string together all of the steps:
```
from mirage.seed_image import catalog_seed_image
from mirage.dark import dark_prep
from mirage.ramp_generator import obs_generator

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
```

### Convenience Functions

**`imaging_pipeline.py`** - wrapper around the three steps needed to create an
imaging mode simulated exposure. This also works for moving target
simulations. Example use shown in the imaging notebook.

**`wfss_pipeline.py`** - wrapper around the steps needed to create an WFSS
simulated exposure. Example use shown in the WFSS notebook.

**`yaml_generator.py`** - Beginning with an Astronomer's Proposal Tool (APT) file,
create the yaml files necessary to simulate the entire proposal. Example use
shown in the imaging and WFSS notebooks.