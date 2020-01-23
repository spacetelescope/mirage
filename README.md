# MIRaGe = Multi Instrument Ramp Generator

[![Build Status](https://travis-ci.org/spacetelescope/mirage.svg?branch=master)](https://travis-ci.org/spacetelescope/mirage)
[![License](https://img.shields.io/pypi/l/Django.svg)](https://github.com/spacetelescope/mirage/blob/master/LICENSE.txt)
[![Python](https://img.shields.io/badge/Python-3.6-blue.svg)](https://www.python.org/)
[![STScI](https://img.shields.io/badge/powered%20by-STScI-blue.svg?colorA=707170&colorB=3e8ddd&style=flat)](http://www.stsci.edu)
[![DOI](https://zenodo.org/badge/109982775.svg)](https://zenodo.org/badge/latestdoi/109982775)


This repository contains code that can be used to generate
simulated NIRCam, NIRISS, or FGS data. These data can be in one
of two formats:

`raw` - No calibrations applied. Detector level effects such as non-linearity,
superbias, etc are still present.

`linearized` - Detector level effects have been removed, and data have been
linearized, but are still in ramp format, where multiple non-destructive
reads of the detector are present.

## Installation and Documentation

The [main documentation for Mirage](https://mirage-data-simulator.readthedocs.io/en/latest/) is located on ReadTheDocs.
[Detailed installation instructions](https://mirage-data-simulator.readthedocs.io/en/latest/install.html) can be found there.


## Examples

See the notebooks in the `examples` subdirectory. There are notebooks
for imaging simulations, WFSS simulations, moving target
(non-sidereal) simulations, and simulations of OTE commissioning.


### Citation
If you find this package useful, please consider citing the Zenodo record using the DOI badge above.
Please find additional citation instructions in [CITATION](CITATION).


## Contributing

Prior to contributing to the `mirage` development, please review our [style guide](https://github.com/spacetelescope/mirage/blob/master/style_guide/style_guide.md).

The following is a bare bones example of a best work flow for contributing to the project:

1. Create a fork off of the `spacetelescope` `mirage` repository.
2. Make a local clone of your fork.
3. Ensure your personal fork is pointing `upstream` properly.
4. Create a branch on that personal fork.
5. Make your software changes.
6. Push that branch to your personal GitHub repository (i.e. `origin`).
7. On the `spacetelescope` `mirage` repository, create a pull request that merges the branch into `spacetelescope:master`.
8. Assign a reviewer from the team for the pull request.
9. Iterate with the reviewer over any needed changes until the reviewer accepts and merges your branch.
10. Delete your local copy of your branch.


## Code of Conduct

Users and contributors to the `mirage` repository should adhere to the [Code of Conduct](https://github.com/spacetelescope/mirage/blob/master/CODE_OF_CONDUCT.md).  Any issues or violations pertaining to the Code of Conduct should be brought to the attention of a `mirage` team member or to `conduct@stsci.edu`.


## Questions

For any questions about the `mirage` project or its software or documentation, please [open an Issue](https://github.com/spacetelescope/mirage/issues).


## Current Development Team
- Bryan Hilbert [@bhilbert4](https://github.com/bhilbert4)
- Johannes Sahlmann [@Johannes-Sahlmann](https://github.com/johannes-sahlmann)
- Lauren Chambers [@laurenmarietta](https://github.com/laurenmarietta)
- Kevin Volk [@KevinVolkSTScI](https://github.com/KevinVolkSTScI)
- Shannon Osborne [@shanosborne](https://github.com/shanosborne)
- Marshall Perrin [@mperrin](https://github.com/mperrin)


## Acknowledgments:
`Mirage` is based on a NIRISS data simulator originally written by Kevin Volk.
