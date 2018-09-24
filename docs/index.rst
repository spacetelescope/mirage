.. FALCO documentation master file, created by
   sphinx-quickstart on Thu Jun 21 20:18:02 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

MIRAGE: Multi-Instrument RAmp GEnerator
===============================================

Mirage is a Python package that creates realistic simulated data products for JWST's NIRCam, NIRISS, and FGS instruments. The output files are in a format that is as close as possible to that which real JWST data will have, which allows Mirage data to be used for data analysis software development and testing, including the official `JWST calibration pipeline <https://jwst-pipeline.readthedocs.io/en/latest/>`_.

Overview
--------
Mirage creates high-fidelity simulated exposures for NIRCam's imaging and wide field slitless spectroscopy (WFSS) modes, NIRISS's imaging, WFSS, and aperture masking interferometery (AMI) modes, and FGS's imaging mode. It supports sidereal as well as non-sidereal tracking (i.e. sources can be made to move across the field of view within an observation). Astronomical sources to add to the scene are specified through the use of source catalogs.

The software is designed such that a single yaml file can be used as input. This file contains various instrument and readout parameters, as well as lists reference files necessary for the production of the simulated data. Details on the contents of the yaml file are given in section 10.

The simulator can be broadly divided into three stages. The first stage is the creation of the “seed image”. This is a noiseless count rate image that contains signal only from the astronomical sources to be simulated. Seed images include instrument distortion effects, so that given RA, Dec values are properly converted to pixel x,y values, with the exception of tangent plane projection, which will be added soon. Mirage currently contains two methods for the construction of seed images:

1. Through the use of source catalog files
2. Extraction from a fits file containing a distortion-free image (e.g. HUDF, GOODS, CANDELS, etc)

For WFSS observations, multiple seed images and optional input spectra are then fed into the disperser software, which simulates the effects of the grism by dispersing the signal from the astronomical sources, creating a “grism seed image”.

The second stage is the preparation of the dark current exposure to use for the simulation. The input dark current exposure is reorganized into the requested readout pattern and number of groups and cropped to the requested subarray size. Detector non-linearity effects are then removed using the initial steps of the JWST calibration pipeline.

By using actual dark current exposures from ground testing, Mirage is able to capture many effects which are specific to the instrument and detector being simulated. For example, the 1/f noise, bias structure, and hot pixel population.

The final stage involves the combination of the seed image and the dark current in order to produce the output exposure. The seed image is expanded into integrations with groups that follow the requested readout pattern. Other effects are also added at this stage, including cosmic rays, interpixel capacitance (IPC) and crosstalk effects.

In addition, the Mirage package includes ancillary convenience functions that allow the user to translate an `Astronomer’s Proposal Tool (APT) <https://jwst-docs.stsci.edu/display/JPP/JWST+Astronomers+Proposal+Tool%2C+APT>`_ proposal into a series of input yaml files for the simulator.


Limitations
-----------
Mirage currently does not simulate coronagraphic observations.

When adding PSFs to simulated data, MIRAGE uses the PSFs in the user-provided library. The libraries currently provided in the Mirage reference files collection contain fits files that are 301x301 pixels. This limited size keeps the file sizes small, but cuts off the wings of the PSFs. For very bright sources, the truncated wings will be visible in the data. In addition, when MIRAGE normalizes a PSF to the requested brightness, signal that should be in the truncated wings will instead be in the central 301x301 pixels of the PSF. The total brightness of the source will be as expected, but the signal will be slightly more concentrated in the center of the PSF than it should be.

When creating data for a target that is moving through the field of view, (e.g. a non-sidereal target in a sidereal observation, or vice versa) the velocity of the target is constant. Speed and direction cannot change within an exposure.

Source brightnesses are constant. There is currently no way to simulate a variable brightness source other than by creating multiple exposures of the same source and varying the source brightness from one exposure to the next.

Tangent plane projection is currently not performed when translating from detector pixel x,y values to RA, Dec values. This leads to small errors in the calculated RA, Dec values. This will be corrected in a future version of MIRAGE.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   install.rst
   seed_images.rst
   catalogs.rst
   dark_preparation.rst
   observation_generator.rst
   yaml_generator.rst
   example_yaml.rst
   quickstart.ipynb
   api.rst

.. admonition:: Getting Help

   For help installing or running Mirage, please `open an issue on the Mirage github page <https://github.com/spacetelescope/mirage/issues>`_.

Contributors
------------
Mirage is based on early NIRISS simulator software written by Kevin Volk. It has been developed by a group of core contributors from STScI:

- `Bryan Hilbert <https://github.com/bhilbert4>`_
- `Kevin Volk <https://github.com/KevinVolkSTScI>`_
- `Lauren Chambers <https://github.com/laurenmarietta>`_
- `Johannes Sahlmann <https://github.com/Johannes-Sahlmann>`_
- `Shannon Osborne <https://github.com/shanosborne>`_
- `Marshall Perrin <https://github.com/mperrin>`_
- `Nor Pirzkal <https://github.com/npirzkal>`_


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
