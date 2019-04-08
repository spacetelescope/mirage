.. FALCO documentation master file, created by
   sphinx-quickstart on Thu Jun 21 20:18:02 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

MIRAGE: Multi-Instrument RAmp GEnerator
===============================================

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   install.rst

   catalogs.rst
   catalog_creation.rst

   three_steps.rst
   seed_images.rst
   dark_preparation.rst
   observation_generator.rst

   yaml_generator.rst
   input_yaml_parameters.rst
   example_yaml.rst

   wfss_simulations.rst
   quickstart.ipynb
   api.rst


Overview
--------

Mirage is a Python package that creates realistic simulated data products for JWST's NIRCam, NIRISS, and FGS instruments. The output files are in a format that is as close as possible to that which real JWST data will have, which allows Mirage data to be used for data analysis software development and testing, including the official `JWST calibration pipeline <https://jwst-pipeline.readthedocs.io/en/latest/>`_.

Mirage creates high-fidelity simulated exposures for NIRCam's imaging and wide field slitless spectroscopy (WFSS) modes, NIRISS's imaging, WFSS, and aperture masking interferometery (AMI) modes, and FGS's imaging mode. It supports sidereal as well as non-sidereal tracking (i.e. sources can be made to move across the field of view within an observation).

Astronomical sources to add to the scene are specified through the use of source catalogs.

The software is designed such that a single yaml file can be used as input. This file contains various instrument and readout parameters, as well as lists reference files necessary for the production of the simulated data. Details on the contents of the yaml file are given on the :ref:`Input Yaml File Parameters <input_yaml_file_parameters>` page. One yaml file contains all the information necessary to produce a simulated exposure from a single detector.

For software developement and testing, users may wish to simulate all data from a given JWST proposal. In order to do this, Mirage requires one yaml file for each exposure/detector combination. To facilitate this, the Mirage package includes ancillary convenience functions that allow the user to translate an Astronomer's Proposal Tool (`APT <https://jwst-docs.stsci.edu/display/JPP/JWST+Astronomers+Proposal+Tool%2C+APT>`_) proposal into a series of input yaml files. These tools are described in detail on the :ref:`Simulating Observations from an APT File <from_apt>` page.

Basic Workflow
--------------

The basic order of steps necessary to produce simulated data with Mirage is:

1. Create the appropriate :ref:`Source Catalogs <catalogs>`. Mirage contains :ref:`tools to aid in the creation of source catalogs <catalog_generation>`.
2. Create the needed :ref:`yaml files <example_yaml>` to descirbe the observations. An easy way to do this is to start with an APT proposal and use Mirage's :ref:`yaml file generator <from_apt>` tool.
3. Run the appropriate module. Most often this will be the imaging simulator **add link to imaging notebook** or the :ref:`wfss simulator <wfss_data>`. **change link to wfss notebook** To produce only noiseless, background-free countrate images of the scene, the :ref:`seed image generator module <seed_images>` can be used.

See the :ref:`Examples <examples>` page for basic examples of each of these situations. There are also several notebooks in the Mirage repository showing more in-depth examples of these steps. These include the:

::

    `Catalog generation notebook <>`_
    `Imaging simulator example notebook <>`_
    `WFSS simulator example notebook <>`_

and for more advanced use:

::

    `Commissionning data generation notebook <>`_
    `Moving Target notebook <>`_



Limitations
-----------
Mirage currently does not simulate coronagraphic observations directly. However the user can input an image of an occulted PSF using the :ref:`extended source catalog <extended>`.

When adding PSFs to simulated data, MIRAGE uses the PSFs in the user-provided library. The libraries currently provided in the Mirage reference files collection contain fits files that are 301x301 pixels. This limited size keeps the file sizes small, but cuts off the wings of the PSFs. For very bright sources, the truncated wings will be visible in the data. In addition, when MIRAGE normalizes a PSF to the requested brightness, signal that should be in the truncated wings will instead be in the central 301x301 pixels of the PSF. The total brightness of the source will be as expected, but the signal will be slightly more concentrated in the center of the PSF than it should be.

When creating data for a target that is moving through the field of view, (e.g. a non-sidereal target in a sidereal observation, or vice versa) the velocity of the target is constant. Speed and direction cannot change within an exposure.

Source brightnesses are constant. There is currently no way to simulate a variable brightness source other than by creating multiple exposures of the same source and varying the source brightness from one exposure to the next.

Tangent plane projection is currently not performed when translating from detector pixel x,y values to RA, Dec values. This leads to small errors in the calculated RA, Dec values. This will be corrected in a future version of MIRAGE.


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
