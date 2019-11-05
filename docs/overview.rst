.. _overview:

Overview
========

Mirage is a Python package that creates realistic simulated data products for JWST's NIRCam, NIRISS, and FGS instruments. The output files are in a format that is as close as possible to that which real JWST data will have, which allows Mirage data to be used for data analysis software development and testing, including the official `JWST calibration pipeline <https://jwst-pipeline.readthedocs.io/en/latest/>`_.

Mirage creates simulated exposures for NIRCam's imaging and wide field slitless spectroscopy (WFSS) modes, NIRISS's imaging, WFSS, and aperture masking interferometery (AMI) modes, and FGS's imaging mode. It supports sidereal as well as non-sidereal tracking (i.e. sources can be made to move across the field of view within an observation).

Astronomical sources to add to the scene are specified through the use of source catalogs.

The software is designed such that a single :ref:`yaml file <example_yaml>` can be used as input to create a single simulated exposure from one detector. This file contains various instrument and readout parameters, as well as lists reference files necessary for the production of the simulated data. Details on the contents of the yaml file are given on the :ref:`Example Yaml <example_yaml>` page. One yaml file contains all the information necessary to produce a simulated exposure from a single detector.

For software developement and testing, users may wish to simulate all data from a given JWST proposal. In order to do this, Mirage requires one yaml file for each exposure/detector combination. To facilitate this, the Mirage package includes ancillary convenience functions that allow the user to translate an Astronomer's Proposal Tool (`APT <https://jwst-docs.stsci.edu/display/JPP/JWST+Astronomers+Proposal+Tool%2C+APT>`_) proposal into a series of input yaml files. These tools are described in detail on the :ref:`Simulating Observations from an APT File <from_apt>` page.

Basic Workflow
--------------

The basic order of steps necessary to produce simulated data with Mirage is:

1. Create the appropriate :ref:`Source Catalogs <catalogs>`.
Mirage contains :ref:`tools to aid in the creation of source catalogs <catalog_generation>`.

2. Create the needed :ref:`yaml files <example_yaml>` to descirbe the observations.
An easy way to do this is to start with an APT proposal and use Mirage's :ref:`yaml file generator <from_apt>` tool.

3. Create the simluated data.
Run the appropriate module. Most often this will be the `imaging simulator <https://github.com/spacetelescope/mirage/blob/master/examples/Imaging_simulator_use_examples.ipynb>`_ or the `wfss simulator <https://github.com/spacetelescope/mirage/blob/master/examples/NIRISS_WFSS_data_creation_example.ipynb>`_. To produce only noiseless, background-free countrate images of the scene, the :ref:`seed image generator module <seed_images>` can be used.

See the `notebooks <https://github.com/spacetelescope/mirage/tree/master/examples>`_ in the Mirage repository for basic examples of each of these situations. These include the:

.. parsed-literal::

    `Catalog generation notebook <https://github.com/spacetelescope/mirage/blob/master/examples/Catalog_Generation_Tools.ipynb>`_
    `Imaging simulator example notebook <https://github.com/spacetelescope/mirage/blob/master/examples/Imaging_simulator_use_examples.ipynb>`_
    `WFSS simulator example notebook <https://github.com/spacetelescope/mirage/blob/master/examples/NIRISS_WFSS_data_creation_example.ipynb>`_
    `Moving Target notebook <https://github.com/spacetelescope/mirage/blob/master/examples/MovingTarget_simulator_use_examples.ipynb>`_

and for use with APT programs:

.. parsed-literal::

    `APT program data generation notebook <https://github.com/spacetelescope/mirage/blob/master/examples/APTProgram_simulator_use_examples.ipynb>`_


Under the Hood
--------------

Broadly speaking, Mirage creates simulated data using three main steps. These include the creation of a noiseless "seed image" of the scene, which contains signal only from astronomical sources. Secondly, a dark current exposure is modified to fit the requested readout pattern. And finally, the seed image and dark current exposures are combined, and noise/cosmic rays are added. See the :ref:`Three Stages <stages>` page for details.


Limitations
-----------
Mirage currently does not simulate coronagraphic observations directly. However the user can input an image of an occulted PSF using the :ref:`extended source catalog <extended_obj>`.

When adding PSFs to simulated data, MIRAGE uses the PSFs from a pre-computed library. The libraries currently provided in the Mirage reference files collection contain fits files that are several hundred pixels on a side. This limited size keeps the file sizes small, but cuts off the wings of the PSFs. For very bright sources, the truncated wings will be visible in the data. In addition, when MIRAGE normalizes a PSF to the requested brightness, signal that should be in the truncated wings will instead be in the pixels of the truncated PSF. The total brightness of the source will be as expected, but the signal will be slightly more concentrated in the center of the PSF than it should be.

When creating data for a target that is moving through the field of view, (e.g. a non-sidereal target in a sidereal observation, or vice versa) the velocity of the target is constant. Speed and direction cannot change within an exposure.

Source brightnesses are constant. There is currently no way to simulate a variable brightness source other than by creating multiple exposures of the same source and varying the source brightness from one exposure to the next.

The *yaml_generator* function currently is only able to parse a subset of all available APT templates. See the :ref:`yaml generator <from_apt>` page for details. If your proposal contains observations that use unsupported APT templates, the easiest work-around at the moment is to make a copy of your APT file and strip out the unsupported observations.

.. admonition:: Getting Help

   For help installing or running Mirage, or if something is not clear in this documentation, please `open an issue on the Mirage github page <https://github.com/spacetelescope/mirage/issues>`_.
