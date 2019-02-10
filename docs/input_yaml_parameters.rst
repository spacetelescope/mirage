.. _input_yaml_file_parameters:

Input yaml file parameters
==========================

The yaml file used as the primary input to mirage contains a large number of telescope and instrument settings, reference files, and catalogs. The :ref:`example yaml<example_yaml>` shows the format of the input file. Here we decribed in detail each of the input lines. The inputs are grouped by section:

1. :ref:`Inst section <inst_section>`
2. :ref:`Readout section <readout_section>`
3. :ref:`Reffiles section <reffiles_section>`
4. :ref:`nonlin section <nonlin_section>`
5. :ref:`cosmicRay section <cosmicray_section>`
6. :ref:`simSignals section <simsignals_section>`
7. :ref:`Telescope section <telescope_section>`
8. :ref:`newRamp section <newramp_section>`
9. :ref:`Output section <output_section>`

.. _inst_section:

Inst secton
-----------

This section of the input yaml file contains information about the instrument being simulated.

instrument
++++++++++

The name of the JWST instrument to be simulated. The simulator will only function if ‘NIRCam’, ‘NIRISS’, or ‘FGS’ is placed in this field.

mode
++++

The observing mode to be simulated. There are three valid options for this field. “imaging” will create imaging data, “wfss” will produce wide field slitless spectroscopic data. The other accepted input is "ami" when simulating NIRISS, although this mode is functionally identical to the use of "imaging".

use_JWST_pipeline
+++++++++++++++++

True/False. Set to False if you wish to proceed without using any JWST pipeline functions. In this case, the input dark current exposure must already be linearized, as the pipeline is used for the linearization process. True is recommneded.

.. _readout_section:

Readout section
---------------

This section of the yaml file contains inputs describing the details of the exposure, including the readout pattern, filter, subarray, etc to use.

readpatt
++++++++

This is the name of the readout timing pattern used for the output simulated exposure. Examples for NIRCam include RAPID, BRIGHT1, BRIGHT2, and DEEP8. Each pattern averages and skips a predefined number of frames when constructing each group of an integration. The list of possible readout patterns and their definitions is provided by an ascii file specified in the **readpattdefs** parameter in the **Reffiles** section of the input file. A more detailed description of readout patterns is given in the detector readout pages for `NIRCam <https://jwst-docs.stsci.edu/display/JTI/JWST+Field+of+View>`_,  `NIRISS <https://jwst-docs.stsci.edu/display/JTI/JWST+Field+of+View>`_, and `FGS <https://jwst-docs.stsci.edu/display/JTI/JWST+Field+of+View>`_.

ngroup
++++++

This parameter lists the number of groups comprising the output integration.

nint
++++

The number of integrations in the output exposure. Each integration is composed of **ngroup** groups. Note that currently, any observation containing a moving target (non-sidereal observation with trailed sidereal objects, or vice versa) cannot have an nint value greater than 1. **(IS THIS STILL TRUE?)**

resets_bet_ints
+++++++++++++++

The number of detector resets between integrations within a single exposure. For all instruments, this should be set to 1.

array_name
++++++++++

This is the name of the aperture used for the simulated data. Generally, this is composed of the name of the detector combined with the name of the subarray used. For example, a full frame observation using NIRCam's A1 detector has an **array_name** of 'NRCA1_FULL', while a full frame NIRISS observation will have an array_name of ‘NIS_CEN’. The list of possible array_name values are given in the **subarray_defs** input file described below. The **array_name** is used to identify several other characteristics of the simulated data, including the detector to use, as well as the proper array dimensions and location on the detector.

filter
++++++

The name of the filter wheel element to use for the simulated data. (e.g. F444W). The filter is used when scaling astronomical sources from the requested brightness in magnitudes to counts on the detector. For NIRCam simulations, the filter name is also used to determine whether the simulated data are to be produced using a shortwave or longwave detector. Lists of instrument filters can be found on the `NIRCam <https://jwst-docs.stsci.edu/display/JTI/JWST+Field+of+View>`_,  `NIRISS <https://jwst-docs.stsci.edu/display/JTI/JWST+Field+of+View>`_, and `FGS <https://jwst-docs.stsci.edu/display/JTI/JWST+Field+of+View>`_ filter pages.

pupil
+++++

The name of the pupil wheel element to use for the simulated data. Some filters for both NIRCam and NIRISS reside in their respective pupil wheels. Therefore this entry is checked when deciding upon scaling factors for simulated sources. Pupil wheel elements are desribed in the `NIRCam <https://jwst-docs.stsci.edu/display/JTI/JWST+Field+of+View>`_,  `NIRISS <https://jwst-docs.stsci.edu/display/JTI/JWST+Field+of+View>`_, and `FGS <https://jwst-docs.stsci.edu/display/JTI/JWST+Field+of+View>`_ pupil wheel pages.

.. _reffiles_section:

Reffiles
--------

This section of the input file lists the various reference files needed for the various steps of the simulator to run.

dark
++++

The name of the raw dark current file that will be used as the basis for the simulated exposure. This file must be in raw format, such that no JWST calibration pipeline steps have been applied to the data. If an already-linearized dark current integration is to be used, that file name should be placed in the **linearized_darkfile** field below. Note that the **linearized_darkfile** entry will take precedence. Only if that is set to __None__ will the file listed in this field be used.

The dark current integration must have a readout pattern of either RAPID/NISRAPID/FGSRAPID or a value identical to that of the integration to be simulated. RAPID/NISRAPID/FGSRAPID data keep every readout frame with no averaging. From this, any other readout pattern can be simulated by averaging and skipping the appropriate frames. Other readout patterns cannot be translated in this way as their data are already averaged or missing some frames. However if simulating, for example a BRIGHT2 integration, then the input dark current integration can be a BRIGHT2 integration, as no translation is necessary in this case.

If a translation between RAPID and another readout pattern is necessary, then frames will be averaged/skipped as necessary. If the input dark current integration does not contain enough frames to be translated into the requested number of output groups, then the script creates enough additional dark current frames to make up the difference. These additional frames are created by making a copy of an appropriate number of existing initial dark current frames, and adding their signals to that in the final dark current frame. Note that this can lead to apparent double cosmic rays in pixels where a cosmic ray appeared in the dark current integration.

.. hint::
	This input can only be used if **use_JWST_pipeline** is set to True.

.. hint::
	The collection of reference files associated with Mirage contains a small library of raw dark current exposures that can be used.

linearized_darkfile
+++++++++++++++++++

The name of a linearized dark current integration to use as input for the simulated data. This file should contain a dark integration that has been processed through the superbias subtraction, reference pixel subtraction, and linearity steps of the JWST calibration pipeline. The resulting linearized signal must be saved in an extension with the name 'SCI'. Also, the subtracted signal from the superbias and reference pixels must be saved in an extension called 'SBANDREFPIX'. This output will be produced and saved for a given dark current file by Mirage.

Using this input rather than the uncalibrated dark above can save significant computing time, especially in the case of creating many output exposures.

.. hint::
	This input can be used for **use_JWST_pipeline** set to True or False.

.. hint::
	The collection of reference files associated with Mirage contains a small library of linearized dark current products that can be used.

badpixmask
++++++++++

If a linearized dark current file is to be used and a linearized output file is requested, this optional bad pixel mask can be used to populate the data quality array in the output simulated data file. The file must be in the format for JWST bad pixel masks that is used by the JWST calibration pipeline.

**ADD LINK HERE TO FILE FORMAT**

.. hint::
	The collection of reference files associated with Mirage contains a library of bad pixel masks that can be used.

superbias
+++++++++

The superbias reference file for the detector of the simulation. This file must match the format of the JWST pipeline superbias reference file. **ADD LINK HERE** If the input dark current integration is a raw file then this superbias file is used to subtract the superbias from the dark. If the input dark is already linearized, this superbias file is not used.

.. hint::
	The collection of reference files associated with Mirage contains a library of superbias files that can be used.

linearity
+++++++++

Name of the reference file containing the linearity correction coefficients. This file must be in the format expected by the JWST calibration pipeline. ** ADD LINK HERE** If the input dark current integration is raw, the coefficients contained in this file are used to linearize the dark current after subtracting the superbias and reference pixel signal. These coefficients are also used to "unlinearize" the final simulated exposure if a raw simulated observation is requested.

In addition, the coefficients in this file are used to linearize the values in the saturation reference file, such that saturated signals in the linear simulated exposure can be found.

.. hint::
	The collection of reference files associated with Mirage contains a library of linearity coefficient files that can be used.

.. _saturation:

saturation
++++++++++

Name of the reference file containing a map of the saturation signal level for all pixels. If the input dark current integration is raw, this file is used by the calibration pipeline to flag saturated pixels in the dark current integration prior to linearizing.

This saturation map, after being linearized, is also used to search for saturated signal values in the combined dark current/simulated source exposure prior to unlinearizing.

.. hint::
	The collection of reference files associated with Mirage contains a library of saturation map files that can be used.

gain
++++

Name of the file containing the gain map appropriate for the detector being used. The gain is used to translate the cosmic rays, which are in units of electrons, to units of ADU prior to adding them to the simulated data.

.. hint::
	The collection of reference files associated with Mirage contains a library of gain map files that can be used.

pixelflat
+++++++++

Name of the pixel flat file to use. Once the simulated integration is created, the result is multiplied by the pixel flat. This is done to un-flatten the image.


illumflat
+++++++++

Name of the illumination flat to use. Once the simulated integration is created, the result is multiplied by the illumination flat.


astrometric
+++++++++++

Name of the astrometric distortion reference file to use for including the effects of distortion in the simulated data.  This file is used to translate input source locations between RA and Dec coordinates and pixel x and y coordinates, and vice versa. This file must be in asdf format and match that expected by the calibration pipeline. **ADD LINK HERE**

.. hint::
	The collection of reference files associated with Mirage contains a library of distortion reference files that can be used.

ipc
+++

File containing the interpixel capacitance (IPC) kernel to apply to the simulated data in order to introduce IPC effects. After all simulated objects have been added to a count rate image, the image is convolved with the IPC kernel. The IPC file must be a fits file with the IPC kernel located in the first (rather than 0th) extension. Typical JWST IPC reference file kernels are a 3x3 array, but Mirage supports kernels of any odd-numbered size, as well as 4-dimensional kernels, where there is a separate 2-dimensional kernel for each pixel. In order to introduce, rather than remove, IPC effects, the kernel must be normalized and have a value in the central pixel which is less than 1.0. This is the inverse of the kernel used in the JWST calibration pipeline IPC removal step, where the central pixel has a value greater than 1.0, and negative values in surrounding pixels. For the simulator, the user can specify a calibration pipeline-formatted kernel file, and then set the **invertIPC** flag below to True, in which case the kernel will be inverted before using.

.. hint::
	The collection of reference files associated with Mirage contains a library of IPC kernel files that can be used.

invertIPC
+++++++++

If set to True, the IPC kernel supplied through the ipc entry is inverted before convolving with the signal rate image. JWST IPC kernel reference files contain the kernel necessary to remove IPC from the data. Therefore these kernels must be inverted before they can add IPC effects to the data in the simulator.


pixelAreaMap
++++++++++++

Fits file containing the pixel area map for the detector to be simulated. If provided, the pixel area map is multiplied into the seed image at a point when the seed image contains only extended sources. Point sources do not have the pixel area map applied to them. **DESCRIBE IN MORE DETAIL WHAT'S GOING ON** The pixel area map file must be in the format of the JWST pixel area map reference file. **ADD LINK HERE**

.. hint::
	The collection of reference files associated with Mirage contains a library of pixel area map files that can be used.

subarray_defs
+++++++++++++

Name of a whitespace-delimited ascii file that lists all of the possible supported subarray apertures. This file is provided with the MIRAGE repository, in the config subdirectory.

.. hint::
	To use the subarray definition files packaged with Mirage, set this to **config** in the input yaml file. This is the default when creating yaml files from an APT file using the :ref:`yaml generator <yaml_generator>`

For each subarray, the file must list the full aperture name (e.g. NRCA1_FULL) as well as the corresponding name used in proposal planning (e.g. FULL), **??as well as the detector on which that subarray is used, and the necessary filter for the subarray (for the special GRISM entries).??**  and the number of amplifiers used in the readout.

readpattdefs
++++++++++++

Ascii file which gives the definitions of the possible readout patterns for the instrument. For each readout pattern, the number of frames averaged to create each group (nframe) and the number of frames skipped beteren each group (nskip) must be specified, as well as the maximum number of allowed groups. For a given readout pattern the simulator will search the entries in this file in order to determine the proper nframe and nskip values to use. The current lists of acceptable NIRCam and NIRISS readout patterns are given on the NIRCam  and NIRISS  detector readouts webpages. These files for all instruments are provided with the MIRAGE repository, in the config subdirectory.

.. hint::
	To use the readout pattern definition files packaged with Mirage, set this to **config** in the input yaml file. This is the default when creating yaml files from an APT file using the :ref:`yaml generator <yaml_generator>`

crosstalk
+++++++++

Ascii file containing crosstalk coefficients. Crosstalk is only applied to data read out through more than one amplifer. The file contains one row for each detector. Each row contains all of the coefficients necessary to fully describe crosstalk. This file is contained in the MIRAGE repository, in the **config** subdirectory.

.. hint::
	To use the crosstalk coefficient files packaged with Mirage, set this to **config** in the input yaml file. This is the default when creating yaml files from an APT file using the :ref:`yaml generator <yaml_generator>`

filtpupilcombo
++++++++++++++

Name of an ascii file containing a list of the filter and pupil wheel elements in place when requesting simulated data for a given filter. This information is used to apply the appropriate conversion between magnitudes and counts when reading in source catalogs. This flux calibration is also added to the header of the seed image, as it is used when seed images are dispersed during the simulation of WFSS data. This file is present in the config subdirectory of the MIRAGE repository.

.. hint::
	To use the filter and pupil wheel definition files packaged with Mirage, set this to **config** in the input yaml file. This is the default when creating yaml files from an APT file using the :ref:`yaml generator <yaml_generator>`

flux_cal
++++++++

Ascii file that lists flux conversion factors and the pivot wavelength associated with each filter. Conversion factors include ABMAG, STMAG, and VEGAMAG to counts per second, as well as FLAM (erg s^(-1) 〖cm〗^(-2) Å^(-1)) and FNU (erg s^(-1) 〖cm〗^(-2) 〖Hz〗^(-1)) to counts per second. This file is used when producing seed images to be fed into the grism disperser code, as well as for translating catalog sources from magnitudes to counts per second. This file is provided with the MIRAGE repository, in the config subdirectory.

.. hint::
	To use the flux calibration files packaged with Mirage, set this to **config** in the input yaml file. This is the default when creating yaml files from an APT file using the :ref:`yaml generator <yaml_generator>`

.. _nonlin_section:

Nonlin
------

The following input fields describe how non-linearity is treated in the input and simulated data.

limit
+++++

Signal limit, in units of ADU, above which the linearity correction is not applied. Pixels with signals above this limit are considered saturated. This single value across the entire detector is only used if a :ref:`saturation reference file <saturation>` is not provided.

accuracy
++++++++

When introducing non-linearity back into the linear data, the Newton-Raphson method is used to essentially run the JWST calibration pipline’s linearity correction step in reverse. The value of this accuracy parameter is the threshold below which the solution is considered to have converged. For example, an accuracy threshold of 0.000001 means that the unlinearization is considered complete when the ratio of the signal values from one iteration to the next is less than 1.000001.

maxiter
+++++++

The maximum number of iterations of the Newton-Raphson method to use when introducing non-linearity back into the data before declaring failure. Default is 10.

robberto
++++++++

If set to False, the simulator assumes that the non-linearity correction function and coefficients match those used in the JWST calibration pipeline. If set to True, the script assumes an alternate linearity function, as defined in Robberto (`2010 <https://jwst.stsci.edu/files/live/sites/jwst/files/home/instrumentation/technical%20documents/JWST-STScI-002163.pdf>`_ , `2011 <https://jwst.stsci.edu/files/live/sites/jwst/files/home/instrumentation/technical%20documents/JWST-STScI-002344.pdf>`_). Currently, no coefficients for the latter method exist, implying this parameter should be set to False

.. _cosmicray_section:

cosmicRay
---------

Input parameters in this section describe how cosmic rays are added to the simulated data.

path
++++

Path of the location of the cosmic ray library to use. The code was developed around the cosmic ray library produced by Robberto (`2009 <https://jwst.stsci.edu/files/live/sites/jwst/files/home/instrumentation/technical%20documents/JWST-STScI-001928.pdf>`_). This library is included in the collection of `reference files <>`_ associated with Mirage. After extracting the library from the tar file, set this path to point to the top level directory of the cosmic ray library.

library
+++++++

Specification of which cosmic ray library to choose cosmic rays from. Options are SUNMIN, SUNMAX, FLARE, each of which assumes a different cosmic ray rate. Details on the three types of libraries are given in Robberto (`2009 <https://jwst.stsci.edu/files/live/sites/jwst/files/home/instrumentation/technical%20documents/JWST-STScI-001928.pdf>`_).

scale
+++++

Scaling factor to apply to the cosmic ray rate. For example, to simulate cosmic rays at a rate twice as high as that in SUNMIN, set library to SUNMIN and scale to 2.0

suffix
++++++

Filename suffix of the cosmic ray library files. The code was developed around files with the suffix of ‘IPC_NIRCam_XX’ where XX is the detector (e.g. B5) for NIRCam, ‘IPC_NIRISS_NIS’ for NIRISS, and ‘IPC_FGS_GUIDERy’ where y is 1 or 2, for FGS. These cosmic ray files are included in Mirage's `reference file collection <>`_. This field will be automatically populated with the correct suffix when creating yaml files using the :ref:`yaml generator <yaml_generator>`.

seed
++++

Random number generator seed to use when selecting cosmic rays to add.

.. _simsignals_section:

simSignals
----------

This section of the input file describes how sources and other signals are added to the simulated data.

.. _ptsrc:

pointsource
+++++++++++

Name of an ascii catalog file listing point sources to add to the simulated image. An example :ref:`point source <point_source>` catalog is provided on the :ref:`Catalogs page <catalogs>`.

psfpath
+++++++

Path name to the PSF library to be used for adding point sources to the data. The code was developed around a PSF library constructed using WebbPSF (Perrin, 2014). This PSF library is included in the collection of Mirage `reference files <>`_ . Once that package is downloaded and the data files extracted from the tar file, set this field to point to the top-level directory of the PSF library.

psfbasename
+++++++++++

Basename of the files in the PSF library. When using the default libraries that are distributed with Mirage, this should be set to the name of the instrument.

psfpixfrac
++++++++++

It is assumed that the PSF library contains a grid of PSFs centered at various sub-pixel locations. This parameter specifies the resolution of this grid. For example, if the library contains PSFs centered at every 0.25 pixels across a pixel in x and y, then this field should be set to 0.25. In the current collection of Mirage `reference files <>`_ the PSF library for NIRCam uses a resolution of 0.25, while those for NIRISS and FGS have a resolution of 0.1 pixels.

psfwfe
++++++

PSF wavefront error value to use when choosing PSF files from the PSF library. The current PSF libraries distributed with the Mirage `reference files <>`_ have two options for wavefront error: “predicted” and “requirements”. These two values represent the predicted in-flight wavefront errors, and the maximum allowed wavefront errors, respectively.

psfwfegroup
+++++++++++

The current PSF library contains 5 different realizations for each filter/wavefront error-specified PSF. In this field, place the realization number to use. With 5 realizations present in the library, this field can have a value of 0 through 4.

.. _galaxylistfile:

galaxyListFile
++++++++++++++

Similar to the :ref:`pointsource <ptsrc>` entry, this is an ascii catalog file containing a list of the galaxies to simulate in the data. See the :ref:`galaxies <galaxies>` entry on the :ref:`catalogs <catalogs>` page for an example of this file.

.. _extendedlist:

extended
++++++++

Name of an ascii file containing a list of "extended images" to add to the simulated data. These are stamp image of sources, contained in small fits files. These stamp images are read in, scaled to the requested magnitude, and added to the seed image.  This is a way to add objects other than point sources or 2D Sersic profiles to the data. The :ref:`extended catalog <extended>` section of the :ref:`catalogs <catalogs>` page shows an example extended source catalog.

extendedScale
+++++++++++++

Multiplicative factor by which to scale the data in the extended image file before adding to the simulated data. The extended image is multiplied by this factor **if the magnitude is set to None in the extended catalog file**.

PSFConvolveExtended
+++++++++++++++++++

True/False. Convolve the extended image with the appropriate instrumental PSF prior to adding to the output image.

movingTargetList
++++++++++++++++

Similar to the :ref:`point source <ptsrc>` list file, this is a file containing a list of targets to treat as moving (non-sidereal) targets.  These sources will move through the field of view as the exposure progresses. This is the list to use if you wish to insert an asteroid or KBO that is moving through the field of view of your observation. See the :ref:`moving point source <moving_point_source>` section on the :ref:`Catalogs <catalogs>` page for an example.

movingTargetSersic
++++++++++++++++++

Similar to the :ref:`galaxy target list file <galaxylistfile>`, this file contains a list of galaxies (2D Sersic profiles) to be used as moving targets. These sources will move through the background of the simulated data. This may be useful for inserting a resolved moon/asteroid into the scene. An example file is shown in the :ref:`Moving Sersic <moving_sersic>` section of the :ref:`Catalogs <catalogs>` page.

movingTargetExtended
++++++++++++++++++++

Similar to the :ref:`extended <extendedlist>` target list, this is an ascii file listing extended targets to move through the background of the image. A description and example of this file are shown in the :ref:`Moving Extended <moving_extended>` section of the :ref:`Catalogs <catalogs>` page.

movingTargetConvolveExtended
++++++++++++++++++++++++++++

Set this input to True if you wish to convolve the images listed in **movingTargetExtended** with the instrumental PSF prior to adding them to the simulated data.

movingTargetToTrack
+++++++++++++++++++

This ascii catalog file is used for what are traditionally (in HST jargon) called 'moving targets'.  Targets listed in this file are treated as non-sidereal targets that JWST will track during the simulated observation. In this case, the target listed in this file will appear static in the output data, but all other sources (e.g. those listed in :ref:`pointSource <ptsrc>`, :ref:`galaxyListFile <galaxylistfile>`, and :ref:`extended <extendedlist>`) will all appear trailed through the data. A description and example of the file are shown in the :ref:`Non-sidereal Source <nonsidereal>` section on the :ref:`Catalogs <catalogs>` page.

.. _zodiacal:

zodiacal
++++++++

Name of a file containing a 2 dimensional count rate image of zodiacal light. This file is read in, scaled by the :ref:`zodiscale <zodiscale>` value, and added to the seed image. Leave as None to skip this step. Note that the :ref:`bkgdrate <bkgdrate>` input parameter, when set to “high”, “medium”, or “low”, will return a background rate image that includes the contribution from zodiacal light, in which case this step should be set to None. The behaviors of this step and the scattered step below are very basic, and identical. There are no requirements on what the count rate images in these files must look like.

.. _zodiscale:

zodiscale
+++++++++

Scaling factor to multiply the :ref:`zodiacal light count rate image <zodiacal>` by prior to adding to the output data.

.. _scattered:

scattered
+++++++++

Scattered light count rate image file. This file is assumed to contain a 2-dimensional array of signals in units of ADU per second. The file is read in, scaled by the :ref:`scatteredscale <scatteredscale>` value, and added to the seed image. Leave as None to skip this step.

.. _scatteredscale:

scatteredscale
++++++++++++++

Scaling factor to multiply the :ref:`scattered light count rate image <scattered>` by prior to adding to the seed image.

.. _bkgdrate:

bkgdrate
++++++++

Constant (across all pixels) background count rate to add to the output data. The value can be a number, in which case it is assumed to have units of counts per pixel per second. Alternately, the value can be “high”, “medium”, or “low”. If one of these options is used, the simulator uses the `jwst_backgrounds <https://github.com/spacetelescope/jwst_backgrounds>`_ repository to calculate the background rate to apply to the simulated data. The package calculates the background signal at the requested pointing on the sky for each night over the course of a year and creates a histogram of these values. If the requested background is "low" then the returned background level is equal to that of the 10th percentile in the histogram. A "medium" background corresponds to the 50th percentile value, and "high" is the 90th percentile value. In this case, the returned background rate includes contributions from zodiacal light and telescope thermal emission.

Note that background rates associated with the "low", "medium", and "high" values are calculated in the same way as when they are used in the `JWST ETC <https://jwst.etc.stsci.edu/>`_.

poissonseed
+++++++++++

Random number generator seed used for Poisson simulation

.. _telescope_section:

Telescope
---------

Inputs in this section of the yaml file describe the telescope pointing to use for the simulation.

ra
+++

Right ascension of the observation. This will be the RA at the reference location on the detector being used for the simulation. The reference location varies with the requested subarray, but is generally in the center of the field of view. This input can be a string "HH:MM:SS.sss", or a float in decimal degrees.

dec
+++

Declination of the observation. This will be the Dec at the reference location on the detector. The reference location varies with the requested subarray, but is generally in the center of the field of view. This input can be a string "DD:MM:SS.sss" or a float in decimal degrees.

rotation
++++++++

Rotation of the y-axis in degrees East of North. Currently this rotation is defined around the reference location of the chosen subarray.

.. _newramp_section:

newRamp
-------

This section of the input file lists JWST calibration pipeline-style configuration files that may be needed when preparing the simulated data. Copies of all of these configuration files are included in the ‘config’ subdirectory of the MIRAGE repository. Therefore, unless you wish to use your own set of configuration files, you can set these fields all to 'config'. This is the default behavior when creating yaml files via the :ref:`yaml generator <yaml_generator>`.

.. hint::
	In order to create your own set of pipeline configuration files, use the shell command:

	> collect_pipeline_cfg /your/destination/directory

dq_configfile
+++++++++++++

Name of the JWST calibration pipeline configuration file to be used in the dq_init step when it is run on the raw dark current integration.

sat_configfile
++++++++++++++

Name of the JWST calibration pipeline configuration file to be used in the saturation step when it is run on the raw dark current integration.

superbias_configfile
++++++++++++++++++++

Name of the JWST calibration pipeline configuration file to be used in the superbias step when it is run on the raw dark current integration.

refpix_configfile
+++++++++++++++++

Name of the JWST calibration pipeline configuration file to be used in the reference pixel subtraction step when it is run on the raw dark current integration.

.. hint::
    If you choose to use your own reference pixel correction configuration file, we recommend setting the **odd_even_rows** entry to False, as this correction is not typically performed on NIRCam, NISISS, or FGS data.

linear_configfile
+++++++++++++++++

Name of the JWST calibration pipeline configuration file to be used in the linearity correction step when it is run on the raw dark current integration.

.. _output_section:

Output
------

This section of the yaml file contains information about the output file, such as filename and location. In addition, this section contains a large number of fields that describe how this particular exposure fits within an observing program/proposal. This information is not used during the creation of the simulated data, but is placed in the header of the output file in order to be consistent with the contents of real JWST data files. In addition, `level 3 of the JWST calibration pipeline <https://jwst-pipeline.readthedocs.io/en/stable/jwst/pipeline/description.html#pipelines>`_, which is used to combine multiple exposures into mosaic images, does require some of this information. The easiest way to correctly populate this information in the simulator yaml files is to :ref:`create the yaml files from an APT file via yaml_generator.py<from_apt>`, in which case the fields are all populated automatically.

file
++++

Filename of the output simulated file (e.g. jw42424024002_01101_00001_nrcb5_uncal.fits). If the linearized ramp is requested as output in the :ref:`datatype<datatype>` field, it will be saved with ‘uncal’ replaced with ‘linear’ in the filename or if ‘uncal’ is not present, ‘linear’ will simply be appended to the filename.  If the raw ramp is requested as output, the given filename will be used with no changes.

We recommend using filenames that end in 'uncal.fits' in order to be consistent with `JWST file naming conventions <https://jwst-docs.stsci.edu/display/JDAT/File+Naming+Conventions+and+Data+Products>`_. The filename is constructed from various pieces of information, including the program ID and visit number. If you wish to use this convention for the output filenames, the easiest way to accomplish this is to :ref:`create the yaml files from an APT file <from_apt>`, in which case the filenames will be generated automatically.

directory
+++++++++

The directory into which the output simulated data will be placed.

.. _datatype:

datatype
++++++++

List of the data format(s) of the output files. Options include:
“linear”, where the output files will contain linearized signals with the superbias and reference pixel signals removed. Bad pixels will also be flagged if a bad pixel file is specified. These files are ready to be run through the jump detection and ramp fitting steps of the JWST calibration pipeline. “raw”, where the output files will be in an uncalibrated state. These files are ready to be run through the entirety of the calibration pipeline, beginning with `calwebb_detector1 <https://jwst-pipeline.readthedocs.io/en/stable/jwst/pipeline/description.html#pipelines>`_.
“linear,raw”, where both the raw and linearized versions of the output files will be saved.

format
++++++

Format of the output file. Currently, only ‘DMS’ is supported, indicating that the fits file format, as well as header keywords, match those expected by the JWST calibration pipeline.

save_intermediates
++++++++++++++++++

True/False.  If True, intermediate products are saved to disk. These products are listed in the table below.

+------------+-----------------------------------------+----------------------------------------------------+
| Module     |  Suffix Appended to Output Filename     | Description                                        |
+============+=========================================+====================================================+
| Seed Image | _pointsources.list                      | Ascii file listing point source x,y                |
| Generator  |                                         | and RA, Dec positions as well as magnitude         |
|            |                                         | and count rate.                                    |
|            +-----------------------------------------+----------------------------------------------------+
|            | _galaxySources.list                     | Ascii file listing galaxy source x,y               |
|            |                                         | and RA, Dec positions, morphology parameters,      |
|            |                                         | magnitudes, and count rates.                       |
|            +-----------------------------------------+----------------------------------------------------+
|            | _extendedsources.list                   | Ascii file listing extended source x,y and RA,     |
|            |                                         | Dec positions as well as magnitude and count rate. |
|            +-----------------------------------------+----------------------------------------------------+
|            | _pointSourceRateImage_elec_per_sec.fits | Count rate image containing only added point       |
|            |                                         | sources                                            |
|            +-----------------------------------------+----------------------------------------------------+
|            | _galaxyRateImage_elec_per_sec.fits      | Count rate image containing only added galaxies    |
|            +-----------------------------------------+----------------------------------------------------+
|            | _extendedObject_elec_per_sec.fits       | Count rate image containing only extended objects  |
|            +-----------------------------------------+----------------------------------------------------+
|            | _AddedSources_elec_per_sec.fits	       | Count rate image containing all added sources      |
+------------+-----------------------------------------+----------------------------------------------------+
| Observation| _doNonLin_accuracy.fits                 | Final accuracy map from the process where the      |
| Generator  |                                         | linearized simulated exposure was “unlinearized”   |
|            +-----------------------------------------+----------------------------------------------------+
|            | _xtalk_correction_image.fits            | Image of the crosstalk signal added to the exposure|
|            +-----------------------------------------+----------------------------------------------------+
|            | _cosmicrays.list                        | Ascii file containing location and magnitude of    |
|            |                                         | added cosmic rays                                  |
+------------+-----------------------------------------+----------------------------------------------------+

grism_source_image
++++++++++++++++++

True/False. If True, the size of the output image is enlarged from the requested array size by a multiplicative factor in the x and y dimensions. For NIRCam this factor is √2, while it NIRISS it is 1.134. This extra area is required if the image is passed to the grism disperser software. In this case, the disperser software is able to include sources which fall just outside the nominal field of view but whose dispersed spectra fall into the nominal field of view.

program_number
++++++++++++++

The proposal ID number. This is placed in the header of the output file in order to match the contents of real observation files.

title
+++++

The title of the proposal. This placed in the header of the output file in order to match the contents of real observation files.

PI_Name
+++++++

Name of the proposal PI. This is placed in the header of the output file in order to match the contents of real observation files.

proposal_category
+++++++++++++++++

Proposal category (e.g. GO, GTO). This is placed in the header of the output file in order to match the contents of real observation files.

science_category
++++++++++++++++

Science category of the proposal, as defined in the APT file. This is placed in the header of the output file in order to match the contents of real observation files.

observation_number
++++++++++++++++++

The observation number containing the output exposure, as defined in the program’s APT file. This is placed in the header of the output file in order to match the contents of real observation files.

observation_label
+++++++++++++++++
The observation label in the APT file under which the output exposure appears. This is placed in the header of the output file in order to match the contents of real observation files.

visit_number
++++++++++++

The visit number, as defined in the APT file, within which the output exposure appears. This is placed in the header of the output file in order to match the contents of real observation files.

visit_group
+++++++++++

The visit group, as defined in the APT file, within which the output exposure appears. This is placed in the header of the output file in order to match the contents of real observation files.

visit_id
++++++++

The visit identifier of the exposure. This can be created by combining the program ID, visit number, and observation number. This is placed in the header of the output file in order to match the contents of real observation files.

sequence_id
+++++++++++

The parallel sequence identifier denotes whether the data were acquired during parallel observations, and with which instrument. Set to 0 for non-parallel observations, 1 for a parallel sequence using the primary instrument, or 2-5 for one of the non-prime instruments.

activity_id
+++++++++++

The activity identifier of the exposure is a base-36 number that is unique to each exposure in a proposal. This is placed in the header of the output file in order to match the contents of real observation files.

exposure_number
+++++++++++++++

A five-character number used to identify the exposure within the current activity.

obs_id
++++++

The observation ID is constructed from several of the other parameters. OBS_ID = 'V' + program_number + observation_id + visit_id + 'P' + parallel-program number + parallel-observation number + visit_group + parallel sequence identifier + activity_identifier.

date_obs
++++++++

UTC date of the start of the exposure with format yyyy-mm-dd.

time_obs
++++++++

UTC time of the start of the exposure with format hh:mm:ss.ssssss.

obs_template
++++++++++++

The name of the observation template used for the exposure (e.g. NIRCam Imaging, NIRCam Time Series)

primary_dither_type
+++++++++++++++++++

Name of the primary dither pattern in use when the data were obtained. For details, see the documentation pages on dither patterns for `NIRCam <https://jwst-docs.stsci.edu/display/JTI/NIRCam+Primary+Dithers>`_, and `NIRISS <https://jwst-docs.stsci.edu/display/JTI/NIRISS+Dithers>`_. (e.g. INTRAMODULEX, INTRASCA).

total_primary_dither_positions
++++++++++++++++++++++++++++++

Total number of primary dither positions in the observation.

primary_dither_position
+++++++++++++++++++++++

Primary dither position number of the exposure being simulated.

subpix_dither_type
++++++++++++++++++

Name of the subpixel dither pattern used for these data. Details on subpixel dither patterns can be found on the `NIRCam subpixel dither patterns page <https://jwst-docs.stsci.edu/display/JTI/NIRCam+Subpixel+Dithers>`_.

total_subpix_dither_positions
+++++++++++++++++++++++++++++

Total number of subpixel dither positions for this observation.

subpix_dither_position
++++++++++++++++++++++

The subpixel dither position number corresponding to the current exposure.

xoffset
+++++++

Offset in the x direction, in arcseconds, of the pointing used for the current exposure relative to the starting position of the dither pattern. This is used to populate header values only. It is not used to determine the pointing when creating the simulated data.

yoffset
+++++++

Offset in the y direction, in arcseconds, of the pointing used for the current exposure relative to the starting position of the dither pattern. This is used to populate header values only. It is not used to determine the pointing when creating the simulated data.



