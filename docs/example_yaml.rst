.. _example_yaml:

Example yaml Input File
=======================

Below is an example yaml input file for *Mirage*. The yaml file used as the primary input to mirage contains a large number of telescope and instrument settings, reference files, and catalogs. The inputs are grouped by section:

1. :ref:`Inst section <Inst>`
2. :ref:`Readout section <Readout>`
3. :ref:`Reffiles section <Reffiles>`
4. :ref:`nonlin section <nonlin>`
5. :ref:`cosmicRay section <cosmicRay>`
6. :ref:`simSignals section <simSignals>`
7. :ref:`Telescope section <Telescope>`
8. :ref:`newRamp section <newRamp>`
9. :ref:`Output section <Output>`

.. For more information on the individual input paramters, see the :ref:`Input Yaml Parameters <input_yaml_file_parameters>` page.


.. parsed-literal::

	Inst_:
	  instrument_: NIRCam          #Instrument name
	  mode_: imaging               #Observation mode (e.g. imaging, WFSS, moving_target)
	  use_JWST_pipeline_: False    #Use pipeline in data transformations

	Readout_:
	  readpatt_: DEEP8         #Readout pattern (RAPID, BRIGHT2, etc) overrides nframe,nskip unless it is not recognized
	  ngroup_: 6               #Number of groups in integration
	  nint_: 1                 #Number of integrations per exposure
	  array_name_: NRCB5_FULL  #Name of array (FULL, SUB160, SUB64P, etc)
	  filter_: F250M           #Filter of simulated data (F090W, F322W2, etc)
	  pupil_: CLEAR            #Pupil element for simulated data (CLEAR, GRISMC, etc)

	Reffiles_:                   #Set to None or leave blank if you wish to skip that step
	  dark_: None                #Dark current integration used as the base
	  linearized_darkfile_: $MIRAGE_DATA/nircam/darks/linearized/B5/Linearized_Dark_and_SBRefpix_NRCNRCBLONG-DARK-60090141241_1_490_SE_2016-01-09T02h46m50_uncal.fits # Linearized dark ramp to use as input. Supercedes dark above
	  badpixmask_: crds          # If linearized dark is used, populate output DQ extensions using this file
	  superbias_: crds           #Superbias file. Set to None or leave blank if not using
	  linearity_: crds           #linearity correction coefficients
	  saturation_: crds          #well depth reference files
	  gain_: crds                #Gain map
	  pixelflat_: None
	  illumflat_: None           #Illumination flat field file
	  astrometric_: crds         #Astrometric distortion file (asdf)
	  ipc_: crds                 #File containing IPC kernel to apply
	  invertIPC_: True           #Invert the IPC kernel before the convolution. True or False. Use True if the kernel is designed for the removal of IPC effects, like the JWST reference files are.
	  occult_: None              #Occulting spots correction image
	  pixelAreaMap_: crds        #Pixel area map for the detector. Used to introduce distortion into the output ramp.
	  subarray_defs_:   config   #File that contains a list of all possible subarray names and coordinates
	  readpattdefs_:    config   #File that contains a list of all possible readout pattern names and associated NFRAME/NSKIP values
	  crosstalk_:       config   #File containing crosstalk coefficients
	  filtpupilcombo_:  config   #File that lists the filter wheel element / pupil wheel element combinations. Used only in writing output file
	  filter_wheel_positions_: config  #File that lists the filter wheel element / pupil wheel element combinations. Used only in writing output file
	  flux_cal_:        config   #File that lists flux conversion factor and pivot wavelength for each filter. Only used when making direct image outputs to be fed into the grism disperser code.

	nonlin_:
	  limit_: 60000.0        #Upper singal limit to which nonlinearity is applied (ADU)
	  accuracy_: 0.000001    #Non-linearity accuracy threshold
	  maxiter_: 10           #Maximum number of iterations to use when applying non-linearity
	  robberto_:  False      #Use Massimo Robberto type non-linearity coefficients

	cosmicRay_:
	  path_: $MIRAGE_DATA/nircam/cosmic_ray_library/    #Path to CR library
	  library_: SUNMIN    								#Type of cosmic rayenvironment (SUNMAX, SUNMIN, FLARE)
	  scale_: 1.5     									#Cosmic ray rate scaling factor
	  suffix_: IPC_NIRCam_B5    					    #Suffix of library file names
	  seed_: 2956411739      							#Seed for random number generator

	simSignals_:
	  pointsource_: my_point_sources.cat               #File containing a list of point sources to add (x,y locations and magnitudes)
	  psfpath_: $MIRAGE_DATA/nircam/gridded_psf_library/   #Path to PSF library
	  gridded_psf_library_row_padding_: 4              # Number of outer rows and columns to avoid when evaluating library. RECOMMEND 4.
  	  psf_wing_threshold_file_: config                 # File defining PSF sizes versus magnitude
  	  add_psf_wings_: True                             # Whether or not to place the core of the psf from the gridded library into an image of the wings before adding.
	  psfwfe_: predicted                               #PSF WFE value ("predicted" or "requirements")
	  psfwfegroup_: 0                                  #WFE realization group (0 to 4)
	  galaxyListFile_: my_galaxies_catalog.list
	  extended_: None                                 #Extended emission count rate image file name
	  extendedscale_: 1.0                             #Scaling factor for extended emission image
	  extendedCenter_: 1024,1024                      #x,y pixel location at which to place the extended image if it is smaller than the output array size
	  PSFConvolveExtended_: True                      #Convolve the extended image with the PSF before adding to the output image (True or False)
	  movingTargetList_: None                         #Name of file containing a list of point source moving targets (e.g. KBOs, asteroids) to add.
	  movingTargetSersic_: None                       #ascii file containing a list of 2D sersic profiles to have moving through the field
	  movingTargetExtended_: None                     #ascii file containing a list of stamp images to add as moving targets (planets, moons, etc)
	  movingTargetConvolveExtended_: True             #convolve the extended moving targets with PSF before adding.
	  movingTargetToTrack_: None                      #File containing a single moving target which JWST will track during observation (e.g. a planet, moon, KBO, asteroid)	This file will only be used if mode is set to "moving_target"
	  zodiacal_:  None                                #Zodiacal light count rate image file
	  zodiscale_:  1.0                                #Zodi scaling factor
	  scattered_:  None                               #Scattered light count rate image file
	  scatteredscale_: 1.0                            #Scattered light scaling factor
	  bkgdrate_: medium                               #Constant background count rate (ADU/sec/pixel in an undispersed image) or "high","medium","low" similar to what is used in the ETC
	  poissonseed_: 2012872553                        #Random number generator seed for Poisson simulation)
	  photonyield_: True                              #Apply photon yield in simulation
	  pymethod_: True                                 #Use double Poisson simulation for photon yield
	  expand_catalog_for_segments_: False             # Expand catalog for 18 segments and use distinct PSFs
	  use_dateobs_for_background_: False              # Use date_obs value to determine background. If False, bkgdrate is used.

	Telescope_:
	  ra_: 53.1                     #RA of simulated pointing
	  dec_: -27.8                   #Dec of simulated pointing
	  rotation_: 0.0                #y axis rotation (degrees E of N)
	  tracking_: sidereal           #sidereal or non-sidereal

	newRamp_:
	  dq_configfile_: config          #config file used by JWST pipeline
	  sat_configfile_: config         #config file used by JWST pipeline
	  superbias_configfile_: config   #config file used by JWST pipeline
	  refpix_configfile_: config      #config file used by JWST pipeline
	  linear_configfile_: config      #config file used by JWST pipeline

	Output_:
	  file_: jw42424024002_01101_00001_nrcb5_uncal.fits   # Output filename
	  directory_: ./                                # Directory in which to place output files
	  datatype_: linear,raw                         # Type of data to save. 'linear' for linearized ramp. 'raw' for raw ramp. 'linear,raw' for both
	  format_: DMS                                  # Output file format Options: DMS, SSR(not yet implemented)
	  save_intermediates_: False                    # Save intermediate products separately (point source image, etc)
	  grism_source_image_: False                    # Create an image to be dispersed?
	  unsigned_: True                               # Output unsigned integers? (0-65535 if true. -32768 to 32768 if false)
	  dmsOrient_: True                              # Output in DMS orientation (vs. fitswriter orientation).
	  program_number_: 42424                        # Program Number
	  title_: Supernovae and Black Holes Near Hyperspatial Bypasses   #Program title
	  PI_Name_: Doug Adams                          # Proposal PI Name
	  Proposal_category_: GO                        # Proposal category
	  Science_category_: Cosmology                  # Science category
	  target_name_: TARG1                           # Name of target
	  target_ra_: 53.1001                           # RA of the target, from APT file.
	  target_dec_: -27.799                          # Dec of the target, from APT file.
	  observation_number_: '002'                    # Observation Number
	  observation_label_: Obs2                      # User-generated observation Label
	  visit_number_: '024'                          # Visit Number
	  visit_group_: '01'                            # Visit Group
	  visit_id_: '42424024002'                      # Visit ID
	  sequence_id_: '1'                             # Sequence ID
	  activity_id_: '01'                            # Activity ID. Increment with each exposure.
	  exposure_number_: '00001'                     # Exposure Number
	  obs_id_: 'V42424024002P0000000001101'         # Observation ID number
	  date_obs_: '2019-10-15'                       # Date of observation
	  time_obs_: '06:29:11.852'                     # Time of observation
	  obs_template_: 'NIRCam Imaging'               # Observation template
	  primary_dither_type_: NONE                    # Primary dither pattern name
	  total_primary_dither_positions_: 1            # Total number of primary dither positions
	  primary_dither_position_: 1                   # Primary dither position number
	  subpix_dither_type_: 2-POINT-MEDIUM-WITH-NIRISS  #Subpixel dither pattern name
	  total_subpix_dither_positions_: 2             # Total number of subpixel dither positions
	  subpix_dither_position_: 2                    # Subpixel dither position number
	  xoffset_: 344.284                             # Dither pointing offset in x (arcsec)
	  yoffset_: 466.768                             # Dither pointing offset in y (arcsec)


.. _inst:

Instrument secton
-----------------

This section of the input yaml file contains information about the instrument being simulated.

.. _instrument:

Instrument Name
+++++++++++++++

*Inst:instrument*

The name of the JWST instrument to be simulated. The simulator will only function if ‘NIRCam’, ‘NIRISS’, or ‘FGS’ is placed in this field.

.. _mode:

Observing mode
++++++++++++++

*Inst:mode*

The observing mode to be simulated. There are three valid options for this field. “imaging” will create imaging data, “wfss” will produce wide field slitless spectroscopic data. The other accepted input is "ami" when simulating NIRISS, although this mode is functionally identical to the use of "imaging".


.. _use_JWST_pipeline:

Create data using JWST pipeline
+++++++++++++++++++++++++++++++

*Inst:use_JWST_pipeline*

True/False. Set to False if you wish to proceed without using any JWST pipeline functions. In this case, the input dark current exposure must already be linearized, as the pipeline is used for the linearization process. True is recommneded.

.. _Readout:

Readout section
---------------

This section of the yaml file contains inputs describing the details of the exposure, including the readout pattern, filter, subarray, etc to use.


.. _readpatt:

Readout pattern
+++++++++++++++

*Readout:readpatt*

This is the name of the readout timing pattern used for the output simulated exposure. Examples for NIRCam include RAPID, BRIGHT1, BRIGHT2, and DEEP8. Each pattern averages and skips a predefined number of frames when constructing each group of an integration. The list of possible readout patterns and their definitions is provided by an ascii file specified in the **readpattdefs** parameter in the **Reffiles** section of the input file. A more detailed description of readout patterns is given in the detector readout pages for `NIRCam <https://jwst-docs.stsci.edu/display/JTI/JWST+Field+of+View>`_,  `NIRISS <https://jwst-docs.stsci.edu/display/JTI/JWST+Field+of+View>`_, and `FGS <https://jwst-docs.stsci.edu/display/JTI/JWST+Field+of+View>`_.

.. _ngroup:

Number of groups per integration
++++++++++++++++++++++++++++++++

*Readout:ngroup*


This parameter lists the number of groups comprising each output integration.


.. _nint:

Number of integrations per exposure
+++++++++++++++++++++++++++++++++++

*Readout:nint*

The number of integrations in the output exposure. Each integration is composed of **ngroup** groups. Note that currently, any observation containing a moving target (non-sidereal observation with trailed sidereal objects, or vice versa) cannot have an nint value greater than 1. **(IS THIS STILL TRUE?)**

.. _resets_bet_ints:

Number of detector resets between integrations
++++++++++++++++++++++++++++++++++++++++++++++

*Readout:resets_bet_ints*

The number of detector resets between integrations within a single exposure. For all instruments, this should be set to 1.

.. _array_name:

Array Name
++++++++++

*Readout:array_name*

This is the name of the aperture used for the simulated data. Generally, this is composed of the name of the detector combined with the name of the subarray used. For example, a full frame observation using NIRCam's A1 detector has an **array_name** of 'NRCA1_FULL', while a full frame NIRISS observation will have an array_name of ‘NIS_CEN’. The list of possible array_name values are given in the **subarray_defs** input file described below. The **array_name** is used to identify several other characteristics of the simulated data, including the detector to use, as well as the proper array dimensions and location on the detector.

.. _filter:

Filter
++++++

*Readout:filter*

The name of the filter wheel element to use for the simulated data. (e.g. F444W). The filter is used when scaling astronomical sources from the requested brightness in magnitudes to counts on the detector. For NIRCam simulations, the filter name is also used to determine whether the simulated data are to be produced using a shortwave or longwave detector. Lists of instrument filters can be found on the `NIRCam <https://jwst-docs.stsci.edu/display/JTI/JWST+Field+of+View>`_,  `NIRISS <https://jwst-docs.stsci.edu/display/JTI/JWST+Field+of+View>`_, and `FGS <https://jwst-docs.stsci.edu/display/JTI/JWST+Field+of+View>`_ filter pages.

.. _pupil:

Pupil
+++++

*Readout:pupil*

The name of the pupil wheel element to use for the simulated data. Some filters for both NIRCam and NIRISS reside in their respective pupil wheels. Therefore this entry is checked when deciding upon scaling factors for simulated sources. Pupil wheel elements are desribed in the `NIRCam <https://jwst-docs.stsci.edu/display/JTI/JWST+Field+of+View>`_,  `NIRISS <https://jwst-docs.stsci.edu/display/JTI/JWST+Field+of+View>`_, and `FGS <https://jwst-docs.stsci.edu/display/JTI/JWST+Field+of+View>`_ pupil wheel pages.

.. _Reffiles:

Reffiles section
----------------

This section of the input file lists the various reference files needed for the various steps of the simulator to run.

.. _dark:

Dark current exposure
+++++++++++++++++++++

*Reffiles:dark*

The name of the raw dark current file that will be used as the basis for the simulated exposure. This file must be in raw format, such that no JWST calibration pipeline steps have been applied to the data. If an already-linearized dark current integration is to be used, that file name should be placed in the **linearized_darkfile** field below. Note that the **linearized_darkfile** entry will take precedence. Only if that is set to __None__ will the file listed in this field be used.

The dark current integration must have a readout pattern of either RAPID/NISRAPID/FGSRAPID or a value identical to that of the integration to be simulated. RAPID/NISRAPID/FGSRAPID data keep every readout frame with no averaging. From this, any other readout pattern can be simulated by averaging and skipping the appropriate frames. Other readout patterns cannot be translated in this way as their data are already averaged or missing some frames. However if simulating, for example a BRIGHT2 integration, then the input dark current integration can be a BRIGHT2 integration, as no translation is necessary in this case.

If a translation between RAPID and another readout pattern is necessary, then frames will be averaged/skipped as necessary. If the input dark current integration does not contain enough frames to be translated into the requested number of output groups, then the script creates enough additional dark current frames to make up the difference. These additional frames are created by making a copy of an appropriate number of existing initial dark current frames, and adding their signals to that in the final dark current frame. Note that this can lead to apparent double cosmic rays in pixels where a cosmic ray appeared in the dark current integration.

.. hint::
	This input can only be used if **use_JWST_pipeline** is set to True.

.. hint::
	The collection of reference files associated with Mirage contains a small library of raw dark current exposures that can be used.

.. _linearized_darkfile:

Linearized dark current exposure
++++++++++++++++++++++++++++++++

*Reffiles:linearized_darkfile*

The name of a linearized dark current integration to use as input for the simulated data. This file should contain a dark integration that has been processed through the superbias subtraction, reference pixel subtraction, and linearity steps of the JWST calibration pipeline. The resulting linearized signal must be saved in an extension with the name 'SCI'. Also, the subtracted signal from the superbias and reference pixels must be saved in an extension called 'SBANDREFPIX'. This output will be produced and saved for a given dark current file by Mirage.

Using this input rather than the uncalibrated dark above can save significant computing time, especially in the case of creating many output exposures.

.. hint::
	This input can be used for **use_JWST_pipeline** set to True or False.

.. hint::
	The collection of :ref:`reference files <reference_files>` associated with Mirage contains a small library of linearized dark current products that can be used.

.. _badpixmask:

Bad pixel mask
++++++++++++++

*Reffiles:badpixmask*

If a linearized dark current file is to be used and a linearized output file is requested, this optional bad pixel mask can be used to populate the data quality array in the output simulated data file. The file must be in the `format for JWST bad pixel masks <https://jwst-pipeline.readthedocs.io/en/stable/jwst/dq_init/reference_files.html>`_ that is used by the JWST calibration pipeline.

.. hint::
	Setting this entry equal to 'crds' will cause Mirage to query the Calibration Reference Database System (CRDS) for the appropriate file, and download that file if it is not already present in your CRDS cache.

.. _superbias:

Superbias
+++++++++

*Reffiles:superbias*

The superbias reference file for the detector of the simulation. This file must match the `format of the JWST pipeline superbias reference file <https://jwst-pipeline.readthedocs.io/en/stable/jwst/superbias/reference_files.html>`_. If the input dark current integration is a raw file then this superbias file is used to subtract the superbias from the dark. If the input dark is already linearized, this superbias file is not used.

.. hint::
	Setting this entry equal to 'crds' will cause Mirage to query the Calibration Reference Database System (CRDS) for the appropriate file, and download that file if it is not already present in your CRDS cache.

.. _linearity:

Linearity correction coefficients
+++++++++++++++++++++++++++++++++

*Reffiles:linearity*

Name of the reference file containing the linearity correction coefficients. This file must be in the `format expected by the JWST calibration pipeline <https://jwst-pipeline.readthedocs.io/en/stable/jwst/linearity/reference_files.html>`_. If the input dark current integration is raw, the coefficients contained in this file are used to linearize the dark current after subtracting the superbias and reference pixel signal. These coefficients are also used to "unlinearize" the final simulated exposure if a raw simulated observation is requested.

In addition, the coefficients in this file are used to linearize the values in the saturation reference file, such that saturated signals in the linear simulated exposure can be found.

.. hint::
	Setting this entry equal to 'crds' will cause Mirage to query the Calibration Reference Database System (CRDS) for the appropriate file, and download that file if it is not already present in your CRDS cache.

.. _saturation:

Saturation
++++++++++

*Reffiles:saturaiton*

Name of the reference file containing a map of the saturation signal level for all pixels. If the input dark current integration is raw, this file is used by the calibration pipeline to flag saturated pixels in the dark current integration prior to linearizing. The `format of this file <https://jwst-pipeline.readthedocs.io/en/stable/jwst/saturation/reference_files.html>`_ must match that used in the saturation flagging step of the JWST calibration pipeline.

This saturation map, after being linearized, is also used to search for saturated signal values in the combined dark current/simulated source exposure prior to unlinearizing.

.. hint::
	Setting this entry equal to 'crds' will cause Mirage to query the Calibration Reference Database System (CRDS) for the appropriate file, and download that file if it is not already present in your CRDS cache.

.. _gain:

Gain
++++

*Reffiles:gain*

Name of the file containing the gain map appropriate for the detector being used. The gain is used to translate the cosmic rays, which are in units of electrons, to units of ADU prior to adding them to the simulated data. The `format of the gain file <https://jwst-pipeline.readthedocs.io/en/stable/jwst/references_general/gain_reffile.html#gain-reffile>`_ must match that used by the JWST calibration pipeline.

.. hint::
	Setting this entry equal to 'crds' will cause Mirage to query the Calibration Reference Database System (CRDS) for the appropriate file, and download that file if it is not already present in your CRDS cache.

.. _pixelflat:

Pixel-to-pixel flat field image
+++++++++++++++++++++++++++++++

*Reffiles:pixelflat*

Name of the pixel flat file to use. Once the simulated integration is created, the result is multiplied by the pixel flat. This is done to un-flatten the image.


.. _illumflat:

Illumination flat (L-flat)
++++++++++++++++++++++++++

*Reffiles:illumflat*

Name of the illumination flat to use. Once the simulated integration is created, the result is multiplied by the illumination flat.


.. _astrometric:

Astrometric distortion file
+++++++++++++++++++++++++++

*Reffiles:astrometric*

Name of the astrometric distortion reference file to use for including the effects of distortion in the simulated data.  This file is used to translate input source locations between RA and Dec coordinates and pixel x and y coordinates, and vice versa. This file must be in `asdf format and match that expected by the calibration pipeline <https://jwst-pipeline.readthedocs.io/en/stable/jwst/references_general/distortion_reffile.html#distortion-reference-file>`_.

.. hint::
	Setting this entry equal to 'crds' will cause Mirage to query the Calibration Reference Database System (CRDS) for the appropriate file, and download that file if it is not already present in your CRDS cache.

.. _ipc:

Interpixel capacitance (IPC)
++++++++++++++++++++++++++++

*Reffiles:ipc*

File containing the interpixel capacitance (IPC) kernel to apply to the simulated data in order to introduce IPC effects. After all simulated objects have been added to a count rate image, the image is convolved with the IPC kernel. The IPC file must be a fits file with the IPC kernel located in the first (rather than 0th) extension. Typical JWST IPC reference file kernels are a 3x3 array, but Mirage supports kernels of any odd-numbered size, as well as 4-dimensional kernels, where there is a separate 2-dimensional kernel for each pixel. In order to introduce, rather than remove, IPC effects, the kernel must be normalized and have a value in the central pixel which is less than 1.0. This is the inverse of the kernel used in the JWST calibration pipeline IPC removal step, where the central pixel has a value greater than 1.0, and negative values in surrounding pixels. For the simulator, the user can specify a `JWST calibration pipeline-formatted kernel file <https://jwst-pipeline.readthedocs.io/en/stable/jwst/ipc/reference_files.html>`_, and then set the **invertIPC** flag below to True, in which case the kernel will be inverted before using.

.. hint::
	Setting this entry equal to 'crds' will cause Mirage to query the Calibration Reference Database System (CRDS) for the appropriate file, and download that file if it is not already present in your CRDS cache.

.. _invertIPC:

Invert IPC
++++++++++

*Reffiles:invertIPC*

If set to True, the IPC kernel supplied through the ipc entry is inverted before convolving with the signal rate image. JWST IPC kernel reference files contain the kernel necessary to remove IPC from the data. Therefore these kernels must be inverted before they can add IPC effects to the data in the simulator.

.. _occult:

Occulting spot image
++++++++++++++++++++

*Reffiles:occult*

This feature is not yet supported and should be set to **None**.

.. _pixelAreaMap:

Pixel area map
++++++++++++++

*Reffiles:pixelAreaMap*

Fits file containing the pixel area map for the detector to be simulated. If provided, the pixel area map is multiplied into the seed image at a point when the seed image contains only extended sources. Point sources have the pixel area map applied to them at the time the PSF libraries were created via `webbpsf <https://webbpsf.readthedocs.io/en/stable/>`_. The pixel area map file must be in the format of the `JWST pixel area map reference file <https://jwst-pipeline.readthedocs.io/en/stable/jwst/photom/reference_files.html#area-reference-file>`_.

.. hint::
	Setting this entry equal to 'crds' will cause Mirage to query the Calibration Reference Database System (CRDS) for the appropriate file, and download that file if it is not already present in your CRDS cache.

.. _subarray_defs:

Subarray definition file
++++++++++++++++++++++++

Reffiles:subarray_defs*

Name of a whitespace-delimited ascii file that lists all of the possible supported subarray apertures. This file is provided with the MIRAGE repository, in the `config <https://github.com/spacetelescope/mirage/tree/master/mirage/config>`_ subdirectory.

.. hint::
	To use the subarray definition files packaged with Mirage, set this to **config** in the input yaml file. This is the default when creating yaml files from an APT file using the :ref:`yaml generator <yaml_generator>`

For each subarray, the file must list the full aperture name (e.g. NRCA1_FULL) as well as the corresponding name used in proposal planning (e.g. FULL), as well as the number of amplifiers used to read out each aperture.

.. _readpattdefs:

Readout pattern definition file
+++++++++++++++++++++++++++++++

*Reffiles:readpattdefs*

Ascii file which gives the definitions of the possible readout patterns for the instrument. For each readout pattern, the number of frames averaged to create each group (nframe) and the number of frames skipped beteren each group (nskip) must be specified, as well as the maximum number of allowed groups. For a given readout pattern the simulator will search the entries in this file in order to determine the proper nframe and nskip values to use. The current lists of acceptable NIRCam and NIRISS readout patterns are given on the NIRCam  and NIRISS  detector readouts webpages. These files for all instruments are provided with the MIRAGE repository, in the `config <https://github.com/spacetelescope/mirage/tree/master/mirage/config>`_ subdirectory.

.. hint::
	To use the readout pattern definition files packaged with Mirage, set this to **config** in the input yaml file. This is the default when creating yaml files from an APT file using the :ref:`yaml generator <yaml_generator>`

.. _crosstalk:

Crosstalk
+++++++++

*Reffiles:crosstalk*

Ascii file containing crosstalk coefficients. Crosstalk is only applied to data read out through more than one amplifer. The file contains one row for each detector. Each row contains all of the coefficients necessary to fully describe crosstalk. This file is contained in the MIRAGE repository, in the `config <https://github.com/spacetelescope/mirage/tree/master/mirage/config>`_ subdirectory.

.. hint::
	To use the crosstalk coefficient files packaged with Mirage, set this to **config** in the input yaml file. This is the default when creating yaml files from an APT file using the :ref:`yaml generator <yaml_generator>`

.. _filtpupilcombo:

Allowed filter/pupil combinations
+++++++++++++++++++++++++++++++++

*Reffiles:filtpupilcombo*

Name of an ascii file containing a list of the filter and pupil wheel elements in place when requesting simulated data for a given filter. This information is used to apply the appropriate conversion between magnitudes and counts when reading in source catalogs. This flux calibration is also added to the header of the seed image, as it is used when seed images are dispersed during the simulation of WFSS data. This file is present in the `config <https://github.com/spacetelescope/mirage/tree/master/mirage/config>`_ subdirectory of the MIRAGE repository.

.. hint::
	To use the filter and pupil wheel definition files packaged with Mirage, set this to **config** in the input yaml file. This is the default when creating yaml files from an APT file using the :ref:`yaml generator <yaml_generator>`

.. _filter_wheel_positions:

Filter/Pupil wheel resolver positions for each optical element
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

*Reffiles:filter_wheel_positions*

Name of an ascii file containing a list of all filter wheel and pupil wheel elements, along with the nominal wheel resolver positions for each. These values are in degrees. This information is passed directly to the header keywords FWCPOS and PWCPOS in the simulated data FITS files. This information is needed to compute the dispersion solution for NIRISS WFSS. Currently the header keywords are only populated for NIRISS observations.

.. hint::
	To use the filter and pupil wheel position files packaged with Mirage, set this to **config** in the input yaml file. This is the default when creating yaml files from an APT file using the :ref:`yaml generator <yaml_generator>`

.. _flux_cal:

Flux calibration
++++++++++++++++

*Reffiles:flux_cal*

Ascii file that lists flux conversion factors and the pivot wavelength associated with each filter. Conversion factors include ABMAG, STMAG, and VEGAMAG to counts per second, as well as FLAM (erg s :sup:`-1` cm :sup:`-2` Å :sup:`-1` and FNU (erg s :sup:`-1` cm :sup:`-2` Hz :sup:`-1` to counts per second. This file is used when producing seed images to be fed into the grism disperser code, as well as for translating catalog sources from magnitudes to counts per second. This file is provided with the MIRAGE repository, in the `config <https://github.com/spacetelescope/mirage/tree/master/mirage/config>`_ subdirectory.

.. hint::
	To use the flux calibration files packaged with Mirage, set this to **config** in the input yaml file. This is the default when creating yaml files from an APT file using the :ref:`yaml generator <yaml_generator>`

.. _nonlin:

Nonlin section
--------------

The following input fields describe how non-linearity is treated in the input and simulated data.

.. _limit:

Limiting Signal
+++++++++++++++

*nonlin:limit*

Signal limit, in units of ADU, above which the linearity correction is not applied. Pixels with signals above this limit are considered saturated. This single value across the entire detector is only used if a :ref:`saturation reference file <saturation>` is not provided.

.. _accuracy:

Accuracy
++++++++

*nonlin:accuracy*

When introducing non-linearity back into the linear data, the Newton-Raphson method is used to essentially run the JWST calibration pipline’s linearity correction step in reverse. The value of this accuracy parameter is the threshold below which the solution is considered to have converged. For example, an accuracy threshold of 0.000001 means that the unlinearization is considered complete when the ratio of the signal values from one iteration to the next is less than 1.000001.

.. _maxiter:

Maximum number of iterations
++++++++++++++++++++++++++++

*nonlin:maxiter*

The maximum number of iterations of the Newton-Raphson method to use when introducing non-linearity back into the data before declaring failure. Default is 10.

.. _robberto:

Robberto
++++++++

*nonlin:robberto*

If set to False, the simulator assumes that the non-linearity correction function and coefficients match those used in the JWST calibration pipeline. If set to True, the script assumes an alternate linearity function, as defined in Robberto (`2010 <https://jwst.stsci.edu/files/live/sites/jwst/files/home/instrumentation/technical%20documents/JWST-STScI-002163.pdf>`_ , `2011 <https://jwst.stsci.edu/files/live/sites/jwst/files/home/instrumentation/technical%20documents/JWST-STScI-002344.pdf>`_). **Currently, no coefficients for the latter method exist, implying this parameter should be set to False.**

.. _cosmicRay:

Cosmic ray section
------------------

Input parameters in this section describe how cosmic rays are added to the simulated data.

.. _path:

Path to cosmic ray library
++++++++++++++++++++++++++

*cosmicRay:path*

Path of the location of the cosmic ray library to use. The code was developed around the cosmic ray library produced by Robberto (`2009 <https://jwst.stsci.edu/files/live/sites/jwst/files/home/instrumentation/technical%20documents/JWST-STScI-001928.pdf>`_). This library is included in the collection of `reference files <reference_files>`_ associated with Mirage. After extracting the library from the tar file, set this path to point to the top level directory of the cosmic ray library.

.. _library:

Library
+++++++

*cosmicRay:library*

Specification of which cosmic ray library to choose cosmic rays from. Options are SUNMIN, SUNMAX, FLARE, each of which assumes a different cosmic ray rate. Details on the three types of libraries are given in Robberto (`2009 <https://jwst.stsci.edu/files/live/sites/jwst/files/home/instrumentation/technical%20documents/JWST-STScI-001928.pdf>`_).

.. _scale:

Scaling value for rate
++++++++++++++++++++++

*cosmicRay:scale*

Scaling factor to apply to the cosmic ray rate. For example, to simulate cosmic rays at a rate twice as high as that in SUNMIN, set library to SUNMIN and scale to 2.0

.. _suffix:

Suffix
++++++

*cosmicRay:suffix*

Filename suffix of the cosmic ray library files. The code was developed around files with the suffix of ‘IPC_NIRCam_XX’ where XX is the detector (e.g. B5) for NIRCam, ‘IPC_NIRISS_NIS’ for NIRISS, and ‘IPC_FGS_GUIDERy’ where y is 1 or 2, for FGS. These cosmic ray files are included in Mirage's `reference file collection <reference_files>`_. This field will be automatically populated with the correct suffix when creating yaml files using the :ref:`yaml generator <yaml_generator>`.

.. _seed:

Seed for random number generator
++++++++++++++++++++++++++++++++

*cosmicRay:seed*

Random number generator seed to use when selecting cosmic rays to add.

.. _simsignals:

simSignals section
------------------

This section of the input file describes how sources and other signals are added to the simulated data.

.. _pointsource:

Point source catalog file
+++++++++++++++++++++++++

*simSignals:pointsource*

Name of an ascii catalog file listing point sources to add to the simulated image. An example :ref:`point source <point_source>` catalog is provided on the :ref:`Catalogs page <catalogs>`.

.. _psfpath:

PSF library path
++++++++++++++++

*simSignals:psfpath*

Path name to the PSF library to be used for adding point sources to the data. The code was developed around a PSF library constructed using WebbPSF (Perrin, 2014). This PSF library is included in the collection of Mirage `reference files <reference_files>`_ . Once that package is downloaded and the data files extracted from the tar file, set this field to point to the top-level directory of the PSF library.

.. _gridded_psf_library_row_padding:

Gridded PSF Library Row Padding
+++++++++++++++++++++++++++++++

The number of outer rows and columns to crop when evaluating the PSF library. This is done to avoid edge effects that can sometimes be
present in the evaluated PSF. Recommended and default value is 4.

.. _psf_wing_threshold_file:

PSF Wing Threshold File
+++++++++++++++++++++++

Ascii file that defines the overall size of the PSF (in pixels) versus magnitude. Through this file, the user can tune the size of the PSFs in the
simulated data. If it is important for your science to see far out into the wings, you can enable that here. These files are located in the ``config``
directory of the repo. There is one file per instrument. The default value for this keyword is ``config``. In this case, Mirage will know to look
for the file in the ``config`` directory.

.. _add_psf_wings:

Add PSF Wings
+++++++++++++

Boolean value stating whether or not to place the core of the psf from the gridded library into an image of the wings before adding.


.. _psfwfe:

PSF library wavefront error
+++++++++++++++++++++++++++

*simSignals:psfwfe*

PSF wavefront error value to use when choosing PSF files from the PSF library. The current PSF libraries distributed with the Mirage `reference files <reference_files>`_ have two options for wavefront error: “predicted” and “requirements”. These two values represent the predicted in-flight wavefront errors, and the maximum allowed wavefront errors, respectively.

.. _psfwfegroup:

PSF realization number
++++++++++++++++++++++

*simSignals:psfwfegroup*

The current PSF library contains 5 different realizations for each filter/wavefront error-specified PSF. In this field, place the realization number to use. With 5 realizations present in the library, this field can have a value of 0 through 4.

.. _galaxyListFile:

Galaxy source catalog file
++++++++++++++++++++++++++

*simSignals:galaxyListFile*

Similar to the :ref:`pointsource <pointsource>` entry, this is an ascii catalog file containing a list of the galaxies to simulate in the data. See the :ref:`galaxies <galaxies>` entry on the :ref:`catalogs <catalogs>` page for an example of this file.

.. _extendedlist:

.. _extended:

Extended source catalog file
++++++++++++++++++++++++++++

*simSignals:extended*

Name of an ascii file containing a list of "extended images" to add to the simulated data. These are stamp image of sources, contained in small fits files. These stamp images are read in, scaled to the requested magnitude, and added to the seed image.  This is a way to add objects other than point sources or 2D Sersic profiles to the data. The :ref:`extended catalog <extended>` section of the :ref:`catalogs <catalogs>` page shows an example extended source catalog.

.. _extendedscale:

Extended source scaling factor
++++++++++++++++++++++++++++++

*simSignals:extendedScale*

Multiplicative factor by which to scale the data in the extended image file before adding to the simulated data. The extended image is multiplied by this factor **if the magnitude is set to None in the extended catalog file**.

.. _extendedCenter:

Extended source center location
+++++++++++++++++++++++++++++++

*simSignals:extendedCenter*

In the case where a single extended source is provided, this entry can be set to the (x,y) pixel location at which to place the center of the exteded image. This functionality is largely replaced by specifying the RA, Dec or x, y of the extended image in the :ref:`extended source catalog file <extended>`.

.. _PSFConvolveExtended:

Convolve extended sources with PSF
++++++++++++++++++++++++++++++++++

*simSignals:PSFConvolveExtended*

True/False. Convolve the extended image with the appropriate instrumental PSF prior to adding to the output image.

.. _movingTargetList:

Moving target source catalog file
+++++++++++++++++++++++++++++++++

*simSignals:movingTargetList*

Similar to the :ref:`point source <pointsource>` list file, this is a file containing a list of targets to treat as moving (non-sidereal) targets.  These sources will move through the field of view as the exposure progresses. This is the list to use if you wish to insert an asteroid or KBO that is moving through the field of view of your observation. See the :ref:`moving point source <moving_point_source>` section on the :ref:`Catalogs <catalogs>` page for an example.

.. _movingTargetSersic:

2D Sersic profile moving target catalog file
++++++++++++++++++++++++++++++++++++++++++++

*simSignals:movingTargetSersic*

Similar to the :ref:`galaxy target list file <galaxyListFile>`, this file contains a list of galaxies (2D Sersic profiles) to be used as moving targets. These sources will move through the background of the simulated data. This may be useful for inserting a resolved moon/asteroid into the scene. An example file is shown in the :ref:`Moving Sersic <moving_sersic>` section of the :ref:`Catalogs <catalogs>` page.

.. _movingTargetExtended:

Moving extended source catalog file
+++++++++++++++++++++++++++++++++++

*simSignals:movingTargetExtended*

Similar to the :ref:`extended <extended>` target list, this is an ascii file listing extended targets to move through the background of the image. A description and example of this file are shown in the :ref:`Moving Extended <moving_extended>` section of the :ref:`Catalogs <catalogs>` page.

.. _movingTargetConvolveExtended:

Convolve moving extended targets with PSF
+++++++++++++++++++++++++++++++++++++++++

*simSignals:movingTargetConvolveExtended*

Set this input to True if you wish to convolve the images listed in **movingTargetExtended** with the instrumental PSF prior to adding them to the simulated data.

.. _movingTargetToTrack:

Tracked non-sidereal target catalog file
++++++++++++++++++++++++++++++++++++++++

*simSignals:movingTargetToTrack*

This ascii catalog file is used for what are traditionally (in HST jargon) called 'moving targets'.  Targets listed in this file are treated as non-sidereal targets that JWST will track during the simulated observation. In this case, the target listed in this file will appear static in the output data, but all other sources (e.g. those listed in :ref:`pointSource <pointsource>`, :ref:`galaxyListFile <galaxyListFile>`, and :ref:`extended <extended>`) will all appear trailed through the data. A description and example of the file are shown in the :ref:`Non-sidereal Source <nonsidereal>` section on the :ref:`Catalogs <catalogs>` page.

.. _zodiacal:

Zodiacal light
++++++++++++++

*simSignals:zodiacal*

This keyword has been depricated in favor of obtaining the zodiacal light from the `JWST backgrounds package <https://github.com/spacetelescope/jwst_backgrounds>`_.

Name of a file containing a 2 dimensional count rate image of zodiacal light. This file is read in, scaled by the :ref:`zodiscale <zodiscale>` value, and added to the seed image. Leave as None to skip this step. The behaviors of this step and the scattered step below are very basic, and identical. There are no requirements on what the count rate images in these files must look like.

.. tip::

    Note that the :ref:`bkgdrate <bkgdrate>` input parameter, when set to “high”, “medium”, or “low”, will return a background rate image that includes the contribution from zodiacal light, in which case this step should be set to None.


.. _zodiscale:

Scaling factor for zodiacal light image
+++++++++++++++++++++++++++++++++++++++

*simSignals:zodiscale*

Scaling factor to multiply the :ref:`zodiacal light count rate image <zodiacal>` by prior to adding to the output data.

.. _scattered:

Scattered light image
+++++++++++++++++++++

*simSignals:scattered*

This keyword is currently not supported.

Scattered light count rate image file. This file is assumed to contain a 2-dimensional array of signals in units of ADU per second. The file is read in, scaled by the :ref:`scatteredscale <scatteredscale>` value, and added to the seed image. Leave as None to skip this step.

.. _scatteredscale:

Scattered light scaling factor
++++++++++++++++++++++++++++++

*simSignals:scatteredscale*

Scaling factor to multiply the :ref:`scattered light count rate image <scattered>` by prior to adding to the seed image.

.. _bkgdrate:

Background signal
+++++++++++++++++

*simSignals:bkgdrate*

This entry, in combination with the :ref:`use_dateobs_for_background <use_dateobs_for_background>` and :ref:`date_obs <date_obs>` parameters, controls the background signal that is added to simulations. The text below describes the way Mirage interprets the various input options:


**Imaging Mode (both NIRCam and NIRISS)**

- Number: The input value is assumed to be in units of ADU/pixel/second. This constant background value is placed in all pixels.
- "low", “medium”, or “high”. If one of these options is used, the simulator uses the `jwst_backgrounds <https://github.com/spacetelescope/jwst_backgrounds>`_ repository to calculate the background rate to apply to the simulated data. The package calculates the background signal at the requested pointing on the sky for each night over the course of a year and creates a histogram of these values. If the requested background is "low" then the returned background level is equal to that of the 10th percentile in the histogram. A "medium" background corresponds to the 50th percentile value, and "high" is the 90th percentile value. In this case, the returned background rate includes contributions from zodiacal light and telescope thermal emission.
- :ref:`use_dateobs_for_background <use_dateobs_for_background>` set to True: (NOTE: currently the bkgdrate value must be set to "low", "medium", or "high" when using this option. If it is set to a number, then that number will be used and use_dateobs_for_background will be ignored.) This is similar to the "low", “medium”, “high” case above, but instead of calculating the background based on a percetile of the distribution of background values, Mirage will select the background value associated with the date in the :ref:`date_obs <date_obs>` parameter.


**WFSS Mode**

NIRCam

- Number: Not supported
- "low", “medium”, or “high”. Similar to the imaging case above. In this case, the background spectrum matching the percentile value is kept. This is fed into the disperser software, which generates a 2D background image.
- :ref:`use_dateobs_for_background <use_dateobs_for_background>` set to True. The background spectrum for the date in the :ref:`date_obs <date_obs>` parameter is fed into the disperser, which generates a 2D background image.

NIRISS

- Number: The input number is assumed to be the desired background value in ADU/pixels/second in the **undispersed view** of the scene. To get the background value in the dispersed image, this number is multiplied by the throughput of the NIRISS grism, which is about 80%. The dispersed background image, which is in the collection of Mirage reference files, is then scaled such that the mean value is equal to the calculated dispersed background value.
- "low", “medium”, or “high”. Same as in the imaging case above. The calculated backrgound value will be multiplied by the throughput of the NIRISS grism, which is about 80%.
- :ref:`use_dateobs_for_background <use_dateobs_for_background>`. Not supported

Note that background rates associated with the "low", "medium", and "high" values are calculated in the same way as when they are used in the `JWST ETC <https://jwst.etc.stsci.edu/>`_.

.. _poissonseed:

Seed value for poisson noise generator
++++++++++++++++++++++++++++++++++++++

*simSignals:poissonseed*

Random number generator seed used for Poisson simulation

.. _photonyield:

Photon Yield
++++++++++++

*simSignals:photonyield*

This keyword is currently not used. T/F. Set this to **True** to include the effects of photon yield in the simulation outputs.

.. _pymethod:

Photon yield method
+++++++++++++++++++

*simSignals:pymethod*

This keyword is currently not used. T/F. Whether or not to use the double photon method when applying photon yield.

.. _expand_catalog_for_segments:

Expand catalog for segments
+++++++++++++++++++++++++++

*simSignals:expand_catalog_for_segments*

This entry controls whether Mirage will look for a separate point source library for each of the mirror segments on the telescope. This
mode is only used for certain wavefront sensing and control observations and should normally be set to False.

.. _use_dateobs_for_background:

Use date_obs for background
+++++++++++++++++++++++++++

*simSignals:use_dateobs_for_background*

This entry controls the way the background signal for the observation is calculated. If it is True, then the background value will be created by extracting the background spectrum assoicated with :ref:`date_obs <date_obs>` from the `jwst_backgrounds <https://github.com/spacetelescope/jwst_backgrounds>`_ package. If False, the background will be determined by calculating the background value at a certain percentile of the collection of backgrounds for the given pointing over 365 days. If :ref:`bkgdrate <bkgdrate>` is "low", "medium", "high", then the percentiles used are 10th, 50th, and 90th, respectively. If it is a float, that value (in ADU/sec/pixel) will be added to all pixels.

.. _Telescope:

Telescope section
-----------------

Inputs in this section of the yaml file describe the telescope pointing to use for the simulation.

.. _ra:

Right Ascension
+++++++++++++++

*Telescope:ra*

Right ascension of the observation. This will be the RA at the reference location on the detector being used for the simulation. The reference location varies with the requested subarray, but is generally in the center of the field of view. This input can be a string "HH:MM:SS.sss", or a float in decimal degrees.

.. _dec:

Declination
+++++++++++

*Telescope:dec*

Declination of the observation. This will be the Dec at the reference location on the detector. The reference location varies with the requested subarray, but is generally in the center of the field of view. This input can be a string "DD:MM:SS.sss" or a float in decimal degrees.

.. _rotation:

Rotation
++++++++

*Telescope:rotation*

Rotation of the y-axis in degrees East of North. Currently this rotation is defined around the reference location of the chosen subarray.

.. _tracking:

Telescope tracking
++++++++++++++++++

*Telescope:tracking*

Either 'sidereal' or 'non-sidereal' depending on the type of exposure. If it is set to non-sidereal then the exposure will be created as if JWST is
tracking on the source in the :ref:`movingTargetToTrack <movingTargetToTrack>` catalog. Sources in the :ref:`pointsource <pointsource>`, :ref:`galaxyListFile <galaxyListFile>`, and :ref:`extended <extended>` catalogs will trail across the field of view over the course of the exposure.

.. _newRamp:

newRamp section
---------------

This section of the input file lists JWST calibration pipeline-style configuration files that may be needed when preparing the simulated data. Copies of all of these configuration files are included in the ‘config’ subdirectory of the MIRAGE repository. Therefore, unless you wish to use your own set of configuration files, you can set these fields all to 'config'. This is the default behavior when creating yaml files via the :ref:`yaml generator <yaml_generator>`.

.. hint::
	In order to create your own set of pipeline configuration files, use the shell command:

	> collect_pipeline_cfg /your/destination/directory

.. _dq_configfile:

DQ step configuration file
++++++++++++++++++++++++++

*newRamp:dq_configfile*

Name of the JWST calibration pipeline configuration file to be used in the dq_init step when it is run on the raw dark current integration.


.. _sat_configfile:

Saturation step configuration file
++++++++++++++++++++++++++++++++++

*newRamp:sat_configfile*

Name of the JWST calibration pipeline configuration file to be used in the saturation step when it is run on the raw dark current integration.

.. _superbias_configfile:

Superbias step configuration file
+++++++++++++++++++++++++++++++++

*newRamp:superbias_configfile*

Name of the JWST calibration pipeline configuration file to be used in the superbias step when it is run on the raw dark current integration.

.. _refpix_configfile:

Reference pixel subtraction configuration file
++++++++++++++++++++++++++++++++++++++++++++++

*newRamp:refpix_configfile*

Name of the JWST calibration pipeline configuration file to be used in the reference pixel subtraction step when it is run on the raw dark current integration.

.. hint::
    If you choose to use your own reference pixel correction configuration file, we recommend setting the **odd_even_rows** entry to False, as this correction is not typically performed on NIRCam, NISISS, or FGS data.

.. _linear_configfile:

Linearity step configuration file
+++++++++++++++++++++++++++++++++

*newRamp:linear_configfile*

Name of the JWST calibration pipeline configuration file to be used in the linearity correction step when it is run on the raw dark current integration.

.. _output:

Output section
--------------

This section of the yaml file contains information about the output file, such as filename and location. In addition, this section contains a large number of fields that describe how this particular exposure fits within an observing program/proposal. This information is not used during the creation of the simulated data, but is placed in the header of the output file in order to be consistent with the contents of real JWST data files. In addition, `level 3 of the JWST calibration pipeline <https://jwst-pipeline.readthedocs.io/en/stable/jwst/pipeline/description.html#pipelines>`_, which is used to combine multiple exposures into mosaic images, does require some of this information. The easiest way to correctly populate this information in the simulator yaml files is to :ref:`create the yaml files from an APT file via yaml_generator.py<from_apt>`, in which case the fields are all populated automatically.

.. _file:

Output filename
+++++++++++++++

*Output:file*

Filename of the output simulated file (e.g. jw42424024002_01101_00001_nrcb5_uncal.fits). If the linearized ramp is requested as output in the :ref:`datatype<datatype>` field, it will be saved with ‘uncal’ replaced with ‘linear’ in the filename or if ‘uncal’ is not present, ‘linear’ will simply be appended to the filename.  If the raw ramp is requested as output, the given filename will be used with no changes.

We recommend using filenames that end in 'uncal.fits' in order to be consistent with `JWST file naming conventions <https://jwst-docs.stsci.edu/display/JDAT/File+Naming+Conventions+and+Data+Products>`_. The filename is constructed from various pieces of information, including the program ID and visit number. If you wish to use this convention for the output filenames, the easiest way to accomplish this is to :ref:`create the yaml files from an APT file <from_apt>`, in which case the filenames will be generated automatically.

.. _directory:

Output directory
++++++++++++++++

*Output:directory*

The directory into which the output simulated data will be placed.

.. _datatype:

Data type
+++++++++

*Output:datatype*

List of the data format(s) of the output files. Options include:
“linear”, where the output files will contain linearized signals with the superbias and reference pixel signals removed. Bad pixels will also be flagged if a bad pixel file is specified. These files are ready to be run through the jump detection and ramp fitting steps of the JWST calibration pipeline. “raw”, where the output files will be in an uncalibrated state. These files are ready to be run through the entirety of the calibration pipeline, beginning with `calwebb_detector1 <https://jwst-pipeline.readthedocs.io/en/stable/jwst/pipeline/description.html#pipelines>`_.
“linear,raw”, where both the raw and linearized versions of the output files will be saved.

.. _format:

Data format
+++++++++++

*Output:format*

Format of the output file. Currently, only ‘DMS’ is supported, indicating that the fits file format, as well as header keywords, match those expected by the JWST calibration pipeline.

.. _save_intermediates:

Save intermediate outputs
+++++++++++++++++++++++++

*Output:save_intermediates*

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



Grism output image
++++++++++++++++++

.. _grism_source_image:

*Output:grism_source_image*

True/False. If True, the size of the output image is enlarged from the requested array size by a multiplicative factor in the x and y dimensions. For NIRCam this factor is √2, while it NIRISS it is 1.134. This extra area is required if the image is passed to the grism disperser software. In this case, the disperser software is able to include sources which fall just outside the nominal field of view but whose dispersed spectra fall into the nominal field of view.

.. _unsigned:

Outputs in unsigned integers
++++++++++++++++++++++++++++

*Output:unsigned*

T/F. If True, output signal values for raw data will be in units of unsigned integers. This matches the output of real JWST data.

.. _dmsOrient:

Output data in DMS orientation
++++++++++++++++++++++++++++++

T/F. If True, data will be output in DMS orientation, as opposed to raw FITSwriter orientation. JWST data will be in DMS orientation.

.. _program_number:

Program number
++++++++++++++

*Output:program_number*

The proposal ID number. This is placed in the header of the output file in order to match the contents of real observation files.

.. _title:

Proposal title
++++++++++++++

*Output:title*

The title of the proposal. This placed in the header of the output file in order to match the contents of real observation files.

.. _PI_Name:

PI name
+++++++

*Output:PI_Name*

Name of the proposal PI. This is placed in the header of the output file in order to match the contents of real observation files.

.. _Proposal_category:

Proposal category
+++++++++++++++++

*Output:proposal_category*

Proposal category (e.g. GO, GTO). This is placed in the header of the output file in order to match the contents of real observation files.

.. _Science_category:

Science category
++++++++++++++++

*Output:science_category*

Science category of the proposal, as defined in the APT file. This is placed in the header of the output file in order to match the contents of real observation files.

.. _target_name:

Target Name
+++++++++++

*Output:target_name*

Name of the target. For yaml files constructed from an APT file, this is the name of the target as input by the user. This value will be propagated into the TARGPROP keyword in the simulated data FITS files.

.. _target_ra:

Target RA
+++++++++

*Output:target_ra*

RA of the target. For yaml files constructed from an APT file, this is the RA of the target as input by the user, translated to units of degrees. This value will be propagated into the TARG_RA keyword in the simulated data FITS files.

.. _target_dec:

Target Dec
++++++++++

*Output:target_dec*

Declination of the target. For yaml files constructed from an APT file, this is the declination of the target as input by the user, translated to units of degrees. This value will be propagated into the TARG_DEC keyword in the simulated data FITS files.

.. _observation_number:

Observation number
++++++++++++++++++

*Output:observation_number*

The observation number containing the output exposure, as defined in the program’s APT file. This is placed in the header of the output file in order to match the contents of real observation files.

.. _observation_label:

Observation label
+++++++++++++++++

*Output:observation_label*

The observation label in the APT file under which the output exposure appears. This is placed in the header of the output file in order to match the contents of real observation files.

.. _visit_number:

Visit number
++++++++++++

*Output:visit_number*

The visit number, as defined in the APT file, within which the output exposure appears. This is placed in the header of the output file in order to match the contents of real observation files.

.. _visit_group:

Visit group number
++++++++++++++++++

*Output:visit_group*

The visit group, as defined in the APT file, within which the output exposure appears. This is placed in the header of the output file in order to match the contents of real observation files.

.. _visit_id:

Visit ID number
+++++++++++++++

*Output:visit_id*

The visit identifier of the exposure. This can be created by combining the program ID, visit number, and observation number. This is placed in the header of the output file in order to match the contents of real observation files.

.. _sequence_id:

Sequence ID
+++++++++++

*Output:sequence_id*

The parallel sequence identifier denotes whether the data were acquired during parallel observations, and with which instrument. Set to 0 for non-parallel observations, 1 for a parallel sequence using the primary instrument, or 2-5 for one of the non-prime instruments.

.. _activity_id:

Activity ID
+++++++++++

*Output:activity_id*

The activity identifier of the exposure is a base-36 number that is unique to each exposure in a proposal. This is placed in the header of the output file in order to match the contents of real observation files.

.. _exposure_number:

Exposure Number
+++++++++++++++

*Output:exposure_number*

A five-character number used to identify the exposure within the current activity.

.. _obs_id:

Observation ID
++++++++++++++

*Output:obs_id*

The observation ID is constructed from several of the other parameters. OBS_ID = 'V' + program_number + observation_id + visit_id + 'P' + parallel-program number + parallel-observation number + visit_group + parallel sequence identifier + activity_identifier.

.. _date_obs:

Observation date
++++++++++++++++

*Output:date_obs*

UTC date of the start of the exposure with format yyyy-mm-dd.

.. _time_obs:

Observation time
++++++++++++++++

*Output:time_obs*

UTC time of the start of the exposure with format hh:mm:ss.ssssss.

.. _obs_template:

Observation template
++++++++++++++++++++

*Output:obs_template*

The name of the observation template used for the exposure (e.g. NIRCam Imaging, NIRCam Time Series)

.. _primary_dither_type:

Primary dither type
+++++++++++++++++++

*Output:primary_dither_type*

Name of the primary dither pattern in use when the data were obtained. For details, see the documentation pages on dither patterns for `NIRCam <https://jwst-docs.stsci.edu/display/JTI/NIRCam+Primary+Dithers>`_, and `NIRISS <https://jwst-docs.stsci.edu/display/JTI/NIRISS+Dithers>`_. (e.g. INTRAMODULEX, INTRASCA).

.. _total_primary_dither_positions:

Number of primary dither positions
++++++++++++++++++++++++++++++++++

*Output:total_primary_dither_positions*

Total number of primary dither positions in the observation.

.. _primary_dither_position:

Primary dither position
+++++++++++++++++++++++

*Output:primary_dither_position*

Primary dither position number of the exposure being simulated.

.. _subpix_dither_type:

Subpixel dither type
++++++++++++++++++++

*Output:subpix_dither_type*

Name of the subpixel dither pattern used for these data. Details on subpixel dither patterns can be found on the `NIRCam subpixel dither patterns page <https://jwst-docs.stsci.edu/display/JTI/NIRCam+Subpixel+Dithers>`_.

.. _total_subpix_dither_positions:

Number of subpixel dither positions
+++++++++++++++++++++++++++++++++++

*Output:total_subpix_dither_positions*

Total number of subpixel dither positions for this observation.

.. _subpix_dither_position:

Subpixel dither position
++++++++++++++++++++++++

*Output:subpix_dither_position*

The subpixel dither position number corresponding to the current exposure.

.. _xoffset:

X offset
++++++++

*Output:xoffset*

Offset in the x direction, in arcseconds, of the pointing used for the current exposure relative to the starting position of the dither pattern. This is used to populate header values only. It is not used to determine the pointing when creating the simulated data.

.. _yoffset:

Y offset
++++++++

*Output:yoffset*

Offset in the y direction, in arcseconds, of the pointing used for the current exposure relative to the starting position of the dither pattern. This is used to populate header values only. It is not used to determine the pointing when creating the simulated data.

