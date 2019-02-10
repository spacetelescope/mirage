.. _example_yaml:

Example yaml Input File
=======================

Below is an example yaml input file for `Mirage`. For more information on the individual input paramters, see the :ref:`Input Yaml Parameters <input_yaml_file_parameters>` page.


::

	Inst:
	  instrument: NIRCam          #Instrument name
	  mode: imaging               #Observation mode (e.g. imaging, WFSS, moving_target)
	  use_JWST_pipeline: False    #Use pipeline in data transformations

	Readout:
	  readpatt: DEEP8         #Readout pattern (RAPID, BRIGHT2, etc) overrides nframe,nskip unless it is not recognized
	  ngroup: 6               #Number of groups in integration
	  nint: 1                 #Number of integrations per exposure
	  array_name: NRCB5_FULL  #Name of array (FULL, SUB160, SUB64P, etc)
	  filter: F250M           #Filter of simulated data (F090W, F322W2, etc)
	  pupil: CLEAR            #Pupil element for simulated data (CLEAR, GRISMC, etc)

	Reffiles:                 #Set to None or leave blank if you wish to skip that step
	  dark: None              #Dark current integration used as the base
	  linearized_darkfile: $MIRAGE_DATA/nircam/darks/linearized/B5/Linearized_Dark_and_SBRefpix_NRCNRCBLONG-DARK-60090141241_1_490_SE_2016-01-09T02h46m50_uncal.fits # Linearized dark ramp to use as input. Supercedes dark above
	  badpixmask: $MIRAGE_DATA/nircam/reference_files/badpix/NRCB5_17161_BPM_ISIMCV3_2016-01-21_ssbspmask_DMSorient.fits # If linearized dark is used, populate output DQ extensions using this file
	  superbias: $MIRAGE_DATA/nircam/reference_files/superbias/NRCB5_superbias_from_list_of_biasfiles.list.fits  #Superbias file. Set to None or leave blank if not using
	  linearity: $MIRAGE_DATA/nircam/reference_files/linearity/NRCBLONG_17161_LinearityCoeff_ADU0_2016-05-22_ssblinearity_v2_DMSorient.fits    #linearity correction coefficients
	  saturation: $MIRAGE_DATA/nircam/reference_files/saturation/NRCB5_17161_WellDepthADU_2016-03-10_ssbsaturation_wfact_DMSorient.fits    #well depth reference files
	  gain: $MIRAGE_DATA/nircam/reference_files/gain/NRCB5_17161_Gain_ISIMCV3_2016-02-25_ssbgain_DMSorient.fits #Gain map
	  pixelflat: None
	  illumflat: None                               #Illumination flat field file
	  astrometric: $MIRAGE_DATA/nircam/reference_files/distortion/NRCB5_FULL_distortion.asdf  #Astrometric distortion file (asdf)
	  ipc: $MIRAGE_DATA/nircam/reference_files/ipc/NRCB5_17161_IPCDeconvolutionKernel_2016-03-18_ssbipc_DMSorient.fits #File containing IPC kernel to apply
	  invertIPC: True       #Invert the IPC kernel before the convolution. True or False. Use True if the kernel is designed for the removal of IPC effects, like the JWST reference files are.
	  occult: None                                    #Occulting spots correction image
	  pixelAreaMap: $MIRAGE_DATA/nircam/reference_files/pam/NIRCam_B5_PAM_imaging.fits #Pixel area map for the detector. Used to introduce distortion into the output ramp.
	  subarray_defs:   config   #File that contains a list of all possible subarray names and coordinates
	  readpattdefs:    config   #File that contains a list of all possible readout pattern names and associated NFRAME/NSKIP values
	  crosstalk:       config   #File containing crosstalk coefficients
	  filtpupilcombo:  config   #File that lists the filter wheel element / pupil wheel element combinations. Used only in writing output file
	  flux_cal:        config   #File that lists flux conversion factor and pivot wavelength for each filter. Only used when making direct image outputs to be fed into the grism disperser code.

	nonlin:
	  limit: 60000.0        #Upper singal limit to which nonlinearity is applied (ADU)
	  accuracy: 0.000001    #Non-linearity accuracy threshold
	  maxiter: 10           #Maximum number of iterations to use when applying non-linearity
	  robberto:  False      #Use Massimo Robberto type non-linearity coefficients

	cosmicRay:
	  path: $MIRAGE_DATA/nircam/cosmic_ray_library/     #Path to CR library
	  library: SUNMIN    								#Type of cosmic rayenvironment (SUNMAX, SUNMIN, FLARE)
	  scale: 1.5     									#Cosmic ray rate scaling factor
	  suffix: IPC_NIRCam_B5    							#Suffix of library file names
	  seed: 2956411739      							#Seed for random number generator

	simSignals:
	  pointsource: my_point_sources.cat               #File containing a list of point sources to add (x,y locations and magnitudes)
	  psfpath: $MIRAGE_DATA/nircam/webbpsf_library/   #Path to PSF library
	  psfbasename: nircam                             #Basename of the files in the psf library
	  psfpixfrac: 0.25                                #Fraction of a pixel between entries in PSF library (e.g. 0.25 = files for PSF centered at 0.25 pixel intervals within pixel)
	  psfwfe: predicted                               #PSF WFE value ("predicted" or "requirements")
	  psfwfegroup: 0                                  #WFE realization group (0 to 4)
	  galaxyListFile: my_galaxies_catalog.list
	  extended: None                                 #Extended emission count rate image file name
	  extendedscale: 1.0                             #Scaling factor for extended emission image
	  extendedCenter: 1024,1024                      #x,y pixel location at which to place the extended image if it is smaller than the output array size
	  PSFConvolveExtended: True                      #Convolve the extended image with the PSF before adding to the output image (True or False)
	  movingTargetList: None                         #Name of file containing a list of point source moving targets (e.g. KBOs, asteroids) to add.
	  movingTargetSersic: None                       #ascii file containing a list of 2D sersic profiles to have moving through the field
	  movingTargetExtended: None                     #ascii file containing a list of stamp images to add as moving targets (planets, moons, etc)
	  movingTargetConvolveExtended: True             #convolve the extended moving targets with PSF before adding.
	  movingTargetToTrack: None                      #File containing a single moving target which JWST will track during observation (e.g. a planet, moon, KBO, asteroid)	This file will only be used if mode is set to "moving_target"
	  zodiacal:  None                                #Zodiacal light count rate image file
	  zodiscale:  1.0                                #Zodi scaling factor
	  scattered:  None                               #Scattered light count rate image file
	  scatteredscale: 1.0                            #Scattered light scaling factor
	  bkgdrate: 0.0                                  #Constant background count rate (electrons/sec/pixel)
	  poissonseed: 2012872553                        #Random number generator seed for Poisson simulation)
	  photonyield: True                              #Apply photon yield in simulation
	  pymethod: True                                 #Use double Poisson simulation for photon yield

	Telescope:
	  ra: 53.1                     #RA of simulated pointing
	  dec: -27.8                   #Dec of simulated pointing
	  rotation: 0.0                #y axis rotation (degrees E of N)

	newRamp:
	  dq_configfile: config          #config file used by JWST pipeline
	  sat_configfile: config         #config file used by JWST pipeline
	  superbias_configfile: config   #config file used by JWST pipeline
	  refpix_configfile: config      #config file used by JWST pipeline
	  linear_configfile: config      #config file used by JWST pipeline

	Output:
	  file: jw42424024002_0112o_NRCB5_uncal.fits   # Output filename
	  directory: ./   							   # Directory in which to place output files
	  datatype: linear,raw 						   # Type of data to save. 'linear' for linearized ramp. 'raw' for raw ramp. 'linear,raw' for both
	  format: DMS          						   # Output file format Options: DMS, SSR(not yet implemented)
	  save_intermediates: False   				   # Save intermediate products separately (point source image, etc)
	  grism_source_image: False   				   # Create an image to be dispersed?
	  unsigned: True   							   # Output unsigned integers? (0-65535 if true. -32768 to 32768 if false)
	  dmsOrient: True    						   # Output in DMS orientation (vs. fitswriter orientation).
	  program_number: 42424    					   # Program Number
	  title: Supernovae and Black Holes Near Hyperspatial Bypasses   #Program title
	  PI_Name: Doug Adams  						   # Proposal PI Name
	  Proposal_category: GO  					   # Proposal category
	  Science_category: Cosmology  				   # Science category
	  observation_number: '002'    				   # Observation Number
	  observation_label: Obs2    				   # User-generated observation Label
	  visit_number: '024'    					   # Visit Number
	  visit_group: '01'    						   # Visit Group
	  visit_id: '42424024002'    				   # Visit ID
	  sequence_id: '2'    						   # Sequence ID
	  activity_id: '2o'    						   # Activity ID. Increment with each exposure.
	  exposure_number: '00001'    				   # Exposure Number
	  obs_id: 'V42424024002P000000000112o'   	   # Observation ID number
	  date_obs: '2019-10-15'  					   # Date of observation
	  time_obs: '06:29:11.852'  				   # Time of observation
	  obs_template: 'NIRCam Imaging'  			   # Observation template
	  primary_dither_type: NONE  				   # Primary dither pattern name
	  total_primary_dither_positions: 1  		   # Total number of primary dither positions
	  primary_dither_position: 1  				   # Primary dither position number
	  subpix_dither_type: 2-POINT-MEDIUM-WITH-NIRISS  #Subpixel dither pattern name
	  total_subpix_dither_positions: 2  		   # Total number of subpixel dither positions
	  subpix_dither_position: 2  				   # Subpixel dither position number
	  xoffset: 344.284  						   # Dither pointing offset in x (arcsec)
	  yoffset: 466.768  						   # Dither pointing offset in y (arcsec)