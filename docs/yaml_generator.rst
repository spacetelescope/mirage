.. _from_apt:

Simulating Observations from an APT File
========================================

The easiest way to simulate a large amount of obsrvtions using Mirage is to begin with the `APT <https://jwst-docs.stsci.edu/display/JPP/JWST+Astronomers+Proposal+Tool%2C+APT>`_ file of a JWST proposal. Mirage contains functions that can (for imaging and WFSS modes) create the collection of yaml files necessary to simulate all of the observations contained in the APT file.

Export XML and Poining files from APT
-------------------------------------
More specifically, open the proposal within APT and under the File menu, choose Export, and select the xml and pointing files. Only one file at a time can be exported, so the process must be repeated for each.

Create an observation list file
-------------------------------
One other input file is needed in addition to those from the APT proposal. This is an observation list file. This file contains a list of the source catalogs associated with each of the observations in the proposal in both the longwave and shortwave channels. This file is necessary for the catalog names to be populated witin all of the yaml files. Through this file, you can also specify the date as well as the roll angle (V3 position angle) for each observation. In this way, you can simulate different epochs in your observations, with associated changes in telescope roll angle. An example of an observation list file is shown below.

::

	# Example observation list for the APT -> simulator code.
	#
	# Observation names must match up with those in the xml
	# file from APT.
	# Observation dates should be in the format: YYYY-MM-DD
	# Roll angle is PAV3, which is available from the
	# JWST General Target Visibility Tool (GTVT)
	# Roll angles and dates are not checked against the
	# GTVT and do not need to agree with the tool's outputs
	#
	# For multiple observations that you want taken back-to-back,
	# use the same date for all

	Observation1:
  	  Name: 'Target #1, Epoch 1'
  	  Date: 2019-10-14
  	  PAV3: 0.
  	  FilterConfig1:
  	    SW:
  	      Filter: F200W
  	      PointSourceCatalog: ptsrc_f200w.cat
  	      GalaxyCatalog: galaxies_f200w.cat
  	      ExtendedCatalog: None
  	      ExtendedScale: 1.0
  	      ExtendedCenter: 1024,1024
  	      MovingTargetList: None
  	      MovingTargetSersic: None
  	      MovingTargetExtended: None
  	      MovingTargetConvolveExtended: True
  	      MovingTargetToTrack: None
  	      BackgroundRate: medium
  	    LW:
  	      Filter: F444W
  	      PointSourceCatalog: ptsrc_f444w.cat
  	      GalaxyCatalog: galaxies_f444w.cat
  	      ExtendedCatalog: None
  	      ExtendedScale: 1.0
  	      ExtendedCenter: 1024,1024
  	      MovingTargetList: None
  	      MovingTargetSersic: None
  	      MovingTargetExtended: None
  	      MovingTargetConvolveExtended: True
  	      MovingTargetToTrack: None
  	      BackgroundRate: medium
  	  FilterConfig2:
  	    SW:
  	      Filter: F210M
  	      PointSourceCatalog: ptsrc_f210m.cat
  	      GalaxyCatalog: galaxies_f210m.cat
  	      ExtendedCatalog: None
  	      ExtendedScale: 1.0
  	      ExtendedCenter: 1024,1024
  	      MovingTargetList: None
  	      MovingTargetSersic: None
  	      MovingTargetExtended: None
  	      MovingTargetConvolveExtended: True
  	      MovingTargetToTrack: None
  	      BackgroundRate: medium
  	    LW:
  	      Filter: F460M
  	      PointSourceCatalog: ptsrc_f460m.cat
  	      GalaxyCatalog: galaxies_f460m.cat
  	      ExtendedCatalog: None
  	      ExtendedScale: 1.0
  	      ExtendedCenter: 1024,1024
  	      MovingTargetList: None
  	      MovingTargetSersic: None
  	      MovingTargetExtended: None
  	      MovingTargetConvolveExtended: True
  	      MovingTargetToTrack: None
  	      BackgroundRate: medium

	Observation2:
	  Name: 'Target #1, Epoch 2'
	  Date: 2020-04-14
	  PAV3: 178.
	  FilterConfig1:
	    SW:
	      Filter: F200W
	      PointSourceCatalog: ptsrc_f200w.cat
	      GalaxyCatalog: galaxies_f200w.cat
	      ExtendedCatalog: None
	      ExtendedScale: 1.0
	      ExtendedCenter: 1024,1024
	      MovingTargetList: None
	      MovingTargetSersic: None
	      MovingTargetExtended: None
	      MovingTargetConvolveExtended: True
	      MovingTargetToTrack: None
	      BackgroundRate: medium
	    LW:
	      Filter: F444W
	      PointSourceCatalog: ptsrc_f444w.cat
	      GalaxyCatalog: galaxies_f444w.cat
	      ExtendedCatalog: None
	      ExtendedScale: 1.0
	      ExtendedCenter: 1024,1024
	      MovingTargetList: None
	      MovingTargetSersic: None
	      MovingTargetExtended: None
	      MovingTargetConvolveExtended: True
	      MovingTargetToTrack: None
	      BackgroundRate: medium
	  FilterConfig2:
	    SW:
	      Filter: F210M
	      PointSourceCatalog: ptsrc_f210m.cat
	      GalaxyCatalog: galaxies_f210m.cat
	      ExtendedCatalog: None
	      ExtendedScale: 1.0
	      ExtendedCenter: 1024,1024
	      MovingTargetList: None
	      MovingTargetSersic: None
	      MovingTargetExtended: None
	      MovingTargetConvolveExtended: True
	      MovingTargetToTrack: None
	      BackgroundRate: medium
	    LW:
	      Filter: F460M
	      PointSourceCatalog: ptsrc_f460m.cat
	      GalaxyCatalog: galaxies_f460m.cat
	      ExtendedCatalog: None
	      ExtendedScale: 1.0
	      ExtendedCenter: 1024,1024
	      MovingTargetList: None
	      MovingTargetSersic: None
	      MovingTargetExtended: None
	      MovingTargetConvolveExtended: True
	      MovingTargetToTrack: None
	      BackgroundRate: medium


.. _yaml_generator:

Run the yaml generator
----------------------

With the XML, pointing, and observation list files in hand, Mirage's *yaml_generator.py* module can be called to create the associated yaml files.

::

	from mirage.yaml import yaml_generator

	# Create a series of data simluator input yaml files
	# from APT files

	yam = yaml_generator.SimInput()
	yam.input_xml = 'example_imaging_program.xml'
	yam.pointing_file = 'example_imaging_program.pointing'
	yam.output_dir = './'
	yam.simdata_output_dir = './'
	yam.observation_table = 'observation_list.yaml'
	yam.use_JWST_pipeline = True
	yam.use_linearized_darks = False
	yam.datatype = 'linear'
	yam.reffile_setup()
	yam.create_inputs()

The outptut from this will be the collection of yaml files needed to run Mirage and create all of the observation files. An :ref:`example yaml file <example_yaml>` shows all of the parameters necessary when simulating an exposure.

Run Mirage
----------

The collection of yaml files can then be fed into Mirage one at a time.

::

	from glob import glob
	from mirage import imaging_simulator

	yaml_files = glob('*.yaml')
	for yfile in yaml_files:
	    im = imaging_simulator.ImgSim()
	    im.paramfile = yfile
	    im.create()




