.. _from_apt:

Simulating Observations from an APT File
========================================

The easiest way to simulate a large amount of obsrvtions using Mirage is to begin with the `APT <https://jwst-docs.stsci.edu/display/JPP/JWST+Astronomers+Proposal+Tool%2C+APT>`_ file of a JWST proposal. Mirage contains functions that can (for imaging and WFSS modes) create the collection of yaml files necessary to simulate all of the observations contained in the APT file.

Export XML and Poining files from APT
-------------------------------------
More specifically, open the proposal within APT and under the File menu, choose Export, and select the xml and pointing files. Only one file at a time can be exported, so the process must be repeated for each.

.. _additional_yaml_generator_inputs:

Additional Yaml Generator Inputs
--------------------------------

In addition to the xml and pointing files, there are several user-inputs needed by the yaml generator.

Source catalogs
+++++++++++++++

The user must specify source catalogs for each instrument and observation. This is done using a nested dictionary following the form:

::

    catalogs = {'nircam': {'sw': {ptsrc: [list of catalogs], 'galaxy': [list of catalogs]}
                           'lw' : {ptsrc: [list of catalogs], 'galaxy': [list of catalogs]}
                'niriss': {'ptsrc': [list of catalogs], 'galaxy': [list of catalogs]}
               }

For each `list of catalogs` instance above, the user may enter a list of catalogs, where there is one catalog for each observation in the APT file where that instrument is used. Alternatively, the user may enter a string containing the name of a single catalog. In this case, the same catalog is used for all observations.

Background Specification
++++++++++++++++++++++++

The user must also supply a nested dictionary containing the background levels to use for each instrument and observation. Allowed values for the background parameter are described in :ref:`bkgdrate <bkgdrate>` section of the Input Yaml Parameters page.

::

    backgrounds = {'nircam': {'sw': [low, medium], 'lw': [low, medium]},
                   'niriss': [low, medium]
                  }

Similar to the case above for the catalog inputs, the given background values can be a list, in which case there must be one entry for each observation, or a string, in which case the value is applied to all observations.

Other Parameters
++++++++++++++++

::

     parameter_defaults['PAV3'] = [110, 110, 20, 20]



.. _yaml_generator:

Run the Yaml Generator
----------------------

With the XML and pointing files in hand, and additional inputs defined above, Mirage's *yaml_generator.py* module can be called to create the associated yaml files. In this case, we specify that the yaml files produced by the *yaml_generator* will be output within a subdirectory of the working directory called *yaml_files*. We also define that when Mirage is run, the simulated data will be output into a subdirectory called *simulated_data*. This information will be placed into the constructed yaml files.

::

	from mirage.yaml import yaml_generator

	# Create a series of data simluator input yaml files
	# from APT files

	yam = yaml_generator.SimInput(xml_file, pointing_file, catalogs=catalogs, verbose=True,
                                  output_dir='yaml_files', simdata_output_dir='simulated_data',
                                  parameter_defaults=params, datatype='raw')
    yam.use_linearized_darks = True
    yam.create_inputs()

The outptut from this will be the collection of yaml files needed to run Mirage and create all of the simulated observation files. An :ref:`example yaml file <example_yaml>` shows all of the parameters necessary when simulating an exposure.

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




