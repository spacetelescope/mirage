.. _from_apt:

Simulating Observations from an APT File
========================================

The easiest way to simulate a large amount of obsrvtions using Mirage is to begin with the `APT <https://jwst-docs.stsci.edu/display/JPP/JWST+Astronomers+Proposal+Tool%2C+APT>`_ file of a JWST proposal. Mirage contains functions that can (for imaging and WFSS modes) create the collection of yaml files necessary to simulate all of the observations contained in the APT file.

Examples of this functionality are shown in the `notebooks <https://github.com/spacetelescope/mirage/tree/master/examples>`_ in the *Mirage* repository. Below are step-by-step instructions explaining how to create *Mirage* input yaml files from an APT file.

Export XML and Pointing files from APT
-------------------------------------
Open the proposal within APT and under the File menu, choose Export, and select the xml and pointing files. Only one file at a time can be exported, so the process must be repeated for each.

.. _additional_yaml_generator_inputs:

Additional Yaml Generator Inputs
--------------------------------

In addition to the xml and pointing files, there are several user-inputs needed by the yaml generator.

Source catalogs
+++++++++++++++

The user must specify source catalogs for each instrument and observation. This is done using a nested dictionary following the form shown below. Catalogs for NIRCam are broken down by channel (SW, LW). In each case, the entry can be a list of catalogs, in which case there must be one catalog in the list for each observation in the program that uses the instrument. For example, in an APT file containing 5 observations, if NIRCam is used in observations 1, 2, and 4, and NIRISS in observations 3 and 5, then the NIRCam catalog lists should have a list of three catalogs, and the NIRISS entry should have a list of 2 catalogs. The code assumes that the catalogs are listed in the same order of the observations in the APT file.

Alternatively, the user may provide a string containing a single catalog name, rather than a list. In this case, *Mirage* assumes that the same catalog will be used in all observations.

::

    catalogs = {'nircam': {'sw': {ptsrc: [list of catalogs], 'galaxy': [list of catalogs]}
                           'lw' : {ptsrc: [list of catalogs], 'galaxy': [list of catalogs]}
                'niriss': {'ptsrc': [list of catalogs], 'galaxy': [list of catalogs]}
               }


Note that currently this format is used for point source catalogs only. Other types of catalogs are limited to a single entry in a separate parameter dictionary as shown in the :ref:`Other Paramteres section <other_params>`. Catalog entries will be made more consistent in a future update. In the meantime, the easiest way to get the proper (non-point-source) catalogs into the yaml files is to use a single catalog for all observations and instruments. This is easily done by adding all necessary magnitude columns into a single catalog, or by combining existing *Mirage* catalogs. See the `Catalog Generation notebook <https://github.com/spacetelescope/mirage/blob/master/examples/Catalog_Generation_Tools.ipynb>`_ for examples.


.. Background Specification
.. ++++++++++++++++++++++++

.. The user must also supply a nested dictionary containing the background levels to use for each instrument and observation. Allowed values for the background parameter are described in :ref:`bkgdrate <bkgdrate>` section of the Input Yaml Parameters page.

.. ::

..     backgrounds = {'nircam': {'sw': [low, medium], 'lw': [low, medium]},
                   'niriss': [low, medium]
                  }

.. Similar to the case above for the catalog inputs, the given background values can be a list, in which case there must be one entry for each observation, or a string, in which case the value is applied to all observations.

.. _other_params:

Other Parameters
++++++++++++++++

There are currently a number of parameters which can only be set using the ``parameter_defaults`` keyword in the yaml_generator.  These parameters are currently limited in that only a single value is accepted for each, and applied to all observations in a proposal. An update to allow different values to different observations will be made soon. For the moment, the list of parameters that can be set with this method are shown below. The values are the defaults used by *Mirage* if the user does not enter their own values.

::

    param_values['Date'] = '2019-07-04'
    param_values['PAV3'] = '111.'
    param_values['GalaxyCatalog'] = 'None'
    param_values['ExtendedCatalog'] = 'None'
    param_values['ExtendedScale'] = '1.0'
    param_values['ExtendedCenter'] = '1024,1024'
    param_values['MovingTargetList'] = 'None'
    param_values['MovingTargetSersic'] = 'None'
    param_values['MovingTargetExtended'] = 'None'
    param_values['MovingTargetConvolveExtended'] = 'True'
    param_values['MovingTargetToTrack'] = 'None'
    param_values['BackgroundRate_sw'] = 'low'
    param_values['BackgroundRate_lw'] = 'low'
    param_values['BackgroundRate'] = '0.5'



.. _yaml_generator:

Run the Yaml Generator
----------------------

With the XML and pointing files in hand, and additional inputs defined above, Mirage's *yaml_generator.py* module can be called to create the associated yaml files. We specify a location for the oputput yaml files using the ``output_dir`` keyword. We also define the directory into which the final simulated data will be placed, using the ``simulated_output_dir`` keyword. This information will be placed into the constructed yaml files.

Setting the ``use_linearized_darks`` option to True will cause the *yaml_generator* to look for linearized dark current files to use with the simulations. These files may be present in the collection of *Mirage* :ref:`reference files <reference_files>`. If linearized darks are not present, leaving this option as False will cause `Mirage` to use raw dark current ramps as inputs.

Note that the point source catalogs and other parameters described above are inputs as well.

::

	  from mirage.yaml import yaml_generator

	  yam = yaml_generator.SimInput(xml_file, pointing_file, catalogs=catalogs, verbose=True,
                                  output_dir='/location/to/place/yaml_files',
                                  simdata_output_dir='/location/to/place/simulated_data',
                                  parameter_defaults=param_values, datatype='raw')
    yam.use_linearized_darks = True
    yam.create_inputs()


The outptut from this will be the collection of yaml files needed to run Mirage and create all of the simulated observation files. An :ref:`example yaml file <example_yaml>` shows all of the parameters necessary when simulating an exposure.

See the Imaging and WFSS notebooks in the `Mirage` repository for examples of *yaml_generator* use.

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




