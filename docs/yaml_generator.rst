.. _from_apt:

Simulating Observations from an APT File
========================================

The easiest way to simulate a large amount of obsrvtions using Mirage is to begin with the `APT <https://jwst-docs.stsci.edu/display/JPP/JWST+Astronomers+Proposal+Tool%2C+APT>`_ file of a JWST proposal. Mirage contains functions that can (for imaging and WFSS modes) create the collection of yaml files necessary to simulate all of the observations contained in the APT file. Specifically, the *yaml_generator.py* module will take the APT files, along with several other inputs, and will output yaml files that Mirage can use to create all the observations contained in the APT files.

Examples of this functionality are shown in the `notebooks <https://github.com/spacetelescope/mirage/tree/master/examples>`_ in the *Mirage* repository. Below are step-by-step instructions explaining how to create *Mirage* input yaml files from an APT file.

Export XML and Pointing files from APT
--------------------------------------
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

.. _override_reffiles:

JWST Calibration Reference Files
++++++++++++++++++++++++++++++++

Mirage makes use of a handful of the `reference file types <https://jwst-pipeline.readthedocs.io/en/stable/jwst/introduction.html#reference-files>`_ used by the JWST calibration pipeline. This includes the `bad pixel mask <https://jwst-pipeline.readthedocs.io/en/stable/jwst/dq_init/reference_files.html#mask-reffile>`_, `saturation level map <https://jwst-pipeline.readthedocs.io/en/stable/jwst/saturation/reference_files.html#saturation-reffile>`_, `superbias <https://jwst-pipeline.readthedocs.io/en/stable/jwst/superbias/reference_files.html#superbias-reffile>`_, `gain <https://jwst-pipeline.readthedocs.io/en/stable/jwst/references_general/gain_reffile.html#gain-reffile>`_, `interpixel capacitance <https://jwst-pipeline.readthedocs.io/en/stable/jwst/ipc/reference_files.html#ipc-reffile>`_, `linearity correction  <https://jwst-pipeline.readthedocs.io/en/stable/jwst/linearity/reference_files.html#linearity-reffile>`_, and `distortion correction <https://jwst-pipeline.readthedocs.io/en/stable/jwst/references_general/distortion_reffile.html#distortion-reffile>`_ files.

Mirage relies on the `CRDS <https://hst-crds.stsci.edu/static/users_guide/index.html>`_ package from STScI to identify the appropriate reference files for a given exposure. These files are automatically downloaded to the user's machine at one of two times:

1) When running *yaml_generator.py* if the ``reffile_defaults`` option is set to 'crds_full_names' (see :ref:`Run the Yaml Generator <yaml_generator>` for more details).
2) When running any of the three main parts of Mirage: the :ref:`catalog seed generator <source_catalogs>`, :ref:`dark preparation <dark_prep>`, or :ref:`observation generator <obs_generator>`.

.. tip::

    In order to specify a location on your machine to store the downloaded reference files, you must have the CRDS_PATH environment variable set to that location. If this environment variable is not set, Mirage will default to use $HOME/crds_cache/

It is also possible to specify that Mirage use reference files other than those downloaded from CRDS. In order to do this, you must supply a dictionary of filenames. Due to the varying number of selection criteria needed to uniquely identify the reference file that matches up with a particular exposure, this dictionary is composed of multiple levels of nested dictionaries. Not all possibilities are required. You may specify reference files only for the particular observing modes you are interested in. Any modes in your APT file that are not contained within the dictionary will revert to using the reference files identified by CRDS. Below is a dictionary showing all nesting required for all reference files. Note that the dictionary structure is instument dependent since a reference file type for different instruments does not necessarily have the same selection criteria.

.. tip::
    If you choose to provide your own reference files, it is best to use these same reference files when running the JWST calibration pipeline on the simualted data files produced by Mirage. Not doing this can lead to systematic errors in your calibrated data.

Here is a view of the dictionary structure required when specifying reference files. All keys are assumed to be lower case strings. The easiest way to see the set of allowed values for a particular property (e.g. exposure_type), go to the `JWST Keyword Dictionary <https://mast.stsci.edu/portal/Mashup/Clients/jwkeywords/index.html>`_, and search for the appropriate keyword. Exposure type is EXP_TYPE, filter is FILTER, pupil is PUPIL, detctor_name is DETECTOR, readpattern is READPATT.

::

    override = {'nircam': {'superbias':  {detector_name: {readpattern: 'reffile_name.fits'}},
                           'linearity':  {detector_name: 'reffile_name.fits'},
                           'saturation': {detector_name: 'reffile_name.fits'},
                           'gain':       {detector_name: 'reffile_name.fits'},
                           'distortion': {detector_name: {filter: {exposure_type: 'reffile_name.asdf'}}},
                           'ipc':        {detector_name: 'reffile_name.fits'},
                           'area':       {detector_name: {filter: {pupil: {exposure_type: 'reffile_name.asdf'}}}},
                           'badpixmask': {detector_name: 'reffile_name.fits'}
                           },
                'niriss': {'superbias':  {readpattern: 'reffile_name.fits'},
                           'linearity':  'reffile_name.fits',
                           'saturation': 'reffile_name.fits',
                           'gain':       'reffile_name.fits',
                           'distortion': {pupil: {exposure_type: 'reffile_name.fits'}},
                           'ipc':        'reffile_name.fits',
                           'area':       {filter: {pupil: {exposure_type: 'reffile_name.asdf'}}},
                           'badpixmask': 'reffile_name.fits'
                           },
                'fgs':    {'superbias':  {detector_name: {readpattern: 'reffile_name.fits'}},
                           'linearity':  {detector_name: 'reffile_name.fits'},
                           'saturation': {detector_name: 'reffile_name.fits'},
                           'gain':       {detector_name: 'reffile_name.fits'},
                           'distortion': {detector_name: {exposure_type: 'reffile_name.fits'}},
                           'ipc':        {detector_name: 'reffile_name.fits'},
                           'area':       {detector_name: 'reffile_name.asdf'},
                           'badpixmask': {detector_name: {exposure_type: 'reffile_name.fits'}}
                           }


Here we show an example dictionary for a particular set of observations.

::

    override = {'nircam': {'superbias':  {'nrcb5': {'bright2': 'my_reffiles/my_superbias_for_b5.fits',
                                                    'rapid': 'my_reffiles/my_superbias_for_b5.fits',
                                                    'shallow4': 'my_reffiles/my_superbias_for_b5.fits'
                                                    },
                                          'nrcb4': {'rapid': 'my_reffiles/my_superbias_for_b4.fits'}
                                          },
                           'linearity':  {'nrcb5': 'my_reffiles/my_linearity_for_b5.fits',
                                          'nrcb4': 'my_reffiles/my_linearity_for_b4.fits'},
                           'saturation': {'nrcb5': 'my_reffiles/my_saturation_for_b5.fits',
                                          'nrcb4': 'my_reffiles/my_saturation_for_b4.fits'},
                           'gain':       {'nrcb5': 'my_reffiles/my_gain_for_b5.fits',
                                          'nrcb4': 'my_reffiles/my_gain_for_b4.fits'},
                           'distortion': {'nrcb5': {'f322w2': {'nrc_image': 'my_reffiles/my_distortion_for_b5.asdf'}},
                                          'nrcb4': {'f444w':  {'nrc_image': 'my_reffiles/my_distortion_for_b4.asdf'}}},
                           'ipc':        {'nrcb5': 'my_reffiles/my_ipc_for_b5.fits',
                                          'nrcb4': 'my_reffiles/my_ipc_for_b4.fits'},
                           'area':       {'nrcb5': {'f322w2': {'clear': {'nrc_image': 'my_reffiles/my_pam_for_b5.asdf'}}},
                                          'nrcb4': {'f444w':  {'clear': {'nrc_image': 'my_reffiles/my_pam_for_b4.asdf'}}}},
                           'badpixmask': {'nrcb5': 'my_reffiles/my_bpm_for_b5.fits',
                                          'nrcb4': 'my_reffiles/my_bpm_for_b4.fits'},
                            }
                }


.. _yaml_generator:

Run the Yaml Generator
----------------------

With the XML and pointing files in hand, and additional inputs defined above, Mirage's *yaml_generator.py* module can be called to create the associated yaml files. We specify a location for the oputput yaml files using the ``output_dir`` keyword. We also define the directory into which the final simulated data will be placed, using the ``simulated_output_dir`` keyword. This information will be placed into the constructed yaml files.

Setting the ``use_linearized_darks`` option to True will cause the *yaml_generator* to look for linearized dark current files to use with the simulations. These files may be present in the collection of *Mirage* :ref:`reference files <reference_files>`. If linearized darks are not present, leaving this option as False will cause `Mirage` to use raw dark current ramps as inputs.

The ``reffile_defaults`` keyword can have one of two values, which induce slightly different behavior. The best value to use depends upon your use case.

``crds`` (default) - This option will place the string 'crds' in the yaml file entries for CRDS reference files. When Mirage (i.e. the seed generator, dark prep, or observation generator) is run, it will query CRDS for the best reference files to use, and download those files to your CRDS_PATH directory if they are not already present. This has the advantage that your yaml files will always have Mirage use the latest, best reference files whenever they are used.

``crds_full_name`` - This option will cause *yaml_generator.py* to query CRDS, which will identify and download the best reference files if they are not already in your CRDS_PATH. It will then place the names of these best reference files into the yaml files being created. This creates yaml files that will always use the same reference files whenever they are run, meaning the outputs should be consistent every time. In this case, if new reference files are delivered to CRDS, the yaml files, and therefore Mirage, will not know that information.

Set ``reffile_overrides`` equal to the name of your nested reference file dictionary, if present.

Set ``parameter_defaults`` equal to the dictionary of parameter values to use.


::

    from mirage.yaml import yaml_generator

    yam = yaml_generator.SimInput(xml_file, pointing_file, catalogs=catalogs, verbose=True,
                                  output_dir='/location/to/place/yaml_files',
                                  simdata_output_dir='/location/to/place/simulated_data',
                                  parameter_defaults=param_values, datatype='raw',
                                  reffile_defaults='crds', reffile_overrides=reffile_overrides)
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




