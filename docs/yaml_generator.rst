.. _from_apt:

Simulating Observations from an APT File
========================================

The easiest way to simulate a large amount of obsrvtions using Mirage is to begin with the `APT <https://jwst-docs.stsci.edu/display/JPP/JWST+Astronomers+Proposal+Tool%2C+APT>`_ file of a JWST proposal. Mirage contains functions that can (for imaging and WFSS modes) create the collection of yaml files necessary to simulate all of the observations contained in the APT file. Specifically, the *yaml_generator.py* module will take the APT files, along with several other inputs, and will output yaml files that Mirage can use to create all the observations contained in the APT files.

Examples of this functionality are shown in the `notebooks <https://github.com/spacetelescope/mirage/tree/master/examples>`_ in the *Mirage* repository. Below are step-by-step instructions explaining how to create *Mirage* input yaml files from an APT file.


.. important::

    Currently Mirage is able to parse the following APT templates. Mirage will fail if your APT file contains observations using a template not on this list.

    +--------------------------------+
    |         NIRCAM                 |
    +================================+
    |  - NircamImaging               |
    |  - NircamWfss                  |
    |  - NircamEngineeringImaging    |
    |  - NircamDark                  |
    |  - NircamTimeSeries            |
    |  - NircamGrismTimeSeries       |
    +--------------------------------+

    +-------------------------------+
    |         NIRISS                |
    +===============================+
    |  - NirissExternalCalibration  |
    |  - NirissWfss                 |
    |  - NirissAmi                  |
    |  - NirissDark                 |
    +-------------------------------+

    +-------------------------------+
    |          FGS                  |
    +===============================+
    |  - FgsExternalCalibration     |
    +-------------------------------+

    +-------------------------------+
    |     Other Instruments         |
    | (for parallel observations)   |
    +===============================+
    |  - NirspecImaging             |
    |  - NirspecInternalLamp        |
    |  - MiriMRS                    |
    +-------------------------------+

    +-------------------------------+
    | Wavefront Sensing and Control |
    +===============================+
    |  - WfscCommissioning          |
    |  - WfscGlobalAlignment        |
    |  - WfscCoarsePhasing          |
    |  - WfscFinePhasing            |
    +-------------------------------+


Export XML and Pointing files from APT
--------------------------------------
Open the proposal within APT and under the File menu, choose Export, and select the xml and pointing files. Only one file at a time can be exported, so the process must be repeated for each.

.. _additional_yaml_generator_inputs:

Additional Yaml Generator Inputs
--------------------------------

In addition to the xml and pointing files, there are several user-inputs needed by the yaml generator. Examples shown below are python statements.

.. _yam_gen_cat_inputs:

Source catalogs
+++++++++++++++

The user must specify source catalogs for each instrument and target in the proposal. This is done using one of several possible input formats, as shown below.

If you do not wish to specify catalogs, and want to leave all catalog inputs as 'none', set your catalogs variable to None, or omit completely:

::

    catalogs = None

To use the same catalogs for all instruments in a given target, create a nested dictionary where the top-level keys exactly match the target names in the APT file. In the second level of the dictionary,
specify catalog filenames for each type of catalog. Note that you only need to include the catalog types that you want to use for your data. All 7 catalog types are shown here for completeness, but
in most cases you will only specify 1-3 catalog types.

::

    catalogs = {'NGC4242': {'point_source': 'target1_ptsrc.cat',
                            'galaxy': 'target1_galaxies.cat',
                            'extended': 'target1_ext.cat',
                            'moving_pointsource': 'target1_mt_ptsrc.cat',
                            'moving_sersic': 'target1_mt_gal.cat',
                            'moving_extended': 'target1_mt_ext.cat',
                            'moving_target_to_track': 'target1_mt_track.cat'
                            },
                'LMC': {'point_source': 'lmc_ptsrc.cat',
                        'galaxy': 'lmc_galaxy.cat',
                        'extended': 'lmc_ex.cat',
                        'moving_pointsource': 'lmc_mt_ptsrc.cat',
                        'moving_sersic': 'lmc_mt_gal.cat',
                        'moving_extended': 'lmc_mt_ext.cat',
                        'moving_target_to_track': 'lmc_mt_track.cat'
                          }
                }

For more fine-grained control, you can also specify catalogs for each instrument and target. This adds another level to the input dictionary. Again, note that the example below shows all three instruments
for completeness. You only need to include the instruments used in the proposal. Instrument and catalog type dictionary keys are case-insensitive. Target name keys must match exactly those in the APT file.
The dictionary currently does not break out separate catalogs for NIRCam shortwave and longwave channels. Note that Mirage contains :ref:`catalog tools <catalog_generation>` that can be used to combine catalogs.
See the `Catalog Generation notebook <https://github.com/spacetelescope/mirage/blob/master/examples/Catalog_Generation_Tools.ipynb>`_ for examples.

::

    catalogs = {'NGC1234': {'nircam': {'point_source': 'ngc1234_ptsrc_nrc.cat',
                                       'galaxy': 'ngc1234_galaxy_nrc.cat',
                                       'extended': 'ngc1234_ex_nrc.cat',
                                       'moving_pointsource': 'ngc1234_mt_ptsrc_nrc.cat',
                                       'moving_sersic': 'ngc1234_mt_gal_nrc.cat',
                                       'moving_extended': 'ngc1234_mt_ext_nrc.cat',
                                       'moving_target_to_track': 'ngc1234_mt_track_nrc.cat'
                                       },
                            'niriss': {'point_source': 'ngc1234_ptsrc_nis.cat',
                                      'galaxy': 'ngc1234_galaxy_nis.cat',
                                      'extended': 'ngc1234_ex_nis.cat',
                                      'moving_pointsource': 'ngc1234_mt_ptsrc_nis.cat',
                                      'moving_sersic': 'ngc1234_mt_gal_nis.cat',
                                      'moving_extended': 'ngc1234_mt_ext_nis.cat',
                                      'moving_target_to_track': 'ngc1234_mt_track_nis.cat'
                                      }
                           },
                'LMC': {'nircam': {'point_source': 'lmc_ptsrc_nrc.cat',
                                     'galaxy': 'lmc_galaxy_nrc.cat',
                                     'extended': 'lmc_ex_nrc.cat',
                                     'moving_pointsource': 'lmc_mt_ptsrc_nrc.cat',
                                     'moving_sersic': 'lmc_mt_gal_nrc.cat',
                                     'moving_extended': 'lmc_mt_ext_nrc.cat',
                                     'moving_target_to_track': 'lmc_mt_track_nrc.cat'
                                     },
                        'niriss': {'point_source': 'lmc_ptsrc_nis.cat',
                                     'galaxy': 'lmc_galaxy_nis.cat',
                                     'extended': 'lmc_ex_nis.cat',
                                     'moving_pointsource': 'lmc_mt_ptsrc_nis.cat',
                                     'moving_sersic': 'lmc_mt_gal_nis.cat',
                                     'moving_extended': 'lmc_mt_ext_nis.cat',
                                     'moving_target_to_track': 'lmc_mt_track_nis.cat'
                                     }
                          },
                }

.. _yam_gen_background_inputs:

Background Specification
++++++++++++++++++++++++

Users may also supply information on the background levels to use for each instrument and observation. Allowed values for the background parameter are described in :ref:`bkgdrate <bkgdrate>` section of the Input Yaml Parameters page.

If *use_dateobs_for_background* is True, then the background for all exposures will be determined based on the value of the observation date (discussed below). The `JWST backgrounds package <https://github.com/spacetelescope/jwst_backgrounds>`_ will be used, and the background associated with the requested pointing and date will be used.

If *use_dateobs_for_background* is False, then the *background* parameter controls the calculation of the background. As with the catalogs above, there are several formats that can be used to supply the background information. Examples are shown below.

To use the default value for background ("low"), either omit this parameter altogether, or set background equal to None.

::

    background = None

To specify a single value for background to be used across all observations and instruments, you may supply a single string ('high', 'medium', or 'low') or a single number. The strings correspond to the background levels from the `JWST ETC <https://jwst.etc.stsci.edu/>`_ and are detailed in the :ref:`bkgdrate <bkgdrate>` section of the Input Yaml Parameters page.
If a single number is provided, it is interpreted as the background signal in units of DN/sec/pixel. This signal will be uniform across all pixels.

::

    background = 'high'
    background = 22.2

In order to use a different background in each observation of the proposal (if for example, your proposal will be broken into multiple epochs observed at different times of year), you can use a dictionary. The keys are the observation
numbers in the proposal file. Note that these are three-character strings. The values can then be the same strings or numbers described above.

::

    background = {'001': 'high', '002': 'medium', '003': 22.3}

For finer control, you can use a nested dictionary to specify the background signal in each instrument and observation. The top level keys of the dictionary are the observation numbers from the APT file (again, as 3-character strings).
The values are then instrument names, which are keys into the second level of the dictionary. The values for these keys can be strings or numbers, as above, or in the case of NIRCam, a further nested dictionary that breaks out the
background level by `channel <https://jwst-docs.stsci.edu/near-infrared-camera/nircam-overview#NIRCamOverview-Channels>`_ (shortwave detectors versus longwave detectors).

::

    background = {'001': {'nircam': {'sw': 0.2, 'lw': 0.3},
                          'niriss': 0.4},
                          'fgs': 0.2},
                  '002': {'nircam': {'sw': 'medium', 'lw': 'high'},
                          'niriss': 'low'},
                          'fgs': 'high'},
                  '003': {'nircam': {'sw': 0.75, 'lw': 'high'},
                          'niriss': 0.2}}
                          'fgs': 0.1}}


.. _yam_gen_pav3_inputs:

Roll Angle
++++++++++

Another optional user input to the *yaml_generator* is the `roll angle <https://jwst-docs.stsci.edu/observatory-functionality/jwst-position-angles-ranges-and-offsets#JWSTPositionAngles,Ranges,andOffsets-Referenceangledefinitions>`_ of the telescope. This is often referred to as PAV3 (or V3PA) as it is the position angle of the V3 axis in degrees east of north. We include this as a tunable parameter in Mirage to allow users to explore different orientations for their data, including the use of different roll angles for different observations within their proposals in order to simulate epochs.

In order to use the Mirage default value (roll angle = 0), simply do not provide a roll angle input to the yaml_generator, or explicitly set it to None.

::

    pav3 = None

To specify a single roll angle to be used in all observations, supply a single number.

::

    pav3 = 34.5

In order to simulate epochs and break up your observations, supply a dictionary where the keys are the (3-character string) observation numbers from your APT file, and the values are the roll angles to use for those observations.

::

    pav3 = {'001': 34.5, '002': 154.5, '003': 37.8}


.. _yam_gen_date_inputs:

Observation Dates
+++++++++++++++++

You may also specify the observation date for each observation in your APT file. This may be used along with roll angle to help define epochs in your observations, or simply
to associate a given dataset with a date. **Note that Mirage does not pay attention to dates in any way** other than to save them into the *date-obs* header keyword in the output
files. Mirage does not check that a given roll angle and pointing are physically realizable on a given date. It is up to you to provide realistic values for these paramters
if they are important to you. The `JWST Target Visibility Tools <http://www.stsci.edu/jwst/science-planning/proposal-planning-toolbox/target-visibility-tools>`_ (TVT) are
useful for this. Note that in all cases below, Mirage will use the entered date (along with a default time) as the starting time of the first exposure in the observation.
Mirage keeps track of exposure times and makes some guesses about overheads, and increments the observation time and date for each exposure.

To use the Mirage default for observation date (arbitrarily set to 2021-10-04), you can either not supply any date information, or explicitly use None.

::

    dates = None

To use a single date for all observations, you can give a date string.

::

    dates = '2022-5-25'

To specify a different date for each observation, use a dictionary where the keys are the (3-character string) observation numbers from your APT file, and the values are
the date strings for each.

::

    dates = {'001': '2022-06-25', '002': '2022-11-15', '003': '2023-03-14'}


.. _yam_gen_cr_inputs:

Cosmic Ray Rates
++++++++++++++++

You may also customize the cosmic ray rates applied to Mirage's outputs. There are two aspects of the cosmic ray behavior that can be controlled. The first is the name
of the library of cosmic ray stamp images to use, and the second is a scaling factor that can be applied to that library. The three library options, in are **SUNMAX**,
**SUNMIN**, and **FLARE**. Each library contains a different collection of cosmic ray images, and each has a default cosmic ray rate (cosmic rays per pixel per second)
associated with it. The SUNMIN and SUNMAX labels refer to the solar activity, and the galactic cosmic ray contribution at L2 is reduced at solar maximum compared to solar
minimum.  The FLARE case is for the largest solar flare event on record and corresponds to conditions under which JWST would presumably not be
operating. The table below give the cosmic ray probabilities for the three libraries. The cosmic ray libraries and default probabilties were taken from
`Robberto 2009 <http://www.stsci.edu/files/live/sites/www/files/home/jwst/documentation/technical-documents/_documents/JWST-STScI-001928.pdf>`_.

+-----------+------------------------+
| *Library* |*Cosmic Ray Probability*|
+-----------+------------------------+
|  SUNMAX   |      5.762e-06         |
+-----------+------------------------+
|  SUNMIN   |      1.587e-05         |
+-----------+------------------------+
|  FLARE    |      0.0098729         |
+-----------+------------------------+

The second configurable aspect of the cosmic ray rate is a scaling factor. This is a multiplicative factor that will be applied to the probability from the selected
library in order to determine the final cosmic ray probability.

To use Mirage's default values of the SUNMAX library and a scaling factor of 1.0, simply do not provide any input, or explicitly set the cosmic ray rate to None.

::

    cr = None

To specify a different library and scale from the default, and apply those to all observations in your proposal, provide a dictionary with 'library' and 'scale' keys
set to your desired values.

::

    cr = {'library': 'FLARE', 'scale': 44.0}

In order to use a different cosmic ray library and scaling factor for each observation, create a nested dictionary where the top-level keys are the (3-character string)
observation numbers from your APT file. Each entry should then contain a dictionary with 'library' and 'scale' values.

::

    cr = {'001': {'library': 'FLARE', 'scale': 1.2},
          '002': {'library': 'SUNMIN', 'scale': 5.5},
          '003': {'library': 'SUNMAX', 'scale': 0.1}}


.. _override_reffiles:

JWST Calibration Reference Files
++++++++++++++++++++++++++++++++

Mirage makes use of a handful of the `reference file types <https://jwst-pipeline.readthedocs.io/en/stable/jwst/introduction.html#reference-files>`_ used by the JWST calibration pipeline. This includes the `bad pixel mask <https://jwst-pipeline.readthedocs.io/en/stable/jwst/dq_init/reference_files.html#mask-reffile>`_, `saturation level map <https://jwst-pipeline.readthedocs.io/en/stable/jwst/saturation/reference_files.html#saturation-reffile>`_, `superbias <https://jwst-pipeline.readthedocs.io/en/stable/jwst/superbias/reference_files.html#superbias-reffile>`_, `gain <https://jwst-pipeline.readthedocs.io/en/stable/jwst/references_general/gain_reffile.html#gain-reffile>`_, `interpixel capacitance <https://jwst-pipeline.readthedocs.io/en/stable/jwst/ipc/reference_files.html#ipc-reffile>`_, `linearity correction  <https://jwst-pipeline.readthedocs.io/en/stable/jwst/linearity/reference_files.html#linearity-reffile>`_, `distortion correction <https://jwst-pipeline.readthedocs.io/en/stable/jwst/references_general/distortion_reffile.html#distortion-reffile>`_ and `pixel to pixel flat field <https://jwst-pipeline.readthedocs.io/en/stable/jwst/flatfield/reference_files.html#flat-reference-file>`_ files.

Mirage relies on the `CRDS <https://hst-crds.stsci.edu/static/users_guide/index.html>`_ package from STScI to identify the appropriate reference files for a given exposure. These files are automatically downloaded to the user's machine at one of two times:

1) When running *yaml_generator.py* if the ``reffile_defaults`` option is set to 'crds_full_names' (see :ref:`Run the Yaml Generator <yaml_generator>` for more details).
2) When running any of the three main parts of Mirage: the :ref:`catalog seed generator <source_catalogs>`, :ref:`dark preparation <dark_prep>`, or :ref:`observation generator <obs_generator>`.

.. tip::

    In order to specify a location on your machine to store the downloaded reference files, you must have the CRDS_PATH environment variable set to that location. If this environment variable is not set, Mirage will default to use $HOME/crds_cache/

.. important::

    Due to a limitation of the **CRDS** package, the CRDS_PATH and CRDS_SERVER_URL environment variables must be set BEFORE importing the **CRDS** package. If you are running Mirage using code that imports **CRDS** *or any other packages that import CRDS* (such as Mirage's **dark_prep** module or the **jwst** package, which contains the calibration pipeline) prior to running Mirage, you should explicitly set the environment variables before importing those packages. If you do not, you will get the following error:

    CRDS - ERROR -  (FATAL) CRDS server connection and cache load FAILED.  Cannot continue.
    CRDS - ERROR -  See `https://hst-crds.stsci.edu/docs/cmdline_bestrefs/ <https://hst-crds.stsci.edu/docs/cmdline_bestrefs/>`_ or `https://jwst-crds.stsci.edu/docs/cmdline_bestrefs/ <https://jwst-crds.stsci.edu/docs/cmdline_bestrefs/>`_
    CRDS - ERROR -  for more information on configuring CRDS,  particularly CRDS_PATH and CRDS_SERVER_URL. : [Errno 2] No such file or directory: '$HOME/crds_cache/config/jwst/server_config'

    If you wish to set the environment variables in your code, simply add lines such as these prior to importing **jwst** or **CRDS**:

    os.environ["CRDS_PATH"] = '{}/crds_cache'.format(os.environ.get('HOME'))
    os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"

    **If your code does not import any packages that rely on the CRDS package, then you may safely neglect setting the two environment variables, and Mirage will set them for you prior to importing CRDS.**


It is also possible to specify that Mirage use reference files other than those downloaded from CRDS. In order to do this, you must supply a dictionary of filenames. Due to the varying number of selection criteria needed to uniquely identify the reference file that matches up with a particular exposure, this dictionary is composed of multiple levels of nested dictionaries. Not all possibilities are required. You may specify reference files only for the particular observing modes you are interested in. Any modes in your APT file that are not contained within the dictionary will revert to using the reference files identified by CRDS. Below is a dictionary showing all nesting required for all reference files. Note that the dictionary structure is instrument dependent since a reference file type for different instruments does not necessarily have the same selection criteria.

.. important::
    If you choose to provide your own reference files, it is best to use these same reference files when running the JWST calibration pipeline on the simulated data files produced by Mirage. Not doing this can lead to systematic errors in your calibrated data.

Here is a view of the dictionary structure required when specifying reference files. The easiest way to see the set of allowed values for a particular property (e.g. exposure_type), go to the `JWST Keyword Dictionary <https://mast.stsci.edu/portal/Mashup/Clients/jwkeywords/index.html>`_, and search for the appropriate keyword. Exposure type is EXP_TYPE, filter is FILTER, pupil is PUPIL, detctor_name is DETECTOR, readpattern is READPATT.

::

    override = {'nircam': {'superbias':  {detector_name: {readpattern: 'reffile_name.fits'}},
                           'linearity':  {detector_name: 'reffile_name.fits'},
                           'saturation': {detector_name: 'reffile_name.fits'},
                           'gain':       {detector_name: 'reffile_name.fits'},
                           'distortion': {detector_name: {filter: {exposure_type: 'reffile_name.asdf'}}},
                           'ipc':        {detector_name: 'reffile_name.fits'},
                           'area':       {detector_name: {filter: {pupil: {exposure_type: 'reffile_name.asdf'}}}},
                           'badpixmask': {detector_name: 'reffile_name.fits'},
                           'pixelflat':  {detector_name: {filter: {pupil: 'reffile_name.fits'}}}
                           },
                'niriss': {'superbias':  {readpattern: 'reffile_name.fits'},
                           'linearity':  'reffile_name.fits',
                           'saturation': 'reffile_name.fits',
                           'gain':       'reffile_name.fits',
                           'distortion': {pupil: {exposure_type: 'reffile_name.fits'}},
                           'ipc':        'reffile_name.fits',
                           'area':       {filter: {pupil: {exposure_type: 'reffile_name.asdf'}}},
                           'badpixmask': 'reffile_name.fits',
                           'pixelflat':  {filter: {pupil: 'reffile_name.fits'}}
                           },
                'fgs':    {'superbias':  {detector_name: {readpattern: 'reffile_name.fits'}},
                           'linearity':  {detector_name: 'reffile_name.fits'},
                           'saturation': {detector_name: 'reffile_name.fits'},
                           'gain':       {detector_name: 'reffile_name.fits'},
                           'distortion': {detector_name: {exposure_type: 'reffile_name.fits'}},
                           'ipc':        {detector_name: 'reffile_name.fits'},
                           'area':       {detector_name: 'reffile_name.asdf'},
                           'badpixmask': {detector_name: {exposure_type: 'reffile_name.fits'}},
                           'pixelflat':  {detector_name: {exposure_type: 'reffile_name.fits'}}
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
                           'pixelflat':  {'nrcb5': {'f322w2': {'clear': 'my_favorites/lw_flat.fits',
                                                               'grismr': 'my_favorites/lwR_flat.fits',
                                                               'grismc': 'my_favorites/lwC_flat.fits'
                                                               }
                                                    },
                                          'nrcb4': {'f070w': {'clear': 'my_SW_favs/sw_flat.fits'}}
                                          }
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
                                  cosmic_rays=crs, background=background, roll_angle=pav3,
                                  dates=dates, datatype='raw', use_dateobs_for_background=False,
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




