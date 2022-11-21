2.3.0
=====

APT
---

Bugfix in the APT xml reader for MOS programs with no dithers. (#806)


CRDS Reference Files
--------------------

Remove references to the pixel area map from code and documentation, as this is no longer used. (#790)


Packaging
---------

Update dependency list. Drop support for python 3.7. (#783)

Remove python 3.6 environment file. (#784)


Catalogs
--------

Update the import location of Irsa from astroquery.irsa to the new location of astroquery.ipac.irsa (#772)


Dark Current
------------

Correct bug that was causing metadata from raw darks to not be read in correctly. Allow input raw darks with no ZEROFRAME extension as long as they were collected with a RAPID readout pattern. (#762)

Remove JWST calibration pipeline config files, as they are non longer supported by the pipeline. (#791)


Moving Targets
--------------

Updated the MOVING_TARGET_POSITION table from detector -> sci. (#757)

Fix a bug that was resulting in moving target sources that were brighter than they should have been. (#807)


Galaxy and Extended Source Position Angles
------------------------------------------

Fixed a bug that was causing incorrect position angles of galaxy and extended sources when PAV3 was non-zero. (#767)

Fix galaxy/extended source rotation test (#778)


NIRISS SOSS
-----------

Stores the input star and planet data (including the orbital parameters) to the output FITS files so users can make direct comparisons. The data is stored in a new 'STAR' and (when applicable) 'PLANET' extensions with the units and orbital params saved as header keywords. (#773)

Moves hardcoded max_frames variable into SossSim.create keyword args so users can adjust memory allocation. (#775)

Allow a user to generate a custom wavelength solution before generating SOSS simulations (#787)


WFSS
----

Save "seed_cube" pickle files of point source, galaxy, and extended sources. These are then picked up by the disperser. This will allow overlapping sources to be treated correctly when everything is disersed. (#779)


2.2.1
=====

NIRISS Ghosts
-------------

Update to new version of NIRISS ghost gap file. (#751)

Remove the NIRISS ghost gap file from the Mirage repository. Download the file from the niriss_ghost repository whenever a
NIRISS simulation including ghosts is run. (#753)


Moving Targets
--------------

Populate the MOVING_TARGET_POSITION table extension in the output files. (#749)


Output File Headers
-------------------

Populate the PPS_APER header keyword in saved uncal and linear files. This information comes from the APT pointing file. (#752)

Populate the NOUTPUTS header keyword in saved uncal and linear files. This is needed in order to find the correct CRDS reference files in some cases. (#756)



2.2.0
=====


APT Reader
----------

Update NIRSpec MOS APT xml reader to make use of ConfigurationPointing fields, in order to determine order and number of exposures. (#707)

Be sure that the number of primary dithers as read from the xml file is an integer, in order to correctly populate the new NRIMDTPT
header keyword. (#718)

Set the observation mode to "wfss" for NIRISS observations that use the External Calibration template and use a grism. (#737)


Documentation
-------------

Remove mentions of outdated mode name from documentation. "moving_target" is not a mode. Instead, refer to "non-sidereal" tracking. (#740)

Use a newer ReadTheDocs build image in order to avoid build errors associated with new version of jupyter-core. (#745)


Example Notebooks
-----------------

Add missing function definition in imaging mode notebook. (#735)


NIRCam Coronagraphic Modes
--------------------------

Correctly populate the array_name entry of the input yaml file for NIRCam coronagraphic observations. This aperture name is needed when
running coronagraphic simuations.  While Mirage still does not directly support coronagraphic simuations, extended images can be used to
create a coronagraphic scene. (#728)


NIRCam TSO Simulations
----------------------

Fix the location of the target in the TSO apertures. The location of the dispersed source is dictated by the intermediate aperture
location in APT, rather than the imaging aperture. (#722)


NIRISS Simulations
------------------

Add F158M to the list of filters for which ghosts can be added. Update the ghost positions for WFSS simualations,
to better reflect ground testing results. (#714)


Non-sidereal target simulations
-------------------------------

Add PSF convolution keyword when making an extended source image with a moving target simulation. (#706)

Updates to the non-sidereal catalog creation method such that users can input the type of source to use for the tracked target. Remove
the movingTargetConvolveExtended entry from the yaml file since it was not being used. Convolution of any extended taget (sidereal,
non-sidereal, tracked, untracked) with the PSF are all controlled by the PSFConvolveExtended yaml entry. Tracked, non-sidereal sources
are now considered opaque. That is, background sources that fall behind the tracked non-sidereal target will not be added to the seed
image. (#743)


Outputs
-------

Following changes in jwst version 1.2.3, replace the NDITHPTS header keyword with NRIMDTPT. (#705)

Correct the file-splitting calculations when creating data with a non-RAPID readout pattern. (#731)


Packaging
---------

Remove secondary dependencies from environment files and requirements.txt. (#733)


PSF Libraries
-------------

Fix bug that was preventing Mirage from finding the correct gridded PSF library for WLP4 simulations. (#724)

In cases where the total signal in a gridded PSF library is greater than 1.0 but less than 1.1, renormalize the PSF to have a total
signal of 1.0. PSFs with total signals greater than 1.1 will raise an exception. The renormalization is done because in some cases,
WebbPSF will create a PSF with a total signal just slightly greater than 1.0. (#725)


Source Catalogs
---------------

Remove CLEAR and CLEARP from the list of filters read in from the pointing file during source catalog creation. (#710)

Fix bug that was causing inconsistencies in the crossmatching of sources between the 2MASS/Gaia/WISE catalogs. (#726)

Fix a bug that was placing NIRISS filter/pupil names in the wrong entry of the exposure dictionary when using for_proposal()
to create source catalogs. (#730)

Subarray Apertures
------------------

Correct the subarray name associated with the NRCB5_MASKLWB aperture name in APT. (#744)


SOSS Mode Simulations
---------------------

Prevent the SOSS simulator from looking for a non-existing 2nd order with F277W. (#700)

Update SOSS simulator to use astropy's BlackBody model rather than the deprecated BlackBody1D. (#720)


Transmission Image
------------------

Create a pom transmission image composed of all 1's for NIRISS AMI and NIRCam imaging TSO simulations, as is done now for regular imaging simulations. (#701)




2.1.0
=====

Flux calibration
----------------

Updated the gain and filter-based zeropoint values for NIRCam observations. These values were changed in order to stay consistent with new gain and photom reference files delivered to CRDS. The new gain values are ~20% lower than the previous versions, which means the zeropoints have changed by the same amount. (#686)


Logging
-------

Fix logging error so that in NIRISS simulations where ghosts are requested, if the filter used does not support the addition of ghosts, this fact is logged only once, rather than for each source. (#667)


NIRISS Simulations
------------------

Fix bug that was preventing NIRISS External Calibration observations from being read in by the yaml_generator when they were in parallel rather than prime. (#679)

When creating "pom" mode simulations, create a seed image that is larger than the detector. (#681)


Reference Files
---------------

Fix the name of the directory into which the NIRCam longwave linearized darks are saved. Prevously the files were being saved into e.g. BLONG, while later, Mirage was looking for B5. (#692)

Allow the name of the transmission file to be None or 'None' when working in imaging mode. (#695)


Seed Image
----------

Enable the use of IRAC mosaic images as input seed images (#670)


SOSS Simulations
----------------

Allow SOSS simulations to be made using a parameter file or keyword arguments (#671)


Source Catalogs
---------------

Updates to the for_proposal() catalog generation function to populate ra_ref and dec_ref dictionary entries. (#696)


Timer
-----

Stop the timer before moving on to the next source in the source catalog when the current source is completely off the detector. Prior to this, sources that were outside the detector would sometimes raise a timer error when Mirage attempted to start the timer for the next source while that for the current source was still running. (#669)


TSO Simulations
---------------

Update the batman package to version 2.4.8 in order to avoid a problem with version 2.4.7. (#680)


V3 Position Angle
-----------------

Remove the skipped_obs_from_xml keyword value in the call to get_pointing_info within default_obs_v3pa_on_date() and all_obs_v3pa_on_date(). (#674)


2.0.2
=====

NIRISS AMI
----------

Updated the PSF normalization check threshold value to match the new value from an updated NRM mask calculation. Prior to this, Mirage was stopping AMI simulations because the normalization threshold was too low. (#664)


2.0.1
=====

NIRISS AMI
----------

Updated the PSF normalization check to lower the expected total PSF signal in the gridded PSF library for cases where
the NRM is in the beam, as well as imaging cases where CLEARP is used. (#662)


2.0.0
=====


APT xml reader
--------------

Correct the name of the NIRISS grism names in the xml reader for cases where input grisms is BOTH (#491)

Update the xml reader to correctly parse the ShouldDither entry For NIRISS WFSS observations. This entry was introduced in APT 2020.2. (#514)

Fix bug in APT reader where, in the case of an observation with a non-CLEAR pupil wheel value,  the “filter” value in the resulting dictionary was being populated with the full string from the APT file (e.g. “F162M+F150W2”). The filter and pupil wheel value are now correctly separated. (#535)

Adjust read_apt_xml to allow the use of SMALL-GRID-DITHER in the APT file for NIRCam imaging mode (#577)

Add new function dedicated to reading in NIRISS AMI observations (#609)

Add function to read NIRCam and MIRI coronagraphy observation templates from the APT xml file (#615)

Set the pupil value for TA exposures to be ‘NRM’ and set the number of dithers to 1 when it is “None” when reading NIRISS AMI xml files (#627)

Add support for all possible prime/parallel observation templates. Also update such that the APT reader will be able to successfully skip over unsupported observation templates, in order to successfully read APT files that have a mix of supported and unsupported templates. (#637)

Add function to read xml from NIRISS External Calibration APT template (#647)



Computation Time/Efficiency
---------------------------

Add a timer module, to help alert the user to the estimated time remaining when creating a simulation (#519)



Configuration Files
-------------------

Fix typos in NIRISS filter names in filter list file. F159M -> F158M, and F580< -> F480M (#500)

Update the names of the throughput files in the config directory to work with changes associated in #535 (#536)

With the weak lens work done recently, the filter/pupil wheel pairing file in the config directory is no longer used. Remove this file. (#543)



Dark current
------------

Update dark_prep.py to allow the input of non-RAPID dark current exposures with no zeroframe extension. In this case, Mirage will construct
an approximate zeroframe extension and add it in to the exposure. This situation should only occur in the case where older ground-testing
darks that have been converted from FITS Writer format are used. (#470)

Update the simulation wrapper modules (e.g. imaging_simulator) such that if the output from a previous run of dark_prep is already present, the user can input that filename and skip the dark_prep step. (#522)

Adds a keyword to dark_prep that allows users to disable file splitting. This is useful primarily for the creation of linearized dark current files to be used in future calls to Mirage. (#549)

Set the output files from dark_prep to contain arrays of 32-bit floats. Output previously contained 64-bit floats, which was doubling file sizes without adding useful information. Output from the calibration pipeline is also 32-bit floats. (#550)

Outputs from dark_prep step now contain the names of any CRDS reference files used in their creation. In addition, for cases where the calibration pipeline is run when creating a linearized dark, the direct output from the pipeline is no longer saved, as it is not useful for future Mirage runs. The output dark_prep_object is still saved and can be used in future Mirage runs. (#551)



Documentation
-------------

Add an example call to create_catalog.galaxy_background() in the documentation (#503)

Update the workflow instructions for Mirage in the README file. (#559)

Fix incorrect units specified in the doctoring of the magnitude_to_countrate  function. Returned results are ADU/sec (#563)

Update the documentation on source catalog creation to use the updated filter_name and magnitude_system keywords, rather than the deprecated filter and mag_sys (#618)

Add the add_ghosts and PSFConvolveGhosts entries to the example yaml file in the documentation (#626)

Update installation instructions to show how to install the master branch from github without having to clone the repository (#656)



FGS Simulations
---------------

Fix a bug that was causing darks for Guider1 and Guider2 to be mixed in cases of FGS exposures with multiple integrations. (#557)

Add a missing import statement for FGS1 dark search string to the yaml generator (#561)



Flat Fielding
-------------

Separate flat fields and POM transmission files. This allows for the correct application of the flat field to a simulation, at the end of the process. Previously the flat field was applied to the seed image, because some flat fields contained e.g. occulters, which the disperser software needs to know about. It was also leading to multiple applications of the flat field in some cases. Mirage now uses “POM transmission files” to add the effects of occulters to the seed images. (#523)



Galaxy Sources
--------------

Fix bug affecting the position angle of galaxies and extended sources (#480)

Add checks to be sure that the galaxies added via create_catalog.galaxy_background() have realistic radius, ellipticity, and sersic index values (#504)

Updates to galaxy stamp creation and scaling. The new strategy for calculating the size of the galaxy stamp fixes a problem where stamps were previously too large. This cuts down computation time. (#516)



Installation and Reference Files
--------------------------------

Update jwst installation instructions (outdated by future PR though) (#468)

Update installation instructions for Mac OS X 10.14 Mojave. (#475)

Fix bug where both linearized and raw dark reference files are necessary (#478)

Fixed bug in reference file downloader that was preventing raw darks from being downloaded when the user asked for raw+linearized darks (#479)

Reorganize reference file setup such that Grism-related reference files must be cloned from the appropriate GitHub repositories, rather than carrying a copy of them within the library of Mirage reference files (#510)

Update jwst installation instuctions - outdated with more recent changes (#545)

Allow download of a single dark from reference file collection, to allow users to get started more quickly and test Mirage (#579)

Update installation instructions to reflect several dependencies that are now available on Pypi (jwst, grismconf, nircam_gsim) (#581)

Add environment files for python 3.7, 3.8, and 3.9. Change to install most packages via pip rather than conda. Update documentation to indicate that support for python 3.6 will be going away soon. (#620)

Expand all directory names in downloader.py so that all are absolute paths (#625)

Reference files related to SOSS mode support added to the reference file downloader script. (#654)



Logging
-------

Add logging to Mirage. The default is to continue printing messages to the screen, but through the logging module. A log file is also produced. Log files are saved to a mirage_logs subdirectory under the directory containing the simulated data files. In the case of a crash, the mirage_latest.log file in the working directory will contain all of the latest information. (#565)

Fix small typo in one call to the logger in dark_prep (#570)

Fix error in logging statement in function to create Besancon source catalog (#587)

When running a NIRISS simulation and asking for optical ghosts, if the filter/pupil pair does not support the addition of ghosts, then log this fact only once. (#636)



NIRISS Simulations
------------------

For NIRISS POM mode observations, save the oversized seed image to a fits file. Fix a bug where point sources outside the detector but within the oversized region were not being populated in the seed image (#493)

Correct a bug in the conversion of magnitudes to count rates for NIRISS AMI simulations, as well as imaging simulations that use filters that are in the filter wheel (as opposed to the pupil wheel). (#527)

Fix a bug that was preventing the selection of the appropriate gridded PSF library for NIRISS NRM simulations (#529)

Add optional optical ghosts when creating NIRISS exopsures (#597)



Non-sidereal Simulations
------------------------

Fix a bug in non-sidereal exposures where slowly moving targets (<1”/hr) were not being added to the scene. Also, a bug in the scaling applied to all non-sidereal sources was fixed, where previously the scaling was too bright. (#555)

Add the option of an ephemeris_file column in the source catalogs for non-sidereal targets. Mirage can now read in a given Horizons-formatted ephemeris file, and calculate the location of the source for each frame of an observation. The option for users to supply constant velocities in arcseconds/hour or pixels/hour remains. (#564)

Update yaml_generator to properly populate input yaml file entries for non-sidereal observations. (#590)



Output Files
------------

Populate APERNAME keyword in headers of output files. This keyword is not used by the jwst calibration pipeline later, but was requested by
people working on WFSC simulations for their data analyses. (#467)

Implement file splitting for imaging mode observations (#506)

Correctly populate the EXPTYPE fits header keyword in NIRISS AMI simulated data (#541)

Fix a bug in observation generation in cases where the seed image was in multiple file segments, but the dark was in a single file. (#571)

Fix the calculation used to populate the DURATION fits header keyword. Small tweaks to correct the EFFINTTM and EFFEXPTM values. (#576)

Update read_apt_xml to allow for several new string values for PrimaryDither (e.g. 4TIGHT). Also, adjustments were made to the information added to the NUMDTHPT (integer) and NDITHPTS (string) header keywords (#578)

Increase the file splitting threshold value to more closely match that used by DMS. The new threshold is the equivalent of 160 full frame reads. (#631)



Repository
----------

Use dependabot to track dependencies for Mirage (#558)

Add requirements.txt so that dependabot can use it (#560)

Change repository from using Travis to Github actions CI (#634)



SOSS Mode
---------

Add support for NIRISS SOSS mode simulations. This was done by integrating the awesimsoss package within Mirage. (#599)

Fix broadcasting error that was preventing SOSS simulator from working. (#658)



Source Catalogs
---------------

Add convenience functions for computing PA_V3 angle for a given target on a given date. (Note that this requires JWST_GTVT, which is not a Mirage dependency. (#494)

Allow weak lens+filter combinations in source catalog creations as well as simulations (#495)

For any provided source catalog, be sure that Mirage will produce an empty seed image in the case where no sources are present on the detector. (#496)

Fix a bug where jwst_gtvt was failing in cases where the user provided the date as a string. (#537)

For WFSS mode simulations, Mirage will now raise an error if a source catalog has a “magnitude” column rather than a more specific column name (e.g. “nircam_f444w_clear_magnitude”) (#580)

An “index” column is now required in input source catalogs. Mirage will check all source catalogs to be sure that index numbers do not overlap between them, and will raise an error if they do. In catalog generation, users can now specify the starting index to use, in order to easily create multiple catalogs with non overlapping index numbers. (#582)

Allow the results returned from the GAIA catalog search to be bytes or strings, in order to support a change in astropy version 4.2. (#619)

Update get_all_catalogs() to specify that the magnitude system is VEGAMAG. Previously no system was specified and Mirage was defaulting to ABMAG. (#622)

Switch WISE query to use ALLWISE source catalog by default. Allow users to specify using WISE All Sky catalog if desired. (#649)

For observations using the weak lens, the magnitude column name in the source catalog must be e.g. nircam_<filter>_wlp8_magnitude. Previously if this name was not found, Mirage would fall back to look for nircam_<filter>_magnitude, but this would ignore the throughput of the weak lens, which is significant. (#650)



Testing
-------

Skip tests related to 2MASS queries (#469)

Increase timeout limit for Travis tests (#474)

witch the build status badge on Github to Github Actions (#646)



Time Series Observation Simulations
-----------------------------------

Fix a bug for TSO observations where the reference pixels were being improperly masked, which was allowing sources to be present. (#497)

Make background sources optional in grism TSO observations. Previously Mirage would crash if no background source catalog was given. (#512)

Fix a bug that was preventing the addition of 2D dispersed background in grism TSO simulations where there were no background sources present. (#517)

Add a more clear error if someone provides catalog_seed_image() with a yaml file for a WFSS or TSO grism observation, but has grism_source_image within that yaml file set to False. (#539)

Update the location where Mirage looks for grism-related config files. The files (which are copied from NIRCAM_GRISM and NIRISS_GRISM repositories) are now assumed to be in $MIRAGE_DATA/<inst_name>/GRISM_<inst_name>/current/ (#621)

Allow user to supply a 2D array of lightcurves. This means that users will not be limited to creating lightcurves with Batman. (#632)

Update the example TSO notebook to create a source catalog column name nircam_f182m_wlp8_magnitude rather than the nircam_f182m_magnitude that was previously used. (#650)



Wavefront Sensing and Control
-----------------------------

Implemented an optional telescope boresight offset and fixed the tilt scaling issue seen between SW and LW data when using segment_psfs. (#462)

Updated the unstacked mirror and nonnominal PSF notebooks with bug fixes and improvements to support upcoming rehearsal. New notebook dealing with
unstacked mirrors added to the repo. Updates made to get_catalog.py, psf_selection.py, yaml_generator.py and catalog_seed_image.py to support the
unique PSF libraries used by WFSC simulations. (#463)

Allow use of strings for jitter input in the webbpsf call when generating PSF libraries for WFSC. A recent update (beyond 0.9.0) to webbpsf allows
for this input. (#464)

Parallelize calculations to create PSFs for mirror segments. (#473)

Update xml reader and yaml generator to include FGS exposures and PSFs for WFSC Global Alignment observations. (#488)

Update the observation dictionary in the yaml_generator to specify FGS apertures for the FGS exposures in WFSC Global Alignment observations (#492)

In WFSC observations where different PSFs are used for each mirror segment, correctly normalize the PSFs based on the area of the primary in each segment. (#546)

Add basic support for DHS simulations in order to allow simulations for OTE Coarse Phasing. In this case, the DHS sources can only be input as extended sources. (#572)

Enables boresight offsets for full-pupil images. (#595)

Populate subpixel dither type in input yaml files for the case of an WFSC observation when using APT 2020.5. This entry has been removed from APT outputs in this version. (#602)

Check for aperture overrides before calculating starting times, to fix a bug affecting FGS WFSC observations with aperture overrides to NIRCam apertures. (#603)

Add new, more direct way of calculating the position of segment PSFs, making use of a set of optional FITS keywords that record the piston, tip, tilt hexike values directly from within the WebbPSF calculation (#607)



WFSS Simulations
----------------

Save the dispersed background image in a WFSS observation to a file. (#490)

Skip rescaling the spectra for sources that are not present in the ascii source catalog. This is to help support the case where sources in the hdf5 file are spread amongst multiple ascii catalog files. (#616)

Update the location where Mirage looks for grism-related config files. The files (which are copied from NIRCAM_GRISM and NIRISS_GRISM repositories) are now assumed to be in $MIRAGE_DATA/<inst_name>/GRISM_<inst_name>/current/ (#621)

Fix spelling error in get_1d_background_spectrum() that was affecting NIRISS. (#643)

Fix a bug that was not passing user-input segmentation map threshold values through to be used when generating the segmentation map. (#655)



Yaml file updates
-----------------

Fix bug in yaml_generator that was creating incorrect filenames in cases where observation numbers in the APT file were not monotonically increasing. (#482)

Add a user-settable parameter that controls the signal rate threshold for adding pixels to the segmentation map. The default value is 0.031 ADU/sec, based on tests with WFSS exposures. (#507)

Make CRDS-hosted reference file entries in the input yaml files optional. Any entry not present in a yaml file will be set to ‘crds’, in which case Mirage will query CRDS to find the appropriate file. (#513)

Update yaml_generator to populate the “tracking” parameter with “non-sidereal” when a non-sidereal target is specified. (#531)

In cases where the input yaml file contains a colon (e.g in the observation name), Mirage creates a copy of the file and removes the colon so that it can be correctly read in. This PR fixes a bug that was only allowing yaml files in the current directory to go through this process. (#532)

Remove the limits on the number of allowed groups per integration. These rules are fully contained and enforced in APT. Better to rely on those than on the simplified case that was used by Mirage. (#567)

Update yaml_generator to properly populate input yaml file entries for non-sidereal observations. (#590)

Compare version of the PRD in the environment to that used to create the input APT file. (#594)



1.3.3
=====

APT Pointing File
-----------------

Bug fix such that the only Target Acquisition observations that are read in are those for NIRCam TSO observations.

Header Keywords
---------------

Corrected schema to populate the XOFFSET and YOFFSET header keywords (#454)

Reference Files
---------------

Fix bug in downloader that was preventing NIRISS darks from being downloaded (#450)


1.3.2
=====

Gain
----

Added a missing import statement for MEAN_GAIN_VALUES in the grism_tso_simulator

Segmentation Map
----------------

Fixed a bug that was causing create_seed_image to crash when updating the segmentation map for extended sources

Grism TSO plots
---------------

Removed call to an unused module in the TSO example notebook. This call was causing the notebook's plotting function to fail


1.3.1
=====

Dependencies
------------

Added batman-package as a dependency. This is used when creating TSO data.


1.3
===

Installation
------------

setup.py has been modified to support installation via pip and Pypi. Installation documentation has been updated to describe the new process.


Gain Values
-----------

Update observation_generator.py, wfss_simulator.py, grism_tso_simualtor.py to use the mean gain value stored in utils/constants.py rather than the values in the gain reference file when translating the dispersed seed image from units of e-/sec to ADU/sec.

Flat Field
----------

Seed images are now multiplied by the flat field reference file rather than the pixel area map reference file in order to get the surface brightnesses correct. Or more simply, since one of the steps in the JWST calibration pipeline is to divide by the flat field, we must multiply by the flat field when creating the data. See #430. For imaging/time series modes, the flat field is multiplied into the seed image. For NIRCam WFSS mode, the flat field is multiplied in to the dispsersed seed image. For NIRISS WFSS, the flat field is multiplied in to the seed image prior to dispersing. This is because the flat field reference file contains both pixel-to-pixel differences in respsonse, as well as images of the occulting spots, which are in the optical train. Ideally the occulting spots would be multiplied into the seed image prior to dispersing, and then the pixel-to-pixel flat would be multiplied into the seed image after dispersing. Unfortunately these two effects are mixed in the flat field reference file and cannot be separated. This will have some implications for calibrated data products.


Backgrounds
-----------

Update the calculation of background signals to better match the values calcualted by the ETC. Values generally are within 10% of those from the ETC, although there are some filters/pointings/levels where the values differ by up to 20%. #430

For the purposes of calculating the background signal in NIRISS WFSS simulations, the system throughput is simply set to 80% of the throughput in the imaging mode (with appropriate filter).

Fixed bug in the Grism TSO simulator where the background signal was being added twice.

Changed the code so that the Grism TSO simulator works in the case where no background source catalogs are provided. In this case, a dummy background point source catalog is generated, as the calculation and addition of background is done using the background sources.

calculate_background function moved into backgrounds.py so that it can be more easily used by modules other than catalog_seed_image.

Besancon Model Query
--------------------

Code relating to the production of Besancon model source catalogs has been updated to reflect the new workflow for querying and retrieving data. This is now a 2-step process. Users must create an account on the Besancon model website. Queries can then be submitting using the `catalogs.create_catalog.besancon` function. The user must then wait for an email which contains a link to download the resulting catalog. Conversion of this catalog to Mirage-format can then proceed. See the `Catalog_Generation_Tools.ipynb` notebook for details.

Non-sidereal
------------

Segmentation map addition bug corrected. Example notebook and input yaml file updated.


1.2.2
=====

Versioning
----------

Update package versioning to be done with setuptools-scm rather than relic.


1.2.1
=====

TSO Modes
---------

- Updated documentation on readthedocs with information on TSO mode work


1.2
===

TSO Modes
---------

- Add the ability to simulate both grism and imaging time series observations for NIRCam. Example notebook included.


1.1.5
=====

PSF Selection
-------------

- Fix bug in PSF library selection code for observations using one of NIRCam's filters present in the pupil wheel. The bug was preventing the correct library file from being found. (#420)


1.1.4
=====

WCS keywords
------------

- Correct the input RA and Dec used to calculate the values of the PC matrix. Remove the calculation of CRVAL1,2 from set_telescope_pointing.py since it is already done in observation_generator.py (#419)


1.1.3
=====

Yaml Generator
--------------

- Update generator to produce yaml files only for the detectors used with a given aperture. e.g. SUB400P with the NIRCam B module only uses NIRCam B1 and B5 detectors. With this update,
yaml files will only be produced for B1 and B5, whereas previously yaml files were generated for all 5 B module detectors. This change only affects NIRCam.


1.1.2
=====

WFSS
----

- Update functionality for rescaling input spectra to desired magnitude in given instrument/filter. Rescaling is now done via synphot's renormalize() function in the prpoper photon-weighted units. (#412)

Catalogs
--------

- Change photometric system in catalog output from 2MASS query from ABmag to Vegamag (#415)

Seed Image
----------

- Remove filter substring from seed image output file name in the case of FGS simulations (#415)


1.1.1
=====

WFSS
----

- Update background scaling calcultions. NIRISS scales pre-existing background image. NIRCam creates image from jwst_background-provided date or level [#399]
