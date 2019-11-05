.. _reference_files:

Reference Files
===============

**NOTE to users internal to STScI: a collection of Mirage reference files is already present on the STScI network, so there is no need to download them. Please email Bryan Hilbert for the path to use for your MIRAGE_DATA environment variable.**

In order to produce data that is as realistic as possible, Mirage is accompanied by a set of reference files that are used to construct the simulated data.


After installing Mirage, these reference files can be downloaded using the *downloader* module. As the collection of reference files is quite large, users have the option of downloading only certain subsets of the files. For example, if you will only be simulating data from one instrument, you can download only the reference files needed for that instrument. You can also choose to download only certain types of reference files in a given call to the downloader. In this way, you can break up the download into multiple, smaller calls to the downloader if desired. The basic commands to download reference files are shown below. The downloader first creates a list of files to download given the input parameters. Before attempting to download a particular file, the script first checks if the file is already present in the requested local directory. If the file is present, the download of that file is skipped. Details on calling the downloader are provided below.

::

  from mirage.reference_files import downloader
  download_path = '/path/into/which/files/are/downlaoded/'
  downloader.download_reffiles(download_path, instrument='all', dark_type='linearized', skip_darks=False, skip_cosmic_rays=False, skip_psfs=False, skip_grism=False)

The ``instrument`` keyword controls which subset of reference files are downloaded. You can give it the name of a single instrument, a string containing a comma-separated list of instruments, or the string ``all``, which will download reference files for NIRCam, NIRISS, and FGS.

The ``dark_type`` keyword controls which dark current exposures are downloaded. *Mirage* requires linearized dark current exposures when creating simulated data. A user may provide raw dark current files, which *Mirage* will linearize on the fly, or linearized dark current files, which will save processing time. Set this keyword to ``linearized`` to download only the linearized versions of the darks, ``raw`` to download only the raw versions, or ``both`` for both. If omitted by the user, the script will default to downloading only the linearized darks. Note that the darks are by far the largest reference files in the collection (3GB per file, with 5 files per NIRCam detector, 20 files for NIRISS, and 8 files for FGS) and will take the most time to download.

If True, the ``skip_dark`` parameter will cause the script not to download the dark current files for the given instrument. Similarly, the ``skip_cosmic_rays`` and ``skip_psfs`` parameters, if True, will cause the script to skip downloading the cosmic ray library and PSF library, respectively, for the indicated instruments. The default for all three of these parameters is False.

The ``skip_grism`` parameter controls whether the reference files associated with the grisms (needed for WFSS simulations) are downloaded. If False, the grism data will be downloaded.

When called, the function will download the appropriate files from the `STScI Box repository <https://stsci.app.box.com/folder/69205492331>`_, unzip the files, and create the directory structure Mirage expects. It will then remind you to point your MIRAGE_DATA environment variable to the top-level location of these files, so that Mirage knows where to find them. You
may wish to add this definition to your .bashrc or .cshrc file.

For example:

::

	export MIRAGE_DATA="/my_files/jwst/simulations/mirage_data"

CRDS Environment Variables
--------------------------

In addition to the MIRAGE_DATA environment variable, there are two environment variables required by the `Calibration References Data System <https://hst-crds.stsci.edu/static/users_guide/overview.html>`_ (CRDS) that should be set. These will be used when Mirage queries CRDS for the appropriate JWST calibration reference files for the observations being simulated. While setting these environment variables before running Mirage is not strictly required, it is recommended in order to avoid confusion when CRDS is queried. The two environment variables are `CRDS_PATH <https://hst-crds.stsci.edu/static/users_guide/environment.html?#user-local-crds-path>`_ and `CRDS_SERVER_URL <https://hst-crds.stsci.edu/static/users_guide/environment.html?#jwst-ops-server>`_. CRDS_PATH should be set to the directory where you would like JWST calibration reference files to be stored. If not set, Mirage will set this to the CRDS default location of $HOME/crds_cache. CRDS_SERVER_URL must be set to the value below.

::

  export CRDS_PATH=$HOME/crds_cache
  export CRDS_SERVER_URL=https://jwst-crds.stsci.edu


Reference File Contents
-----------------------

There are three main groups of reference files for each instrument: dark current exposures, a PSF library, a cosmic ray library. Note that users can run Mirage with alternate reference files. For example, if a specific science case requires larger PSFs, the user can create a new PSF library and replace the downloaded PSF library, or simply update the ``psf_path`` entry in their Mirage input yaml files to point to the new library. Similarly, to use a non-standard JWST calibration reference file, the appropriate entry in the yaml input file can be changed to point to the new reference file.