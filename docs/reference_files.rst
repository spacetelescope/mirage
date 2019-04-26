.. _reference_files:

Reference Files
===============

**NOTE to users internal to STScI: a collection of Mirage reference files is already present on the STScI network, so there is no need to download them. Please email Bryan Hilbert for the path to use for your MIRAGE_DATA environment variable.**

In order to produce data that is as realistic as possible, Mirage is accompanied by a set of reference files that are used to construct the simulated data.

..
   After installing Mirage, these reference files can be downloaded using the *reference_files* module. As the collection of reference files is quite large, users have the option of downloading only certain subsets of the files. For example, if you will only be simulating data from one instrument, you can download only the reference files needed for that instrument. The basic commands to download reference files are shown below.

  ::

      from mirage.reference_files import downloader
      downloader.download_reffiles(download_path, instrument='all', psf_version='subpixel', dark_type='linearized')

  The ``instrument`` keyword controls which subset of reference files are downloaded. You can give it the name of a single instrument, a string containing a comma-separated list of instruments, or ``all``, which will download reference files for NIRCam, NIRISS, and FGS.

  The ``psf_version`` keyword controls the PSF libraries that are downloaded. The current version of *Mirage* uses libraries composed of many PSFs at various sub-pixel locations. These libraries do not include the effects of distortion. Future PSF libraries used by *Mirage* will be composed of sub-sampled PSFs at various locations across the detector, and will include distortion. In order to download the current PSF libraries, the ``psf_version`` keyword should be set to ``subpixel``, or omitted, in which case the script will default to retrieving the current libraries.

  The ``psf_version`` keyword controls the PSF libraries that are downloaded. Older versions of *Mirage* used libraries composed of many PSFs at various sub-pixel locations. These libraries did not include the effects of distortion. The current PSF libraries used by *Mirage* are composed of sub-sampled PSFs at various locations across the detector, and do include distortion. In order to download the current PSF libraries, the ``psf_version`` keyword should be set to ``gridded``, or omitted, in which case the script will default to retrieving the current libraries.

  The ``dark_type`` keyword controls which dark current exposures are downloaded. *Mirage* requires linearized dark current exposures when creating simulated data. A user may provide raw dark current files, which *Mirage* will linearize on the fly, or linearized dark current files, which will save processing time. Set this keyword to ``linearized`` to download only the linearized versions of the darks, ``raw`` to download only the raw versions, or ``both`` for both. If omitted by the user, the script will default to downloading only the linearized darks. Note that the darks are by far the largest reference files in the collection (3GB per file, with 5 files per NIRCam detector, 20 for NIRISS, and 8 for FGS) and will take the most time to download.

  When called, the function will download the appropriate files from the `STScI Box repository <https://stsci.app.box.com/folder/69205492331>`_, unzip the files, and create the directory structure Mirage expects. It will then remind you to point your MIRAGE_DATA environment variable to the top-level location of these files, so that Mirage knows where to find them.

  For example:

A download function for obtaining a copy of the reference file collection is coming soon. Once you have access to the reference files, you must define the **MIRAGE_DATA** environment variable so that *Mirage* knows where to find the files.

::

	export MIRAGE_DATA="/my_files/jwst/simulations/mirage_data"


Reference File Contents
-----------------------

There are four main groups of reference files for each instrument: dark current exposures, a PSF library, a cosmic ray library, and a set of reference files used by the `JWST calibration pipeline <https://jwst-pipeline.readthedocs.io/en/stable/jwst/introduction.html#reference-files>`_ for calibrating data. Note that users are able to run Mirage with alternate reference files. For example, if a specific science case requires larger PSFs, the user can create a new PSF library and replace the downloaded PSF library, or simply update the ``psf_path`` entry in their Mirage input yaml files to point to the new library.