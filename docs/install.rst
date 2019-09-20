Installing MIRAGE
=================

The easiest way to get a working installation of Mirage is to use the YAML file in the root directory, ``environment.yml``, to create a MIRaGe-specific conda environment.

You can then use this environment whenever you need to use MIRaGe, or build upon it for any larger projects that require MIRaGe.

Set up a conda environment
--------------------------
First, be sure that your conda installation is up to date. Exit any currently activated environment, if necessary.::

    conda deactivate
    conda update conda

Get a copy of the MIRaGe repository::

    git clone https://github.com/spacetelescope/mirage.git

Move into the mirage directory and use the environment file to create a new environment. Give your environment a name you will recognize (here we choose ``mirage``)::

    cd mirage
    conda env create -f environment.yml -n mirage

Once the environment is built, activate it.::

    conda activate mirage

Move out of the mirage directory before cloning and installing any of the depenencies.::

    cd ../


Install Development-Version Dependencies
----------------------------------------

Unfortunately, there are conflicts with some of MIRaGe's dependencies that have not been addressed in the latest pip and conda versions. So, in order for MIRaGe to work, we need to download and install the development versions of those packages. (This will no longer be necessary once both packages make new releases.)

If you do not have a copy of the ``poppy`` or ``webbpsf`` git repositories, clone them::

    git clone git://github.com/spacetelescope/poppy.git
    git clone git://github.com/spacetelescope/webbpsf.git

If you already have a local copy of these packages, first make sure you have the most up-to-date version. Here we assume that you have previously cloned
the repository and that the ``origin`` remote points to the package. This is git's default behavior. If in doubt, simply delete the local copies and re-clone
using the commands above.::

    cd poppy
    git fetch origin master
    git pull origin master

    cd ../webbpsf
    git fetch origin master
    git pull origin master
    cd ../

Then, install both ``poppy`` and ``webbpsf``::

    pip install -e poppy
    pip install -e webbpsf


Install MIRaGe
--------------

Next, install MIRaGe::

    pip install -e mirage

This ``pip`` command will also always install all required dependencies (though in this case, they were already installed when we built our conda environment). If you want to know what dependencies MIRaGe requires,, the list of packages can
be viewed in the ``install_requires`` part of this repository's `setup.py file <../setup.py>`_.

Install Supporting Packages
---------------------------

If you plan to use MIRaGe to simulate **Wide Field Slitless (WFSS)** observations, you will also need to download two supporting packages:

- `NIRCam_Gsim <https://github.com/npirzkal/NIRCAM_Gsim>`_
- `GRISMCONF <https://github.com/npirzkal/GRISMCONF>`_

For each of these packages, download the repository if you do not have a local copy::

    git clone https://github.com/npirzkal/NIRCAM_Gsim.git
    git clone https://github.com/npirzkal/GRISMCONF.git

Or, if you already have a local copy, make sure you have the most recent version. Here we assume that you have previously cloned
the repository and that the ``origin`` remote points to the package. This is git's default behavior. If in doubt, simply delete the local copies and re-clone
using the commands above.::

    cd ../NIRCAM_Gsim
    git fetch upstream master
    git pull upstream master
    cd ../

    cd ../GRISMCONF
    git fetch upstream master
    git pull upstream master
    cd ../

Then install::

    pip install -e NIRCAM_Gsim
    pip install -e GRISMCONF

.. _ref_file_collection:

Reference Files and MIRAGE_DATA Environment Variable
----------------------------------------------------

In addition to the code itself, there is a set of reference files that accompany Mirage, and are necessary for Mirage to function. These
files include dark current ramps and cosmic ray and PSF libraries.

Instructions for downloading the reference files are provided on the :ref:`reference files <reference_files>` page.



