Installing MIRAGE
=================
There are two aspects to Mirage installation. First, the software itself must be installed. Once this is complete, there is a set of reference files which
must be downloaded. The preferred installation method is via :ref:`Pypi <pypi>`, as this is the latest stable version of the software.


.. attention::
    **For those running Mac OSX 10.14:**

    Some users have reported errors when installing Mirage on their machines running Mac OSX 10.14. If you see installation failures for the synphot or batman packages, they are most likely related to the OpenMP library. See the section below on :ref:`Troubleshooting for Mac OSX 10.14 installtion <osx1014>`.


.. _pypi:

Install from Pypi
-----------------

Mirage is now hosted on `Pypi <https://pypi.org/project/mirage/>`_. To install the latest stable version of Mirage, use the commands below. In this example, we create
a conda environment called "mirage" and then install the software into that environment. After installing Mirage, there are three packages which must be installed separately.
These are: **jwst**, which is the JWST calibration pipeline software, and two packages that help to create dispersed data using the grisms.

::

    conda create -n mirage python=3.6 -y
    conda activate mirage
    pip install healpy==1.12.5
    pip install mirage
    pip install git+https://github.com/npirzkal/GRISMCONF#egg=grismconf
    pip install git+https://github.com/npirzkal/NIRCAM_Gsim#egg=nircam_gsim
    pip install git+https://github.com/spacetelescope/jwst#0.16.2

.. tip::
    Some of Mirage's dependencies rely on `Healpy <https://healpy.readthedocs.io/en/latest/>`_,. Healpy has released different wheels for different versions of Mac OSX. For example, healpy version 1.12.5
    works for MacOSX 10.13 (High Sierra). If the version of healpy above does not work for your system, you may need to install a different version.

.. tip::
    This method installs `webbpsf <https://webbpsf.readthedocs.io/en/latest/>`_ via pip. In this case, you must also `manually download the collection of webbpsf data files <https://webbpsf.readthedocs.io/en/latest/installation.html#installing-the-required-data-files>`_ If you install webbpsf via conda, the data files are downloaded and installed for you.


Install the Development Version
-------------------------------

For those wishing to contribute to the code base, you can install Mirage by cloning and installing the repository. This is only
recommended for those looking to help with development. In general, those wishing only to use Mirage should install the latest stable version from :ref:`Pypi <pypi>`.


Clone the Mirage repository::

    git clone https://github.com/spacetelescope/mirage.git

Installation can then be done via pip, which uses setup.py, or using the conda environment file that is included in the package.

To install using pip and setup.py:
Create and activate a new environment. In this example we call the environment "mirage". Then move into the mirage directory, and install Mirage into the new environment::

    conda create -n mirage python=3.6 -y
    conda activate mirage
    cd mirage
    pip install healpy==1.12.5
    pip install .
    pip install git+https://github.com/npirzkal/GRISMCONF#egg=grismconf
    pip install git+https://github.com/npirzkal/NIRCAM_Gsim#egg=nircam_gsim
    pip install git+https://github.com/spacetelescope/jwst@0.16.2

.. tip::
    Some of Mirage's dependencies rely on `Healpy <https://healpy.readthedocs.io/en/latest/>`_,. Healpy has released different wheels for different versions of Mac OSX. For example, healpy version 1.12.5
    works for MacOSX 10.13 (High Sierra). If the version of healpy above does not work for your system, you may need to install a different version.

.. tip::
    This method installs `webbpsf <https://webbpsf.readthedocs.io/en/latest/>`_ via pip. In this case, you must also `manually download the collection of webbpsf data files <https://webbpsf.readthedocs.io/en/latest/installation.html#installing-the-required-data-files>`_ If you install webbpsf via conda, the data files are downloaded and installed for you.

.. _env_file_install:

**Or, to install using the environment file, again creating an environment called "mirage"**::

    cd mirage
    conda env create -f environment.yml --name mirage python=3.6
    conda activate mirage
    pip install .

.. tip::
    For this latter case, packages are installed via conda. For `webbpsf <https://webbpsf.readthedocs.io/en/latest/installation.html#requirements-installation>`_, this means the data files will be downloaded and installed with the software itself. No manual installation of the data files is necessary.


.. _osx1014:

Troubleshooting for Mac OSX 10.14 installtion
---------------------------------------------

If you have installation errors on your machine running 10.14 (Mojave), try these solutions.

Synphot
+++++++

If the synphot package fails to build, try installing via conda using the conda-forge channel. Do this before installing Mirage, using the command:

    - conda install synphot -c conda-forge

Batman
++++++

If the `Batman <https://github.com/lkreidberg/batman>`_ package fails to build, the work-around is more complex. Mirage uses the Batman package when simulating imaging and grism Time Series Observations (TSO).

The installation errors are related to supporting Batman's ability to run calculations in parallel. There are two options for modifying the installation, which are described in this `Batman issue on github <https://github.com/lkreidberg/batman/issues/32https://github.com/lkreidberg/batman/issues/32>`_

    1. If you do want to make use of parallel processing (or simply want to try the less invasive installation fix), you must install LLVM and OpenMP on your machine prior to installing Mirage. See this `StackOverflow issue <https://stackoverflow.com/questions/43555410/enable-openmp-support-in-clang-in-mac-os-x-sierra-mojave>`_ for details. If you successfully install these, then you should be able to install Mirage following the instructions in the sections above.


    2. If you do not wish to use parallel processing within Batman, or the option above fails, then you can modify Batman such that it does not use parallel processing. This involves modifying the Batman and Mirage *setup.py* files and install using those. Clone the `Batman <https://github.com/lkreidberg/batman>`_ package, open its *setup.py* file, and remove "-fopenmp". Then you must clone Mirage and remove Batman from Mirage's *environment.yml* and *setup.py* files. Then create the environment using *environment.yml*, pip install the local copy of Batman, and pip install the local copy of Mirage.

    ::

        cd mirage
        conda env create -f environment.yml --name mirage python=3.6
        conda activate mirage
        pip install .
        cd ../batman
        pip install .

    3. If you are having installtion problems and will not be creating TSO simulations, you could skip Batman installation altogether. In this case you will still need to clone Mirage and remove Batman from the *environment.yml* and *setup.py* files. Then :ref:`install Mirage via the environment file <env_file_install>`.


.. _ref_file_collection:

Reference Files and MIRAGE_DATA Environment Variable
----------------------------------------------------

In addition to the code itself, there is a set of reference files that accompany Mirage, and are necessary for Mirage to function. These
files include dark current ramps and cosmic ray and PSF libraries.

Instructions for downloading the reference files are provided on the :ref:`reference files <reference_files>` page.



