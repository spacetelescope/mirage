Installing MIRAGE
=================
There are two aspects to Mirage installation. First, the software itself must be installed. Once this is complete, there is a set of reference files which
must be downloaded. The preferred installation method is via :ref:`Pypi <pypi>`, as this is the latest stable version of the software.


.. attention::
    **For those running Mac OSX 10.14:**

    There are several extra steps that must be taken when installing Mirage on a machine running Mac OSX 10.14 (Mojave). These changes are relalted to the option of running calculations in the `Batman <https://github.com/lkreidberg/batman>`_ package in parallel. There are two options for this modified installation, which are described in this `Batman issue on github <https://github.com/lkreidberg/batman/issues/32https://github.com/lkreidberg/batman/issues/32>`_

    1. (The less invasive method) If you do want to make use of parallel processing, you must install LLVM and OpenMP on your machine prior to installing Mirage as described below. See this `StackOverflow issue <https://stackoverflow.com/questions/43555410/enable-openmp-support-in-clang-in-mac-os-x-sierra-mojave>`_ for details.

    2. If you do not wish to use parallel processing within Batman, then you must clone the `Batman <https://github.com/lkreidberg/batman>`_ package, open its *setup.py* file, and remove "-fopenmp". In this case, it is easiest to then :ref:`install Mirage via the environment file <env_file_install>`. Before creating the environment, remove Batman from the environment file. Then create the environment, and then pip install the local copy of Batman.

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
    pip install git+https://github.com/spacetelescope/jwst#0.15.1

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
    pip install git+https://github.com/spacetelescope/jwst@0.14.2

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

.. _ref_file_collection:

Reference Files and MIRAGE_DATA Environment Variable
----------------------------------------------------

In addition to the code itself, there is a set of reference files that accompany Mirage, and are necessary for Mirage to function. These
files include dark current ramps and cosmic ray and PSF libraries.

Instructions for downloading the reference files are provided on the :ref:`reference files <reference_files>` page.



