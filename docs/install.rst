Installing MIRAGE
=================
There are two aspects to Mirage installation. First, the software itself must be installed. Once this is complete, there is a set of reference files which
must be downloaded. The preferred installation method is via :ref:`Pypi <pypi>`, as this is the latest stable version of the software.

.. tip::
    Mirage currently supports python 3.11 and 3.12.

.. _pypi:

Install from Pypi
-----------------

Mirage is hosted on `Pypi <https://pypi.org/project/mirage/>`_. To install the latest stable version of Mirage, use the commands below. In this example, we create
a conda environment called "mirage" and then install the software into that environment.

::

    conda create -n mirage python=3.11 -y
    conda activate mirage
    pip install mirage


Install the Development Version
-------------------------------

For those wishing to use the most up-to-date version of Mirage, you can install directly from github. In general, the version in the master branch should work, although that is not guaranteed. By installing this version, you will have access to updates made since the latest release on Pypi.

The installation procedure is nearly identical to the case of installing from Pypi above. The only difference is that the install command points to github.

::

    conda create -n mirage python=3.11 -y
    conda activate mirage
    pip install git+https://github.com/spacetelescope/mirage


Installation for Developers
---------------------------

For those wishing to contribute to the code base, you can install Mirage by cloning and installing the repository. This is only
recommended for those looking to help with development. In general, those wishing only to use Mirage should install the latest stable version from :ref:`Pypi <pypi>`.


Clone the Mirage repository::

    git clone https://github.com/spacetelescope/mirage.git

Installation can then be done via pip, which uses setup.py, or using the conda environment file that is included in the package.

To install using pip and setup.py:
Create and activate a new environment. In this example we call the environment "mirage". Then move into the mirage directory, and install Mirage into the new environment::

    conda create -n mirage python=3.11 -y
    conda activate mirage
    cd mirage
    pip install .

.. _env_file_install:

Install via Environment File
----------------------------

The Mirage repository also contains environment files, which can be used to create an environment with proper versions of all of Mirage's dependencies. After cloning the Mirage repository, the environment file (located within the top-level directory) can be used via the following commands. The *name* keyword is used to specify that the name of the environment. You can name the environment anything you like.

Create a python 3.11 environment using the environment file, activate the environment, and install mirage::

    cd mirage
    conda env create -f environment_python_3.11.yml
    conda activate mirage_py3.11
    pip install .


There is also an environment file that can be used to create python 3.12 environment::

    cd mirage
    conda env create -f environment_python_3.12.yml
    conda activate mirage_py3.12
    pip install .


.. _batman

Batman versions
---------------

Python 3.12
+++++++++++

As of October 2025, the Batman package does not work when using python 3.12. If you require the Batman package (used to create Time Series data), you should use a python 3.11 environment.


Intel-based Mac installation
++++++++++++++++++++++++++++

As of Oct 2025, the most recent version of the Batman package is 2.5.3. However, this version does not have a wheel that works for Intel-based Macs. It is possible to revert back to version 2.5.2, however this version requires numpy < 2.0, which conflicts with several other dependencies, including the jwst pipeline. Therefore, for those on Intel-based Macs, we recommend omitting the Batman package when creating the environment. If you are on an Intel-based Mac and require the Batman package, you can either create your environment using the environment file for python 3.11 and then manually downgrade numpy, jwst, and Batman, or use the environment_python_3.11_intel_mac.yml to create the environment.

Build failure
+++++++++++++

This is an issue from previous versions of Mac OXS and may not longer be relevent. We keep it here just in case it is still an issue for anyone.

If the `Batman <https://github.com/lkreidberg/batman>`_ package fails to build, the work-around is more complex. Mirage uses the Batman package when simulating imaging and grism Time Series Observations (TSO).

The installation errors are related to supporting Batman's ability to run calculations in parallel. There are two options for modifying the installation, which are described in this `Batman issue on github <https://github.com/lkreidberg/batman/issues/32https://github.com/lkreidberg/batman/issues/32>`_

    1. If you do want to make use of parallel processing (or simply want to try the less invasive installation fix), you must install LLVM and OpenMP on your machine prior to installing Mirage. See this `StackOverflow issue <https://stackoverflow.com/questions/43555410/enable-openmp-support-in-clang-in-mac-os-x-sierra-mojave>`_ for details. If you successfully install these, then you should be able to install Mirage following the instructions in the sections above.


    2. If you do not wish to use parallel processing within Batman, or the option above fails, then you can modify Batman such that it does not use parallel processing. This involves modifying the Batman and Mirage *setup.py* files and install using those. Clone the `Batman <https://github.com/lkreidberg/batman>`_ package, open its *setup.py* file, and remove "-fopenmp". Then you must clone Mirage and remove Batman from Mirage's *environment.yml* and *setup.py* files. Then create the environment using *environment.yml*, pip install the local copy of Batman, and pip install the local copy of Mirage.

    ::

        cd mirage
        conda env create -f environment_python_3.11.yml
        conda activate mirage-py3.11
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



