Installing MIRAGE
=================

The easiest way to get a working installation of Mirage is to use one of the conda environment files
in the repository to create a conda environment. There are two environment files in the package:

1. `mirage_environment_linux.yml`
2. `mirage_environment_osx.yml`



Create an environment using the environment file
------------------------------------------------
First, be sure that your conda installation is up to date. Exit the current environment if necessary.

::

	source deactivate
	conda update conda

Next, use the appropriate environment file for your operating system to create a new environment.

::

    conda env create -f mirage_environment_osx.yml

Activate the environment
------------------------

::

    source activate mirage_osx


Install Mirage
--------------

In the top-level Mirage directory:

::

    python setup.py install


Install Supporting Packages
---------------------------

`NIRCam_Gsim <https://github.com/npirzkal/NIRCAM_Gsim>`_ -- Support for Wide Field Slitless observations

`GRISMCONF <https://github.com/npirzkal/GRISMCONF>`_ -- Support for Wide Field Slitless observations

For each of these packages, clone or download the repository, cd into the top-level directory, and then install using:

::

    python setup.py install

.. _reference_files:

Reference Files and MIRAGE_DATA Environment Variable
----------------------------------------------------

There is a set of reference files that accompany Mirage, and are necessary for Mirage to function. These
files include dark current ramps, grism throughput curves, cosmic ray and PSF libraries, as well as more standard
JWST reference files, such as superbias images, linearity correction coefficients, etc.

These reference files can be downloaded at **SOMEWHERE**. And you can set them up **SOMEHOW**

Once the reference files are present, create a new environment variable called MIRAGE_DATA, and point
it to the top level directory of the reference files.

::

	export MIRAGE_DATA="/my_files/jwst/simulations/mirage_data"


