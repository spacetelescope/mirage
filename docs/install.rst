Installing MIRAGE
=================

The easiest way to get a working installation of Mirage is to use the YAML file in the root directory, ``environment.yml``, to create a MIRaGe-specific conda environment.

You can then use this environment whenever you need to use MIRaGe, or build upon it for any larger projects that require MIRaGe.

Set up a conda environment
--------------------------
First, be sure that your conda installation is up to date. Exit any currently activated environment, if necessary.::

    conda deactivate
    conda update conda

Next, use the environment file to create a new environment. Give your environment a name you will recognize (here we chose just ``mirage``).::

    conda env create -f environment.yml -n mirage

Once the environment is built, activate it.::

    conda activate mirage


Install Development-Version Dependencies
----------------------------------------

Unfortunately, there are conflicts with some of MIRaGe's dependencies that have not been addressed in the latest pip and conda versions. So, in order for MIRaGe to work, we need to download and install the development versions of those packages. (This will no longer be necessary once both packages make new releases.)

To do this, first clone the git repositories for both the ``poppy`` and ``webbpsf`` packages somewhere on your computer::

    git clone git://github.com/spacetelescope/poppy.git
    git clone git://github.com/spacetelescope/webbpsf.git

Then, install both ``poppy`` and ``webbpsf``::

    pip install -e /path/to/poppy/
    pip install -e /path/to/webbpsf/


Install MIRaGe
--------------

Finally, to install the `mirage` package itself, first clone the repository:::

    git clone https://github.com/spacetelescope/mirage.git

Then, install the package:::

    cd mirage
    pip install .

This ``pip`` command will also always install all required dependencies (though in this case, they were already installed when we built our conda environment). If you want to know what dependencies MIRaGe requires,, the list of packages can
be viewed in the ``install_requires`` part of this repository's `setup.py file <../setup.py>`_.

Install Supporting Packages
---------------------------

If you plan to use MIRaGe to simulate **Wide Field Slitless (WFSS)** observations, you will also need to download two supporting packages:

- `NIRCam_Gsim <https://github.com/npirzkal/NIRCAM_Gsim>`_
- `GRISMCONF <https://github.com/npirzkal/GRISMCONF>`_

For each of these packages, clone or download the linked repository, cd into the top-level directory, and then install using:::

    pip install .

.. _reference_files:

Reference Files and MIRAGE_DATA Environment Variable
----------------------------------------------------

In addition to the code itself, there is a set of reference files that accompany Mirage, and are necessary for Mirage to function. These
files include dark current ramps, grism throughput curves, cosmic ray and PSF libraries, as well as more standard
JWST reference files, such as superbias images, linearity correction coefficients, etc.

Instructions for downloading the reference files are provided on the :ref:`reference files <reference_files>` page.



