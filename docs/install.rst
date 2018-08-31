Installing MIRAGE
=================

The easiest way to get a working installation of Mirage is to use the `mirage_environment.yml` file in the repository to create a conda environment.

Create an environment from the environment file
-----------------------------------------------

::

    conda env create -f mirage_environment.yml

Activate the environment
------------------------

::

    source activate mirage


Install Mirage
--------------

In the top-level Mirage directory:
::
    python setup.py install

Install supporting packages
---------------------------
`NIRCam_Gsim <https://github.com/npirzkal/NIRCAM_Gsim>`_ -- Support for Wide Field Slitless observations

`GRISMCONF <https://github.com/npirzkal/GRISMCONF>`_ -- Support for Wide Field Slittless observations

For each of these packages, clone or download the repository, cd into the top-level directory, and then install using:
::
    python setup.py install



