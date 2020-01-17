Installing MIRAGE
=================
There are two aspects to Mirage installation. First, the software itself must be installed. Once this is complete, there is a set of reference files which
must be downloaded.


Install from Pypi
-----------------

Mirage is now hosted on `Pypi <https://pypi.org/project/mirage/>`_. To install the latest stable version of Mirage, use the commands below. In this example, we create
a conda environment called "mirage" and then install the software into that environment.

::

    conda create -n mirage python=3.6 -y
    conda activate mirage
    pip install mirage

.. tip::
    `Healpy <https://healpy.readthedocs.io/en/latest/>`_, upon which Mirage relies has released different wheels for different versions of Mac OSX. Mirage's setup file currently specifies healpy version 1.12.5,
    which works for MacOSX 10.13 (High Sierra). If this version of healpy does not work for your system, you may need to install a different version of healpy. If this
    is the case, replacing the third command above with `pip install mirage healpy==x.y.z`, where x.y.z is the version number you require, will install that version of healpy.


Install the Development Version
-------------------------------

For those wishing to contribute to the code base, you can install Mirage by cloning and installing the repository. This is only
recommended for those looking to help with development. In general, those wishing only to use Mirage should install the latest stable version from Pypi.


Clone the Mirage repository::

    git clone https://github.com/spacetelescope/mirage.git

Create and activate a new environment. In this example we call the environment "mirage". Then move into the mirage directory, and install Mirage into the new environment::

    conda create -n mirage python=3.6 -y
    conda activate mirage
    cd mirage
    pip install -e .


.. _ref_file_collection:

Reference Files and MIRAGE_DATA Environment Variable
----------------------------------------------------

In addition to the code itself, there is a set of reference files that accompany Mirage, and are necessary for Mirage to function. These
files include dark current ramps and cosmic ray and PSF libraries.

Instructions for downloading the reference files are provided on the :ref:`reference files <reference_files>` page.



