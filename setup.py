#! /usr/bin/env python

from setuptools import setup
from setuptools import find_packages

setup(
    name='mirage',
    version = '1.1',
    description = 'Create simulated JWST data',
    long_description = ('A tool to create simulated NIRCam, NIRISS,'
                        'and FGS exposures'
                        'using an input dark current exposure and a'
                        'noiseless seed image, which can be produced'
                        'from source catalogs. Data can optionally be'
                        'dispersed as well, to simulate wide field'
                        'slitless data files.'),
    author = 'Bryan Hilbert',
    author_email = 'hilbert@stsci.edu',
    keywords = ['astronomy'],
    classifiers = ['Programming Language :: Python'],
    packages = find_packages(exclude=["examples","catalogs"]),
    install_requires = [],
    include_package_data = True
    )
