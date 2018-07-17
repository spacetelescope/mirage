#! /usr/bin/env python

from setuptools import setup
from setuptools import find_packages
from setuptools.command.test import test as TestCommand
import subprocess
import sys

class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = ['tests']
        self.test_suite = True

    def run_tests(self):
        # import here, cause outside the eggs aren't loaded
        import pytest
        errno = pytest.main(self.test_args)
        sys.exit(errno)



# make sure jwst is available
try:
    import jwst
except ImportError:
    try:
        subprocess.check_call(['git', 'clone',
                               'https://github.com/spacetelescope/jwst.git'])
        sys.path.insert(1, 'jwst')
        # import jwst
    except subprocess.CalledProcessError as e:
        print(e)
        exit(1)



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
    author = 'STScI (Hilbert, Volk, Chambers, Sahlmann et al.)',
    author_email = 'hilbert@stsci.edu',
    url='https://github.com/spacetelescope/mirage',
    keywords = ['astronomy'],
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    packages = find_packages(exclude=["examples","catalogs"]),
    package_data={'mirage': ['tests/test_data/*/*.yaml',
                             'tests/test_data/*/*.list',
                             'tests/test_data/*/*.xml',
                             'tests/test_data/*/*.pointing',
                                ]},
    install_requires=[
        'astropy>=1.2',
        'numpy>=1.9',
        'matplotlib>=1.4.3',
        'lxml>=3.6.4',
        'asdf>=1.2.0',
        # 'scipy>=0.17',
        # 'jwst>=0.9.0',
    ],
    include_package_data = True,
    cmdclass={
        'test': PyTest,
        # 'build_sphinx': BuildSphinx
    },
    )
