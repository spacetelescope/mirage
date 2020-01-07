#!/usr/bin/env python
import os
import pkgutil
import subprocess
import sys
from setuptools import setup, find_packages, Extension, Command
from setuptools.command.test import test as TestCommand


# allows you to build sphinx docs from the package
# main directory with "python setup.py build_sphinx"

try:
    from sphinx.cmd.build import build_main
    from sphinx.setup_command import BuildDoc

    class BuildSphinx(BuildDoc):
        """Build Sphinx documentation after compiling C source files"""

        description = 'Build Sphinx documentation'

        def initialize_options(self):
            BuildDoc.initialize_options(self)

        def finalize_options(self):
            BuildDoc.finalize_options(self)

        def run(self):
            build_cmd = self.reinitialize_command('build_ext')
            build_cmd.inplace = 1
            self.run_command('build_ext')
            build_main(['-b', 'html', './docs', './docs/_build/html'])

except ImportError:
    class BuildSphinx(Command):
        user_options = []

        def initialize_options(self):
            pass

        def finalize_options(self):
            pass

        def run(self):
            print('!\n! Sphinx is not installed!\n!', file=sys.stderr)
            exit(1)


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
    description='Create simulated JWST data',
    long_description=('A tool to create simulated NIRCam, NIRISS,'
                      'and FGS exposures'
                      'using an input dark current exposure and a'
                      'noiseless seed image, which can be produced'
                      'from source catalogs. Data can optionally be'
                      'dispersed as well, to simulate wide field'
                      'slitless data files.'),
    author='STScI (Hilbert, Volk, Chambers, Sahlmann et al.)',
    author_email='hilbert@stsci.edu',
    url='https://github.com/spacetelescope/mirage',
    keywords=['astronomy'],
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    packages=find_packages(exclude=["examples"]),
    package_data={'mirage': ['tests/test_data/*/*.yaml',
                             'tests/test_data/*/*.list',
                             'tests/test_data/*/*.xml',
                             'tests/test_data/*/*.pointing',
                             'config/*.*']
                  },
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    install_requires=[
        'asdf>=1.2.0',
        'astropy>=3.2.1',
        'jwst-backgrounds>=1.1.1',
        'lxml>=3.6.4',
        'matplotlib>=1.4.3',
        'numpy',
        'photutils>=0.4.0',
        'pysiaf>=0.6.1'
        'scipy>=0.17',
    ],
    include_package_data=True,
    cmdclass={
        'build_sphinx': BuildSphinx
    },)
