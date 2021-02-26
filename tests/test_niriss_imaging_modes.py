"""System test of mirage/NIRISS for regular and NRM imaging.

Authors
-------
    - Kevin Volk

Use
---
    >>> pytest -s test_niriss_imaging_modes.py

Description of the test:
------------------------

The test runs mirage for the same NIRISs imaging scene in regular imaging and
in the NRM imaging mode.  There is only one star in the scene of the same
magnitude in both instances.  The code reads the two pointsource output list
files to get the count rates for the two cases and verifies that the proper
scaling has been used.  The count rate ratio should be exactly 0.15/0.84
between the NRM case and the regular imaging case.  Due to limitations on the
precision of the output format, the values differ by a fraction amount of
about 1.e-08 in my test.  Here the threhsold for agreement is a deviaiton of
less than 1.e-06.

"""

import os
import pytest
import numpy

from mirage import imaging_simulator

os.environ['TEST_DATA'] = os.path.join(os.path.dirname(__file__), 'test_data/NIRISS')

# Determine if tests are being run on Github Actions CI
ON_GITHUB = '/home/runner' in os.path.expanduser('~')


@pytest.mark.skipif(ON_GITHUB,
                   reason="Cannot access mirage data in the central storage directory from Github Actions CI.")
def test_niriss_imaging():
    nis = imaging_simulator.ImgSim(offline=True)
    nis.paramfile = os.path.join(os.path.dirname(__file__), 'test_data/NIRISS/niriss_imaging_test.yaml')
    nis.create()
    nis.paramfile = os.path.join(os.path.dirname(__file__), 'test_data/NIRISS/niriss_nrm_test.yaml')
    nis.create()
    value1 = numpy.loadtxt('V88888024002P000000000112o_NIS_F480M_uncal_pointsources.list',usecols=[8,])
    value2 = numpy.loadtxt('V88888024002P000000000112o_NIS_NRM_F480M_uncal_pointsources.list',usecols=[8,])
    fluxratio = value2 / value1

    # The 0.15 factor for the NRM and the 0.84 factor from the CLEARP element
    # are now baked into the PSF from WebbPSF.
    #targetratio = 0.15 / 0.84
    targetratio = 1.0
    deviation = abs(fluxratio/targetratio - 1.)
    assert deviation < 1.e-06
    # clean up the output files in the test directory from this test.
    os.system('/bin/rm V88888024002P000000000112o*')
