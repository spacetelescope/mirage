#! /usr/bin/env python

"""Download the reference files needed to run Mirage. Extract and unzip
files, and place into the proper directory. Inform the user how to set
the MIRAGE_DATA encironment variable.
"""
import os
import requests
import shutil
import tarfile
import gzip

from mirage.utils.utils import ensure_dir_exists


# Cosmic ray libraries
NIRCAM_CR_LIBRARY_URL = ['https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/cosmic_ray_library/mirage_nircam_cr_library.tar.gz']
NIRISS_CR_LIBRARY_URL = ['https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/cosmic_ray_library/mirage_niriss_cr_library.tar.gz']
FGS_CR_LIBRARY_URL = ['https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/fgs/cosmic_ray_library/mirage_fgs_cr_library.tar.gz']

# Gridded PSF libraries
NIRCAM_GRIDDED_PSF_URLS = ['https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/psf_libraries/gridded_psf_libraries/nircam_psf_wings_library.tar.gz',
                           'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/psf_libraries/gridded_psf_libraries/nircam_Amodule_det12_gridded_psf_library.tar.gz',
                           'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/psf_libraries/gridded_psf_libraries/nircam_Amodule_det345_gridded_psf_library.tar.gz',
                           'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/psf_libraries/gridded_psf_libraries/nircam_Bmodule_det12_gridded_psf_library.tar.gz',
                           'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/psf_libraries/gridded_psf_libraries/nircam_Bmodule_det345_gridded_psf_library.tar.gz',
                           ]
NIRISS_GRIDDED_PSF_URLS = ['https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/psf_libraries/gridded_psf_libraries/niriss_gridded_psf_library.tar.gz']
FGS_GRIDDED_PSF_URLS = ['https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/fgs/psf_libraries/gridded_psf_libraries/fgs_gridded_psf_library.tar.gz']

# Grism-related files
NIRCAM_GRISM_FILES = ['https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/grism/nircam_grism_library.tar.gz']
NIRISS_GRISM_FILES = ['https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/grism/niriss_grism_library.tar.gz']

# Dark current files
NIRCAM_RAW_DARK_URLS = ['https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/A1/NRCNRCA1-DARK-60082202011_1_481_SE_2016-01-09T00h03m58_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/A1/NRCNRCA1-DARK-60090213141_1_481_SE_2016-01-09T02h53m12_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/A1/NRCNRCA1-DARK-60090604481_1_481_SE_2016-01-09T06h52m47_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/A1/NRCNRCA1-DARK-60091005411_1_481_SE_2016-01-09T10h56m36_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/A1/NRCNRCA1-DARK-60091434481_1_481_SE_2016-01-09T15h50m45_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/A2/NRCNRCA2-DARK-60082224241_1_482_SE_2016-01-09T00h10m36_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/A2/NRCNRCA2-DARK-60090235001_1_482_SE_2016-01-09T04h17m03_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/A2/NRCNRCA2-DARK-60090635511_1_482_SE_2016-01-09T07h05m19_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/A2/NRCNRCA2-DARK-60091030561_1_482_SE_2016-01-09T11h03m17_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/A2/NRCNRCA2-DARK-60091457131_1_482_SE_2016-01-09T15h50m45_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/A3/NRCNRCA3-DARK-60082245481_1_483_SE_2016-01-09T00h04m26_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/A3/NRCNRCA3-DARK-60090321241_1_483_SE_2016-01-09T04h17m10_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/A3/NRCNRCA3-DARK-60090656591_1_483_SE_2016-01-09T07h31m27_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/A3/NRCNRCA3-DARK-60091052561_1_483_SE_2016-01-09T11h28m06_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/A3/NRCNRCA3-DARK-60091522581_1_483_SE_2016-01-09T16h30m34_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/A4/NRCNRCA4-DARK-60082307391_1_484_SE_2016-01-09T00h04m08_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/A4/NRCNRCA4-DARK-60090259591_1_484_SE_2016-01-09T04h16m50_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/A4/NRCNRCA4-DARK-60090720041_1_484_SE_2016-01-09T07h58m26_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/A4/NRCNRCA4-DARK-60091117441_1_484_SE_2016-01-09T11h52m23_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/A4/NRCNRCA4-DARK-60091548131_1_484_SE_2016-01-09T16h30m50_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/A5/NRCNRCALONG-DARK-60082329041_1_485_SE_2016-01-09T00h04m16_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/A5/NRCNRCALONG-DARK-60090344021_1_485_SE_2016-01-09T04h16m42_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/A5/NRCNRCALONG-DARK-60090746381_1_485_SE_2016-01-09T08h21m48_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/A5/NRCNRCALONG-DARK-60091140151_1_485_SE_2016-01-09T14h23m49_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/A5/NRCNRCALONG-DARK-60091611271_1_485_SE_2016-01-09T17h16m35_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/B1/NRCNRCB1-DARK-60082356471_1_486_SE_2016-01-09T02h47m00_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/B1/NRCNRCB1-DARK-60090405201_1_486_SE_2016-01-09T05h33m56_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/B1/NRCNRCB1-DARK-60090807581_1_486_SE_2016-01-09T08h48m11_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/B1/NRCNRCB1-DARK-60091205311_1_486_SE_2016-01-09T14h30m08_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/B1/NRCNRCB1-DARK-60091636021_1_486_SE_2016-01-09T17h16m13_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/B2/NRCNRCB2-DARK-60090021181_1_487_SE_2016-01-09T02h51m54_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/B2/NRCNRCB2-DARK-60090427541_1_487_SE_2016-01-09T05h33m14_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/B2/NRCNRCB2-DARK-60090830131_1_487_SE_2016-01-09T08h59m50_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/B2/NRCNRCB2-DARK-60091230011_1_487_SE_2016-01-09T14h23m47_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/B2/NRCNRCB2-DARK-60091735131_1_487_SE_2016-01-09T18h09m45_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/B3/NRCNRCB3-DARK-60090043151_1_488_SE_2016-01-09T02h53m21_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/B3/NRCNRCB3-DARK-60090451471_1_488_SE_2016-01-09T05h33m25_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/B3/NRCNRCB3-DARK-60090852451_1_488_SE_2016-01-09T09h35m03_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/B3/NRCNRCB3-DARK-60091254111_1_488_SE_2016-01-09T14h23m58_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/B3/NRCNRCB3-DARK-60091757401_1_488_SE_2016-01-09T18h40m55_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/B4/NRCNRCB4-DARK-60090118261_1_489_SE_2016-01-09T02h46m53_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/B4/NRCNRCB4-DARK-60090513431_1_489_SE_2016-01-09T05h57m50_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/B4/NRCNRCB4-DARK-60090914351_1_489_SE_2016-01-09T09h52m02_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/B4/NRCNRCB4-DARK-60091316411_1_489_SE_2016-01-09T14h23m38_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/B4/NRCNRCB4-DARK-60091822061_1_489_SE_2016-01-09T18h53m02_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/B5/NRCNRCBLONG-DARK-60090141241_1_490_SE_2016-01-09T02h46m50_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/B5/NRCNRCBLONG-DARK-60090535381_1_490_SE_2016-01-09T06h17m51_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/B5/NRCNRCBLONG-DARK-60090939281_1_490_SE_2016-01-09T10h22m25_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/B5/NRCNRCBLONG-DARK-60091338491_1_490_SE_2016-01-09T15h46m43_level1b_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/raw/B5/NRCNRCBLONG-DARK-60101408431_1_490_SE_2016-01-10T15h01m09_level1b_uncal.fits.gz']
NIRCAM_LINEARIZED_DARK_URLS = ['https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/A1/Linearized_Dark_and_SBRefpix_NRCNRCA1-DARK-60082202011_1_481_SE_2016-01-09T00h03m58_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/A1/Linearized_Dark_and_SBRefpix_NRCNRCA1-DARK-60090213141_1_481_SE_2016-01-09T02h53m12_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/A1/Linearized_Dark_and_SBRefpix_NRCNRCA1-DARK-60090604481_1_481_SE_2016-01-09T06h52m47_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/A1/Linearized_Dark_and_SBRefpix_NRCNRCA1-DARK-60091005411_1_481_SE_2016-01-09T10h56m36_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/A1/Linearized_Dark_and_SBRefpix_NRCNRCA1-DARK-60091434481_1_481_SE_2016-01-09T15h50m45_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/A2/Linearized_Dark_and_SBRefpix_NRCNRCA2-DARK-60082224241_1_482_SE_2016-01-09T00h10m36_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/A2/Linearized_Dark_and_SBRefpix_NRCNRCA2-DARK-60090235001_1_482_SE_2016-01-09T04h17m03_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/A2/Linearized_Dark_and_SBRefpix_NRCNRCA2-DARK-60090635511_1_482_SE_2016-01-09T07h05m19_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/A2/Linearized_Dark_and_SBRefpix_NRCNRCA2-DARK-60091030561_1_482_SE_2016-01-09T11h03m17_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/A2/Linearized_Dark_and_SBRefpix_NRCNRCA2-DARK-60091457131_1_482_SE_2016-01-09T15h50m45_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/A3/Linearized_Dark_and_SBRefpix_NRCNRCA3-DARK-60082245481_1_483_SE_2016-01-09T00h04m26_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/A3/Linearized_Dark_and_SBRefpix_NRCNRCA3-DARK-60090321241_1_483_SE_2016-01-09T04h17m10_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/A3/Linearized_Dark_and_SBRefpix_NRCNRCA3-DARK-60090656591_1_483_SE_2016-01-09T07h31m27_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/A3/Linearized_Dark_and_SBRefpix_NRCNRCA3-DARK-60091052561_1_483_SE_2016-01-09T11h28m06_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/A3/Linearized_Dark_and_SBRefpix_NRCNRCA3-DARK-60091522581_1_483_SE_2016-01-09T16h30m34_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/A4/Linearized_Dark_and_SBRefpix_NRCNRCA4-DARK-60082307391_1_484_SE_2016-01-09T00h04m08_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/A4/Linearized_Dark_and_SBRefpix_NRCNRCA4-DARK-60090259591_1_484_SE_2016-01-09T04h16m50_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/A4/Linearized_Dark_and_SBRefpix_NRCNRCA4-DARK-60090720041_1_484_SE_2016-01-09T07h58m26_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/A4/Linearized_Dark_and_SBRefpix_NRCNRCA4-DARK-60091117441_1_484_SE_2016-01-09T11h52m23_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/A4/Linearized_Dark_and_SBRefpix_NRCNRCA4-DARK-60091548131_1_484_SE_2016-01-09T16h30m50_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/A5/Linearized_Dark_and_SBRefpix_NRCNRCALONG-DARK-60082329041_1_485_SE_2016-01-09T00h04m16_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/A5/Linearized_Dark_and_SBRefpix_NRCNRCALONG-DARK-60090344021_1_485_SE_2016-01-09T04h16m42_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/A5/Linearized_Dark_and_SBRefpix_NRCNRCALONG-DARK-60090746381_1_485_SE_2016-01-09T08h21m48_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/A5/Linearized_Dark_and_SBRefpix_NRCNRCALONG-DARK-60091140151_1_485_SE_2016-01-09T14h23m49_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/A5/Linearized_Dark_and_SBRefpix_NRCNRCALONG-DARK-60091611271_1_485_SE_2016-01-09T17h16m35_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/B1/Linearized_Dark_and_SBRefpix_NRCNRCB1-DARK-60082356471_1_486_SE_2016-01-09T02h47m00_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/B1/Linearized_Dark_and_SBRefpix_NRCNRCB1-DARK-60090405201_1_486_SE_2016-01-09T05h33m56_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/B1/Linearized_Dark_and_SBRefpix_NRCNRCB1-DARK-60090807581_1_486_SE_2016-01-09T08h48m11_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/B1/Linearized_Dark_and_SBRefpix_NRCNRCB1-DARK-60091205311_1_486_SE_2016-01-09T14h30m08_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/B1/Linearized_Dark_and_SBRefpix_NRCNRCB1-DARK-60091636021_1_486_SE_2016-01-09T17h16m13_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/B2/Linearized_Dark_and_SBRefpix_NRCNRCB2-DARK-60090021181_1_487_SE_2016-01-09T02h51m54_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/B2/Linearized_Dark_and_SBRefpix_NRCNRCB2-DARK-60090427541_1_487_SE_2016-01-09T05h33m14_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/B2/Linearized_Dark_and_SBRefpix_NRCNRCB2-DARK-60090830131_1_487_SE_2016-01-09T08h59m50_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/B2/Linearized_Dark_and_SBRefpix_NRCNRCB2-DARK-60091230011_1_487_SE_2016-01-09T14h23m47_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/B2/Linearized_Dark_and_SBRefpix_NRCNRCB2-DARK-60091735131_1_487_SE_2016-01-09T18h09m45_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/B3/Linearized_Dark_and_SBRefpix_NRCNRCB3-DARK-60090043151_1_488_SE_2016-01-09T02h53m21_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/B3/Linearized_Dark_and_SBRefpix_NRCNRCB3-DARK-60090451471_1_488_SE_2016-01-09T05h33m25_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/B3/Linearized_Dark_and_SBRefpix_NRCNRCB3-DARK-60090852451_1_488_SE_2016-01-09T09h35m03_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/B3/Linearized_Dark_and_SBRefpix_NRCNRCB3-DARK-60091254111_1_488_SE_2016-01-09T14h23m58_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/B3/Linearized_Dark_and_SBRefpix_NRCNRCB3-DARK-60091757401_1_488_SE_2016-01-09T18h40m55_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/B4/Linearized_Dark_and_SBRefpix_NRCNRCB4-DARK-60090118261_1_489_SE_2016-01-09T02h46m53_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/B4/Linearized_Dark_and_SBRefpix_NRCNRCB4-DARK-60090513431_1_489_SE_2016-01-09T05h57m50_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/B4/Linearized_Dark_and_SBRefpix_NRCNRCB4-DARK-60090914351_1_489_SE_2016-01-09T09h52m02_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/B4/Linearized_Dark_and_SBRefpix_NRCNRCB4-DARK-60091316411_1_489_SE_2016-01-09T14h23m38_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/B4/Linearized_Dark_and_SBRefpix_NRCNRCB4-DARK-60091822061_1_489_SE_2016-01-09T18h53m02_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/B5/Linearized_Dark_and_SBRefpix_NRCNRCBLONG-DARK-60090141241_1_490_SE_2016-01-09T02h46m50_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/B5/Linearized_Dark_and_SBRefpix_NRCNRCBLONG-DARK-60090535381_1_490_SE_2016-01-09T06h17m51_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/B5/Linearized_Dark_and_SBRefpix_NRCNRCBLONG-DARK-60090939281_1_490_SE_2016-01-09T10h22m25_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/B5/Linearized_Dark_and_SBRefpix_NRCNRCBLONG-DARK-60091338491_1_490_SE_2016-01-09T15h46m43_uncal.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/darks/linearized/B5/Linearized_Dark_and_SBRefpix_NRCNRCBLONG-DARK-60101408431_1_490_SE_2016-01-10T15h01m09_uncal.fits.gz']

NIRISS_RAW_DARK_URLS = ['https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/raw/NISNIRISSDARK-153451235_11_496_SE_2015-12-11T16h05m20_dms_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/raw/NISNIRISSDARK-153451235_12_496_SE_2015-12-11T16h23m51_dms_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/raw/NISNIRISSDARK-153451235_13_496_SE_2015-12-11T16h42m52_dms_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/raw/NISNIRISSDARK-153451235_14_496_SE_2015-12-11T17h01m50_dms_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/raw/NISNIRISSDARK-153451235_15_496_SE_2015-12-11T17h20m40_dms_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/raw/NISNIRISSDARK-153451235_16_496_SE_2015-12-11T17h40m30_dms_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/raw/NISNIRISSDARK-153451235_17_496_SE_2015-12-11T17h59m52_dms_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/raw/NISNIRISSDARK-153451235_18_496_SE_2015-12-11T18h16m31_dms_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/raw/NISNIRISSDARK-153451235_19_496_SE_2015-12-11T18h36m32_dms_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/raw/NISNIRISSDARK-153451235_20_496_SE_2015-12-11T18h53m52_dms_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/raw/NISNIRISSDARK-172500017_13_496_SE_2017-09-07T04h48m22_dms_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/raw/NISNIRISSDARK-172500017_14_496_SE_2017-09-07T05h06m42_dms_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/raw/NISNIRISSDARK-172500017_15_496_SE_2017-09-07T05h28m22_dms_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/raw/NISNIRISSDARK-172500017_16_496_SE_2017-09-07T05h47m42_dms_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/raw/NISNIRISSDARK-172500017_17_496_SE_2017-09-07T06h09m02_dms_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/raw/NISNIRISSDARK-172500017_18_496_SE_2017-09-07T06h29m12_dms_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/raw/NISNIRISSDARK-172500017_19_496_SE_2017-09-07T06h49m52_dms_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/raw/NISNIRISSDARK-172500017_20_496_SE_2017-09-07T07h09m22_dms_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/raw/NISNIRISSDARK-172500017_21_496_SE_2017-09-07T07h29m52_dms_uncal.fits.gz',
                        'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/raw/NISNIRISSDARK-172500017_22_496_SE_2017-09-07T07h50m32_dms_uncal.fits.gz']
NIRISS_LINEARIZED_DARK_URLS = ['https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/linearized/Linearized_Dark_and_SBRefpix_NISNIRISSDARK-153451235_11_496_SE_2015-12-11T16h05m20_dms_uncal_linear_dark_prep_object.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/linearized/Linearized_Dark_and_SBRefpix_NISNIRISSDARK-153451235_12_496_SE_2015-12-11T16h23m51_dms_uncal_linear_dark_prep_object.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/linearized/Linearized_Dark_and_SBRefpix_NISNIRISSDARK-153451235_13_496_SE_2015-12-11T16h42m52_dms_uncal_linear_dark_prep_object.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/linearized/Linearized_Dark_and_SBRefpix_NISNIRISSDARK-153451235_14_496_SE_2015-12-11T17h01m50_dms_uncal_linear_dark_prep_object.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/linearized/Linearized_Dark_and_SBRefpix_NISNIRISSDARK-153451235_15_496_SE_2015-12-11T17h20m40_dms_uncal_linear_dark_prep_object.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/linearized/Linearized_Dark_and_SBRefpix_NISNIRISSDARK-153451235_16_496_SE_2015-12-11T17h40m30_dms_uncal_linear_dark_prep_object.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/linearized/Linearized_Dark_and_SBRefpix_NISNIRISSDARK-153451235_17_496_SE_2015-12-11T17h59m52_dms_uncal_linear_dark_prep_object.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/linearized/Linearized_Dark_and_SBRefpix_NISNIRISSDARK-153451235_18_496_SE_2015-12-11T18h16m31_dms_uncal_linear_dark_prep_object.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/linearized/Linearized_Dark_and_SBRefpix_NISNIRISSDARK-153451235_19_496_SE_2015-12-11T18h36m32_dms_uncal_linear_dark_prep_object.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/linearized/Linearized_Dark_and_SBRefpix_NISNIRISSDARK-153451235_20_496_SE_2015-12-11T18h53m52_dms_uncal_linear_dark_prep_object.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/linearized/Linearized_Dark_and_SBRefpix_NISNIRISSDARK-172500017_13_496_SE_2017-09-07T04h48m22_dms_uncal_linear_dark_prep_object.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/linearized/Linearized_Dark_and_SBRefpix_NISNIRISSDARK-172500017_14_496_SE_2017-09-07T05h06m42_dms_uncal_linear_dark_prep_object.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/linearized/Linearized_Dark_and_SBRefpix_NISNIRISSDARK-172500017_15_496_SE_2017-09-07T05h28m22_dms_uncal_linear_dark_prep_object.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/linearized/Linearized_Dark_and_SBRefpix_NISNIRISSDARK-172500017_16_496_SE_2017-09-07T05h47m42_dms_uncal_linear_dark_prep_object.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/linearized/Linearized_Dark_and_SBRefpix_NISNIRISSDARK-172500017_17_496_SE_2017-09-07T06h09m02_dms_uncal_linear_dark_prep_object.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/linearized/Linearized_Dark_and_SBRefpix_NISNIRISSDARK-172500017_18_496_SE_2017-09-07T06h29m12_dms_uncal_linear_dark_prep_object.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/linearized/Linearized_Dark_and_SBRefpix_NISNIRISSDARK-172500017_19_496_SE_2017-09-07T06h49m52_dms_uncal_linear_dark_prep_object.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/linearized/Linearized_Dark_and_SBRefpix_NISNIRISSDARK-172500017_20_496_SE_2017-09-07T07h09m22_dms_uncal_linear_dark_prep_object.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/linearized/Linearized_Dark_and_SBRefpix_NISNIRISSDARK-172500017_21_496_SE_2017-09-07T07h29m52_dms_uncal_linear_dark_prep_object.fits.gz',
                               'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/niriss/darks/linearized/Linearized_Dark_and_SBRefpix_NISNIRISSDARK-172500017_22_496_SE_2017-09-07T07h50m32_dms_uncal_linear_dark_prep_object.fits.gz']

FGS_RAW_DARK_URLS = ['https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/fgs/darks/raw/29722_1x88_FGSF03512-D-NR-G2-5339214947_1_498_SE_2015-12-05T22h27m19_dms_uncal.fits.gz',
                     'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/fgs/darks/raw/29782_1x88_FGSF03872-PAR-5340074326_1_498_SE_2015-12-06T12h22m47_dms_uncal.fits.gz',
                     'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/fgs/darks/raw/29813_1x88_FGSF037221-MR-2-5340161743_1_498_SE_2015-12-06T16h45m10_dms_uncal.fits.gz',
                     'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/fgs/darks/raw/30632_1x88_FGSF03511-D-NR-G1-5346180117_1_497_SE_2015-12-12T19h00m12_dms_uncal.fits.gz',
                     'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/fgs/darks/raw/30670_1x88_FGSF03511-D-NR-G2-5346181816_1_498_SE_2015-12-12T21h31m01_dms_uncal.fits.gz',
                     'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/fgs/darks/raw/30742_1x88_FGSF03871-PAR-5347035139_1_498_SE_2015-12-13T05h23m30_dms_uncal.fits.gz',
                     'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/fgs/darks/raw/30749_1x88_FGSF03881-PAR-5347043800_1_497_SE_2015-12-13T09h02m01_dms_uncal.fits.gz',
                     'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/fgs/darks/raw/30829_1x88_FGSF037111-G1NRNC-5347151640_1_497_SE_2015-12-13T16h28m38_dms_uncal.fits.gz']
FGS_LINEARIZED_DARK_URLS = ['https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/fgs/darks/linearized/29722_1x88_FGSF03512-D-NR-G2-5339214947_1_498_SE_2015-12-05T22h27m19_dms_uncal_linearized.fits.gz',
                            'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/fgs/darks/linearized/29782_1x88_FGSF03872-PAR-5340074326_1_498_SE_2015-12-06T12h22m47_dms_uncal_linearized.fits.gz',
                            'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/fgs/darks/linearized/29813_1x88_FGSF037221-MR-2-5340161743_1_498_SE_2015-12-06T16h45m10_dms_uncal_linearized.fits.gz',
                            'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/fgs/darks/linearized/30632_1x88_FGSF03511-D-NR-G1-5346180117_1_497_SE_2015-12-12T19h00m12_dms_uncal_linearized.fits.gz',
                            'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/fgs/darks/linearized/30670_1x88_FGSF03511-D-NR-G2-5346181816_1_498_SE_2015-12-12T21h31m01_dms_uncal_linearized.fits.gz',
                            'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/fgs/darks/linearized/30742_1x88_FGSF03871-PAR-5347035139_1_498_SE_2015-12-13T05h23m30_dms_uncal_linearized.fits.gz',
                            'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/fgs/darks/linearized/30749_1x88_FGSF03881-PAR-5347043800_1_497_SE_2015-12-13T09h02m01_dms_uncal_linearized.fits.gz',
                            'https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/fgs/darks/linearized/30829_1x88_FGSF037111-G1NRNC-5347151640_1_497_SE_2015-12-13T16h28m38_dms_uncal_linearized.fits.gz']

TEMP_DISTORTION_REFERENCE_FILES = ['https://data.science.stsci.edu/redirect/JWST/jwst-simulations/mirage_reference_files/nircam/reference_files/nircam_distortion_files.tar.gz']

DISK_USAGE = {'nircam': {'crs': 1.1, 'psfs': 23, 'raw_darks': 79, 'lin_darks': 319, 'grism': 0.005},
              'niriss': {'crs': 0.26, 'psfs': 0.87, 'raw_darks': 31, 'lin_darks': 121, 'grism': 0.167},
              'fgs': {'crs': 0.31, 'psfs': .04, 'raw_darks': 11, 'lin_darks': 39}}

def download_file(url, file_name, output_directory='./'):
    """Download into the current working directory the
    file from Box given the direct URL

    Parameters
    ----------
    url : str
        URL to the file to be downloaded

    Returns
    -------
    download_filename : str
        Name of the downloaded file
    """
    download_filename = os.path.join(output_directory, file_name)

    # Only download the file if it doesn't already exist
    if not os.path.isfile(download_filename):
        print('Downloading: {}'.format(file_name))
        with requests.get(url, stream=True) as response:
            if response.status_code != 200:
                raise RuntimeError("Wrong URL - {}".format(url))
            with open(download_filename, 'wb') as f:
                for chunk in response.iter_content(chunk_size=2048):
                    if chunk:
                        f.write(chunk)
        print('Download of {} complete.'.format(file_name))
    else:
        print('{} already exists. Skipping download.'.format(file_name))
    return download_filename


def download_reffiles(directory, instrument='all', dark_type='linearized',
                      skip_darks=False, skip_cosmic_rays=False, skip_psfs=False,
                      skip_grism=False):
    """Download tarred and gzipped reference files. Expand, unzip and
    organize into the necessary directory structure such that Mirage
    can use them.

    Parameters
    ----------
    directory : str
        Directory into which the reference files are placed. This will
        be the directory set to the MIRAGE_DATA environment variable

    instrument : str
        Instruments for which to download data. Single string with
        comma-separated instrument names.

        If ``all``: download files for NIRCam, NIRISS, and FGS.

        If the name of an individual instrument (e.g. 'nircam'),
        download only the data for that instrument.

        If a list of instrument names, (e.g. 'nircam, niriss') download
        all data for those instruments.

    dark_type : str
        Type of dark current files to download. Options are:
        'linearized': download linearize dark current ramps
        'raw': download raw dark current ramps
        'both': download both raw and linearized dark current ramps

    skip_darks : bool
        If False (default), download the requested darks. If True,
        do not download the darks

    skip_comsic_rays : bool
        If False (default), download the requested cosmic ray libraries.
        If True, do not download the cosmic ray library.

    skip_psfs : bool
        If False (default), download the requested PSF libraries.
        If True, do not download the libraries.

    skip_grism : bool
        If False (default), download the grism-related reference files
        If True, do not download the grism data.
    """
    # Be sure the input instrument is a list
    file_list = get_file_list(instrument.lower(), dark_type.lower(),
                              skip_darks=skip_darks,
                              skip_cosmic_rays=skip_cosmic_rays,
                              skip_psfs=skip_psfs, skip_grism=skip_grism)

    # Download everything first
    for file_url in file_list:
        filename = os.path.split(file_url)[-1]
        local_file = os.path.join(directory, filename)
        download_file(file_url, filename, directory)

    # Now untar/organize. This way if the download is interrupted, it can
    # pick up where it left off, since no downloaded files will have been
    # moved yet.
    for file_url in file_list:
        filename = os.path.split(file_url)[-1]
        local_file = os.path.join(directory, filename)
        # Unzip and untar file
        if 'tar.gz' in local_file:
            print('Unzipping/extracting {}'.format(filename))
            file_object = tarfile.open(name=local_file, mode='r:gz')
            file_object.extractall(path=directory)
        else:
            # Darks need to be unzipped into the correct directory
            if 'linearized' in filename.lower():
                cal = 'linearized'
            else:
                cal = 'raw'

            # Determine directory
            if 'NRCNRC' in filename:
                det_str = filename.split('NRCNRC')[1].split('-')[0]
                if 'LONG' in det_str:
                    det_str.replace('LONG', '5')
                darks_dir = os.path.join(directory, 'mirage_data', 'nircam', 'darks')
                sub_directory = os.path.join(darks_dir, cal, det_str)
            elif 'NIRISS' in filename.lower():
                darks_dir = os.path.join(directory, 'mirage_data', 'niriss', 'darks')
                sub_directory = os.path.join(darks_dir, cal)
                #sub_directory = os.path.join(directory, 'mirage_data', 'niriss', 'darks', cal)
            elif 'FGS' in filename:
                darks_dir = os.path.join(directory, 'mirage_data', 'fgs', 'darks')
                sub_directory = os.path.join(darks_dir, cal)

            # Create the directory if it does not yet exist
            ensure_dir_exists(darks_dir)
            ensure_dir_exists(os.path.join(darks_dir, cal))
            ensure_dir_exists(sub_directory)

            final_location = os.path.join(sub_directory, filename)

            # Move the zipped file into the correct subdirectory
            if not os.path.isfile(final_location):
                print('Moving {} to {}'.format(filename, sub_directory))
                shutil.move(local_file, final_location)

            # Unzip
            unzipped_filename = final_location.replace('.gz', '')
            if not os.path.isfile(unzipped_filename):
                unzip_file(final_location)
            else:
                print('Unzipped file {} already exists. Skipping unzip.'.format(unzipped_filename))

    full_dir = os.path.abspath(directory)
    print(('Mirage reference files downloaded and extracted. \nBefore '
           'using Mirage, be sure to set the MIRAGE_DATA environment '
           'variable to point to:\n{}/mirage_data'.format(full_dir)))
    print('In bash: ')
    print('export MIRAGE_DATA="{}"'.format(os.path.join(full_dir, 'mirage_data')))


def get_file_list(instruments, dark_current, skip_darks=False, skip_cosmic_rays=False,
                  skip_psfs=False, skip_grism=False):
    """Collect the list of URLs corresponding to the Mirage reference
    files to be downloaded

    Parameters
    ----------
    instruments : list
        List of instrument names for which to download data

    dark_current : str
        Type of dark current files to download. Options are:
        'linearized': download linearize dark current ramps
        'raw': download raw dark current ramps
        'both': download both raw and linearized dark current ramps

    skip_darks : bool
        If False (default), include the requested darks. If True,
        do not include the darks

    skip_comsic_rays : bool
        If False (default), include the requested cosmic ray libraries.
        If True, do not include the cosmic ray library.

    skip_psfs : bool
        If False (default), include the requested PSF libraries.
        If True, do not include the libraries.

    skip_grism : bool
        If False (default), include the grism-related reference files.
        If True, do not include the grism files.

    Returns
    -------
    urls : list
        List of tuples. Each tuple contains:
        (URL for downloading, name of file)
    """
    urls = []
    total_download_size = 0.
    instrument_names = [name.strip().lower() for name in instruments.split(',')]

    if 'all' in instrument_names:
        instrument_names = ['nircam', 'niriss', 'fgs']

    for instrument_name in instrument_names:
        # NIRCam
        if instrument_name.lower() == 'nircam':

            if not skip_cosmic_rays:
                urls.extend(NIRCAM_CR_LIBRARY_URL)
                added_size = DISK_USAGE['nircam']['crs']
                total_download_size += added_size
                print('Size of NIRCam cosmic ray library file: {} Gb'.format(added_size))

            if not skip_psfs:
                urls.extend(NIRCAM_GRIDDED_PSF_URLS)
                added_size = DISK_USAGE['nircam']['psfs']
                total_download_size += added_size
                print('Size of NIRCam PSF library file: {} Gb'.format(added_size))

            if not skip_darks:
                if dark_current in ['linearized', 'both']:
                    urls.extend(NIRCAM_LINEARIZED_DARK_URLS)
                    added_size = DISK_USAGE['nircam']['lin_darks']
                    total_download_size += added_size
                    print('Size of NIRCam linearized dark files: {} Gb'.format(added_size))
                elif dark_current in ['raw', 'both']:
                    urls.extend(NIRCAM_RAW_DARK_URLS)
                    added_size = DISK_USAGE['nircam']['raw_darks']
                    total_download_size += added_size
                    print('Size of NIRCam raw dark files: {} Gb'.format(added_size))

            if not skip_grism:
                # Get grism-related files
                urls.extend(NIRCAM_GRISM_FILES)
                added_size = DISK_USAGE['nircam']['grism']
                total_download_size += added_size
                print('Size of NIRCam grism files: {} Gb'.format(added_size))

            # Get the temporary distortion reference files with the
            # correct coefficients
            urls.extend(TEMP_DISTORTION_REFERENCE_FILES)

        # NIRISS
        elif instrument_name.lower() == 'niriss':
            if not skip_cosmic_rays:
                urls.extend(NIRISS_CR_LIBRARY_URL)
                added_size = DISK_USAGE['niriss']['crs']
                total_download_size += added_size
                print('Size of NIRISS cosmic ray library file: {} Gb'.format(added_size))

            if not skip_psfs:
                urls.extend(NIRISS_GRIDDED_PSF_URLS)
                added_size = DISK_USAGE['niriss']['psfs']
                total_download_size += added_size
                print('Size of NIRISS PSF library file: {} Gb'.format(added_size))

            if not skip_darks:
                if dark_current in ['linearized', 'both']:
                    urls.extend(NIRISS_LINEARIZED_DARK_URLS)
                    added_size = DISK_USAGE['niriss']['lin_darks']
                    total_download_size += added_size
                    print('Size of NIRISS linearized dark files: {} Gb'.format(added_size))
                elif dark_current in ['raw', 'both']:
                    urls.extend(NIRISS_RAW_DARK_URLS)
                    added_size = DISK_USAGE['niriss']['raw_darks']
                    total_download_size += added_size
                    print('Size of NIRISS raw dark files: {} Gb'.format(added_size))

            if not skip_grism:
                # Get grism-related files
                urls.extend(NIRISS_GRISM_FILES)
                added_size = DISK_USAGE['niriss']['grism']
                total_download_size += added_size
                print('Size of NIRISS grism files: {} Gb'.format(added_size))

        # FGS
        elif instrument_name.lower() == 'fgs':
            if not skip_cosmic_rays:
                urls.extend(FGS_CR_LIBRARY_URL)
                added_size = DISK_USAGE['fgs']['crs']
                total_download_size += added_size
                print('Size of FGS cosmic ray library file: {} Gb'.format(added_size))

            if not skip_psfs:
                urls.extend(FGS_GRIDDED_PSF_URLS)
                added_size = DISK_USAGE['fgs']['psfs']
                total_download_size += added_size
                print('Size of FGS PSF library file: {} Gb'.format(added_size))

            if not skip_darks:
                if dark_current in ['linearized', 'both']:
                    urls.extend(FGS_LINEARIZED_DARK_URLS)
                    added_size = DISK_USAGE['fgs']['lin_darks']
                    total_download_size += added_size
                    print('Size of FGS linearized dark files: {} Gb'.format(added_size))
                elif dark_current in ['raw', 'both']:
                    urls.extend(FGS_RAW_DARK_URLS)
                    added_size = DISK_USAGE['fgs']['raw_darks']
                    total_download_size += added_size
                    print('Size of FGS raw dark files: {} Gb'.format(added_size))

    print("Total size of files to be downloaded: {} Gb".format(total_download_size))
    return urls


def unzip_file(filename):
    """Unzip a file

    Parameters
    ----------
    filename : str
        Name of file to unzip

    dir_name : str
        Directory into which the file is unzipped
    """
    if not os.path.isfile(filename):
        print('File {} does not exist.'.format(filename))
    print('Unzipping {}'.format(filename))
    with gzip.open(filename, 'rb') as file_in:
        unzipped_name = filename.replace('.gz', '')
        with open(unzipped_name, 'wb') as file_out:
            shutil.copyfileobj(file_in, file_out)


#if __name__ == '__main__':
#      params = parse_args()
#
#      download_reffiles(directory, instrument=instrument, dark_type=dark_type,
#                        skip_darks=skip_darks, skip_cosmic_rays=skip_cosmic_rays,
#                        skip_psfs=skip_psfs)
