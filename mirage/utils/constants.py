"""Various universal constants and values that are useful across the MIRAGE package.

Use
---
    ::
        from mirage.utils import constants
        inst_abb = constants.instrument_abbreviations
"""
import os
import pkg_resources

from astropy.time import Time
import astropy.units as u
import numpy as np

instrument_abbreviations = {'nircam': 'NRC',
                            'fgs': 'FGS',
                            'niriss': 'NIS',
                            'nirspec': 'NRS',
                            'miri': 'MIR'}

EXPTYPES = {"nircam": {"imaging": "NRC_IMAGE", "ts_imaging": "NRC_TSIMAGE",
                       "wfss": "NRC_WFSS", "ts_grism": "NRC_TSGRISM"},
            "niriss": {"imaging": "NIS_IMAGE", "ami": "NIS_AMI", "pom": "NIS_IMAGE",
                       "wfss": "NIS_WFSS", "soss": "NIS_SOSS"},
            "fgs": {"imaging": "FGS_IMAGE"}}

# Number of detector resets prior to the start of an exposure
NUM_RESETS_BEFORE_EXP = {"nircam": {"full": 0, "sub": 1},
                         "niriss": {"full": 0, "sub": 1},
                         "fgs": {"full": 0, "sub": 0}
                         }

# Numer of detector resets between integrations in an exposure
NUM_RESETS_BEFORE_INT = {"nircam": 1, "niriss": 1, "fgs":1}

# SEARCH STRINGS TO USE FOR FGS DARKS (needed because darks for the two
# detectors are mixed in a single directory)
FGS1_DARK_SEARCH_STRING = '*_497_*fits'
FGS2_DARK_SEARCH_STRING = '*_498_*fits'

# Supported NIRISS filters
NIRISS_FILTER_WHEEL_FILTERS = ['F277W', 'F356W', 'F380M', 'F430M', 'F444W', 'F480M']
NIRISS_PUPIL_WHEEL_FILTERS = ['F090W', 'F115W', 'F158M', 'F140M', 'F150W', 'F200W']
NIRISS_FILTERS = NIRISS_FILTER_WHEEL_FILTERS + NIRISS_PUPIL_WHEEL_FILTERS #+ \
                 #['{}/CLEARP'.format(element) for element in NIRISS_FILTER_WHEEL_FILTERS] + \
                 #['{}/CLEAR'.format(element) for element in NIRISS_FILTER_WHEEL_FILTERS] + \
                 #['{}/CLEAR'.format(element) for element in NIRISS_PUPIL_WHEEL_FILTERS] + \
                 #['CLEARP/{}'.format(element) for element in NIRISS_FILTER_WHEEL_FILTERS] + \
                 #['CLEAR/{}'.format(element) for element in NIRISS_FILTER_WHEEL_FILTERS] + \
                 #['CLEAR/{}'.format(element) for element in NIRISS_PUPIL_WHEEL_FILTERS]


# Supported NIRCam filters
NIRCAM_PUPIL_WHEEL_FILTERS = ['F162M', 'F164N', 'F323N', 'F405N', 'F466N', 'F470N']
NIRCAM_GO_PW_FILTER_PAIRINGS = {'F162M': 'F150W2', 'F164N': 'F150W2', 'F323N': 'F322W2',
                                'F405N': 'F444W', 'F466N': 'F444W', 'F470N': 'F444W'}
NIRCAM_WL8_CROSSING_FILTERS = ['F070W', 'F090W', 'F115W', 'F140M', 'F150W2', 'F150W',
                               'F162M', 'F164N', 'F182M', 'F187N', 'F200W', 'F210M',
                               'F212N', 'WLP4']
NIRCAM_CLEAR_CROSSING_FILTERS = ['F070W', 'F090W', 'F115W', 'F140M', 'F150W', 'F150W2',
                                 'F182M', 'F187N', 'F200W', 'F210M', 'F212N', 'WLP4',
                                 'F250M', 'F277W', 'F300M', 'F322W2', 'F335M', 'F356W',
                                 'F360M', 'F410M', 'F430M', 'F444W', 'F460M', 'F480M']
NIRCAM_2_FILTER_CROSSES = ['WLP4/WLP8', 'WLP4/WLM8', 'F150W2/F162M', 'F150W2/F164N',
                           'F322W2/F323N', 'F444W/F405N', 'F444W/F466N', 'F444W/F470N']
NIRCAM_UNSUPPORTED_PUPIL_VALUES = ['GDHS0', 'GDHS60', 'FLAT', 'MASKBAR', 'MASKIPR', 'MASKRND', 'PINHOLES']

NIRCAM_FILTERS = NIRCAM_CLEAR_CROSSING_FILTERS + NIRCAM_2_FILTER_CROSSES + NIRCAM_PUPIL_WHEEL_FILTERS + \
                ['{}/CLEAR'.format(element) for element in NIRCAM_CLEAR_CROSSING_FILTERS] + \
                ['CLEAR/{}'.format(element) for element in NIRCAM_CLEAR_CROSSING_FILTERS] + \
                ['{}/WLP8'.format(element) for element in NIRCAM_WL8_CROSSING_FILTERS] + \
                ['WLP8/{}'.format(element) for element in NIRCAM_WL8_CROSSING_FILTERS] + \
                ['{}/WLM8'.format(element) for element in NIRCAM_WL8_CROSSING_FILTERS] + \
                ['WLM8/{}'.format(element) for element in NIRCAM_WL8_CROSSING_FILTERS]

NIRCAM_SW_GRISMTS_APERTURES = ['NRCA1_GRISMTS256', 'NRCA1_GRISMTS128', 'NRCA1_GRISMTS64', 'NRCA1_GRISMTS',
                               'NRCA3_GRISMTS256', 'NRCA3_GRISMTS128', 'NRCA3_GRISMTS64', 'NRCA3_GRISMTS']

NIRCAM_LW_GRISMTS_APERTURES = ['NRCA5_GRISM256_F277W', 'NRCA5_GRISM128_F277W', 'NRCA5_GRISM64_F277W', 'NRCA5_GRISM_F277W',
                               'NRCA5_GRISM256_F322W2', 'NRCA5_GRISM128_F322W2', 'NRCA5_GRISM64_F322W2', 'NRCA5_GRISM_F322W2',
                               'NRCA5_GRISM256_F356W', 'NRCA5_GRISM128_F356W', 'NRCA5_GRISM64_F356W', 'NRCA5_GRISM_F356W',
                               'NRCA5_GRISM256_F444W', 'NRCA5_GRISM128_F444W', 'NRCA5_GRISM64_F444W', 'NRCA5_GRISM_F444W']

# For consistency with NIRCam/NIRISS filters
FGS_FILTERS = ["GUIDER1", "GUIDER2"]

FLAMBDA_CGS_UNITS = u.erg / u.second / u.cm / u.cm / u.AA
FLAMBDA_MKS_UNITS = u.watt / u.meter**2 / u.micron
FNU_CGS_UNITS = u.erg / u.second / u.cm / u.cm / u.Hz
FNU_MKS_UNITS = u.watt / u.m**2 / u.Hz

CATALOG_YAML_ENTRIES = ['pointsource', 'galaxyListFile', 'extended', 'movingTargetList',
                        'movingTargetSersic', 'movingTargetExtended', 'movingTargetToTrack',
                        'tso_imaging_catalog', 'tso_grism_catalog']

IMAGING_ALLOWED_CATALOGS = ['pointsource', 'galaxyListFile', 'extended', 'movingTargetList',
                            'movingTargetSersic', 'movingTargetExtended', 'movingTargetToTrack']

WFSS_ALLOWED_CATALOGS = ['pointsource', 'galaxyListFile', 'extended']

TS_IMAGING_ALLOWED_CATALOGS = ['pointsource', 'galaxyListFile', 'extended', 'tso_imaging_catalog']
TS_GRISM_ALLOWED_CATALOGS = ['pointsource', 'galaxyListFile', 'extended', 'tso_grism_catalog']

TSO_MODES = ['ts_imaging', 'ts_grism']

# Upper limit to the size of a seed image or dark current array. Arrays
# containing more pixels than this limit will be split into segment files.
FILE_SPLITTING_LIMIT = 160. * 2048. * 2048

CRDS_FILE_TYPES = {'badpixmask': 'mask',
                   'astrometric': 'distortion',
                   'gain': 'gain',
                   'ipc': 'ipc',
                   'linearity': 'linearity',
                   'pixelAreaMap': 'area',
                   'transmission': 'transmission',
                   'saturation': 'saturation',
                   'superbias': 'superbias',
                   'pixelflat': 'flat',
                   'photom': 'photom'}

MEAN_GAIN_VALUES = {'nircam': {'swa': 2.05, 'swb': 2.05, 'lwa': 1.82, 'lwb': 1.82,
                               'nrca1': 2.05, 'nrca2': 2.05, 'nrca3': 2.05, 'nrca4': 2.05,
                               'nrcb1': 2.05, 'nrcb2': 2.05, 'nrcb3': 2.05, 'nrcb4': 2.05,
                               'nrca5': 1.82, 'nrcb5': 1.82},
                    'niriss': 1.611902,
                    'fgs': {'guider1': 1.9712874, 'guider2': 1.6897644}
                    }

# Factor by which throughput is decreased when the Grism is used compared
# to when it is not in the beam
NIRISS_GRISM_THROUGHPUT_FACTOR = 0.8

# Factors by which total throughput is reduced and gridded PSF library
# signal (output from WebbPSF) is reduced. Used in normalization_check()
# in psf_selection.py when performing a sanity check on the input PSF signal
NIRISS_NRM_PSF_THROUGHPUT_REDUCTION = 0.17
NIRISS_CLEARP_PSF_THROUGHPUT_REDUCTION = 0.84
NIRCAM_WLP12_PSF_THROUGHPUT_REDUCTION = 0.5

# PSF normalization check - default min and max allowed values for the
# total signal in the PSF divided by the square of the oversampling factor.
# In most cases this value will be close to 1.0. In cases where NIRISS
# NRM or CLEARP are used, these values will be reduced by the factors above.
PSF_NORM_MAX = 1.0
PSF_NORM_MIN = 0.8

# Minimum signal rate for a pixel to be included in the segmentation map.
SEGMENTATION_MIN_SIGNAL_RATE = 0.031  # ADU/sec
SUPPORTED_SEGMENTATION_THRESHOLD_UNITS = ['adu/s', 'adu/sec', 'e/s', 'e/sec', 'mjy/str', 'mjy/sr', 'erg/cm2/a', 'erg/cm2/hz']

# For use in converting background MJy/sr to e-/sec
PRIMARY_MIRROR_AREA = 25.326 * u.meter * u.meter
PLANCK = 6.62607004e-34  * u.meter * u.meter * u.kg / u.second

# Fraction of the total Sersic
SERSIC_FRACTIONAL_SIGNAL = 0.9995

# Configuration file for defining the logs
LOG_CONFIG_FILENAME = 'logging_config.yaml'

# Standard log file names. These are needed because we need to have a
# log file name before we know the yaml/xml filename
STANDARD_LOGFILE_NAME = 'mirage_latest.log'

# Default stamp image to use for ghost sources resulting from point sources
MODULE_PATH = pkg_resources.resource_filename('mirage', '')
CONFIG_DIR = os.path.join(MODULE_PATH, 'config')
NIRISS_GHOST_GAP_FILE = os.path.join(CONFIG_DIR, 'niriss_ghost_gap_summary.txt')
NIRISS_GHOST_GAP_URL = 'https://raw.githubusercontent.com/spacetelescope/niriss_ghost/main/niriss_ghost/niriss_ghost_gap_summary.txt'
DEFAULT_NIRISS_PTSRC_GHOST_FILE = os.path.join(os.path.expandvars('$MIRAGE_DATA'), 'niriss/ghosts/', 'niriss_ghost_cen.fits')

# Vega spectrum to use when translating/normalizing input spectra to a given vegamag
# (in spectra_from_catalog.rescale_normalized_spectra). This needs to be consistent with the
# vega spectrum used to generate zeropoints in the zeropoints files in the config dirctory!
# For the June 2021 update to the zeropoints (based on updated NIRCam gain values), we used
# alpha_lyr_stis_010.fits from synphot.
VEGA_SPECTRUM = 'http://ssb.stsci.edu/cdbs/calspec/alpha_lyr_stis_010.fits'

# The base time to use when populating the group times and moving target locations tables.
# Times are reported relative to this basetime
TABLE_BASETIME = Time('2000-01-01T00:00:00')