"""Various universal constants and values that are useful across the MIRAGE package.

Use
---
    ::
        from mirage.utils import constants
        inst_abb = constants.instrument_abbreviations
"""
import astropy.units as u
import numpy as np

instrument_abbreviations = {'nircam': 'NRC',
                            'fgs': 'FGS',
                            'niriss': 'NIS',
                            'nirspec': 'NRS',
                            'miri': 'MIR'}

EXPTYPES = {"nircam": {"imaging": "NRC_IMAGE", "ts_imaging": "NRC_TSIMAGE",
                       "wfss": "NRC_WFSS", "ts_grism": "NRC_TSGRISM"},
            "niriss": {"imaging": "NIS_IMAGE", "ami": "NIS_IMAGE", "pom": "NIS_IMAGE",
                       "wfss": "NIS_WFSS"},
            "fgs": {"imaging": "FGS_IMAGE"}}

NIRISS_FILTER_WHEEL_FILTERS = ['F277W', 'F356W', 'F380M', 'F430M', 'F444W', 'F480M']
NIRISS_PUPIL_WHEEL_FILTERS = ['F090W', 'F115W', 'F158M', 'F140M', 'F150W', 'F200W']

FLAMBDA_CGS_UNITS = u.erg / u.second / u.cm / u.cm / u.AA
FLAMBDA_MKS_UNITS = u.watt / u.meter**2 / u.micron
FNU_CGS_UNITS = u.erg / u.second / u.cm / u.cm / u.Hz
FNU_MKS_UNITS = u.watt / u.m**2 / u.Hz

CATALOG_YAML_ENTRIES = ['pointsource', 'galaxyListFile', 'extended', 'movingTargetList',
                        'movingTargetSersic', 'movingTargetExtended', 'movingTargetToTrack',
                        'tso_imaging_catalog', 'tso_grism_catalog']

TSO_MODES = ['ts_imaging', 'ts_grism']

# Upper limit to the size of a seed image or dark current array. Arrays
# containing more pixels than this limit will be split into segment files.
FILE_SPLITTING_LIMIT = 38. * 2048. * 2048

CRDS_FILE_TYPES = {'badpixmask': 'mask',
                   'astrometric': 'distortion',
                   'gain': 'gain',
                   'ipc': 'ipc',
                   'linearity': 'linearity',
                   'pixelAreaMap': 'area',
                   'saturation': 'saturation',
                   'superbias': 'superbias',
                   'pixelflat': 'flat'}

# Apertures that use fewer than all detectors in the module
LIMITED_DETECTOR_APERTURES = {'SUBGRISM64': ['NRCA1', 'NRCA3', 'NRCA5', 'NRCB4', 'NRCB2', 'NRCB5'],
                              'SUBGRISM128': ['NRCA1', 'NRCA3', 'NRCA5', 'NRCB4', 'NRCB2', 'NRCB5'],
                              'SUBGRISM256': ['NRCA1', 'NRCA3', 'NRCA5', 'NRCB4', 'NRCB2', 'NRCB5'],
                              'SUB32TATSGRISM': ['NRCA5', 'NRCB5'],
                              'DHSPIL': ['NRCA3', 'NRCB4'],
                              'DHSPIL_SUB96': ['NRCA3', 'NRCB4'],
                              'DHSPIL_WEDGES': ['NRCA3', 'NRCB4'],
                              'FP1': ['NRCA3', 'NRCB4'],
                              'FP1_SUB8': ['NRCA3', 'NRCB4'],
                              'FP1_SUB64': ['NRCA3', 'NRCB4'],
                              'FP2MIMF': ['NRCA3', 'NRCB4'],
                              'FP3MIMF': ['NRCA1', 'NRCB2'],
                              'FP4MIMF': ['NRCA2', 'NRCB1'],
                              'FP5MIMF': ['NRCA4', 'NRCB3'],
                              'SUB64P': ['NRCA3', 'NRCA5', 'NRCB1', 'NRCB5'],
                              'SUB160P': ['NRCA3', 'NRCA5', 'NRCB1', 'NRCB5'],
                              'SUB400P': ['NRCA3', 'NRCA5', 'NRCB1', 'NRCB5'],
                              }


def grism_factor(instrument_name):
    """Return the factor by which the field of view is expanded when
    creating grism simulations compared to direct image simulations

    Parameters
    ----------
    instrument_name : str
        JWST instrument name

    Returns
    -------
    factor : float
        Multiplicative factor by which the fov is enlarged
    """
    if instrument_name.lower() == 'nircam':
        return np.sqrt(2.)
    elif instrument_name.lower() == 'niriss':
        return 2322./2048.
