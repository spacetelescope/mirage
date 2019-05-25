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

NIRISS_FILTER_WHEEL_FILTERS = ['F277W', 'F356W', 'F380M', 'F430M', 'F444W', 'F480M']
NIRISS_PUPIL_WHEEL_FILTERS = ['F090W', 'F115W', 'F158M', 'F140M', 'F150W', 'F200W']

FLAMBDA_CGS_UNITS = u.erg / u.second / u.cm / u.cm / u.AA
FLAMBDA_MKS_UNITS = u.watt / u.meter**2 / u.micron
FNU_CGS_UNITS = u.erg / u.second / u.cm / u.cm / u.Hz
FNU_MKS_UNITS = u.watt / u.m**2 / u.Hz

CATALOG_YAML_ENTRIES = ['pointsource', 'galaxyListFile', 'extended', 'movingTargetList',
                        'movingTargetSersic', 'movingTargetExtended', 'movingTargetToTrack']

TSO_MODES = ['ts_imaging', 'ts_wfss']

# Upper limit to the size of a seed image or dark current array. Arrays
# containing more pixels than this limit will be split into segment files.
FILE_SPLITTING_LIMIT = 38. * 2048. * 2048


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
