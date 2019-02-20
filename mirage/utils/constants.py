"""Various universal constants and values that are useful across the MIRAGE package.

Use
---
    ::
        from mirage.utils import constants
        inst_abb = constants.instrument_abbreviations
"""
import astropy.units as u

instrument_abbreviations = {'nircam': 'NRC',
                            'fgs': 'FGS',
                            'niriss': 'NIS',
                            'nirspec': 'NRS',
                            'miri': 'MIR'}

FLAMBDA_UNITS = u.erg / u.second / u.cm / u.cm / u.AA
FNU_UNITS = u.erg / u.second / u.cm / u.cm / u.Hz