"""Various universal constants and values that are useful across the MIRAGE package.

Use
---
    ::
        from mirage.utils import constants
        inst_abb = constants.instrument_abbreviations
"""

instrument_abbreviations = {'nircam': 'NRC',
                            'fgs': 'FGS',
                            'niriss': 'NIS',
                            'nirspec': 'NRS',
                            'miri': 'MIR'}

NIRISS_FILTER_WHEEL_FILTERS = ['F277W', 'F356W', 'F380M', 'F430M', 'F444W', 'F480M']
NIRISS_PUPIL_WHEEL_FILTERS = ['F090W', 'F115W', 'F158M', 'F140M', 'F150W', 'F200W']