#! /usr/bin/env python

"""
This module contains functions used to indentify and download reference files
from CRDS and place them in the expected location for Mirage reference files.

This module uses the crds software package (LINK HERE) which is installed when
the JWST calibration pipeline package is installed. Reference files are
identified by supplying some basic metadata from the exposure being calibrated.

See https://hst-crds.stsci.edu/static/users_guide/library_use.html#crds-getreferences
for a description of the function used for this task.

Author
------

    - Bryan Hilbert


Use
---

    This module can be used as such:
    ::
        from mirage.reference_files import crds
        params = {'INSTRUME': 'NIRCAM', 'DETECTOR': 'NRCA1'}
        reffiles = crds.get_reffiles(params)
"""

import crds
import datetime

from mirage.utils.constants import EXPTYPES

def dict_from_yaml(yaml_dict):
    """Create a dictionary to be used as input to the CRDS getreferences
    function from the nested dictionary created when a standard Mirage
    input yaml file is read in.

    Parameters
    ----------
    yaml_dict : dict
        Nested dictionary from reading in yaml file

    Returns
    -------
    crds_dict : dict
        Dictionary of information necessary to select refernce files
        via getreferences().
    """
    crds_dict = {}
    instrument = yaml_dict['Inst']['instrument'].upper()
    crds_dict['INSTRUME'] = instrument
    crds_dict['READPATT'] = yaml_dict['Readout']['readpatt'].upper()

    # Currently, all reference files that use SUBARRAY as a selection
    # criteria contain SUBARRAY = 'GENERIC', meaning that SUBARRAY
    # actually isn't important. So let's just set it to FULL here.
    crds_dict['SUBARRAY'] = 'FULL'

    # Use the current date and time in order to get the most recent
    # reference file
    crds_dict['DATE-OBS'] = datetime.date.today().isoformat()
    current_date = datetime.datetime.now()
    crds_dict['TIME-OBS'] = current_date.time().isoformat()

    array_name = yaml_dict['Readout']['array_name']
    crds_dict['DETECTOR'] = array_name.split('_')[0].upper()

    crds_dict['EXP_TYPE'] = EXPTYPES[instrument.lower()][yaml_dict['Inst']['mode'].lower()]

    # This assumes that filter and pupil names match up with reality,
    # as opposed to the more user-friendly scheme of allowing any
    # filter to be in the filter field.
    crds_dict['FILTER'] = yaml_dict['Readout']['filter']
    crds_dict['PUPIL'] = yaml_dict['Readout']['pupil']
    return crds_dict


def get_reffiles(yaml_input):
    """
    """
    dict_for_crds = dict_from_yaml(yaml_input)
    reffile_mapping = crds.getreferences(dict_for_crds)
    return reffile_mapping




