"""Various utility functions.

Authors
-------
    - Lauren Chambers
    - Bryan Hilbert

Use
---
    This module can be imported as such:

    >>> import utils
    siaf_files = get_siaf()

Dependencies
------------
    The user must have a configuration file named ``siaf_config.json``
    placed in mirage/config/ directory.
"""

import copy
import json
import os
import pkg_resources
import re
import yaml

from astropy.io import ascii as asc
import numpy as np
from scipy.stats import sigmaclip

from mirage.utils.constants import CRDS_FILE_TYPES, NIRISS_FILTER_WHEEL_FILTERS, NIRISS_PUPIL_WHEEL_FILTERS


def append_dictionary(base_dictionary, added_dictionary, braid=False):
    """Append the content of added_dictionary key-by-key to the base_dictionary.

    This assumes that the keys refer to lists.

    Parameters
    ----------
    base_dictionary : dict
    added_dictionary : dict
    braid : bool
        If true, the elements of added_dictionary are added in alternating sequence.
        This is used to synchronize parallel observations with the pointing file.

    Returns
    -------
    new_dictionary : dict
        Dictionary where every key holds a list of lists

    """
    new_dictionary = copy.deepcopy(base_dictionary)

    # extract an arbitrary key name
    first_key = [key for i, key in enumerate(base_dictionary.keys()) if i == 0][0]

    # Insert keys from added_dictionary that are not yet present in base_dictionary
    for key in added_dictionary.keys():
        if key not in base_dictionary.keys():
            new_dictionary[key] = ['None'] * len(base_dictionary[first_key])

    # Append the items
    for key in new_dictionary.keys():
        if key not in added_dictionary.keys():
            continue
        # print('{} {}'.format(key, new_dictionary[key]))
        if len(new_dictionary[key]) == 0:
            new_dictionary[key] = added_dictionary[key]
        else:
            if braid:
                # solution from https://stackoverflow.com/questions/3678869/pythonic-way-to-combine-two-lists-in-an-alternating-fashion
                new_dictionary[key] = [sub[i] for i in range(len(added_dictionary[key])) for sub in
                                       [new_dictionary[key], added_dictionary[key]]]
            else:
                new_dictionary[key] = new_dictionary[key] + added_dictionary[key]

    return new_dictionary


def calc_frame_time(instrument, aperture, xdim, ydim, amps):
    """Calculate the readout time for a single frame
    of a given size and number of amplifiers. Note that for
    NIRISS and FGS, the fast readout direction is opposite to
    that in NIRCam, so we switch xdim and ydim so that we can
    keep a single equation.

    Parameters:
    -----------
    instrument : str
        Name of the instrument being simulated

    aperture : str
        Name of aperture being simulated (e.g "NRCA1_FULL")
        Currently this is only used to check for the FGS
        ACQ1 aperture, which uses a unique value of colpad
        below.

    xdim : int
        Number of columns in the frame

    ydim : int
        Number of rows in the frame

    amps : int
        Number of amplifiers used to read out the frame

    Returns:
    --------
    frametime : float
        Readout time in seconds for the frame
    """
    instrument = instrument.lower()
    if instrument == "nircam":
        colpad = 12

        # Fullframe
        if amps == 4:
            rowpad = 1
            fullpad = 1
        else:
            # All subarrays
            rowpad = 2
            fullpad = 0

            if ((xdim <= 8) & (ydim <= 8)):
                # The smallest subarray
                rowpad = 3

    elif instrument == "niriss":
        # Reverse x and y since NIRISS's fast readout direction is
        # opposite of NIRCam's
        tmpx = copy.deepcopy(xdim)
        xdim = copy.deepcopy(ydim)
        ydim = tmpx

        colpad = 12

        # Fullframe
        if amps == 4:
            rowpad = 1
            fullpad = 1
        else:
            rowpad = 2
            fullpad = 0

    elif instrument == 'fgs':
        # Reverse x and y since FGS's fast readout direction is
        # opposite of NIRCam's
        tmpx = copy.deepcopy(xdim)
        xdim = copy.deepcopy(ydim)
        ydim = tmpx

        colpad = 12
        fullpad = 0

        if ((xdim == 2048) & (ydim == 2048)):
            rowpad = 1
            fullpad = 1
        else:
            rowpad = 2

        if ((xdim <= 32) & (ydim <= 32)):
            colpad = 6
            rowpad = 1

    return ((1.0 * xdim / amps + colpad) * (ydim + rowpad) + fullpad) * 1.e-5


def check_niriss_filter(oldfilter, oldpupil):
    """
    This is a utility function that checks the FILTER and PUPIL parameters read in from the .yaml file and makes sure
    that for NIRISS the filter and pupil names are correct, releaving the user of the need to remember which of the 12
    filters are in the filter wheel and which are in the pupil wheel.  If the user puts the correct values for filter and
    pupil nothing is done, but the user can just put the filter name in the filter parameter and CLEAR in the pupil
    parameter, in which case this routine sorts out which value actually goes where.  Hence for example one can have a
    parameter file with FILTER = F277W and PUPIL = CLEAR or FILTER = F090W and PUPIL = CLEAR for imaging and either will
    produce the desired result although the actual pairings would be F277W/CLEARP and CLEAR/F090W respectively.

    The code also corrects CLEARP to CLEAR if needed.

    Parameters
    ----------

    oldfilter :  str
        Assumed to be the self.params['Readout']['filter'] value from the .yaml file

    oldpupil :  str
        Assumed to be the self.params['Readout']['pupil'] value from the .yaml file

    Returns
    -------

    newfilter : str
        The proper FILTER parameter for the requested NIRISS filter name

    newpupil : str
        The proper PUPIL parameter for the requested NIRISS filter name

    Note that this routine should only be called when the instrument is NIRISS.
    """
    str1 = oldfilter
    str2 = oldpupil

    if oldfilter in NIRISS_PUPIL_WHEEL_FILTERS:
        newfilter = str2
        newpupil = str1
        if newfilter == 'CLEARP':
            newfilter = 'CLEAR'
    elif oldfilter in NIRISS_FILTER_WHEEL_FILTERS:
        if oldpupil == 'CLEAR':
            newpupil = 'CLEARP'
        else:
            newpupil = oldpupil
        newfilter = oldfilter
    else:
        # Case where oldfilter is an optic other than a filter (e.g. a grism)
        newfilter = oldfilter
        newpupil = oldpupil
    return newfilter, newpupil


def crop_to_subarray(data, bounds):
    """
    Crop the given full frame array down to the appropriate
    subarray size and location based on the requested subarray
    name.

    Parameters
    ----------
    data : numpy.ndarray
        Full frame image or ramp. (x,y) = (2048, 2048)
        May be 2D, 3D, or 4D

    bounds : list
        4-element list containing the full frame indices that
        define the position of the subarray.
        [xstart, ystart, xend, yend]

    Returns
    -------
    data : numpy.ndarray
        Input array cropped in x and y dimensions
    """
    dimensions = len(data.shape)
    yl, xl = data.shape[-2:]

    valid = [False, False, False, False]
    valid = [(b >= 0 and b < xl) for b in bounds[0:3:2]]
    validy = [(b >= 0 and b < yl) for b in bounds[1:4:2]]
    valid.extend(validy)

    if all(valid):
        if dimensions == 2:
            return data[bounds[1]:bounds[3] + 1, bounds[0]:bounds[2] + 1]
        elif dimensions == 3:
            return data[:, bounds[1]:bounds[3] + 1, bounds[0]:bounds[2] + 1]
        elif dimensions == 4:
            return data[:, :, bounds[1]:bounds[3] + 1, bounds[0]:bounds[2] + 1]
        else:
            raise ValueError(("In crop_to_subarray, input array is not 2, 3, or 4D."))
    else:
            raise ValueError(("WARNING: subarray bounds are outside the "
                              "dimensions of the input array."))


def ensure_dir_exists(fullpath):
    """Creates dirs from ``fullpath`` if they do not already exist.
    """
    if not os.path.exists(fullpath):
        os.makedirs(fullpath)


def expand_environment_variable(variable_name, offline=False):
    """ Expand an environment variable and check that the directory exists

    Parameters
    ----------
    variable_name : str
        Environment variable name

    Returns
    -------
    variable_directory : str
        Path pointed to by the environment variable
    """
    variable_directory = os.environ.get(variable_name)
    if variable_directory is None:
        raise ValueError(("{} environment variable is not set. "
                          "This must be set to the base directory "
                          "containing the darks, cosmic ray, PSF, etc "
                          "input files needed for the simulation. "
                          "These files must be downloaded separately "
                          "from the Mirage package.".format(variable_name)))
    if not offline:
        if not os.path.isdir(variable_directory):
            raise FileNotFoundError(("The directory contained in the {} "
                                     "environment variable: {} does not exist or "
                                     "is not accessible.".format(variable_name, variable_directory)))
    return variable_directory


def full_paths(params, module_path, crds_dictionary, offline=False):
    """Expand the relevant input paths from the input yaml file to be full
    paths.

    Parameters
    ----------
    params : dict
        Nested dictionary read in from an input yaml file

    module_path : str
        Path to the Mirage module

    crds_dictionary : dict
        Dictionary of basic observation metadata that is used by CRDS to
        select appropriate reference files.

    offline : bool
        Set to True in order to skip downloading reference files from CRDS
        and instead only return the dictionary of reference file names.

    Returns
    -------
    params : dict
        Modified nested dictionary including full paths
    """
    # Place import statement here in order to avoid circular import
    # with crd_tools
    from mirage.reference_files import crds_tools
    pathdict = {'Reffiles': ['dark', 'linearized_darkfile', 'badpixmask',
                             'superbias', 'linearity', 'saturation', 'gain',
                             'pixelflat', 'illumflat', 'ipc', 'astrometric',
                             'crosstalk', 'occult', 'pixelAreaMap',
                             'subarray_defs', 'filtpupilcombo',
                             'flux_cal', 'readpattdefs', 'filter_throughput',
                             'filter_wheel_positions'],
                'simSignals': ['pointsource', 'psfpath', 'galaxyListFile',
                               'extended', 'movingTargetList', 'movingTargetSersic',
                               'movingTargetExtended', 'movingTargetToTrack',
                               'psf_wing_threshold_file', 'zodiacal', 'scattered'],
                'cosmicRay': ['path'],
                'newRamp': ['dq_configfile', 'sat_configfile', 'superbias_configfile',
                            'refpix_configfile', 'linear_configfile'],
                'Output': ['file', 'directory']}

    all_config_files = {'nircam': {'Reffiles-subarray_defs': 'NIRCam_subarray_definitions.list',
                                   'Reffiles-flux_cal': 'NIRCam_zeropoints.list',
                                   'Reffiles-crosstalk': 'xtalk20150303g0.errorcut.txt',
                                   'Reffiles-readpattdefs': 'nircam_read_pattern_definitions.list',
                                   'Reffiles-filter_throughput': 'placeholder.txt',
                                   'Reffiles-filtpupilcombo': 'nircam_filter_pupil_pairings.list',
                                   'Reffiles-filter_wheel_positions': 'nircam_filter_and_pupil_wheel_positions.txt',
                                   'simSignals-psf_wing_threshold_file': 'nircam_psf_wing_rate_thresholds.txt',
                                   'newRamp-dq_configfile': 'dq_init.cfg',
                                   'newRamp-sat_configfile': 'saturation.cfg',
                                   'newRamp-superbias_configfile': 'superbias.cfg',
                                   'newRamp-refpix_configfile': 'refpix.cfg',
                                   'newRamp-linear_configfile': 'linearity.cfg'},
                        'niriss': {'Reffiles-subarray_defs': 'niriss_subarrays.list',
                                   'Reffiles-flux_cal': 'niriss_zeropoints.list',
                                   'Reffiles-crosstalk': 'niriss_xtalk_zeros.txt',
                                   'Reffiles-readpattdefs': 'niriss_readout_pattern.txt',
                                   'Reffiles-filter_throughput': 'placeholder.txt',
                                   'Reffiles-filtpupilcombo': 'niriss_dual_wheel_list.txt',
                                   'Reffiles-filter_wheel_positions': 'niriss_filter_and_pupil_wheel_positions.txt',
                                   'simSignals-psf_wing_threshold_file': 'niriss_psf_wing_rate_thresholds.txt',
                                   'newRamp-dq_configfile': 'dq_init.cfg',
                                   'newRamp-sat_configfile': 'saturation.cfg',
                                   'newRamp-superbias_configfile': 'superbias.cfg',
                                   'newRamp-refpix_configfile': 'refpix.cfg',
                                   'newRamp-linear_configfile': 'linearity.cfg'},
                        'fgs': {'Reffiles-subarray_defs': 'guider_subarrays.list',
                                'Reffiles-flux_cal': 'guider_zeropoints.list',
                                'Reffiles-crosstalk': 'guider_xtalk_zeros.txt',
                                'Reffiles-readpattdefs': 'guider_readout_pattern.txt',
                                'Reffiles-filter_throughput': 'placeholder.txt',
                                'Reffiles-filtpupilcombo': 'guider_filter_dummy.list',
                                'Reffiles-filter_wheel_positions': 'dummy.txt',
                                'simSignals-psf_wing_threshold_file': 'fgs_psf_wing_rate_thresholds.txt',
                                'newRamp-dq_configfile': 'dq_init.cfg',
                                'newRamp-sat_configfile': 'saturation.cfg',
                                'newRamp-superbias_configfile': 'superbias.cfg',
                                'newRamp-refpix_configfile': 'refpix.cfg',
                                'newRamp-linear_configfile': 'linearity.cfg'}}
    config_files = all_config_files[params['Inst']['instrument'].lower()]

    for key1 in pathdict:
        for key2 in pathdict[key1]:
            if params[key1][key2].lower() not in ['none', 'config', 'crds']:
                params[key1][key2] = os.path.abspath(os.path.expandvars(params[key1][key2]))
            elif params[key1][key2].lower() == 'config':
                cfile = config_files['{}-{}'.format(key1, key2)]
                fpath = os.path.join(module_path, 'config', cfile)
                params[key1][key2] = fpath
                print("'config' specified: Using {} for {}:{} input file".format(fpath, key1, key2))
            elif key2 in CRDS_FILE_TYPES.keys() and params[key1][key2].lower() == 'crds':
                # 'crds' set for one of the reference files means to
                # search for the best reference file currently in CRDS
                # and download if it is not already present in
                # self.crds_datadir
                mapping = crds_tools.get_reffiles(crds_dictionary, [CRDS_FILE_TYPES[key2]], download=not offline)
                params[key1][key2] = mapping[CRDS_FILE_TYPES[key2]]
                print("From CRDS, found {} as the {} reference file.".format(mapping[CRDS_FILE_TYPES[key2]], key2))

                # If we grab the IPC reference file from CRDS, then we need to invert it prior to use
                if key2 == 'ipc':
                    ipc_path, ipc_file = os.path.split(mapping['ipc'])
                    inverted_kernel_file = os.path.join(ipc_path, 'Kernel_to_add_IPC_effects_from_{}'.format(ipc_file))
                    # Check to see if the version of the IPC file with an
                    # inverted kernel already exists.
                    if os.path.isfile(inverted_kernel_file):
                        print("Found an existing inverted kernel for this IPC file: {}".format(inverted_kernel_file))
                        params[key1][key2] = inverted_kernel_file
                        params['Reffiles']['invertIPC'] = False
                    else:
                        print("No inverted IPC kernel file found. Kernel will be inverted prior to use.")
                        params['Reffiles']['invertIPC'] = True
    return params


def get_all_reffiles(param_dict):
    """
    """
    # Place import statement here in order to avoid circular import
    # with crd_tools
    from mirage.reference_files import crds_tools

    mapping = crds_tools.get_reffiles(param_dict, list(CRDS_FILE_TYPES.values()))
    return mapping


def get_siaf():
    '''Return a dictionary that holds the contents of the SIAF config
    file.

    Returns
    -------
    siaf_files : dict
        A dictionary that holds the contents of the config file.
    '''
    utils_dir = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    package_dir = os.path.dirname(utils_dir)
    config_file = os.path.join(package_dir, 'config', 'siaf_config.json')

    with open(config_file, 'r') as config_file:
        siaf_files = json.load(config_file)

    return siaf_files


def get_aperture_definition(aperture_name, instrument):
    '''Parses the SIAF to get the definition of a given aperture
    '''

    siaf_file = get_siaf()[instrument]
    siaf_table = asc.read(siaf_file, header_start=1)
    row_ind = list(siaf_table['AperName']).index(aperture_name)
    siaf_row = siaf_table[row_ind]

    detector = aperture_name[3:5]  # Won't work for ALL, A, or B case

    if re.search(r"_F\d\d\d[MWN]", aperture_name):
        filt = aperture_name.split('_')[-1]
    else:
        filt = 'ANY'

    refpix_x = siaf_row['XSciRef']
    refpix_y = siaf_row['YSciRef']
    refpix_v2 = siaf_row['V2Ref']
    refpix_v3 = siaf_row['V3Ref']

    x_size = siaf_row['XSciSize']
    y_size = siaf_row['YSciSize']

    xstart = refpix_x - x_size / 2
    xend = refpix_x + x_size / 2
    ystart = refpix_y - y_size / 2
    yend = refpix_y + y_size / 2

    # Need to add num_amps

    aperture_definition = [aperture_name, detector, filt, xstart, ystart, xend,
                           yend, refpix_x, refpix_y, refpix_v2, refpix_v3]
    return aperture_definition


def get_frame_count_info(numints, numgroups, numframes, numskips, numresets):
    """Calculate information on the number of frames per group and
    per integration

    Parameters
    ----------
    numints : int
        Number of integraitons in the exposure

    numgroups : int
        Number of groups per integration

    numframes : int
        Number of frames averaged together to create a group

    numskips : int
        Number of skipped frames per group

    numresets : int
        Number of detector resets between integrations

    Returns
    -------
    frames_per_group
    """
    frames_per_group = numframes + numskips
    frames_per_integration = numgroups * frames_per_group
    total_frames = numgroups * frames_per_group

    if numints > 1:
        # Frames for all integrations
        total_frames *= numints
        # Add the resets for all but the first integration
        total_frames += (numresets * (numints - 1))

    return frames_per_group, frames_per_integration, total_frames


def get_filter_throughput_file(instrument, filter_name, nircam_module=None, fgs_detector=None):
    """Locate the filter throughput file in the config directory that
    corresponds to the given instrument/module/filter

    Parameters
    ----------
    instrument : str
        Name of JWST instrument.

    filter_name : str
        Name of the filter. Ignored in the case of FGS.

    nircam_module : str
        Name of module. Ignored except in the case of NIRCam

    fgs_detector : str
        Name of FGS detector (e.g. 'GUIDER1'). Ignored in the case of other
        instruments

    Returns
    -------
    throughput_file : str
        Name of ascii file containing the filter throughput curve
    """
    modpath = pkg_resources.resource_filename('mirage', '')
    config_path = os.path.join(modpath, 'config')

    instrument = instrument.lower()
    if instrument == 'nircam':
        tp_file = '{}_nircam_plus_ote_throughput_mod{}_sorted.txt'.format(filter_name.upper(),
                                                                          nircam_module.lower())
    elif instrument == 'niriss':
        tp_file = '{}_niriss_throughput1.txt'.format(filter_name.lower())
    elif instrument == 'fgs':
        det_number = fgs_detector[-1]
        tp_file = 'guider{}_throughput_py.txt'.format(det_number)
    else:
        raise ValueError('ERROR: invalid instrument: {}'.format(instrument))

    throughput_file = os.path.join(config_path, tp_file)
    return throughput_file


def get_subarray_info(params, subarray_table):
    """Find aperture-specific information from the subarray information config file

    Parameters:
    -----------
    params : dict
        Nested dictionary containing Mirage parameters.

    subarray_table : astropy.table.Table
        Table containing subarray definition information. Output from
        read_subarray_definition_file

    Returns:
    --------
    params : dict
        Updated nested dictionary
    """
    if params['Readout']['array_name'] in subarray_table['AperName']:
        mtch = params['Readout']['array_name'] == subarray_table['AperName']
        namps = subarray_table['num_amps'].data[mtch][0]
        if namps != 0:
            params['Readout']['namp'] = int(namps)
        else:
            try:
                if ((params['Readout']['namp'] == 1) or
                   (params['Readout']['namp'] == 4)):
                    print(("CAUTION: Aperture {} can be used with either "
                           "a 1-amp".format(subarray_table['AperName'].data[mtch][0])))
                    print("or a 4-amp readout. The difference is a factor of 4 in")
                    print(("readout time. You have requested {} amps."
                           .format(params['Readout']['namp'])))
                else:
                    raise ValueError(("WARNING: {} requires the number of amps to be 1 or 4. Please set "
                                      "'Readout':'namp' in the input yaml file to one of these values."
                                      .format(params['Readout']['array_name'])))
            except KeyError:
                raise KeyError(("WARNING: 'Readout':'namp' not present in input yaml file. "
                                "{} aperture requires the number of amps to be 1 or 4. Please set "
                                "'Readout':'namp' in the input yaml file to one of these values."
                                .format(params['Readout']['array_name'])))
    else:
        raise ValueError(("WARNING: subarray name {} not found in the "
                          "subarray dictionary {}."
                          .format(params['Readout']['array_name'],
                                  params['Reffiles']['subarray_defs'])))
    return params


def magnitude_to_countrate(observation_mode, magsys, mag, photfnu=None, photflam=None,
                           vegamag_zeropoint=None):
    """Convert a given source magnitude into count rate

    Parameters
    ----------
    observation_mode : str
        e.g. 'imaging', 'wfss'

    magsys : str
        Magnitude system of the input magnitudes. Allowed values are:
        'abmag', 'stmag', 'vegamag'

    mag : float or list
        Magnitude value(s) to transform

    photfnu : float
        Photfnu value that relates count rate and flux density. Only used
        in ABMAG conversions

    photflam : float
        Photflam value that relates count rate and flux density. Only used
        in STMAG conversions

    vegamag_zeropoint : float
        Filter offset in the VEGAMAG (or A0V mag) system. Only used in
        VEGAMAG conversions.


    Returns
    -------
    count_rate : float or list
        Count rate (e/s) corresponding to the input magnutude(s)

    """
    # For NIRISS AMI mode, the count rate values calculated need to be
    # scaled by a factor 0.15/0.84 = 0.17857.  The 0.15 value is the
    # throughput of the NRM, while the 0.84 value is the throughput of the
    # imaging CLEARP element that is in place in the pupil wheel for the
    # normal imaging observations.
    if observation_mode in ['ami', 'AMI']:
        count_scale = 0.15 / 0.84
    else:
        count_scale = 1.
    if magsys.lower() == 'abmag':
        try:
            return count_scale * (10**((mag + 48.599934378) / -2.5) / photfnu)
        except:
            raise ValueError(("AB mag to countrate conversion failed."
                              "magnitude = {}, photfnu = {}".format(mag, photfnu)))
    if magsys.lower() == 'vegamag':
        try:
            return count_scale * (10**((vegamag_zeropoint - mag) / 2.5))
        except:
            raise ValueError(("Vega mag to countrate conversion failed."
                              "magnitude = {}".format(mag)))
    if magsys.lower() == 'stmag':
        try:
            return count_scale * (10**((mag + 21.099934378) / -2.5) / photflam)
        except:
            raise ValueError(("ST mag to countrate conversion failed."
                              "magnitude = {}, photflam = {}".format(mag, photflam)))


def sigma_clipped_mean_value_of_image(array, sigma_value):
    """Return the sigma-clipped mean value of a 2D array

    Parameters
    ----------
    array : numpy.ndarray
        array of values

    sigma_value : int
        Number of sigma to use as clipping threshold. Clipping will be
        symmetric above and below the mean

    Returns
    -------
    meanval : float
        3-sigma iteratively sigma-clipped mean value of ``array``
    """
    good = np.isfinite(array)
    clipped, lo, hi = sigmaclip(array[good], low=sigma_value, high=sigma_value)
    meanval = clipped.mean()
    return meanval


def parse_RA_Dec(ra_string, dec_string):
    """Convert input RA and Dec strings to floats

    Parameters
    ----------

    ra_string : str
        String containing RA information. Can be in the form 10:12:13.2
        or 10h12m13.2s

    dec_string : str
        String containing Declination information. Can be in the form
        10:12:2.4 or 10d12m2.4s

    Returns
    -------

    ra_degrees : float
        Right Ascention value in degrees

    dec_degrees : float
        Declination value in degrees
    """
    try:
        ra_string = ra_string.lower()
        ra_string = ra_string.replace("h", ":")
        ra_string = ra_string.replace("m", ":")
        ra_string = ra_string.replace("s", "")
        ra_string = re.sub(r"\s+", "", ra_string)
        dec_string = dec_string.lower()
        dec_string = dec_string.replace("d", ":")
        dec_string = dec_string.replace("m", ":")
        dec_string = dec_string.replace("s", "")
        dec_string = re.sub(r"\s+", "", dec_string)

        values = ra_string.split(":")
        ra_degrees = 15.*(int(values[0]) + int(values[1])/60. + float(values[2])/3600.)

        values = dec_string.split(":")
        if "-" in values[0]:
            sign = -1
            values[0] = values[0].replace("-", "")
        else:
            sign = +1

        dec_degrees = sign*(int(values[0]) + int(values[1])/60. + float(values[2])/3600.)
        return ra_degrees, dec_degrees
    except:
        raise ValueError("Error parsing RA, Dec strings: {} {}".format(ra_string, dec_string))


def read_pattern_check(parameters):
    # Check the readout pattern that's entered and set nframe and nskip
    # accordingly
    parameters['Readout']['readpatt'] = parameters['Readout']['readpatt'].upper()

    # Read in readout pattern definition file
    # and make sure the possible readout patterns are in upper case
    readpatterns = asc.read(parameters['Reffiles']['readpattdefs'])
    readpatterns['name'] = [s.upper() for s in readpatterns['name']]

    # If the requested readout pattern is in the table of options,
    # then adopt the appropriate nframe and nskip
    if parameters['Readout']['readpatt'] in readpatterns['name']:
        mtch = parameters['Readout']['readpatt'] == readpatterns['name']
        parameters['Readout']['nframe'] = int(readpatterns['nframe'][mtch].data[0])
        parameters['Readout']['nskip'] = int(readpatterns['nskip'][mtch].data[0])
        print(('Requested readout pattern {} is valid. '
               'Using the nframe = {} and nskip = {}'
               .format(parameters['Readout']['readpatt'],
                       parameters['Readout']['nframe'],
                       parameters['Readout']['nskip'])))

        print('in read_pattern_check, nframe and nskip are:')
        print(parameters['Readout']['nframe'], parameters['Readout']['nskip'])
        print(type(parameters['Readout']['nframe']), type(parameters['Readout']['nskip']))
        print('maxiter', parameters['nonlin']['maxiter'], type(parameters['nonlin']['maxiter']))

    else:
        # If the read pattern is not present in the definition file
        # then quit.
        raise ValueError(("WARNING: the {} readout pattern is not defined in {}."
                          .format(parameters['Readout']['readpatt'],
                                  parameters['Reffiles']['readpattdefs'])))
    return parameters


def read_subarray_definition_file(filename):
    """Read in the file that contains a list of subarray names and related information

    Parameters:
    -----------
    filename : str
        Name of the ascii file containing the table of subarray information

    Returns:
    --------
    data : astropy.table.Table
        Table containing subarray information
    """
    try:
        data = asc.read(filename, data_start=1, header_start=0)
    except:
        raise RuntimeError(("Error: could not read in subarray definitions file: {}"
                            .format(filename)))
    return data


def read_yaml(filename):
    """Read the contents of a yaml file into a nested dictionary

    Parameters
    ----------
    filename : str
        Name of yaml file to be read in

    Returns
    -------
    data : dict
        Nested dictionary of file contents
    """
    try:
        with open(filename, 'r') as f:
            data = yaml.load(f, Loader=yaml.FullLoader)
    except FileNotFoundError as e:
            print(e)
    return data


def write_yaml(data, filename):
    """Write a nested dictionary to a yaml file

    Parameters
    ----------
    data : dict
        Nested dictionary of paramters

    filename : str
        Output filename
    """
    with open(filename, 'w') as output:
        yaml.dump(data, output, default_flow_style=False)
