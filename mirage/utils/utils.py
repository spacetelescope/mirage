"""Various utility functions.

Authors
-------
    - Lauren Chambers

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
import re

from astropy.io import ascii as asc


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
        xs = xdim
        ys = ydim
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
        xs = ydim
        ys = xdim
        colpad = 12

        # Fullframe
        if amps == 4:
            rowpad = 1
            fullpad = 1
        else:
            rowpad = 2
            fullpad = 0

    elif instrument == 'fgs':
        xs = ydim
        ys = xdim
        colpad = 6
        if 'acq1' in aperture.lower():
            colpad = 12
        rowpad = 1
        if amps == 4:
            fullpad = 1
        else:
            fullpad = 0

    return ((1.0 * xs / amps + colpad) * (ys + rowpad) + fullpad) * 1.e-5


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
            params['Readout']['namp'] = namps
        else:
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
