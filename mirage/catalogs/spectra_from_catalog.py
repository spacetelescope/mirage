#! /usr/bin.env python

"""This module reads in a Mirage-formatted source catalog and for each
source, interpolates the provided magnitudes in order to produce a
continuum spectrum. The collection of all sources' spectra are then
saved to an hdf5 file in the proper format to be used by the
NIRCAM_Gsim disperser software.

Authors
-------

    - Bryan Hilbert

Use
---

    This script is intended to be executed as such:

    ::

        from mirage.catalogs import spectra_from_catalog
        spectra_dict = spectra_from_catalog.calculate_flambda('my_catalog_file.cat')
"""
import os
import pkg_resources

from astropy.io import ascii
import astropy.units as u

import hdf5_catalog

MODULE_PATH = pkg_resources.resource_filename('mirage', '')
CONFIG_PATH = os.path.join(MODULE_PATH, 'config')
ZEROPOINT_FILES = {'niriss': os.path.join(CONFIG_PATH, 'niriss_zeropoints.list'),
                   'nircam': os.path.join(CONFIG_PATH, 'NIRCam_zeropoints.list'),
                   'fgs': os.path.join(CONFIG_PATH, 'guider_zeropoints.list')}


def calculate_flambda(source_catalog, magnitude_system, outfile=None):
    """Calculate and add f_lambda columns to an existing source catalog"""
    updated_catalog, filter_information = add_flam_columns(source_catalog, magnitude_system)

    # Create catalog output name if none is given
    if outfile is None:
        cat_dir, cat_file = os.path.split(catalog_file)
        index = cat_file.rindex('.')
        suffix = cat_file[index:]
        outbase = cat_file[0: index] + '_with_flambda' + suffix
        outfile = os.path.join(cat_dir, outbase)

    updated_catalog.write(outfile, format='ascii', overwrite=True)
    print('Catalog updated with f_lambda columns, saved to: {}'.format(outfile))
    return updated_catalog, filter_information


def create_spectra(catalog_with_flambda, filter_info):
    """Create object spectra from an astropy Table source catalog

    Parameters
    ----------
    table : astropy.table.Table
        Source catalog containing f_lambda columns

        OR MAKE THIS THE OVERALL WRAPPER FUNCTION...

    Returns
    -------
    spectra : dict
        Object spectra in a dictionary. Keys are index numbers, values are
        dictionaries containing 'wavelengths' and 'fluxes' keys. Those values
        are numpy arrays of values. These can optionally have astropy units
        attached to them.
    """
    # Calculate f_lambda values and add as columns to the catalog
    #catalog_with_flambda = calculate_flambda(input_catalog)

    flambda_cols = [col for col in catalog_with_flambda.colnames if 'flambda' in col]
    filter_name = [colname.split[1] for colname in flambda_cols]
    print('above: what about fgs? are mag cols for fgs "fgs_guider1_magnitude?" or just "fgs_magnitude"?')
    instrument = [colname.split[0] for colname in flambda_cols]

    pivots = []
    for inst_setup in filter_info:
        f_lambda, zeropoint, pivot = filter_info[inst_setup]
        pivots.append(pivot)

    sorted_indexes = np.argsort(np.array(pivots))
    wavelengths = pivots[sorted_indexes]
    spectra = OrderedDictionary({})
    for idx in range(len(catalog_with_flambda)):
        source = catalog_with_flambda[idx]
        index = source['index']
        flux = source[flambda_cols]
        sorted_fluxes = flux[sorted_indexes]
        spectra[index] = {'wavelengths': wavelengths, 'fluxes': sorted_indexes}

    return spectra


def overall_wrapper(catalog_file, input_spectra=None, input_spectra_file=None):
    """Overall wrapper function

    Parameters
    ----------
    source_catalog : str
        Name of Mirage-formatted source catalog file

    input_spectra : dict
        Dictionary containing spectra for some/all targets

    input_specctra_file : str
        Name of an hdf5 file containing spectra for some/all targets

    Returns
    -------
    ??
    """
    all_input_spectra = {}
    # Read in input spectra from file, add to all_input_spectra dictionary
    if input_spectra_file is not None:
        spectra_from_file = hdf5_catalog.open(input_spectra_file)
        all_input_spectra = {**all_input_spectra, **spectra_from_file}

    # If both an input spectra file and dictionary are given, combine into
    # a single dictionary. Running in this order means that dictionary
    # inputs override inputs from a file.
    if input_spectra is not None:
        all_input_spectra = {**all_input_spectra, **input_spectra}

    # Read in input catalog
    ascii_catalog, mag_sys = read_catalog(catalog_file)
    catalog, filter_info = calculate_flambda(ascii_catalog, mag_sys)

    # Check to see if input spectra are normalized and if so, scale using
    # magnitudes from catalog_file
    for key in all_input_spectra:
        check_if_normalized()
        rescale_if_necessary()

    # For sources in catalog_file but not in all_input_spectra, use the
    # magnitudes in the catalog_file, interpolate/extrapolate as necessary
    # and create continuum spectra
    for line in ascii_catalog:
        match = line['index'] == all_input_spectra.keys
        if not any(match):
            print('Source {} is in ascii catalog but not input spectra.'.format(line['index']))
            continuum = create_spectra(line)
            all_input_spectra[line['index']] = continuum

    # For convenience, reorder the sources by index number
    spectra = OrderedDictionary({})
    bleck()

    # Save the source spectra in an hdf5 file
    hdf5_catalog.save(spectra, spectra_filename, wavelength_unit='micron', flux_unit='flam')
    print('Spectra catalog file saved to {}'.format(spectra_filename))




 when creating spectrum, you can mix nircam module X and niriss and guider N, but you
 cannot mix nircam mod A and nircam mod B nor guider 1 and guider 2. ACTUALLY, MAYBE YOU
 CAN MIX EVERYTHING, AS LONG AS THE MAGNITUDES ARE CORRECT...






def read_catalog(filename):
    """Read in the Mirage-formatted catalog

    Parameters
    ----------
    filename : str
        Name of catalog file

    Returns
    -------
    catalog : astropy.table.Table
        Catalog contents

    mag_sys : str
        Magnitude system (e.g. 'abmag', 'stmag', 'vegamag')
    """
    catalog = ascii.read(filename)

    # Check to be sure index column is present
    if 'index' not in catalog.colnames:
        raise ValueError(('WARNING: input catalog must have an index column '
                          'so that object spectra can be associated with the '
                          'correct object in the seed image/segmentation map.'))
    min_index_value = min(catalog['index'].data)
    if min_index_value < 1:
        raise ValueError(("WARNING: input catalog cannot have index values less"
                          "than 1, for the purposes of creating the segmentation map."))

    # Get the magnitude system
    mag_sys = [l.lower() for l in catalog.meta['comments'][0:4] if 'mag' in l][0]
    if mag_sys not in ['abmag', 'stmag', 'vegamag']:
        raise ValueError(('WARNING: unrecognized magnitude system {}. Must be'
                          'one of "abmag", "stmag", or "vegamag".').format(mag_sys))
    return catalog, mag_sys


def add_flam_columns(cat, mag_sys):
    """Convert magnitudes to flux densities in units of f_lambda and add
    to catalog.

    Parameters
    ----------
    cat : astropy.table.Table
        Catalog table (output from read_catalog)

    mag_sys : str
        Magnitude system of the data (e.g. 'abmag')

    Returns
    -------
    cat : astropy.table.Table
        Modified table with "<instrument>_<filter>_flam" columns added
    """
    # Get all of the magnitude columns
    magnitude_columns = [col for col in cat.colnames if 'magnitude' in col]

    # Organize magnitude columns by instrument to minimize file i/o when
    # getting photflam and zeropoint values
    nircam_mag_cols = []
    niriss_mag_cols = []
    fgs_mag_cols = []
    for col in magnitude_columns:
        if 'nircam' in col.lower():
            nircam_mag_cols.append(col)
        elif 'niriss' in col.lower():
            niriss_mag_cols.append(col)
        elif 'fgs' in col.lower():
            fgs_mag_cols.append(col)

    # Get PHOTFLAM and filter zeropoint values
    parameters = {}
    if len(nircam_mag_cols) > 0:
        nrc_params = get_filter_info(nircam_mag_cols, mag_sys)
        parameters = {**parameters, **nrc_params}
    if len(niriss_mag_cols) > 0:
        nrs_params = get_filter_info(niriss_mag_cols, mag_sys)
        parameters = {**parameters, **nrs_params}
    if len(fgs_mag_cols) > 0:
        fgs_params = get_filter_info(fgs_mag_cols, mag_sys)
        parameters = {**parameters, **fgs_params}

    # Create a corresponding f_lambda column for each magnitude column
    for key in parameters:
        magnitude_values = cat[key].data
        flam_values = convert_to_flam(key, magnitude_values, parameters[key])
        new_column_name = key.replace('magnitude', 'flam')
        cat[new_column_name] = flam_values
    return cat, parameters


def get_filter_info(column_names, magsys):
    """Given a list of catalog columns names (e.g. 'nircam_f480m_magnitude')
    get the corresponding PHOTFLAM and filter zeropoint values

    Parameters
    ----------
    column_names : list
        List of Mirage catalog column names

    magsys : str
        Magnitude system (e.g. 'abmag')

    Returns
    -------
    info : dict
        Dictionary of values
        (e.g. info['nircam_f480m_magnitude'] = (photflam, zeropoint))
    """
    # Get the correct zeropoint file name
    instrument = column_names[0].split('_')[0].lower()
    zp_file = ZEROPOINT_FILES[instrument]

    # Read in the file
    zp_table = ascii.read(zp_file)

    magsys = magsys.upper()
    info = {}

    if instrument == 'nircam':
        for entry in column_names:
            filter_name = entry.split('_')[1].upper()
            match = zp_table['Filter'] == filter_name
            modules = zp_table['Module'].data[match]
            zp = zp_table[magsys].data[match]
            photflam = zp_table['PHOTFLAM'].data[match]
            pivot = zp_table['Pivot_wave']
            for mod, flam, z, pvt in zip(modules, photflam, zp, pivot):
                new_col_name = entry + '_mod{}'.format(mod)
                info[new_col_name] = (flam, z, pvt)

    elif instrument == 'niriss':
        for entry in column_names:
            filter_name = entry.split('_')[1].upper()
            match = zp_table['Filter'] == filter_name
            zp = zp_table[magsys].data[match][0]
            photflam = zp_table['PHOTFLAM'].data[match][0]
            pivot = zp_table['Pivot_wave']
            info[entry] = (photflam, zp, pivot)

    elif instrument == 'fgs':
        for line in zp_table:
            detector = line['Detector']
            zp = line[magsys]
            photflam = line['PHOTFLAM']
            pivot = line['Pivot_wave']
            new_col_name = column_names[0] + '_{}'.format(detector.lower())
            info[new_col_name] = (photflam, zp, pivot)

    return info


def convert_to_flam(colname, magnitudes, param_tuple):
    """Convert the magnitude values for a given magnitude column into
    units of f_lambda.

    Parameters
    ----------
    colname : str
        Name of magnitude column (e.g. 'nircam_f480m_magnitude_modA')

    magnitudes : list
        List of magnitude values for the column

    param_tuple : tup
        Tuple of (photflam, zeropoint) for the given filter

    Returns
    -------
    flambda : list
        List of f_lambda values corresponding to the list of input magnitudes
    """
    see mag_to_countrate in catalog_seed_image
    photflam, zeropoint, pivot = param_tuple
    f_lam = 10**((magnitudes - zeropoint)/-2.5) * photflam
    return f_lam
