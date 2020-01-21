#! /usr/bin/env python

"""
Tools for generating Mirage-compatible catalogs from surveys
Note that once you have mirage-formatted catalogs, there is an "add_catalog" function
in catalog_generator.py that can combine catalogs
"""

from collections import OrderedDict
import copy
import math
import numpy as np
import os
import pkg_resources

from astropy.coordinates import SkyCoord, Galactic
from astropy.io import ascii
from astropy.table import Table, join
import astropy.units as u
from astroquery.gaia import Gaia
from astroquery.irsa import Irsa
from pysiaf.utils.projection import deproject_from_tangent_plane

from mirage.apt.apt_inputs import get_filters
from mirage.catalogs.catalog_generator import PointSourceCatalog, GalaxyCatalog, \
    ExtendedCatalog, MovingPointSourceCatalog, MovingExtendedCatalog, \
    MovingSersicCatalog
from mirage.utils.utils import ensure_dir_exists


def create_basic_exposure_list(xml_file, pointing_file):
    """Create an exposure table from an APT (xml and pointing) file.
    This is a shortened version of what is done in apt_inputs, in order to
    collect the information needed to generate source catalogs.

    Parameters
    ----------
    xml_file : str
        Name of xml file exported from a given APT proposal

    pointing_file : str
        Name of pointing file exported from a given APT proposal

    Returns
    -------
    apt.exposure_tab : dict
        Dictionary containing observation information from the xml and
        pointing files
    """
    from mirage.yaml import yaml_generator
    from mirage.apt import apt_inputs

    info = yaml_generator.SimInput(input_xml=xml_file, pointing_file=pointing_file,
                                   offline=True)

    apt = apt_inputs.AptInput(input_xml=xml_file, pointing_file=pointing_file)
    apt.observation_list_file = info.observation_list_file
    apt.apt_xml_dict = info.apt_xml_dict
    apt.output_dir = info.output_dir
    apt.create_input_table()
    return apt.exposure_tab


def for_proposal(xml_filename, pointing_filename, point_source=True, extragalactic=True,
                 catalog_splitting_threshold=0.12, besancon_catalog_file=None,
                 ra_column_name='RAJ2000', dec_column_name='DECJ2000', out_dir=None,
                 save_catalogs=True, galaxy_seed=None):
    """
    Given a pointing dictionary from an APT file, generate source catalogs
    that cover all of the coordinates specifired.

    Parameters
    ----------
    xml_filename : str
        Name of the xml file of the APT proposal

    pointing_filename : str
        Name of the pointing file associated with the APT proposal

    point_source : bool
        If True, a point source catalog is created

    extragalactic : bool
        If True, a catalog of background galaxies is created

    catalog_splitting_threshold : float
        Maximum distance in degrees between targets that will be contained
        in a given catalog. Sources farther than this distance will be placed
        into a separate catalog.

    besancon_catalog_file : str
        Name of ascii catalog containing background stars. The code was
        developed around this being a catalog output by the Besaoncon
        model (via ``besancon()``), but it can be any ascii catalog
        of sources as long as it contains the columns specified in
        the ``catalog`` parameter of ``johnson_catalog_to_mirage_catalog()``
        If None, the Besancon step will be skipped and a catalog will be
        built using only GAIA/2MASS/WISE

    ra_column_name : str
        Name of the column within ``besancon_catalog_file`` containing the
        right ascension of the sources.

    dec_column_name : str
        Name of the column within ``besancon_catalog_file`` containing the
        declination of the sources

    out_dir : str
        Directory in which to save catalog files

    save_catalogs : bool
        If True, save the catalogs to ascii files, in addition to returning
        them. If False, the catalog objects are returned, but not saved.

    galaxy_seed : int
        Seed to use in the random number generator used in galaxy_background

    Returns
    -------
    ptsrc_catalog_list : list
        List of Mirage point source catalog objects

    galaxy_catalog_list : list
        List of Mirage galaxy source catalog objects

    ptsrc_catalog_names : list
        List of filenames of the saved point source catalogs

    galaxy_catalog_names : list
        List of filenames of the saved galaxy catalogs
    """
    pointing_dictionary = create_basic_exposure_list(xml_filename, pointing_filename)
    instrument_filter_dict = get_filters(pointing_dictionary)

    threshold = catalog_splitting_threshold * u.deg
    ra_apertures = np.array(pointing_dictionary['ra_ref'] * u.deg)
    dec_apertures = np.array(pointing_dictionary['dec_ref'] * u.deg)

    ra_targets = np.array([np.float(num) for num in pointing_dictionary['ra']] * u.deg)
    dec_targets = np.array([np.float(num) for num in pointing_dictionary['dec']] * u.deg)
    mapped = np.array([False] * len(ra_targets))
    index = np.arange(len(ra_targets))

    targets = SkyCoord(ra=pointing_dictionary['ra'], dec=pointing_dictionary['dec'], frame='icrs', unit=u.deg)
    apertures = []
    for ra, dec in zip(pointing_dictionary['ra_ref'], pointing_dictionary['dec_ref']):
        apertures.append(SkyCoord(ra=ra, dec=dec, frame='icrs', unit=u.deg))

    observation_ids = np.array(pointing_dictionary['ObservationID'])
    unique_observation_ids = list(set(observation_ids))
    mapped_observations = {}
    ptsrc_catalog_mapping = {}
    ptsrc_catalog_list = []
    ptsrc_catalog_names = []
    galaxy_catalog_mapping = {}
    galaxy_catalog_list = []
    galaxy_catalog_names = []
    for observation in unique_observation_ids:
        # Skip observations that have already been examined because their
        # targets match those in a previous observation
        mapped_observation = False
        for values in mapped_observations.values():
            if observation in values:
                mapped_observation = True
            else:
                mapped_observation = False
        if mapped_observation:
            continue

        match = np.where(observation_ids == observation)[0]
        target_ra = ra_targets[match[0]]
        target_dec = dec_targets[match[0]]
        target_pos = targets[match[0]]

        dithers = ((ra_targets == target_ra) & (dec_targets == target_dec))
        all_dithers = np.array(copy.deepcopy(dithers))
        aperture_position = [(r, d) for r, d in zip(np.array(pointing_dictionary['ra_ref'])[dithers],
                                                    np.array(pointing_dictionary['dec_ref'])[dithers])]
        all_obs = observation_ids[dithers]
        mapped_observations[str(observation)] = list(set(all_obs))

        # Create a list of relative SkyCoords
        targ_pos_relative = [pos.transform_to(targets[match[0]].skyoffset_frame()) for pos in targets[dithers]]

        separation_distances = target_pos.separation(targets)
        nearby_targets = ((separation_distances < threshold) & (mapped == False))

        if any(nearby_targets):
            nearby_target_ra = np.array(pointing_dictionary['ra'])[nearby_targets]
            nearby_target_dec = np.array(pointing_dictionary['dec'])[nearby_targets]
            near_positions = np.array([(ra_val, dec_val) for ra_val, dec_val in zip(nearby_target_ra,
                                                                                    nearby_target_dec)])
            unique_nearby = np.unique(near_positions, axis=0)
            for near_coords in unique_nearby:
                near_ra = float(near_coords[0])
                near_dec = float(near_coords[1])
                tmp_dither = ((ra_targets == near_ra) & (dec_targets == near_dec))
                all_dithers = all_dithers | np.array(tmp_dither)
                tmp_aperture_position = [(r, d) for r, d in zip(np.array(pointing_dictionary['ra_ref'])[tmp_dither],
                                                                np.array(pointing_dictionary['dec_ref'])[tmp_dither])]
                aperture_position.extend(tmp_aperture_position)

                # Add to list of relative SkyCoords
                relative = [pos.transform_to(targets[match[0]].skyoffset_frame()) for pos in targets[tmp_dither]]
                targ_pos_relative.extend(relative)

                tmp_obs = observation_ids[tmp_dither]
                temp = mapped_observations[str(observation)]
                temp.extend(list(set(tmp_obs)))
                mapped_observations[str(observation)] = list(set(temp))

        longitudes = [pos.lon.arcsec for pos in targ_pos_relative]
        latitudes = [pos.lat.arcsec for pos in targ_pos_relative]
        ralist = np.array([pos.ra.degree for pos in targets[all_dithers]])
        declist = np.array([pos.dec.degree for pos in targets[all_dithers]])

        pad = 0.062 * 1024 * 1.5
        delta_lon = np.max(longitudes) - np.min(longitudes) + pad
        delta_lat = np.max(latitudes) - np.min(latitudes) + pad
        full_width = max(delta_lon, delta_lat)

        if ((np.min(ralist) <= (catalog_splitting_threshold)) & (np.max(ralist) >= (360.-catalog_splitting_threshold))):
            ralist[ralist >= 360.-catalog_splitting_threshold] -= 360.

        mean_ra = np.mean(ralist)
        mean_dec = np.mean(declist)

        # Generate a string listing which observations the catalog covers
        for_obs = mapped_observations[str(observation)]
        for_obs_str = ''
        for value in for_obs:
            for_obs_str += str(value)+'_'
        for_obs_str = for_obs_str[0:-1]
        xml_base = os.path.basename(xml_filename).split('.xml')[0]
        if out_dir is None:
            out_dir = os.path.dirname(xml_filename)

        if point_source:
            for i, instrument in enumerate(instrument_filter_dict):
                print('\n--- Creating {} point source catalog ---'.format(instrument))
                filter_list = instrument_filter_dict[instrument]
                tmp_cat, tmp_filters = get_all_catalogs(mean_ra, mean_dec, full_width,
                                                        besancon_catalog_file=besancon_catalog_file,
                                                        instrument=instrument, filters=filter_list,
                                                        ra_column_name=ra_column_name, dec_column_name=dec_column_name)
                if i == 0:
                    ptsrc_cat = copy.deepcopy(tmp_cat)
                else:
                    ptsrc_cat = combine_catalogs(ptsrc_cat, tmp_cat)

            if save_catalogs:
                ptsrc_catalog_name = 'ptsrc_for_{}_observations_{}.cat'.format(xml_base, for_obs_str)

                # Populate dictionary listing ptsrc catalog name associated with each
                # observation number
                ptsrc_catalog_mapping[str(observation)] = ptsrc_catalog_name
                for value in for_obs:
                    ptsrc_catalog_mapping[str(value)] = ptsrc_catalog_name

                ensure_dir_exists(out_dir)
                full_catalog_path = os.path.join(out_dir, ptsrc_catalog_name)
                ptsrc_cat.save(full_catalog_path)
                print('\nPOINT SOURCE CATALOG SAVED: {}'.format(full_catalog_path))
                ptsrc_catalog_names.append(full_catalog_path)

            ptsrc_catalog_list.append(ptsrc_cat)
        else:
            ptsrc_cat = None

        if extragalactic:
            for i, instrument in enumerate(instrument_filter_dict):
                print('\n--- Creating {} extragalactic catalog ---'.format(instrument))
                filter_list = instrument_filter_dict[instrument]
                tmp_cat, tmp_seed = galaxy_background(mean_ra, mean_dec, 0., full_width, instrument,
                                                      filter_list, boxflag=False, brightlimit=14.0,
                                                      seed=galaxy_seed)

                if i == 0:
                    galaxy_cat = copy.deepcopy(tmp_cat)
                else:
                    galaxy_cat = combine_catalogs(galaxy_cat, tmp_cat)


            if save_catalogs:
                gal_catalog_name = 'galaxies_for_{}_observations_{}.cat'.format(xml_base, for_obs_str)

                # Populate dictionary listing galaxy catalog name associated with each
                # observation number
                galaxy_catalog_mapping[str(observation)] = gal_catalog_name
                for value in for_obs:
                    galaxy_catalog_mapping[str(value)] = gal_catalog_name

                full_catalog_path = os.path.join(out_dir, gal_catalog_name)
                galaxy_cat.save(full_catalog_path)
                print('\nGALAXY CATALOG SAVED: {}'.format(full_catalog_path))
                galaxy_catalog_names.append(full_catalog_path)

            galaxy_catalog_list.append(galaxy_cat)
        else:
            galaxy_cat = None

    return (ptsrc_catalog_list, galaxy_catalog_list, ptsrc_catalog_names, galaxy_catalog_names,
            ptsrc_catalog_mapping, galaxy_catalog_mapping)


def query_2MASS_ptsrc_catalog(ra, dec, box_width):
    """
    Query the 2MASS All-Sky Point Source Catalog in a square region around
    the RA and Dec provided. Box width must be in units of arcseconds

    Parameters
    ----------
    ra : float or str
        Right ascention of the center of the catalog. Can be decimal degrees
        or HMS string

    dec : float or str
        Declination of the center of the catalog. Can be decimal degrees of
        DMS string

    box_width : float
        Width of the box in arcseconds containing the catalog.

    Returns
    -------
    query_table : astropy.table.Table
        Catalog composed of MaskedColumns containing 2MASS sources

    magnitude_column_names : list
        List of column header names corresponding to columns containing
        source magnitude
    """
    # Don't artificially limit how many sources are returned
    Irsa.ROW_LIMIT = -1

    ra_dec_string = "{}  {}".format(ra, dec)
    query_table = Irsa.query_region(ra_dec_string, catalog='fp_psc', spatial='Box',
                                    width=box_width * u.arcsec)

    # Exclude any entries with missing RA or Dec values
    radec_mask = filter_bad_ra_dec(query_table)
    query_table = query_table[radec_mask]

    # Column names of interest
    magnitude_column_names = ['j_m', 'h_m', 'k_m']
    return query_table, magnitude_column_names


def get_wise_ptsrc_catalog(ra, dec, box_width):
    """Wrapper around WISE query and creation of mirage-formatted catalog

    Parameters
    ----------
    ra : float or str
        Right ascention of the center of the catalog. Can be decimal degrees
        or HMS string

    dec : float or str
        Declination of the center of the catalog. Can be decimal degrees of
        DMS string
    """
    wise_cat, wise_mag_cols = query_WISE_ptsrc_catalog(ra, dec, box_width)
    wise_mirage = mirage_ptsrc_catalog_from_table(wise_cat, 'WISE', wise_mag_cols)
    return wise_mirage, wise_cat


def query_WISE_ptsrc_catalog(ra, dec, box_width):
    """Query the WISE All-Sky Point Source Catalog in a square region around the RA and Dec
    provided. Box width must be in units of arcseconds

    Parameters
    ----------
    ra : float or str
        Right ascention of the center of the catalog. Can be decimal degrees
        or HMS string

    dec : float or str
        Declination of the center of the catalog. Can be decimal degrees of
        DMS string

    box_width : float
        Width of the box in arcseconds containing the catalog.

    Returns
    -------
    query_table : astropy.table.Table
        Catalog composed of MaskedColumns containing WISE sources

    magnitude_column_names : list
        List of column header names corresponding to columns containing source magnitude
    """
    # Don't artificially limit how many sources are returned
    Irsa.ROW_LIMIT = -1

    ra_dec_string = "{}  {}".format(ra, dec)
    query_table = Irsa.query_region(ra_dec_string, catalog='allsky_4band_p3as_psd', spatial='Box',
                                    width=box_width * u.arcsec)

    # Exclude any entries with missing RA or Dec values
    radec_mask = filter_bad_ra_dec(query_table)
    query_table = query_table[radec_mask]

    # Column names of interest
    magnitude_column_names = ['w1mpro', 'w2mpro', 'w3mpro', 'w4mpro']
    return query_table, magnitude_column_names


def mirage_ptsrc_catalog_from_table(table, instrument, mag_colnames, magnitude_system='vegamag'):
    """Create a mirage-formatted point source catalog from an input astropy
    table (e.g. from one of the query functions), along with the magnitude
    column names of interest

    Parameters
    ----------
    table : astropy.table.Table
        Source catalog from (e.g.) Gaia or 2MASS query

    instrument : str
        Unique identifier for where data came from (e.g. '2MASS', 'WISE', 'nircam_f200w')

    mag_colnames : list
        List of strings corresponding to the columns in 'table' that contain magnitude values

    magnitude_system : str
        This is the label for the magnitude system, 'vegamag', 'abmag', or 'stmag'.
    """
    cat = PointSourceCatalog(ra=table['ra'].data.data, dec=table['dec'].data.data)

    for magcol in mag_colnames:
        data = table[magcol].filled().data
        cat.add_magnitude_column(data, instrument=instrument, filter_name=magcol,
                                 magnitude_system=magnitude_system)
    return cat


def get_2MASS_ptsrc_catalog(ra, dec, box_width):
    """Wrapper around 2MASS query and creation of mirage-formatted catalog

    Parameters
    ----------
    ra : float or str
        Right ascention of the center of the catalog. Can be decimal degrees or HMS string

    dec : float or str
        Declination of the center of the catalog. Can be decimal degrees of DMS string

    box_width : float
        Width of the box in arcseconds containing the catalog.

    Returns
    -------
    twomass : tup
        2-element tuple. The first is a PointSourceCatalog object containing
        the 2MASS source information. The second is an astropy.table.Table
        object containing the exact query results from 2MASS. Note that these
        outputs will contain only JHK magnitudes for the returned sources.
    """
    twomass_cat, twomass_mag_cols = query_2MASS_ptsrc_catalog(ra, dec, box_width)
    twomass_mirage = mirage_ptsrc_catalog_from_table(twomass_cat, '2MASS', twomass_mag_cols)
    return twomass_mirage, twomass_cat


def twoMASS_plus_background(ra, dec, box_width, kmag_limits=(17, 29), email='', seed=None):
    """Convenience function to create a catalog from 2MASS and add a population of
    fainter stars. In this case, cut down the magnitude limits for the call to the
    Besancon model so that we don't end up with a double population of bright stars

    Parameters
    ----------
    ra : float or str
        Right ascention of the center of the catalog. Can be decimal degrees
        or HMS string

    dec : float or str
        Declination of the center of the catalog. Can be decimal degrees of
        DMS string

    box_width : float
        Width of the box in arcseconds containing the catalog.

    kmag_limits : tup
        Tuple of minimum and maximum K magnitude values to use in Besancon
        model query for background stars

    email : str
        A valid email address is necessary to perform a Besancon query

    seed : int
        Seed to use in the random number generator when choosing RA and
        Dec values for Besancon sources.
    """
    two_mass, twomass_cat = get_2MASS_ptsrc_catalog(ra, dec, box_width)
    background, background_cat = besancon(ra, dec, box_width, coords='ra_dec', email=email, seed=seed)
    two_mass.add_catalog(background)
    return two_mass


def get_all_catalogs(ra, dec, box_width, besancon_catalog_file=None, instrument='NIRISS', filters=[],
                     ra_column_name='RAJ2000', dec_column_name='DECJ2000'):
    """
    This is a driver function to query the GAIA/2MASS/WISE catalogues
    plus the Besancon model and combine these into a single JWST source list.
    The routine cuts off the bright end of the Besancon model to match the input
    2MASS catalogue, so as to avoid getting too many bright stars in the
    output catalogue.
    For the observed values, the code interpolates the photon-weighted mean
    flux density values in wavelength to get the NIRCam/NIRISS/Guider magnitudes.
    This process is nominally OK when there is extinction, but only when the
    spectrum is relatively smooth (so it is suspect for stars with very deep
    moelcular bands over the near-infrared wavelengths).
    For the Bescancon values, the VJHKL magnitudes are matched to Kurucz models
    and then the extinction is applied.  This process is limited by the range
    of Kurucz models considered., and to some degree by the simple extinction
    model used in the Besancon code.

    Parameters
    ----------
        ra : float or str
            Right ascension of the target field in degrees or HMS

        dec : float or str
            Declination of the target field in degrees or DMS

        box_width : float
            Size of the (square) target field in arc-seconds

        besancon_catalog_file : str
            Name of ascii catalog containing background stars. The code was
            developed around this being a catalog output by the Besaoncon
            model (via ``besancon()``), but it can be any ascii catalog
            of sources as long as it contains the columns specified in
            the ``catalog`` parameter of ``johnson_catalog_to_mirage_catalog()``
            If None, the Besancon step will be skipped and a catalog will be
            built using only GAIA/2MASS/WISE

        instrument : str
            One of "all", "NIRISS", "NIRCam", or "Guider"

        filters : list
            Either an empty list (which gives all filters) or a list
            of filter names (i.e. F090W) to be calculated.

        ra_column_name : str
            Name of the column within ``besancon_catalog_file`` containing the
            right ascension of the sources.

        dec_column_name : str
            Name of the column within ``besancon_catalog_file`` containing the
            declination of the sources

    Returns
    -------
        source_list : mirage.catalogs.create_catalogs.PointSourceCatalog
            A table with the filter magnitudes

        filter_names : list
            A list of the filter name header strings for writing to
            an output file.
    """
    if isinstance(ra, str):
        pos = SkyCoord(ra, dec, frame='icrs')
        outra = pos.ra.deg
        outdec = pos.dec.deg
    else:
        outra = ra
        outdec = dec
    filter_names = make_filter_names(instrument, filters)

    gaia_cat, gaia_mag_cols, gaia_2mass, gaia_2mass_crossref, gaia_wise, \
        gaia_wise_crossref = query_GAIA_ptsrc_catalog(outra, outdec, box_width)
    twomass_cat, twomass_cols = query_2MASS_ptsrc_catalog(outra, outdec, box_width)
    wise_cat, wise_cols = query_WISE_ptsrc_catalog(outra, outdec, box_width)

    if besancon_catalog_file is not None:
        filter_dict = {instrument: filters}
        besancon_jwst = johnson_catalog_to_mirage_catalog(besancon_catalog_file, filter_dict, ra_column_name=ra_column_name,
                                                          dec_column_name=dec_column_name, magnitude_system='vegamag')

    # Combine data from GAIA/2MASS/WISE to create single catalog with JWST filters
    observed_jwst = combine_and_interpolate(gaia_cat, gaia_2mass, gaia_2mass_crossref, gaia_wise,
                                            gaia_wise_crossref, twomass_cat, wise_cat, instrument, filters)

    if besancon_catalog_file is not None:
        print('Adding %d sources from Besancon to %d sources from the catalogues.' % (len(besancon_jwst.ra),
                                                                                      len(observed_jwst.ra)))
        source_list = combine_catalogs(observed_jwst, besancon_jwst)
    else:
        source_list = observed_jwst
    return source_list, filter_names


def transform_johnson_to_jwst(johnson_cat, filter_names):
    """
    Given the output from a Besancon model, transform to JWST magnitudes by
    making a best match to the standard BOSZ model colours for the VJHKL
    magnitudes, then taking the BOSZ results for the JWST filters.  Apply
    any ISM extinction after the matching.
    One possible improvement would be to interpolate between the best matches
    rather than just taking the best one.

    Parameters
    ----------
    johnson_cat : astropy.table.Table
        Table containing sources with magnitudes in Johnson filters.
        V, J, H, and K are required, along with Av (ISM extinction)

    filter_names : list
        The list of the requried output NIRCam/NIRISS/Guider filter names.

    Returns
    -------
    johnson_cat : astropy.table.Table
        Modified table with new columns for JWST filter-based magnitudes
    """
    standard_magnitudes, standard_values, standard_filters, standard_labels = read_standard_magnitudes()
    nstars = len(johnson_cat['V'].data)
    nfilters = len(filter_names)

    # out_magnitudes is a two-dimensional float numpy array with the star
    # number in the first dimension and the estimated JWST magnitudes
    # in the second dimension. A value of None is returned if the matching
    # fails.  This should not happen with the regular inputs.
    out_magnitudes = np.zeros((nstars, nfilters), dtype=np.float32)
    inds = crossmatch_filter_names(filter_names, standard_filters)
    in_magnitudes = np.zeros((4), dtype=np.float32)

    in_filters = ['Johnson V', 'Johnson J', 'Johnson H', 'Johnson K']
    for loop in range(nstars):
        # Exclude any bad values returned by the Besancon query
        if johnson_cat['V'].data[loop] < 90:
            in_magnitudes[0] = johnson_cat['V'].data[loop]
            in_magnitudes[1] = johnson_cat['J'].data[loop]
            in_magnitudes[2] = johnson_cat['H'].data[loop]
            in_magnitudes[3] = johnson_cat['K'].data[loop]
            newmags = match_model_magnitudes(in_magnitudes, in_filters, standard_magnitudes, standard_values,
                                             standard_filters, standard_labels)
            if newmags is None:
                return None
            newmags = johnson_cat['Av'].data[loop] * standard_values[3, :] + newmags
            out_magnitudes[loop, :] = newmags[inds]
        else:
            out_magnitudes[loop, :] = np.repeat(99., len(inds))

    # Create new columns for the JWST filter-based magnitudes
    for magnitudes, filter_name in zip(out_magnitudes.transpose(), filter_names):
        johnson_cat[filter_name] = magnitudes

    return johnson_cat


def catalog_colors_to_vjhk(color_catalog):
    """Given a Table containing a catalog of sources with V magnitudes as
    well as colors, generate J, H, and K magnitudes for all sources.

    Parameters
    ----------

    color_catalog : astropy.table.Table
        Table of sources. Should have columns named: 'K', 'V-K', 'J-H',
        and 'J-K'

    Returns
    -------

    color_catalog : astropy.table.Table
        Modified Table with 'J', 'H', 'K' columns added
    """
    # Calculate magnitudes
    kmags = color_catalog['K'].data
    vmags = (color_catalog['K'] + color_catalog['V-K']).data
    jmags = kmags + color_catalog['J-K'].data
    hmags = jmags - color_catalog['J-H'].data

    # Add to table
    color_catalog['J'] = jmags
    color_catalog['H'] = hmags
    color_catalog['V'] = kmags
    return color_catalog


def johnson_catalog_to_mirage_catalog(catalog_file, filters, ra_column_name='RAJ2000', dec_column_name='DECJ2000',
                                      magnitude_system='abmag', output_file=None):
    """Create a Mirage-formatted catalog containing sources from a query of the Besancon model

    Parameters
    ----------
    catalog_file : str
        Name of ascii file containing the input catalog. The catalog must have
        RA and Dec columns (whose names can be specified using ``ra_column_name``
        and ``dec_column_name``), as well as a ``K`` column containing source
        magnitudes in the K band, and an ``Av`` column containing ISM extinction.
        It must also have either: 1) ``J``, ``H``,  and ``V`` columns with
        appropriate magnitudes, or 2) ``V-K``, ``J-K``, and ``J-H`` columns.

    filters : dict
        Dictionary with keys equal to JWST instrument names, and values
        that are lists of filter names. Besancon sources will have magnitudes
        transformed into thiese filters

    ra_column_name : str
        Name of the column in the input catalog that contains the Right Ascension
        values for the sources

    dec_column_name : str
        Name of the column in the input catalog that contains the Declination
        values for the sources

    magnitude_system : str
        Magnitude system of the values in the catalog. Options are 'abmag', 'stmag',
        and 'vegamag'

    output_file : str
        If not None, save the catalog to a file with this name

    Returns
    -------
    transformed_besancon_cat : mirage.catalogs.create_catalog.PointSourceCatalog
        Catalog containing Besancon model stars with magnitudes in requested
        JWST filters
    """
    # Check the input magnitude system
    if magnitude_system not in ['abmag', 'stmag', 'vegamag']:
        raise ValueError(("ERROR: magnitude_system for {} must be one of: 'abmag', 'stmag', 'vegamag'.")
                         .format(catalog_file))

    # Create list of filter names from filters dictionary
    all_filters = []
    for instrument in filters:
        filt_list = make_filter_names(instrument.lower(), filters[instrument])
        all_filters.extend(filt_list)

    # Read in the input catalog
    catalog = ascii.read(catalog_file)

    # Quick check to be sure required columns are present
    req_cols = ['K', 'Av', ra_column_name, dec_column_name]
    for colname in req_cols:
        if colname not in catalog.colnames:
            raise ValueError(('ERROR: Required column {} is missing from {}.'.format(colname, catalog_file)))

    # Check which columns are present to see if colors need to be translated
    # into magnitudes
    if 'J' in catalog.colnames and 'H' in catalog.colnames and 'V' in catalog.colnames:
        calculate_magnitudes_from_colors = False
    else:
        if 'V-K' in catalog.colnames and 'J-K' in catalog.colnames and 'J-H' in catalog.colnames:
            calculate_magnitudes_from_colors = True
        else:
            raise ValueError(("ERROR: {} must contain either 'J', 'H', and 'K' columns, or 'V-K', 'J-K' and "
                              "'J-H' columns").format(catalog_file))

    # If the input catalog contains only V magnitudes and colors, calculate
    # the magnitudes in the other Johnson filters
    if calculate_magnitudes_from_colors:
        catalog = catalog_colors_to_vjhk(catalog)

    # Translate Johnson-based filter magnitudes into JWST filter-based
    # magnitudes for the requested filters.
    catalog = transform_johnson_to_jwst(catalog, all_filters)

    # Extract the relevant columns and create a Mirage-formatted point
    # source catalog
    mirage_cat = PointSourceCatalog(ra=catalog[ra_column_name].data,
                                    dec=catalog[dec_column_name].data)

    for filt in all_filters:
        instrument, filter_name, _ = filt.split('_')
        mirage_cat.add_magnitude_column(catalog[filt], instrument=instrument, filter_name=filter_name)

    # Save to output file if requested
    if output_file is not None:
        mirage_cat.save(output_file)

    return mirage_cat


def crossmatch_filter_names(filter_names, standard_filters):
    """Return a list of filter names matching those of standard filters

    Parameters
    ----------
    filter_names : list
        List of input filter names

    standard_filters : list
        List of recognized filter names

    Returns
    -------
    inds : list
        List of matching filter names
    """
    inds = []
    for loop in range(len(filter_names)):
        for n1 in range(len(standard_filters)):
            if filter_names[loop] in standard_filters[n1]:
                inds.append(n1)
    return inds


def match_model_magnitudes(in_magnitudes, in_filters, standard_magnitudes,
                           standard_values, standard_filters, standard_labels):
    """
    This code attempts to make the best match between a set of input magnitudes
    and a set of BOSZ simulated magnitudes.  It is assumed that the input
    magnitudes are not reddened.  The match is for the smallest root-mean-square
    deviation between the input magnitudes and the BOSZ model magnitudes, with
    an overall offset factor applied to the latter.

    Parameters
    ----------
    in_magnitudes : numpy.ndarray
        The list of (A0V) magnitude values to match.

    in_filters : list
        The labels for the magnitudes.

    standard_magntudes : numpy,ndarray
        2D array of standard simulated magnitude values.

    standard_values : numpy.ndarray
        2D array of other filter values (wavelengths,
        zero magnitude flux density values, extinction)

    standard_filters : list
        The list of the standard magnitude filter names

    standard_labels : list
        The labels for the input stellar atmosphere models
        used to calculate the standard magnitudes.

    Returns
    -------
    out_magnitudes : numpy.ndarray or None
        1D array of the full set of estimated magnitudes from the model
        matching, or None if a problem occurs.
    """
    inds = crossmatch_filter_names(in_filters, standard_filters)
    nmatch = float(len(inds))
    if nmatch != len(in_filters):
        print('Error in matching the requested filters for model matching.')
        return None

    subset = np.copy(standard_magnitudes[:, inds])
    del1 = subset - in_magnitudes
    offset = np.mean(del1, axis=1)
    offset_exp = np.expand_dims(offset, axis=1)
    offset_stack = np.repeat(offset_exp, len(in_magnitudes), axis=1)
    delm = (subset - offset_stack - in_magnitudes)
    rms = np.sqrt(np.sum(delm * delm, axis=1) / nmatch)
    smallest = np.where(rms == np.min(rms))[0][0]
    out_magnitudes = standard_magnitudes[smallest, :] - offset[smallest]

    return out_magnitudes


def read_standard_magnitudes():
    """
    The code reads a file magslist_bosz_normal_mirage1.new to get the simulated
    magnitudes.  This file gives simulated magnitudes for the Johnson VJHKL
    filters, the 2MASS filters, the WISE filters, the GAIA filters, and the
    JWST filters.  Also some filter specific values are read from the header
    lines of the file.

    Returns
    -------
    standard_magntudes : numpy.ndarray
        2D array of standard simulated magnitude values.

    standard_values : numpy.ndarray
        2D array of other filter values (wavelengths,
        zero magnitude flux density values, extinction)

    standard_filters : list
        The list of the standard magnitude filter names

    standard_labels : list
        The labels for the input stellar atmosphere models
        used to calculate the standard magnitudes.
    """
    # read in the values needed to transform the Besancon model magnitudes
    #
    module_path = pkg_resources.resource_filename('mirage', '')
    standard_mag_file = os.path.join(module_path, 'config/magslist_bosz_normal_mirage.new')
    with open(standard_mag_file, 'r') as infile:
        lines = infile.readlines()

    standard_magnitudes = np.loadtxt(standard_mag_file, comments='#')
    # standard_values holds the wavelengths (microns), zero magnitude flux
    # density values (W/m^2/micron and Jy) and the relative ISM extinction
    standard_values = np.zeros((4, 58), dtype=np.float32)
    # The following list is manually produced, but must match the order
    # of filters in the input file, where the names are a bit different.
    # Note that for the GAIA g filter the trailing space is needed to
    # allow the code to differentiate the G, BP, and RP filters.
    standard_filters = ['Johnson V', 'Johnson J', 'Johnson H', 'Johnson K',
                        '2MASS J', '2MASS H', '2MASS Ks', 'Johnson L',
                        'WISE W1', 'WISE W2', 'WISE W3', 'WISE W4', 'GAIA g ',
                        'GAIA gbp', 'GAIA grp',
                        'niriss_f090w_magnitude', 'niriss_f115w_magnitude',
                        'niriss_f140m_magnitude', 'niriss_f150w_magnitude',
                        'niriss_f158m_magnitude', 'niriss_f200w_magnitude',
                        'niriss_f277w_magnitude', 'niriss_f356w_magnitude',
                        'niriss_f380m_magnitude', 'niriss_f430m_magnitude',
                        'niriss_f444w_magnitude', 'niriss_f480m_magnitude',
                        'fgs_guider1_magnitude', 'fgs_guider2_magnitude',
                        'nircam_f070w_magnitude', 'nircam_f090w_magnitude',
                        'nircam_f115w_magnitude', 'nircam_f140m_magnitude',
                        'nircam_f150w_magnitude', 'nircam_f150w2_magnitude',
                        'nircam_f162m_magnitude', 'nircam_f164n_magnitude',
                        'nircam_f182m_magnitude', 'nircam_f187n_magnitude',
                        'nircam_f200w_magnitude', 'nircam_f210m_magnitude',
                        'nircam_f212n_magnitude', 'nircam_f250m_magnitude',
                        'nircam_f277w_magnitude', 'nircam_f300m_magnitude',
                        'nircam_f322w2_magnitude', 'nircam_f323n_magnitude',
                        'nircam_f335m_magnitude', 'nircam_f356w_magnitude',
                        'nircam_f360m_magnitude', 'nircam_f405n_magnitude',
                        'nircam_f410m_magnitude', 'nircam_f430m_magnitude',
                        'nircam_f444w_magnitude', 'nircam_f460m_magnitude',
                        'nircam_f466n_magnitude', 'nircam_f470n_magnitude',
                        'nircam_f480m_magnitude']
    standard_labels = []
    n1 = 0
    for line in lines:
        line = line.strip('\n')
        if '#' in line[0:1]:
            values = line.split('#')
            if len(values) == 3:
                v1 = values[-1].split()
                for loop in range(4):
                    standard_values[loop, n1] = float(v1[loop])
                n1 = n1 + 1
        else:
            values = line.split('#')
            standard_labels.append(values[-1])
    return standard_magnitudes, standard_values, standard_filters, standard_labels


def combine_and_interpolate(gaia_cat, gaia_2mass, gaia_2mass_crossref, gaia_wise,
                            gaia_wise_crossref, twomass_cat, wise_cat, instrument, filter_names):
    """
    This function combines GAIA/2MASS/WISE photometry to estimate JWST filter
    magnitudes.  The algorithm depends a bit on what magnitudes are available.

    Note it is implicitly assumed that the GAIA/2MASS/WISE data are for the
    same area of sky.
    The code gives different results depending on what magnitudes are available.
    If only the GAIA g magnitude is available, this is used for all filters.
    If only the GAIA g/BP/RP magnitudes are available, the BP - RP magnitude
    value is used to predict the infrared magnitudes based on standard stellar
    atmosphere simulations.
    If any 2MASS or WISE magnitudes (excluding upper limits) are available the
    JWST magnitudes are found by wavelength interpolation using the pivot
    wavelengths of the different filters including those for GAIA/2MASS/WISE.

    Parameters
    ----------
    gaia_cat : astropy.table.Table
        contains GAIA DR2 search results in table form

    gaia_2mass : astropy.table.Table
        contains 2MASS catalogue values from the GAIA DR2 archive

    gaia_2mass_crossref : astropy.table.Table
        contains GAIA/2MASS cross-references from the GAIA DR2 archive

    gaia_wise : astropy.table.Table
        contains WISE catalogue values from the GAIA DR2 archive
\
    gaia_wise_crossref : astropy.table.Table
        contains GAIA/WISE cross-references from the GAIA DR2 archive

    twomass_cat : astropy.table.Table
        contains 2MASS data from IPAC in table form

    wise_cat : astropy.table.Table
        contains WISE data from IPAC in table form

    instrument : str
        Name of the instrument for which filter magnitudes are
        needed:  "NIRcam", "NIRISS", "Guider", or "All"

    filter_names : list or None
        List of names  of the filters to select for the instrument.
        If the list is empty, or the value is None, all filters are
        selected.

    Returns
    -------
    outcat : mirage.catalogs.create_catalog.PointSourceCatalog
        This is the catalog of positions/magnitudes.
    """
    standard_magnitudes, standard_values, standard_filters, standard_labels = read_standard_magnitudes()
    nfilters = len(filter_names)
    ngaia = len(gaia_cat['ra'])
    n2mass1 = len(gaia_2mass['ra'])
    nwise1 = len(gaia_wise_crossref['ra'])
    n2mass2 = len(twomass_cat['ra'])
    nwise2 = len(wise_cat['ra'])
    nout = ngaia + n2mass2 + nwise2
    in_magnitudes = np.zeros((nout, 10), dtype=np.float32) + 10000.0
    raout = np.zeros((nout), dtype=np.float32)
    decout = np.zeros((nout), dtype=np.float32)
    # magnitudes Gaia bp, g, rp; 2MASS J, H, Ks; WISE W1, W2, W3, W4
    in_filters = ['GAIA gbp', 'GAIA g ', 'GAIA grp', '2MASS J', '2MASS H',
                  '2MASS Ks', 'WISE W1', 'WISE W2', 'WISE W3', 'WISE W4']
    inds = crossmatch_filter_names(in_filters, standard_filters)
    if len(inds) != len(in_filters):
        print('Error matching the filters to the standard set.')
        return None
    # first populate the gaia sources, with cross-references
    in_magnitudes[0:ngaia, 1] = gaia_cat['phot_g_mean_mag']
    raout[0:ngaia] = gaia_cat['ra']
    decout[0:ngaia] = gaia_cat['dec']
    ngaia2masscr = twomass_crossmatch(gaia_cat, gaia_2mass, gaia_2mass_crossref, twomass_cat)
    twomassflag = [True] * n2mass2
    for n1 in range(n2mass2):
        for loop in range(len(gaia_2mass['ra'])):
            if gaia_2mass['designation'][loop] == twomass_cat['designation'][n1]:
                if ngaia2masscr[n1] >= 0:
                    twomassflag[n1] = False
    matchwise, gaiawiseinds, twomasswiseinds = wise_crossmatch(gaia_cat, gaia_wise, gaia_wise_crossref, wise_cat, twomass_cat)
    wisekeys = ['w1sigmpro', 'w2sigmpro', 'w3sigmpro', 'w4sigmpro']

    # Set invalid values to NaN
    try:
        gaia_bp_mags = gaia_cat['phot_bp_mean_mag'].filled(np.nan)
        gaia_rp_mags = gaia_cat['phot_rp_mean_mag'].filled(np.nan)
    except:
        gaia_bp_mags = gaia_cat['phot_bp_mean_mag']
        gaia_rp_mags = gaia_cat['phot_rp_mean_mag']

    for loop in range(ngaia):
        try:
            in_magnitudes[loop, 0] = gaia_bp_mags[loop]
            in_magnitudes[loop, 2] = gaia_rp_mags[loop]
        except:
            pass
        # see if there is a 2MASS match
        for n1 in range(n2mass1):
            if loop == ngaia2masscr[n1]:
                in_magnitudes[loop, 3] = gaia_2mass['j_m'][n1]
                in_magnitudes[loop, 4] = gaia_2mass['h_m'][n1]
                in_magnitudes[loop, 5] = gaia_2mass['ks_m'][n1]
                for l1 in range(3):
                    if gaia_2mass['ph_qual'][n1][l1] == 'U':
                        in_magnitudes[loop, 3+l1] = 10000.
        # see if there is a WISE match
        for n2 in range(len(matchwise)):
            if matchwise[n2] and (gaiawiseinds[n2] == loop):
                in_magnitudes[loop, 6] = wise_cat['w1mpro'][n2]
                in_magnitudes[loop, 7] = wise_cat['w2mpro'][n2]
                in_magnitudes[loop, 8] = wise_cat['w3mpro'][n2]
                in_magnitudes[loop, 9] = wise_cat['w4mpro'][n2]
                for l1 in range(4):
                    if not isinstance(wise_cat[wisekeys[l1]][n2], float):
                        in_magnitudes[loop, 6+l1] = 10000.
    # Add in any 2MASS sources with no GAIA match
    n1 = 0
    noff = ngaia
    for loop in range(n2mass2):
        if twomassflag[loop]:
            raout[noff+n1] = twomass_cat['ra'][loop]
            decout[noff+n1] = twomass_cat['dec'][loop]
            in_magnitudes[noff+n1, 3] = twomass_cat['j_m'][loop]
            in_magnitudes[noff+n1, 4] = twomass_cat['h_m'][loop]
            in_magnitudes[noff+n1, 5] = twomass_cat['k_m'][loop]
            for l1 in range(3):
                if twomass_cat['ph_qual'][loop][l1] == 'U':
                    in_magnitudes[noff+n1, 3+l1] = 10000.
            # Check to see if there is a WISE cross-match
            for l1 in range(len(twomasswiseinds)):
                if twomasswiseinds[l1] == loop:
                    in_magnitudes[noff+n1, 6] = wise_cat['w1mpro'][l1]
                    in_magnitudes[noff+n1, 7] = wise_cat['w2mpro'][l1]
                    in_magnitudes[noff+n1, 8] = wise_cat['w3mpro'][l1]
                    in_magnitudes[noff+n1, 9] = wise_cat['w4mpro'][l1]
                    for l2 in range(4):
                        if not isinstance(wise_cat[wisekeys[l2]][l1], float):
                            in_magnitudes[noff+n1, 6+l2] = 10000.
            n1 = n1 + 1
    # Finally, add in WISE sources that have not been cross-matched to GAIA
    # or 2MASS.
    noff = ngaia+n1
    n1 = 0
    for loop in range(nwise2):
        if (not matchwise[loop]) and (twomasswiseinds[loop] < 0):
            raout[noff+n1] = wise_cat['ra'][loop]
            decout[noff+n1] = wise_cat['dec'][loop]
            in_magnitudes[noff+n1, 6] = wise_cat['w1mpro'][loop]
            in_magnitudes[noff+n1, 7] = wise_cat['w2mpro'][loop]
            in_magnitudes[noff+n1, 8] = wise_cat['w3mpro'][loop]
            in_magnitudes[noff+n1, 9] = wise_cat['w3mpro'][loop]
            for l1 in range(4):
                if not isinstance(wise_cat[wisekeys[l1]][loop], float):
                    in_magnitudes[noff+n1, 6+l1] = 10000.
            n1 = n1+1
    # Now, convert to JWST magnitudes either by transformation (for sources
    # with GAIA G/BP/RP magnitudes) or by interpolation (all other
    # cases).
    out_filter_names = make_filter_names(instrument, filter_names)
    inds = crossmatch_filter_names(in_filters, standard_filters)
    in_wavelengths = np.squeeze(np.copy(standard_values[0, inds]))
    inds = crossmatch_filter_names(out_filter_names, standard_filters)
    if len(inds) < 1:
        return None
    out_wavelengths = np.squeeze(np.copy(standard_values[0, inds]))
    if len(inds) == 1:
        out_wavelengths = np.zeros((1), dtype=np.float32)+out_wavelengths
    nfinal = noff + n1

    out_magnitudes = np.zeros((nout, len(out_filter_names)),
                              dtype=np.float32)
    for loop in range(nfinal):
        values = interpolate_magnitudes(in_wavelengths, in_magnitudes[loop, :],
                                        out_wavelengths, out_filter_names)
        out_magnitudes[loop, :] = np.copy(values)
    raout = np.copy(raout[0:nfinal])
    decout = np.copy(decout[0:nfinal])
    out_magnitudes = np.copy(out_magnitudes[0:nfinal, :])
    outcat = PointSourceCatalog(ra=raout, dec=decout)
    n1 = 0
    for filter in out_filter_names:
        values = filter.split('_')
        if len(values) == 3:
            inst = values[0].upper()
            fil1 = values[1].upper()
        else:
            inst = values[0].upper()
            fil1 = ''
        outcat.add_magnitude_column(np.squeeze(out_magnitudes[:, n1]), instrument=inst,
                                    filter_name=fil1, magnitude_system='vegamag')
        n1 = n1+1
    return outcat


def twomass_crossmatch(gaia_cat, gaia_2mass, gaia_2mass_crossref, twomass_cat):
    """
    Take the GAIA to 2MASS cross references and make sure that there is only
    one GAIA source cross-matched to a given 2MASS source in the table.

    Parameters
    ----------
    gaia_cat : astropy.table.Table
        contains the GAIA DR2 catalogue values from the GAIA archive

    gaia_2mass : astropy.table.Table
        contains 2MASS catalogue values from the GAIA DR2 archive

    gaia_2mass_crossref : astropy.table.Table
        contains GAIA/2MASS cross-references from the GAIA DR2 archive

    twomass_cat : astropy.table.Table
        contains the 2MASS catalogue values from IPAC

    Returns
    -------
    ngaia2masscr : numpy.ndarray
                    an integer array of cross match indexes, giving for
                    each 2MASS source the index number of the associated
                    GAIA source in the main GAIA table, or a value of -10
                    where there is no match
    """
    ntable1 = len(gaia_cat['designation'])
    ntable2 = len(gaia_2mass['ra'])
    ntable3 = len(gaia_2mass_crossref['designation'])
    ntable4 = len(twomass_cat['ra'])
    ngaia2mass = np.zeros((ntable2), dtype=np.int16) - 10
    ngaia2masscr = np.zeros((ntable4), dtype=np.int16) - 10
    for loop in range(ntable2):
        # find the number of entries of each 2MASS source in the cross references
        nmatch = 0
        namematch = []
        match1 = []
        for l1 in range(ntable3):
            if gaia_2mass['designation'][loop] == gaia_2mass_crossref['designation_2'][l1]:
                nmatch = nmatch + 1
                namematch.append(gaia_2mass_crossref['designation'][l1])
                match1.append(loop)
        # Find the matching GAIA sources and select the one with the best
        # magnitude match within a radius of 0.3 arc-seconds of the 2MASS
        # position.
        magkeys = ['j_m', 'h_m', 'ks_m']
        if nmatch > 0:
            mindelm = 10000.0
            ncross = -10
            for l1 in range(nmatch):
                for l2 in range(ntable1):
                    gmag = 0.
                    irmag = -10000.0
                    if gaia_cat['designation'][l2] == namematch[l1]:
                        ra1 = gaia_cat['ra'][l2]
                        dec1 = gaia_cat['dec'][l2]
                        ra2 = gaia_2mass['ra'][match1[l1]]
                        dec2 = gaia_2mass['dec'][match1[l1]]
                        p1 = SkyCoord(ra1*u.deg, dec1*u.deg)
                        p2 = SkyCoord(ra2*u.deg, dec2*u.deg)
                        if p2.separation(p1).arcsec < 0.3:
                            gmag = gaia_cat['phot_g_mean_mag'][l2]
                            # select 2MASS magnitude: first ph_qual = A or if none
                            # is of quality A the first ph_qual = B or if none is
                            # of quality A or B then the first non U value.
                            for l3 in range(3):
                                if (irmag < -100.) and (gaia_2mass['ph_qual'][match1[l1]].decode()[l3:l3+1] == "A"):
                                    irmag = gaia_2mass[magkeys[l3]][match1[l1]]
                            for l3 in range(3):
                                if (irmag < -100.) and (gaia_2mass['ph_qual'][match1[l1]].decode()[l3:l3+1] == "B"):
                                    irmag = gaia_2mass[magkeys[l3]][match1[l1]]
                            for l3 in range(3):
                                if (irmag < -100.) and (gaia_2mass['ph_qual'][match1[l1]].decode()[l3:l3+1] != "U"):
                                    irmag = gaia_2mass[magkeys[l3]][match1[l1]]
                            delm = gmag - irmag
                            if (delm > -1.2) and (delm < 30.0):
                                if delm < mindelm:
                                    ncross = l2
                                    mindelm = delm
                ngaia2mass[loop] = ncross
    # Now locate the 2MASS sources in the IPAC 2MASS table, and put in the
    # index values.
    for loop in range(ntable4):
        for n1 in range(ntable2):
            if twomass_cat['designation'][loop] == gaia_2mass['designation'][n1]:
                ngaia2masscr[loop] = ngaia2mass[n1]
    return ngaia2masscr


def wise_crossmatch(gaia_cat, gaia_wise, gaia_wise_crossref, wise_cat, twomass_cat):
    """
    Relate the GAIA/WISE designations to the WISE catalogue designations, since the names
    change a little between the different catalogues.  Return the boolean list of matches and
    the index values in wise_cat.

    Parameters
    ----------
    gaia_cat : astropy.table.Table
        contains the GAIA DR2 catalogue values from the GAIA archive

    gaia_wise : astropy.table.Table
        contains WISE catalogue values from the GAIA DR2 archive

    gaia_wise_crossref : astropy.table.Table
        contains GAIA/WISE cross-references from the GAIA DR2 archive

    wise_cat : astropy.table.Table
        contains WISE data from IPAC in table form

    twomass_cat : astropy.table.Table
        contains 2MASS data from IPAC in table form

    gaia2masscr : list
        list of integers of the cross-reference indexes from 2MASS to GAIA

    Returns
    -------
    matchwise : list
        boolean list of length equal to wise_cat with True if there is a
        cross-match with GAIA

    gaiawiseinds :  list
        list of integer index values from wise_cat to gaia_cat (i.e. the
        GAIA number to which the WISE source corresponds)

    twomasswiseinds : list
        list of integer index values from wise_cat to twomass_cat (i.e.
        the 2MASS number to which the WISE source corresponds)
    """
    num_entries = len(wise_cat['ra'])
    num_gaia = len(gaia_cat['ra'])
    matchwise = [False] * num_entries
    gaiawiseinds = [-1] * num_entries
    twomasswiseinds = [-1] * num_entries
    ra1 = np.copy(wise_cat['ra'])
    dec1 = np.copy(wise_cat['dec'])
    ra2 = np.copy(gaia_wise['ra'])
    dec2 = np.copy(gaia_wise['dec'])
    ra3 = np.copy(gaia_wise_crossref['ra'])
    dec3 = np.copy(gaia_wise_crossref['dec'])
    sc1 = SkyCoord(ra=ra1*u.degree, dec=dec1*u.degree)
    sc2 = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)
    sc3 = SkyCoord(ra=ra3*u.degree, dec=dec3*u.degree)
    n2mass = len(twomass_cat['ra'])
    # look at the WISE data and find the sources with listed 2MASS counterparts
    for loop in range(num_entries):
        if not np.isnan(wise_cat['h_m_2mass'][loop]):
            for n1 in range(n2mass):
                if (abs(twomass_cat['j_m'][n1] - wise_cat['j_m_2mass'][loop]) < 0.001) and \
                   (abs(twomass_cat['h_m'][n1] - wise_cat['h_m_2mass'][loop]) < 0.001) and \
                   (abs(twomass_cat['k_m'][n1] - wise_cat['k_m_2mass'][loop]) < 0.001):
                    twomasswiseinds[loop] = n1
                    break
    # match WISE to gaia_wise by position
    idx, d2d, d3d = sc3.match_to_catalog_sky(sc1)
    for loop in range(len(idx)):
        if (d2d[loop].arcsec) < 0.4:
            matchwise[idx[loop]] = True
            for n2 in range(num_gaia):
                if gaia_cat['designation'][n2] == gaia_wise_crossref['designation'][loop]:
                    gaiawiseinds[idx[loop]] = n2
                    break
    return matchwise, gaiawiseinds, twomasswiseinds


def interpolate_magnitudes(wl1, mag1, wl2, filternames):
    """
    Given an input set of magnitudes and associated wavelengths, interpolate
    these in wavelength to get approximate JWST magnitudes.  It is assumed that
    the first 3 magnitudes in the input vector mag1 are the GAIA BP, g, and
    RP filters and that the rest are long wavelength magnitudes.  If only the
    GAIA g filter magnitude is available it is used for all output magnitudes.
    If only GAIA BP/g/RP magnitudes are available the GAIA BP/RP magnitudes
    are transformed to JWST magnitudes based on stellar model values.  In the
    case where there is some infrared data, the numpy interp function is used
    to interpolate the magnitudes.  In the inputs filters without data are
    assigned magnitude > 100, and these are not used in the interpolation.

    As is normal for numpy interp, wavelength values outside the defined
    range in the input magnitudes get nearest magnitude value.

    The input magnitudes must be A0V magnitudes not AB or ST magnitudes.

    Parameters
    ----------
    wl1 : numpy.ndarray
        The pivot wavelengths, in microns, for the input filters.  The
        values need to be sorted before passing to the routine.

    mag1 : numpy.ndarray
        The associated magnitudes (A0V by assumption).  Values > 100.
        indicate "no data".

    wl2 : numpy.ndarray
        The pivot wavelengths, in microns, for the output filters

    filternames : list
        The names of the output filters, used when the GAIA blue/red
        magnitudes are available but no near-infrared magnitudes
        are available.

    Returns
    -------
    out_magnitudes : numpy.ndarray
        The output interpolated magnitudes corresponding to the
        wavelengths.
    """
    inds = np.isnan(mag1)
    mag1[inds] = 10000.
    nout = len(wl2)
    outmags = wl2*0.+10000.
    # Case 1:  All dummy values, return all magnitudes = 10000.0 (this should
    #          not happen)
    if np.min(mag1) > 100.:
        return outmags
    # Case 2,  Only GAIA magnitudes.  Either transform from the GAIA BP and RP
    #          colour to the JWST magnitudes or assume a default colour value
    #          for the star, equivalent to a star of type K4V.
    if np.min(mag1[3:]) > 100.:
        if (mag1[0] > 100.) or (mag1[2] > 100.):
            # Where the BP and RP magnitudes are not available, make colours
            # matching a K4V star (assumed T=4500, log(g)=5.0)
            inmags = np.zeros((2), dtype=np.float32)
            inmags[0] = mag1[1] + 0.5923
            inmags[1] = mag1[1] - 0.7217
        else:
            inmags = np.copy(mag1[[0, 2]])
        standard_magnitudes, standard_values, standard_filters, standard_labels = read_standard_magnitudes()
        in_filters = ['GAIA gbp', 'GAIA grp']
        newmags = match_model_magnitudes(inmags, in_filters, standard_magnitudes, standard_values,
                                         standard_filters, standard_labels)
        inds = crossmatch_filter_names(filternames, standard_filters)
        outmags = np.copy(newmags[inds])
        return outmags
    # Case 3, some infrared magnitudes are available, interpolate good values
    # (magnitude = 10000 for bad values)
    inds = np.where(mag1 < 100.)
    inmags = mag1[inds]
    inwl = wl1[inds]
    outmags = np.interp(wl2, inwl, inmags)
    return outmags


def make_filter_names(instrument, filters):
    """
    Given as input the instrument name and the list of filters needed, this
    routine generates the list of output header values for the Mirage input
    file.

    Parameters
    ----------
    insrument : str
        'NIRCam', 'NIRISS', or 'Guider"

    filters : list
        List of the filters needed (names such as F090W) or either None or
        an empty list to select all filters)

    Returns
    -------
    headerstrs : list
        The header string list for the output Mirage magnitudes line;
        if concatentated with spaces, it is the list of magnitudes header
        values that needs to be written to the source list file.
    """
    instrument_names = ['FGS', 'NIRCam', 'NIRISS']
    guider_filters = ['guider1', 'guider2']
    guider_filter_names = ['fgs_guider1_magnitude', 'fgs_guider2_magnitude']
    niriss_filters = ['F090W', 'F115W', 'F140M', 'F150W', 'F158M', 'F200W',
                      'F277W', 'F356W', 'F380M', 'F430M', 'F444W', 'F480M']
    niriss_filter_names = ['niriss_{}_magnitude'.format(filt.lower()) for filt in niriss_filters]
    nircam_filters = ['F070W', 'F090W', 'F115W', 'F140M', 'F150W', 'F150W2',
                      'F162M', 'F164N', 'F182M', 'F187N', 'F200W', 'F210M',
                      'F212N', 'WLP4',  'F250M', 'F277W', 'F300M', 'F322W2', 'F323N',
                      'F335M', 'F356W', 'F360M', 'F405N', 'F410M', 'F430M',
                      'F444W', 'F460M', 'F466N', 'F470N', 'F480M']
    nircam_filter_names = ['nircam_{}_magnitude'.format(filt.lower()) for filt in nircam_filters]
    names1 = [guider_filters, nircam_filters, niriss_filters]
    names2 = [guider_filter_names, nircam_filter_names, niriss_filter_names]
    headerstrs = []
    for loop in range(len(instrument_names)):
        if (instrument_names[loop].lower() == instrument.lower()) or (instrument.lower() == 'all'):
            headerstrs = add_filter_names(headerstrs, names1[loop], names2[loop], filters)
    return headerstrs


def add_filter_names(headerlist, filter_names, filter_labels, filters):
    """
    Add a set of filter header labels (i.e. niriss_f090w_magnitude for example)
    to a list, by matching filter names.

    Parameters
    ----------
    headerlist : list
        An existing (possibly empty) list to hold the header string for the
        output magnitudes

    filter_names : list
        The list of available filter names to match to

    filter_labels : list
        The corresponding list of filter labels

    filters : list
        The list of filter names to match, or an empty list or None to get
        all available filter labels

    Returns
    -------
    headerlist : list
        The revised list of labels with the filter labels requested appended
    """
    try:
        n1 = len(filters)
    except:
        n1 = 0
    if (filters is None) or (n1 == 0):
        for loop in range(len(filter_labels)):
            headerlist.append(filter_labels[loop])
    if n1 > 0:
        for loop in range(n1):
            for k in range(len(filter_names)):
                if filters[loop].lower() == filter_names[k].lower():
                    headerlist.append(filter_labels[k])
    return headerlist


def combine_catalogs(cat1, cat2, magnitude_fill_value=99.):
    """Combine two Mirage catalog objects. Catalogs must be of the same
    type (e.g. PointSourceCatalog), and have the same values for position
    units (RA, Dec or x, y) velocity units (arcsec/hour vs pixels/hour),
    and galaxy radius units (arcsec vs pixels), if present. This uses the
    astropy.table join method. A source common to both tables will be
    identified and combined into a single row. Sources that are not common
    to both input catalogs will have a placeholder value for magnitude
    columns where there is no data.

    Parameters
    ----------
    cat1 : mirage.catalogs.catalog_generator.XXCatalog
        First catalog to be joined

    cat2 : mirage.catalogs.catalog_generator.XXCatalog
        Second catalog to be joined

    Returns
    -------
    new_cat : mirage.catalogs.catalog_generator.XXCatalog
        Combined catalog
    """
    # Be sure that location coordinate systems are the same between the
    # two catalogs
    # Make sure the catalogs are the same type
    if type(cat1) != type(cat2):
        raise TypeError("Catalogs are different types. Cannot be combined.")

    if cat1.location_units != cat2.location_units:
        raise ValueError('Coordinate mismatch in catalogs to combine.')

    # Join catalog tables. Set fill value for all magnitude columns
    combined = join(cat1.table, cat2.table, join_type='outer')
    mag_cols = [col for col in combined.colnames if 'magnitude' in col]
    for col in mag_cols:
        combined[col].fill_value = magnitude_fill_value
    combined = combined.filled()

    # NOTE that the order of this if/elif statement matters, as different
    # catalog classes inherit from each other.

    # --------------Moving Galaxies---------------------------------------
    if isinstance(cat1, MovingSersicCatalog):
        if cat1.velocity_units != cat2.velocity_units:
            raise ValueError('Velocity unit mismatch in catalogs to combine.')
        if cat1.radius_units != cat2.radius_units:
            raise ValueError('Radius unit mismatch in catalogs to combine.')

        if cat1.location_units == 'position_RA_Dec':
            if cat1.velocity_units == 'velocity_RA_Dec':
                new_cat = MovingSersicCatalog(ra=combined['x_or_RA'].data,
                                              dec=combined['y_or_Dec'].data,
                                              ra_velocity=combined['ra_velocity'].data,
                                              dec_velocity=combined['dec_velocity'].data,
                                              ellipticity=combined['ellipticity'].data,
                                              radius=combined['radius'].data,
                                              sersic_index=combined['sersic_index'].data,
                                              position_angle=combined['pos_angle'].data,
                                              radius_units=cat1.radius_units)
            else:
                new_cat = MovingSersicCatalog(ra=combined['x_or_RA'].data,
                                              dec=combined['y_or_Dec'].data,
                                              x_velocity=combined['ra_velocity'].data,
                                              y_velocity=combined['dec_velocity'].data,
                                              ellipticity=combined['ellipticity'].data,
                                              radius=combined['radius'].data,
                                              sersic_index=combined['sersic_index'].data,
                                              position_angle=combined['pos_angle'].data,
                                              radius_units=cat1.radius_units)
        else:
            if cat1.velocity_units == 'velocity_RA_Dec':
                new_cat = MovingSersicCatalog(x=combined['x_or_RA'].data,
                                              y=combined['y_or_Dec'].data,
                                              ra_velocity=combined['ra_velocity'].data,
                                              dec_velocity=combined['dec_velocity'].data,
                                              ellipticity=combined['ellipticity'].data,
                                              radius=combined['radius'].data,
                                              sersic_index=combined['sersic_index'].data,
                                              position_angle=combined['pos_angle'].data,
                                              radius_units=cat1.radius_units)
            else:
                new_cat = MovingSersicCatalog(x=combined['x_or_RA'].data,
                                              y=combined['y_or_Dec'].data,
                                              x_velocity=combined['ra_velocity'].data,
                                              y_velocity=combined['dec_velocity'].data,
                                              ellipticity=combined['ellipticity'].data,
                                              radius=combined['radius'].data,
                                              sersic_index=combined['sersic_index'].data,
                                              position_angle=combined['pos_angle'].data,
                                              radius_units=cat1.radius_units)

    # --------------Moving Extended Sources-------------------------------
    elif isinstance(cat1, MovingExtendedCatalog):
        if cat1.velocity_units != cat2.velocity_units:
            raise ValueError('Velocity unit mismatch in catalogs to combine.')

        if cat1.location_units == 'position_RA_Dec':
            if cat1.velocity_units == 'velocity_RA_Dec':
                new_cat = MovingExtendedCatalog(ra=combined['x_or_RA'].data,
                                                dec=combined['y_or_Dec'].data,
                                                ra_velocity=combined['ra_velocity'].data,
                                                dec_velocity=combined['dec_velocity'].data,
                                                filenames=combined['filename'].data,
                                                position_angle=combined['pos_angle'].data)
            else:
                new_cat = MovingExtendedCatalog(ra=combined['x_or_RA'].data,
                                                dec=combined['y_or_Dec'].data,
                                                x_velocity=combined['ra_velocity'].data,
                                                y_velocity=combined['dec_velocity'].data,
                                                filenames=combined['filename'].data,
                                                position_angle=combined['pos_angle'].data)
        else:
            if cat1.velocity_units == 'velocity_RA_Dec':
                new_cat = MovingExtendedCatalog(x=combined['x_or_RA'].data,
                                                y=combined['y_or_Dec'].data,
                                                ra_velocity=combined['ra_velocity'].data,
                                                dec_velocity=combined['dec_velocity'].data,
                                                filenames=combined['filename'].data,
                                                position_angle=combined['pos_angle'].data)
            else:
                new_cat = MovingExtendedCatalog(x=combined['x_or_RA'].data,
                                                y=combined['y_or_Dec'].data,
                                                x_velocity=combined['ra_velocity'].data,
                                                y_velocity=combined['dec_velocity'].data,
                                                filenames=combined['filename'].data,
                                                position_angle=combined['pos_angle'].data)

    # --------------------Galaxies------------------------------------
    elif isinstance(cat1, GalaxyCatalog):
        if cat1.radius_units != cat2.radius_units:
            raise ValueError('Radius unit mismatch in catalogs to combine.')

        if cat1.location_units == 'position_RA_Dec':
            new_cat = GalaxyCatalog(ra=combined['x_or_RA'].data,
                                    dec=combined['y_or_Dec'].data,
                                    ellipticity=combined['ellipticity'].data,
                                    radius=combined['radius'].data,
                                    sersic_index=combined['sersic_index'].data,
                                    position_angle=combined['pos_angle'].data,
                                    radius_units=cat1.radius_units)
        else:
            new_cat = GalaxyCatalog(x=combined['x_or_RA'].data,
                                    y=combined['y_or_Dec'].data,
                                    ellipticity=combined['ellipticity'].data,
                                    radius=combined['radius'].data,
                                    sersic_index=combined['sersic_index'].data,
                                    position_angle=combined['pos_angle'].data,
                                    radius_units=cat1.radius_units)

    # ------------------Extended Sources-------------------------------
    elif isinstance(cat1, ExtendedCatalog):
        if cat1.location_units == 'position_RA_Dec':
            new_cat = ExtendedCatalog(ra=combined['x_or_RA'].data,
                                      dec=combined['y_or_Dec'].data,
                                      filenames=combined['filename'].data,
                                      position_angle=combined['pos_angle'].data)
        else:
            new_cat = ExtendedCatalog(x=combined['x_or_RA'].data,
                                      y=combined['y_or_Dec'].data,
                                      filenames=combined['filename'].data,
                                      position_angle=combined['pos_angle'].data)

    # -------------Moving Point Sources--------------------------------
    elif isinstance(cat1, MovingPointSourceCatalog):
        if cat1.velocity_units != cat2.velocity_units:
            raise ValueError('Velocity unit mismatch in catalogs to combine.')

        if cat1.location_units == 'position_RA_Dec':
            if cat1.velocity_units == 'velocity_RA_Dec':
                new_cat = MovingPointSourceCatalog(ra=combined['x_or_RA'].data,
                                                   dec=combined['y_or_Dec'].data,
                                                   ra_velocity=combined['ra_velocity'].data,
                                                   dec_velocity=combined['dec_velocity'].data)
            else:
                new_cat = MovingPointSourceCatalog(ra=combined['x_or_RA'].data,
                                                   dec=combined['y_or_Dec'].data,
                                                   x_velocity=combined['ra_velocity'].data,
                                                   y_velocity=combined['dec_velocity'].data)
        else:
            if cat1.location_units == 'velocity_RA_Dec':
                new_cat = MovingPointSourceCatalog(x=combined['x_or_RA'].data,
                                                   y=combined['y_or_Dec'].data,
                                                   ra_velocity=combined['ra_velocity'].data,
                                                   dec_velocity=combined['dec_velocity'].data)
            else:
                new_cat = MovingPointSourceCatalog(x=combined['x_or_RA'].data,
                                                   y=combined['y_or_Dec'].data,
                                                   x_velocity=combined['ra_velocity'].data,
                                                   y_velocity=combined['dec_velocity'].data)

    # --------------------Point Sources-------------------------------
    elif isinstance(cat1, PointSourceCatalog):
        # Create new catalog object and populate
        if cat1.location_units == 'position_RA_Dec':
            new_cat = PointSourceCatalog(ra=combined['x_or_RA'].data,
                                         dec=combined['y_or_Dec'].data)
        else:
            new_cat = PointSourceCatalog(x=combined['x_or_RA'].data,
                                         y=combined['y_or_Dec'].data)


    # -------------Add magnitude columns-------------------------------
    mag_cols = [colname for colname in combined.colnames if 'magnitude' in colname]
    for col in mag_cols:
        elements = col.split('_')
        instrument = elements[0]
        minus_inst = col.split('{}_'.format(instrument))[-1]
        filter_name = minus_inst.split('_magnitude')[0]
        new_cat.add_magnitude_column(combined[col].data, instrument=instrument,
                                     filter_name=filter_name)

    return new_cat


def combine_catalogs_v0(observed_jwst, besancon_jwst):
    """
    Replaced by combine_catalogs. Do we need to keep this?
    This code takes two input PointSourceCatalog objects and returns the
    combined PointSourceCatalog.  The two catalogs have to be in the same
    magnitude units and have the same set of filters.
    Input values:
    observed_jwst:    (mirage.catalogs.catalog_generator.PointSourceCatalog)
                      Catalog object one
    besancon_jwst:    (mirage.catalogs.catalog_generator.PointSourceCatalog)
                      Catalog object two
    Return value:
    outcat:           (mirage.catalogs.catalog_generator.PointSourceCatalog)
                      A new catalog object combining the two input catalogs
    """
    keys1 = list(observed_jwst.magnitudes.keys())
    keys2 = list(besancon_jwst.magnitudes.keys())
    besanconinds = []
    for key in keys1:
        for loop in range(len(keys2)):
            if key == keys2[loop]:
                besanconinds.append(loop)
    if len(keys1) != len(besanconinds):
        print('Magnitude mismatch in catalogs to combine.  Will return None.')
        return None
    if observed_jwst.location_units != besancon_jwst.location_units:
        print('Coordinate mismatch in catalogs to combine.  Will return None.')
        return None
    ra1 = observed_jwst.ra
    dec1 = observed_jwst.dec
    ra2 = besancon_jwst.ra
    dec2 = besancon_jwst.dec
    raout = np.concatenate((ra1, ra2))
    decout = np.concatenate((dec1, dec2))

    outcat = PointSourceCatalog(ra=raout, dec=decout)
    #outcat.location_units = observed_jwst.location_units
    for key in keys1:
        mag1 = observed_jwst.magnitudes[key][1]
        mag2 = besancon_jwst.magnitudes[key][1]
        magout = np.concatenate((mag1, mag2))
        values = key.split('_')
        instrument = values[0]
        filter = values[1]
        outcat.add_magnitude_column(magout, magnitude_system=observed_jwst.magnitudes[key][0],
                                    instrument=instrument, filter_name=filter)
    return outcat


def get_gaia_ptsrc_catalog(ra, dec, box_width):
    """Wrapper around Gaia query and creation of mirage-formatted catalog"""
    gaia_cat, gaia_mag_cols, gaia_2mass, gaia_2mass_crossref, gaia_wise, \
        gaia_wise_crossref = query_GAIA_ptsrc_catalog(ra, dec, box_width)
    gaia_mirage = mirage_ptsrc_catalog_from_table(gaia_cat, 'gaia', gaia_mag_cols)
    return gaia_mirage, gaia_cat, gaia_2mass_crossref, gaia_wise_crossref


def query_GAIA_ptsrc_catalog(ra, dec, box_width):
    """
    This code is adapted from gaia_crossreference.py by Johannes Sahlmann.  It
    queries the GAIA DR2 archive for sources within a given square region of
    the sky and rerurns the catalogue along withe the 2MASS and WISE
    cross-references for use in combining the catalogues to get the infrared
    magnitudes for the sources that are detected in the other telescopes.

    The GAIA main table has all the GAIA values, but the other tables have only
    specific subsets of the data values to save space.  In the "crossref"
    tables the 'designation' is the GAIA name while 'designation_2' is the
    2MASS or WISE name.

    Parameters
    ----------
    ra : float
        Right ascension of the target field in degrees

    dec : float
        Declination of the target field in degrees

    box_width : float
        Width of the (square) sky area, in arc-seconds

    Returns
    -------
    gaia_cat : astropy.table.Table
        The gaia DR2 magnitudes and other data

    gaia_mag_cols : list
        List of the GAIA magnitude column names

    gaia_2mass : astropy.table.Table
        The 2MASS values as returned from the GAIA archive

    gaia_2mass_crossref : astropy.table.Table
        The cross-reference list with 2MASS sources

    gaia_wise : astropy.table.Table
        The WISE values as returned from the GAIA archive

    gaia_wise_crossref : astropy.table.Table
        The cross-reference list with WISE sources
    """
    data = OrderedDict()
    data['gaia'] = OrderedDict()
    data['tmass'] = OrderedDict()
    data['wise'] = OrderedDict()
    data['tmass_crossmatch'] = OrderedDict()
    data['wise_crossmatch'] = OrderedDict()
    # convert box width to degrees for the GAIA query
    boxwidth = box_width/3600.
    data['gaia']['query'] = """SELECT * FROM gaiadr2.gaia_source AS gaia
                        WHERE 1=CONTAINS(POINT('ICRS',gaia.ra,gaia.dec), BOX('ICRS',{}, {}, {}, {}))
                        """.format(ra, dec, boxwidth, boxwidth)

    data['tmass']['query'] = """SELECT ra,dec,ph_qual,j_m,h_m,ks_m,designation FROM gaiadr1.tmass_original_valid AS tmass
                        WHERE 1=CONTAINS(POINT('ICRS',tmass.ra,tmass.dec), BOX('ICRS',{}, {}, {}, {}))
                        """.format(ra, dec, boxwidth, boxwidth)

    data['tmass_crossmatch']['query'] = """SELECT field.ra,field.dec,field.designation,tmass.designation from
            (SELECT gaia.*
            FROM gaiadr2.gaia_source AS gaia
            WHERE 1=CONTAINS(POINT('ICRS',gaia.ra,gaia.dec), BOX('ICRS',{}, {}, {}, {})))
            AS field
            INNER JOIN gaiadr2.tmass_best_neighbour AS xmatch
                ON field.source_id = xmatch.source_id
            INNER JOIN gaiadr1.tmass_original_valid AS tmass
                ON tmass.tmass_oid = xmatch.tmass_oid
        """.format(ra, dec, boxwidth, boxwidth)

    data['wise']['query'] = """SELECT ra,dec,ph_qual,w1mpro,w2mpro,w3mpro,w4mpro,designation FROM gaiadr1.allwise_original_valid AS wise
                        WHERE 1=CONTAINS(POINT('ICRS',wise.ra,wise.dec), BOX('ICRS',{}, {}, {}, {}))
                        """.format(ra, dec, boxwidth, boxwidth)

    data['wise_crossmatch']['query'] = """SELECT field.ra,field.dec,field.designation,allwise.designation from
            (SELECT gaia.*
            FROM gaiadr2.gaia_source AS gaia
            WHERE 1=CONTAINS(POINT('ICRS',gaia.ra,gaia.dec), BOX('ICRS',{}, {}, {}, {})))
            AS field
            INNER JOIN gaiadr2.allwise_best_neighbour AS xmatch
                ON field.source_id = xmatch.source_id
            INNER JOIN gaiadr1.allwise_original_valid AS allwise
                ON allwise.designation = xmatch.original_ext_source_id
        """.format(ra, dec, boxwidth, boxwidth)

    outvalues = {}
    print('Searching the GAIA DR2 catalog')
    for key in data.keys():
        job = Gaia.launch_job_async(data[key]['query'], dump_to_file=False)
        table = job.get_results()
        outvalues[key] = table
        print('Retrieved {} sources for catalog {}'.format(len(table), key))
    gaia_mag_cols = ['phot_g_mean_mag', 'phot_bp_mean_mag', 'phot_rp_mean_mag']
    return outvalues['gaia'], gaia_mag_cols, outvalues['tmass'], outvalues['tmass_crossmatch'], outvalues['wise'], outvalues['wise_crossmatch']


def besancon(ra, dec, box_width, username='', kmag_limits=(13, 29)):
    """
    This routine calls a server to get a Besancon star count model over a given
    small sky area at a defined position.  For documentation of the Besancon
    model see the web site: http://model.obs-besancon.fr/
    The star count model returns V, J, H, K, and L magnitudes plus other
    information.  Here the magnitudes and the extinction are of most interest
    in producing the simulated NIRISS/NIRCam/Guider magnitudes.
    Note that the Besancon model uses a simplified exinction model, and so
    this is carried over into the Mirage point source lists.
    An email address is required for the besancon call.

    Parameters
    ----------
        ra : float
            Right ascension or galactic longitude for the sky
            position where the model is to be calculated, in degrees.
            Which it is depends on the "coords" parameter.

        dec : float
            Declinatinon or galactic latitude for the sky
            position where the model is to be calculated, in degrees.
            Which it is depends on the "coords" parameter.

        box_width : float
            Size of the (square) sky area to be simulated, in arc-seconds

        coords : str
            If "ra_dec", the default value, then the input values are assumed to
            be RA/Dec in degrees.  Otherwise the values are assumed to be (l,b)
            in degrees.

        kmag_limits : tup
            The range of allowable K magnitudes for the model, given
            as a tuple (min, max).  The default is (13,29).  The
            bright limit will generally be set by 2MASS completeness
            limit.  The 2MASS faint limit for any given sky position
            is roughly magnitude 15 to 16 in general.  As the 2MASS
            completeness limit varies with position and the
            completeness limit is above the faint limit the Besancon
            bright limit is taken as magnitude 14 by default.  Note
            that for the JWST instruments the 2MASS sources will
            saturate in full frame imaging in many cases.
    """
    from astropy import units as u

    # Specified coordinates. Will need to convert to galactic long and lat
    # when calling model
    ra = ra * u.deg
    dec = dec * u.deg
    box_width = box_width * u.arcsec

    min_ra = ra - box_width / 2
    max_ra = ra + box_width / 2
    min_dec = dec - box_width / 2
    max_dec = dec + box_width / 2

    # Define the list of colors to return
    colors = 'V-K,J-K,J-H,J-L'

    # Band minimum and maximum values correspond to the 9 Johnson-Cousins
    # filters used by the Besancon model: V, B, U, R, I, J, H, K, L
    # in that order.
    band_min = '{},-99.0,-99.0,-99.0,-99.0,-99.0,-99.0,-99.0,-99.0'.format(kmag_limits[0])
    band_max = '{},99.0,99.0,99.0,99.0,99.0,99.0,99.0,99.0'.format(kmag_limits[1])

    # Query the model
    path = os.path.dirname(__file__)
    client = os.path.join(path, 'galmod_client.py')
    command = (('python {} --url "https://model.obs-besancon.fr/ws/" --user {} '
               '--create -p KLEH 2 -p Coor1_min {} -p Coor2_min {} -p Coor1_max {} -p Coor2_max {} '
               '-p ref_filter K -p acol {} -p band_min {} -p band_max {} --run')
               .format(client, username, min_ra.value, min_dec.value, max_ra.value, max_dec.value,
                       colors, band_min, band_max))
    print('Running command: ', command)
    os.system(command)


def crop_besancon(ra, dec, box_width, catalog_file, ra_column_name='RAJ2000', dec_column_name='DECJ2000'):
    """Given a file containing an ascii source catalog, this function
    will crop the catalog such that it contains only sources within the
    provided RA/Dec range.

    Parameters
    ----------
    ra : float
        Right Ascension, in degrees, of the center of the area to keep

    dec : float
        Declination, in degerees, of the center of the area to keep

    box_width : float
        Full width, in arcseconds, of the area to be kept. We assume a
        square box.

    catalog_file : str
        Name of ascii file containing catalog

    ra_column_name : str
        Name of the column in the catalog containing the RA values

    dec_column_name : str
        Name of the column in the catalog containing the Dec values

    Returns
    -------
    cropped_catalog : astropy.table.Table
        Table containing cropped catalog
    """
    # Define min and max RA, Dec values to be kept
    ra = ra * u.deg
    dec = dec * u.deg
    box_width = box_width * u.arcsec

    min_ra = ra - box_width / 2
    max_ra = ra + box_width / 2
    min_dec = dec - box_width / 2
    max_dec = dec + box_width / 2

    # Read in input catalog
    catalog = ascii.read(catalog_file)

    # Check for requested columns
    if ra_column_name not in catalog.colnames:
        raise ValueError(('ERROR: {} does not contain a {} column.').format(catalog_file, ra_column_name))

    if dec_column_name not in catalog.colnames:
        raise ValueError(('ERROR: {} does not contain a {} column.').format(catalog_file, dec_column_name))

    # Keep only sources within the range of RA and Dec
    good = np.where((catalog[ra_column_name].data >= min_ra) & (catalog[ra_column_name].data <= max_ra) &
                    (catalog[dec_column_name].data >= min_dec) & (catalog[dec_column_name].data <= max_dec))[0]

    cropped_catalog = Table()
    for colname in catalog.colnames:
        cropped_catalog[colname] = catalog[colname][good]

    return cropped_catalog


def galactic_plane(box_width, instrument, filter_list, besancon_catalog_file,
                   ra_column_name='RAJ2000', dec_column_name='DECJ2000'):
    """Convenience function to create a typical scene looking into the disk of
    the Milky Way, using the besancon function.

    Parameters
    ----------
    box_width : float
        Size of the (square) sky area to be simulated, in arc-seconds

    instrument : str
        Name of instrument for the siumulation

    filter_list : list
        List of filters to use to generate the catalog. (e.g ['F480M', 'F444W'])

    besancon_catalog_file : str
        Name of ascii catalog containing background stars. The code was
        developed around this being a catalog output by the Besaoncon
        model (via ``besancon()``), but it can be any ascii catalog
        of sources as long as it contains the columns specified in
        the ``catalog`` parameter of ``johnson_catalog_to_mirage_catalog()``
        If None, the Besancon step will be skipped and a catalog will be
        built using only GAIA/2MASS/WISE

    ra_column_name : str
        Name of the column within ``besancon_catalog_file`` containing the
        right ascension of the sources.

    dec_column_name : str
        Name of the column within ``besancon_catalog_file`` containing the
        declination of the sources

    Returns
    -------
    cat : mirage.catalogs.create_catalog.PointSourceCatalog
        Catalog with representative stars pulled from Besancon model
    """
    galactic_longitude = 45.0 * u.deg  # deg
    galactic_latitude = 0.0 * u.deg  # deg
    coord = SkyCoord(galactic_longitude, galactic_latitude, frame=Galactic)

    cat, column_filter_list = get_all_catalogs(coord.icrs.ra.value, coord.icrs.dec.value, box_width,
                                               besancon_catalog_file=besancon_catalog_file, instrument=instrument,
                                               filters=filter_list, ra_column_name=ra_column_name,
                                               dec_column_name=dec_column_name)
    return cat


def out_of_galactic_plane(box_width, instrument, filter_list, besancon_catalog_file,
                          ra_column_name='RAJ2000', dec_column_name='DECJ2000'):
    """Convenience function to create typical scene looking out of the plane of
    the Milky Way by querying the Besancon model

    Parameters
    ----------
    box_width : float
        Size of the (square) sky area to be simulated, in arc-seconds

    instrument : str
        Name of instrument for the siumulation

    filter_list : list
        List of filters to use to generate the catalog. (e.g ['F480M', 'F444W'])

    besancon_catalog_file : str
        Name of ascii catalog containing background stars. The code was
        developed around this being a catalog output by the Besaoncon
        model (via ``besancon()``), but it can be any ascii catalog
        of sources as long as it contains the columns specified in
        the ``catalog`` parameter of ``johnson_catalog_to_mirage_catalog()``
        If None, the Besancon step will be skipped and a catalog will be
        built using only GAIA/2MASS/WISE

    ra_column_name : str
        Name of the column within ``besancon_catalog_file`` containing the
        right ascension of the sources.

    dec_column_name : str
        Name of the column within ``besancon_catalog_file`` containing the
        declination of the sources

    Returns
    -------
    cat : mirage.catalogs.create_catalog.PointSourceCatalog
        Catalog with representative stars pulled from Besancon model
    """
    galactic_longitude = 45.0 * u.deg  # deg
    galactic_latitude = 85.0 * u.deg  # deg
    coord = SkyCoord(galactic_longitude, galactic_latitude, frame=Galactic)

    cat, column_filter_list = get_all_catalogs(coord.icrs.ra.value, coord.icrs.dec.value, box_width,
                                               besancon_catalog_file=besancon_catalog_file, instrument=instrument,
                                               filters=filter_list, ra_column_name=ra_column_name,
                                               dec_column_name=dec_column_name)
    return cat


def galactic_bulge(box_width, instrument, filter_list, besancon_catalog_file,
                   ra_column_name='RAJ2000', dec_column_name='DECJ2000'):
    """Convenience function to create typical scene looking into bulge of
    the Milky Way

    Parameters
    ----------
    box_width : float
        Size of the (square) sky area to be simulated, in arc-seconds

    instrument : str
        Name of instrument for the siumulation

    filter_list : list
        List of filters to use to generate the catalog. (e.g ['F480M', 'F444W'])

    besancon_catalog_file : str
        Name of ascii catalog containing background stars. The code was
        developed around this being a catalog output by the Besaoncon
        model (via ``besancon()``), but it can be any ascii catalog
        of sources as long as it contains the columns specified in
        the ``catalog`` parameter of ``johnson_catalog_to_mirage_catalog()``
        If None, the Besancon step will be skipped and a catalog will be
        built using only GAIA/2MASS/WISE

    ra_column_name : str
        Name of the column within ``besancon_catalog_file`` containing the
        right ascension of the sources.

    dec_column_name : str
        Name of the column within ``besancon_catalog_file`` containing the
        declination of the sources

    Returns
    -------
    cat : mirage.catalogs.create_catalog.PointSourceCatalog
        Catalog with representative stars pulled from Besancon model
    """
    galactic_longitude = 0.0 * u.deg  # deg
    galactic_latitude = 5.0 * u.deg  # deg
    coord = SkyCoord(galactic_longitude, galactic_latitude, frame=Galactic)

    cat, column_filter_list = get_all_catalogs(coord.icrs.ra.value, coord.icrs.dec.value, box_width,
                                               besancon_catalog_file=besancon_catalog_file, instrument=instrument,
                                               filters=filter_list, ra_column_name=ra_column_name,
                                               dec_column_name=dec_column_name)
    return cat


def filter_bad_ra_dec(table_data):
    """Remove targets with bad RA and/or Dec values from the input table.
    Use the column masks to find which entries have bad RA or Dec values.
    These will be excluded from the Mirage catalog.

    Parameters
    ----------
    table_data : astropy.table.Table
        Input table from e.g. 2MASS query

    Returns
    -------
    position_mask : np.ndarray
        1D boolean array. True for good sources, False for bad.
    """
    ra_mask = ~table_data['ra'].data.mask
    dec_mask = ~table_data['dec'].data.mask
    position_mask = ra_mask & dec_mask
    return position_mask


def generate_ra_dec(number_of_stars, ra_min, ra_max, dec_min, dec_max, seed=None):
    """
    Generate a list of random RA, Dec values in a square region.  Note that
    this assumes a small sky area so that the change in the sky area per
    a degree of right ascension is negligible.  This routine will break down
    at the north or south celestial poles.
    The np uniform random number generator is used to make the source
    postions.

    Parameters
    ----------
        number_of_stars: int
            The number of stars to generate positions for

        ra_min : float
            The minimum RA value of the area, in degrees

        ra_max : float
            The maximum RA value of the area, in degrees

        dec_min : float
            The minimum Dec value of the area, in degrees

        dec_max : float
            The minimum Dec value of the area, in degrees

        seed : int
            Optional seed for random number generator

    Returns
    -------
        ra_list : numpy.ndarray
            The list of output RA values in degrees.

        dec_list : numpy.ndarray
            The list of output Dec values in degrees.
    """
    delta_ra = ra_max - ra_min
    delta_dec = dec_max - dec_min

    # If a seed is provided, seed the generator
    if seed is not None:
        np.random.seed(seed)

    # Generate random numbers
    numbers = np.random.random(2 * number_of_stars)

    # Create RA and Dec values
    ra_list = numbers[0: number_of_stars] * delta_ra + ra_min
    dec_list = numbers[number_of_stars:] * delta_dec + dec_min
    return ra_list, dec_list


def galaxy_background(ra0, dec0, v3rotangle, box_width, instrument, filters,
                      boxflag=True, brightlimit=14.0, seed=None):
    """
    Given a sky position (ra0,dec0), and a V3 rotation angle (v3rotangle) this
    routine makes a fake galaxy background for a square region of sky or a circle
    provided that the area is smaller than the GODDS-S field.  The fake galaxies
    are randomly distributed over the area and have random changes to the input
    Sersic parameters, so the field is different at every call.

    Parameters
    ----------
    ra0 : float
        The sky position right ascension in decimal degrees

    dec0 : float
        The sky position declimation in decimal degrees

    v3rotangle : float
        The v3 rotation angle of the y coordinate in the image, in
        degrees E of N

    box_width : float
        The width of the square on the sky in arc-seconds, or the
        radius of the circle on the sky in arc-seconds; must be a
        value larger and 1.0 and smaller than 779.04 for a square
        or 439.52 for a circle.

    boxflag : boolean
        Flag for whether the region is a square (if True) or a
        circle (if False).

    instrument : str
        One of "All", "NIRCam", "NIRISS", or "Guider".

    filters : list
        A list of the required filters.  If set to None or an empty
        list all filters are used.

    brightlimit : float
        The bright limit AB magnitude for the galaxies in the output
        catalogue, applied in the F200W filter (close to K-band).

    seed : integer or None
        If a value is given, it is used to seed the random number
        generator so that a scene can be reproduced.  If the value
        is None the seed is set randomly.  If a non-integer value is
        given it is converted to integer.

    Returns
    -------
    galaxy_cat : mirage.catalogs.catalog_generate.GalaxyCatalog
        This is the output list of galaxy properties.  If there is
        an error, None is returned.

    seedvalue : integer
        The seed value used with numpy.random to generate the values.
    """
    # The following is the area of the GOODS-S field catalogue from Gabe Brammer
    # in square arc-seconds
    goodss_area = 606909.
    if boxflag:
        outarea = box_width*box_width
    else:
        outarea = math.pi*box_width*box_width
    if outarea >= goodss_area:
        print('Error: requested sky area is too large.  Values will not be produced.')
        return None, None
    if seed is None:
        seedvalue = int(950397468.*np.random.random())
    else:
        if not isinstance(seed, int):
            seedvalue = int(abs(seed))
        else:
            seedvalue = seed
    np.random.seed(seedvalue)
    threshold = outarea/goodss_area
    filter_names = make_filter_names(instrument, filters)
    nfilters = len(filter_names)
    if nfilters < 1:
        print('Error matching filters to standard list.  Inputs are:')
        print('Instrument: ', instrument)
        print('Filter names: ', filters)
        return None, None
    # add 8 to these indexes to get the columns in the GODDS-S catalogue file
    #
    # Note: NIRCam filters are used in proxy for the NIRISS long wavelength
    # filters.  The broad NIRCam F150W2 filter is used as a proxy for the
    # Guiders.
    filterinds = {'niriss_f090w_magnitude': 0, 'niriss_f115w_magnitude': 1,
                  'niriss_f150w_magnitude': 2, 'niriss_f200w_magnitude': 3,
                  'niriss_f140m_magnitude': 4, 'niriss_f158m_magnitude': 5,
                  'nircam_f070w_magnitude': 6, 'nircam_f090w_magnitude': 7,
                  'nircam_f115w_magnitude': 8, 'nircam_f150w_magnitude': 9,
                  'nircam_f200w_magnitude': 10, 'nircam_f150w2_magnitude': 11,
                  'nircam_f140m_magnitude': 12, 'nircam_f162m_magnitude': 13,
                  'nircam_f182m_magnitude': 14, 'nircam_f210m_magnitude': 15,
                  'nircam_f164n_magnitude': 16, 'nircam_f187n_magnitude': 17,
                  'nircam_f212n_magnitude': 18, 'nircam_f277w_magnitude': 19,
                  'nircam_f356w_magnitude': 20, 'nircam_f444w_magnitude': 21,
                  'nircam_f322w2_magnitude': 22, 'nircam_f250m_magnitude': 23,
                  'nircam_f300m_magnitude': 24, 'nircam_f335m_magnitude': 25,
                  'nircam_f360m_magnitude': 26, 'nircam_f410m_magnitude': 27,
                  'nircam_f430m_magnitude': 28, 'nircam_f460m_magnitude': 29,
                  'nircam_f480m_magnitude': 30, 'nircam_f323n_magnitude': 31,
                  'nircam_f405n_magnitude': 32, 'nircam_f466n_magnitude': 33,
                  'nircam_f470n_magnitude': 34, 'niriss_f277w_magnitude': 19,
                  'niriss_f356w_magnitude': 20, 'niriss_f380m_magnitude': 26,
                  'niriss_f430m_magnitude': 28, 'niriss_f444w_magnitude': 21,
                  'niriss_f480m_magnitude': 30, 'fgs_guider1_magnitude': 11,
                  'fgs_guider2_magnitude': 11}
    module_path = pkg_resources.resource_filename('mirage', '')
    catalog_file = os.path.join(module_path, 'config/goodss_3dhst.v4.1.jwst_galfit.cat')
    catalog_values = np.loadtxt(catalog_file, comments='#')
    outinds = np.zeros((nfilters), dtype=np.int16)
    try:
        loop = 0
        for filter in filter_names:
            outinds[loop] = filterinds[filter] + 8
            loop = loop+1
    except:
        print('Error matching filter %s to standard list.' % (filter))
        return None, None
    # The following variables hold the Sersic profile index values
    # (radius [arc-seconds], sersic index, ellipticity, position angle)
    # and the assocated uncertaintie index values
    sersicinds = [59, 61, 63, 65]
    sersicerrorinds = [60, 62, 64, 66]
    ncat = catalog_values.shape[0]
    select = np.random.random(ncat)
    magselect = np.copy(catalog_values[:, filterinds['niriss_f200w_magnitude']+8])
    outputinds = []
    for loop in range(ncat):
        if (magselect[loop] >= brightlimit) and (select[loop] < threshold):
            outputinds.append(loop)
    nout = len(outputinds)
    if boxflag:
        delx0 = (-0.5+np.random.random(nout))*box_width/3600.
        dely0 = (-0.5+np.random.random(nout))*box_width/3600.
        radius = np.sqrt(delx0*delx0+dely0*dely0)
        angle = np.arctan2(delx0, dely0)+v3rotangle*math.pi/180.
    else:
        radius = box_width*np.sqrt(np.random.random(nout))/3600.
        angle = 2.*math.pi*np.random.random(nout)
    delx = radius*np.cos(angle)
    dely = radius*np.sin(angle)
    raout = delx * 0.
    decout = dely * 0.
    for loop in range(len(delx)):
        raout[loop], decout[loop] = deproject_from_tangent_plane(delx[loop], dely[loop], ra0, dec0)
    rot1 = 360.*np.random.random(nout)-180.
    rout = np.copy(catalog_values[outputinds, sersicinds[0]])
    drout = np.copy(catalog_values[outputinds, sersicerrorinds[0]])
    rout = rout+2.*drout*np.random.normal(0., 1., nout)
    rout[rout < 0.01] = 0.01
    elout = np.copy(catalog_values[outputinds, sersicinds[2]])
    delout = np.copy(catalog_values[outputinds, sersicerrorinds[2]])
    elout = elout+delout*np.random.normal(0., 1., nout)
    elout[elout > 0.98] = 0.98
    sindout = np.copy(catalog_values[outputinds, sersicinds[1]])
    dsindout = np.copy(catalog_values[outputinds, sersicinds[1]])
    sindout = sindout+dsindout*np.random.normal(0., 1., nout)
    sindout[sindout < 0.1] = 0.1
    paout = np.copy(catalog_values[outputinds, sersicinds[3]])
    dpaout = np.copy(catalog_values[outputinds, sersicinds[3]])
    paout = paout+dpaout*np.random.normal(0., 1., nout)
    for loop in range(len(paout)):
        if paout[loop] < -180.:
            paout[loop] = paout[loop]+360.
        if paout[loop] > 180.:
            paout[loop] = paout[loop]-360.
    galaxy_cat = GalaxyCatalog(ra=raout, dec=decout, ellipticity=elout,
                               radius=rout, sersic_index=sindout,
                               position_angle=paout, radius_units='arcsec')
    for loop in range(len(filter_names)):
        mag1 = catalog_values[outputinds, outinds[loop]]
        dmag1 = -0.2*np.random.random(nout)+0.1
        mag1 = mag1 + dmag1
        if 'niriss' in filter_names[loop]:
            inst1 = 'NIRISS'
            filter_name = filter_names[loop].split('_')[1].upper()
        elif 'nircam' in filter_names[loop]:
            inst1 = 'NIRCam'
            filter_name = filter_names[loop].split('_')[1].upper()
        else:
            inst1 = 'FGS'
            filter_name = filter_names[loop].split('_')[1].upper()
        galaxy_cat.add_magnitude_column(mag1, instrument=inst1,
                                        filter_name=filter_name,
                                        magnitude_system='abmag')
    return galaxy_cat, seedvalue
