#! /usr/bin/env python

"""
Tools for generating Mirage-compatible catalogs from surveys

Note that once you have mirage-formatted catalogs, there is an "add_catalog" function
in catalog_generator.py that can combine catalogs

"""

from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.gaia import Gaia
from astroquery.irsa import Irsa
import numpy as np

from mirage.catalogs.catalog_generator import PointSourceCatalog, GalaxyCatalog


def ptsrc_for_proposal(pointing_dictionary, catalog_splitting_threshold=1., email=''):
    """
    NOT YET FUNCTIONAL

    Given a pointing dictionary from an APT file, generate source catalogs
    that cover all of the coordinates specifired.

    Parameters
    ----------
    pointing_dictionary : dict
        Output from miarge.apt.apt_inputs.ra_dec_update()

    catalog_splitting_threshold : float
        Maximum distance in degrees between two pointings where a single catalog will contain
        both pointings

    Returns
    -------
    something
    """
    threshold = catalog_splitting_threshold * u.deg
    ra_apertures = pointing_dictionary['ra_ref'] * u.deg
    dec_apertures = pointing_dictionary['dec_ref'] * u.deg
    ra_target = pointing_dictionary['ra'] * u.deg
    dec_target = pointing_dictionary['dec'] * u.deg
    mapped = np.array([False] * len(ra_target))
    index = np.arange(len(ra_target))

    for indx, target_ra, target_dec in zip(index, ra_target, dec_target):
        if mapped[indx] is False:
            dithers = ra_target == target_ra & dec_target == target_dec
            min_ra = np.min(ra_apertures[dithers])
            max_ra = np.max(ra_apertures[dithers])
            min_dec = np.min(dec_apertures[dithers])
            max_dec = np.max(dec_apertures[dithers])

            nearby_targets = separation(target_ra, target_dec to ra_target, dec_target) < catalog_splitting_threshold &
                            mapped is False
            for close_target in nearby_targets:
                dithers = ra_target == ra_target[close_target] & dec_target == dec_target[close_target]
                min_ra = np.min([ra_apertures[dithers]].append(min_ra))
                max_ra = np.max([ra_apertures[dithers]].append(max_ra))
                min_dec = np.min([dec_apertures[dithers]].append(min_dec))
                max_dec = np.max([dec_apertures[dithers]].append(max_dec))
            mapped[nearby_targets] = True

            # Pad min and max RA and Dec values since so far we have values only at the reference
            # location. Add at least half a (LW) detector width * sqrt(2).
            pad = 0.062 * 1024 * 1.5
            min_ra -= pad
            max_ra += pad
            min_dec -= pad
            max_dec += pad
            mean_ra = (max_ra + min_ra) / 2.
            mean_dec = (max_dec + min_dec) / 2.
            delta_ra = max_ra - min_ra
            delta_dec = max_dec - min_dec
            full_width = np.max([delta_ra, delta_dec])

            # Create catalog
            cat = twoMASS_plus_background(mean_ra, mean_dec, full_width, email=email)
            cat_filename = '???????'
            cat_dir = '???????'
            cat.save(os.path.join(cat_dir, cat_filename))


def query_gaia(ra, dec, box_width, frame='icrs'):
    """Query the Gaia archive in a square region around the RA and Dec provided.

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

    cat : catalog_generator:PointSourceCatalog object
        Catalog containing Gaia sources
    """
    if isinstance(ra, float):
        ra_units = u.deg
    elif isinstance(ra, str):
        ra_units = u.hourangle
    else:
        raise ValueError(("WARNING: RA must be a decimal degrees or a string (e.g. '12h23m34.5s', "
                          "'12:23:34.5'"))
    if isinstance(dec, float):
        dec_units = u.deg
    elif isinstance(dec, str):
        dec_units = u.hourangle
    else:
        raise ValueError(("WARNING: Dec must be a decimal degrees or a string (e.g. '12d23m34.5s', "
                          "'12:23:34.5'"))

    coord = SkyCoord(ra=ra, dec=dec, unit=(ra_units, dec_units), frame=frame)
    width = u.Quantity(box_width, u.arcsec)
    height = u.Quantity(box_width, u.arcsec)
    gaia_cat = Gaia.query_object_async(coordinate=coord, width=width, height=height)

    magnitude_column_names = ['phot_g_mean_mag', 'phot_bp_mean_mag', 'phot_rp_mean_mag']
    # Place the results into a Mirage-formatted catalog
    #cat = PointSourceCatalog(ra=gaia_cat['ra'].data.data, dec=gaia_cat['dec'].data.data)
    #for key in ['phot_g_mean_mag', 'phot_bp_mean_mag', 'phot_rp_mean_mag']:
    #    filter_name = key.split('_')[1]
    #    cat.add_magnitude_column(gaia_cat[key].data.data, instrument='Gaia', filter_name=filter_name)

    return gaia_cat, magnitude_column_names


def query_2MASS_ptsrc_catalog(ra, dec, box_width):
    """Query the 2MASS All-Sky Point Source Catalog in a square region around the RA and Dec
    provided. Box width must be in units of arcseconds

    Parameters
    ----------
    something
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
    #cat = PointSourceCatalog(ra=query_table['ra'].data.data.data[radec_mask],
    #                         dec=query_table['dec'].data.data.data[radec_mask])
    #
    # Add the J, H, and K magnitudes as they may be useful for magnitude conversions later
    # Add the values that have had fill_values applied. The fill_value is 1e20.
    #for key in ['j_m', 'h_m', 'k_m']:
    #    data = query_table[key].filled().data
    #    cat.add_magnitude_column(data, instrument='2MASS', filter_name=key)

    return query_table, magnitude_column_names


def mirage_ptsrc_catalog_from_table(table, instrument, mag_colnames):
    """Create a mirage-formatted point source catalog from an input
    table (e.g. from one of the query functions), along with the magnitude
    column names of interest

    Parameters
    ----------

    instrument : str
        Unique identifier for where data came from (e.g. '2MASS', 'WISE', 'nircam_f200w')
    """
    cat = PointSourceCatalog(ra=table['ra'].data.data, dec=table['dec'].data.data)

    for magcol in mag_colnames:
        data = table[key].filled().data
        cat.add_magnitude_column(data, instrument=instrument, filter_name=key)
    return cat


def get_2MASS_ptsrc_catalog(ra, dec, box_width):
    """Wrapper around 2MASS query and creation of mirage-formatted catalog"""
    twomass_cat, twomass_mag_cols = query_2MASS_ptsrc_catalog(ra, dec, box_width)
    twomass_mirage = mirage_ptsrc_catalog_from_table(twomass_cat, '2MASS', twomass_mag_cols)
    return twomass_mirage


def twoMASS_plus_background(ra, dec, box_width, kmag_limits=(17, 29), email=''):
    """Convenience function to create a catalog from 2MASS and add a population of
    fainter stars. In this case, cut down the magnitude limits for the call to the
    Besancon model so that we don't end up with a double population of bright stars
    """
    two_mass = get_2MASS_ptsrc_catalog(ra, dec, box_width)
    background = besancon(ra, dec, box_width, coords='ra_dec', email=email)
    two_mass.add_catalog(background)
    return two_mass


def get_gaia_ptsrc_catalog(ra, dec, box_width):
    """Wrapper around Gaia query and creation of mirage-formatted catalog"""
    gaia_cat, gaia_mag_cols = query_2MASS_ptsrc_catalog(ra, dec, box_width)
    gaia_mirage = mirage_ptsrc_catalog_from_table(gaia_cat, 'gaia', gaia_mag_cols)
    return gaia_mirage


def besancon(ra, dec, box_width, coords='ra_dec', email='', kmag_limits=(10, 29)):
    """test it out
    coords keyword can be 'ra_dec', which is default, or
    'galactic' in which case ra and dec are interpreted as galactic longitude and
    latitude

    email is required for the besancon call. kind of annoying.

    """
    from astroquery.besancon import Besancon
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    # Specified coordinates. Will need to convert to galactic long and lat
    # when calling model
    ra = ra * u.deg
    dec = dec * u.deg
    box_width = box_width * u.arcsec

    if coords == 'ra_dec':
        location = SkyCoord(ra=ra, dec=dec, frame='icrs')
        coord1 = location.galactic.l.value
        coord2 = location.galactic.b.value
    elif coords == 'galactic':
        coord1 = ra.value
        coord2 = dec.value

    # Area of region to search (model expects area in square degrees)
    area = box_width * box_width
    area = area.to(u.deg * u.deg)

    # Query the model
    model = Besancon.query(coord1, coord2, smallfield=True, area=area.value,
                           colors_limits={"J-H": (-99, 99), "J-K": (-99, 99),
                                          "J-L": (-99, 99), "V-K": (-99, 99)},
                           mag_limits={'K': kmag_limits}, retrieve_file=True, email=email)

    # Calculate magnitudes in given bands
    k_mags = (model['V'] - model['V-K']).data
    j_mags = k_mags + model['J-K'].data
    h_mags = j_mags - model['J-H'].data
    l_mags = j_mags - model['J-L'].data

    # Since these are theoretical stars generated by a model, we need to provide RA and Dec values.
    # The model is run in 'smallfield' mode, which assumes a constant density of stars across the given
    # area. So let's just select RA and Dec values at random across the fov.
    half_width = box_width * 0.5
    min_ra = ra - half_width
    max_ra = ra + half_width
    min_dec = dec - half_width
    max_dec = dec + half_width
    ra_values, dec_values = generate_ra_dec(len(k_mags), min_ra, max_ra, min_dec, max_dec)

    # Create the catalog object
    cat = PointSourceCatalog(ra=ra_values, dec=dec_values)

    # Add the J, H, K and L magnitudes as they may be useful for magnitude conversions later
    cat.add_magnitude_column(j_mags, instrument='Besancon', filter_name='j')
    cat.add_magnitude_column(h_mags, instrument='Besancon', filter_name='h')
    cat.add_magnitude_column(k_mags, instrument='Besancon', filter_name='k')
    cat.add_magnitude_column(l_mags, instrument='Besancon', filter_name='l')
    return cat


def galactic_plane(box_width, email=''):
    """Convenience function to create a typical scene looking into the disk of
    the Milky Way, using the besancon function

    Parameters
    ----------
    box_width : float

    Returns
    -------
    cat : obj
        mirage.catalogs.create_catalog.PointSourceCatalog

    RA and Dec values of various features

    Center of MW
    center_ra = 17h45.6m
    center_dec = -28.94 deg

    anti-center-ra = 5h45.6m
    anti-center-dec = 28.94deg

    Galactic Poles
    north_pole_ra = 12h51.4m
    north_pole_dec = 27.13deg

    south_pole_ra = 0h51.4m
    south_pole_dec = -27.13deg
    """
    representative_galactic_longitude = 45.0  # deg
    representative_galactic_latitude = 0.0  # deg

    cat = besancon(representative_galactic_longitude, representative_galactic_latitude,
                   box_width, coords='galactic', email=email)
    return cat


def out_of_plane(box_width, email=''):
    """Convenience function to create typical scene looking out of the plane of
    the Milky Way

    Parameters
    ----------
    box_width : float

    Returns
    -------
    cat : obj
        mirage.catalogs.create_catalog.PointSourceCatalog
    """
    representative_galactic_longitude = 45.0  # deg
    representative_galactic_latitude = 85.0  # deg

    cat = besancon(representative_galactic_longitude, representative_galactic_latitude,
                   box_width, coords='galactic', email=email)
    return cat


def galactic_bulge(box_width, email=''):
    """Convenience function to create typical scene looking into bulge of
    the Milky Way

    Parameters
    ----------
    box_width : float

    Returns
    -------
    cat : obj
        mirage.catalogs.create_catalog.PointSourceCatalog
    """
    #Look up Besancon limitations. Model breaks down somewhere close to the
    #galactic core.
    representative_galactic_longitude = 0.  # ? deg
    representative_galactic_latitude = 5.0  # ? deg

    cat = besancon(representative_galactic_longitude, representative_galactic_latitude,
                   box_width, coords='galactic', email=email)
    return cat

#def from_luminosity_function(self, luminosity_function):
#    more customizable


def filter_bad_ra_dec(table_data):
    """Use the column masks to find which entries have bad RA or Dec values.
    These will be excluded from the Mirage catalog

    Parameters
    ----------
    something

    Returns
    -------

    position_mask : numpy.ndarray
        1D boolean array. True for good sources, False for bad.
    """
    #ra_data = table_data['ra'].data.data
    ra_mask = ~table_data['ra'].data.mask
    #dec_data = table_data['dec'].data.data
    dec_mask = ~table_data['dec'].data.mask
    position_mask = ra_mask & dec_mask
    return position_mask


def generate_ra_dec(number_of_stars, ra_min, ra_max, dec_min, dec_max):
    """Generate a list of random RA, Dec values"""
    delta_ra = ra_max - ra_min
    ra_list = np.random.random(number_of_stars) * delta_ra + ra_min
    delta_dec = dec_max - dec_min
    dec_list = np.random.random(number_of_stars) * delta_dec + dec_min
    return ra_list, dec_list




