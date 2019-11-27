'''Create necessary catalog files for targets in a given APT proposal.

Extracts the RA and Dec of all observed targets from the APT .pointing
file, and creates a catalog file that lists the sources within a
defined radius of each of these targets in mirage/catalogs/.
This module uses astroquery and Vizier to query 2MASS and WISE for
short-wave (F212N) and long-wave (F480M) catalogs, respectively, and
writes .list files.

Authors
-------
    - Lauren Chambers

Use
---
    This module can be executed in a Python shell as such:

    ::
        from mirage.catalogs import get_catalog
        sw_cats, lw_cats = get_catalog.get_all_catalogs(pointing_file, prop_id)

    Required arguments:
        ``pointing file`` - filepath to the APT pointing file
        ``prop_id`` - APT proposal ID number

References
----------
    Example of a query with Vizier:
        http://astroquery.readthedocs.io/en/latest/vizier/vizier.html
'''
import os

from astropy.io import ascii as asc
from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier

from ..apt import apt_inputs

SCRIPTS_DIR = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
PACKAGE_DIR = os.path.dirname(SCRIPTS_DIR)
CATALOG_DIR = os.path.join(PACKAGE_DIR, 'catalogs')


def get_target_coords(pointing_file, prop_id):
    '''Parse APT file to extract RA and Dec from each observation

    Parameters
    ----------
    pointing_file : str
        Path to APT .pointing file
    prop_id : str
        APT proposal ID number

    Returns
    -------
    target_coords : list
        List of astropy SkyCoord objects for all targets
    '''
    APT = apt_inputs.AptInput()
    pointing = APT.get_pointing_info(pointing_file, prop_id)

    # Split up by observation, finding all RA and Decs in pointing file
    obs_nums = []
    coords = []
    for ra, dec, obs_id in zip(pointing['ra'], pointing['dec'], pointing['obs_num']):
        if obs_id not in obs_nums:
            obs_nums.append(obs_id)
            coords.append((ra, dec))
        else:
            # If there are different RA/Decs within an observation.... eek
            if (ra != coords[-1][0]) or (dec != coords[-1][1]):
                print('** WARNING: Different RA/Decs within observation {} **'.format(obs_id))

    target_coords = []
    for ra, dec in coords:
        target_coords.append(SkyCoord(ra, dec, unit=(u.deg, u.deg)))

    for tup in set([(t.ra.deg, t.dec.deg) for t in target_coords]):
        print('Target coordinates: ', tup[0], tup[1])

    return target_coords


def get_sw_catalog(target_coords, search_radius=15 * u.arcmin):
    '''Query 2MASS catalog and write out shortwave catalog file for each target

    Parameters
    ----------
    target_coords : list
        List of astropy SkyCoord objects for all targets

    Returns
    -------
    catalog_filenames_sw : list
        List of filenames for all shortwave catalogs

    References
    ----------
    Catalog: http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=II/246
    '''
    catalog_filenames_sw = []

    # For each target
    for i_obs, t in enumerate(target_coords):
        # Generate the catalog filename
        catalog_filename_sw = os.path.join(CATALOG_DIR,
                                           '2MASS_RA{:.2f}deg_Dec{:.2f}deg.list'.format(t.ra.degree,
                                                                                        t.dec.degree))

        # Check that the shortwave catalog has not already been created
        if os.path.isfile(catalog_filename_sw):
            catalog_filenames_sw.append(catalog_filename_sw)
            if catalog_filename_sw not in catalog_filenames_sw[:i_obs]:
                print('Shortwave catalog file {} already exists. Will not overwrite.'.\
                      format(os.path.basename(catalog_filename_sw)))
            continue

        # If not, query shortwave sources from the 2MASS catalog from Vizier
        v = Vizier(catalog='II/246/out', columns=['_RAJ2000', '_DEJ2000', 'Kmag'])
        v.ROW_LIMIT = -1
        result = v.query_region(t, radius=search_radius)
        queried_catalog_sw = result['II/246/out']

        print('Queried {} 2MASS objects within {} {} of RA, Dec ({:.2f}, {:.2f}).'.
              format(len(queried_catalog_sw), search_radius.value,
                     str(search_radius.unit), t.ra.deg, t.dec.deg))

        # Rename table columns to match mirage's expectations
        queried_catalog_sw.rename_column('_RAJ2000', 'x_or_RA')
        queried_catalog_sw.rename_column('_DEJ2000', 'y_or_Dec')
        queried_catalog_sw.rename_column('Kmag', 'magnitude')

        # Save shortwave catalog and header comments to file
        comments = [
            '',
            'vegamag',
            'Catalog for ramp_simulator.py. created from 2MASS catalog',
            'Matching NIRCam SW F212N filter',
            'Magnitudes are Kmag (2-3 microns).'
        ]
        queried_catalog_sw.meta['comments'] = comments
        asc.write(queried_catalog_sw, output=catalog_filename_sw, overwrite=True)
        print('Successfully saved shortwave catalog to {}.'.format(catalog_filename_sw))

        catalog_filenames_sw.append(catalog_filename_sw)

    return catalog_filenames_sw


def get_lw_catalog(target_coords, search_radius=15 * u.arcmin):
    '''Query WISE catalog and write out longwave catalog file for each target

    Parameters
    ----------
    target_coords : list
        List of astropy SkyCoord objects for all targets

    Returns
    -------
    catalog_filenames_lw : list
        List of filenames for all longwave catalogs

    References
    ----------
    Catalog: http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=II/328
    '''
    catalog_filenames_lw = []

    # For each target
    for i_obs, t in enumerate(target_coords):
        # Generate the catalog filename
        catalog_filename_lw = os.path.join(CATALOG_DIR,
                                           'WISE_RA{:.2f}deg_Dec{:.2f}deg.list'.format(t.ra.degree,
                                                                                       t.dec.degree))

        # Check that the shortwave catalog has not already been created
        if os.path.isfile(catalog_filename_lw):
            catalog_filenames_lw.append(catalog_filename_lw)
            if catalog_filename_lw not in catalog_filenames_lw[:i_obs]:
                print('Longwave catalog file {} already exists. Will not overwrite.'.\
                      format(os.path.basename(catalog_filename_lw)))
            continue

        # If not, query longwave sources from the WISE catalog from Vizier
        v = Vizier(catalog='II/328/allwise', columns=['RAJ2000', 'DEJ2000', 'W2mag'])
        v.ROW_LIMIT = -1
        result = v.query_region(t, radius=search_radius)
        queried_catalog_lw = result['II/328/allwise']

        print('Queried {} WISE objects within {} {} of RA, Dec ({:.2f}, {:.2f}).'.
              format(len(queried_catalog_lw), search_radius.value,
                     str(search_radius.unit), t.ra.deg, t.dec.deg))

        # Rename table columns to match nircam simulator expectations
        queried_catalog_lw.rename_column('RAJ2000', 'x_or_RA')
        queried_catalog_lw.rename_column('DEJ2000', 'y_or_Dec')
        queried_catalog_lw.rename_column('W2mag', 'magnitude')

        # Save longwave catalog and header comments to file
        comments = [
            '',
            'abmag',
            'Catalog for ramp_simulator.py. created from WISE catalog',
            'Matching NIRCam LW F460M & F480M filters',
            'Magnitudes are W2 (4-8 microns).'
        ]
        queried_catalog_lw.meta['comments'] = comments
        asc.write(queried_catalog_lw, output=catalog_filename_lw, overwrite=True)
        print('Successfully saved longwave catalog to {}.'.format(catalog_filename_lw))

        catalog_filenames_lw.append(catalog_filename_lw)

    return catalog_filenames_lw


def get_all_catalogs(pointing_file, prop_id):
    '''Query WISE and 2MASS catalogs and write out long- and shortwave
    catalog files for each target in the provided APT proposal

    Parameters
    ----------
    pointing_file : str
        Path to APT .pointing file
    prop_id : str
        APT proposal ID number

    Returns
    -------
    tuple
        Lists of shortwave and longwave catalog filenames
    '''
    target_coords = get_target_coords(pointing_file, prop_id)
    catalog_filenames_sw = get_sw_catalog(target_coords)
    catalog_filenames_lw = get_lw_catalog(target_coords)

    return target_coords, catalog_filenames_sw, catalog_filenames_lw
