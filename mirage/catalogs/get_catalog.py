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
import copy
import os

from astropy.io import ascii as asc
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Column
from astroquery.vizier import Vizier
import numpy as np
from numpy import ma

from ..apt import apt_inputs
from ..apt.read_apt_xml import ReadAPTXML
from ..utils.utils import ensure_dir_exists


def get_target_coords(xml_filename):
    '''Parse APT xml file to extract RA and Dec from each target

    Parameters
    ----------
    xml_filename : str
        Path to APT .xml file

    Returns
    -------
    target_coords : list
        List of astropy SkyCoord objects for all targets
    '''
    xml_info = ReadAPTXML()
    xml_info.read_xml(xml_filename)
    targ_dict = xml_info.target_info
    target_list = []
    target_coords = []
    for key in targ_dict:
        target_list.append(key)
        ra_str, dec_str = targ_dict[key]
        target_coords.append(SkyCoord('{} {}'.format(ra_str, dec_str), unit=(u.hourangle, u.deg)))

    return target_list, target_coords


def get_sw_catalog(target_coords, output_directory, search_radius=15 * u.arcmin):
    '''Query 2MASS catalog and write out shortwave catalog file for each target

    Parameters
    ----------
    target_coords : list
        List of astropy SkyCoord objects for all targets

    Returns
    -------
    catalog_filenames_sw : list
        List of filenames for all shortwave catalogs

    output_directory : str
        Name of directory into which the catalog files are written

    search_radius : astropy.units.quantity
        Search radius to use in query

    References
    ----------
    Catalog: http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=II/246
    '''
    catalog_filenames_sw = []

    # For each target
    for i_obs, t in enumerate(target_coords):
        # Generate the catalog filename
        catalog_filename_sw = os.path.join(output_directory,
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
        queried_catalog_sw.rename_column('Kmag', 'nircam_f212_nmagnitude')

        # Add index column
        col_ind = Column(np.arange(1, len(queried_catalog_sw)+1), name='index')
        queried_catalog_sw.add_column(col_ind, index=0)

        # Save shortwave catalog and header comments to file
        comments = [
            '',
            'vegamag',
            'Catalog for Mirage. created from 2MASS catalog',
            'Matching NIRCam SW F212N filter',
            'Magnitudes are Kmag (2-3 microns).'
        ]
        queried_catalog_sw.meta['comments'] = comments
        asc.write(queried_catalog_sw, output=catalog_filename_sw, overwrite=True)
        print('Successfully saved shortwave catalog to {}.'.format(catalog_filename_sw))

        catalog_filenames_sw.append(catalog_filename_sw)

    return catalog_filenames_sw


def get_catalog(target_coords, output_directory, search_radius=15 * u.arcmin):
    '''Query WISE and 2MASS catalogs and write out catalog file for each target

    Parameters
    ----------
    target_coords : list
        List of astropy SkyCoord objects for all targets

    output_directory : str
        Name of directory into which the catalog files are written

    search_radius : astropy.units.quantity
        Search radius to use in query

    Returns
    -------
    catalog_filenames_lw : list
        List of filenames for all longwave catalogs

    References
    ----------
    Catalog: http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=II/328
    '''
    catalog_filenames = []

    # For each target
    for i_obs, t in enumerate(target_coords):
        # Generate the catalog filename
        catalog_filename = os.path.join(output_directory,
                                        'WISE_and_2MASS_RA{:.2f}deg_Dec{:.2f}deg.list'.format(t.ra.degree,
                                                                                              t.dec.degree))

        # Check that the shortwave catalog has not already been created
        if os.path.isfile(catalog_filename):
            catalog_filenames.append(catalog_filename)
            if catalog_filename not in catalog_filenames[:i_obs]:
                print('Longwave catalog file {} already exists. Will not overwrite.'.\
                      format(os.path.basename(catalog_filename)))
            continue

        # If not, query longwave sources from the WISE catalog from Vizier
        v = Vizier(catalog='II/328/allwise', columns=['RAJ2000', 'DEJ2000', 'W2mag', 'Kmag'])
        v.ROW_LIMIT = -1
        result = v.query_region(t, radius=search_radius)
        queried_catalog = copy.deepcopy(result['II/328/allwise'])

        print('Queried {} WISE objects within {} {} of RA, Dec ({:.2f}, {:.2f}).'.
              format(len(queried_catalog), search_radius.value,
                     str(search_radius.unit), t.ra.deg, t.dec.deg))

        # Deal with masked values
        queried_catalog['Kmag'].fill_value = 99.
        queried_catalog['W2mag'].fill_value = 99.
        queried_catalog = queried_catalog.filled()

        # Rename table columns to match nircam simulator expectations
        queried_catalog.rename_column('RAJ2000', 'x_or_RA')
        queried_catalog.rename_column('DEJ2000', 'y_or_Dec')
        queried_catalog.rename_column('W2mag', 'nircam_f480m_magnitude')
        queried_catalog.rename_column('Kmag', 'nircam_f212n_magnitude')

        # Duplicate F480M column for F460M
        f460m = Column(queried_catalog['nircam_f480m_magnitude'].data, name='nircam_f460m_magnitude')
        queried_catalog.add_column(f460m, index=3)

        # Add index column
        col_ind = Column(np.arange(1, len(queried_catalog)+1), name='index')
        queried_catalog.add_column(col_ind, index=0)

        # Save longwave catalog and header comments to file
        comments = [
            '',
            'vegamag',
            'Catalog for Mirage. Created from WISE W2 band and 2MASS Ks band',
            'Matching NIRCam F460M & F480M, and F212N filters respectively',
            'No interpolation done to change magnitudes from the original bandpasses to NIRCam bandpasses'
        ]
        queried_catalog.meta['comments'] = comments

        asc.write(queried_catalog, output=catalog_filename, overwrite=True)
        print('Successfully saved catalog to {}.'.format(catalog_filename))

        catalog_filenames.append(catalog_filename)

    return catalog_filenames


def get_all_catalogs(xml_file, out_dir='./'):
    '''Query WISE and 2MASS catalogs and write out
    catalog files for each target in the provided APT proposal

    Parameters
    ----------
    xml_file : str
        Path to APT .xml file

    prop_id : str
        APT proposal ID number

    out_dir : str
        Output directory for catalog files

    Returns
    -------
    tuple
        Lists of shortwave and longwave catalog filenames
    '''
    target_list, target_coords = get_target_coords(xml_file)

    ensure_dir_exists(out_dir)
    #catalog_filenames_sw = get_sw_catalog(target_coords, out_dir)
    catalog_filenames = get_catalog(target_coords, out_dir)

    return target_list, target_coords, catalog_filenames
