#! /usr/bin/env python

"""This module contains functions for creating Mirage-formatted source
catalogs

maybe i should change self.ra to self._ra and same with self.dec, to encourage
people to use the get_ra() and get_dec() functions

Right now create_table returns a table, which is fine, but would also be nice
to be able to do:

tab = PointSourceCatalog(stuff)
tab_ap_format = tab.create_table()
tab.save() - and have tab.save save the tab_ap_format table


Maybe:
tab = PointSourceCatalog(stuff)
tab.create_table()  # creates tab.astropy_table or something like that
tab.save()

"""

import copy
import os

from astropy.io import ascii
from astropy.table import Table, Column
import numpy as np


class PointSourceCatalog():
    def __init__(self, ra=[], dec=[], x=[], y=[]):
        """Initialize the point source catalog. Users can enter lists of RA and Dec values
        or x and y values for source potisions

        Parameters
        ----------
        something
        """
        # Make sure we are working with numpy arrays
        ra = np.array(ra)
        dec = np.array(dec)
        x = np.array(x)
        y = np.array(y)

        if len(ra) != len(dec):
            raise ValueError(("WARNING: inconsistent length between RA and Dec inputs."))
        if len(x) != len(y):
            raise ValueError(("WARNING: inconsistent length between x and y inputs."))
        if len(ra) > 0 and len(x) > 0:
            raise ValueError(("WARNING: Provide either RA and Dec values, or x and y values."))

        if len(ra) > 0:
            self._ra = ra
            self._dec = dec
        else:
            self._ra = x
            self._dec = y
        self.magnitudes = {}

        # Determine the units for the location fields. All that Mirage needs to know is whether
        # the units are pixels or not. Degrees vs hour angle is not important at this point.
        if len(ra) > 0:
            self.location_units = 'position_RA_Dec'
        if len(x) > 0:
            self.location_units = 'position_pixels'

    def add_catalog(self, catalog_to_add, magnitude_fill_value=99.):
        """Add a catalog to the current catalog instance"""
        # The the source positions in the two catalogs have different units, then the catalogs
        # can't be combined.
        if self.location_units != catalog_to_add.location_units:
            raise ValueError("WARNING: Sources in the catalogs do not have matching units (RA/Dec or x/y)")

        # Check the magnitude system from each and make sure they match
        current_mag_labels = list(self.magnitudes.keys())
        new_mag_labels = list(catalog_to_add.magnitudes.keys())
        mag_sys = self.magnitudes[current_mag_labels[0]][0]
        if mag_sys != catalog_to_add.magnitudes[new_mag_labels[0]][0]:
            print("WARNING: Magnitude systems of the two catalogs do not match. Cannot combine.")

        # Get the length of the two catalogs
        current_length = len(self._ra)
        new_length = len(catalog_to_add._ra)

        # Combine location columns
        self._ra = np.append(self._ra, catalog_to_add._ra)
        self._dec = np.append(self._dec, catalog_to_add._dec)

        # Now we need to compare magnitude columns. Columns common to both catalogs can be
        # combined. Columns not common will have to have fill values added so that everything
        # is the same length
        orig_current_mag_labels = copy.deepcopy(current_mag_labels)
        current_mag_labels.extend(new_mag_labels)
        mag_label_set = set(current_mag_labels)
        for label in mag_label_set:
            if ((label in orig_current_mag_labels) and (label not in new_mag_labels)):
                fill = [magnitude_fill_value] * new_length
                self.magnitudes[label][1] = np.append(self.magnitudes[label][1], fill)
            if ((label not in orig_current_mag_labels) and (label in new_mag_labels)):
                fill = [magnitude_fill_value] * current_length
                self.magnitudes[label] = [mag_sys, np.append(fill, catalog_to_add.magnitudes[label][1])]
            if ((label in orig_current_mag_labels) and (label in new_mag_labels)):
                self.magnitudes[label][1] = np.append(self.magnitudes[label][1],
                                                      catalog_to_add.magnitudes[label][1])

        # Update the catalog table if it already exists, so that it is consistent
        self.create_table()

    def add_magnitude_column(self, magnitude_list, magnitude_system='abmag', instrument='', filter_name=''):
        """Add a list of magnitudes to the catalog

        Parameters
        ----------
        some_stuff
        """
        # Force magnitude list to be a numpy array
        magnitude_list = np.array(magnitude_list)

        # Make sure instrument and filter are allowed values
        instrument = instrument.lower()
        filter_name = filter_name.lower()
        self.filter_check(instrument, filter_name)

        # Create the string for the column header
        if instrument == '' or filter_name == '':
            header = 'magnitude'
        else:
            if 'guider' in instrument.lower():
                header = '{}_magnitude'.format(instrument)
            else:
                header = '{}_{}_magnitude'.format(instrument, filter_name)

        # Get a list of magnitude_system for the existing catalog. No mixing of
        # magnitude systems is currently allowed.
        keys = self.magnitudes.keys()
        if len(keys) > 0:
            current_mag_sys = self.magnitudes[list(keys)[0]][0]
        else:
            current_mag_sys = magnitude_system

        # Make sure this is not a duplicate column
        if header in self.magnitudes.keys():
            raise ValueError(("WARNING: {} entry already exists in catalog. No duplicate entries allowed."
                              .format(header)))
        # Make sure the magnitude system is consistent
        elif magnitude_system != current_mag_sys:
                raise ValueError(("WARNING: only a single magnitude system is allowed per catalog. "
                                  "Current catalog is using {}, while the new entry is {}."
                                  .format(current_mag_sys, magnitude_system)))
        else:
            self.magnitudes[header] = [magnitude_system, magnitude_list]

        # Update the catalog table
        self.create_table()

    def filter_check(self, inst_name, filt_name):
        """Make sure the requested instrument/filter pair is valid
        """
        if inst_name == 'nircam':
            filter_file = 'nircam_filter_pupil_pairings.list'
        elif inst_name == 'niriss':
            filter_file = 'niriss_dual_wheel_list.txt'

        if inst_name in ['nircam', 'niriss']:
            filter_file_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../config/',
                                                            filter_file))
            filter_table = ascii.read(filter_file_path)
            if filt_name.upper() not in filter_table['filter']:
                raise ValueError("WARNING: {} is not a valid filter for {}.".format(filt_name,
                                                                                    inst_name.upper()))
        if inst_name == 'fgs':
            if filt_name not in ['', 'na']:
                raise ValueError("WARNING: {} is not a valid filter for FGS.".format(filt_name))

    def instrument_check(self, inst_name):
        """Only certain instrument values are allowed
        NOT USED - so that we can return J, H, K from searches
        """
        if inst_name not in ['nircam', 'niriss', 'fgs']:
            raise ValueError("WARNING: {} is not a valid instrument.".format(inst_name))

    def dec(self):
        """Return Dec values from catalog"""
        if self.location_units == 'position_RA_Dec':
            return self._dec
        else:
            return []

    def get_magnitudes(self, key):
        """Return the magnitude list (but not mag_sys) associated with the given key

        Parameters
        ----------
        justone
        """
        try:
            return self.magnitudes[key][1]
        except KeyError:
            print("WARNING: No {} magnitude column present.".format(key))

    def ra(self):
        """Return RA values from catalog"""
        if self.location_units == 'position_RA_Dec':
            return self._ra
        else:
            return []

    def x(self):
        """Return x values from catalog"""
        if self.location_units == 'position_pixels':
            return self._ra
        else:
            return []

    def y(self):
        """Return y values from catalog"""
        if self.location_units == 'position_pixels':
            return self._dec
        else:
            return []

    def create_table(self):
        """Create an astropy table containing the catalog
        """
        tab = create_basic_table(self._ra, self._dec, self.magnitudes, self.location_units)

        # Make sure there are at least 4 comment lines at the top
        self.table = pad_table_comments(tab)

    def save(self, output_name):
        """Write out the catalog to an ascii file

        Parameters
        ----------
        something
        """
        self.table.write(output_name, format='ascii', overwrite=True)


class GalaxyCatalog(PointSourceCatalog):
    def __init__(self, ra=[], dec=[], x=[], y=[], ellipticity=[], radius=[], sersic_index=[],
                 position_angle=[], radius_units='arcsec'):
        """Create instance of a galaxy catalog
        """
        # Add location information
        PointSourceCatalog.__init__(self, ra=ra, dec=dec, x=x, y=y)

        # Add galaxy-specific information
        self.morphology = {'radius': radius, 'ellipticity': ellipticity, 'sersic_index': sersic_index,
                           'pos_angle': position_angle}

        if radius_units not in ['arcsec', 'pixels']:
            raise ValueError("WARNING: Galaxy radii must be in units of 'arcsec' or 'pixels'.")
        self.radius_units = radius_units

    def create_table(self):
        """Create an astropy table containing the catalog
        """
        tab = create_basic_table(self._ra, self._dec, self.magnitudes, self.location_units)

        # Add morphology columns
        for key in self.morphology:
            values = self.morphology[key]
            col = Column(values, name=key)
            tab.add_column(col, index=3)

        # Add magnitude system and position units as metadata
        tab.meta['comments'].append("radius_{}".format(self.radius_units))

        # Make sure there are at least 4 comment lines at the top
        self.table = pad_table_comments(tab)


class ExtendedCatalog(PointSourceCatalog):
    def __init__(self, filenames=[], ra=[], dec=[], x=[], y=[], position_angle=[]):
        # Add location information
        PointSourceCatalog.__init__(self, ra=ra, dec=dec, x=x, y=y)

        # Add extended source-specific information
        self.filenames = filenames
        self.pos_angle = position_angle

    def create_table(self):
        """Create an astropy table containing the catalog
        """
        tab = create_basic_table(self._ra, self._dec, self.magnitudes, self.location_units)

        # Add the filename column
        file_col = Column(self.filenames, name='filename')
        tab.add_column(file_col, index=3)

        # Add the position angle column
        pa_col = Column(self.pos_angle, name='pos_angle')
        tab.add_column(pa_col)

        # Make sure there are at least 4 comment lines at the top
        self.table = pad_table_comments(tab)


class MovingPointSourceCatalog(PointSourceCatalog):
    def __init__(self, ra=[], dec=[], x=[], y=[], ra_velocity=[], dec_velocity=[], x_velocity=[],
                 y_velocity=[]):
        # Add location information
        PointSourceCatalog.__init__(self, ra=ra, dec=dec, x=x, y=y)

        if len(ra_velocity) > 0 and len(x_velocity) > 0:
            raise ValueError(("WARNING: Provide either RA and Dec velocities, or x and y "
                              "velocities, but not both."))
        # Object velocities
        if len(ra_velocity) > 0:
            self.ra_velocity = ra_velocity
            self.dec_velocity = dec_velocity
            self.velocity_units = 'velocity_RA_Dec'
        else:
            self.ra_velocity = x_velocity
            self.dec_velocity = y_velocity
            self.velocity_units = 'velocity_pixels'

    def create_table(self):
        tab = create_basic_velocity_table(self._ra, self._dec, self.magnitudes, self.location_units,
                                          self.ra_velocity, self.dec_velocity, self.velocity_units)

        # Make sure there are at least 4 comment lines at the top
        self.table = pad_table_comments(tab)


class MovingSersicCatalog(GalaxyCatalog, MovingPointSourceCatalog):
    def __init__(self, ra=[], dec=[], x=[], y=[], ra_velocity=[], dec_velocity=[], x_velocity=[],
                 y_velocity=[], ellipticity=[], radius=[], sersic_index=[], position_angle=[],
                 radius_units='arcsec'):
        # Add location information
        MovingPointSourceCatalog.__init__(self, ra=ra, dec=dec, x=x, y=y, ra_velocity=ra_velocity,
                                          dec_velocity=dec_velocity, x_velocity=x_velocity,
                                          y_velocity=y_velocity)
        GalaxyCatalog.__init__(self, ra=ra, dec=dec, x=x, y=y, ellipticity=ellipticity, radius=radius,
                               sersic_index=sersic_index, position_angle=position_angle,
                               radius_units='arcsec')

    def create_table(self):
        tab = create_basic_velocity_table(self._ra, self._dec, self.magnitudes, self.location_units,
                                          self.ra_velocity, self.dec_velocity, self.velocity_units)
        # Add morphology columns
        for key in self.morphology:
            values = self.morphology[key]
            col = Column(values, name=key)
            tab.add_column(col, index=3)

        # Add magnitude system and position units as metadata
        tab.meta['comments'].append("radius_{}".format(self.radius_units))

        # Make sure there are at least 4 comment lines at the top
        self.table = pad_table_comments(tab)


class MovingExtendedCatalog(ExtendedCatalog, MovingPointSourceCatalog):
    def __init__(self, ra=[], dec=[], x=[], y=[], ra_velocity=[], dec_velocity=[], x_velocity=[],
                 y_velocity=[], filenames=[], position_angle=[]):
        # Add location information
        MovingPointSourceCatalog.__init__(self, ra=ra, dec=dec, x=x, y=y, ra_velocity=ra_velocity,
                                          dec_velocity=dec_velocity, x_velocity=x_velocity,
                                          y_velocity=y_velocity)
        ExtendedCatalog.__init__(self, ra=ra, dec=dec, x=x, y=y, filenames=filenames,
                                 position_angle=position_angle)

    def create_table(self):
        tab = create_basic_velocity_table(self._ra, self._dec, self.magnitudes, self.location_units,
                                          self.ra_velocity, self.dec_velocity, self.velocity_units)
        # Add filename column
        file_col = Column(self.filenames, name='filename')
        tab.add_column(file_col, index=1)

        # Add position_angle column
        pa_col = Column(self.pos_angle, name='pos_angle')
        tab.add_column(pa_col, index=6)

        # Make sure there are at least 4 comment lines at the top
        self.table = pad_table_comments(tab)


class NonSiderealCatalog(MovingPointSourceCatalog):
    def __init__(self, ra=[], dec=[], x=[], y=[], ra_velocity=[], dec_velocity=[], x_velocity=[],
                 y_velocity=[], object_type=[]):
        # Add location information
        MovingPointSourceCatalog.__init__(self, ra=ra, dec=dec, x=x, y=y, ra_velocity=ra_velocity,
                                          dec_velocity=dec_velocity, x_velocity=x_velocity,
                                          y_velocity=y_velocity)

        # List of the type of sources in the non-sidereal catalog
        valid_objects = ['pointSource', 'galaxies', 'extended']
        for obj in object_type:
            if obj not in valid_objects:
                raise ValueError(("WARNING: object_type list members must be one of: {}"
                                  .format(obj, valid_objects)))
        self.object_type = object_type

    def create_table(self):
        """Create an astropy table containing the catalog
        """
        tab = create_basic_velocity_table(self._ra, self._dec, self.magnitudes, self.location_units,
                                          self.ra_velocity, self.dec_velocity, self.velocity_units)
        obj_col = Column(self.object_type, name='object')
        tab.add_column(obj_col, index=1)

        # Make sure there are at least 4 comment lines at the top
        self.table = pad_table_comments(tab)


def add_velocity_columns(input_table, ra_velocities, dec_velocities, velocity_units):
    """
    DO IT
    """
    ra_vel_col = Column(ra_velocities, name='x_or_RA_velocity')
    dec_vel_col = Column(dec_velocities, name='y_or_Dec_velocity')
    input_table.add_columns([ra_vel_col, dec_vel_col], indexes=[3, 3])
    input_table.meta['comments'].append(velocity_units)
    return input_table

def combine_catalogs(catalog1, catalog2):
    """Combine two catalogs into one
    input two catalog objects
    """


def create_basic_table(ra_values, dec_values, magnitudes, location_units):
    """Create astropy table containing the basic catalog info
    NOTE THAT THIS IS OUTSIDE CLASSES
    """
    basic_table = Table()

    # Add index, filename, RA, Dec or x, y columns
    index_col = Column(np.arange(1, len(ra_values)+1), name='index')
    ra_col = Column(ra_values, name='x_or_RA')
    dec_col = Column(dec_values, name='y_or_Dec')
    basic_table.add_columns([index_col, ra_col, dec_col])

    # Add magnitude columns
    for key in magnitudes:
        mag_values = magnitudes[key][1]
        mag_sys = magnitudes[key][0]
        mag_column = Column(mag_values, name=key)
        basic_table.add_column(mag_column)

    # Add magnitude system and position units as metadata
    basic_table.meta['comments'] = [location_units, mag_sys]
    return basic_table


def create_basic_velocity_table(ra_values, dec_values, magnitudes, location_units,
                                ra_velocities, dec_velocities, velocity_units):
    tab = create_basic_table(ra_values, dec_values, magnitudes, location_units)
    tab = add_velocity_columns(tab, ra_velocities, dec_velocities, velocity_units)
    return tab


def pad_table_comments(input_table):
    """do it"""
    if len(input_table.meta['comments']) < 4:
        pad = [''] * (4 - len(input_table.meta['comments']))
        input_table.meta['comments'].extend(pad)
    return input_table
