#! /usr/bin/env python

"""This module contains functions for creating Mirage-formatted source
catalogs"""

import os

from astropy.io import ascii
from astropy.table import Table, Column


class PointSourceCatalog():
    def __init__(self, ra=[], dec=[], x=[], y=[]):
        """Initialize the point source catalog. Users can enter lists of RA and Dec values
        or x and y values for source potisions

        Parameters
        ----------
        something
        """
        if len(ra) != len(dec):
            raise ValueError(("WARNING: inconsistent length between RA and Dec inputs."))
        if len(x) != len(y):
            raise ValueError(("WARNING: inconsistent length between x and y inputs."))
        if len(ra) > 0 and len(x) > 0:
            raise ValueError(("WARNING: Provide either both RA and Dec values, or both x and y values."))

        self.ra = ra
        self.dec = dec
        self.x = x
        self.y = y
        self.magnitudes = {}

        # Determine the units for the location fields. All that Mirage needs to know is whether
        # the units are pixels or not. Degrees vs hour angle is not important at this point.
        if len(ra) > 0:
            self.location_untis = 'degrees'
        if len(x) > 0:
            self.location_units = 'pixels'

    def add_magnitude_column(self, magnitude_list, magnitude_system='abmag', instrument='', filter_name=''):
        """Add a list of magnitudes to the catalog

        Parameters
        ----------
        some_stuff
        """
        # Make sure instrument and filter are allowed values
        instrument = instrument.lower()
        filter_name = filter_name.lower()
        #self.instrument_check(instrument)
        self.filter_check(instrument, filter_name)

        # Create the string for the column header
        if instrument == '' or filter_name == '':
            header = 'magnitude'
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
            self.magnitudes[header] = (magnitude_system, magnitude_list)

    def filter_check(self, inst_name, filt_name):
        """Make sure the requested instrument/filter pair is valid
        """
        if inst_name == 'nircam':
            filter_file = 'nircam_filter_pupil_pairings.list'
        elif inst_name == 'niriss':
            filter_file = 'niriss_dual_wheel_list.txt'
        filter_file_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../config/', filter_file))

        if inst_name in ['nircam', 'niriss']:
            filter_table = ascii.read(filter_file_path)
            if filt_name.upper() not in filter_table['filter']:
                raise ValueError("WARNING: {} is not a valid filter for {}.".format(filt_name,
                                                                                    inst_name.upper()))
        if inst_name == 'fgs':
            if filt_name not in ['', 'na']:
                raise ValueError("WARNING: {} is not a valid filter for FGS.".format(filt_name))

    def instrument_check(self, inst_name):
        """Only certain instrument values are allowed
        """
        if inst_name not in ['nircam', 'niriss', 'fgs']:
            raise ValueError("WARNING: {} is not a valid instrument.".format(inst_name))

    def get_dec(self):
        """Return Dec values from catalog"""
        return self.dec

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

    def get_ra(self):
        """Return RA values from catalog"""
        return self.ra

    def get_x(self):
        """Return x values from catalog"""
        return self.x

    def get_y(self):
        """Return y values from catalog"""
        return self.y

    def create_table(self):
        """Create an astropy table containing the catalog
        """
        tab = Table()
        convert self.magnitudes[list(keys)[0]][0] and self.location_units to metadata
        for key in self.magnitudes:
            mag_values = self.magnitudes[key][1]
            mag_column = Column(mag_values, name=key)
            tab.add_column(mag_column)

    def save(self, catalog_table, output_name):
        """Write out the catalog to an ascii file

        Parameters
        ----------
        something
        """
        ascii.write(catalog_table, output_name, overwrite=True)


class GalaxyCatalog(PointSourceCatalog):
    def __init__(self, ra=[], dec=[], x=[], y=[], ellipticity=[], radius=[], sersic_index=[]):
        """Create instance of a galaxy catalog
        """
        # Add location information
        PointSourceCatalog.__init__(self, ra=ra, dec=dec, x=x, y=y)

        # Add galaxy-specific information
        self.radius = radius
        self.ellipticity = ellipticity
        self.sersic = sersic_index
