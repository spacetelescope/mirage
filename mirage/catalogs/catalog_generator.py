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

TSO_GRISM_INDEX = 99999

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
            self._location_units = 'position_RA_Dec'
        if len(x) > 0:
            self._location_units = 'position_pixels'

    def add_catalog(self, catalog_to_add, magnitude_fill_value=99.):
        """Add a catalog to the current catalog instance"""
        # If the the source positions in the two catalogs have different units, then the catalogs
        # can't be combined.
        if self._location_units != catalog_to_add.location_units:
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
            if instrument.lower() in ['fgs1', 'fgs2']:
                header = '{}_magnitude'.format(instrument.lower())
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
            if filt_name.lower() not in ['guider1', 'guider2']:
                raise ValueError("WARNING: {} is not a valid filter for FGS.".format(filt_name))

    def instrument_check(self, inst_name):
        """Only certain instrument values are allowed
        NOT USED - so that we can return J, H, K from searches
        """
        if inst_name not in ['nircam', 'niriss', 'fgs']:
            raise ValueError("WARNING: {} is not a valid instrument.".format(inst_name))

    @property
    def dec(self):
        """Return Dec values from catalog"""
        if self._location_units == 'position_RA_Dec':
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

    @property
    def location_units(self):
        """Return RA values from catalog"""
        return self._location_units

    @property
    def ra(self):
        """Return RA values from catalog"""
        if self._location_units == 'position_RA_Dec':
            return self._ra
        else:
            return []

    @property
    def x(self):
        """Return x values from catalog"""
        if self._location_units == 'position_pixels':
            return self._ra
        else:
            return []

    @property
    def y(self):
        """Return y values from catalog"""
        if self._location_units == 'position_pixels':
            return self._dec
        else:
            return []

    def create_table(self):
        """Create an astropy table containing the catalog
        """
        tab = create_basic_table(self._ra, self._dec, self.magnitudes, self._location_units)

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
        self._morphology = {'radius': radius, 'ellipticity': ellipticity, 'sersic_index': sersic_index,
                            'pos_angle': position_angle}

        self._ellipticity = ellipticity
        self._pos_angle = position_angle
        self._radius = radius
        self._sersic_index = sersic_index

        if radius_units not in ['arcsec', 'pixels']:
            raise ValueError("WARNING: Galaxy radii must be in units of 'arcsec' or 'pixels'.")
        self._radius_units = radius_units

    def create_table(self):
        """Create an astropy table containing the catalog
        """
        tab = create_basic_table(self._ra, self._dec, self.magnitudes, self.location_units)

        # Add morphology columns
        for key in self._morphology:
            values = self._morphology[key]
            col = Column(values, name=key)
            tab.add_column(col, index=3)

        # Add magnitude system and position units as metadata
        tab.meta['comments'].append("radius_{}".format(self.radius_units))

        # Make sure there are at least 4 comment lines at the top
        self.table = pad_table_comments(tab)

    @property
    def ellipticity(self):
        """Return Dec values from catalog"""
        return self._ellipticity

    @property
    def position_angle(self):
        """Return Dec values from catalog"""
        return self._pos_angle

    @property
    def radius(self):
        """Return Dec values from catalog"""
        return self._radius

    @property
    def sersic_index(self):
        """Return Dec values from catalog"""
        return self._sersic_index

    @property
    def radius_units(self):
        """Return Dec values from catalog"""
        return self._radius_units

    @property
    def morphology(self):
        """Return Dec values from catalog"""
        return self._morphology


class ExtendedCatalog(PointSourceCatalog):
    def __init__(self, filenames=[], ra=[], dec=[], x=[], y=[], position_angle=[]):
        # Add location information
        PointSourceCatalog.__init__(self, ra=ra, dec=dec, x=x, y=y)

        # Add extended source-specific information
        self._filenames = filenames
        self._pos_angle = position_angle

    def create_table(self):
        """Create an astropy table containing the catalog
        """
        tab = create_basic_table(self._ra, self._dec, self.magnitudes, self._location_units)

        # Add the filename column
        file_col = Column(self._filenames, name='filename')
        tab.add_column(file_col, index=3)

        # Add the position angle column
        pa_col = Column(self._pos_angle, name='pos_angle')
        tab.add_column(pa_col)

        # Make sure there are at least 4 comment lines at the top
        self.table = pad_table_comments(tab)

    @property
    def filename(self):
        """Return Dec values from catalog"""
        return self._filenames

    @property
    def position_angle(self):
        """Return Dec values from catalog"""
        return self._pos_angle


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
            self._ra_velocity = ra_velocity
            self._dec_velocity = dec_velocity
            self._velocity_units = 'velocity_RA_Dec'
        else:
            self._ra_velocity = x_velocity
            self._dec_velocity = y_velocity
            self._velocity_units = 'velocity_pixels'

    def create_table(self):
        tab = create_basic_velocity_table(self._ra, self._dec, self.magnitudes, self._location_units,
                                          self._ra_velocity, self._dec_velocity, self._velocity_units)

        # Make sure there are at least 4 comment lines at the top
        self.table = pad_table_comments(tab)

    @property
    def dec_velocity(self):
        """Return Dec values from catalog"""
        return self._dec_velocity

    @property
    def ra_velocity(self):
        """Return Dec values from catalog"""
        return self._ra_velocity

    @property
    def velocity_units(self):
        """Return Dec values from catalog"""
        return self._velocity_units


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
        tab = create_basic_velocity_table(self._ra, self._dec, self.magnitudes, self._location_units,
                                          self._ra_velocity, self._dec_velocity, self._velocity_units)
        # Add morphology columns
        for key in self._morphology:
            values = self._morphology[key]
            col = Column(values, name=key)
            tab.add_column(col, index=3)

        # Add magnitude system and position units as metadata
        tab.meta['comments'].append("radius_{}".format(self._radius_units))

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
        tab = create_basic_velocity_table(self._ra, self._dec, self.magnitudes, self._location_units,
                                          self._ra_velocity, self._dec_velocity, self._velocity_units)
        # Add filename column
        file_col = Column(self._filenames, name='filename')
        tab.add_column(file_col, index=1)

        # Add position_angle column
        pa_col = Column(self._pos_angle, name='pos_angle')
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
        tab = create_basic_velocity_table(self._ra, self._dec, self.magnitudes, self._location_units,
                                          self._ra_velocity, self._dec_velocity, self._velocity_units)
        obj_col = Column(self.object_type, name='object')
        tab.add_column(obj_col, index=1)

        # Make sure there are at least 4 comment lines at the top
        self.table = pad_table_comments(tab)


class ImagingTSOCatalog(PointSourceCatalog):
    def __init__(self, ra=[], dec=[], x=[], y=[], lightcurve_file=[]):
        """Create instance of an imaging TSO catalog
        """
        # Add location information
        PointSourceCatalog.__init__(self, ra=ra, dec=dec, x=x, y=y)

        # Add TSO-specific information
        self._lightcurve_file = lightcurve_file

    def create_table(self):
        """Create an astropy table containing the catalog
        """
        tab = create_basic_table(self._ra, self._dec, self.magnitudes, self.location_units,
                                 minimum_index=TSO_GRISM_INDEX)

        # Add the lightcurve filename column
        lc_col = Column(self._lightcurve_file, name='lightcurve_file')
        tab.add_column(lc_col, index=3)

        # Make sure there are at least 4 comment lines at the top
        self.table = pad_table_comments(tab)

    @property
    def lightcurve_file(self):
        """Return lightcurve filenames"""
        return self._lightcurve_file


class GrismTSOCatalog(PointSourceCatalog):
    def __init__(self, ra=[], dec=[], x=[], y=[], semimajor_axis=[], orbital_inclination=[],
                 eccentricity=[], orbital_period=[], longitude_of_periastron=[], limb_dark_model=[],
                 limb_dark_coeffs=[], time_units=[], start_time=[], end_time=[], inferior_conj=[],
                 transmission_spectrum=[]):
        """Create instance of a grism TSO catalog
        """
        # Add location information
        PointSourceCatalog.__init__(self, ra=ra, dec=dec, x=x, y=y)

        # Add TSO-specific information
        self._semimajor_axis = np.array(semimajor_axis)
        self._orbital_inclination = np.array(orbital_inclination)
        self._eccentricity = np.array(eccentricity)
        self._orbital_period = np.array(orbital_period)
        self._longitude_of_periastron = np.array(longitude_of_periastron)
        self._limb_dark_model = np.array(limb_dark_model)
        self._time_units = np.array(time_units)
        self._start_time = np.array(start_time)
        self._end_time = np.array(end_time)
        self._inferior_conj = np.array(inferior_conj)
        self._transmission_spectrum = np.array(transmission_spectrum)

        # Limb darkening coefficients should be a string of comma-separated
        # numbers
        str_limb_dark_coeffs = []
        for entry in limb_dark_coeffs:
            str_limb_dark_coeffs.append(str(entry).replace('[', '').replace(']', ''))
        self._limb_dark_coeffs = np.array(str_limb_dark_coeffs)

    def create_table(self):
        """Create an astropy table containing the catalog
        """
        tab = create_basic_table(self._ra, self._dec, self.magnitudes, self.location_units,
                                 minimum_index=TSO_GRISM_INDEX)

        # Add the lightcurve filename column
        semimajor_col = Column(self._semimajor_axis, name='Semimajor_axis_in_stellar_radii')
        tab.add_column(semimajor_col, index=3)

        orb_inc_col = Column(self._orbital_inclination, name='Orbital_inclination_deg')
        tab.add_column(orb_inc_col, index=4)

        ecc_col = Column(self._eccentricity, name='Eccentricity')
        tab.add_column(ecc_col, index=5)

        pastron_col = Column(self._longitude_of_periastron, name='Longitude_of_periastron')
        tab.add_column(pastron_col, index=6)

        limb_model_col = Column(self._limb_dark_model, name='Limb_darkening_model')
        tab.add_column(limb_model_col, index=7)

        limb_coeff_col = Column(self._limb_dark_coeffs, name='Limb_darkening_coeffs')
        tab.add_column(limb_coeff_col, index=8)

        time_unit_col = Column(self._time_units, name='Time_units')
        tab.add_column(time_unit_col, index=9)

        start_col = Column(self._start_time, name='Start_time')
        tab.add_column(start_col, index=10)

        end_col = Column(self._end_time, name='End_time')
        tab.add_column(end_col, index=11)

        inf_conj_col = Column(self._inferior_conj, name='Time_of_inferior_conjunction')
        tab.add_column(inf_conj_col, index=12)

        orb_per_col = Column(self._orbital_period, name='Orbital_period')
        tab.add_column(orb_per_col, index=13)

        spec_col = Column(self._transmission_spectrum, name='Transmission_spectrum')
        tab.add_column(spec_col, index=14)

        # Make sure there are at least 4 comment lines at the top
        self.table = pad_table_comments(tab)

    @property
    def semimajor_axis(self):
        return self._semimajor_axis

    @property
    def orbital_inclination(self):
        return self._orbital_inclination

    @property
    def eccentricity(self):
        return self._eccentricity

    @property
    def longitude_of_periastron(self):
        return self._longitude_of_periastron

    @property
    def limb_dark_model(self):
        return self._limb_dark_model

    @property
    def limb_dark_coeffs(self):
        return self._limb_dark_coeffs

    @property
    def time_units(self):
        return self._time_units

    @property
    def start_time(self):
        return self._start_time

    @property
    def end_time(self):
        return self._end_time

    @property
    def inferior_conj(self):
        return self._inferior_conj

    @property
    def orbital_period(self):
        return self._orbital_period

    @property
    def transmission_spectrum(self):
        return self._transmission_spectrum


def add_velocity_columns(input_table, ra_velocities, dec_velocities, velocity_units):
    """
    DO IT
    """
    ra_vel_col = Column(ra_velocities, name='x_or_RA_velocity')
    dec_vel_col = Column(dec_velocities, name='y_or_Dec_velocity')
    input_table.add_columns([ra_vel_col, dec_vel_col], indexes=[3, 3])
    input_table.meta['comments'].append(velocity_units)
    return input_table


def cat_from_file(filename, catalog_type='point_source'):
    """Read in a Mirage-formatted ascii catalog file and place into an
    instance of a catalog object

    Parameters
    ----------
    filename : str
        Name of ascii catalog to be read in

    catalog_type : str
        Type of source catalog. Allowed values are:
        point_source, galaxy, extended, moving_point_source, moving_sersic,
        moving_extended, non_sidereal

    Returns
    -------
    catalog_object : mirage.catalogs.catalog_generator object
    """
    catalog_type = catalog_type.lower()
    allowed_types = ['point_source', 'galaxy', 'extended', 'moving_point_source',
                     'moving_sersic', 'moving_extended', 'non_sidereal']
    if catalog_type not in allowed_types:
        raise ValueError(("Input catalog type {} is not one of the allowed types: {}"
                         .format(catalog_type, allowed_types)))

    cat_table = ascii.read(filename)

    if 'position_pixels' in cat_table.meta['comments'][0:4]:
        xpos = 'x'
        ypos = 'y'
    else:
        xpos = 'ra'
        ypos = 'dec'

    magnitude_columns = [col for col in cat_table.colnames if 'magnitude' in col]

    # Point Source catalog
    if catalog_type == 'point_source':
        catalog_object = PointSourceCatalog(**{xpos: cat_table['x_or_RA']}, **{ypos: cat_table['y_or_Dec']})
        for magcol in magnitude_columns:
            instrument, filtername = get_inst_filter_from_colname(magcol)
            catalog_object.add_magnitude_column(cat_table[magcol], instrument=instrument,
                                                filter_name=filtername)

    # Galaxy catalog
    elif catalog_type == 'galaxy':
        radius_units = 'arcsec'
        if 'radius_pixels' in cat_table.meta['comments'][0:4]:
            radius_units = 'pixels'
        catalog_object = GalaxyCatalog(**{xpos: cat_table['x_or_RA']}, **{ypos: cat_table['y_or_Dec']},
                                       ellipticity=cat_table['ellipticity'], radius=cat_table['radius'],
                                       sersic_index=cat_table['sersic_index'], position_angle=cat_table['pos_angle'],
                                       radius_units=radius_units)
        for magcol in magnitude_columns:
            instrument, filtername = get_inst_filter_from_colname(magcol)
            catalog_object.add_magnitude_column(cat_table[magcol], instrument=instrument,
                                                filter_name=filtername)

    # Extended catalog
    elif catalog_type == 'extended':
        catalog_object = ExtendedCatalog(**{xpos: cat_table['x_or_RA']}, **{ypos: cat_table['y_or_Dec']},
                                         filenames=cat_table['filename'], position_angle=cat_table['pos_angle'])
        for magcol in magnitude_columns:
            instrument, filtername = get_inst_filter_from_colname(magcol)
            catalog_object.add_magnitude_column(cat_table[magcol], instrument=instrument,
                                                filter_name=filtername)

    # Moving point source catalog
    elif catalog_type == 'moving_point_source':
        velocity_units = 'velocity_RA_Dec'
        xvel = 'ra_velocity'
        yvel = 'dec_velocity'
        if 'velocity_pixels' in cat_table.meta['comments'][0:4]:
            velocity_units = 'velocity_pixels'
            xvel = 'x_velocity'
            yvel = 'y_velocity'
        catalog_object = MovingPointSourceCatalog(**{xpos: cat_table['x_or_RA']},
                                                  **{ypos: cat_table['y_or_Dec']},
                                                  **{xvel: cat_table['x_or_RA_velocity']},
                                                  **{yvel: cat_table['y_or_Dec_velocity']})
        for magcol in magnitude_columns:
            instrument, filtername = get_inst_filter_from_colname(magcol)
            catalog_object.add_magnitude_column(cat_table[magcol], instrument=instrument,
                                                filter_name=filtername)

    # Moving sersic catalog
    elif catalog_type == 'moving_sersic':
        velocity_units = 'velocity_RA_Dec'
        xvel = 'ra_velocity'
        yvel = 'dec_velocity'
        if 'velocity_pixels' in cat_table.meta['comments'][0:4]:
            velocity_units = 'velocity_pixels'
            xvel = 'x_velocity'
            yvel = 'y_velocity'
        radius_units = 'arcsec'
        if 'radius_pixels' in cat_table.meta['comments'][0:4]:
            radius_units = 'pixels'
        catalog_object = MovingSersicCatalog(**{xpos: cat_table['x_or_RA']},
                                             **{ypos: cat_table['y_or_Dec']},
                                             **{xvel: cat_table['x_or_RA_velocity']},
                                             **{yvel: cat_table['y_or_Dec_velocity']},
                                             ellipticity=cat_table['ellipticity'],
                                             radius=cat_table['radius'],
                                             sersic_index=cat_table['sersic_index'],
                                             position_angle=cat_table['pos_angle'],
                                             radius_units=radius_units)
        for magcol in magnitude_columns:
            instrument, filtername = get_inst_filter_from_colname(magcol)
            catalog_object.add_magnitude_column(cat_table[magcol], instrument=instrument,
                                                filter_name=filtername)

    # Moving extended catalog
    elif catalog_type == 'moving_extended':
        velocity_units = 'velocity_RA_Dec'
        xvel = 'ra_velocity'
        yvel = 'dec_velocity'
        if 'velocity_pixels' in cat_table.meta['comments'][0:4]:
            velocity_units = 'velocity_pixels'
            xvel = 'x_velocity'
            yvel = 'y_velocity'
        catalog_object = MovingExtendedCatalog(**{xpos: cat_table['x_or_RA']},
                                               **{ypos: cat_table['y_or_Dec']},
                                               **{xvel: cat_table['x_or_RA_velocity']},
                                               **{yvel: cat_table['y_or_Dec_velocity']},
                                               filenames=cat_table['filename'],
                                               position_angle=cat_table['pos_angle'])
        for magcol in magnitude_columns:
            instrument, filtername = get_inst_filter_from_colname(magcol)
            catalog_object.add_magnitude_column(cat_table[magcol], instrument=instrument,
                                                filter_name=filtername)

    # Non-sidereal catalog
    elif catalog_type == 'non_sidereal':
        velocity_units = 'velocity_RA_Dec'
        xvel = 'ra_velocity'
        yvel = 'dec_velocity'
        if 'velocity_pixels' in cat_table.meta['comments'][0:4]:
            velocity_units = 'velocity_pixels'
            xvel = 'x_velocity'
            yvel = 'y_velocity'
        catalog_object = NonSiderealCatalog(**{xpos: cat_table['x_or_RA']},
                                            **{ypos: cat_table['y_or_Dec']},
                                            **{xvel: cat_table['x_or_RA_velocity']},
                                            **{yvel: cat_table['y_or_Dec_velocity']},
                                            object_type=cat_table['object'])
        for magcol in magnitude_columns:
            instrument, filtername = get_inst_filter_from_colname(magcol)
            catalog_object.add_magnitude_column(cat_table[magcol], instrument=instrument,
                                                filter_name=filtername)
    return catalog_object


def create_basic_table(ra_values, dec_values, magnitudes, location_units, minimum_index=1):
    """Create astropy table containing the basic catalog info
    NOTE THAT THIS IS OUTSIDE CLASSES
    """
    basic_table = Table()

    # Add index, filename, RA, Dec or x, y columns
    index_col = Column(np.arange(minimum_index, minimum_index + len(ra_values)), name='index')
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


def get_inst_filter_from_colname(column_name):
    """Get the instrument and filter name from a magnitude column string

    Parameters
    ----------
    column_name : str
        Magnitude column header from catalog (e.g. 'nircam_f444w_magnitude')

    Returns
    -------
    instrument : str
        Instrument name

    filten_name : str
        Filter name
    """
    parts = column_name.split('_')
    instrument = parts[0].lower()
    if instrument == 'guider':
        filter_name = ''
    else:
        filter_name = parts[1].lower()
    return instrument, filter_name


def pad_table_comments(input_table):
    """do it"""
    if len(input_table.meta['comments']) < 4:
        pad = [''] * (4 - len(input_table.meta['comments']))
        input_table.meta['comments'].extend(pad)
    return input_table
