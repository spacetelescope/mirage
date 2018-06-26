#! /usr/bin/env python

'''
Functions for creating catalogs needed as inputs to
the nircam_simulator
'''

from astropy.io import ascii

class Catalog():
    def __init__(self):
        self.something = something

    def write_pointsource(self,ra,dec,mag,abmag=True,pixpos=False):
        stuff

    def write_galaxy(self,ra,dec,rad,ellip,pa,sersic,mag,abmag=True,pixpos=False):
        stuff

    def write_extended(self,something):
        stuff

    def write_moving_target_to_track(self,objtype,ra,dec,ravel,decvel,mag,abmag=True,pixpos=False,pixvel=False):
        stuff

    def write_moving_pointsource(self,ra,dec,mag,ravel,decvel,abmag=True,pixpos=False,pixvel=False):
        stuff

    def write_moving_galaxies(self,ra,dec,rad,ellip,pa,sersic,mag,ravel,decvel,abmag=True,pixpos=False,pixvel=False):
        stuff

    def write_moving_extended(self,file,ra,dec,mag,pa,ravel,decvel,abmag=True,pixpos=False,pixvel=False):
        stuff


    def add_radec_xy_note(self):
        return ("Position can be x,y or RA,Dec. If x,y, place "
                "the phrase 'position_pixels' in the top line of the file."
                "If RA,Dec, values can be floats or RA, Dec strings of the "
                "form RA: 'HH:MM:SS' and Dec: 'DD:MM:SS'")

    def add_velocity_note(self):
        return ("Velocity can be in units of pix/hour or arcsec/hour. "
                "If using pix/hour, place 'velocity_pixels' in the second "
                "line of the file. Note that if using velocities of pix/hour, "
                "the results will not be strictly correct because in reality "
                "distortion will cause objects' velocities to vary in pixels/hour. "
                "Velocities in arcsec/hour will be constant.")

    def add_mtt_sourcetype_note(self):
        return ("An object value containing 'point' will be interpreted "
                "as a point source. Anything containing 'sersic' will "
                "create a 2D sersic profile. Any other value will be "
                "interpreted as an extended source.")
        
    def add_mtt_velocity_note(self):
        return ("x_or_RA_velocity is the proper motion of the target in "
                "units of arcsec (or pixels) per year Y_or_Dec_velocity "
                "is the proper motion of the target in units of arcsec "
                "(or pixels) per year if the units are pixels per year, "
                "include 'velocity pixels' in line 2 above.")
    
    def add_gal_radius_note(self):
        return ("radius can be in units of pixels or arcseconds. Place "
                "'radius_pixels' at top of file to specify radii in pixels.")
