#! /usr/bin/env python

"""
This module contains utilities used for reading and writing various files
used or created by Mirage.

Author
------

    - Bryan Hilbert

Use
---

    This module can be imported and used as such:

    ::

        from mirage.utils import file_io
        gain, gain_head = file_io.read_gain_file(filename)
"""
import logging
import os

from astropy.io import ascii, fits
import astropy.units as q
import numpy as np

from mirage.logging import logging_functions
from mirage.utils.constants import LOG_CONFIG_FILENAME, STANDARD_LOGFILE_NAME


classdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../'))
log_config_file = os.path.join(classdir, 'logging', LOG_CONFIG_FILENAME)
logging_functions.create_logger(log_config_file, STANDARD_LOGFILE_NAME)


def read_file_spectrum(file, wave_units=q.um, flux_units=q.erg/q.s/q.cm**2/q.AA, survey=None):
    """Create a spectrum from an ASCII, XML, or FITS file

    Parameters
    ----------
    file: str
        The path to the file
    wave_units: str, astropy.units.quantity.Quantity
        The wavelength units
    flux_units: str, astropy.units.quantity.Quantity, None
        The flux units
    survey: str
        The name of the survey, ['SDSS']

    Returns
    -------
    list
        The [wavelength, flux, error] of the file spectrum with units
    """
    # Read the fits data...
    if file.endswith('.fits'):

        if file.endswith('.fits'):
            data, head = fits.getdata(file, header=True)

        elif survey == 'SDSS':
            raw, head = fits.getdata(file, header=True)
            flux_units = 1E-17 * q.erg / q.s / q.cm ** 2 / q.AA
            wave_units = q.AA
            log_w = head['COEFF0'] + head['COEFF1'] * np.arange(len(raw.flux))
            data = [10 ** log_w, raw.flux, raw.ivar]

        # Check if it is a recarray
        elif isinstance(raw, fits.fitsrec.FITS_rec):

            # Check if it's an SDSS spectrum
            raw = fits.getdata(file, ext=ext)
            data = raw['WAVELENGTH'], raw['FLUX'], raw['ERROR']

        # Otherwise just an array
        else:
            print("Sorry, I cannot read the file at", file)

    # ...or the ascii data...
    elif file.endswith('.txt'):
        data = np.genfromtxt(file, unpack=True)

    # ...or the VO Table
    elif file.endswith('.xml'):
        vot = vo.parse_single_table(file)
        data = np.array([list(i) for i in vot.array]).T

    else:
        raise IOError('The file needs to be ASCII, XML, or FITS.')

    # Make sure units are astropy quantities
    if isinstance(wave_units, str):
        wave_units = q.Unit(wave_units)
    if isinstance(flux_units, str):
        flux_units = q.Unit(flux_units)

    # Sanity check for wave_units
    if data[0].min() > 100 and wave_units == q.um:
        print("WARNING: Your wavelength range ({} - {}) looks like Angstroms. Are you sure it's {}?".format(
            data[0].min(), data[0].max(), wave_units))

    # Apply units
    wave = data[0] * wave_units
    flux = data[1] * (flux_units or 1.)
    if len(data) > 2:
        unc = data[2] * (flux_units or 1.)
    else:
        unc = None

    return [wave, flux, unc]


def read_filter_throughput(filename):
    """Read in the ascii file containing a filter throughput curve

    Parameters
    ----------
    filename : str
        Name of ascii file containing throughput info

    Returns
    -------
    (wavelengths, throughput) : tup
        Tuple of 1D numpy arrays containing the wavelength and throughput
        values from the file
    """
    tab = ascii.read(filename)
    return tab['Wavelength_microns'].data, tab['Throughput'].data


def read_gain_file(filename):
    """
    Read in CRDS-formatted gain reference file

    Paramters
    ---------
    filename : str
        Name of (fits) gain reference file

    Returns
    -------
    image : numpy.ndarray
        2D array gain map

    header : dict
        Information contained in the header of the file
    """
    logger = logging.getLogger('mirage.utils.file_io.read_gain_file')

    try:
        with fits.open(filename) as h:
            image = h[1].data
            header = h[0].header
    except (FileNotFoundError, IOError) as e:
        logger.error(e)

    mngain = np.nanmedian(image)

    # Set pixels with a gain value of 0 equal to mean
    image[image == 0] = mngain

    # Set any pixels with non-finite values equal to mean
    image[~np.isfinite(image)] = mngain
    return image, header


def readMTFile(filename):
    """
    Read in moving target list file

    Arguments:
    ----------
    filename : str
        name of moving target catalog file

    Returns:
    --------
    returns : obj
        Table containing moving target entries
        pixelflag (boolean) -- If true, locations are in units of
            pixels. If false, locations are RA, Dec
        pixelvelflag (boolean) -- If true, moving target velocities
            are in units of pixels/hour. If false, arcsec/hour
        magsys -- magnitude system of the moving target magnitudes
    """
    mtlist = ascii.read(filename, comment='#')

    # Convert all relevant columns to floats
    for col in mtlist.colnames:
        if mtlist[col].dtype in ['int64', 'int']:
            mtlist[col] = mtlist[col].data * 1.

    # Check to see whether the position is in x,y or ra,dec
    pixelflag = False
    try:
        if 'position_pixels' in mtlist.meta['comments'][0:4]:
            pixelflag = True
    except:
        pass

    # If present, check whether the velocity entries are pix/sec
    # or arcsec/sec.
    pixelvelflag = False
    try:
        if 'velocity_pixels' in mtlist.meta['comments'][0:4]:
            pixelvelflag = True
    except:
        pass

    # If present, check whether the radius entries (for galaxies)
    # are in arcsec or pixels. If in arcsec, change to pixels
    if 'radius' in mtlist.colnames:
        if 'radius_pixels' not in mtlist.meta['comments'][0:4]:
            mtlist['radius'] /= self.siaf.XSciScale

    # If galaxies are present, change position angle from degrees
    # to radians
    if 'pos_angle' in mtlist.colnames:
        mtlist['pos_angle'] = mtlist['pos_angle'] * np.pi / 180.

    # Check to see if magnitude system is specified in comments
    # If not, assume AB mags
    msys = 'abmag'

    condition = ('stmag' in mtlist.meta['comments'][0:4]) | ('vegamag' in mtlist.meta['comments'][0:4])
    if condition:
        msys = [l for l in mtlist.meta['comments'][0:4] if 'mag' in l][0]
        msys = msys.lower()

    return mtlist, pixelflag, pixelvelflag, msys.lower()
