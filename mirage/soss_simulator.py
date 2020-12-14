#! /usr/bin/env python

"""
A module to generate simulated 2D time-series SOSS data

Authors: Joe Filippazzo, Kevin Volk, Nestor Espinoza, Jonathan Fraine, Michael Wolfe
"""

from copy import copy
import datetime
from functools import partial, wraps
from multiprocessing.pool import ThreadPool
from multiprocessing import cpu_count
from pkg_resources import resource_filename
import time
import warnings
import os
import datetime

import astropy.units as q
import astropy.constants as ac
from astropy.io import fits
from astropy.modeling.models import BlackBody1D, Voigt1D, Gaussian1D, Lorentz1D
from astropy.modeling.blackbody import FLAM
import astropy.table as at
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord
from astroquery.simbad import Simbad
import batman
# from bokeh.plotting import figure, show
# import bokeh.palettes as pal
from hotsoss import utils as hu, plotting, locate_trace
import numpy as np

from mirage import wfss_simulator
from mirage.catalogs import model_atmospheres
from mirage.seed_image import catalog_seed_image
from mirage.dark import dark_prep
from mirage.logging import logging_functions
from mirage.ramp_generator import obs_generator
from mirage.reference_files import crds_tools
from mirage.utils import read_fits
from mirage.utils.constants import CATALOG_YAML_ENTRIES, MEAN_GAIN_VALUES, LOG_CONFIG_FILENAME, STANDARD_LOGFILE_NAME
from mirage.utils.file_splitting import find_file_splits, SplitFileMetaData
from mirage.utils import utils, file_io, backgrounds
from mirage.psf import soss_trace
from mirage.utils.timer import Timer
from mirage.yaml import yaml_update


classpath = os.path.dirname(__file__)
log_config_file = os.path.join(classpath, 'logging', LOG_CONFIG_FILENAME)
logging_functions.create_logger(log_config_file, STANDARD_LOGFILE_NAME)

warnings.simplefilter('ignore')

SUB_SLICE = {'SUBSTRIP96': slice(1792, 1888), 'SUBSTRIP256': slice(1792, 2048), 'FULL': slice(0, 2048)}
SUB_DIMS = {'SUBSTRIP96': (96, 2048), 'SUBSTRIP256': (256, 2048), 'FULL': (2048, 2048)}


def check_psf_files():
    """Function to run on import to verify that the PSF files have been precomputed"""
    if not os.path.isfile(resource_filename('awesimsoss', 'files/SOSS_CLEAR_PSF_order1_1.npy')):
        print("Looks like you haven't generated the SOSS PSFs yet, which are required to produce simulations.")
        print("This takes about 10 minutes but you will only need to do it this one time.")
        compute = input("Would you like to do it now? [y] ")

        if compute is None or compute.lower() in ['y', 'yes']:
            multip = input("Do you want to run this using multiprocessing? [y]")
            if multip is None or multip.lower() in ['y', 'yes']:
                soss_trace.nuke_psfs(mprocessing=True)
            else:
                soss_trace.nuke_psfs(mprocessing=False)


def run_required(func):
    """A wrapper to check that the simulation has been run before a method can be executed"""
    @wraps(func)
    def _run_required(*args, **kwargs):
        """Check that the 'tso' attribute is not None"""
        if args[0].tso_order1_ideal is None:
            print("No simulation found! Please run the 'simulate' method first.")

        else:
            return func(*args, **kwargs)

    return _run_required


# TODO: Generate NIRISS SOSS PSFs the same way the other instruments do it,
#  with `mirage.psf.psf_selection.get_gridded_psf_library`. Right now it's just
#  pointing to awesimsoss/files directory
# TODO: Enable YAML paramfile as input in addition to args

check_psf_files()


class SossSim():
    """
    Generate NIRISS SOSS time series observations
    """
    def __init__(self, ngrps, nints, star=None, planet=None, tmodel=None, snr=700,
                 filter='CLEAR', subarray='SUBSTRIP256', orders=[1, 2],
                 obs_date=None, target='New Target', title=None, verbose=True):
        """
        Initialize the TSO object and do all pre-calculations

        Parameters
        ----------
        ngrps: int
            The number of groups per integration
        nints: int
            The number of integrations for the exposure
        star: sequence
            The wavelength and flux of the star
        planet: sequence
            The wavelength and transmission of the planet
        snr: float
            The signal-to-noise
        filter: str
            The name of the filter to use, ['CLEAR', 'F277W']
        subarray: str
            The name of the subarray to use, ['SUBSTRIP256', 'SUBSTRIP96', 'FULL']
        orders: int, list
            The orders to simulate, [1], [1, 2], [1, 2, 3]
        obs_date: str, datetime.datetime, astropy.time.Time
            The datetime of the start of the observation
        target: str (optional)
            The name of the target
        title: str (optionl)
            A title for the simulation
        verbose: bool
            Print status updates throughout calculation

        Example
        -------
        # Imports
        import numpy as np
        from awesimsoss import TSO, STAR_DATA
        import astropy.units as q
        from pkg_resources import resource_filename
        star = np.genfromtxt(resource_filename('awesimsoss', 'files/scaled_spectrum.txt'), unpack=True)
        star1D = [star[0]*q.um, (star[1]*q.W/q.m**2/q.um).to(q.erg/q.s/q.cm**2/q.AA)]

        # Initialize simulation
        tso = TSO(ngrps=3, nints=10, star=star1D)
        """
        # Metadata
        self.verbose = verbose
        self.target = target
        self.title = title or '{} Simulation'.format(self.target)

        # Set static values
        self._star = None

        # Set NISRAPID values from PRD Datamodes
        self.readpatt = 'NISRAPID'
        self.groupgap = 0
        self.nframes = self.nsample = 1
        self.nresets = self.nresets1 = self.nresets2 = 1
        self.dropframes1 = self.dropframes3 = 0

        # Set instance attributes for the exposure
        self.obs_date = obs_date or Time.now()
        self.ngrps = ngrps
        self.nints = nints
        self.total_groups = self.ngrps * self.nints

        self.orders = orders
        self.filter = filter
        self.header = ''
        self.snr = snr
        self.model_grid = 'ACES'
        self.subarray = subarray

        # Get correct reference files
        params = {"INSTRUME": "NIRISS",
                  "READPATT": "NIS",
                  "EXP_TYPE": "NIS_SOSS",
                  "DETECTOR": "NIS",
                  "PUPIL": "GR700XD",
                  "DATE-OBS": "2020-07-28",
                  "TIME-OBS": "00:00:00",
                  "INSTRUMENT": "NIRISS",
                  "FILTER": self.filter,
                  "SUBARRAY": self.subarray}
        self.refs = crds_tools.get_reffiles(params, ['photom'])

        # Set instance attributes for the target
        self.lines = at.Table(names=('name', 'profile', 'x_0', 'amp', 'fwhm', 'flux'), dtype=('S20', 'S20', float, float, 'O', 'O'))
        self.star = star
        self.tmodel = tmodel
        self.ld_coeffs = np.zeros((3, 2048, 2))
        self.planet = planet

    def add_line(self, x_0, amplitude, fwhm, profile='lorentz', name='Line I'):
        """
        Add an emission or absorption line to the spectrum

        Parameters
        ----------
        x_0: astropy.units.quantity.Quantity
            The rest wavelength of the line
        amplitude: astropy.units.quantity.Quantity
            The amplitude of the line relative to the continuum,
            with negative value for absorption and positive for emission
        fwhm: astropy.units.quantity.Quantity, sequence
            The full-width-half-max(s) of the line ('voigt' requires 2)
        profile: str
            The profile to use, ['voigt', 'lorentz', 'gaussian']
        name: str
            A name for the line
        """
        # Check the profile
        profiles = {'voigt': Voigt1D, 'gaussian': Gaussian1D, 'lorentz': Lorentz1D}
        if profile not in profiles:
            raise ValueError("'{}' profile not supported. Please select from {}".format(profile, list(profiles.keys())))

        # Select the profile
        prof = profiles[profile]

        # Convert to match star units and remove units
        x_0 = x_0.to(self.star[0].unit).value
        amplitude = amplitude.to(self.star[1].unit).value

        # Generate the line function
        if profile == 'voigt':
            if len(fwhm) != 2:
                raise TypeError("fwhm must be sequence of two values for Voigt profile.")
            else:
                fwhm_L, fwhm_G = [fw.to(self.star[0].unit).value for fw in fwhm]
            func = prof(amplitude_L=amplitude, x_0=x_0, fwhm_L=fwhm_L, fwhm_G=fwhm_G)

        elif profile == 'lorentz':
            fwhm = fwhm.to(self.star[0].unit).value
            func = prof(amplitude=amplitude, x_0=x_0, fwhm=fwhm)

        elif profile == 'gaussian':
            fwhm = fwhm.to(self.star[0].unit).value
            func = prof(amplitude=amplitude, mean=x_0, stddev=fwhm / 2.355)

        # Evaluate the profile
        line = func(self.star[0].value) * self.star[1].unit

        # Add the line to the line list
        self.lines.add_row([name, profile, x_0, amplitude, fwhm, line])

    @run_required
    def add_noise(self, c_pink=9.6, u_pink=3.2, bias_amp=5358.87, bias_offset=20944.06, acn=2.0, pca0_amp=0., rd_noise=12.95, pedestal_drift=18.3, gain=1.612, dc_seed=914075, noise_seed=2879328, zodi_scale=1., skip=[]):
        """
        Run the signal data through the ramp generator
        """
        self.message('Adding noise to TSO...')
        start = time.time()

        # Run signal through the ramp generator here:
        tso = self.tso_ideal.copy()

        # Put into 4D
        self.tso = tso.reshape(self.dims)

        self.message('Noise model finished: {} {}'.format(round(time.time() - start, 3), 's'))

    def export(self, outfile):
        """
        Export the simulated data to a JWST pipeline ingestible FITS file

        Parameters
        ----------
        outfile: str
            The path of the output file
        """
        if not outfile.endswith('_uncal.fits'):
            raise ValueError("Filename must end with '_uncal.fits'")

        # Make a RampModel
        data = copy(self.tso) if self.tso is not None else np.ones((1, 1, self.nrows, self.ncols))
        mod = ju.jwst_ramp_model(data=data, groupdq=np.zeros_like(data), pixeldq=np.zeros((self.nrows, self.ncols)), err=np.zeros_like(data))
        pix = hu.subarray_specs(self.subarray)

        # Set meta data values for header keywords
        mod.meta.telescope = 'JWST'
        mod.meta.instrument.name = 'NIRISS'
        mod.meta.instrument.detector = 'NIS'
        mod.meta.instrument.filter = self.filter
        mod.meta.instrument.pupil = 'GR700XD'
        mod.meta.exposure.type = 'NIS_SOSS'
        mod.meta.exposure.nints = self.nints
        mod.meta.exposure.ngroups = self.ngrps
        mod.meta.exposure.nframes = self.nframes
        mod.meta.exposure.readpatt = self.readpatt
        mod.meta.exposure.groupgap = self.groupgap
        mod.meta.exposure.frame_time = self.frame_time
        mod.meta.exposure.group_time = self.group_time
        mod.meta.exposure.exposure_time = self.exposure_time
        mod.meta.exposure.duration = self.duration
        mod.meta.exposure.nresets_at_start = self.nresets1
        mod.meta.exposure.nresets_between_ints = self.nresets2
        mod.meta.subarray.name = self.subarray
        mod.meta.subarray.xsize = data.shape[3]
        mod.meta.subarray.ysize = data.shape[2]
        mod.meta.subarray.xstart = pix.get('xloc', 1)
        mod.meta.subarray.ystart = pix.get('yloc', 1)
        mod.meta.subarray.fastaxis = -2
        mod.meta.subarray.slowaxis = -1
        mod.meta.observation.date = self.obs_date.iso.split()[0]
        mod.meta.observation.time = self.obs_date.iso.split()[1]
        mod.meta.target.ra = self.ra
        mod.meta.target.dec = self.dec
        mod.meta.target.source_type = 'POINT'

        # Save the file
        mod.save(outfile, overwrite=True)

        # Save input data
        with fits.open(outfile) as hdul:

            # Save input star data
            hdul.append(fits.ImageHDU(data=np.array([i.value for i in self.star], dtype=np.float64), name='STAR'))
            hdul['STAR'].header.set('FUNITS', str(self.star[1].unit))
            hdul['STAR'].header.set('WUNITS', str(self.star[0].unit))

            # Save input planet data
            if self.planet is not None:
                hdul.append(fits.ImageHDU(data=np.asarray(self.planet, dtype=np.float64), name='PLANET'))
                for param, val in self.tmodel.__dict__.items():
                    if isinstance(val, (float, int, str)):
                        hdul['PLANET'].header.set(param.upper()[:8], val)
                    elif isinstance(val, np.ndarray) and len(val) == 1:
                        hdul['PLANET'].header.set(param.upper(), val[0])
                    elif isinstance(val, type(None)):
                        hdul['PLANET'].header.set(param.upper(), '')
                    elif param == 'u':
                        for n, v in enumerate(val):
                            hdul['PLANET'].header.set('U{}'.format(n + 1), v)
                    else:
                        print(param, val, type(val))

            # Write to file
            hdul.writeto(outfile, overwrite=True)

        print('File saved as', outfile)

    @property
    def filter(self):
        """Getter for the filter"""
        return self._filter

    @filter.setter
    def filter(self, filt):
        """Setter for the filter

        Properties
        ----------
        filt: str
            The name of the filter to use,
            ['CLEAR', 'F277W']
        """
        # Valid filters
        filts = ['CLEAR', 'F277W']

        # Check the value
        if not isinstance(filt, str) or filt.upper() not in filts:
            raise ValueError("'{}' not a supported filter. Try {}".format(filt, filts))

        # Set it
        filt = filt.upper()
        self._filter = filt

        # If F277W, set orders to 1 to speed up calculation
        if filt == 'F277W':
            self.orders = [1]

        # Update the results
        self._reset_data()

    @property
    def info(self):
        """Summary table for the observation settings"""
        # Pull out relevant attributes
        track = ['_ncols', '_nrows', '_nints', '_ngrps', '_nresets', '_subarray', '_filter', '_obs_date', '_orders', 'ld_profile', '_target', 'title', 'ra', 'dec']
        settings = {key.strip('_'): val for key, val in self.__dict__.items() if key in track}
        return settings

    def message(self, message_text):
        """
        Print message

        Parameters
        ----------
        message_text: str
            The message to print
        """
        if self.verbose:
            print(message_text)

    @property
    def ncols(self):
        """Getter for the number of columns"""
        return self._ncols

    @ncols.setter
    def ncols(self, err):
        """Error when trying to change the number of columns
        """
        raise TypeError("The number of columns is fixed by setting the 'subarray' attribute.")

    @property
    def ngrps(self):
        """Getter for the number of groups"""
        return self._ngrps

    @ngrps.setter
    def ngrps(self, ngrp_val):
        """Setter for the number of groups

        Properties
        ----------
        ngrp_val: int
            The number of groups
        """
        # Check the value
        if not isinstance(ngrp_val, int):
            raise TypeError("The number of groups must be an integer")

        # Set it
        self._ngrps = ngrp_val

        # Update the results
        self._reset_data()
        self._reset_time()

    @property
    def nints(self):
        """Getter for the number of integrations"""
        return self._nints

    @nints.setter
    def nints(self, nint_val):
        """Setter for the number of integrations

        Properties
        ----------
        nint_val: int
            The number of integrations
        """
        # Check the value
        if not isinstance(nint_val, int):
            raise TypeError("The number of integrations must be an integer")

        # Set it
        self._nints = nint_val

        # Update the results
        self._reset_data()
        self._reset_time()

    def noise_report(self, plot=True, exclude=[]):
        """
        Table and plot of contributions from different noise sources

        Parameters
        ----------
        plot: bool
            Plot a figure of the noise contributions
        exclude: list
            A list of the sources to exclude from the plot
        """
        # Make empty table of inventory
        inv = at.Table()
        inv['nint'] = [int(i) for j in (np.ones((self.ngrps, self.nints)) * np.arange(1, self.nints + 1)).T for i in j]
        inv['ngrp'] = [int(i) for j in np.ones((self.nints, self.ngrps)) * np.arange(1, self.ngrps + 1) for i in j]

        # Add signal
        self.noise_model.noise_sources['signal'] = list(np.nanmean(self.tso_ideal, axis=(2, 3)))

        # Add the data
        cols = []
        for col, data in self.noise_model.noise_sources.items():
            inv[col] = [i for j in data for i in j]
            cols.append(col)

        # Print the table
        inv.pprint(max_width=-1, max_lines=-1)

        # Plot it
        if plot:
            fig = figure(x_axis_label='Group', y_axis_label='Electrons/ADU', width=1000, height=500)
            palette = pal.viridis(len(inv.colnames) - 2)
            x = np.arange(self.ngrps * self.nints) + 1
            for n, col in enumerate(cols):
                if col not in exclude:
                    fig.line(x, inv[col], color=palette[n], line_width=2, legend=col)

            # Legend
            fig.legend.click_policy = 'hide'

            show(fig)

    @property
    def nresets(self):
        """Getter for the number of resets"""
        return self._nresets

    @nresets.setter
    def nresets(self, nreset_val):
        """Setter for the number of resets

        Properties
        ----------
        nreset_val: int
            The number of resets
        """
        # Check the value
        if not isinstance(nreset_val, int) or nreset_val < 1:
            raise TypeError("The number of resets must be an integer greater than 0")

        # Set it
        self._nresets = nreset_val

        # Update the time (data shape doesn't change)
        self._reset_time()

    @property
    def nrows(self):
        """Getter for the number of rows"""
        return self._nrows

    @nrows.setter
    def nrows(self, err):
        """Error when trying to change the number of rows
        """
        raise TypeError("The number of rows is fixed by setting the 'subarray' attribute.")

    @property
    def obs_date(self):
        """Getter for observation start date"""
        return self._obs_date

    @obs_date.setter
    def obs_date(self, obsdate):
        """Setter for observation start date

        Properties
        ----------
        obsdate: str, datetime.datetime, astropy.time.Time
            The datetime of the start of the observation
        """
        # Acceptible time formats
        good_dtypes = str, datetime.datetime, Time

        # Try to make an astropy.time object
        if not isinstance(obsdate, good_dtypes):
            raise ValueError("'{}' not a supported obs_date. Try a dtype of {}".format(obsdate, good_dtypes))

        # Set the transit midpoint
        self._obs_date = Time(obsdate)

        # Reset the data and time arrays
        self._reset_data()
        self._reset_time()

    @property
    def orders(self):
        """Getter for the orders"""
        return self._orders

    @orders.setter
    def orders(self, ords):
        """Setter for the orders

        Properties
        ----------
        ords: list
            The orders to simulate, [1, 2, 3]
        """
        # Valid order lists
        orderlist = [[1], [1, 2], [1, 2, 3]]

        # Check the value and set single order to list
        if isinstance(ords, int):
            ords = [ords]
        if not all([o in [1, 2, 3] for o in ords]):
            raise ValueError("'{}' is not a valid list of orders. Try {}".format(ords, orderlist))

        # Set it
        self._orders = ords

        # Update the results
        self._reset_data()

    @property
    def planet(self):
        """Getter for the stellar data"""
        return self._planet

    @planet.setter
    def planet(self, spectrum):
        """Setter for the planetary data

        Parameters
        ----------
        spectrum: sequence
            The [W, F] or [W, F, E] of the planet to simulate
        """
        # Check if the planet has been set
        if spectrum is None:
            self._planet = None

        else:

            # Check planet is a sequence of length 2 or 3
            if not isinstance(spectrum, (list, tuple)) or not len(spectrum) in [2, 3]:
                raise ValueError(type(spectrum), ': Planet input must be a sequence of [W, F] or [W, F, E]')

            # Check the units
            if not spectrum[0].unit.is_equivalent(q.um):
                raise ValueError(spectrum[0].unit, ': Wavelength must be in units of distance')

            # Check the transmission spectrum is less than 1
            if not all(spectrum[1] < 1):
                raise ValueError('{} - {}: Transmission must be between 0 and 1'.format(min(spectrum[1]), max(spectrum[1])))

            # Check the wavelength range
            spec_min = np.nanmin(spectrum[0][spectrum[0] > 0.])
            spec_max = np.nanmax(spectrum[0][spectrum[0] > 0.])
            sim_min = np.nanmin(self.wave[self.wave > 0.]) * q.um
            sim_max = np.nanmax(self.wave[self.wave > 0.]) * q.um
            if spec_min > sim_min or spec_max < sim_max:
                print("Wavelength range of input spectrum ({} - {} um) does not cover the {} - {} um range needed for a complete simulation. Interpolation will be used at the edges.".format(spec_min, spec_max, sim_min, sim_max))

            # Good to go
            self._planet = spectrum

    @run_required
    def plot(self, idx=0, scale='linear', order=None, noise=True, traces=False, saturation=0.8, draw=True):
        """
        Plot a TSO frame

        Parameters
        ----------
        idx: int
            The frame index to plot
        scale: str
            Plot scale, ['linear', 'log']
        order: int (optional)
            The order to isolate
        noise: bool
            Plot with the noise model
        traces: bool
            Plot the traces used to generate the frame
        saturation: float
            The fraction of full well defined as saturation
        draw: bool
            Render the figure instead of returning it
        """
        # Get the data cube
        tso = self._select_data(order, noise)

        # Set the plot args
        wavecal = self.wave
        title = '{} - Frame {}'.format(self.title, idx)
        coeffs = locate_trace.trace_polynomial() if traces else None

        # Plot the frame
        fig = plotting.plot_frames(data=tso, idx=idx, scale=scale, trace_coeffs=coeffs, saturation=saturation, title=title, wavecal=wavecal)

        if draw:
            show(fig)
        else:
            return fig

    @run_required
    def plot_ramp(self, order=None, noise=True, draw=True):
        """
        Plot the total flux on each frame to display the ramp

        Parameters
        ----------
        order: sequence
            The order to isolate
        noise: bool
            Plot with the noise model
        draw: bool
            Render the figure instead of returning it
        """
        # Get the data cube
        tso = self._select_data(order, noise)

        # Make the figure
        fig = plotting.plot_ramp(tso)

        if draw:
            show(fig)
        else:
            return fig

    def _reset_data(self):
        """Reset the results to all zeros"""
        # Check that all the appropriate values have been initialized
        if all([i in self.info for i in ['nints', 'ngrps', 'nrows', 'ncols']]):

            # Update the dimensions
            self.dims = (self.nints, self.ngrps, self.nrows, self.ncols)
            self.dims3 = (self.nints * self.ngrps, self.nrows, self.ncols)

            # Reset the results
            for arr in ['tso'] + ['tso_order{}_ideal'.format(n) for n in self.orders]:
                setattr(self, arr, None)

    def _reset_time(self):
        """Reset the time axis based on the observation settings"""
        # Check that all the appropriate values have been initialized
        if all([i in self.info for i in ['subarray', 'nints', 'ngrps', 'obs_date', 'nresets']]):

            # Get frame time based on the subarray
            self.frame_time = self.subarray_specs.get('tfrm')
            self.group_time = (self.groupgap + self.nframes) * self.frame_time

            # The indexes of valid groups, skipping resets
            grp_idx = np.concatenate([np.arange(self.nresets, self.nresets + self.ngrps) + n * (self.nresets + self.ngrps) for n in range(self.nints)])

            # The time increments
            dt = TimeDelta(self.frame_time, format='sec')

            # Integration time of each frame in seconds
            self.inttime = np.tile((dt * grp_idx).sec[:self.ngrps], self.nints)

            # Datetime of each frame
            self.time = self.obs_date + (dt * grp_idx)

            # Exposure duration
            # self.duration = TimeDelta(self.time.max() - self.obs_date, format='sec')
            self.exposure_time = self.frame_time * ( (self.ngrps * self.nframes + (self.ngrps - 1) * self.groupgap + self.dropframes1) * self.nints)
            self.duration = self.exposure_time + self.frame_time * (self.dropframes3 * self.nints + self.nresets1 + self.nresets2 * (self.nints - 1))

    def _reset_psfs(self):
        """Scale the psf for each detector column to the flux from the 1D spectrum"""
        # Check that all the appropriate values have been initialized
        if all([i in self.info for i in ['filter', 'subarray']]) and self.star is not None:

            for order in self.orders:

                # Get the wavelength map
                wave = self.avg_wave[order - 1]

                # Get relative spectral response for the order
                photom = fits.getdata(self.refs['photom'])
                throughput = photom[(photom['order'] == order) & (photom['filter'] == self.filter) & (photom['pupil'] == 'GR700XD')]
                ph_wave = throughput.wavelength[throughput.wavelength > 0][1:-2]
                ph_resp = throughput.relresponse[throughput.wavelength > 0][1:-2]
                response = np.interp(wave, ph_wave, ph_resp)

                # Add spectral lines if necessary
                for line in self.lines:
                    self.star[1] += line['flux']

                # Convert response in [mJy/ADU/s] to [Flam/ADU/s] then invert so
                # that we can convert the flux at each wavelegth into [ADU/s]
                response = self.frame_time / (response * q.mJy * ac.c / (wave * q.um)**2).to(self.star[1].unit)
                flux = np.interp(wave, self.star[0].value, self.star[1].value, left=0, right=0) * self.star[1].unit * response
                cube = soss_trace.SOSS_psf_cube(filt=self.filter, order=order, subarray=self.subarray) * flux[:, None, None]
                setattr(self, 'order{}_response'.format(order), response)
                setattr(self, 'order{}_psfs'.format(order), cube)

    @run_required
    def _select_data(self, order, noise, reshape=True):
        """
        Select the data given the order and noise args

        Parameters
        ----------
        order: int (optional)
            The order to use, [1, 2, 3]
        noise: bool
            Include noise model
        reshape: bool
            Reshape to 3 dimensions

        Returns
        -------
        np.ndarray
            The selected data
        """
        if order in [1, 2]:
            tso = copy(getattr(self, 'tso_order{}_ideal'.format(order)))
        else:
            if noise:
                tso = copy(self.tso)
            else:
                tso = copy(self.tso_ideal)

        # Reshape data
        if reshape:
            tso.shape = self.dims3

        return tso

    def simulate(self, ld_coeffs=None, n_jobs=-1, params=None, supersample_factor=None, **kwargs):
        """
        Generate the simulated 4D ramp data given the initialized TSO object

        Parameters
        ----------
        ld_coeffs: array-like (optional)
            A 3D array that assigns limb darkening coefficients to each pixel, i.e. wavelength
        ld_profile: str (optional)
            The limb darkening profile to use
        n_jobs: int
            The number of cores to use in multiprocessing

        Example
        -------
        # Run simulation of star only
        tso.simulate()

        # Simulate star with transiting exoplanet by including transmission spectrum and orbital params
        import batman
        from hotsoss import PLANET_DATA
        params = batman.TransitParams()
        params.t0 = 0.                                      # time of inferior conjunction
        params.per = 5.7214742                              # orbital period (days)
        params.a = 3.5                                      # semi-major axis (in units of stellar radii)
        params.inc = 89.8                                   # orbital inclination (in degrees)
        params.ecc = 0.                                     # eccentricity
        params.w = 90.                                      # longitude of periastron (in degrees)
        params.limb_dark = 'quadratic'                      # limb darkening profile to use
        params.u = [0.1, 0.1]                               # limb darkening coefficients
        params.rp = 1.                                      # planet radius (placeholder)
        tmodel = batman.TransitModel(params, tso.time.jd)
        tmodel.teff = 3500                                  # effective temperature of the host star
        tmodel.logg = 5                                     # log surface gravity of the host star
        tmodel.feh = 0                                      # metallicity of the host star
        tso.simulate(planet=PLANET_DATA, tmodel=tmodel)
        """
        # Check that there is star data
        if self.star is None:
            print("No star to simulate! Please set the self.star attribute!")
            return

        # Check kwargs for updated attrs
        for key, val in kwargs.items():
            setattr(self, key, val)

        # Clear out old simulation
        self._reset_data()

        # Reset relative response function
        self._reset_psfs()

        # Get correct reference files
        # self.refs = ju.get_references(self.subarray, self.filter)

        # Start timer
        begin = time.time()

        # Set the number of cores for multiprocessing
        max_cores = cpu_count()
        if n_jobs == -1 or n_jobs > max_cores:
            n_jobs = max_cores

        # Chunk along the time axis so results can be dumped into a file and then deleted
        max_frames = 50
        nints_per_chunk = max_frames // self.ngrps
        nframes_per_chunk = self.ngrps * nints_per_chunk

        # Chunk the time arrays
        time_chunks = [self.time[i:i + nframes_per_chunk] for i in range(0, self.total_groups, nframes_per_chunk)]
        inttime_chunks = [self.inttime[i:i + nframes_per_chunk] for i in range(0, self.total_groups, nframes_per_chunk)]
        n_chunks = len(time_chunks)

        # Iterate over chunks
        self.message('Simulating {target} in {title}\nConfiguration: {subarray} + {filter}\nGroups: {ngrps}, Integrations: {nints}\n'.format(**self.info))
        for chunk, (time_chunk, inttime_chunk) in enumerate(zip(time_chunks, inttime_chunks)):

            # Run multiprocessing to generate lightcurves
            self.message('Constructing frames for chunk {}/{}...'.format(chunk + 1, n_chunks))
            start = time.time()

            # Re-define lightcurve model for the current chunk
            if params is not None:
                if supersample_factor is None:
                    c_tmodel = batman.TransitModel(params, time_chunk.jd)
                else:
                    frame_days = (self.frame_time * q.s).to(q.d).value
                    c_tmodel = batman.TransitModel(params, time_chunk.jd, supersample_factor=supersample_factor, exp_time=frame_days)
            else:
                c_tmodel = self.tmodel

            # Generate simulation for each order
            for order in self.orders:

                # Get the wavelength map
                wave = self.avg_wave[order - 1]

                # Get the psf cube and filter response function
                psfs = getattr(self, 'order{}_psfs'.format(order))

                # Get limb darkening coeffs and make into a list
                ld_coeffs = self.ld_coeffs[order - 1]
                ld_coeffs = list(map(list, ld_coeffs))

                # Set the radius at the given wavelength from the transmission
                # spectrum (Rp/R*)**2... or an array of ones
                if self.planet is not None:
                    tdepth = np.interp(wave, self.planet[0].to(q.um).value, self.planet[1])
                else:
                    tdepth = np.ones_like(wave)
                tdepth[tdepth < 0] = np.nan
                self.rp = np.sqrt(tdepth)

                # Generate the lightcurves at each wavelength
                pool = ThreadPool(n_jobs)
                func = partial(soss_trace.psf_lightcurve, time=time_chunk, tmodel=c_tmodel)
                data = list(zip(psfs, ld_coeffs, self.rp))
                lightcurves = np.asarray(pool.starmap(func, data), dtype=np.float16)
                pool.close()
                pool.join()
                del pool

                # Reshape to make frames
                lightcurves = lightcurves.swapaxes(0, 1)

                # Multiply by the integration time to convert to [ADU]
                lightcurves *= inttime_chunk[:, None, None, None]

                # Make the 2048*N lightcurves into N frames
                pool = ThreadPool(n_jobs)
                frames = np.asarray(pool.map(soss_trace.make_frame, lightcurves))
                pool.close()
                pool.join()
                del pool

                # Add it to the individual order
                order_name = 'tso_order{}_ideal'.format(order)
                if getattr(self, order_name) is None:
                    setattr(self, order_name, frames)
                else:
                    setattr(self, order_name, np.concatenate([getattr(self, order_name), frames]))

                # Clear memory
                del frames, lightcurves, psfs, wave

            self.message('Chunk {}/{} finished: {} {}'.format(chunk + 1, n_chunks, round(time.time() - start, 3), 's'))

        # Trim SUBSTRIP256 array if SUBSTRIP96
        if self.subarray == 'SUBSTRIP96':
            for order in self.orders:
                order_name = 'tso_order{}_ideal'.format(order)
                setattr(self, order_name, getattr(self, order_name)[:, :96, :])

        # Expand SUBSTRIP256 array if FULL frame
        if self.subarray == 'FULL':
            for order in self.orders:
                order_name = 'tso_order{}_ideal'.format(order)
                full = np.zeros((self.total_groups, 2048, 2048))
                full[:, -256:, :] = getattr(self, order_name)
                setattr(self, order_name, full)
                del full

        # Reshape into (nints, ngrps, y, x)
        for order in self.orders:
            order_name = 'tso_order{}_ideal'.format(order)
            setattr(self, order_name, getattr(self, order_name).reshape(self.dims).astype(np.float64))

        # Make ramps and add noise to the observations
        self.add_noise(**kwargs)

        # Simulate reference pixels
        # self.tso = ju.add_refpix(self.tso)

        self.message('\nTotal time: {} {}'.format(round(time.time() - begin, 3), 's'))

    @property
    def star(self):
        """Getter for the stellar data"""
        return self._star

    @star.setter
    def star(self, spectrum):
        """Setter for the stellar data

        Parameters
        ----------
        spectrum: sequence
            The [W, F] or [W, F, E] of the star to simulate
        """
        # Check if the star has been set
        if spectrum is None:
            self.message("No star to simulate! Please set the self.star attribute!")
            self._star = None

        else:

            # Check star is a sequence of length 2 or 3
            if not isinstance(spectrum, (list, tuple)) or not len(spectrum) in [2, 3]:
                raise ValueError(type(spectrum), ': Star input must be a sequence of [W, F] or [W, F, E]')

            # Check star has units
            if not all([isinstance(i, q.quantity.Quantity) for i in spectrum]):
                types = ', '.join([str(type(i)) for i in spectrum])
                raise ValueError('[{}]: Spectrum must be in astropy units'.format(types))

            # Check the units
            if not spectrum[0].unit.is_equivalent(q.um):
                raise ValueError(spectrum[0].unit, ': Wavelength must be in units of distance')

            if not all([i.unit.is_equivalent(q.erg / q.s / q.cm**2 / q.AA) for i in spectrum[1:]]):
                raise ValueError(spectrum[1].unit, ': Flux density must be in units of F_lambda')

            # Check the wavelength range
            spec_min = np.nanmin(spectrum[0][spectrum[0] > 0.])
            spec_max = np.nanmax(spectrum[0][spectrum[0] > 0.])
            sim_min = np.nanmin(self.wave[self.wave > 0.]) * q.um
            sim_max = np.nanmax(self.wave[self.wave > 0.]) * q.um
            if spec_min > sim_min or spec_max < sim_max:
                print("Wavelength range of input spectrum ({} - {} um) does not cover the {} - {} um range needed for a complete simulation. Interpolation will be used at the edges.".format(spec_min, spec_max, sim_min, sim_max))

            # Good to go
            self._star = spectrum

    @property
    def subarray(self):
        """Getter for the subarray"""
        return self._subarray

    @subarray.setter
    def subarray(self, subarr):
        """Setter for the subarray

        Properties
        ----------
        subarr: str
            The name of the subarray to use,
            ['SUBSTRIP256', 'SUBSTRIP96', 'FULL']
        """
        subs = ['SUBSTRIP256', 'SUBSTRIP96', 'FULL']

        # Check the value
        if subarr not in subs:
            raise ValueError("'{}' not a supported subarray. Try {}".format(subarr, subs))

        # Set the subarray
        self._subarray = subarr
        self.subarray_specs = hu.subarray_specs(subarr)
        self.row_slice = SUB_SLICE[subarr]

        # Set the dependent quantities
        self._ncols = 2048
        self._nrows = self.subarray_specs.get('y')
        self.wave = hu.wave_solutions(subarr)
        self.avg_wave = np.mean(self.wave, axis=1)
        self.coeffs = locate_trace.trace_polynomial(subarray=subarr)

        # Reset the data and time arrays
        self._reset_data()
        self._reset_time()

    @property
    def target(self):
        """Getter for target name"""
        return self._target

    @target.setter
    def target(self, name):
        """Setter for target name and coordinates

        Properties
        ----------
        tmid: str
            The transit midpoint
        """
        # Check the name
        if not isinstance(name, str):
            raise TypeError("Target name must be a string.")

        # Set the subarray
        self._target = name
        self.ra = 1.23456
        self.dec = 2.34567

        # Query Simbad for target RA and Dec
        if self.target != 'New Target':

            try:
                rec = Simbad.query_object(self.target)
                coords = SkyCoord(ra=rec[0]['RA'], dec=rec[0]['DEC'], unit=(q.hour, q.degree), frame='icrs')
                self.ra = coords.ra.degree
                self.dec = coords.dec.degree
                if self.verbose:
                    print("Coordinates {} {} for '{}' found in Simbad!".format(self.ra, self.dec, self.target))
            except TypeError:
                if self.verbose:
                    print("Could not resolve target '{}' in Simbad. Using ra={}, dec={}.".format(self.target, self.ra, self.dec))
                    print("Set coordinates manually by updating 'ra' and 'dec' attributes.")

    @property
    def tmodel(self):
        """Getter for the transit model"""
        return self._tmodel

    @tmodel.setter
    def tmodel(self, model, time_unit='days', model_grid='ACES'):
        """Setter for the transit model

        Parameters
        ----------
        model: batman.transitmodel.TransitModel
            The transit model
        time_unit: string
            The units of model.t, ['seconds', 'minutes', 'hours', 'days']
        """
        # Check if the transit model has been set
        if model is None:
            self._tmodel = None

        else:

            # Check transit model type
            mod_type = str(type(model))
            if not mod_type == "<class 'batman.transitmodel.TransitModel'>":
                raise TypeError("{}: Transit model must be of type batman.transitmodel.TransitModel".format(mod_type))

            # Check time units
            time_units = {'seconds': 86400., 'minutes': 1440., 'hours': 24., 'days': 1.}
            if time_unit not in time_units:
                raise ValueError("{}: time_unit must be {}".format(time_unit, time_units.keys()))

            # Convert seconds to days in order to match the Period and T0 parameters
            model.t /= time_units[time_unit]

            # Get the stellar parameters
            params = [model.teff, model.logg, model.feh]

            # Update the transit model
            self._tmodel = model

            # Update ld_coeffs
            self.ld_coeffs = [soss_trace.generate_SOSS_ldcs(self.avg_wave[order - 1], model.limb_dark, params, model_grid = self.model_grid) for order in self.orders]

    @property
    def tso_ideal(self):
        """Getter for TSO data without noise"""
        if self.tso_order1_ideal is None:
            return None

        if 2 in self.orders:
            return np.sum([self.tso_order1_ideal, self.tso_order2_ideal], axis=0)

        else:
            return self.tso_order1_ideal


class SossSpecSim(SossSim):
    """Generate a SossSim object from a default spectrum"""
    def __init__(self, ngrps=2, nints=2, filter='CLEAR', subarray='SUBSTRIP256', run=True, add_planet=False, **kwargs):
        """Get the test data and load the object

        Parameters
        ----------
        ngrps: int
            The number of groups per integration
        nints: int
            The number of integrations for the exposure
        filter: str
            The name of the filter to use, ['CLEAR', 'F277W']
        subarray: str
            The name of the subarray to use, ['SUBSTRIP256', 'SUBSTRIP96', 'FULL']
        run: bool
            Run the simulation after initialization
        add_planet: bool
            Add a transiting exoplanet
        """
        # Initialize base class
        super().__init__(ngrps=ngrps, nints=nints, star=hu.STAR_DATA, subarray=subarray, filter=filter, **kwargs)

        # Add planet
        if add_planet:
            self.planet = hu.PLANET_DATA
            self.tmodel = hu.transit_params(self.time.jd)

        # Run the simulation
        if run:
            self.simulate()


class SossBlackbodySim(SossSim):
    """Generate a SossSim object with a blackbody spectrum"""
    def __init__(self, ngrps=2, nints=2, teff=1800, filter='CLEAR', subarray='SUBSTRIP256', run=True, add_planet=False, scale=1., **kwargs):
        """Get the test data and load the object

        Parameters
        ---------
        ngrps: int
            The number of groups per integration
        nints: int
            The number of integrations for the exposure
        teff: int
            The effective temperature [K] of the test source
        filter: str
            The name of the filter to use, ['CLEAR', 'F277W']
        subarray: str
            The name of the subarray to use, ['SUBSTRIP256', 'SUBSTRIP96', 'FULL']
        run: bool
            Run the simulation after initialization
        add_planet: bool
            Add a transiting exoplanet
        scale: int, float
            Scale the flux by the given factor
        """
        # Generate a blackbody at the given temperature
        bb = BlackBody1D(temperature=teff * q.K)
        wav = np.linspace(0.5, 2.9, 1000) * q.um
        flux = bb(wav).to(FLAM, q.spectral_density(wav)) * 1E-8 * scale

        # Initialize base class
        super().__init__(ngrps=ngrps, nints=nints, star=[wav, flux], subarray=subarray, filter=filter, **kwargs)

        # Add planet
        if add_planet:
            self.planet = hu.PLANET_DATA
            self.tmodel = hu.transit_params(self.time.jd)

        # Run the simulation
        if run:
            self.simulate()


class SossModelSim(SossSim):
    """Generate a SossSim object with a theoretical ATLAS or PHOENIX stellar spectrum of choice"""
    def __init__(self, ngrps=2, nints=2, teff=5700.0, logg=4.0, feh=0.0, alpha=0.0, jmag=9.0, stellar_model='ATLAS', filter='CLEAR', subarray='SUBSTRIP256', run=True, add_planet=False, scale=1., **kwargs):
        """Get the test data and load the object

        Parameters
        ---------
        ngrps: int
            The number of groups per integration
        nints: int
            The number of integrations for the exposure
        teff: double
            The effective temperature [K] of the stellar source
        logg: double
            The log-gravity of the stellar source
        feh: double
            The [Fe/H] of the stellar source
        alpha: double
            The alpha enhancement of the stellar source
        jmag: double
            The J magnitude of the source. This will be used to scale the model stellar flux to Earth-values.
        stellar_model: str
            The stellar model grid to use. Can either be 'ATLAS' or 'PHOENIX'. Default is 'ATLAS'
        filter: str
            The name of the filter to use, ['CLEAR', 'F277W']
        subarray: str
            The name of the subarray to use, ['SUBSTRIP256', 'SUBSTRIP96', 'FULL']
        run: bool
            Run the simulation after initialization
        add_planet: bool
            Add a transiting exoplanet
        scale: int, float
            Scale the flux by the given factor
        """
        # Retrieve PHOENIX or ATLAS stellar models:
        if stellar_model.lower() == 'phoenix':
            w, f = model_atmospheres.get_phoenix_model(feh, alpha, teff, logg)
        elif stellar_model.lower() == 'atlas':
            w, f = model_atmospheres.get_atlas_model(feh, teff, logg)

        # Now scale model spectrum to user-input J-band:
        f = model_atmospheres.scale_spectrum(w, f, jmag)

        # Initialize base class
        super().__init__(ngrps=ngrps, nints=nints, star=[w, f], subarray=subarray, filter=filter, **kwargs)

        # Add planet
        if add_planet:
            self.planet = hu.PLANET_DATA
            self.tmodel = hu.transit_params(self.time.jd)

        # Run the simulation
        if run:
            self.simulate()
