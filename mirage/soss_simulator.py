#! /usr/bin/env python

"""
A module to generate simulated 2D time-series SOSS data

Authors: Joe Filippazzo, Kevin Volk, Nestor Espinoza, Jonathan Fraine, Michael Wolfe
"""

from copy import copy
from functools import partial, wraps
from multiprocessing.pool import ThreadPool
from multiprocessing import cpu_count
from pkg_resources import resource_filename
import logging
import time
import warnings
import os
import datetime
import yaml

import astropy.units as q
import astropy.constants as ac
from astropy.io import fits
from astropy.modeling.models import BlackBody, Voigt1D, Gaussian1D, Lorentz1D
import astropy.table as at
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord
from astroquery.simbad import Simbad
import batman
from bokeh.plotting import show
from hotsoss import utils as hu, plotting, locate_trace
import numpy as np
import synphot as sphot

from mirage.seed_image import save_seed, segmentation_map
from mirage.dark import dark_prep
from mirage.logging import logging_functions
from mirage.ramp_generator import obs_generator
from mirage.reference_files import crds_tools
from mirage.utils.constants import FLAMBDA_CGS_UNITS, LOG_CONFIG_FILENAME, STANDARD_LOGFILE_NAME
from mirage.utils import utils, file_io
from mirage.psf import soss_trace
from mirage.yaml import yaml_generator


classpath = os.path.dirname(__file__)
log_config_file = os.path.join(classpath, 'logging', LOG_CONFIG_FILENAME)
logging_functions.create_logger(log_config_file, STANDARD_LOGFILE_NAME)

warnings.simplefilter('ignore')

SUB_SLICE = {'SUBSTRIP96': slice(1802, 1898), 'SUBSTRIP256': slice(1792, 2048), 'FULL': slice(0, 2048)}
SUB_DIMS = {'SUBSTRIP96': (96, 2048), 'SUBSTRIP256': (256, 2048), 'FULL': (2048, 2048)}


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


class SossSim():
    """
    Generate NIRISS SOSS time series observations
    """
    def __init__(self, ngrps=2, nints=1, star=None, planet=None, tmodel=None, filter='CLEAR',
                 subarray='SUBSTRIP256', orders=[1, 2], paramfile=None, obs_date=None, target='New Target',
                 title=None, offline=True, test=False, override_dark=None, verbose=True):
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
        filter: str
            The name of the filter to use, ['CLEAR', 'F277W']
        subarray: str
            The name of the subarray to use, ['SUBSTRIP256', 'SUBSTRIP96', 'FULL']
        orders: int, list
            The orders to simulate, [1], [1, 2], [1, 2, 3]
        paramfile: str
            The path to a parameter file for the observation
        obs_date: str, datetime.datetime, astropy.time.Time
            The datetime of the start of the observation
        target: str (optional)
            The name of the target
        title: str (optionl)
            A title for the simulation
        offline: bool
            Offline mode
        test: bool
            Skip simulation
        verbose: bool
            Print status updates throughout calculation

        Example
        -------
        from mirage.soss_simulator import SossModelSim
        tso = SossModelSim(teff=3500, ngrps=2, nints=2)
        """
        # Metadata
        self.verbose = verbose
        self.target = target
        self.title = title or '{} Simulation'.format(self.target)
        self.offline = offline
        self.test = test
        self.params = None
        self.logger = logging.getLogger('mirage.soss_simulator')

        # Set default reference file parameters
        self._star = None
        self.ref_params = {"INSTRUME": "NIRISS", "READPATT": "NISRAPID", "EXP_TYPE": "NIS_SOSS", "DETECTOR": "NIS", "PUPIL": "GR700XD", "DATE-OBS": "2020-07-28", "TIME-OBS": "00:00:00", "INSTRUMENT": "NIRISS"}

        # Additional parmeters
        self.orders = orders
        self.groupgap = 0
        self.nframes = self.nsample = 1
        self.nresets = self.nresets1 = self.nresets2 = 1
        self.dropframes1 = self.dropframes3 = self.nskip = 0
        self.model_grid = 'ACES'
        self.override_dark = override_dark
        self.obs_datetime = obs_date or Time.now()
        self.ngrps = ngrps
        self.nints = nints
        self.filter = filter
        self.subarray = subarray
        self.readpatt = 'NISRAPID'
        self.paramfile = paramfile

        # Set instance attributes for the target
        self.lines = at.Table(names=('name', 'profile', 'x_0', 'amp', 'fwhm', 'flux'), dtype=('S20', 'S20', float, float, 'O', 'O'))
        self.star = star
        self.tmodel = tmodel
        self.ld_coeffs = np.zeros((3, 2048, 2))
        self.planet_radius = np.ones_like(self.avg_wave)
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
    def add_noise(self):
        """
        Run the signal data through the ramp generator
        """
        self.message("Starting noise generator...")
        start = time.time()

        # Generate segmentation map with correct dimensions
        self.segmap = segmentation_map.SegMap()
        self.segmap.ydim = self.nrows

        # Prepare dark current exposure if needed.
        if self.override_dark is None:
            self.logger.info('Running dark prep')
            d = dark_prep.DarkPrep(offline=self.offline)
            d.paramfile = self.paramfile
            d.prepare()
            use_darks = d.dark_files
        else:
            self.logger.info('\noverride_dark is set. Skipping call to dark_prep and using these files instead.')
            use_darks = self.override_dark

        # Make a seed files split like the dark file
        if isinstance(use_darks, str):
            use_darks = [use_darks]

        nint = 0
        nfiles = len(use_darks)
        self.tso = np.empty_like(self.tso_ideal)
        for n, dfile in enumerate(use_darks):
            dhead = fits.getheader(dfile)
            dnint = dhead['NINTS']

            # Save the raw signal as a seed image
            seed_seg = self.tso_ideal[nint:nint + dnint, :, :, :]
            seedfile, seedinfo = save_seed.save(seed_seg, self.paramfile, self.params, True, False, 1., 2048, (self.nrows, self.ncols), {'xoffset': 0, 'yoffset': 0}, 1, frametime=self.frame_time)

            # Combine into final observation
            self.logger.info('Running observation generator for segment {}/{}'.format(n, nfiles))
            obs = obs_generator.Observation(offline=self.offline)
            obs.linDark = dfile
            obs.seed = seed_seg
            obs.segmap = self.segmap
            obs.seedheader = seedinfo
            obs.paramfile = self.paramfile
            obs.create()

            if nfiles > 1:
                os.system('mv {} {}'.format(seedfile, seedfile.replace('_seed_image.fits', '_seg{}_part001_seed_image.fits'.format(str(n + 1).zfill(3)))))

            # Save ramp to tso attribute
            self.tso[nint:nint + dnint, :, :, :] = obs.raw_outramp

            # Pickup where last segment left off
            nint += dnint

        # Log it
        self.logger.info('\nSOSS simulator complete')
        self.message('Noise model finished: {} {}'.format(round(time.time() - start, 3), 's'))

    def create(self, n_jobs=-1, noise=True, override_dark=None, **kwargs):
        """
        Create the simulated 4D ramp data given the initialized TSO object

        Parameters
        ----------
        n_jobs: int
            The number of cores to use in multiprocessing
        tparams: dict
            The transit parameters of the system
        supersample_factor: int
            The

        Example
        -------
        # Run simulation of star only
        tso.create()

        # Simulate star with transiting exoplanet by including transmission spectrum and orbital params
        import batman
        from hotsoss import PLANET_DATA
        tparams = batman.TransitParams()
        tparams.t0 = 0.                                      # time of inferior conjunction
        tparams.per = 5.7214742                              # orbital period (days)
        tparams.a = 3.5                                      # semi-major axis (in units of stellar radii)
        tparams.inc = 89.8                                   # orbital inclination (in degrees)
        tparams.ecc = 0.                                     # eccentricity
        tparams.w = 90.                                      # longitude of periastron (in degrees)
        tparams.limb_dark = 'quadratic'                      # limb darkening profile to use
        tparams.u = [0.1, 0.1]                               # limb darkening coefficients
        tparams.rp = 1.                                     # planet radius (placeholder)

        tso.tmodel = tparams
        tso.planet = PLANET_DATA
        tso.create()
        """
        # Override the dark file
        self.override_dark = override_dark

        # Check that there is star data
        if self.star is None:
            print("No star to simulate! Please set the self.star attribute!")
            return

        # Check kwargs for updated attrs
        for key, val in kwargs.items():
            setattr(self, key, val)

        # Clear out old simulation
        self._reset_data()

        # Reset time parameters
        self._reset_time()

        # Reset relative response function
        self._reset_psfs()

        # Logging
        self.logger.info('\n\nRunning soss_simulator....\n')
        self.logger.info('Using parameter file: ')
        self.logger.info('{}'.format(self.paramfile))
        begin = time.time()

        # Make a simulation of all ones to save time
        if self.test:

            self.message("Generating test simulation of shape {}...".format(self.dims))
            self.tso_order1_ideal = self.tso_order2_ideal = np.ones((self.nints * self.ngrps, 256, 2048))

        # Make a true simulation
        else:
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

                # Get transit model for the current time chunk
                c_tmodel = None if self.tmodel is None else batman.TransitModel(self.tmodel, time_chunk.jd)

                # Generate simulation for each order
                for order in self.orders:

                    # Get the psf cube and filter response function
                    psfs = getattr(self, 'order{}_psfs'.format(order))

                    # Make a transit model for each radius and limb darkening coefficients
                    # (Dumb but multiprocessing requires it)
                    tmodels = [None] * self.ncols
                    if self.planet is not None:
                        for idx, (radius, ldc) in enumerate(zip(self.planet_radius[order - 1], self.ld_coeffs[order - 1])):
                            tmod = copy(c_tmodel)
                            tmod.rp = radius
                            tmod.u = ldc
                            tmodels[idx] = tmod
                            del tmod

                    # Generate the lightcurves at each wavelength
                    pool = ThreadPool(n_jobs)
                    func = partial(soss_trace.psf_lightcurve, time=time_chunk)
                    data = list(zip(psfs, tmodels))
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
                    del frames, lightcurves, psfs

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
        if noise:
            self.add_noise(**kwargs)
        else:
            self.tso = self.tso_ideal

        self.message('\nTotal time: {} {}'.format(round(time.time() - begin, 3), 's'))

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
        track = ['_ncols', '_nrows', '_nints', '_ngrps', '_nresets', '_subarray', '_filter', '_obs_datetime', '_orders', 'ld_profile', '_target', 'title', 'ra', 'dec']
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
        """Error when trying to change the number of columns"""
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
        """Error when trying to change the number of rows"""
        raise TypeError("The number of rows is fixed by setting the 'subarray' attribute.")

    @property
    def obs_datetime(self):
        """Getter for observation start date"""
        return self._obs_datetime

    @obs_datetime.setter
    def obs_datetime(self, obsdate):
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
            raise ValueError("'{}' not a supported obs_datetime. Try a dtype of {}".format(obsdate, good_dtypes))

        # Set the transit midpoint
        self._obs_datetime = Time(obsdate)

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
    def paramfile(self):
        """Getter for the paramfile"""
        return self._paramfile

    @paramfile.setter
    def paramfile(self, pfile):
        """Setter for the parameter file

        Parameters
        ----------
        pfile: str
            The path to the parameter file
        """
        # Populate params attribute if no file given
        if pfile is None:

            # Read template file
            pfile = resource_filename('mirage', 'tests/test_data/NIRISS/niriss_soss_substrip256_clear.yaml').replace('mirage/mirage', 'mirage')
            self._read_parameter_file(pfile)

            # Populate params dict
            self.params['Readout']['nint'] = self.nints
            self.params['Readout']['ngroup'] = self.ngrps
            self.params['Readout']['array_name'] = 'NIS_' + self.subarray.replace('FULL', 'SOSSFULL')
            self.params['Readout']['filter'] = self.filter
            self.params['Readout']['readpatt'] = self.readpatt
            self.params['Readout']['nframe'] = self.nframes
            self.params['Readout']['nskip'] = self.nskip
            self.params['Output']['target_name'] = self.target
            self.params['Output']['target_ra'] = self.ra
            self.params['Output']['target_dec'] = self.dec
            self.params['Output']['obs_id'] = self.title
            self.params['Output']['file'] = '{}_NIS_SOSS_{}_{}.fits'.format(self.title.replace(' ', '-'), self.filter, self.subarray)
            self.params['Output']['directory'] = '.'
            self.params['Output']['date_obs'], self.params['Output']['time_obs'] = self.obs_datetime.to_value('iso').split()
            self.params['Output']['paramfile'] = pfile = yaml_generator.yaml_from_params(self.params)

        else:

            # Check that it's a valid file
            if not os.path.isfile(pfile):
                raise IOError("{}: No such file.".format(pfile))

            # Set paramfile and read contents
            self._read_parameter_file(pfile)
            self.params['Output']['paramfile'] = pfile

            # Override observation parameters
            self._nints = self.params['Readout']['nint']
            self._ngrps = self.params['Readout']['ngroup']
            self._subarray = self.params['Readout']['array_name'].replace('NIS_', '').replace('SOSS', '')
            self._filter = self.params['Readout']['filter']
            self._obs_datetime = Time('{0[date_obs]} {0[time_obs]}'.format(self.params['Output']))
            self._readpatt = self.params['Readout']['readpatt']
            self._nrows, self._ncols = SUB_DIMS[self.subarray]
            self._target = self.params['Output']['target_name']
            self.ra = self.params['Output']['target_ra']
            self.dec = self.params['Output']['target_dec']
            self.title = self.params['Output']['obs_id']

        # Save paramfile
        self._paramfile = pfile

        # Only first order if F277W
        if self.filter == 'F277W':
            self.orders = [1]

        # Set the dependent quantities
        self.wave = hu.wave_solutions(self.subarray)
        self.avg_wave = np.mean(self.wave, axis=1)
        self.coeffs = locate_trace.trace_polynomial(subarray=self.subarray)

        # Reset data, time and psfs
        self._reset_data()
        self._reset_time()
        self._reset_psfs()

    @property
    def total_groups(self):
        """Getter for total groups"""
        return self.ngrps * self.nints

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
            self.planet_radius = np.ones_like(self.avg_wave)

        else:

            # Read file if path is provided
            if isinstance(spectrum, str):
                if os.path.exists(spectrum):
                    spectrum = file_io.read_file_spectrum(spectrum, flux_units=None)
                else:
                    raise IOError("{}: No file at that location.".format(spectrum))

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

            # Set the radius at the given wavelength from the transmission spectrum (Rp/R*)**2
            for order, wave in enumerate(self.avg_wave):
                tdepth = np.interp(wave, self.planet[0].to(q.um).value, self.planet[1])
                tdepth[tdepth < 0] = np.nan
                self.planet_radius[order] = np.sqrt(tdepth)

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

    def _read_parameter_file(self, pfile):
        """Read in the yaml parameter file"""
        try:
            with open(pfile, 'r') as infile:
                self.params = yaml.safe_load(infile)

        except FileNotFoundError:
            print("Unable to open {}".format(pfile))

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
        if all([i in self.info for i in ['subarray', 'nints', 'ngrps', 'obs_datetime', 'nresets']]):

            # Get frame time based on the subarray
            self.frame_time = utils.calc_frame_time('NIRISS', None, self.ncols, self.nrows, 1)
            self.group_time = (self.groupgap + self.nframes) * self.frame_time

            # The indexes of valid groups, skipping resets
            grp_idx = np.concatenate([np.arange(self.nresets, self.nresets + self.ngrps) + n * (self.nresets + self.ngrps) for n in range(self.nints)])

            # The time increments
            dt = TimeDelta(self.frame_time, format='sec')

            # Integration time of each frame in seconds
            self.inttime = np.tile((dt * grp_idx).sec[:self.ngrps], self.nints)

            # Datetime of each frame
            self.time = self.obs_datetime + (dt * grp_idx)

            # Exposure duration
            # self.duration = TimeDelta(self.time.max() - self.obs_datetime, format='sec')
            self.exposure_time = self.frame_time * ( (self.ngrps * self.nframes + (self.ngrps - 1) * self.groupgap + self.dropframes1) * self.nints)
            self.duration = self.exposure_time + self.frame_time * (self.dropframes3 * self.nints + self.nresets1 + self.nresets2 * (self.nints - 1))

            # Update params dict
            if self.params is not None:
                self.params['Output']['date_obs'], self.params['Output']['time_obs'] = self.obs_datetime.to_value('iso').split()

    def _reset_psfs(self):
        """Scale the psf for each detector column to the flux from the 1D spectrum"""
        # Check that all the appropriate values have been initialized
        if all([i in self.info for i in ['filter', 'subarray']]) and self.star is not None:

            for order in self.orders:

                # Get the wavelength map
                wave = self.avg_wave[order - 1]

                # Get relative spectral response for the order
                photom = fits.getdata(crds_tools.get_reffiles(self.ref_params, ['photom'])['photom'])
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

    @property
    def star(self):
        """Getter for the stellar data"""
        return self._star

    @star.setter
    def star(self, spectrum):
        """Setter for the stellar data

        Parameters
        ----------
        spectrum: sequence, str
            The [W, F] or [W, F, E] of the star to simulate
        """
        # Check if the star has been set
        if spectrum is None:
            self.message("No star to simulate! Please set the self.star attribute!")
            self._star = None

        else:

            # Read file if path is provided
            if isinstance(spectrum, str):
                if os.path.exists(spectrum):
                    spectrum = file_io.read_file_spectrum(spectrum)
                else:
                    raise IOError("{}: No file at that location.".format(spectrum))

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
                print("Wavelength range of input spectrum ({} - {}) does not cover the {} - {} range needed for a complete simulation. Interpolation will be used at the edges.".format(spec_min, spec_max, sim_min, sim_max))

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
        self.row_slice = SUB_SLICE[subarr]
        self._nrows, self._ncols = SUB_DIMS[subarr]

        # Set the dependent quantities
        self.wave = hu.wave_solutions(subarr)
        self.avg_wave = np.mean(self.wave, axis=1)
        self.coeffs = locate_trace.trace_polynomial(subarray=subarr)

        # Get correct reference files
        # self.ref_params['SUBARRAY'] = subarr

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
            self.create()


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
        bb = BlackBody(temperature=teff * q.K)
        wav = np.linspace(0.5, 2.9, 1000) * q.um
        flux = (bb(wav) * q.sr / bb.bolometric_flux.value).to(FLAMBDA_CGS_UNITS, q.spectral_density(wav)) * 1E-8 * scale

        # Initialize base class
        super().__init__(ngrps=ngrps, nints=nints, star=[wav, flux], subarray=subarray, filter=filter, **kwargs)

        # Add planet
        if add_planet:
            self.planet = hu.PLANET_DATA
            self.tmodel = hu.transit_params(self.time.jd)

        # Run the simulation
        if run:
            self.create()


class SossModelSim(SossSim):
    """Generate a SossSim object with a theoretical ATLAS or PHOENIX stellar spectrum of choice"""
    def __init__(self, ngrps=2, nints=2, teff=5700.0, logg=4.0, feh=0.0, alpha=0.0, jmag=9.0, models='ck04models', filter='CLEAR', subarray='SUBSTRIP256', run=True, add_planet=False, scale=1., **kwargs):
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
        # # Retrieve PHOENIX or ATLAS stellar models:
        # if stellar_model.lower() == 'phoenix':
        #     w, f = model_atmospheres.get_phoenix_model(feh, alpha, teff, logg)
        # elif stellar_model.lower() == 'atlas':
        #     w, f = model_atmospheres.get_atlas_model(feh, teff, logg)
        #
        # # Now scale model spectrum to user-input J-band:
        # f = model_atmospheres.scale_spectrum(w, f, jmag)

        # TODO: Retrieve stellar model (https://synphot.readthedocs.io/en/latest/#id13)
        # TODO: and trash the whole model_atmospheres.py module
        spectrum = sphot.icat(models, teff, feh, logg) # TODO: Where is icat function?
        bandpass = sphot.spectrum.SpectralElement.from_filter('2mass_j')
        obs = sphot.observation.Observation(spectrum, bandpass)

        sp_norm = spectrum.normalize(jmag, band=obs.bandpass)
        w = sp_norm.waveset
        f = sp_norm(w)

        # Initialize base class
        super().__init__(ngrps=ngrps, nints=nints, star=[w, f], subarray=subarray, filter=filter, **kwargs)

        # Add planet
        if add_planet:
            self.planet = hu.PLANET_DATA
            self.tmodel = hu.transit_params(self.time.jd)

        # Run the simulation
        if run:
            self.create()


class SossSeedSim(SossSim):
    """
    Generate a SossSim object from a 4D seed image
    """
    def __init__(self, seed, filter='CLEAR', paramfile=None, noise=True, **kwargs):
        """
        Parameters
        ----------
        seed: np.ndarray
            4D array of data
        """
        # Ingest seed file if possible
        if isinstance(seed, str):
            seed = fits.getdata(seed)

        # Get shape
        nints, ngrps, nrows, ncols = seed.shape

        # Determine subarray
        if nrows == 256 and ncols == 2048:
            subarray = 'SUBSTRIP256'
        elif nrows == 96 and ncols == 2048:
            subarray == 'SUBSTRIP96'
        elif nrows == 2048 and ncols == 2048:
            subarray == 'FULL'
        else:
            raise ValueError("{}: Axes 2 and 3 don't match a valid SOSS subarray. Try {}".format(seed.shape, SUB_DIMS))

        # Initialize base class
        super().__init__(ngrps=ngrps, nints=nints, star=None, subarray=subarray, filter=filter, paramfile=paramfile, **kwargs)

        # Set the ideal (noiseless) simulation
        self.tso_order1_ideal = seed
        self.tso_order2_ideal = np.zeros_like(seed)

        # Run noise and ramp generator
        if noise:
            self.add_noise()