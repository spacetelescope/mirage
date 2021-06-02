"""
A module to generate simulated 2D time-series SOSS data

Authors: Joe Filippazzo
"""

import os
from pkg_resources import resource_filename
import multiprocessing
import time
from functools import partial
import warnings

import numpy as np
from astropy.io import fits
from bokeh.plotting import figure, show
from hotsoss import utils, locate_trace
from scipy.interpolate import interp1d
from scipy.ndimage.interpolation import rotate
from scipy.interpolate import interp2d, RectBivariateSpline

warnings.simplefilter('ignore')

# Set the MIRAGE_DATA environment variable
if os.environ.get('MIRAGE_DATA') is None:
    PSF_DIR = None
else:
    PSF_DIR = os.path.join(os.environ['MIRAGE_DATA'], 'niriss/soss_psfs/')


def calculate_psf_tilts():
    """
    Calculate the tilt of the psf at the center of each column
    using all binned pixels in the given wavelength calibration file
    for both orders and save to file
    """
    for order in [1, 2]:

        # Get the file
        psf_file = os.path.join(PSF_DIR, 'SOSS_PSF_tilt_order{}.npy'.format(order))

        # Dimensions
        subarray = 'SUBSTRIP256'
        X = range(2048)
        Y = range(256)

        # Get the wave map
        wave_map = utils.wave_solutions(subarray, order).astype(float)

        # Get the y-coordinate of the trace polynomial in this column
        # (center of the trace)
        coeffs = locate_trace.trace_polynomial(subarray=subarray, order=order)
        trace = np.polyval(coeffs, X)

        # Interpolate to get the wavelength value at the center
        wave = interp2d(X, Y, wave_map)

        # Get the wavelength of the trace center in each column
        trace_wave = []
        for x, y in zip(X, trace):
            trace_wave.append(wave(x, y)[0])

        # For each column wavelength (defined by the wavelength at
        # the trace center) define an isowavelength contour
        angles = []
        for n, x in enumerate(X):

            w = trace_wave[x]

            # Edge cases
            try:
                w0 = trace_wave[x-1]
            except IndexError:
                w0 = 0

            try:
                w1 = trace_wave[x+1]
            except IndexError:
                w1 = 10

            # Define the width of the wavelength bin as half-way
            # between neighboring points
            dw0 = np.mean([w0, w])
            dw1 = np.mean([w1, w])

            # Get the coordinates of all the pixels in that range
            yy, xx = np.where(np.logical_and(wave_map >= dw0, wave_map < dw1))

            # Find the angle between the vertical and the tilted wavelength bin
            if len(xx) >= 1:
                angle = get_angle([xx[-1], yy[-1]], [x, trace[x]])
            else:
                angle = 0

            # Don't flip them upside down
            angle = angle % 180

            # Add to the array
            angles.append(angle)

        # Save the file
        np.save(psf_file, np.array(angles))
        print('Angles saved to', psf_file)


def nuke_psfs(tilts=True, raw=True, final=True, mprocessing=True):
    """Generate all the psf cubes from scratch"""
    # Calculate the psf tilts
    if tilts:
        calculate_psf_tilts()

    for filt in ['CLEAR', 'F277W']:

        # Calculate the raw psfs from WebbPSF
        if raw:
            generate_SOSS_psfs(filt)

        # Generate the rotated and interpolated psfs ready for trace assembly
        if final:
            SOSS_psf_cube(filt=filt, generate=True, mprocessing=mprocessing)


def generate_SOSS_ldcs(wavelengths, ld_profile, params, model_grid='ACES', subarray='SUBSTRIP256', n_bins=100):
    """
    Generate a lookup table of limb darkening coefficients for full
    SOSS wavelength range

    Parameters
    ----------
    wavelengths: sequence
        The wavelengths at which to calculate the LDCs
    ld_profile: str
        A limb darkening profile name supported by
        `ExoCTK.ldc.ldcfit.ld_profile()`
    params: sequence
        The stellar parameters [Teff, logg, FeH]
    model_grid: modelgrid.ModelGrid
        The grid of stellar intensity models to calculate LDCs from
    subarray: str
        The name of the subarray to use, ['SUBSTRIP96', 'SUBSTRIP256', 'FULL']
    n_bins: int
        The number of bins to break up the grism into

    Returns
    -------
    np.ndarray
        An array of limb darkening coefficients for each wavelength

    Example
    -------
    from mirage.psf import soss_trace as st
    lookup = st.generate_SOSS_ldcs(np.linspace(1., 2., 3), 'quadratic', [3300, 4.5, 0])
    """
    try:

        from exoctk import modelgrid
        from exoctk.limb_darkening import limb_darkening_fit as lf
        from svo_filters import svo

        # Break the bandpass up into n_bins pieces
        bandpass = svo.Filter('NIRISS.GR700XD.1', n_bins=n_bins, verbose=False)

        # Calculate the LDCs
        ldcs = lf.LDC(model_grid=model_grid)
        ldcs.calculate(params[0], params[1], params[2], ld_profile, mu_min=0.08, bandpass=bandpass, verbose=False)

        # Interpolate the LDCs to the desired wavelengths
        # TODO: Propagate errors
        coeff_cols = [col for col in ldcs.results.colnames if col.startswith('c') and len(col) == 2]
        coeff_errs = [err for err in ldcs.results.colnames if err.startswith('e') and len(err) == 2]
        coeffs = [[np.interp(wav, list(ldcs.results['wave_eff']), list(ldcs.results[c])) for c in coeff_cols] for wav in wavelengths]
        coeffs = np.array(coeffs)

        del ldcs

    except Exception as exc:

        print(exc)
        print('There was a problem computing those limb darkening coefficients. Using all zeros.')

        n_coeffs = 1 if ld_profile in ['uniform', 'linear'] else 3 if ld_profile == '3-parameter' else 4 if ld_profile == '4-parameter' else 2
        coeffs = np.zeros((len(wavelengths), n_coeffs))

    return coeffs


def generate_SOSS_psfs(filt):
    """
    Gnerate a cube of the psf at 100 wavelengths from the min to the max wavelength

    Parameters
    ----------
    filt: str
        The filter to use, ['CLEAR', 'F277W']
    """
    try:

        import webbpsf

        # Get the file
        file = os.path.join(PSF_DIR, 'SOSS_{}_PSF.fits'.format(filt))

        # Get the NIRISS class from webbpsf and set the filter
        ns = webbpsf.NIRISS()
        ns.filter = filt
        ns.pupil_mask = 'GR700XD'

        # Get the min and max wavelengths
        wavelengths = utils.wave_solutions('SUBSTRIP256').flatten()
        wave_min = np.max([ns.SHORT_WAVELENGTH_MIN * 1E6, np.min(wavelengths[wavelengths > 0])])
        wave_max = np.min([ns.LONG_WAVELENGTH_MAX * 1E6, np.max(wavelengths[wavelengths > 0])])

        # webbpsf.calc_datacube can only handle 100 but that's sufficient
        W = np.linspace(wave_min, wave_max, 100) * 1E-6

        # Calculate the psfs
        print("Generating SOSS psfs. This takes about 8 minutes...")
        start = time.time()
        PSF = ns.calc_datacube(W, oversample=1)[0].data
        print("Finished in", time.time()-start)

        # Make the HDUList
        psfhdu = fits.PrimaryHDU(data=PSF)
        wavhdu = fits.ImageHDU(data=W * 1E6, name='WAV')
        hdulist = fits.HDUList([psfhdu, wavhdu])

        # Write the file
        hdulist.writeto(file, overwrite=True)
        hdulist.close()

    except (ImportError, OSError, IOError):

        print("Could not import `webbpsf` package. Functionality limited. Generating dummy file.")

def get_angle(pf, p0=np.array([0, 0]), pi=None):
    """Compute angle (in degrees) for pf-p0-pi corner

    Parameters
    ----------
    pf: sequence
        The coordinates of a point on the rotated vector
    p0: sequence
        The coordinates of the pivot
    pi: sequence
        The coordinates of the fixed vector

    Returns
    -------
    float
        The angle in degrees
    """
    if pi is None:
        pi = p0 + np.array([0, 1])
    v0 = np.array(pf) - np.array(p0)
    v1 = np.array(pi) - np.array(p0)

    angle = np.math.atan2(np.linalg.det([v0, v1]), np.dot(v0, v1))
    angle = np.degrees(angle)

    return angle


def get_SOSS_psf(wavelength, filt='CLEAR', psfs=None, cutoff=0.005, plot=False):
    """
    Retrieve the SOSS psf for the given wavelength,
    scale the total flux to 1, and set pixels below
    cutoff value to zero

    Parameters
    ----------
    wavelength: float
        The wavelength to retrieve [um]
    filt: str
        The filter to use, ['CLEAR', 'F277W']
    psfs: numpy.interp1d object (optional)
        The interpolator
    plot: bool
        Plot the psf

    Returns
    -------
    np.ndarray
        The 2D psf for the input wavelength
    """
    if psfs is None:

        # Get the file
        file = os.path.join(PSF_DIR, 'SOSS_{}_PSF.fits'.format(filt))

        # Load the SOSS psf cube
        cube = fits.getdata(file).swapaxes(-1, -2)
        wave = fits.getdata(file, ext=1)

        # Initilize interpolator
        psfs = interp1d(wave, cube, axis=0, kind=3)

    # Check the wavelength
    if wavelength < psfs.x[0]:
        wavelength = psfs.x[0]

    if wavelength > psfs.x[-1]:
        wavelength = psfs.x[-1]

    # Interpolate and scale psf
    psf = psfs(wavelength)
    psf *= 1. / np.sum(psf)

    # Remove background
    # psf[psf < cutoff] = 0

    if plot:

        fig = figure()
        fig.image([psf], x=0, y=0, dw=psf.shape[0], dh=psf.shape[1])
        show(fig)

    else:
        return psf


def make_frame(psfs):
    """
    Generate a frame from an array of psfs

    Parameters
    ----------
    psfs: sequence
        An array of psfs of shape (2048, 76, 76)

    Returns
    -------
    np.ndarray
        An array of the SOSS psf at 2048 wavelengths for each order
    """
    # Empty frame
    frame = np.zeros((256, 2124))

    # Add each psf
    for n, psf in enumerate(psfs):
        frame[:, n:n + 76] += psf

    return frame[:, 38:-38]


def psf_lightcurve(psf, tmodel=None, time=None):
    """
    Generate a lightcurve for a (76, 76) psf of a given wavelength

    Parameters
    ----------
    psf: sequencs
        The flux-scaled psf for the given wavelength
    tmodel: batman.transitmodel.TransitModel
        The transit model of the planet
    time: sequence
        The time axis for the TSO

    Returns
    -------
    sequence
        A 1D array of the lightcurve with the same length as *t*

    Example 1
    ---------
    # No planet
    import numpy as np
    from mirage.psf import soss_trace as st
    psf = np.ones((76, 76))
    time = np.linspace(-0.2, 0.2, 200)
    lc = st.psf_lightcurve(psf, time)

    Example 2
    ---------
    # With a planet
    import batman
    import numpy as np
    import astropy.units as q
    from mirage.psf import soss_trace as st
    params = batman.TransitParams()
    params.t0 = 0.                                # time of inferior conjunction
    params.per = 5.7214742                        # orbital period (days)
    params.a = 0.0558*q.AU.to(q.R_sun)*0.66       # semi-major axis (in units of stellar radii)
    params.inc = 89.8                             # orbital inclination (in degrees)
    params.ecc = 0.                               # eccentricity
    params.w = 90.                                # longitude of periastron (in degrees)
    params.limb_dark = 'quadratic'                # limb darkening profile to use
    params.u = [1, 1]                             # limb darkening coefficients
    params.rp = 0.15                              # radius of the planet
    psf = np.ones((76, 76))
    time = np.linspace(-0.2, 0.2, 200)
    tmodel = batman.TransitModel(params, time)
    lc = st.psf_lightcurve(psf, tmodel, time)
    """
    # Expand to shape of time axis
    flux = np.tile(psf, (len(time), 1, 1))

    # If there is a transiting planet...
    if tmodel is not None:

        # Generate the light curve for this pixel
        lightcurve = tmodel.light_curve(tmodel)

        # Scale the flux with the lightcurve
        flux *= lightcurve[:, None, None]

    return flux


def psf_tilts(order):
    """
    Get the psf tilts for the given order

    Parameters
    ----------
    order: int
        The order to use, [1, 2]

    Returns
    -------
    np.ndarray
        The angle from the vertical of the psf in each of the 2048 columns
    """
    if order not in [1, 2]:
        raise ValueError('Only orders 1 and 2 are supported.')

    # Get the file
    psf_file = os.path.join(PSF_DIR, 'SOSS_PSF_tilt_order{}.npy'.format(order))

    if not os.path.exists(psf_file):
        calculate_psf_tilts()

    return np.load(psf_file)


def put_psf_on_subarray(psf, y, frame_height=256):
    """Make a 2D SOSS trace from a sequence of psfs and trace center locations

    Parameters
    ----------
    psf: sequence
        The 2D psf
    y: float
        The grid y value to place the center of the psf
    grid: sequence
        The [x, y] grid ranges

    Returns
    -------
    np.ndarray
        The 2D frame with the interpolated psf
    """
    # Create spline generator
    dim = psf.shape[0]
    mid = (dim - 1.0) / 2.0
    arr = np.arange(dim, dtype=np.float)
    spline = RectBivariateSpline(arr, arr, psf.T, kx=3, ky=3, s=0)

    # Create output frame, shifted as necessary
    yg, xg = np.indices((frame_height, dim), dtype=np.float64)
    yg += mid - y

    # Resample onto the subarray
    frame = spline.ev(xg, yg)

    # Fill resampled points with zeros
    extrapol = (((xg < -0.5) | (xg >= dim - 0.5)) | ((yg < -0.5) | (yg >= dim - 0.5)))
    frame[extrapol] = 0

    return frame


def SOSS_psf_cube(filt='CLEAR', order=1, subarray='SUBSTRIP256', generate=False, mprocessing=True):
    """
    Generate/retrieve a data cube of shape (3, 2048, 76, 76) which is a
    76x76 pixel psf for 2048 wavelengths for each trace order. The PSFs
    are scaled to unity and rotated to reproduce the trace tilt at each
    wavelength then placed on the desired subarray.

    Parameters
    ----------
    filt: str
        The filter to use, ['CLEAR', 'F277W']
    order: int
        The trace order
    subarray: str
        The subarray to use, ['SUBSTRIP96', 'SUBSTRIP256', 'FULL']
    generate: bool
        Generate a new cube
    mprocessing: bool
        Use multiprocessing

    Returns
    -------
    np.ndarray
        An array of the SOSS psf at 2048 wavelengths for each order
    """
    if generate:

        print('Coffee time! This takes about 5 minutes.')

        # Get the wavelengths
        wavelengths = np.mean(utils.wave_solutions(subarray), axis=1)[:2 if filt == 'CLEAR' else 1]
        coeffs = locate_trace.trace_polynomial(subarray)

        # Get the file
        psf_file = os.path.join(PSF_DIR, 'SOSS_{}_PSF.fits'.format(filt))

        # Load the SOSS psf cube
        cube = fits.getdata(psf_file).swapaxes(-1, -2)
        wave = fits.getdata(psf_file, ext=1)

        # Initilize interpolator
        psfs = interp1d(wave, cube, axis=0, kind=3)
        trace_cols = np.arange(2048)

        # Run datacube
        for n, wavelength in enumerate(wavelengths):

            # Evaluate the trace polynomial in each column to get the y-position of the trace center
            trace_centers = np.polyval(coeffs[n], trace_cols)

            # Don't calculate order2 for F277W or order 3 for either
            if (n == 1 and filt.lower() == 'f277w') or n == 2:
                pass

            else:

                # Get the psf for each column
                print('Calculating order {} SOSS psfs for {} filter...'.format(n + 1, filt))
                start = time.time()
                func = partial(get_SOSS_psf, filt=filt, psfs=psfs)

                if mprocessing:
                    pool = multiprocessing.Pool(8)
                    raw_psfs = np.array(pool.map(func, wavelength))
                    pool.close()
                    pool.join()
                    del pool
                else:
                    raw_psfs = []
                    for i in range(len(wavelength)):
                        raw_psfs.append(func(wavelength[i]))
                    raw_psfs = np.array(raw_psfs)

                print('Finished in {} seconds.'.format(time.time()-start))

                # Rotate the psfs
                print('Rotating order {} SOSS psfs for {} filter...'.format(n + 1, filt))
                start = time.time()
                func = partial(rotate, reshape=False)

                # Get the PSF tilt at each column
                angles = psf_tilts(order)

                if mprocessing:
                    pool = multiprocessing.Pool(8)
                    rotated_psfs = np.array(pool.starmap(func, zip(raw_psfs, angles)))
                    pool.close()
                    pool.join()
                    del pool
                else:
                    rotated_psfs = []
                    for rp, ang in zip(raw_psfs, angles):
                        rotated_psfs.append(func(rp, ang))
                    rotated_psfs = np.array(rotated_psfs)

                print('Finished in {} seconds.'.format(time.time()-start))

                # Scale psfs to 1
                rotated_psfs = np.abs(rotated_psfs)
                scale = np.nansum(rotated_psfs, axis=(1, 2))[:, None, None]
                rotated_psfs = rotated_psfs / scale

                # Split it into 4 chunks to be below Github file size limit
                chunks = rotated_psfs.reshape(4, 512, 76, 76)
                for N, chunk in enumerate(chunks):

                    idx0 = N * 512
                    idx1 = idx0 + 512
                    centers = trace_centers[idx0:idx1]

                    # Interpolate the psfs onto the subarray
                    print('Interpolating chunk {}/4 for order {} SOSS psfs for {} filter onto subarray...'.format(N + 1, n + 1, filt))
                    start = time.time()
                    func = put_psf_on_subarray

                    if mprocessing:
                        pool = multiprocessing.Pool(8)
                        data = zip(chunk, centers)
                        subarray_psfs = pool.starmap(func, data)
                        pool.close()
                        pool.join()
                        del pool
                    else:
                        subarray_psfs = []
                        for ch, ce in zip(chunk, centers):
                            subarray_psfs.append(func(ch, ce))

                    print('Finished in {} seconds.'.format(time.time()-start))

                    # Get the filepath
                    file = os.path.join(PSF_DIR, 'SOSS_{}_PSF_order{}_{}.npy'.format(filt, n+1, N+1))

                    # Delete the file if it exists
                    if os.path.isfile(file):
                        os.system('rm {}'.format(file))

                    # Write the data
                    np.save(file, np.array(subarray_psfs))

                    print('Data saved to', file)

    else:

        if PSF_DIR is None:
            print("No PSF files detected. Using all ones.")
            subarr = 256 if subarray == 'SUBSTRIP256' else 96 if subarray == 'SUBSTRIP256' else 2048
            return np.ones((2048, subarr, 76))

        else:

            # Get the chunked data and concatenate
            full_data = []
            for chunk in [1, 2, 3, 4]:
                file = os.path.join(PSF_DIR, 'SOSS_{}_PSF_order{}_{}.npy'.format(filt, order, chunk))
                full_data.append(np.load(file))

            return np.concatenate(full_data, axis=0)
