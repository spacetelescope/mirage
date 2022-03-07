#!/usr/stsci/pyssg/Python-2.7/bin/python
"""
Module to download and normalize stellar model atmospheres for soss_simulator.SossModelSim class

Author: Nestor Espinoza
"""

import numpy as np
import os
from pkg_resources import resource_filename
from scipy import interpolate, ndimage
import shutil
import urllib.request as request

import astropy.units as q
from astropy.io import fits
from contextlib import closing


def closest_value(input_value, possible_values):
    """
    This function calculates, given an input_value and an array of possible_values,
    the closest value to input_value in the array.
    Parameters
    ----------
    input_value: double
         Input value to compare against possible_values.
    possible_values: np.ndarray
         Array of possible values to compare against input_value.
    Returns
    -------
    double
        Closest value on possible_values to input_value.
    """
    distance = np.abs(possible_values - input_value)
    idx = np.where(distance == np.min(distance))[0]

    return possible_values[idx[0]]


def get_atlas_folder(feh):
    """
    Given input metallicity, this function defines the first part of the URL that will define what
    file to download from the STScI website.
    Parameters
    ----------
    feh: np.double
         [Fe/H] of the desired spectrum.
    Returns
    -------
    string
        URL of ATLAS models closer to the input metallicity.
    """
    # Define closest possible metallicity from ATLAS models:
    model_metallicity = closest_value(feh, np.array([-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.2, 0.5]))
    met_sign = 'm'

    # Define the sign before the filename, obtain absolute value if needed:
    if model_metallicity >= 0.0:
        met_sign = 'p'
    else:
        model_metallicity = np.abs(model_metallicity)

    model_metallicity = ''.join(str(model_metallicity).split('.'))
    fname = 'https://archive.stsci.edu/hlsps/reference-atlases/cdbs/grid/ck04models/ck{0:}{1:}/'.format(met_sign,
                                                                                                        model_metallicity)

    return fname


def get_phoenix_folder(feh, alpha):
    """
    Given input metallicity and alpha-enhancement, this function defines the first part of the URL that will define what
    file to download from the PHOENIX site.
    Parameters
    ----------
    feh: np.double
         [Fe/H] of the desired spectrum.
    alpha: np.double
         Alpha-enhancement of the desired spectrum.
    Returns
    -------
    string
        FTP URL of PHOENIX file with the closest properties to the input properties.
    """
    # Define closest possible metallicity from PHOENIX models:
    model_metallicity = closest_value(feh, np.array([-4.0, -3.0, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0]))

    # Same for alpha-enhancement:
    model_alpha = closest_value(alpha, np.array([-0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.20]))
    met_sign, alpha_sign = '-', '-'

    # Define the sign before the filename, obtain absolute value if needed:
    if model_metallicity > 0.0:
        met_sign = '+'
    else:
        model_metallicity = np.abs(model_metallicity)
    if model_alpha > 0.0:
        alpha_sign = '+'
    else:
        model_alpha = np.abs(model_alpha)

    # Create the folder name
    if alpha == 0.0:
        fname = 'ftp://phoenix.astro.physik.uni-goettingen.de/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z{0:}{1:.1f}/'.format(
            met_sign, model_metallicity)
    else:
        fname = 'ftp://phoenix.astro.physik.uni-goettingen.de/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z{0:}{1:.1f}.Alpha={2:}{3:.2f}/'.format(
            met_sign, model_metallicity, alpha_sign, model_alpha)

    return fname


def download(url, fname):
    """
    Download files from ftp server at url in filename fname. Obtained/modified from jfs here: https://stackoverflow.com/questions/11768214/python-download-a-file-over-an-ftp-server
    Parameters
    ----------
    url: string
        URL pointing to the file to download.
    fname: string
        Output filename of the file to download.
    """
    with closing(request.urlopen(url)) as r:
        with open(fname, 'wb') as f:
            shutil.copyfileobj(r, f)


def get_vega():
    """
    This functions reads in the spectrum of Vega (Alpha Lyr) from CALSPEC
    Returns
    -------
    astropy.units.quantity.Quantity
        Wavelength in um of Vega spectrum.
    astropy.units.quantity.Quantity
        Flux in erg/s/cm2/A of Vega spectrum.
    """
    data = fits.getdata(resource_filename('mirage', 'reference_files/alpha_lyr_stis_009.fits'), header=False)

    # Wavelength is in Angstroms, convert to microns to match the get_phoenix_model function.
    # Flux is in Flambda (same as Phoenix; i.e., erg/s/cm2/A):
    return (data['WAVELENGTH'] * q.angstrom).to(q.um), data['FLUX'] * (q.erg / q.s / q.cm ** 2 / q.AA)


def read_phoenix_list(phoenix_model_list):
    """
    This function extracts filenames, effective temperatures and log-gs given a filename that contains a list of PHOENIX model filenames.
    Parameters
    ----------
    phoenix_model_list: string
         Filename of file containing, on each row, a PHOENIX model filename.
    Returns
    -------
    np.ndarray
        Array of PHOENIX model filenames
    np.ndarray
        Array containing the effective temperatures (K) of each PHOENIX model filename.
    np.nadarray
        Array containing the log-g (cgs) of each PHOENIX model filename.
    """
    fin = open(phoenix_model_list, 'r')
    fnames = np.array([])
    teffs = np.array([])
    loggs = np.array([])

    while True:
        line = fin.readline()

        if line != '':
            fname = line.split()[-1]
            teff, logg = fname.split('-')[:2]
            fnames = np.append(fnames, fname)
            teffs = np.append(teffs, np.double(teff[3:]))
            loggs = np.append(loggs, np.double(logg))

        else:
            break

    return fnames, teffs, loggs


def get_phoenix_model(feh, alpha, teff, logg):
    """
    This function gets you the closest PHOENIX high-resolution model to the input stellar parameters from the Goettingen website (ftp://phoenix.astro.physik.uni-goettingen.de).
    Parameters
    ----------
    feh: np.double
         [Fe/H] of the desired spectrum.
    alpha: np.double
         Alpha-enhancement of the desired spectrum.
    teff: np.double
         Effective temperature (K) of the desired spectrum.
    logg: np.double
         Log-gravity (cgs) of the desired spectrum.
    Returns
    -------
    np.ndarray
        Wavelength in um of the closest spectrum to input properties.
    np.ndarray
        Surface flux in f-lambda of the closest spectrum to input properties in units of erg/s/cm**2/angstroms.
    """
    # First get grid corresponding to input Fe/H and alpha:
    url_folder = get_phoenix_folder(feh, alpha)

    # Now define details for filenames and folders. First, extract metallicity and alpha-enhancement in
    # the PHOENIX filename format (and get rid of the "Z" in, e.g., "Z-1.0.Alpha=-0.20"):
    phoenix_met_and_alpha = url_folder.split('/')[-2][1:]

    # Define folders where we will save (1) all stellar model data and (2) all phoenix models:
    stellarmodels_folder_path = resource_filename('mirage', 'reference_files/stellarmodels/')
    phoenix_folder_path = resource_filename('mirage', 'reference_files/stellarmodels/phoenix/')
    model_folder_path = resource_filename('mirage', 'reference_files/stellarmodels/phoenix/' + phoenix_met_and_alpha + '/')

    # Check if we even have stellarmodels folder created. Create it if not:
    if not os.path.exists(stellarmodels_folder_path):
        os.mkdir(stellarmodels_folder_path)

    # Same for phoenix folder:
    if not os.path.exists(phoenix_folder_path):
        os.mkdir(phoenix_folder_path)

    # Check if the current metallicity-alpha folder exists as well:
    if not os.path.exists(model_folder_path):
        os.mkdir(model_folder_path)

    # Check if we have the PHOENIX wavelength solution. If not, download it:
    if not os.path.exists(phoenix_folder_path + 'wavsol.fits'):
        download('ftp://phoenix.astro.physik.uni-goettingen.de/HiResFITS/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits',
                      phoenix_folder_path + 'wavsol.fits')

    # Extract wavelength solution:
    wavelengths = fits.getdata(phoenix_folder_path + 'wavsol.fits')

    # Now, figure out the closest model to the input stellar parameters. For this, first figure out the range of teff and logg
    # for the current metallicity and alpha. For this, either retrieve from the system or download the full list of PHOENIX models
    # for the current metallicity and alpha. If not already here, save it on the system:
    phoenix_model_list = model_folder_path + 'model_list.txt'
    if not os.path.exists(phoenix_model_list):
        download(url_folder, phoenix_model_list)

    # Extract information from this list:
    model_names, possible_teffs, possible_loggs = read_phoenix_list(model_folder_path + 'model_list.txt')

    # Search the closest to the input teff:
    phoenix_teff = closest_value(teff, possible_teffs)

    # Raise a warning in case the found teff is outside the PHOENIX model range, give some
    # guidance on how to proceed:
    if np.abs(phoenix_teff - teff) > 200.:
        print(
            '\t Warning: the input stellar effective temperature is outside the {0:}-{1:} K model range of PHOENIX models for {2:}.'.format(
                np.min(possible_teffs),
                np.max(possible_teffs), phoenix_met_and_alpha))

        if 'Alpha' in phoenix_met_and_alpha:
            print(
                '\t Modelling using a {0:} K model. Using models without alpha-enhancement (alpha = 0.0), which range from 2300 to 12000 K would perhaps help find more suitable temperature models.'.format(
                    phoenix_teff))

        else:
            print('\t Modelling using a {0:} K model.'.format(phoenix_teff))

    # Same excercise for logg, given the teffs:
    idx_logg = np.where(np.abs(phoenix_teff - possible_teffs) == 0.)[0]
    phoenix_logg = closest_value(logg, possible_loggs[idx_logg])

    # Select final model:
    idx = np.where((np.abs(phoenix_teff - possible_teffs) == 0.) & (np.abs(possible_loggs == phoenix_logg)))[0]
    phoenix_model, phoenix_logg = model_names[idx][0], possible_loggs[idx][0]

    # Raise warning for logg as well:
    if np.abs(phoenix_logg - logg) > 0.5:
        print(
        '\t Warning: the input stellar log-gravity is outside the {0:}-{1:} model range of PHOENIX models for {2:} and Teff {3:}.'.format(
            np.min(possible_loggs[idx_logg]), np.max(possible_loggs[idx_logg]), phoenix_met_and_alpha, phoenix_teff))

    # Check if we already have the downloaded model. If not, download the corresponding file:
    if not os.path.exists(model_folder_path + phoenix_model):
        print('\t PHOENIX stellar models for {0:} not found in {1:}. Downloading...'.format(phoenix_met_and_alpha,
                                                                                            model_folder_path))
        download(url_folder + phoenix_model, model_folder_path + phoenix_model)

    # Once we have the file, simply extract the data:
    print('\t Using the {0:} PHOENIX model (Teff {1:}, logg {2:}).'.format(phoenix_model, phoenix_teff, phoenix_logg))
    flux = fits.getdata(model_folder_path + phoenix_model, header=False)

    # Change units in order to match what is expected by the TSO modules:
    wav = (wavelengths * q.angstrom).to(q.um)
    flux = (flux * (q.erg / q.s / q.cm ** 2 / q.cm)).to(q.erg / q.s / q.cm ** 2 / q.AA)

    return wav, flux


def get_atlas_model(feh, teff, logg):
    """
    This function gets you the closest ATLAS9 Castelli and Kurucz model to the input stellar parameters from the STScI website
    (http://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/castelli-and-kurucz-atlas).
    Parameters
    ----------
    feh: np.double
         [Fe/H] of the desired spectrum.
    teff: np.double
         Effective temperature (K) of the desired spectrum.
    logg: np.double
         Log-gravity (cgs) of the desired spectrum.
    Returns
    -------
    np.ndarray
        Wavelength in um of the closest spectrum to input properties.
    np.ndarray
        Surface flux in f-lambda of the closest spectrum to input properties in units of erg/s/cm**2/angstroms.
    """
    # First get grid corresponding to input Fe/H:
    url_folder = get_atlas_folder(feh)

    # Now define details for filenames and folders. Extract foldername with the metallicity info from the url_folder:
    atlas_met = url_folder.split('/')[-2]

    # Define folders where we will save (1) all stellar model data and (2) all atlas models:
    stellarmodels_folder_path = resource_filename('mirage', 'reference_files/stellarmodels/')
    atlas_folder_path = resource_filename('mirage', 'reference_files/stellarmodels/atlas/')
    model_folder_path = resource_filename('mirage', 'reference_files/stellarmodels/atlas/' + atlas_met + '/')

    # Check if we even have stellarmodels folder created. Create it if not:
    if not os.path.exists(stellarmodels_folder_path):
        os.mkdir(stellarmodels_folder_path)

    # Same for phoenix folder:
    if not os.path.exists(atlas_folder_path):
        os.mkdir(atlas_folder_path)

    # Check if the current metallicity-alpha folder exists as well:
    if not os.path.exists(model_folder_path):
        os.mkdir(model_folder_path)

    # Define possible teff and loggs (thankfully, this is easier for ATLAS models):
    possible_teffs, possible_loggs = np.append(np.arange(3500, 13250, 250), np.arange(14000, 51000, 1000)), np.arange(
        0.0, 5.5, 0.5)

    # Check the closest teff and logg to input ones:
    atlas_teff = closest_value(teff, possible_teffs)
    atlas_logg = closest_value(logg, possible_loggs)

    # Raise a warning in case the found teff is outside the ATLAS model range, give some
    # guidance on how to proceed:
    if np.abs(atlas_teff - teff) > 200.:
        print(
            '\t Warning: the input stellar effective temperature is outside the {0:}-{1:} K model range of ATLAS models for {2:}.'.format(
                np.min(possible_teffs),
                np.max(possible_teffs), atlas_met))
        print('\t Modelling using a {0:} K model.'.format(atlas_teff))

    # Now, if not already in the system, download the model corresponding to the chosen teff:
    atlas_fname = model_folder_path + atlas_met + '_{0:}.fits'.format(atlas_teff)
    if not os.path.exists(atlas_fname):
        download(url_folder + atlas_met + '_{0:}.fits'.format(atlas_teff), atlas_fname)

    # Read the file:
    d = fits.getdata(atlas_fname)

    # This variable will save non-zero logg at the given temperatures. Only useful to report back to the user and/or input logg
    # doesn't have data:
    real_possible_loggs = np.array([])

    # Check if the closest requested logg has any data. If not, check all possible loggs for non-zero data, and select the closest
    # to the input logg that has data:
    s_logg = 'g' + ''.join('{0:.1f}'.format(atlas_logg).split('.'))
    if np.count_nonzero(d[s_logg]) != 0:
        w, f = d['WAVELENGTH'], d[s_logg]
    else:
        real_possible_loggs = np.array([])
        for loggs in possible_loggs:
            s_logg = 'g' + ''.join('{0:.1f}'.format(loggs).split('.'))
            if np.count_nonzero(d[s_logg]) != 0:
                real_possible_loggs = np.append(real_possible_loggs, loggs)
        atlas_logg = closest_value(logg, real_possible_loggs)
        s_logg = 'g' + ''.join('{0:.1f}'.format(atlas_logg).split('.'))
        w, f = d['WAVELENGTH'], d[s_logg]

    # Raise warning for logg as well:
    if np.abs(atlas_logg - logg) > 0.5:

        # If real_possible_loggs is empty, calculate it:
        if len(real_possible_loggs) == 0:
            for loggs in possible_loggs:
                s_logg = 'g' + ''.join('{0:.1f}'.format(loggs).split('.'))
                if np.count_nonzero(d[s_logg]) != 0:
                    real_possible_loggs = np.append(real_possible_loggs, loggs)

        print(
        '\t Warning: the input stellar log-gravity is outside the {0:}-{1:} model range of ATLAS models for {2:} and Teff {3:}.'.format(
            np.min(real_possible_loggs), np.max(real_possible_loggs), atlas_met, atlas_teff))

    # Change units in order to match what is expected by the TSO modules:
    wav = (w * q.angstrom).to(q.um)
    flux = f * q.erg / q.s / q.cm ** 2 / q.AA
    return wav, flux


def get_resolution(w, f):
    """
    This function returns the (w) wavelength (median) resolution of input spectra (f)
    Parameters
    ----------
    w: np.ndarray
        Wavelengths of the spectrum
    f: np.ndarray
        Value at the given wavelength (can be flux, transmission, etc.)
    Returns
    -------
    np.double
        The median resolution of the spectrum.
    """
    eff_wav = np.sum(w * f) / np.sum(f)
    delta_wav = np.median(np.abs(np.diff(w)))

    return eff_wav / delta_wav


def spec_integral(input_w, input_f, wT, TT):
    """
    This function computes the integral of lambda*f*T divided by the integral of lambda*T, where
    lambda is the wavelength, f the flux (in f-lambda) and T the transmission function. The input
    stellar spectrum is given by wavelength w and flux f. The input filter response wavelengths
    are given by wT and transmission curve by TT. It is assumed both w and wT are in the same
    wavelength units.
    Parameters
    ----------
    input_w: np.ndarray
        Wavelengths of the input spectrum
    input_f: np.ndarray
        Flux (in f-lambda) of the input spectrum
    wT: np.ndarray
        Wavelength of the input transmission function
    TT: np.ndarray
        Spectral response function of the transmission function
    Returns
    -------
    np.double
        Value of the integral (over dlambda) of lambda*f*T divided by the integral (over dlambda) of lambda*T.
    """
    # If resolution of input spectra in the wavelength range of the response function
    # is higher than it, degrade it to match the transmission function resolution. First,
    # check that resolution of input spectra is indeed higher than the one of the
    # transmisssion. Resolution of input transmission first:
    min_wav, max_wav = np.min(wT), np.max(wT)
    resT = get_resolution(wT, TT)

    # Resolution of input spectra in the same wavelength range:
    idx = np.where((input_w >= min_wav - 10) & (input_w <= max_wav + 10))[0]
    res = get_resolution(input_w[idx], input_f[idx])

    # If input spetrum resolution is larger, degrade:
    if res > resT:

        # This can be way quicker if we just take the gaussian weight *at* the evaluated
        # points in the interpolation. TODO: make faster.
        f = ndimage.gaussian_filter(input_f[idx], int(np.double(len(idx)) / np.double(len(wT))))
        w = input_w[idx]
    else:
        w, f = input_w, input_f

    interp_spectra = interpolate.interp1d(w, f, fill_value='extrapolate')
    numerator = np.trapz(wT * interp_spectra(wT) * TT, x=wT)
    denominator = np.trapz(wT * TT, x=wT)

    return numerator / denominator


def scale_spectrum(w, f, jmag):
    """
    This function scales an input spectrum to a given jmag. This follows eq. (8) in Casagrande et al. (2014, MNRAS, 444, 392).
    Parameters
    ----------
    w: np.ndarray
        Wavelengths of the spectrum in microns.
    f: np.ndarray
        Flux of the spectrum in erg/s/cm2/A.
    jmag: np.double
        2MASS J-magnitude to which we wish to re-scale the spectrum.
    Returns
    -------
    np.ndarray
        Rescaled spectrum at wavelength w.
    """
    # Get filter response (note wT is in microns):
    wT, TT = np.loadtxt(resource_filename('mirage', 'reference_files/jband_transmission.dat'), unpack=True, usecols=(0, 1))

    # Get spectrum of vega:
    w_vega, f_vega = get_vega()

    # Use those two to get the absolute flux calibration for Vega (left-most term in equation (9) in Casagrande et al., 2014).
    # Multiply wavelengths by 1e4 as they are in microns (i.e., transform back to angstroms both wavelength ranges):
    vega_weighted_flux = spec_integral(np.array(w_vega.to(q.AA)), np.array(f_vega), wT * 1e4, TT)

    # J-band zero-point is thus (maginutde of Vega, m_*, obtained from Table 1 in Casagrande et al, 2014):
    ZP = -0.001 + 2.5 * np.log10(vega_weighted_flux)

    # Now compute (inverse?) bolometric correction for target star. For this, compute same integral as for vega, but for target:
    target_weighted_flux = spec_integral(np.array(w) * 1e4, f, np.array(wT) * 1e4, TT)

    # Get scaling factor for target spectrum (this ommits any extinction):
    scaling_factor = 10 ** (-((jmag + 2.5 * np.log10(target_weighted_flux) - ZP) / 2.5))

    return f * scaling_factor