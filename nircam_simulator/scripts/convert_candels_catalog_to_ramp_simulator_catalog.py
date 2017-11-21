# ! /usr/bin/env python

'''
convert a CANDELS source catalog (downloaded from MAST)
to the right format to be used in the ramp simulator
'''

from astropy.io import ascii
from astropy.table import Table
import numpy as np
import os

file = '../catalogs/hlsp_candels_hst_wfc3_goodss-tot-multiband_f160w_v1-1photom_cat.txt'
# file = '../catalogs/hlsp_candels_hst_wfc3_egs-tot-multiband_f160w_v1-1photom_cat.txt'
# file = '../source_catalogs/candels/CANDELS.GOODSS.F160W.v1_1.photom.cat.txt'
catdir = '../catalogs/'
catpath, catname = os.path.split(file)
catname = catdir + catname
# in output catalog save only one source out of each
# crop factor number of lines in the original
crop_factor = 10 # 100 for GOODS-S. Use 10 for candels catalog in north

galaxy_crop_factor = 40 #  120 for goods-s. use 40-for candels catalog in north

# for file ../source_catalogs/hlsp_candels_hst_wfc3_egs-tot-multiband_f160w_v1-1photom_cat.txt, the central RA and Dec I'm cropping around are:
# RA: 214.9
# Dec: 2.88
# Want to keep enough coverage for all 8 detectors at any roll angle,
#  + /- 0.1 degrees in each dimension should be plenty.

# # for GOODS-S catalog
# ra_min = 53.05
# ra_max = 53.20
# dec_min = -27.84
# dec_max = -27.76

# Trying to do OTE-01 ?? 11/13
ra_diff = 0.15
dec_diff = .08

ra_cent = 9.79
dec_cent = 63.25

ra_min = ra_cent - 0.5 * ra_diff
ra_max = ra_cent + 0.5 * ra_diff
dec_min = dec_cent - 0.5 * dec_diff
dec_max = dec_cent + 0.5 * dec_diff

# files listing zeropoints and Flam so that we can convert fluxes to
# magnitudes later
photoma = '/grp/jwst/wit/nircam/reference_files/photom/NIRCam_zeropoints_modA.txt'
photomb = '/grp/jwst/wit/nircam/reference_files/photom/NIRCam_zeropoints_modB.txt'

# read in catalog
dat = ascii.read(file)
print('RA: {} - {}'.format(min(dat['RA']), max(dat['RA'])))

# columns to keep
# for candels catalog in north:
# keeplong = ['RA', 'DEC', 'WIRCAM_J_FLUX', 'WIRCAM_H_FLUX', 'WIRCAM_K_FLUX', 'IRAC_CH1_FLUX', 'IRAC_CH2_FLUX']

# for goods-s
keeplong = ['RA', 'DEC', 'WFC3_F160W_FLUX', 'ISAAC_KS_FLUX', 'HAWKI_KS_FLUX', 'IRAC_CH1_FLUX', 'IRAC_CH2_FLUX']
keepshort = ['RA', 'DEC', 'WFC3_F125W_FLUX']
keepall = [keeplong, keepshort]

filterslong = ['F250M', 'F300M', 'F335M', 'F410M', 'F480M']
filtersshort = ['F182M']
filtersall = [filterslong, filtersshort]

effwavlong = np.array([2.503, 2.989, 3.362, 4.082, 4.874]) * 1.e4
effwavshort = np.array([1.845]) * 1.e4
effwavall = [effwavlong, effwavshort]

i=0
for keep, filters, effwav in zip(keepall, filtersall, effwavall):
    dat2 = dat[keep]

    good = np.where((dat2['RA'].data >= ra_min) & (dat2['RA'].data <= ra_max) & (dat2['DEC'].data >= dec_min) & (dat2['DEC'].data <= dec_max))
    dat2 = dat2[good]

    print("After cropping by RA, Dec, number of remaining entries: {}".format(len(dat2)))

    # deal with bad entries
    for col in dat2.colnames:
        if col not in ['RA', 'DEC']:
            column = dat2[col]
            bad = column <= 0
            column[bad] = 0.0000005

            bad = np.isnan(column)
            column[bad] = 0.0000005

    # crop using ra_min, ra_max, dec_min, dec_max
    good = np.where((dat2['RA'] >= ra_min) & (dat2['RA'] <= ra_max) & (dat2['DEC'] < dec_min) & (dat2['DEC'] <= dec_max))[0]


# all input fluxes are in microjansky. Let's convert to vegamag

# flam units: erg s-1 cm-2 A-1
# microjansky: 10-29 erg s-1 cm-2 Hz-1
# c = 2.99792458e18 # A/s

    # read in file with zeropoints
    zp = ascii.read(photomb)

    # filters = ['F250M', 'F300M', 'F335M', 'F410M', 'F480M']
    # effwav = np.array([2.503, 2.989, 3.362, 4.082, 4.874]) * 1.e4
    fluxes = keep[2:]

    dat3 = Table()
    dat3['x_or_RA'] = dat2['RA']
    dat3['y_or_Dec'] = dat2['DEC']
    meta0 = ''
    meta1 = ''
    meta2 = 'catalog for ramp_simulator.py. created from candels catalog'
    meta3 = 'Magnitudes are converted from input flux densities. See'
    meta4 = 'convert_candels_catalog_to_ramp_simulator_catalog.py for details'
    dat3.meta['comments'] = [meta0, meta1, meta2, meta3]


    for flux, filter, wave in zip(fluxes, filters, effwav):
        mtch = zp['Filter'] == filter
        vega = zp['VEGAMAG'][mtch]
        flam = zp['FLAM'][mtch]

        # convert microjansky to units of flam
        obj_flam = 2.99792458E-5 * (dat2[flux] / 1.e6) / wave**2
        print(filter, wave, flux)
        # print(dat2[flux][0])
        # print(obj_flam[0])
        # print(wave)

        # we can't just convert filter to filter because of
        # bandwidth differences, etc etc. So find the calculated
        # flambda for the first source, and scale it so that it'll
        # come out to be, say 20th magnitude. Then scale all other
        # sources by the same factor
        mag0 = 2.5 * np.log10(flam/obj_flam[0])

        deltamag = 20. - mag0
        fluxratio = 10.**(deltamag/(0-2.5))

        # scale all fluxes
        obj_flam *= fluxratio

        # now convert flam to vegamags
        mag = 2.5 * np.log10(flam/obj_flam)
        dat3['magnitude'] = mag

        # print(mag0)
        # print(deltamag, fluxratio)
        # print(flam, obj_flam[0])
        print(np.min(mag), np.max(mag))
        # stop

        if crop_factor > 1:
            good = np.arange(0, len(dat3), crop_factor)
            dat4 = dat3[good]

            # let's use some of the cropped entries to make a galaxy catalog
            goodgal = np.setdiff1d(np.arange(0, len(dat3)), good)
            dat_gal1 = dat3[goodgal]
        else:
            dat4 = dat3

        dat4.meta['comments'] = [meta0, meta1, meta2, meta3]


        # save to a catalog file. one per filter
        outname = catname + '_' + filter + '_cropby' + str(crop_factor) + '_ra_' + \
                  str(ra_min) + '_' + str(ra_max) + '_dec_' + str(dec_min) + '_' + \
                  str(dec_max) + '.cat'
        ascii.write(dat4, outname, overwrite=True)# , format='ascii')
        print('Point source output saved to {}'.format(outname))

        if galaxy_crop_factor > 1:
            good2 = np.arange(0, len(dat_gal1), galaxy_crop_factor)
            dat_gal2 = dat_gal1[good2]
        else:
            dat_gal2 = dat_gal1

        dat_gal2.meta['comments'] = [meta0, meta1, meta2, meta3]
        dat_gal2['magnitude'] -= 1

        # Now we need to add columns for radius, ellipticity, pos_angle, sersic_index
        # and we should also bump up the magnitudes a bit
        if i == 0:
            nrow = len(dat_gal2)
            radii = np.random.random_sample(nrow) * 0.75
            ellip = np.random.random_sample(nrow)
            pos = np.random.random_sample(nrow) * 359.9
            sersic = np.random.random_sample(nrow) * 4.


        dat_gal3 = Table([dat_gal2['x_or_RA'].data, dat_gal2['y_or_Dec'].data, radii, ellip, pos, sersic, dat_gal2['magnitude'].data], names=['x_or_RA', 'y_or_Dec', 'radius', 'ellipticity', 'pos_angle', 'sersic_index', 'magnitude'])
        dat_gal3.meta['comments'] = [meta0, meta1, meta1, meta2, meta3]
        goutname = catname + '_toGALAXYcat_' + filter + '_cropby' + str(galaxy_crop_factor) + '_ra_' + str(ra_min) + '_' + str(ra_max) + '_dec_' + str(dec_min) + '_' + str(dec_max) + '.cat'
        ascii.write(dat_gal3, goutname, overwrite=True)
        print('Galaxy source output saved to {}'.format(goutname))
        i = i + 1
