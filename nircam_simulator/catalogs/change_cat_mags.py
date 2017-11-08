#! /usr/bin/env python

'''
Randomly change the magnitudes in input catalogs, to
help seed image testing and development
'''

#infile = 'seed_im_from_catalog_test_ptsrc_catalog.list'
infile = 'seed_im_from_catalog_test_galaxies_catalog.list'
outfile = infile[0:-5] + '_filter2.list'

from astropy.io import ascii
from numpy.random import normal

cat = ascii.read(infile)

mags = cat['magnitude'].data
add = normal(mags,0.5)
cat['magnitude'] = add

ascii.write(cat,outfile,overwrite=True)

