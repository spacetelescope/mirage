#! /usr/bin/env python

'''
Randomly change the magnitudes in input catalogs, to
help seed image testing and development
'''

from astropy.io import ascii
from numpy.random import normal

#infile = 'seed_im_from_catalog_test_ptsrc_catalog.list'
#infile = 'seed_im_from_catalog_test_galaxies_catalog.list'
#outfile = infile[0:-5] + '_filter2.list'
infile = 'seed_im_from_catalog_test_ptsrc_bright_catalog.list'
outfile = infile

cat = ascii.read(infile)

mags = cat['magnitude'].data
add = mags - 2.
#add = normal(mags,0.5)
cat['magnitude'] = add

ascii.write(cat,outfile,overwrite=True)
