#! /usr/bin/env python

"""
Create a version of the xtalk file with all zeros for coefficients
"""

from astropy.io import ascii

infile = "niriss_xtalk.txt"

tab = ascii.read(infile)

for cname in tab.colnames:
    if cname not in ['Det', 'DetID']:
        tab[cname] *= 0.

ascii.write(tab,'niriss_xtalk_zeros.txt')
