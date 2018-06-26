# SIAF file - I got an excel file from Colin that I saved as a csv

from astropy.io import ascii
import numpy as np

siaf_file = '/ifs/jwst/wit/witserv/data4/nrc/hilbert/simulated_data/NIRCam_SIAF_2017-03-28.csv'

# Read in SIAF
siaf = ascii.read(siaf_file,header_start=1)

# Find the desired row of the SIAF based on
# aperture name (e.g. "NRCA1_FULL")
match = siaf['AperName'] == ap_name
if np.any(match) == False:
    print("Aperture name {} not found in input CSV file.".
          format(aperture))
    sys.exit()

# Extract the row from the SIAF
siaf_row = siaf[match]

# Turn coefficients into functions
sci2idlx, sci2idly, from_units, to_units = read_siaf_table.get_siaf_transform(siaf_row,ap_name,'science','ideal',5)

# Translate detector x,y to x_ideal,y_ideal
# x and y can be single numbers, or lists of numbers
x = [1024,1032]
y = [1024,1000]
x_ideal = sci2idlx(x,y)
y_ideal = sci2idly(x,y)

