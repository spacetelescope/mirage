"""Based on SIAF transforms PATHS
NIRCAM
------
science --> ideal --> V2,V3
V2,V3 --> ideal --> science
"""
import sys

import numpy as np
from astropy.modeling.models import Polynomial2D, Shift


#t = ascii.read("NIRCam_SIAF_2016-09-29.csv",header_start=1)

def get_siaf_transform(row, aperture, from_system, to_system, degree):
    """
    This reads in the file with transformations that the TEL team
    is using to construct the SIAF file. These transformations are
    defined as two polynomials (one in x, and one in y) transforming coordinates
    from one coordinate system to another.
    Parameters
    ----------
    aperture: str
        Name of aperture on NIRCam, composed of the detector name followed
        by an underscore and the subarray name. (e.g. "NRCA1_FULL",
        "NRCB5_SUB160")
    from_system : str
        Starting system (e.g. "science", "ideal")
    to_system : str
        Ending coordinate system (e.g. "science" , "ideal")
    degree : int
        Degree of polynomial
    Returns
    -------
    x_model : astropy.modeling.Model
        Correction in x
    y_model : astropy.modeling.Model
        Correction in y
    from_units : str
        Units in the starting system
    to_units : str
        Units in the ending system
    Examples
    --------
    >>> get_siaf_transform('NRCA1_FULL', "science", "ideal", 5)
    """
    #read in csv file of coefficients
    #t = ascii.read("NIRCam_SIAF_2016-09-29.csv",header_start=1)
    #t = ascii.read(coefffile,header_start=1)

    #from_system and to_system are very limited. Can only be "ideal" for the
    #distortion-free coords, and "science" for distorted coords
    from_system = from_system.lower()
    if from_system not in ['ideal','science']:
        print("Requested from_system of {} not recognized.".format(from_system))
        sys.exit()

    to_system = to_system.lower()
    if to_system not in ['ideal','science']:
        print("Requested to_system of {} not recognized.".format(to_system))
        sys.exit()

    #Generate the string corresponding to the requested coefficient labels
    if from_system == 'ideal' and to_system == 'science':
        label = 'Idl2Sci'
        from_units = 'arcsec'
        to_units = 'distorted pixels'
    elif from_system == 'science' and to_system == 'ideal':
        label = 'Sci2Idl'
        from_units = 'distorted pixels'
        to_units = 'arcsec'

    #Find the row that matches the requested aperture
    #match = t['AperName'] == aperture
    #if np.any(match) == False:
    #    print("Aperture name {} not found in input CSV file.".format(aperture))
    #    sys.exit()
    #row = t[match]

    #Get the coefficients for "science" to "ideal" transformation (and back)
    #"science" is distorted pixels. "ideal" is undistorted arcsec from the
    #the reference pixel location.

    #we need the parity value, which describes the relationship between the
    #v2,v3 coordinate system and the science pixel coordinate system
    parity = row['VIdlParity'].data[0]

    #Then create the model for the transformation
    X_cols = [c for c in row.colnames if label+'X' in c]
    #for kw in X_cols:
    #    ele = row[kw]
    #    ele *= parity
    #    row[kw] = ele
    x_coeffs = row[X_cols]
    X_model = to_model(x_coeffs,degree)

    Y_cols = [c for c in row.colnames if label+'Y' in c]
    y_coeffs = row[Y_cols]
    Y_model = to_model(y_coeffs,degree)

    return X_model, Y_model, from_units, to_units


def to_model(coeffs, degree=5):
    """
    Creates an astropy.modeling.Model object
    Parameters
    ----------
    coeffs : array like
        Coefficients from the ISIM transformations file.
    degree : int
        Degree of polynomial.
        Default is 5 as in the ISIM file but many of the polynomials are of
        a smaller degree.
    Returns
    -------
    poly : astropy.modeling.Polynomial2D
        Polynomial model transforming one coordinate (x or y) between two systems.
    """

    #map Colin's coefficients into the order expected by Polynomial2D
    c = {}
    for cname in coeffs.colnames:
        siaf_i = int(cname[-2])
        siaf_j = int(cname[-1])
        name = 'c{0}_{1}'.format(siaf_i-siaf_j,siaf_j)
        c[name] = coeffs[cname].data[0]

    #0,0 coefficient should not be used, according to Colin's TR
    #JWST-STScI-001550
    c['c0_0'] = 0

    return Polynomial2D(degree, **c)


def get_siaf_v2v3_transform(row,aperture,from_system='v2v3',to_system='v2v3'):
    """
    Generate transformation model to go to/from V2/V3 from
    undistorted angular distnaces from the reference pixel ("ideal")
    """
    #read in csv file of coefficients
    #t = ascii.read("NIRCam_SIAF_2016-09-29.csv",header_start=1)
    #t = ascii.read(coefffile,header_start=1)

    from_system = from_system.lower()
    to_system = to_system.lower()

    if from_system != 'v2v3' and to_system != 'v2v3':
        print("WARNING, either from_system or to_system must be 'v2v3'")
        sys.exit()

    #Find the row that matches the requested aperture
    #match = t['AperName'] == aperture
    #if np.any(match) == False:
    #    print("Aperture name {} not found in input CSV file.".format(aperture))
    #    sys.exit()
    #row = t[match]

    #Then create the model for the transformation
    parity = row['VIdlParity'].data[0]
    v3_ideal_y_angle = row['V3IdlYAngle'].data[0] * np.pi / 180.

    #print("parity and angle are {}, {}".format(parity,v3_ideal_y_angle))

    X_model, Y_model = v2v3_model(from_system,to_system,parity,v3_ideal_y_angle)

    return X_model, Y_model


def v2v3_model(from_sys, to_sys, par, angle):
    """
    Creates an astropy.modeling.Model object
    for the undistorted ("ideal") to V2V3 coordinate translation
    """
    if from_sys != 'v2v3' and to_sys != 'v2v3':
        print("This function is designed to generate the transformation either to or from V2V3.")
        sys.exit()

    #cast the transform functions as 1st order polynomials
    xc = {}
    yc = {}
    if to_sys == 'v2v3':
        xc['c1_0'] = par * np.cos(angle)
        xc['c0_1'] = np.sin(angle)
        yc['c1_0'] = (0.-par) * np.sin(angle)
        yc['c0_1'] = np.cos(angle)

    if from_sys == 'v2v3':
        xc['c1_0'] = par * np.cos(angle)
        xc['c0_1'] = par * (0. - np.sin(angle))
        yc['c1_0'] = np.sin(angle)
        yc['c0_1'] = np.cos(angle)

    #0,0 coeff should never be used.
    xc['c0_0'] = 0
    yc['c0_0'] = 0

    #print("coeffs for v2v3 transform:")
    #for key in xc:
    #    print("{} {}".format(key,xc[key]))
    #sys.exit()

    xmodel = Polynomial2D(1, **xc)
    ymodel = Polynomial2D(1, **yc)

    return xmodel, ymodel


def get_refpix(row):
    #Return the reference location within the given aperture
    xref = row['XSciRef'].data[0]
    yref = row['YSciRef'].data[0]
    return Shift(-xref) & Shift(-yref)


def get_v2v3ref(row):
    #Return v2 and v3 at the reference location
    #These are arcsec in the SIAF file. Convert to degrees
    v2ref = row['V2Ref'].data[0]
    v3ref = row['V3Ref'].data[0]
    return Shift(v2ref) & Shift(v3ref)
