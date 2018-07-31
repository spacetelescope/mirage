'''A collection of basic routines for doing rotation calculations
Changed to express quaternion as vector before scalar part(a, p, )
'''

from math import radians, cos, sin, degrees, asin, atan2, sqrt, acos

import numpy as np


def unit(ra, dec):
    ''' Converts vector expressed in Euler angles to unit vector components.
    ra and dec in degrees
    Can be used for V2V3 after converting from arcsec to degrees)'''

    rar = radians(ra)
    decr = radians(dec)
    u = np.array([cos(rar)*cos(decr), sin(rar)*cos(decr), sin(decr)])
    return u


def radec(u):
    '''convert unit vector to Euler angles
    u is an array or list of length 3'''

    if len(u) != 3:
        print('Not a vector')
        return
    norm = sqrt(u[0]**2 + u[1]**2 + u[2]**2) # Works for list or array
    dec = degrees(asin(u[2]/norm))
    ra = degrees(atan2(u[1], u[0])) # atan2 puts it in the correct quadrant
    if ra < 0.0: ra += 360.0    # Astronomers prefer the range 0 to 360 degrees
    return (ra, dec)


def v2v3(u):
    '''Convert unit vector to v2v3'''
    if len(u) != 3:
        print ('Not a vector')
        return
    norm = sqrt(u[0]**2 + u[1]**2 + u[2]**2) # Works for list or array
    v2 = 3600*degrees(atan2(u[1], u[0])) # atan2 puts it in the correct quadrant
    v3 = 3600*degrees(asin(u[2]/norm))
    return (v2, v3)


def rotate(axis, angle):
    '''Fundamental rotation matrices.
    Rotate by angle measured in degrees, about axis 1 2 or 3'''
    if axis not in list(range(1, 4)):
        print ('Axis must be in range 1 to 3')
        return
    r = np.zeros((3, 3))
    ax0 = axis-1 #Allow for zero offset numbering
    theta = radians(angle)
    r[ax0, ax0] = 1.0
    ax1 = (ax0+1) % 3
    ax2 = (ax0+2) % 3
    r[ax1, ax1] = cos(theta)
    r[ax2, ax2] = cos(theta)
    r[ax1, ax2] = -sin(theta)
    r[ax2, ax1] = sin(theta)
    return r


def rv(v2, v3):
    '''Rotate from v2, v3 position to V1 axis'''
    v2d = v2/3600.0 #convert from arcsec to degrees
    v3d = v3/3600.0
    mv2 = rotate(3, -v2d)
    mv3 = rotate(2, v3d)
    rv = np.dot(mv3, mv2)
    return rv


def slew(v2t, v3t, v2a, v3a):
    ''' Calculate matrix which slews from target (v2t, v3t)
    to aperture position (v2a, v3a)'''
    v2td = v2t/3600.0
    v3td = v3t/3600.0
    v2ad = v2a/3600.0
    v3ad = v3a/3600.0
    r1 = rotate(3, -v2td)
    r2 = rotate(2, v3td-v3ad)
    r3 = rotate(3, v2ad)
    # Combine r3 r2 r1
    mv = np.dot(r2, r1)
    mv = np.dot(r3, mv)
    return mv


def attitude(v2, v3, ra, dec, pa):
    '''This will make a rotation matrix which rotates a unit vector representing a v2, v3 position
    to a unit vector representing an RA, Dec pointing with an assigned position angle
    Described in JWST-STScI-001550, SM-12, section 6.1'''

    # v2, v3 in arcsec, ra, dec and position angle in degrees
    v2d = v2/3600.0
    v3d = v3/3600.0

    # Get separate rotation matrices
    mv2 = rotate(3, -v2d)
    mv3 = rotate(2, v3d)
    mra = rotate(3, ra)
    mdec = rotate(2, -dec)
    mpa = rotate(1, -pa)

    # Combine as mra*mdec*mpa*mv3*mv2
    m = np.dot(mv3, mv2)
    m = np.dot(mpa, m)
    m = np.dot(mdec, m)
    m = np.dot(mra, m)

    return m


def pointing(attitude, v2, v3):
    '''Using the attitude matrix to calculate where any v2v3 position points on the sky'''
    v2d = v2/3600.0
    v3d = v3/3600.0
    v = unit(v2d, v3d)
    w = np.dot(attitude, v)
    rd = radec(w)
    return rd # tuple containing ra and dec


def getv2v3(attitude, ra, dec):
    ''' Using the inverse of attitude matrix
    find v2, v3 position of any RA and DEC'''
    urd = unit(ra, dec)
    invAtt = np.linalg.inv(attitude)
    uv = np.dot(invAtt, urd)
    (v2, v3) = radec(uv)
    if v2 > 180.0: v2 = v2-360.0
    v2 = 3600.0*v2
    v3 = 3600.0*v3
    return (v2, v3)


def posangle(attitude, v2, v3):
    ''' Using the attitude matrix find the V3 angle at arbitrary v2, v3
    This is the angle measured from North to V3 in an anti-clockwise direction
    i.e. North to East
    Formulae from JWST-STScI-001550, SM-12, section 6.2
    Subtract 1 from each index in the text to allow for python zero indexing'''

    A = attitude  # Synonym to simplify typing
    v2r = radians(v2/3600.0)
    v3r = radians(v3/3600.0)
    x = -(A[2, 0]*cos(v2r) + A[2, 1]*sin(v2r))*sin(v3r) + A[2, 2]*cos(v3r)
    y = (A[0, 0]*A[1, 2] - A[1, 0]*A[0, 2])*cos(v2r) + (A[0, 1]*A[1, 2] - A[1, 1]*A[0, 2])*sin(v2r)
    pa = degrees(atan2(y, x))
    return pa


def rodrigues(attitude):
    '''Interpret rotation matrix as a single rotation by angle phi around unit length axis
    Return axis, angle and matching quaternion'''

    A = attitude # Synonym for clarity and to save typing
    cos_phi = 0.5*(A[0, 0] + A[1, 1] + A[2, 2] - 1.0)
    phi = acos(cos_phi)
    axis = np.array([A[2, 1]-A[1, 2], A[0, 2]-A[2, 0], A[1, 0]-A[0, 1]])/(2.0*sin(phi))
    # Make corresponding quaternion
    q = np.hstack(( axis*sin(phi/2.0), [cos(phi/2.0)]))
    phi = degrees(phi)
    return (axis, phi, q)


def cross(a, b):
    """cross product of two vectors"""
    c = np.array([a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]])
    return c


def axial(ax, phi, u):
    ''' Apply direct rotation to a vector using Rodrigues' formula
    ax is unit axis vector  phi is rotation angle in degrees
    u is initial vector'''
    rphi = radians(phi)
    v = u*cos(rphi) + cross(ax, u)*sin(rphi) + ax*np.dot(ax, u)*(1-cos(rphi))
    return v
