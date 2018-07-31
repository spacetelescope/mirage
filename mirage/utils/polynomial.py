#!/usr/stsci/pyssg/Python-2.7/bin/python

"""
polynomial.py

Created by Colin Cox on 2013-04-29.
Based on surfit.pro
"""

from math import *

import scipy
from scipy import linalg
import matplotlib.pyplot as P


def triangle(A, order):
    '''Print coefficients in triangular layout'''
    k = 0
    for i in range(order+1):
        for j in range(i+1):
            #print('%12.5e' %A[k], end=' ')
            print('{:f.3}'.format(A[k]))
            k += 1
        print()


def triangulate(A, order):
    '''Convert linear array to 2-D array with triangular coefficient layout'''
    AT = scipy.zeros((order+1, order+1))
    k = 0
    for i in range(order+1):
        for j in range(i+1):
            AT[i, j] = A[k]
            k += 1
    return AT


def flatten(A, order):
    '''Convert triangular layout to linear array'''
    terms = (order+1)*(order+2)//2
    AF = scipy.zeros(terms)
    k = 0
    for i in range(order+1):
        for j in range(i+1):
            AF[k] = A[i, j]
            k += 1
    return AF


def reorder(A, B, verbose=False) :
    '''Reorder Sabatke coefficients to my convention'''
    order = 5
    terms = (order+1)*(order+2)//2
    Aarray = scipy.zeros((order+1, order+1))
    Barray = scipy.zeros((order+1, order+1))
    #Carray = scipy.zeros((order+1, order+1))
    #Darray = scipy.zeros((order+1, order+1))
    k1 = 0
    for i in range(order+1):
        for j in range(order+1-i):
            Aarray[j, i] = A[k1]
            Barray[j, i] = B[k1]
            #Carray[j, i] = C[k1]
            #Darray[j, i] = D[k1]
            k1 += 1

    A2 = scipy.zeros((terms))
    B2 = scipy.zeros((terms))
    #C2 = scipy.zeros((terms))
    #D2 = scipy.zeros((terms))
    k2 = 0
    for i in range(order+1):
        for j in range(i+1):
            A2[k2] = Aarray[j, i-j]
            B2[k2] = Barray[j, i-j]
            #C2[k2] = Carray[j, i-j]
            #D2[k2] = Darray[j, i-j]
            k2 += 1

    if verbose:
        print('A')
        triangle(A2, order)
        print('\nB')
        triangle(B2, order)
        #print '\nC'
        #polyfit(C2, order)
        #print '\nD'
        #triangle(D2, order)

    return (A2, B2) #, C2, D2)


def poly(a, x, y, order):
    pol = 0.0
    k = 0 # index for coefficients
    for i in range(order+1):
        for j in range(i+1):
            pol = pol + a[k]*x**(i-j)*y**j
            k+=1
    return pol


def dpdx(a, x, y, order): # Differential wrt x
    dpdx = 0.0
    k = 1 # index for coefficients
    for i in range(1, order+1):
        for j in range(i+1):
            if i-j > 0: dpdx = dpdx + (i-j)*a[k]*x**(i-j-1)*y**j
            k+=1
    return dpdx


def dpdy(a, x, y, order): # Differential wrt y
    dpdy = 0.0
    k = 1 # index for coefficients
    for i in range(1, order+1):
        for j in range(i+1):
            if j > 0: dpdy = dpdy + j*a[k]*x**(i-j)*y**(j-1)
            k+=1
    return dpdy


def jacob(a, b, x, y, order):
    '''Calculation of Jacobean, or relative area'''
    j = dpdx(a, x, y, order)*dpdy(b, x, y, order) - dpdx(b, x, y, order)*dpdy(a, x, y, order)
    j = scipy.fabs(j)
    return j


def invert(a, b, u, v, n, verbose=False):
    '''Given that order n polynomials of (x, y) have the result (u, v), find (x, y)
    Newton Raphson method in two dimensions'''
    tol = 1.0e-10
    err = 1.0
    #Initial guesses - Linear approximation
    det = a[1]*b[2] - a[2]*b[1]
    x0 = (b[2]*u-a[2]*v)/det
    y0 = (-b[1]*u + a[1]*v)/det
    if verbose: print('Initial guesses', x0, y0)
    x = x0
    y = y0
    X = scipy.array([x, y])
    iter = 0
    while err > tol:
        f1 = scipy.array([poly(a, x, y, n)-u, poly(b, x, y, n)-v])
        j = scipy.array([[dpdx(a, x, y, n), dpdy(a, x, y, n)], [dpdx(b, x, y, n), dpdy(b, x, y, n)]])
        invj = scipy.linalg.inv(j)
        X = X - scipy.dot(invj, f1)
        if verbose: print('[X1, Y1]', X)
        x1 = X[0]
        y1 = X[1]
        err = hypot(x-x1, y-y1)
        if verbose: print('Error %10.2e' % err)
        [x, y] = [x1, y1]
        iter += 1

    return  (x, y, err, iter)


def choose(n, r):
    '''The number of ways of choosing r items from n'''
    if n<0 or r <0:
        print('Negative values not allowed')
        return 0
    if r > n:
        print('r must not be greater than n')
        return 0

    combin = 1
    if r > n/2: r1 = n-r
    else: r1 = r
    for k in range(r1): combin = combin*(n-k)//(k+1)
    return combin


def ShiftCoeffs(a, xshift, yshift, order, verbose=False):
    '''Calculate coefficients of polynomial when shifted to new origin'''

    # First place in triangular layout
    at = scipy.zeros([order+1, order+1])
    atshift = scipy.zeros([order+1, order+1])
    ashift = scipy.zeros([len(a)]) # Copy shape of a
    k = 0
    for p in range(order+1):
        for q in range(p+1):
            at[p, q] = a[k]
            k+=1

    # Apply shift
    for p in range(order+1):
        for q in range(p+1):
            if verbose: print("A'%1d%1d" %(p, q))
            for i in range(p, order+1):
                for j in range(q, i+1-(p-q)):
                    f = choose(j, q)*choose(i-j, p-q)
                    atshift[p, q] = atshift[p, q] + f*xshift**((i-j)-(p-q))*yshift**(j-q)*at[i, j]
                    if verbose: print('%2d A(%1d, %1d) x^%1d y^%1d' %(f, i , j, i-j-(p-q), (j-q)))
            if verbose: print()

    # Put back in linear layout
    k = 0
    for p in range(order+1):
        for q in range(p+1):
            ashift[k] = atshift[p, q]
            k+=1

    return ashift


def polyfit(u, x, y, order):
    '''Fit polynomial to a set of u values on an x, y grid
    u is a function u(x, y) being a polynomial of the form
    u = a[i, j] x**(i-j) y**j. x and y can be on a grid or be arbitrary values'''

    # First set up x and y powers for each coefficient
    px = []
    py = []
    for i in range(order+1):
        for j in range(i+1):
            px.append(i-j)
            py.append(j)
    terms = len(px)
    #print terms, ' terms for order ', order
    #print px
    #print py

    # Make up matrix and vector
    vector = scipy.zeros((terms))
    mat = scipy.zeros((terms, terms))
    for i in range(terms):
        vector[i] = (u*x**px[i]*y**py[i]).sum()
        for j in range(terms):
            mat[i, j] =  (x**px[i]*y**py[i]*x**px[j]*y**py[j]).sum()

    #print 'Vector', vector
    #print 'Matrix'
    #print mat
    imat = linalg.inv(mat)
    #print 'Inverse'
    #print imat
    # Check that inversion worked
    #print scipy.dot(mat, imat)
    coeffs = scipy.dot(imat, vector)
    return coeffs


def polyfit2(u, x, y, order):
    '''Fit polynomial to a set of u values on an x, y grid
    u is a function u(x, y) being a polynomial of the form
    u = a[i, j]x**(i-j)y**j. x and y can be on a grid or be arbitrary values
    This version uses solve instead of matrix inversion'''

    # First set up x and y powers for each coefficient
    px = []
    py = []
    for i in range(order+1):
        for j in range(i+1):
            px.append(i-j)
            py.append(j)
    terms = len(px)
    #print terms, ' terms for order ', order
    #print px
    #print py

    # Make up matrix and vector
    vector = scipy.zeros((terms))
    mat = scipy.zeros((terms, terms))
    for i in range(terms):
        vector[i] = (u*x**px[i]*y**py[i]).sum() # Summing over all x, y
        for j in range(terms):
            mat[i, j] =  (x**px[i]*y**py[i]*x**px[j]*y**py[j]).sum()

    coeffs = linalg.solve(mat, vector)
    return coeffs


def testpoly():
    [x, y] = scipy.mgrid[0:10, 0:10]
    #print 'X'
    #print x
    #print 'Y'
    #print y
    u = scipy.zeros((10, 10))
    v = scipy.zeros((10, 10))
    # Random polynomials
    a0 = scipy.random.rand(1)
    a1 = 0.1*(scipy.random.rand(2)-0.5)
    a2 = 0.01*(scipy.random.rand(3)-0.5)
    a = scipy.concatenate((a0, a1))
    a = scipy.concatenate((a, a2))
    a[2] = 0.01*a[2]
    print('A coefficients')
    print(a)
    b0 = scipy.random.rand(1)
    b1 = 0.1*(scipy.random.rand(2)-0.5)
    b2 = 0.01*(scipy.random.rand(3)-0.5)
    b = scipy.concatenate((b0, b1))
    b = scipy.concatenate((b, b2))
    b[1] = 0.01*b[1]
    print('B coeffcicients')
    print(b)
    for i in range(10):
        for j in range(10):
            u[i, j] = poly(a, x[i, j], y[i, j], 2) #+ scipy.random.normal(0.0, 0.01)
            v[i, j] = poly(b, x[i, j], y[i, j], 2)  #+ scipy.random.normal(0.0, 0.01)
    #print z
    s1 = polyFit2(u, x, y, 2)
    s2 = polyFit2(v, x, y, 2)
    print('S1', s1)
    print('S2', s2)
    uc = poly(s1, x, y, 2)
    vc = poly(s2, x, y, 2)

    P.figure(1)
    P.clf()
    P.grid(True)
    P.plot(u, v, 'gx')
    P.plot(uc, vc, 'r+')


def FlipX(A, order):
    """Change sign of all coefficients with odd x power"""
    terms = (order+1)*(order+2)//2
    AF = scipy.zeros(terms)
    k = 0
    for i in range(order+1):
        for j in range (i+1):
            AF[k] = (-1)**(i-j)*A[k]
            k += 1
    return  AF


def FlipY(A, order):
    """Change sign of all coefficients with odd y power"""
    terms = (order+1)*(order+2)//2
    AF = scipy.zeros(terms)
    k = 0
    for i in range(order+1):
        for j in range (i+1):
            AF[k] = (-1)**(j)*A[k]
            k += 1
    return AF


def FlipXY(A, order):
    "Change sign for coeffs where sum of x and y powers is odd"
    terms = (order+1)*(order+2)//2
    AF = scipy.zeros(terms)
    k = 0
    for i in range(order+1):
        for j in range (i+1):
            AF[k] = (-1)**(i)*A[k]
            k += 1
    return AF


def RotateCoeffs(a, theta, order, verbose=False):
    '''Rotate axes of coefficients by theta degrees'''
    c = cos(radians(theta))
    s = sin(radians(theta))

    # First place in triangular layout
    at = scipy.zeros([order+1, order+1])
    k = 0
    for m in range(order+1):
        for n in range(m+1):
            at[m, n] = a[k]
            k+=1

    # Apply rotation
    atrotate = scipy.zeros([order+1, order+1])
    arotate = scipy.zeros([len(a)]) # Copy shape of a
    for m in range(order+1):
        for n in range(m+1):
            for mu in range(0, m-n+1):
                for j in range(m-n-mu, m-mu+1):
                    factor = (-1)**(m-n-mu)*choose(m-j, mu)*choose(j, m-n-mu)
                    cosSin = c**(j+2*mu-m+n)*s**(2*m-2*mu-j-n)
                    atrotate[m, n] = atrotate[m, n] + factor*cosSin*at[m, j]
                    if verbose: print(m, n, j, factor, 'cos^', j+2*mu-m+n, 'sin^', 2*m-2*mu-j-n, ' A', m, j)
    # Put back in linear layout
    k = 0
    for m in range(order+1):
        for n in range(m+1):
            arotate[k] = atrotate[m, n]
            k+=1

    return arotate


def TransCoeffs(A, a, b, c, d, order, verbose=False):
    '''Transform polynomial coefficients to allow for
    xp = a*x + b*y
    yp = c*x + d*y'''

    A1 = scipy.zeros((order+1, order+1))
    A2 = scipy.zeros((order+1, order+1))
    ncoeffs = (order+1)*(order+2)//2
    if verbose: print(ncoeffs, 'coefficients for order', order)
    AT = scipy.zeros((ncoeffs))

    # First place A in triangular layout
    k = 0
    for i in range(order+1):
        for j in range(i+1):
            A1[i, j] = A[k]
            k+=1

    for m in range(order+1):
        for n in range(m+1):
            if verbose: print('\nM, N', m, n)
            for mu in range(m-n+1):
                for j in range(m-n-mu, m-mu+1):
                    if verbose: print('J, MU', j, mu)
                    if verbose: print('Choose', m-j, mu, 'and', j, m-n-mu)
                    factor = choose(m-j, mu)*choose(j, m-n-mu)
                    A2[m, n] += factor*a**mu*b**(m-j-mu)*c**(m-n-mu)*d**(mu+j-m+n)*A1[m, j]
                    if verbose: print(m, j, ' Factor', factor)
    # Restore A2 to flat layout in AT
    k = 0
    for m in range(order+1):
        for n in range(m+1):
            AT[k] = A2[m, n]
            k += 1
    return AT


def TwoStep(A, B, a, b, order):
    """Change coefficients when
    xp = a[0] + a[1].x + a[2].y
    yp = b[0] + b[1].x + b[2].y"""
    terms = (order+1)*(order+2)//2
    A2 = scipy.zeros((order+1, order+1))
    B2 = scipy.zeros((order+1, order+1))

    k=0
    for i in range(order+1):
        for j in range(i+1):
            for alpha in range(i-j+1):
                for beta in range(i-j-alpha+1):
                    f1 = choose(i-j, alpha)*choose(i-j-alpha, beta)*a[0]**(i-j-alpha-beta)*a[1]**alpha*a[2]**beta
                    for gamma in range(j+1):
                        for delta in range(j-gamma+1):
                            f2 = choose(j, gamma)*choose(j-gamma, delta)*b[0]**(j-gamma-delta)*b[1]**gamma*b[2]**delta
                            A2[alpha+beta+gamma+delta, beta+delta] += A[k]*f1*f2
                            B2[alpha+beta+gamma+delta, beta+delta] += B[k]*f1*f2
            k += 1

    # Flatten A@ and B2
    k = 0
    Aflat = scipy.zeros(terms)
    Bflat = scipy.zeros(terms)
    for i in range(order+1):
        for j in range(i+1):
            Aflat[k] = A2[i, j]
            Bflat[k] = B2[i, j]
            k += 1
    return (Aflat, Bflat)


def TestTwoStep():
    A = scipy.array([10.0, 2.0, 0.1, 0.01, -0.02, 0.03])
    B = scipy.array([4.0, 1.8, 0.2, 0.02, 0.03, -0.02])
    a = scipy.array([1.0, 0.5, 0.1])
    b = scipy.array([2.0, 0.2, 0.6])
    print('\nA')
    triangle(A, 2)
    print('B')
    triangle(B, 2)
    print('a\n', a)
    print('b\n', b)
    (A2, B2) = TwoStep(A, B, a, b, 2)
    print('\nA2')
    triangle(A2, 2)
    print('B2')
    triangle(B2, 2)

    # Now do a test calculation
    (x, y) = (10, 5)
    xp = a[0] + a[1]*x + a[2]*y
    yp = b[0] + b[1]*x + b[2]*y
    print('x, y', x, y)
    print('xp, yp', xp, yp)

    u = poly(A, xp, yp, 2)
    v = poly(B, xp, yp, 2)
    up = poly(A2, x, y, 2)
    vp = poly(B2, x, y, 2)
    print('Two step', u, v)
    print('One step', up, vp)
    return
