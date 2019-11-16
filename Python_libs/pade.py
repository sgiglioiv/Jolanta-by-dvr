import numpy
import sys
from scipy import sqrt
from scipy.optimize import minimize


def genpade2_via_lstsq(nq, np, ns, xs, ys, rcond=1e-15, return_lists=False):
    """ 
    compute the coeffients of a generalized Pade np/nq/ns approximant using numpy.linalg.lstsq
    the lstsq function will take any number of equations, so this will
    work for any len(xs)=len(ys) 
    returned are the Pade coefficients Ps, Qs, and Ss in highest-power first order
    Ps, Qs, and Ss are returns as numpy polynomials: poly1d package 
    """
    N_coef = ns + np + nq + 2
    M_data = len(xs)
    A=numpy.zeros((M_data,N_coef))
    b=numpy.zeros(M_data)
    for k_data in range(M_data):
        i_coef = 0
        E=ys[k_data]
        E2=E*E
        #
        for iq in range(nq,1,-1): # if nq <= 2 this range is automatically empty
            A[k_data,i_coef] = E2*xs[k_data]**iq
            i_coef += 1
        if nq > 0:
            A[k_data,i_coef] = E2*xs[k_data]
            i_coef += 1
        #
        for ip in range(np,1,-1): # if np <= 2 this range is automatically empty
            A[k_data,i_coef] = E*xs[k_data]**ip
            i_coef += 1
        if np > 0:
            A[k_data,i_coef] = E*xs[k_data]
            i_coef += 1
        A[k_data,i_coef] = E
        i_coef += 1
        #
        for ks in range(ns,1,-1): # if ns <= 2 this range is automatically empty
            A[k_data,i_coef] = xs[k_data]**ks
            i_coef += 1
        if ns > 0:
            A[k_data,i_coef] = xs[k_data]
            i_coef += 1
        A[k_data,i_coef] = 1.0
        #
        b[k_data] = -E2 
    coefs, residual, rank, s = numpy.linalg.lstsq(A,b,rcond=rcond)
    Qs = numpy.array(list(coefs[:nq]) + [1])
    Ps = coefs[nq:nq+np+1]
    Ss = coefs[-ns-1:]
    if return_lists:
        return Qs, Ps, Ss
    else:
        return numpy.poly1d(Qs), numpy.poly1d(Ps), numpy.poly1d(Ss)

#
# return the two roots, and possibly the derivative of E+
# of a generalized quadratic: A(x)*E^2 + B(x)*E + C(x) = 0
# where A, B, and C are polynomials in x (numpy.poly1d)
#
def E_minus(x, A, B, C):
    """ minus root """
    return -0.5 * (B(x) - sqrt(B(x)*B(x)-4*A(x)*C(x)) ) / A(x)
    
def E_plus(x, A, B, C, der=0):
    """ plus root """
    a=A(x)
    b=B(x)
    c=C(x)
    s=sqrt(b*b - 4*a*c)
    if a == 0:
        E = -b/c
    else:
        E = -0.5 * (b + s) / a
    if der == 0:
        return E
    else:
        ap = A.deriv()(x)
        bp = B.deriv()(x)
        cp = C.deriv()(x)
        sp = (b*bp - 2*ap*c - 2*a*cp) / s
        dEdx = (ap/a*(b-s) + sp - bp)  /  (2*a)
        return E, dEdx
    
def E_lower(x, A, B, C):
    """ lower branch of a generalized Pade approximant """
    if A(x) < 0:
        return -0.5 * (B(x) - numpy.sqrt(B(x)*B(x)-4*A(x)*C(x)) ) / A(x)
    else:
        return -0.5 * (B(x) + numpy.sqrt(B(x)*B(x)-4*A(x)*C(x)) ) / A(x) 
    
def E_upper(x, A, B, C):
    """ upper branch of a generalized Pade approximant """
    if A(x) >= 0:
        return -0.5 * (B(x) - numpy.sqrt(B(x)*B(x)-4*A(x)*C(x)) ) / A(x)
    else:
        return -0.5 * (B(x) + numpy.sqrt(B(x)*B(x)-4*A(x)*C(x)) ) / A(x)

"""

         Standard Pade             

"""

    
def pade_via_lstsq(np, nq, xs, ys, rcond=1e-15):
    """ 
    compute the coeffients of a Pade np/np approximant using numpy.linalg.lstsq
    the lstsq function will take any number of equations, so this should
    work for any len(xs)=len(ys) 
    returned are the Pade coefficients ps and qs in highest-power first order
    """
    N_coef = np + nq + 1
    M_data = len(xs)
    A=numpy.zeros((M_data,N_coef))
    b=numpy.zeros(M_data)
    for k_data in range(M_data):
        i_coef = 0
        for ip in range(np,1,-1): # if np >= 2 this range is automatically empty
            A[k_data,i_coef] = xs[k_data]**ip
            i_coef += 1
        A[k_data,i_coef] = xs[k_data]
        i_coef += 1
        A[k_data,i_coef] = 1.0
        i_coef += 1
        for iq in range(nq,1,-1): # if nq >= 2 this range is automatically empty
            A[k_data,i_coef] = -ys[k_data] * xs[k_data]**iq
            i_coef += 1
        A[k_data,i_coef] = -ys[k_data] * xs[k_data]
        b[k_data] = ys[k_data]
    coefs, residual, rank, s = numpy.linalg.lstsq(A,b,rcond=rcond)
    ps = coefs[:np+1]
    qs = numpy.array(list(coefs[np+1:np+nq+1]) + [1])
    return ps, qs



def eval_pade(x, ps, qs):
    """ evaluate a standard Pade approximant """
    return numpy.polyval(ps,x)/numpy.polyval(qs,x)


