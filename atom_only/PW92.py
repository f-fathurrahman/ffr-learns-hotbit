import numpy as np
from math import sqrt, log, pi

def exc_PW92(n,der=0):
    """ Exchange-correlation with electron density n. """
    SMALL = 1E-90
    if n < SMALL:
        return 0.0
    else:
        return e_x_PW92(n, der=der) + e_corr_PW92(n,der=der)

def e_x_PW92(n,der=0):
    """ Exchange. """
    if der==0:
        return -3.0/4*(3*n/pi)**(1.0/3)
    elif der==1:
        return -3.0/(4*pi)*(3*n/pi)**(-2.0/3)


def e_corr_PW92(n,der=0):
    """ Correlation energy. """
    a1 = 0.21370
    c0 = 0.031091
    c1 = 0.046644
    b1 = 1.0/2.0/c0*np.exp(-c1/2.0/c0)
    b2 = 2*c0*b1**2
    b3 = 1.6382
    b4 = 0.49294
    rs = (3.0/(4*pi*n))**(1.0/3)
    aux = 2*c0*( b1*sqrt(rs)+b2*rs+b3*rs**(3.0/2)+b4*rs**2 )
    if der==0:
        return -2*c0*(1+a1*rs)*log(1+aux**-1)
    elif der==1:
        return ( -2*c0*a1*log(1+aux**-1) \
               -2*c0*(1+a1*rs)*(1+aux**-1)**-1*(-aux**-2)\
               *2*c0*(b1/(2*sqrt(rs))+b2+3*b3*sqrt(rs)/2+2*b4*rs) )*( -(4*pi*n**2*rs**2)**-1 )

def vxc_PW92(n):
    """ Exchange-correlation potential (functional derivative of exc). """
    eps = 1E-9*n  # ???????
    SMALL = 1E-90
    if n < SMALL:
        return 0.0
    else:
        return exc_PW92(n) + n*exc_PW92(n,der=1)

