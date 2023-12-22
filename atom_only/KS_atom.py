import numpy as np
from PW92 import *
from GLOBAL import *

def V_nuclear(Z, r):
    return -Z/r


def guess_density(Z, rgrid):
    """ Guess initial density. """
    r2 = 0.02*Z # radius at which density has dropped to half; improve this!
    dens = np.exp( -rgrid/(r2/np.log(2)) )
    dens = dens/radialgrid_integrate( dens, rgrid,use_dV=True)*nel
    #pl.plot(self.rgrid,dens)
    return dens

def radialgrid_integrate(f, rgrid, use_dV=False):
    """
    Integrate function f (given with N grid points).
    int_rmin^rmax f*dr (use_dv=False) or int_rmin^rmax*f dV (use_dV=True)
    """
    N = len(rgrid)
    dr = rgrid[1:N] - rgrid[0:N-1]
    r0 = rgrid[0:N-1] + dr/2
    # first dV is sphere (treat separately), others are shells
    dV = 4*np.pi*r0**2 * dr
    dV *= (4*np.pi*rmax**3/3)/sum(dV)
    if use_dV:
        return ( (f[0:N-1]+f[1:N])*dV).sum()*0.5
    else:
        return ( (f[0:N-1]+f[1:N])*dr).sum()*0.5


def calculate_Hartree_potential():
    """
    Calculate Hartree potential.

    Everything is very sensitive to the way this is calculated.
    If you can think of how to improve this, please tell me!

    """
    #dV = self.grid.get_dvolumes()
    
    N = len(rgrid)

    dr = rgrid[1:N] - rgrid[0:N-1]
    r0 = rgrid[0:N-1] + dr/2
    dV = 4*np.pi*r0**2 * dr
    dV *= (4*np.pi*rmax**3/3)/sum(dV)
    #
    n0 = 0.5*( dens[1:] + dens[:-1] )
    n0 *= nel/sum(n0*dV)

    lo, hi = np.zeros(N), np.zeros(N)
    lo[0] = 0.0
    for i in range(1,N):
        lo[i] = lo[i-1] + dV[i-1]*n0[i-1]

    hi[-1] = 0.0
    for i in range(N-2,-1,-1):
        hi[i] = hi[i+1] + n0[i]*dV[i]/r0[i]

    for i in range(N):
        VHartree[i] = lo[i]/rgrid[i] + hi[i]

def calculate_veff():
    """ Calculate effective potential. """
    for i in range(N):
        vxc[i] = vxc_PW92(dens[i])
    #
    return nucl + VHartree + vxc + conf


def construct_coefficients(l, eps):
    c = 137.036
    c2 = np.ones(N)
    if scalarrel == False:
        c0 = -2*( 0.5*l*(l+1) + rgrid**2*(veff-eps) )
        c1 = -np.ones(N)
    else:
        print('SHOULD NOT PASS HERE !!!')
        # from Paolo Giannozzi: Notes on pseudopotential generation
        #ScR_mass = array([1 + 0.5*(eps-V)/c**2 for V in self.veff])
        #c0 = -l*(l+1) - 2*ScR_mass*self.rgrid**2*(self.veff-eps) - self.dveff*self.rgrid/(2*ScR_mass*c**2)
        #c1 = self.rgrid*self.dveff/(2*ScR_mass*c**2) - 1
    return c0, c1, c2

