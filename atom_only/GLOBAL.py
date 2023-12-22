from atoms_data import *
import numpy as np
from copy import copy

angular_momenta = ['s','p','d','f','g','h','i','j','k','l']
def orbit_transform(nl,string):
    """ Transform orbitals into strings<->tuples, e.g. (2,1)<->'2p'. """
    if string==True and type(nl)==type(''):
        return nl #'2p'->'2p'
    elif string==True:
        return '%i%s' %(nl[0],angular_momenta[nl[1]]) #(2,1)->'2p'
    elif string==False and type(nl)==type((2,1)):
        return nl      #(2,1)->(2,1)
    elif string==False:
        l=angular_momenta.index(nl[1])
        n=int(nl[0])
        return (n,l)  # '2p'->(2,1)


def list_states(maxn, maxl, occu):
    """ List all potential states {(n,l,'nl')}. """
    states = []
    for l in range(maxl+1):
        for n in range(1,maxn+1):
            nl = orbit_transform((n,l),string=True)
            if nl in occu:
                states.append((n,l,nl))
    return states

symbol = "C"

configuration = {}
valence = []
confinement = None
xc = "PW92"
convergence = {"density":1E-7,"energies":1E-7}
scalarrel = False
rmax = 100.0
nodegpts = 500
mix =0.2
itmax = 200
timing = False
verbose = False
txt = None
restart = None
write = None

# element data
data = copy( atoms_data[symbol] )

Z = data["Z"]
if valence == []:
    valence = copy( atoms_data[symbol]["valence_orbitals"] )

print("Symbol = ", symbol)
print("Valence = ", valence)

# more specific
occu = copy( atoms_data[symbol]["configuration"] )
nel_neutral = Z

assert sum(occu.values()) == nel_neutral

occu.update( configuration )
print("Occupations = ", occu)

nel = sum(occu.values())
print("nel = ", nel)

charge = nel_neutral - nel
print("charge = ", charge)


"""
if self.confinement==None:
    self.confinement_potential=ConfinementPotential("none")
else:
    self.confinement_potential=ConfinementPotential(**self.confinement)
"""


conf=None
nucl=None
exc=None


# technical stuff
maxl = 9
maxn = 9
plotr = {}
unlg = {}
Rnlg = {}
unl_fct = {}
Rnl_fct = {}
veff_fct = None
total_energy = 0.0

maxnodes = max( [n-l-1 for n,l,nl in list_states(maxn,maxl,occu)] )

rmin, rmax, N=( 1E-2/Z, rmax, (maxnodes+1)*nodegpts )

if scalarrel:
    print("Using scalar relativistic corrections.")

print("max %i nodes, %i grid points" %(maxnodes,N))
print("rmin = ", rmin, " rmax = ", rmax)

xgrid = np.linspace(0,np.log(rmax/rmin),N)
rgrid = rmin*np.exp(xgrid)

#grid = RadialGrid(rgrid)


veff = np.zeros(N)
VHartree = np.zeros(N)
vxc = np.zeros(N)
exc = np.zeros(N)

# make confinement and nuclear potentials
#conf = array([self.confinement_potential(r) for r in self.rgrid])
conf = np.zeros(N)

#nucl = array([self.V_nuclear(r) for r in self.rgrid])
nucl = np.zeros(N)

dens = np.zeros(N)