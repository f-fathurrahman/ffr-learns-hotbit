from math import log, pi, sqrt
from copy import copy

import numpy as np

from box.data import data as box_data
from hotbit.io.fortran import fortran_readline

fileobj = "param/PTBP_Ni_N_H/Ni-Ni.skf" # should be homonuclear species
symbol = "Ni"

if isinstance(fileobj, str):
    fileobj = open(fileobj)

# First line contains grid spacing and number of grid points, ignore
fileobj.readline()

# Self-energies, spin-polarization energies, Hubbard-U
eself = [ 0.0 ]*3
U = [ 0.0 ]*3
q = [ 0.0 ]*3
eself[0], eself[1], eself[2], espin, \
    U[0], U[1], U[2], q[0], q[1], q[2], *_ = fortran_readline(fileobj)

# Initialize data from database
data = copy(box_data[symbol])

data['epsilon'] = { }
for orbital in data['valence_orbitals']:
    if 's' in orbital:
        data['epsilon'][orbital] = eself[2]
    elif 'p' in orbital:
        data['epsilon'][orbital] = eself[1]
    elif 'd' in orbital:
        data['epsilon'][orbital] = eself[0]

# Apparently, only the last U is used
data['U'] = U[2]
data['FWHM'] = sqrt(8*log(2.0)/pi)/data['U']

energies = []            
for orbital in data['valence_orbitals']:            
    eps = data['epsilon'][orbital]
    if   's' in orbital: n=1
    elif 'p' in orbital: n=3
    elif 'd' in orbital: n=5
    energies.extend( [eps]*n )                
data['onsite_energies'] = energies
data['nr_basis_orbitals'] = len(energies)
data['valence_energies'] = np.array(energies, dtype=float)

data['comment'] = None

# .skf files do not contain basis information
#return data, None
