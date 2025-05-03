from math import pi

import numpy as np

from hotbit.io.fortran import fortran_readline

fileobj = "param/PTBP_Ni_N_H/Ni-H.skf"
symboli = "Ni"
symbolj = "H"

# Read Hamitonian and overlap data from DFTB-style .skf file.
# 
# Parameters:
# -----------
# fileobj:   filename of file-object to read from
# symboli:   chemical symbol of the first element
# symbolj:   chemical symbol of the second element

if isinstance(fileobj, str):
    fileobj = open(fileobj)

# Homoatomic interactions also contain element data
if symboli == symbolj:
    # First line contains grid spacing and number of grid points
    l = fortran_readline(fileobj)
    dx = l[0]
    n = l[1]
    n = int(n)

    # Contains self-energies, spin-polarization energies, Hubbard-U, ...
    l = fileobj.readline()
else:
    # First line contains grid spacing and number of grid points
    dx, n = fortran_readline(fileobj)
    n = int(n)

# The n information is sometimes wrong, better count while reading
#x = dx*np.arange(0, n)

HS = [ ]
l = fileobj.readline()
while l and l.strip() != 'Spline':
    if l.strip() == 'Spline':
        if i != n-1:
            raise RuntimeError('Spline keyword encountered reading tables '
                                'for %s-%s before reaching the end of the '
                                'HS table.' % ( symboli, symbolj ))
    else:
        HS += [ fortran_readline(l) ]

    l = fileobj.readline()

# don't care if not spline...
#if not l:
#    raise RuntimeError('Premature end-of-file: Keyword "Spline" not found '
#                       'for %s-%s.' % ( symboli, symbolj ))

HS = np.array(HS)
x = dx*np.arange(0, HS.shape[0])

#return x[0:HS.shape[0]], np.array(HS)
#return x, np.array(HS)