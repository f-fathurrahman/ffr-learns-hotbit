from math import pi

import numpy as np

from box.data import data as box_data
from hotbit.io.fortran import fortran_readline

fileobj = "param/PTBP_Ni_N_H/Ni-Ni.skf" # should be homonuclear species
symbol = "Ni"


# Read repulsion from DFTB-style .skf file.
# The repulsion in the .skf-file consists of an exponential and a spline
# part. Both will be converted to a single table.
#
# Parameters:
# -----------
# fileobj:   filename of file-object to read from
# rep_x0:    minimum value of the repulsion table
# rep_dx:    step size for discretization

# default values
rep_x0 = 0.1
rep_dx = 0.005

if isinstance(fileobj, str):
    fileobj = open(fileobj)

l = fileobj.readline()
while l and l.strip() != 'Spline':
    l = fileobj.readline()

if not l:
    # no spline
    fileobj.seek(0)
    dx, n = fortran_readline(fileobj) #ignore
    items = fortran_readline(fileobj)
    if len(items)<20: # in homonuclear tables one element info came before repulsion
        items = fortran_readline(fileobj) 
    cutoff = items[9]
    c = [0,0] + items[1:9]         
    x = np.linspace(rep_x0, cutoff, int((cutoff-rep_x0)/rep_dx)+1)
    y = np.zeros_like(x)
    for i in range(2,10):
        y = y + c[i]*(cutoff-x)**i
    
    #raise RuntimeError('Could not find "Spline" keyword when reading '
    #                   'repulsion.')

else:
    n, cutoff = fortran_readline(fileobj)
    n = int(n)

    # Coefficients for the exponential
    c1, c2, c3 = fortran_readline(fileobj)

    x1, x2, splc1, splc2, splc3, splc4 = fortran_readline(fileobj)

    x = np.linspace(rep_x0, cutoff, int((cutoff-rep_x0)/rep_dx)+1)
    y = c3 + np.exp(c2-c1*x)

    i0 = np.searchsorted(x, x1)+1
    for j in range(n-1):
        if j > 0:
            last_x2 = x2
            x1, x2, splc1, splc2, splc3, splc4 = fortran_readline(fileobj)
            assert x1 == last_x2
        i1 = np.searchsorted(x, x2)+1
        y[i0:i1] = \
            splc1 + \
            splc2 * (x[i0:i1]-x1) + \
            splc3 * (x[i0:i1]-x1)**2 + \
            splc4 * (x[i0:i1]-x1)**3
        i0 = i1

    # The last entry is a fifth-order polynomial
    last_x2 = x2
    x1, x2, splc1, splc2, splc3, splc4, splc5, splc6 = fortran_readline(fileobj)
    assert x1 == last_x2

    i1 = np.searchsorted(x, x2)+1
    y[i0:i1] = \
        splc1 + \
        splc2 * (x[i0:i1]-x1) + \
        splc3 * (x[i0:i1]-x1)**2 + \
        splc4 * (x[i0:i1]-x1)**3 + \
        splc5 * (x[i0:i1]-x1)**4 + \
        splc6 * (x[i0:i1]-x1)**5

# return x, y
