from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import splrep, splev

from KS_atom import *
from GLOBAL import *

exec(open("shoot.py").read())
exec(open("solve_eigenstates.py").read())

exec(open("calculate_density.py").read())
exec(open("calculate_energies.py").read())

exec(open("SplineFunction.py").read())
exec(open("FastSplineFunction.py").read())
exec(open("Function.py").read())

"""
Solve the self-consistent potential.
"""
print("\nStart iteration...")

enl = {}
d_enl = {}
for n,l,nl in list_states(maxn,maxl,occu):
    enl[nl] = 0.0
    d_enl[nl] = 0.0

dens = guess_density(Z, rgrid)
nucl = V_nuclear(Z, rgrid)
veff = nucl + conf

calculate_Hartree_potential()

# SCF start here

print("mix = ", mix)
print("itmax = ", itmax)

for it in range(itmax):

    veff = mix*calculate_veff() + (1-mix)*veff

    d_enl_max, itmax_eigen = solve_eigenstates(it)

    # Evaluate convergence in electron density
    dens0 = dens.copy()
    dens = calculate_density()
    diff = radialgrid_integrate( np.abs(dens-dens0), rgrid, use_dV=True )

    if diff < convergence["density"] and d_enl_max < convergence["energies"] and it > 5:
        break

    calculate_Hartree_potential()

    print("%3d %18.10e" % (it, diff))

    if it == itmax-1:
        raise RuntimeError("Density not converged in %i iterations" %(it+1))


calculate_energies()
print("converged in %i iterations" % it)
print("%9.4f electrons, should be %9.4f" %
      (radialgrid_integrate(dens,rgrid,use_dV=True), nel))


for n,l,nl in list_states(maxn,maxl,occu):
    Rnl_fct[nl] = Function("spline", rgrid, Rnlg[nl])
    unl_fct[nl] = Function("spline", rgrid, unlg[nl])

plt.clf()
plt.plot( rgrid, Rnlg["1s"], marker="o", label="Rnlg" )
plt.plot( rgrid, unlg["1s"], marker="o", label="unlg" )
plt.xlim(0,10)
plt.legend()
plt.grid()
plt.savefig("Radial_1s.png")


print("Program ended normally")
