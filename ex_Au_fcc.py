from hotbit import Hotbit
from ase.lattice.cubic import FaceCenteredCubic
import numpy as np

pi = np.pi

d = 4.08
atoms = FaceCenteredCubic(
    'Au',
    latticeconstant=d,
    directions=((0,1,1),(1,0,1),(1,1,0)),
    align=False
)
                          
calc = Hotbit(SCC=False, kpts=(8,8,8), txt="-")
atoms.set_calculator(calc)

Etot = atoms.get_potential_energy()
print("Etot = ", Etot)
