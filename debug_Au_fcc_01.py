from my_hotbit import MyHotbit
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
                          
calc = MyHotbit(SCC=False, kpts=(8,8,8), txt="LOG_Au_fcc")

atoms.set_calculator(calc) # initialize will be called here

Etot = atoms.get_potential_energy()
print("Etot = ", Etot)
