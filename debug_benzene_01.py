from ase import *
from my_hotbit import MyHotbit
from ase.build import molecule

atoms = molecule('C6H6')
atoms.set_cell([10.0, 10.0, 10.0])

calc = MyHotbit(SCC=False, txt="-")
atoms.set_calculator(calc)

Etot = atoms.get_potential_energy()
print("Etot = ", Etot)

forces = atoms.get_forces()
print("forces = ", forces)