from my_ase import *
from hotbit import *
from my_ase.build import molecule

atoms = molecule('C6H6')
atoms.set_cell([10.0, 10.0, 10.0])

# This will only set some parameters, no significant computations
calc = Hotbit(SCC=False, width=0.05, txt=None)
el = Elements(calc, atoms)
el.set_atoms(atoms)
print(el.greetings()) 