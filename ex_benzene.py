from ase import *
from hotbit import *
from ase.build import molecule

atoms = molecule('C6H6')
atoms.set_cell([10.0, 10.0, 10.0])

#calc = Hotbit(SCC=True, width=0.05, txt=None)
#calc = Hotbit(SCC=True, width=0.05, txt=None, 
#    mixer={'name':'Anderson','mixing_constant':0.1, 'memory':4}
#)
calc = Hotbit(SCC=True, width=0.05, txt='LOG_benzene', 
    mixer={'name': 'Pulay', 'mixing_constant': 0.2, 'memory':4}
)
#calc = Hotbit(SCC=False, width=0.05, txt=None)

print("\nENTER Atoms.set_calculator")
atoms.set_calculator(calc)
print("EXIT  Atoms.set_calculator")

Etot = atoms.get_potential_energy()
print("Etot = ", Etot)