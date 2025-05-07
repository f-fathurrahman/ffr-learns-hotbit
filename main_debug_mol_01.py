from ase import *
from my_hotbit import MyHotbit
from ase.build import molecule

atoms = molecule('C6H6')
atoms.set_cell([10.0, 10.0, 10.0])

calc = MyHotbit(SCC=False, txt="-")
atoms.set_calculator(calc)

# This is called in solve_ground_state
calc.el.update_geometry(atoms)

from my_occupations import Occupations

if calc.st.nk == None:
    physical = calc.st.calc.get('physical_k') # call directly calc?
    calc.st.nk, calc.st.k, calc.st.kl, calc.st.wk = calc.st.setup_k_sampling(
        calc.st.calc.get('kpts'),
        physical=physical,
        rs=calc.st.calc.get('rs')
    )
    width = calc.st.calc.get('width')
    calc.st.occu = Occupations(calc.st.calc.el.get_number_of_electrons(), width, calc.st.wk)

from debug_get_matrices_01 import debug_get_matrices
H0, S, dH0, dS = debug_get_matrices(calc.ia)


"""
Etot = atoms.get_potential_energy()
print("Etot = ", Etot)

forces = atoms.get_forces()
print("forces = ", forces)
"""