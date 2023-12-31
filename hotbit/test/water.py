from hotbit import Hotbit 
from my_ase import *
from box import Atoms
from box.md import quench
from hotbit.test.misc import default_param
#import cgitb; cgitb.enable()

h2o=read('H2O.xyz')
h2o.center(vacuum=5)

calc=Hotbit(verbose=True,SCC=True,verbose_SCC=False,txt='h2o.cal',**default_param)
h2o.set_calculator(calc)
print(h2o.get_potential_energy())
#print h2o.get_forces()

calc.timer.summary()
print(calc.timer.get_timings())
