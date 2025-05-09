# Copyright (C) 2008 NSC Jyvaskyla
# Please see the accompanying LICENSE file for further information.

"""
    ASE-calculator interface for HOTBIT.
    (Hybrid Open-source Tight-Binding Tool)

"""
from __future__ import print_function
import os
import glob
import sys
import numpy as np

from my_ase.units import Bohr, Hartree

from hotbit.auxil import k_to_kappa_points
from box.timing import Timer

from my_elements import Elements

from my_interactions import Interactions

# This is probably can be ignored for the moment
from my_environment import Environment

from my_pairpotential import PairPotential

from my_repulsion import Repulsion

from my_states import States

from my_grids import Grids

from hotbit.version import hotbit_version
from hotbit.analysis import MullikenAnalysis
from hotbit.analysis import MullikenBondAnalysis
from hotbit.analysis import DensityOfStates
from hotbit.output import Output
from hotbit.vdw import setup_vdw

from box.mix import broaden
import box.mix as mix
from time import time

hbar = 0.02342178268


class MyHotbit(Output):
    
    def __init__(self,parameters=None,
                      elements=None,
                      tables=None,
                      verbose=False,
                      charge=0.0,
                      SCC=True,
                      kpts=(1,1,1),
                      rs='kappa',
                      physical_k=True,
                      maxiter=50,
                      gamma_cut=None,
                      txt=None,
                      verbose_SCC=False,
                      width=0.02,
                      mixer=None,
                      coulomb_solver=None,
                      charge_density='Gaussian',
                      vdw=False,
                      vdw_parameters=None,
                      internal={}):

        print("\n<div> ENTER Hotbit.__init__\n")

        import os

        # Convert gamma_cut from Bohr to angstrom
        if gamma_cut != None:
            print("gamma_cut = ", gamma_cut)
            gamma_cut = gamma_cut/Bohr
            print("gamma_cut = ", gamma_cut)

        self.__dict__ = {
            'parameters':parameters,
            'elements':elements,
            'tables':tables,
            'verbose':verbose,
            'charge':charge,
            'width':width/Hartree,
            'SCC':SCC,
            'kpts':kpts,
            'rs':rs,
            'physical_k':physical_k,
            'maxiter':maxiter,
            'gamma_cut':gamma_cut,
            'vdw':vdw,
            'vdw_parameters':vdw_parameters,
            'txt':txt,
            'verbose_SCC':verbose_SCC,
            'mixer':mixer,
            'coulomb_solver':coulomb_solver,
            'charge_density':charge_density,
            'internal':internal
        }

        if parameters != None:
            os.environ['HOTBIT_PARAMETERS'] = parameters

        self.init = False
        self.notes = []
        self.dry_run = '--dry-run' in sys.argv
        internal0 = {
            'sepsilon':0.,               # add this to the diagonal of S to avoid LAPACK error in diagonalization
            'tol_imaginary_e': 1E-13,    # tolerance for imaginary band energy
            'tol_mulliken':1E-5,         # tolerance for mulliken charge sum deviation from integer
            'tol_eigenvector_norm':1E-6, # tolerance for eigenvector norm for eigensolver
            'symop_range':5              # range for the number of symmetry operations in all symmetries
        }              
        internal0.update(internal)
        for key in internal0:
            self.set(key,internal0[key])
        #self.set_text(self.txt)
        #self.timer=Timer('Hotbit',txt=self.get_output())
        print("\n</div> EXIT Hotbit.__init__\n")


    def __del__(self):

        if self.get('SCC'):
            try:
                print(self.st.solver.get_iteration_info())
                self.txt.flush()
            except:
                pass
        if len(self.notes)>0:
            print('Notes and warnings:')
            for note in self.notes:
                print(note)
        if self.init:
            self.timer.summary()
            Output.__del__(self)


    def write_electronic_data(self,filename,keys=None):

        data = {}
        data['N'] = self.el.N
        data['norb'] = self.st.norb
        data['charge'] = self.get('charge')
        data['nelectrons'] = self.el.get_number_of_electrons()
        data['erep'] = self.rep.get_repulsive_energy()
        data['ecoul'] = self.get_coulomb_energy(self.el.atoms)
        data['ebs'] = self.get_band_structure_energy(self.el.atoms)
        data['epot'] = self.get_potential_energy(self.el.atoms)
        data['forces'] = self.get_forces(self.el.atoms)
        data['symbols'] = self.el.symbols
        data['e'] = self.st.e
        data['occ'] = self.st.f
        data['nk'] = self.st.nk
        data['k'] = self.st.k
        data['wk'] = self.st.wk
        data['dq'] = self.st.mulliken()
        data['gap'], data['gap_prob'] = self.get_energy_gap()
        data['dose'], data['dos'] = self.get_density_of_states(False)

        for key in list(data.keys()):
            if keys!=None and key not in keys:
                del data[key]
        import pickle
        f = open(filename, 'w')
        pickle.dump(data,f)
        f.close()


    def set(self,key,value):
        if key == 'txt':
            self.set_text(value)
        elif self.init==True and key not in ['charge']:
            raise AssertionError('Parameters cannot be set after initialization.')
        else:
            self.__dict__[key]=value


    def get_atoms(self):
        atoms = self.el.atoms.copy()
        atoms.set_calculator(self)
        return atoms


    def add_note(self,note):
        self.notes.append(note)


    def greetings(self):
        from time import asctime
        from os import uname
        from os.path import abspath, curdir
        from os import environ

        self.version = hotbit_version
        #print('\n\n\n\n\n')
        #
        print('\n')
        #
        print(' _           _    _     _ _')
        print('| |__   ___ | |_ | |__ |_| |_')
        print('|  _ \ / _ \|  _||  _ \| |  _|')
        print('| | | | ( ) | |_ | ( ) | | |_')
        print('|_| |_|\___/ \__|\____/|_|\__|  ver.',self.version)
        print('Distributed under GNU GPL; see %s' %environ.get('HOTBIT_DIR')+'/LICENSE')
        print('Date:',asctime())
        dat=uname()
        print('Nodename:', dat[1])
        print('Arch:', dat[4])
        print('Dir:', abspath(curdir))
        print('System:', self.el.get_name())
        print('       Charge=%4.1f' % self.charge)
        print('       Container', self.el.container_info())
        print('Symmetry operations (if any):')
        rs = self.get('rs')
        kpts = self.get('kpts')
        M = self.el.get_number_of_transformations()
        for i in range(3):
            print('       %i: pbc=' %i, self.el.atoms.get_pbc()[i], end=' ')
            if type(kpts)==type([]):
                print(', %s-points=%i, M=%.f' %(rs,len(kpts),M[i]))
            else:
                print(', %s-points=%i, M=%.f' %(rs,kpts[i],M[i]))
        print('Electronic temperature:', self.width*Hartree,'eV')
        mixer = self.st.solver.mixer
        print('Mixer:', mixer.get('name'), 'with memory =', mixer.get('memory'), ', mixing parameter =', mixer.get('beta'))
        print(self.el.greetings())
        print(self.ia.greetings())
        print(self.rep.greetings())
        if self.pp.exists():
            print(self.pp.greetings())


    def out(self,text):
        print(text)
        self.txt.flush()


    def set_text(self,txt):
        """ Set up the output file. """
        if txt=='-' or txt=='null':
            self.txt = open('/dev/null','w')
        elif hasattr(txt, 'write'):
            self.txt = txt
        elif txt is None:
            from sys import stdout
            self.txt=stdout
        else:
            self.txt=open(txt,'a')
        # check if the output of timer must be changed also
        if 'timer' in self.__dict__:
            self.timer.txt = self.get_output()


    def get(self,arg=None):

        if arg==None:
            return self.__dict__
        else:
            return self.__dict__[arg]


    def memory_estimate(self):

        if self.st.nk>1:
            number = 16. #complex
        else:
            number = 8. #real
        M = self.st.nk*self.st.norb**2*number
        #     H   S   dH0   dS    wf  H1  dH   rho rhoe
        mem = M + M + 3*M + 3*M + M + M + 3*M + M + M
        print('Memory consumption estimate: > %.2f GB' %(mem/1E9))
        self.txt.flush()
        if self.dry_run:
            raise SystemExit


    def solve_ground_state(self,atoms):

        print("\n<div> ENTER Hotbit.solve_ground_state")
        if not self.init:
            assert type(atoms) != type(None)
            print("Doing initialization")
            self._initialize(atoms)
        #
        #
        if type(atoms) == type(None):
            pass
        elif self.calculation_required(atoms, 'ground state'):
            #
            self.el.update_geometry(atoms)
            t0 = time()
            self.st.solve()
            self.el.set_solved('ground state')
            t1 = time()
            self.flags['Mulliken'] = False
            self.flags['DOS'] = False
            self.flags['bonds'] = False
            if self.verbose:
                print("Solved in %0.2f seconds" % (t1-t0), file=self.get_output())
            #if self.get('SCC'):
            #    atoms.set_charges(-self.st.get_dq())
        else:
            print("Will not do anything")
            print("Because calculation_required = ", self.calculation_required(atoms, 'ground state'))
            pass
        print("\n</div> EXIT Hotbit.solve_ground_state\n")


    def _initialize(self, atoms):

        print("\n<div> ENTER Hotbit._initialize\n")
        if not self.init:
            self.set_text(self.txt)
            self.timer = Timer('Hotbit',txt=self.get_output())
            self.start_timing('initialization')
            self.el = Elements(self, atoms)
            self.ia = Interactions(self)
            self.st = States(self)
            self.rep = Repulsion(self)
            self.pp = PairPotential(self)
            if self.get('vdw'):
                if self.get('vdw_parameters') is not None:
                    self.el.update_vdw(self.get('vdw_parameters'))
                setup_vdw(self)
            self.env = Environment(self)
            pbc = atoms.get_pbc()
            # FIXME: gamma_cut -stuff
            #if self.get('SCC') and np.any(pbc) and self.get('gamma_cut')==None:
            #    raise NotImplementedError('SCC not implemented for periodic systems yet (see parameter gamma_cut).')
            if np.any(pbc) and abs(self.get('charge'))>0.0 and self.get('SCC'):
                raise AssertionError('Charged system cannot be periodic.')
            self.flush()
            self.flags = {}
            self.flags['Mulliken'] = False
            self.flags['DOS'] = False
            self.flags['bonds'] = False
            self.flags['grid'] = False
            self.stop_timing('initialization')
        self.el.set_atoms(atoms)

        if not self.init:
            print("Will set init to true")
            self.init = True
            print("Will call Hotbit.greetings")
            self.greetings()
        print("\n</div> EXIT Hotbit._initialize\n")


    def calculation_required(self,atoms,quantities):
        """ Check if a calculation is required.

        Check if the quantities in the quantities list have already been calculated
        for the atomic configuration atoms. The quantities can be one or more of:
        'ground state', 'energy', 'forces', 'magmoms', and 'stress'.
        """
        return self.el.calculation_required(atoms,quantities)


    def get_potential_energy(self, atoms, force_consistent=False):
        """ Return the potential energy of present system. """
        print("\n<div> ENTER Hotbit.get_potential_energy\n")
        #
        if force_consistent:
            raise NotImplementedError
        if self.calculation_required(atoms,['energy']):
            #
            self.solve_ground_state(atoms)
            #
            self.start_timing('energy')
            #
            ebs = self.get_band_structure_energy(atoms)
            ecoul = self.get_coulomb_energy(atoms)
            erep = self.rep.get_repulsive_energy()
            epp = self.pp.get_energy()
            #
            self.epot = ebs + ecoul + erep + epp - self.el.efree*Hartree
            #
            self.stop_timing('energy')
            self.el.set_solved('energy')
        print("\n</div> EXIT  Hotbit.get_potential_energy\n")
        return self.epot.copy()


    def get_forces(self,atoms):
        """
        Return forces (in eV/Angstrom)

        Ftot = F(band structure) + F(coulomb) + F(repulsion).
        """
        if self.calculation_required(atoms,['forces']):
            self.solve_ground_state(atoms)
            self.start_timing('forces')
            fbs = self.st.get_band_structure_forces()
            frep = self.rep.get_repulsive_forces()
            fcoul = self.st.es.gamma_forces() #zero for non-SCC
            fpp = self.pp.get_forces()
            self.stop_timing('forces')
            self.f = (fbs+frep+fcoul+fpp)*(Hartree/Bohr)
            self.el.set_solved('forces')
        return self.f.copy()


    def get_band_energies(self, kpts=None, shift=True, rs='kappa', h1=False):

        if kpts is None:
            e = self.st.e * Hartree
        else:
            if rs=='k':
                klist = k_to_kappa_points(kpts,self.el.atoms)
            elif rs=='kappa':
                klist = kpts
            e = self.st.get_band_energies(klist,h1)*Hartree

        if shift:
            return e-self.get_fermi_level()
        else:
            return e


    def get_stress(self,atoms):
        self.solve_ground_state(atoms)
        # TODO: ASE needs an array from this method, would it be proper to
        # somehow inform that the stresses are not calculated?
        return np.zeros((6,))


    def get_charge(self):
        """ Return system's total charge. """
        return self.get('charge')


    def get_eigenvalues(self):
        return self.st.get_eigenvalues()*Hartree


    def get_energy_gap(self):

        eigs = (self.get_eigenvalues() - self.get_fermi_level()).flatten()
        occ = self.get_occupations().flatten()
        ehi, elo = 1E10, -1E10
        for e,f in zip(eigs,occ):
            if elo < e <= 0.0:
                elo = e
                flo = f
            elif 0.0 < e < ehi:
                ehi = e
                fhi = f
        return ehi-elo, (flo-fhi)/2


    def get_state_indices(self, state):
        """
        Return the k-point index and band index of given state.

        parameters:
        -----------
        state:    'HOMO', or 'LUMO'

                  HOMO is the first state below Fermi-level.
                  LUMO is the first state above Fermi-level.
        """
        eigs = (self.get_eigenvalues() - self.get_fermi_level()).flatten()
        if state=='HOMO':
            k,a = np.unravel_index(np.ma.masked_array(eigs,eigs>0.0).argmax(),(self.st.nk,self.st.norb))
        if state=='LUMO':
            k,a = np.unravel_index(np.ma.masked_array(eigs,eigs<0.0).argmin(),(self.st.nk,self.st.norb))
        return k,a


    def get_occupations(self):
        #self.solve_ground_state(atoms)
        return self.st.get_occupations()


    def get_band_structure_energy(self,atoms):
        if self.calculation_required(atoms, ['ebs']):
            self.solve_ground_state(atoms)
            self.ebs = self.st.get_band_structure_energy()*Hartree
            self.el.set_solved('ebs')
        return self.ebs


    def get_coulomb_energy(self,atoms):
        if self.calculation_required(atoms,['ecoul']):
            self.solve_ground_state(atoms)
            self.ecoul = self.st.es.coulomb_energy()*Hartree
            self.st
        return self.ecoul


    # some not implemented ASE-assumed methods
    def get_fermi_level(self):
        """
        Return the Fermi-energy (chemical potential) in eV.
        """
        return self.st.occu.get_mu() * Hartree


    def set_atoms(self,atoms):
        """ Initialize the calculator for given atomic system. """
        print("\n<div> ENTER Hotbit.set_atoms\n")
        if self.init==True and atoms.get_chemical_symbols()!=self.el.atoms.get_chemical_symbols():
            raise RuntimeError('Calculator initialized for %s. Create new calculator for %s.'
                               %(self.el.get_name(),mix.parse_name_for_atoms(atoms)))
        else:
            self._initialize(atoms)
        print("\n</div> EXIT Hotbit.set_atoms\n")


    def get_occupation_numbers(self,kpt=0):
        """ Return occupation numbers for given k-point index. """
        return self.st.f[kpt].copy()


    def get_number_of_bands(self):
        """ Return the total number of orbitals. """
        return self.st.norb


    def start_timing(self, label):
        self.timer.start(label)


    def stop_timing(self, label):
        self.timer.stop(label)


    #
    #    various analysis methods
    #
    def get_dielectric_function(self,width=0.05,cutoff=None,N=400):
        """
        Return the imaginary part of the dielectric function for non-SCC.

        Note: Uses approximation that requires that the orientation of
              neighboring unit cells does not change much.
              (Exact for Bravais lattice.)

        See, e.g., Marder, Condensed Matter Physics, or
        Popov New J. Phys 6, 17 (2004)

        parameters:
        -----------
        width:     energy broadening in eV
        cutoff:    cutoff energy in eV
        N:         number of points in energy grid

        return:
        -------
        e[:], d[:,0:2]
        """
        self.start_timing('dielectric function')
        width = width/Hartree
        otol = 0.05 # tolerance for occupations
        if cutoff==None:
            cutoff = 1E10
        else:
            cutoff = cutoff/Hartree

        st = self.st
        nk, e, f, wk = st.nk, st.e, st.f, st.wk
        ex, wt = [], []
        for k in range(nk):
            wf = st.wf[k]
            wfc = wf.conjugate()
            dS = st.dS[k].transpose((0,2,1))
            ek = e[k]
            fk = f[k]
            kweight = wk[k]
            # electron excitation ka-->kb; restrict the search:
            bmin = list(fk<2-otol).index(True)
            amin = list(ek>ek[bmin]-cutoff).index(True)
            amax = list(fk<otol).index(True)
            for a in range(amin,amax+1):
                bmax = list(ek>ek[a]+cutoff).index(True)
                for b in range(max(a+1,bmin),bmax+1):
                    de = ek[b]-ek[a]
                    df = fk[a]-fk[b]
                    if df<otol:
                        continue
                    # P = < ka | P | kb >
                    P = 1j*hbar*np.dot(wfc[a],np.dot(dS,wf[b]))
                    ex.append( de )
                    wt.append( kweight*df*np.abs(P)**2 )

        ex, wt = np.array(ex), np.array(wt)
        cutoff = min( ex.max(),cutoff )
        y = np.zeros((N,3))
        for d in range(3):
            # Lorenzian should be used, but long tail would bring divergence at zero energy
            x,y[:,d] = broaden( ex,wt[:,d],width,'gaussian',N=N,a=width,b=cutoff )
            y[:,d] = y[:,d]/x**2
        const = (4*np.pi**2/hbar)
        self.stop_timing('dielectric function')
        return x*Hartree, y*const #y also in eV, Ang

    #
    #   grid stuff
    #
    def set_grid(self,h=0.2,cutoff=3.0):
        if self.calculation_required(self.el.atoms,['energy']):
            raise AssertionError('Electronic structure is not solved yet!')
        if self.flags['grid']==False:
            self.gd = Grids(self,h,cutoff)
            self.flags['grid']=True


    def get_grid_basis_orbital(self,I,otype,k=0,pad=True):
        """
        Return basis orbital on grid.

        parameters:
        ===========
        I:     atom index
        otype: orbital type ('s','px','py',...)
        k:     k-point index (basis functions are really the extended
               Bloch functions for periodic systems)
        pad:   padded edges in the array
        """
        if self.flags['grid']==False:
            raise AssertionError('Grid needs to be set first by method "set_grid".')
        return self.gd.get_grid_basis_orbital(I,otype,k,pad)


    def get_grid_wf(self,a,k=0,pad=True):
        """
        Return eigenfunction on a grid.

        parameters:
        ===========
        a:     state (band) index
        k:     k-vector index
        pad:   padded edges
        """
        if self.flags['grid']==False:
            raise AssertionError('Grid needs to be set first by method "set_grid".')
        return self.gd.get_grid_wf(a,k,pad)


    def get_grid_wf_density(self,a,k=0,pad=True):
        """
        Return eigenfunction density.

        Density is not normalized; accurate quantitative analysis
        on this density are best avoided.

        parameters:
        ===========
        a:     state (band) index
        k:     k-vector index
        pad:   padded edges
        """
        if self.flags['grid']==False:
            raise AssertionError('Grid needs to be set first by method "set_grid".')
        return self.gd.get_grid_wf_density(a,k,pad)


    def get_grid_density(self,pad=True):
        """
        Return electron density on grid.

        Do not perform accurate analysis on this density.
        Integrated density differs from the total number of electrons.
        Bader analysis inaccurate.

        parameters:
        pad:      padded edges
        """
        if self.flags['grid']==False:
            raise AssertionError('Grid needs to be set first by method "set_grid".')
        return self.gd.get_grid_density(pad)


    def get_grid_LDOS(self,bias=None,window=None,pad=True):
        """
        Return electron density over selected states around the Fermi-level.

        parameters:
        -----------
        bias:      bias voltage (eV) with respect to Fermi-level.
                   Negative means probing occupied states.
        window:    2-tuple for lower and upper bounds wrt. Fermi-level
        pad:       padded edges
        """
        if self.flags['grid']==False:
            raise AssertionError('Grid needs to be set first by method "set_grid".')
        return self.gd.get_grid_LDOS(bias,window,pad)


    #
    # Mulliken population analysis tools
    #
    def _init_mulliken(self):
        """ Initialize Mulliken analysis. """
        if self.calculation_required(self.el.atoms,['energy']):
            raise AssertionError('Electronic structure is not solved yet!')
        if self.flags['Mulliken']==False:
            self.MA = MullikenAnalysis(self)
            self.flags['Mulliken']=True

    def get_dq(self,atoms=None):
        """ Return atoms' excess Mulliken populations.

        The total populations subtracted by
        the numbers of valence electrons.

        """
        self.solve_ground_state(atoms)
        return self.st.get_dq()

    def get_charges(self,atoms=None):
        """ Return atoms' electric charges (Mulliken). """
        return -self.get_dq(atoms)


    def get_atom_mulliken(self,I):
        """
        Return Mulliken population for atom I.

        This is the total population, without the number
        of valence electrons subtracted.

        parameters:
        ===========
        I:        atom index
        """
        self._init_mulliken()
        return self.MA.get_atom_mulliken(I)


    def get_basis_mulliken(self,mu):
        """
        Return Mulliken population of given basis state.

        parameters:
        ===========
        mu:     orbital index (see Elements' methods for indices)
        """
        self._init_mulliken()
        return self.MA.get_basis_mulliken(mu)


    def get_basis_wf_mulliken(self,mu,k,a,wk=True):
        """
        Return Mulliken population for given basis state and wavefunction.

        parameters:
        ===========
        mu:     basis state index
        k:      k-vector index
        a:      eigenstate index
        wk:     include k-point weight in the population?
        """
        self._init_mulliken()
        return self.MA.get_basis_wf_mulliken(mu,k,a,wk)


    def get_atom_wf_mulliken(self,I,k,a,wk=True):
        """
        Return Mulliken population for given atom and wavefunction.

        parameters:
        ===========
        I:      atom index (if None, return an array for all atoms)
        k:      k-vector index
        a:      eigenstate index
        wk:     embed k-point weight in population
        """
        self._init_mulliken()
        return self.MA.get_atom_wf_mulliken(I,k,a,wk)


    def get_atom_wf_all_orbital_mulliken(self,I,k,a):
        """
        Return orbitals' Mulliken populations for given atom and wavefunction.

        parameters:
        ===========
        I:      atom index (returned array size = number of orbitals on I)
        k:      k-vector index
        a:      eigenstate index
        """
        self._init_mulliken()
        return self.MA.get_atom_wf_all_orbital_mulliken(I,k,a)


    def get_atom_wf_all_angmom_mulliken(self,I,k,a,wk=True):
        """
        Return atom's Mulliken populations for all angmom for given wavefunction.

        parameters:
        ===========
        I:        atom index
        k:        k-vector index
        a:        eigenstate index
        wk:       embed k-point weight into population

        return: array (length 3) containing s,p and d-populations
        """
        self._init_mulliken()
        return self.MA.get_atom_wf_all_angmom_mulliken(I,k,a,wk)


    #
    #  Densities of states methods
    #
    def _init_DOS(self):
        """ Initialize Density of states analysis. """
        if self.calculation_required(self.el.atoms,['energy']):
            raise AssertionError('Electronic structure is not solved yet!')
        if self.flags['DOS']==False:
            self.DOS = DensityOfStates(self)
            self.flags['DOS']=True


    def get_local_density_of_states(self,projected=False,width=0.05,window=None,npts=501):
        """
        Return state density for all atoms as a function of energy.

        parameters:
        ===========
        projected: return local density of states projected for
                   angular momenta 0,1 and 2 (s,p and d)
        width:     energy broadening (in eV)
        window:    energy window around Fermi-energy; 2-tuple (eV)
        npts:      number of grid points for energy

        return:    projected==False:
                        energy grid, ldos[atom,grid]
                   projected==True:
                        energy grid,
                        ldos[atom, grid],
                        pldos[atom, angmom, grid]
        """
        self._init_DOS()
        return self.DOS.get_local_density_of_states(projected,width,window,npts)


    def get_density_of_states(self,broaden=False,projected=False,occu=False,width=0.05,window=None,npts=501):
        """
        Return the full density of states.

        Sum of states over k-points. Zero is the Fermi-level.
        Spin-degeneracy is NOT counted.

        parameters:
        ===========
        broaden:     * If True, return broadened DOS in regular grid
                       in given energy window.
                     * If False, return energies of all states, followed
                       by their k-point weights.
        projected:   project DOS for angular momenta
        occu:        for not broadened case, return also state occupations
        width:       Gaussian broadening (eV)
        window:      energy window around Fermi-energy; 2-tuple (eV)
        npts:        number of data points in output

        return:      * if projected: e[:],dos[:],pdos[l,:] (angmom l=0,1,2)
                     * if not projected: e[:],dos[:]
                       * if broaden: e[:] is on regular grid, otherwise e[:] are
                         eigenvalues and dos[...] corresponding weights
                     * if occu: e[:],dos[:],occu[:]

        """
        self._init_DOS()
        return self.DOS.get_density_of_states(broaden,projected,occu,width,window,npts)



    # Bonding analysis
    def _init_bonds(self):
        """ Initialize Mulliken bonding analysis. """
        if self.calculation_required(self.el.atoms,['energy']):
            raise AssertionError('Electronic structure is not solved yet!')
        if self.flags['bonds']==False:
            self.bonds = MullikenBondAnalysis(self)
            self.flags['bonds']=True


    def get_atom_energy(self,I=None):
        """
        Return the energy of atom I (in eV).

        Warning: bonding & atom energy analysis less clear for
        systems where orbitals overlap with own periodic images.

        parameters:
        ===========
        I:         atom index. If None, return all atoms' energies
                   as an array.
        """
        self._init_bonds()
        return self.bonds.get_atom_energy(I)



    def get_mayer_bond_order(self,i,j):
        """
        Return Mayer bond-order between two atoms.

        Warning: bonding & atom energy analysis less clear for
        systems where orbitals overlap with own periodic images.

        parameters:
        ===========
        I:        first atom index
        J:        second atom index
        """
        self._init_bonds()
        return self.bonds.get_mayer_bond_order(i,j)


    def get_promotion_energy(self,I=None):
        """
        Return atom's promotion energy (in eV).

        Defined as:
            E_prom,I = sum_(mu in I) [q_(mu) - q_(mu)^0] epsilon_mu

        parameters:
        ===========
        I:         atom index. If None, return all atoms' energies
                   as an array.
        """
        self._init_bonds()
        return self.bonds.get_promotion_energy(I)


    def get_bond_energy(self,i,j):
        """
        Return the absolute bond energy between atoms (in eV).

        Warning: bonding & atom energy analysis less clear for
        systems where orbitals overlap with own periodic images.

        parameters:
        ===========
        i,j:     atom indices
        """
        self._init_bonds()
        return self.bonds.get_bond_energy(i,j)


    def get_atom_and_bond_energy(self,i=None):
        """
        Return given atom's contribution to cohesion.

        parameters:
        ===========
        i:    atom index. If None, return all atoms' energies
              as an array.
        """
        self._init_bonds()
        return self.bonds.get_atom_and_bond_energy(i)


    def get_covalent_energy(self,mode='default',i=None,j=None,width=None,window=None,npts=501):
        """
        Return covalent bond energies in different modes. (eV)

        ecov is described in
        Bornsen, Meyer, Grotheer, Fahnle, J. Phys.:Condens. Matter 11, L287 (1999) and
        Koskinen, Makinen Comput. Mat. Sci. 47, 237 (2009)



        parameters:
        ===========
        mode:    'default' total covalent energy
                 'orbitals' covalent energy for orbital pairs
                 'atoms' covalent energy for atom pairs
                 'angmom' covalent energy for angular momentum components
        i,j:     atom or orbital indices, or angular momentum pairs
        width:   * energy broadening (in eV) for ecov
                 * if None, return energy eigenvalues and corresponding
                   covalent energies in arrays, directly
        window:  energy window (in eV wrt Fermi-level) for broadened ecov
        npts:    number of points in energy grid (only with broadening)

        return:
        =======
        x,y:     * if width==None, x is list of energy eigenvalues (including k-points)
                   and y covalent energies of those eigenstates
                 * if width!=None, x is energy grid for ecov.
                 * energies (both energy grid and ecov) are in eV.

        Note: energies are always shifted so that Fermi-level is at zero.
              Occupations are not otherwise take into account (while k-point weights are)
        """
        self._init_bonds()
        return self.bonds.get_covalent_energy(mode,i,j,width,window,npts)


    def add_pair_potential(self,i,j,v,eVA=True):
        """
        Add pair interaction potential function for elements or atoms

        parameters:
        ===========
        i,j:    * atom indices, if integers (0,1,2,...)
                * elements, if strings ('C','H',...)
        v:      Pair potential function.
                Only one potential per element and atom pair allowed.
                Syntax:  v(r,der=0), v(r=None) returning the
                interaction range in Bohr or Angstrom.
        eVA:    True for v in eV and Angstrom
                False for v in Hartree and Bohr
        """
        self.pp.add_pair_potential(i,j,v,eVA)





### Helper functions

def database_from_path(path):
    if path is None:
        path = '.'

    if not type(path) == list:
        path = [ path ]

    fns = [ ]
    for p in path:
        fns += glob.glob('%s/*.elm' % p)

    elements = { }
    tables = { }

    if len(fns) > 0:
        for fn in fns:
            i0 = fn.rfind('/')
            i1 = fn.rfind('.')
            el1 = fn[i0+1:i1]

            elements[el1] = fn

            for fn2 in fns:
                i0 = fn.rfind('/')
                i1 = fn.rfind('.')
                el2 = fn2[i0+1:i1]

                for p in path:
                    if os.path.exists('%s/%s_%s.par' % ( p, el1, el2 )):
                        tables['%s%s' % ( el1, el2 )] = \
                                      '%s/%s_%s.par' % ( p, el1, el2 )
                    else:
                        if os.path.exists('%s/%s_%s.par' % ( p, el2, el1 )):
                            tables['%s%s' % ( el1, el2 )] = \
                                          '%s/%s_%s.par' % ( p, el2, el1 )

    else:
        fns = [ ]
        for p in path:
            fns = glob.glob('%s/*-*.skf' % p)

        if len(fns) > 0:
            for fn in fns:
                i0 = fn.rfind('/')
                i1 = fn.rfind('-')
                i2 = fn.rfind('.')
                el1 = fn[i0+1:i1]
                el2 = fn[i1+1:i2]

                if el1 == el2:
                    elements[el1] = fn
                tables['%s%s' % ( el1, el2 )] = fn
        else:
            raise RuntimeError('No Slater-Koster database found in directory '
                               '%s.' % path)

    return { 'elements': elements, 'tables': tables }
