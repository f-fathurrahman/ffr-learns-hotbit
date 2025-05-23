# Copyright (C) 2008 NSC Jyvaskyla
# Please see the accompanying LICENSE file for further information.

from box import mix
import numpy as np
from hotbit import auxil
from box.interpolation import MultipleSplineFunction
from weakref import proxy
from copy import copy
from os import path
from hotbit.io import read_HS
from math import cos, sin, sqrt

from _hotbit import fast_slako_transformations

dot = np.dot
array = np.array
norm = np.linalg.norm
outer = np.outer
zeros = np.zeros

aux={'s':0,'p':1,'d':2}
itypes={'ss':['s'],'sp':['s'],'ps':['s'],'sd':['s'],'ds':['s'],
        'pp':['s','p'],'pd':['s','p'],'dp':['s','p'],'dd':['s','p','d']}
integrals={'dds':0,'ddp':1,'ddd':2,'pds':3,'pdp':4,'pps':5,'ppp':6,\
           'sds':7,'sps':8,'sss':9,'dps':10,'dpp':11,'dss':12,'pss':13}





class Interactions:

    def __init__(self, calc):

        print("\n<div> ENTER Interactions constructor\n")

        from os import environ
        from os.path import isfile

        tables = copy(calc.get('tables'))
        print("tables = ", tables)

        present = calc.el.get_present()
        print("present = ", present)

        default = environ.get('HOTBIT_PARAMETERS')
        # current = path.abspath('.')
        self.nullpar = path.join(default,'null.par')

        files = {}
        for k1 in present:
            for k2 in present:
                files[k1+k2] = None

        # set customized files
        if tables != None:
            for key in tables:
                if key == 'rest':
                    continue
                #
                e1, e2 = auxil.separate_symbols(key)
                file = tables[key]
                if file == None:
                    file = self.nullpar
                #if file[-4:]!='.par':
                #    file+='.par'
                if not isfile(file):
                    raise RuntimeError('Custom parameter file "%s" for %s-%s interaction not found.' %(file,e1,e2))
                file = path.abspath(file)
                files[e1+e2] = file
                if not e2+e1 in tables:
                    files[e2+e1] = file


        # set interactions from default place
        if tables==None or tables!=None and 'rest' in tables and tables['rest']=='default':
            for e1 in present:
                for e2 in present:
                    if files[e1+e2] != None:
                        continue
                    def12 = path.join(default, '%s_%s.par' %(e1,e2))
                    def21 = path.join(default, '%s_%s.par' %(e2,e1))

                    if isfile(def12):
                        file = def12
                    elif isfile(def21):
                        file = def21
                    else:
                        raise RuntimeError('Default parameter files "%s" or "%s" for %s-%s interaction not found.' %(def12,def21,e1,e2))
                    file = path.abspath(file)
                    files[e1+e2] = file
                    files[e2+e1] = file
        #
        print("Files: ")
        for kf,ff in files.items():
            print(f"pair: {kf} file: {ff}")

        self.files = files
        self.calc = proxy(calc)
        self.present = present
        self.max_cut = 0.0 # maximum interaction range in Bohrs
        self.read_tables()
        self.first = True

        print("\n</div>EXIT Interactions constructor\n")



    def __del__(self):
        pass

    def get_files(self):
        """ Return the list of Slater-Koster table files."""
        return self.files

    def get_cutoffs(self):
        """ Return element pair cutoff dictionary. """
        return self.cut

    def greetings(self):
        """ Return the interaction specifications """
        txt='Interactions:\n'
        for i,s1 in enumerate(self.present):
            for s2 in self.present:
                file = self.files[s1+s2]
                if file==None:
                    txt+='  %s%s: None\n' %(s1,s2)
                else:
                    txt+='  %s%s in %s\n' %(s1,s2,file)
                    doc=mix.find_value(file,'slako_comment',fmt='strings',default=['no slako doc'])
                    for line in doc:
                        txt+='    *'+line.lstrip()+'\n'
        return txt


    def read_tables(self):
        """
        Read par-file tables. Files are tabulated with |ket> having >= angular momentum,
        so one has to copy values for interaction the other way round.
        """
        #
        print("\n<div> ENTER Interactions.read_tables\n")
        #
        self.h = {}
        self.s = {}
        self.cut = {}
        self.maxh = {}
#        self.kill_radii={}
        #
        for si in self.present:
            for sj in self.present:
                if self.files[si+sj] is None:
                    raise RuntimeError('No parametrization specified for %s-%s interaction.' % ( si, sj ))
                if self.files[sj+si] is None:
                    raise RuntimeError('No parametrization specified for %s-%s interaction.' % ( sj, si ))
                #
                x_ij, table_ij = read_HS(self.files[si+sj], si, sj)
                self.cut[si+sj] = x_ij[-1]
                print("Reading file: ", self.files[si+sj])
                #
                x_ji, table_ji = read_HS(self.files[sj+si], sj, si)
                self.cut[sj+si] = x_ji[-1]
                print("Reading file: ", self.files[sj+si])
                #
                self.max_cut = max(self.max_cut,self.cut[si+sj])
                self.maxh[si+sj] = max( [max(np.abs(table_ij[:,i])) for i in range(1,11)] )

                ei, ej = self.calc.el.elements[si], self.calc.el.elements[sj]
                valence_i, valence_j = ei.get_valence_orbitals(), ej.get_valence_orbitals()

                pair = si + sj
                self.h[pair] = MultipleSplineFunction(x_ij)
                self.s[pair] = MultipleSplineFunction(x_ji)
                for vi in valence_i:
                    for vj in valence_j:
                        li, lj = vi[1], vj[1]
                        # for given valence orbitals, go through all possible integral types (sigma,pi,...)
                        for itype in itypes[li+lj]:
                            table = '%s(%s)-%s(%s)-%s' %(si,li,sj,lj,itype)
                            short = '%s%s%s' %(li,lj,itype)
                            # self.h['C(p)-H(s)-sigma']=...
                            if short[0:2] == 'ps' or short[0:2] == 'ds' or short[0:2] == 'dp':
                                # this is tabulated in other table; switch order -> parity factor
                                parity = (-1)**( aux[li]+aux[lj] )
                                index = integrals[short[1]+short[0]+short[2]]
                                self.h[pair].add_function(table_ji[:,index]*parity,table,integrals[short])
                                self.s[pair].add_function(table_ji[:,index+10]*parity,table,integrals[short])
                            else:
                                index=integrals[short]
                                self.h[pair].add_function(table_ij[:,index],table,integrals[short])
                                self.s[pair].add_function(table_ij[:,index+10],table,integrals[short])

        # cutoffs for atom pair indices
        N = self.calc.el.N
        self.hscut = np.zeros((N,N),float)
        for i,si in enumerate(self.calc.el.symbols):
            for j,sj in enumerate(self.calc.el.symbols):
                self.hscut[i,j] = self.cut[si+sj]
        self.calc.el.set_cutoffs(self.cut)
        #
        print("\n</div> EXIT Interactions.read_tables\n")
        #

    def get_tables(self, si, sj):
        return self.h[si+sj], self.s[si+sj]


    def plot_table(self,e1,e2,der=0):
        """ Plot SlaKo table for given elements. """
        import pylab as pl
        R=np.linspace(0,self.cut[e1+e2],1000)
        n=self.h[e1+e2].get_number_of_functions()
        nx=max(n/2, 1)
        ny=n/nx
        for i in range(n):
            pl.subplot(nx,ny,i+1)
            pl.plot(R,[self.h[e1+e2](r,der=der)[i] for r in R],label='H')
            pl.plot(R,[self.s[e1+e2](r,der=der)[i] for r in R],label='S')
            pl.xlim(0,R[-1])
            pl.ylim(-self.maxh[e1+e2],self.maxh[e1+e2])
            #pl.title('%s-%s integral %s (der=%i)' %(e1,e2,self.integrals[i],der))
            pl.title('%s-%s (der=%i)' %(e1,e2,der))
            pl.xlabel('r (Bohr)')
            pl.ylabel('H (Hartree) and S')
        pl.show()


    def rotation_transformation(self,n):
        '''
        Return the transpose of 9x9 rotation transformation matrix for given symmetry operation.

        For example, below transformation for a rotation around z-axis:

            | 1   0   0   0    0    0    0   0   0 |   s  } s-orbital
            | 0  ca  sa   0    0    0    0   0   0 |   px
            | 0 -sa  ca   0    0    0    0   0   0 |   py } p-orbitals         ca = cos(a)
            | 0   0   0   1    0    0    0   0   0 |   pz                      sa = sin(a)
         D= | 0   0   0   0   c2a   0    0 -s2a  0 |   dxy                     c2a = cos(2a)
            | 0   0   0   0    0   ca   -sa  0   0 |   dyz                     s2a = sin(2a)
            | 0   0   0   0    0   sa   ca   0   0 |   dzx } d-orbitals
            | 0   0   0   0  s2a    0    0 c2a   0 |   dx2-y2
            | 0   0   0   0    0    0    0   0   1 |   d3z2-r2

        @param n: 3-tuple of symmetry transformation
        '''
        
        print("Calling Interactions.rotation_transformation")

        R = self.calc.el.rotation(n)

        if np.all(abs(R.diagonal()-1)<1E-12):
            #no rotation
            return np.eye(9)

        elif np.any(np.array(self.calc.el.get_property_lists(['no']))>4) and abs(R[2,2]-1)>1E-12:
            # These are more complex transformations for d-orbitals. Use only if needed, for they
            # may take more time.
            if abs(R[2,2]+1)<1E-12:
                # only xy-plane reflection (orbitals antisymm. wrt z-axis change sign)
                D = np.eye(9)
                D[3,3] = -1
                D[5,5] = -1
                D[6,6] = -1
                return D
            else:
                # This is the general transformation with full expressions:
                rot = self.calc.el.rotation(n,angles=True)
                theta,phi,angle = rot
                D = orbital_transformations_from_angles((theta,phi,angle))
                return D.transpose()

        else:
            # rotation around z-axis, bit more simple transformations than the general case
            ca, sa = R[0,0], -R[0,1]
            c2a, s2a = 2*ca**2-1, 2*sa*ca
            D = np.diag((1.0,ca,ca,1.0,c2a,ca,ca,c2a,1.0))

            D[1,2] = sa
            D[2,1] = -sa
            D[5,6] = -sa
            D[6,5] = sa
            D[4,7] = -s2a
            D[7,4] = s2a

            D[1:4,1:4] = R[:,:].transpose()
            return D.transpose()


    def get_phases(self):
        """ Return phases for symmetry operations and k-points.

        phases[n,k]
        """
        return np.array(self.phases)


    def get_matrices(self, kpts=None):
        """ Hamiltonian and overlap matrices. """
        
        print("\n<div> ENTER Interactions.get_matrices")

        timing = False
        el = self.calc.el
        states = self.calc.st

        seps = self.calc.get('sepsilon')
        orbs = el.orbitals()
        norb = len(orbs)

        if kpts is None:
            nk = states.nk
            ks = states.k
        else:
            ks       = np.asarray(kpts)
            ks.shape = (-1, 3)
            nk       = ks.shape[0]

        H0  = np.zeros( (nk,norb,norb), complex )
        S   = np.zeros( (nk,norb,norb), complex )
        dH0 = np.zeros( (nk,norb,norb,3), complex )
        dS  = np.zeros( (nk,norb,norb,3), complex )

        # FFR: fixed size????
        h, s, dh, ds = zeros((14,)), zeros((14,)), zeros((14,3)), zeros((14,3))

        phases = []
        DTn = []
        Rot = []
        for n in range(len(el.ntuples)):
            nt = el.ntuples[n]
            phases.append( np.array([np.exp(1j*np.dot(nt,k))
                                     for k in ks]) )
            DTn.append( self.rotation_transformation(nt) )
            Rot.append( self.calc.el.rotation(nt) )
        self.phases = phases

        lst = el.get_property_lists(['i', 's', 'no', 'o1'])
        Rijn, dijn = self.calc.el.get_distances()
        # defined in elements.py
        # 'i'=index; 's'=symbol; 'no'=number of orbitals; 'o1'= first orbital
        for i,si,noi,o1i in lst:
            a, b = o1i, o1i+noi # block indices?
            # on-site energies only for n==0
            for orb in el.orbitals(i):
                ind = orb['index']
                H0[:,ind,ind] = orb['energy']
                S[:,ind,ind]  = 1.0 + seps
            #
            for j,sj,noj,o1j in lst[i:]:
                c, d = o1j, o1j+noj
                print(f"Loop (i={i},j={j}) (si={si},sj={sj}) (noi={noi},noj={noj}) (o1i={o1i},o1j={o1j})")
                htable = self.h[si+sj]
                stable = self.s[si+sj]
                ij_interact = False
                r1, r2 = htable.get_range()
                #
                for n, (rij,dij) in enumerate(zip(Rijn[:,i,j],dijn[:,i,j])):
                    if i==j and n==0:
                        continue
                    nt = el.ntuples[n]
                    h.fill(0)
                    s.fill(0)
                    dh.fill(0)
                    ds.fill(0)

                    rijh  = rij/dij
                    if dij < 0.1:
                        print(nt)
                        raise AssertionError('Distance between atoms %i and %i is only %.4f Bohr' %(i,j,dij) )
                    #
                    assert dij > 0.1
                    #
                    if not r1 <= dij <= r2:
                        continue
                    ij_interact = True

                    # interpolate Slater-Koster tables and derivatives
                    hij, dhij = htable(dij)
                    sij, dsij = stable(dij)

                    indices = htable.get_indices()
                    h[indices], s[indices] = hij, sij
                    dh[indices], ds[indices] = outer(dhij,rijh), outer(dsij,rijh)

                    # make the Slater-Koster transformations
                    ht, st, dht, dst = \
                        fast_slako_transformations(rijh,dij,noi,noj,h,s,dh,ds)

                    # Here we do the MEL transformation;
                    # H'_ij = sum_k H_ik * D_kj^T  ( |j> = sum_k D_jk |k> )
                    DT = DTn[n]
                    ht = dot( ht,DT[0:noj,0:noj] )
                    st = dot( st,DT[0:noj,0:noj] )
                    dht = dot( dht.transpose((2,0,1)),DT[0:noj,0:noj] ).transpose((1,2,0))
                    dst = dot( dst.transpose((2,0,1)),DT[0:noj,0:noj] ).transpose((1,2,0))
                    phase = phases[n]
                    hblock  = outer(phase,ht.flatten()).reshape(nk,noi,noj)
                    sblock  = outer(phase,st.flatten()).reshape(nk,noi,noj)
                    dhblock = outer(phase,-dht.flatten()).reshape(nk,noi,noj,3)
                    dsblock = outer(phase,-dst.flatten()).reshape(nk,noi,noj,3)

                    H0[ :,a:b,c:d]   += hblock
                    S[  :,a:b,c:d]   += sblock
                    dH0[:,a:b,c:d,:] += dhblock
                    dS[ :,a:b,c:d,:] += dsblock

                    if i!=j and ij_interact:
                        # construct the other derivatives wrt. atom j.
                        Rot = self.calc.el.Rot[n]
                        dht2 = dot( dht,Rot )
                        dst2 = dot( dst,Rot )
                        dh2block = outer(phase,dht2.flatten()).reshape(nk,noi,noj,3)
                        ds2block = outer(phase,dst2.flatten()).reshape(nk,noi,noj,3)
                        dH0[:,c:d,a:b,:] += dh2block.transpose((0,2,1,3)).conjugate()
                        dS[ :,c:d,a:b,:] += ds2block.transpose((0,2,1,3)).conjugate()

                if i != j and ij_interact:
                    # Hermitian (and anti-Hermitian) conjugates;
                    # only if matrix block non-zero
                    # ( H(k) and S(k) can be (anti)symmetrized as a whole )
                    # TODO: symmetrization should be done afterwards
                    H0[ :,c:d,a:b]   =  H0[ :,a:b,c:d].transpose((0,2,1)).conjugate()
                    S[  :,c:d,a:b]   =  S[  :,a:b,c:d].transpose((0,2,1)).conjugate()

        if kpts is None:
            self.H0, self.S, self.dH0, self.dS = H0, S, dH0, dS

        # print info about Hamiltonian
        if self.first:
            nonzero = sum( abs(S[0].flatten())>1E-15 )
            self.calc.out('Hamiltonian ~%.3f %% filled.' %(nonzero*100.0/norb**2) )
            self.first = False

        print("\n</div> EXIT Interactions.get_matrices\n")
        return H0, S, dH0, dS


    def get_cutoff(self):
        """ Maximum cutoff. """
        return self.max_cut




def orbital_transformations_from_angles(angles):
    """
    Return the unitary transformation matrix, given rotation angles

    theta, phi, angle = angles, where theta, phi is the direction of rotation axis, and angle
    is the rotated angle

    Transformation equations are direct output from SAGE calculation (20.5 2013)
    """
    theta, phi, angle = angles
    D = np.zeros((9,9))
    # s-orbital
    D[0,0] = 1.0
    pi = np.pi
    ca, sa, ct, st, cp, sp = cos(angle), sin(angle), cos(theta), sin(theta), cos(phi), sin(phi)
    ca2, sa2, ct2, st2, cp2, sp2 = ca**2, sa**2, ct**2, st**2, cp**2, sp**2
    ca1 = 1-ca
    # p-orbitals (transforms the same as rotation matrix itself)
    D[1,1] = (ca1*st2*cp2 + ca)
    D[1,2] = (ca1*sp*st2*cp + sa*ct)
    D[1,3] = (ca1*cp*ct - sa*sp)*st
    D[2,1] = (ca1*sp*st2*cp - sa*ct)
    D[2,2] = (ca1*sp2*st2 + ca)
    D[2,3] = (ca1*sp*ct + sa*cp)*st
    D[3,1] = (ca1*cp*ct + sa*sp)*st
    D[3,2] = (ca1*sp*ct - sa*cp)*st
    D[3,3] = (ca1*ct2 + ca)

    # d-orbitals
    D[4,4] = (2*(pi + pi*ca2 - 2*pi*ca)*sp2*st**4*cp2 - pi*sa2*ct2 - ((pi*ca2 - pi*ca)*sp2 + (pi*ca2 - pi*ca)*cp2)*st2 + pi*ca2)/pi
    D[4,5] = -(((pi - pi*ca)*sa*sp**3 - (pi - pi*ca)*sa*sp*cp2 - 2*(pi + pi*ca2 - 2*pi*ca)*sp2*cp*ct)*st**3 - ((pi - pi*ca)*sa*sp*ct2 - pi*sa*sp*ca + (pi*sa2 - pi*ca2 + pi*ca)*cp*ct)*st)/pi
    D[4,6] = -(((pi - pi*ca)*sa*sp2*cp - (pi - pi*ca)*sa*cp**3 - 2*(pi + pi*ca2 - 2*pi*ca)*sp*cp2*ct)*st**3 + ((pi - pi*ca)*sa*cp*ct2 - pi*sa*ca*cp - (pi*sa2 - pi*ca2 + pi*ca)*sp*ct)*st)/pi
    D[4,7] = -(((pi + pi*ca2 - 2*pi*ca)*sp**3*cp - (pi + pi*ca2 - 2*pi*ca)*sp*cp**3)*st**4 + ((pi - pi*ca)*sa*sp2 + (pi - pi*ca)*sa*cp2)*st2*ct + 2*pi*sa*ca*ct)/pi
    D[4,8] = -1./3*(((pi + pi*ca2 - 2*pi*ca)*sp**3*cp + (pi + pi*ca2 - 2*pi*ca)*sp*cp**3)*st**4 - (2*(pi + pi*ca2 - 2*pi*ca)*sp*cp*ct2 - 2*(pi*sa2 - pi*ca2 + pi*ca)*sp*cp - 3*((pi - pi*ca)*sa*sp2 - (pi - pi*ca)*sa*cp2)*ct)*st2)*sqrt(3)/pi
    D[5,4] = (((pi - pi*ca)*sa*sp**3 - (pi - pi*ca)*sa*sp*cp2 + 2*(pi + pi*ca2 - 2*pi*ca)*sp2*cp*ct)*st**3 - ((pi - pi*ca)*sa*sp*ct2 - pi*sa*sp*ca - (pi*sa2 - pi*ca2 + pi*ca)*cp*ct)*st)/pi
    D[5,5] = -((pi*ca2 - pi*ca)*ct2 - (2*(pi + pi*ca2 - 2*pi*ca)*sp2*ct2 - pi*sa2*cp2 - (pi*ca2 - pi*ca)*sp2)*st2 - pi*ca2)/pi
    D[5,6] = -((pi - pi*ca)*sa*ct**3 + pi*sa*ca*ct - (2*(pi + pi*ca2 - 2*pi*ca)*sp*cp*ct2 + (pi*sa2 - pi*ca2 + pi*ca)*sp*cp + ((pi - pi*ca)*sa*sp2 + (pi - pi*ca)*sa*cp2)*ct)*st2)/pi
    D[5,7] = ((2*(pi - pi*ca)*sa*sp2*cp - ((pi + pi*ca2 - 2*pi*ca)*sp**3 - (pi + pi*ca2 - 2*pi*ca)*sp*cp2)*ct)*st**3 - ((pi - pi*ca)*sa*cp*ct2 - pi*sa*ca*cp + (pi*sa2 - pi*ca2 + pi*ca)*sp*ct)*st)/pi
    D[5,8] = -1./3*(((pi + pi*ca2 - 2*pi*ca)*sp**3 + (pi + pi*ca2 - 2*pi*ca)*sp*cp2)*st**3*ct - (3*(pi - pi*ca)*sa*cp*ct2 + 2*(pi + pi*ca2 - 2*pi*ca)*sp*ct**3 + 3*pi*sa*ca*cp + (pi*sa2 - pi*ca2 + pi*ca)*sp*ct)*st)*sqrt(3)/pi
    D[6,4] = (((pi - pi*ca)*sa*sp2*cp - (pi - pi*ca)*sa*cp**3 + 2*(pi + pi*ca2 - 2*pi*ca)*sp*cp2*ct)*st**3 + ((pi - pi*ca)*sa*cp*ct2 - pi*sa*ca*cp + (pi*sa2 - pi*ca2 + pi*ca)*sp*ct)*st)/pi
    D[6,5] = ((pi - pi*ca)*sa*ct**3 + pi*sa*ca*ct + (2*(pi + pi*ca2 - 2*pi*ca)*sp*cp*ct2 + (pi*sa2 - pi*ca2 + pi*ca)*sp*cp - ((pi - pi*ca)*sa*sp2 + (pi - pi*ca)*sa*cp2)*ct)*st2)/pi
    D[6,6] = -((pi*ca2 - pi*ca)*ct2 - (2*(pi + pi*ca2 - 2*pi*ca)*cp2*ct2 - pi*sa2*sp2 - (pi*ca2 - pi*ca)*cp2)*st2 - pi*ca2)/pi
    D[6,7] = ((2*(pi - pi*ca)*sa*sp*cp2 - ((pi + pi*ca2 - 2*pi*ca)*sp2*cp - (pi + pi*ca2 - 2*pi*ca)*cp**3)*ct)*st**3 - ((pi - pi*ca)*sa*sp*ct2 - pi*sa*sp*ca - (pi*sa2 - pi*ca2 + pi*ca)*cp*ct)*st)/pi
    D[6,8] = -1./3*(((pi + pi*ca2 - 2*pi*ca)*sp2*cp + (pi + pi*ca2 - 2*pi*ca)*cp**3)*st**3*ct + (3*(pi - pi*ca)*sa*sp*ct2 - 2*(pi + pi*ca2 - 2*pi*ca)*cp*ct**3 + 3*pi*sa*sp*ca - (pi*sa2 - pi*ca2 + pi*ca)*cp*ct)*st)*sqrt(3)/pi
    D[7,4] = -(((pi + pi*ca2 - 2*pi*ca)*sp**3*cp - (pi + pi*ca2 - 2*pi*ca)*sp*cp**3)*st**4 - ((pi - pi*ca)*sa*sp2 + (pi - pi*ca)*sa*cp2)*st2*ct - 2*pi*sa*ca*ct)/pi
    D[7,5] = -((2*(pi - pi*ca)*sa*sp2*cp + ((pi + pi*ca2 - 2*pi*ca)*sp**3 - (pi + pi*ca2 - 2*pi*ca)*sp*cp2)*ct)*st**3 - ((pi - pi*ca)*sa*cp*ct2 - pi*sa*ca*cp - (pi*sa2 - pi*ca2 + pi*ca)*sp*ct)*st)/pi
    D[7,6] = -((2*(pi - pi*ca)*sa*sp*cp2 + ((pi + pi*ca2 - 2*pi*ca)*sp2*cp - (pi + pi*ca2 - 2*pi*ca)*cp**3)*ct)*st**3 - ((pi - pi*ca)*sa*sp*ct2 - pi*sa*sp*ca + (pi*sa2 - pi*ca2 + pi*ca)*cp*ct)*st)/pi
    D[7,7] = 1./2*(((pi + pi*ca2 - 2*pi*ca)*sp**4 - 2*(pi + pi*ca2 - 2*pi*ca)*sp2*cp2 + (pi + pi*ca2 - 2*pi*ca)*cp**4)*st**4 - 2*pi*sa2*ct2 - 2*((pi*ca2 - pi*ca)*sp2 + (pi*ca2 - pi*ca)*cp2)*st2 + 2*pi*ca2)/pi
    D[7,8] = 1./6*(((pi + pi*ca2 - 2*pi*ca)*sp**4 - (pi + pi*ca2 - 2*pi*ca)*cp**4)*st**4 - 2*(6*(pi - pi*ca)*sa*sp*cp*ct + ((pi + pi*ca2 - 2*pi*ca)*sp2 - (pi + pi*ca2 - 2*pi*ca)*cp2)*ct2 - (pi*sa2 - pi*ca2 + pi*ca)*sp2 + (pi*sa2 - pi*ca2 + pi*ca)*cp2)*st2)*sqrt(3)/pi
    D[8,4] = -1./3*(((pi + pi*ca2 - 2*pi*ca)*sp**3*cp + (pi + pi*ca2 - 2*pi*ca)*sp*cp**3)*st**4 - (2*(pi + pi*ca2 - 2*pi*ca)*sp*cp*ct2 - 2*(pi*sa2 - pi*ca2 + pi*ca)*sp*cp + 3*((pi - pi*ca)*sa*sp2 - (pi - pi*ca)*sa*cp2)*ct)*st2)*sqrt(3)/pi
    D[8,5] = -1./3*(((pi + pi*ca2 - 2*pi*ca)*sp**3 + (pi + pi*ca2 - 2*pi*ca)*sp*cp2)*st**3*ct + (3*(pi - pi*ca)*sa*cp*ct2 - 2*(pi + pi*ca2 - 2*pi*ca)*sp*ct**3 + 3*pi*sa*ca*cp - (pi*sa2 - pi*ca2 + pi*ca)*sp*ct)*st)*sqrt(3)/pi
    D[8,6] = -1./3*(((pi + pi*ca2 - 2*pi*ca)*sp2*cp + (pi + pi*ca2 - 2*pi*ca)*cp**3)*st**3*ct - (3*(pi - pi*ca)*sa*sp*ct2 + 2*(pi + pi*ca2 - 2*pi*ca)*cp*ct**3 + 3*pi*sa*sp*ca + (pi*sa2 - pi*ca2 + pi*ca)*cp*ct)*st)*sqrt(3)/pi
    D[8,7] = 1./6*(((pi + pi*ca2 - 2*pi*ca)*sp**4 - (pi + pi*ca2 - 2*pi*ca)*cp**4)*st**4 + 2*(6*(pi - pi*ca)*sa*sp*cp*ct - ((pi + pi*ca2 - 2*pi*ca)*sp2 - (pi + pi*ca2 - 2*pi*ca)*cp2)*ct2 + (pi*sa2 - pi*ca2 + pi*ca)*sp2 - (pi*sa2 - pi*ca2 + pi*ca)*cp2)*st2)*sqrt(3)/pi
    D[8,8] = 1./6*(((pi + pi*ca2 - 2*pi*ca)*sp**4 + 2*(pi + pi*ca2 - 2*pi*ca)*sp2*cp2 + (pi + pi*ca2 - 2*pi*ca)*cp**4)*st**4 + 4*(pi + pi*ca2 - 2*pi*ca)*ct**4 + 2*(pi*sa2 - 4*pi*ca2 + 4*pi*ca)*ct2 - 2*(2*((pi + pi*ca2 - 2*pi*ca)*sp2 + (pi + pi*ca2 - 2*pi*ca)*cp2)*ct2 + (2*pi*sa2 + pi*ca2 - pi*ca)*sp2 + (2*pi*sa2 + pi*ca2 - pi*ca)*cp2)*st2 + 6*pi*ca2)/pi
    return D

def simple_table_notation(table):
    a,b,i=table[2:].split('-')
    return a[-2]+b[-2]+i[0]



s3=np.sqrt(3.0)

def slako_transformations(rhat,dist,noi,noj,h,s,dh,ds):
    """
    Apply Slater-Koster transformation rules to orbitals iorbs and orbitals jorbs,
    where rhat is vector i->j and table gives the values for given tabulated
    matrix elements. Convention: orbital name starts with s,p,d,...
    """
    l,m,n=rhat
    ll,mm,nn=rhat**2
    dl=(np.array([1,0,0])-l*rhat)/dist
    dm=(np.array([0,1,0])-m*rhat)/dist
    dn=(np.array([0,0,1])-n*rhat)/dist
    dll, dmm, dnn = 2*l*dl, 2*m*dm, 2*n*dn

    mat=np.zeros((9,9,14))
    ind=np.zeros((9,9,14),dtype=int)
    der=np.zeros((9,9,14,3))
    cnt=np.zeros((9,9),dtype=int)+1
    cnt[1:,1:]=2
    cnt[4:,4:]=3

    ht=np.zeros((noi,noj))
    st=np.zeros((noi,noj))
    dht=np.zeros((noi,noj,3))
    dst=np.zeros((noi,noj,3))
    mxorb=max(noi,noj)

    mat[0,0,0]=1  #ss
    der[0,0,0]=0
    ind[0,0,0]=9

    if mxorb>=2:   #sp
        mat[0,1,0]=l
        der[0,1,0,:]=dl
        ind[0,1,0]=8

        mat[0,2,0]=m
        der[0,2,0,:]=dm
        ind[0,2,0]=8

        mat[0,3,0]=n
        der[0,3,0,:]=dn
        ind[0,3,0]=8

    if mxorb>=5:   #sd
        mat[0,4,0]=s3*l*m
        der[0,4,0,:]=s3*(dl*m+l*dm)
        ind[0,4,0]=7

        mat[0,5,0]=s3*m*n
        der[0,5,0,:]=s3*(dm*n+m*dn)
        ind[0,5,0]=7

        mat[0,6,0]=s3*n*l
        der[0,6,0,:]=s3*(dn*l+n*dl)
        ind[0,6,0]=7

        mat[0,7,0]=0.5*s3*(ll-mm)
        der[0,7,0,:]=0.5*s3*(dll-dmm)
        ind[0,7,0]=7

        mat[0,8,0]=nn-0.5*(ll+mm)
        der[0,8,0,:]=dnn-0.5*(dll+dmm)
        ind[0,8,0]=7

    if mxorb>=2: #pp
        mat[1,1,0:2]=[ll, 1-ll]
        der[1,1,0:2,:]=[dll, -dll]
        ind[1,1,0:2]=[5,6]

        mat[1,2,0:2]=[l*m, -l*m]
        der[1,2,0:2,:]=[dl*m+l*dm, -(dl*m+l*dm)]
        ind[1,2,0:2]=[5,6]

        mat[1,3,0:2]=[l*n, -l*n]
        der[1,3,0:2,:]=[dl*n+l*dn, -(dl*n+l*dn)]
        ind[1,3,0:2]=[5,6]

    if mxorb>=5: #pd
        mat[1,4,0:2]=[s3*ll*m, m*(1-2*ll)]
        der[1,4,0:2,:]=[s3*(dll*m+ll*dm), dm*(1-2*ll)+m*(-2*dll)]
        ind[1,4,0:2]=[3,4]

        mat[1,5,0:2]=[s3*l*m*n, -2*l*m*n]
        der[1,5,0:2,:]=[s3*(dl*m*n+l*dm*n+l*m*dn), -2*(dl*m*n+l*dm*n+l*m*dn)]
        ind[1,5,0:2]=[3,4]

        mat[1,6,0:2]=[s3*ll*n, n*(1-2*ll)]
        der[1,6,0:2,:]=[s3*(dll*n+ll*dn), dn*(1-2*ll)+n*(-2*dll)]
        ind[1,6,0:2]=[3,4]

        mat[1,7,0:2]=[0.5*s3*l*(ll-mm), l*(1-ll+mm)]
        der[1,7,0:2,:]=[0.5*s3*(dl*(ll-mm)+l*(dll-dmm)), dl*(1-ll+mm)+l*(-dll+dmm)]
        ind[1,7,0:2]=[3,4]

        mat[1,8,0:2]=[l*(nn-0.5*(ll+mm)), -s3*l*nn]
        der[1,8,0:2,:]=[dl*(nn-0.5*(ll+mm))+l*(dnn-0.5*(dll+dmm)), -s3*(dl*nn+l*dnn)]
        ind[1,8,0:2]=[3,4]

    if mxorb>=2:
        mat[2,2,0:2]=[mm, 1-mm]
        der[2,2,0:2,:]=[dmm, -dmm]
        ind[2,2,0:2]=[5,6]

        mat[2,3,0:2]=[m*n, -m*n]
        der[2,3,0:2,:]=[dm*n+m*dn, -(dm*n+m*dn)]
        ind[2,3,0:2]=[5,6]

    if mxorb>=5:
        mat[2,4,0:2]=[s3*mm*l, l*(1-2*mm)]
        der[2,4,0:2,:]=[s3*(dmm*l+mm*dl), dl*(1-2*mm)+l*(-2*dmm)]
        ind[2,4,0:2]=[3,4]

        mat[2,5,0:2]=[s3*mm*n, n*(1-2*mm)]
        der[2,5,0:2,:]=[s3*(dmm*n+mm*dn), dn*(1-2*mm)+n*(-2*dmm)]
        ind[2,5,0:2]=[3,4]

        mat[2,6,0:2]=[s3*m*n*l, -2*m*n*l]
        der[2,6,0:2,:]=[s3*(dm*n*l+m*dn*l+m*n*dl), -2*(dm*n*l+m*dn*l+m*n*dl)]
        ind[2,6,0:2]=[3,4]

        mat[2,7,0:2]=[0.5*s3*m*(ll-mm), -m*(1+ll-mm)]
        der[2,7,0:2,:]=[0.5*s3*(dm*(ll-mm)+m*(dll-dmm)), -(dm*(1+ll-mm)+m*(dll-dmm))]
        ind[2,7,0:2]=[3,4]

        mat[2,8,0:2]=[m*(nn-0.5*(ll+mm)), -s3*m*nn]
        der[2,8,0:2,:]=[dm*(nn-0.5*(ll+mm))+m*(dnn-0.5*(dll+dmm)), -s3*(dm*nn+m*dnn)]
        ind[2,8,0:2]=[3,4]

    if mxorb>=2:
        mat[3,3,0:2]=[nn, 1-nn]
        der[3,3,0:2,:]=[dnn, -dnn]
        ind[3,3,0:2]=[5,6]

    if mxorb>=5:
        mat[3,4,0:2]=[s3*l*m*n, -2*m*n*l]
        der[3,4,0:2,:]=[s3*(dl*m*n+l*dm*n+l*m*dn), -2*(dm*n*l+m*dn*l+m*n*dl)]
        ind[3,4,0:2]=[3,4]

        mat[3,5,0:2]=[s3*nn*m, m*(1-2*nn)]
        der[3,5,0:2,:]=[s3*(dnn*m+nn*dm), dm*(1-2*nn)+m*(-2*dnn)]
        ind[3,5,0:2]=[3,4]

        mat[3,6,0:2]=[s3*nn*l, l*(1-2*nn)]
        der[3,6,0:2,:]=[s3*(dnn*l+nn*dl), dl*(1-2*nn)+l*(-2*dnn)]
        ind[3,6,0:2]=[3,4]

        mat[3,7,0:2]=[0.5*s3*n*(ll-mm), -n*(ll-mm)]
        der[3,7,0:2,:]=[0.5*s3*(dn*(ll-mm)+n*(dll-dmm)), -(dn*(ll-mm)+n*(dll-dmm))]
        ind[3,7,0:2]=[3,4]

        mat[3,8,0:2]=[n*(nn-0.5*(ll+mm)), s3*n*(ll+mm)]
        der[3,8,0:2,:]=[dn*(nn-0.5*(ll+mm))+n*(dnn-0.5*(dll+dmm)), s3*(dn*(ll+mm)+n*(dll+dmm))]
        ind[3,8,0:2]=[3,4]

    if mxorb>=5:
        mat[4,4,0:3]=[3*ll*mm, ll+mm-4*ll*mm, nn+ll*mm]
        der[4,4,0:3,:]=[3*(dll*mm+ll*dmm), dll+dmm-4*(dll*mm+ll*dmm), dnn+(dll*mm+ll*dmm)]
        ind[4,4,0:3]=[0,1,2]

        mat[4,5,0:3]= [3*l*mm*n, l*n*(1-4*mm), l*n*(mm-1)]
        der[4,5,0:3,:]= [3*(dl*mm*n+l*dmm*n+l*mm*dn), dl*n*(1-4*mm)+l*dn*(1-4*mm)+l*n*(-4*dmm), dl*n*(mm-1)+l*dn*(mm-1)+l*n*(dmm)]
        ind[4,5,0:3]=[0,1,2]

        mat[4,6,0:3]=[3*ll*m*n, m*n*(1-4*ll), m*n*(ll-1)]
        der[4,6,0:3,:]=[3*(dll*m*n+ll*dm*n+ll*m*dn), dm*n*(1-4*ll)+m*dn*(1-4*ll)+m*n*(-4*dll), dm*n*(ll-1)+m*dn*(ll-1)+m*n*(dll)]
        ind[4,6,0:3]=[0,1,2]

        mat[4,7,0:3]=[1.5*l*m*(ll-mm), 2*l*m*(mm-ll), 0.5*l*m*(ll-mm)]
        der[4,7,0:3,:]=[1.5*(dl*m*(ll-mm)+l*dm*(ll-mm)+l*m*(dll-dmm)),\
                    2*(dl*m*(mm-ll)+l*dm*(mm-ll)+l*m*(dmm-dll)),\
                    0.5*(dl*m*(ll-mm)+l*dm*(ll-mm)+l*m*(dll-dmm))]
        ind[4,7,0:3]=[0,1,2]

        mat[4,8,0:3]=[s3*l*m*(nn-0.5*(ll+mm)), - 2*s3*l*m*nn, 0.5*s3*l*m*(1+nn)]
        der[4,8,0:3,:]=[s3*( dl*m*(nn-0.5*(ll+mm))+l*dm*(nn-0.5*(ll+mm))+l*m*(dnn-0.5*(dll+dmm)) ),\
                    -2*s3*(dl*m*nn+l*dm*nn+l*m*dnn),\
                    0.5*s3*( dl*m*(1+nn)+l*dm*(1+nn)+l*m*(dnn) )]
        ind[4,8,0:3]=[0,1,2]

        mat[5,5,0:3]=[3*mm*nn,  (mm+nn-4*mm*nn), (ll+mm*nn)]
        der[5,5,0:3,:]=[3*(dmm*nn+mm*dnn), (dmm+dnn-4*(dmm*nn+mm*dnn)),  (dll+dmm*nn+mm*dnn)]
        ind[5,5,0:3]=[0,1,2]

        mat[5,6,0:3]=[3*m*nn*l, m*l*(1-4*nn), m*l*(nn-1)]
        der[5,6,0:3,:]=[3*(dm*nn*l+m*dnn*l+m*nn*dl),\
                    dm*l*(1-4*nn)+m*dl*(1-4*nn)+m*l*(-4*dnn),\
                    dm*l*(nn-1)+m*dl*(nn-1)+m*l*(dnn)]
        ind[5,6,0:3]=[0,1,2]

        mat[5,7,0:3]=[1.5*m*n*(ll-mm), - m*n*(1+2*(ll-mm)), m*n*(1+0.5*(ll-mm))]
        der[5,7,0:3,:]=[1.5*( dm*n*(ll-mm)+m*dn*(ll-mm)+m*n*(dll-dmm) ),\
                    - ( dm*n*(1+2*(ll-mm))+m*dn*(1+2*(ll-mm))+m*n*(2*dll-2*dmm) ),\
                    dm*n*(1+0.5*(ll-mm))+m*dn*(1+0.5*(ll-mm))+m*n*(0.5*(dll-dmm))]
        ind[5,7,0:3]=[0,1,2]

        mat[5,8,0:3]=[s3*m*n*(nn-0.5*(ll+mm)), s3*m*n*(ll+mm-nn), -0.5*s3*m*n*(ll+mm)]
        der[5,8,0:3,:]=[s3*( dm*n*(nn-0.5*(ll+mm)) + m*dn*(nn-0.5*(ll+mm))+m*n*(dnn-0.5*(dll+dmm)) ),\
                    s3*( dm*n*(ll+mm-nn)+m*dn*(ll+mm-nn)+m*n*(dll+dmm-dnn) ),\
                    - 0.5*s3*( dm*n*(ll+mm)+m*dn*(ll+mm)+m*n*(dll+dmm) )]
        ind[5,8,0:3]=[0,1,2]

        mat[6,6,0:3]=[3*nn*ll, (nn+ll-4*nn*ll), (mm+nn*ll)]
        der[6,6,0:3,:]=[3*(dnn*ll+nn*dll), dnn+dll-4*(dnn*ll+nn*dll), (dmm+dnn*ll+nn*dll)]
        ind[6,6,0:3]=[0,1,2]

        mat[6,7,0:3]=[1.5*n*l*(ll-mm), n*l*(1-2*(ll-mm)), - n*l*(1-0.5*(ll-mm))]
        der[6,7,0:3,:]=[1.5*( dn*l*(ll-mm)+n*dl*(ll-mm)+n*l*(dll-dmm) ),\
                    dn*l*(1-2*(ll-mm))+n*dl*(1-2*(ll-mm))+n*l*(-2*(dll-dmm)),\
                    -( dn*l*(1-0.5*(ll-mm))+n*dl*(1-0.5*(ll-mm))+n*l*(-0.5*(dll-dmm)) )]
        ind[6,7,0:3]=[0,1,2]

        mat[6,8,0:3]=[s3*l*n*(nn-0.5*(ll+mm)), s3*l*n*(ll+mm-nn), - 0.5*s3*l*n*(ll+mm)]
        der[6,8,0:3,:]=[s3*( dl*n*(nn-0.5*(ll+mm))+l*dn*(nn-0.5*(ll+mm))+l*n*(dnn-0.5*(dll+dmm)) ),\
                    s3*( dl*n*(ll+mm-nn)+l*dn*(ll+mm-nn)+l*n*(dll+dmm-dnn) ),\
                    - 0.5*s3*( dl*n*(ll+mm)+l*dn*(ll+mm)+l*n*(dll+dmm) )]
        ind[6,8,0:3]=[0,1,2]

        mat[7,7,0:3]=[0.75*(ll-mm)**2, (ll+mm-(ll-mm)**2), (nn+0.25*(ll-mm)**2)]
        der[7,7,0:3,:]=[0.75*2*(ll-mm)*(dll-dmm), (dll+dmm-2*(ll-mm)*(dll-dmm)), (dnn+0.25*2*(ll-mm)*(dll-dmm))]
        ind[7,7,0:3]=[0,1,2]

        mat[7,8,0:3]=[0.5*s3*(ll-mm)*(nn-0.5*(ll+mm)), s3*nn*(mm-ll), 0.25*s3*(1+nn)*(ll-mm)]
        der[7,8,0:3,:]=[0.5*s3*( (dll-dmm)*(nn-0.5*(ll+mm))+(ll-mm)*(dnn-0.5*(dll+dmm)) ),\
                    s3*( dnn*(mm-ll)+nn*(dmm-dll) ),\
                    0.25*s3*( dnn*(ll-mm)+(1+nn)*(dll-dmm) )]
        ind[7,8,0:3]=[0,1,2]

        mat[8,8,0:3]=[(nn-0.5*(ll+mm))**2, 3*nn*(ll+mm), 0.75*(ll+mm)**2]
        der[8,8,0:3,:]=[2*(nn-0.5*(ll+mm))*(dnn-0.5*(dll+dmm)),\
                    3*( dnn*(ll+mm)+nn*(dll+dmm) ),\
                    0.75*2*(ll+mm)*(dll+dmm)]
        ind[8,8,0:3]=[0,1,2]

    # use the same rules for orbitals when they are reversed (pd ->dp)...
    for a in range(9):
        for b in range(a+1,9):
            mat[b,a,:]=mat[a,b,:]
            der[b,a,:,:]=der[a,b,:,:]
            ind[b,a,:]=ind[a,b,:]

    # ...but use different indices from table
    #pd 3:5-->10:12
    #sd 7->12
    #sp 8->13
    ind[1,0,0]=13
    ind[2,0,0]=13
    ind[3,0,0]=13
    ind[4,0,0]=12
    ind[5,0,0]=12
    ind[6,0,0]=12
    ind[7,0,0]=12
    ind[8,0,0]=12
    ind[4,1,0:2]=[10,11]
    ind[5,1,0:2]=[10,11]
    ind[6,1,0:2]=[10,11]
    ind[7,1,0:2]=[10,11]
    ind[8,1,0:2]=[10,11]
    ind[4,2,0:2]=[10,11]
    ind[5,2,0:2]=[10,11]
    ind[6,2,0:2]=[10,11]
    ind[7,2,0:2]=[10,11]
    ind[8,2,0:2]=[10,11]
    ind[4,3,0:2]=[10,11]
    ind[5,3,0:2]=[10,11]
    ind[6,3,0:2]=[10,11]
    ind[7,3,0:2]=[10,11]
    ind[8,3,0:2]=[10,11]

    for i in range(noi):
        for j in range(noj):
            ht[i,j]=sum( [mat[i,j,k]*h[ind[i,j,k]] for k in range(cnt[i,j])] )
            st[i,j]=sum( [mat[i,j,k]*s[ind[i,j,k]] for k in range(cnt[i,j])] )
            for a in range(3):
                dht[i,j,a]=sum( [mat[i,j,k]*dh[ind[i,j,k],a]+der[i,j,k,a]*h[ind[i,j,k]] for k in range(cnt[i,j])] )
                dst[i,j,a]=sum( [mat[i,j,k]*ds[ind[i,j,k],a]+der[i,j,k,a]*s[ind[i,j,k]] for k in range(cnt[i,j])] )
    return ht, st, dht, dst
