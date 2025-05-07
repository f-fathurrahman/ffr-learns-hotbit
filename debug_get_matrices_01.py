import numpy as np
import matplotlib.pyplot as plt

from _hotbit import fast_slako_transformations

def debug_get_matrices(my_ia, kpts=None):
    
    print("\n<div> ENTER debug_get_matrices")

    el = my_ia.calc.el
    states = my_ia.calc.st

    seps = my_ia.calc.get('sepsilon')
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
    print(f"norb = {norb}")

    # ffr: fixed size???? 
    # ffr: Allocate memory for submatrices
    h = np.zeros((14,))
    s = np.zeros((14,))
    dh = np.zeros((14,3))
    ds = np.zeros((14,3))

    phases = []
    DTn = []
    Rot = []
    for n in range(len(el.ntuples)):
        nt = el.ntuples[n]
        print(f"n={n} nt={nt}")
        #
        phases.append( np.array([np.exp(1j*np.dot(nt,k)) for k in ks]) )
        DTn.append( my_ia.rotation_transformation(nt) )
        Rot.append( my_ia.calc.el.rotation(nt) )
    my_ia.phases = phases

    print("phases = ", phases)
    print("DTn = ", DTn)
    print("Rot = ", Rot)

    lst = el.get_property_lists(['i', 's', 'no', 'o1'])
    Rijn, dijn = my_ia.calc.el.get_distances()
    # defined in elements.py
    # 'i'=index; 's'=symbol; 'no'=number of orbitals; 'o1'= first orbital
    icalc_matrix = 0
    for i,si,noi,o1i in lst:
        a = o1i # start index
        b = o1i + noi # stop index
        #
        # on-site energies only for n==0
        for orb in el.orbitals(i):
            ind = orb['index']
            H0[:,ind,ind] = orb['energy']
            S[:,ind,ind]  = 1.0 + seps
        #
        for j, sj, noj, o1j in lst[i:]:
            c = o1j
            d = o1j + noj
            #
            print(f"\nLoop (i={i},j={j}) (si={si},sj={sj}) (noi={noi},noj={noj}) (o1i={o1i},o1j={o1j})")
            print(f"indices: {a}-{b} {c}-{d}")
            #
            htable = my_ia.h[si+sj]
            stable = my_ia.s[si+sj]
            ij_interact = False
            r1, r2 = htable.get_range()
            #
            for n, (rij,dij) in enumerate(zip(Rijn[:,i,j], dijn[:,i,j])):
                #
                print(f"n = {n} rij={rij} dij={dij}")
                #
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
                print("htable indices = ", indices)
                #
                h[indices], s[indices] = hij, sij
                dh[indices], ds[indices] = np.outer(dhij,rijh), np.outer(dsij,rijh)

                # make the Slater-Koster transformations
                ht, st, dht, dst = fast_slako_transformations(rijh, dij, noi, noj, h, s, dh, ds)
                print(f"rijh={rijh}")
                print(f"h.shape={h.shape}")
                print(f"dh.shape={dh.shape}")

                # Here we do the MEL transformation;
                # H'_ij = sum_k H_ik * D_kj^T  ( |j> = sum_k D_jk |k> )
                DT = DTn[n]
                ht = np.dot( ht,DT[0:noj,0:noj] )
                st = np.dot( st,DT[0:noj,0:noj] )
                dht = np.dot( dht.transpose((2,0,1)),DT[0:noj,0:noj] ).transpose((1,2,0))
                dst = np.dot( dst.transpose((2,0,1)),DT[0:noj,0:noj] ).transpose((1,2,0))
                phase = phases[n]
                hblock  = np.outer(phase,ht.flatten()).reshape(nk,noi,noj)
                sblock  = np.outer(phase,st.flatten()).reshape(nk,noi,noj)
                dhblock = np.outer(phase,-dht.flatten()).reshape(nk,noi,noj,3)
                dsblock = np.outer(phase,-dst.flatten()).reshape(nk,noi,noj,3)

                H0[ :,a:b,c:d]   += hblock
                S[  :,a:b,c:d]   += sblock
                dH0[:,a:b,c:d,:] += dhblock
                dS[ :,a:b,c:d,:] += dsblock

                icalc_matrix += 1

                if i!=j and ij_interact:
                    # construct the other derivatives wrt. atom j.
                    Rot = my_ia.calc.el.Rot[n]
                    dht2 = np.dot( dht,Rot )
                    dst2 = np.dot( dst,Rot )
                    dh2block = np.outer(phase,dht2.flatten()).reshape(nk,noi,noj,3)
                    ds2block = np.outer(phase,dst2.flatten()).reshape(nk,noi,noj,3)
                    dH0[:,c:d,a:b,:] += dh2block.transpose((0,2,1,3)).conjugate()
                    dS[ :,c:d,a:b,:] += ds2block.transpose((0,2,1,3)).conjugate()

            if i != j and ij_interact:
                # Hermitian (and anti-Hermitian) conjugates;
                # only if matrix block non-zero
                # ( H(k) and S(k) can be (anti)symmetrized as a whole )
                # TODO: symmetrization should be done afterwards
                H0[ :,c:d,a:b]   =  H0[ :,a:b,c:d].transpose((0,2,1)).conjugate()
                S[  :,c:d,a:b]   =  S[  :,a:b,c:d].transpose((0,2,1)).conjugate()
        
            submatAffected = np.zeros(H0[0,:,:].shape)
            submatAffected[a:b,c:d] = 1.0
            plt.clf()
            #plt.matshow(H0[0,:,:].real)
            #plt.imshow(np.abs(H0[0,:,:]), vmin=0.0, vmax=0.1)
            #plt.imshow(np.abs(H0[0,:,:]))
            plt.imshow(submatAffected)
            plt.colorbar()
            plt.title(f"(si={si},sj={sj}) (o1i={o1i},o1j={o1j})")
            plt.savefig(f"IMG_Ham_{i}_{j}.png", dpi=150)


    if kpts is None:
        my_ia.H0, my_ia.S, my_ia.dH0, my_ia.dS = H0, S, dH0, dS

    # print info about Hamiltonian
    if my_ia.first:
        nonzero = sum( abs(S[0].flatten())>1E-15 )
        my_ia.calc.out('Hamiltonian ~%.3f %% filled.' %(nonzero*100.0/norb**2) )
        my_ia.first = False

    print("\n</div> EXIT debug_get_matrices\n")
    return H0, S, dH0, dS