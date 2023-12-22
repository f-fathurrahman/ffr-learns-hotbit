def calculate_energies():
    """
    Calculate energy contributions.
    """
    bs_energy = 0.0
    for n,l,nl in list_states(maxn, maxl, occu):
        bs_energy += occu[nl]*enl[nl]

    #exc = array([self.xcf.exc(self.dens[i]) for i in xrange(self.N)])
    for i in range(N):
        exc[i] = exc_PW92(dens[i])

    Hartree_energy = radialgrid_integrate( VHartree*dens, rgrid, use_dV=True )/2
    vxc_energy = radialgrid_integrate( vxc*dens, rgrid, use_dV=True )
    exc_energy = radialgrid_integrate( exc*dens, rgrid, use_dV=True )
    confinement_energy = radialgrid_integrate( conf*dens, rgrid, use_dV=True )
    total_energy = bs_energy- Hartree_energy - vxc_energy + exc_energy

    print('\n\nEnergetics:')
    print('-------------')
    print('\nsingle-particle energies')
    print('------------------------')
    for n,l,nl in list_states(maxn, maxl, occu):
        print('%s, energy %22.15f' % (nl,enl[nl]))

    print('\nvalence orbital energies')
    print('--------------------------')
    for nl in atoms_data[symbol]['valence_orbitals']:
        print('%s, energy %22.15f' %(nl, enl[nl]))

    print('\n')
    print('total energies:')
    print('---------------')
    print('sum of eigenvalues:     %22.15f' % bs_energy)
    print('Hartree energy:         %22.15f' % Hartree_energy)
    print('vxc correction:         %22.15f' % vxc_energy)
    print('exchange + corr energy: %22.15f' % exc_energy)
    print('----------------------------')
    print('total energy:           %22.15f\n\n' % total_energy)
