def calculate_density():
    """ Calculate the radial electron density.; sum_nl |Rnl(r)|**2/(4*pi) """
    # ffr: change dens to rho
    rho = np.zeros_like(rgrid)
    for n,l,nl in list_states(maxn, maxl, occu):
        rho += occu[nl] * unlg[nl]**2

    # Nel: local var nel
    Nel = radialgrid_integrate(rho, rgrid)
    if abs( Nel - nel ) > 1E-10:
        raise RuntimeError('Integrated density %.3g, number of electrons %.3g' % (Nel,nel) )
    rho = rho/(4*np.pi*rgrid**2)

    return rho
