import matplotlib.pyplot as plt

def solve_eigenstates(iteration,itmax_eigen=100):
    """
    Solve the eigenstates for given effective potential.

    u""(r) - 2*(v_eff(r)+l*(l+1)/(2r**2)-e)*u(r)=0
    ( u""(r) + c0(r)*u(r) = 0 )

    r=r0*exp(x) --> (to get equally spaced integration mesh)

    u""(x) - u"(x) + c0(x(r))*u(r) = 0
    """

    dx = xgrid[1] - xgrid[0]
    c2 = np.ones(N)
    c1 = -np.ones(N)
    d_enl_max = 0.0
    itmax_eigen = 0

    for n, l, nl in list_states(maxn,maxl,occu):
        nodes_nl = n-l-1
        #
        if iteration==0:
            eps = -1.0*Z**2/n**2
        else:
            eps=enl[nl]

        if iteration <= 3:
            delta=0.5*Z**2/n**2  #previous!!!!!!!!!!
        else:
            delta = d_enl[nl]

        direction = "none"
        epsmax = veff[-1] - l*(l+1)/(2*rgrid[-1]**2)
        it = 0
        u = np.zeros(N)
        hist = []

        while True:
            eps0 = eps
            c0, c1, c2 = construct_coefficients(l, eps)

            # boundary conditions for integration from analytic behaviour (unscaled)
            # u(r)~r**(l+1)   r->0
            # u(r)~exp( -sqrt(c0(r)) ) (set u[-1]=1 and use expansion to avoid overflows)
            u[0:2] = rgrid[0:2]**(l+1)

            if not( c0[-2]<0 and c0[-1]<0 ):
                plt.plot(c0)
                plt.show()

            assert c0[-2]<0 and c0[-1]<0

            u, nodes, A, ctp = shoot(u,dx,c2,c1,c0,N)
            it+=1
            norm = radialgrid_integrate(u**2, rgrid)
            u = u/sqrt(norm)

            if nodes>nodes_nl:
                # decrease energy
                if direction=="up": delta/=2
                eps -= delta
                direction = "down"
            elif nodes < nodes_nl:
                # increase energy
                if direction=="down": delta/=2
                eps+=delta
                direction="up"
            elif nodes==nodes_nl:
                shift=-0.5*A/(rgrid[ctp]*norm)
                if abs(shift)<1E-8: #convergence
                    break
                if shift>0:
                    direction="up"
                elif shift<0:
                    direction="down"
                eps+=shift
            if eps>epsmax:
                eps=0.5*(epsmax+eps0)
            hist.append(eps)

            if it > 100:
                print("Epsilon history for %s" % nl)
                for h in hist:
                    print(h)
                print("nl=%s, eps=%f" %(nl,eps) )
                print("max epsilon",epsmax)
                raise RuntimeError("Eigensolver out of iterations. Atom not stable?")

        itmax_eigen=max(it,itmax_eigen)
        unlg[nl] = u
        Rnlg[nl] = unlg[nl]/rgrid
        d_enl[nl] = abs(eps-enl[nl])
        d_enl_max = max( d_enl_max, d_enl[nl])
        enl[nl] = eps
        #print("-- state %s, %i eigensolver iterations, e=%9.5f, de=%9.5f" %(nl, it, enl[nl], d_enl[nl]))

        assert nodes==nodes_nl
        assert u[1] > 0.0
    #
    return d_enl_max, itmax_eigen
