class XC_PW92:
    def __init__(self):
        """ The Perdew-Wang 1992 LDA exchange-correlation functional. """
        self.small=1E-90
        self.a1 = 0.21370
        self.c0 = 0.031091
        self.c1 = 0.046644
        self.b1 = 1.0/2.0/self.c0*np.exp(-self.c1/2.0/self.c0)
        self.b2 = 2*self.c0*self.b1**2
        self.b3 = 1.6382
        self.b4 = 0.49294

    def exc(self,n,der=0):
        """ Exchange-correlation with electron density n. """
        if n<self.small:
            return 0.0
        else:
            return self.e_x(n,der=der)+self.e_corr(n,der=der)

    def e_x(self,n,der=0):
        """ Exchange. """
        if der==0:
            return -3.0/4*(3*n/pi)**(1.0/3)
        elif der==1:
            return -3.0/(4*pi)*(3*n/pi)**(-2.0/3)

    def e_corr(self,n,der=0):
        """ Correlation energy. """
        rs = (3.0/(4*pi*n))**(1.0/3)
        aux=2*self.c0*( self.b1*sqrt(rs)+self.b2*rs+self.b3*rs**(3.0/2)+self.b4*rs**2 )
        if der==0:
            return -2*self.c0*(1+self.a1*rs)*log(1+aux**-1)
        elif der==1:
            return ( -2*self.c0*self.a1*log(1+aux**-1) \
                   -2*self.c0*(1+self.a1*rs)*(1+aux**-1)**-1*(-aux**-2)\
                   *2*self.c0*(self.b1/(2*sqrt(rs))+self.b2+3*self.b3*sqrt(rs)/2+2*self.b4*rs) )*( -(4*pi*n**2*rs**2)**-1 )

    def vxc(self,n):
        """ Exchange-correlation potential (functional derivative of exc). """
        eps=1E-9*n
        if n<self.small:
            return 0.0
        else:
            return self.exc(n)+n*self.exc(n,der=1)

            
