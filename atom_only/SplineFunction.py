class SplineFunction:
    def __init__(self,x,y,k=3,s=0,name=None):
        """ Simple B-spline function; order k is cubic by default.

        Parameters:
        -----------
        x:  x-grid
        y:  y-values for given x-grid points
        k:  order of spline (cubic by default)
        s:  smoothness parameters (means exact reproduction of x,y values)
        name: name of the function

        """
        if s==-1:
            s=len(x)-np.sqrt(2.0*len(x))
        self.tck = splrep(x,y,s=s,k=k)
        self.x = x
        self.y = y
        self.a = x[0]
        self.b = x[-1]
        self.M = len(y)
        self.name = name

    def __call__(self,x,der=0):
        """ Return der'th derivative of f(x)

        Return zero if x beyond the original grid range.
        """
        if isinstance(x, np.ndarray):
            return np.where(x < self.x[0],
                            np.zeros(len(x)),
                            np.where(x > self.x[-1],
                                     np.zeros(len(x)),
                                     splev(x, self.tck, der=der)
                                     )
                            )
        else:
            if x < self.x[0] or x > self.x[-1]:
                return 0.0
            else:
                return splev(x, self.tck, der=der)

    def get_name(self):
        """ Return the name of the function. """
        return self.name

    def get_range(self):
        return (self.x[0],self.x[-1])

    def limits(self):
        return self.get_range()

    def solve(self,y,a=None,b=None):
        """ Solve x for f(x)=y, where x in [a,b]. """
        if a==None: a=self.x[0]
        if b==None: b=self.x[-1]
        assert a<b
        return brentq(lambda x:self(x)-y,a=a,b=b)

    def integrate(self,a,b):
        """
        Integrate given function within [a,b]
        """
        return splint(a,b,self.tck)

    def plot(self,return_pylab=False,der=0,filename=None):
        """
        Plot the function with matplolib
        """
        p=pl
        p.clf()
        from numpy import linspace
        X=linspace(self.a,self.b,self.M*10)
        Y=[self(x,der=der) for x in X]
        p.plot(X,Y)
        if der==0:
            p.scatter(self.x,self.y)
        f1=SplineFunction(self.x,self.y,k=1)
        p.plot(X,[f1(x,der=der) for x in X])
        if return_pylab:
            return p
        elif filename!=None:
            p.savefig(filename)
        else:
            p.show()
        p.close()

    def max_deviation_from_linear(self):
        """
        For given spline (default cubic), return maximum difference
        wrt linear (k=0) interpolation.
        """
        from numpy import array,linspace,abs
        min_dx=min( array(self.x[1:])-array(self.x[0:-1]) )
        M=10*(self.b-self.a)/min_dx
        X=linspace(self.a,self.b,M)
        f1=SplineFunction(self.x,self.y,k=1)
        return max( array([abs(f1(x)-self(x)) for x in X]) )

    def smoothness(self):
        """
        Return a measure for the ''smoothness'' of interpolation.

        It is measured by max_deviation_from_linear/average|dy|.
        Smooth interpolations should have <<1.
        If ~1, then interpolation deviates as much as is the variations
        in y and interpolation is not smooth.
        """
        from numpy import abs,average,array
        avg=average( abs(array(self.y[1:])-array(self.y[0:-1])) )
        return self.max_deviation_from_linear()/avg
