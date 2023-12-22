class FastSplineFunction:
    def __init__(self,x,y,grid,**kwargs):
        """ Initialize second derivatives.

        Parameters:
        -----------
        x: x-grid
        y: y-values on given x-grid
        grid: type of grid
            'linear':   x[i]=xmin + i/(N-1)*(xmax-xmin) (i~(N-1)*(x-xmin)/(xmax-xmin))
            'power':  x[i]=xmin + (i/(N-1))**p*(xmax-xmin)
        **kwargs: parameters for grid (can vary) (xmin,xmax,p,...)
            (N comes from length of x)
        """


        self.x=x
        self.y=y
        self.a=x[0]
        self.b=x[-1]
        d2=np.zeros_like(y)
        u=np.zeros_like(y)
        n=len(x)
        qn       = 0.0
        un       = 0.0

        # Solve the second derivatives
        for i in range(1,n-1):
            sig  = (x[i]-x[i-1])/(x[i+1]-x[i-1])
            p    = sig*d2[i-1]+2
            d2[i]= (sig-1)/p
            u[i] = ( 6*((y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]))/(x[i+1]-x[i-1])-sig*u[i-1] )/p

        d2[n-1]=(un-qn*u[n-2])/(qn*d2[n-2]+1)
        for i in range(n-1)[::-1]:
            d2[i]=d2[i]*d2[i+1]+u[i]
        self.n=n
        self.d2=d2
        self.grid=grid
        self.xmin, self.xmax=x[0], x[-1]
        if grid=='linear':
            #self.xmin,self.xmax=kwargs['xmin'],kwargs['xmax']
            self.find_lower_index=self.find_lower_index_linear
        if grid=='power':
            self.p=kwargs['p']
            self.find_lower_index=self.find_lower_index_power

    def get_range(self):
        return self.x[0],self.x[-1]

    def find_lower_index_linear(self,x):
        """ For given x, return i such that x_i<=x<x_i+1 """
        return int( (x-self.xmin)/(self.xmax-self.xmin)*(self.n-1) )

    def find_lower_index_power(self,x):
        """ For given x, return i such that x_i<=x<x_i+1 """
        return int( ((x-self.xmin)/(self.xmax-self.xmin))**(1.0/self.p)*(self.n-1) )

    def __call__(self,x,der=0):
        """! If x is outside the interpolation range, return 0."""

        if x<self.x[0] or x>=self.x[-1]:
            return 0.0

        lo=self.find_lower_index(x)
        hi=lo+1

        h=self.x[hi]-self.x[lo]
        assert self.x[lo]<=x<self.x[hi]
        a, b=(self.x[hi]-x)/h, (x-self.x[lo])/h

        if der==0:
            return a*self.y[lo] + b*self.y[hi] + ((a**3-a)*self.d2[lo]+(b**3-b)*self.d2[hi])*(h**2)/6
        elif der==1:
            return (self.y[hi]-self.y[lo])/h - (3*a**2-1)/6*h*self.d2[lo] + (3*b**2-1)/6*h*self.d2[hi]
        elif der==2:
            return a*self.d2[lo] + b*self.d2[hi]
