class Function:

    def __init__(self,mode,*args,**kwargs):
        if mode == 'spline':
            self.f=SplineFunction(*args,**kwargs)
        elif mode == 'string':
            raise NotImplementedError('todo')
            self.args=args
            self.kwargs=kwargs
        elif mode == 'fastspline':
            self.f=FastSplineFunction(*args,**kwargs)
        else:
            raise NotImplementedError('todo')

    def __call__(self,x,der=0):
        return self.f(x,der)

    def plot(self,der=0,a=None,b=None,npoints=1000,filename=None,return_pylab=False):
        """ Plot the function with matplolib. """
        #import pylab as pl
        #from numpy import linspace
        a0,b0 = self.f.get_range()
        lower = [a0,a][a!=None]
        upper = [b0,b][b!=None]
        X = np.linspace(lower,upper,npoints)
        Y = [self(x,der=der) for x in X]
        plt.plot(X,Y)
        if return_pylab:
            return plt
        elif filename!=None:
            plt.savefig(filename)
        else:
            plt.show()
        plt.close()
