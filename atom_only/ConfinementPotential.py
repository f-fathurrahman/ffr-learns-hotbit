class ConfinementPotential:
    def __init__(self,mode,**kwargs):
        self.mode=mode
        if mode=='none':
            self.f=self.none #lambda r:0.0
            self.comment='none'
        elif mode=='quadratic':
            self.r0=kwargs['r0']
            self.f=self.quadratic #lambda r:(r/self.r0)**2
            self.comment='quadratic r0=%.3f' %self.r0
        elif mode=='general':
            self.r0=kwargs['r0']
            self.s=kwargs['s']
            self.f=self.general #lambda r:(r/self.r0)**s
            self.comment='general r0=%.3f s=%.3f' %self.r0, self.s
        else:
            raise NotImplementedError('implement new confinements')

    def get_comment(self):
        return self.comment

    def none(self,r):
        return 0.0

    def quadratic(self,r):
        return (r/self.r0)**2

    def general(self,r):
        return (r/self.r0)**self.s

    def __call__(self,r):
        return self.f(r)
