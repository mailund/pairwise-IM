'''

Concrete parameterized IM models.


'''

from likelihood import Likelihood
from scipy.optimize import fmin

class FourParamIM(object):
    def __init__(self, ts, loci_20, loci_11, loci_02):
        self.ts = ts
        self.lik = Likelihood(ts, loci_20, loci_11, loci_02)
    
    def log_likelihood(self, c1, c2, m12, m21):
        '''Compute the log likelihood for the four parameter model.'''
        c1s = [c1] * len(self.ts)
        c2s = [c2] * len(self.ts)
        m12s = [m12] * len(self.ts)
        m21s = [m21] * len(self.ts)
        return self.lik.log_likelihood(c1s, c2s, m12s, m21s)
    
    def __call__(self, c1, c2, m12, m21):
        return self.log_likelihood(c1, c2, m12, m21)
        
    def mle_parameters(self, c1, c2, m12, m21):
        ''' Estimate MLE parameters.  Provided parameters are used as start values.'''
        
        def function(x):
            c1, c2, m12, m21 = x
            return -self.log_likelihood(c1, c2, m12, m21)
        
        return fmin(function, [c1,c2,m12,m21])
