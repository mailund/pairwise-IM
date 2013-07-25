'''

Concrete parameterized IM models.


'''

from likelihood import Likelihood
from scipy.optimize import fmin
from scipy import linspace

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

class EpochsIM(object):
    def __init__(self, epsilon, tsplits, n, loci_20, loci_11, loci_02):
        ''' Build a model with n intervals per epochs, epochs separated by tsplits.'''
        self.no_epochs = len(tsplits)
        self.no_intervals = n
        ts = list(linspace(epsilon, tsplits[0], num=n))
        for i in xrange(1,len(tsplits)):
            ts.extend(linspace(tsplits[i-1]+epsilon, tsplits[i], num=n))
        self.lik = Likelihood(ts, loci_20, loci_11, loci_02)
    
    def log_likelihood(self, c1s, c2s, m12s, m21s):
        '''Compute the log likelihood for the four parameter model.'''
        assert len(c1s) == self.no_epochs
        assert len(c2s) == self.no_epochs
        assert len(m12s) == self.no_epochs
        assert len(m21s) == self.no_epochs
        
        def rep(lst, n):
            result = []
            for elm in lst:
                result.extend([elm] * n)
            return result
            
        IM_c1s = rep(c1s, self.no_intervals)
        IM_c2s = rep(c2s, self.no_intervals)
        IM_m12s = rep(m12s, self.no_intervals)
        IM_m21s = rep(m21s, self.no_intervals)
        
        return self.lik.log_likelihood(IM_c1s, IM_c2s, IM_m12s, IM_m21s)
    
    def __call__(self, c1s, c2s, m12s, m21s):
        return self.log_likelihood(c1s, c2s, m12s, m21s)
    
    def _wrap_params(self, c1s, c2s, m12s, m21s):
        return c1s + c2s + m12s + m21s
    
    def _unwrap_params(self, x):
        c1s = x[:self.no_epochs]
        c2s = x[self.no_epochs:2*self.no_epochs]
        m12s = x[2*self.no_epochs:3*self.no_epochs]
        m21s = x[3*self.no_epochs:]
        return c1s, c2s, m12s, m21s
    
    def mle_parameters(self, c1s, c2s, m12s, m21s):
        ''' Estimate MLE parameters.  Provided parameters are used as start values.'''
        
        def function(x):
            if min(x) < 0:
                return float('Inf')
            c1, c2, m12, m21 = self._unwrap_params(x)
            return -self.log_likelihood(c1, c2, m12, m21)
        
        def callback(x):
            c1s,c2s,m12s,m21s = self._unwrap_params(x)
            print 'c1s:', c1s
            print 'c2s:', c2s
            print 'm12s:', m12s
            print 'm21s:', m21s
            print
        
        return self._unwrap_params(fmin(function, 
                                        self._wrap_params(c1s,c2s,m12s,m21s),
                                        maxiter=1000, callback=callback))
