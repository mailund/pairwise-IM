'''

Representation of loci. They must somehow be able to give a likelihood for a coalescence time, P(D | t), however that is done, which will then be combined with the coalescence probabilities computed by the IMSystem.

'''

from numpy import zeros

class Locus(object):
    '''Abstract interface for loci.
    
    Depending on the model, the way loci computes P(D | t) will vary, but this
    is the interface they must implement.
    '''
    
    def __init__(self):
        ''' Create a locus that can compute likelihoods for the discritized
        time space.'''
        pass
                
    def likelihood(self, ts):
        ''' Return a vector of probabilities P(D | t[i]) for i an index over intervals
        where P(D | t[i]) should be the likelihood of the locus data if the coalescence
        time for the locus was in t[i-1] to t[i] (with t[-1] == 0).'''
        pass
        


class IndicatorLocus(Locus):
    '''This is a dummy locus that assumes that we actually *know* the coalescence time
    and we thus put 1 on the right interval and 0 elsewhere.  For testing purposes.'''
    
    def __init__(self, coal_time):
        super(IndicatorLocus,self).__init__()
        self.coal_time = coal_time
        
    def _get_coal_interval(self, ts):
        if self.coal_time <= ts[0]:
            # we coalesce in first interval
            return 0
        elif self.coal_time > ts[-1]:
            return len(ts)
        else:
            for i in xrange(0,len(ts)):
                if ts[i-1] < self.coal_time <= ts[i]:
                    return i
        assert False, "We shouldn't reach this point."
    
    def likelihood(self, ts):
        probabilities = zeros(len(ts) + 1)
        probabilities[self._get_coal_interval(ts)] = 1.0
        return probabilities

        
if __name__ == '__main__':
    from scipy import linspace
    ts = linspace(0.1,4, num=5)
    
    l1 = IndicatorLocus(0.01)
    l2 = IndicatorLocus(0.11)
    l3 = IndicatorLocus(4.10)
    
    print ts
    print l1.likelihood(ts)
    print l2.likelihood(ts)
    print l3.likelihood(ts)
