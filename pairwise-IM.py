from numpy import matrix
from scipy.linalg import expm

## Constants used as indices in rate and transition matrices
LINEAGES_IN_SEP_POPS = 0
LINEAGES_IN_POP_1 = 1
LINEAGES_IN_POP_2 = 2
COALESCED = 3
NOT_COALESCED = [0,1,2]

def make_rate_matrix(c1, c2, m12, m21):
    '''Create a rate matrix based on coalescence rates c1 and c2 and
    migration rates m12 and m21.'''
    Q = matrix(
        [
        # State 1: lineages in different populations
        [-(m12+m21), m21, m12, 0],
        # State 2: both lineages in population 1
        [2*m12, -(2*m12+c1), 0, c1],
        # State 3: both lineages in population 2
        [2*m21, 0, -(2*m21+c2), c2],
        # State 4: coalesced (catches both populations; absorbing)
        [0, 0, 0, 0]
        ])
    return Q

class IMSystem:
    '''Wrapping a two-population isolation-with-migration system.'''
    
    def __init__(self, ts, c1s, c2s, m12s, m21s):
        '''Build the system based on end-points of time intervals, ts,
        coalescence rates c1s and c2s and migration rates m12s and m21s.
        '''
        self.ts = ts
        self.c1s = c1s
        self.c2s = c2s
        self.m12s = m12s
        self.m21s = m21s
        
        self.no_intervals = len(ts)
        assert len(self.c1s) == self.no_intervals
        assert len(self.c2s) == self.no_intervals
        assert len(self.m12s) == self.no_intervals
        assert len(self.m21s) == self.no_intervals
        
        self.Qs = [make_rate_matrix(self.c1s[i],self.c2s[i],self.m12s[i],self.m21s[i])
                   for i in xrange(self.no_intervals)]
        self.Ps = [None] * self.no_intervals
        self.Ps[0] = matrix(expm(self.Qs[0] * self.ts[0]))
        for i in xrange(1,self.no_intervals):
            self.Ps[i] = matrix(expm(self.Qs[i] * (self.ts[i]-self.ts[i-1])))

    def coalescence_distribution(self):
        '''Returns the (discritized) coalescence distribution for the time
        intervals. Implicitly the time interval from the last ts till infinity
        is assumed to carry the probability mass that gets this to sum to 1.'''
        pdm_20 = [0] * self.no_intervals
        pdm_11 = [0] * self.no_intervals
        pdm_02 = [0] * self.no_intervals
        
        pdm_20[0] = self.Ps[0][LINEAGES_IN_POP_1,COALESCED]
        pdm_11[0] = self.Ps[0][LINEAGES_IN_SEP_POPS,COALESCED]
        pdm_02[0] = self.Ps[0][LINEAGES_IN_POP_2,COALESCED]
        prevP = self.Ps[0]
        
        for i in xrange(self.no_intervals):
            for state in NOT_COALESCED:
                pdm_20[i] += prevP[LINEAGES_IN_POP_1,state] * self.Ps[i][state,COALESCED]
                pdm_11[i] += prevP[LINEAGES_IN_SEP_POPS,state] * self.Ps[i][state,COALESCED]
                pdm_02[i] += prevP[LINEAGES_IN_POP_2,state] * self.Ps[i][state,COALESCED]
                prevP = prevP * self.Ps[i]
                
        return (pdm_20,pdm_11,pdm_02)

from scipy import linspace
ts = linspace(0.1,4)
c1s = [1] * len(ts)
c2s = [2] * len(ts)
m12s = [0.1] * len(ts)
m21s = [0.1] * len(ts)

im = IMSystem(ts, c1s, c2s, m12s, m21s)
pdm_20,pdm_11,pdm_02 = im.coalescence_distribution()

from matplotlib import pyplot
pyplot.plot(im.ts,pdm_20)
pyplot.plot(im.ts,pdm_11)
pyplot.plot(im.ts,pdm_02)
pyplot.axis([0, max(ts), 0, max([max(pdm_20),max(pdm_11),max(pdm_02)])])
pyplot.show()
