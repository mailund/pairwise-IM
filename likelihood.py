'''

Code for computing the likelihood of a model given a set of loci.


'''

from IMSystem import IMSystem
from scipy import log

class Likelihood:

    def __init__(self, ts, loci_20, loci_11, loci_02):
        '''Build a model based on a time discritization and a set of loci...'''
        self.ts = ts
        self.loci_20 = loci_20
        self.loci_11 = loci_11
        self.loci_02 = loci_02
        self.loci_20_likelihoods = [l.likelihood(ts) for l in self.loci_20]
        self.loci_11_likelihoods = [l.likelihood(ts) for l in self.loci_11]
        self.loci_02_likelihoods = [l.likelihood(ts) for l in self.loci_02]
        
    def log_likelihood(self, c1s, c2s, m12s, m21s):
        model = IMSystem(self.ts, c1s, c2s, m12s, m21s)
        pdm_20,pdm_11,pdm_02 = model.coalescence_distribution()

        logL = 0.0
        for l in self.loci_20_likelihoods:
            logL += log(sum(l * pdm_20))
        for l in self.loci_11_likelihoods:
            logL += log(sum(l * pdm_11))
        for l in self.loci_02_likelihoods:
            logL += log(sum(l * pdm_02))

        return logL


if __name__ == '__main__':
    from scipy import linspace
    ts = linspace(0.1,4)

    c1s = [1] * len(ts)
    c2s = [2] * len(ts)
    m12s = [0.1] * len(ts)
    m21s = [0.2] * len(ts)

    from loci import IndicatorLocus
    l1 = IndicatorLocus(0.01)
    l2 = IndicatorLocus(0.11)
    l3 = IndicatorLocus(4.10)

    loci_20 = [l1,l2,l3]
    loci_11 = [l1,l2,l3]
    loci_02 = [l1,l2,l3]
    
    lik = Likelihood(ts, loci_20, loci_11, loci_02)
    print lik.log_likelihood(c1s, c2s, m12s, m21s)

