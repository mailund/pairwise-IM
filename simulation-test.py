from loci import IndicatorLocus
from likelihood import Likelihood

class FourParamIM(object):
    def __init__(self, ts, loci_20, loci_11, loci_02):
        self.ts = ts
        self.lik = Likelihood(ts, loci_20, loci_11, loci_02)
    
    def __call__(self, c1, c2, m12, m21):
        c1s = [c1] * len(self.ts)
        c2s = [c2] * len(self.ts)
        m12s = [m12] * len(self.ts)
        m21s = [m21] * len(self.ts)
        return self.lik.log_likelihood(c1s, c2s, m12s, m21s)


loci_20 = [IndicatorLocus(ct) 
           for ct in map(float,open('simulations/coal.2-0.txt').read().split())]
loci_11 = [IndicatorLocus(ct) 
           for ct in map(float,open('simulations/coal.1-1.txt').read().split())]
loci_02 = [IndicatorLocus(ct) 
           for ct in map(float,open('simulations/coal.0-2.txt').read().split())]
           

from scipy import linspace
ts = linspace(0.1,10)

model = FourParamIM(ts, loci_20, loci_11, loci_02)
cs = linspace(0.1,1.5)
c1_curve = [model(c1,0.5,0.1,0.2) for c1 in cs]
c2_curve = [model(1,c2,0.1,0.2) for c2 in cs]

from matplotlib import pyplot as plt

plt.plot(cs, c1_curve)
plt.plot(cs, c2_curve)
#plt.show()

ms = linspace(0.0,1.0)
m12_curve = [model(1,0.5,m12,0.2) for m12 in ms]
m21_curve = [model(1,0.5,0.1,m21) for m21 in ms]

plt.plot(ms, m12_curve)
plt.plot(ms, m21_curve)
plt.show()
