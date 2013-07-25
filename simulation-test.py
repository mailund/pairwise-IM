from loci import IndicatorLocus
from models import FourParamIM, EpochsIM



loci_20 = [IndicatorLocus(ct) 
           for ct in map(float,open('simulations/coal.2-0.txt').read().split())]
loci_11 = [IndicatorLocus(ct) 
           for ct in map(float,open('simulations/coal.1-1.txt').read().split())]
loci_02 = [IndicatorLocus(ct) 
           for ct in map(float,open('simulations/coal.0-2.txt').read().split())]
           

from scipy import linspace
ts = linspace(0.1,10)

#model = FourParamIM(ts, loci_20, loci_11, loci_02)
#print model.mle_parameters(1.0, 0.5, 0.1, 0.2)

model = EpochsIM(1e-5, [1,2,5], 5, loci_20, loci_11, loci_02)
c1s, c2s, m12s, m21s = [1.0]*3, [0.5]*3, [0.1]*3, [0.2]*3
print model.log_likelihood(c1s, c2s, m12s, m21s)
wrapped = model._wrap_params(c1s, c2s, m12s, m21s)
print wrapped
print model._unwrap_params(wrapped)
print

print model.mle_parameters(c1s,c2s,m12s,m21s)

#cs = linspace(0.1,1.5)
#c1_curve = [model(c1,0.5,0.1,0.2) for c1 in cs]
#c2_curve = [model(1,c2,0.1,0.2) for c2 in cs]

# from matplotlib import pyplot as plt
# 
# plt.plot(cs, c1_curve)
# plt.plot(cs, c2_curve)
# #plt.show()
# 
# ms = linspace(0.0,1.0)
# m12_curve = [model(1,0.5,m12,0.2) for m12 in ms]
# m21_curve = [model(1,0.5,0.1,m21) for m21 in ms]
# 
# plt.plot(ms, m12_curve)
# plt.plot(ms, m21_curve)
# plt.show()
