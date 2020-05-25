import numpy as np
import matplotlib.pyplot as plt
import itertools
from mpl_toolkits.mplot3d import Axes3D

def nSamples(n):
    return n*(n+1)/2

#samples should be an m x n matrix: each row is a sample, each column is a coupling value
nOps = 2 #number of operators being considered: important so that we know how to construct/interpret the coupling matrices

samples = np.array(  [[1.0, 0.0, 0.0], [1.0, 2.0, 0.0], [1.0, 2.0, 50.0], [1.0, 2.0, -50.0], [1.0, -2.0, 50.0], [1.0, -2.0, -50.0] ]  )
#samples = {1.0, 0.0, 0.0}, {1.0, 2.0, 0.0}, {1.0, 2.0, 50.0}, {1.0, 2.0, -50.0}, {1.0, -2.0, 50.0}, {1.0, -2.0, -50.0},

# the n values of c_i that define the 'input' distributions: c_sm, c_i1...c_in
#in this input format, each row corresponds to the couplings for a given input sample
# ordering SM^2, SM*C1...SM*Cn, Ci*Cj...

print("n prediction required is  = " + str(nSamples(nOps+1)))

#cInputAr is generated automatically from 'samples'
nTerms = 1 + (2*nOps) + len(list(itertools.combinations(samples[0],2))) - nOps

print("n terms is  = " + str(nTerms))

cInputAr = np.array([])

for sample in samples:
    couplings = np.array([])
    combs = list(itertools.combinations(sample,2))
    #SM term
    couplings = np.append(couplings, sample[0]**2)
    #'SM--Oi' interference terms
    for comb in range(0, nOps):
        couplings = np.append(couplings, combs[comb][0]*combs[comb][1])
    #'squared' terms
    for ci in range(1, nOps+1):
        couplings = np.append(couplings, sample[ci]**2)
    #'cross' terms
    for comb in range(nOps, len(combs)):
        couplings = np.append(couplings, combs[comb][0]*combs[comb][1])
    cInputAr = np.append(cInputAr, couplings , axis=0)

cInputAr = np.reshape(cInputAr,(nTerms,nTerms))
#cInputAr = np.array([[1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [1.0, 2.0, 0.0, 4.0, 0.0, 0.0], [1.0, 2.0, 50.0, 4.0, 2500.0, 100.0],[1.0, 2.0, -50.0, 4.0, 2500.0, -100.0], [1.0, -2.0, 50.0, 4.0, 2500.0, -100.0], [1.0, -2.0, -50.0, 4.0, 2500.0, 100.0]])
cInputAr = np.transpose(cInputAr)
inWeightsAr = np.linalg.inv(cInputAr)
Tin = np.array([ [20.46], [34.54], [158.9], [151.6], [135.6], [132.2]])
TinT = Tin.transpose()

# 'sanity check' example
ciAr = np.array([])
cjAr = np.array([])
#xs = np.array([])

def calcXS(ci, cj):
    ciTargetAr = np.array([[1.0], [ci], [cj], [ci**2], [cj**2], [ci*cj]  ])
    #ciTargetAr = np.array([])
    morphWeightsAr = inWeightsAr.dot(ciTargetAr)
    Tout = TinT.dot(morphWeightsAr)
    #print("ci  = " + str(ci) + ", cj = " + str(cj))
    return Tout

print("val prediction for ci = 3.65, cj = 17.12 = " + str(calcXS(3.65, 17.12)))
print("MadGraph prediction for ci = 3.65, cj = 17.12 = 66.55")

vCalcXS = np.vectorize(calcXS)

ciAr = np.arange(-8, 5, 0.5)
cjAr = np.arange(-40, 40, 0.5)
X, Y = np.meshgrid(ciAr, cjAr)

plt.figure()

xs = vCalcXS(X,Y)

print("xs shape = " + str(xs.shape))
cols = np.arange(0, 85, 1)
contour_filled = plt.contourf(ciAr, cjAr, xs, levels=cols)
cbar = plt.colorbar(contour_filled)
cbar.set_label(r"$\sigma_{t\bar{t}}$ [pb]", rotation=90)

plt.xlabel(r"$c_tG$")
plt.ylabel(r"$c^{3,8}_Qq$")

plt.show()

#valPtsX = [-3.31, 5.72]
#valPtsY = [11.17, 77.03]

#inPtsX = [0.0, 2.0, -2.0]
#inPtsY = [20.76, 34.38, 12.71]

#fig, ax = plt.subplots()
#ax.set_xlabel('ctg')
#ax.set_ylabel('cross section [pb]')

#plt.plot(x,y, label='morphing model')
#plt.plot(inPtsX, inPtsY, "b+", markersize=7, label='input predictions')
#plt.plot(valPtsX, valPtsY, "r+", markersize=7,label='validation predictions')

#legend = plt.legend(loc='upper center', shadow=True, fontsize='x-large')

#plt.show()
