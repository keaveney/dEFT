import numpy as np
import emcee
import time

start = time.time()

bins = [-1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5]

coefficients = ["ctg", "ctw", "ctphi", "ctb", "cphit", "cphiQ1", "cphiQ3"]

preds = {'SM': np.array([3.0,3.0,3.0, 3.0,3.0,3.0, 3.0]),
    'ctg-': np.array([2.0,3.0,3.0,3.0,3.0,3.0,3.0]),
    'ctg+': np.array([4.0,3.0,3.0,3.0,3.0,3.0,3.0]),
    'ctw-': np.array([3.0,2.0,3.0,3.0,3.0,3.0,3.0]),
    'ctw+': np.array([3.0,4.0,3.0,3.0,3.0,3.0,3.0]),
    'ctphi-': np.array([3.0,3.0,2.0,3.0,3.0,3.0,3.0]),
    'ctphi+': np.array([3.0,3.0,4.0,3.0,3.0,3.0,3.0]),
    'ctb-': np.array([3.0,3.0,3.0,2.0,3.0,3.0,3.0]),
    'ctb+': np.array([3.0,3.0,3.0,4.0,3.0,3.0,3.0]),
    'cphit-': np.array([3.0,3.0,3.0,3.0,2.0,3.0,3.0]),
    'cphit+': np.array([3.0,3.0,3.0,3.0,4.0,3.0,3.0]),
    'cphiQ1-': np.array([3.0,3.0,3.0,3.0,3.0,2.0,3.0]),
    'cphiQ1+': np.array([3.0,3.0,3.0,3.0,3.0,4.0,3.0]),
    'cphiQ3-': np.array([3.0,3.0,3.0,3.0,3.0,3.0,2.0]),
    'cphiQ3+': np.array([3.0,3.0,3.0,3.0,3.0,3.0,4.0])
}

data = [5.0, 5.0, 3.0, 3.0, 3.0, 3.0, 3.0]
#icov = [[1.5, -0.49, -0.5], [-0.49, 1.5,  -0.5], [-0.50, -0.5,1.5]]
cov = [[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],[0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],[0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]]
icov = np.linalg.inv(cov)

pred_basis = {}
        
#isolate SM--c_i term
# SM--c_i interference term


def make_pred(c):
    pred = [0.0,0.0,0.0,0.0,0.0,0.0,0.0]
    for ci in range(0, len(coefficients)):
        basis_c_sm_name = coefficients[ci] + "_sm"
        linear_contrib = (c[ci]*pred_basis[basis_c_sm_name])
        pred += linear_contrib
    pred = pred + preds['SM']
    return pred

print "SM pred is  = " + str(make_pred([0.0,0.0,0.0,0.0,0.0,0.0,0.0]))

def lnprob(c, data, invcov):
    pred = make_pred(c)
    diff = pred - data
    #print" c = " +str(c) + " pred = " + str(pred) + " data = " + str(data) + " diff = " + str(diff) + " ll = " + str(-np.dot(diff,np.dot(invcov,diff))/2.0)
    return -np.dot(diff,np.dot(invcov,diff))/2.0

ndim = 7

nwalkers = 200

p0 = np.random.rand(ndim * nwalkers).reshape((nwalkers, ndim))

sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[data, icov])

pos, prob, state = sampler.run_mcmc(p0, 100)

#print "p0 = " + str(p0) + " pos = " + str(pos)

sampler.reset()
sampler.run_mcmc(pos, 500)
samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

import matplotlib.pyplot as pl

#for i in range(ndim):
# pl.figure()
#pl.hist(sampler.flatchain[:,i], 300, color="k", histtype="step")
#pl.title("Dimension {0:d}".format(i))
#pl.show()

# red dashes, blue squares and green triangles
#plt.plot(t, t, 'r--', t, t**2, 'bs', t, t**3, 'g^')
#pl.plot(bins, data, 'ro')
#pl.errorbar(bins, data, 'ro', xerr=0.2, yerr=0.4)

pl.figure()

print "len xs = " + str(len(bins)) +  " len ys = " + str(make_pred([1.0,0.0,0.0,0.0,0.0,0.0,0.0]))

pl.errorbar(bins, make_pred([1.0,0.0,0.0,0.0,0.0,0.0,0.0]), xerr=0.0, yerr=0.0, label='ctg=1.0')
pl.errorbar(bins, make_pred([0.0,1.0,0.0,0.0,0.0,0.0,0.0]), xerr=0.0, yerr=0.0, label='ctw=1.0')
pl.errorbar(bins, make_pred([0.0,0.0,1.0,0.0,0.0,0.0,0.0]), xerr=0.0, yerr=0.0, label='ctphi=1.0')
pl.errorbar(bins, make_pred([0.0,0.0,0.0,1.0,0.0,0.0,0.0]), xerr=0.0, yerr=0.0, label='ctb=1.0')
pl.errorbar(bins, make_pred([0.0,0.0,0.0,0.0,1.0,0.0,0.0]), xerr=0.0, yerr=0.0, label='cphit=1.0')
pl.errorbar(bins, make_pred([0.0,0.0,0.0,0.0,0.0,1.0,0.0]), xerr=0.0, yerr=0.0, label='ctphiQ1=1.0')
pl.errorbar(bins, make_pred([0.0,0.0,0.0,0.0,0.0,0.0,1.0]), xerr=0.0, yerr=0.0, label='ctphiQ3=1.0')
pl.errorbar(bins, preds['SM'], xerr=0.0, yerr=0.0, label='SM')
pl.errorbar(bins, data, fmt="o",xerr=0.25, yerr=0.05, label='Data')

pl.axis([-1.0, 1.0, -5.0, 10.0])

pl.legend()
pl.show()

import corner
fig = corner.corner(samples, labels=["$ctg$", "$ctw$", "$ctphi$", "$ctb$", "$cphit$","$ctphiQ1$","$ctphiQ3$"],
                    truths=[0.0, 0.0, 0.0,0.0, 0.0, 0.0,0.0])

fig.show()
fig.savefig("triangle.png")

end = time.time()
print "Total elapsed wall time  = " +  str(end - start)
