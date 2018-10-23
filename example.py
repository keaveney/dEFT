from tools import splash
from tools import pred_builder
from asciimatics.screen import Screen
import numpy as np
import itertools
import matplotlib.pyplot as pl
import corner
import emcee
import time

#intro graphics
Screen.wrapper(splash.deft_splash)

#start time
start = time.time()

######################################################
###############   DEFINE THE DATA   ##################
######################################################
#bins
bins = np.array([-1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5])

#data
data = np.array([5.0, 5.0, 3.0, 3.0, 3.0, 3.0, 3.0])

#covariance matrix and inverse
cov = [[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],[0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],[0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]]
icov = np.linalg.inv(cov)

######################################################
###############   DEFINE THE MODEL   #################
######################################################

#list of coeffieicnts/operators to be considered
coefficients = ["ctg", "ctw", "ctphi", "ctb", "cphit", "cphiQ1", "cphiQ3"]

#maximum power of c_i/Lambda^2 to be considered
# this can only be 1 or 2
# 1: Lambda^{-2} - interference between diagrams with single dim6 operator and SM diagrams
# 2: Lambda^{-2} + Lambda^{-4} - the above with the 'squared' terms (diagrams with single dim6 operator squared)
#                                  and 'cross' terms (interefernece between pairs of diagrams with single dim6 operator each)
max_coeff_power = 1

#the tricky part - list of predictions for data and for +/- X for each coefficient
# if max_coeff_power = 1, and more than one operator is being considred,
# this dictonary also needs predictions for (ci,cj) = (X,X)
# for all pairs of operators

preds = {
    'SM': np.array([3.0,3.0,3.0, 3.0,3.0,3.0, 3.0]),
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

#define prediction builder function for convenience
pb = pred_builder.pred_builder()
pred_basis = pb.make_basis(preds,coefficients,max_coeff_power)

######################################################
###############   DEFINE THE FIT/SCAN   ##############
######################################################

#define define the objective function that will be minimised/maximised
def lnprob(c, data, icov):
    pred = pb.make_pred(c)
    diff = pred - data
    #print" c = " +str(c) + " pred = " + str(pred) + " data = " + str(data) + " diff = " + str(diff) + " ll = " + str(-np.dot(diff,np.dot(invcov,diff))/2.0)
    return -np.dot(diff,np.dot(icov,diff))/2.0

#define define the objective function that will be minimised/maximiseds
nwalkers = 100
ndim = len(coefficients)

p0 = np.random.rand(ndim * nwalkers).reshape((nwalkers, ndim))

sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[data, icov])

pos, prob, state = sampler.run_mcmc(p0, 100)

sampler.reset()

sampler.run_mcmc(pos, 500)

samples = sampler.chain[:, 50:, :].reshape((-1, ndim))


######################################################
###############   MAKE RESULTS PLOTS  ################
######################################################
pl.figure()

pl.errorbar(bins, pb.make_pred([1.0,0.0,0.0,0.0,0.0,0.0,0.0]), xerr=0.0, yerr=0.0, label='ctg=1.0')
pl.errorbar(bins, pb.make_pred([0.0,1.0,0.0,0.0,0.0,0.0,0.0]), xerr=0.0, yerr=0.0, label='ctw=1.0')
pl.errorbar(bins, pb.make_pred([0.0,0.0,1.0,0.0,0.0,0.0,0.0]), xerr=0.0, yerr=0.0, label='ctphi=1.0')
pl.errorbar(bins, pb.make_pred([0.0,0.0,0.0,1.0,0.0,0.0,0.0]), xerr=0.0, yerr=0.0, label='ctb=1.0')
pl.errorbar(bins, pb.make_pred([0.0,0.0,0.0,0.0,1.0,0.0,0.0]), xerr=0.0, yerr=0.0, label='cphit=1.0')
pl.errorbar(bins, pb.make_pred([0.0,0.0,0.0,0.0,0.0,1.0,0.0]), xerr=0.0, yerr=0.0, label='ctphiQ1=1.0')
pl.errorbar(bins, pb.make_pred([0.0,0.0,0.0,0.0,0.0,0.0,1.0]), xerr=0.0, yerr=0.0, label='ctphiQ3=1.0')
pl.errorbar(bins, preds['SM'], xerr=0.0, yerr=0.0, label='SM')
pl.errorbar(bins, data, fmt="o",xerr=0.25, yerr=0.05, label='Data')

pl.axis([-1.0, 1.0, -5.0, 10.0])
pl.legend()
pl.show()

fig = corner.corner(samples, labels=["$ctg$", "$ctw$", "$ctphi$", "$ctb$", "$cphit$","$ctphiQ1$","$ctphiQ3$"],
                    truths=[0.0, 0.0, 0.0,0.0, 0.0, 0.0,0.0])

fig.show()
fig.savefig("triangle.png")


#end time
end = time.time()
print "Total elapsed wall time  = " +  str(end - start)
