from tools.predBuilder import predBuilder
from tools.summaryPlotter import summaryPlotter
import numpy as np
import emcee
import time
import json
import sys
import matplotlib.pyplot as pl
import pandas as pd
from tools.configReader import configReader

from multiprocessing import Pool

start = time.time()

######################################################
###############   READ CONFIG   ######################
######################################################
filename = sys.argv[1]
config = configReader()
config.init(filename)
pb = predBuilder()

#could pb.init just be called with 'config' as argument?
#pb.init(config.predictions,config.coefficients,config.params["config"]["model"]["max_coeff_power"],config.params["config"]["model"]["c_i_benchmark"],config.cross_terms,config.params["config"]["model"]["inclusive_k_factor"])

#self, nOps, inputCoeffValues, preds

if(config.params["config"]["model"]["input"] == "numpy"):
    pb.init(len(config.params["config"]["model"]["prior_limits"]), config.params["config"]["model"]["samples"], np.asarray(config.params["config"]["model"]["predictions"]))


#print("testing pb init = " + str(config.params["config"]["model"]["samples"]) + " --- " + str(np.asarray(config.params["config"]["model"]["predictions"])))

#pb.init(len(config.params["config"]["model"]["prior_limits"]), config.params["config"]["model"]["samples"], config.params["config"]["model"]["predictions"])

######################################################
###############   DEFINE LIKELIHOOD   ################
######################################################
def lnprior(c):
    lnp = 0.0
    #print("ln prior: " + str(len(c)) + "  " + str(config.prior_limits))
    #print("ln test: " + str(config.coefficients[0]))
    for scan_ind in range(0, len(c)):
        if ((c[scan_ind] < config.prior_limits[config.coefficients[scan_ind]][0]) | (c[scan_ind] > config.prior_limits[config.coefficients[scan_ind]][1])):
            lnp = -np.inf
    return lnp

def lnprob(c, data, icov):
    #print("lnprob::icov = " + str(icov))
    pred = pb.makePred(c)
    diff = pred - data
    #ll =  (-np.dot(diff,diff)) + (lnprior(c))
    ll =  (-np.dot(diff,np.dot(icov,diff))/2.0) + (lnprior(c))
    #ll =  (-np.dot(diff,np.dot(icov,diff)))
    return ll

##############################################################
###############  VALIDATION OF MORPHING MODEL  ###############
##############################################################
"""val = np.array([-2.19, 11.4, -9.99, -0.1, 6.5, -13.34])
print("val prediction for [-2.19, 11.4, -9.99, -0.1, 6.5, -13.34] = " + str(pb.makePred(val)))
print("MadGraph prediction for [-2.19, 11.4, -9.99, -0.1, 6.5, -13.34] = 1.027")


val = np.array([0.5, 5.9, 11.0, 18.0, 5.3, -0.6])
print("val prediction for [0.5, 5.9, 11.0, 18.0, 5.3, -0.6]= " + str(pb.makePred(val)))
print("MadGraph prediction for [0.5, 5.9, 11.0, 18.0, 5.3, -0.6] = 1.089")

val = np.array([-9.4, -3.0, 4.0, 61.0, -23.0, 111.0])
print("val prediction for [-9.4, -3.0, 4.0, 61.0, -23.0, 111.0] = " + str(pb.makePred(val)))
print("MadGraph prediction for [-9.4, -3.0, 4.0, 61.0, -23.0, 111.0] = 12.8")

val = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
print("val prediction for [0.0, 0.0, 0.0, 0.0, 0.0, 0.0] = " + str(pb.makePred(val)))
print("MadGraph prediction for [0.0, 0.0, 0.0, 0.0, 0.0, 0.0] = 0.5612")
"""

######################################################
###############   RUN THE FIT   ######################
######################################################
nWalkers = config.n_walkers
#ndim = int((len(config.predictions) - 1.0)/2.0)
ndim = int(len(config.params["config"]["model"]["prior_limits"]))
nBurnIn = config.n_burnin
nTotal = config.n_total
p0 = np.random.rand(ndim * nWalkers).reshape(nWalkers, ndim)
#p0 = [np.zeros(ndim) + 1e-4*np.random.randn(ndim) for i in range(nWalkers)]
#p0 = np.zeros(ndim * nWalkers).reshape(nWalkers, ndim)

print("icov = " + str(config.icov))

with Pool() as pool:
    sampler = emcee.EnsembleSampler(nWalkers, ndim, lnprob,  pool=pool, args=[config.params["config"]["data"]["central_values"], config.icov])
    pos, prob, state = sampler.run_mcmc(p0, nBurnIn,  progress=True)
    #print "pos " + str(pos)
    sampler.reset()
    sampler.run_mcmc(pos,nTotal,  progress=True)

samples = sampler.chain.reshape((-1, ndim))
if isinstance(samples, np.ndarray):
    print("is numpy array")

mcmc_params = np.mean(sampler.flatchain,axis=0)
mcmc_params_cov = np.cov(np.transpose(sampler.flatchain))

######################################################
###############   PLOT RESULTS   ####################
######################################################
#print "best fit vals = " + str(mcmc_params)
#print "uncertainties = " + str(np.sqrt(np.diag(mcmc_params_cov)))
#print "Summarising: pb =  "  +  str(ndim) + " samples=  " + str(samples)

sp = summaryPlotter()
sp.summarise(config, pb, sampler, samples, mcmc_params)
end = time.time()
#print "Total elapsed wall time  = " +  str(end - start)
