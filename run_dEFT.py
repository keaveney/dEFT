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

start = time.time()

######################################################
###############   READ CONFIG   ######################
######################################################
filename = sys.argv[1]
config = configReader()
config.init(filename)
pb = predBuilder()
pb.init(config.predictions,config.coefficients,config.params["config"]["model"]["max_coeff_power"],config.params["config"]["model"]["c_i_benchmark"],config.cross_terms,config.params["config"]["model"]["inclusive_k_factor"])


######################################################
###############   DEFINE LIKELIHOOD   ################
######################################################
def lnprior(c):
    lnp = 0.0
    for scan_ind in c:
        for coeff in config.coefficients:
            if ((scan_ind < config.prior_limits[coeff][0]) | (scan_ind > config.prior_limits[coeff][1])):
                lnp = -np.inf
    return lnp
def lnprob(c, data, icov):
    pred = pb.make_pred(c)
    diff = pred - data
    #ll =  (-np.dot(diff,np.dot(icov,diff))/2.0) + (lnprior(c))
    ll =  (-np.dot(diff,np.dot(icov,diff)))
    return ll


######################################################
###############   RUN THE FIT   ######################
######################################################
nWalkers = config.n_walkers
ndim = int((len(config.predictions) - 1.0)/2.0)
nBurnIn = config.n_burnin
nTotal = config.n_total
p0 = np.random.rand(ndim * nWalkers).reshape(nWalkers, ndim)
#p0 = [np.zeros(ndim) + 1e-4*np.random.randn(ndim) for i in range(nWalkers)]
#p0 = np.zeros(ndim * nWalkers).reshape(nWalkers, ndim)
sampler = emcee.EnsembleSampler(nWalkers, ndim, lnprob, args=[config.params["config"]["data"]["central_values"], config.icov])
pos, prob, state = sampler.run_mcmc(p0, nBurnIn)
#print "pos " + str(pos)
sampler.reset()
sampler.run_mcmc(pos,nTotal)
samples = sampler.chain.reshape((-1, ndim))
if isinstance(samples, np.ndarray):
    print "is numpy array"

print "summaryPlotter 6 = "
print str(samples)

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
