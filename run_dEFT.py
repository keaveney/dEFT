from tools import splash
from tools.predBuilder import predBuilder
from tools.summaryPlotter import summaryPlotter
from asciimatics.screen import Screen
import numpy as np
import emcee
import time
import json
import sys
import matplotlib.pyplot as pl
from tools.configReader import configReader

#intro graphics
Screen.wrapper(splash.deft_splash)

#start time
start = time.time()
filename = sys.argv[1]
config = configReader()
config.init(filename)

pb = predBuilder()
#pred_basis = predBuilder.make_basis(config.predictions,config.coefficients,config.params["config"]["model"]["max_coeff_power"],config.params["config"]["model"]["c_i_benchmark"])

pb.init(config.predictions,config.coefficients,config.params["config"]["model"]["max_coeff_power"],config.params["config"]["model"]["c_i_benchmark"],config.cross_terms)

print "pred basis = " + str(pb.pred_basis)

def lnprior(c):
    #one could play with implementing non-flat priors here
    # method here so far seems inefficient and/or wrong
    lnp = 0.0
    for scan_ind in c:
        for coeff in config.coefficients:
            if ((scan_ind < config.prior_limits[coeff][0]) | (scan_ind > config.prior_limits[coeff][1])):
                lnp = -np.inf
    #print "scan ind = " + str(scan_ind) + " lnp = " + str(lnp)
    return lnp

def lnprob(c, data, icov):
    pred = pb.make_pred(c)
    diff = pred - data
    ll =  (-np.dot(diff,np.dot(icov,diff))/2.0) + (lnprior(c))
    return ll

######################################################
###############   DEFINE THE FIT/SCAN   ##############
######################################################
nWalkers = config.n_walkers
ndim = len(config.coefficients)
nBurnIn = config.n_burnin
nTotal = config.n_total

p0 = np.random.rand(ndim * nWalkers).reshape(nWalkers, ndim)

sampler = emcee.EnsembleSampler(nWalkers, ndim, lnprob, args=[config.params["config"]["data"]["central_values"], config.icov])

pos, prob, state = sampler.run_mcmc(p0, nBurnIn)

sampler.reset()

sampler.run_mcmc(pos,nTotal)

samples = sampler.chain.reshape((-1, ndim))

mcmc_params = np.mean(sampler.flatchain,axis=0)
mcmc_params_cov = np.cov(np.transpose(sampler.flatchain))

print "best fit vals = " + str(mcmc_params)
#print "uncertainties = " + str(np.sqrt(np.diag(mcmc_params_cov)))

sp = summaryPlotter()

sp.summarise(config, pb, samples)


#end time
end = time.time()
print "Total elapsed wall time  = " +  str(end - start)
