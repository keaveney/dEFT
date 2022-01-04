from tools.predBuilder import predBuilder
from tools.summaryPlotter import summaryPlotter
from tools.modelValidator import modelValidator

import numpy as np
import emcee
import time
import json
import sys
import math
from tools.configReader import configReader
from multiprocessing import Pool
import matplotlib.pyplot as pl

np.random.seed(42) # hard code random seed for reproducibility, could be configurable


start = time.time()

###################################################
####### READ CONFIG, BUILD MORPHING MODEL #########
###################################################
filename = sys.argv[1]
config = configReader()
config.init(filename)
pb = predBuilder()

if(config.params["config"]["model"]["input"] == "numpy"):
    pb.initRM(len(config.prior_limits), config.samples, config.predictions, config.kfac)
    print("SM pred == " +  str(pb.makeRMPred(np.zeros(len(config.prior_limits)))) + " ")


##########################################################
########  VALIDATE OF MORPHING MODEL (OPTIONAL)  #########
##########################################################

if (len(sys.argv) > 2):
    filenameTest = sys.argv[2]
    configTest = configReader()
    configTest.init(filenameTest)
    
    mv = modelValidator()
    mv.validate(configTest, pb)
    
#######################################
######### ESTIMATE POSTERIOR ##########
#######################################

def lnprior(c):
    lnp = 0.0
    for scan_ind in range(0, len(c)):
        if ((c[scan_ind] < config.prior_limits[config.coefficients[scan_ind]][0]) | (c[scan_ind] > config.prior_limits[config.coefficients[scan_ind]][1])):
            lnp = -np.inf
    return lnp
    
def lnprob(c, data, icov):
    pred =  pb.makeRMPred(c)
    diff = pred - data
    ll =  (-0.5 * np.dot(diff,np.dot(icov,diff))) + (lnprior(c))

    return ll

if (len(sys.argv) <= 2):
    nWalkers = config.n_walkers
    ndim = int(len(config.params["config"]["model"]["prior_limits"]))
    nBurnIn = config.n_burnin
    nTotal = config.n_total
    p0 = [np.zeros(ndim) + 1e-4*np.random.randn(ndim) for i in range(nWalkers)]
    
    """
    multiprocessing is broken in python > 3.6, needs investiagtion, reverting to serial for now
    with MPIPool() as pool:
    with Pool(processes=4) as pool:
    with Pool() as pool:
        if not pool.is_master():
            pool.wait()
            sys.exit(0)
    sampler = emcee.EnsembleSampler(nWalkers, ndim, lnprob,  pool=pool, args=[config.params["config"]["data"]["central_values"], config.icov])
    """
    
    np.random.seed(42) # hard code random seed for reproducibility, could be configurable
    sampler = emcee.EnsembleSampler(nWalkers, ndim, lnprob, args=[config.params["config"]["data"]["central_values"], config.icov])
    pos, prob, state = sampler.run_mcmc(p0, nBurnIn, progress=True)
    sampler.reset()
    sampler.run_mcmc(pos,nTotal, progress=True)
    samples = sampler.chain.reshape((-1, ndim))

    sp = summaryPlotter()
    sp.summarise(config, pb, sampler, samples)
    end = time.time()
    print("Total elapsed wall time  = " +  str(int(end - start)) + " seconds")
