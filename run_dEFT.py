from tools.predBuilder import predBuilder
from tools.summaryPlotter import summaryPlotter
from tools.modelValidator import modelValidator

from schwimmbad import MPIPool

from scipy.optimize import minimize
import numpy as np
import emcee
import time
import json
import sys
import matplotlib.pyplot as pl
import pandas as pd
from tools.configReader import configReader
from multiprocessing import Pool
import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d import Axes3D

start = time.time()
######################################################
###############   READ CONFIG   ######################
######################################################
filename = sys.argv[1]
config = configReader()
config.init(filename)
pb = predBuilder()

if(config.params["config"]["model"]["input"] == "numpy"):
    pb.initRM(len(config.params["config"]["model"]["prior_limits"]), np.asarray(config.params["config"]["model"]["samples"]), np.asarray(config.params["config"]["model"]["predictions"]))


##############################################################
###############  VALIDATION OF MORPHING MODEL  ###############
##############################################################

if (len(sys.argv) > 2):

    filenameTest = sys.argv[2]
    configTest = configReader()
    configTest.init(filenameTest)
    
    mv = modelValidator()
    mv.validate(configTest, pb)
    

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
    #print("lnprob::c = " + str(c))
    #print("lnprob::data = " + str(data))
    #print("lnprob::icov = " + str(icov))
    #pred = pb.makePred(c)
    pred =  pb.makeRMPred(c)
    #print("lnprob::pred = " + str(pred))

    diff = pred - data
    #ll =  (-np.dot(diff,diff)) + (lnprior(c))
    #ll =  (-np.dot(diff,np.dot(icov,diff))/2.0) + (lnprior(c))
    ll =  (-np.dot(diff,np.dot(icov,diff))/2.0) 
    return ll

######################################################
###############   RUN THE FIT   ######################
######################################################
if (len(sys.argv) <= 2):
    nWalkers = config.n_walkers
    #ndim = int((len(config.predictions) - 1.0)/2.0)
    ndim = int(len(config.params["config"]["model"]["prior_limits"]))
    nBurnIn = config.n_burnin
    nTotal = config.n_total
    #p0 = np.random.rand(ndim * nWalkers).reshape(nWalkers, ndim)
    p0 = [np.zeros(ndim) + 1e-4*np.random.randn(ndim) for i in range(nWalkers)]
    #print("p0 = " + str(p0))
    #p0 = np.zeros(ndim * nWalkers).reshape(nWalkers, ndim)

    print("RUNNING THE FIT 1")

    with Pool(processes=40) as pool:
#    with MPIPool() as pool:
#        if not pool.is_master():
#            pool.wait()
#            sys.exit(0)
        sampler = emcee.EnsembleSampler(nWalkers, ndim, lnprob, pool=pool, args=[config.params["config"]["data"]["central_values"], config.icov])

#    sampler = emcee.EnsembleSampler(nWalkers, ndim, lnprob,  pool=pool, args=[config.params["config"]["data"]["central_values"], config.icov])
        pos, prob, state = sampler.run_mcmc(p0, nBurnIn)
    #print "pos " + str(pos)
        sampler.reset()
        sampler.run_mcmc(pos,nTotal)

    samples = sampler.chain.reshape((-1, ndim))
    if isinstance(samples, np.ndarray):
        print("is numpy array")

    sp = summaryPlotter()
    sp.summarise(config, pb, sampler, samples)
    end = time.time()
    #print "Total elapsed wall time  = " +  str(end - start)
