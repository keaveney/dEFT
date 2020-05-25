from tools.predBuilder import predBuilder
from tools.summaryPlotter import summaryPlotter
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
if (len(sys.argv) > 2):
    filenameTest = sys.argv[2]
    configTest = configReader()
    configTest.init(filenameTest)
    
    testSamples = configTest.params["config"]["model"]["samples"]
    testPreds = np.asarray(configTest.params["config"]["model"]["predictions"])

    print("##############################################################")
    print("###############  VALIDATION OF MORPHING MODEL  ###############")
    print("##############################################################")
    
    ndDistances = []
    morphPreds = ([])
    diffs = []
    testPredStrs =  ([])
    iDistances = np.zeros((len(testPreds),len(config.params["config"]["model"]["prior_limits"]) ))
    
    if (len(config.params["config"]["model"]["prior_limits"]) == 2): # a special validation of diffs vs (ci, cj) can be visualised if nOPs=2
        ciP = ([])
        cjP = ([])
        absDiff = ([])
    
    for p in range(0, len(testPreds)):
        print("coefficient values at test point =        " + str(testSamples[p]))
        print("prediction at test point =                " + str(testPreds[p]))
        testPred = testPreds[p]
        ciPoint = np.delete(testSamples[p], 0)
        morphPred = pb.makePred(ciPoint)
        morphPreds = np.append(morphPreds,morphPred)
        print("morphing model prediction at test point = " +str(morphPred))
        testPredStr = str(ciPoint)
        testPredStrs = np.append(testPredStrs, testPredStr)
        
        if (len(ciPoint) == 2):
            ciP = np.append(ciP,ciPoint[0])
            cjP = np.append(cjP,ciPoint[1])
            absDiff = np.append(absDiff, (morphPred - testPred))

        smPoint = np.zeros(len(config.params["config"]["model"]["prior_limits"]))
        # find Euclidean distance between test pred and SM in ci space
        ndDistance = np.linalg.norm(ciPoint-smPoint)
        print("ND distance (SM -> test point)          = " +str(ndDistance))
        print(" ")
        ndDistances.append(ndDistance)
        diff = 0.0
        
        for c in range(0, len(ciPoint)):
            iDistances[p][c] = ciPoint[c]
        
        # find sum of squared differences between morph pred and test pred
        for bin in range(0, len(morphPred)):
            diff = diff + ((morphPred[bin] - testPred[bin])**2)
            
        diff = np.abs(morphPred)
        diffs.append(diff)
        
        pl.figure()
        pl.errorbar(config.x_vals, testPred, xerr=0.0, yerr=0.0, label="benchmark prediction", fmt="o")
        pl.errorbar(config.x_vals, morphPred, xerr=0.0, yerr=0.0, label="morphing model", fmt="r+")
        #ax = pl.gca()
        #labely = ax.set_xlabel(xlabel, fontsize = 18)
        #ax.xaxis.set_label_coords(0.85, -0.065)
        #labely = ax.set_ylabel(ylabel, fontsize = 18)
        #ax.yaxis.set_label_coords(-0.037, 0.83)
        pl.legend(loc=2)
        plotfilename = str(config.params["config"]["run_name"]) + "_results/" + str(config.params["config"]["run_name"]) + "_MorphingValidation_" + str(p)  +".png"
        pl.savefig(plotfilename)
        pl.close()
        
    ax = pl.axes(projection='3d')
    ax.scatter3D(ciP, cjP, absDiff, cmap='Greens')
    pl.savefig("diff.png")

    pl.figure()
    
    itr = np.arange(len(morphPreds))
    
    #fig = pl.errorbar(x, bestFits, yerr=[marginUncsDown, marginUncsUp], fmt='o')
    fig = pl.errorbar(itr, morphPreds, fmt='o', label="morph model")
    fig = pl.errorbar(itr, testPreds, fmt='o', label="test predictions")
    pl.gcf().subplots_adjust(bottom=0.15)
    pl.legend(loc=2)
    pl.xticks(range(len(testPredStrs)), testPredStrs, rotation='45')
    pl.savefig(config.params["config"]["run_name"] + "_results/" + config.params["config"]["run_name"] + "_MorphingValidation_summary.png")
    pl.close()

    for c in range(0,len(ciPoint)):
        ciVals = iDistances[:, c]

        pl.figure()
        pl.scatter(ciVals, diffs)
        pl.legend(loc=2)
        plotfilename = str(config.params["config"]["run_name"]) + "_results/" + str(config.params["config"]["run_name"]) + "_SquaredDiff_vs_" + str(c) + "_iDistance.png"
        pl.savefig(plotfilename)
        pl.close()

    pl.figure()
    pl.scatter(ndDistances,diffs)
    pl.legend(loc=2)
    plotfilename = str(config.params["config"]["run_name"]) + "_results/" + str(config.params["config"]["run_name"]) + "_SquaredDiff_vs_nDDistance.png"
    pl.savefig(plotfilename)
    pl.close()

    pl.figure()
    #pl.errorbar(config.x_vals, pb.makePred([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]),xerr=0.0, yerr=0.0, label="morphing model")
    pl.errorbar(config.x_vals, testPreds[0], xerr=0.0, yerr=0.0, label="benchmark prediction", fmt="o")
    pl.legend(loc=2)
    plotfilename = str(config.params["config"]["run_name"]) + "_results/" + str(config.params["config"]["run_name"]) + "_MorphingValidation_SM.png"
    pl.savefig(plotfilename)
    pl.close()

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
    #p0 = np.zeros(ndim * nWalkers).reshape(nWalkers, ndim)

    ###############   FIND ML RESULT ######################
    nll = lambda *args: lnprob(*args)
    initial = p0
    #soln = minimize(nll, initial, args=(config.data, config.icov), method='SLSQP', options={'maxiter': 100, 'ftol': 1e-01, 'iprint': 1, 'disp': False, 'eps': 1.4901161193847656e-08})
    #print("ML soln = " + str(soln.x))

    with Pool(processes=6) as pool:
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
