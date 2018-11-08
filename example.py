from tools import splash
from tools import pred_builder
from asciimatics.screen import Screen
import numpy as np
import itertools
import matplotlib.pyplot as pl
import corner
import emcee
import time
import json
from config import *
import json
import sys
from configReader import configReader

#intro graphics
Screen.wrapper(splash.deft_splash)

#start time
start = time.time()


filename = sys.argv[1]

config = configReader()

config.init(filename)

print config.params["config"]["data"]["bins"]


def make_basis(preds,coefficients,max_coeff_power,c_i_benchmark):
    pred_basis = {}
    for c in coefficients:
        ci_m_name = c + "-"
        ci_p_name = c + "+"
        # sm--ci interference contribution
        sigma_sm_ci = (np.subtract(preds[ci_p_name], preds[ci_m_name]) / (2*c_i_benchmark))
        # ci--ci squared contribution
        if max_coeff_power > 1.0:
            sigma_ci_ci = (np.add(preds[ci_p_name], preds[ci_m_name]) - (2.0*preds['SM']) ) / (2*(c_i_benchmark**2))
            for d in coefficients:
                # TODO ci--cj interference contribution
                if c != d:
                    ci_cj_m_name = c + "-" + d + "-"
                    ci_cj_p_name = c + "+" + d + "+"
                    if max_coeff_power > 1.0:
                        sigma_cj_cj = (np.add(preds[ci_cj_m_name], preds[ci_cj_p_name]) - (2.0*preds['SM']) ) / (2*(c_i_benchmark**2))
                        sigma_ci_cj = (np.add(preds[ci_cj_m_name], preds[ci_cj_p_name]) - (2.0*preds['SM'])  - (2.0*(sigma_ci_ci**2)) - (2.0*(sigma_cj_cj**2)) ) / (2*(c_i_benchmark**2))
                    else:
                        sigma_ci_ci = np.zeros(len(preds['SM']))
                    basis_ci_cj_name = c + "_" + d
                    pred_basis[basis_ci_cj_name] = sigma_ci_cj
        else:
            sigma_ci_ci = np.zeros(len(preds['SM']))
        
        basis_ci_sm_name = c + "_sm"
        basis_ci_ci_name = c + "_" + c

        pred_basis[basis_ci_sm_name] = sigma_sm_ci
        pred_basis[basis_ci_ci_name] = sigma_ci_ci

    return pred_basis

pred_basis = make_basis(config.predictions,config.coefficients,config.params["config"]["model"]["max_coeff_power"],config.params["config"]["model"]["c_i_benchmark"] )

def make_pred(c):
    pred = np.zeros(nbins)
    for ci in range(0, len(config.coefficients)):
        basis_ci_sm_name = config.coefficients[ci] + "_sm"
        basis_ci_ci_name = config.coefficients[ci] + "_" + config.coefficients[ci]

        ci_sm_contrib = (c[ci]*pred_basis[basis_ci_sm_name])
        ci_ci_contrib = ((c[ci]**2)*pred_basis[basis_ci_ci_name])

        pred += ci_sm_contrib
    #pred += ci_ci_contrib
                       
    pred = pred + config.predictions['SM']
    return pred

def lnprior(c):
    #one could play with implementing non-flat priors here
    lnp = 0.0
    for c in config.coefficients:
        if ((c < config.prior_limits[c][0]) | (c > config.prior_limits[c][1])):
            prob = -np.inf
    return lnp

def lnprob(c, data, icov):
    pred = make_pred(c)
    diff = pred - data
    return (-np.dot(diff,np.dot(icov,diff))/2.0) + lnprior(c)

######################################################
###############   DEFINE THE FIT/SCAN   ##############
######################################################
nWalkers = 100
ndim = 7
nBurnIn = 50
nTotal = 300

p0 = np.random.rand(ndim * nWalkers).reshape((nWalkers, ndim))
sampler = emcee.EnsembleSampler(nWalkers, ndim, lnprob, args=[data, icov])
pos, prob, state = sampler.run_mcmc(p0, nBurnIn)
sampler.reset()
sampler.run_mcmc(pos,nTotal)
samples = sampler.chain.reshape((-1, ndim))

mcmc_params = np.mean(sampler.flatchain,axis=0)
mcmc_params_cov = np.cov(np.transpose(sampler.flatchain))

print "best fit vals = " + str(mcmc_params)
print "uncertainties = " + str(np.sqrt(np.diag(mcmc_params_cov)))

######################################################
###############   MAKE RESULTS PLOTS  ################
######################################################
# raw predictions and data
pl.figure()
for c in range(0,len(config.coefficients)):
    vals = [0.0,0.0,0.0,0.0,0.0,0.0,0.0]
    vals[c] = 1.0
    label_string = config.coefficients[c] + "=1.0"
    pl.errorbar(bins, make_pred(vals), xerr=0.0, yerr=0.0, label=label_string)
pl.errorbar(config.params["config"]["data"]["bins"], config.predictions['SM'], xerr=0.0, yerr=0.0, label='SM')
pl.errorbar(config.params["config"]["data"]["bins"], config.params["config"]["data"]["central_values"], fmt="o",xerr=0.25, yerr=0.05, label='Data')
pl.axis([config.params["config"]["data"]["bins"][0]-0.25, config.params["config"]["data"]["bins"][nbins-1]+0.25, -5.0, 10.0])
pl.legend()
#pl.show()

#best fit prediction and data
#TODO


#corner plot
#fig = corner.corner(samples, labels=["$ctg$", "$ctw$", "$ctphi$", "$ctb$", "$cphit$","$ctphiQ1$","$ctphiQ3$"],
#                    truths=[0.0, 0.0, 0.0,0.0, 0.0, 0.0,0.0])

fig = corner.corner(samples, labels=["$ctg$", "$ctw$", "$ctphi$", "$ctb$", "$cphit$", "$ctphiQ1$", "$ctphiQ3$"],
                    quantiles=[0.05, 0.95],
                       range=[1.0,1.0,1.0,1.0,1.0,1.0,1.0], truths=[0.0, 0.0, 0.0,0.0, 0.0, 0.0,0.0],
                       show_titles=True, title_kwargs={"fontsize": 18})

fig.savefig("triangle.png")


#end time
end = time.time()
print "Total elapsed wall time  = " +  str(end - start)
