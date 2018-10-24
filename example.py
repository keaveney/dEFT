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


#intro graphics
Screen.wrapper(splash.deft_splash)

#start time
start = time.time()

######################################################
###############   DEFINE THE DATA   ##################
######################################################
#bins
bins = np.array([-1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5])
nbins = len(bins)

#data
data = (4.0)*(np.array([7.0, 7.0, 3.0, 3.0, 3.0, 7.0, 7.0]))

#covariance matrix and inverse
cov = [[1.0, 0.7, 0.0, 0.0, 0.0, 0.0, 0.0], [0.7, 1.0, 0.5, 0.0, 0.0, 0.0, 0.0], [0.0, 0.5, 1.0, 0.5, 0.0, 0.0, 0.0], [0.0, 0.0, 0.5, 1.0, 0.5, 0.0, 0.0],[0.0, 0.0, 0.0, 0.5, 1.0, 0.5, 0.0],[0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.5], [0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0]]
icov = np.linalg.inv(cov)

######################################################
###############   DEFINE THE MODEL   #################
######################################################

#list of coeffieicnts/operators to be considered
coefficients = ["ctg", "ctw", "ctphi", "ctb", "cphit", "cphiQ1", "cphiQ3"]
ndim = len(coefficients)


preds =  {
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


#just generate some dummy cross predictions
for c in coefficients:
    ci_m_name = c + "-"
    ci_p_name = c + "+"
    for d in coefficients:
        if c != d:
            cj_m_name = c + "-"
            cj_p_name = c + "+"
            ci_cj_m_name = c + "-" + d + "-"
            ci_cj_p_name = c + "+" + d + "+"
            #wild guess
            preds[ci_cj_m_name] = preds[ci_m_name] + preds[cj_m_name] - (0.8*preds['SM'])
            preds[ci_cj_p_name] = preds[ci_p_name] + preds[cj_p_name] - (0.8*preds['SM'])


priors =  {
        'ctg': np.array([-3.0,3.0]),
        'ctw': np.array([-3.0,3.0]),
        'ctphi': np.array([-3.0,3.0]),
        'ctb': np.array([-3.0,3.0]),
        'cphit': np.array([-3.0,3.0]),
        'cphiQ1': np.array([-3.0,3.0]),
        'cphiQ3': np.array([-3.0,3.0])
}

#maximum power of c_i/Lambda^2 to be considered
# this can only be 1 or 2
# 1: Lambda^{-2} - interference between diagrams with single dim6 operator and SM diagrams
# 2: Lambda^{-2} + Lambda^{-4} - the above with the 'squared' terms (diagrams with single dim6 operator squared)
#                                  and 'cross' terms (interefernece between pairs of diagrams with single dim6 operator each)
max_coeff_power = 2

#the tricky part - list of predictions for data and for +/- X for each coefficient
# if max_coeff_power = 1, and more than one operator is being considred,
# this dictonary also needs predictions for (ci,cj) = (X,X)
# for all pairs of operators

# for now assume all input predictions correspond to c_i = +/- 2.0
c_i_benchmark = 2.0


def make_basis(preds,coefficients,max_coeff_power):
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

pred_basis = make_basis(preds,coefficients,max_coeff_power)

def make_pred(c):
    pred = np.zeros(nbins)
    for ci in range(0, len(coefficients)):
        basis_ci_sm_name = coefficients[ci] + "_sm"
        basis_ci_ci_name = coefficients[ci] + "_" + coefficients[ci]

        ci_sm_contrib = (c[ci]*pred_basis[basis_ci_sm_name])
        ci_ci_contrib = ((c[ci]**2)*pred_basis[basis_ci_ci_name])

        pred += ci_sm_contrib
        pred += ci_ci_contrib
                       
    pred = pred + preds['SM']
    return pred

def lnprior(c):
    #one could play with implementing non-flat priors here
    lnp = 0.0
    for c in coefficients:
        if ((c < priors[c][0]) | (c > priors[c][1])):
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
for c in range(0,len(coefficients)):
    vals = [0.0,0.0,0.0,0.0,0.0,0.0,0.0]
    vals[c] = 1.0
    label_string = coefficients[c] + "=1.0"
    pl.errorbar(bins, make_pred(vals), xerr=0.0, yerr=0.0, label=label_string)
pl.errorbar(bins, preds['SM'], xerr=0.0, yerr=0.0, label='SM')
pl.errorbar(bins, data, fmt="o",xerr=0.25, yerr=0.05, label='Data')
pl.axis([bins[0]-0.25, bins[nbins-1]+0.25, -5.0, 10.0])
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
