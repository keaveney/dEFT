import numpy as np
import matplotlib.pyplot as pl
from scipy.stats import chi2, norm
import scipy


class modelValidator:

    def chiSq(self, testPred, testUncs, morphPred):
        testPredArr = np.array(testPred)
        testUncsArr = np.array(testUncs)
        morphPredArr = np.array(morphPred)
        
        return np.sum(np.square(testPredArr - morphPredArr)/(np.square(testUncsArr)))
        
    def validate(self, configTest, pb):

        nbins =  len(np.asarray(configTest.params["config"]["model"]["predictions"][0]))
        testSamples = configTest.params["config"]["model"]["samples"]
        testPreds = np.asarray(configTest.params["config"]["model"]["predictions"])
        #uncertainties = np.asarray(configTest.params["config"]["model"]["uncertainties"])
        
        print("##############################################################")
        print("###############  VALIDATION OF MORPHING MODEL  ###############")
        print("##############################################################")
        
        ndDistances = []
        chiSqs = []
        morphPreds = ([])

        testPredStrs =  ([])
        iDistances = np.zeros((len(testPreds),len(configTest.params["config"]["model"]["prior_limits"]) ))
        
        if (len(configTest.params["config"]["model"]["prior_limits"]) == 2): # a special validation of diffs vs (ci, cj) can be visualised if nOPs=2
            ciP = ([])
            cjP = ([])
            absDiff = ([])
        
        for p in range(0, len(testPreds)):
            print("coefficient values at test point =        " + str(testSamples[p]))
            print("prediction at test point =                " + str(testPreds[p]))
            testPred = testPreds[p]
            #uncs = uncertainties[p]

            ciPoint = np.delete(testSamples[p], 0)
            morphPred = pb.makeRMPred(ciPoint)
            morphPreds = np.append(morphPreds,morphPred)

            print("morphing model prediction at test point = " +str(morphPred))
            
            #chisq = self.chiSq(testPred,uncs,morphPred)
            #chiSqs = np.append(chiSqs, chisq)
            #print("chi2 = " + str(chisq))
            
            testPredStr = str(ciPoint)
            testPredStrs = np.append(testPredStrs, testPredStr)
 
            fig = pl.figure()
            ax = fig.add_subplot(111)

            pl.errorbar(configTest.x_vals, testPred, xerr=0.0, label="benchmark prediction", fmt="o")
            pl.errorbar(configTest.x_vals, morphPred, xerr=0.0, yerr=0.0, label="morphing model", fmt="r+")

            pos = np.max(morphPred)/2.0

            #txtLabel0 = "$\chi^{2} = $" + str(round(chisq, 1))
            #ax.text(300, pos, txtLabel0, fontsize=15)
            
            #txtLabel1 = "$\chi^{2}/ndf = $" + str(round(chisq/float(nbins), 1))
            #ax.text(300, pos/1.5, txtLabel1, fontsize=15)
            
            #ax = pl.gca()
            #labely = ax.set_xlabel(xlabel, fontsize = 18)
            #ax.xaxis.set_label_coords(0.85, -0.065)
            #labely = ax.set_ylabel(ylabel, fontsize = 18)
            #ax.yaxis.set_label_coords(-0.037, 0.83)
            pl.legend(loc=2)
            plotfilename = str(configTest.params["config"]["run_name"]) + "_results/" + str(configTest.params["config"]["run_name"]) + "_MorphingValidation_" + str(p)  +".png"
            pl.savefig(plotfilename)
            pl.close()
            
        pl.figure()
        
        itr = np.arange(len(morphPreds))
        
        print ("morph preds = " + str(len(morphPreds)))
        print ("test preds = " + str(len(testPreds)))
        
        residuals = (testPreds.flatten()- morphPreds)/(testPreds.flatten())
        _, bins, _ = pl.hist(residuals, 15, density=1, alpha=0.5)

        mu, sigma = scipy.stats.norm.fit(residuals)
        
        print ("% residual = " + str(100.0*(residuals)))
        print ("mean/sigma residuala = " + str(mu) + ", " + str(sigma))

        #fig = pl.errorbar(x, bestFits, yerr=[marginUncsDown, marginUncsUp], fmt='o')
        #fig = pl.errorbar(itr, morphPreds, fmt='o', label="morph model")
        #fig = pl.errorbar(itr, testPreds, fmt='o', label="test predictions")
        pl.gcf().subplots_adjust(bottom=0.18)
        axes = pl.gca()
        ymax = np.amax(morphPreds)*1.3
        axes.set_ylim([0.0,ymax])
        pl.legend(loc=2)
        pl.xticks(range(len(testPredStrs)), testPredStrs, rotation='45')
        pl.savefig(configTest.params["config"]["run_name"] + "_results/" + configTest.params["config"]["run_name"] + "_MorphingValidation_summary.png")
        pl.close()

        pl.figure()
        pl.errorbar(configTest.x_vals, testPreds[0], xerr=0.0, yerr=0.0, label="benchmark prediction", fmt="o")
        pl.legend(loc=2)
        plotfilename = str(configTest.params["config"]["run_name"]) + "_results/" + str(configTest.params["config"]["run_name"]) + "_MorphingValidation_SM.png"
        pl.savefig(plotfilename)
        pl.close()
        
        pl.figure()
        fig, ax = pl.subplots(1, 1)
        #x = np.linspace(norm.ppf(0.01, 0.0, 0.01), norm.ppf(0.99, 0.0, 0.01), 100)
        #ax.plot(x, norm.pdf(x, 0.0, 0.01),  'r-', lw=5, alpha=0.6, label='norm pdf')
        x = np.linspace(norm.ppf(0.0001, mu, sigma), norm.ppf(0.9999, mu, sigma), 100)
        pl.hist(residuals, bins=bins, label="test points", density=True)
        best_fit_line = scipy.stats.norm.pdf(x, mu, sigma)
        pl.plot(x, best_fit_line)
        pl.legend(loc=2)
        plotfilename = str(configTest.params["config"]["run_name"]) + "_results/" + str(configTest.params["config"]["run_name"]) + "_residuals.png"
        pl.savefig(plotfilename)
        pl.close()
        
        pl.figure()
        fig, ax = pl.subplots(1, 1)
        #pl.hist(chiSqs, bins=12, label="test points", density=True)
        df = nbins
        x = np.linspace(chi2.ppf(0.01, df), chi2.ppf(0.99, df), 100)
        ax.plot(x, chi2.pdf(x, df),  'r-', lw=5, alpha=0.6, label='chi2 pdf')
        pl.legend(loc=2)
        plotfilename = str(configTest.params["config"]["run_name"]) + "_results/" + str(configTest.params["config"]["run_name"]) + "_testChiSquares.png"
        pl.savefig(plotfilename)
        pl.close()
        
        pl.figure()
        fig, ax = pl.subplots(1, 1)
        #pl.hist(chiSqs/nbins, bins=9, label="test points", density=True)
        pl.legend(loc=2)
        plotfilename = str(configTest.params["config"]["run_name"]) + "_results/" + str(configTest.params["config"]["run_name"]) + "_testChiSquaresNDF.png"
        pl.savefig(plotfilename)
        pl.close()
