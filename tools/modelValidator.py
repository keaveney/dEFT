import numpy as np
import matplotlib.pyplot as pl


class modelValidator:

    def validate(self, configTest, pb):

        testSamples = configTest.params["config"]["model"]["samples"]
        testPreds = np.asarray(configTest.params["config"]["model"]["predictions"])
        
        print("##############################################################")
        print("###############  VALIDATION OF MORPHING MODEL  ###############")
        print("##############################################################")
        
        ndDistances = []
        morphPreds = ([])
        diffs = []
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
            ciPoint = np.delete(testSamples[p], 0)
            morphPred = pb.makeRMPred(ciPoint)
            morphPreds = np.append(morphPreds,morphPred)
            print("morphing model prediction at test point = " +str(morphPred))
            testPredStr = str(ciPoint)
            testPredStrs = np.append(testPredStrs, testPredStr)
            
            if (len(ciPoint) == 2):
                ciP = np.append(ciP,ciPoint[0])
                cjP = np.append(cjP,ciPoint[1])
                absDiff = np.append(absDiff, (morphPred - testPred))

            smPoint = np.zeros(len(configTest.params["config"]["model"]["prior_limits"]))
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
            pl.errorbar(configTest.x_vals, testPred, xerr=0.0, yerr=0.0, label="benchmark prediction", fmt="o")
            pl.errorbar(configTest.x_vals, morphPred, xerr=0.0, yerr=0.0, label="morphing model", fmt="r+")
            #ax = pl.gca()
            #labely = ax.set_xlabel(xlabel, fontsize = 18)
            #ax.xaxis.set_label_coords(0.85, -0.065)
            #labely = ax.set_ylabel(ylabel, fontsize = 18)
            #ax.yaxis.set_label_coords(-0.037, 0.83)
            pl.legend(loc=2)
            plotfilename = str(configTest.params["config"]["run_name"]) + "_results/" + str(configTest.params["config"]["run_name"]) + "_MorphingValidation_" + str(p)  +".png"
            pl.savefig(plotfilename)
            pl.close()
            
        #ax = pl.axes(projection='3d')
        #ax.scatter3D(ciP, cjP, absDiff, cmap='Greens')
        #pl.savefig("diff.png")

        pl.figure()
        
        itr = np.arange(len(morphPreds))
        
        print ("morph preds = " + str(len(morphPreds)))
        print ("test preds = " + str(len(testPreds)))

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

        for c in range(0,len(ciPoint)):
            ciVals = iDistances[:, c]

            pl.figure()
            #pl.scatter(ciVals, diffs)
            pl.legend(loc=2)
            plotfilename = str(configTest.params["config"]["run_name"]) + "_results/" + str(configTest.params["config"]["run_name"]) + "_SquaredDiff_vs_" + str(c) + "_iDistance.png"
            pl.savefig(plotfilename)
            pl.close()

        pl.figure()
        #pl.scatter(ndDistances,diffs)
        pl.legend(loc=2)
        plotfilename = str(configTest.params["config"]["run_name"]) + "_results/" + str(configTest.params["config"]["run_name"]) + "_SquaredDiff_vs_nDDistance.png"
        pl.savefig(plotfilename)
        pl.close()

        pl.figure()
        #pl.errorbar(configTest.x_vals, pb.makeRMPred([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]),xerr=0.0, yerr=0.0, label="morphing model")
        pl.errorbar(configTest.x_vals, testPreds[0], xerr=0.0, yerr=0.0, label="benchmark prediction", fmt="o")
        pl.legend(loc=2)
        plotfilename = str(configTest.params["config"]["run_name"]) + "_results/" + str(configTest.params["config"]["run_name"]) + "_MorphingValidation_SM.png"
        pl.savefig(plotfilename)
        pl.close()
