import numpy as np
import itertools
import sklearn

# linear regression for multioutput
from sklearn.datasets import make_regression
from sklearn.linear_model import LinearRegression

from itertools import combinations

class predBuilder:

        def nSamples(self,n):
            return int((n+1)*(n+2)/2)
            
        def triangularNumber(self, n):
            tN = int(n*(n+1)/2)
            return tN
            
        def makeCoeff(self, ci):
            X = np.array([])
            for row in ci:
                X = np.append(X, row**2)
                combs = combinations(list(row), 2)
                num_rows = np.shape(ci)[0]
                for comb in combs:
                    X = np.append(X, comb[0]*comb[1])
    
            X = X.reshape(num_rows, self.triangularNumber(len(ci[0])))

            return X
            
        def makeCoeffPoint(self, ci):
            X = np.array([])
            X = np.append(X, ci**2)
            combs = combinations(list(ci), 2)
            for comb in combs:
                X = np.append(X, comb[0]*comb[1])
    
            X = X.reshape(1, self.triangularNumber(len(ci)))
        
            return X
            
        def initRM(self, nOps, samples, preds):
            # declare "model" object for scikit-learn (probably outside of init() method)
            # read in train points and predictions from json input
            # run 'fit()'
            # need method to convert ci point to 'couplings' vector?
            
            if (len(preds) <= self.nSamples(nOps)):
                print("ERROR: constraining " + str(nOps) + " coefficients requires at least " + str(self.nSamples(nOps)) + " samples,  but " + str(len(preds)) + " are provided")

            self.nOps = nOps
            cInputAr = np.array([])
            

            #convert to vector of "couplings" to make the morphing a linear system
            cInputAr = self.makeCoeff(samples)
            
            # define model
            model = LinearRegression()

            # fit model
            model.fit(cInputAr, preds)
            
            #print("n test points = " + str(cInputAr[0]) + "n predictions = " + str(len(preds[0])) )
            
            c = np.array([1.0, 24.920763888374044, 6.988302790442056, 0.24257571778552145, -24.381823700819854, 2.6365359630424203, 5.850613285030119, 7.7411147706635095])
            test = self.makeCoeffPoint(c)

            self.model = model
            
        def makeRMPred(self,c):

            testc = np.array([1.0, 24.920763888374044, 6.988302790442056, 0.24257571778552145, -24.381823700819854, 2.6365359630424203, 5.850613285030119, 7.7411147706635095])
            test = self.makeCoeffPoint(testc)

            c = np.append(1.0, c)
            cInputAr = self.makeCoeffPoint(c)
            pred = self.model.predict(cInputAr)
            
            return pred[0]
        
        def init(self, nOps, samples, preds):
        
            if (len(preds) <= self.nSamples(nOps)):
                print("ERROR: constraining " + str(nOps) + " coefficients requires at least " + str(self.nSamples(nOps)) + " samples,  but " + str(len(preds)) + " are provided")

            self.nOps = nOps
            cInputAr = np.array([])
            
            for sample in samples:
                couplings = np.array([])
                combs = list(itertools.combinations(sample,2))
                print("sample = " + str(sample))

                #SM term
                couplings = np.append(couplings, sample[0]**2)
                #'SM--Oi' interference terms
                for comb in range(0, nOps):
                    couplings = np.append(couplings, combs[comb][0]*combs[comb][1])
                #'squared' terms
                for ci in range(1, nOps+1):
                    couplings = np.append(couplings, sample[ci]**2)
                #'cross' terms
                for comb in range(nOps, len(combs)):
                    couplings = np.append(couplings, combs[comb][0]*combs[comb][1])
                    
                cInputAr = np.append(cInputAr, couplings , axis=0)

            cInputAr = np.reshape(cInputAr,(len(samples),len(samples)))
            cInputAr = np.transpose(cInputAr)
            print("cInputAr = " + str(cInputAr))

            cInputAr.tofile("cInputAr.txt", sep=" ", format="%s")

            print("cInputAr det = " + str(np.linalg.det(cInputAr)))

            inWeightsAr = np.linalg.inv(cInputAr)
            print("inWeightsAr = " + str(inWeightsAr))

            self.inWeightsAr = inWeightsAr
            inPreds = preds.transpose()
            self.inPreds = inPreds
            inPreds.tofile("inPreds.txt", sep=" ", format="%s")


