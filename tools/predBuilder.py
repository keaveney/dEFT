import numpy as np
import itertools
from sklearn.linear_model import LinearRegression

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
                combs = itertools.combinations(list(row), 2)
                num_rows = np.shape(ci)[0]
                for comb in combs:
                    X = np.append(X, comb[0]*comb[1])
    
            X = X.reshape(num_rows, self.triangularNumber(len(ci[0])))

            return X
            
        def makeCoeffPoint(self, ci):
            X = np.array([])
            X = np.append(X, ci**2)
            combs = itertools.combinations(list(ci), 2)
            for comb in combs:
                X = np.append(X, comb[0]*comb[1])
    
            X = X.reshape(1, self.triangularNumber(len(ci)))
        
            return X
            
        def initRM(self, nOps, samples, preds):
  
            if (len(preds) <= self.nSamples(nOps)):
                raise TypeError("morphing with " + str(nOps) + " coefficients requires at least " + str(self.nSamples(nOps)) + " samples,  but only " + str(len(preds)) + " are provided")
            
            self.nOps = nOps
            cInputAr = np.array([])
            
            #convert to vector of coefficient factors to linearise the morphing
            cInputAr = self.makeCoeff(samples)
            
            # define model
            model = LinearRegression()

            # fit model
            model.fit(cInputAr, preds)

            self.model = model
            
        def makeRMPred(self,c):

            c = np.append(1.0, c)
            cInputAr = self.makeCoeffPoint(c)
            pred = self.model.predict(cInputAr)
            
            return pred[0]
        
        def init(self, nOps, samples, preds):
            # 'exact' morphing, soon to be defunct
        
            if (len(preds) <= self.nSamples(nOps)):
                raise TypeError("morphing with " + str(nOps) + " coefficients requires at least " + str(self.nSamples(nOps)) + " samples,  but only " + str(len(preds)) + " are provided")

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
            cInputAr.tofile("cInputAr.txt", sep=" ", format="%s")

            inWeightsAr = np.linalg.inv(cInputAr)

            self.inWeightsAr = inWeightsAr
            inPreds = preds.transpose()
            self.inPreds = inPreds
            inPreds.tofile("inPreds.txt", sep=" ", format="%s")


