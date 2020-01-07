import numpy as np
import itertools


class predBuilder:
        def nSamples(self,n):
            return (n+1)*(n+2)/2
        
        def init(self, nOps, samples, preds):
            if (len(preds) != self.nSamples(nOps)):
                print("ERROR: constraining " + str(nOps) + " coefficients requires " + str(self.nSamples(nOps)) + " samples,  but " + str(len(preds)) + " are provided")

            self.nOps = nOps
            cInputAr = np.array([])
            
            for sample in samples:
                couplings = np.array([])
                combs = list(itertools.combinations(sample,2))

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

            inWeightsAr = np.linalg.inv(cInputAr)
            self.inWeightsAr = inWeightsAr
            inPreds = preds.transpose()
            self.inPreds = inPreds

        def makePred(self,c):
            #ciTargetAr = np.array([[1.0], [ci], [cj], [ci**2], [cj**2], [ci*cj] ])
            ciTargetAr = np.array([])
            
            #SM term
            ciTargetAr = np.append(ciTargetAr, 1.0 )
            #linear 'SM--Oi' interference terms
            for ci in c:
                ciTargetAr = np.append(ciTargetAr, ci)

            #squared 'SM--Oi' interference terms
            for ci in c:
                ciTargetAr = np.append(ciTargetAr, ci**2)

            combs = list(itertools.combinations(c,2))

            for comb in range(0, len(combs)):
                ciTargetAr = np.append(ciTargetAr, combs[comb][0]*combs[comb][1])

            morphWeightsAr = self.inWeightsAr.dot(ciTargetAr)
            pred = self.inPreds.dot(morphWeightsAr)

            return pred
