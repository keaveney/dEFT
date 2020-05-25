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
                    
                print("couplings = " + str(couplings))
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

        def makePred(self,c):
            #ciTargetAr = np.array([[1.0], [ci], [cj], [ci**2], [cj**2], [ci*cj] ])
            ciTargetAr = np.array([])
            
            #SM term
            ciTargetAr = np.append(ciTargetAr, 1.0)
            
            #linear 'SM--Oi' interference terms
            for ci in c:
                ciTargetAr = np.append(ciTargetAr, ci)

            #squared 'Oi' interference terms
            for ci in c:
                ciTargetAr = np.append(ciTargetAr, ci**2)
            
            combs = list(itertools.combinations(c,2))
            
            #print("In makePred(): c = " +  str(c))
            #print("In makePred(): pred = " +  str(combs))

            #cross 'Oi--Oj' interference terms
            for comb in range(0, len(combs)):
                ciTargetAr = np.append(ciTargetAr, combs[comb][0]*combs[comb][1])

            #NEED TO EXPLICITLY CHECK IF THE STRUCTURE OF THE 'G' MATRIX MATCHES
            #THAT OF THE FIXED 'A' MATRIX.

            ciTargetAr.tofile("ciTargetAr.txt", sep=" ", format="%s")
            #print("In makePred(): ciTargetAr " + str(ciTargetAr))
            morphWeightsAr = self.inWeightsAr.dot(ciTargetAr)
            #print("In makePred(): inpred = " +  str(self.inPreds))
            #print("In makePred(): morphWeightsAr = " +  str(morphWeightsAr))

            pred = self.inPreds.dot(morphWeightsAr)

            #print("In makePred(): pred = " +  str(pred))

            return pred
