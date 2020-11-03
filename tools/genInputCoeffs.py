#script to generate suitable sets of coefficient values
# for predictions for an observable in an EFT
import numpy as np
import random
from predBuilder import predBuilder

#ops = ["4", "5", "6", "8", "10"]

# 10 = ctW, 16 = ctG R, 6 = cPhit, 8 = cPhitB, 4 = cPhiQM, 5 = cPhiQt

#ops = ["6", "8", "10"]

#ops = ["4", "5", "6", "8", "10", "16"]

#ops = ["22", "23", "24"]

#cHWB:3 cW,:6, cHB:10, cHt:15, ctZ:22, ctW:23, ctG:24,

ops = ["3", "6", "10", "15", "22", "23", "24"]

nOps = len(ops)

pb = predBuilder()

def genRandomCoeffSets(nOps):
    boostFact = 3.0
    nSamples = 2*pb.nSamples(nOps)
    nSamples = 10
    coeffSets = []
    
    for s in range(0, int(nSamples)):
        coeffSet = []
        coeffSet.append(1.0)

        for c in range(0, nOps):
            signs = [-1.0, 1.0]
            coeffSign = random.choice(signs)
            coeffMag = float(random.uniform(0,30))
            randCoeff = coeffSign*coeffMag
            coeffSet.append(randCoeff)

        coeffSets.append(coeffSet)
    #coeffSets = coeffSets.reshape(int(nSamples), int(nSamples))
    return coeffSets

def genRandomPreds(nOps):
    
    nSamples = pb.nSamples(nOps)
    #preds = np.array([])
    preds = []
    
    for s in range(0, int(nSamples)):
        #pred = np.array([])
        pred = []

        #can make loop here for binned predictions
        #pred = np.append(pred, random.uniform(-100.0,100.0))
        #preds = np.append(preds, pred)
        pred.append(random.uniform(-50.0,50.0))
        preds.append(pred)

    #np.reshape(preds, (int(nSamples),1))
    return preds

def writeProcCard(randCoeffSet):

    #print("import dim6top_LO_UFO")
    print("import SMEFTatNLO-NLO")
    print("define p = p b b~")
    print("define tp = t t~")
    print("define w = w+ w-")
    print("define l+ = e+ mu+")
    print("define l- = e- mu-")
    print("define vl = ve vm")
    print("define vl~ = ve~ vm~")
    print("define lept = l+ l- vl vl~")
    print("generate p p > tp w z QCD=1 QED=2 NP=1 [QCD]")
    #print("generate p p > tp w z FCNC=0 DIM6=1")
    #print("generate p p > t t~ FCNC=0 DIM6=1, (t > w+ b DIM6=0, w+ > lept lept DIM6=0),(t~ > w- b~ DIM6=0, w- > lept lept DIM6=0)")
    print("output /tmp/twz_train_7ops/")
        
    for sample in range(0, len(randCoeffSet)):
        print("launch /tmp/twz_train_7ops/")
        print("order=LO")
        print("fixed_order=ON")
        #if sample ==0:
            #print("order=LO")
            #print("fixed_order=ON")
            #print("madspin=OFF")
            #print("shower=OFF")
            #print("reweight=OFF")
            #print("set nevents=5000")
        for c in range(1, len(randCoeffSet[sample])):
            print("set DIM62F " + str(ops[c-1]) + " " + str(randCoeffSet[sample][c]))

randCoeffSet = genRandomCoeffSets(nOps)
randPreds = genRandomPreds(nOps)

print("random coeffs = " + str(randCoeffSet) + " random preds = " + str(randPreds) )

#pb.init(nOps,randCoeffSet, randPreds)

writeProcCard(randCoeffSet)

