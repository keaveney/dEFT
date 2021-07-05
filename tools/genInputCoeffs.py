#script to generate suitable sets of coefficient values
# for predictions for an observable in an EFT
import numpy as np
import random
import sys
from predBuilder import predBuilder

#ops = ["4", "5", "6", "8", "10"]

# 10 = ctW, 16 = ctG R, 6 = cPhit, 8 = cPhitB, 4 = cPhiQM, 5 = cPhiQt

#ops = ["6", "8", "10"]

#ops = ["4", "5", "6", "8", "10", "16"]

#ops = ["22", "23", "24"]

#cHWB:3 cW,:6, cHB:10, cHt:15, ctZ:22, ctW:23, ctG:24,

nTasks  = 1
nSamples = 30

#ops = ["3", "6", "10", "15", "22", "23", "24"]

# SMEFTatNLO list
ops = ["12", "13", "15", "22", "23", "24"]

nOps = len(ops)

pb = predBuilder()

def genRandomCoeffSets(nOps):
    boostFact = 3.0
    #nSamples = 2*pb.nSamples(nOps)
    coeffSets = []

    for s in range(0, int(nSamples)):
        coeffSet = []
        coeffSet.append(1.0)

        for c in range(0, nOps):
            signs = [-1.0, 1.0]
            coeffSign = random.choice(signs)
            coeffMag = float(random.uniform(0,2))
            #coeffMag = float(random.uniform(0,30))
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
        #pred.append(random.uniform(-50.0,50.0))
        pred.append(random.uniform(-5.0,5.0))
        preds.append(pred)

    #np.reshape(preds, (int(nSamples),1))
    return preds

def writeProcCard(randCoeffSet):

    genFileName = "TWZ_6ops_gen.txt"
    genFile = open(genFileName, "w")

    #genFile.write("import dim6top_LO_UFO")
    genFile.write("import SMEFTatNLO-NLO")
    genFile.write("define p = p b b~")
    genFile.write("define tp = t t~")
    genFile.write("define w = w+ w-")
    genFile.write("define l+ = e+ mu+")
    genFile.write("define l- = e- mu-")
    genFile.write("define vl = ve vm")
    genFile.write("define vl~ = ve~ vm~")
    genFile.write("define lept = l+ l- vl vl~")
    genFile.write("generate p p > tp w z QCD=1 QED=2 NP=1 [QCD]")
    #genFile.write("generate p p > tp w z FCNC=0 DIM6=1")
    #genFile.write("generate p p > t t~ FCNC=0 DIM6=1, (t > w+ b DIM6=0, w+ > lept lept DIM6=0),(t~ > w- b~ DIM6=0, w- > lept lept DIM6=0)")
    #genFile.write("output /tmp/twz_train_7ops/")
    genFile.write("output /scratch/james/twz_train_6ops/")
    genFile.close()
    
    randCoeffSetAr = np.array(randCoeffSet)
    print("SHAPE = " + str(randCoeffSetAr.shape))

    #randCoeffSetAr = np.reshape(randCoeffSetAr, (nTasks,-1))
    print("SHAPE = " + str(randCoeffSetAr.shape[0]))
   
    pointsFileName = "pointsFile.txt"
    np.savetxt(pointsFileName,randCoeffSetAr, delimiter=',', newline='\n' )
    
    if ((nSamples % nTasks) != 0):
        sys.exit("nPoints must be an integer multiple of nTasks... exiting")
    
    for task in range(0,nTasks):
        trainFileName = "TWZ_6ops_train_" + str(task) + ".txt"
        #trainFileName = "TZQ_7ops_train_" + str(task) + ".txt"
        trainFile = open(trainFileName, "w")
    
        nSamplesPerTask = int(nSamples/nTasks)
        print("samples per task  = " +str(nSamplesPerTask))

        for point in range(0, nSamplesPerTask):
            trainFile.write("launch /scratch/james/twz_train_6ops_lo/\n")
            #trainFile.write("launch /scratch/james/tzq_train/\n")
            trainFile.write("order=LO\n")
            trainFile.write("fixed_order=ON\n")
            
            pointIndex = (task * ((nSamples/nTasks))) + point
            print("pointIndex = " +str(int(pointIndex)))
            for c in range(1, len(randCoeffSet[int(pointIndex)])):
                trainFile.write("set DIM62F " + str(ops[c-1]) + " " + str(randCoeffSet[int(pointIndex)][c])+ "\n")
        trainFile.close()

    for sample in range(0, len(randCoeffSet)):
        #print("launch /tmp/twz_train_7ops/")
        print("launch /scratch/james/twz_train_6ops_lo/")
        #print("launch /scratch/james/tzq_train/")
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

