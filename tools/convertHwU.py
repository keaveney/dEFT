import numpy as np
from array import *
import math
import sys
import matplotlib.pyplot as pl
from matplotlib import gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable


startRun = 1
endRun   = 150

#startRun = 226
#endRun   = 243


#startRun = 163
#endRun   = 180

nBins = 5
nOps = 6
obs = "ptZ"

procDirName = "twz_train_6ops_lo"
procFileStem = "TWZ_6ops_train_"

def convert(hwufile):
    #print "tools/yoda2array: converting yodafile to numpy array"
    ignore_strings = ["# ID", "# xlow", "Total   ", "Underflow", "Overflow"]
    read = "false"
    central_values = []
    final_array = []
    
    for run in range(startRun,(endRun+1)):
        central_values = []
        if (run < 10):
            run_string = "run_0" + str(run) + "_LO/MADatNLO.HwU"
        else:
            run_string = "run_"  + str(run) + "_LO/MADatNLO.HwU"
            
        file_string = hwufile + run_string
        with open(file_string, 'r') as f:
            for line in f:
                if ( (obs in line)):
                    read = "true"
                if (read == "true"):
                    if ((len(line.split("   ")) == 15)):
                        #print("LINE " + str(line))
                        #bin_width = float(line.split("\t")[1]) - float(line.split("\t")[0])
                        central_values.append(float(line.split("   ")[2]) )
                        #print(line.split("   ")[2][1:]  )
                if ("<\histogram>" in line):
                    read = "false"
        print(str(central_values) + ",")
        final_array = np.append(final_array, central_values)
    final_array.reshape(((endRun-startRun)+1), nBins)
    return final_array
    
preds = convert("../analyses/twz_train_6ops_lo/Events/")
print(repr(preds.reshape(((endRun-startRun)+1),nBins)))

def convertRun(runfile):
    #print "tools/yoda2array: converting yodafile to numpy array"
    read = "false"
    central_values = np.array([])
    final_array = []

    #print("converting " + str(runfile))
    
    central_values = []
            
    with open(runfile, 'r') as f:
        for line in f:
            if ( ("fixed_order=ON" in line)):
                read = "true"
                central_values = np.append(central_values,float(1.00000001))
            if (read == "true"):
                if (len(line.split(" ")) >3 ):
                    central_values = np.append(central_values,float(line.split(" ")[3]) )
                    #central_values.append(float(line.split(" ")[3]) )
            if ("#launch /tmp/" + procDirName + "/" in line):
                read = "false"
        #print(central_values)
        #central_values.reshape(endRun, nOps)
    return central_values
    
totalPreds = np.array([])

for p in range(0, 5):
    fileName = procFileStem + str(p) + ".txt"
    preds = convertRun(fileName)
    totalPreds = np.concatenate((totalPreds, preds), axis=None)

totalPreds = totalPreds.reshape(((endRun-startRun)+1),(nOps+1))

print("shap totalPreds = " + str(totalPreds.shape))

np.set_printoptions(threshold=sys.maxsize)
#print(np.array2string(totalPreds, separator=','))

#print(repr(totalPreds.reshape(endRun,(nOps+1))))

#for pred in totalPreds:
#    print(str(repr(pred)) + ",")






