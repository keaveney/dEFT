import numpy as np
from array import *
import math

def convert(hwufile):
    #print "tools/yoda2array: converting yodafile to numpy array"
    ignore_strings = ["# ID", "# xlow", "Total   ", "Underflow", "Overflow"]
    read = "false"
    central_values = []
    final_array = []

    print("converting " + str(hwufile))
    
    for run in range(191,201):
        central_values = []
        if (run < 10):
            run_string = "run_0" + str(run) + "_LO/MADatNLO.HwU"
        else:
            run_string = "run_"  + str(run) + "_LO/MADatNLO.HwU"
            
        print(run_string)
        file_string = hwufile + run_string
        with open(file_string, 'r') as f:
            for line in f:
                if ( ("<histogram> 3" in line)):
                    read = "true"
                if (read == "true"):
                    if ((len(line.split("   ")) == 15)):
                        #print("LINE " + str(line))
                        #bin_width = float(line.split("\t")[1]) - float(line.split("\t")[0])
                        central_values.append(float(line.split("   ")[2]) )
                        #print(line.split("   ")[2][1:]  )
                if ("<\histogram>" in line):
                    read = "false"
        print(central_values)
        final_array = np.append(final_array, central_values)
    final_array.reshape(10, 3)
    return final_array
    
preds = convert("../analyses/Events/")
print(repr(preds.reshape(10,3)))


def convertRun(runfile):
    #print "tools/yoda2array: converting yodafile to numpy array"
    read = "false"
    central_values = np.array([])
    final_array = []

    print("converting " + str(runfile))
    
    central_values = []
            
    with open(runfile, 'r') as f:
        for line in f:
            if ( ("fixed_order=ON" in line)):
                read = "true"
            if (read == "true"):
                if (len(line.split(" ")) >3 ):
                    central_values = np.append(central_values,float(line.split(" ")[3]) )
                    #central_values.append(float(line.split(" ")[3]) )
            if ("#launch /tmp/twz_train/" in line):
                read = "false"
        print(central_values)
        central_values.reshape(20, 3)
    return central_values
    
#preds = convertRun("../run_twz_SMEFTATNLO_train.txt")
#print(preds.reshape(20,3))
