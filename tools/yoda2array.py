""" module to convert yoda histograms in yoda files into numpy arrays
    so that the user can directly suppply the yoda file as input into a dEFT
    fit
"""

import numpy as np

def convert(yodafile, histname):
    print "tools/yoda2array: converting yodafile to numpy array"
    ignore_strings = ["# ID", "# xlow", "Total   ", "Underflow", "Overflow"]
    read = "false"
    central_values = []

    with open(yodafile, 'r') as f:
        for line in f:
            if ( ("BEGIN" in line) & (histname in line)):
                read = "true"
            if (read == "true"):
                if ((len(line.split("\t")) == 7) & (line.split("\t")[0] not in ignore_strings) ):
                    central_values.append(float(line.split("\t")[2]))
            if ("END" in line):
                read = "false"
    return  np.asarray(central_values)




