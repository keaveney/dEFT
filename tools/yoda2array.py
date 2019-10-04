""" module to convert yoda histograms in yoda files into numpy arrays
    so that the user can directly suppply the yoda file as input into a dEFT
    fit
"""
import numpy as np
from ROOT import TFile, TH1D
from array import *
import math

def convert(yodafile, histname):
    print "tools/yoda2array: converting yodafile to numpy array"
    ignore_strings = ["# ID", "# xlow", "Total   ", "Underflow", "Overflow"]
    read = "false"
    central_values = []

    #print "converting " + str(yodafile)
    with open(yodafile, 'r') as f:
        for line in f:
            if ( ("BEGIN" in line) & (histname in line)):
                read = "true" 
            if (read == "true"):
                if ((len(line.split("\t")) == 7) & (line.split("\t")[0] not in ignore_strings) ):
                    print "LINE " + str(line)
                    bin_width = float(line.split("\t")[1]) - float(line.split("\t")[0])
                    central_values.append(float(line.split("\t")[2]) / bin_width )
            if ("END" in line):
                read = "false"
    return  np.asarray(central_values)

def convert2root(yodafile, histname):
    print "tools/yoda2array: converting yodafile to root file"
    ignore_strings = ["# ID", "# xlow", "Total   ", "Underflow", "Overflow"]
    read = "false"
    central_values = []
    uncertainties = []
    bin_edges = []
    last_line = ""
    n_bin_line = 0

    with open(yodafile, 'r') as f:
        for line in f:
            #print "RAW-LINE = " + line  + " histname  "+ str(histname)
            if ( ("BEGIN" in line) & (histname in line)):
                read = "true"
            if (read == "true"):
                if ((len(line.split("\t")) == 7) & (line.split("\t")[0] not in ignore_strings) ):
                    n_bin_line = n_bin_line + 1
                    if(n_bin_line) == 1:
                         bin_edges.append(float(line.split("\t")[0]))
                         bin_edges.append(float(line.split("\t")[1]))
                    else:
                         bin_edges.append(float(line.split("\t")[1]))
                    bin_width = float(line.split("\t")[1]) - float(line.split("\t")[0])
                    central_values.append(float(line.split("\t")[2]) / bin_width )
                    uncertainties.append(float(math.sqrt(float(line.split("\t")[3]))) / bin_width )

            if ("END" in line):
                read = "false"

        bin_edges = sorted(bin_edges)
        edgesArray = array('d',bin_edges)
        #print "edges = " + str(edgesArray) 
        histo = TH1D("", "", len(bin_edges)-1, edgesArray)
        # fill histo
        for bin in range(0, len(central_values)):
            histo.SetBinContent(bin+1, central_values[bin])
            histo.SetBinError(bin+1, uncertainties[bin])
            #print "error bar is " + str(uncertainties[bin])
    return  histo
