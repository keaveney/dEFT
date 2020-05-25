import matplotlib.pyplot as pl
import matplotlib.image as image
import numpy as np
import corner
import os
import sys
import json

filenames = []
filenames.append(sys.argv[1])
filenames.append(sys.argv[2])

print("filenames = " + str(filenames))

div = 0.0
nf = 0
colours  = ["r", "b", "g"]

pl.figure()
for f in filenames:
    with open(f, 'r') as fs:
        lbl = str(f.split("/")[1].split(".")[0])
        print("label = " + str(f.split("/")[1].split(".")[0]))
        summ = json.load(fs)
        print(str(summ))
        label = f.split(".")[1]
        newX = [xElem + div for xElem in summ["x"]]
        div = div + 0.2
        fmtStr = "o" + colours[nf]
        nf = nf + 1
        fig = pl.errorbar(newX, summ["bestFit"], yerr=[summ["UncsDown"], summ["UncsUp"]], fmt=fmtStr, label=lbl)

axes = pl.gca()
axes.set_ylim([-12.0, 12.0])
pl.xticks(range(len(summ["labels"])), summ["labels"], rotation='45')
pl.legend(loc=2)
pl.savefig("test_bestFits.png")
pl.close()
