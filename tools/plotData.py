import numpy as np
from array import *
import math
import sys
import matplotlib.pyplot as pl
from matplotlib import gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable

nBins = 7
nOps = 7
obs = "ptZ"

procDirName = "tzq_train"
procFileStem = "TZQ_7ops_train_"

def plotHWU(hwufile):
    central_values = []
    file_string = hwufile
    
    lumi = 3000000 # pb^{-1}
    eff = 0.4
    BR = 0.0115
    
    scale = lumi*eff*BR
    
    with open(file_string, 'r') as f:
        for line in f:
            if (obs in line):
                read = "true"
            if (read == "true"):
                print("N " + str(len(line.split("   "))))

                if ((len(line.split("   ")) == 11)):
                    print("LINE " + str(line))
                    #bin_width = float(line.split("\t")[1]) - float(line.split("\t")[0])
                    central_values.append(float(line.split("   ")[2]) )
                    #print(line.split("   ")[2][1:]  )
                if ("<\histogram>" in line):
                    read = "false"
    print(str(central_values) + ",")
    
    central_values = [0.05300887, 0.03509872, 0.0088587,  0.00106592, 0.00062776]
    dsig =  [(x*10.0) for x in central_values]
    
    final_values = [x * scale for x in central_values]
    stat_unc = np.sqrt(final_values)
    rel_stat_unc = (stat_unc/final_values)

    xvals = np.linspace(50, 450, 5)

    print("xvals " + str(len(xvals)))
    print("final_values " + str(len(final_values)))
    print("sum of central values " + str(sum(central_values)))
    print("rel stat unc  " + str(rel_stat_unc))
    
    #fill special arrays for error bands
    bin_edges = np.linspace(0, 500, 6)
    bin_edges_for_errors = np.array([])
    final_values_for_band = np.array([])

    stat_unc_for_band = np.array([])
    stat_unc_for_ratio = np.array([])

    for bin in range(0, len(xvals)):
        bin_edges_for_errors = np.append(xvals[len(xvals) -  bin - 1] + 50.0, bin_edges_for_errors)
        bin_edges_for_errors = np.append(xvals[len(xvals) -  bin - 1] - 50.0, bin_edges_for_errors)
        final_values_for_band = np.append(final_values[len(xvals) -  bin - 1], final_values_for_band)
        final_values_for_band = np.append(final_values[len(xvals) -  bin - 1], final_values_for_band)
        stat_unc_for_band = np.append(stat_unc[len(xvals) -  bin - 1], stat_unc_for_band)
        stat_unc_for_band = np.append(stat_unc[len(xvals) -  bin - 1], stat_unc_for_band)
        stat_unc_for_ratio = np.append( ((stat_unc[len(xvals) -  bin - 1])/(final_values[len(xvals) -  bin - 1]))*100.0 , stat_unc_for_ratio)
        stat_unc_for_ratio = np.append( ((stat_unc[len(xvals) -  bin - 1])/(final_values[len(xvals) -  bin - 1]))*100.0 , stat_unc_for_ratio)

    print("befe = " + str(bin_edges_for_errors))
    print("fvefe = " + str(final_values_for_band  ))
    print("uefe = " + str(stat_unc_for_band ))
    zeros_for_ratio =  np.zeros(len(xvals)*2)

    fig = pl.figure()
    spec = gridspec.GridSpec(ncols=1, nrows=2, height_ratios=[3, 1], hspace=0.0)
    ax0 = fig.add_subplot(spec[0])
    ax0.set_ylabel(r'$N_{tWZ}$', fontsize = 16)
    ax0.yaxis.set_label_coords(-0.075, 0.83)
    ax0.errorbar(xvals, final_values, fmt="",  ls='none', xerr=50.)
    ax0.fill_between(bin_edges_for_errors, final_values_for_band-stat_unc_for_band, final_values_for_band+stat_unc_for_band, facecolor='#F0F8FF', alpha=1.0, edgecolor='none')
    
    #ax1 = ax0.twinx()
    #ax1.set_ylabel(r'$\frac{d \sigma}{d p^{Z}_{T}}\; [fb/GeV]$', fontsize = 14)
    #ax1.yaxis.set_label_coords(1.06, 0.79)
    #ax1.tick_params(axis='y')
    #ax1.errorbar(xvals, dsig, fmt="",  ls='none', alpha=0.0, xerr=50.0)
    
    #ratio plot
    ax2 = fig.add_subplot(spec[1])
    ax2.errorbar(xvals, np.zeros(len(xvals)), fmt="",  ls='none', xerr=50.0, label=r'$N_{tWZ}$')
    ax2.fill_between(bin_edges_for_errors, zeros_for_ratio-stat_unc_for_ratio, zeros_for_ratio+stat_unc_for_ratio, facecolor='#F0F8FF', alpha=1.0, edgecolor='none')

    ax2.xaxis.set_label_coords(0.85, -0.25)
    labelx = ax2.set_xlabel(r'$p^{Z}_{T} \;[GeV]$', fontsize = 16)
    labely = ax2.set_ylabel(r'rel. unc. [%]', fontsize = 10)
    
    fig.savefig("sm.png")
    pl.close()
    
    return rel_stat_unc*central_values

stat_unc = plotHWU("sm_pred")

def plotDataCov(stat_unc):
    
    cov_mat = np.diag(stat_unc**2)
    diag = np.array([])
    
    cov1 = 0.05

    for bin in range(1, len(stat_unc)):
        diag = np.append(diag,((stat_unc[bin])*(stat_unc[bin-1])*cov1))

    cov_mat = cov_mat + np.diag(diag, -1) + np.diag(diag, 1)

    print(repr(cov_mat))

    pl.figure()
    im = pl.matshow(cov_mat, cmap=pl.cm.Blues, extent=[0.0,500.0,0.0, 500.0])
    pl.xlabel(r'$p^{Z}_{T} \;[GeV]$')
    pl.ylabel(r'$p^{Z}_{T} \;[GeV]$')

    ax = pl.gca()

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    pl.colorbar(im, cax=cax)
    ax.tick_params(axis="x", bottom=True, top=False, labelbottom=True, labeltop=False)
    ax.xaxis.set_label_coords(0.85, -0.08)
    ax.yaxis.set_label_coords(-0.12, 0.9)

    pl.savefig("data_cov.png")
    pl.close()
    
    return
    
plotDataCov(stat_unc)

