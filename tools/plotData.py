import numpy as np
from array import *
import math
import sys
import matplotlib.pyplot as pl
from matplotlib import gridspec
import matplotlib.colors as colors

from mpl_toolkits.axes_grid1 import make_axes_locatable

obs = "ptZ"
procDirName = "tzq_train"
procFileStem = "TZQ_7ops_train_"

def cmsCorr():

    pt_ll = np.array([0.019132, 0.0529744, 0.0859099, 0.123766, 0.108097, 0.0344017, 0.00142478])
    #unc vec of data
    pt_ll_unc = np.array([0.00154984, 0.00370626, 0.00561757, 0.00776465, 0.00639847, 0.00209706, 0.000108188])
    
    #cov matrix of data
    pt_ll_cov = np.array([[2.53276e-06, 5.05439e-06, 7.52378e-06, 1.01732e-05, 8.23324e-06, 2.55454e-06, 1.1711e-07],
    [5.05439e-06, 1.40578e-05, 1.93376e-05, 2.61039e-05, 2.15583e-05, 6.90085e-06, 3.10456e-07],
    [7.52378e-06, 1.93376e-05, 3.18233e-05, 4.26285e-05, 3.49307e-05, 1.1103e-05,5.20878e-07],
    [1.01732e-05, 2.61039e-05, 4.26285e-05, 6.06592e-05, 4.8833e-05, 1.56136e-05, 7.25042e-07],
    [8.23324e-06, 2.15583e-05, 3.49307e-05, 4.8833e-05, 4.10895e-05, 1.31334e-05, 6.2368e-07],
    [2.55454e-06, 6.90085e-06, 1.1103e-05, 1.56136e-05,1.31334e-05,4.43529e-06, 2.02826e-07],
    [1.1711e-07, 3.10456e-07, 5.20878e-07, 7.25042e-07, 6.2368e-07, 2.02826e-07,1.20047e-08]])
    
    outer_prod = np.outer(pt_ll_unc,pt_ll_unc)
    corr  = np.true_divide(pt_ll_cov, outer_prod)
    
    return corr[:5,:5]

def plotHWU(hwufile):
    central_values = []
    file_string = hwufile
    
    lumi = 3000000 # pb^{-1}
    eff = 0.4
    #BR = 0.0115
    BR = 0.026
    flatSys = 0.05
    
    scale = lumi*eff*BR
    
    with open(file_string, 'r') as f:
        for line in f:
            if (obs in line):
                read = "true"
            if (read == "true"):
                #print("N " + str(len(line.split("   "))))

                if ((len(line.split("   ")) == 11)):
                    #print("LINE " + str(line))
                    #bin_width = float(line.split("\t")[1]) - float(line.split("\t")[0])
                    central_values.append(float(line.split("   ")[2]) )
                    #print(line.split("   ")[2][1:]  )
                if ("<\histogram>" in line):
                    read = "false"
    #print(str(central_values) + ",")
    
    #just overwrite with SM pred from regMorph model = in units of pb (cross section, not divided by bin width)
    central_values = [0.05432116, 0.03366438, 0.01031508, 0.00257308, 0.00099607]
    
    #this needs a k-factor to bring it to the NLO predcition
    lo_xs = np.sum(central_values)
    nlo_xs = 0.160
    print("lo_xs = " + str(lo_xs))
    k_fac = nlo_xs/lo_xs
    print("incl. kfac = " + str(k_fac))

    central_values = [x*k_fac for x in central_values]
    scaled_xs = np.sum(central_values)
    print("scaled_xs = " + str(scaled_xs))

    diff_xs = [x/100.0 for x in central_values]
        
    final_values = [x * scale for x in central_values]
    stat_unc = np.sqrt(final_values)
    
    sys_unc = [x * flatSys for x in final_values]
    diff_xs_sys_unc = [x * flatSys for x in diff_xs]

    rel_stat_unc = (stat_unc/final_values)
    
    stat_unc_diffxs = rel_stat_unc * diff_xs
    #print("stat_unc_diffxs" + str(stat_unc_diffxs))

    xvals = np.linspace(50, 450, 5)

    
    #fill special arrays for error bands
    bin_edges = np.linspace(0, 500, 6)
    bin_edges_for_errors = np.array([])
    
    final_values_for_band = np.array([])
    final_diffxs_values_for_band = np.array([])

    stat_unc_for_band = np.array([])
    stat_unc_for_ratio = np.array([])
    
    diffxs_stat_unc_for_band = np.array([])
    diffxs_stat_unc_for_ratio = np.array([])
    
    diffxs_sys_unc_for_band = np.array([])
    diffxs_sys_unc_for_ratio = np.array([])
    
    sys_unc_for_band = np.array([])
    sys_unc_for_ratio = np.array([])

    for bin in range(0, len(xvals)):
        bin_edges_for_errors = np.append(xvals[len(xvals) -  bin - 1] + 50.0, bin_edges_for_errors)
        bin_edges_for_errors = np.append(xvals[len(xvals) -  bin - 1] - 50.0, bin_edges_for_errors)
        final_values_for_band = np.append(final_values[len(xvals) -  bin - 1], final_values_for_band)
        final_values_for_band = np.append(final_values[len(xvals) -  bin - 1], final_values_for_band)
        final_diffxs_values_for_band = np.append(diff_xs[len(xvals) -  bin - 1], final_diffxs_values_for_band)
        final_diffxs_values_for_band = np.append(diff_xs[len(xvals) -  bin - 1], final_diffxs_values_for_band)
        
        stat_unc_for_band = np.append(stat_unc[len(xvals) -  bin - 1], stat_unc_for_band)
        stat_unc_for_band = np.append(stat_unc[len(xvals) -  bin - 1], stat_unc_for_band)
        sys_unc_for_band = np.append(sys_unc[len(xvals) -  bin - 1], sys_unc_for_band)
        sys_unc_for_band = np.append(sys_unc[len(xvals) -  bin - 1], sys_unc_for_band)
        
        diffxs_stat_unc_for_band = np.append(stat_unc_diffxs[len(xvals) -  bin - 1], diffxs_stat_unc_for_band)
        diffxs_stat_unc_for_band = np.append(stat_unc_diffxs[len(xvals) -  bin - 1], diffxs_stat_unc_for_band)
        diffxs_sys_unc_for_band = np.append(diff_xs_sys_unc[len(xvals) -  bin - 1], diffxs_sys_unc_for_band)
        diffxs_sys_unc_for_band = np.append(diff_xs_sys_unc[len(xvals) -  bin - 1], diffxs_sys_unc_for_band)
                
        stat_unc_for_ratio = np.append( ((stat_unc[len(xvals) -  bin - 1])/(final_values[len(xvals) -  bin - 1]))*100.0 , stat_unc_for_ratio)
        stat_unc_for_ratio = np.append( ((stat_unc[len(xvals) -  bin - 1])/(final_values[len(xvals) -  bin - 1]))*100.0 , stat_unc_for_ratio)
        sys_unc_for_ratio = np.append( ((sys_unc[len(xvals) -  bin - 1])/(final_values[len(xvals) -  bin - 1]))*100.0 , sys_unc_for_ratio)
        sys_unc_for_ratio = np.append( ((sys_unc[len(xvals) -  bin - 1])/(final_values[len(xvals) -  bin - 1]))*100.0 , sys_unc_for_ratio)
        
        diffxs_stat_unc_for_ratio = np.append( ((stat_unc_diffxs[len(xvals) -  bin - 1])/(diff_xs[len(xvals) -  bin - 1]))*100.0 , diffxs_stat_unc_for_ratio)
        diffxs_stat_unc_for_ratio = np.append( ((stat_unc_diffxs[len(xvals) -  bin - 1])/(diff_xs[len(xvals) -  bin - 1]))*100.0 , diffxs_stat_unc_for_ratio)
        diffxs_sys_unc_for_ratio = np.append( ((diff_xs_sys_unc[len(xvals) -  bin - 1])/(diff_xs[len(xvals) -  bin - 1]))*100.0 , diffxs_sys_unc_for_ratio)
        diffxs_sys_unc_for_ratio = np.append( ((diff_xs_sys_unc[len(xvals) -  bin - 1])/(diff_xs[len(xvals) -  bin - 1]))*100.0 , diffxs_sys_unc_for_ratio)
        
    #print("befe = " + str(bin_edges_for_errors))
    #print("fvefe = " + str(final_values_for_band  ))
    #print("uefe = " + str(stat_unc_for_band ))
    zeros_for_ratio =  np.zeros(len(xvals)*2)

    fig = pl.figure()
    spec = gridspec.GridSpec(ncols=1, nrows=2, height_ratios=[3, 1], hspace=0.0)
    ax0 = fig.add_subplot(spec[0])
    ax0.set_ylabel(r'$N_{tWZ}$', fontsize = 16)
    ax0.yaxis.set_label_coords(-0.082, 0.83)
    ax0.set_xlim(0.0, 500.0)
    ax0.set_xticklabels([])

    ax0.errorbar(xvals, final_values, fmt="",  ls='none', xerr=50., label ='central values')
    ax0.fill_between(bin_edges_for_errors, final_values_for_band-sys_unc_for_band, final_values_for_band+sys_unc_for_band, facecolor='orange', alpha=0.3, edgecolor='none', label ='syst. unc.')
    ax0.fill_between(bin_edges_for_errors, final_values_for_band-stat_unc_for_band, final_values_for_band+stat_unc_for_band, facecolor='#1f77b4', alpha=0.3, edgecolor='none', label ='stat. unc.')
    pl.legend()

    #ratio plot
    ax2 = fig.add_subplot(spec[1])
    ax2.errorbar(xvals, np.zeros(len(xvals)), fmt="",  ls='none', xerr=50.0, label=r'$N_{tWZ}$')
    ax2.fill_between(bin_edges_for_errors, zeros_for_ratio-sys_unc_for_ratio, zeros_for_ratio+sys_unc_for_ratio, facecolor='orange', alpha=0.3, edgecolor='none')
    ax2.fill_between(bin_edges_for_errors, zeros_for_ratio-stat_unc_for_ratio, zeros_for_ratio+stat_unc_for_ratio, facecolor='#1f77b4', alpha=0.3, edgecolor='none')

    ax2.set_xlim(0.0, 500.0)
    ax2.xaxis.set_label_coords(0.85, -0.25)
    labelx = ax2.set_xlabel(r'$p^{Z}_{T} \;[GeV]$', fontsize = 16)
    labely = ax2.set_ylabel(r'rel. unc. [%]', fontsize = 10)
    
    fig.savefig("twz_expected.png")
    pl.close()
    
    ### diff xsection plot
    fig = pl.figure()
    spec = gridspec.GridSpec(ncols=1, nrows=2, height_ratios=[3, 1], hspace=0.0)
    ax0 = fig.add_subplot(spec[0])
    ax0.set_ylabel(r'$\frac{d\sigma_{tWZ}}{p^{Z}_{T}}\;[pb\; GeV^{-1}] $', fontsize = 15)
    ax0.yaxis.set_label_coords(-0.115, 0.75)
    ax0.set_xlim(0.0, 500.0)
    ax0.set_xticklabels([])

    ax0.errorbar(xvals, diff_xs, fmt="",  ls='none', xerr=50., label ='central values')
    ax0.fill_between(bin_edges_for_errors, final_diffxs_values_for_band-diffxs_sys_unc_for_band, final_diffxs_values_for_band+diffxs_sys_unc_for_band, facecolor='orange', alpha=0.3, edgecolor='none', label ='syst. unc.')
    ax0.fill_between(bin_edges_for_errors, final_diffxs_values_for_band-diffxs_stat_unc_for_band, final_diffxs_values_for_band+diffxs_stat_unc_for_band, facecolor='#1f77b4', alpha=0.3, edgecolor='none', label ='stat. unc.')
    #ax0.ticklabel_format(style='sci')
    ax0.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3), useOffset=False)

    pl.legend()

    #ratio plot
    ax2 = fig.add_subplot(spec[1])
    ax2.errorbar(xvals, np.zeros(len(xvals)), fmt="",  ls='none', xerr=50.0, label=r'$N_{tWZ}$')
    ax2.fill_between(bin_edges_for_errors, zeros_for_ratio-diffxs_sys_unc_for_ratio, zeros_for_ratio+sys_unc_for_ratio, facecolor='orange', alpha=0.3, edgecolor='none')
    ax2.fill_between(bin_edges_for_errors, zeros_for_ratio-diffxs_stat_unc_for_ratio, zeros_for_ratio+diffxs_stat_unc_for_ratio, facecolor='#1f77b4', alpha=0.3, edgecolor='none')

    ax2.set_xlim(0.0, 500.0)
    ax2.xaxis.set_label_coords(0.82, -0.25)

    labelx = ax2.set_xlabel(r'$p^{Z}_{T} \;[GeV]$', fontsize = 15)
    labely = ax2.set_ylabel(r'rel. unc. [%]', fontsize = 12)
    
    fig.tight_layout()

    fig.savefig("twz_diff_xs.pdf")
    
    ### diff xsection plot with operator effects
    fig = pl.figure(figsize=(5, 9))
    pl.rcParams['image.cmap'] = 'Pastel1'
    spec = gridspec.GridSpec(ncols=1, nrows=8, height_ratios=[2.5, 1, 1, 1, 1, 1, 1, 1], hspace=0.0)
    ax0 = fig.add_subplot(spec[0])
    ax0.set_ylabel(r'$\frac{d\sigma_{tWZ}}{p^{Z}_{T}}\;[pb\; GeV^{-1}] $', fontsize = 15)
    ax0.yaxis.set_label_coords(-0.06, 0.6)
    ax0.set_xlim(0.0, 500.0)
    ax0.set_xticklabels([])

    ax0.errorbar(xvals, diff_xs, fmt="",  ls='none', xerr=50., label ='central values')
    ax0.fill_between(bin_edges_for_errors, final_diffxs_values_for_band-diffxs_sys_unc_for_band, final_diffxs_values_for_band+diffxs_sys_unc_for_band, facecolor='orange', alpha=0.3, edgecolor='none', label ='syst. unc.')
    ax0.fill_between(bin_edges_for_errors, final_diffxs_values_for_band-diffxs_stat_unc_for_band, final_diffxs_values_for_band+diffxs_stat_unc_for_band, facecolor='#1f77b4', alpha=0.3, edgecolor='none', label ='stat. unc.')
    #ax0.ticklabel_format(style='sci')
    ax0.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3), useOffset=False)

    pl.legend()

    #ratio plot
    ax2 = fig.add_subplot(spec[1])
    ax2.errorbar(xvals, np.zeros(len(xvals)), fmt="",  ls='none', xerr=50.0, label=r'$N_{tWZ}$')
    ax2.fill_between(bin_edges_for_errors, zeros_for_ratio-diffxs_sys_unc_for_ratio, zeros_for_ratio+sys_unc_for_ratio, facecolor='orange', alpha=0.3, edgecolor='none')
    ax2.fill_between(bin_edges_for_errors, zeros_for_ratio-diffxs_stat_unc_for_ratio, zeros_for_ratio+diffxs_stat_unc_for_ratio, facecolor='#1f77b4', alpha=0.3, edgecolor='none')
    ax2.set_xlim(0.0, 500.0)
    ax2.set_xticklabels([])

    smeft_preds = np.array( [ [0.05394412, 0.03393586, 0.01046979, 0.00303195, 0.00095483],
    [0.0677107,  0.04120467, 0.01246974, 0.00375727, 0.0012688],
    [0.05572474, 0.03449813, 0.01049137, 0.00320375, 0.00104812],
    [0.05441161, 0.03512581, 0.01124541, 0.00374836, 0.0014506 ],
    [0.05289766, 0.03561033, 0.01316591, 0.00538783, 0.00258019],
    [0.06113224, 0.03949186, 0.01325608, 0.00492008, 0.00209033]])
    
    smeft_preds = np.true_divide(smeft_preds, 100.0)
    
    tex_labels = ["$c^{3}_\\phi Q}$", "$c^{M}_\\phi Q}$","$c_{tZ}$","$c_{t W}$","$c_{tG}$"]
    
    colors = ['b' ,'g', 'r', 'c', 'm', 'y']

    for p in range(1, len(smeft_preds)):
        ax2 = fig.add_subplot(spec[p+1])
        ax2.annotate(tex_labels[p-1], # this is the text
                 (50.0, 2.5), # these are the coordinates to position the label
                 textcoords="offset points", # how to position the text
                 xytext=(0,10), # distance from text to points (x,y)
                 ha='center')
        ax2.errorbar(xvals, ( 100.0*((smeft_preds[p] - smeft_preds[0]) /smeft_preds[0])), fmt="", alpha=0.8, ls='none', xerr=50.0, label=tex_labels[p-1])
        if (p <= 2):
            ax2.set_ylim(0.0, 40.0)
        elif (p > 2):
            ax2.set_ylim(-0.5, 180.0)
            #ax2.set_yticklabels([0,50])

        ax2.set_xlim(0.0, 500.0)
        ax2.set_xticklabels([])

        #ax2.yaxis.set_ticks([0,1,3,4,5,6])

    ax2.set_xticklabels([0,100,200,300,400,500])
    ax2.set_xlim(0.0, 500.0)
    ax2.xaxis.set_label_coords(0.82, -0.4)

    labelx = ax2.set_xlabel(r'$p^{Z}_{T} \;[GeV]$', fontsize = 15)
    labely = ax2.set_ylabel(r'(SMEFT-SM)/SM', fontsize = 12)
    ax2.yaxis.set_label_coords(-0.115, 3.2)

    fig.tight_layout()

    fig.savefig("twz_diff_xs_smeft.pdf")
    
    
    pl.close()
    
    sys_unc_for_cov = [x * flatSys for x in central_values]
    
    print("central values  = " + str(central_values))
    print("final values  = " + str(final_values))
    print("sqrt(final values)  = " + str(np.sqrt(final_values)))
    stat_unc_central = [ (np.sqrt(final_values[x])/final_values[x])*(central_values[x]) for x in range(0,len(central_values))]
    print("stat unc (central values)  = " + str(stat_unc_central))
    sys_unc_central = [x * flatSys for x in central_values]
    print("syst(central values)  = " + str(sys_unc_central))
    total_unc = [np.sqrt(stat_unc_central[x]**2 + sys_unc_central[x]**2 ) for x in range(0,len(central_values))]
    print("unc (central values)  = " + str(total_unc))

    #total_unc = [ np.sqrt((rel_stat_unc[x]*central_values[x])**2 + (sys_unc_for_cov[x])**2) for x in range(0,5)]
    total_unc_diff_xs = [ np.sqrt(diffxs_sys_unc_for_band[x]**2 + diffxs_stat_unc_for_band[x]**2) for x in range(0,5)]

    return np.array(total_unc)

unc = plotHWU("sm_pred")

def plotDataCov(unc):
    
    #get estimated correlation matrix of CMS measurement
    corr = cmsCorr()
    
    print("CMS CORR:")
    print(repr(corr))

    ####### SIMPLE METHOD ########
    #cov_mat = np.diag(unc**2)
    #diag = np.array([])
    
    #cov1 = 0.05
    #cov1 = 0.5

    #for bin in range(1, len(unc)):
    #    diag = np.append(diag,((unc[bin])*(unc[bin-1])*cov1))
    #cov_mat = cov_mat + np.diag(diag, -1) + np.diag(diag, 1)
    ###############################

    outer_prod = np.outer(unc, unc)
    cov_mat = np.multiply(corr, outer_prod)

    print("total covariance matrix:")
    print(repr(cov_mat))

    fig = pl.figure()
#    vmin=0, vmax=99
    im = pl.matshow(cov_mat, cmap=pl.cm.Blues, norm=colors.LogNorm(vmin=2e-9, vmax=2e-5), extent=[0.0, 500.0, 500.0,0.0])
#    im = pl.matshow(cov_mat, cmap=pl.cm.Blues, vmin=2e-9, vmax=2e-5,extent=[0.0, 500.0, 500.0,0.0])

    pl.xlabel(r'$p^{Z}_{T} \;[GeV]$', fontsize = 15)
    pl.ylabel(r'$p^{Z}_{T} \;[GeV]$', labelpad=30, fontsize = 15)

    ax = pl.gca()

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = pl.colorbar(im, cax=cax)
    cbar.set_label('covariance', x=0.7, y=0.83, fontsize = 15)

    ax.tick_params(axis="x", bottom=True, top=False, labelbottom=True, labeltop=False)
    ax.xaxis.set_label_coords(0.85, -0.08)
    ax.yaxis.set_label_coords(-0.12, 0.85)

    pl.savefig("data_cov.pdf",bbox_inches='tight')
    pl.close()
    
    return
    
plotDataCov(unc)



    
    
