import matplotlib.pyplot as pl
import matplotlib.image as image
import matplotlib.ticker as mticker
from matplotlib import gridspec

import numpy as np
import corner
import os
import json

#from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('text')

#rc('font',**{'family':'serif','serif':['Palatino']})
#pl.rcParams['pdf.fonttype'] = 42


class summaryPlotter:
    def summarise(self, config, pb, sampler, samples):
                    
            #make directory to hold results of this run
            cmd = "mkdir -p " + str(config.params["config"]["run_name"]) + "_results"
            os.system(cmd)
            
            mcmc_params = np.mean(sampler.flatchain,axis=0)
            mcmc_params_cov = np.cov(np.transpose(sampler.flatchain))
        
            data_label = "Data" + " (" + config.run_name + ")"
            max_val = (1.5)*(max(config.params["config"]["data"]["central_values"]))
            min_val = (0.0)*(min(config.params["config"]["data"]["central_values"]))
            xlabel = config.observable
            ylabel = "d$\sigma_{tt}$/dX"

            #make plots of prediction at max point of posterior versus data
            coefficients = list(config.coefficients)
            pl.figure()
            label_string_bestfit = "best-fit: ("
            for c in range(0, len(config.coefficients)):
                if (c == (len(config.coefficients) - 1)):
                    label_string_bestfit = label_string_bestfit + coefficients[c] + " = " + '%.1f' % mcmc_params[c] + ")"
                else:
                    label_string_bestfit = label_string_bestfit + coefficients[c] + " = " + '%.1f' % mcmc_params[c] + ", "

            print(" xvals = " + str(config.x_vals)  +  "  pb.makeRMPred(mcmc_params) "   +  str(pb.makeRMPred(mcmc_params))  )
            pl.errorbar(config.x_vals, pb.makeRMPred(mcmc_params), fmt="m", xerr=0.0, yerr=0.0, label=label_string_bestfit)
            pl.errorbar(config.x_vals, config.params["config"]["data"]["central_values"], fmt="o",xerr=0.25, yerr=0.05, label=data_label)
            pl.axis([config.x_vals[0]-0.25, config.x_vals[len(config.x_vals)-1]+0.25, min_val, max_val])
            ax = pl.gca()
            labely = ax.set_xlabel(xlabel, fontsize = 18)
            labely = ax.set_ylabel(ylabel, fontsize = 18)
            ax.xaxis.set_label_coords(0.85, -0.065)
            ax.yaxis.set_label_coords(-0.037, 0.83)
            pl.legend(loc=2)
            pl.savefig(config.params["config"]["run_name"] + "_results/" + config.params["config"]["run_name"] + "_bestfit" + "_predictions.png")
            pl.close()

            #  make "corner" plot
            labels = []
            ranges = []
            #ranges  = np.full((len(config.params["config"]["model"]["prior_limits"].keys()),1), 0.95)
            for c in config.params["config"]["model"]["prior_limits"].keys():
                label = "$" + c + "$"
                labels.append(label)
                ranges.append(config.params["config"]["model"]["prior_limits"][c])

            fig = corner.corner(samples, labels=config.tex_labels,
                                label_kwargs={"fontsize": 18},
                                range=ranges,
                                #quantiles=[0.16, 0.84],
                                #levels=(1-np.exp(-0.5),),
                                truths=np.zeros(len(labels)),
                                show_titles=True, title_kwargs={"fontsize": 18})

            plotfilename = str(config.params["config"]["run_name"]) + "_results/" + str(config.params["config"]["run_name"]) + ".png"
            logo = image.imread('logo/dEFT_logo.png')

            ax0 = fig.add_subplot(999)
            ax0.axis('off')
            img = ax0.imshow(logo)

            fig.savefig(plotfilename)
            pl.close()
            
            bestFits = []
            marginUncsUp = []
            marginUncsDown = []
            marginUncsUp95 = []
            marginUncsDown95 = []
            x = []
    
            for i in range(len(labels)):
                mcmc = np.percentile(samples[:, i], [16, 50, 84])
                q = np.diff(mcmc)
                txt = "\mathrm{{{3}}} = {0:.3f}_{{-{1:.3f}}}^{{{2:.3f}}}"
                txt = txt.format(mcmc[1], q[0], q[1], labels[i])
                x.append(i)
                bestFits.append(mcmc[1])
                marginUncsUp.append(q[1])
                marginUncsDown.append(q[0])
                
                #q_16, q_50, q_84 = quantile(x, [0.16, 0.5, 0.84],weights=weights)
                
                mcmc95 = np.percentile(samples[:, i], [2.5, 50.0, 97.5])
                #print(" np.percentile(samples[:, i], [2.5, 50.0, 97.5]) =   " + str(np.percentile(samples[:, i], [2.5, 50.0, 97.5]))   )
                q = np.diff(mcmc95)
                #print("q =   " + str(q)  )

                marginUncsUp95.append(q[1])
                marginUncsDown95.append(q[0])

            pl.figure()
            pl.tight_layout()
            pl.gcf().subplots_adjust(bottom=0.15)

            #marginUncsDown95ForPlot = np.array(marginUncsDown95) - np.array(bestFits)
            #marginUncsUp95ForPlot = np.array(marginUncsUp95) - np.array(bestFits)
            
            marginUncsDown95ForPlot = np.array(marginUncsDown95)
            marginUncsUp95ForPlot = np.array(marginUncsUp95) 

            fig = pl.errorbar(x, np.zeros(len(labels)), yerr=[marginUncsDown95ForPlot, marginUncsUp95ForPlot], fmt='none', linestyle=None, elinewidth=18, label=r'all bins')
                        
            #marginUncsDown95ForPlot_validty = np.array([2.485571981823523569e+00, 3.814826737191868933e+00, 1.719548151262410673e+00, 2.106877931308388519e+00])
            #marginUncsUp95ForPlot_validty = np.array([1.750224475309016592e+00, 3.712112260384360152e+00, 1.859996135276679396e+00, 1.954559834197093338e+00])

            marginUncsDown95ForPlot_validty = np.array([2.497898126870446767e+00, 4.027889810657558911e+00, 1.820171638994192431e+00, 2.250182560665191733e+00])
            marginUncsUp95ForPlot_validty = np.array([1.955803247345035745e+00, 3.794527477091643064e+00, 1.965576305035904747e+00, 2.064912976805433065e+00])

            fig = pl.errorbar(x, np.zeros(len(labels)), yerr=[marginUncsDown95ForPlot_validty,marginUncsUp95ForPlot_validty], fmt='none', linestyle=None, alpha=0.3, elinewidth=18, label=r'last bin omitted')
            
            rangeUp = np.max(marginUncsUp95ForPlot_validty)*(1.2)
            rangeDown = (-1.0)*(np.max(marginUncsDown95ForPlot_validty))*(1.2)

            pl.ylim([rangeDown, rangeUp])

            pl.legend(loc=1, fontsize = 16)
                
            pl.xticks(range(len(labels)), config.tex_labels, rotation='45', fontsize = 25)
            pl.yticks(fontsize = 22)
            pl.locator_params(axis='y', nbins=6)
            plotTitle ="95% credible intervals"
            pl.title(plotTitle,fontsize = 28 )

            ax = pl.gca()
            ax.tick_params(axis='x', pad=-3)
            ax.yaxis.set_label_coords(-0.085, 0.75)
            labely = ax.set_ylabel(r"$c_{i}/\Lambda^{2}$ [$TeV^{-2}$]", fontsize = 22)
            #pl.legend(loc=1, fontsize = 21)
            pl.savefig(config.params["config"]["run_name"] + "_results/" + config.params["config"]["run_name"] + "_bestFits.png")
            pl.close()
            
            #Write medians/CIs to text file
            results_file_name = config.params["config"]["run_name"] + "_results/" + config.params["config"]["run_name"] + "_postFit" + "_results.txt"
            np.savetxt(results_file_name, (bestFits,marginUncsDown,marginUncsUp))   # x,y,z equal sized 1D arrays
            
            results_file_name = config.params["config"]["run_name"] + "_results/" + config.params["config"]["run_name"] + "_postFit" + "_results_95.txt"
            np.savetxt(results_file_name, (np.zeros(len(labels)),marginUncsDown95ForPlot,marginUncsUp95ForPlot))   # x,y,z equal sized 1D arrays
            
            ############################################
            ######## Corner with 1-d CL overlaid  ######
            ############################################
                        
            fig_overlay = corner.corner(samples, labels=config.tex_labels,
                    label_kwargs={"fontsize": 22},
                    labelpad=-0.12,
                    range=ranges,
                    color='k',
                    smooth=(0.8,0.8),
#                    quantiles=[0.16, 0.84],
                    #quantiles=[0.025, 0.975],
#                    levels=(1-np.exp(-0.5),),
                    #levels=(1 - np.exp( -2.0 ),),
                    truths=np.zeros(len(labels)),
                    show_titles=True,
                    title_kwargs={"fontsize": 15},
                    hist2d_kwargs={"fill_contours": True, "plot_density": True})
                    
            resplot = image.imread(config.params["config"]["run_name"] + "_results/" + config.params["config"]["run_name"] + "_bestFits.png")

            ax0 = fig_overlay.add_subplot(322)
            ax0.axis('off')
            img = ax0.imshow(resplot)
            logo = image.imread('logo/dEFT_logo.png')

            ax0 = fig_overlay.add_subplot(5,5,15)
            ax0.axis('off')
            img = ax0.imshow(logo)
            
            fig_overlay.tight_layout()

            plotfilename = str(config.params["config"]["run_name"] + "_results/" + config.params["config"]["run_name"] + "_overlay.pdf")

            fig_overlay.savefig(plotfilename)
            pl.close()
            
            # make plots of best fit prediction versus data
            pl.figure()
            label_string_bestfit = "best-fit prediction"

            pl.errorbar(config.x_vals, pb.makeRMPred(bestFits), fmt="m", xerr=0.0, yerr=0.0, label=label_string_bestfit)
            pl.errorbar(config.x_vals, config.params["config"]["data"]["central_values"], fmt="o",xerr=0.25, yerr=np.sqrt(np.diagonal(config.params["config"]["data"]["covariance_matrix"])), label=data_label)
            pl.axis([config.x_vals[0]-0.25, config.x_vals[len(config.x_vals)-1]+0.25, min_val, max_val])
            ax = pl.gca()
            labely = ax.set_xlabel(xlabel, fontsize = 18)
            labely = ax.set_ylabel(ylabel, fontsize = 18)
            ax.xaxis.set_label_coords(0.85, -0.065)
            ax.yaxis.set_label_coords(-0.037, 0.83)
            pl.legend(loc=1)
            pl.savefig(config.params["config"]["run_name"] + "_results/" + config.params["config"]["run_name"] + "_bestfit" + "_predictions.png")
            pl.close()
            
            #make fit summary json
            fitSummary = {}
            fitSummary["bestFit"] = bestFits
            fitSummary["bestFit"] = bestFits
            fitSummary["UncsUp"] = marginUncsUp
            fitSummary["UncsDown"] = marginUncsDown
            fitSummary["labels"] = labels
            fitSummary["x"] = x

            with open(config.params["config"]["run_name"] + ".json", 'w') as fs:
                json.dump(fitSummary, fs)

            pl.figure()

            # best-fit prediction and random samples compared to data
            uncX = np.absolute((np.diff(config.bins)))/2.0
            dataUncY = np.sqrt(np.diag(config.cov))
            pl.errorbar(config.x_vals, config.params["config"]["data"]["central_values"], fmt="o",xerr=uncX, yerr=dataUncY, label=data_label)

            inds = np.random.randint(len(samples), size=499)
            for ind in inds:
                sample = samples[ind]
                pl.plot(config.x_vals, pb.makeRMPred(sample), "C1",alpha=0.02)

            pl.plot(config.x_vals, pb.makeRMPred(samples[42]), "C1", label="500 random samples" ,alpha=0.1)

            ax = pl.gca()
            labelx = ax.set_xlabel(xlabel, fontsize = 18)
            labely = ax.set_ylabel(ylabel, fontsize = 18)
            ax.xaxis.set_label_coords(0.85, -0.065)
            ax.yaxis.set_label_coords(-0.037, 0.83)
            pl.legend(loc=1)
            pl.savefig(config.params["config"]["run_name"] + "_results/" + config.params["config"]["run_name"] + "_postFit" + "_predictions.png")
            pl.close()

            pl.figure()
            fig, axes = pl.subplots(len(labels), figsize=(10, 7), sharex=True)
            mc_samples = sampler.get_chain()
            if (len(labels)==1):
                axes.plot(mc_samples[:, :, i], "k", alpha=0.3)
                axes.set_xlim(0, len(mc_samples))
                axes.set_ylabel(labels[i])
                axes.yaxis.set_label_coords(-0.1, 0.5)
            else:
                for i in range(len(labels)):
                    ax = axes[i]
                    ax.plot(mc_samples[:, :, i], "k", alpha=0.3)
                    ax.set_xlim(0, len(mc_samples))
                    ax.set_ylabel(labels[i])
                    ax.yaxis.set_label_coords(-0.1, 0.5)

                axes[-1].set_xlabel("step number");

            pl.savefig(config.params["config"]["run_name"] + "_results/" + config.params["config"]["run_name"] + "_walkerPaths.png")
            pl.close()
                        
            #covariance matrix of coefficicents
            fig, ax = pl.subplots()
            #im = ax.imshow(mcmc_params_cov, cmap=pl.cm.Blues)
            
            ax.set_xticks(np.arange(len(config.tex_labels)))
            ax.set_yticks(np.arange(len(config.tex_labels)))
            # ... and label them with the respective list entries
            ax.set_xticklabels(config.tex_labels)
            ax.set_yticklabels(config.tex_labels)

            # fixing yticks with matplotlib.ticker "FixedLocator"
            #ticks_loc = ax.get_yticks().tolist()
            #ticks_loc = config.tex_labels
            
            #ax.yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
            #label_format = '{:,.0f}'

            #pl.imshow(mcmc_params_cov, cmap=pl.cm.Blues)
            #ax.set_xticks(np.arange(len(config.tex_labels)))

            #ax.set_yticklabels(ticks_loc)
            
            pl.savefig(config.params["config"]["run_name"] + "_results/" + config.params["config"]["run_name"] + "mcmc_params_cov.png")
            pl.close()
        
            #Write SM pred to text file for validation
            sm_file_name = config.params["config"]["run_name"] + "_results/" + config.params["config"]["run_name"] + "_postFit" + "_sm_pred.txt"
            f = open(sm_file_name, "w")
            sm_pred = str(repr(pb.makeRMPred(np.zeros(len(config.params["config"]["model"]["prior_limits"])))) )
            f.write(sm_pred)
            f.close()

            pl.figure()
            #pl.errorbar(config.x_vals, pb.makeRMPred([0.0,0.0,0.0,0.0,0.0,0.0]), fmt="",  ls='none', xerr=50., label ='SM')
            id = np.identity(len(config.params["config"]["model"]["prior_limits"]))
            for op in range(0, len(id)):
                pl.errorbar(config.x_vals, pb.makeRMPred(id[op]), fmt="",  ls='none', xerr=50., label=config.tex_labels[op])
                pl.plot(config.x_vals, pb.makeRMPred(id[op]), label=config.tex_labels[op])
            pl.legend()
            pl.savefig(config.params["config"]["run_name"] + "_results/" + config.params["config"]["run_name"] + "_benchmarks.png")
            pl.close()


            id = np.identity(len(config.params["config"]["model"]["prior_limits"]))
            for op in range(0, len(id)):
                fig = pl.figure()
                spec = gridspec.GridSpec(ncols=1, nrows=2, height_ratios=[3, 1], hspace=0.0)
                ax0 = fig.add_subplot(spec[0])
                ax0.set_ylabel(r'$\frac{d\sigma_{tWZ}}{p^{Z}_{T}}\;[pb/GeV] $', fontsize = 15)
                ax0.yaxis.set_label_coords(-0.115, 0.75)
                ax0.set_xlim(0.0, 500.0)
                ax0.set_xticklabels([])

                sm_pred =  pb.makeRMPred(np.zeros(len(config.params["config"]["model"]["prior_limits"])))
                ax0.errorbar(config.x_vals, sm_pred, fmt="",  ls='none', xerr=50., label ='SM')
            
                plot = ax0.errorbar(config.x_vals, pb.makeRMPred(id[op]), fmt="", ls='none', xerr=50., label=config.tex_labels[op])
                clr = plot[0].get_color()
                #eb = ax0.errorbar(config.x_vals, pb.makeRMPred((-1.0)*(id[op])), ecolor=clr, ls='--', xerr=50., label=config.tex_labels[op]+"-")
                #eb[-1][0].set_linestyle('--')

                pl.legend()

                #ratio plot
                ax2 = fig.add_subplot(spec[1])
                #eb = ax2.errorbar(config.x_vals, np.zeros(len(config.x_vals)), fmt="",  ls='none', xerr=50.0, label=r'')
                neg_rat =  (pb.makeRMPred((-1.0)*(id[op]))) / (sm_pred)
                pos_rat =  (  (pb.makeRMPred(id[op]) - (sm_pred) ) / (sm_pred))
                
                print("pos rat argument = " + str(id[op])  )
                
                #eb = ax2.errorbar(config.x_vals, neg_rat, fmt="",  ls='--', xerr=50.0, label=r'')
                #clr = eb[0].get_color()

                eb = ax2.errorbar(config.x_vals, pos_rat, fmt="", ecolor=clr, ls='none', xerr=50.0, label=r'')

                #eb[-1][0].set_linestyle('--')

                #rel_effect_p = (  (pb.makeRMPred(id[op]) - (sm_pred) ) / (sm_pred))
                #rel_effect_n = (100.0*(  (pb.makeRMPred((-1.0)*(id[op])) - (sm_pred) ) / (sm_pred)) )
                #plot =  ax2.errorbar(config.x_vals, rel_effect_p, fmt="",  ls='none', xerr=50., label=config.tex_labels[op])
                #clr = plot[0].get_color()
                #eb = ax2.errorbar(config.x_vals, rel_effect_n, fmt="",  ls='none', xerr=50., ecolor=clr, label=config.tex_labels[op])
                #eb[-1][0].set_linestyle('--')

            #ax2.fill_between(bin_edges_for_errors, zeros_for_ratio-diffxs_sys_unc_for_ratio, zeros_for_ratio+sys_unc_for_ratio, facecolor='orange', alpha=0.3, edgecolor='none')
            #ax2.fill_between(bin_edges_for_errors, zeros_for_ratio-diffxs_stat_unc_for_ratio, zeros_for_ratio+diffxs_stat_unc_for_ratio, facecolor='#1f77b4', alpha=0.3, edgecolor='none')

                #rel_effects = np.append(rel_effect_n, rel_effect_p)
                ax2.set_xlim(0.0, 500.0)
                #ax2.set_ylim((1.2*np.min(rel_effects)), (1.2*np.max(rel_effects)))

                ax2.xaxis.set_label_coords(0.82, -0.25)
                start, end = ax.get_xlim()
                ax2.yaxis.set_ticks(np.arange(0.75, 3.5, 0.5))

                labelx = ax2.set_xlabel(r'$p^{Z}_{T} \;[GeV]$', fontsize = 15)
                labely = ax2.set_ylabel(r'SMEFT/SM', fontsize = 12)
    
                fig.tight_layout()

                fig.savefig(config.params["config"]["run_name"] + "_results/" + config.params["config"]["run_name"] + "_benchmarks_"+str(op)+".png")
                pl.close()
