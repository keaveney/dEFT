import matplotlib.pyplot as pl
import matplotlib.image as image

import numpy as np
import corner
import os
import json

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


class summaryPlotter:
    def summarise(self, config, pb, sampler, samples):
        
            mcmc_params = np.mean(sampler.flatchain,axis=0)
            mcmc_params_cov = np.cov(np.transpose(sampler.flatchain))
            
            cmd = "mkdir -p " + str(config.params["config"]["run_name"]) + "_results"
            os.system(cmd)
        
            #print "pred len " + str (config.params["config"]["data"]["bins"])
            #print "sm pred len  " + str (config.predictions['SM'])
            data_label = "Data" + " (" + config.run_name + ")"
            max_val = (1.5)*(max(config.params["config"]["data"]["central_values"]))
            min_val = (0.0)*(min(config.params["config"]["data"]["central_values"]))
            xlabel = config.observable
            ylabel = "d$\sigma_{tt}$/dX"
                
            #zero step, compare explicit predictions with output of the morphing model
            #at a few test points to validate the morhing model.
            #This expects an extra json file containing the explicit predictions to
            #be present with the same name and a '_test' tag.
            

            #first make plots of basis predictions versus data
            coefficients = list(config.coefficients)
            for c in range(0,len(coefficients)):
                valsp = np.zeros(len(coefficients))
                valsn = np.zeros(len(coefficients))
                valsp[c] = 1.0
                valsn[c] = -1.0
                label_stringp = coefficients[c] + "/$\Lambda^{2}$" + " = 1.0 " + "[TeV$^{-2}$]"
                label_stringn = coefficients[c] + "/$\Lambda^{2}$" + " = -1.0 " + "[TeV$^{-2}$]"
                pl.figure()
                pl.errorbar(config.x_vals, pb.makeRMPred(valsp), xerr=0.0, yerr=0.0, label=label_stringp)
                pl.errorbar(config.x_vals, pb.makeRMPred(valsn), xerr=0.0, yerr=0.0, label=label_stringn)
                pl.errorbar(config.x_vals, config.params["config"]["model"]["predictions"][0], xerr=0.0, yerr=0.0, label='SM')
                pl.errorbar(config.x_vals, config.params["config"]["data"]["central_values"], fmt="o",xerr=0.25, yerr=0.05, label=data_label)
                pl.axis([config.x_vals[0]-0.25, config.x_vals[len(config.x_vals)-1]+0.25, min_val, max_val])
                ax = pl.gca()
                labely = ax.set_xlabel(xlabel, fontsize = 18)
                ax.xaxis.set_label_coords(0.85, -0.065)
                labely = ax.set_ylabel(ylabel, fontsize = 18)
                ax.yaxis.set_label_coords(-0.037, 0.83)
                pl.legend(loc=2)
                plotfilename = str(config.params["config"]["run_name"]) + "_results/" + str(config.params["config"]["run_name"]) + "_" + coefficients[c] + "_predictions.png"
                pl.savefig(plotfilename)
                pl.close()

            
            #second make plots of best fit prediction versus data
            pl.figure()
            label_string_bestfit = "best-fit: ("
            for c in range(0, len(config.coefficients)):
                if (c == (len(config.coefficients) - 1)):
                    label_string_bestfit = label_string_bestfit + coefficients[c] + " = " + '%.3f' % mcmc_params[c] + ")"
                else:
                    label_string_bestfit = label_string_bestfit + coefficients[c] + " = " + '%.3f' % mcmc_params[c] + ", "

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

            #third  make "corner" plots
            labels = []
            ranges  = []
            for c in config.params["config"]["model"]["prior_limits"].keys():
                label = "$" + c + "$"
                labels.append(label)
                #ranges.append(1.0)
                ranges.append(config.params["config"]["model"]["prior_limits"][c])

            
            #print "samples " + str(samples)
            #df = pd.DataFrame.from_records(sampler.get_blobs(flat=True, discard=100, thin=30))
    
    #        fig = corner.corner(x, quantiles=(0.16, 0.84), levels=(1-np.exp(-0.5),))

            texLabels = ["$c_{HWB}$","$c_{W}$","$c_{HB}$","$c_{Ht}$","$c_{tZ}$","$c_{t W}$","$c_{tG}$"]
            
            #fig = corner.corner(samples, plot_contours=False, labels=texLabels,
            #                    range=ranges, truths=np.zeros(len(labels)),
            #                    show_titles=True, title_kwargs={"fontsize": 18})

            fig = corner.corner(samples, labels=texLabels,
                                range=ranges,
                                quantiles=[0.16, 0.5, 0.84],
                                truths=np.zeros(len(labels)),
                                show_titles=True, title_kwargs={"fontsize": 18})

            plotfilename = str(config.params["config"]["run_name"]) + "_results/" + str(config.params["config"]["run_name"]) + ".png"
            logo = image.imread('logo/dEFT_logo.png')

            ax0 = fig.add_subplot(999)
            ax0.axis('off')
            img = ax0.imshow(logo)

            fig.savefig(plotfilename)
            pl.close()

            fig = pl.figure(figsize=(6, 3.2))
            ax = fig.add_subplot(111)
            ax.set_title('colorMap')
            #pl.imshow(config.cov, interpolation='nearest')
            ax.set_aspect('equal')

            cax = fig.add_axes([0.12, 0.1, 0.78, 0.8])
            cax.get_xaxis().set_visible(False)
            cax.get_yaxis().set_visible(False)
            cax.patch.set_alpha(0)

            cax.set_frame_on(False)
            #pl.colorbar(orientation='vertical')
            #pl.show()

            bestFits = []
            marginUncsUp = []
            marginUncsDown = []
            x = []
    
            for i in range(len(labels)):
                mcmc = np.percentile(samples[:, i], [16, 50, 84])
                q = np.diff(mcmc)
                txt = "\mathrm{{{3}}} = {0:.3f}_{{-{1:.3f}}}^{{{2:.3f}}}"
                txt = txt.format(mcmc[1], q[0], q[1], labels[i])
                x.append(i)
                bestFits.append(mcmc[1])
                marginUncsUp.append(q[0])
                marginUncsDown.append(q[0])
                
            pl.figure()
            fig = pl.errorbar(x, bestFits, yerr=[marginUncsDown, marginUncsUp], fmt='o')
            pl.xticks(range(len(labels)), labels, rotation='45')
            pl.savefig(config.params["config"]["run_name"] + "_results/" + config.params["config"]["run_name"] + "_bestFits.png")
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
            # best-fit prediction and CI compared to data
            uncX = np.absolute((np.diff(config.bins)))/2.0
            dataUncY = np.sqrt(np.diag(config.cov))
            pl.errorbar(config.x_vals, config.params["config"]["data"]["central_values"], fmt="o",xerr=uncX, yerr=dataUncY, label=data_label)
            pl.errorbar(config.x_vals, pb.makeRMPred(bestFits), fmt="m", xerr=0.0, yerr=0.0, label=label_string_bestfit)
            
            inds = np.random.randint(len(samples), size=100)
            for ind in inds:
                sample = samples[ind]
                pl.plot(config.x_vals, pb.makeRMPred(sample), "C1", alpha=0.1)
                #pl.plot(config.x_vals, np.dot(np.vander(config.x_vals, 2), sample[:2]), "C1", alpha=0.1)
                #pl.plot(config.x_vals, sample[:2], "C1", alpha=0.1)

            #pl.axis([config.x_vals[0]-0.25, config.x_vals[len(config.x_vals)-1]+0.25, min_val*0.5, max_val*1.5])
            #pl.axis([config.x_vals[0]-0.25, config.x_vals[len(config.x_vals)-1]+0.25, 0.0, 10.0])

            #pl.plot(x0, m_true * x0 + b_true, "k", label="truth")
            pl.savefig(config.params["config"]["run_name"] + "_results/" + config.params["config"]["run_name"] + "_postFit" + "_predictions.png")
            pl.close()

            pl.figure()
            fig, axes = pl.subplots(len(labels), figsize=(10, 7), sharex=True)
            samples = sampler.get_chain()
            for i in range(len(labels)):
                ax = axes[i]
                ax.plot(samples[:, :, i], "k", alpha=0.3)
                ax.set_xlim(0, len(samples))
                ax.set_ylabel(labels[i])
                ax.yaxis.set_label_coords(-0.1, 0.5)

            axes[-1].set_xlabel("step number");

            pl.savefig(config.params["config"]["run_name"] + "_results/" + config.params["config"]["run_name"] + "_walkerPaths.png")
            pl.close()
            

            ######################################################
            ###############   PLOT RESULTS   ####################
            ######################################################
            print("best fit vals = " + str(mcmc_params))
            print("uncertainties = " + str(np.sqrt(np.diag(mcmc_params_cov))))
            print("cov mat = " + str(mcmc_params_cov))

            pl.figure()

            pl.matshow(mcmc_params_cov, cmap=pl.cm.Blues)
            pl.savefig(config.params["config"]["run_name"] + "_results/" + config.params["config"]["run_name"] + "mcmc_params_cov.png")
            pl.close()



