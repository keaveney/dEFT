import matplotlib.pyplot as pl
import matplotlib.image as image

import numpy as np
import corner

class summaryPlotter:
    def summarise(self, config, pb, sampler, samples, mcmc_params):
            #print "pred len " + str (config.params["config"]["data"]["bins"])
            #print "sm pred len  " + str (config.predictions['SM'])
            data_label = "Data" + " (" + config.run_name + ")"
            max_val = (1.5)*(max(config.params["config"]["data"]["central_values"]))
            min_val = (0.0)*(min(config.params["config"]["data"]["central_values"]))
            xlabel = config.observable
            ylabel = "d$\sigma_{tt}$/dX"

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
                pl.errorbar(config.x_vals, pb.makePred(valsp), xerr=0.0, yerr=0.0, label=label_stringp)
                pl.errorbar(config.x_vals, pb.makePred(valsn), xerr=0.0, yerr=0.0, label=label_stringn)
                pl.errorbar(config.x_vals, config.params["config"]["model"]["predictions"][0], xerr=0.0, yerr=0.0, label='SM')
                pl.errorbar(config.x_vals, config.params["config"]["data"]["central_values"], fmt="o",xerr=0.25, yerr=0.05, label=data_label)
                pl.axis([config.x_vals[0]-0.25, config.x_vals[len(config.x_vals)-1]+0.25, min_val, max_val])
                ax = pl.gca()
                labely = ax.set_xlabel(xlabel, fontsize = 18)
                ax.xaxis.set_label_coords(0.85, -0.065)
                labely = ax.set_ylabel(ylabel, fontsize = 18)
                ax.yaxis.set_label_coords(-0.037, 0.83)
                pl.legend(loc=2)
                plotfilename = str(config.params["config"]["run_name"] + "_" + coefficients[c] + "_predictions.png")
                pl.savefig(plotfilename)
            
            #second make plots of best fit prediction versus data
            pl.figure()
            label_string_bestfit = "best-fit: ("
            for c in range(0, len(config.coefficients)):
                if (c == (len(config.coefficients) - 1)):
                    label_string_bestfit = label_string_bestfit + coefficients[c] + " = " + '%.3f' % mcmc_params[c] + ")"
                else:
                    label_string_bestfit = label_string_bestfit + coefficients[c] + " = " + '%.3f' % mcmc_params[c] + ", "

            pl.errorbar(config.x_vals, pb.makePred(mcmc_params), fmt="m", xerr=0.0, yerr=0.0, label=label_string_bestfit)
            pl.errorbar(config.x_vals, config.params["config"]["data"]["central_values"], fmt="o",xerr=0.25, yerr=0.05, label=data_label)
            pl.axis([config.x_vals[0]-0.25, config.x_vals[len(config.x_vals)-1]+0.25, min_val, max_val])
            ax = pl.gca()
            labely = ax.set_xlabel(xlabel, fontsize = 18)
            labely = ax.set_ylabel(ylabel, fontsize = 18)
            ax.xaxis.set_label_coords(0.85, -0.065)
            ax.yaxis.set_label_coords(-0.037, 0.83)
            pl.legend(loc=2)
            pl.savefig(config.params["config"]["run_name"] + "_bestfit" + "_predictions.png")
            #third  make "corner" plots
            labels = []
            ranges  = []
            for c in config.params["config"]["model"]["prior_limits"].keys():
                label = "$" + c + "$"
                labels.append(label)
                ranges.append(1.0)
            
            #print "samples " + str(samples)
            #df = pd.DataFrame.from_records(sampler.get_blobs(flat=True, discard=100, thin=30))
            
            fig = corner.corner(samples, labels=labels,
                                quantiles=[0.32, 0.68, 0.05, 0.95],
                                range=ranges, truths=np.zeros(len(labels)),
                                show_titles=True, title_kwargs={"fontsize": 18})

            plotfilename = config.params["config"]["run_name"] + ".png"
            logo = image.imread('logo/dEFT_logo.png')

            ax0 = fig.add_subplot(999)
            ax0.axis('off')
            img = ax0.imshow(logo)

            fig.savefig(plotfilename)

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
            



