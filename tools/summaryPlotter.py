import matplotlib.pyplot as pl
import numpy as np
import corner



class summaryPlotter:
    def summarise(self, config, pb, samples):
            pl.figure()

            print "pred len " + str (config.params["config"]["data"]["bins"])
            print "sm pred len  " + str (config.predictions['SM'])

            valsp = np.zeros(len(config.coefficients))
            valsn = np.zeros(len(config.coefficients))
            for c in range(0,len(config.coefficients)):
                valsp[c] = 1.0
                valsn[c] = -1.0
                label_stringp = config.coefficients[c] + "=1.0"
                label_stringn = config.coefficients[c] + "=-1.0"
            pl.errorbar(config.x_vals, pb.make_pred(valsp), xerr=0.0, yerr=0.0, label=label_stringp)
            pl.errorbar(config.x_vals, pb.make_pred(valsn), xerr=0.0, yerr=0.0, label=label_stringn)
            pl.errorbar(config.x_vals, config.predictions['SM'], xerr=0.0, yerr=0.0, label='SM')
            pl.errorbar(config.x_vals, config.params["config"]["data"]["central_values"], fmt="o",xerr=0.25, yerr=0.05, label='Data')
            pl.axis([config.x_vals[0]-0.25, config.x_vals[len(config.x_vals)-1]+0.25, 0.0, 7.0])
            pl.legend()

            plotfilename = config.params["config"]["run_name"] +"_predictions.png"

            pl.savefig(plotfilename)

            pl.show()

            #best fit prediction and data
            #TODO

            #corner plot
            #fig = corner.corner(samples, labels=["$ctg$", "$ctw$", "$ctphi$", "$ctb$", "$cphit$","$ctphiQ1$","$ctphiQ3$"],
            #                    truths=[0.0, 0.0, 0.0,0.0, 0.0, 0.0,0.0])

            #fig = corner.corner(samples, labels=["$ctg$", "$ctw$", "$ctphi$", "$ctb$", "$cphit$", "$ctphiQ1$", "$ctphiQ3$"],
            #                    quantiles=[0.05, 0.95],
            #                       range=[1.0,1.0,1.0,1.0,1.0,1.0,1.0], truths=[0.0, 0.0, 0.0,0.0, 0.0, 0.0,0.0],
            #                       show_titles=True, title_kwargs={"fontsize": 18})

            #read in slabels from config

            fig = corner.corner(samples, labels=config.labels,
                                quantiles=[0.05, 0.95],
                                range=[1.0], truths=np.zeros(len(config.coefficients)),
                                show_titles=True, title_kwargs={"fontsize": 18})

            plotfilename = config.params["config"]["run_name"] + ".png"

            fig.savefig(plotfilename)

            fig = pl.figure(figsize=(6, 3.2))
            ax = fig.add_subplot(111)
            ax.set_title('colorMap')
            pl.imshow(config.cov, interpolation='nearest')
            ax.set_aspect('equal')

            cax = fig.add_axes([0.12, 0.1, 0.78, 0.8])
            cax.get_xaxis().set_visible(False)
            cax.get_yaxis().set_visible(False)
            cax.patch.set_alpha(0)
            cax.set_frame_on(False)
            pl.colorbar(orientation='vertical')
            pl.show()



