from ROOT import *
import ROOT
f = TFile("root_inputs/nlo_ctG_p2_lhe_nom_pt_nom_xxx.root")
h = f.Get("TTbarSpinDensityMatrix/snake_spinCorr_part_cut_0")

sigma_CtG_0_nlo = 687.14
sigma_CtG_neg2_nlo = 404.48
sigma_CtG_pos2_nlo = 1270.204
acceptance = 1.0
branching_ratio = 0.046
k_factor = 1.21
n_bases =22.0

fid_xsec = (sigma_CtG_pos2_nlo*k_factor*acceptance*n_bases)
integral = h.Integral()
h.Scale(fid_xsec/integral)

print "[",
for bin in range(0,h.GetNbinsX()):
    #print str(    h.GetBinCenter(bin+1)  -  ( (h.GetBinWidth(bin+1))/2.0)  ) + ",",
    print str(    h.GetBinContent(bin+1)    ) + ",",
print "]"


f2 = TFile("root_inputs/unfolded_data_181109.root")

h2 = f2.Get("TTbarSpinDensityMatrix/snake_spinCorr_absolute_data_b")

#print "[",
#for bin in range(0,h2.GetNbinsX()):
#    print str(    h2.GetBinContent(bin+1)  / (h2.GetBinWidth(bin+1))  ) + ",",
#print "]"

f3 = TFile("root_inputs/covmat_181109_Systematics_AllVars_final.root")
cov = f3.Get("TotalStatSystCovMatrix_AllVar_rebinnedB")

#print "[",
#for bini in range(0,cov.GetNbinsX()):
# print "[",
#for binj in range(0,cov.GetNbinsY()):
#    cov_ij = cov.GetBinContent(bini+1, binj+1)
#    print str(cov_ij) + ",",
#print "]"
#print "]"


"""
    #TOP-17-014
f = TFile("/Users/keaveney/EFT_Fitter/EFT_Fitter/files/CtG_2_nominal_v10.root")
h = f.Get("CMS_dilepton_diff/ll_delphi_abs")

#for bin in range(0,h.GetNbinsX()):
#h.SetBinContent(bin+1, (h.GetBinContent(bin+1)/h.GetBinWidth(bin+1)))

sigma_CtG_0_nlo = 687.14
sigma_CtG_neg2_nlo = 404.48
sigma_CtG_pos2_nlo = 1270.204
acceptance = 0.28371
branching_ratio = 0.046
k_factor = 1.21

fid_xsec = (sigma_CtG_pos2_nlo*k_factor*branching_ratio*acceptance)
integral = h.Integral()
h.Scale(fid_xsec/integral)

print "[",
for bin in range(0,h.GetNbinsX()):
    print str(    h.GetBinContent(bin+1)  / (h.GetBinWidth(bin+1))  ) + ",",
print "]"

f2 = TFile("/Users/keaveney/EFT_Fitter/EFT_Fitter/files/Nov1/particle/absolute/results/DiffXS_HypLLBarDPhi_source.root")

g = f2.Get("data")
x, y = ROOT.Double(0), ROOT.Double(0)

print "[",
for bin in range(0,h.GetNbinsX()):
    g.GetPoint(bin, x, y)
    print str(y) + ",",
print "]"

f3 = TFile("/Users/keaveney/EFT_Fitter/EFT_Fitter/files/Nov1/particle/absolute/covariance/HypLLBarDPhi_totCovEnvXSMtrxFile.root")

cov = f3.Get("inv_cov")

print "[",
for bini in range(0,cov.GetNbinsX()):
    print "[",
    for binj in range(0,cov.GetNbinsY()):
        cov_ij = cov.GetBinContent(bini+1, binj+1)
        print str(cov_ij) + ",",
    print "]"
print "]"

"""
