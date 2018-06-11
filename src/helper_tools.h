//This header file contains a few simple tools that were needed for the EFT analysis
//but ended up being useful for the charge asymmetry analysis and creation
//of pval tables and plots for the TOP-17-014 paper.
//Therefore they are now decoupled from those class and live here
// as standalone methods that can be intgreated into a similar analysis.
#include <iostream>
#include <sstream>
#include <fstream>
#include "TH1F.h"
#include "TLatex.h"
#include "TMatrixD.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooPlot.h>
#include <RooMCStudy.h>
#include <RooBinning.h>
using namespace std;

TH2D* make_covariance_matrix(std::string, std::string, std::string);
TH1* make_poisson_toy(TH1*,TH1*, int, double);
double calculate_test_statistic(TH1F*, TH1F*, TGraphAsymmErrors*, std::string);
TCanvas* make_results_canvas(TH1F*,vector<double>);
vector<double> get_CL_boundaries(TGraphErrors* gr);

vector<double> get_bestfit_CL_boundaries(TGraphErrors* gr){
    
    double min_relchi =999999.9, best_fit, cl_68_lo, cl_95_lo , cl_68_hi, cl_95_hi;
    bool  hi_68 =true , lo_68= true , lo_95 = true, hi_95 =true;
    double x, y;
    vector<double> boundaries;
    boundaries.clear();
    
    for(Int_t i=0; i < gr->GetN(); i++){
        gr->GetPoint(i,x,y);
        
        if(y < min_relchi){
            min_relchi = y;
            best_fit = x;
        }
        
        if(y <= 3.84 && lo_95){
            cl_95_lo = x;
            lo_95 = false;
            cout <<"found lo 95"<< endl;
        }
        if(y <= 1.0 && lo_68 && !(lo_95)){
            cl_68_lo = x;
            lo_68 = false;
            cout <<"found lo 68"<< endl;
        }
        if(y >= 1.0 && hi_68 && !lo_95 && !lo_68){
            cl_68_hi = x;
            hi_68 = false;
            cout <<"found hi 68" <<endl;
        }
        if(y >= 3.84 && hi_95 && !lo_95 && !lo_68 && !hi_68){
            cl_95_hi = x;
            hi_95 = false;
            cout <<"found hi 95"<< endl;
        }
    }
    
    cout <<"boundaries: "<< cl_95_lo  <<" "<< cl_68_lo << " " << best_fit << "  "<< cl_68_hi <<" "<<cl_95_hi << endl;
    cout <<"uncertainties : "<< fabs(best_fit - cl_95_lo )  <<" "<< fabs(best_fit - cl_68_lo ) << " " << best_fit << "  "<< fabs(best_fit - cl_68_hi )  <<" "<< fabs(best_fit - cl_95_hi )  << endl;

    if ((fabs(best_fit - cl_68_lo ) !=  fabs(best_fit - cl_68_hi )) || (fabs(best_fit - cl_95_lo ) !=  fabs(best_fit - cl_95_hi ))) std::cout <<" boudaries not symmetric around best-fit value -  proceed with caution "<< std::endl;
    
    boundaries.push_back(cl_95_lo);
    boundaries.push_back(cl_68_lo);
    boundaries.push_back(best_fit);
    boundaries.push_back(cl_68_hi);
    boundaries.push_back(cl_95_hi);
    return boundaries;

}

TCanvas* make_results_canvas(TH1F* base_histo, vector<double> boundaries){
    
    base_histo->GetYaxis()->SetRangeUser(0.0, 5.0);
    base_histo->GetXaxis()->SetRangeUser(-0.5, 0.5);
    base_histo->GetYaxis()->SetTitle("#Delta #chi^{2}");
    base_histo->GetXaxis()->SetTitle("C_{tG}/#Lambda^{2} [TeV^{-2}]");
    base_histo->GetYaxis()->SetLabelSize(0.04);
    base_histo->GetXaxis()->SetLabelSize(0.04);
    base_histo->GetXaxis()->SetTitleSize(0.04);
    base_histo->GetYaxis()->SetTitleSize(0.05);
    base_histo->GetXaxis()->SetTitleOffset(1.1);
    base_histo->GetYaxis()->SetTitleOffset(0.85);
    
    base_histo->GetXaxis()->SetLabelOffset(0.01);
    base_histo->GetYaxis()->SetLabelOffset(0.01);
    base_histo->GetYaxis()->SetNdivisions(8);
    base_histo->GetXaxis()->SetNdivisions(12);

    
    TCanvas * all_relscans_c = new TCanvas("","",800,900);
    //all_relscans_c->SetTickx();
    //all_relscans_c->SetTicky();
    
    base_histo->Draw();
    all_relscans_c->SetTopMargin(0.1);
    all_relscans_c->SetBottomMargin(0.15);
    all_relscans_c->SetLeftMargin(0.0935);
    all_relscans_c->SetRightMargin(0.09);

    float H = all_relscans_c->GetWh();
    float W = all_relscans_c->GetWw();
    float l = all_relscans_c->GetLeftMargin();
    float t = all_relscans_c->GetTopMargin();
    float r = all_relscans_c->GetRightMargin();
    float b = all_relscans_c->GetBottomMargin();
    float extraOverCmsTextSize  = 0.8;
    
    TString cmsText, extraText, lumiText;
    cmsText += "CMS";
    extraText += "Preliminary";
    lumiText += "35.9 fb^{-1} (13 TeV)";
    
    TLatex latex;
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextSize(0.48*t);
    latex.SetTextColor(kBlack);
    latex.SetTextFont(61);
    latex.SetTextAlign(31);
    latex.DrawLatex(0.295,0.912,cmsText);
    
    latex.SetTextFont(52);
    latex.SetTextSize(0.48*t*extraOverCmsTextSize);
    latex.DrawLatex(0.495,0.912,extraText);
    
    latex.SetTextFont(42);
    latex.SetTextSize(0.46*t);
    latex.DrawLatex(0.9,0.912,lumiText);
    
    
    return all_relscans_c;

}





double calculate_test_statistic(TH1F* data, TH1F* pred, TGraphAsymmErrors * gr_errors, std::string cov_file){
    
    //std::cout <<"***CALCULATING TEST STATISTIC***   "<<  cov_file  << std::endl;
    
    double test_statistic = 0.0;
    int nbins  = data->GetSize() - 2;
    bool schmitt_fit = true;
    std::vector<double> deltas;
    std::vector<double> errors;
    
    double corr_coff = 1.0;
    bool data_errors_only = true;
    double chisq = 0.0;
    double chi2_bin = 0.0;
    double theory_error =0.0;
    double total_error =0.0;
    
    //Retrieve covariance matrix from text file
    TFile *cov_rootfile = new TFile(cov_file.c_str());
    TH2D * cov = (TH2D*)cov_rootfile->Get("inv_cov");
    
    //TH2D * cov;
    //cov = make_covariance_matrix(cov_file.c_str(),data);
    
    deltas.clear();
    errors.clear();
    
    double cumul_chi2 = 0;
    //std::cout <<" "<< endl;


    for (int i=1;i <=  data->GetSize() -2; i++){
        
        double data_bin = data->GetBinContent(i);
        double pred_bin = pred->GetBinContent(i);
        double delta = data_bin - pred_bin;
        deltas.push_back(delta);
        
        cumul_chi2 = cumul_chi2 + (delta*delta) / ( data->GetBinError(i) * data->GetBinError(i));
        
        if (data_errors_only){
            
   // if (cov_file == "files/Nov1/particle/normalised/covariance/HypJetMultpt30_totCovEnvXSMtrxFile.root") std::cout <<"data " <<data_bin <<", pred_bin "<< pred_bin <<" delta " << delta<< " delta/error = "<< delta/data->GetBinError(i) << " delta_sq " << delta*delta <<", data error squared  = "<< data->GetBinError(i)* data->GetBinError(i) <<  " chi2 in bin =  " <<  (delta*delta)/(data->GetBinError(i)*data->GetBinError(i))<< " cumul. chi2 = " << cumul_chi2 <<"\n";
            errors.push_back(data->GetBinError(i));
            // errors.push_back(pred_bin);
            
        }else{
            //implmement quadrature addition of data and pred errors
            
            if (data_bin > pred_bin ){
                theory_error = gr_errors->GetErrorYhigh(i-1);//accroding to aMC@NLO manual, the scale uncertainty of the NLO predicitons is 10%
            }else if (data_bin < pred_bin ){
                theory_error = gr_errors->GetErrorYlow(i-1);
            }
            total_error =  pow(  pow( (theory_error +  data->GetBinError(i) ) ,2.0), 0.5);
            errors.push_back( total_error);
        }
    }
    

    //TMatrixD mat_inv_cov =  TMatrixD(matrix.size(), matrix.size());
    //TMatrixD mat_deltas =  TMatrixD(matrix.size(), matrix.size());
    
    TMatrixD mat_cov = TMatrixD(nbins,nbins);
    TMatrixD mat_inv_cov(nbins,nbins);
    TMatrixD mat_deltas_tp(1, nbins);
    TMatrixD mat_deltas(nbins, 1);
    TMatrixD mat_chi(1, 1);

    
    for (int i=1;i<=nbins; i++) {
      //  std::cout <<"  "<< std::endl;
        mat_deltas[i-1][0] = deltas[i-1];
        mat_deltas_tp[0][i-1] = deltas[i-1];
        
        for (int j=1;j<=nbins; j++){
            
            if (i==j){
                mat_cov[i-1][j-1] = (errors[i-1]*errors[j-1]);
            }
            else if ((i<4) &&(j<4) && (fabs(i-j) ==1.0)){
                mat_cov[i-1][j-1] = (errors[i-1]*errors[j-1])*(0.5);
            }
            else{
                mat_cov[i-1][j-1] = 0.0;
            }
            
            double corr_coff = cov->GetBinContent(i,j);
            if (schmitt_fit){
                mat_inv_cov[i-1][j-1] = corr_coff;
                chi2_bin = deltas[i-1]*deltas[j-1]*corr_coff;
                chisq += chi2_bin;
//if (cov_file == "files/Nov1/parton/normalised/covariance/HypAntiToppT_totCovEnvXSMtrxFile.root")std::cout <<"data "<<data->GetBinContent(i) <<" data error "<<errors[i-1] <<"  delta/unc  "<< deltas[i-1]/errors[i-1]  <<"  data error sq "<< errors[i-1]*errors[j-1]  <<" delta i "<< deltas[i-1]<< " delta j " <<  deltas[j-1] <<" corr coff "<< corr_coff  <<" chi2 bin "<<chi2_bin <<" chi2 "<< chisq << "\n";
                }
            else {
                if(i == j){
                    //std::cout<<"Calc chi2 " << "bin  "<< i << " delta =  " << deltas[i] <<  ", error  = "<< errors[i]  <<"\n";
                    double delta_sq = deltas[i-1]*deltas[i-1];
                    //double chi2_bin = delta_sq/errors[i-1];
                    //                    double chi2_bin = delta_sq/(errors[i-1]*errors[i-1]);
                    chi2_bin = delta_sq/(errors[i-1]*errors[i-1]);
                    //  std::cout <<"chi2 "<< chisq << "\n";
                    chisq += chi2_bin;
                }
            }
        }
    }
    
    mat_cov.Invert();
    
    if (cov_file == "files/Nov1/particle/normalised/covariance/HypJetMultpt30_totCovEnvXSMtrxFile.root"){
        
    mat_chi = mat_deltas_tp * mat_cov * mat_deltas;
    //std::cout <<"chi2_xcheck "<< mat_deltas_tp.Print() << endl;
   // mat_deltas_tp.Print();
   // mat_inv_cov.Print();
   // mat_deltas.Print();
    mat_chi.Print();
    mat_cov.Print();

    }
    
    cov_rootfile->Close();
        
    return chisq;
}



TH2D* make_covariance_matrix(std::string filename, std::string data_filename, std::string varname){
    
    //std::cout <<" make_covariance_matrix "<< filename  <<std::endl;

    //first fill vectors from text file
    std::ifstream readFile(filename.c_str());
    std::vector<std::vector<std::string>> matrix;
    std::vector<std::string> tokens;
    std::string line;
    std::string delimiter = " ";
    
    while(getline(readFile,line)){
        tokens.clear();
        std::istringstream iss(line);
        size_t pos = 0;
        std::string token;
    
        while ((pos = line.find(delimiter)) != std::string::npos){
            token = line.substr(0, pos);
            //std::cout <<" mini-loop token "<< token << std::endl;
            line.erase(0, pos + delimiter.length());
            tokens.push_back(token);
        }
        //tokens.push_back(line);
        matrix.push_back(tokens);
    }
    readFile.close();
    
    //std::cout <<" rows "<< matrix.size() <<"  columns  "<< matrix[0].size() <<std::endl;

    TFile * data_file = new TFile(data_filename.c_str());
    TMatrixD mat =  TMatrixD(matrix.size(), matrix.size());
    TH1F * h_template = (TH1F*)data_file->Get("mc");
    TGraphAsymmErrors * data = (TGraphAsymmErrors*)data_file->Get("data");
    
    int nbins = h_template->GetNbinsX();
    double lower_range = h_template->GetBinLowEdge(1);
    double upper_range = h_template->GetBinLowEdge(nbins) + h_template->GetBinWidth(nbins);

    TH2D* cov = new TH2D("cov","cov", matrix.size(), 1.0, matrix.size() +1.0 , matrix.size(), 1.0, matrix.size() + 1.0);
    TH2D* inv_cov = new TH2D("inv_cov","inv_cov", matrix.size(), 1.0, matrix.size() +1.0, matrix.size(), 1.0, matrix.size() + 1.0);
    
    cov->GetXaxis()->SetNdivisions(matrix.size());
    cov->GetYaxis()->SetNdivisions(matrix.size());
    inv_cov->GetXaxis()->SetNdivisions(matrix.size());
    inv_cov->GetYaxis()->SetNdivisions(matrix.size());
    
    cov->GetXaxis()->CenterLabels();
    cov->GetYaxis()->CenterLabels();
    inv_cov->GetXaxis()->CenterLabels();
    inv_cov->GetYaxis()->CenterLabels();
    
    
    double point_i_x,point_i_y, point_j_x, point_j_y;
    
    for (int x = 0; x < matrix.size(); x++){
        data->GetPoint(x, point_i_x, point_i_y);
        
        for (int y = 0; y < matrix[x].size(); y++){
            //std::cout <<" xm "<< matrix[x][y] <<std::endl;
            std::string::size_type sz;
            std::string cov_elem_str = matrix[x][y];
            
            data->GetPoint(y, point_j_x, point_j_y);
            double cov_elem = std::stod(cov_elem_str, &sz);
            //double cov_elem = cov_elem_rel*point_i_y*point_j_y;
            //cov_elem = matrix[x][y];
            cov->SetBinContent(x+1, y+1, cov_elem);
            mat[x][y] = cov_elem;
           //
            //std::cout <<" x = "<< x <<"  y "<< y <<"  "<< cov_elem_str <<std::endl;
        }
    }
    
    //invert covariance matrix and write
    Double_t det2;
    mat.Invert(&det2);
    
    for (int x = 0; x < matrix.size(); x++){
        for (int y = 0; y < matrix[x].size(); y++){
            inv_cov->SetBinContent(x+1, y+1, mat[x][y]);
        }
    }
    
    std::string rootfile_name = filename.substr(0, filename.length() - 3) + "root";
    std::string pdffile_name = filename.substr(0, filename.length() - 3) + "pdf";
    std::string axis_title = varname + "   bin number";

    cov->SetTitle("");
    cov->SetXTitle(axis_title.c_str());
    cov->SetYTitle(axis_title.c_str());
    
    TCanvas * c = new TCanvas("","",800,800);
    gStyle->SetOptStat(00000);
    cov->Draw("COLZTEXT");

    c->SaveAs(pdffile_name.c_str());

    TFile * file = new TFile(rootfile_name.c_str(), "RECREATE");
    cov->Write();
    inv_cov->Write();
    mat.Write();
    file->Close();
    data_file->Close();

    
    return inv_cov;
}


TH1* make_poisson_toy(TH1* hh, TH1* data, int nevents, double data_integral){
    
    //cout <<"making poisson toy " << endl;
    
    std::string opt = "manual";
    
    //double n_toys = 523850.625;
    
    //This function make a Poisson toy for a given prediction at particle level
    // corresponding to the same number of events as observed in data.
    
    //Histos to contain raw distributions for prediction and data.
    
    // TH1F* hh_t = new TH1F("","", hh->GetNbinsX() , 0.0, 3.15);
    // TH1F* hh_data = new TH1F("","", hh->GetNbinsX() , 0.0, 3.15);
    
    //TH1F* hh_t = (TH1F*)hh->Clone();
    //TH1F* hh_data = (TH1F*)hh->Clone();
    
    TH1 * hh_f;
    hh_f =(TH1F*)hh->Clone();
    hh_f->Reset();
    
    
    //hh_t->Reset();
    //hh_data->Reset();
    double bin_height_pred;
    double bin_height_data;
    TRandom3 * gRand = new TRandom3();

    if (opt == "manual"){
        for (int bin = 1; bin <= hh->GetNbinsX(); bin++){
            double mean = hh->GetBinContent(bin);
            double width = data->GetBinError(bin);
            double toy_bin_height = gRand->Gaus(mean,width);
            double toy_bin_error = width;
            //  cout <<"TOY =  "<< "mean " << mean <<" width  "<<width<< " toy height " <<  toy_bin_height <<endl;
            hh_f->SetBinContent(bin, toy_bin_height);
            hh_f->SetBinError(bin, toy_bin_error);
            //hh_t->SetBinContent(bin, bin_height_pred);
            //hh_data->SetBinContent(bin, bin_height_data);
        }
    }
    else if (opt=="roofit"){
        //Multiply back by bin widths to get cross sections as bin heights
        for (int bin = 1; bin <= hh->GetNbinsX(); bin++){
            if (bin == hh->GetNbinsX() ){
                bin_height_pred = (0.24)*hh->GetBinContent(bin);
                bin_height_data = (0.24)*data->GetBinContent(bin);
            }else{
                bin_height_pred = (0.25)*hh->GetBinContent(bin);
                bin_height_data = (0.25)*data->GetBinContent(bin);
            }
            //hh_t->SetBinContent(bin, bin_height_pred);
            //     hh_data->SetBinContent(bin, bin_height_data);
        }
        
        
        //now scale by integrated lumi
        //double scaling = hh->Integral() / hh_t->Integral();
        
        //    double int_lumi = 37500000.0;
        //double int_lumi = 3750.0;
        
        double int_lumi = 6000.0;
        
        //hh_t->Scale(int_lumi);
        
        //TH1F* h_forPDF = new TH1F("","", hh->GetNbinsX() , -1.0,1.0);
        TH1F* h_forPDF = (TH1F*)hh->Clone();
        h_forPDF->Reset();
        
        
        //fill histo with raw event counts for this prediction
        // for (int bin = 1; bin <= hh->GetNbinsX(); bin++){
        //     h_forPDF->SetBinContent(bin, int_lumi *  hh_t->GetBinContent(bin));
        // }
        
        
        //some neccessary Roofit declarations
        RooRealVar x("x","x",0,3.15);
        x.setBins(hh->GetNbinsX());
        
        
        //convert raw prediction into pdf
        RooDataHist dh("hh_t","hh_t",x,RooFit::Import(*h_forPDF)) ;
        RooHistPdf histpdf1("histpdf1","histpdf1",x,dh,0) ;
        
        //generate toy event count
        TRandom * r = new TRandom();
        double n_expected = r->Poisson(h_forPDF->Integral());
        
        //generate binned toys according the PDF with Poisson errors
        RooDataHist* toy_data = histpdf1.generateBinned(x,n_expected);
        
        //convert it back to a TH1
        hh_f = toy_data->createHistogram("toy", x);
        
        
        double bbin;
        double bbin_error;
        
        //convert back to 'cross section result' histo
        //divide by bin width
        for (int bin = 1; bin <= hh->GetNbinsX(); bin++){
            
            if (bin == hh->GetNbinsX() ){
                bbin = hh_f->GetBinContent(bin)/(0.24);
                // bbin_error = hh_f->GetBinError(bin)/(0.24);
                bbin_error = bbin *(0.01);
                
            }else{
                bbin = hh_f->GetBinContent(bin)/(0.25);
                // bbin_error = hh_f->GetBinError(bin)/(0.25);
                bbin_error = bbin *(0.01);
                
            }
            
            hh_f->SetBinContent(bin,bbin);
            hh_f->SetBinError(bin,bbin_error);
        }
        
        //divide by lumi
        hh_f->Scale(1.0/int_lumi);
        
        
        
        // double scaling_2 = data_integral/hh_f->Integral();
        double scaling_2 = hh->Integral()/hh_f->Integral();
        
        //hh_f->Scale(scaling_2);
        
    }
    //cout <<" toy integral  = "<<  hh_f->Integral()<< endl;
    
    return hh_f;
}


