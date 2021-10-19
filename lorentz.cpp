#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <TH1F.h>
#include <TF1.h>
#include <math.h>
#include <TMath.h>
#include <TTree.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <THistPainter.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TMinuit.h>
#include <TStyle.h>
#include <TText.h>
#include <TPaveStats.h>
#include <RooLandau.h>
#include "timestamp.cpp"
#include "lorentz_fit.cpp"
using namespace std;
void lorentz()
{
    TCanvas* c1         = new TCanvas("c1");
    TString run_time    = timestamp();
    TFile *myfile       = new TFile("/home/deshan/Documents/Internships/M2/New_env_mu/Q_corr.root", "READ");
    TFile *outfile      = new TFile("L_fit_" + run_time + ".root", "RECREATE");
    TH1F *hist          = (TH1F*)myfile->Get("pe_mm_80");
    hist->Draw("HIST");
    Double_t N, gamma, x0;
    Double_t a0, a1, a2, a3, a4, a5;
    TF1 *fit_lorentz = new TF1 ("fit_lorentz", "RooLandau::RooLandau(x, [1],[2])", 0, 10);
    //TF1 *fit_lorentz = new TF1 ("fit_lorentz", "lorentz_fit", 0, 10, 2);
    //TF1 *fit_lorentz = new TF1 ("fit_lorentz", "gaus + landau", 0, 10);
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);
    c1->Update();
    fit_lorentz->SetNpx(30000);
    //fit_lorentz->SetParameters(1., 0., 1.);
    //fit_lorentz->SetParameters(80000., 0.5);
    //fit_lorentz->FixParameter(0, 10000);
    //fit_lorentz->FixParameter(1, 0.4891);
    //fit_lorentz->FixParameter(2, 0.4891);
    //fit_lorentz->SetParLimits(0, 12000, 14000);
    //fit_lorentz->SetParLimits(1, -0.55, -0.5);
    //fit_lorentz->SetParNames("N", "a");
    //fit_lorentz->SetParNames("N", "#Gamma", "x_{0}");
    for(int i=0; i<10; ++i)
    {
        hist->Fit("fit_lorentz","RE");
        N = fit_lorentz->GetParameter(0);   gamma = fit_lorentz->GetParameter(1);     x0 = fit_lorentz->GetParameter(2); 
        fit_lorentz->SetParameter(0,N);     fit_lorentz->SetParameter(1,gamma);       fit_lorentz->SetParameter(2,x0); 
        cout<<"Chisq/NDOF = "<<fit_lorentz->GetChisquare()/fit_lorentz->GetNDF()<<endl;
        if( ((gMinuit->fCstatu.Contains("CONVERGED")) || (gMinuit->fCstatu.Contains("SUCCESSFUL"))) && (fit_lorentz->GetChisquare()/fit_lorentz->GetNDF()<5) )  {fit_lorentz->SetLineColor(3);  break;}
    } 
    fit_lorentz->Draw("SAME");
    outfile->cd();
    outfile->Write();
}
