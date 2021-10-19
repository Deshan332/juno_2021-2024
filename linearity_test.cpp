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
#include "TGrErr_getter.cpp"
#include "hist_getter.cpp"
#include "timestamp.cpp"
#include "Q_inj_fit.cpp"
using namespace std;
void linearity_test()
{
    TCanvas* c1         = new TCanvas("c1");
    TString run_time    = timestamp();
    TFile *myfile       = new TFile("/home/deshan/Documents/Internships/M2/lin_test_data/2021-04-09/att_c3L_2021-04-09_15_48.root", "READ");
    TFile *outfile      = new TFile(TString("linear_test_") + TString(run_time) + TString(".root"), "RECREATE");
    TTree *mytree       = (TTree*)myfile->Get("TT_ROB_data");
    UShort_t charge[64];
    Int_t att           = 24;
    Int_t count         = 0;
    Int_t nEntries      = 100000;
    Int_t k             = 0;
    Int_t stat1         = 0;
    Int_t stat2         = 0;
    Double_t test1      = 0;
    Double_t test2      = 0;
    Double_t a0, a00, a1, a2, a3, a4, a5, b;
    Double_t ATTS[24]   = {46, 44, 42, 40, 38, 36, 34, 32, 30, 28, 26, 24, 22, 20, 18, 16, 14, 12, 10, 8, 6, 4, 2, 0};
    Int_t channel       = 14;        //channel number starting from 1
    vector<Double_t> mean, mean_err, charge_inj, charge_inj_err;
    Double_t expo,charge_max, alpha=1;  //alpha = mean_err correction for chisq test
    charge_max = 0;
    TString histname, hist_title;
    mytree->SetBranchAddress("charge",  &charge[0]);
    TH1F *myhist[att];
    TGraphErrors *mygraph;
    TH1F *ped1          = hist_getter(TString("/home/deshan/Documents/Internships/M2/lin_test_data/2021-04-09/ped_c3L_2021-04-09_15_41_hist.root"), channel);
    TH1F *ped2          = hist_getter(TString("/home/deshan/Documents/Internships/M2/lin_test_data/2021-04-09/ped_c3L_2021-04-09_16_20_hist.root"), channel);
    Double_t PED_mean   = ( ped1->GetMean() + ped2->GetMean() )/2;
    TString labels[3]   = {TString("Charge Injected (pC)"), TString("Charge collected (ADC)"), TString("Graph of Charge collected vs. Charge injected for FEB no. 808, card 3L, channel ") + TString(Form("%d",channel))};        //array storing labels of graphs  xlabel, ylabel, title
    Double_t ddata[5]   = {1.2, 0., 10.5, 450, 0.};       // array storing double values needed for graphs        //{marker size, x_LL, x_UL, y_max, y_min}
    Int_t idata[4]      = {4, 20, 2, 7};                   // array storing int values needed for graphs           //{marker color, marker size, linewidth, line color}
    Int_t Ntot = mytree->GetEntries();
    if(Ntot/nEntries != att){cout<<"ERROR! Total no. of in the TTree does not match att*nEntries!"<<endl;}
    c1->Print( TString("lin_test_") + TString(run_time) + TString(".pdf[") );  
    for(int i=0; i<att; ++i)
    {
        expo        = -(ATTS[i]/20.0)+1 ;
        histname    = TString("Attenuation = ") +TString(Form("%.0f",ATTS[i])) +TString("dB");
        hist_title  = TString("Attenuation = ") +TString(Form("%.0f",ATTS[i])) + TString("dB, Charge injected = ") + TString(Form("%.2f",pow(10, expo)/0.16) + TString("PE"));
        myhist[i]   = new TH1F(histname, histname, 400, 0, 400);
        for(int j=i*nEntries; j<(i+1)*nEntries; ++j)
        {
            mytree->GetEntry(j);
            myhist[i]->Fill(charge[channel-1]); 
        }
        myhist[i]->SetTitle(hist_title);
        myhist[i]->GetXaxis()->SetTitle("Charge(ADC)");
        myhist[i]->GetYaxis()->SetTitle("Counts");
        myhist[i]->Draw("HIST");
        c1->Print( TString("lin_test_") + TString(run_time) + TString(".pdf") );
        outfile->cd();
        myhist[i]->Write();
        mean.push_back(myhist[i]->GetMean()-PED_mean);
        mean_err.push_back(alpha*myhist[i]->GetMeanError());
        charge_inj.push_back(pow(10, expo));
        charge_inj_err.push_back(0);
        if(myhist[i]->GetMean()-PED_mean > charge_max) {charge_max = myhist[i]->GetMean()-PED_mean;}
    }
    struct graph_data gdata = {charge_inj, charge_inj_err, mean, mean_err};
    mygraph = TGrErr_getter(gdata, ddata, idata, labels); 
    mygraph->Draw("AP");
    //TF1 *fit_pol3 = new TF1("fit_pol3", "[0] + [1]*x + [2]*x^2 + [3]*x^3", 0, 7); 
    TF1 *fit_pol3 = new TF1 ("fit_pol3", Q_inj_fit, 0, 10, 5);
    //TF1 *fit_pol3 = new TF1("fit_pol3", "([0])/( 1 + pow( (x/[1]),[2] ) ) + [3]  ", 0, 10); 
    fit_pol3->SetParameters(0., 0., 50., 1., 2.);
    fit_pol3->FixParameter(0,0);
    fit_pol3->FixParameter(4,2);
    gStyle->SetOptStat(111);
    gStyle->SetOptFit(111);
    gStyle->SetStatX(0.5);
    gStyle->SetStatY(0.9);
    fit_pol3->SetNpx(30000);
    c1->Update();
    fit_pol3->SetParNames("a_{0}", "a_{00}", "a_{1}", "a_{5}", "b"); 
    for(int x=0; x<10; ++x)
    {
        mygraph->Fit("fit_pol3","RE");
        a0 = fit_pol3->GetParameter(0);  a00 = fit_pol3->GetParameter(1);       a1 = fit_pol3->GetParameter(2); a5 = fit_pol3->GetParameter(3); b = fit_pol3->GetParameter(4);
        fit_pol3->SetParameter(0,a0);      fit_pol3->SetParameter(1,a00);          fit_pol3->SetParameter(2,a1);   fit_pol3->SetParameter(3,a5);  fit_pol3->SetParameter(4,b);
        cout<<"Chisq/NDOF = "<<fit_pol3->GetChisquare()/fit_pol3->GetNDF()<<endl;
        if( (gMinuit->fCstatu.Contains("CONVERGED")) || (gMinuit->fCstatu.Contains("SUCCESSFUL"))  )  {fit_pol3->SetLineColor(3);  break;}
    }
    fit_pol3->Draw("SAME");
    c1->Print( TString("lin_test_") + TString(run_time) + TString(".pdf") );
    outfile->cd();
    mygraph->Write();
    c1->Print( TString("lin_test_") + TString(run_time) + TString(".pdf]") );
    outfile->Close();
    delete outfile;
}