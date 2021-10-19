#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <TH1F.h>
#include <TF1.h>
#include <math.h>
#include <TMath.h>
#include<TGraph.h>
#include <TTree.h>
#include<TGraphErrors.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TMinuit.h>
#include <TStyle.h>
#include <TLatex.h>
using namespace std;

void mu_hist()
{
    TFile *myfile = new TFile("Cards_Analysis_test.root", "READ");
    TTree *mytree = (TTree*)myfile->Get("par_tree");
    Double_t mu, mu_err, Chisq;
    Int_t status, bin, NB=7000; 
    mytree->SetBranchAddress("mu",        &mu);
    mytree->SetBranchAddress("mu_err",    &mu_err);
    mytree->SetBranchAddress("Chisq",     &Chisq);
    mytree->SetBranchAddress("status",    &status);
    Int_t Ntot = mytree->GetEntries();
    Int_t counts[NB];
    TH1F *myhist = new TH1F("myhist", "Graph of #delta#mu/#mu vs. mu",NB, 0, 7);
    for(Int_t i=0; i<Ntot; ++i)
    {
        mytree->GetEntry(i);
        if( (Chisq<5) && (status==1) && (mu_err/mu<0.1))
        {
            bin = myhist->Fill(mu, mu_err/mu);
            if(mu_err/mu != 0) {counts[bin] += 1;}
        }
    }
    Int_t Nbins = myhist->GetXaxis()->GetNbins();
    for(Int_t j=0; j<Nbins; ++j)
    {
        myhist->SetBinContent(j, myhist->GetBinContent(j)/counts[j] );
        if(counts[j] == 0)  {myhist->SetBinContent(j, 0);}
    }
    myhist->Rebin(2);
    myhist->Draw("HIST");
}