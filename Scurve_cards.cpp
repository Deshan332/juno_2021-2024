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
#include <TApplication.h>
#include "TGrErr_getter.cpp"
#include "Q_function.cpp"

using namespace std;
int Scurve_cards()
{
    Int_t FEB_idc, FEB_ID=0, cat_idc, cat_ID, test_ID, test_idc, data_uidc, data_UID, CH, status, SID, sidc;
    Double_t mu, mu_err, mu_mean, mu_sig, sigma, sig_err, sig_mean, sig_sig, THRESH, THRESH_mean, ChiSq, ChiSq_mean, ChiSq_sig;
    TString filename = "/home/pdeshans/Documents/PhD_files/from_M2_cpy/Scurve_ThirdPE_sum.root", tree_name = "TB_Scurve_par";
    TFile *myfile       = new TFile(filename); //input TFile
    TTree *TB_pars      = (TTree*)myfile->Get(tree_name); //input TTree

    TB_pars->SetBranchAddress("FEB_ID",       &FEB_ID);
    TB_pars->SetBranchAddress("cat_ID",       &cat_ID);
    TB_pars->SetBranchAddress("data_UID",     &data_UID); 
    TB_pars->SetBranchAddress("test_ID",      &test_ID); 
    TB_pars->SetBranchAddress("CH",           &CH); 
    TB_pars->SetBranchAddress("SID",          &SID); 
    TB_pars->SetBranchAddress("mu",           &mu); 
    TB_pars->SetBranchAddress("mu_err",       &mu_err); 
    TB_pars->SetBranchAddress("sigma",        &sigma); 
    TB_pars->SetBranchAddress("sig_err",      &sig_err); 
    TB_pars->SetBranchAddress("ChiSq",        &ChiSq); 
    TB_pars->SetBranchAddress("THRESH",        &THRESH);
    TB_pars->SetBranchAddress("status",       &status); 
    
    TFile *outfile      = new TFile("TB_cards_ThirdPE.root", "RECREATE"); //output TFile with histos + TTree
    TTree *TB_cards     = new TTree("TB_cards","TestBenchDB ThirdPE data card by card analysis"); //output TTree
 
    TB_cards->Branch("FEB_idc",      &FEB_idc,   "FEB_idc/I"); 
    TB_cards->Branch("cat_idc",      &cat_idc,   "cat_idc/I"); 
    TB_cards->Branch("data_uidc",    &data_uidc, "data_uidc/I"); 
    TB_cards->Branch("test_idc",     &test_idc,  "test_idc/I"); 
    TB_cards->Branch("sidc",         &sidc,      "sidc/I");
    TB_cards->Branch("mu_mean",      &mu_mean,   "mu_mean/D"); 
    TB_cards->Branch("mu_sig",       &mu_sig,    "mu_sig/D");  
    TB_cards->Branch("sig_mean",     &sig_mean,  "sig_mean/D");
    TB_cards->Branch("sig_sig",      &sig_sig,   "sig_sig/D");  
    TB_cards->Branch("THRESH_mean",  &THRESH_mean,  "THRESH_mean/D");
    TB_cards->Branch("ChiSq_mean",   &ChiSq_mean,"ChiSq_mean/D"); 
   
    TString cat_label[8]    = {TString("For TT"), TString("Good spares"), TString("Medium spares"), TString("Bad Spares"), TString("Broken cards"), TString("Non existent cards"), TString("For Telescope"), TString("Other") };
    TH1F *hist_cards[2000][7];
    
    for(int i=0; i<TB_pars->GetEntries(); ++i)
    {
        TB_pars->GetEntry(i);  
        if(CH==0)
        {
            cout<<"Opening new FEB id "<<FEB_ID<<":"<<endl;
            TString hist_tag = TString("card") + TString(Form("%d",FEB_ID));
            hist_cards[FEB_ID][0] = new TH1F(hist_tag + TString("_par_") + TString(Form("%d",0)), hist_tag + TString("_par_") + TString(Form("%d",0)), Int_t(TB_pars->GetMaximum("mu")      - TB_pars->GetMinimum("mu"))    + 20, TB_pars->GetMinimum("mu"),        TB_pars->GetMaximum("mu") );  
            hist_cards[FEB_ID][1] = new TH1F(hist_tag + TString("_par_") + TString(Form("%d",1)), hist_tag + TString("_par_") + TString(Form("%d",1)), Int_t(TB_pars->GetMaximum("sigma")   - TB_pars->GetMinimum("sigma")) + 20, TB_pars->GetMinimum("sigma"),     TB_pars->GetMaximum("sigma") );  
            hist_cards[FEB_ID][2] = new TH1F(hist_tag + TString("_par_") + TString(Form("%d",2)), hist_tag + TString("_par_") + TString(Form("%d",2)), Int_t(TB_pars->GetMaximum("ChiSq")   - TB_pars->GetMinimum("ChiSq")) + 20, TB_pars->GetMinimum("ChiSq"),     TB_pars->GetMaximum("ChiSq") );           
        }
        hist_cards[FEB_ID][0]-> Fill(mu);
        hist_cards[FEB_ID][1]-> Fill(sigma);
        hist_cards[FEB_ID][2]-> Fill(ChiSq);
        if(CH==63)
        {
            mu_mean     = hist_cards[FEB_ID][0]->GetMean();
            mu_sig      = hist_cards[FEB_ID][0]->GetRMS();
            sig_mean    = hist_cards[FEB_ID][1]->GetMean();
            sig_sig     = hist_cards[FEB_ID][1]->GetRMS();
            ChiSq_mean  = hist_cards[FEB_ID][2]->GetMean();
            ChiSq_sig   = hist_cards[FEB_ID][2]->GetRMS();
            FEB_idc = FEB_ID;    cat_idc = cat_ID;    data_uidc = data_UID;     sidc = SID;     test_idc = test_ID; THRESH_mean = THRESH; 
            for(int k=0; k<3; ++k)  hist_cards[FEB_ID][k]->Write();
            TB_cards->Fill();
            cout<<"FEB id "<<FEB_ID<<" results have been written."<<endl;
        }
    }
    TB_cards->Write();
    outfile->Close();
    delete outfile;
    cout<<"TFile closed."<<endl;
    return 0;
}

