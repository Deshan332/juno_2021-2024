#include <iostream>
#include <string>
#include <fstream>
#include <TF1.h>
#include <math.h>
#include <TMath.h>
#include <TTree.h>
#include <TFile.h>
#include <TStyle.h>
#include <TText.h>
#include <TList.h>
using namespace std;

void Qinj_ThirdPE_combine()
{
    Int_t FEB_ID, cat_ID, data_UID, CH;
    Double_t a0, a00, a00_err, a1, a1_err, a2, a2_err, a3, a4, a5, a5_err, ChiSq_Qinj, mu, mu_err, sigma, sig_err, ChiSq_ThPE;
    TFile *Qinj_file    = new TFile("/home/pdeshans/Documents/PhD_files/from_M2_cpy/Qinj_fitdata/After_constrain/test_b2.root", "READ");
    TTree *Qinj_tree    = (TTree*)Qinj_file->Get("TB_lin_par");
    
    Qinj_tree->SetBranchAddress("FEB_ID",       &FEB_ID); 
    Qinj_tree->SetBranchAddress("cat_ID",       &cat_ID); 
    Qinj_tree->SetBranchAddress("data_UID",     &data_UID); 
    Qinj_tree->SetBranchAddress("CH",           &CH); 
    Qinj_tree->SetBranchAddress("a00",          &a00); 
    Qinj_tree->SetBranchAddress("a00_err",      &a00_err);
    Qinj_tree->SetBranchAddress("a1",           &a1); 
    Qinj_tree->SetBranchAddress("a1_err",       &a1_err);
    Qinj_tree->SetBranchAddress("a2",           &a2); 
    Qinj_tree->SetBranchAddress("a3",           &a3); 
    Qinj_tree->SetBranchAddress("a4",           &a4); 
    Qinj_tree->SetBranchAddress("a5",           &a5); 
    Qinj_tree->SetBranchAddress("a5_err",       &a5_err);
    Qinj_tree->SetBranchAddress("ChiSq",        &ChiSq_Qinj);

    TFile *ThPE_file    = new TFile("/home/pdeshans/Documents/PhD_files/from_M2_cpy/ThirdPE_fitdata/corrected_for_variable_plateau_27_01_2022/Scurve_ThirdPE_sum.root", "READ");  
    TTree *ThPE_tree    = (TTree*)ThPE_file->Get("TB_Scurve_par_datetime");
  
    ThPE_tree->SetBranchAddress("mu",           &mu); 
    ThPE_tree->SetBranchAddress("mu_err",       &mu_err);
    ThPE_tree->SetBranchAddress("sigma",        &sigma);  
    ThPE_tree->SetBranchAddress("sig_err",      &sig_err);
    ThPE_tree->SetBranchAddress("ChiSq",        &ChiSq_ThPE); 

    TFile *outfile          = new TFile("Qinj_ThPE_combined.root","recreate");
    TTree *TB_combined      = new TTree("TB_combined","TestBenchDB Qinj data and ThPE data combined");
    TB_combined->Branch("FEB_ID",       &FEB_ID,    "FEB_ID/I");
    TB_combined->Branch("cat_ID",       &cat_ID,    "cat_ID/I");
    TB_combined->Branch("data_UID",     &data_UID,  "data_UID/I"); 
    TB_combined->Branch("CH",           &CH,        "CH/I"); 
    TB_combined->Branch("a00",          &a00,       "a00/D"); 
    TB_combined->Branch("a00_err",      &a00_err,   "a00_err/D");
    TB_combined->Branch("a1",           &a1,        "a1/D"); 
    TB_combined->Branch("a1_err",       &a1_err,    "a1_err/D");
    TB_combined->Branch("a2",           &a2,        "a2/D"); 
    TB_combined->Branch("a3",           &a3,        "a3/D"); 
    TB_combined->Branch("a4",           &a4,        "a4/D"); 
    TB_combined->Branch("a5",           &a5,        "a5/D"); 
    TB_combined->Branch("a5_err",       &a5_err,    "a5_err/D");
    TB_combined->Branch("ChiSq_Qinj",   &ChiSq_Qinj,"ChiSq_Qinj/D"); 
    TB_combined->Branch("mu",           &mu,        "mu/D");
    TB_combined->Branch("mu_err",       &mu_err,    "mu_err/D");  
    TB_combined->Branch("sigma",        &sigma,     "sigma/D");
    TB_combined->Branch("sig_err",      &sig_err,   "sig_err/D");
    TB_combined->Branch("ChiSq_ThPE",   &ChiSq_ThPE,"ChiSq_ThPE/D");

    for(int i=0; i<Qinj_tree->GetEntries(); ++i)
    {
        Qinj_file->cd();
        Qinj_tree->GetEntry(i);
        ThPE_file->cd();
        ThPE_tree->GetEntry(i);
        outfile->cd();
        TB_combined->Fill();
    }

    TB_combined->Write();
    outfile->Close();
    delete outfile;

}