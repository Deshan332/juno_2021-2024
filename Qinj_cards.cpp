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
#include "Q_inj_fit.cpp"

using namespace std;
int Qinj_cards()
{
    TCanvas* c1         = new TCanvas("c1");
    Int_t FEB_idc, FEB_ID=0, cat_idc, cat_ID, data_uidc, data_UID, CH;
    Double_t b, a0, a00, a00_err, a00_mean, a00_sig, a1, a1_err, a1_mean, a1_sig, a2, a2_err, a2_mean, a2_sig, a3, a3_err, a3_mean, a3_sig, a4, a4_err, a4_mean, a4_sig, a5, a5_err, a5_mean, a5_sig, ChiSq, ChiSq_mean, ChiSq_sig;
    TString filename = "/home/pdeshans/Documents/PhD_files/from_M2_cpy/Qinj_fitdata/test_b2.root", tree_name = "TB_lin_par";
    TFile *myfile       = new TFile(filename, "READ"); //input TFile
    TTree *TB_pars      = (TTree*)myfile->Get(tree_name); //input TTree
    TB_pars->SetBranchAddress("FEB_ID",       &FEB_ID);
    TB_pars->SetBranchAddress("cat_ID",       &cat_ID);
    TB_pars->SetBranchAddress("data_UID",     &data_UID); 
    TB_pars->SetBranchAddress("CH",           &CH); 
    TB_pars->SetBranchAddress("a0",           &a0); 
    TB_pars->SetBranchAddress("a00",          &a00); //j=0
    TB_pars->SetBranchAddress("a00_err",      &a00_err); 
    TB_pars->SetBranchAddress("a1",           &a1); //j=1
    TB_pars->SetBranchAddress("a1_err",       &a1_err); 
    TB_pars->SetBranchAddress("a2",           &a2); //j=2
    TB_pars->SetBranchAddress("a3",           &a3); //j=3
    TB_pars->SetBranchAddress("a4",           &a4);//j=4
    TB_pars->SetBranchAddress("a5",           &a5); //j=5
    TB_pars->SetBranchAddress("a5_err",       &a5_err); 
    TB_pars->SetBranchAddress("b",            &b);
    TB_pars->SetBranchAddress("ChiSq",        &ChiSq); //j=6
    
    TFile *outfile      = new TFile("TB_cards_b2.root", "RECREATE"); //output TFile with histos + TTree
    TTree *TB_cards     = new TTree("TB_cards","TestBenchDB Qinj data card by card analysis"); //output TTree
    
    TB_cards->Branch("FEB_idc",      &FEB_idc,   "FEB_idc/I"); 
    TB_cards->Branch("cat_idc",      &cat_idc,   "cat_idc/I"); 
    TB_cards->Branch("data_uidc",    &data_uidc, "data_uidc/I"); 
    TB_cards->Branch("a00_mean",     &a00_mean,  "a00_mean/D");
    TB_cards->Branch("a00_sig",      &a00_sig,   "a00_sig/D"); 
    TB_cards->Branch("a1_mean",      &a1_mean,   "a1_mean/D"); 
    TB_cards->Branch("a1_sig",       &a1_sig,    "a1_sig/D");  
    TB_cards->Branch("a2_mean",      &a2_mean,   "a2_mean/D");
    TB_cards->Branch("a2_sig",       &a2_sig,    "a2_sig/D");  
    TB_cards->Branch("a3_mean",      &a3_mean,   "a3_mean/D");
    TB_cards->Branch("a3_sig",       &a3_sig,    "a3_sig/D");  
    TB_cards->Branch("a4_mean",      &a4_mean,   "a4_mean/D");
    TB_cards->Branch("a4_sig",       &a4_sig,    "a4_sig/D");  
    TB_cards->Branch("a5_mean",      &a5_mean,   "a5_mean/D");
    TB_cards->Branch("a5_sig",       &a5_sig,    "a5_sig/D");   
    TB_cards->Branch("ChiSq_mean",   &ChiSq_mean,"ChiSq_mean/D"); 
    TB_pars->Branch("ChiSq_sig",     &ChiSq_sig, "ChiSq_sig/D");   

    TString cat_label[8]    = {TString("For TT"), TString("Good spares"), TString("Medium spares"), TString("Bad Spares"), TString("Broken cards"), TString("Non existent cards"), TString("For Telescope"), TString("Other") };
    TH1F *hist_cards[2000][7];
    
    for(int i=0; i<TB_pars->GetEntries(); ++i)
    {
        TB_pars->GetEntry(i);  
        if(CH==0)
        {
            cout<<"Opening new FEB id "<<FEB_ID<<":"<<endl;
            TString hist_tag = TString("card") + TString(Form("%d",FEB_ID));
            hist_cards[FEB_ID][0] = new TH1F(hist_tag + TString("_par_") + TString(Form("%d",0)), hist_tag + TString("_par_") + TString(Form("%d",0)), Int_t(TB_pars->GetMaximum("a00")     - TB_pars->GetMinimum("a00"))   + 20, TB_pars->GetMinimum("a00"),   TB_pars->GetMaximum("a00") );  
            hist_cards[FEB_ID][1] = new TH1F(hist_tag + TString("_par_") + TString(Form("%d",1)), hist_tag + TString("_par_") + TString(Form("%d",1)), Int_t(TB_pars->GetMaximum("a1")      - TB_pars->GetMinimum("a1"))    + 20, TB_pars->GetMinimum("a1"),    TB_pars->GetMaximum("a1") );  
            hist_cards[FEB_ID][2] = new TH1F(hist_tag + TString("_par_") + TString(Form("%d",2)), hist_tag + TString("_par_") + TString(Form("%d",2)), Int_t(TB_pars->GetMaximum("a2")      - TB_pars->GetMinimum("a2"))    + 20, TB_pars->GetMinimum("a2"),    TB_pars->GetMaximum("a2") );  
            hist_cards[FEB_ID][3] = new TH1F(hist_tag + TString("_par_") + TString(Form("%d",3)), hist_tag + TString("_par_") + TString(Form("%d",3)), Int_t(TB_pars->GetMaximum("a3")      - TB_pars->GetMinimum("a3"))    + 20, TB_pars->GetMinimum("a3"),    TB_pars->GetMaximum("a3") );  
            hist_cards[FEB_ID][4] = new TH1F(hist_tag + TString("_par_") + TString(Form("%d",4)), hist_tag + TString("_par_") + TString(Form("%d",4)), Int_t(TB_pars->GetMaximum("a4")      - TB_pars->GetMinimum("a4"))    + 20, TB_pars->GetMinimum("a4"),    TB_pars->GetMaximum("a4") );  
            hist_cards[FEB_ID][5] = new TH1F(hist_tag + TString("_par_") + TString(Form("%d",5)), hist_tag + TString("_par_") + TString(Form("%d",5)), Int_t(TB_pars->GetMaximum("a5")      - TB_pars->GetMinimum("a5"))    + 20, TB_pars->GetMinimum("a5"),    TB_pars->GetMaximum("a5") );  
            hist_cards[FEB_ID][6] = new TH1F(hist_tag + TString("_par_") + TString(Form("%d",6)), hist_tag + TString("_par_") + TString(Form("%d",6)), Int_t(TB_pars->GetMaximum("ChiSq")   - TB_pars->GetMinimum("ChiSq")) + 20, TB_pars->GetMinimum("ChiSq"), TB_pars->GetMaximum("ChiSq") );           
        }
        hist_cards[FEB_ID][0]-> Fill(a00);
        hist_cards[FEB_ID][1]-> Fill(a1);
        hist_cards[FEB_ID][2]-> Fill(a2);
        hist_cards[FEB_ID][3]-> Fill(a3);
        hist_cards[FEB_ID][4]-> Fill(a4);
        hist_cards[FEB_ID][5]-> Fill(a5);
        hist_cards[FEB_ID][6]-> Fill(ChiSq);
        if(CH==63)
        {
            a00_mean    = hist_cards[FEB_ID][0]->GetMean();
            a00_sig     = hist_cards[FEB_ID][0]->GetRMS();
            a1_mean     = hist_cards[FEB_ID][1]->GetMean();
            a1_sig      = hist_cards[FEB_ID][1]->GetRMS();
            a2_mean     = hist_cards[FEB_ID][2]->GetMean();
            a2_sig      = hist_cards[FEB_ID][2]->GetRMS();
            a3_mean     = hist_cards[FEB_ID][3]->GetMean();
            a3_sig      = hist_cards[FEB_ID][3]->GetRMS();
            a4_mean     = hist_cards[FEB_ID][4]->GetMean();
            a4_sig      = hist_cards[FEB_ID][4]->GetRMS();
            a5_mean     = hist_cards[FEB_ID][5]->GetMean();
            a5_sig      = hist_cards[FEB_ID][5]->GetRMS();
            ChiSq_mean  = hist_cards[FEB_ID][6]->GetMean();
            ChiSq_sig   = hist_cards[FEB_ID][6]->GetRMS();
            FEB_idc = FEB_ID;    cat_idc = cat_ID;    data_uidc = data_UID; 
            for(int k=0; k<7; ++k)  hist_cards[FEB_ID][k]->Write();
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

