#include <iostream>
#include <string>
#include <fstream>
#include <istream>
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
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
using namespace std;

int CB_match()
{
    TFile *infile   = TFile::Open("cbv4Tel_cCB3_2021-10-27_18_12.root");
    TFile *outfile  = new TFile("CB_data_sum.root", "RECREATE");
    
    TTree *ROB_FR   = (TTree*)infile->Get("TT_CB_ROBFR_data");
    TTree *CB_TS    = (TTree*)infile->Get("TT_CB_TS");
    TTree *datanoTS = new TTree("datanoTS", "data but no TS");
    TTree *TSnodata = new TTree("TSnodata", "TS but no data");
    TTree *dataTS   = new TTree("dataTS",   "Good ones");
    TTree *RB_dT    = new TTree("RB_dT",   "Time steps in RB");
    TTree *CB_dT    = new TTree("CB_dT",   "Time steps in CB");
    TTree *CB_per_RB_dT    = new TTree("CB_per_RB_dT",   "Time steps in CB by ROB for CB data");
    TTree *RB_per_RB_dT    = new TTree("RB_per_RB_dT",   "Time steps in CB by ROB for ROB data");

    Int_t fNanoSec_RB, fNanoSec_CB, FNanoSec, i=0, n=0, m=0, min=0, max=0, matchcount=0, temp=0, check[8] = {0};
    UShort_t ROB_id_RB, ROB_id_CB, ROB_no, flag_RB, flag_CB, Flag;
    Int_t  fSec_RB, fSec_CB, FSec, offset = 0, t_before=0, delta_T=0, dT_by_RB=0, entry_no_CB=0, entry_no_RB=0, delta_ROB[8]={0}, delta_CB[8]={0}, t_ROB[8] = {0}, t_CB[8]={0};
    int one;

    ROB_FR->SetBranchAddress("ROB_id",       &ROB_id_RB);
    ROB_FR->SetBranchAddress("CB_flag",      &flag_RB);
    ROB_FR->SetBranchAddress("fSec",         &fSec_RB); 
    ROB_FR->SetBranchAddress("fNanoSec",     &fNanoSec_RB); 

    CB_TS->SetBranchAddress("ROB_id",       &ROB_id_CB);
    CB_TS->SetBranchAddress("flag",         &flag_CB);
    CB_TS->SetBranchAddress("fSec",         &fSec_CB); 
    CB_TS->SetBranchAddress("fNanoSec",     &fNanoSec_CB); 

    datanoTS->Branch("entry_no_RB",         &entry_no_RB,       "entry_no_RB/I");
    datanoTS->Branch("ROB_no",              &ROB_no,            "ROB_no/s");
    datanoTS->Branch("Flag",                &Flag,              "Flag/s");
    datanoTS->Branch("FSec",                &FSec,              "FSec/I");
    datanoTS->Branch("FNanoSec",            &FNanoSec,          "FNanoSec/I");

    TSnodata->Branch("entry_no_CB",         &entry_no_CB,       "entry_no_CB/I");
    TSnodata->Branch("ROB_no",              &ROB_no,            "ROB_no/s");
    TSnodata->Branch("Flag",                &Flag,              "Flag/s");
    TSnodata->Branch("FSec",                &FSec,              "FSec/I");
    TSnodata->Branch("FNanoSec",            &FNanoSec,          "FNanoSec/I");
 
    dataTS->Branch("entry_no_RB",       &entry_no_RB,       "entry_no_RB/I");
    dataTS->Branch("entry_no_CB",       &entry_no_CB,       "entry_no_CB/I");
    dataTS->Branch("ROB_no",            &ROB_no,            "ROB_no/s");
    dataTS->Branch("Flag",              &Flag,              "Flag/s");
    dataTS->Branch("FSec",              &FSec,              "FSec/I");
    dataTS->Branch("FNanoSec",          &FNanoSec,          "FNanoSec/I");


    RB_dT->Branch("ROB_no",         &ROB_no,            "ROB_no/s");
    RB_dT->Branch("Flag",           &Flag,              "Flag/s");
    RB_dT->Branch("FSec",           &FSec,              "FSec/I");
    RB_dT->Branch("FNanoSec",       &FNanoSec,          "FNanoSec/I");
    RB_dT->Branch("delta_T",        &delta_T,           "delta_T/I");


    CB_dT->Branch("ROB_no",         &ROB_no,            "ROB_no/s");
    CB_dT->Branch("Flag",           &Flag,              "Flag/s");
    CB_dT->Branch("FSec",           &FSec,              "FSec/I");
    CB_dT->Branch("FNanoSec",       &FNanoSec,          "FNanoSec/I");
    CB_dT->Branch("delta_T",        &delta_T,           "delta_T/I");

    CB_per_RB_dT->Branch("ROB_no",         &ROB_no,            "ROB_no/s");
    CB_per_RB_dT->Branch("Flag",           &Flag,              "Flag/s");
    CB_per_RB_dT->Branch("FSec",           &FSec,              "FSec/I");
    CB_per_RB_dT->Branch("FNanoSec",       &FNanoSec,          "FNanoSec/I");
    CB_per_RB_dT->Branch("delta_T",        &dT_by_RB,          "dT_by_RB/I");

    RB_per_RB_dT->Branch("ROB_no",         &ROB_no,            "ROB_no/s");
    RB_per_RB_dT->Branch("Flag",           &Flag,              "Flag/s");
    RB_per_RB_dT->Branch("FSec",           &FSec,              "FSec/I");
    RB_per_RB_dT->Branch("FNanoSec",       &FNanoSec,          "FNanoSec/I");
    RB_per_RB_dT->Branch("delta_T",        &dT_by_RB,          "dT_by_RB/I");

    for(Int_t i=0; i<ROB_FR->GetEntries();  ++i) 
    {
        ROB_FR->GetEntry(i);
        //checking entries with TS!=0
        if((fNanoSec_RB>0)&&(fSec_RB>0))
        {
            if(i==0)                     t_before = fNanoSec_RB;                 //for entry no. 0
            if(check[ROB_id_RB]==0)     
            {
                t_ROB[ROB_id_RB]    = fNanoSec_RB;         //for the 1st event for each ROB, the corresponding TS is the earliest TS (so that delta_ROB = 0 for 1st entry)
                check[ROB_id_RB]    = 1;
                cout<<"1st entry for ROB "<<ROB_id_RB<<" detected at entry no. "<<i<<" on ROB data."<<endl;
            }
            if(i>0) 
            {
                delta_T                 = fNanoSec_RB-t_before; 
                delta_ROB[ROB_id_RB]    = fNanoSec_RB-t_ROB[ROB_id_RB];
                t_before = fNanoSec_RB;         t_ROB[ROB_id_RB] = fNanoSec_RB;
                ROB_no = ROB_id_RB; Flag = flag_RB; FSec = fSec_RB; FNanoSec = fNanoSec_RB; entry_no_RB = i;
                RB_dT->Fill();
                dT_by_RB = delta_ROB[ROB_id_RB];
                RB_per_RB_dT->Fill();
            }
        }

        min = i+offset-100;  max = i+offset+100;
        if(min< 0) min=0;    if(max>CB_TS->GetEntries()) max = CB_TS->GetEntries();

        for(Int_t j=min; j<max;  ++j) 
        {
           temp=0;
           CB_TS->GetEntry(j);
           if((ROB_id_CB==ROB_id_RB) && (flag_CB==flag_RB) && (fSec_CB==fSec_RB) && (fNanoSec_CB==fNanoSec_RB) ) //finding matching entries between 2 TTrees
           {
               offset = j-i;
               //cout<<"offset = "<<offset<<endl;
               ROB_no = ROB_id_RB; Flag = flag_RB; FSec = fSec_RB; FNanoSec = fNanoSec_RB;  entry_no_RB = i;    entry_no_CB = j;
               dataTS->Fill();
               temp=1;
           }
           if(temp==1) break;
        }
        if(temp==0) //filter entries in ROB that aren't seen in CB, and fill them in datanoTS
        {
            ROB_no = ROB_id_RB; Flag = flag_RB; FSec = fSec_RB; FNanoSec = fNanoSec_RB; entry_no_RB = i;
            datanoTS->Fill(); 
        }
    }

    offset=0;   delta_T=0;  t_before=0; //resetting offset, delta_T, t_before counters
    for(Int_t p=0; p<8; ++p)            //resetting check and t_ROB
    {
        check[p]    = 0;
        t_ROB[p]    = 0;
    }       
    for(Int_t j=0; j<CB_TS->GetEntries();  ++j) 
    {
        CB_TS->GetEntry(j);
        if((fNanoSec_CB>0)&&(fSec_CB>0))
        {
            if(j==0)                     t_before = fNanoSec_CB;                 //for entry no. 0
            if(check[ROB_id_CB]==0)     
            {
                t_CB[ROB_id_CB]     = fNanoSec_CB;         //for the 1st event for each ROB, the corresponding TS is the earliest TS (so that delta_ROB = 0 for 1st entry)
                check[ROB_id_CB]    = 1;
                cout<<"1st entry for ROB "<<ROB_id_CB<<" detected at entry no. "<<j<<" on CB data."<<endl;
            }
            if(j>0) 
            {
                delta_T                = fNanoSec_CB-t_before; 
                delta_CB[ROB_id_CB]    = fNanoSec_CB-t_CB[ROB_id_CB];
                t_before = fNanoSec_CB;         t_CB[ROB_id_CB] = fNanoSec_CB;
                ROB_no = ROB_id_CB; Flag = flag_CB; FSec = fSec_CB; FNanoSec = fNanoSec_CB; entry_no_CB = j;
                CB_dT->Fill();
                dT_by_RB = delta_CB[ROB_id_CB];
                CB_per_RB_dT->Fill();
            }
        }

        min = j-offset-100;  max = j-offset+100;
        if(min< 0) min=0;    if(max>ROB_FR->GetEntries()) max = ROB_FR->GetEntries();
        for(Int_t i=min; i<max;  ++i) 
        {
           temp=0;
           ROB_FR->GetEntry(i);
           if((ROB_id_RB==ROB_id_CB) && (flag_RB==flag_CB) && (fSec_RB==fSec_CB) && (fNanoSec_RB==fNanoSec_CB) ) //filter entries in CB that aren't seen in ROB, and fill them in TSnodata
           {
               offset = j-i;
               //cout<<"offset type 2 = "<<offset<<endl;
               temp=1;
               break;
           }
        }
        if(temp==0) 
        {
            ROB_no = ROB_id_CB; Flag = flag_CB; FSec = fSec_CB; FNanoSec = fNanoSec_CB; entry_no_CB = j;  
            TSnodata->Fill();
        }
        
    }
    
    TSnodata->Write();
    datanoTS->Write();
    dataTS->Write();
    RB_dT->Write();
    CB_dT->Write();
    RB_per_RB_dT->Write();
    CB_per_RB_dT->Write();
    outfile->Close();
    delete outfile;
    return 0;

}