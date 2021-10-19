#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <TMath.h>
#include <TGraph.h>
#include <TAxis.h>
#include<TGraphErrors.h>
#include<THistPainter.h>
#include <TFile.h>
#include <TStyle.h>
#include <TTree.h>
using namespace std;
int TSaver_Scurve()
{
    TString filename = "Scurve_data_20_10_2021.txt";
    ifstream read;
    read.open(filename);
    if (!read)
    {
      cout << "the file " + filename + " does not exist!" << endl;
      return 0;
    } 

    TFile * Outfile     =   new TFile("Scurve_tree.root", "RECREATE"); 
    TTree * TB_Scurve   =   new TTree("TB_Scurve","TestBenchDB Scurve data");
    
    Int_t FEB_id, cat_id, data_uid, test_id, ch, ch_inj, sid, DAC, counts;//TBranches
    Double_t thresh;
    TString category;
    string Th;
    TB_Scurve->Branch("FEB_id",       &FEB_id,    "FEB_id/I");
    TB_Scurve->Branch("cat_id",       &cat_id,    "cat_id/I");
    TB_Scurve->Branch("data_uid",     &data_uid,  "data_uid/I"); 
    TB_Scurve->Branch("test_id",      &test_id,   "test_id/I"); 
    TB_Scurve->Branch("ch",           &ch,        "ch/I"); 
    TB_Scurve->Branch("ch_inj",       &ch_inj,    "ch_inj/I"); 
    TB_Scurve->Branch("thresh",       &thresh,    "thresh/D"); 
    TB_Scurve->Branch("sid",          &sid,       "sid/I");
    TB_Scurve->Branch("DAC",          &DAC,       "DAC/I"); 
    TB_Scurve->Branch("counts",       &counts,    "counts/I");

    
    read.ignore(1000000000000, '\n');//ignore 1st line of text file
    while (read >> FEB_id >> cat_id >> data_uid >> test_id >> ch >> ch_inj >> Th >> sid >> DAC >> counts) 
    {
       if(Th =="NULL") continue;
       else
       {
           thresh = stod(Th);
           TB_Scurve->Fill();
       }
    }
    read.close();    
    Outfile->Write();
    Outfile->Close();
    delete Outfile;
    return 0;
}
