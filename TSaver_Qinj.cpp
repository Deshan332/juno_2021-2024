#include <iostream>
#include <string>
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
int TSaver_Qinj(TString filename)
{
    ifstream read;
    read.open(filename);
    if (!read)
    {
      cout << "the file " + filename + " does not exist!" << endl;
      return 0;
    } 

    TFile * Outfile     =   new TFile("Q_inj_tree.root", "RECREATE"); 
    TTree * TB_Qinj     =   new TTree("TB_Qinj","TestBenchDB Qinj data");
    
    Int_t FEB_id, cat_id, data_uid, ch, QinjNum, n_entries, max, min;//TBranches
    Double_t mean, std;
    TString category;
    TB_Qinj->Branch("FEB_id",       &FEB_id,    "FEB_id/I");
    TB_Qinj->Branch("cat_id",       &cat_id,    "cat_id/I");
    TB_Qinj->Branch("data_uid",     &data_uid,  "data_uid/I"); 
    TB_Qinj->Branch("ch",           &ch,        "ch/I"); 
    TB_Qinj->Branch("QinjNum",      &QinjNum,   "QinjNum/I"); 
    TB_Qinj->Branch("mean",         &mean,      "mean/D"); 
    TB_Qinj->Branch("std",          &std,       "std/D");
    TB_Qinj->Branch("n_entries",    &n_entries, "n_entries/I"); 
    TB_Qinj->Branch("max",          &max,       "max/I");
    TB_Qinj->Branch("min",          &min,       "min/I"); 
    
    read.ignore(1000000000000, '\n');
    while (read >> FEB_id >> cat_id >> data_uid >> ch >> QinjNum >> mean >> std >> n_entries >> max >> min) 
      {TB_Qinj->Fill();}
    read.close();    
    Outfile->Write();
    Outfile->Close();
    delete Outfile;

}

int main(int argc, char ** argv)
{
  if(argc<2) return 1;
  TString fname = argv[1];
  TApplication app("app", &argc, argv);
  TSaver_Qinj(fname);
  app.Run();
  return 0;
}