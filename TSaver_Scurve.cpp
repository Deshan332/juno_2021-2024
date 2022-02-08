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
#include <TString.h>
#include <TDatime.h>
using namespace std;
int TSaver_Scurve()
{
    TString filename = "Scurve_data_26_11_2021.txt";
    ifstream read;
    read.open(filename);
    if (!read)
    {
      cout << "the file " + filename + " does not exist!" << endl;
      return 0;
    } 

    TFile * Outfile     =   new TFile("Scurve_tree_datetime.root", "RECREATE"); 
    TTree * TB_Scurve   =   new TTree("TB_Scurve","TestBenchDB Scurve data");
    
    Int_t FEB_id, cat_id, data_uid, test_id, ch, ch_inj, sid, DAC, counts, TStamp;//TBranches
    Double_t thresh;
    Int_t date, time;
    TString category, d, t, th;
    string Th, dd, tt, TS;
    TDatime date_time;
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
    TB_Scurve->Branch("date",         &date,      "date/I");
    TB_Scurve->Branch("time",         &time,      "time/I");
    //TB_Scurve->Branch("TStamp",       &TStamp,    "TStamp/I");
    //TB_Scurve->Branch("date_time",    &date_time);

    
    read.ignore(1000000000000, '\n');//ignore 1st line of text file
    while (read >> FEB_id >> cat_id >> data_uid >> test_id >> ch >> ch_inj >> Th >> sid >> DAC >> counts >> dd >> tt) 
    {
       if(Th =="NULL") continue;
       else
       {
           
           th = TString(Th);
           d  = TString(dd);
           t  = TString(tt);
           d.ReplaceAll("-","");
           t.ReplaceAll(":","");
           thresh = th.Atof();
           date   = d.Atoi();
           time   = t.Atoi();
           //date_time.Set(TS.c_str());
           //TStamp = date_time.Convert();
           //cout<<TStamp<<endl;
           TB_Scurve->Fill();
       }
    }
    read.close();    
    Outfile->Write();
    Outfile->Close();
    delete Outfile;
    return 0;
}
