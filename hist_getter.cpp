#include <iostream>
#include <string>
#include <fstream>
#include <TF1.h>
#include <math.h>
#include <TMath.h>
#include <TFile.h>
using namespace std;
TH1F* hist_getter(TString hist_path, int channel_num)
{
    TFile *filename = TFile::Open(hist_path);
    TString hist_name;
    if(channel_num <10)  {hist_name = TString("ch0") + Form("%d",channel_num);}
    if(channel_num >=10) {hist_name = TString("ch") + Form("%d",channel_num);}
    TH1F *hname = (TH1F*)filename->Get(hist_name);
    hname->SetDirectory(0);
    filename->Close();
    return hname;
}