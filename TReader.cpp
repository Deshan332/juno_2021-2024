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
void TReader()
{
    TFile *myfile       = new TFile("/home/deshan/Documents/Internships/M2/New_env_mu/Cards_rerun_tree_2021-07-02_12:47.root", "READ");
    TTree *mytree       = (TTree*)myfile->Get("par_tree");
    ofstream ff ("gain_data.txt");
    Int_t module = 0;
    Int_t LRhand = 0;
    Int_t channel = 0;
    Double_t Q_1 = 0.;
    Int_t status = 0, temp = 0, N = 0;
    double a0 = 0., a00 = 1.488, a1 = 69.94, a5 = 0.5384, b = 2.;
    mytree->SetBranchAddress("module",  &module);
    mytree->SetBranchAddress("LRhand",  &LRhand);
    mytree->SetBranchAddress("channel",  &channel);
    mytree->SetBranchAddress("Q_1",  &Q_1);
    mytree->SetBranchAddress("status",  &status);
    N = mytree->GetEntries();
    for(int i=0; i<N; ++i)
    {
        mytree->GetEntry(i);
        if(LRhand==1){temp = 1;}    if(LRhand==0){temp = 0;}
        ff<<module<<" "<<temp<<" "<<channel-1<<" "<<Q_1<<" "<<status<<" "<<a0<<" "<<a00<<" "<<a1<<" "<<a5<<" "<<b<<endl;
        //ff<<module<<" "<<temp<<" "<<channel-1<<" "<<Q_1<<" "<<status<<endl;
    }
    ff.close();
}