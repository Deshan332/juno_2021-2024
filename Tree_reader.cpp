#include <iostream>
#include <string>
#include <fstream>
#include <TF1.h>
#include <math.h>
#include <TMath.h>
#include <TTree.h>
#include <THistPainter.h>
#include <TFile.h>
void Tree_reader()
{
    ofstream ff ("gaindata_t1.txt");
    TFile *filename     = TFile::Open("/home/deshan/Documents/Internships/M2/New_env_mu/Cards_rerun_tree_2021-07-02_12:47.root","READ");
    TTree *par_tree     = (TTree*)filename->Get("par_tree");
    double Q_1;
    int module, LRhand, channel, status;
    par_tree->SetBranchAddress("module",&module);
    par_tree->SetBranchAddress("LRhand",&LRhand);
    par_tree->SetBranchAddress("channel",&channel);
    par_tree->SetBranchAddress("Q_1",&Q_1);
    par_tree->SetBranchAddress("status",&status);
    int N = par_tree->GetEntries();
    for(int i=0; i<N; ++i)
    {
        par_tree->GetEntry(i);
        ff << module << "  "<< LRhand << "  " << channel << "  " << Q_1 << "  " << status <<endl; 
    }
    ff.close();
}