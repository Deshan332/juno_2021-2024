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

int Merge2trees()
{
    TList *list     = new TList();
    
    TFile *infile1  = new TFile("TB_Scurve_ThirdPE_pars.root", "READ");
    TTree *tree1    = (TTree*)infile1->Get("TB_Scurve_par");
    list->Add(tree1);
    
    TFile *infile2  = new TFile("TB_Scurve_ThirdPE_pars2.root", "READ");  
    TTree *tree2    = (TTree*)infile2->Get("TB_Scurve_par");
    list->Add(tree2);

    TFile *outfile  = new TFile("Scurve_ThirdPE_sum.root","recreate");
    TTree *newtree = TTree::MergeTrees(list);

    newtree->Write();
    outfile->Close();
    delete outfile;
    return 0;

}