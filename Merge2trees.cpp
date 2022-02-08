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
    TString file1, file2, tr1, tr2;

    cout<<"Input the ROOT file no. 1"<<endl;
    cin>>file1;
    cout<<"Input the corresponding TTree name"<<endl;
    cin>>tr1;

    cout<<"Input the ROOT file no. 2"<<endl;
    cin>>file2;
    cout<<"Input the corresponding TTree name"<<endl;
    cin>>tr2;

    TFile *infile1  = new TFile(file1, "READ");
    TTree *tree1    = (TTree*)infile1->Get(tr1);
    list->Add(tree1);
    
    TFile *infile2  = new TFile(file2, "READ");  
    TTree *tree2    = (TTree*)infile2->Get(tr2);
    list->Add(tree2);

    TFile *outfile  = new TFile("Scurve_ThirdPE_sum.root","recreate");
    TTree *newtree = TTree::MergeTrees(list);

    newtree->Write();
    outfile->Close();
    delete outfile;
    return 0;

}