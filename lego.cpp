#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>


void lego()
{
    TCanvas* c1 = new TCanvas("c1");
    TFile *myfile       = new TFile("/home/deshan/Documents/Internships/M2/New_env_mu/HV_Tree_2021-04-26_04:01.root", "READ");
    TTree *mytree       = (TTree*)myfile->Get("par_tree");
    TH2F *LEGO = new TH2F("2D", "MA-PMT Gain Distribution", 8,0.5,8.5,8,0.5,8.5);
    Int_t channel = 0, nEntries=0;
    Double_t Q_1 = 0;
    mytree->SetBranchAddress("Q_1",  &Q_1);
    for(int y=0; y<8; ++y)
        {
            for(int x=0; x<8; ++x)
                {
                    mytree->GetEntry(8*y+x);
                    LEGO->Fill(x+1,y+1,Q_1);
                }
        }
    LEGO->Draw("lego");
}