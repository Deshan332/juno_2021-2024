#include <iostream>
#include <string>
#include <fstream>
#include <TH1F.h>
#include <TF1.h>
#include <TFile.h>
using namespace std;
TF1* Ped_fitter(TH1F* hist, Double_t *par, Double_t *parerror)
{
    TF1 *Ped_fit = new TF1 ("Ped_fit", Gaus_fit, par[1]-2, par[1]+2, 3);
    Ped_fit->SetParameter(0,par[0]);    //N_0
    Ped_fit->SetParameter(1,par[1]);    //Q_0
    Ped_fit->SetParameter(2,par[2]);    //s_0
    Ped_fit->SetNpx(100);
    hist->Fit("Ped_fit", "Rq");
    for(int j=0; j<3; ++j)
    {
        par[j] = Ped_fit->GetParameter(j);
        parerror[j] = Ped_fit->GetParError(j);
    }
    par[8] = Ped_fit->GetChisquare()/Ped_fit->GetNDF();
    return Ped_fit;
}