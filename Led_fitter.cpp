#include <iostream>
#include <string>
#include <fstream>
#include <TH1F.h>
#include <TF1.h>
#include <math.h>
#include <TMath.h>
#include <TFile.h>
using namespace std;
TF1* Led_fitter(TH1F* hist, Double_t *par, Double_t *par_limits, Double_t *parerror, Double_t LL, Double_t UL)
{
    TF1 *Led_fit = new TF1 ("Led_fit", Bellamy, LL, UL,8);
    TString parnames[9] = {TString("N_{0}"), TString("Q_{0}"), TString("Q_{1}"), TString("#sigma_{0}"), TString("#sigma_{1}"), TString("w"), TString("#alpha"), TString("#mu"), TString("#Chi^{2}/NDOF") };
    for(int k=0; k<8; ++k)
        {
            Led_fit->SetParameter(k,par[k]);
            Led_fit->SetParLimits(k,par_limits[2*k],par_limits[2*k+1]);
            Led_fit->SetParName(k,parnames[k]);  
        }
    Led_fit->SetNpx(3000);
    hist->Fit("Led_fit","Rq","",LL,UL);
    for(int j=0; j<8; ++j)
    {
        par[j] = Led_fit->GetParameter(j);
        parerror[j] = Led_fit->GetParError(j);
    }
    par[8] = Led_fit->GetChisquare()/Led_fit->GetNDF();
    return Led_fit;
}