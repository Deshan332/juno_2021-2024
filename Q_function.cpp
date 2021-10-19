#include <math.h>
#include <TMath.h>
Double_t Q_function(Double_t *Xarr, Double_t *par)
{
    Double_t x      = *Xarr;
    Double_t N_0    = par[0];
    Double_t mu     = par[1];
    Double_t sigma  = par[2];
    Double_t y      = (x-mu)/sigma;
    Double_t Q      = 0.5*N_0*TMath::Erfc(y/sqrt(2));
    return Q;
}