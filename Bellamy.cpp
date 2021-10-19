#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <math.h>
#include <TMath.h>
using namespace std;
Double_t Bellamy(Double_t *X, Double_t *par)
{
    Double_t x   =*X;
    Double_t N_0 = par[0];
    Double_t Q_0 = par[1];
    Double_t Q_1 = par[2];
    Double_t s_0 = par[3];
    Double_t s_1 = par[4];
    Double_t w = par[5];
    Double_t alpha = par[6];
    Double_t mu = par[7];
    Double_t Q_sh = w/alpha;
    Double_t theta; if(x>=Q_0){theta=1.;}   else {theta=0.;}
    Double_t arg1_ped = -1*pow((x-Q_0),2)/ ( 2*pow(s_0,2) );
    Double_t arg2_ped = -1*alpha*(x-Q_0);
    Double_t Ped_fn = ( (1-w)*exp(arg1_ped)/(s_0*sqrt(2*M_PI)) + w*theta*alpha*exp(arg2_ped) )*exp(-1.*mu);
    Double_t Bel =0.;
    Double_t n=0.;
    Double_t arg_sig;
    for(int i=1; i<101; ++i)
    {
        n = 1.*i;
        arg_sig = -1*pow( (x- Q_0 - Q_sh - n*Q_1),2 )/( 2*n*pow(s_1,2) );        
        Bel += ( TMath::Poisson(n,mu) * exp(arg_sig) / (s_1*sqrt(2*M_PI*n)) );
    }    
    Bel= N_0*(Bel+Ped_fn);
    if(s_1>0.95*Q_1)    {Bel=0;}
    return Bel;
}