#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <TMath.h>
using namespace std;
Double_t Gaus_fit(Double_t *X, Double_t *par)
{
    Double_t x   = X[0];
    Double_t N_0 = par[0];
    Double_t Q_0 = par[1];
    Double_t s_0 = par[2];
    Double_t G   = 0;
    if(s_0 != 0)    {G = N_0/( s_0*sqrt(2*M_PI) ) * exp(-pow( (x-Q_0),2 ) / ( 2 * pow(s_0,2) ) );}
    else            {cout<<"Error! division by 0 is undefined!"<<endl;}
    return G;
}