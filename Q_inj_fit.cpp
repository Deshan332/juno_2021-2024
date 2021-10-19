#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <math.h>
#include <TMath.h>
using namespace std;
Double_t Q_inj_fit(Double_t *X, Double_t *par)
{
    Double_t x  =*X;
    Double_t a0 = par[0];
    Double_t a00= par[1];
    Double_t a1 = par[2];
    Double_t a2;
    Double_t a3;
    Double_t a4;
    Double_t a5 = par[3];
    Double_t b  = par[4];
    Double_t val= 0;
    a4 = a1*b*(exp(a5) -1)/(a5* (a0-a00+a1*b) );
    a3 = a5*pow(b,-1*a4);
    a2 = a1*b*exp(a5)/(a4*a5);
    if(x<=b)        {val = a0 + a1*x;}
    else if (x>=b)  {val = a00 + a2*( 1 - exp( -1*a3*pow(x, a4) ) ); }
    return val;
}