#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <math.h>
#include <TMath.h>
using namespace std;
Double_t lorentz_fit(Double_t *X, Double_t *par)
{
    Double_t x  =*X;
    Double_t N      = par[0];
    Double_t gamma  = par[1];
    //Double_t a      = par[1];
    Double_t x0     = par[2];
    //Double_t val    = 0;
    //Double_t den    = pow((x-x0),2) + pow((gamma/2),2);
    //val             = gamma/(2*M_PI*den);
    return (0.5*N*gamma/M_PI) / TMath::Max(1.e-10, (x-x0)*(x-x0) + 0.25*gamma*gamma);//lorentz
    //return N*(sqrt(2/M_PI)*x*x*exp(-x*x/(2*a*a))/(a*a*a)) ;//M-B
    //return (N/(x*x0*sqrt(2*M_PI)))*exp(-1*pow((log(x)-gamma),2)/(2*x0*x0)); //log-normal
    //return TMath::gaus(0) + TMath::Landau(3); //gaus + landau
}