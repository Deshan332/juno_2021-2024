#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <TH1F.h>
#include <TF1.h>
#include <math.h>
#include <TMath.h>
#include<TGraph.h>
#include<TGraphErrors.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TMinuit.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TText.h>
#include <TPaveStats.h>
#include "TGrErr_getter.cpp"
#include "Q_function.cpp"
using namespace std;

void s_fitter()
{
    TCanvas* c1     = new TCanvas("c1");
    TFile *outfile  = new TFile("s_test_2021_04_14.root", "RECREATE");
    vector<Double_t> x, x_err, y, y_err;
    Int_t count=0;
    Double_t X, X_err, Y, Y_err, Nmax, status =0;
    TGraphErrors* s_curve;
    ifstream ff ("DAC_data_2021_04_14_3L_ch14_1_3_PE.txt");
    if(!ff) {cout<<"No input file detected in the directory."<<endl;    exit(0);}
    Double_t LL = 380;
    Double_t UL = 480;
    while(ff >> X >> Y)
    {
        if(count==0){Nmax = Y;  ++count;}
        x.push_back(X);
        y.push_back(Y/Nmax);
        x_err.push_back(0);
        Y_err   =   (Y/Nmax)*sqrt( (1/Y) + (1/Nmax) );
        y_err.push_back(Y_err);
    }
    ff.close();
    struct graph_data gdata = {x, x_err, y, y_err};
    Double_t par[3]         = {1., 400., 1.};
    Double_t par_limits[6]  = {0.95*par[0], 1.05*par[0], 380., 480., 0., 10.};
    Double_t ddata[5]       = {1, 380., 480., 1.1, -0.1};       // array storing double values needed for graphs        //{marker size, x_LL, x_UL, y_max, y_min}
    Int_t idata[4]          = {4, 20, 2, 7};                // array storing int values needed for graphs           //{marker color, marker size, linewidth, line color}
    TString labels[3]       = {TString("DAC (a.u.)"), TString("Trigger rate / N_{max}"), TString("Card 3L, channel 14, DAC threshold determination for 1/3PE") };
    s_curve                 = TGrErr_getter(gdata, ddata, idata, labels);   
    s_curve->Draw("ALP"); 
    TF1 *Q_fit = new TF1 ("Q_fit", Q_function, LL, UL,3);
    TString parnames[3] = {TString("N_{0}"), TString("#mu"), TString("#sigma_{0}")};
    
    for(int k=0; k<3; ++k)
        {
            Q_fit->SetParameter(k,par[k]);
            Q_fit->SetParLimits(k,par_limits[2*k],par_limits[2*k+1]);
            Q_fit->SetParName(k,parnames[k]);  
        }
    Q_fit->SetNpx(3000);
    for(int i=0; i<40; ++i)
    {
        Q_fit->FixParameter(0, 1);
        s_curve->Fit("Q_fit","Rq","",LL,UL);
        for(int j=0; j<3; ++j)
        {
            par[j] = Q_fit->GetParameter(j);
            Q_fit->SetParameter(j,par[j]);
        }
        if( (gMinuit->fCstatu.Contains("CONVERGED"))  && (Q_fit->GetChisquare()/Q_fit->GetNDF() < 5) )  {status =1;    Q_fit->SetLineColor(2);  break;}
        else    {status =0;    Q_fit->SetLineColor(1);}
    }
    gStyle->SetOptStat(111);
    gStyle->SetOptFit(111);
    s_curve->Write();
    outfile->Write();
    cout<<"Fit "<<gMinuit->fCstatu<<" with Chisq/NDOF = "<<Q_fit->GetChisquare()/Q_fit->GetNDF()<<endl;
    outfile->Close();
}