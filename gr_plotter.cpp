#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <TCanvas.h>
#include <TMath.h>
#include<TGraph.h>
#include<TGraphErrors.h>
#include<THistPainter.h>
#include <TFile.h>
#include <TStyle.h>
#include "TGrErr_getter.cpp"
using namespace std;
void gr_plotter()
{
    TCanvas* c1     = new TCanvas("c1");
    TFile *outfile  = new TFile("HD_10_1PE.root", "RECREATE");
    TGraphErrors* sshaper;
    Double_t X, X_err, Y, Y_err, Nmax, status =0;
    Double_t corr = 2.2;        //correction factor for cable delay between two signals
    vector<Double_t> x, x_err, y, y_err;
    ifstream ff ("HD_10_1PE.txt");
    if(!ff) {cout<<"No input file detected in the directory."<<endl;    exit(0);}
    while(ff >> X >> Y>>Y_err)
    {
        x.push_back(X-corr);
        y.push_back(Y);
        x_err.push_back(0);
        y_err.push_back(Y_err);
    }
    struct graph_data grdata = {x, x_err, y, y_err};        //structs needed to feed into TGrErr_getter
    TString labels[3]   = {TString("#Deltat(ns)"), TString("Mean charge (ADC)"), TString("Graph of #Deltat vs. Mean charge for card 3L (#Deltat = delay on ch64 - delay on ch14), Hold delay = 10, injecting 1PE")};        //array storing labels of graphs
    Double_t ddata[5]   = {1.2, -100., 310., 150., 0.};       // array storing double values needed for graphs        //{marker size, x_LL, x_UL, y_max, y_min}
    Int_t idata[4]      = {2, 20, 2, 1};                    // array storing int values needed for graphs           //{marker color, marker size, linewidth, line color}   
    sshaper             = TGrErr_getter(grdata, ddata, idata, labels);
    gStyle->SetPadGridX(kTRUE);
    gStyle->SetPadGridY(kTRUE);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    sshaper->Draw("ALP");
    sshaper->Write();
    outfile->Close();
    delete outfile;
}