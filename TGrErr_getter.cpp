#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <TMath.h>
#include <TGraph.h>
#include <TAxis.h>
#include<TGraphErrors.h>
#include<THistPainter.h>
#include <TFile.h>
#include <TStyle.h>
using namespace std;
struct graph_data
{
    vector <Double_t> x;
    vector <Double_t> x_err;
    vector <Double_t> y;
    vector <Double_t> y_err;
};

TGraphErrors* TGrErr_getter(graph_data &gr_data, Double_t *double_data, Int_t *int_data, TString *labels)
{
    TGraphErrors* par_graph = new TGraphErrors(gr_data.x.size(),  &(gr_data.x[0]), &(gr_data.y[0]), &(gr_data.x_err[0]), &(gr_data.y_err[0]) );
    par_graph->SetMarkerSize(double_data[0]);
    par_graph->SetMarkerColor(int_data[0]);
    par_graph->SetMarkerStyle(int_data[1]);
    par_graph->SetLineWidth(int_data[2]);
    par_graph->SetLineColor(int_data[3]);
    par_graph->GetXaxis()->SetTitle(labels[0]);
    par_graph->GetYaxis()->SetTitle(labels[1]);
    par_graph->SetTitle(labels[2]);
    par_graph->GetXaxis()->SetRangeUser(double_data[1],double_data[2]);
    par_graph->SetMaximum(double_data[3]);
    par_graph->SetMinimum(double_data[4]);
    return par_graph;
}