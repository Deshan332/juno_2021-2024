#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <TH1F.h>
#include <TF1.h>
#include <math.h>
#include <TMath.h>
#include <TTree.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <THistPainter.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TMinuit.h>
#include <TStyle.h>
#include <TText.h>
#include <TPaveStats.h>
#include <TApplication.h>
#include "TGrErr_getter.cpp"
#include "Q_inj_fit.cpp"

using namespace std;
int Q_inj_fft()
{
    TCanvas* c1         = new TCanvas("c1");
    Int_t FEB_id=0, FEB_ID, cat_id, cat_ID, data_uid, data_UID, ch, CH, QinjNum, n_entries, max, min, cats = 7,  count=0, temp=0, i=0, fcount = 0, n=0, f_id=1, status;
    Double_t mean, std, a0, a00, a00_err, a1, a1_err, a2, a3, a3_err, a4, a5, a5_err, b, PED=0, ChiSq;
    TString filename = "Q_inj_tree.root", tree_name = "TB_Qinj";
    TFile *myfile       = new TFile(filename, "READ"); //input TFile
    TTree *mytree       = (TTree*)myfile->Get(tree_name);
    mytree ->SetBranchAddress("FEB_id",       &FEB_id);
    mytree ->SetBranchAddress("cat_id",       &cat_id);
    mytree ->SetBranchAddress("data_uid",     &data_uid); 
    mytree ->SetBranchAddress("ch",           &ch); 
    mytree ->SetBranchAddress("QinjNum",      &QinjNum); 
    mytree ->SetBranchAddress("mean",         &mean); 
    mytree ->SetBranchAddress("std",          &std);
    mytree ->SetBranchAddress("n_entries",    &n_entries); 
    mytree ->SetBranchAddress("max",          &max);
    mytree ->SetBranchAddress("min",          &min); 

    TFile *outfile[20]; 
    TFile *treefile     = new TFile("test_b2.root", "RECREATE"); 
    TTree *TB_pars      = new TTree("TB_lin_par","TestBenchDB Qinj data linearity fit parameters");
    
    TB_pars->Branch("FEB_ID",       &FEB_ID,    "FEB_ID/I");
    TB_pars->Branch("cat_ID",       &cat_ID,    "cat_ID/I");
    TB_pars->Branch("data_UID",     &data_UID,  "data_UID/I"); 
    TB_pars->Branch("CH",           &CH,        "CH/I"); 
    TB_pars->Branch("a0",           &a0,        "a0/D"); 
    TB_pars->Branch("a00",          &a00,       "a00/D");
    TB_pars->Branch("a00_err",      &a00_err,   "a00_err/D"); 
    TB_pars->Branch("a1",           &a1,        "a1/D");
    TB_pars->Branch("a1_err",       &a1_err,    "a1_err/D"); 
    TB_pars->Branch("a2",           &a2,        "a2/D");
    TB_pars->Branch("a3",           &a3,        "a3/D");
    TB_pars->Branch("a4",           &a4,        "a4/D");
    TB_pars->Branch("a5",           &a5,        "a5/D");
    TB_pars->Branch("a5_err",       &a5_err,    "a5_err/D"); 
    TB_pars->Branch("b",            &b,         "b/D");
    TB_pars->Branch("ChiSq",        &ChiSq,     "ChiSq/D");  
    TB_pars->Branch("status",       &status,    "status/I"); 

    vector<Double_t> ch_pC, ch_pC_err, ch_ADC, ch_ADC_err; 
    TGraphErrors *mygraph;
    TString labels[3]       = {TString("Charge Injected (pC)"), TString("Charge collected (ADC)"), TString("")};        //array storing labels of graphs  xlabel, ylabel, title
    Double_t charges[5] {0.056, 0.158, 0.5, 2., 10.}; //x charges in pC
    for(int a=0; a<5; ++a){ch_pC.push_back(charges[a]);     ch_pC_err.push_back(0);}    //x and x error
    Double_t ddata[5]       = {1.2, 0., 10.5, 450, 0.};       // array storing double values needed for graphs        //{marker size, x_LL, x_UL, y_max, y_min}
    Int_t idata[4]          = {4, 20, 2, 7};                  // array storing int values needed for graphs           //{marker color, marker size, linewidth, line color}
    TString cat_label[8]    = {TString("For TT"), TString("Good spares"), TString("Medium spares"), TString("Bad Spares"), TString("Broken cards"), TString("Non existent cards"), TString("For Telescope"), TString("Other") };
    //cout<<"here"<<endl;
    
    for(int i=0; i<mytree->GetEntries(); ++i)
    {
        if(f_id==1) //f_id is the TFile creation index. f_id =1 ---> call to Create a new TFile
        {
            TString file_tag    = TString("testb2_") + TString(Form("%d",n)) + TString(".root");
            outfile[n]          = new TFile(file_tag, "RECREATE"); //n is the index of TFile count.
            f_id = 0;
            cout<<"New TFile created. (n = "<<n<<")"<<endl;
        }
        mytree->GetEntry(i);
        FEB_ID = FEB_id;    cat_ID = cat_id;    data_UID = data_uid;    CH = ch;     
        
        if(cat_id==1)
        {
            if(FEB_id!=temp)    {++count;  temp = FEB_id;} //count keeps the count of concerned FEB type only (per batch)

            if(QinjNum==0)      {ch_ADC.clear();    ch_ADC_err.clear();    PED=mean;}
            else
            {
                ch_ADC.push_back(mean - PED); //y in ADC
                ch_ADC_err.push_back(std);
            }

            if(QinjNum==5)
            {
                labels[2] = TString("Graph of charge collected (ADC) vs. charge injected (pC) for FEB no. ") + TString(Form("%d",FEB_id)) + TString(", channel ") + TString(Form("%d",ch)) + TString(" (Category: ") + cat_label[cat_id-1] + TString(")");
                struct graph_data gdata = {ch_pC, ch_pC_err, ch_ADC, ch_ADC_err};
                mygraph = TGrErr_getter(gdata, ddata, idata, labels); 
                mygraph->SetName(TString("FEB") + TString(Form("%d", FEB_id)) + TString("ch") + TString(Form("%d",ch)));
                mygraph->Draw("AP");

                TF1 *fit_pol3 = new TF1 ("fit_pol3", Q_inj_fit, 0, 10, 5);
                fit_pol3->SetParameters(0., 0., 50., 1., 2.);
                fit_pol3->FixParameter(0,0);
                fit_pol3->FixParameter(4,2);
                gStyle->SetOptStat(111);
                gStyle->SetOptFit(111);
                gStyle->SetStatX(0.5);
                gStyle->SetStatY(0.9);
                fit_pol3->SetNpx(30000);
                c1->Update();
                fit_pol3->SetParNames("a_{0}", "a_{00}", "a_{1}", "a_{5}", "b");
                fit_pol3->SetParLimits(3,0,10);  
                fit_pol3->SetParLimits(1,-20,20);
                for(int x=0; x<10; ++x)
                {
                    mygraph->Fit("fit_pol3","REq");
                    a0 = fit_pol3->GetParameter(0);  a00 = fit_pol3->GetParameter(1);       a1 = fit_pol3->GetParameter(2); a5 = fit_pol3->GetParameter(3); b = fit_pol3->GetParameter(4);
                    fit_pol3->SetParameter(0,a0);      fit_pol3->SetParameter(1,a00);          fit_pol3->SetParameter(2,a1);   fit_pol3->SetParameter(3,a5);  fit_pol3->SetParameter(4,b);
                    ChiSq = fit_pol3->GetChisquare()/fit_pol3->GetNDF();
                    if( (gMinuit->fCstatu.Contains("CONVERGED")) || (gMinuit->fCstatu.Contains("SUCCESSFUL"))  )  {fit_pol3->SetLineColor(3); status =1;  break;}
                    else status = 0;
                }
                a00_err = fit_pol3->GetParError(1);  a1_err = fit_pol3->GetParError(2);  a5_err = fit_pol3->GetParError(3);   
                a4 = a1*b*(exp(a5) -1)/(a5* (a0-a00+a1*b) );
                a2 = a1*b*exp(a5)/(a4*a5);
                a3 = a5*pow(b,-1*a4);
                
                //cout<<"Chisq/NDOF = "<<ChiSq<<endl;
                fit_pol3->Draw("SAME");
                treefile->cd();
                TB_pars->Fill();
                outfile[n]->cd();
                mygraph->Write();
            }
        }

        if(count==100)
        {
            f_id=1;
            outfile[n]->Close();
            delete outfile[n];
            cout<<"TFile closed."<<endl;
            ++n;  
            count=0;  
        }
    }
    outfile[n]->Close();
    delete outfile[n];
    cout<<"TFile closed."<<endl;
    treefile->cd();
    TB_pars->Write();
    treefile->Close();
    cout<<"treefile closed."<<endl;
    delete treefile;
    return 0;
}

