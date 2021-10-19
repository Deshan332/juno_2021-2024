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
#include "TGrErr_getter.cpp"
#include "Q_inj_fit.cpp"
#include "timestamp.cpp"
using namespace std;
void Q_test()
{
    TCanvas* c1         = new TCanvas("c1");
    TString run_time    = timestamp();
    TFile *outfile      = new TFile(TString("linear_test_TB") + TString(run_time) + TString(".root"), "RECREATE");
    ofstream ff ("gain_par.txt");
    TGraphErrors *mygraph[8][64];
    TH1F *hname[5][8][64];
    TH1F *pname[5][8][64];
    TH1F *g_hist= new TH1F("g", "Distribution of pC->ADC linear conversion factors for the PMT channels of the Muon Telescope", 120, 50., 80.);
    TH1F *a2_hist= new TH1F("a2", "Distribution of a_{2} for the PMT channels of the Muon Telescope", 400, 250., 450.);
    TH1F *a3_hist= new TH1F("a3", "Distribution of a_{3} for the PMT channels of the Muon Telescope", 50, 0., 1.);
    g_hist->GetXaxis()->SetTitle("a_{1} (ADC/pC)");
    g_hist->GetYaxis()->SetTitle("Counts");
    a2_hist->GetXaxis()->SetTitle("a_{2} (ADC)");
    a2_hist->GetYaxis()->SetTitle("Counts");
    TFile *filename[8];
    TString hist_name;
    TString card_num[8] = {"5002", "5008", "5003", "5009", "5007", "5011", "4004", "4003"};
    TString ped_name;
    TString labels[3]   = {TString("Charge Injected (pC)"), TString("Charge collected (ADC)"), TString("") };        //array storing labels of graphs  xlabel, ylabel, title
    Double_t ddata[5]   = {1.2, 0., 10.5, 450, 0.};       // array storing double values needed for graphs        //{marker size, x_LL, x_UL, y_max, y_min}
    Int_t idata[4]      = {4, 20, 2, 7};                   // array storing int values needed for graphs           //{marker color, marker size, linewidth, line color}
    Double_t charges[5]      = {0.056, 0.158, 0.5, 2., 10.};
    Double_t charge_err[5]   = {0., 0., 0., 0., 0.};
    Double_t a0, a00, a1, a2, a3, a4, a5, b;
    vector<Double_t> mean[8][64], mean_err[8][64], charge_inj[8][64], charge_inj_err[8][64];
    struct graph_data gdata[8][64];
    TF1 *fit_pol3[8][64];
    c1->Print( TString("lin_test_TB") + TString(run_time) + TString(".pdf[") );  
    for(int card=0; card<8; ++card)
    {   
        filename[card]     = TFile::Open("/home/deshan/Documents/Internships/M2/TestBenchDB_MT_data/Test08_QinjTest_hist_" + card_num[card] + ".root","READ");

        for(int ch=0; ch<64; ++ch)
        {
            for(int Q_num=0; Q_num<5; ++Q_num)
            {
                if(ch <10)  
                {
                    hist_name = TString("h_ch0") + Form("%d",ch) + TString("_Q #") + TString(Form("%d",Q_num+1));
                    ped_name  = TString("h_ch0") + Form("%d",ch) + TString("_ped");
                }
                if(ch >=10) 
                {
                    hist_name = TString("h_ch")  + Form("%d",ch) + TString("_Q #") + TString(Form("%d",Q_num+1));
                    ped_name  = TString("h_ch") + Form("%d",ch) + TString("_ped");
                }
                pname[Q_num][card][ch] = (TH1F*)filename[card]->Get(ped_name);
                hname[Q_num][card][ch] = (TH1F*)filename[card]->Get(hist_name);
                mean[card][ch].push_back( hname[Q_num][card][ch]->GetMean() - pname[Q_num][card][ch]->GetMean() );
                mean_err[card][ch].push_back( sqrt( pow( hname[Q_num][card][ch]->GetMeanError() , 2 ) + pow( pname[Q_num][card][ch]->GetMeanError() , 2 ) ) );
                charge_inj[card][ch].push_back(charges[Q_num]);
                charge_inj_err[card][ch].push_back(0);
            }
            labels[2]   = TString("Graph of Charge collected vs. Charge injected for FEB no." + card_num[card] + ", channel ") + TString(Form("%d",ch+1));
            gdata[card][ch]   = {charge_inj[card][ch], charge_inj_err[card][ch], mean[card][ch], mean_err[card][ch]};
            mygraph[card][ch] = TGrErr_getter(gdata[card][ch], ddata, idata, labels); 
            mygraph[card][ch]->Draw("AP");
            fit_pol3[card][ch] = new TF1 ("fit_pol3", Q_inj_fit, 0, 10, 5);
            fit_pol3[card][ch]->SetParameters(0., 0., 50., 1., 2.);
            fit_pol3[card][ch]->FixParameter(0,0);
            fit_pol3[card][ch]->FixParameter(4,2);
            gStyle->SetOptStat(111);
            gStyle->SetOptFit(111);
            gStyle->SetStatX(0.5);
            gStyle->SetStatY(0.9);
            fit_pol3[card][ch]->SetNpx(30000);
            c1->Update();
            fit_pol3[card][ch]->SetParNames("a_{0}", "a_{00}", "a_{1}", "a_{5}", "b"); 
            for(int x=0; x<10; ++x)
            {
                mygraph[card][ch]->Fit("fit_pol3","RqE");
                a0 = fit_pol3[card][ch]->GetParameter(0);   a00 = fit_pol3[card][ch]->GetParameter(1);      a1 = fit_pol3[card][ch]->GetParameter(2);   a5 = fit_pol3[card][ch]->GetParameter(3);   b = fit_pol3[card][ch]->GetParameter(4);
                fit_pol3[card][ch]->SetParameter(0,a0);     fit_pol3[card][ch]->SetParameter(1,a00);        fit_pol3[card][ch]->SetParameter(2,a1);     fit_pol3[card][ch]->SetParameter(3,a5);     fit_pol3[card][ch]->SetParameter(4,b);
                cout<<"Chisq/NDOF = "<<fit_pol3[card][ch]->GetChisquare()/fit_pol3[card][ch]->GetNDF()<<endl;
                if( (gMinuit->fCstatu.Contains("CONVERGED")) || (gMinuit->fCstatu.Contains("SUCCESSFUL"))  )  {fit_pol3[card][ch]->SetLineColor(3);  break;}
            }
            fit_pol3[card][ch]->Draw("SAME");
            ff<<fit_pol3[card][ch]->GetParameter(0)<<" "<<fit_pol3[card][ch]->GetParameter(1)<<"  "<<fit_pol3[card][ch]->GetParameter(2)<<"  "<<fit_pol3[card][ch]->GetParameter(3)<<"  "<<fit_pol3[card][ch]->GetParameter(4)<<endl;
            a0 = fit_pol3[card][ch]->GetParameter(0);   a00 = fit_pol3[card][ch]->GetParameter(1);      a1 = fit_pol3[card][ch]->GetParameter(2);   a5 = fit_pol3[card][ch]->GetParameter(3);   b = fit_pol3[card][ch]->GetParameter(4);
            a4 = a1*b*(exp(a5)-1)/(a5*(a0-a00+a1*b));
            a3 = a5*pow(b, -a4);
            a2 = a1*b*exp(a5)/(a4*a5); cout<<a2<<endl;
            g_hist->Fill(a1);
            a2_hist->Fill(a2);
            a3_hist->Fill(a3);
            c1->Print( TString("lin_test_TB") + TString(run_time) + TString(".pdf") );
            outfile->cd();
            mygraph[card][ch]->Write();
        }
        filename[card]->Close();
    }
    g_hist->Draw("hist");
    c1->Print( TString("lin_test_TB") + TString(run_time) + TString(".pdf") );
    a2_hist->Draw("hist");
    c1->Print( TString("lin_test_TB") + TString(run_time) + TString(".pdf") );
    a3_hist->Draw("hist");
    c1->Print( TString("lin_test_TB") + TString(run_time) + TString(".pdf") );
    outfile->cd();
    g_hist->Write();
    a2_hist->Write();
    a3_hist->Write();
    c1->Print( TString("lin_test_TB") + TString(run_time) + TString(".pdf]") );
    outfile->Close();
    delete outfile;
}