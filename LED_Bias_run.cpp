#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <time.h>
#include <TH1F.h>
#include <TF1.h>
#include <math.h>
#include <TMath.h>
#include<TGraph.h>
#include<TTree.h>
#include<TGraphErrors.h>
#include<THistPainter.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TMinuit.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TText.h>
#include <TPaveStats.h>
#include "hist_getter.cpp"
#include "Gaus_fit.cpp"
#include "Bellamy.cpp"
#include "Ped_fitter.cpp"
#include "Led_fitter.cpp"
#include "TGrErr_getter.cpp"
#include "timestamp.cpp"
using namespace std;
void LED_Bias_run()
{
    TCanvas* c1 = new TCanvas("c1");
    TString run_time = timestamp();
    TFile *outfile  = new TFile(TString("LB_plots_") + TString(run_time) + TString(".root"), "RECREATE");
    ofstream ff ("gaindata_norm.txt");
    c1->SetWindowSize(1400,800);     
    c1->SetCanvasSize(1400,800); 
    c1->Update();
    c1->SetLogy();
    gStyle->SetStatW(0.1);
    int B=1;   //B= no. of led bias values. Change when analyzing different numbers of led biases
    TString histname, filename_LED, filename_PED, handedness;
    TString SUCCESS, STATUS, BV;
    TString cname;
    Double_t UL = 250;          // Upper limit of charge fit range
    Double_t LL = 27.5;         // lower limit of charge fit range
    TString time_PED        =   TString("2021-03-24_16_53");        //insert initial ped timestamp here. edit for different data.
    TString time_PED_last   =   TString("2021-03-24_16_53");        //insert final ped timestamp here. edit for different data.
    TString time_LED[12]    =   {TString("2021-03-24_16_35"),TString("2021-03-24_16:35"), TString("2021-03-24_16:26"), TString("2021-03-24_16:16"), TString("2021-03-24_16:07"), TString("2021-03-24_15:57"), TString("2021-03-24_14:52"), TString("2021-03-24_15:01"), TString("2021-03-24_15:15"), TString("2021-03-24_15:24"), TString("2021-03-24_15:37"), TString("2021-03-24_15:47")}; //timestamps of led functions depending on LED intensity in ascending order. edit depending on changes in B.
    TH1F* hist_PED[B][4][2][64];               //hist_PED[bias][mod][hand][ch]  --- PED histos
    TH1F* hist_LED[B][4][2][64];               //hist_LED[bias][mod][hand][ch]  --- LED histos
    TF1* Ped_fit[B][4][2][64];                 //Ped_fit[bias][mod][hand][ch]   --- PED fits
    TF1* Led_fit[B][4][2][64];                                           //---LED fits
    Double_t par_ped[10];                       //arrays storing ped fit details
    Double_t par_ped_err[10];
    Double_t par_led[10];
    Double_t par_led_err[10], N_ped=0;
    int iterations=0;
    vector<Double_t> mu_arr[B][4][2], mu_err_arr[12][4][2], Q_1_arr[12][4][2], Q_1_err_arr[12][4][2], Chisq_arr[12][4][2], chan[12][4][2], Chisq_err_arr[12][4][2], chan_err[12][4][2];  //vectors used to plot graphs for every card
    TGraphErrors* MU[B][4][2];         //the TGraphErrors definitions 
    TGraphErrors* GAIN[B][4][2];
    TGraphErrors* CHI[B][4][2];
    TFile *f = new TFile( TString("LB_Tree_") + TString(run_time) + TString(".root"), "RECREATE");       //TFile saving TTree
    TString parnames[9]  = {TString("N_{0}"),    TString("Q_{0}"),       TString("Q_{1}"),       TString("#sigma_{0}"),      TString("#sigma_{1}"),      TString("w"),   TString("#alpha"),  TString("#mu"),     TString("#Chi^{2}/NDOF") };  //names of parameters fed to the stat box
    TString par_axis[9]  = {TString("N_{0}"),    TString("Q_{0} (ADC)"), TString("Q_{1} (ADC)"), TString("#sigma_{0} (ADC)"),TString("#sigma_{1} (ADC)"),TString("w"),   TString("#alpha"),  TString("#mu"),     TString("#Chi^{2}/NDOF") };  // y-axis labels
    TString LED_Bias[12] = {TString("-400"),        TString("-440"),        TString("-480"),            TString("-520"),            TString("-560"),TString("-600"),    TString("-640"),    TString("-680"),            TString("-720"),    TString("-760"),    TString("-800"), TString("")};  //led biases. used to make graph titles and keep track of running led bias when running. edit depending on changes in B.
    TString title;      //title of graphs
    c1->Print( TString("LB_analysis_") + TString(run_time) + TString(".pdf[") );     //pdf storing all fits and graphs
    TTree *par_tree = new TTree("par_tree","Fit_Cards");        //TTree
    Int_t status=0, module, LRhand, channel, range_passed, bin_max, BIAS, ped_val, ped_stat=0;//TBranches
    Double_t N_0, Q_0, Q_1, s_0, s_1, w, alpha, mu, Chisq, Q_0_err, Q_1_err, s_0_err, s_1_err, w_err, alpha_err, mu_err, bin_max_x;
    par_tree->Branch("BIAS",    &BIAS,      "BIAS/I");
    par_tree->Branch("module",  &module,    "module/I");
    par_tree->Branch("LRhand",  &LRhand,    "LRhand/I"); 
    par_tree->Branch("channel", &channel,   "channel/I"); 
    par_tree->Branch("N_0",     &N_0,       "N_0/D"); 
    par_tree->Branch("Q_0",     &Q_0,       "Q_0/D");
    par_tree->Branch("Q_1",     &Q_1,       "Q_1/D"); 
    par_tree->Branch("s_0",     &s_0,       "s_0/D"); 
    par_tree->Branch("s_1",     &s_1,       "s_1/D"); 
    par_tree->Branch("w",       &w,         "w/D"); 
    par_tree->Branch("alpha",   &alpha,     "alpha/D"); 
    par_tree->Branch("mu",      &mu,        "mu/D");
    par_tree->Branch("Chisq",   &Chisq,     "Chisq/D");
    par_tree->Branch("Q_0_err",     &Q_0_err,       "Q_0_err/D");
    par_tree->Branch("Q_1_err",     &Q_1_err,       "Q_1_err/D"); 
    par_tree->Branch("s_0_err",     &s_0_err,       "s_0_err/D"); 
    par_tree->Branch("s_1_err",     &s_1_err,       "s_1_err/D"); 
    par_tree->Branch("w_err",       &w_err,         "w_err/D"); 
    par_tree->Branch("alpha_err",   &alpha_err,     "alpha_err/D"); 
    par_tree->Branch("mu_err",      &mu_err,        "mu_err/D");
    par_tree->Branch("status",      &status,        "status/I");   //status of convergence  
    par_tree->Branch("range_passed",      &range_passed,        "range_passed/I");  //keeps track of whether the LL of fit range exceeds ped position. 1 if not, 0 if passed.
    for(int bias=0; bias<B; ++bias)
    {
        for(int mod=0; mod<1; ++mod)            //start fitting in a loop, mod for module {0,1,2,3}
        {
            for(int hand=0; hand<2; ++hand)     //hand for L or R
            {            
                if( (mod==0) && (hand==0) ) {continue;} //skip 0R. Remove when analyzing all cards.
                //if( (mod==1) && (hand==1) ) {continue;} //skip 1R. Remove when analyzing all cards.
                for(int ch=0; ch<64; ++ch)      //ch for channel 1-64
                {
                    if(hand==1){handedness=TString("L");} 
                    if(hand==0){handedness=TString("R");}
                    c1->SetRightMargin(0.2);
                    filename_PED = TString("ped_c") + TString(Form("%d",mod)) + handedness + TString("_") + time_PED + TString("_hist.root");   //ped histo file
                    filename_LED = TString("led") + LED_Bias[bias] +TString("_c") + TString(Form("%d",mod)) + handedness + TString("_") + time_LED[bias] + TString("_hist.root");       //led histo file
                    hist_PED[bias][mod][hand][ch] = hist_getter(TString("/home/pdeshans/Documents/M2/led_data_new/") + filename_PED, ch+1);       //retrieve histos from file using hist_getter
                    hist_LED[bias][mod][hand][ch] = hist_getter(TString("/home/pdeshans/Documents/M2/led_data_new/") + filename_LED, ch+1);
                    hist_LED[bias][mod][hand][ch]->GetXaxis()->SetRangeUser(0,UL);
                    c1->Update();
                    par_ped[0] = hist_PED[bias][mod][hand][ch]->Integral();     //initializations for ped fit
                    par_ped[1] = hist_PED[bias][mod][hand][ch]->GetMean();
                    par_ped[2] = hist_PED[bias][mod][hand][ch]->GetRMS();
                    cout<<"***************************"<<endl;
                    cname = TString("Analyzing Card c") + Form("%d",mod) + handedness + TString(", channel ") + TString(Form("%d",ch+1)) + TString(", LED bias = ") + LED_Bias[bias] + TString("mV");
                    if(B==1){cname = TString("Analyzing Card c") + Form("%d",mod) + handedness + TString(", channel ") + TString(Form("%d",ch+1)) ;}
                    cname.ReplaceAll("-","");
                    cout<<cname<<endl;
                    cout<<"Initializations:"<<endl;
                    cout<<"N_0       Q_0       s_0"<<endl;
                    for(int i=0; i<3; ++i)  {cout<<par_ped[i]<<"   ";}  //print initializations
                    cout<<"  "<<endl;
                    cout<<"***************************"<<endl;
                    iterations=0;
                    gStyle->SetOptStat(1111);
                    gStyle->SetOptFit(1);
                    c1->SetLogy();
                    c1->SetGridx(0);
                    c1->SetGridy(0);
                    title = TString("Card ") + TString(Form("%d",mod)) + handedness + TString(", LED Bias = ")  +  LED_Bias[bias] + TString("mV");
                    if(B==1){title = TString("card C") + TString(Form("%d",mod)) + handedness ;}
                    title.ReplaceAll("-","");
                    hist_LED[bias][mod][hand][ch]->SetTitle(title);
                    gStyle->SetStatX(0.98);
                    gStyle->SetStatY(0.8);
                    hist_LED[bias][mod][hand][ch]->Draw("HIST");
                    cout<<"Fitting the PED file "<<filename_PED<<" ..."<<endl;
                    Ped_fit[bias][mod][hand][ch] = Ped_fitter(hist_PED[bias][mod][hand][ch], par_ped, par_ped_err);    //the first PED fit, using initial ped file
                    for(int j=0; j<100; ++j) //repeat ped fit 100 times until converged
                    {
                        if( (gMinuit->fCstatu.Contains("CONVERGED")) && (Ped_fit[bias][mod][hand][ch]->GetChisquare()/Ped_fit[bias][mod][hand][ch]->GetNDF() < 5) )  {Ped_fit[bias][mod][hand][ch]->SetLineColor(3);    ped_stat=1; break;}
                        else {Ped_fit[bias][mod][hand][ch] = Ped_fitter(hist_PED[bias][mod][hand][ch], par_ped, par_ped_err);   Ped_fit[bias][mod][hand][ch]->SetLineColor(7);}
                        iterations=j;
                    }
                    if(ped_stat==0)     //if ped did not converge, try the other ped file for 100 iterations max
                    {
                        filename_PED = TString("ped_c") + TString(Form("%d",mod)) + handedness + TString("_") + time_PED_last + TString("_hist.root");
                        cout<<"Fitting the PED file "<<filename_PED<<" ..."<<endl;
                        hist_PED[bias][mod][hand][ch] = hist_getter(TString("/home/pdeshans/Documents/M2/led_data_new/") + filename_PED, ch+1);
                        for(int j=0; j<100; ++j) //repeat ped fit 100 times until converged
                        {
                            if( (gMinuit->fCstatu.Contains("CONVERGED")) && (Ped_fit[bias][mod][hand][ch]->GetChisquare()/Ped_fit[bias][mod][hand][ch]->GetNDF() < 5) )  {Ped_fit[bias][mod][hand][ch]->SetLineColor(3);    ped_stat=1; break;}
                            else {Ped_fit[bias][mod][hand][ch] = Ped_fitter(hist_PED[bias][mod][hand][ch], par_ped, par_ped_err);   Ped_fit[bias][mod][hand][ch]->SetLineColor(7);}
                            iterations=j;
                        }
                    } 
                    cout<<"PED_fit status: "<< gMinuit->fCstatu<<" after "<<iterations+1<<" iterations."<<endl;
                    cout<<"Fit Results:"<<endl;
                    cout<<"N_0            Q_0            s_0"<<endl;
                    cout<<par_ped[0]<<"   "<<par_ped[1]<<"   "<<par_ped[2]<<endl;       //give ped fit stats
                    iterations=0;
                    par_led[0] = hist_LED[bias][mod][hand][ch]->Integral();           //N_0   //initialization of led fit parameters
                    par_led[1] = par_ped[1];                                    //Q_0 (from pedestal fit)
                    par_led[2] = 7.;    //hist1->GetMean();                     //Q_1 (initialized to the mean of the led distribution)
                    par_led[3] = par_ped[2];                                    //s_0 (from pedestal fit)
                    par_led[4] = 5;    //hist1->GetRMS();                      //s_1 (initialized to SD of led distribution)
                    par_led[5] = 0.05;                                          //w 
                    par_led[6] = 0.1;                                            //alpha
                    for(int y=0; y<3; ++y)      //calculate no. of entries in 3 bins surrounding the pedestal ---- 1st estimate of mu
                    {
                        N_ped += hist_LED[bias][mod][hand][ch]->GetBinContent(hist_LED[bias][mod][hand][ch]->FindBin(par_led[1])-1 + y);
                    }
                    par_led[7] = TMath::Log(par_led[0]/N_ped);              //mu   = ln (Nped/Ntot)
                    if(par_led[7]<0.1){par_led[7]=0.5;}                     // if mu<0.1 reassign to 0.5
                    Double_t par_limits[16]={0.95*par_led[0],1.05*par_led[0],par_ped[1]-3,par_ped[1]+3,0.,30.,0.,2.,0.,12.,0.,1.,0.,2.,0.,4}; //limits for led fit parameters. [0] & [1] for p0, [2] & [3] for p1 etc.
                    if(bias>6){par_limits[15] = 6.5;}       //lift the mu limit for high LED biases
                    cout<<"************************************************************"<<endl;
                    cout<<"Initializations:"<<endl;
                    cout<<"N_0       Q_0       Q_1   s_0     s_1   w   alpha   mu"<<endl;
                    for(int i=0; i<8; ++i)  {cout<<par_led[i]<<"   ";}
                    cout<<"  "<<endl;
                    cout<<"************************************************************"<<endl;
                    ped_val = hist_LED[bias][mod][hand][ch]->GetXaxis()->FindBin(par_ped[1]);           //find the bin containing pedestal
                    status=0;       //parameter used to check if fit converged and chisq/ndof < 5. 0 if FALSE, 1 if TRUE
                    LL = ped_val - 3;       //set lower limit of fit 3 bins below the pedestal
                    for(int x=0; x<14; ++x) //fit the LED for 14 rounds shifting the LL up by 0.5 ADC each round until convergence with chisq/ndf < 5 is obtained
                    {
                        LL=LL+0.5;
                        cout<<"Fitting on the range ["<<LL<<","<<UL<<"] :"<<endl;
                        for(int k=0; k<20; ++k) //led fit 20 times until converged or chisq/NDOF < 5
                        {
                            par_led[2] = 7 + 0.5*k;
                            iterations = iterations + 1;
                            Led_fit[bias][mod][hand][ch] = Led_fitter(hist_LED[bias][mod][hand][ch], par_led, par_limits, par_led_err, LL, UL);     //fit to bellamy using Led_fitter
                            cout<<"Chisq/NDOF for LED_fit = "<<Led_fit[bias][mod][hand][ch]->GetChisquare()/Led_fit[bias][mod][hand][ch]->GetNDF()<<" , "<<gMinuit->fCstatu<<endl;
                            if( ( (gMinuit->fCstatu.Contains("CONVERGED")) || (gMinuit->fCstatu.Contains("SUCCESSFUL")) ) && (Led_fit[bias][mod][hand][ch]->GetChisquare()/Led_fit[bias][mod][hand][ch]->GetNDF()<3) )  {break;}
                        }
                        if( ( (gMinuit->fCstatu.Contains("CONVERGED")) || (gMinuit->fCstatu.Contains("SUCCESSFUL")) ) && (Led_fit[bias][mod][hand][ch]->GetChisquare()/Led_fit[bias][mod][hand][ch]->GetNDF()<3) )  {break;}
                    }
                    if ((gMinuit->fCstatu.Contains("CONVERGED")) || (gMinuit->fCstatu.Contains("SUCCESSFUL")) ){Led_fit[bias][mod][hand][ch]->SetLineColor(3);   status=1; } //set color and status value depending on fit status
                    else {Led_fit[bias][mod][hand][ch]->SetLineColor(2);  status=0; }
                    cout<<"LED_fit status: "<< gMinuit->fCstatu<<" after "<<iterations<<" iterations."<<endl;
                    cout<<"************************************************************"<<endl;
                    cout<<"Fit Results:"<<endl;
                    cout<<"Chisq/NDOF for LED_fit = "<<Led_fit[bias][mod][hand][ch]->GetChisquare()/Led_fit[bias][mod][hand][ch]->GetNDF()<<endl;
                    cout<<"N_0       Q_0       Q_1      s_0      s_1      w      alpha       mu"<<endl;
                    for(int i=0; i<8; ++i)  {cout<<par_led[i]<<"   ";}          //print fit results
                    cout<<"  "<<endl;
                    gStyle->SetOptStat(1111);
                    gStyle->SetOptFit(1);
                    int bin=0;
                    f->cd();        //change directory to TFile we need to store the TTree
                    BIAS    =   bias;
                    module  =   mod;    //assign values to TBranch values
                    LRhand  =   hand;
                    channel =   ch+1;
                    N_0     =   par_led[0];
                    Q_0     =   par_led[1];
                    Q_1     =   par_led[2];     
                    s_0     =   par_led[3];
                    s_1     =   par_led[4];
                    w       =   par_led[5];
                    alpha   =   par_led[6];
                    mu      =   par_led[7];
                    Chisq   =   Led_fit[bias][mod][hand][ch]->GetChisquare()/Led_fit[bias][mod][hand][ch]->GetNDF();
                    Q_0_err =   par_led_err[1];
                    Q_1_err =   par_led_err[2];
                    s_0_err =   par_led_err[3]; 
                    s_1_err =   par_led_err[4];
                    w_err   =   par_led_err[5];
                    alpha_err   =   par_led_err[6];
                    mu_err  =   par_led_err[7];
                    if(LL>Q_0-0.1)  {range_passed=0;}   //if last LL is beyond the pedestal, mark that it has passed. 0 if TRUE.
                    else        {range_passed=1;}
                    if( (Chisq<5) && (status==1) )      // push values to vectors used to draw graphs only if converged and chisq/NDF < 5
                    {
                        mu_arr[bias][mod][hand].push_back(mu);
                        mu_err_arr[bias][mod][hand].push_back(mu_err);
                        Q_1_arr[bias][mod][hand].push_back(Q_1);
                        Q_1_err_arr[bias][mod][hand].push_back(Q_1_err);
                        Chisq_arr[bias][mod][hand].push_back(Chisq);
                        Chisq_err_arr[bias][mod][hand].push_back(0.);
                        chan[bias][mod][hand].push_back(ch+1);
                        chan_err[bias][mod][hand].push_back(0.);
                    }
                    ff << bias <<"  "<< mod << "  "<< hand << "  " << ch << "  " << par_led[2] << "  " << status <<endl;          
                    c1->Update();
                    c1->SetLogy();
                    Led_fit[bias][mod][hand][ch]->Draw("SAMES");    //draw led fit on histo
                    if(LL>Q_0-0.1)    // if range has passed, draw the ped from ped fit
                    {
                        cout<<"LED_fit converged beyond the Pedestal. Invoking results from PED fit..."<<endl;
                        cout<<"************************************************************"<<endl;
                        Ped_fit[bias][mod][hand][ch]->Draw("SAME");
                        //c1->Update();   //commenting this out stops the stat box from updating to new fit parameters
                        Q_0 =   par_ped[1];
                        s_0 =   par_ped[2];
                        Q_0_err =   par_ped_err[1];
                        s_0_err =   par_ped_err[2];
                        TPaveStats *ps = (TPaveStats*)c1->GetPrimitive("stats");        //update stat box with new stats from ped fit
                        ps->SetName("mystats");
                        TList *listoflines = ps->GetListOfLines();
                        TText *ped = ps->GetLineWith("Q_{0}");
                        TText *sig = ps->GetLineWith("#sigma_{0}");
                        listoflines->Remove(ped);
                        listoflines->Remove(sig);
                        TString Q = TString("Q_{0} = ") + TString(Form("%.3f", Q_0)) + TString(" +/- ") + TString(Form("%.3f",Q_0_err));
                        TString S = TString("#sigma_{0} = ") + TString(Form("%.5f", s_0)) + TString(" +/- ") + TString(Form("%.5f",s_0_err));
                        TLatex *PED = new TLatex(0,0,Q);
                        TLatex *SIG = new TLatex(0,0,S);
                        PED ->SetTextFont(42);
                        PED ->SetTextSize(0.025);
                        PED ->SetTextColor(kRed);
                        SIG ->SetTextFont(42);
                        SIG ->SetTextSize(0.025);
                        SIG ->SetTextColor(kRed);
                        listoflines->Add(PED);
                        listoflines->Add(SIG);
                        hist_LED[bias][mod][hand][ch]->SetStats(0);
                    }
                    par_tree->Fill();           //fill tree
                    c1->Update();
                    outfile->cd();
                    hist_LED[bias][mod][hand][ch]->Write();
                    //Led_fit[bias][mod][hand][ch]->Write();
                    c1->Print( TString("LB_analysis_") + TString(run_time) + TString(".pdf") );
                }
                BV = TString(" , LED Bias = ") + LED_Bias[bias] + TString("mV");    BV.ReplaceAll("-","");      //TString indicating LED bias in title
                TString labels[3]   = {TString("Channel"), TString("#mu"), TString("Graph of #mu vs. Channel for card ") + TString(Form("%d",mod)) + handedness + BV };        //array storing labels of graphs
                Double_t ddata[5]   = {1.2, 0., 65., 7., 0.};       // array storing double values needed for graphs        //{marker size, x_LL, x_UL, y_max, y_min}
                Int_t idata[4]      = {2, 20, 2, 1};                // array storing int values needed for graphs           //{marker color, marker size, linewidth, line color}
                struct graph_data mudata = {chan[bias][mod][hand], chan_err[bias][mod][hand], mu_arr[bias][mod][hand], mu_err_arr[bias][mod][hand]};        //structs needed to feed into TGrErr_getter
                MU[bias][mod][hand]     = TGrErr_getter(mudata, ddata, idata, labels);  //get TGraphErrors for mu vs. channel
                c1->SetLogy(0);
                c1->SetGridx();
                c1->SetGridy();
                c1->SetRightMargin(0.1);
                MU[bias][mod][hand]->Draw("ALP");
                outfile->cd();
                MU[bias][mod][hand]->Write();
                c1->Print( TString("LB_analysis_") + TString(run_time) + TString(".pdf") );

                labels[1]   = TString("Q_{1} (ADC)");   labels[2] = TString("Graph of Q_{1} vs. Channel for card C") + TString(Form("%d",mod)) + handedness + BV; 
                ddata[3]    = 25.; ddata[4] = 0;
                idata[0]    = 4;   idata[3] = 7;
                struct graph_data gaindata  = {chan[bias][mod][hand], chan_err[bias][mod][hand], Q_1_arr[bias][mod][hand], Q_1_err_arr[bias][mod][hand]};
                GAIN[bias][mod][hand]       = TGrErr_getter(gaindata, ddata, idata, labels); //get TGraphErrors for gain vs. channel
                c1->SetGridx();
                c1->SetGridy();
                GAIN[bias][mod][hand]->Draw("ALP");
                outfile->cd();
                GAIN[bias][mod][hand]->Write();
                c1->Print( TString("LB_analysis_") + TString(run_time) + TString(".pdf") );

                labels[1]   = TString("#Chi^{2}/NDF");   labels[2] = TString("Graph of #Chi^{2}/NDF vs. Channel for card C") + TString(Form("%d",mod)) + handedness + BV; 
                ddata[3]    = 5.;  ddata[4] = 0;
                idata[0]    = 1;   idata[3] = 3;
                struct graph_data chidata   = {chan[bias][mod][hand], chan_err[bias][mod][hand], Chisq_arr[bias][mod][hand], Chisq_err_arr[bias][mod][hand]};
                CHI[bias][mod][hand]        = TGrErr_getter(chidata, ddata, idata, labels); //get TGraphErrors for Chisq/NDF vs. channel
                c1->SetGridx();
                c1->SetGridy();
                CHI[bias][mod][hand]->Draw("ALP");
                outfile->cd();
                CHI[bias][mod][hand]->Write();
                c1->Print( TString("LB_analysis_") + TString(run_time) + TString(".pdf") );
            }
        }
    }
    c1->Print( TString("LB_analysis_") + TString(run_time) + TString(".pdf]") );
    outfile->cd();
    outfile->Write();
    outfile->Close();
    ff.close();
    delete outfile;
    f->cd();
    f->Write();
    f->Close();
    delete f;
}