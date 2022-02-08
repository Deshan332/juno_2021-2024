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
#include "Q_function.cpp"

using namespace std;
int TB_Scurve_fit()
{
    TCanvas* c1         = new TCanvas("c1");
    Int_t p=0, FEB_id=601, FEB_ID, FEB_tmp = 0, cat_id, test_id, test_ID, cat_ID, data_uid, data_UID, ch, CH, ch_tmp=0, ch_inj, FEB_count=0, 
    ch_INJ, sid, SID, DAC, counts, counts_max=0, cats = 7,  entrycount=0, temp=0, s_temp=0, i=0, fcount = 0, n=0, f_id=1, status, v_size=0,
    cat_tmp=0, data_tmp=0, test_tmp=0, CH_tmp=0, ch_INJ_tmp=0, THRESH_tmp=0, SID_tmp=0, cts=0, counts_tmp=0, spikes=0, diffs=0, date, time;
    Double_t thresh, THRESH, ChiSq, UL, LL, N_0, N_0_err, mu, mu_err, sigma, sig_err;
    TString filename = "Scurve_tree_datetime.root", tree_name = "TB_Scurve";
    TFile *myfile           = new TFile(filename, "READ"); //input TFile
    TTree *TB_Scurve        = (TTree*)myfile->Get(tree_name);
    TB_Scurve->SetBranchAddress("FEB_id",       &FEB_id);
    TB_Scurve->SetBranchAddress("cat_id",       &cat_id);
    TB_Scurve->SetBranchAddress("data_uid",     &data_uid); 
    TB_Scurve->SetBranchAddress("test_id",      &test_id); 
    TB_Scurve->SetBranchAddress("ch",           &ch); 
    TB_Scurve->SetBranchAddress("ch_inj",       &ch_inj); 
    TB_Scurve->SetBranchAddress("thresh",       &thresh); 
    TB_Scurve->SetBranchAddress("sid",          &sid);
    TB_Scurve->SetBranchAddress("DAC",          &DAC); 
    TB_Scurve->SetBranchAddress("counts",       &counts);
    TB_Scurve->SetBranchAddress("date",         &date); 
    TB_Scurve->SetBranchAddress("time",         &time);
    
    UL = TB_Scurve->GetMaximum("DAC");
    LL = TB_Scurve->GetMinimum("DAC");
    TFile *outfile[30]; 
    TFile *treefile     = new TFile("TB_Scurve_ThirdPE_pars_datetime_2.root", "RECREATE"); 
    TTree *TB_pars      = new TTree("TB_Scurve_par_datetime","TestBenchDB Scurve NOPE fit parameters");
    
    TB_pars->Branch("FEB_ID",       &FEB_ID,    "FEB_ID/I");
    TB_pars->Branch("cat_ID",       &cat_ID,    "cat_ID/I");
    TB_pars->Branch("data_UID",     &data_UID,  "data_UID/I");
    TB_pars->Branch("test_ID",      &test_ID,   "test_ID/I"); 
    TB_pars->Branch("CH",           &CH,        "CH/I"); 
    TB_pars->Branch("SID",          &SID,       "SID/I");
    TB_pars->Branch("N_0",          &N_0,       "N_0/D"); 
    TB_pars->Branch("mu",           &mu,        "mu/D");
    TB_pars->Branch("ch_INJ",       &ch_INJ,    "ch_INJ/I"); //if ch_inj = -1, it corresponds to Scurve_Noinj. else Scurve_ThirdPE.
    TB_pars->Branch("THRESH",       &THRESH,    "THRESH/D");
    TB_pars->Branch("sigma",        &sigma,     "sigma/D"); 
    TB_pars->Branch("N_0_err",      &N_0_err,   "N_0_err/D");
    TB_pars->Branch("mu_err",       &mu_err,    "mu_err/D"); 
    TB_pars->Branch("sig_err",      &sig_err,   "sig_err/D");
    TB_pars->Branch("ChiSq",        &ChiSq,     "ChiSq/D");  
    TB_pars->Branch("status",       &status,    "status/I"); 
    TB_pars->Branch("spikes",       &spikes,    "spikes/I");
    TB_pars->Branch("diffs",        &diffs,     "diffs/I");  


    vector<Double_t> hit_counts, hit_counts_err, sc_DAC, sc_DAC_err; //plot is counts vs. DAC
    TGraphErrors *mygraph;
    TString labels[3]       = {TString("Threshold (DAC)"), TString("Counts (Normalized to unity)"), TString("")};        //array storing labels of graphs  xlabel, ylabel, title
    Double_t ddata[5]       = {1, LL, UL, 1.1, -0.1};       // array storing double values needed for graphs        //{marker size, x_LL, x_UL, y_max, y_min}
    Int_t idata[4]          = {4, 20, 2, 7};                    // array storing int values needed for graphs           //{marker color, marker size, linewidth, line color}
    TString cat_label[8]    = {TString("For TT"), TString("Good spares"), TString("Medium spares"), TString("Bad Spares"), TString("Broken cards"), TString("Non existent cards"), TString("For Telescope"), TString("Other") };
    
    for(int i=0; i<TB_Scurve->GetEntries(); ++i)
    {
        TB_Scurve->GetEntry(i);
        if((cat_id==1)&&(ch_inj>-1)&&(FEB_id>600))
        {
            if((date<20200518)&&(time<151500))  counts_max = 5000;
            else if(date<20200525)              counts_max = 5130;
            else                                counts_max = 5170;
            FEB_tmp = FEB_ID;   cat_tmp = cat_ID;   data_tmp = data_UID;    test_tmp = test_ID; CH_tmp = CH;    SID_tmp = SID;  THRESH_tmp = THRESH;    ch_INJ_tmp = ch_INJ;  //this set holds the data from previous read                           
            if(p==0){s_temp=sid; FEB_tmp = FEB_id; p=10;}
            if(FEB_id!=temp)    
            {
                ++FEB_count;  //FEB_count keeps the count of concerned FEB type only (per batch)
                temp = FEB_id;
                cout<<"FEB_id = "<<FEB_id<<endl;
            } 

            if((sid!=s_temp)&&(entrycount==1))    //before switching to a new scurve, do the fit for the current one 
            {
                v_size = hit_counts.size();
                cat_ID =  cat_tmp;
                labels[2] = TString("Graph of no. of TRT counts vs. DAC threshold for FEB no. ") + TString(Form("%d",FEB_tmp)) + TString(", channel ") + TString(Form("%d",CH_tmp)) + TString(" (Injecting 0PE, Category: ") + cat_label[cat_ID-1] + TString(")");
                    
                for(int m=0; m<v_size; ++m) 
                {
                    //hit_counts[m] = hit_counts[m]/counts_max;    //normallize counts to unity 
                    hit_counts_err.push_back ( sqrt ( ( (hit_counts[m]+1)*(hit_counts[m]+2))/((counts_max+2)*(counts_max+3) ) - ( pow(hit_counts[m]+1 , 2)/pow(counts_max+2 , 2) ) ) );
                    //hit_counts_err.push_back( sqrt( hit_counts[m] * (1 - hit_counts[m])/counts_max ) );
                    hit_counts[m] = hit_counts[m]/counts_max;           // now normalizing to the highest global TRT
                }

                struct graph_data gdata = {sc_DAC, sc_DAC_err, hit_counts, hit_counts_err};         //the graph
                Double_t par[3]         = {1., 360., 1.};       //initial parametrization
                Double_t par_limits[6]  = {0.95*par[0], 1.05*par[0], LL, UL, 0., 10.}; //parameter bounds
                mygraph = TGrErr_getter(gdata, ddata, idata, labels);

                mygraph->SetName(TString("FEB") + TString(Form("%d", FEB_tmp)) + TString("ch") + TString(Form("%d",CH_tmp)));
                mygraph->Draw("AP");
                    
                TF1 *Q_fit = new TF1 ("Q_fit", Q_function, LL, UL,3);                       //fit method
                TString parnames[3] = {TString("N_{0}"), TString("#mu"), TString("#sigma_{0}")};

                for(int k=0; k<3; ++k)
                {
                    Q_fit->SetParameter(k,par[k]);
                    Q_fit->SetParLimits(k,par_limits[2*k],par_limits[2*k+1]);
                    Q_fit->SetParName(k,parnames[k]);  
                }
                Q_fit->SetNpx(30000);
                for(int p=0; p<10; ++p)
                {
                    Q_fit->FixParameter(0, 1); //commented this line after normalizing to highest global TRT
                    mygraph->Fit("Q_fit","RqE","",LL,UL);

                    for(int j=0; j<3; ++j)
                    {
                        par[j] = Q_fit->GetParameter(j);
                        Q_fit->SetParameter(j,par[j]);
                    }
                    if( (gMinuit->fCstatu.Contains("CONVERGED"))  && (Q_fit->GetChisquare()/Q_fit->GetNDF() < 5) )  {status =1;    Q_fit->SetLineColor(2);   break;}
                    else    {status =0;    Q_fit->SetLineColor(1);}
                }
                gStyle->SetOptStat(111);
                gStyle->SetOptFit(111);
                Q_fit->Draw("SAME");
                N_0     = Q_fit->GetParameter(0);       mu      = Q_fit->GetParameter(1);   sigma   = Q_fit->GetParameter(2);
                N_0_err = Q_fit->GetParError(0);        mu_err  = Q_fit->GetParError(1);    sig_err = Q_fit->GetParError(2);
                ChiSq   = Q_fit->GetChisquare()/Q_fit->GetNDF();
                FEB_ID = FEB_tmp;    cat_ID = cat_tmp;    data_UID = data_tmp;    test_ID = test_tmp;  CH = CH_tmp;        SID=SID_tmp;        THRESH = THRESH_tmp;        ch_INJ = ch_INJ_tmp; //set the correct values from memory to write to TTree
                hit_counts.clear();    
                sc_DAC.clear();    
                treefile->cd();
                TB_pars->Fill();
                outfile[n]->cd();
                mygraph->Write();
                counts_max = counts;
                s_temp=sid;
                FEB_ID = FEB_id;    cat_ID = cat_id;    data_UID = data_uid;    test_ID = test_id;  CH = ch;    SID=sid;    THRESH = thresh;    ch_INJ = ch_inj; //assign values from current set 
            }
            if(cts==0) //save first line to memory 
                {
                    FEB_ID = FEB_id;    cat_ID = cat_id;    data_UID = data_uid;    test_ID = test_id;  CH = ch;    SID = sid;    THRESH = thresh;    ch_INJ = ch_inj;   cts=1;
                }  
            if(ch==0) spikes=0;
            s_temp=sid; entrycount=1;
            if(counts<=counts_max)
            {
                hit_counts.push_back(counts); //y 
                //cout<<counts<< "  "<<counts_tmp<<endl;
                if((DAC<452) && (DAC>448) && (counts>counts_tmp)) 
                {
                    spikes=DAC; diffs = counts-counts_tmp;
                    cout<<"Spike found at FEB no. "<<FEB_id<<", channel no. "<<ch<<" at "<<DAC<<" DAC"<< endl;
                }
                counts_tmp = counts;
                //hit_counts_err.push_back(0); //y_err 
                sc_DAC.push_back(DAC); //x in DAC
                sc_DAC_err.push_back(0); //x_err
            }
            //if(counts>counts_max) counts_max = counts;  
        }

        if(FEB_count==100)
        {
            f_id=1;
            outfile[n]->Close();
            delete outfile[n];
            cout<<"TFile closed."<<endl;
            ++n;  
            FEB_count=0;  
        }
        if(f_id==1) //f_id is the TFile creation index. f_id =1 ---> call to Create a new TFile
        {
            TString file_tag    = TString("TB_Scurve_ThirdPE_") + TString(Form("%d",FEB_id)) + TString(".root");
            outfile[n]          = new TFile(file_tag, "RECREATE"); //n is the index of TFile count.
            f_id = 0;
            cout<<"New TFile created. (n = "<<n<<")"<<endl;
        }

    }

    //to get the last fit


    v_size = hit_counts.size();
    cat_ID =  cat_tmp;
    labels[2] = TString("Graph of no. of TRT counts vs. DAC threshold for FEB no. ") + TString(Form("%d",FEB_tmp)) + TString(", channel ") + TString(Form("%d",CH_tmp)) + TString(" (Injecting 0PE, Category: ") + cat_label[cat_ID-1] + TString(")");
                    
    for(int m=0; m<v_size; ++m) 
    {
        hit_counts[m] = hit_counts[m]/counts_max;    //normallize counts to unity
    }

    struct graph_data gdata1 = {sc_DAC, sc_DAC_err, hit_counts, hit_counts_err};         //the graph
    Double_t par1[3]         = {1., 360., 1.};       //initial parametrization
    Double_t par_limits1[6]  = {0.95*par1[0], 1.05*par1[0], LL, UL, 0., 10.}; //parameter bounds
    mygraph = TGrErr_getter(gdata1, ddata, idata, labels);

    mygraph->SetName(TString("FEB") + TString(Form("%d", FEB_tmp)) + TString("ch") + TString(Form("%d",CH_tmp)));
    mygraph->Draw("AP");
                    
    TF1 *Q_fit1 = new TF1 ("Q_fit", Q_function, LL, UL,3);                       //fit method
    TString parnames1[3] = {TString("N_{0}"), TString("#mu"), TString("#sigma_{0}")};

    for(int k=0; k<3; ++k)
    {
        Q_fit1->SetParameter(k,par1[k]);
        Q_fit1->SetParLimits(k,par_limits1[2*k],par_limits1[2*k+1]);
        Q_fit1->SetParName(k,parnames1[k]);  
    }
        Q_fit1->SetNpx(30000);
        for(int p=0; p<10; ++p)
        {
            Q_fit1->FixParameter(0, 1);
            mygraph->Fit("Q_fit","RqE","",LL,UL);

            for(int j=0; j<3; ++j)
            {
                par1[j] = Q_fit1->GetParameter(j);
                Q_fit1->SetParameter(j,par1[j]);
            }
            if( (gMinuit->fCstatu.Contains("CONVERGED"))  && (Q_fit1->GetChisquare()/Q_fit1->GetNDF() < 5) )  {status =1;    Q_fit1->SetLineColor(2);   break;}
            else    {status =0;    Q_fit1->SetLineColor(1);}
        }
        gStyle->SetOptStat(111);
        gStyle->SetOptFit(111);
        Q_fit1->Draw("SAME");
        if(ch>62) cout<<"fitted FEB "<<FEB_tmp<<" channel "<<CH_tmp<<endl; 
        N_0     = Q_fit1->GetParameter(0);       mu      = Q_fit1->GetParameter(1);   sigma   = Q_fit1->GetParameter(2);
        N_0_err = Q_fit1->GetParError(0);        mu_err  = Q_fit1->GetParError(1);    sig_err = Q_fit1->GetParError(2);
        ChiSq   = Q_fit1->GetChisquare()/Q_fit1->GetNDF();
        FEB_ID = FEB_tmp;    cat_ID = cat_tmp;    data_UID = data_tmp;    test_ID = test_tmp;  CH = CH_tmp;        SID=SID_tmp;        THRESH = THRESH_tmp;        ch_INJ = ch_INJ_tmp; //set the correct values from memory to write to TTree
        hit_counts.clear();    
        sc_DAC.clear();    
        treefile->cd();
        TB_pars->Fill();
        outfile[n]->cd();
        mygraph->Write();
        counts_max = counts;
        s_temp=sid;


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

