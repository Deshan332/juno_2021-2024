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
#include <TTreeIndex.h>
using namespace std;

int Rank_cards()
{
    Int_t FEB_id, FEB_ID, FEB_idc, FEB_IDTh, FEB_idcTh, CH, CHTh, rank_ch=0, rank_card=0, temp_id=0, status;
    Double_t a1, a1_mean, a1_r, a2, a2_mean, a2_r, ChiSq, ChiSq_mean, 
    ChiSqTh, ChiSq_meanTh, mu, mu_mean, mu_r, sigma, sig_mean, sig_r, score_ch, score_card, scores_ch[1300], scores_card[1300];
    
    TString Qinj_file = "/home/pdeshans/Documents/PhD_files/from_M2_cpy/Qinj_fitdata/After_constrain/test_b2.root",     Qinj_file_tree = "TB_lin_par";
    TString Qinj_card = "/home/pdeshans/Documents/PhD_files/from_M2_cpy/Qinj_fitdata/After_constrain/TB_cards_b2.root", Qinj_card_tree = "TB_cards";
    TString ThPE_file = "/home/pdeshans/Documents/PhD_files/from_M2_cpy/ThirdPE_fitdata/corrected_for_variable_plateau_27_01_2022/Scurve_ThirdPE_sum.root",       ThPE_file_tree = "TB_Scurve_par_datetime";
    TString ThPE_card = "/home/pdeshans/Documents/PhD_files/from_M2_cpy/ThirdPE_fitdata/corrected_for_variable_plateau_27_01_2022/TB_cards_ThirdPE.root",         ThPE_card_tree = "TB_cards";
    
    TFile *Qinj_Tfile       = new TFile(Qinj_file); //input TFile with Qinj data by channels
    TTree *Qinj_pars        = (TTree*)Qinj_Tfile->Get(Qinj_file_tree);
    Qinj_pars->SetBranchAddress("FEB_ID",       &FEB_ID); 
    Qinj_pars->SetBranchAddress("CH",           &CH); 
    Qinj_pars->SetBranchAddress("a1",           &a1); 
    Qinj_pars->SetBranchAddress("a2",           &a2); 
    Qinj_pars->SetBranchAddress("ChiSq",        &ChiSq);
    Qinj_pars->SetBranchAddress("status",       &status);

    TFile *Qinj_Tcard       = new TFile(Qinj_card); //input Tfile with Qinj data by cards 
    TTree *Qinj_card_pars   = (TTree*)Qinj_Tcard->Get(Qinj_card_tree);
    Qinj_card_pars->SetBranchAddress("FEB_idc",      &FEB_idc);  
    Qinj_card_pars->SetBranchAddress("a1_mean",      &a1_mean);   
    Qinj_card_pars->SetBranchAddress("a2_mean",      &a2_mean);   
    Qinj_card_pars->SetBranchAddress("ChiSq_mean",   &ChiSq_mean); 

    Qinj_card_pars->Draw("a1_mean>>a1_his", "", "goff");
    TH1F *a1_his = (TH1F*)gDirectory->Get("a1_his");
    Double_t a1_g   = a1_his->GetMean(); //mean of means
    Qinj_card_pars->Draw("a2_mean>>a2_his", "", "goff");
    TH1F *a2_his = (TH1F*)gDirectory->Get("a2_his");
    Double_t a2_g   = a2_his->GetMean();    //mean of means
  
    TFile *ThPE_Tfile       = new TFile(ThPE_file); //input TFile with ThirdPE data by channels
    TTree *ThPE_pars        = (TTree*)ThPE_Tfile->Get(ThPE_file_tree);
    ThPE_pars->SetBranchAddress("FEB_ID",       &FEB_IDTh); 
    ThPE_pars->SetBranchAddress("CH",           &CHTh); 
    ThPE_pars->SetBranchAddress("mu",           &mu); 
    ThPE_pars->SetBranchAddress("sigma",        &sigma);  
    ThPE_pars->SetBranchAddress("ChiSq",        &ChiSqTh); 

    TFile *ThPE_Tcard       = new TFile(ThPE_card);  //input Tfile with ThPE data by cards 
    TTree *ThPE_card_pars   = (TTree*)ThPE_Tcard->Get(ThPE_card_tree); 
    ThPE_card_pars->SetBranchAddress("FEB_idc",      &FEB_idcTh); 
    ThPE_card_pars->SetBranchAddress("mu_mean",        &mu_mean); 
    ThPE_card_pars->SetBranchAddress("sig_mean",       &sig_mean);
    ThPE_card_pars->SetBranchAddress("ChiSq_mean",   &ChiSq_meanTh); 

    ThPE_card_pars->Draw("mu_mean>>mu_his", "", "goff");
    TH1F *mu_his = (TH1F*)gDirectory->Get("mu_his");
    Double_t mu_g   = mu_his->GetMean();    //mean of means
    ThPE_card_pars->Draw("sig_mean>>sig_his", "", "goff");
    TH1F *sig_his = (TH1F*)gDirectory->Get("sig_his");
    Double_t sig_g  = sig_his->GetMean();  //mean of means

    TFile *outfile      = new TFile("Rank_cards_v3.root", "RECREATE"); //output TFile with rankings //v1: when score_ch is divided by 64 after sqrt 
    TTree *rank         = new TTree("Ranks","TestBenchDB TT card rankings"); //output TTree         //v2: when score_ch is not divided by 64
    rank->Branch("FEB_id",      &FEB_id,        "FEB_id/I");                                        //v0: when score_ch was divided by 64 each time it was incremented
    rank->Branch("a1_r",        &a1_r,          "a1_r/D"); 
    rank->Branch("a2_r",        &a2_r,          "a2_r/D"); 
    rank->Branch("mu_r",        &mu_r,          "mu_r/D"); 
    rank->Branch("sig_r",       &sig_r,         "sig_r/D"); 
    rank->Branch("score_ch",    &score_ch,      "score_ch/D");
    rank->Branch("score_card",  &score_card,    "score_card/D"); 
    //rank->Branch("rank_no",  &rank_no,  "rank_no/I");

    for(int i=0; i<Qinj_card_pars->GetEntries(); ++i) //1160  //Qinj_card_pars->GetEntries();
    {
        Qinj_Tcard->cd();
        Qinj_card_pars->GetEntry(i);
        ThPE_Tcard->cd();
        ThPE_card_pars->GetEntry(i);
        score_ch = 0;
        score_card = 0;
        for(int j=0; j<64; ++j)
        {
            Qinj_Tfile->cd();
            Qinj_pars->GetEntry(i*64 + j);
            ThPE_Tfile->cd();
            ThPE_pars->GetEntry(i*64 + j);
            score_ch += ( pow((a1-a1_mean)/a1_mean , 2) + pow((a2-a2_mean)/a2_mean , 2) + pow((mu-mu_mean)/mu_mean , 2) + pow((sigma-sig_mean)/sig_mean , 2) )/64;
            
            if(status == 0) {score_ch += 0.1;   score_card += 0.1;} //demerit points for bad channels which don't fit
            if( (sigma <=1) || (mu==360) || (mu==520) ) {score_ch += 0.1;  score_card += 0.1;}
        }
        //if( (FEB_ID == 566) || (FEB_ID == 964) || (FEB_ID==867) || (FEB_ID=94) ||(FEB_ID==1092) ) {score_ch += 0.1;   score_card += 0.1;} 
        if(FEB_ID == 653) {score_ch += 1.;   score_card += 1.;}//too large a2
        if(FEB_ID == 895) {score_ch += 1.;   score_card += 1.;} //too large sigma
        if(FEB_ID == 566) {score_ch += 1.;   score_card += 1.;}////awkward "dips" on almost all channels in the NoQinj/ThirdPE test s-curves
        if(FEB_ID == 964) {score_ch += 1.;   score_card += 1.;} ////awkward "dips" on almost all channels in the NoQinj/ThirdPE test s-curves
        if(FEB_ID == 867) {score_ch += 1.;   score_card += 1.;}////awkward "dips" on almost all channels in the NoQinj/ThirdPE test s-curves
        if(FEB_ID == 94) {score_ch += 1.;   score_card += 1.;} ////awkward "dips" on almost all channels in the NoQinj/ThirdPE test s-curves
        if(FEB_ID == 1092) {score_ch += 1.;   score_card += 1.;}////awkward "dips" on almost all channels in the NoQinj/ThirdPE test s-curves

        score_card  += pow( pow((a1_mean-a1_g)/a1_g , 2) + pow((a2_mean-a2_g)/a2_g , 2) + pow((mu_mean-mu_g)/mu_g , 2) + pow((sig_mean-sig_g)/sig_g , 2) , 0.5);
        score_ch += pow(score_ch, 0.5);    FEB_id = FEB_ID;    a1_r = a1_mean; a2_r = a2_mean; mu_r = mu_mean;    sig_r = sig_mean;  //score_card = scores_card[FEB_ID];
        outfile->cd();
        rank->Fill();
    }
    outfile->cd();
    
    Int_t nEntries = rank->GetEntries();
    rank->Draw("score_ch:score_card", "", "goff");
    Int_t *index = new Int_t[nEntries];
    TMath::Sort(nEntries,rank->GetV1(),index);

    //TTree *rank_sort_ch = (TTree*)rank->CloneTree(0);
    TTree *rank_sort_ch         = new TTree("Ranks_sorted","TestBenchDB TT card rankings sorted by ch & card distributions");
    rank_sort_ch ->Branch("FEB_id",      &FEB_id,        "FEB_id/I");  
    rank_sort_ch ->Branch("score_ch",    &score_ch,      "score_ch/D");
    rank_sort_ch ->Branch("rank_ch",     &rank_ch,       "rank_ch/I"); 

	for (Int_t i=0; i<nEntries; i++) 
    {
		rank_ch = nEntries - i;
        rank->GetEntry(index[i]);
		rank_sort_ch->Fill();
	}
    //rank_sort_ch->Write();
	delete [] index;

    Int_t *index1 = new Int_t[nEntries];
    TMath::Sort(nEntries,rank->GetV2(),index1);

    //TTree *rank_sort_card = (TTree*)rank->CloneTree(0);
    TTree *rank_sort_card         = new TTree("Ranks_sort_card","TestBenchDB TT card rankings sorted by cards distribution");
    rank_sort_card ->Branch("FEB_id",      &FEB_id,        "FEB_id/I"); 
    rank_sort_card ->Branch("score_card",  &score_card,    "score_card/D");
    rank_sort_card ->Branch("rank_card",   &rank_card,     "rank_card/I");

	for (Int_t i=0; i<nEntries; i++) 
    {
		rank_card = nEntries - i;
        rank->GetEntry(index1[i]);
		rank_sort_card->Fill();
	}
    //rank_sort_card->Write();
	delete [] index1;

    TBranch *newBranch1 = rank_sort_ch-> Branch("score_card",  &score_card,    "score_card/D");
    TBranch *newBranch2 = rank_sort_ch-> Branch("rank_card",   &rank_card,     "rank_card/I");
    for(Int_t i=0; i<nEntries; ++i)
    {
        rank_sort_ch->GetEntry(i);
        temp_id = FEB_id;
        for(Int_t j=0; j<nEntries; ++j)
        {
            rank_sort_card->GetEntry(j);
            if(temp_id==FEB_id)
            {
                newBranch1->Fill();
                newBranch2->Fill();
                break;
            }
        }
    }
    rank_sort_ch->Write("",TObject::kOverwrite);
    rank->Write();
    delete outfile;
    cout<<"TFile closed."<<endl;
    return 0;
}