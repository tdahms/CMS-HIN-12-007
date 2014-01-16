#include <iostream>
#include <string>

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1.h"


void make_Jpsi_psi2s_histos(std::string fname="~/Jpsi_Histos_181611-183013_150invmub_mass2-5.root",bool isHI=true)
{
  gStyle->SetOptStat("emri");
  gStyle->SetOptFit(1);
  
  TFile *inf;
  inf = new TFile(fname.c_str(),"READ");
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->cd();
  TH1F *hmass_os_0100 = new TH1F("hmass_os_0100","h1mass_os_0100;m_{inv} (#mu^{+}#mu^{-}) [GeV/c^{2}];counts",80,2.6,4.2);
  TH1F *hmass_ls_0100 = new TH1F("hmass_ls_0100","h1mass_ls_0100;m_{inv} (#mu^{#pm}#mu^{#pm}) [GeV/c^{2}];counts",80,2.6,4.2);
  TH1F *hmass_su_0100 = new TH1F("hmass_su_0100","h1mass_su_0100;m_{inv} (#mu^{+}#mu^{-}) [GeV/c^{2}];counts",80,2.6,4.2);

  hmass_os_0100->Sumw2();
  hmass_os_0100->SetMarkerSize(1.0);
  hmass_os_0100->SetMarkerStyle(24);
  hmass_os_0100->SetMarkerColor(kBlack);
  hmass_os_0100->SetLineColor(kBlack);

  hmass_ls_0100->Sumw2();
  hmass_ls_0100->SetMarkerSize(1.0);
  hmass_ls_0100->SetMarkerStyle(24);
  hmass_ls_0100->SetMarkerColor(kRed);
  hmass_ls_0100->SetLineColor(kRed);

  hmass_su_0100->Sumw2();
  hmass_su_0100->SetMarkerSize(1.0);
  hmass_su_0100->SetMarkerStyle(20);
  hmass_su_0100->SetMarkerColor(kBlue);
  hmass_su_0100->SetLineColor(kBlue);


  TH1F *hmass_os_020 = (TH1F*) hmass_os_0100->Clone("hmass_os_020");
  TH1F *hmass_ls_020 = (TH1F*) hmass_ls_0100->Clone("hmass_ls_020");
  TH1F *hmass_su_020 = (TH1F*) hmass_su_0100->Clone("hmass_su_020");

  TH1F *hmass_os_2040 = (TH1F*) hmass_os_0100->Clone("hmass_os_2040");
  TH1F *hmass_ls_2040 = (TH1F*) hmass_ls_0100->Clone("hmass_ls_2040");
  TH1F *hmass_su_2040 = (TH1F*) hmass_su_0100->Clone("hmass_su_2040");

  TH1F *hmass_os_40100 = (TH1F*) hmass_os_0100->Clone("hmass_os_40100");
  TH1F *hmass_ls_40100 = (TH1F*) hmass_ls_0100->Clone("hmass_ls_40100");
  TH1F *hmass_su_40100 = (TH1F*) hmass_su_0100->Clone("hmass_su_40100");

  if (isHI) {
    // 0-100%
    ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>hmass_os_0100","Reco_QQ_sign==0&&Reco_QQ_type==0&&Reco_QQ_4mom.Pt()>3.0&&Reco_QQ_4mom.Pt()<6.5&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Centrality<40","e");
    ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>hmass_ls_0100","Reco_QQ_sign!=0&&Reco_QQ_type==0&&Reco_QQ_4mom.Pt()>3.0&&Reco_QQ_4mom.Pt()<6.5&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Centrality<40","e");
    hmass_su_0100->Add(hmass_os_0100,hmass_ls_0100,1,-1);

    // 0-20%
    ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>hmass_os_020","Reco_QQ_sign==0&&Reco_QQ_type==0&&Reco_QQ_4mom.Pt()>3.0&&Reco_QQ_4mom.Pt()<6.5&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Centrality<8","e");
    ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>hmass_ls_020","Reco_QQ_sign!=0&&Reco_QQ_type==0&&Reco_QQ_4mom.Pt()>3.0&&Reco_QQ_4mom.Pt()<6.5&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Centrality<8","e");
    hmass_su_020->Add(hmass_os_020,hmass_ls_020,1,-1);

    // 20-40%
    ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>hmass_os_2040","Reco_QQ_sign==0&&Reco_QQ_type==0&&Reco_QQ_4mom.Pt()>3.0&&Reco_QQ_4mom.Pt()<6.5&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Centrality>=8&&Centrality<16","e");
    ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>hmass_ls_2040","Reco_QQ_sign!=0&&Reco_QQ_type==0&&Reco_QQ_4mom.Pt()>3.0&&Reco_QQ_4mom.Pt()<6.5&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Centrality>=8&&Centrality<16","e");
    hmass_su_2040->Add(hmass_os_2040,hmass_ls_2040,1,-1);

    // 40-100%
    ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>hmass_os_40100","Reco_QQ_sign==0&&Reco_QQ_type==0&&Reco_QQ_4mom.Pt()>3.0&&Reco_QQ_4mom.Pt()<6.5&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Centrality>=16&&Centrality<40","e");
    ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>hmass_ls_40100","Reco_QQ_sign!=0&&Reco_QQ_type==0&&Reco_QQ_4mom.Pt()>3.0&&Reco_QQ_4mom.Pt()<6.5&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Centrality>=16&&Centrality<40","e");
    hmass_su_40100->Add(hmass_os_40100,hmass_ls_40100,1,-1);
  }
  else {
    // pp
    ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>hmass_os_0100","Reco_QQ_sign==0&&Reco_QQ_type==0&&Reco_QQ_4mom.Pt()>3.0&&Reco_QQ_4mom.Pt()<6.5&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4","e");
    ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>hmass_ls_0100","Reco_QQ_sign!=0&&Reco_QQ_type==0&&Reco_QQ_4mom.Pt()>3.0&&Reco_QQ_4mom.Pt()<6.5&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4","e");
    hmass_su_0100->Add(hmass_os_0100,hmass_ls_0100,1,-1);
  }

  hmass_os_0100->Draw();
  hmass_ls_0100->Draw("histsame");
  hmass_su_0100->Draw("same");


  TFile *outf;
  if (isHI)
    outf = new TFile("histos_pbpb_pt3065_y1624.root","CREATE");
  else
    outf = new TFile("histos_pp_pt3065_y1624.root","CREATE");
  hmass_os_0100->Write();
  hmass_ls_0100->Write();
  hmass_su_0100->Write();

  hmass_os_020->Write();
  hmass_ls_020->Write();
  hmass_su_020->Write();

  hmass_os_2040->Write();
  hmass_ls_2040->Write();
  hmass_su_2040->Write();

  hmass_os_40100->Write();
  hmass_ls_40100->Write();
  hmass_su_40100->Write();

  outf->Close();

  return;
}
