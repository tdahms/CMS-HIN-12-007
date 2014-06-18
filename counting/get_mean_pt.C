#include <iostream>
#include <string>

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1.h"
#include "TMath.h"


void get_mean_pt(bool isPbPb=true,bool savePlots=false)
{
  std::string infname;
  if (isPbPb)
    infname="/data/CMS/data/PbPb/All_Histos_cmssw445p5_RegIt_EvtPlane_hStats.root";
  else
    infname="/data/CMS/data/pp/All_v2.24_Histos_Runs_211739-211831_GlbGlb_woPileUpRej_muLessPV.root";

  TFile *inf;
  inf = new TFile(infname.c_str(),"READ");
  TCanvas *c0 = new TCanvas("c0","c0");
  TCanvas *c1 = new TCanvas("c1","c1");
  TCanvas *c2 = new TCanvas("c2","c2");
  c2->Divide(2,2);

  // TCanvas *c2 = new TCanvas("c2","c2");
  // c2->Divide(2,2);

  // TCanvas *c3 = new TCanvas("c3","c3");
  // c3->Divide(2,1);
  c2->cd(1);

  TH1F *hmassMid = new TH1F("hmassMid","hmassMid;m_{#mu#mu} (GeV/c^{2});counts",160,1.8,5);
  TH1F *hmassMida = new TH1F("hmassMida","hmassMida;m_{#mu#mu} (GeV/c^{2});counts",160,1.8,5);
  TH1F *hmassFwd = new TH1F("hmassFwd","hmassFwd;m_{#mu#mu} (GeV/c^{2});counts",160,1.8,5);
  TH1F *hmassFwda = new TH1F("hmassFwda","hmassFwda;m_{#mu#mu} (GeV/c^{2});counts",160,1.8,5);
  TH1F *hmassFwdP = new TH1F("hmassFwdP","hmassFwdP;p_{#mu^{+}#mu^{+}} (GeV/c^{2});counts",160,1.8,5);
  TH1F *hmassFwdM = new TH1F("hmassFwdM","hmassFwdM;p_{#mu^{-}#mu^{-}} (GeV/c^{2});counts",160,1.8,5);
  hmassMid->Sumw2();
  hmassMida->Sumw2();
  hmassFwd->Sumw2();
  hmassFwda->Sumw2();
  hmassFwdP->Sumw2();
  hmassFwdM->Sumw2();

  hmassMida->SetMarkerStyle(24);
  hmassFwdP->SetMarkerStyle(24);
  hmassFwdM->SetMarkerStyle(25);

  hmassMida->SetMarkerColor(2);
  hmassFwdP->SetMarkerColor(2);
  hmassFwdM->SetMarkerColor(4);
  hmassFwda->SetLineColor(kGreen+2);
  hmassFwda->SetLineWidth(2);


  TH1F *h0 = new TH1F("h0","h0;p_{T} (#mu^{+}#mu^{-}) [GeV/c];counts",300,0.0,30.0);
  TH1F *h1 = new TH1F("h1","h1;p_{T} (#mu^{+}#mu^{-}) [GeV/c];counts",300,0.0,30.0);
  TH1F *h2 = new TH1F("h2","h2;p_{T} (#mu^{+}#mu^{-}) [GeV/c];counts",300,0.0,30.0);
  TH1F *h3 = new TH1F("h3","h3;p_{T} (#mu^{+}#mu^{-}) [GeV/c];counts",300,0.0,30.0);

  TH1F *h4 = new TH1F("h4","h4;p_{T} (#mu^{+}#mu^{-}) [GeV/c];counts",300,0.0,30.0);
  TH1F *h5 = new TH1F("h5","h5;p_{T} (#mu^{+}#mu^{-}) [GeV/c];counts",300,0.0,30.0);
  TH1F *h6 = new TH1F("h6","h6;p_{T} (#mu^{+}#mu^{-}) [GeV/c];counts",300,0.0,30.0);
  TH1F *h7 = new TH1F("h7","h7;p_{T} (#mu^{+}#mu^{-}) [GeV/c];counts",300,0.0,30.0);

  TH1F *h8 = new TH1F("h8","h8;p_{T} (#mu^{+}#mu^{-}) [GeV/c];counts",300,0.0,30.0);
  TH1F *h9 = new TH1F("h9","h9;p_{T} (#mu^{+}#mu^{-}) [GeV/c];counts",300,0.0,30.0);


  TH1F *h0a = new TH1F("h0a","h0a;p_{T} (#mu^{#pm}#mu^{#pm}) [GeV/c];counts",300,0.0,30.0);
  TH1F *h1a = new TH1F("h1a","h1a;p_{T} (#mu^{#pm}#mu^{#pm}) [GeV/c];counts",300,0.0,30.0);
  TH1F *h2a = new TH1F("h2a","h2a;p_{T} (#mu^{#pm}#mu^{#pm}) [GeV/c];counts",300,0.0,30.0);
  TH1F *h3a = new TH1F("h3a","h3a;p_{T} (#mu^{#pm}#mu^{#pm}) [GeV/c];counts",300,0.0,30.0);

  TH1F *h4a = new TH1F("h4a","h4a;p_{T} (#mu^{+}#mu^{-}) [GeV/c];counts",300,0.0,30.0);
  TH1F *h5a = new TH1F("h5a","h5a;p_{T} (#mu^{+}#mu^{-}) [GeV/c];counts",300,0.0,30.0);
  TH1F *h6a = new TH1F("h6a","h6a;p_{T} (#mu^{+}#mu^{-}) [GeV/c];counts",300,0.0,30.0);
  TH1F *h7a = new TH1F("h7a","h7a;p_{T} (#mu^{+}#mu^{-}) [GeV/c];counts",300,0.0,30.0);

  TH1F *h8a = new TH1F("h8a","h8a;p_{T} (#mu^{+}#mu^{-}) [GeV/c];counts",300,0.0,30.0);
  TH1F *h9a = new TH1F("h9a","h9a;p_{T} (#mu^{+}#mu^{-}) [GeV/c];counts",300,0.0,30.0);

  h0->Sumw2();
  h1->Sumw2();
  h2->Sumw2();
  h3->Sumw2();

  h4->Sumw2();
  h5->Sumw2();
  h6->Sumw2();
  h7->Sumw2();

  h8->Sumw2();
  h9->Sumw2();

  h0a->Sumw2();
  h1a->Sumw2();
  h2a->Sumw2();
  h3a->Sumw2();

  h4a->Sumw2();
  h5a->Sumw2();
  h6a->Sumw2();
  h7a->Sumw2();

  h8a->Sumw2();
  h9a->Sumw2();

  h0a->SetMarkerStyle(24);
  h1a->SetMarkerStyle(24);

  h0a->SetMarkerColor(2);
  h1a->SetMarkerColor(2);

  if (isPbPb) {
    c0->cd(1);
    ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>hmassMid","(HLTriggers&1)==1&&(Reco_QQ_trig&1)==1&&Reco_QQ_sign==0&&Reco_QQ_type==0&&abs(Reco_QQ_4mom.Rapidity())<1.6&&Reco_QQ_4mom.Pt()>6.5&&Reco_QQ_4mom.Pt()<30","e");
    ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>hmassMida","(HLTriggers&1)==1&&(Reco_QQ_trig&1)==1&&Reco_QQ_sign!=0&&Reco_QQ_type==0&&abs(Reco_QQ_4mom.Rapidity())<1.6&&Reco_QQ_4mom.Pt()>6.5&&Reco_QQ_4mom.Pt()<30","esame");
    c1->cd(1);
    ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>hmassFwd","(HLTriggers&1)==1&&(Reco_QQ_trig&1)==1&&Reco_QQ_sign==0&&Reco_QQ_type==0&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Reco_QQ_4mom.Pt()>3&&Reco_QQ_4mom.Pt()<30&&Centrality<8","e");
    ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>hmassFwdP","(HLTriggers&1)==1&&(Reco_QQ_trig&1)==1&&Reco_QQ_sign==1&&Reco_QQ_type==0&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Reco_QQ_4mom.Pt()>3&&Reco_QQ_4mom.Pt()<30&&Centrality<8","esame");
    ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>hmassFwdM","(HLTriggers&1)==1&&(Reco_QQ_trig&1)==1&&Reco_QQ_sign==2&&Reco_QQ_type==0&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Reco_QQ_4mom.Pt()>3&&Reco_QQ_4mom.Pt()<30&&Centrality<8","esame");
    c2->cd(1);
    ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.Pt()>>h0","(HLTriggers&1)==1&&(Reco_QQ_trig&1)==1&&Reco_QQ_sign==0&&Reco_QQ_type==0&&Reco_QQ_4mom.M()>=3.6&&Reco_QQ_4mom.M()<3.76&&abs(Reco_QQ_4mom.Rapidity())<1.6","e");
    ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.Pt()>>h0a","(HLTriggers&1)==1&&(Reco_QQ_trig&1)==1&&Reco_QQ_sign!=0&&Reco_QQ_type==0&&Reco_QQ_4mom.M()>=3.6&&Reco_QQ_4mom.M()<3.76&&abs(Reco_QQ_4mom.Rapidity())<1.6","esame");
    c2->cd(2);
    ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.Pt()>>h1","(HLTriggers&1)==1&&(Reco_QQ_trig&1)==1&&Reco_QQ_sign==0&&Reco_QQ_type==0&&Reco_QQ_4mom.M()>=3.6&&Reco_QQ_4mom.M()<3.76&&abs(Reco_QQ_4mom.Rapidity())>=1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4","e");
    ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.Pt()>>h1a","(HLTriggers&1)==1&&(Reco_QQ_trig&1)==1&&Reco_QQ_sign!=0&&Reco_QQ_type==0&&Reco_QQ_4mom.M()>=3.6&&Reco_QQ_4mom.M()<3.76&&abs(Reco_QQ_4mom.Rapidity())>=1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4","esame");
    // c2->cd(1);
    // ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.Pt()>>h0","(HLTriggers&1)==1&&(Reco_QQ_trig&1)==1&&Reco_QQ_sign==0&&Reco_QQ_type==0&&Reco_QQ_4mom.M()>=3.0&&Reco_QQ_4mom.M()<3.2&&abs(Reco_QQ_4mom.Rapidity())<1.6&&Centrality>=16&&Centrality<40","e");
    // ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.Pt()>>h0a","(HLTriggers&1)==1&&(Reco_QQ_trig&1)==1&&Reco_QQ_sign!=0&&Reco_QQ_type==0&&Reco_QQ_4mom.M()>=3.0&&Reco_QQ_4mom.M()<3.2&&abs(Reco_QQ_4mom.Rapidity())<1.6&&Centrality>=16&&Centrality<40","esame");
    // c2->cd(2);
    // ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.Pt()>>h1","(HLTriggers&1)==1&&(Reco_QQ_trig&1)==1&&Reco_QQ_sign==0&&Reco_QQ_type==0&&Reco_QQ_4mom.M()>=3.0&&Reco_QQ_4mom.M()<3.2&&abs(Reco_QQ_4mom.Rapidity())>=1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Centrality>=16&&Centrality<40","e");
    // ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.Pt()>>h1a","(HLTriggers&1)==1&&(Reco_QQ_trig&1)==1&&Reco_QQ_sign!=0&&Reco_QQ_type==0&&Reco_QQ_4mom.M()>=3.0&&Reco_QQ_4mom.M()<3.2&&abs(Reco_QQ_4mom.Rapidity())>=1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Centrality>=16&&Centrality<40","esame");
  }
  else {
    c0->cd(1);
    ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>hmassMid","(HLTriggers&2)==2&&(Reco_QQ_trig&2)==2&&Reco_QQ_sign==0&&Reco_QQ_type==0&&abs(Reco_QQ_4mom.Rapidity())<1.6","e");
    ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>hmassMida","(HLTriggers&2)==2&&(Reco_QQ_trig&2)==2&&Reco_QQ_sign!=0&&Reco_QQ_type==0&&abs(Reco_QQ_4mom.Rapidity())<1.6","esame");
    c1->cd(1);
    ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>hmassFwd","(HLTriggers&2)==2&&(Reco_QQ_trig&2)==2&&Reco_QQ_sign==0&&Reco_QQ_type==0&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4","e");
    ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>hmassFwda","(HLTriggers&2)==2&&(Reco_QQ_trig&2)==2&&Reco_QQ_sign!=0&&Reco_QQ_type==0&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4","esame");

    c2->cd(1);
    ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.Pt()>>h0","(HLTriggers&2)==2&&(Reco_QQ_trig&2)==2&&Reco_QQ_sign==0&&Reco_QQ_type==0&&Reco_QQ_4mom.M()>=3.6&&Reco_QQ_4mom.M()<3.76&&abs(Reco_QQ_4mom.Rapidity())<1.6","e");
    ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.Pt()>>h0a","(HLTriggers&2)==2&&(Reco_QQ_trig&2)==2&&Reco_QQ_sign!=0&&Reco_QQ_type==0&&Reco_QQ_4mom.M()>=3.6&&Reco_QQ_4mom.M()<3.76&&abs(Reco_QQ_4mom.Rapidity())<1.6","esame");
    c2->cd(2);
    ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.Pt()>>h1","(HLTriggers&2)==2&&(Reco_QQ_trig&2)==2&&Reco_QQ_sign==0&&Reco_QQ_type==0&&Reco_QQ_4mom.M()>=3.6&&Reco_QQ_4mom.M()<3.76&&abs(Reco_QQ_4mom.Rapidity())>=1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4","e");
    ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.Pt()>>h1a","(HLTriggers&2)==2&&(Reco_QQ_trig&2)==2&&Reco_QQ_sign!=0&&Reco_QQ_type==0&&Reco_QQ_4mom.M()>=3.6&&Reco_QQ_4mom.M()<3.76&&abs(Reco_QQ_4mom.Rapidity())>=1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4","esame");
    // c2->cd(1);
    // ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.Pt()>>h0","(HLTriggers&2)==2&&(Reco_QQ_trig&2)==2&&Reco_QQ_sign==0&&Reco_QQ_type==0&&Reco_QQ_4mom.M()>=3.0&&Reco_QQ_4mom.M()<3.2&&abs(Reco_QQ_4mom.Rapidity())<1.6","e");
    // ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.Pt()>>h0a","(HLTriggers&2)==2&&(Reco_QQ_trig&2)==2&&Reco_QQ_sign!=0&&Reco_QQ_type==0&&Reco_QQ_4mom.M()>=3.0&&Reco_QQ_4mom.M()<3.2&&abs(Reco_QQ_4mom.Rapidity())<1.6","esame");
    // c2->cd(2);
    // ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.Pt()>>h1","(HLTriggers&2)==2&&(Reco_QQ_trig&2)==2&&Reco_QQ_sign==0&&Reco_QQ_type==0&&Reco_QQ_4mom.M()>=3.0&&Reco_QQ_4mom.M()<3.2&&abs(Reco_QQ_4mom.Rapidity())>=1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4","e");
    // ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.Pt()>>h1a","(HLTriggers&2)==2&&(Reco_QQ_trig&2)==2&&Reco_QQ_sign!=0&&Reco_QQ_type==0&&Reco_QQ_4mom.M()>=3.0&&Reco_QQ_4mom.M()<3.2&&abs(Reco_QQ_4mom.Rapidity())>=1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4","esame");
  }

  for (int i=1;i<hmassFwdP->GetNbinsX();++i) {
    hmassFwda->SetBinContent(i,2*sqrt(hmassFwdP->GetBinContent(i)*hmassFwdM->GetBinContent(i)));
    hmassFwda->SetBinError(i,sqrt(hmassFwda->GetBinContent(i)));
  }

  c0->cd();
  c0->SetLogy();
  hmassMid->GetXaxis()->CenterTitle(1);
  if (!isPbPb)
    hmassMid->GetYaxis()->SetRangeUser(1,1e4);
  hmassMid->Draw();
  hmassMida->Draw("same");
  c1->cd();
  c1->SetLogy();
  hmassFwd->GetXaxis()->CenterTitle(1);
  if (!isPbPb)
    hmassFwd->GetYaxis()->SetRangeUser(10,6e3);
  hmassFwd->Draw();
  hmassFwda->Draw("same");


  if (savePlots) {
    if (isPbPb) {
      c0->SaveAs("invmass_os_ss_comp_mid_pp.pdf");
      c1->SaveAs("invmass_os_ss_comp_fwd_pp.pdf");
    }
    else {
      c0->SaveAs("invmass_os_ss_comp_mid_pp.pdf");
      c1->SaveAs("invmass_os_ss_comp_fwd_pp.pdf");
    }
  }
  h0->Add(h0a,-1.0);
  h1->Add(h1a,-1.0);
  // h2->Add(h2a,-1.0);
  // h3->Add(h3a,-1.0);
  // h4->Add(h4a,-1.0);
  // h5->Add(h5a,-1.0);
  // h6->Add(h6a,-1.0);
  // h7->Add(h7a,-1.0);
  // h8->Add(h8a,-1.0);
  // h9->Add(h9a,-1.0);

  c2->cd(3);
  h0->Draw();
  c2->cd(4);
  h1->Draw();
  // c2->cd(3);
  // h2->Draw();
  // c2->cd(4);
  // h3->Draw();

  // c2->cd(1);
  // h4->Draw();
  // c2->cd(2);
  // h5->Draw();
  // c2->cd(3);
  // h6->Draw();
  // c2->cd(4);
  // h7->Draw();

  // c3->cd(1);
  // h8->Draw();
  // c3->cd(2);
  // h9->Draw();

  std::cout << "<pT> \t (RMS)" << std::endl;
  // 0-100%
  // All eta
  // pt 0-30
  // pt 0-6.5
  // pt 6.5-10
  // pt 10-30

  // Barrel
  // pt 6.5-30

  // Overlap
  // pt 2.0-30

  // EndCap
  // pt 0.0-30

  if (isPbPb)
    std::cout << "PbPb" << std::endl;
  else
    std::cout << "PbPb" << std::endl;
  double epsilon=0.0001;
  int ptMin = h0->GetXaxis()->FindBin(0.0+epsilon);
  int ptMax = h0->GetXaxis()->FindBin(30.0-epsilon);
  // h0->GetXaxis()->SetRange(ptMin,ptMax);
  // std::cout << h0->GetMean() << " \t " << h0->GetRMS() << std::endl;

  // ptMin = h0->GetXaxis()->FindBin(0.0+epsilon);
  // ptMax = h0->GetXaxis()->FindBin(6.5-epsilon);
  // h0->GetXaxis()->SetRange(ptMin,ptMax);
  // std::cout << h0->GetMean() << " \t " << h0->GetRMS() << std::endl;

  // ptMin = h0->GetXaxis()->FindBin(6.5+epsilon);
  // ptMax = h0->GetXaxis()->FindBin(10.0-epsilon);
  // h0->GetXaxis()->SetRange(ptMin,ptMax);
  // std::cout << h0->GetMean() << " \t " << h0->GetRMS() << std::endl;

  std::cout << "Barrel (|y|<1.6 & 6.5<pT<30 GeV/c)" << std::endl;
  ptMin = h0->GetXaxis()->FindBin(6.5+epsilon);
  ptMax = h0->GetXaxis()->FindBin(30.0-epsilon);
  h0->GetXaxis()->SetRange(ptMin,ptMax);
  std::cout << h0->GetMean() << " \t " << h0->GetRMS() << std::endl;


  std::cout << "EndCap (1.6<|y|<2.4 & 3<pT<30 GeV/c)" << std::endl;
  ptMin = h1->GetXaxis()->FindBin(3.0+epsilon);
  ptMax = h1->GetXaxis()->FindBin(30.0-epsilon);
  h1->GetXaxis()->SetRange(ptMin,ptMax);
  std::cout << h1->GetMean() << " \t " << h1->GetRMS() << std::endl;

  // 0-10%
  // 10-20%
  // 20-60%
  // 40-100%
  // All eta
  // pt 0-30
  // pt 6.5-30

  // std::cout << "0-10" << std::endl;
  // ptMin = h4->GetXaxis()->FindBin(0.0+epsilon);
  // ptMax = h4->GetXaxis()->FindBin(30.0-epsilon);
  // h4->GetXaxis()->SetRange(ptMin,ptMax);
  // std::cout << h4->GetMean() << " \t " << h4->GetRMS() << std::endl;

  // ptMin = h4->GetXaxis()->FindBin(6.5+epsilon);
  // ptMax = h4->GetXaxis()->FindBin(30.0-epsilon);
  // h4->GetXaxis()->SetRange(ptMin,ptMax);
  // std::cout << h4->GetMean() << " \t " << h4->GetRMS() << std::endl;

  // std::cout << "10-20" << std::endl;
  // ptMin = h5->GetXaxis()->FindBin(0.0+epsilon);
  // ptMax = h5->GetXaxis()->FindBin(30.0-epsilon);
  // h5->GetXaxis()->SetRange(ptMin,ptMax);
  // std::cout << h5->GetMean() << " \t " << h5->GetRMS() << std::endl;

  // ptMin = h5->GetXaxis()->FindBin(6.5+epsilon);
  // ptMax = h5->GetXaxis()->FindBin(30.0-epsilon);
  // h5->GetXaxis()->SetRange(ptMin,ptMax);
  // std::cout << h5->GetMean() << " \t " << h5->GetRMS() << std::endl;

  // std::cout << "20-40" << std::endl;
  // ptMin = h6->GetXaxis()->FindBin(0.0+epsilon);
  // ptMax = h6->GetXaxis()->FindBin(30.0-epsilon);
  // h6->GetXaxis()->SetRange(ptMin,ptMax);
  // std::cout << h6->GetMean() << " \t " << h6->GetRMS() << std::endl;

  // ptMin = h6->GetXaxis()->FindBin(6.5+epsilon);
  // ptMax = h6->GetXaxis()->FindBin(30.0-epsilon);
  // h6->GetXaxis()->SetRange(ptMin,ptMax);
  // std::cout << h6->GetMean() << " \t " << h6->GetRMS() << std::endl;

  // std::cout << "40-100" << std::endl;
  // ptMin = h7->GetXaxis()->FindBin(0.0+epsilon);
  // ptMax = h7->GetXaxis()->FindBin(30.0-epsilon);
  // h7->GetXaxis()->SetRange(ptMin,ptMax);
  // std::cout << h7->GetMean() << " \t " << h7->GetRMS() << std::endl;

  // ptMin = h7->GetXaxis()->FindBin(6.5+epsilon);
  // ptMax = h7->GetXaxis()->FindBin(30.0-epsilon);
  // h7->GetXaxis()->SetRange(ptMin,ptMax);
  // std::cout << h7->GetMean() << " \t " << h7->GetRMS() << std::endl;


  // 0-100%
  // All eta
  // pt 0-30
  // pt 6.5-30

  // EndCap
  // pt 0-6.5
  // pt 6.5-30

  // 0-20%
  // 20-100%
  // All eta
  // pt 0-30
  // pt 6.5-30

  // std::cout << "Non-prompt binning" << std::endl;
  // std::cout << "All" << std::endl;
  // ptMin = h0->GetXaxis()->FindBin(0.0+epsilon);
  // ptMax = h0->GetXaxis()->FindBin(30.0-epsilon);
  // h0->GetXaxis()->SetRange(ptMin,ptMax);
  // std::cout << h0->GetMean() << " \t " << h0->GetRMS() << std::endl;

  // ptMin = h0->GetXaxis()->FindBin(6.5+epsilon);
  // ptMax = h0->GetXaxis()->FindBin(30.0-epsilon);
  // h0->GetXaxis()->SetRange(ptMin,ptMax);
  // std::cout << h0->GetMean() << " \t " << h0->GetRMS() << std::endl;


  // std::cout << "EndCap" << std::endl;
  // ptMin = h3->GetXaxis()->FindBin(0.0+epsilon);
  // ptMax = h3->GetXaxis()->FindBin(6.5-epsilon);
  // h3->GetXaxis()->SetRange(ptMin,ptMax);
  // std::cout << h3->GetMean() << " \t " << h3->GetRMS() << std::endl;

  // ptMin = h3->GetXaxis()->FindBin(6.5+epsilon);
  // ptMax = h3->GetXaxis()->FindBin(30.0-epsilon);
  // h3->GetXaxis()->SetRange(ptMin,ptMax);
  // std::cout << h3->GetMean() << " \t " << h3->GetRMS() << std::endl;

  // std::cout << "All" << std::endl;
  // std::cout << "0-20" << std::endl;
  // ptMin = h8->GetXaxis()->FindBin(0.0+epsilon);
  // ptMax = h8->GetXaxis()->FindBin(30.0-epsilon);
  // h8->GetXaxis()->SetRange(ptMin,ptMax);
  // std::cout << h8->GetMean() << " \t " << h8->GetRMS() << std::endl;

  // ptMin = h8->GetXaxis()->FindBin(6.5+epsilon);
  // ptMax = h8->GetXaxis()->FindBin(30.0-epsilon);
  // h8->GetXaxis()->SetRange(ptMin,ptMax);
  // std::cout << h8->GetMean() << " \t " << h8->GetRMS() << std::endl;

  // std::cout << "20-100" << std::endl;
  // ptMin = h9->GetXaxis()->FindBin(0.0+epsilon);
  // ptMax = h9->GetXaxis()->FindBin(30.0-epsilon);
  // h9->GetXaxis()->SetRange(ptMin,ptMax);
  // std::cout << h9->GetMean() << " \t " << h9->GetRMS() << std::endl;

  // ptMin = h9->GetXaxis()->FindBin(6.5+epsilon);
  // ptMax = h9->GetXaxis()->FindBin(30.0-epsilon);
  // h9->GetXaxis()->SetRange(ptMin,ptMax);
  // std::cout << h9->GetMean() << " \t " << h9->GetRMS() << std::endl;


  return;
}
