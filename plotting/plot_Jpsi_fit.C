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


double crystal_ball_expo(double *x, double *par)
{
  double N = par[0];

  double mean = par[1];
  double sigma = par[2];
  double alpha = par[3];
  double n = par[4];

  double A = pow(n/fabs(alpha),n)*exp(-pow(alpha,2)/2);
  double B = n/fabs(alpha) - fabs(alpha);

  // factor 0.050 for bin width to get proper integrated yield
  if ((x[0]-mean)/sigma>-alpha)
    return 0.0250*N*TMath::Gaus(x[0],mean,sigma,1)+TMath::Exp(par[5]-x[0]*par[6]);
  else
    return 0.0250*N/(sqrt(TMath::TwoPi())*sigma)*A*pow(B-(x[0]-mean)/sigma,-n)+TMath::Exp(par[5]-x[0]*par[6]);
}

double two_crystal_ball_expo(double *x, double *par)
{
  double N1 = par[0];

  double mean = par[1];
  double sigma = par[2];
  double alpha = par[3];
  double n = par[4];


  double N2 = par[5];

  double A = pow(n/fabs(alpha),n)*exp(-pow(alpha,2)/2);
  double B = n/fabs(alpha) - fabs(alpha);

  // factor 0.050 for bin width to get proper integrated yield
  if ((x[0]-mean)/sigma>-alpha)
    return 0.0250*(N1*TMath::Gaus(x[0],mean,sigma,1) + N2*TMath::Gaus(x[0],mean*3.686/3.0969,sigma*3.686/3.0969,1))+TMath::Exp(par[6]-x[0]*par[7]);
  else
    return 0.0250*(N1/(sqrt(TMath::TwoPi())*sigma)*A*pow(B-(x[0]-mean)/sigma,-n) + N2/(sqrt(TMath::TwoPi())*sigma*3.686/3.0969)*A*pow(B-(x[0]-mean*3.686/3.0969)/sigma*3.686/3.0969,-n))+TMath::Exp(par[6]-x[0]*par[7]);
}


void plot_Jpsi_fit(std::string fname="Histos_OniaSkim_ReReco_v13_h3_sum.root", bool isHI=true, bool fitPsiP=false)
{
  gStyle->SetOptStat("emri");
  gStyle->SetOptFit(1);

  TF1 *f1 = new TF1("f1","0.0250*[0]*TMath::Gaus(x,[1],[2],1)",2.6,3.5);
  TF1 *f2 = new TF1("f2","0.0250*[0]*TMath::Gaus(x,[1],[2],1)+TMath::Exp([3]-x*[4])",2.6,3.5);
  f1->SetParNames("Yield","Mean","Sigma");
  f2->SetParNames("Yield","Mean","Sigma","a","b");
  f1->SetNpx(10000);
  f2->SetNpx(10000);
  TF1 *f4 = new TF1("f4","0.0250*[0]*TMath::Gaus(x,[1],[2],1)+TMath::Exp([3]-x*[4])",2.6,3.5);
  f4->SetParNames("Yield (#psi')","Mean","Sigma","a","b");
  f4->SetNpx(10000);

  // TF1 *f5 = new TF1("f5","0.0250*([0]*TMath::Gaus(x,[1],[2],1)+[3]*TMath::Gaus(x,[4],[5],1))+TMath::Exp([6]-x*[7])",2.6,3.5);
  // f5->SetParNames("Yield (J/#psi)","Mean (J/#psi)","Sigma (J/#psi)","Yield (#psi')","Mean (#psi')","Sigma (#psi')","a","b");
  // f5->SetNpx(10000);
  // f5->SetLineWidth(2);
  TF1 *f5 = new TF1("f5","0.0250*([0]*TMath::Gaus(x,[1],[2],1)+[3]*TMath::Gaus(x,[1]*3.686/3.0969,[2]*3.686/3.0969,1))+TMath::Exp([6]-x*[7])",2.6,3.5);
  f5->SetParNames("Yield (J/#psi)","Mean (J/#psi)","Sigma (J/#psi)","Yield (#psi')","a","b");
  f5->SetNpx(10000);
  f5->SetLineWidth(1);

  TF1 *f5a = new TF1("f5",two_crystal_ball_expo,2.6,3.5,8);
  f5a->SetParNames("Yield (J/#psi)","Mean (J/#psi)","Sigma (J/#psi)","#alpha","n","Yield (#psi')","a","b");
  f5a->SetNpx(10000);
  f5a->SetLineWidth(1);


  TF1 *f6 = new TF1("f6",crystal_ball_expo,2.6,3.5,8);
  f6->SetParNames("Yield (J/#psi)","Mean","Sigma","#alpha","n","a","b");
  f6->SetNpx(10000);
  f6->SetLineWidth(1);

  TF1 *f7 = new TF1("f6",crystal_ball_expo,3.2,4.0,8);
  f7->SetParNames("Yield (#psi')","Mean","Sigma","#alpha","n","a","b");
  f7->SetNpx(10000);
  f7->SetLineWidth(1);


  TF1 *f3 = new TF1("f3","TMath::Exp([0]-x*[1])",2.6,3.5);
  f3->SetParNames("a","b");
  
  TFile *inf;
  inf = new TFile(fname.c_str(),"READ");
  TCanvas *c1 = new TCanvas("c1","c1");
  if (!isHI)
    c1->SetLogy();
  c1->cd();
  TH1F *h1 = new TH1F("h1","h1;m_{inv} (#mu^{#pm}#mu^{#pm}) [GeV/c^{2}];counts",480,2.0,14.0);
  TH1F *h2 = new TH1F("h2","h2;m_{inv} (#mu^{+}#mu^{-}) [GeV/c^{2}];counts",480,2.0,14.0);
  TH1F *h3 = new TH1F("h3","h3;m_{inv} (#mu^{#pm}#mu^{#pm}) [GeV/c^{2}];counts",480,2.0,14.0);
  TH1F *h4 = new TH1F("h4","h4;m_{inv} (#mu^{+}#mu^{-}) [GeV/c^{2}];counts",480,2.0,14.0);
  TH1F *h5 = new TH1F("h5","h5;m_{inv} (#mu^{#pm}#mu^{#pm}) [GeV/c^{2}];counts",480,2.0,14.0);
  TH1F *h6 = new TH1F("h6","h6;m_{inv} (#mu^{+}#mu^{-}) [GeV/c^{2}];counts",480,2.0,14.0);
  TH1F *h7 = new TH1F("h7","h7;m_{inv} (#mu^{#pm}#mu^{#pm}) [GeV/c^{2}];counts",480,2.0,14.0);
  TH1F *h8 = new TH1F("h8","h8;m_{inv} (#mu^{+}#mu^{-}) [GeV/c^{2}];counts",480,2.0,14.0);

  TH1F *hPt = new TH1F("hPt","hPt;p_{T} (#mu^{+}#mu^{-}) [GeV/c];counts",300,0.0,30.0);
  TH1F *hPt2 = new TH1F("hPt2","hPt2;p_{T} (#mu^{#pm}#mu^{#pm}) [GeV/c];counts",300,0.0,30.0);

  TH1F *hRap = new TH1F("hRap","hRap;|y| (#mu^{+}#mu^{-});counts",24,0.0,2.4);
  TH1F *hRap2 = new TH1F("hRap2","hRap2;|y| (#mu^{#pm}#mu^{#pm});counts",24,0.0,2.4);

  h1->Sumw2();
  h1->SetMarkerSize(1.0);
  h1->SetMarkerStyle(24);
  h1->SetMarkerColor(kBlue);
  h1->SetLineColor(kBlue);
  h2->Sumw2();
  h2->SetMarkerSize(1.0);
  h2->SetMarkerStyle(20);
  h2->SetMarkerColor(kRed);
  h2->SetLineColor(kRed);
  h3->Sumw2();
  h3->SetMarkerStyle(24);
  h3->SetMarkerColor(kBlue);
  h3->SetLineColor(kBlue);
  h4->Sumw2();
  h4->SetMarkerStyle(20);
  h4->SetMarkerColor(kRed);
  h4->SetLineColor(kRed);
  h5->Sumw2();
  h5->SetMarkerStyle(24);
  h5->SetMarkerColor(kBlue);
  h5->SetLineColor(kBlue);
  h6->Sumw2();
  h6->SetMarkerStyle(20);
  h6->SetMarkerColor(kRed);
  h6->SetLineColor(kRed);
  h7->Sumw2();
  h7->SetMarkerStyle(24);
  h7->SetMarkerColor(kBlue);
  h7->SetLineColor(kBlue);
  h8->Sumw2();
  h8->SetMarkerStyle(20);
  h8->SetMarkerColor(kRed);
  h8->SetLineColor(kRed);
  // h1->Rebin();
  // h2->Rebin();
  // h3->Rebin();
  hPt->Sumw2();
  hPt2->Sumw2();
  hRap->Sumw2();
  hRap2->Sumw2();

  // default
  // 38305
  ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>h2","Reco_QQ_sign==0&&Reco_QQ_type==0&&Reco_QQ_VtxProb>0.01&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Centrality>=4","e");//&&abs(Reco_QQ_mupl_4mom.Eta())<2.1&&abs(Reco_QQ_mumi_4mom.Eta())<2.1&&Reco_QQ_4mom.Pt()>6.5
  ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>h1","Reco_QQ_sign!=0&&Reco_QQ_type==0&&Reco_QQ_VtxProb>0.01&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Centrality>=4","e");//&&abs(Reco_QQ_mupl_4mom.Eta())<2.1&&abs(Reco_QQ_mumi_4mom.Eta())<2.1&&Reco_QQ_4mom.Pt()>6.5
  // h1->Scale(0.5);

  h3->Add(h2,h1,1,-1);

  // atlas cuts
  // 0-10%
//  ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>h2","HLTriggers&1==1&&Reco_QQ_trig&1==1&&Reco_QQ_sign==0&&Reco_QQ_type==0&&Reco_QQ_4mom.Pt()>3&&Centrality<4","e");
//  ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>h1","HLTriggers&1==1&&Reco_QQ_trig&1==1&&Reco_QQ_sign!=0&&Reco_QQ_type==0&&Reco_QQ_4mom.Pt()>3&&Centrality<4","e");
  //10-20%
//  ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>h2","HLTriggers&1==1&&Reco_QQ_trig&1==1&&Reco_QQ_sign==0&&Reco_QQ_type==0&&Reco_QQ_4mom.Pt()>3&&Centrality>=4&&Centrality<8","e");
//  ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>h1","HLTriggers&1==1&&Reco_QQ_trig&1==1&&Reco_QQ_sign!=0&&Reco_QQ_type==0&&Reco_QQ_4mom.Pt()>3&&Centrality>=4&&Centrality<8","e");
//  //20-40%
//  ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>h2","HLTriggers&1==1&&Reco_QQ_trig&1==1&&Reco_QQ_sign==0&&Reco_QQ_type==0&&Reco_QQ_4mom.Pt()>3&&Centrality>=8&&Centrality<16","e");
//  ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>h1","HLTriggers&1==1&&Reco_QQ_trig&1==1&&Reco_QQ_sign!=0&&Reco_QQ_type==0&&Reco_QQ_4mom.Pt()>3&&Centrality>=8&&Centrality<16","e");
//  //40-80%
  // ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>h2","HLTriggers&1==1&&Reco_QQ_trig&1==1&&Reco_QQ_sign==0&&Reco_QQ_type==0&&Reco_QQ_4mom.Pt()>3&&Centrality>=16&&Centrality<32","e");
  // ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>h1","HLTriggers&1==1&&Reco_QQ_trig&1==1&&Reco_QQ_sign!=0&&Reco_QQ_type==0&&Reco_QQ_4mom.Pt()>3&&Centrality>=16&&Centrality<32","e");

  /*
  ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.Pt()>>hPt","(HLTriggers&1)==1&&(Reco_QQ_trig&1)==1&&Reco_QQ_sign==0&&Reco_QQ_type==0&&Reco_QQ_VtxProb>0.01&&Reco_QQ_4mom.M()>=2.9&&Reco_QQ_4mom.M()<3.3&&abs(Reco_QQ_4mom.Rapidity())>=1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Centrality>=0&&Centrality<40","e");//&&Centrality<8&&abs(Reco_QQ_4mom.Rapidity())>=1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4
  ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.Pt()>>hPt2","(HLTriggers&1)==1&&(Reco_QQ_trig&1)==1&&Reco_QQ_sign!=0&&Reco_QQ_type==0&&Reco_QQ_VtxProb>0.01&&Reco_QQ_4mom.M()>=2.9&&Reco_QQ_4mom.M()<3.3&&abs(Reco_QQ_4mom.Rapidity())>=1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Centrality>=0&&Centrality<40","e");//&&Centrality<8&&abs(Reco_QQ_4mom.Rapidity())>=1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4

  ((TTree*) gROOT->FindObject("myTree"))->Draw("abs(Reco_QQ_4mom.Rapidity())>>hRap","(HLTriggers&1)==1&&(Reco_QQ_trig&1)==1&&Reco_QQ_sign==0&&Reco_QQ_type==0&&Reco_QQ_VtxProb>0.01&&Reco_QQ_4mom.M()>=2.9&&Reco_QQ_4mom.M()<3.3&&Reco_QQ_4mom.Pt()>6.5&&Reco_QQ_4mom.Pt()<30&&abs(Reco_QQ_4mom.Rapidity())<2.4","e");//&&Centrality<8&&abs(Reco_QQ_4mom.Rapidity())>=1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4
  ((TTree*) gROOT->FindObject("myTree"))->Draw("abs(Reco_QQ_4mom.Rapidity())>>hRap2","(HLTriggers&1)==1&&(Reco_QQ_trig&1)==1&&Reco_QQ_sign!=0&&Reco_QQ_type==0&&Reco_QQ_VtxProb>0.01&&Reco_QQ_4mom.M()>=2.9&&Reco_QQ_4mom.M()<3.3&&Reco_QQ_4mom.Pt()>6.5&&Reco_QQ_4mom.Pt()<30&&abs(Reco_QQ_4mom.Rapidity())<2.4","e");//&&Centrality<8&&abs(Reco_QQ_4mom.Rapidity())>=1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4
  */
  if (!isHI) 
    h2->GetYaxis()->SetRangeUser(2.0,20000.0);
  h2->GetXaxis()->SetRangeUser(2.0,3.99);

  f1->SetParameters(1000,3.1,0.050);

  h2->Fit(f1,"ME","",2.6,3.4);
  f2->SetParameter(0, f1->GetParameter(0));
  f2->SetParameter(1, f1->GetParameter(1));
  f2->SetParameter(2, f1->GetParameter(2));
  f2->SetParameter(3, 1.0);
  f2->SetParameter(4, 0.1);

  h2->Fit(f2,"ME","",2.6,3.5);

  if (isHI) {
    f6->SetParameter(0, f2->GetParameter(0));
    f6->SetParameter(1, f2->GetParameter(1));
    f6->SetParameter(2, f2->GetParameter(2));
    f6->SetParameter(3, 1.5);
    f6->SetParameter(4, 2.0);
    f6->SetParameter(5, f2->GetParameter(3));
    f6->SetParameter(6, f2->GetParameter(4));

    h2->Fit(f6,"RME","");

    if (fitPsiP) {
      f4->SetParameter(0, 100);
      f4->SetParameter(1, 3.686);
      f4->SetParameter(2, f2->GetParameter(2)*3.69/3.1);
      f4->SetParameter(3, f2->GetParameter(3));
      f4->SetParameter(4, f2->GetParameter(4));

      h2->Fit(f4,"ME","",3.4,4.0);
    }

    // f5->SetParameter(0, f2->GetParameter(0));
    // f5->SetParameter(1, f2->GetParameter(1));
    // f5->SetParameter(2, f2->GetParameter(2));
    // f5->SetParameter(3, f4->GetParameter(0));
    // //    f5->SetParameter(4, f4->GetParameter(1));
    // //    f5->SetParameter(5, f4->GetParameter(2));
    // //    f5->SetParameter(6, f2->GetParameter(3));
    // //    f5->SetParameter(7, f2->GetParameter(4));
    // f5->SetParameter(4, f2->GetParameter(3));
    // f5->SetParameter(5, f2->GetParameter(4));

    if (fitPsiP) {
      f7->SetParameter(0, f4->GetParameter(0));
      f7->SetParameter(1, f4->GetParameter(1));
      f7->SetParameter(2, f4->GetParameter(2));
      f7->SetParameter(3, f6->GetParameter(3));
      f7->SetParameter(4, f6->GetParameter(4));
      f7->SetParameter(5, f4->GetParameter(3));
      f7->SetParameter(6, f4->GetParameter(4));

      h2->Fit(f6,"RME","");
      h2->Fit(f7,"RME+","");
      /*
      f5a->SetParameter(0, f6->GetParameter(0));
      f5a->SetParameter(1, f6->GetParameter(1));
      f5a->SetParameter(2, f6->GetParameter(2));
      f5a->SetParameter(3, f6->GetParameter(3));
      f5a->SetParameter(4, f6->GetParameter(4));
      f5a->SetParameter(5, f4->GetParameter(0));
      f5a->SetParameter(6, f6->GetParameter(5));
      f5a->SetParameter(7, f6->GetParameter(6));

      h2->Fit(f5a,"ME","",2.6,4.0);
      */
    }
  }
  else {
    f4->SetParameter(0, 100);
    f4->SetParameter(1, 3.686);
    f4->SetParameter(2, f2->GetParameter(2));
    f4->SetParameter(3, f2->GetParameter(3));
    f4->SetParameter(4, f2->GetParameter(4));

    h2->Fit(f4,"ME","",3.4,4.0);

    f5->SetParameter(0, f2->GetParameter(0));
    f5->SetParameter(1, f2->GetParameter(1));
    f5->SetParameter(2, f2->GetParameter(2));
    f5->SetParameter(3, f4->GetParameter(0));
    f5->SetParameter(4, f4->GetParameter(1));
    f5->SetParameter(5, f4->GetParameter(2));
    f5->SetParameter(6, f2->GetParameter(3));
    f5->SetParameter(7, f2->GetParameter(4));

    h2->Fit(f5,"ME","",2.0,4.0);
  }
  h2->Draw();
  h1->Draw("same");

  f3->SetParameters(0, f2->GetParameter(3));
  f3->SetParameters(1, f2->GetParameter(4));

  double signal = f2->Integral(f2->GetParameter(1)-2*f2->GetParameter(3), f2->GetParameter(1)+2*f2->GetParameter(3));
  double bg = f3->Integral(f2->GetParameter(1)-2*f2->GetParameter(3), f2->GetParameter(1)+2*f2->GetParameter(3))/0.050;

  std::cout << "S/B = " << signal << " / " << bg << " = " << signal/bg << std::endl;

  
  TCanvas *c2 = new TCanvas("c2","c2");
  c2->SetLogy();
  c2->cd();
  hPt->Add(hPt2,-1.0);
  hPt->Draw();

  // All+Barrel
  hPt->GetXaxis()->SetRangeUser(6.5,30.0);
  std::cout << "<pT> = "<< hPt->GetMean() << " +/- " << hPt->GetRMS() << " (RMS)" << std::endl;
  hPt->GetXaxis()->SetRangeUser(6.5,10.0);
  std::cout << "<pT> = "<< hPt->GetMean() << " +/- " << hPt->GetRMS() << " (RMS)" << std::endl;
  hPt->GetXaxis()->SetRangeUser(10.0,30.0);
  std::cout << "<pT> = "<< hPt->GetMean() << " +/- " << hPt->GetRMS() << " (RMS)" << std::endl;
  hPt->GetXaxis()->SetRangeUser(0.0,30.0);

  // Overlap
  // hPt->GetXaxis()->SetRangeUser(6.5,30.0);
  // std::cout << "<pT> = "<< hPt->GetMean() << " +/- " << hPt->GetRMS() << " (RMS)" << std::endl;
  // hPt->GetXaxis()->SetRangeUser(5.5,30.0);
  // std::cout << "<pT> = "<< hPt->GetMean() << " +/- " << hPt->GetRMS() << " (RMS)" << std::endl;

  // EndCap
  // hPt->GetXaxis()->SetRangeUser(6.5,30.0);
  // std::cout << "<pT> = "<< hPt->GetMean() << " +/- " << hPt->GetRMS() << " (RMS)" << std::endl;
  // hPt->GetXaxis()->SetRangeUser(3.0,30.0);
  // std::cout << "<pT> = "<< hPt->GetMean() << " +/- " << hPt->GetRMS() << " (RMS)" << std::endl;

  double error;
  double errorLS;

  std::cout << "Yields" << std::endl;
  hPt2->IntegralAndError(66,300,errorLS);
  hPt->IntegralAndError(66,300,error);
  std::cout << "6.5 < pT < 30.0: "<< hPt->IntegralAndError(66,300,error) << " +/- " << error << " +/- " << errorLS << std::endl;
  hPt2->IntegralAndError(66,100,errorLS);
  std::cout << "6.5 <= pT < 10.0: "<< hPt->IntegralAndError(66,100,error) << " +/- " << error << " +/- " << errorLS << std::endl;
  hPt2->IntegralAndError(101,300,errorLS);
  std::cout << "10.0 <= pT < 30.0: "<< hPt->IntegralAndError(101,300,error) << " +/- " << error << " +/- " << errorLS << std::endl;

  hPt2->IntegralAndError(1,65,errorLS);
  std::cout << "pT < 6.5: "<< hPt->IntegralAndError(1,65,error) << " +/- " << error << " +/- " << errorLS << std::endl;
  hPt2->IntegralAndError(31,300,errorLS);
  std::cout << "3.0 <= pT < 30.0: "<< hPt->IntegralAndError(31,300,error) << " +/- " << error << " +/- " << errorLS << std::endl;
  hPt2->IntegralAndError(56,300,errorLS);
  std::cout << "5.5 <= pT < 30.0: "<< hPt->IntegralAndError(56,300,error) << " +/- " << error << " +/- " << errorLS << std::endl;


  // cross check
  std::cout << "Yields" << std::endl;
  std::cout << "3-3.2: "<< h3->IntegralAndError(21,24,error) << " +/- " << error << std::endl;



  TCanvas *c4 = new TCanvas("c4","c4");
  c4->SetLogy();
  c4->cd();
  hRap->Add(hRap2,-1.0);
  hRap->Draw();


  hRap->GetXaxis()->SetRangeUser(0.0,2.4);
  std::cout << "<|y|> = "<< hRap->GetMean() << " +/- " << hRap->GetRMS() << " (RMS)" << std::endl;
  hRap->GetXaxis()->SetRangeUser(0.0,1.2);
  std::cout << "<|y|> = "<< hRap->GetMean() << " +/- " << hRap->GetRMS() << " (RMS)" << std::endl;
  hRap->GetXaxis()->SetRangeUser(1.2,1.6);
  std::cout << "<|y|> = "<< hRap->GetMean() << " +/- " << hRap->GetRMS() << " (RMS)" << std::endl;
  hRap->GetXaxis()->SetRangeUser(1.6,2.4);
  std::cout << "<|y|> = "<< hRap->GetMean() << " +/- " << hRap->GetRMS() << " (RMS)" << std::endl;
  hRap->GetXaxis()->SetRangeUser(0.0,2.4);


  return;
}
