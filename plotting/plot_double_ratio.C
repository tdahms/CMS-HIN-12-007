#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TPaveStats.h"
#include "TLatex.h"
#include "TLegend.h"


void plot_double_ratio(bool isPaper=false, bool allRapidity=false, bool plotGlobalPP = false, bool savePlots=false)
{
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);

  // 0-20%, 20-40%, 40-100%
  double Npart[3] = {308.3853, 158.5653, 32.7694};
  double Npart2[3] = {318.3853, 158.5653, 22.7694};
  double Npart_err[3] = {0,0,0};
  double Npart_sys[3] = {10,10,10};


  // 20140225
  double ratio_fwd_pt3065[1] = {10.131};
  double ratio_fwd_pt3065_err[1] = {1.802};
  double ratio_fwd_pt3065_sys[1] = {2.773};
  double ratio_fwd_pt3065_glb[1] = {6.698};
    
  double ratio_fwd_pt6530[1] = {1.803};
  double ratio_fwd_pt6530_err[1] = {0.598};
  double ratio_fwd_pt6530_sys[1] = {0.479};
  double ratio_fwd_pt6530_glb[1] = {0.957};

  double ratio_fwd_pt3030[3] = {2.231860, 0.851082, 0.833084};
  double ratio_fwd_pt3030_err[3] = {0.501988, 0.441254, 0.368925};
  double ratio_fwd_pt3030_sys[3] = {0.559082, 0.821231, 0.362263};
  double ratio_fwd_pt3030_sysHi[3] = {0.559082, 0.000000, 0.362263};
  double ratio_fwd_pt3030_sysLo[3] = {0.312271, 0.821231, 0.030393};
  double ratio_fwd_pt3030_glb[3] = {0.212159, 0.080903, 0.079192};
  double ratio_fwd_pt3030_glbHi[3] = {0.213386, 0.081371, 0.079650};
  double ratio_fwd_pt3030_glbLo[3] = {0.399951, 0.152515, 0.149289};

  double ratio_mid_pt6530[3] = {0.430459, 0.684042, 0.00};
  double ratio_mid_pt6530_err[3] = {0.203423, 0.201516, 0.161972};
  double ratio_mid_pt6530_sys[3] = {0.168448, 0.052673, 0.036599};
  double ratio_mid_pt6530_sysHi[3] = {0.002913, 0.052673, 0.036599};
  double ratio_mid_pt6530_sysLo[3] = {0.168448, 0.006997, 0.027120};
  double ratio_mid_pt6530_glb[3] = {0.029885, 0.047490, 0.000631};
  double ratio_mid_pt6530_glbHi[3] = {0.030242, 0.048058, 0.000638};
  double ratio_mid_pt6530_glbLo[3] = {0.025392, 0.040350, 0.000536};

  double ratio_all_pt6530[3] = {0.751, 0.525, 0.344};
  double ratio_all_pt6530_err[3] = {0.177, 0.172, 0.162};
  //  double ratio_all_pt6530_sys[3] = {0.031, 0.024, 0.251};
  double ratio_all_pt6530_sys[3] = {0.068, 0.044, 0.050};
  double ratio_all_pt6530_glb[3] = {0.145, 0.101, 0.066};
  
  TGraphErrors *gr_fwd_pt3065 = new TGraphErrors(3,Npart, ratio_fwd_pt3065, Npart_err, ratio_fwd_pt3065_err);
  TGraphErrors *gr_fwd_pt3065_sys = new TGraphErrors(3,Npart, ratio_fwd_pt3065, Npart_sys, ratio_fwd_pt3065_sys);
  TGraphErrors *gr_fwd_pt3065_glb = new TGraphErrors(3,Npart2, ratio_fwd_pt3065, Npart_sys, ratio_fwd_pt3065_glb);
  TGraphErrors *gr_fwd_pt3065P = new TGraphErrors(3,Npart, ratio_fwd_pt3065, Npart_err, Npart_err);
  gr_fwd_pt3065->SetName("gr_fwd_pt3065");
  gr_fwd_pt3065_sys->SetName("gr_fwd_pt3065_sys");
  gr_fwd_pt3065_glb->SetName("gr_fwd_pt3065_glb");
  gr_fwd_pt3065P->SetName("gr_fwd_pt3065P");
									               
  TGraphErrors *gr_fwd_pt3030 = new TGraphErrors(3,Npart, ratio_fwd_pt3030, Npart_err, ratio_fwd_pt3030_err);
  // TGraphErrors *gr_fwd_pt3030_sys = new TGraphErrors(3,Npart, ratio_fwd_pt3030, Npart_sys, ratio_fwd_pt3030_sys);
  // TGraphErrors *gr_fwd_pt3030_glb = new TGraphErrors(3,Npart2, ratio_fwd_pt3030, Npart_sys, ratio_fwd_pt3030_glb);
  TGraphAsymmErrors *gr_fwd_pt3030_sys = new TGraphAsymmErrors(3,Npart, ratio_fwd_pt3030, Npart_sys, Npart_sys, ratio_fwd_pt3030_sysLo, ratio_fwd_pt3030_sysHi);
  TGraphAsymmErrors *gr_fwd_pt3030_glb = new TGraphAsymmErrors(3,Npart2, ratio_fwd_pt3030, Npart_sys, Npart_sys, ratio_fwd_pt3030_glbLo, ratio_fwd_pt3030_glbHi);
  TGraphErrors *gr_fwd_pt3030P = new TGraphErrors(3,Npart, ratio_fwd_pt3030, Npart_err, Npart_err);
  gr_fwd_pt3030->SetName("gr_fwd_pt3030");
  gr_fwd_pt3030_sys->SetName("gr_fwd_pt3030_sys");
  gr_fwd_pt3030_glb->SetName("gr_fwd_pt3030_glb");
  gr_fwd_pt3030P->SetName("gr_fwd_pt3030P");
				               
  TGraphErrors *gr_fwd_pt6530 = new TGraphErrors(3,Npart, ratio_fwd_pt6530, Npart_err, ratio_fwd_pt6530_err);
  TGraphErrors *gr_fwd_pt6530_sys = new TGraphErrors(3,Npart, ratio_fwd_pt6530, Npart_sys, ratio_fwd_pt6530_sys);
  TGraphErrors *gr_fwd_pt6530_glb = new TGraphErrors(3,Npart2, ratio_fwd_pt6530, Npart_sys, ratio_fwd_pt6530_glb);
  TGraphErrors *gr_fwd_pt6530P = new TGraphErrors(3,Npart, ratio_fwd_pt6530, Npart_err, Npart_err);
  gr_fwd_pt6530->SetName("gr_fwd_pt6530");
  gr_fwd_pt6530_sys->SetName("gr_fwd_pt6530_sys");
  gr_fwd_pt6530_sys->SetName("gr_fwd_pt6530_glb");
  gr_fwd_pt6530P->SetName("gr_fwd_pt6530P");
									               
  TGraphErrors *gr_mid_pt6530 = new TGraphErrors(3,Npart, ratio_mid_pt6530, Npart_err, ratio_mid_pt6530_err);
  // TGraphErrors *gr_mid_pt6530_sys = new TGraphErrors(3,Npart, ratio_mid_pt6530, Npart_sys, ratio_mid_pt6530_sys);
  // TGraphErrors *gr_mid_pt6530_glb = new TGraphErrors(3,Npart2, ratio_mid_pt6530, Npart_sys, ratio_mid_pt6530_glb);
  TGraphAsymmErrors *gr_mid_pt6530_sys = new TGraphAsymmErrors(3,Npart, ratio_mid_pt6530, Npart_sys, Npart_sys, ratio_mid_pt6530_sysLo, ratio_mid_pt6530_sysHi);
  TGraphAsymmErrors *gr_mid_pt6530_glb = new TGraphAsymmErrors(3,Npart2, ratio_mid_pt6530, Npart_sys, Npart_sys, ratio_mid_pt6530_glbLo, ratio_mid_pt6530_glbHi);
  TGraphErrors *gr_mid_pt6530P = new TGraphErrors(3,Npart, ratio_mid_pt6530, Npart_err, Npart_err);
  gr_mid_pt6530->SetName("gr_mid_pt6530");
  gr_mid_pt6530_sys->SetName("gr_mid_pt6530_sys");
  gr_mid_pt6530_sys->SetName("gr_mid_pt6530_glb");
  gr_mid_pt6530P->SetName("gr_mid_pt6530P");

  TGraphErrors *gr_all_pt6530 = new TGraphErrors(3,Npart, ratio_all_pt6530, Npart_err, ratio_all_pt6530_err);
  TGraphErrors *gr_all_pt6530_sys = new TGraphErrors(3,Npart, ratio_all_pt6530, Npart_sys, ratio_all_pt6530_sys);
  TGraphErrors *gr_all_pt6530_glb = new TGraphErrors(3,Npart2, ratio_all_pt6530, Npart_sys, ratio_all_pt6530_glb);
  TGraphErrors *gr_all_pt6530P = new TGraphErrors(3,Npart, ratio_all_pt6530, Npart_err, Npart_err);
  gr_all_pt6530->SetName("gr_all_pt6530");
  gr_all_pt6530_sys->SetName("gr_all_pt6530_sys");
  gr_all_pt6530_sys->SetName("gr_all_pt6530_glb");
  gr_all_pt6530P->SetName("gr_all_pt6530P");

  gr_fwd_pt3065->SetMarkerSize(1.2);
  gr_fwd_pt3030->SetMarkerSize(1.2);
  gr_fwd_pt6530->SetMarkerSize(1.2);
  gr_mid_pt6530->SetMarkerSize(2.0);
  gr_all_pt6530->SetMarkerSize(2.0);

  gr_fwd_pt3065P->SetMarkerSize(1.2);
  gr_fwd_pt3030P->SetMarkerSize(1.2);
  gr_fwd_pt6530P->SetMarkerSize(1.2);
  gr_mid_pt6530P->SetMarkerSize(2.0);
  gr_all_pt6530P->SetMarkerSize(2.0);

  if (isPaper) {
    gr_fwd_pt3065->SetMarkerColor(kRed+2);
    gr_fwd_pt3065->SetLineColor(kRed+2);
    gr_fwd_pt3030->SetMarkerColor(kRed+2);
    gr_fwd_pt3030->SetLineColor(kRed+2);
    gr_fwd_pt6530->SetMarkerColor(kBlue+2);
    gr_fwd_pt6530->SetLineColor(kBlue+2);

    gr_fwd_pt3065_sys->SetLineColor(kRed+2);
    gr_fwd_pt3030_sys->SetLineColor(kRed+2);
    gr_fwd_pt6530_sys->SetLineColor(kBlue+2);

    gr_fwd_pt3065_glb->SetFillColor(kRed+2);
    gr_fwd_pt3030_glb->SetFillColor(kRed+2);
    gr_fwd_pt6530_glb->SetFillColor(kBlue+2);
  }
  else {
    gr_fwd_pt3065->SetMarkerColor(kRed);
    gr_fwd_pt3065->SetLineColor(kRed);
    gr_fwd_pt3030->SetMarkerColor(kRed);
    gr_fwd_pt3030->SetLineColor(kRed);
    gr_fwd_pt6530->SetMarkerColor(kBlue);
    gr_fwd_pt6530->SetLineColor(kBlue);

    gr_fwd_pt3065_sys->SetLineColor(kRed);
    gr_fwd_pt3030_sys->SetLineColor(kRed);
    gr_fwd_pt6530_sys->SetLineColor(kBlue);

    gr_fwd_pt3065_glb->SetFillColor(kRed);
    gr_fwd_pt3030_glb->SetFillColor(kRed);
    gr_fwd_pt6530_glb->SetFillColor(kBlue);
  }

  gr_fwd_pt3065_sys->SetFillColor(kRed-9);
  gr_fwd_pt3030_sys->SetFillColor(kRed-9);
  gr_fwd_pt6530_sys->SetFillColor(kBlue-9);

  gr_mid_pt6530->SetMarkerColor(kGreen+2);
  gr_mid_pt6530->SetLineColor(kGreen+2);
  gr_mid_pt6530_sys->SetLineColor(kGreen+2);
  gr_mid_pt6530_sys->SetFillColor(kGreen-9);
  gr_mid_pt6530_glb->SetFillColor(kGreen+2);

  gr_all_pt6530->SetMarkerColor(kGreen+2);
  gr_all_pt6530->SetLineColor(kGreen+2);
  gr_all_pt6530_sys->SetLineColor(kGreen+2);
  gr_all_pt6530_sys->SetFillColor(kGreen-9);
  gr_all_pt6530_glb->SetFillColor(kGreen+2);

  gr_fwd_pt6530->SetMarkerStyle(21);
  gr_mid_pt6530->SetMarkerStyle(33);
  gr_all_pt6530->SetMarkerStyle(33);

  gr_fwd_pt3065P->SetMarkerStyle(24);
  gr_fwd_pt3030P->SetMarkerStyle(24);
  gr_fwd_pt6530P->SetMarkerStyle(25);
  gr_mid_pt6530P->SetMarkerStyle(27);
  gr_all_pt6530P->SetMarkerStyle(27);

  gr_fwd_pt3065_glb->SetLineColor(kRed+2);
  gr_fwd_pt3030_glb->SetLineColor(kRed+2);
  gr_fwd_pt6530_glb->SetLineColor(kBlue+2);
  gr_mid_pt6530_glb->SetLineColor(kGreen+2);
  gr_all_pt6530_glb->SetLineColor(kGreen+2);

  gr_fwd_pt3065_glb->SetFillStyle(0);
  gr_fwd_pt3030_glb->SetFillStyle(0);
  gr_fwd_pt6530_glb->SetFillStyle(0);
  gr_mid_pt6530_glb->SetFillStyle(0);
  gr_all_pt6530_glb->SetFillStyle(0);

  TGraphErrors *gr_fwd_pt3065_glb2 = (TGraphErrors*) gr_fwd_pt3065_glb->Clone("gr_fwd_pt3065_glb2");
  TGraphErrors *gr_fwd_pt3030_glb2 = (TGraphErrors*) gr_fwd_pt3030_glb->Clone("gr_fwd_pt3030_glb2");
  TGraphErrors *gr_fwd_pt6530_glb2 = (TGraphErrors*) gr_fwd_pt6530_glb->Clone("gr_fwd_pt6530_glb2");
  TGraphErrors *gr_mid_pt6530_glb2 = (TGraphErrors*) gr_mid_pt6530_glb->Clone("gr_mid_pt6530_glb2");
  TGraphErrors *gr_all_pt6530_glb2 = (TGraphErrors*) gr_all_pt6530_glb->Clone("gr_all_pt6530_glb2");

  gr_fwd_pt3065_glb2->SetFillStyle(1001);
  gr_fwd_pt3030_glb2->SetFillStyle(1001);
  gr_fwd_pt6530_glb2->SetFillStyle(1001);
  gr_mid_pt6530_glb2->SetFillStyle(1001);
  gr_all_pt6530_glb2->SetFillStyle(1001);

  gr_fwd_pt3065_glb2->SetFillColor(kGray);
  gr_fwd_pt3030_glb2->SetFillColor(kGray);
  gr_fwd_pt6530_glb2->SetFillColor(kGray);
  gr_mid_pt6530_glb2->SetFillColor(kGray);
  gr_all_pt6530_glb2->SetFillColor(kGray);

  double pp_all_pt6530_err = 0.194;//update
  double pp_mid_pt6530_err = 0.210;//update
  double pp_fwd_pt3030_err = 0.492;
  double pp_fwd_pt3065_err = 0.670;//update
  double pp_fwd_pt6530_err = 0.503;//update

  TBox *pp_fwd_pt3065 = new TBox(0.0,1.0-pp_fwd_pt3065_err,7.0,1.0+pp_fwd_pt3065_err);
  pp_fwd_pt3065->SetFillColor(kRed-9);
  TBox *pp_fwd_pt3030 = new TBox(0.0,1.0-pp_fwd_pt3030_err,10.0,1.0+pp_fwd_pt3030_err);
  pp_fwd_pt3030->SetFillColor(kRed-9);
  TBox *pp_fwd_pt6530 = new TBox(7.0,1.0-pp_fwd_pt6530_err,14.0,1.0+pp_fwd_pt6530_err);
  pp_fwd_pt6530->SetFillColor(kBlue-9);
  TBox *pp_mid_pt6530 = new TBox(14.0,1.0-pp_mid_pt6530_err,21.0,1.0+pp_mid_pt6530_err);
  pp_mid_pt6530->SetFillColor(kGreen-9);

  TBox *pp_fwd2_pt6530 = new TBox(0.0,1.0-pp_fwd_pt6530_err,7.0,1.0+pp_fwd_pt6530_err);
  pp_fwd2_pt6530->SetFillColor(kBlue-9);
  TBox *pp_mid2_pt6530 = new TBox(7.0,1.0-pp_mid_pt6530_err,14.0,1.0+pp_mid_pt6530_err);
  pp_mid2_pt6530->SetFillColor(kGreen-9);

  TBox *pp_all_pt6530;
  if (allRapidity) {
    pp_all_pt6530 = new TBox(0.0,1.0-pp_all_pt6530_err,10.0,1.0+pp_all_pt6530_err);
  }
  else {
    pp_all_pt6530 = new TBox(10.0,1.0-pp_all_pt6530_err,20.0,1.0+pp_all_pt6530_err);
  }
  pp_all_pt6530->SetFillColor(kGreen-9);



  TF1 *f1 = new TF1("f1","1",0,400);
  f1->SetLineWidth(1);
  f1->GetXaxis()->SetTitle("N_{part}");
  f1->GetYaxis()->SetTitle("[#psi(2S) #/J/#psi]_{PbPb} #/[#psi(2S) #/J/#psi]_{pp}");
  f1->GetYaxis()->SetTitleOffset(1.1);
  if (allRapidity)
    f1->GetYaxis()->SetRangeUser(0.0,1.5);
  else
    f1->GetYaxis()->SetRangeUser(0.0,12.0);
  f1->GetXaxis()->CenterTitle(kTRUE);

  TCanvas *c1 = new TCanvas("c1","c1");
  f1->Draw();

  if (allRapidity) {
    if (plotGlobalPP)
      pp_all_pt6530->Draw();
    else
      gr_all_pt6530_glb->Draw("3");
    f1->Draw("same");
    gr_all_pt6530_sys->Draw("5");
    gr_all_pt6530->Draw("P");
    gr_all_pt6530P->Draw("PX");
  }
  else {
    if (plotGlobalPP) {
      pp_fwd_pt3065->Draw();
      pp_fwd_pt6530->Draw();
      pp_mid_pt6530->Draw();
    }
    else {
      gr_fwd_pt3065_glb2->Draw("3");
      gr_fwd_pt6530_glb2->Draw("3");
      gr_mid_pt6530_glb2->Draw("3");
    }

    f1->Draw("same");

    gr_fwd_pt3065_sys->Draw("5");
    gr_fwd_pt6530_sys->Draw("5");
    gr_mid_pt6530_sys->Draw("5");

    gr_fwd_pt3065->Draw("P");
    gr_fwd_pt6530->Draw("P");
    gr_mid_pt6530->Draw("P");

    gr_fwd_pt3065P->Draw("PX");
    gr_fwd_pt6530P->Draw("PX");
    gr_mid_pt6530P->Draw("PX");

    gr_fwd_pt3065_glb->Draw("3");
    gr_fwd_pt6530_glb->Draw("3");
    gr_mid_pt6530_glb->Draw("3");
  }

  TLegend *leg1;
  if (allRapidity)
    leg1 = new TLegend(0.15,0.78,0.72,0.88);
  else
    leg1 = new TLegend(0.15,0.68,0.72,0.88);

  leg1->SetFillStyle(0);
  leg1->SetFillColor(0);
  leg1->SetBorderSize(0);
  leg1->SetMargin(0.15);
  leg1->SetTextSize(0.035);
  leg1->SetHeader("PbPb  #sqrt{s_{NN}} = 2.76 TeV");
  if (allRapidity)
    leg1->AddEntry(gr_all_pt6530,"6.5 <  p_{T} < 30 GeV/c, |y| < 2.4","P");
  else {
    leg1->AddEntry(gr_fwd_pt3065,"3 < p_{T} < 6.5 GeV/c, 1.6 < |y| < 2.4","P");
    leg1->AddEntry(gr_fwd_pt6530,"6.5 < p_{T} < 30 GeV/c, 1.6 < |y| < 2.4","P");
    leg1->AddEntry(gr_mid_pt6530,"6.5 < p_{T} < 30 GeV/c, |y| < 1.6","P");
  }
  leg1->Draw();
  
  TLatex *pre = new TLatex(18,11,"CMS Preliminary");
  pre->SetTextFont(42);
  pre->SetTextSize(0.05);
  if (allRapidity)
    pre->SetY(1.375);

  //  pre->Draw();

  gPad->RedrawAxis();


  if (savePlots) {
    if (allRapidity) {
      c1->SaveAs("double_ratio_ppGlobal_all.pdf");
      c1->SaveAs("double_ratio_ppGlobal_all.png");
    }
    else {
      c1->SaveAs("double_ratio_ppGlobal.pdf");
      c1->SaveAs("double_ratio_ppGlobal.png");
    }
  }



  TF1 *f2 = new TF1("f2","1",0,400);
  f2->SetLineWidth(1);
  f2->GetXaxis()->SetTitle("N_{part}");
  f2->GetYaxis()->SetTitle("[#psi(2S) #/J/#psi]_{PbPb} #/[#psi(2S) #/J/#psi]_{pp}");
  f2->GetYaxis()->SetTitleOffset(1.1);
  f2->GetYaxis()->SetRangeUser(0.0,4.0);
  f2->GetXaxis()->CenterTitle(kTRUE);

  TCanvas *c2 = new TCanvas("c2","c2");
  f2->Draw();

  if (plotGlobalPP) {
    pp_fwd2_pt6530->Draw();
    pp_mid2_pt6530->Draw();
  }
  else {
    gr_fwd_pt6530_glb2->Draw("3");
    gr_mid_pt6530_glb2->Draw("3");
  }
  f2->Draw("same");

  gr_fwd_pt6530_sys->Draw("5");
  gr_mid_pt6530_sys->Draw("5");

  gr_fwd_pt6530->Draw("P");
  gr_mid_pt6530->Draw("P");

  gr_fwd_pt6530P->Draw("PX");
  gr_mid_pt6530P->Draw("PX");

  gr_fwd_pt6530_glb->Draw("3");
  gr_mid_pt6530_glb->Draw("3");

  TLegend *leg2 = new TLegend(0.15,0.58,0.72,0.88);
  //leg1 = new TLegend(0.15,0.68,0.72,0.88);

  leg2->SetFillStyle(0);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetMargin(0.15);
  leg2->SetTextSize(0.035);
  leg2->SetHeader("PbPb  #sqrt{s_{NN}} = 2.76 TeV");
  leg2->AddEntry(gr_all_pt6530,"6.5 <  p_{T} < 30 GeV/c, |y| < 1.6","P");
  leg2->AddEntry(gr_fwd_pt6530,"6.5 < p_{T} < 30 GeV/c, 1.6 < |y| < 2.4","P");
  leg2->Draw();
  
  TLatex *pre2 = new TLatex(18,3.667,"CMS Preliminary");
  pre2->SetTextFont(42);
  pre2->SetTextSize(0.05);
  pre2->Draw();

  gPad->RedrawAxis();

  if (savePlots) {
    c2->SaveAs("double_ratio_ppGlobal_highpt.pdf");
    c2->SaveAs("double_ratio_ppGlobal_highpt.png");
  }

  TF1 *f3 = new TF1("f3","1",0,400);
  f3->SetLineWidth(1);
  f3->GetXaxis()->SetTitle("N_{part}");
  f3->GetYaxis()->SetTitle("[ #psi (2S) #/J/#psi ]_{PbPb} #/[ #psi (2S) #/J/#psi ]_{pp}");
  f3->GetYaxis()->SetTitleOffset(1.1);
  f3->GetYaxis()->SetRangeUser(0.0,3.0);
  f3->GetXaxis()->CenterTitle(kTRUE);

  TF1 *f3a = new TF1("f3a","1",0,400);
  f3a->SetLineWidth(1);
  f3a->GetXaxis()->SetTitle("N_{part}");
  f3a->GetYaxis()->SetTitle("[ #psi (2S) #/J/#psi ]_{PbPb} #/[ #psi (2S) #/J/#psi ]_{pp}");
  f3a->GetYaxis()->SetTitleOffset(1.1);
  f3a->GetYaxis()->SetRangeUser(0.0,3);
  f3a->GetXaxis()->CenterTitle(kTRUE);

  //  TCanvas *c3 = new TCanvas("c3","c3",600,800);
  TCanvas *c3 = new TCanvas("c3","c3");//,1000,600);
  //  c3->Divide(1,2);

  // c3->GetPad(1)->SetPad(0.0,0.5,1.0,1.0);
  // c3->GetPad(2)->SetPad(0.0,0.0,1.0,0.5);
  // c3->GetPad(1)->SetBottomMargin(0.0);
  // c3->GetPad(2)->SetTopMargin(0.0);

  // c3->GetPad(1)->SetPad(0.0,0.0,0.5,1.0);
  // c3->GetPad(2)->SetPad(0.5,0.0,1.0,1.0);
  c3->cd(1);
  //  c3->SetLogy();
  f3->Draw();

  if (plotGlobalPP)
    pp_fwd_pt3030->Draw();
  else
    gr_fwd_pt3030_glb2->Draw("3");

  f3->Draw("same");

  gr_fwd_pt3030_sys->Draw("5");
  gr_fwd_pt3030->Draw("P");
  gr_fwd_pt3030P->Draw("PX");
  gr_fwd_pt3030_glb->Draw("3");

  TLegend *leg3;
  if (plotGlobalPP)
    leg3 = new TLegend(0.15,0.77,0.72,0.88);
  else
    leg3 = new TLegend(0.15,0.78,0.72,0.94);

  //leg1 = new TLegend(0.15,0.68,0.72,0.88);

  leg3->SetFillStyle(0);
  leg3->SetFillColor(0);
  leg3->SetBorderSize(0);
  leg3->SetMargin(0.15);
  //  leg3->SetTextSize(0.04);
  leg3->SetTextSize(0.035);
  leg3->SetHeader("PbPb  #sqrt{s_{NN}} = 2.76 TeV");
  leg3->AddEntry(gr_fwd_pt3030,"3 < p_{T} < 30 GeV/c, 1.6 < |y| < 2.4","P");
  leg3->AddEntry(gr_fwd_pt3030_glb2,"pp uncertainty (global)","F");
  leg3->Draw();
  
  //  pre->SetTextSize(0.06);
  //  pre->SetY(9.1667);
  pre->SetY(0.4);
  pre->SetX(210);
  //  pre->Draw();

  TLine *l1 = new TLine(400,1.0,490,7.1);
  //  l1->Draw("same");
  TLine *l2 = new TLine(400,0.0,490.0,0.0);
  //  l2->Draw("same");

  gPad->RedrawAxis();


  TCanvas *c4 = new TCanvas("c4","c4");//,1000,600);
  //  c3->cd(2);
  f3a->Draw();
  if (plotGlobalPP)
    pp_mid_pt6530->Draw();
  else
    gr_mid_pt6530_glb2->Draw("3");
  f3a->Draw("same");

  gr_mid_pt6530_sys->Draw("5");
  gr_mid_pt6530->Draw("P");
  gr_mid_pt6530P->Draw("PX");
  gr_mid_pt6530_glb->Draw("3");


  TLegend *leg3a;
  if (plotGlobalPP)
    leg3a = new TLegend(0.15,0.77,0.72,0.88);
  else
    leg3a = new TLegend(0.15,0.78,0.72,0.94);

  //leg1 = new TLegend(0.15,0.68,0.72,0.88);

  leg3a->SetFillStyle(0);
  leg3a->SetFillColor(0);
  leg3a->SetBorderSize(0);
  leg3a->SetMargin(0.15);
  //  leg3a->SetTextSize(0.04);
  leg3a->SetTextSize(0.035);
  leg3a->SetHeader("PbPb  #sqrt{s_{NN}} = 2.76 TeV");

  leg3a->AddEntry(gr_mid_pt6530,"6.5 < p_{T} < 30 GeV/c, |y| < 1.6","P");
  leg3a->AddEntry(gr_mid_pt6530_glb2,"pp uncertainty (global)","F");
  leg3a->Draw();

  TLatex *pre3 = (TLatex*)pre->Clone();
  pre3->SetY(0.056);
  //  pre3->Draw();

  TLine *l1a = new TLine(-90,0.15,0.0,1.0);
  //  l1a->Draw("same");

  TLine *l2a = new TLine(-90,0.0,-18.0,0.0);
  //  l2a->Draw("same");


  gPad->RedrawAxis();

  if (savePlots) {
    c3->SaveAs("double_ratio_ppGlobal_fwd_pp2013_PbPbRegit_Asym.pdf");
    c3->SaveAs("double_ratio_ppGlobal_fwd_pp2013_PbPbRegit_Asym.png");
    c4->SaveAs("double_ratio_ppGlobal_mid_pp2013_PbPbRegit_Asym.pdf");
    c4->SaveAs("double_ratio_ppGlobal_mid_pp2013_PbPbRegit_Asym.png");
  }


  return;
}
