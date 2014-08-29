#include <iostream>
#include <fstream>
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
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TPaveStats.h"
#include "TLatex.h"
#include "TLegend.h"

TMultiGraph *read_alice();
TGraphErrors *gr_alice_pt03P;
TGraphErrors *gr_alice_pt38;
TGraphAsymmErrors *gr_alice_pt38_bar;


TMultiGraph *read_pas();
TGraphErrors *gr_pas_fwd_pt3030;
TGraphErrors *gr_pas_fwd_pt3030P;
TGraphErrors *gr_pas_mid_pt6530;
TGraphErrors *gr_pas_fwd_pt3030_glb;

void plot_double_ratio(bool isPaper=false, bool plotGlobalPP = false, bool overlay=true, bool showSplit=false, bool overlayAlice=false, bool savePlots=false)
{
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  
  bool preliminary = true;
  bool overlayPAS = false;

  double split = 0.86;
  double zoom = split/(1-split);

  // 0-20%, 20-40%, 40-100%
  double Npart[3] = {308.3853, 158.5653, 32.7694};
  double Npart2[3] = {313.3853, 163.5653, 37.7694};
  double Npart_err[3] = {0,0,0};
  double Npart_sys[3] = {10,10,10};

  double ratio_fwd_pt3030[3] = {};
  double ratio_fwd_pt3030_err[3] = {};
  double ratio_fwd_pt3030_sys[3] = {};
  double ratio_fwd_pt3030_glb[3] = {};

  // centrality integrated points
  double ratio_fwd_pt3030_mb[1] = {};
  double ratio_fwd_pt3030_mb_err[1] = {};
  double ratio_fwd_pt3030_mb_sys[1] = {};

  double Npart3[2] = {355.4, 261.4};
  double ratio_fwd_pt3030_split[2] = {};
  double ratio_fwd_pt3030_split_err[2] = {};

  std::ifstream ifs("fwd.txt");
  int i=0;
  if (ifs.is_open()) {
    while (ifs.good()) {
      ifs >> ratio_fwd_pt3030[i] >> ratio_fwd_pt3030_err[i] >> ratio_fwd_pt3030_sys[i] >> ratio_fwd_pt3030_glb[i];
      i++;
    }
    ifs.close();
  }

  ifs.open("fwd_mb.txt");
  i=0;
  if (ifs.is_open()) {
    while (ifs.good()) {
      ifs >> ratio_fwd_pt3030_mb[i] >> ratio_fwd_pt3030_mb_err[i] >> ratio_fwd_pt3030_mb_sys[i];
      i++;
    }
    ifs.close();
  }

  ifs.open("fwd_split.txt");
  i=0;
  if (ifs.is_open()) {
    while (ifs.good()) {
      ifs >> ratio_fwd_pt3030_split[i] >> ratio_fwd_pt3030_split_err[i];
      i++;
    }
    ifs.close();
  }

  double ratio_mid_pt6530[3] =     {};
  double ratio_mid_pt6530_err[3] = {};
  double ratio_mid_pt6530_sys[3] = {};
  double ratio_mid_pt6530_glb[3] = {};

  // centrality integrated points
  double ratio_mid_pt6530_mb[1] = {};
  double ratio_mid_pt6530_mb_err[1] = {};
  double ratio_mid_pt6530_mb_sys[1] = {};

  ifs.open("mid.txt");
  i=0;
  if (ifs.is_open()) {
    while (ifs.good()) {
      ifs >> ratio_mid_pt6530[i] >> ratio_mid_pt6530_err[i] >> ratio_mid_pt6530_sys[i] >> ratio_mid_pt6530_glb[i];
      i++;
    }
    ifs.close();
  }

  ifs.open("mid_mb.txt");
  i=0;
  if (ifs.is_open()) {
    while (ifs.good()) {
      ifs >> ratio_mid_pt6530_mb[i] >> ratio_mid_pt6530_mb_err[i] >> ratio_mid_pt6530_mb_sys[i];
      i++;
    }
    ifs.close();
  }
  
  TGraphErrors *gr_fwd_pt3030;
  TGraphErrors *gr_fwd_pt3030_sys;
  TGraphErrors *gr_fwd_pt3030_glb;
  TGraphErrors *gr_fwd_pt3030P;

  if (overlay) {
    gr_fwd_pt3030 = new TGraphErrors(3,Npart2, ratio_fwd_pt3030, Npart_err, ratio_fwd_pt3030_err);
    gr_fwd_pt3030_sys = new TGraphErrors(3,Npart2, ratio_fwd_pt3030, Npart_sys, ratio_fwd_pt3030_sys);
    gr_fwd_pt3030_glb = new TGraphErrors(3,Npart2, ratio_fwd_pt3030, Npart_sys, ratio_fwd_pt3030_glb);
    gr_fwd_pt3030P = new TGraphErrors(3,Npart2, ratio_fwd_pt3030, Npart_err, Npart_err);
  }
  else {
    gr_fwd_pt3030 = new TGraphErrors(3,Npart, ratio_fwd_pt3030, Npart_err, ratio_fwd_pt3030_err);
    gr_fwd_pt3030_sys = new TGraphErrors(3,Npart, ratio_fwd_pt3030, Npart_sys, ratio_fwd_pt3030_sys);
    gr_fwd_pt3030_glb = new TGraphErrors(3,Npart, ratio_fwd_pt3030, Npart_sys, ratio_fwd_pt3030_glb);
    // TGraphAsymmErrors *gr_fwd_pt3030_sys = new TGraphAsymmErrors(3,Npart, ratio_fwd_pt3030, Npart_sys, Npart_sys, ratio_fwd_pt3030_sysLo, ratio_fwd_pt3030_sysHi);
    // TGraphAsymmErrors *gr_fwd_pt3030_glb = new TGraphAsymmErrors(3,Npart, ratio_fwd_pt3030, Npart_sys, Npart_sys, ratio_fwd_pt3030_glbLo, ratio_fwd_pt3030_glbHi);
    gr_fwd_pt3030P = new TGraphErrors(3,Npart, ratio_fwd_pt3030, Npart_err, Npart_err);
  }
  gr_fwd_pt3030->SetName("gr_fwd_pt3030");
  gr_fwd_pt3030_sys->SetName("gr_fwd_pt3030_sys");
  gr_fwd_pt3030_glb->SetName("gr_fwd_pt3030_glb");
  gr_fwd_pt3030P->SetName("gr_fwd_pt3030P");				               


  TGraphErrors *gr_fwd_pt3030_split = new TGraphErrors(2,Npart3,ratio_fwd_pt3030_split, Npart_sys, ratio_fwd_pt3030_split_err);
  gr_fwd_pt3030_split->SetName("gr_fwd_pt3030_split");
  gr_fwd_pt3030_split->SetMarkerColor(kRed);
  gr_fwd_pt3030_split->SetMarkerStyle(24);
  gr_fwd_pt3030_split->SetMarkerSize(1.2);

  TGraphErrors *gr_mid_pt6530 = new TGraphErrors(3,Npart, ratio_mid_pt6530, Npart_err, ratio_mid_pt6530_err);
  TGraphErrors *gr_mid_pt6530_sys = new TGraphErrors(3,Npart, ratio_mid_pt6530, Npart_sys, ratio_mid_pt6530_sys);
  TGraphErrors *gr_mid_pt6530_glb = new TGraphErrors(3,Npart, ratio_mid_pt6530, Npart_sys, ratio_mid_pt6530_glb);
  // TGraphAsymmErrors *gr_mid_pt6530_sys = new TGraphAsymmErrors(3,Npart, ratio_mid_pt6530, Npart_sys, Npart_sys, ratio_mid_pt6530_sysLo, ratio_mid_pt6530_sysHi);
  // TGraphAsymmErrors *gr_mid_pt6530_glb = new TGraphAsymmErrors(3,Npart, ratio_mid_pt6530, Npart_sys, Npart_sys, ratio_mid_pt6530_glbLo, ratio_mid_pt6530_glbHi);
  TGraphErrors *gr_mid_pt6530P = new TGraphErrors(3,Npart, ratio_mid_pt6530, Npart_err, Npart_err);
  gr_mid_pt6530->SetName("gr_mid_pt6530");
  gr_mid_pt6530_sys->SetName("gr_mid_pt6530_sys");
  gr_mid_pt6530_sys->SetName("gr_mid_pt6530_glb");
  gr_mid_pt6530P->SetName("gr_mid_pt6530P");

  gr_fwd_pt3030->SetMarkerSize(1.2);
  gr_mid_pt6530->SetMarkerSize(1.2);

  gr_fwd_pt3030P->SetMarkerSize(1.2);
  gr_mid_pt6530P->SetMarkerSize(1.2);

  gr_mid_pt6530->RemovePoint(2);
  gr_mid_pt6530_sys->RemovePoint(2);
  gr_mid_pt6530_glb->RemovePoint(2);
  gr_mid_pt6530P->RemovePoint(2);

  double ratio_mid_pt6530_limit[1] = {0.452};
  double ratio_mid_pt6530_arrow[1] = {0.36};
  double Npart_limit[1] = {32.7694};
  TGraphAsymmErrors *gr_mid_pt6530_limit = new TGraphAsymmErrors(1, Npart_limit, ratio_mid_pt6530_limit, Npart_err, Npart_err, ratio_mid_pt6530_arrow, Npart_err);
  TGraphAsymmErrors *gr_mid_pt6530_bar = new TGraphAsymmErrors(1, Npart_limit, ratio_mid_pt6530_limit, Npart_sys, Npart_sys, Npart_err, Npart_err);
  gr_mid_pt6530_limit->SetName("gr_mid_pt6530_limit");
  gr_mid_pt6530_bar->SetName("gr_mid_pt6530_bar");
  gr_mid_pt6530_limit->SetMarkerStyle(1);
  gr_mid_pt6530_limit->SetMarkerSize(2.0);
  gr_mid_pt6530_bar->SetMarkerSize(0.0);
  gr_mid_pt6530_limit->SetLineWidth(1);
  gr_mid_pt6530_bar->SetLineWidth(2);

  if (isPaper) {
    gr_fwd_pt3030->SetMarkerColor(kRed+2);
    gr_fwd_pt3030->SetLineColor(kRed+2);
    gr_fwd_pt3030_sys->SetLineColor(kRed+2);
    gr_fwd_pt3030_glb->SetLineColor(kRed+2);

    gr_mid_pt6530->SetMarkerColor(kBlue+1);
    gr_mid_pt6530->SetLineColor(kBlue+1);
    gr_mid_pt6530_sys->SetLineColor(kBlue+1);
    gr_mid_pt6530_glb->SetLineColor(kBlue+1);

    gr_mid_pt6530_limit->SetMarkerColor(kBlue+1);
    gr_mid_pt6530_limit->SetLineColor(kBlue+1);
    gr_mid_pt6530_bar->SetMarkerColor(kBlue+1);
    gr_mid_pt6530_bar->SetLineColor(kBlue+1);
  }
  else {
    gr_fwd_pt3030->SetMarkerColor(kRed);
    gr_fwd_pt3030->SetLineColor(kRed);
    gr_fwd_pt3030_sys->SetLineColor(kRed);
    gr_fwd_pt3030_glb->SetLineColor(kRed);

    gr_mid_pt6530->SetMarkerColor(kBlue);
    gr_mid_pt6530->SetLineColor(kBlue);
    gr_mid_pt6530_sys->SetLineColor(kBlue);
    gr_mid_pt6530_glb->SetLineColor(kBlue);

    gr_mid_pt6530_limit->SetMarkerColor(kBlue);
    gr_mid_pt6530_limit->SetLineColor(kBlue);
    gr_mid_pt6530_bar->SetMarkerColor(kBlue);
    gr_mid_pt6530_bar->SetLineColor(kBlue);
  }

  gr_fwd_pt3030_sys->SetFillColor(kRed-9);
  gr_fwd_pt3030_glb->SetFillColor(kRed-9);
  gr_mid_pt6530_sys->SetFillColor(kAzure-9);
  gr_mid_pt6530_glb->SetFillColor(kAzure-9);

  gr_mid_pt6530->SetMarkerStyle(21);

  gr_fwd_pt3030P->SetMarkerStyle(24);
  gr_mid_pt6530P->SetMarkerStyle(25);

  gr_fwd_pt3030_glb->SetFillStyle(0);
  gr_mid_pt6530_glb->SetFillStyle(0);

  TGraphErrors *gr_fwd_pt3030_glb2 = (TGraphErrors*) gr_fwd_pt3030_glb->Clone("gr_fwd_pt3030_glb2");
  TGraphErrors *gr_mid_pt6530_glb2 = (TGraphErrors*) gr_mid_pt6530_glb->Clone("gr_mid_pt6530_glb2");

  gr_fwd_pt3030_glb2->SetFillStyle(1001);
  gr_mid_pt6530_glb2->SetFillStyle(1001);

  gr_fwd_pt3030_glb2->SetFillColor(kGray);
  gr_mid_pt6530_glb2->SetFillColor(kGray);

  double pp_mid_pt6530_err = 0.057;
  double pp_fwd_pt3030_err = 0.066;

  TBox *pp_fwd_pt3030 = new TBox(0.4,1.0-pp_fwd_pt3030_err,10.0,1.0+pp_fwd_pt3030_err);
  pp_fwd_pt3030->SetFillColor(kRed-9);
  TBox *pp_mid_pt6530 = new TBox(0.4,1.0-pp_mid_pt6530_err,10.0,1.0+pp_mid_pt6530_err);
  pp_mid_pt6530->SetFillColor(kAzure-9);

  TBox *pp_mid2_pt6530 = new TBox(10.0,1.0-pp_mid_pt6530_err,20.0,1.0+pp_mid_pt6530_err);
  pp_mid2_pt6530->SetFillColor(kAzure-9);


  // centrality integrated points
  double xval[1] = {0.5*416/zoom};
  TGraphErrors *gr_fwd_pt3030_mb = new TGraphErrors(1,xval,ratio_fwd_pt3030_mb, Npart_err, ratio_fwd_pt3030_mb_err);
  TGraphErrors *gr_fwd_pt3030_mb_sys = new TGraphErrors(1,xval,ratio_fwd_pt3030_mb, Npart_sys, ratio_fwd_pt3030_mb_sys);
  TGraphErrors *gr_fwd_pt3030_mbP = new TGraphErrors(1,xval,ratio_fwd_pt3030_mb, Npart_err, Npart_err);

  gr_fwd_pt3030_mb->SetName("gr_fwd_pt3030_mb");
  gr_fwd_pt3030_mb_sys->SetName("gr_fwd_pt3030_mb_sys");
  gr_fwd_pt3030_mbP->SetName("gr_fwd_pt3030_mbP");				               

  TGraphErrors *gr_mid_pt6530_mb = new TGraphErrors(1,xval,ratio_mid_pt6530_mb, Npart_err, ratio_mid_pt6530_mb_err);
  TGraphErrors *gr_mid_pt6530_mb_sys = new TGraphErrors(1,xval,ratio_mid_pt6530_mb, Npart_sys, ratio_mid_pt6530_mb_sys);
  TGraphErrors *gr_mid_pt6530_mbP = new TGraphErrors(1,xval,ratio_mid_pt6530_mb, Npart_err, Npart_err);

  gr_mid_pt6530_mb->SetName("gr_mid_pt6530_mb");
  gr_mid_pt6530_mb_sys->SetName("gr_mid_pt6530_mb_sys");
  gr_mid_pt6530_mbP->SetName("gr_mid_pt6530_mbP");


  if (isPaper) {
    gr_fwd_pt3030_mb->SetMarkerColor(kRed+2);
    gr_fwd_pt3030_mb->SetLineColor(kRed+2);
    gr_fwd_pt3030_mb_sys->SetLineColor(kRed+2);
    
    gr_mid_pt6530_mb->SetMarkerColor(kBlue+1);
    gr_mid_pt6530_mb->SetLineColor(kBlue+1);
    gr_mid_pt6530_mb_sys->SetLineColor(kBlue+1);
  }
  else {
    gr_fwd_pt3030_mb->SetMarkerColor(kRed);
    gr_fwd_pt3030_mb->SetLineColor(kRed);
    gr_fwd_pt3030_mb_sys->SetLineColor(kRed);
    
    gr_mid_pt6530_mb->SetMarkerColor(kBlue);
    gr_mid_pt6530_mb->SetLineColor(kBlue);
    gr_mid_pt6530_mb_sys->SetLineColor(kBlue);
  }

  gr_fwd_pt3030_mb_sys->SetFillColor(kRed-9);
  gr_mid_pt6530_mb_sys->SetFillColor(kAzure-9);

  gr_fwd_pt3030_mbP->SetMarkerStyle(24);
  gr_mid_pt6530_mbP->SetMarkerStyle(25);
  gr_mid_pt6530_mb->SetMarkerStyle(21);

  gr_fwd_pt3030_mb->SetMarkerSize(1.2);
  gr_mid_pt6530_mb->SetMarkerSize(1.2);
  gr_fwd_pt3030_mbP->SetMarkerSize(1.2);
  gr_mid_pt6530_mbP->SetMarkerSize(1.2);


  TMultiGraph *alice = NULL;
  alice = read_alice();

  TMultiGraph *pas = NULL;
  pas = read_pas();


  TF1 *f3 = new TF1("f3","1",0,416);
  f3->SetLineWidth(1);
  f3->SetNpx(200);
  f3->GetXaxis()->SetTitle("N_{part}");
  f3->GetYaxis()->SetTitle("[ N_{#psi(2S)}/N_{J/#psi} ]_{PbPb} / [ N_{#psi(2S)}/N_{J/#psi} ]_{pp}");
  f3->GetYaxis()->SetTitleOffset(1.15);
  f3->GetYaxis()->SetRangeUser(0.0,3.2);
  f3->GetXaxis()->CenterTitle(kTRUE);

  TF1 *f3a = new TF1("f3a","1",0,416/zoom);
  f3a->SetLineWidth(1);
  f3a->SetNpx(10);
  //  f3a->GetXaxis()->SetTitle("0-100%");
  f3a->GetXaxis()->SetTitle("");
  f3a->GetYaxis()->SetTitle("");
  f3a->GetYaxis()->SetTitleOffset(1.15);
  f3a->GetYaxis()->SetTickLength(zoom*f3->GetYaxis()->GetTickLength());
  //  f3a->GetXaxis()->SetTitleOffset(0.165);
  //  f3a->GetXaxis()->SetTitleSize(zoom*f3->GetXaxis()->GetTitleSize()*0.99);
  //  f3a->GetXaxis()->CenterTitle(kTRUE);
  f3a->GetXaxis()->SetNdivisions(0);
  f3a->GetYaxis()->SetRangeUser(0.0,3.2);

  //  TCanvas *c3 = new TCanvas("c3","c3",600,800);
  TCanvas *c3;
  if (overlay)
    c3 = new TCanvas("c3","c3",650,600);//,1000,600);
  else
    c3 = new TCanvas("c3","c3");//,1000,600);

  if(overlay) {
    c3->Divide(2,1);
    c3->GetPad(1)->SetPad(0.0,0.0,split,1.0);
    c3->GetPad(2)->SetPad(split,0.0,1.0,1.0);
    c3->GetPad(1)->SetRightMargin(0.0);
    c3->GetPad(2)->SetLeftMargin(0.0);
    c3->GetPad(1)->SetLeftMargin(0.13);
  }

  //  c3->Divide(1,2);

  // c3->GetPad(1)->SetPad(0.0,0.5,1.0,1.0);
  // c3->GetPad(2)->SetPad(0.0,0.0,1.0,0.5);
  // c3->GetPad(1)->SetBottomMargin(0.0);
  // c3->GetPad(2)->SetTopMargin(0.0);

  // c3->GetPad(1)->SetPad(0.0,0.0,0.5,1.0);
  // c3->GetPad(2)->SetPad(0.5,0.0,1.0,1.0);
  c3->cd(1);

  if (overlayAlice) {
    f3->GetYaxis()->SetRangeUser(0,6);
    f3a->GetYaxis()->SetRangeUser(0,6);
  }

  if (overlayPAS) {
    f3->GetYaxis()->SetRangeUser(0.0,10);
    f3a->GetYaxis()->SetRangeUser(0.0,10);
  }
  //  c3->SetLogy();
  f3->Draw();

  if (plotGlobalPP) {
    pp_fwd_pt3030->Draw();
    if (overlay)
      pp_mid2_pt6530->Draw();
  }
  else
    gr_fwd_pt3030_glb2->Draw("3");

  if (overlayAlice)
    alice->Draw();

  if (overlayPAS)
    pas->Draw();

  f3->Draw("same");

  gr_fwd_pt3030_sys->Draw("5");
  if (overlay)
    gr_mid_pt6530_sys->Draw("5");
  gr_fwd_pt3030->Draw("P");
  gr_fwd_pt3030P->Draw("PX");
  if (!plotGlobalPP)
    gr_fwd_pt3030_glb->Draw("3");

  if (showSplit)
    gr_fwd_pt3030_split->Draw("P");

  if (overlay) {
    gr_mid_pt6530->Draw("P");
    gr_mid_pt6530P->Draw("PX");
    if (!plotGlobalPP)
      gr_mid_pt6530_glb->Draw("3");

    gr_mid_pt6530_limit->Draw(">");
    gr_mid_pt6530_bar->Draw("Z");
  }

  TLatex *lCMS;
  if (preliminary)
    lCMS = new TLatex(0.16,0.85,"PbPb & pp #sqrt{s_{NN}} = 2.76 TeV");
  else
    lCMS = new TLatex(0.16,0.9,"CMS PbPb & pp #sqrt{s_{NN}} = 2.76 TeV");

  lCMS->SetNDC();
  lCMS->SetTextAlign(11);
  //  lCMS->SetTextSize(0.05);
  lCMS->SetTextSize(0.045);
  TLatex *lPre = new TLatex(0.16,0.90,"CMS Preliminary");
  lPre->SetNDC();
  lPre->SetTextAlign(11);
  lPre->SetTextSize(0.045);
  if(overlay && !preliminary) {
    lCMS->Draw();
    // lPre->Draw();
  }
  TLegend *leg3;
  TLegend *leg3a = NULL;
  if (plotGlobalPP&&!overlay)
    leg3 = new TLegend(0.16,0.77,0.72,0.88);
  else if (overlay) {
    if (overlayAlice) {
      if (preliminary)
	leg3a = new TLegend(0.16,0.54,0.72,0.74);
      else
	leg3a = new TLegend(0.16,0.53,0.72,0.73);
    }
    else if (overlayPAS) {
      if (preliminary)
	leg3a = new TLegend(0.16,0.54,0.72,0.74);
      else
	leg3a = new TLegend(0.16,0.53,0.72,0.73);
    }

    //      leg3 = new TLegend(0.15,0.555,0.72,0.865);
    //    else
    if (preliminary) 
      leg3 = new TLegend(0.16,0.74,0.72,0.94);
    else
      leg3 = new TLegend(0.16,0.735,0.72,0.885);
    //      leg3 = new TLegend(0.16,0.665,0.72,0.865);
  }
  else
    leg3 = new TLegend(0.16,0.78,0.72,0.94);

  //leg1 = new TLegend(0.15,0.68,0.72,0.88);

  leg3->SetFillStyle(0);
  leg3->SetFillColor(0);
  leg3->SetBorderSize(0);
  leg3->SetMargin(0.1);
  //  leg3->SetTextSize(0.04);
  leg3->SetTextSize(0.035);
  if (overlayAlice || overlayPAS) {
    leg3a->SetFillStyle(0);
    leg3a->SetFillColor(0);
    leg3a->SetBorderSize(0);
    leg3a->SetMargin(0.1);
    leg3a->SetTextSize(0.035);
  }
  //  if (!overlay)
  if (preliminary)
    leg3->SetHeader("CMS Preliminary PbPb & pp #sqrt{s_{NN}} = 2.76 TeV");

  leg3->AddEntry(gr_fwd_pt3030,"3 < p_{T} < 30 GeV/c, 1.6 < |y| < 2.4","P");
  if (overlay) {
    if (overlayAlice) {
      leg3->AddEntry(gr_mid_pt6530,"6.5 < p_{T} < 30 GeV/c, |y| < 1.6","P");
      leg3->AddEntry(gr_mid_pt6530_bar,"95% CL","L");
      leg3a->SetHeader("ALICE Preliminary PbPb #sqrt{s_{NN}} = 2.76 TeV & pp #sqrt{s} = 7 TeV");
      leg3a->AddEntry(gr_alice_pt03P,"p_{T} < 3 GeV/c, 2.5 < |y| < 4","P");
      leg3a->AddEntry(gr_alice_pt38,"3 < p_{T} < 8 GeV/c, 2.5 < |y| < 4","P");
      leg3a->AddEntry(gr_alice_pt38_bar,"95% CL","L");//, 3 < p_{T} < 8 GeV/c, 2.5 < |y| < 4
    }
    else if (overlayPAS) {
      leg3->AddEntry(gr_mid_pt6530,"6.5 < p_{T} < 30 GeV/c, |y| < 1.6","P");
      leg3->AddEntry(gr_mid_pt6530_bar,"95% CL","L");
      leg3a->SetHeader("CMS-PAS-HIN-12-007");
      leg3a->AddEntry(gr_pas_fwd_pt3030P,"3 < p_{T} < 30 GeV/c, 1.6 < |y| < 2.4","P");
      leg3a->AddEntry(gr_pas_mid_pt6530,"6.5 < p_{T} < 30 GeV/c, |y| < 1.6","P");
      leg3a->AddEntry(gr_pas_fwd_pt3030_glb,"pp uncertainty (global)","F");
    }
    else {
      leg3->AddEntry(gr_mid_pt6530,"6.5 < p_{T} < 30 GeV/c, |y| < 1.6","P");
      leg3->AddEntry(gr_mid_pt6530_bar,"95% CL","L");//6.5 < p_{T} < 30 GeV/c, |y| < 1.6
    }
  }
  if (!plotGlobalPP)
    leg3->AddEntry(gr_fwd_pt3030_glb2,"pp uncertainty (global)","F");
  leg3->Draw();
  if (overlayAlice||overlayPAS)
    leg3a->Draw();

  //  pre->SetTextSize(0.06);
  //  pre->SetY(9.1667);
  // pre->SetY(0.4);
  // pre->SetX(210);
  //  pre->Draw();

  // TLine *l1 = new TLine(400,1.0,490,7.1);
  // //  l1->Draw("same");
  // TLine *l2 = new TLine(400,0.0,490.0,0.0);
  // //  l2->Draw("same");

  gPad->RedrawAxis();

  TLatex *label = new TLatex(55/zoom,2.675,"#splitline{Cent.}{0-100%}");
  label->SetTextSize(0.23);
  if (overlayAlice)
    label->SetY(5.2);
  label->SetY(8.36);

  if (overlay) {
    c3->cd(2);
    f3a->Draw();
    gr_fwd_pt3030_mb_sys->Draw("5");
    gr_mid_pt6530_mb_sys->Draw("5");
    gr_fwd_pt3030_mb->Draw("P");
    gr_fwd_pt3030_mbP->Draw("PX");
    
    gr_mid_pt6530_mb->Draw("P");
    gr_mid_pt6530_mbP->Draw("PX");
    label->Draw();
  }

  TCanvas *c4=NULL;
  if (!overlay) {
    c4 = new TCanvas("c4","c4");//,1000,600);
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
    if (!plotGlobalPP)
      gr_mid_pt6530_glb->Draw("3");

    gr_mid_pt6530_limit->Draw(">");
    gr_mid_pt6530_bar->Draw("Z");

    TLegend *leg3b;
    if (plotGlobalPP)
      leg3b = new TLegend(0.15,0.74,0.72,0.94);
    else
      leg3b = new TLegend(0.15,0.78,0.72,0.94);

    //leg1 = new TLegend(0.15,0.68,0.72,0.88);

    leg3b->SetFillStyle(0);
    leg3b->SetFillColor(0);
    leg3b->SetBorderSize(0);
    leg3b->SetMargin(0.15);
    //  leg3b->SetTextSize(0.04);
    leg3b->SetTextSize(0.035);
    leg3b->SetHeader("CMS PbPb #sqrt{s_{NN}} = 2.76 TeV");

    leg3b->AddEntry(gr_mid_pt6530,"6.5 < p_{T} < 30 GeV/c, |y| < 1.6","P");
    if (!plotGlobalPP)
      leg3b->AddEntry(gr_mid_pt6530_glb2,"pp uncertainty (global)","F");
    leg3b->Draw();

    // TLatex *pre3 = (TLatex*)pre->Clone();
    // pre3->SetY(0.056);
    //  pre3->Draw();

    // TLine *l1a = new TLine(-90,0.15,0.0,1.0);
    // //  l1a->Draw("same");

    // TLine *l2a = new TLine(-90,0.0,-18.0,0.0);
    // //  l2a->Draw("same");


    gPad->RedrawAxis();
  }
  if (savePlots) {
    if (overlay) {
      if (overlayAlice) {
	c3->SaveAs("double_ratio_compAlice_ppGlobal_pp2013_PbPbRegit.pdf");
	c3->SaveAs("double_ratio_compAlice_ppGlobal_pp2013_PbPbRegit.png");
      }
      else if (overlayPAS) {
	c3->SaveAs("double_ratio_compPAS_ppGlobal_pp2013_PbPbRegit.pdf");
	c3->SaveAs("double_ratio_compPAS_ppGlobal_pp2013_PbPbRegit.png");
      }
      else {
	c3->SaveAs("double_ratio_ppGlobal_pp2013_PbPbRegit.pdf");
	c3->SaveAs("double_ratio_ppGlobal_pp2013_PbPbRegit.png");
      }
    }
    else{
      c3->SaveAs("double_ratio_ppGlobal_fwd_pp2013_PbPbRegit.pdf");
      c3->SaveAs("double_ratio_ppGlobal_fwd_pp2013_PbPbRegit.png");
      c4->SaveAs("double_ratio_ppGlobal_mid_pp2013_PbPbRegit.pdf");
      c4->SaveAs("double_ratio_ppGlobal_mid_pp2013_PbPbRegit.png");
    }
  }


  return;
}

TMultiGraph* read_alice()
{
  TMultiGraph *gr = new TMultiGraph();

  double Npart_alice_err[3] = {0.0, 0.0, 0.0};
  double Npart_alice_sys[3] = {10.0, 10.0, 10.0};

  double Npart_alice_pt03[3] = {16.7901, 69.1358, 157.037};
  double ratio_alice_pt03[3] = {1.9887, 1.06215, 0.94915};
  double ratio_alice_pt03_err[3] = {0.779661, 0.745763, 0.610169};
  double ratio_alice_pt03_sys[3] = {0.282486, 0.700565, 0.293785};
  double ratio_alice_pt03_glb[3] = {0.508475, 0.259887, 0.259887};
  
  double Npart_alice_pt38[2] = {114.568, 312.099};
  double ratio_alice_pt38[2] = {0.632768, 1.74011};
  double ratio_alice_pt38_err[2] = {0.463277, 0.0};
  double ratio_alice_pt38_sys[2] = {0.112994, 0.0};
  double ratio_alice_pt38_glb[2] = {0.169492, 0.0};


  TGraphErrors *gr_alice_pt03 = new TGraphErrors(3, Npart_alice_pt03, ratio_alice_pt03, Npart_alice_err, ratio_alice_pt03_err);
  TGraphErrors *gr_alice_pt03_sys = new TGraphErrors(3, Npart_alice_pt03, ratio_alice_pt03, Npart_alice_sys, ratio_alice_pt03_sys);
  TGraphErrors *gr_alice_pt03_glb = new TGraphErrors(3, Npart_alice_pt03, ratio_alice_pt03, Npart_alice_err, ratio_alice_pt03_glb);
  gr_alice_pt03P = new TGraphErrors(3, Npart_alice_pt03, ratio_alice_pt03, Npart_alice_err, Npart_alice_err);

  gr_alice_pt03->SetName("alice_pt03");
  gr_alice_pt03_sys->SetName("alice_pt03_sys");
  gr_alice_pt03_glb->SetName("alice_pt03_glb");
  gr_alice_pt03P->SetName("alice_pt03P");

  gr_alice_pt03->SetMarkerStyle(33);
  gr_alice_pt03_sys->SetMarkerStyle(27);
  gr_alice_pt03_glb->SetMarkerStyle(27);
  gr_alice_pt03P->SetMarkerStyle(27);
  gr_alice_pt03->SetMarkerSize(2);
  gr_alice_pt03_sys->SetMarkerSize(2);
  gr_alice_pt03_glb->SetMarkerSize(2);
  gr_alice_pt03P->SetMarkerSize(2);
  gr_alice_pt03->SetMarkerColor(kGray);
  gr_alice_pt03_sys->SetMarkerColor(kBlack);
  gr_alice_pt03_glb->SetMarkerColor(kBlack);
  gr_alice_pt03P->SetMarkerColor(kBlack);

  gr_alice_pt03->SetLineColor(kBlack);
  gr_alice_pt03_glb->SetLineColor(kBlack);
  gr_alice_pt03_sys->SetLineColor(kBlack);
  gr_alice_pt03P->SetLineColor(kBlack);

  gr_alice_pt03_glb->SetFillStyle(0);
  gr_alice_pt03_glb->SetFillColor(kBlack);
  gr_alice_pt03_sys->SetFillColor(kGray);

  gr_alice_pt38 = new TGraphErrors(1, Npart_alice_pt38, ratio_alice_pt38, Npart_alice_err, ratio_alice_pt38_err);
  TGraphErrors *gr_alice_pt38_sys = new TGraphErrors(1, Npart_alice_pt38, ratio_alice_pt38, Npart_alice_sys, ratio_alice_pt38_sys);
  TGraphErrors *gr_alice_pt38_glb = new TGraphErrors(1, Npart_alice_pt38, ratio_alice_pt38, Npart_alice_sys, ratio_alice_pt38_glb);

  gr_alice_pt38->SetName("alice_pt38");
  gr_alice_pt38_glb->SetName("alice_pt38_glb");
  gr_alice_pt38_sys->SetName("alice_pt38_sys");
  
  gr_alice_pt38->SetMarkerStyle(33);
  gr_alice_pt38_sys->SetMarkerStyle(33);
  gr_alice_pt38_glb->SetMarkerStyle(33);
  gr_alice_pt38->SetMarkerSize(2);
  gr_alice_pt38_sys->SetMarkerSize(2);
  gr_alice_pt38_glb->SetMarkerSize(2);
  gr_alice_pt38->SetMarkerColor(kBlack);
  gr_alice_pt38_sys->SetMarkerColor(kBlack);
  gr_alice_pt38_glb->SetMarkerColor(kBlack);

  gr_alice_pt38->SetLineColor(kBlack);
  gr_alice_pt38_glb->SetLineColor(kBlack);
  gr_alice_pt38_sys->SetLineColor(kBlack);

  gr_alice_pt38_glb->SetFillColor(kBlack);
  gr_alice_pt38_glb->SetFillStyle(0);
  gr_alice_pt38_sys->SetFillColor(kGray);


  double ratio_alice_pt38_limit[1] = {1.74011};
  double ratio_alice_pt38_arrow[1] = {1.0};
  double Npart_alice_limit[1] = {312.099};
  TGraphAsymmErrors *gr_alice_pt38_limit = new TGraphAsymmErrors(1, Npart_alice_limit, ratio_alice_pt38_limit, Npart_alice_err, Npart_alice_err, ratio_alice_pt38_arrow, Npart_alice_err);
  gr_alice_pt38_bar = new TGraphAsymmErrors(1, Npart_alice_limit, ratio_alice_pt38_limit, Npart_alice_sys, Npart_alice_sys, Npart_alice_err, Npart_alice_err);
  gr_alice_pt38_limit->SetName("gr_alice_pt38_limit");
  gr_alice_pt38_bar->SetName("gr_alice_pt38_bar");
  gr_alice_pt38_limit->SetMarkerStyle(1);
  gr_alice_pt38_limit->SetMarkerSize(2.0);
  gr_alice_pt38_bar->SetMarkerSize(0.0);
  gr_alice_pt38_limit->SetLineWidth(1);
  gr_alice_pt38_bar->SetLineWidth(2);

  gr_alice_pt38_limit->SetLineColor(kBlack);
  gr_alice_pt38_bar->SetLineColor(kBlack);


  gr->Add(gr_alice_pt03_sys,"5");
  gr->Add(gr_alice_pt38_sys,"5");
  gr->Add(gr_alice_pt03_glb,"3");
  gr->Add(gr_alice_pt38_glb,"2");
  gr->Add(gr_alice_pt03,"P");
  gr->Add(gr_alice_pt03P,"P");
  gr->Add(gr_alice_pt38,"P");
  gr->Add(gr_alice_pt38_limit,">");
  gr->Add(gr_alice_pt38_bar,"Z");

  return gr;
}

TMultiGraph* read_pas()
{
  TMultiGraph *gr = new TMultiGraph();

  double Npart_pas[3] = {308.3853, 158.5653, 32.7694};
  double Npart2_pas[3] = {310.3853, 160.5653, 34.7694};

  double Npart_pas_err[3] = {0.0, 0.0, 0.0};
  double Npart_pas_sys[3] = {10.0, 10.0, 10.0};

  double ratio_pas_fwd_pt3030[3] = {5.32, 2.81, 2.54};
  double ratio_pas_fwd_pt3030_err[3] = {1.03, 0.97, 0.87};
  double ratio_pas_fwd_pt3030_sys[3] = {0.79, 0.29, 0.39};
  double ratio_pas_fwd_pt3030_glb[3] = {2.58, 1.36, 1.23};
  
  double ratio_pas_mid_pt6530[3] = {0.49, 0.33, 0.18};
  double ratio_pas_mid_pt6530_err[3] = {0.16, 0.15, 0.15};
  double ratio_pas_mid_pt6530_sys[3] = {0.08, 0.06, 0.07};
  double ratio_pas_mid_pt6530_glb[3] = {0.10, 0.07, 0.04};


  gr_pas_fwd_pt3030 = new TGraphErrors(3, Npart_pas, ratio_pas_fwd_pt3030, Npart_pas_err, ratio_pas_fwd_pt3030_err);
  TGraphErrors *gr_pas_fwd_pt3030_sys = new TGraphErrors(3, Npart_pas, ratio_pas_fwd_pt3030, Npart_pas_sys, ratio_pas_fwd_pt3030_sys);
  gr_pas_fwd_pt3030_glb = new TGraphErrors(3, Npart_pas, ratio_pas_fwd_pt3030, Npart_pas_err, ratio_pas_fwd_pt3030_glb);
  gr_pas_fwd_pt3030P = new TGraphErrors(3, Npart_pas, ratio_pas_fwd_pt3030, Npart_pas_err, Npart_pas_err);

  gr_pas_fwd_pt3030->SetName("pas_fwd_pt3030");
  gr_pas_fwd_pt3030_sys->SetName("pas_fwd_pt3030_sys");
  gr_pas_fwd_pt3030_glb->SetName("pas_fwd_pt3030_glb");
  gr_pas_fwd_pt3030P->SetName("pas_fwd_pt3030P");

  gr_pas_fwd_pt3030->SetMarkerStyle(20);
  gr_pas_fwd_pt3030_sys->SetMarkerStyle(24);
  gr_pas_fwd_pt3030_glb->SetMarkerStyle(24);
  gr_pas_fwd_pt3030P->SetMarkerStyle(24);
  gr_pas_fwd_pt3030->SetMarkerSize(1.2);
  gr_pas_fwd_pt3030_sys->SetMarkerSize(1.2);
  gr_pas_fwd_pt3030_glb->SetMarkerSize(1.2);
  gr_pas_fwd_pt3030P->SetMarkerSize(1.2);
  gr_pas_fwd_pt3030->SetMarkerColor(kGray);
  gr_pas_fwd_pt3030_sys->SetMarkerColor(kBlack);
  gr_pas_fwd_pt3030_glb->SetMarkerColor(kBlack);
  gr_pas_fwd_pt3030P->SetMarkerColor(kBlack);

  gr_pas_fwd_pt3030->SetLineColor(kBlack);
  gr_pas_fwd_pt3030_glb->SetLineColor(kBlack);
  gr_pas_fwd_pt3030_sys->SetLineColor(kBlack);
  gr_pas_fwd_pt3030P->SetLineColor(kBlack);

  gr_pas_fwd_pt3030_glb->SetFillStyle(0);
  gr_pas_fwd_pt3030_glb->SetFillColor(kBlack);
  gr_pas_fwd_pt3030_sys->SetFillColor(kGray);

  gr_pas_mid_pt6530 = new TGraphErrors(3, Npart2_pas, ratio_pas_mid_pt6530, Npart_pas_err, ratio_pas_mid_pt6530_err);
  TGraphErrors *gr_pas_mid_pt6530_sys = new TGraphErrors(3, Npart2_pas, ratio_pas_mid_pt6530, Npart_pas_sys, ratio_pas_mid_pt6530_sys);
  TGraphErrors *gr_pas_mid_pt6530_glb = new TGraphErrors(3, Npart2_pas, ratio_pas_mid_pt6530, Npart_pas_sys, ratio_pas_mid_pt6530_glb);

  gr_pas_mid_pt6530->SetName("pas_mid_pt6530");
  gr_pas_mid_pt6530_glb->SetName("pas_mid_pt6530_glb");
  gr_pas_mid_pt6530_sys->SetName("pas_mid_pt6530_sys");
  
  gr_pas_mid_pt6530->SetMarkerStyle(25);
  gr_pas_mid_pt6530_sys->SetMarkerStyle(25);
  gr_pas_mid_pt6530_glb->SetMarkerStyle(25);
  gr_pas_mid_pt6530->SetMarkerSize(1.2);
  gr_pas_mid_pt6530_sys->SetMarkerSize(1.2);
  gr_pas_mid_pt6530_glb->SetMarkerSize(1.2);
  gr_pas_mid_pt6530->SetMarkerColor(kBlack);
  gr_pas_mid_pt6530_sys->SetMarkerColor(kBlack);
  gr_pas_mid_pt6530_glb->SetMarkerColor(kBlack);

  gr_pas_mid_pt6530->SetLineColor(kBlack);
  gr_pas_mid_pt6530_glb->SetLineColor(kBlack);
  gr_pas_mid_pt6530_sys->SetLineColor(kBlack);

  gr_pas_mid_pt6530_glb->SetFillColor(kBlack);
  gr_pas_mid_pt6530_glb->SetFillStyle(0);
  gr_pas_mid_pt6530_sys->SetFillColor(kGray);


  gr->Add(gr_pas_fwd_pt3030_sys,"5");
  gr->Add(gr_pas_mid_pt6530_sys,"5");
  gr->Add(gr_pas_fwd_pt3030_glb,"3");
  gr->Add(gr_pas_mid_pt6530_glb,"3");
  gr->Add(gr_pas_fwd_pt3030,"P");
  gr->Add(gr_pas_fwd_pt3030P,"P");
  gr->Add(gr_pas_mid_pt6530,"P");

  return gr;
}
