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
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TPaveStats.h"
#include "TLatex.h"
#include "TLegend.h"


void plot_double_ratio(bool isPaper=false, bool plotGlobalPP = false, bool savePlots=false)
{
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);

  // 0-20%, 20-40%, 40-100%
  double Npart[3] = {308.3853, 158.5653, 32.7694};
  double Npart2[3] = {/*318.3853, */308.3853, 158.5653, 32.7694, /*22.7694*/};
  double Npart_err[3] = {0,0,0};
  double Npart_sys[3] = {10,10,10};

  double ratio_fwd_pt3030[3];
  double ratio_fwd_pt3030_err[3];
  double ratio_fwd_pt3030_sys[3];
  double ratio_fwd_pt3030_glb[3];

  std::ifstream ifs("fwd.txt");
  int i=0;
  if (ifs.is_open()) {
    while (ifs.good()) {
      ifs >> ratio_fwd_pt3030[i] >> ratio_fwd_pt3030_err[i] >> ratio_fwd_pt3030_sys[i] >> ratio_fwd_pt3030_glb[i];
      i++;
    }
    ifs.close();
  }

  double ratio_mid_pt6530[3] =     {};
  double ratio_mid_pt6530_err[3] = {};
  double ratio_mid_pt6530_sys[3] = {};
  double ratio_mid_pt6530_glb[3] = {};

  ifs.open("mid.txt");
  i=0;
  if (ifs.is_open()) {
    while (ifs.good()) {
      ifs >> ratio_mid_pt6530[i] >> ratio_mid_pt6530_err[i] >> ratio_mid_pt6530_sys[i] >> ratio_mid_pt6530_glb[i];
      i++;
    }
    ifs.close();
  }
									               
  TGraphErrors *gr_fwd_pt3030 = new TGraphErrors(3,Npart, ratio_fwd_pt3030, Npart_err, ratio_fwd_pt3030_err);
  TGraphErrors *gr_fwd_pt3030_sys = new TGraphErrors(3,Npart, ratio_fwd_pt3030, Npart_sys, ratio_fwd_pt3030_sys);
  TGraphErrors *gr_fwd_pt3030_glb = new TGraphErrors(3,Npart2, ratio_fwd_pt3030, Npart_sys, ratio_fwd_pt3030_glb);
  // TGraphAsymmErrors *gr_fwd_pt3030_sys = new TGraphAsymmErrors(3,Npart, ratio_fwd_pt3030, Npart_sys, Npart_sys, ratio_fwd_pt3030_sysLo, ratio_fwd_pt3030_sysHi);
  // TGraphAsymmErrors *gr_fwd_pt3030_glb = new TGraphAsymmErrors(3,Npart2, ratio_fwd_pt3030, Npart_sys, Npart_sys, ratio_fwd_pt3030_glbLo, ratio_fwd_pt3030_glbHi);
  TGraphErrors *gr_fwd_pt3030P = new TGraphErrors(3,Npart, ratio_fwd_pt3030, Npart_err, Npart_err);
  gr_fwd_pt3030->SetName("gr_fwd_pt3030");
  gr_fwd_pt3030_sys->SetName("gr_fwd_pt3030_sys");
  gr_fwd_pt3030_glb->SetName("gr_fwd_pt3030_glb");
  gr_fwd_pt3030P->SetName("gr_fwd_pt3030P");				               
									               
  TGraphErrors *gr_mid_pt6530 = new TGraphErrors(3,Npart, ratio_mid_pt6530, Npart_err, ratio_mid_pt6530_err);
  TGraphErrors *gr_mid_pt6530_sys = new TGraphErrors(3,Npart, ratio_mid_pt6530, Npart_sys, ratio_mid_pt6530_sys);
  TGraphErrors *gr_mid_pt6530_glb = new TGraphErrors(3,Npart2, ratio_mid_pt6530, Npart_sys, ratio_mid_pt6530_glb);
  // TGraphAsymmErrors *gr_mid_pt6530_sys = new TGraphAsymmErrors(3,Npart, ratio_mid_pt6530, Npart_sys, Npart_sys, ratio_mid_pt6530_sysLo, ratio_mid_pt6530_sysHi);
  // TGraphAsymmErrors *gr_mid_pt6530_glb = new TGraphAsymmErrors(3,Npart2, ratio_mid_pt6530, Npart_sys, Npart_sys, ratio_mid_pt6530_glbLo, ratio_mid_pt6530_glbHi);
  TGraphErrors *gr_mid_pt6530P = new TGraphErrors(3,Npart, ratio_mid_pt6530, Npart_err, Npart_err);
  gr_mid_pt6530->SetName("gr_mid_pt6530");
  gr_mid_pt6530_sys->SetName("gr_mid_pt6530_sys");
  gr_mid_pt6530_sys->SetName("gr_mid_pt6530_glb");
  gr_mid_pt6530P->SetName("gr_mid_pt6530P");

  gr_fwd_pt3030->SetMarkerSize(1.2);
  gr_mid_pt6530->SetMarkerSize(2.0);

  gr_fwd_pt3030P->SetMarkerSize(1.2);
  gr_mid_pt6530P->SetMarkerSize(2.0);

  gr_mid_pt6530->RemovePoint(2);
  gr_mid_pt6530_sys->RemovePoint(2);
  gr_mid_pt6530_glb->RemovePoint(2);
  gr_mid_pt6530P->RemovePoint(2);

  double ratio_mid_pt6530_limit[1] = {0.46};
  double ratio_mid_pt6530_arrow[1] = {0.36};
  double Npart_limit[1] = {32.7694};
  TGraphAsymmErrors *gr_mid_pt6530_limit = new TGraphAsymmErrors(1, Npart_limit, ratio_mid_pt6530_limit, Npart_err, Npart_err, ratio_mid_pt6530_arrow, Npart_err);
  TGraphAsymmErrors *gr_mid_pt6530_bar = new TGraphAsymmErrors(1, Npart_limit, ratio_mid_pt6530_limit, Npart_sys, Npart_sys, ratio_mid_pt6530_arrow, Npart_err);
  gr_mid_pt6530_limit->SetName("gr_mid_pt6530_limit");
  gr_mid_pt6530_bar->SetName("gr_mid_pt6530_bar");
  gr_mid_pt6530_limit->SetMarkerSize(2.0);
  gr_mid_pt6530_bar->SetMarkerSize(0.0);
  gr_mid_pt6530_limit->SetMarkerColor(kGreen+2);
  gr_mid_pt6530_limit->SetLineColor(kGreen+2);
  gr_mid_pt6530_bar->SetMarkerColor(kGreen+2);
  gr_mid_pt6530_bar->SetLineColor(kGreen+2);

  if (isPaper) {
    gr_fwd_pt3030->SetMarkerColor(kRed+2);
    gr_fwd_pt3030->SetLineColor(kRed+2);
    gr_fwd_pt3030_sys->SetLineColor(kRed+2);
    gr_fwd_pt3030_glb->SetFillColor(kRed+2);
  }
  else {
    gr_fwd_pt3030->SetMarkerColor(kRed);
    gr_fwd_pt3030->SetLineColor(kRed);
    gr_fwd_pt3030_sys->SetLineColor(kRed);
    gr_fwd_pt3030_glb->SetFillColor(kRed);
  }

  gr_fwd_pt3030_sys->SetFillColor(kRed-9);

  gr_mid_pt6530->SetMarkerColor(kGreen+2);
  gr_mid_pt6530->SetLineColor(kGreen+2);
  gr_mid_pt6530_sys->SetLineColor(kGreen+2);
  gr_mid_pt6530_sys->SetFillColor(kGreen-9);
  gr_mid_pt6530_glb->SetFillColor(kGreen+2);
  gr_mid_pt6530->SetMarkerStyle(33);

  gr_fwd_pt3030P->SetMarkerStyle(24);
  gr_mid_pt6530P->SetMarkerStyle(27);

  gr_fwd_pt3030_glb->SetLineColor(kRed+2);
  gr_mid_pt6530_glb->SetLineColor(kGreen+2);

  gr_fwd_pt3030_glb->SetFillStyle(0);
  gr_mid_pt6530_glb->SetFillStyle(0);

  TGraphErrors *gr_fwd_pt3030_glb2 = (TGraphErrors*) gr_fwd_pt3030_glb->Clone("gr_fwd_pt3030_glb2");
  TGraphErrors *gr_mid_pt6530_glb2 = (TGraphErrors*) gr_mid_pt6530_glb->Clone("gr_mid_pt6530_glb2");

  gr_fwd_pt3030_glb2->SetFillStyle(1001);
  gr_mid_pt6530_glb2->SetFillStyle(1001);

  gr_fwd_pt3030_glb2->SetFillColor(kGray);
  gr_mid_pt6530_glb2->SetFillColor(kGray);

  double pp_mid_pt6530_err = 0.056;
  double pp_fwd_pt3030_err = 0.064;

  TBox *pp_fwd_pt3030 = new TBox(0.0,1.0-pp_fwd_pt3030_err,10.0,1.0+pp_fwd_pt3030_err);
  pp_fwd_pt3030->SetFillColor(kRed-9);
  TBox *pp_mid_pt6530 = new TBox(0.0,1.0-pp_mid_pt6530_err,10.0,1.0+pp_mid_pt6530_err);
  pp_mid_pt6530->SetFillColor(kGreen-9);

  TBox *pp_mid2_pt6530 = new TBox(7.0,1.0-pp_mid_pt6530_err,14.0,1.0+pp_mid_pt6530_err);
  pp_mid2_pt6530->SetFillColor(kGreen-9);

  TF1 *f3 = new TF1("f3","1",0,400);
  f3->SetLineWidth(1);
  f3->GetXaxis()->SetTitle("N_{part}");
  f3->GetYaxis()->SetTitle("[ #psi (2S) #/J/#psi ]_{PbPb} #/[ #psi (2S) #/J/#psi ]_{pp}");
  f3->GetYaxis()->SetTitleOffset(1.15);
  f3->GetYaxis()->SetRangeUser(0.0,3.0);
  f3->GetXaxis()->CenterTitle(kTRUE);

  TF1 *f3a = new TF1("f3a","1",0,400);
  f3a->SetLineWidth(1);
  f3a->GetXaxis()->SetTitle("N_{part}");
  f3a->GetYaxis()->SetTitle("[ #psi (2S) #/J/#psi ]_{PbPb} #/[ #psi (2S) #/J/#psi ]_{pp}");
  f3a->GetYaxis()->SetTitleOffset(1.15);
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
  if (!plotGlobalPP)
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
  leg3->SetHeader("CMS PbPb #sqrt{s_{NN}} = 2.76 TeV");
  leg3->AddEntry(gr_fwd_pt3030,"3 < p_{T} < 30 GeV/c, 1.6 < |y| < 2.4","P");
  if (!plotGlobalPP)
    leg3->AddEntry(gr_fwd_pt3030_glb2,"pp uncertainty (global)","F");
  leg3->Draw();
  
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
  if (!plotGlobalPP)
    gr_mid_pt6530_glb->Draw("3");

  gr_mid_pt6530_limit->Draw(">");
  gr_mid_pt6530_bar->Draw("Z");

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
  leg3a->SetHeader("CMS PbPb #sqrt{s_{NN}} = 2.76 TeV");

  leg3a->AddEntry(gr_mid_pt6530,"6.5 < p_{T} < 30 GeV/c, |y| < 1.6","P");
  if (!plotGlobalPP)
    leg3a->AddEntry(gr_mid_pt6530_glb2,"pp uncertainty (global)","F");
  leg3a->Draw();

  // TLatex *pre3 = (TLatex*)pre->Clone();
  // pre3->SetY(0.056);
  //  pre3->Draw();

  // TLine *l1a = new TLine(-90,0.15,0.0,1.0);
  // //  l1a->Draw("same");

  // TLine *l2a = new TLine(-90,0.0,-18.0,0.0);
  // //  l2a->Draw("same");


  gPad->RedrawAxis();

  if (savePlots) {
    c3->SaveAs("double_ratio_ppGlobal_fwd_pp2013_PbPbRegit.pdf");
    c3->SaveAs("double_ratio_ppGlobal_fwd_pp2013_PbPbRegit.png");
    c4->SaveAs("double_ratio_ppGlobal_mid_pp2013_PbPbRegit.pdf");
    c4->SaveAs("double_ratio_ppGlobal_mid_pp2013_PbPbRegit.png");
  }


  return;
}
