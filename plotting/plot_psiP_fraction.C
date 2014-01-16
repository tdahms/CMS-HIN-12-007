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


void plot_psiP_fraction(bool isPaper=false, bool savePlots=false)
{
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(2);

  // 0-20%, 20-40%, 40-100%
  double Npart[4] = {308.3853, 158.5653, 32.7694, 2};
  double Npart_err[4] = {0,0,0,0};
  double Npart_sys[4] = {10,10,10,5};
  /* 20120421
  // CB
     double ratio_fwd_pt3065[4] = {0.149, 0.064, 0.049, 0.019};
     double ratio_fwd_pt3065_err[4] = {0.027, 0.025, 0.022, 0.010};

     double ratio_fwd_pt3040[4] = {0.115, 0.059, 0.034, 0.022};
     double ratio_fwd_pt3040_err[4] = {0.019, 0.017, 0.014, 0.008};

     double ratio_fwd_pt6540[4] = {0.054, 0.050, 0.013, 0.030};
     double ratio_fwd_pt6540_err[4] = {0.017, 0.017, 0.014, 0.012};

     double ratio_all_pt6540[4] = {0.035, 0.024, 0.017, 0.046};
     double ratio_all_pt6540_err[4] = {0.007, 0.007, 0.007, 0.008};
  */

  // 20120426
  // CBG
 double ratio_fwd_pt3065[4] = {0.1408, 0.0756, 0.0687, 0.0166};
 double ratio_fwd_pt3065_err[4] = {0.0288, 0.0279, 0.0254, 0.0095};
 // double ratio_fwd_pt3065_sys[4] = {0.0289, 0.0128, 0.0061, 0.0003};
 double ratio_fwd_pt3065_sys[4] = {0.0289, 0.0220, 0.0165, 0.0058};
 
 double ratio_fwd_pt3030[4] = {0.1253, 0.0629, 0.0502, 0.0195};
 double ratio_fwd_pt3030_err[4] = {0.0200, 0.0189, 0.0168, 0.0075};
 // double ratio_fwd_pt3030_sys[4] = {0.0229, 0.0107, 0.0016, 0.0023};
 double ratio_fwd_pt3030_sys[4] = {0.0169, 0.0083, 0.0072, 0.0060};
 
 double ratio_fwd_pt6530[4] = {0.0537, 0.0460, 0.0330, 0.0258};
 double ratio_fwd_pt6530_err[4] = {0.017, 0.018, 0.017, 0.011};
 // double ratio_fwd_pt6530_sys[4] = {0.0068, 0.0057, 0.0053, 0.0077};
 double ratio_fwd_pt6530_sys[4] = {0.0104, 0.0093, 0.0065, 0.0062};
 
 double ratio_mid_pt6530[4] = {0.0243, 0.0177, 0.0084, 0.0481};
 double ratio_mid_pt6530_err[4] = {0.0078, 0.0181, 0.0071, 0.0100};
 // double ratio_mid_pt6530_sys[4] = {0.0003, 0.0012, 0.0110, 0.0010};
 double ratio_mid_pt6530_sys[4] = {0.0038, 0.0019, 0.0052, 0.0015};
 
 double ratio_all_pt6530[4] = {0.0279, 0.0221, 0.0127, 0.0415};
 double ratio_all_pt6530_err[4] = {0.0072, 0.0071, 0.0069, 0.0078};
 // double ratio_all_pt6530_sys[4] = {0.0012, 0.0010, 0.0104, 0.0021};
 double ratio_all_pt6530_sys[4] = {0.0057, 0.0011, 0.0036, 0.0020};
 
 

  // CBG
  // double ratio_fwd_pt3065[4] = {0.142, 0.075, 0.069, 0.017};
  // double ratio_fwd_pt3065_err[4] = {0.030, 0.039, 0.026, 0.009};

  // double ratio_fwd_pt3040[4] = {0.110, 0.064, 0.051, 0.020};
  // double ratio_fwd_pt3040_err[4] = {0.020, 0.019, 0.017, 0.007};

  // double ratio_fwd_pt6540[4] = {0.053, 0.046, 0.033, 0.026};
  // double ratio_fwd_pt6540_err[4] = {0.017, 0.027, 0.015, 0.011};

  // double ratio_all_pt6540[4] = {0.030, 0.023, 0.012, 0.042};
  // double ratio_all_pt6540_err[4] = {0.007, 0.007, 0.007, 0.008};
  

  TGraphErrors *gr_fwd_pt3065 = new TGraphErrors(4,Npart, ratio_fwd_pt3065, Npart_err, ratio_fwd_pt3065_err);
  TGraphErrors *gr_fwd_pt3065_sys = new TGraphErrors(4,Npart, ratio_fwd_pt3065, Npart_sys, ratio_fwd_pt3065_sys);
  TGraphErrors *gr_fwd_pt3065P = new TGraphErrors(4,Npart, ratio_fwd_pt3065, Npart_err, Npart_err);
  gr_fwd_pt3065->SetName("gr_fwd_pt3065");				               
  gr_fwd_pt3065_sys->SetName("gr_fwd_pt3065_sys");				               
  gr_fwd_pt3065P->SetName("gr_fwd_pt3065P");				               
									               
  TGraphErrors *gr_fwd_pt3030 = new TGraphErrors(4,Npart, ratio_fwd_pt3030, Npart_err, ratio_fwd_pt3030_err);
  TGraphErrors *gr_fwd_pt3030_sys = new TGraphErrors(4,Npart, ratio_fwd_pt3030, Npart_sys, ratio_fwd_pt3030_sys);
  TGraphErrors *gr_fwd_pt3030P = new TGraphErrors(4,Npart, ratio_fwd_pt3030, Npart_err, Npart_err);
  gr_fwd_pt3030->SetName("gr_fwd_pt3030");				               
  gr_fwd_pt3030_sys->SetName("gr_fwd_pt3030_sys");				               
  gr_fwd_pt3030P->SetName("gr_fwd_pt3030P");				               
									               
  TGraphErrors *gr_fwd_pt6530 = new TGraphErrors(4,Npart, ratio_fwd_pt6530, Npart_err, ratio_fwd_pt6530_err);
  TGraphErrors *gr_fwd_pt6530_sys = new TGraphErrors(4,Npart, ratio_fwd_pt6530, Npart_sys, ratio_fwd_pt6530_sys);
  TGraphErrors *gr_fwd_pt6530P = new TGraphErrors(4,Npart, ratio_fwd_pt6530, Npart_err, Npart_err);
  gr_fwd_pt6530->SetName("gr_fwd_pt6530");				               
  gr_fwd_pt6530_sys->SetName("gr_fwd_pt6530_sys");				               	
  gr_fwd_pt6530P->SetName("gr_fwd_pt6530P");				               
								               
  TGraphErrors *gr_mid_pt6530 = new TGraphErrors(4,Npart, ratio_mid_pt6530, Npart_err, ratio_mid_pt6530_err);
  TGraphErrors *gr_mid_pt6530_sys = new TGraphErrors(4,Npart, ratio_mid_pt6530, Npart_sys, ratio_mid_pt6530_sys);
  TGraphErrors *gr_mid_pt6530P = new TGraphErrors(4,Npart, ratio_mid_pt6530, Npart_err, Npart_err);
  gr_mid_pt6530->SetName("gr_mid_pt6530");
  gr_mid_pt6530_sys->SetName("gr_mid_pt6530_sys");
  gr_mid_pt6530P->SetName("gr_mid_pt6530P");

  TGraphErrors *gr_all_pt6530 = new TGraphErrors(4,Npart, ratio_all_pt6530, Npart_err, ratio_all_pt6530_err);
  TGraphErrors *gr_all_pt6530_sys = new TGraphErrors(4,Npart, ratio_all_pt6530, Npart_sys, ratio_all_pt6530_sys);
  TGraphErrors *gr_all_pt6530P = new TGraphErrors(4,Npart, ratio_all_pt6530, Npart_err, Npart_err);
  gr_all_pt6530->SetName("gr_all_pt6530");
  gr_all_pt6530_sys->SetName("gr_all_pt6530_sys");
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
    gr_fwd_pt6530->SetMarkerColor(kBlue+2);
    gr_fwd_pt6530->SetLineColor(kBlue+2);
  }
  else {
    gr_fwd_pt3065->SetMarkerColor(kRed);
    gr_fwd_pt3065->SetLineColor(kRed);
    
    gr_fwd_pt6530->SetMarkerColor(kBlue);
    gr_fwd_pt6530->SetLineColor(kBlue);
  }

  gr_mid_pt6530->SetMarkerColor(kGreen+2);
  gr_mid_pt6530->SetLineColor(kGreen+2);

  gr_all_pt6530->SetMarkerColor(kGreen+2);
  gr_all_pt6530->SetLineColor(kGreen+2);

  gr_fwd_pt3065_sys->SetFillColor(kRed-9);
  gr_fwd_pt6530_sys->SetFillColor(kBlue-9);
  gr_mid_pt6530_sys->SetFillColor(kGreen-9);
  gr_all_pt6530_sys->SetFillColor(kGreen-9);


  //  gr_all_pt6530->SetMarkerColor(kRed);
  gr_fwd_pt6530->SetMarkerStyle(21);
  gr_mid_pt6530->SetMarkerStyle(33);
  gr_all_pt6530->SetMarkerStyle(33);

  gr_fwd_pt3065P->SetMarkerStyle(24);
  gr_fwd_pt6530P->SetMarkerStyle(25);
  gr_mid_pt6530P->SetMarkerStyle(27);
  gr_all_pt6530P->SetMarkerStyle(27);


  TF1 *f1 = new TF1("f1","1",-5,400);
  f1->SetLineWidth(1);
  f1->GetXaxis()->SetTitle("N_{part}");
  f1->GetYaxis()->SetTitle("R_{#psi(2S)}");
  f1->GetYaxis()->SetRangeUser(0.0,0.20);
  f1->GetXaxis()->CenterTitle(kTRUE);

  TCanvas *c1 = new TCanvas("c1","c1");
  f1->Draw();
  gr_fwd_pt3065_sys->Draw("2");
  gr_fwd_pt6530_sys->Draw("2");
  gr_mid_pt6530_sys->Draw("2");
  gr_fwd_pt3065->Draw("P");
  gr_fwd_pt6530->Draw("P");
  gr_mid_pt6530->Draw("P");
  gr_fwd_pt3065P->Draw("PX");
  gr_fwd_pt6530P->Draw("PX");
  gr_mid_pt6530P->Draw("PX");


  TLegend *leg1 = new TLegend(0.15,0.68,0.72,0.88);
  leg1->SetFillStyle(0);
  leg1->SetFillColor(0);
  leg1->SetBorderSize(0);
  leg1->SetMargin(0.15);
  leg1->SetTextSize(0.035);
  leg1->SetHeader("PbPb  #sqrt{s_{NN}} = 2.76 TeV");
  leg1->AddEntry(gr_fwd_pt3065,"3 < p_{T} < 6.5 GeV/c, 1.6 < |y| < 2.4","P");
  leg1->AddEntry(gr_fwd_pt6530,"6.5 < p_{T} < 30 GeV/c, 1.6 < |y| < 2.4","P");
  leg1->AddEntry(gr_mid_pt6530,"6.5 <  p_{T} < 30 GeV/c, |y|<1.6","P");
  leg1->Draw();

  TLatex *pre = new TLatex(18,0.183,"CMS Preliminary");
  pre->SetTextFont(42);
  pre->SetTextSize(0.05);
  pre->Draw();

  gPad->RedrawAxis();

  if (savePlots) {
    c1->SaveAs("single_ratio.pdf");
    c1->SaveAs("single_ratio.png");
  }


  // calculate double ratio
  double double_ratio_fwd_pt3065[4];
  double double_ratio_fwd_pt3030[4];
  double double_ratio_fwd_pt6530[4];
  double double_ratio_mid_pt6530[4];
  double double_ratio_all_pt6530[4];

  double double_ratio_fwd_pt3065_err[4];
  double double_ratio_fwd_pt3030_err[4];
  double double_ratio_fwd_pt6530_err[4];
  double double_ratio_mid_pt6530_err[4];
  double double_ratio_all_pt6530_err[4];

  for (int i=0;i<4;++i) {
    double_ratio_fwd_pt3065[i] = ratio_fwd_pt3065[i]/ratio_fwd_pt3065[3];
    double_ratio_fwd_pt3030[i] = ratio_fwd_pt3030[i]/ratio_fwd_pt3030[3];
    double_ratio_fwd_pt6530[i] = ratio_fwd_pt6530[i]/ratio_fwd_pt6530[3];
    double_ratio_mid_pt6530[i] = ratio_mid_pt6530[i]/ratio_mid_pt6530[3];
    double_ratio_all_pt6530[i] = ratio_all_pt6530[i]/ratio_all_pt6530[3];

    double_ratio_fwd_pt3065_err[i] = double_ratio_fwd_pt3065[i] * ratio_fwd_pt3065_err[i]/ratio_fwd_pt3065[i];
    double_ratio_fwd_pt3030_err[i] = double_ratio_fwd_pt3030[i] * ratio_fwd_pt3030_err[i]/ratio_fwd_pt3030[i];
    double_ratio_fwd_pt6530_err[i] = double_ratio_fwd_pt6530[i] * ratio_fwd_pt6530_err[i]/ratio_fwd_pt6530[i];
    double_ratio_mid_pt6530_err[i] = double_ratio_mid_pt6530[i] * ratio_mid_pt6530_err[i]/ratio_mid_pt6530[i];
    double_ratio_all_pt6530_err[i] = double_ratio_all_pt6530[i] * ratio_all_pt6530_err[i]/ratio_all_pt6530[i];
  }

  TGraphErrors *gr_dbl_fwd_pt3065 = new TGraphErrors(4,Npart, double_ratio_fwd_pt3065, Npart_err, double_ratio_fwd_pt3065_err);
  gr_dbl_fwd_pt3065->SetName("gr_dbl_fwd_pt3065");				               
									               
  TGraphErrors *gr_dbl_fwd_pt3030 = new TGraphErrors(4,Npart, double_ratio_fwd_pt3030, Npart_err, double_ratio_fwd_pt3030_err);
  gr_dbl_fwd_pt3030->SetName("gr_dbl_fwd_pt3030");				               
									               
  TGraphErrors *gr_dbl_fwd_pt6530 = new TGraphErrors(4,Npart, double_ratio_fwd_pt6530, Npart_err, double_ratio_fwd_pt6530_err);
  gr_dbl_fwd_pt6530->SetName("gr_dbl_fwd_pt6530");				               
									               
  TGraphErrors *gr_dbl_mid_pt6530 = new TGraphErrors(4,Npart, double_ratio_mid_pt6530, Npart_err, double_ratio_mid_pt6530_err);
  gr_dbl_mid_pt6530->SetName("gr_dbl_mid_pt6530");

  TGraphErrors *gr_dbl_all_pt6530 = new TGraphErrors(4,Npart, double_ratio_all_pt6530, Npart_err, double_ratio_all_pt6530_err);
  gr_dbl_all_pt6530->SetName("gr_dbl_all_pt6530");


  gr_dbl_fwd_pt3065->RemovePoint(3);
  gr_dbl_fwd_pt3030->RemovePoint(3);
  gr_dbl_fwd_pt6530->RemovePoint(3);
  gr_dbl_mid_pt6530->RemovePoint(3);
  gr_dbl_all_pt6530->RemovePoint(3);

  gr_dbl_fwd_pt3065->SetMarkerSize(1.2);
  gr_dbl_fwd_pt3030->SetMarkerSize(1.2);
  gr_dbl_fwd_pt6530->SetMarkerSize(1.2);
  gr_dbl_mid_pt6530->SetMarkerSize(2.0);
  gr_dbl_all_pt6530->SetMarkerSize(2.0);


  if (isPaper) {
    gr_dbl_fwd_pt3065->SetMarkerColor(kRed+2);
    gr_dbl_fwd_pt6530->SetMarkerColor(kBlue+2);
  }
  else {
    gr_dbl_fwd_pt3065->SetMarkerColor(kRed);
    //  gr_dbl_fwd_pt3065->SetLineColor(kRed);
    
    gr_dbl_fwd_pt6530->SetMarkerColor(kBlue);
    //  gr_dbl_fwd_pt6530->SetLineColor(kBlue);
  }

  gr_dbl_mid_pt6530->SetMarkerColor(kGreen+2);
  //  gr_dbl_mid_pt6530->SetLineColor(kGreen+2);
  gr_dbl_all_pt6530->SetMarkerColor(kGreen+2);
  //  gr_dbl_all_pt6530->SetLineColor(kGreen+2);
  //  gr_dbl_all_pt6530->SetMarkerColor(kRed);

  gr_dbl_fwd_pt6530->SetMarkerStyle(21);
  gr_dbl_mid_pt6530->SetMarkerStyle(33);
  gr_dbl_all_pt6530->SetMarkerStyle(33);



  TBox *pp_fwd_pt3065 = new TBox(0.0,1.0-double_ratio_fwd_pt3065_err[3],10.0,1.0+double_ratio_fwd_pt3065_err[3]);
  pp_fwd_pt3065->SetFillColor(kRed-9);
  TBox *pp_fwd_pt3030 = new TBox(0.0,1.0-double_ratio_fwd_pt3030_err[3],7.0,1.0+double_ratio_fwd_pt3030_err[3]);
  TBox *pp_fwd_pt6530 = new TBox(7.0,1.0-double_ratio_fwd_pt6530_err[3],14.0,1.0+double_ratio_fwd_pt6530_err[3]);
  pp_fwd_pt6530->SetFillColor(kBlue-9);
  TBox *pp_mid_pt6530 = new TBox(14.0,1.0-double_ratio_all_pt6530_err[3],21.0,1.0+double_ratio_all_pt6530_err[3]);
  pp_mid_pt6530->SetFillColor(kGreen-9);
  TBox *pp_all_pt6530 = new TBox(20.0,1.0-double_ratio_all_pt6530_err[3],30.0,1.0+double_ratio_all_pt6530_err[3]);
  pp_all_pt6530->SetFillColor(kGreen-9);



  TF1 *f2 = new TF1("f2","1",0,400);
  f2->SetLineWidth(1);
  f2->GetXaxis()->SetTitle("N_{part}");
  f2->GetYaxis()->SetTitle("R_{#psi(2S)}(PbPb)/R_{#psi(2S)}(pp)");
  f2->GetYaxis()->SetTitleOffset(1.1);
  f2->GetYaxis()->SetRangeUser(0.0,12.0);
  f2->GetXaxis()->CenterTitle(kTRUE);


  TCanvas *c2 = new TCanvas("c2","c2");
  f2->Draw();
  pp_fwd_pt3065->Draw();
  //  pp_fwd_pt3030->Draw();
  pp_fwd_pt6530->Draw();
  pp_mid_pt6530->Draw();

  gr_dbl_fwd_pt3065->Draw("P");
  //  gr_dbl_fwd_pt3030->Draw("P");
  gr_dbl_fwd_pt6530->Draw("P");
  gr_dbl_mid_pt6530->Draw("P");
  f2->Draw("same");

  TLegend *leg2 = new TLegend(0.15,0.68,0.72,0.88);
  leg2->SetFillStyle(0);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetMargin(0.15);
  leg2->SetTextSize(0.035);
  leg2->SetHeader("PbPb  #sqrt{s_{NN}} = 2.76 TeV");
  leg2->AddEntry(gr_dbl_fwd_pt3065,"3 < p_{T} < 6.5 GeV/c, 1.6 < |y| < 2.4","P");
  leg2->AddEntry(gr_dbl_fwd_pt6530,"6.5 < p_{T} < 30 GeV/c, 1.6 < |y| < 2.4","P");
  leg2->AddEntry(gr_dbl_mid_pt6530,"6.5 <  p_{T} < 30 GeV/c, |y|<1.6","P");
  leg2->Draw();

  TLatex *pre2 = new TLatex(18,11,"CMS Preliminary");
  pre2->SetTextFont(42);
  pre2->SetTextSize(0.05);
  pre2->Draw();

  gPad->RedrawAxis();


  return;
}
