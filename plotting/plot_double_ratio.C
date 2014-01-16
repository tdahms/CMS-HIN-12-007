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
  /*
  // 20120418
  // CBG              CB
  fwd_pt3065
  0-100%
  5.473   3.506       5.538   1.908
  0-20%
  7.427   3.761       8.292   3.438
  20-40%
  3.266   10.984      3.294   2.069
  40-100%
  2.342   1.805       2.484   1.628

  fwd_pt3040
  0-100%
  3.602   1.175       3.808   1.427
  0-20%
  4.836   1.681       5.152   1.959
  20-40%
  2.548   3.171       2.650   1.195
  40-100%
  1.510   1.054       1.535   0.841

  fwd_pt6540
  0-100%
  1.498   0.674       1.518   0.662
  0-20%
  1.697   0.805       1.771   0.878
  20-40%
  1.642   0.820       1.679   0.841
  40-100%
  0.420   0.396       0.416   0.490

  all_pt6540
  0-100%
  0.634   0.171       0.618   0.146
  0-20%
  0.769   0.154       0.774   0.209
  20-40%
  0.532   0.178       0.530   0.182
  40-100%
  0.361   0.157       0.367   0.158

  // CB
  double ratio_fwd_pt3065[3] = {8.292, 3.294, 2.069};
  double ratio_fwd_pt3065_err[3] = {3.438, 2.069, 1.628};

  double ratio_fwd_pt3040[3] = {5.152, 2.650, 1.535};
  double ratio_fwd_pt3040_err[3] = {1.959, 1.195, 0.841};

  double ratio_fwd_pt6540[3] = {1.771, 1.679, 0.416};
  double ratio_fwd_pt6540_err[3] = {0.878, 0.841, 0.490};

  double ratio_all_pt6540[3] = {0.774, 0.530, 0.367};
  double ratio_all_pt6540_err[3] = {0.209, 0.182, 0.158};
  
  // CBG
  double ratio_fwd_pt3065[3] = {7.427, 3.266, 2.342};
  double ratio_fwd_pt3065_err[3] = {3.761, 10.984, 1.805};

  double ratio_fwd_pt3040[3] = {4.836, 2.548, 1.510};
  double ratio_fwd_pt3040_err[3] = {1.681, 3.171, 1.054};

  double ratio_fwd_pt6540[3] = {1.697, 1.642, 0.420};
  double ratio_fwd_pt6540_err[3] = {0.805, 0.820, 0.396};

  double ratio_all_pt6540[3] = {0.769, 0.532, 0.361};
  double ratio_all_pt6540_err[3] = {0.154, 0.178, 0.157};
  

  // 20120425
  // CBG              CB
  fwd_pt3065
  0-100%
  5.995   4.367       6.071   4.468
  0-20%
  7.834   4.479       8.501   3.698
  20-40%
  3.763   2.419       3.800   1.932
  40-100%
  3.395   2.100       3.601   2.446

  fwd_pt330
  0-100%
  4.1695  1.695       4.461   1.908
  0-20%
  5.521   1.566       5.769   1.755
  20-40%
  2.841   1.373       3.012   1.405
  40-100%
  2.535   1.102       2.636   1.371

  fwd_pt6530
  0-100%
  1.735   0.259       1.815   0.964
  0-20%
  1.845   0.974       1.971   1.097
  20-40%
  1.725   1.011       1.841   1.092
  40-100%
  1.300   0.883       1.324   0.893

  all_pt6530
  0-100%
  0.612   0.181       0.595   0.162
  0-20%
  0.721   0.279       0.724   0.227
  20-40%
  0.560   0.194       0.563   0.214
  40-100%
  0.372   0.140       0.358   0.188
  */


  // 0-20%, 20-40%, 40-100%
  double Npart[3] = {308.3853, 158.5653, 32.7694};
  double Npart2[3] = {318.3853, 158.5653, 22.7694};
  double Npart_err[3] = {0,0,0};
  double Npart_sys[3] = {10,10,10};
  // CBG
  double ratio_fwd_pt3065[3] = {10.131, 5.652, 5.949};
  double ratio_fwd_pt3065_err[3] = {1.802, 1.726, 1.644};
  //  double ratio_fwd_pt3065_sys[3] = {1.784, 0.773, 0.382};
  double ratio_fwd_pt3065_sys[3] = {2.773, 2.594, 2.915};
  double ratio_fwd_pt3065_glb[3] = {6.698, 3.737, 3.933};
  
  // double ratio_fwd_pt3030[3] = {6.420, 3.222, 2.573};
  // double ratio_fwd_pt3030_err[3] = {1.025, 0.970, 0.858};
  // //  double ratio_fwd_pt3030_sys[3] = {1.239, 0.553, 0.102};
  // double ratio_fwd_pt3030_sys[3] = {1.085, 0.508, 0.407};
  // double ratio_fwd_pt3030_glb[3] = {3.159, 1.585, 1.266};

  // Abdulla 20130814
  double ratio_fwd_pt3030[3] = {2.00, 1.03, 0.79};
  double ratio_fwd_pt3030_err[3] = {0.46, 0.44, 0.38};
  //  double ratio_fwd_pt3030_sys[3] = {1.239, 0.553, 0.102};
  double ratio_fwd_pt3030_sys[3] = {0.188, 0.208, 0.128};
  double ratio_fwd_pt3030_glb[3] = {0.126, 0.064, 0.050};
  
  double ratio_fwd_pt6530[3] = {1.803, 1.773, 1.241 };
  double ratio_fwd_pt6530_err[3] = {0.598, 0.649, 0.611};
  //  double ratio_fwd_pt6530_sys[3] = {0.290, 0.220, 0.208};
  double ratio_fwd_pt6530_sys[3] = {0.479, 0.313, 0.302};
  double ratio_fwd_pt6530_glb[3] = {0.957, 0.941, 0.659};

  // CBG
  // double ratio_mid_pt6530[3] = {0.531, 0.374, 0.208};
  // double ratio_mid_pt6530_err[3] = {0.158, 0.154, 0.146};
  // //  double ratio_mid_pt6530_sys[3] = {0.010, 0.026, 0.228};
  // double ratio_mid_pt6530_sys[3] = {0.078, 0.061, 0.192};
  // double ratio_mid_pt6530_glb[3] = {0.110, 0.078, 0.043};
  // CB

  // Abdulla 20130814
  double ratio_mid_pt6530[3] = {0.81, 0.67, 0.29};
  double ratio_mid_pt6530_err[3] = {0.17, 0.17, 0.17};
  //  double ratio_mid_pt6530_sys[3] = {0.010, 0.026, 0.228};
  double ratio_mid_pt6530_sys[3] = {0.106, 0.024, 0.108};
  double ratio_mid_pt6530_glb[3] = {0.042, 0.035, 0.015};
  
  // double ratio_all_pt6530[3] = {0.674, 0.534, 0.305};
  // double ratio_all_pt6530_err[3] = {0.173, 0.172, 0.166};
  // //  double ratio_all_pt6530_sys[3] = {0.031, 0.024, 0.251};
  // double ratio_all_pt6530_sys[3] = {0.138, 0.029, 0.087};
  // double ratio_all_pt6530_glb[3] = {0.131, 0.104, 0.059};
  double ratio_all_pt6530[3] = {0.751, 0.525, 0.344};
  double ratio_all_pt6530_err[3] = {0.177, 0.172, 0.162};
  //  double ratio_all_pt6530_sys[3] = {0.031, 0.024, 0.251};
  double ratio_all_pt6530_sys[3] = {0.068, 0.044, 0.050};
  double ratio_all_pt6530_glb[3] = {0.145, 0.101, 0.066};


  
  /*
  // CBG
  double ratio_fwd_pt3065[3] = {7.427, 3.266, 2.342};
  double ratio_fwd_pt3065_err[3] = {3.761, 10.984, 1.805};

  double ratio_fwd_pt3040[3] = {4.836, 2.548, 1.510};
  double ratio_fwd_pt3040_err[3] = {1.681, 3.171, 1.054};

  double ratio_fwd_pt6540[3] = {1.697, 1.642, 0.420};
  double ratio_fwd_pt6540_err[3] = {0.805, 0.820, 0.396};

  double ratio_all_pt6540[3] = {0.769, 0.532, 0.361};
  double ratio_all_pt6540_err[3] = {0.154, 0.178, 0.157};
  */
  
  TGraphErrors *gr_fwd_pt3065 = new TGraphErrors(3,Npart, ratio_fwd_pt3065, Npart_err, ratio_fwd_pt3065_err);
  TGraphErrors *gr_fwd_pt3065_sys = new TGraphErrors(3,Npart, ratio_fwd_pt3065, Npart_sys, ratio_fwd_pt3065_sys);
  TGraphErrors *gr_fwd_pt3065_glb = new TGraphErrors(3,Npart2, ratio_fwd_pt3065, Npart_sys, ratio_fwd_pt3065_glb);
  TGraphErrors *gr_fwd_pt3065P = new TGraphErrors(3,Npart, ratio_fwd_pt3065, Npart_err, Npart_err);
  gr_fwd_pt3065->SetName("gr_fwd_pt3065");
  gr_fwd_pt3065_sys->SetName("gr_fwd_pt3065_sys");
  gr_fwd_pt3065_glb->SetName("gr_fwd_pt3065_glb");
  gr_fwd_pt3065P->SetName("gr_fwd_pt3065P");
									               
  TGraphErrors *gr_fwd_pt3030 = new TGraphErrors(3,Npart, ratio_fwd_pt3030, Npart_err, ratio_fwd_pt3030_err);
  TGraphErrors *gr_fwd_pt3030_sys = new TGraphErrors(3,Npart, ratio_fwd_pt3030, Npart_sys, ratio_fwd_pt3030_sys);
  TGraphErrors *gr_fwd_pt3030_glb = new TGraphErrors(3,Npart2, ratio_fwd_pt3030, Npart_sys, ratio_fwd_pt3030_glb);
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
  TGraphErrors *gr_mid_pt6530_sys = new TGraphErrors(3,Npart, ratio_mid_pt6530, Npart_sys, ratio_mid_pt6530_sys);
  TGraphErrors *gr_mid_pt6530_glb = new TGraphErrors(3,Npart2, ratio_mid_pt6530, Npart_sys, ratio_mid_pt6530_glb);
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
  f3->GetYaxis()->SetTitle("[#psi(2S) #/J/#psi]_{PbPb} #/[#psi(2S) #/J/#psi]_{pp}");
  f3->GetYaxis()->SetTitleOffset(1.1);
  f3->GetYaxis()->SetRangeUser(0.0,10.0);
  f3->GetXaxis()->CenterTitle(kTRUE);

  TF1 *f3a = new TF1("f3a","1",0,400);
  f3a->SetLineWidth(1);
  f3a->GetXaxis()->SetTitle("N_{part}");
  //  f3a->GetYaxis()->SetTitle("[#psi(2S) #/J/#psi]_{PbPb} #/[#psi(2S) #/J/#psi]_{pp}");
  f3a->GetYaxis()->SetTitleOffset(1.1);
  f3a->GetYaxis()->SetRangeUser(0.0,1.4);
  f3a->GetXaxis()->CenterTitle(kTRUE);

  //  TCanvas *c3 = new TCanvas("c3","c3",600,800);
  TCanvas *c3 = new TCanvas("c3","c3",1000,600);
  c3->Divide(1,2);

  // c3->GetPad(1)->SetPad(0.0,0.5,1.0,1.0);
  // c3->GetPad(2)->SetPad(0.0,0.0,1.0,0.5);
  // c3->GetPad(1)->SetBottomMargin(0.0);
  // c3->GetPad(2)->SetTopMargin(0.0);

  c3->GetPad(1)->SetPad(0.0,0.0,0.5,1.0);
  c3->GetPad(2)->SetPad(0.5,0.0,1.0,1.0);
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
  l1->Draw("same");
  TLine *l2 = new TLine(400,0.0,490.0,0.0);
  l2->Draw("same");

  gPad->RedrawAxis();


  c3->cd(2);
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
  l1a->Draw("same");

  TLine *l2a = new TLine(-90,0.0,-18.0,0.0);
  l2a->Draw("same");


  gPad->RedrawAxis();

  if (savePlots) {
    c3->SaveAs("double_ratio_ppGlobal_midpt_pp2013_PbPbRegit.pdf");
    c3->SaveAs("double_ratio_ppGlobal_midpt_pp2013_PbPbRegit.png");
  }


  return;
  /*
  double ptbins[] = {4.75, 10};
  double ptbins_err[] = {1.25, 3.5};
  double ptSystX[] = {0.75, 0.75};

  // CBG
  double ratio_fwd[3] = {7.488, 1.806};
  double ratio_fwd_err[3] = {1.118, 0.424};
  double ratio_fwd_sys[3] = {1.311, 0.175};


  // CBG
  // double ratio_fwd[3] = {5.473, 1.498};
  // double ratio_fwd_err[3] = {3.506, 0.674};

  TGraphErrors *gr_fwd_pt = new TGraphErrors(2, ptbins, ratio_fwd, ptbins_err, ratio_fwd_err);
  TGraphErrors *gr_fwd_pt_sys = new TGraphErrors(2, ptbins, ratio_fwd, ptSystX, ratio_fwd_sys);
  gr_fwd_pt->SetName("gr_fwd_pt");
  gr_fwd_pt_sys->SetName("gr_fwd_pt_sys");
  gr_fwd_pt_sys->SetFillColor(kGray-2);

  TF1 *f2 = new TF1("f2","1",0,20);
  f2->SetLineWidth(1);
  f2->GetXaxis()->SetTitle("p_{T}");
  f2->GetYaxis()->SetTitle("R_{AA}(#psi(2S))/R_{AA}(J/#psi)");
  f2->GetYaxis()->SetRangeUser(0.0,10.0);
  f2->GetXaxis()->CenterTitle(kTRUE);

  TCanvas *c2 = new TCanvas("c2","c2");
  f2->Draw();
  gr_fwd_pt_sys->Draw("2");
  gr_fwd_pt->Draw("P");

  return;

  */
}
