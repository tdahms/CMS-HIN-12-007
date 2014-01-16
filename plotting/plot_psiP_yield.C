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


void plot_psiP_yield()
{
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);

  // 0-20%, 20-40%, 40-100%
  double Npart[4] = {308.3853, 158.5653, 32.7694, 2};
  double Npart_err[4] = {0,0,0,0};

  double TAA[4] = {18.8317, 6.9361, 0.8482, 1};

  // CB
  double yield_fwd_pt3065[4] = {409.5, 91.7, 42.5, 13.2};
  double yield_fwd_pt3065_err[4] = {82.9, 38.1, 18.9, 6.9};

  double yield_fwd_pt3040[4] = {492.9, 133.5, 44.8, 22.9};
  double yield_fwd_pt3040_err[4] = {86.0, 40.9, 98.8, 8.0};

  double yield_fwd_pt6540[4] = {84.0, 42.6, 6.2, 10.5};
  double yield_fwd_pt6540_err[4] = {71.0, 15.0, 6.6, 4.1};

  double yield_all_pt6540[4] = {181.8, 67.8, 26.0, 49.9};
  double yield_all_pt6540_err[4] = {39.2, 20.6, 10.9, 8.5};
  
  for (int i=0; i<4; ++i) {
    yield_fwd_pt3065[i]/=TAA[i];
    yield_fwd_pt3065_err[i]/=TAA[i];
    yield_fwd_pt3040[i]/=TAA[i];
    yield_fwd_pt3040_err[i]/=TAA[i];
    yield_fwd_pt6540[i]/=TAA[i];
    yield_fwd_pt6540_err[i]/=TAA[i];
    yield_all_pt6540[i]/=TAA[i];
    yield_all_pt6540_err[i]/=TAA[i];
  }


  TGraphErrors *gr_fwd_pt3065 = new TGraphErrors(4,Npart, yield_fwd_pt3065, Npart_err, yield_fwd_pt3065_err);
  gr_fwd_pt3065->SetName("gr_fwd_pt3065");				               
									               
  TGraphErrors *gr_fwd_pt3040 = new TGraphErrors(4,Npart, yield_fwd_pt3040, Npart_err, yield_fwd_pt3040_err);
  gr_fwd_pt3040->SetName("gr_fwd_pt3040");				               
									               
  TGraphErrors *gr_fwd_pt6540 = new TGraphErrors(4,Npart, yield_fwd_pt6540, Npart_err, yield_fwd_pt6540_err);
  gr_fwd_pt6540->SetName("gr_fwd_pt6540");				               
									               
  TGraphErrors *gr_all_pt6540 = new TGraphErrors(4,Npart, yield_all_pt6540, Npart_err, yield_all_pt6540_err);
  gr_all_pt6540->SetName("gr_all_pt6540");


  gr_fwd_pt3065->SetMarkerSize(1.0);
  gr_fwd_pt3040->SetMarkerSize(1.0);
  gr_fwd_pt6540->SetMarkerSize(1.0);
  gr_all_pt6540->SetMarkerSize(1.0);


  gr_fwd_pt3065->SetMarkerColor(kBlue);
  gr_fwd_pt3065->SetLineColor(kBlue);

  gr_fwd_pt6540->SetMarkerColor(kRed);
  gr_fwd_pt6540->SetLineColor(kRed);

  //  gr_all_pt6540->SetMarkerColor(kRed);
  gr_all_pt6540->SetMarkerStyle(24);


  TF1 *f1 = new TF1("f1","0",0,400);
  f1->SetLineWidth(1);
  f1->GetXaxis()->SetTitle("N_{part}");
  f1->GetYaxis()->SetTitle("N_{#psi(2S)}/T_{AA}");
  f1->GetYaxis()->SetRangeUser(0.0,100);
  f1->GetXaxis()->CenterTitle(kTRUE);

  TCanvas *c1 = new TCanvas("c1","c1");
  f1->Draw();
  gr_fwd_pt3065->Draw("P");
  //  gr_fwd_pt3040->Draw("P");
  gr_fwd_pt6540->Draw("P");
  gr_all_pt6540->Draw("P");

  return;
}
