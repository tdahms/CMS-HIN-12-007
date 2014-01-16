#include <iostream>
#include <string>

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1.h"


void count_Jpsi_psi2s_yield(std::string fname="histos_pbpb_pt6530_y024.root",bool isHI=true)
{
  gStyle->SetOptStat("emri");
  gStyle->SetOptFit(1);
  TF1 *f0 = new TF1("f0","0",2.6,4.2);
  f0->SetLineWidth(1);
  f0->SetLineStyle(2);

  TFile *inf;
  inf = new TFile(fname.c_str(),"READ");
  TCanvas *c1 = new TCanvas("c1","c1");
  if (isHI)
    c1->Divide(2,2);
  c1->cd(1);
  ((TH1F*) gROOT->FindObject("hmass_os_0100"))->SetMarkerSize(1.0);
  ((TH1F*) gROOT->FindObject("hmass_os_0100"))->SetMarkerStyle(24);
  ((TH1F*) gROOT->FindObject("hmass_os_0100"))->SetMarkerColor(kBlack);
  ((TH1F*) gROOT->FindObject("hmass_os_0100"))->SetLineColor(kBlack);

  ((TH1F*) gROOT->FindObject("hmass_ls_0100"))->SetMarkerSize(1.0);
  ((TH1F*) gROOT->FindObject("hmass_ls_0100"))->SetMarkerStyle(24);
  ((TH1F*) gROOT->FindObject("hmass_ls_0100"))->SetMarkerColor(kRed);
  ((TH1F*) gROOT->FindObject("hmass_ls_0100"))->SetLineColor(kRed);

  ((TH1F*) gROOT->FindObject("hmass_su_0100"))->SetMarkerSize(1.0);
  ((TH1F*) gROOT->FindObject("hmass_su_0100"))->SetMarkerStyle(20);
  ((TH1F*) gROOT->FindObject("hmass_su_0100"))->SetMarkerColor(kBlue);
  ((TH1F*) gROOT->FindObject("hmass_su_0100"))->SetLineColor(kBlue);

  ((TH1F*) gROOT->FindObject("hmass_os_0100"))->Rebin();
  ((TH1F*) gROOT->FindObject("hmass_ls_0100"))->Rebin();
  ((TH1F*) gROOT->FindObject("hmass_su_0100"))->Rebin();

  ((TH1F*) gROOT->FindObject("hmass_os_0100"))->SetMaximum(((TH1F*) gROOT->FindObject("hmass_os_0100"))->GetMaximum()*0.2);
  ((TH1F*) gROOT->FindObject("hmass_os_0100"))->SetMinimum(((TH1F*) gROOT->FindObject("hmass_os_0100"))->GetMaximum()*(-0.2));

  ((TH1F*) gROOT->FindObject("hmass_os_0100"))->Draw();
  ((TH1F*) gROOT->FindObject("hmass_ls_0100"))->Draw("histsame");
  ((TH1F*) gROOT->FindObject("hmass_su_0100"))->Draw("same");
  f0->Draw("same");

  TLatex *l0100;
  if (isHI) 
    l0100 = new TLatex(3.4,((TH1F*) gROOT->FindObject("hmass_os_0100"))->GetMaximum()*0.80,"0-100%");
  else
    l0100 = new TLatex(3.4,((TH1F*) gROOT->FindObject("hmass_os_0100"))->GetMaximum()*0.80,"pp");
  l0100->SetTextSize(0.05);
  l0100->Draw("same");

  if (isHI) {
    c1->cd(2);
    ((TH1F*) gROOT->FindObject("hmass_os_020"))->SetMarkerSize(1.0);
    ((TH1F*) gROOT->FindObject("hmass_os_020"))->SetMarkerStyle(24);
    ((TH1F*) gROOT->FindObject("hmass_os_020"))->SetMarkerColor(kBlack);
    ((TH1F*) gROOT->FindObject("hmass_os_020"))->SetLineColor(kBlack);

    ((TH1F*) gROOT->FindObject("hmass_ls_020"))->SetMarkerSize(1.0);
    ((TH1F*) gROOT->FindObject("hmass_ls_020"))->SetMarkerStyle(24);
    ((TH1F*) gROOT->FindObject("hmass_ls_020"))->SetMarkerColor(kRed);
    ((TH1F*) gROOT->FindObject("hmass_ls_020"))->SetLineColor(kRed);

    ((TH1F*) gROOT->FindObject("hmass_su_020"))->SetMarkerSize(1.0);
    ((TH1F*) gROOT->FindObject("hmass_su_020"))->SetMarkerStyle(20);
    ((TH1F*) gROOT->FindObject("hmass_su_020"))->SetMarkerColor(kBlue);
    ((TH1F*) gROOT->FindObject("hmass_su_020"))->SetLineColor(kBlue);

    ((TH1F*) gROOT->FindObject("hmass_os_020"))->Rebin();
    ((TH1F*) gROOT->FindObject("hmass_ls_020"))->Rebin();
    ((TH1F*) gROOT->FindObject("hmass_su_020"))->Rebin();

    ((TH1F*) gROOT->FindObject("hmass_os_020"))->SetMaximum(((TH1F*) gROOT->FindObject("hmass_os_020"))->GetMaximum()*0.2);
    ((TH1F*) gROOT->FindObject("hmass_os_020"))->SetMinimum(((TH1F*) gROOT->FindObject("hmass_os_020"))->GetMaximum()*(-0.2));
    ((TH1F*) gROOT->FindObject("hmass_os_020"))->Draw();
    ((TH1F*) gROOT->FindObject("hmass_ls_020"))->Draw("histsame");
    ((TH1F*) gROOT->FindObject("hmass_su_020"))->Draw("same");
    f0->Draw("same");

    TLatex *l020 = new TLatex(3.4,((TH1F*) gROOT->FindObject("hmass_os_020"))->GetMaximum()*0.80,"0-20%");
    l020->SetTextSize(0.05);
    l020->Draw("same");

    c1->cd(3);
    ((TH1F*) gROOT->FindObject("hmass_os_2040"))->SetMarkerSize(1.0);
    ((TH1F*) gROOT->FindObject("hmass_os_2040"))->SetMarkerStyle(24);
    ((TH1F*) gROOT->FindObject("hmass_os_2040"))->SetMarkerColor(kBlack);
    ((TH1F*) gROOT->FindObject("hmass_os_2040"))->SetLineColor(kBlack);

    ((TH1F*) gROOT->FindObject("hmass_ls_2040"))->SetMarkerSize(1.0);
    ((TH1F*) gROOT->FindObject("hmass_ls_2040"))->SetMarkerStyle(24);
    ((TH1F*) gROOT->FindObject("hmass_ls_2040"))->SetMarkerColor(kRed);
    ((TH1F*) gROOT->FindObject("hmass_ls_2040"))->SetLineColor(kRed);

    ((TH1F*) gROOT->FindObject("hmass_su_2040"))->SetMarkerSize(1.0);
    ((TH1F*) gROOT->FindObject("hmass_su_2040"))->SetMarkerStyle(20);
    ((TH1F*) gROOT->FindObject("hmass_su_2040"))->SetMarkerColor(kBlue);
    ((TH1F*) gROOT->FindObject("hmass_su_2040"))->SetLineColor(kBlue);

    ((TH1F*) gROOT->FindObject("hmass_os_2040"))->Rebin();
    ((TH1F*) gROOT->FindObject("hmass_ls_2040"))->Rebin();
    ((TH1F*) gROOT->FindObject("hmass_su_2040"))->Rebin();

    ((TH1F*) gROOT->FindObject("hmass_os_2040"))->SetMaximum(((TH1F*) gROOT->FindObject("hmass_os_2040"))->GetMaximum()*0.1);
    ((TH1F*) gROOT->FindObject("hmass_os_2040"))->SetMinimum(((TH1F*) gROOT->FindObject("hmass_os_2040"))->GetMaximum()*(-0.2));
    ((TH1F*) gROOT->FindObject("hmass_os_2040"))->Draw();
    ((TH1F*) gROOT->FindObject("hmass_ls_2040"))->Draw("histsame");
    ((TH1F*) gROOT->FindObject("hmass_su_2040"))->Draw("same");
    f0->Draw("same");

    TLatex *l2040 = new TLatex(3.4,((TH1F*) gROOT->FindObject("hmass_os_2040"))->GetMaximum()*0.80,"20-40%");
    l2040->SetTextSize(0.05);
    l2040->Draw("same");

    c1->cd(4);
    ((TH1F*) gROOT->FindObject("hmass_os_40100"))->SetMarkerSize(1.0);
    ((TH1F*) gROOT->FindObject("hmass_os_40100"))->SetMarkerStyle(24);
    ((TH1F*) gROOT->FindObject("hmass_os_40100"))->SetMarkerColor(kBlack);
    ((TH1F*) gROOT->FindObject("hmass_os_40100"))->SetLineColor(kBlack);

    ((TH1F*) gROOT->FindObject("hmass_ls_40100"))->SetMarkerSize(1.0);
    ((TH1F*) gROOT->FindObject("hmass_ls_40100"))->SetMarkerStyle(24);
    ((TH1F*) gROOT->FindObject("hmass_ls_40100"))->SetMarkerColor(kRed);
    ((TH1F*) gROOT->FindObject("hmass_ls_40100"))->SetLineColor(kRed);

    ((TH1F*) gROOT->FindObject("hmass_su_40100"))->SetMarkerSize(1.0);
    ((TH1F*) gROOT->FindObject("hmass_su_40100"))->SetMarkerStyle(20);
    ((TH1F*) gROOT->FindObject("hmass_su_40100"))->SetMarkerColor(kBlue);
    ((TH1F*) gROOT->FindObject("hmass_su_40100"))->SetLineColor(kBlue);

    ((TH1F*) gROOT->FindObject("hmass_os_40100"))->Rebin();
    ((TH1F*) gROOT->FindObject("hmass_ls_40100"))->Rebin();
    ((TH1F*) gROOT->FindObject("hmass_su_40100"))->Rebin();

    ((TH1F*) gROOT->FindObject("hmass_os_40100"))->SetMaximum(((TH1F*) gROOT->FindObject("hmass_os_40100"))->GetMaximum()*0.1);
    ((TH1F*) gROOT->FindObject("hmass_os_40100"))->SetMinimum(((TH1F*) gROOT->FindObject("hmass_os_40100"))->GetMaximum()*(-0.2));
    ((TH1F*) gROOT->FindObject("hmass_os_40100"))->Draw();
    ((TH1F*) gROOT->FindObject("hmass_ls_40100"))->Draw("histsame");
    ((TH1F*) gROOT->FindObject("hmass_su_40100"))->Draw("same");
    f0->Draw("same");

    TLatex *l40100 = new TLatex(3.4,((TH1F*) gROOT->FindObject("hmass_os_40100"))->GetMaximum()*0.80,"40-100%");
    l40100->SetTextSize(0.05);
    l40100->Draw("same");
  }

  int min1 = hmass_su_0100->GetXaxis()->FindBin(2.95);
  int max1 = hmass_su_0100->GetXaxis()->FindBin(3.25);

  int min2 = hmass_su_0100->GetXaxis()->FindBin(3.48);
  int max2 = hmass_su_0100->GetXaxis()->FindBin(3.88);

  double err1,err2;

  std::cout << "|y|<2.4\t6.5<p_T<30GeV\t0-100%:\t" << hmass_su_0100->IntegralAndError(min1,max1,err1) << " +/- " << err1 << "\t" << hmass_su_0100->IntegralAndError(min2,max2,err2) << " +/- " << err2 << "\t" << hmass_su_0100->Integral(min2,max2)/hmass_su_0100->Integral(min1,max1) << " +/- " << sqrt(pow(err1/hmass_su_0100->Integral(min1,max1),2) + pow(err2/hmass_su_0100->Integral(min2,max2),2))*hmass_su_0100->Integral(min2,max2)/hmass_su_0100->Integral(min1,max1) << std::endl;

  if (isHI) {
    std::cout << "|y|<2.4\t6.5<p_T<30GeV\t0-20%:\t"
	      << hmass_su_020->IntegralAndError(min1,max1,err1) << " +/- " << err1 << "\t"
	      << hmass_su_020->IntegralAndError(min2,max2,err2) << " +/- " << err2 << "\t"
	      << hmass_su_020->Integral(min2,max2)/hmass_su_020->Integral(min1,max1) << " +/- " 
	      << sqrt(pow(err1/hmass_su_020->Integral(min1,max1),2) + pow(err2/hmass_su_020->Integral(min2,max2),2))*hmass_su_020->Integral(min2,max2)/hmass_su_020->Integral(min1,max1)
	      << std::endl;
    std::cout << "|y|<2.4\t6.5<p_T<30GeV\t20-40%:\t"
	      << hmass_su_2040->IntegralAndError(min1,max1,err1) << " +/- " << err1 << "\t"
	      << hmass_su_2040->IntegralAndError(min2,max2,err2) << " +/- " << err2 << "\t"
	      << hmass_su_2040->Integral(min2,max2)/hmass_su_2040->Integral(min1,max1) << " +/- "
	      << sqrt(pow(err1/hmass_su_2040->Integral(min1,max1),2) + pow(err2/hmass_su_2040->Integral(min2,max2),2))*hmass_su_2040->Integral(min2,max2)/hmass_su_2040->Integral(min1,max1)
	      << std::endl;
    std::cout << "|y|<2.4\t6.5<p_T<30GeV\t40-100%:\t"
	      << hmass_su_40100->IntegralAndError(min1,max1,err1) << " +/- " << err1 << "\t"
	      << hmass_su_40100->IntegralAndError(min2,max2,err2) << " +/- " << err2 << "\t"
	      << hmass_su_40100->Integral(min2,max2)/hmass_su_40100->Integral(min1,max1) << " +/- "
	      << sqrt(pow(err1/hmass_su_40100->Integral(min1,max1),2) + pow(err2/hmass_su_40100->Integral(min2,max2),2))*hmass_su_40100->Integral(min2,max2)/hmass_su_40100->Integral(min1,max1)
	      << std::endl;
  }

  return;
}
