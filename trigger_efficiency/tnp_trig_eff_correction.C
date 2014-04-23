#include <iostream>
#include <string>

#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"

using namespace std;

void tnp_trig_eff_correction(string dataFile="../root_files/Jpsi_Regit_RD_2013_Trg3_Eff.root", string mcFile="../root_files/Jpsi_Regit_MC_2013_Trg3_Eff.root", bool savePlots=false)
{
  gStyle->SetOptFit(0);

  bool isTrg=false;
  bool isTrk=false;
  bool isMuId=false;

  bool isForward=false;
  if (dataFile.find("Trg3")!=string::npos || dataFile.find("MuId3")!=string::npos)
    isForward=true;

  bool isPbPb=true;
  if (dataFile.find("pp")!=string::npos)
    isPbPb=false;

  if (dataFile.find("Trg")!=string::npos)
    isTrg=true;
  else if (dataFile.find("Trk")!=string::npos)
    isTrk=true;
  else if (dataFile.find("MuId")!=string::npos)
    isMuId=true;
  else {
    cout << "No good TnP file provided. Exiting" << endl;
    return;
  }

  TGraphAsymmErrors *grEffRD=NULL;
  TGraphAsymmErrors *grEffMC=NULL;

  TFile *tnpRD = new TFile(dataFile.c_str(),"READ");
  if (isTrg)
    grEffRD = (TGraphAsymmErrors*) tnpRD->Get("Trg_pt_All")->Clone("grEffRD");
  else if (isTrk)
    grEffRD = (TGraphAsymmErrors*) tnpRD->Get("Trk_pt_All")->Clone("grEffRD");
  else if (isMuId)
    grEffRD = (TGraphAsymmErrors*) tnpRD->Get("MuId_pt_All")->Clone("grEffRD");

  //  tnpRD->Close();
  grEffRD->SetMarkerColor(kBlue);
  grEffRD->SetMarkerStyle(21);
  TFile *tnpMC = new TFile(mcFile.c_str(),"READ");
  if (isTrg)
    grEffMC = (TGraphAsymmErrors*) tnpMC->Get("Trg_pt_All")->Clone("grEffMC");
  else if (isTrk)
    grEffMC = (TGraphAsymmErrors*) tnpMC->Get("Trk_pt_All")->Clone("grEffMC");
  else if (isMuId)
    grEffMC = (TGraphAsymmErrors*) tnpMC->Get("MuId_pt_All")->Clone("grEffMC");

  //  tnpMC->Close();
  grEffMC->SetMarkerColor(kRed);
  grEffMC->SetMarkerStyle(20);

  TF1 *f1 = new TF1("f1","[0]*TMath::Erf((x-[1])/[2])",0.0,20);

  f1->SetParNames("eff0","x0","m");
  f1->SetParameters(0.8,0.5,2.5);
  f1->SetLineWidth(2);

  grEffRD->Fit(f1,"RME");
  f1->SetParameters(0.9,0.5,2.5);
  grEffMC->Fit(f1,"WRME");

  grEffRD->GetFunction("f1")->SetLineColor(kBlue);
  grEffMC->GetFunction("f1")->SetLineColor(kRed);

  TCanvas *c1 = new TCanvas("c1","c1");
  grEffRD->Draw("AP");
  grEffMC->Draw("P");

  grEffRD->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  if (isTrg)
    grEffRD->GetYaxis()->SetTitle("Trigger Efficiency");
  else if (isTrk)
    grEffRD->GetYaxis()->SetTitle("Tracking Efficiency");
  else if (isMuId)
    grEffRD->GetYaxis()->SetTitle("Muon ID Efficiency");




  TGraphAsymmErrors *grEffRatio = (TGraphAsymmErrors*) grEffMC->Clone("grEffRatio");
  grEffRatio->GetFunction("f1")->Delete();

  double xMC = 0;
  double yMC = 0;
  double xRD = 0;
  double yRD = 0;

  for (int i=0; i<grEffRatio->GetN(); ++i) {
    grEffMC->GetPoint(i,xMC,yMC);
    grEffRD->GetPoint(i,xRD,yRD);
    grEffRatio->SetPoint(i,xMC,yRD/yMC);

    double errxl = grEffMC->GetErrorXlow(i);
    double errxh = grEffMC->GetErrorXhigh(i);

    // double erryl = yRD/yMC*sqrt( pow(grEffMC->GetErrorYlow(i)/yMC,2) + pow(grEffRD->GetErrorYlow(i)/yRD,2) );
    // double erryh = yRD/yMC*sqrt( pow(grEffMC->GetErrorYhigh(i)/yMC,2) + pow(grEffRD->GetErrorYhigh(i)/yRD,2) );
    // let's ignore the MC uncertainties
    double erryl = grEffRD->GetErrorYlow(i)/yMC;
    double erryh = grEffRD->GetErrorYhigh(i)/yMC;

    grEffRatio->SetPointError(i, errxl, errxh, erryl, erryh);

    cout << "pt bin: " << xMC-errxl  << "--" << xMC+errxh << " weight: " << yRD/yMC << endl;

  }

  TF1 *f2 = new TF1("f2","[0]*TMath::Erf((x-[1])/[2])/TMath::Erf((x-[3])/[4])",1.5,20);
  f2->SetParameter(0, grEffRD->GetFunction("f1")->GetParameter(0)/grEffMC->GetFunction("f1")->GetParameter(0));
  f2->SetParameter(1, grEffRD->GetFunction("f1")->GetParameter(1));
  f2->SetParameter(2, grEffRD->GetFunction("f1")->GetParameter(2));
  f2->SetParameter(3, grEffMC->GetFunction("f1")->GetParameter(1));
  f2->SetParameter(4, grEffMC->GetFunction("f1")->GetParameter(2));

  grEffRatio->Fit(f2,"RME"); // fit to ratio

  TCanvas *c2 = new TCanvas("c2","c2");
  grEffRatio->Draw("AP");
  grEffRatio->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  grEffRatio->GetYaxis()->SetTitle("Data (fit) / MC (fit)");
  grEffRatio->GetYaxis()->SetRangeUser(0.8,2.0);
  f2->SetLineStyle(2);
  f2->SetParameter(0, grEffRD->GetFunction("f1")->GetParameter(0)/grEffMC->GetFunction("f1")->GetParameter(0));
  f2->SetParameter(1, grEffRD->GetFunction("f1")->GetParameter(1));
  f2->SetParameter(2, grEffRD->GetFunction("f1")->GetParameter(2));
  f2->SetParameter(3, grEffMC->GetFunction("f1")->GetParameter(1));
  f2->SetParameter(4, grEffMC->GetFunction("f1")->GetParameter(2));
  f2->Draw("same"); // ratio of fits

  TLegend *leg = new TLegend(0.45,0.75,0.9,0.90);
  leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetMargin(0.15); leg->SetTextSize(0.035);
  leg->AddEntry(grEffRatio,"(Data points)/(MC points)","P");
  leg->AddEntry(f2,"(Data Fit)/(MC Fit)","L");
  leg->AddEntry(grEffRatio->GetFunction("f2"),"(Data/MC) Fit","L");
  leg->Draw();

  if (savePlots) {
    if (isTrg) {
      if (isForward) {
	if (isPbPb) {
	  c1->SaveAs("20140324/tnp_pbpb_trg_eff_fwd.pdf");
	  c2->SaveAs("20140324/tnp_pbpb_trg_eff_fwd_dataMC_Ratio.pdf");
	}
	else {
	  c1->SaveAs("20140324/tnp_pp_trg_eff_fwd.pdf");
	  c2->SaveAs("20140324/tnp_pp_trg_eff_fwd_dataMC_Ratio.pdf");
	}
      }
      else { 
	if (isPbPb) {
	  c1->SaveAs("20140324/tnp_pbpb_trg_eff_mid.pdf");
	  c2->SaveAs("20140324/tnp_pbpb_trg_eff_mid_dataMC_Ratio.pdf");
	}
	else {
	  c1->SaveAs("20140324/tnp_pp_trg_eff_mid.pdf");
	  c2->SaveAs("20140324/tnp_pp_trg_eff_mid_dataMC_Ratio.pdf");
	}
      }
    }
    else if (isTrk) {
	c1->SaveAs("20140324/tnp_trk_eff.pdf");
	c2->SaveAs("20140324/tnp_trk_eff_dataMC_Ratio.pdf");
    }
    else if (isMuId) {
      if (isForward) {
	c1->SaveAs("20140324/tnp_muId_eff_fwd.pdf");
	c2->SaveAs("20140324/tnp_muId_eff_fwd_dataMC_Ratio.pdf");
      }
      else {
	c1->SaveAs("20140324/tnp_muId_eff_mid.pdf");
	c2->SaveAs("20140324/tnp_muId_eff_mid_dataMC_Ratio.pdf");
      }
    }
  }

  return;
}
