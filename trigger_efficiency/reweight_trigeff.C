#include "TFile.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TGraphAsymmErrors.h"


void reweight_trigeff()
{
  TGraphAsymmErrors *grEffRD;
  TGraphAsymmErrors *grEffMC;

  TFile *tnpRD = new TFile("MC/Jpsi_RD_2012_Trg_L1_Eff_cbExpo.root","READ");
  grEffRD = (TGraphAsymmErrors*) Trg_pt_All->Clone("grEffRD");
  tnpRD->Close();
  grEffRD->SetMarkerColor(kBlue);
  grEffRD->SetMarkerStyle(21);
  TFile *tnpMC = new TFile("MC/Jpsi_MC_2012_Trg_L1_Eff.root","READ");
  grEffMC = (TGraphAsymmErrors*) Trg_pt_All->Clone("grEffMC");
  tnpMC->Close();
  grEffMC->SetMarkerColor(kRed);
  grEffMC->SetMarkerStyle(20);

  // double effMC[8] = {0.48, 0,68, 0.87, 0.94, 0.96, 0.97, 0.96, 0.98}
  // double effRD[8] = {0.59, 0,76, 0.91, 0.97, 0.98, 0.99, 0.96, 0.98}

  TF1 *f1 = new TF1("f1","[0]*TMath::Erf((x-[1])/[2])",1.5,20);

  f1->SetParNames("eff0","x0","m");
  f1->SetParameters(1.0,0.5,1);
  f1->SetLineWidth(2);

  grEffRD->Fit(f1,"RME");
  grEffMC->Fit(f1,"RME");

  grEffRD->GetFunction("f1")->SetLineColor(kBlue);
  grEffMC->GetFunction("f1")->SetLineColor(kRed);

  TCanvas *c1 = new TCanvas("c1","c1");
  grEffRD->Draw("AP");
  grEffMC->Draw("P");

  TH1F * hMuPt_psi2s;
  TH1F * hMuPt_Jpsi;

  TFile *infPsi2s = new TFile("SingleMuPt_fwd_Psi2SMC.root","READ");
  hMuPt_psi2s = (TH1F*) hMuPtRec->Clone("hMuPt_psi2s");
  hMuPt_psi2s->SetDirectory(0);
  infPsi2s->Close();
  TH1F *hMuPtTrgRD_psi2s = (TH1F*)hMuPt_psi2s->Clone("hMuPtTrgRD_psi2s");
  TH1F *hMuPtTrgMC_psi2s = (TH1F*)hMuPt_psi2s->Clone("hMuPtTrgMC_psi2s");

  TFile *infJpsi = new TFile("SingleMuPt_fwd_JpsiMC.root","READ");
  hMuPt_Jpsi = (TH1F*) hMuPtRec->Clone("hMuPt_Jpsi");
  hMuPt_Jpsi->SetDirectory(0);
  infJpsi->Close();

  hMuPt_psi2s->SetMarkerColor(kRed);
  hMuPt_Jpsi->SetMarkerColor(kBlue);
  hMuPt_Jpsi->SetMarkerStyle(21);

  cout << "psi(2S) mean pT: " << hMuPt_psi2s->GetMean() << " RMS: " << hMuPt_psi2s->GetRMS() << endl;
  cout << "J/psi mean pT: " << hMuPt_Jpsi->GetMean() << " RMS: " << hMuPt_Jpsi->GetRMS() << endl;

  TH1F *hMuPtTrgRD_Jpsi = (TH1F*)hMuPt_Jpsi->Clone("hMuPtTrgRD_Jpsi");
  TH1F *hMuPtTrgMC_Jpsi = (TH1F*)hMuPt_Jpsi->Clone("hMuPtTrgMC_Jpsi");

  hMuPtTrgRD_psi2s->Reset();
  hMuPtTrgRD_Jpsi->Reset();
  hMuPtTrgMC_psi2s->Reset();
  hMuPtTrgMC_Jpsi->Reset();

  for (int i=1;i<=hMuPtTrgRD_psi2s->GetNbinsX();++i) {
    double pt = hMuPtTrgRD_psi2s->GetXaxis()->GetBinCenter(i);
    hMuPtTrgRD_psi2s->SetBinContent(i,hMuPt_psi2s->GetBinContent(i)*grEffRD->GetFunction("f1")->Eval(pt));
    hMuPtTrgRD_Jpsi->SetBinContent(i,hMuPt_Jpsi->GetBinContent(i)*grEffRD->GetFunction("f1")->Eval(pt));

    hMuPtTrgMC_psi2s->SetBinContent(i,hMuPt_psi2s->GetBinContent(i)*grEffMC->GetFunction("f1")->Eval(pt));
    hMuPtTrgMC_Jpsi->SetBinContent(i,hMuPt_Jpsi->GetBinContent(i)*grEffMC->GetFunction("f1")->Eval(pt));
  }

  double ptbins[10] = {0.0,1.5,2.5,3.0,4.5,6.0,9.0,12.0,15.0,20.0};
  TH1F *hMuPt_psi2s_rebin = (TH1F*) hMuPt_psi2s->Rebin(9,"hMuPt_psi2s_rebin",ptbins);
  TH1F *hMuPt_Jpsi_rebin = (TH1F*) hMuPt_Jpsi->Rebin(9,"hMuPt_Jpsi_rebin",ptbins);
  TH1F *hMuPtTrgRD_psi2s_rebin = (TH1F*) hMuPtTrgRD_psi2s->Rebin(9,"hMuPtTrgRD_psi2s_rebin",ptbins);
  TH1F *hMuPtTrgRD_Jpsi_rebin = (TH1F*) hMuPtTrgRD_Jpsi->Rebin(9,"hMuPtTrgRD_Jpsi_rebin",ptbins);
  TH1F *hMuPtTrgMC_psi2s_rebin = (TH1F*) hMuPtTrgMC_psi2s->Rebin(9,"hMuPtTrgMC_psi2s_rebin",ptbins);
  TH1F *hMuPtTrgMC_Jpsi_rebin = (TH1F*) hMuPtTrgMC_Jpsi->Rebin(9,"hMuPtTrgMC_Jpsi_rebin",ptbins);

  hMuPt_psi2s_rebin->GetXaxis()->SetRangeUser(0,9);
  hMuPt_psi2s_rebin->GetXaxis()->SetTitle("p_{T}^{#mu} [GeV/c]");
  hMuPt_psi2s_rebin->GetYaxis()->SetTitle("counts [arb. units]");


  TCanvas *c2 = new TCanvas("c2","c2");
  hMuPt_psi2s->DrawNormalized("");
  hMuPt_Jpsi->DrawNormalized("same");

  cout << "Trigger efficiency RD (J/psi): " << hMuPtTrgRD_Jpsi->Integral()/hMuPt_Jpsi->Integral() << endl;
  cout << "Trigger efficiency RD (psi(2S)): " << hMuPtTrgRD_psi2s->Integral()/hMuPt_psi2s->Integral() << endl;
  cout << "Trigger efficiency MC (J/psi): " << hMuPtTrgMC_Jpsi->Integral()/hMuPt_Jpsi->Integral() << endl;
  cout << "Trigger efficiency MC (psi(2S)): " << hMuPtTrgMC_psi2s->Integral()/hMuPt_psi2s->Integral() << endl;


  int max = hMuPt_Jpsi->GetXaxis()->FindBin(4.0-0.0001);
  cout << "In 0.0<pt<4 GeV/c:" << endl;
  cout << "Trigger efficiency RD (J/psi): " << hMuPtTrgRD_Jpsi->Integral(1,max)/hMuPt_Jpsi->Integral(1,max) << endl;
  cout << "Trigger efficiency RD (psi(2S)): " << hMuPtTrgRD_psi2s->Integral(1,max)/hMuPt_psi2s->Integral(1,max) << endl;
  cout << "Trigger efficiency MC (J/psi): " << hMuPtTrgMC_Jpsi->Integral(1,max)/hMuPt_Jpsi->Integral(1,max) << endl;
  cout << "Trigger efficiency MC (psi(2S)): " << hMuPtTrgMC_psi2s->Integral(1,max)/hMuPt_psi2s->Integral(1,max) << endl;


  TLegend *leg1 = new TLegend(0.6,0.8,0.9,0.9);
  leg1->SetFillStyle(0);
  leg1->SetFillColor(0);
  leg1->SetBorderSize(0);
  leg1->SetMargin(0.15);
  leg1->SetTextSize(0.035);
  leg1->AddEntry(hMuPt_psi2s_rebin,"#mu from #psi(2S)","P");
  leg1->AddEntry(hMuPt_Jpsi_rebin,"#mu from J/#psi","P");
  leg1->Draw();

  TGraphAsymmErrors *effPsi2S_RD = new TGraphAsymmErrors();
  TGraphAsymmErrors *effJpsi_RD = new TGraphAsymmErrors();
  TGraphAsymmErrors *effPsi2S_MC = new TGraphAsymmErrors();
  TGraphAsymmErrors *effJpsi_MC = new TGraphAsymmErrors();

  effPsi2S_RD->BayesDivide(hMuPtTrgRD_psi2s_rebin,hMuPt_psi2s_rebin);
  effJpsi_RD->BayesDivide(hMuPtTrgRD_Jpsi_rebin,hMuPt_Jpsi_rebin);

  effPsi2S_MC->BayesDivide(hMuPtTrgMC_psi2s_rebin,hMuPt_psi2s_rebin);
  effJpsi_MC->BayesDivide(hMuPtTrgMC_Jpsi_rebin,hMuPt_Jpsi_rebin);

  effPsi2S_RD->SetMarkerColor(kRed);
  effJpsi_RD->SetMarkerColor(kBlue);
  effPsi2S_MC->SetMarkerColor(kRed);
  effJpsi_MC->SetMarkerColor(kBlue);

  effPsi2S_RD->SetMarkerStyle(20);
  effJpsi_RD->SetMarkerStyle(21);
  effPsi2S_MC->SetMarkerStyle(24);
  effJpsi_MC->SetMarkerStyle(25);

  TCanvas *c3 = new TCanvas("c3","c3",800,600);
  c3->Divide(2,1);
  c3->cd(1);
  effPsi2S_RD->Draw("AP");

  effPsi2S_RD->GetXaxis()->SetRangeUser(0,9);
  effPsi2S_RD->GetYaxis()->SetRangeUser(0.3,1.0);
  effPsi2S_RD->GetXaxis()->SetTitle("p_{T}^{#mu} [GeV/c]");
  effPsi2S_RD->GetYaxis()->SetTitle("Trigger efficiency (Data)");

  effJpsi_RD->Draw("P");

  TLegend *leg2 = new TLegend(0.4,0.2,0.9,0.45);
  leg2->SetFillStyle(0);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetMargin(0.15);
  leg2->SetTextSize(0.035);
  leg2->SetHeader("Tag & Probe efficiency (Data)");
  leg2->AddEntry(effPsi2S_RD,"Reweighted with #mu from #psi(2S)","P");
  leg2->AddEntry(effJpsi_RD,"Reweighted with #mu from J/#psi","P");
  leg2->Draw();

  c3->cd(2);
  effPsi2S_MC->Draw("AP");

  effPsi2S_MC->GetXaxis()->SetRangeUser(0,9);
  effPsi2S_MC->GetYaxis()->SetRangeUser(0.3,1.0);
  effPsi2S_MC->GetXaxis()->SetTitle("p_{T}^{#mu} [GeV/c]");
  effPsi2S_MC->GetYaxis()->SetTitle("Trigger efficiency (MC)");

  effJpsi_MC->Draw("P");

  TLegend *leg2a = new TLegend(0.4,0.2,0.9,0.45);
  leg2a->SetFillStyle(0);
  leg2a->SetFillColor(0);
  leg2a->SetBorderSize(0);
  leg2a->SetMargin(0.15);
  leg2a->SetTextSize(0.035);
  leg2a->SetHeader("Tag & Probe efficiency (MC)");
  leg2a->AddEntry(effPsi2S_MC,"Reweighted with #mu from #psi(2S)","P");
  leg2a->AddEntry(effJpsi_MC,"Reweighted with #mu from J/#psi","P");
  leg2a->Draw();

  return;
}
