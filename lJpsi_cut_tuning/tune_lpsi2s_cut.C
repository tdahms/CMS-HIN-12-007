#include <iostream>
#include "TCanvas.h"
#include "TCut.h"
#include "TLegend.h"
#include "TChain.h"
#include "TH1.h"
#include "TLine.h"
#include "TLatex.h"
#include "TRandom.h"

void tune_lpsi2s_cut(float ptmin=0.0, float ptmax=30.0, float ymin=0.0, float ymax=2.4, bool absRapidity=true, double threshold=0.8, double efficiency=0.9, bool savePlot=false)
{
  float lCut=-10.0;
  //  float efficiency = 1.0;
  TLegend *leg = new TLegend(0.42,0.4,0.95,0.6);
  leg->SetHeader("PYTHIA: pp #sqrt{s} = 2.76 TeV");
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetMargin(0.15);
  leg->SetTextSize(0.035);

  //  TCut defaultCut = "Reco_QQ_sign==0&&Gen_QQ_4mom.M()>3.68&&Gen_QQ_4mom.M()<3.70&&Reco_QQ_4mom.M()>3.3&&Reco_QQ_4mom.M()<4.0&&Reco_QQ_ctauTrue>-10";
  TCut defaultCut = "Reco_QQ_sign==0&&Reco_QQ_4mom.M()>3.35&&Reco_QQ_4mom.M()<4.0&&Reco_QQ_ctauTrue>-10";
  TCut ptCut = Form("Reco_QQ_4mom.Pt()>%4.1f&&Reco_QQ_4mom.Pt()<%4.1f",ptmin,ptmax);
  TCut rapCut;
  if (absRapidity)
    rapCut = Form("abs(Reco_QQ_4mom.Rapidity())>%3.1f&&abs(Reco_QQ_4mom.Rapidity())<%3.1f",ymin,ymax);
  else
    rapCut = Form("Reco_QQ_4mom.Rapidity()>%3.1f&&Reco_QQ_4mom.Rapidity()<%3.1f",ymin,ymax);

  unsigned int trigBit=2; // DoubleMu0_HighQ
  TCut trigCut = Form("(HLTriggers&%u)==%u&&(Reco_QQ_trig&%u)=%u",trigBit,trigBit,trigBit,trigBit);

  TString fname;
  fname = Form("20140115/psi2s_pp_eff_%3.1f_Rap_%3.1f-%3.1f_Pt_%3.1f-%3.1f.pdf",efficiency,ymin,ymax,ptmin,ptmax);

  std::cout << fname << std::endl;


  std::cout << "default: " << defaultCut.GetTitle() << std::endl;
  std::cout << "pt cut: " << ptCut.GetTitle() << std::endl;
  std::cout << "rapidity cut: " << rapCut.GetTitle() << std::endl;
  std::cout << "trigger bit: " << trigCut.GetTitle() << std::endl;

  TLatex *lpt;
  if (ptmin==0.0)
    lpt = new TLatex(0.5,0.85,Form("p_{T} < %3.1f GeV/c",ptmax));
  else
    lpt = new TLatex(0.5,0.85,Form("%3.1f < p_{T} < %3.1f GeV/c",ptmin,ptmax));
  TLatex *lrap;
  if (absRapidity){
    if (ymin==0.0)
      lrap = new TLatex(0.5,0.80,Form("|y| < %3.1f",ymax));
    else
      lrap = new TLatex(0.5,0.80,Form("%3.1f < |y| < %3.1f",ymin,ymax));
  }
  else {
    if (ymin==0.0)
      lrap = new TLatex(0.5,0.80,Form("y < %3.1f",ymax));
    else
      lrap = new TLatex(0.5,0.80,Form("%3.1f < y < %3.1f",ymin,ymax));
  }

  lpt->SetNDC(kTRUE);
  lrap->SetNDC(kTRUE);

  TCanvas *c1 = new TCanvas("c1","c1",1000,600);
  c1->Divide(2,1);
  c1->GetPad(1)->SetLogy();
  c1->cd(1);

  //  TFile *_file0 = TFile::Open("/data/CMS/MC/pp/NPMC_Psi2S_GlbGlb_Histos_muLessPV.root");
  TChain *NPTree = new TChain("myTree");
  // PV includes the muons
  //  NPTree->Add("/data/CMS/MC/pp/NPpsi2SMC_Histos_2013pp_GlbGlb_STARTHI53_V28-v1.root");
  // PV without the muons
  NPTree->Add("/data/CMS/MC/pp/NPpsi2SMC_Histos_2013pp_GlbGlb_STARTHI53_V28-v1_allPV.root");

  TH1F *hNPorig = new TH1F("hNPorig","hNPorig;l_{J/#psi} (mm);counts",1500,-5,10);
  hNPorig->Sumw2();
  hNPorig->SetLineColor(2);
  TH1F *hNP = new TH1F("hNP","hNP;l_{#psi(2S)} (mm);counts",1500,-5,10);
  hNP->SetLineColor(2);
  hNP->GetXaxis()->SetRangeUser(-1,3);
  TH1F* hNPCut = (TH1F*) hNP->Clone("hNPCut");
  TH1F* hNPEff = (TH1F*) hNP->Clone("hNPEff");
  hNPCut->GetYaxis()->SetTitle("#varepsilon_{prompt} or 1-#varepsilon_{non-prompt}");
  hNPEff->GetYaxis()->SetTitle("#varepsilon_{non-prompt}");

  //  hNP->SetDirectory(0);
  //  TH1::AddDirectory(kFALSE);

  NPTree->Draw("Reco_QQ_ctau>>hNPorig",defaultCut&&ptCut&&rapCut&&trigCut,"e");

  //  TFile *_file1 = TFile::Open("/data/CMS/MC/pp/PRMC_Psi2S_GlbGlb_Histos_muLessPV.root");
  TChain *PRTree = new TChain("myTree");
  // PV includes the muons
  //  PRTree->Add("/data/CMS/MC/pp/PRpsi2SMC_Histos_2013pp_GlbGlb_STARTHI53_V28-v1.root");
  // PV without the muons
  PRTree->Add("/data/CMS/MC/pp/PRpsi2SMC_Histos_2013pp_GlbGlb_STARTHI53_V28-v1_muLessPV.root");

  TH1F *hPRorig = new TH1F("hPRorig","hPRorig;l_{J/#psi} (mm);counts",1500,-5,10);
  hPRorig->Sumw2();
  hPRorig->SetLineColor(1);

  TH1F *hPR = new TH1F("hPR","hPR;l_{#psi(2S)} (mm);counts",1500,-5,10);
  hPR->SetLineColor(1);

  //  hPR->SetDirectory(0);

  TH1F* hPRCut = (TH1F*) hPR->Clone("hPRCut");

  PRTree->Draw("Reco_QQ_ctau>>hPRorig",defaultCut&&ptCut&&rapCut&&trigCut,"e");

  for (int i=1; i<=hPRorig->GetNbinsX();++i) {
    float lJpsi = hPRorig->GetBinCenter(i);
    float lpsi2s;
    for (int j=0; j<hPRorig->GetBinContent(i); ++j) {
      lpsi2s = lJpsi*3.686109/3.096916+hPRorig->GetBinWidth(i)*(gRandom->Rndm()-0.5);
      hPR->Fill(lpsi2s);
    }

    for (int j=0; j<hNPorig->GetBinContent(i); ++j) {
      lpsi2s = lJpsi*3.686109/3.096916+hNPorig->GetBinWidth(i)*(gRandom->Rndm()-0.5);
      hNP->Fill(lpsi2s);
    }

  }
  // activate Sumw2 after filling to ignore weights and set errors to sqrt(N)
  hNP->Sumw2();
  hPR->Sumw2();

  //  _file0->cd();
  hNP->Scale(1.0/hNP->GetEntries());
  hNP->SetMaximum(0.2);
  hNP->Draw("hist");
  //  _file1->cd();
  hPR->Scale(1.0/hPR->GetEntries());
  hPR->Draw("histsame");
  leg->AddEntry(hPR,"Prompt #psi(2S)","L");
  lpt->Draw();
  lrap->Draw();


  c1->cd(2);
  //  _file0->cd();
  leg->AddEntry(hNP,"Non-prompt #psi(2S)","L");
  bool stillLooking = true;
  int lCutBin=0;
  for (int i=1;i<=hNP->GetNbinsX();++i) {
    hNPCut->SetBinContent(i,hNP->Integral(i,hNP->GetNbinsX()+1)/hNP->Integral(1,hNP->GetNbinsX()+1));
    hNPEff->SetBinContent(i,hNP->Integral(1,i)/hNP->Integral(1,hNP->GetNbinsX()+1));
    hNPCut->SetBinError(i,0.0);
    hNPEff->SetBinError(i,0.0);
    if (hNPCut->GetBinContent(i)<threshold && stillLooking) {
      lCut = hNPCut->GetBinLowEdge(i);
      lCutBin = i;
      stillLooking = false;
    }
  }
  hNPCut->Draw("hist");

  //  _file1->cd();
  for (int i=1;i<=hPR->GetNbinsX();++i) {
    hPRCut->SetBinContent(i,hPR->Integral(1,i)/hPR->Integral(1,hPR->GetNbinsX()+1));
    hPRCut->SetBinError(i,0.0);
  }
  hPRCut->Draw("histsame");

  TH1F *hRatio = (TH1F*) hNPEff->Clone("hRatio");
  TH1F *hRatio_fb1 = (TH1F*) hNPEff->Clone("hRatio_fb1");
  TH1F *hRatio_fb2 = (TH1F*) hNPEff->Clone("hRatio_fb2");
  TH1F *hRatio_fb3 = (TH1F*) hNPEff->Clone("hRatio_fb3");
  TH1F *hRatio_fb5 = (TH1F*) hNPEff->Clone("hRatio_fb5");
  TH1F *hRatio_fb7 = (TH1F*) hNPEff->Clone("hRatio_fb7");
  hRatio->Divide(hPRCut);
  hRatio->SetLineColor(4);
  //  hRatio->Draw("histsame");

  for (int i=1;i<=hRatio_fb1->GetNbinsX();++i) {
    if (hPRCut->GetBinContent(i)>0.0) {
      hRatio_fb1->SetBinContent(i, 0.9*hPRCut->GetBinContent(i)/sqrt(0.1*hNPEff->GetBinContent(i)+0.9*hPRCut->GetBinContent(i)) );
      hRatio_fb2->SetBinContent(i, 0.8*hPRCut->GetBinContent(i)/sqrt(0.2*hNPEff->GetBinContent(i)+0.8*hPRCut->GetBinContent(i)) );
      hRatio_fb3->SetBinContent(i, 0.7*hPRCut->GetBinContent(i)/sqrt(0.3*hNPEff->GetBinContent(i)+0.7*hPRCut->GetBinContent(i)) );
      hRatio_fb5->SetBinContent(i, 0.5*hPRCut->GetBinContent(i)/sqrt(0.5*hNPEff->GetBinContent(i)+0.5*hPRCut->GetBinContent(i)) );
      hRatio_fb7->SetBinContent(i, 0.3*hPRCut->GetBinContent(i)/sqrt(0.7*hNPEff->GetBinContent(i)+0.3*hPRCut->GetBinContent(i)) );
    }
    else {
      hRatio_fb1->SetBinContent(i,0.0);
      hRatio_fb2->SetBinContent(i,0.0);
      hRatio_fb3->SetBinContent(i,0.0);
      hRatio_fb5->SetBinContent(i,0.0);
      hRatio_fb7->SetBinContent(i,0.0);
    } 
    hRatio_fb1->SetBinError(i,0.0);
    hRatio_fb2->SetBinError(i,0.0);
    hRatio_fb3->SetBinError(i,0.0);
    hRatio_fb5->SetBinError(i,0.0);
    hRatio_fb7->SetBinError(i,0.0);
 }
  hRatio_fb1->SetLineColor(1);
  hRatio_fb2->SetLineColor(2);
  hRatio_fb3->SetLineColor(4);
  hRatio_fb5->SetLineColor(kGreen+2);
  hRatio_fb7->SetLineColor(kOrange+2);

  std::cout << "\nFixing rejection:" << std::endl;
  std::cout << "Non-prompt J/psi rejection of " << threshold << " at l_Jpsi = " << lCut << std::endl;
  std::cout << "Prompt J/psi efficiency of " <<  hPRCut->GetBinContent(lCutBin) << " at l_Jpsi = " << lCut << std::endl;

  stillLooking = true;
  for (int i=1;i<=hPRCut->GetNbinsX();++i) {
    if (hPRCut->GetBinContent(i)>efficiency && stillLooking) {
      lCut = hPRCut->GetBinLowEdge(i);
      lCutBin = i;
      break;
    }
  }
  std::cout << "\nFixing efficiency:" << std::endl;
  std::cout << "Non-prompt psi(2S) rejection of " <<  hNPCut->GetBinContent(lCutBin) << " at l_Jpsi = " << lCut << std::endl;
  std::cout << "Prompt psi(2S) efficiency of " <<  efficiency << " at l_Jpsi = " << lCut << std::endl;

  TLine *l1 = new TLine(-1.0,hNPCut->GetBinContent(lCutBin),lCut,hNPCut->GetBinContent(lCutBin));
  l1->SetLineStyle(2);
  l1->SetLineColor(2);
  l1->Draw();

  TLine *l2 = new TLine(lCut,0.0,lCut,1.0);
  l2->SetLineStyle(2);
  l2->Draw();

  TLine *l3 = new TLine(-1.0,efficiency,lCut,efficiency);
  l3->SetLineStyle(2);
  l3->Draw();


  TLatex *lcut = new TLatex(0.5,0.88,Form("Keep only l_{#psi(2S)} < %4.3f mm:",lCut));
  TLatex *leff = new TLatex(0.5,0.80,Form("#varepsilon_{prompt} = %4.3f",efficiency));
  TLatex *lrej = new TLatex(0.5,0.75,Form("1-#varepsilon_{non-prompt} = %4.3f",hNPCut->GetBinContent(lCutBin)));
  lcut->Draw();
  leff->Draw();
  lrej->Draw();

  leg->Draw();

  if (savePlot)
    c1->SaveAs(fname);

  TCanvas *c2 = new TCanvas("c2","c2");
  c2->cd();
  hRatio_fb1->Draw("hist");
  hRatio_fb2->Draw("histsame");
  hRatio_fb3->Draw("histsame");
  hRatio_fb5->Draw("histsame");
  hRatio_fb7->Draw("histsame");

  std::cout << "fb = 0.1: " << hRatio_fb1->GetXaxis()->GetBinLowEdge(hRatio_fb1->GetMaximumBin()) << std::endl;
  std::cout << "fb = 0.2: " << hRatio_fb2->GetXaxis()->GetBinLowEdge(hRatio_fb2->GetMaximumBin()) << std::endl;
  std::cout << "fb = 0.3: " << hRatio_fb3->GetXaxis()->GetBinLowEdge(hRatio_fb3->GetMaximumBin()) << std::endl;
  std::cout << "fb = 0.5: " << hRatio_fb5->GetXaxis()->GetBinLowEdge(hRatio_fb5->GetMaximumBin()) << std::endl;
  std::cout << "fb = 0.7: " << hRatio_fb7->GetXaxis()->GetBinLowEdge(hRatio_fb7->GetMaximumBin()) << std::endl;

  return;
}
