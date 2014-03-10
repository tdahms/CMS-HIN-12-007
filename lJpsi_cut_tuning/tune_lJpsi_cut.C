#include <iostream>
#include "TCanvas.h"
#include "TCut.h"
#include "TLegend.h"
#include "TChain.h"
#include "TH1.h"
#include "TLine.h"
#include "TLatex.h"
#include "TFile.h"

void tune_lJpsi_cut(bool isHI=false, float ptmin=0.0, float ptmax=30.0, float ymin=0.0, float ymax=2.4, bool absRapidity=true, double threshold=0.8, double efficiency=0.9, bool savePlot=false, int centmin=0, int centmax=40)
{
  float lCut=-10.0;
  //  float efficiency = 1.0;
  TLegend *leg = new TLegend(0.42,0.4,0.95,0.6);
  if (isHI)
    leg->SetHeader("PYTHIA: PbPb #sqrt{s_{NN}} = 2.76 TeV");
  else
    leg->SetHeader("PYTHIA: pp #sqrt{s} = 2.76 TeV");
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetMargin(0.15);
  leg->SetTextSize(0.035);

  //  TCut defaultCut = "Reco_QQ_sign==0&&Gen_QQ_4mom.M()>3.09&&Gen_QQ_4mom.M()<3.10&&Reco_QQ_4mom.M()>2.9&&Reco_QQ_4mom.M()<3.2&&Reco_QQ_ctauTrue>-10";
  TCut defaultCut = "Reco_QQ_sign==0&&Reco_QQ_4mom.M()>2.5&&Reco_QQ_4mom.M()<3.35&&Reco_QQ_ctauTrue>-10";
  TCut ptCut = Form("Reco_QQ_4mom.Pt()>%3.1f&&Reco_QQ_4mom.Pt()<%3.1f",ptmin,ptmax);
  TCut rapCut;
  if (absRapidity)
    rapCut = Form("abs(Reco_QQ_4mom.Rapidity())>%3.1f&&abs(Reco_QQ_4mom.Rapidity())<%3.1f",ymin,ymax);
  else
    rapCut = Form("Reco_QQ_4mom.Rapidity()>%3.1f&&Reco_QQ_4mom.Rapidity()<%3.1f",ymin,ymax);
  TCut centCut;
  if (isHI)
    centCut = Form("Centrality>=%d&&Centrality<%d",centmin,centmax);

  unsigned int trigBit;
  if (isHI)
    trigBit=1; // DoubleMu0_HighQ
  else
    trigBit=2; // DoubleMu0_HighQ
  TCut trigCut = Form("(HLTriggers&%u)==%u&&(Reco_QQ_trig&%u)==%u",trigBit,trigBit,trigBit,trigBit);

  TString fname;
  TString outfname;
  if (isHI) {
    //    fname = "PbPb_eff_" + efficiency + "_Rap_" + ymin + "-" + ymax + "_Pt_" + ptmin + "-" + ptmax + ".pdf";
    fname = Form("20140210/Jpsi_PbPb_eff_%3.1f_Rap_%3.1f-%3.1f_Pt_%3.1f-%3.1f_Cent_%d-%d.pdf",efficiency,ymin,ymax,ptmin,ptmax,int(centmin*2.5),int(centmax*2.5));
    outfname = Form("20140210/Jpsi_PbPb_eff_%3.1f_Rap_%3.1f-%3.1f_Pt_%3.1f-%3.1f_Cent_%d-%d.root",efficiency,ymin,ymax,ptmin,ptmax,int(centmin*2.5),int(centmax*2.5));
  } 
  else {
    fname = Form("20140210/Jpsi_pp_eff_%3.1f_Rap_%3.1f-%3.1f_Pt_%3.1f-%3.1f.pdf",efficiency,ymin,ymax,ptmin,ptmax);
    outfname = Form("20140210/Jpsi_pp_eff_%3.1f_Rap_%3.1f-%3.1f_Pt_%3.1f-%3.1f.root",efficiency,ymin,ymax,ptmin,ptmax);
  }
  std::cout << fname << std::endl;


  std::cout << "default: " << defaultCut.GetTitle() << std::endl;
  std::cout << "pt cut: " << ptCut.GetTitle() << std::endl;
  std::cout << "rapidity cut: " << rapCut.GetTitle() << std::endl;
  std::cout << "centrality cut: " << centCut.GetTitle() << std::endl;
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


  TLatex *lcent = NULL;
  if (isHI) {
    lcent = new TLatex(0.5,0.75,Form("Cent. %d-%d%%",int(centmin*2.5),int(centmax*2.5)));
    lcent->SetNDC(kTRUE);
  }

  TCanvas *c1 = new TCanvas("c1","c1",1000,600);
  c1->Divide(2,1);
  c1->GetPad(1)->SetLogy();
  c1->cd(1);

  //  TCut prompt_weight = ;
  //  TCut nonprompt_weight = "((Gen_QQ_4mom.Pt()<3.0)*(1.0/171458.0))||((Gen_QQ_4mom.Pt()>3.0&&Gen_QQ_4mom.Pt()<6.0)*(2.25*0.623/172208.0))||((Gen_QQ_4mom.Pt()>6.0&&Gen_QQ_4mom.Pt()<9.0)*(2.55*0.245/171178.0))||((Gen_QQ_4mom.Pt()>9.0&&Gen_QQ_4mom.Pt()<12.0)*(2.45*0.080/172208.0))||((Gen_QQ_4mom.Pt()>12.0&&Gen_QQ_4mom.Pt()<15.0)*(2.0*0.0299/170670.0))||((Gen_QQ_4mom.Pt()>15.0&&Gen_QQ_4mom.Pt()<30.0)*(3.35*0.0109/169302.0))"


  TChain *NPTree = new TChain("myTree");
  if (isHI) {
    NPTree->Add("/data/CMS/MC/PbPb/bJpsiMuMu_JpsiPt03_Histos_cmssw445p5_RegIt_hStats.root");
    NPTree->Add("/data/CMS/MC/PbPb/bJpsiMuMu_JpsiPt36_Histos_cmssw445p5_RegIt_hStats.root");
    NPTree->Add("/data/CMS/MC/PbPb/bJpsiMuMu_JpsiPt69_Histos_cmssw445p5_RegIt_hStats.root");
    NPTree->Add("/data/CMS/MC/PbPb/bJpsiMuMu_JpsiPt912_Histos_cmssw445p5_RegIt_hStats.root");
    NPTree->Add("/data/CMS/MC/PbPb/bJpsiMuMu_JpsiPt1215_Histos_cmssw445p5_RegIt_hStats.root");
    NPTree->Add("/data/CMS/MC/PbPb/bJpsiMuMu_JpsiPt1530_Histos_cmssw445p5_RegIt_hStats.root");
  }
  else
    // PV includes the muons
    //    NPTree->Add("/data/CMS/MC/pp/NPMC_Histos_2013pp_GlbGlb_STARTHI53_V28-v1.root");
    // PV without the muons
    NPTree->Add("/data/CMS/MC/pp/NPMC_Histos_2013pp_GlbGlb_STARTHI53_V28-v1_muLessPV.root");

  // TFile *_file0;
  // if (isHI) {
  //   _file0 = TFile::Open("/data/CMS/MC/PbPb/NPMC_Histos_2013pp_GlbGlb_STARTHI53_V28-v1.root");
  // else
  //   _file0 = TFile::Open("/data/CMS/MC/pp/NPMC_Histos_2013pp_GlbGlb_STARTHI53_V28-v1.root");

  TH1F *hNP = new TH1F("hNP","hNP;l_{J/#psi} (mm);counts",1500,-5,10);
  hNP->Sumw2();
  hNP->SetLineColor(2);
  hNP->GetXaxis()->SetRangeUser(-1,3);
  TH1F* hNPCut = (TH1F*) hNP->Clone("hNPCut");
  TH1F* hNPEff = (TH1F*) hNP->Clone("hNPEff");
  hNPCut->GetYaxis()->SetTitle("#varepsilon_{prompt} or 1-#varepsilon_{non-prompt}");
  hNPEff->GetYaxis()->SetTitle("#varepsilon_{non-prompt}");

  //  hNP->SetDirectory(0);
  //  TH1::AddDirectory(kFALSE);

  NPTree->Draw("Reco_QQ_ctau>>hNP",defaultCut&&ptCut&&rapCut&&centCut&&trigCut,"e");

  TChain *PRTree = new TChain("myTree");
  if (isHI) {
    PRTree->Add("/data/CMS/MC/PbPb/jpsiMuMu_JpsiPt03_Histos_cmssw445p5_RegIt_hStats.root");
    PRTree->Add("/data/CMS/MC/PbPb/jpsiMuMu_JpsiPt36_Histos_cmssw445p5_RegIt_hStats.root");
    PRTree->Add("/data/CMS/MC/PbPb/jpsiMuMu_JpsiPt69_Histos_cmssw445p5_RegIt_hStats.root");
    PRTree->Add("/data/CMS/MC/PbPb/jpsiMuMu_JpsiPt912_Histos_cmssw445p5_RegIt_hStats.root");
    PRTree->Add("/data/CMS/MC/PbPb/jpsiMuMu_JpsiPt1215_Histos_cmssw445p5_RegIt_hStats.root");
    PRTree->Add("/data/CMS/MC/PbPb/jpsiMuMu_JpsiPt1530_Histos_cmssw445p5_RegIt_hStats.root");
  }
  else
    // PV includes the muons
    //    PRTree->Add("/data/CMS/MC/pp/PRMC_Histos_2013pp_GlbGlb_STARTHI53_V28-v1.root");
    // PV without the muons
    PRTree->Add("/data/CMS/MC/pp/PRMC_Histos_2013pp_GlbGlb_STARTHI53_V28-v1_muLessPV.root");
  // TFile *_file1;
  // if (isHI)
  //   _file1 = TFile::Open("/data/CMS/MC/pp/PRMC_Histos_2013pp_GlbGlb_STARTHI53_V28-v1.root");
  // else
  //   _file1 = TFile::Open("/data/CMS/MC/pp/PRMC_Histos_2013pp_GlbGlb_STARTHI53_V28-v1.root");

  TH1F *hPR = new TH1F("hPR","hPR;l_{J/#psi} (mm);counts",1500,-5,10);
  hPR->Sumw2();
  hPR->SetLineColor(1);

  //  hPR->SetDirectory(0);

  TH1F* hPRCut = (TH1F*) hPR->Clone("hPRCut");

  PRTree->Draw("Reco_QQ_ctau>>hPR",defaultCut&&ptCut&&rapCut&&centCut&&trigCut,"e");

  //  _file0->cd();
  if (!isHI) {
    hNP->Scale(1.0/hNP->GetEntries());
    hPR->Scale(1.0/hPR->GetEntries());
  }

  hNP->SetMaximum(0.2);
  hNP->Draw("hist");
  //  _file1->cd();
  hPR->Draw("histsame");
  leg->AddEntry(hPR,"Prompt J/#psi","L");
  lpt->Draw();
  lrap->Draw();
  if (isHI)
    lcent->Draw();

  c1->cd(2);
  //  _file0->cd();
  leg->AddEntry(hNP,"Non-prompt J/#psi","L");
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
  std::cout << "Non-prompt J/psi rejection of " <<  hNPCut->GetBinContent(lCutBin) << " at l_Jpsi = " << lCut << std::endl;
  std::cout << "Prompt J/psi efficiency of " <<  efficiency << " at l_Jpsi = " << lCut << std::endl;

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


  TLatex *lcut = new TLatex(0.5,0.88,Form("Keep only l_{J/#psi} < %4.3f mm:",lCut));
  TLatex *leff = new TLatex(0.5,0.80,Form("#varepsilon_{prompt} = %4.3f",efficiency));
  TLatex *lrej = new TLatex(0.5,0.75,Form("1-#varepsilon_{non-prompt} = %4.3f",hNPCut->GetBinContent(lCutBin)));
  lcut->Draw();
  leff->Draw();
  lrej->Draw();

  leg->Draw();

  std::cout << "RMS = " << hPR->GetRMS()*1000 << " μm" << std::endl;

  if (savePlot) {
    c1->SaveAs(fname);

    TFile *outf = new TFile(outfname,"RECREATE");
    hPR->Write();
    hNP->Write();
    outf->Close();
  }

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
