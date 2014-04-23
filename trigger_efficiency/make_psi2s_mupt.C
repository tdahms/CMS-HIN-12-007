#include <iostream>
#include "TCanvas.h"
#include "TCut.h"
#include "TChain.h"
#include "TString.h"
#include "TH1.h"
#include "TFile.h"

void make_psi2s_mupt(float ptmin=0.0, float ptmax=30.0, float ymin=0.0, float ymax=2.4, bool absRapidity=true, bool saveFile=false)
{
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->SetLogy();
  TCut defaultCut = "Reco_QQ_sign==0&&Reco_QQ_4mom.M()>3.35&&Reco_QQ_4mom.M()<4.0&&Reco_QQ_ctauTrue>-10";
  TCut ptCut = Form("Reco_QQ_4mom.Pt()>%4.1f&&Reco_QQ_4mom.Pt()<%4.1f",ptmin,ptmax);
  TCut rapCut;
  if (absRapidity)
    rapCut = Form("abs(Reco_QQ_4mom.Rapidity())>%3.1f&&abs(Reco_QQ_4mom.Rapidity())<%3.1f",ymin,ymax);
  else
    rapCut = Form("Reco_QQ_4mom.Rapidity()>%3.1f&&Reco_QQ_4mom.Rapidity()<%3.1f",ymin,ymax);

  unsigned int trigBit=2; // DoubleMu0_HighQ
  TCut trigCut = Form("(HLTriggers&%u)==%u&&(Reco_QQ_trig&%u)==%u",trigBit,trigBit,trigBit,trigBit);

  TString fname;
  fname = Form("20140324/MC_psi2s_pp_mupt_Rap_%3.1f-%3.1f_Pt_%3.1f-%3.1f.root",ymin,ymax,ptmin,ptmax);

  std::cout << fname << std::endl;


  std::cout << "default: " << defaultCut.GetTitle() << std::endl;
  std::cout << "pt cut: " << ptCut.GetTitle() << std::endl;
  std::cout << "rapidity cut: " << rapCut.GetTitle() << std::endl;
  std::cout << "trigger bit: " << trigCut.GetTitle() << std::endl;

  TH1F *hMuPlPtRec = new TH1F("hMuPlPtRec","hMuPlPtRec;p_{T} (#mu^{+}) (GeV/c);Events",200,0,20);
  TH1F *hMuMiPtRec = new TH1F("hMuMiPtRec","hMuMiPtRec;p_{T} (#mu^{-}) (GeV/c);Events",200,0,20);
  TH1F *hMuPtRec = new TH1F("hMuPtRec","hMuPtRec;p_{T} (#mu^{#pm}) (GeV/c);Events",200,0,20);
  hMuPlPtRec->Sumw2();
  hMuMiPtRec->Sumw2();
  hMuPtRec->Sumw2();

  hMuPlPtRec->SetMarkerColor(kRed);
  hMuMiPtRec->SetMarkerColor(kBlue);

  TChain *myTree = new TChain("myTree");
  myTree->Add("../root_files/PRpsi2SMC_Histos_2013pp_GlbGlb_STARTHI53_V28-v1_muLessPV.root");

  myTree->Draw("Reco_QQ_mupl_4mom.Pt()>>hMuPlPtRec",defaultCut&&ptCut&&rapCut&&trigCut,"e");
  myTree->Draw("Reco_QQ_mumi_4mom.Pt()>>hMuMiPtRec",defaultCut&&ptCut&&rapCut&&trigCut,"e");

  hMuPtRec->Add(hMuPlPtRec,hMuMiPtRec);

  hMuPtRec->Draw();
  hMuPlPtRec->Draw("same");
  hMuMiPtRec->Draw("same");

  TFile *outf = NULL;
  if (saveFile) {
    outf = new TFile(fname,"RECREATE");
    hMuPtRec->Write();
    outf->Close();
  }

  return;
}
