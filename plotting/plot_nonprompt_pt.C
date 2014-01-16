void plot_nonprompt_pt()
{
  // TH1F *h0 = new TH1F("h0","h0;p_{T} [GeV/c];counts",500,0,50);
  // TH1F *h1 = new TH1F("h1","h1;p_{T} [GeV/c];counts",500,0,50);
  // TH1F *h2 = new TH1F("h2","h2;p_{T} [GeV/c];counts",500,0,50);
  // TH1F *h3 = new TH1F("h3","h3;p_{T} [GeV/c];counts",500,0,50);
  // TH1F *h4 = new TH1F("h4","h4;p_{T} [GeV/c];counts",500,0,50);
  // TH1F *h5 = new TH1F("h5","h5;p_{T} [GeV/c];counts",500,0,50);
  // TH1F *h6 = new TH1F("h6","h6;p_{T} [GeV/c];counts",500,0,50);
  // h0->Sumw2();
  // h1->Sumw2();
  // h2->Sumw2();
  // h3->Sumw2();
  // h4->Sumw2();
  // h5->Sumw2();
  // h6->Sumw2();
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->SetLogy();

  TCut genCuts = "abs(Gen_QQ_4mom.Rapidity())<2.4";
  TCut muPlusAcc = "abs(Gen_QQ_mupl_4mom.Rapidity()<2.4)&&Gen_QQ_mupl_4mom.P()>2.5";
  TCut muMinusAcc = "abs(Gen_QQ_mumi_4mom.Rapidity()<2.4)&&Gen_QQ_mumi_4mom.P()>2.5";
  TFile *f0 = new TFile("/data/CMS/MC/PbPb/bJpsiMuMu_JpsiPt03_Histos_cmssw445p1_RegIt.root","READ");
  myTree->Draw("Reco_QQ_4mom.Pt()>>h0(500,0,50)",(genCuts&&muPlusAcc&&muMinusAcc&&"Gen_QQ_4mom.Pt()>0.0&&Gen_QQ_4mom.Pt()<3.0")*"(1.0/171458.0)","e");
  //  h0->Scale(1.0/171458.0);
  //  h0->Scale(1.0/171458.0);
  h0->SetLineColor(2);
  TH1F *hpt = (TH1F*)h0->Clone("hpt");
  hpt->SetLineColor(1);
  hpt->SetMarkerStyle(24);
  hpt->SetDirectory(0);

  TFile *f1 = new TFile("/data/CMS/MC/PbPb/bJpsiMuMu_JpsiPt36_Histos_cmssw445p1_RegIt.root","READ");
  myTree->Draw("Reco_QQ_4mom.Pt()>>h1(500,0,50)",(genCuts&&muPlusAcc&&muMinusAcc&&"Gen_QQ_4mom.Pt()>3.0&&Gen_QQ_4mom.Pt()<6.0")*"(2.25*0.623/172208.0)","e");
  //  h1->Scale(2.25*0.623/172208.0);
  h1->SetLineColor(4);
  hpt->Add(h1);

  TFile *f2 = new TFile("/data/CMS/MC/PbPb/bJpsiMuMu_JpsiPt69_Histos_cmssw445p1_RegIt.root","READ");
  myTree->Draw("Reco_QQ_4mom.Pt()>>h2(500,0,50)",(genCuts&&muPlusAcc&&muMinusAcc&&"Gen_QQ_4mom.Pt()>6.0&&Gen_QQ_4mom.Pt()<9.0")*"(2.55*0.245/171178.0)","e");
  //  h2->Scale(2.55*0.245/171178.0);
  h2->SetLineColor(kGreen+2);
  hpt->Add(h2);

  TFile *f3 = new TFile("/data/CMS/MC/PbPb/bJpsiMuMu_JpsiPt912_Histos_cmssw445p1_RegIt.root","READ");
  myTree->Draw("Reco_QQ_4mom.Pt()>>h3(500,0,50)",(genCuts&&muPlusAcc&&muMinusAcc&&"Gen_QQ_4mom.Pt()>9.0&&Gen_QQ_4mom.Pt()<12.0")*"(2.45*0.080/172208.0)","e");
  //  h3->Scale(2.45*0.080/172208.0);
  h3->SetLineColor(kOrange+2);
  hpt->Add(h3);

  TFile *f4 = new TFile("/data/CMS/MC/PbPb/bJpsiMuMu_JpsiPt1215_Histos_cmssw445p1_RegIt.root","READ");
  myTree->Draw("Reco_QQ_4mom.Pt()>>h4(500,0,50)",(genCuts&&muPlusAcc&&muMinusAcc&&"Gen_QQ_4mom.Pt()>12.0&&Gen_QQ_4mom.Pt()<15.0")*"(2.0*0.0299/170670.0)","e");
  //  h4->Scale(2.0*0.0299/170670.0);
  h4->SetLineColor(kMagenta+1);
  hpt->Add(h4);

  TFile *f5 = new TFile("/data/CMS/MC/PbPb/bJpsiMuMu_JpsiPt1530_Histos_cmssw445p1_RegIt.root","READ");
  myTree->Draw("Reco_QQ_4mom.Pt()>>h5(500,0,50)",(genCuts&&muPlusAcc&&muMinusAcc&&"Gen_QQ_4mom.Pt()>15.0&&Gen_QQ_4mom.Pt()<30.0")*"(3.35*0.0109/169302.0)","e");
  //  h5->Scale(3.35*0.0109/169302.0);
  h5->SetLineColor(kCyan-2);
  hpt->Add(h5);

  hpt->SetTitle(";p_{T} [GeV/c];counts");
  hpt->Draw("e");
  f0->cd();
  h0->Draw("histsame");
  f1->cd();
  h1->Draw("histsame");
  f2->cd();
  h2->Draw("histsame");
  f3->cd();
  h3->Draw("histsame");
  f4->cd();
  h4->Draw("histsame");
  f5->cd();
  h5->Draw("histsame");
  return;
}
