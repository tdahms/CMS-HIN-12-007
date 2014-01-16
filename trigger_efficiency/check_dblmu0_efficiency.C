void check_dblmu0_efficiency()
{
  TCanvas *c1 = new TCanvas("c1","c1",1000,600);
  c1->Divide(2,1);
  c1->GetPad(1)->SetLogy(1);
  //  c1->GetPad(2)->SetLogy(1);

  c1->cd(1);
  TFile *_file0 = TFile::Open("/data/CMS/MC/pp/PRMC_Jpsi_GlbGlb_Histos_muLessPV.root");
  
  myTree->Draw("Reco_QQ_4mom.Rapidity()>>hgen_old(100,-2.5,2.5)","Reco_QQ_sign==0&&Reco_QQ_ctauTrue>-10&&(Reco_QQ_trig&1)==1","e");
  myTree->Draw("Reco_QQ_4mom.Rapidity()>>hreco_old(100,-2.5,2.5)","Reco_QQ_sign==0&&Reco_QQ_ctauTrue>-10&&(Reco_QQ_trig&2)==2","esame");
  hgen_old->Scale(1.0/hgen_old->GetEntries());
  hreco_old->Scale(1.0/hgen_old->GetEntries());

  hreco_old->SetMarkerColor(kGreen+2);

  hgen_old->SetTitle(";y_{J/#psi};counts [arb. units]");

  c1->cd(2);
  TGraphAsymmErrors *eff_old = new TGraphAsymmErrors();
  eff_old->BayesDivide(hreco_old,hgen_old);
  eff_old->Draw("AP");
  eff_old->SetTitle(";y_{J/#psi};N(DoubleMu0_HighQ) / N(DoubleMuOpen)");
  eff_old->GetXaxis()->SetRangeUser(-2.5,2.5);

  c1->cd(1);
  TFile *_file1 = TFile::Open("/data/CMS/MC/pp/PRMC_Histos_2013pp_GlbGlb_STARTHI53_V28-v1.root");
  
  myTree->Draw("Reco_QQ_4mom.Rapidity()>>hgen_new(100,-2.5,2.5)","Reco_QQ_sign==0&&Reco_QQ_ctauTrue>-10&&(Reco_QQ_trig&1)==1","esame");
  myTree->Draw("Reco_QQ_4mom.Rapidity()>>hreco_new(100,-2.5,2.5)","Reco_QQ_sign==0&&Reco_QQ_ctauTrue>-10&&(Reco_QQ_trig&2)==2","esame");

  hgen_new->Scale(1.0/hgen_new->GetEntries());
  hreco_new->Scale(1.0/hgen_new->GetEntries());

  hgen_new->SetMarkerStyle(25);
  hgen_new->SetMarkerColor(2);

  hreco_new->SetMarkerStyle(25);
  hreco_new->SetMarkerColor(4);

  hgen_old->SetTitle(";y_{J/#psi};counts [arb. units]");

  c1->cd(2);
  TGraphAsymmErrors *eff_new = new TGraphAsymmErrors();
  eff_new->BayesDivide(hreco_new,hgen_new);
  eff_new->Draw("P");
  eff_new->SetMarkerStyle(24);
  eff_new->SetMarkerColor(4);


  return;
}
