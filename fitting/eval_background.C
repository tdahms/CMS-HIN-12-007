void eval_background(bool isHI=true)
{
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->SetLogy();
  TFile *inf;
  if (isHI) {
    inf = new TFile("/Users/tdahms/CMS/HIN-12-007/Jpsi_Histos_PbPb_RegIT_glbglb.root","READ");
    myTree->Draw("Reco_QQ_4mom.M()>>h1(100,2.2,4.2)","(HLTriggers&1)==1&&(Reco_QQ_trig&1)==1&&Reco_QQ_sign==0&&Reco_QQ_4mom.Pt()>6.5&&Reco_QQ_4mom.Pt()<30.0&&abs(Reco_QQ_4mom.Rapidity())>0.0&&abs(Reco_QQ_4mom.Rapidity())<1.6","e");
    myTree->Draw("Reco_QQ_4mom.M()>>h2(100,2.2,4.2)","(HLTriggers&1)==1&&(Reco_QQ_trig&1)==1&&Reco_QQ_sign!=0&&Reco_QQ_4mom.Pt()>6.5&&Reco_QQ_4mom.Pt()<30.0&&abs(Reco_QQ_4mom.Rapidity())>0.0&&abs(Reco_QQ_4mom.Rapidity())<1.6","e");
  }
  else {
    inf = new TFile("/Users/tdahms/CMS/HIN-12-007/Jpsi_Histos_pp_glbglb_woPileUpRej_muLessPV.root","READ");
    myTree->Draw("Reco_QQ_4mom.M()>>h1(100,2.2,4.2)","(HLTriggers&2)==2&&(Reco_QQ_trig&2)==2&&Reco_QQ_sign==0&&Reco_QQ_4mom.Pt()>6.5&&Reco_QQ_4mom.Pt()<30.0&&abs(Reco_QQ_4mom.Rapidity())>0.0&&abs(Reco_QQ_4mom.Rapidity())<1.6","e");
    myTree->Draw("Reco_QQ_4mom.M()>>h2(100,2.2,4.2)","(HLTriggers&2)==2&&(Reco_QQ_trig&2)==2&&Reco_QQ_sign!=0&&Reco_QQ_4mom.Pt()>6.5&&Reco_QQ_4mom.Pt()<30.0&&abs(Reco_QQ_4mom.Rapidity())>0.0&&abs(Reco_QQ_4mom.Rapidity())<1.6","e");
  }


  h2->SetMarkerColor(kRed);
  h2->SetMarkerStyle(24);

  h1->SetTitle(";m_{#mu#mu} (GeV/c^{2});counts");
  h1->GetXaxis()->SetTitleSize(0.04);
  h1->GetXaxis()->SetTitleOffset(1.3);
  h1->GetXaxis()->CenterTitle(true);
  h1->Draw();
  h2->Draw("same");

  if (isHI)
    c1->SaveAs("PbPb_mass_os_ss_rap0-16_pT65-30_cent0-100.pdf");
  else
    c1->SaveAs("pp_mass_os_ss_rap0-16_pT65-30_cent0-100.pdf");

  return;
}
