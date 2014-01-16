void calc_efficiency(std::string fname = "Psi2s_Histos_MC_PbPb.root", bool applyWeight=true, bool isHI=true, bool isJpsi=false)
{
  TFile *inf = new TFile(fname.c_str(),"READ");

  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Divide(2,1);
  c1->GetPad(1)->SetLogy();
  c1->GetPad(2)->SetLogy();
  c1->cd(1);
  // y<2.4
  TCut gen = "(abs(Gen_QQ_4mom.Rapidity())>=0.0&&abs(Gen_QQ_4mom.Rapidity())<2.4&&Gen_QQ_mupl_acc==1&&Gen_QQ_mumi_acc==1)";
  TCut rec = "(goodEvent&&abs(Reco_QQ_4mom.Rapidity())>=0.0&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Reco_QQ_mupl_acc==1&&Reco_QQ_mumi_acc==1&&Reco_QQ_mupl_cuts==1&&Reco_QQ_mumi_cuts==1)";
  TCut mass = "(Reco_QQ_4mom.M()>3.48&&Reco_QQ_4mom.M()<3.88)";
  if (isJpsi)
    mass = "(Reco_QQ_4mom.M()>2.95&&Reco_QQ_4mom.M()<3.25)";
  TCut trg = "((HLTriggers&1)==1&&(Reco_QQ_trig&1)==1)";
  TCut weight = "1";
  TCut cent = "Centrality>=0&&Centrality<40";
  if (!isHI)
    cent="1";

  // y<1.6
  TCut gen_y016 = "(abs(Gen_QQ_4mom.Rapidity())<1.6)";
  TCut rec_y016 = "(abs(Reco_QQ_4mom.Rapidity())<1.6)";

  // 1.6<y<2.4
  TCut gen_y1624 = "(abs(Gen_QQ_4mom.Rapidity())>=1.6)";
  TCut rec_y1624 = "(abs(Reco_QQ_4mom.Rapidity())>=1.6)";

  if (applyWeight)
    weight = "eventWeight";

  myTree->Draw("Gen_QQ_4mom.M()>>hMAcc(80,2.6,4.2)",(gen+cent)*weight,"e");
  myTree->Draw("Reco_QQ_4mom.M()>>hMRec(80,2.6,4.2)",(gen+cent+rec)*weight,"e");
  myTree->Draw("Reco_QQ_4mom.M()>>hMTrg(80,2.6,4.2)",(gen+cent+rec+trg)*weight,"e");

  myTree->Draw("Gen_QQ_4mom.Pt()>>hPtAcc(60,0,60)",(gen+cent)*weight,"e");
  myTree->Draw("Reco_QQ_4mom.Pt()>>hPtRec(60,0,60)",(gen+cent+mass+rec)*weight,"e");
  myTree->Draw("Reco_QQ_4mom.Pt()>>hPtTrg(60,0,60)",(gen+cent+mass+rec+trg)*weight,"e");

  myTree->Draw("Gen_QQ_4mom.Pt()>>hPtAcc_y016(60,0,60)",(gen+cent+gen_y016)*weight,"e");
  myTree->Draw("Reco_QQ_4mom.Pt()>>hPtRec_y016(60,0,60)",(gen+cent+mass+rec+gen_y016+rec_y016)*weight,"e");
  myTree->Draw("Reco_QQ_4mom.Pt()>>hPtTrg_y016(60,0,60)",(gen+cent+mass+rec+trg+gen_y016+rec_y016)*weight,"e");

  myTree->Draw("Gen_QQ_4mom.Pt()>>hPtAcc_y1624(60,0,60)",(gen+cent+gen_y1624)*weight,"e");
  myTree->Draw("Reco_QQ_4mom.Pt()>>hPtRec_y1624(60,0,60)",(gen+cent+mass+rec+gen_y1624+rec_y1624)*weight,"e");
  myTree->Draw("Reco_QQ_4mom.Pt()>>hPtTrg_y1624(60,0,60)",(gen+cent+mass+rec+trg+gen_y1624+rec_y1624)*weight,"e");

  hMRec->SetMarkerColor(kBlue);
  hMTrg->SetMarkerStyle(24);
  hMTrg->SetMarkerColor(kRed);

  hPtRec->SetMarkerColor(kBlue);
  hPtTrg->SetMarkerStyle(24);
  hPtTrg->SetMarkerColor(kRed);

  hMAcc->Draw();
  hMRec->Draw("same");
  hMTrg->Draw("same");

  double min = hMTrg->GetMaximum()/10000.0;
  hMAcc->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c^{2})");
  hMAcc->GetYaxis()->SetTitle("counts");
  hMAcc->SetMinimum(min);

  c1->cd(2);
  hPtAcc->Draw();
  hPtRec->Draw("same");
  hPtTrg->Draw("same");

  min = hPtTrg->GetMaximum()/10000.0;
  hPtAcc->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hPtAcc->GetYaxis()->SetTitle("counts");
  hPtAcc->SetMinimum(min);
  if (isJpsi)
    hPtAcc->GetXaxis()->SetRangeUser(0.0,30.0-0.001);

  
  int minPt = hPtAcc->GetXaxis()->FindBin(6.5+0.001);
  int maxPt = hPtAcc->GetXaxis()->FindBin(30.0-0.001);
  double errorRec, errorTrg, errorGen;
  cout << "|y|<2.4 && 6.5<pT<30.0 && 0-100%:   Gen: " << hPtAcc->IntegralAndError(minPt,maxPt,errorGen) << " +/- " << errorGen
       << " Reco: " << hPtRec->IntegralAndError(minPt,maxPt,errorRec) << " +/- " << errorRec
       << " eff: " << hPtRec->Integral(minPt,maxPt)/hPtAcc->Integral(minPt,maxPt) << "\t"
       << " Trg: " << hPtTrg->IntegralAndError(minPt,maxPt,errorTrg) << " +/- " << errorTrg
       << " eff: " << hPtTrg->Integral(minPt,maxPt)/hPtAcc->Integral(minPt,maxPt) << "+/-" << sqrt(pow(errorTrg/hPtTrg->Integral(minPt,maxPt),2)+pow(errorGen/hPtAcc->Integral(minPt,maxPt),2))*hPtTrg->Integral(minPt,maxPt)/hPtAcc->Integral(minPt,maxPt) << endl;
  cout << "|y|<1.6 && 6.5<pT<30.0 && 0-100%:   Gen: " << hPtAcc_y016->Integral(minPt,maxPt)
       << " Reco: " << hPtRec_y016->Integral(minPt,maxPt)
       << " eff: " << hPtRec_y016->Integral(minPt,maxPt)/hPtAcc_y016->Integral(minPt,maxPt) << "\t"
       << " Trg: " << hPtTrg_y016->Integral(minPt,maxPt)
       << " eff: " << hPtTrg_y016->Integral(minPt,maxPt)/hPtAcc_y016->Integral(minPt,maxPt) << endl;
  cout << "1.6<|y|<2.4 && 6.5<pT<30.0 && 0-100%:   Gen: " << hPtAcc_y1624->Integral(minPt,maxPt)
       << " Reco: " << hPtRec_y1624->Integral(minPt,maxPt)
       << " eff: " << hPtRec_y1624->Integral(minPt,maxPt)/hPtAcc_y1624->Integral(minPt,maxPt) << "\t"
       << " Trg: " << hPtTrg_y1624->Integral(minPt,maxPt)
       << " eff: " << hPtTrg_y1624->Integral(minPt,maxPt)/hPtAcc_y1624->Integral(minPt,maxPt) << endl;
  minPt = hPtAcc->GetXaxis()->FindBin(3.0+0.001);
  maxPt = hPtAcc->GetXaxis()->FindBin(30.0-0.001);
  cout << "1.6<|y|<2.4 && 3.0<pT<30.0 && 0-100%:   Gen: " << hPtAcc_y1624->Integral(minPt,maxPt)
       << " Reco: " << hPtRec_y1624->Integral(minPt,maxPt)
       << " eff: " << hPtRec_y1624->Integral(minPt,maxPt)/hPtAcc_y1624->Integral(minPt,maxPt) << "\t"
       << " Trg: " << hPtTrg_y1624->Integral(minPt,maxPt)
       << " eff: " << hPtTrg_y1624->Integral(minPt,maxPt)/hPtAcc_y1624->Integral(minPt,maxPt) << endl;
  minPt = hPtAcc->GetXaxis()->FindBin(3.0+0.001);
  maxPt = hPtAcc->GetXaxis()->FindBin(6.5-0.001);
  cout << "1.6<|y|<2.4 && 3.0<pT<6.5 && 0-100%:   Gen: " << hPtAcc_y1624->Integral(minPt,maxPt)
       << " Reco: " << hPtRec_y1624->Integral(minPt,maxPt)
       << " eff: " << hPtRec_y1624->Integral(minPt,maxPt)/hPtAcc_y1624->Integral(minPt,maxPt) << "\t"
       << " Trg: " << hPtTrg_y1624->Integral(minPt,maxPt)
       << " eff: " << hPtTrg_y1624->Integral(minPt,maxPt)/hPtAcc_y1624->Integral(minPt,maxPt) << endl;


  //  c1->cd(3);
  TCanvas *c2 = new TCanvas("c2","c2");
  // TGraphAsymmErrors *effRec = new TGraphAsymmErrors();
  // effRec->SetName("effRec");
  // effRec->SetMarkerColor(kBlue);
  // effRec->BayesDivide(hPtRec,hPtAcc);
  TH1F* effRec = (TH1F*) hPtRec->Clone("effRec");
  effRec->Divide(hPtRec,hPtAcc,1,1,"B");
  effRec->Draw("");
  effRec->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  effRec->GetYaxis()->SetTitle("#varepsilon");
  effRec->GetXaxis()->SetRangeUser(0,29.99);
  effRec->GetYaxis()->SetRangeUser(0,1);

  // TGraphAsymmErrors *effTrg = new TGraphAsymmErrors();
  // effTrg->SetName("effTrg");
  // effTrg->SetMarkerStyle(24);
  // effTrg->SetMarkerColor(kRed);
  // effTrg->BayesDivide(hPtTrg,hPtAcc);
  TH1F* effTrg = (TH1F*) hPtTrg->Clone("effTrg");
  effTrg->Divide(hPtTrg,hPtAcc,1,1,"B");
  effTrg->Draw("same");

  

  TCanvas *c3 = new TCanvas("c3","c3");
  TCut gen_pt330 = "(Gen_QQ_4mom.Pt()>3.0&&Gen_QQ_4mom.Pt()<30)";
  TCut rec_pt330 = "(Reco_QQ_4mom.Pt()>3.0&&Reco_QQ_4mom.Pt()<30)";


  myTree->Draw("Gen_QQ_mupl_4mom.Pt()>>hMuPlPtAcc(400,0,20)",(gen+cent+gen_y1624+gen_pt330)*weight,"e");
  myTree->Draw("Reco_QQ_mupl_4mom.Pt()>>hMuPlPtRec(400,0,20)",(gen+cent+mass+rec+gen_y1624+gen_pt330+rec_y1624+rec_pt330)*weight,"e");
  myTree->Draw("Reco_QQ_mupl_4mom.Pt()>>hMuPlPtTrg(400,0,20)",(gen+cent+mass+rec+trg+gen_y1624+gen_pt330+rec_y1624+rec_pt330)*weight,"e");

  myTree->Draw("Gen_QQ_mumi_4mom.Pt()>>hMuMiPtAcc(400,0,20)",(gen+cent+gen_y1624+gen_pt330)*weight,"e");
  myTree->Draw("Reco_QQ_mumi_4mom.Pt()>>hMuMiPtRec(400,0,20)",(gen+cent+mass+rec+gen_y1624+gen_pt330+rec_y1624+rec_pt330)*weight,"e");
  myTree->Draw("Reco_QQ_mumi_4mom.Pt()>>hMuMiPtTrg(400,0,20)",(gen+cent+mass+rec+trg+gen_y1624+gen_pt330+rec_y1624+rec_pt330)*weight,"e");

  TH1F *hMuPtAcc = new TH1F("hMuPtAcc","hMuPtAcc",400,0,20);
  TH1F *hMuPtRec = new TH1F("hMuPtRec","hMuPtRec",400,0,20);
  TH1F *hMuPtTrg = new TH1F("hMuPtTrg","hMuPtTrg",400,0,20);
  hMuPtAcc->Add(hMuPlPtAcc,hMuMiPtAcc);
  hMuPtRec->Add(hMuPlPtRec,hMuMiPtRec);
  hMuPtTrg->Add(hMuPlPtTrg,hMuMiPtTrg);
  hMuPtRec->SetMarkerColor(kBlue);
  hMuPtTrg->SetMarkerColor(kRed);
  hMuPtTrg->SetMarkerStyle(24);

  hMuPtAcc->Draw();
  hMuPtRec->Draw("same");
  //  hMuPtTrg->Draw();

  // TFile *outf;
  // if (isJpsi)
  //   outf = new TFile("SingleMuPt_fwd_JpsiMC.root","CREATE");
  // else
  //   outf = new TFile("SingleMuPt_fwd_Psi2SMC.root","CREATE");

  // hMuPtAcc->Write();
  // hMuPtRec->Write();
  // hMuPtTrg->Write();
  // outf->Close();

  return;
}
