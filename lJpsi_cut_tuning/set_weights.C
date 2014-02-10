void set_weights()
{
  // prompt
  TFile *f0 = new TFile("/data/CMS/MC/PbPb/jpsiMuMu_JpsiPt03_Histos_cmssw445p5_RegIt.root","UPDATE");
  myTree->SetWeight(1.0/117512.0);
  myTree->AutoSave();
  f0->Close();
  TFile *f1 = new TFile("/data/CMS/MC/PbPb/jpsiMuMu_JpsiPt36_Histos_cmssw445p5_RegIt.root","UPDATE");
  myTree->SetWeight(0.96075/124802.0);
  myTree->AutoSave();
  f1->Close();
  TFile *f2 = new TFile("/data/CMS/MC/PbPb/jpsiMuMu_JpsiPt69_Histos_cmssw445p5_RegIt.root","UPDATE");
  myTree->SetWeight(0.2136/172176.0);
  myTree->AutoSave();
  f2->Close();
  TFile *f3 = new TFile("/data/CMS/MC/PbPb/jpsiMuMu_JpsiPt912_Histos_cmssw445p5_RegIt.root","UPDATE");
  myTree->SetWeight(0.04385/167336.0);
  myTree->AutoSave();
  f3->Close();
  TFile *f4 = new TFile("/data/CMS/MC/PbPb/jpsiMuMu_JpsiPt1215_Histos_cmssw445p5_RegIt.root","UPDATE");
  myTree->SetWeight(0.0112585/143944.0);
  myTree->AutoSave();
  f4->Close();
  TFile *f5 = new TFile("/data/CMS/MC/PbPb/jpsiMuMu_JpsiPt1530_Histos_cmssw445p5_RegIt.root","UPDATE");
  myTree->SetWeight(0.00579/166960.0);
  myTree->AutoSave();
  f5->Close();
  /*
  // non-prompt
  TFile *f0 = new TFile("/data/CMS/MC/PbPb/bJpsiMuMu_JpsiPt03_Histos_cmssw445p5_RegIt.root","UPDATE");
  myTree->SetWeight(1.0/171458.0);
  myTree->AutoSave();
  f0->Close();
  TFile *f1 = new TFile("/data/CMS/MC/PbPb/bJpsiMuMu_JpsiPt36_Histos_cmssw445p5_RegIt.root","UPDATE");
  myTree->SetWeight(2.25*0.623/172208.0);
  myTree->AutoSave();
  f1->Close();
  TFile *f2 = new TFile("/data/CMS/MC/PbPb/bJpsiMuMu_JpsiPt69_Histos_cmssw445p5_RegIt.root","UPDATE");
  myTree->SetWeight(2.55*0.245/171178.0);
  myTree->AutoSave();
  f2->Close();
  TFile *f3 = new TFile("/data/CMS/MC/PbPb/bJpsiMuMu_JpsiPt912_Histos_cmssw445p5_RegIt.root","UPDATE");
  myTree->SetWeight(2.45*0.080/172208.0);
  myTree->AutoSave();
  f3->Close();
  TFile *f4 = new TFile("/data/CMS/MC/PbPb/bJpsiMuMu_JpsiPt1215_Histos_cmssw445p5_RegIt.root","UPDATE");
  myTree->SetWeight(2.0*0.0299/170670.0);
  myTree->AutoSave();
  f4->Close();
  TFile *f5 = new TFile("/data/CMS/MC/PbPb/bJpsiMuMu_JpsiPt1530_Histos_cmssw445p5_RegIt.root","UPDATE");
  myTree->SetWeight(3.35*0.0109/169302.0);
  myTree->AutoSave();
  f5->Close();
  */
  return;
}
