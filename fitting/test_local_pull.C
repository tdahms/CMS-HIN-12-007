void test_local_pull(string infname="20140423_SimFits_M2242_DblMu0_AllCent_WithSyst_final/fracLogCBG_PbPbpol3_HI020_pol2_HI2040_pol1_HI40100_pppol1_rap16-24_pT3-30_centMult_plots.root")
{
  gStyle->SetOptStat("emr");

  TFile *inf = new TFile(infname.c_str(),"READ");
  inf->cd("pp");
  RooHist *pull = (RooHist*) gROOT->FindObject("hpullhist_pp");

  TH1F *hpull = new TH1F("hpull","hpull",100,-10,10);
  double x, y;
  double chi2 = 0;
  int n = 0;
  for (int i=0;i<pull->GetN();++i) {
    pull->GetPoint(i,x,y);
    if (x>3.5&&x<3.9) {
      hpull->Fill(y);
      chi2 += pow(y,2);
      n++;
    }
  }

  TCanvas *c1 = new TCanvas("c1","c1");
  c1->cd();
  hpull->Draw();

  cout << "Local chi2 = " << chi2 << " with dof = " << n-1 << " (p-value = " << TMath::Prob(chi2, n-1) << ")" << endl;

  return;
}
