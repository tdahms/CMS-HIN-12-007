void plot_CLsb(std::string infname="../fitting/Good_CFscan_WithSyst/HI020/Freq_Cls+b_grid_ts2_fracLogCBG_PbPbpol3_HI020_pol2_HI2040_pol1_HI40100_pppol1_rap16-24_pT3-30_centMult_Workspace.root")
{
  TFile *inf = new TFile(infname.c_str(),"READ");
  
  const char * resultName = result_doubleRatio_HI020->GetName();
  string range = "1.6<|y|<2.4 and 3<p_{T}<30 GeV/c";
  TString plotTitle = TString::Format("Frequentist CL Scan for %s",range.c_str());
  RooStats::HypoTestInverterPlot *plot = new RooStats::HypoTestInverterPlot("HTI_Result_Plot",plotTitle,result_doubleRatio_HI020);
  
  TCanvas * c1 = new TCanvas("c1","c1"); 
  c1->SetLogy(false);
  TF1 *f1 = new TF1("f1","0.05",0,4.0);
  f1->SetLineWidth(1);
  f1->SetLineColor(kRed);
  f1->Draw();
  f1->GetYaxis()->SetRangeUser(0,1);
  f1->SetTitle(plotTitle);
  f1->GetXaxis()->SetTitle("[ #psi (2S) #/J/#psi ]_{PbPb} #/[ #psi (2S) #/J/#psi ]_{pp}(0-20%)");
  f1->GetYaxis()->SetTitle("p value");

  TF1 *f2 = new TF1("f2","0.32",0,4.0);
  f2->SetLineWidth(2);
  f2->SetLineStyle(7);
  f2->SetLineColor(kBlue);
  f2->Draw("same");

  plot->Draw("OBS same");  // plot all and Clb
  return;
}
