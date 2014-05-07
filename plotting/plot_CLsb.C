void plot_CLsb(std::string infname="../fitting/20140424_CFscan_WithSyst_Good/HI020/Freq_Cls+b_grid_ts2_fracLogCBG_PbPbpol3_HI020_pol2_HI2040_pol1_HI40100_pppol1_rap16-24_pT3-30_centMult_Workspace.root",double min=0, double max=4)
{
  TFile *inf = new TFile(infname.c_str(),"READ");
  
  const char * resultName = result_doubleRatio_HI020->GetName();
  string range = "1.6<|y|<2.4 and 3<p_{T}<30 GeV/c";
  TString plotTitle = TString::Format("Frequentist CL Scan for %s",range.c_str());
  RooStats::HypoTestInverterPlot *plot = new RooStats::HypoTestInverterPlot("HTI_Result_Plot",plotTitle,result_doubleRatio_HI020);
  
  TCanvas * c1 = new TCanvas("c1","c1"); 
  c1->SetBottomMargin(0.15);
  c1->SetLogy(false);
  TF1 *f1 = new TF1("f1","0.05",min,max);
  f1->SetLineWidth(1);
  f1->SetLineColor(kRed);
  f1->Draw();
  f1->GetYaxis()->SetRangeUser(0,1);
  f1->SetTitle(plotTitle);
  f1->GetXaxis()->SetTitle("[N_{#psi(2S)}/N_{J/#psi}]_{PbPb}^{0-20%} / [N_{#psi(2S)}/N_{J/#psi}]_{pp}");
  f1->GetXaxis()->CenterTitle();
  f1->GetXaxis()->SetTitleOffset(1.2);
  f1->GetYaxis()->SetTitle("p value");

  TF1 *f2 = new TF1("f2","0.32",min,max);
  f2->SetLineWidth(2);
  f2->SetLineStyle(7);
  f2->SetLineColor(kBlue);
  f2->Draw("same");

  TGraphErrors *gr = plot->MakePlot("CLs+b");//Draw("OBS same");  // plot all and Clb
  gr->SetName("gr");
  gr->SetMarkerSize(1.0);
  gr->SetMarkerColor(kRed);
  //  gr->SetLineWidth(2);
  gr->Draw("PL");

  TLegend *leg;
  leg = new TLegend(0.6,0.83,0.93,0.93);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetMargin(0.15);
  //  leg->SetTextSize(0.04);
  leg->SetTextSize(0.035);
  leg->AddEntry(gr,"Observed CLs+b","P");
  leg->Draw();

  return;
}
