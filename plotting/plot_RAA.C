void plot_RAA()
{

  // Jpsi, psi(2S), Y(1S), Y(2S), Y(3S)
  double xRadius[] = {0.5, 0.9, 0.28, 0.56, 0.78};
  double x_err[] =   {0.0, 0.0, 0.0, 0.0, 0.0};
  double x_sys[] =   {0.025, 0.025, 0.025, 0.025, 0.025};

  double RAA[] = {0.30, 0.11, 0.53, 0.13, 0.1};
  double RAA_err[] = {0.03, 0.03, 0.07, 0.04, 0.0};
  double RAA_sys[] = {0.02, 0.18, 0.07, 0.02, 0.0};


  TGraphErrors *gJpsi = new TGraphErrors(1, xRadius, RAA, x_err, RAA_err);
  TGraphErrors *gPsiP = new TGraphErrors(2, xRadius, RAA, x_err, RAA_err);
  TGraphErrors *gUps1 = new TGraphErrors(3, xRadius, RAA, x_err, RAA_err);
  TGraphErrors *gUps2 = new TGraphErrors(4, xRadius, RAA, x_err, RAA_err);
  //  TGraphErrors *gUps3 = new TGraphErrors(5, xRadius, RAA, x_err, RAA_err);

  gJpsi->SetName("gJpsi");
  gPsiP->SetName("gPsiP");
  gUps1->SetName("gUps1");
  gUps2->SetName("gUps2");
  //  gUps3->SetName("gUps3");

  TGraphErrors *gJpsi_sys = new TGraphErrors(1, xRadius, RAA, x_sys, RAA_sys);
  TGraphErrors *gPsiP_sys = new TGraphErrors(2, xRadius, RAA, x_sys, RAA_sys);
  TGraphErrors *gUps1_sys = new TGraphErrors(3, xRadius, RAA, x_sys, RAA_sys);
  TGraphErrors *gUps2_sys = new TGraphErrors(4, xRadius, RAA, x_sys, RAA_sys);
  //  TGraphErrors *gUps3_sys = new TGraphErrors(5, xRadius, RAA, x_sys, RAA_sys);

  gJpsi_sys->SetName("gJpsi_sys");
  gPsiP_sys->SetName("gPsiP_sys");
  gUps1_sys->SetName("gUps1_sys");
  gUps2_sys->SetName("gUps2_sys");
  //  gUps3_sys->SetName("gUps3_sys");


  TGraphErrors *gJpsi_p = new TGraphErrors(1, xRadius, RAA, x_err, RAA_err);
  TGraphErrors *gPsiP_p = new TGraphErrors(2, xRadius, RAA, x_err, RAA_err);
  TGraphErrors *gUps1_p = new TGraphErrors(3, xRadius, RAA, x_err, RAA_err);
  TGraphErrors *gUps2_p = new TGraphErrors(4, xRadius, RAA, x_err, RAA_err);
  //  TGraphErrors *gUps3_p = new TGraphErrors(5, xRadius, RAA, x_err, RAA_err);

  gJpsi_p->SetName("gJpsi_p");
  gPsiP_p->SetName("gPsiP_p");
  gUps1_p->SetName("gUps1_p");
  gUps2_p->SetName("gUps2_p");
  //  gUps3_p->SetName("gUps3_p");



  gPsiP->RemovePoint(0);

  gUps1->RemovePoint(0);
  gUps1->RemovePoint(0);

  gUps2->RemovePoint(0);
  gUps2->RemovePoint(0);
  gUps2->RemovePoint(0);

  // gUps3->RemovePoint(0);
  // gUps3->RemovePoint(0);
  // gUps3->RemovePoint(0);
  // gUps3->RemovePoint(0);

  gPsiP_sys->RemovePoint(0);

  gUps1_sys->RemovePoint(0);
  gUps1_sys->RemovePoint(0);

  gUps2_sys->RemovePoint(0);
  gUps2_sys->RemovePoint(0);
  gUps2_sys->RemovePoint(0);

  // gUps3_sys->RemovePoint(0);
  // gUps3_sys->RemovePoint(0);
  // gUps3_sys->RemovePoint(0);
  // gUps3_sys->RemovePoint(0);

  gPsiP_p->RemovePoint(0);

  gUps1_p->RemovePoint(0);
  gUps1_p->RemovePoint(0);

  gUps2_p->RemovePoint(0);
  gUps2_p->RemovePoint(0);
  gUps2_p->RemovePoint(0);

  // gUps3_p->RemovePoint(0);
  // gUps3_p->RemovePoint(0);
  // gUps3_p->RemovePoint(0);
  // gUps3_p->RemovePoint(0);

  gJpsi_sys->SetFillColor(kRed-9);
  gPsiP_sys->SetFillColor(kOrange-9);
  gUps1_sys->SetFillColor(kGreen-9);
  gUps2_sys->SetFillColor(kBlue-9);
  //  gUps3_sys->SetFillColor(kGray);


  gJpsi->SetMarkerColor(kRed+2);
  gPsiP->SetMarkerColor(kOrange+2);
  gUps1->SetMarkerColor(kGreen+2);
  gUps2->SetMarkerColor(kBlue+2);

  gJpsi->SetMarkerSize(1.2);
  gPsiP->SetMarkerSize(1.2);
  gUps1->SetMarkerSize(1.2);
  gUps2->SetMarkerSize(1.2);

  gJpsi_p->SetMarkerSize(1.2);
  gPsiP_p->SetMarkerSize(1.2);
  gUps1_p->SetMarkerSize(1.2);
  gUps2_p->SetMarkerSize(1.2);

  gJpsi_p->SetMarkerStyle(24);
  gPsiP_p->SetMarkerStyle(24);
  gUps1_p->SetMarkerStyle(24);
  gUps2_p->SetMarkerStyle(24);


  TF1 *f1 = new TF1("f1","1",0,1);

  f1->SetLineWidth(1);
  f1->GetXaxis()->SetTitle("r_{0} (fm)");
  f1->GetYaxis()->SetTitle("R_{AA}");
  f1->GetYaxis()->SetTitle("R_{AA}");
  f1->GetYaxis()->SetRangeUser(0.0,1.2);
  f1->GetXaxis()->CenterTitle(kTRUE);

  TCanvas *c1 = new TCanvas("c1","c1");
  f1->Draw();
  gJpsi_sys->Draw("2");
  gPsiP_sys->Draw("2");
  gUps1_sys->Draw("2");
  gUps2_sys->Draw("2");

  gJpsi->Draw("P");
  gPsiP->Draw("P");
  gUps1->Draw("P");
  gUps2->Draw("P");


  gJpsi_p->Draw("P");
  gPsiP_p->Draw("P");
  gUps1_p->Draw("P");
  gUps2_p->Draw("P");


  TLegend *leg1;
  leg1 = new TLegend(0.15,0.625,0.72,0.875);

  leg1->SetFillStyle(0);
  leg1->SetFillColor(0);
  leg1->SetBorderSize(0);
  leg1->SetMargin(0.15);
  leg1->SetTextSize(0.035);
  leg1->SetHeader("PbPb  #sqrt{s_{NN}} = 2.76 TeV");
  leg1->AddEntry(gJpsi,"prompt J/#psi (6.5 < p_{T} < 30 GeV/c, |y| < 2.4)","P");
  leg1->AddEntry(gPsiP,"incl. #psi(2S) (6.5 < p_{T} < 30 GeV/c, |y| < 1.6)","P");
  leg1->AddEntry(gUps1,"#varUpsilon(1S) (|y| < 2.4)","P");
  leg1->AddEntry(gUps2,"#varUpsilon(2S) (|y| < 2.4)","P");
  //  leg1->AddEntry(gUps3,"#varUpsilon(3S) (|y| < 2.4)","P");
  leg1->Draw();


  TLatex *pre = new TLatex(0.05,1.1,"CMS Preliminary");
  pre->SetTextFont(42);
  pre->SetTextSize(0.05);

  pre->Draw();

  gPad->RedrawAxis();


  return;
}
