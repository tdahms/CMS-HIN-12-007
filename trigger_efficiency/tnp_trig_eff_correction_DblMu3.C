#include "TFile.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TGraphAsymmErrors.h"


void tnp_trig_eff_correction_DblMu3()
{
  //  TGraphAsymmErrors *grEffRD;
  //  TGraphAsymmErrors *grEffMC;


  TGraphAsymmErrors *grEffRD = new TGraphAsymmErrors(8);
  grEffRD->SetPoint(0,2.121837,0.09001519);
  grEffRD->SetPointError(0,0.6218373,0.3781627,0.0179337,0.01957792);
  grEffRD->SetPoint(1,2.74728,0.2420874);
  grEffRD->SetPointError(1,0.2472804,0.2527196,0.02191146,0.02382281);
  grEffRD->SetPoint(2,3.795761,0.7201343);
  grEffRD->SetPointError(2,0.7957612,0.7042388,0.008815487,0.008733398);
  grEffRD->SetPoint(3,5.124785,0.9004093);
  grEffRD->SetPointError(3,0.624785,0.875215,0.007220636,0.01203388);
  grEffRD->SetPoint(4,7.137906,0.9583283);
  grEffRD->SetPointError(4,1.137906,1.862094,0.005572095,0.00520308);
  grEffRD->SetPoint(5,10.11656,0.9744146);
  grEffRD->SetPointError(5,1.116558,1.883442,0.010175,0.01530677);
  grEffRD->SetPoint(6,13.30238,0.972695);
  grEffRD->SetPointError(6,1.302379,1.697621,0.01795625,0.01170048);
  grEffRD->SetPoint(7,16.91028,0.9769509);
  grEffRD->SetPointError(7,1.910281,3.089719,0.02078673,0.01312097);


  TFile *tnpRD = new TFile("MC/Jpsi_RD_2012_Trg_L1_Eff_cbExpo.root","READ");
  axis = (TGraphAsymmErrors*) Trg_pt_All->Clone("axis");
  // tnpRD->Close();
  grEffRD->SetMarkerColor(kBlue);
  grEffRD->SetMarkerStyle(21);

  TGraphAsymmErrors *grEffMC = new TGraphAsymmErrors(8);
  grEffMC->SetPoint(0,2.174589,0.06302429);
  grEffMC->SetPointError(0,0.6745893,0.3254107,0.004583018,0.004804124);
  grEffMC->SetPoint(1,2.765836,0.2272458);
  grEffMC->SetPointError(1,0.2658358,0.2341642,0.007388634,0.007505177);
  grEffMC->SetPoint(2,3.809204,0.6989028);
  grEffMC->SetPointError(2,0.8092038,0.6907962,0.003390834,0.003373203);
  grEffMC->SetPoint(3,5.100813,0.9100785);
  grEffMC->SetPointError(3,0.6008126,0.8991874,0.002932325,0.08992151);
  grEffMC->SetPoint(4,7.101783,0.9396493);
  grEffMC->SetPointError(4,1.101783,1.898217,0.003849575,0.003786857);
  grEffMC->SetPoint(5,10.11028,0.9554703);
  grEffMC->SetPointError(5,1.110284,1.889716,0.007352275,0.00667199);
  grEffMC->SetPoint(6,13.20572,0.9714306);
  grEffMC->SetPointError(6,1.205722,1.794278,0.01131661,0.00935118);
  grEffMC->SetPoint(7,17.06318,0.9757057);
  grEffMC->SetPointError(7,2.063182,2.936818,0.0189386,0.01287194);

  // TFile *tnpMC = new TFile("MC/Jpsi_MC_2012_Trg_L1_Eff.root","READ");
  // grEffMC = (TGraphAsymmErrors*) Trg_pt_All->Clone("grEffMC");
  // tnpMC->Close();
  grEffMC->SetMarkerColor(kRed);
  grEffMC->SetMarkerStyle(20);

  // double effMC[8] = {0.48, 0,68, 0.87, 0.94, 0.96, 0.97, 0.96, 0.98}
  // double effRD[8] = {0.59, 0,76, 0.91, 0.97, 0.98, 0.99, 0.96, 0.98}

  TF1 *f1 = new TF1("f1","[0]*TMath::Erf((x-[1])/[2])",1.5,20);

  f1->SetParNames("eff0","x0","m");
  f1->SetParameters(1.0,0.5,1);
  f1->SetLineWidth(2);

  grEffRD->Fit(f1,"RME");
  grEffMC->Fit(f1,"RME");

  grEffRD->GetFunction("f1")->SetLineColor(kBlue);
  grEffMC->GetFunction("f1")->SetLineColor(kRed);


  TCanvas *c1 = new TCanvas("c1","c1");
  grEffRD->Draw("AP");
  grEffMC->Draw("P");

  grEffRD->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  grEffRD->GetYaxis()->SetTitle("Trigger Efficiency");


  TGraphAsymmErrors *grEffRatio = (TGraphAsymmErrors*) grEffMC->Clone("grEffRatio");

  double xMC = 0;
  double yMC = 0;
  double xRD = 0;
  double yRD = 0;

  for (int i=0; i<grEffRatio->GetN(); ++i) {
    grEffMC->GetPoint(i,xMC,yMC);
    grEffRD->GetPoint(i,xRD,yRD);
    grEffRatio->SetPoint(i,xMC,yRD/yMC);

    double errxl = grEffMC->GetErrorXlow(i);
    double errxh = grEffMC->GetErrorXhigh(i);

    double erryl = yRD/yMC*sqrt( pow(grEffMC->GetErrorYlow(i)/yMC,2) + pow(grEffRD->GetErrorYlow(i)/yRD,2) );
    double erryh = yRD/yMC*sqrt( pow(grEffMC->GetErrorYhigh(i)/yMC,2) + pow(grEffRD->GetErrorYhigh(i)/yRD,2) );

    grEffRatio->SetPointError(i, errxl, errxh, erryl, erryh);

    cout << "pt bin: " << xMC-errxl  << "--" << xMC+errxh << " weight: " << yRD/yMC << endl;

  }

  TF1 *f2 = new TF1("f2","[0]*TMath::Erf((x-[1])/[2])/TMath::Erf((x-[3])/[4])",1.5,20);
  f2->SetParameter(0, grEffRD->GetFunction("f1")->GetParameter(0)/grEffMC->GetFunction("f1")->GetParameter(0));
  f2->SetParameter(1, grEffRD->GetFunction("f1")->GetParameter(1));
  f2->SetParameter(2, grEffRD->GetFunction("f1")->GetParameter(2));
  f2->SetParameter(3, grEffMC->GetFunction("f1")->GetParameter(1));
  f2->SetParameter(4, grEffMC->GetFunction("f1")->GetParameter(2));

  TCanvas *c2 = new TCanvas("c2","c2");
  f2->Draw();
  grEffRatio->Draw("P");
  f2->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  f2->GetYaxis()->SetTitle("Data (fit) / MC (fit)");

  cout << "f2->SetParameters(" << f2->GetParameter(0) << ", " << f2->GetParameter(1) << ", " << f2->GetParameter(2) << ", " << f2->GetParameter(3) << ", " << f2->GetParameter(4) << ");" << endl;

  return;
}
