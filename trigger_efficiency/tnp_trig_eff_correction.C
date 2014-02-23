#include "TFile.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TGraphAsymmErrors.h"


void tnp_trig_eff_correction()
{
  //  TGraphAsymmErrors *grEffRD;
  TGraphAsymmErrors *grEffMC;


  TGraphAsymmErrors *grEffRD = new TGraphAsymmErrors(8);
  grEffRD->SetPoint(0,2.121837,0.5283981);
  grEffRD->SetPointError(0,0.6218373,0.3781627,0.03749582,0.03747286);
  grEffRD->SetPoint(1,2.74728,0.7798603);
  grEffRD->SetPointError(1,0.2472804,0.2527196,0.02377983,0.02311939);
  grEffRD->SetPoint(2,3.795761,0.89572);
  grEffRD->SetPointError(2,0.7957612,0.7042388,0.00602868,0.005937595);
  grEffRD->SetPoint(3,5.124785,0.9653656);
  grEffRD->SetPointError(3,0.624785,0.875215,0.004552644,0.004252621);
  grEffRD->SetPoint(4,7.137906,0.9763212);
  grEffRD->SetPointError(4,1.137906,1.862094,0.004288511,0.003897451);
  grEffRD->SetPoint(5,10.11656,0.985752);
  grEffRD->SetPointError(5,1.116558,1.883442,0.007086854,0.00549465);
  grEffRD->SetPoint(6,13.30238,0.9940162);
  grEffRD->SetPointError(6,1.302379,1.697621,0.01671768,0.005983758);
  grEffRD->SetPoint(7,16.91028,0.9754653);
  grEffRD->SetPointError(7,1.910281,3.089719,0.0295633,0.01524883);


  // TFile *tnpRD = new TFile("MC/Jpsi_RD_2012_Trg_L1_Eff_cbExpo.root","READ");
  // axis = (TGraphAsymmErrors*) Trg_pt_All->Clone("axis");
  // tnpRD->Close();
  grEffRD->SetMarkerColor(kBlue);
  grEffRD->SetMarkerStyle(21);

  TGraphAsymmErrors *grEffMC = new TGraphAsymmErrors(8);
  grEffMC->SetPoint(0,2.174589,0.4886081);
  grEffMC->SetPointError(0,0.6745893,0.3254107,0.009670632,0.009548092);
  grEffMC->SetPoint(1,2.765836,0.6395594);
  grEffMC->SetPointError(1,0.2658358,0.2341642,0.00878076,0.008728074);
  grEffMC->SetPoint(2,3.809204,0.8714975);
  grEffMC->SetPointError(2,0.8092038,0.6907962,0.002470183,0.002445725);
  grEffMC->SetPoint(3,5.100813,0.9544173);
  grEffMC->SetPointError(3,0.6008126,0.8991874,0.002139408,0.002069145);
  grEffMC->SetPoint(4,7.101783,0.959291);
  grEffMC->SetPointError(4,1.101783,1.898217,0.003232583,0.003085069);
  grEffMC->SetPoint(5,10.11028,0.9572457);
  grEffMC->SetPointError(5,1.110284,1.889716,0.007067502,0.006395253);
  grEffMC->SetPoint(6,13.20572,0.9646924);
  grEffMC->SetPointError(6,1.205722,1.794278,0.01269083,0.01040044);
  grEffMC->SetPoint(7,17.06318,0.9769534);
  grEffMC->SetPointError(7,2.063182,2.936818,0.01891466,0.01287603);

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
  grEffRatio->GetFunction("f1")->Delete();

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


  return;
}
