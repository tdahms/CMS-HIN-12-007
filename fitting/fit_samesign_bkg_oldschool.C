#include <iostream>
#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TBox.h"

using namespace std;

void fit_samesign_bkg_oldschool(string infname="mid_pbpb.root")
{
  gStyle->SetOptFit(0);

  TFile *inf = new TFile(infname.c_str(),"READ");
  TH1F *h1 = (TH1F*) inf->Get("h1");
  h1->GetXaxis()->SetTitle("m_{#mu^{#pm}#mu^{#pm}} (GeV/c^{2})");
  h1->GetXaxis()->CenterTitle(1);

  TF1 *f1 = new TF1("f1","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)",1.8,5.0);
  f1->SetParNames("a","b","c","d");
  f1->SetParameters(1,1,1,0);
  f1->SetNpx(200);
  //f1->FixParameter(3,0);
  f1->SetLineWidth(3);
  f1->SetLineColor(kBlue);
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->SetBottomMargin(0.14);

  h1->Draw("");
  //plotting
  TFitResultPtr s1 = h1->Fit(f1,"RMEILS","");
  double R1850 = f1->Integral(3.52,3.84)/0.02;///f1->Integral(1.8,5.0)*h1->Integral(1,160);
  double errR1850 = f1->IntegralError(3.52,3.84,s1->GetParams(), s1->GetCovarianceMatrix().GetMatrixArray())/0.02;
  cout << "M1850: " << R1850 << " +/- " << errR1850 << endl;

  f1->SetLineColor(kRed);
  f1->SetLineStyle(7);
  TFitResultPtr s2 = h1->Fit(f1,"MEILS+","",2.0,4.5);
  double R2045 = f1->Integral(3.52,3.84)/0.02;///f1->Integral(1.8,5.0)*h1->Integral(1,160);
  double errR2045 = f1->IntegralError(3.52,3.84,s2->GetParams(), s2->GetCovarianceMatrix().GetMatrixArray())/0.02;
  cout << "M2045: " << R2045 << " +/- " << errR2045 << endl;

  f1->SetLineColor(kGreen+2);
  f1->SetLineStyle(10);
  TFitResultPtr s3 = h1->Fit(f1,"MEILS+","",2.2,4.2);
  double R2242 = f1->Integral(3.52,3.84)/0.02;///f1->Integral(1.8,5.0)*h1->Integral(1,160);
  double errR2242 = f1->IntegralError(3.52,3.84,s3->GetParams(), s3->GetCovarianceMatrix().GetMatrixArray())/0.02;
  cout << "M2242: " << R2242 << " +/- " << errR2242 << endl;

  double err = 0.0;
  int min = h1->GetXaxis()->FindBin(3.52001);
  int max = h1->GetXaxis()->FindBin(3.83999);
  double Rdata = h1->IntegralAndError(min,max,err);
  cout << "data: " << Rdata << " +/- " << err << " (" << err/Rdata*100 << "%)" << endl;

  cout << R1850 << " +/- " << errR1850 << endl;
  cout << R2045 << " +/- " << errR2045 << endl;
  cout << R2242 << " +/- " << errR2242 << endl;
  
  cout << 1-(R1850/Rdata) << endl;
  cout << 1-(R2045/Rdata) << endl;
  cout << 1-(R2242/Rdata) << endl;

  double maxY = h1->GetMaximum();
  maxY+=1.1*sqrt(maxY);
  TBox *b1 = new TBox(3.52,0,3.84,maxY);
  b1->SetFillColor(kYellow);
  //  b1->Draw();

  b1->Draw();
  h1->Draw("same");
  gPad->RedrawAxis();

  // TLine *l1 = new TLine(3.52,0,3.52,100);
  // TLine *l2 = new TLine(3.84,0,3.84,100);

  // s1.Get()->Print("v");
  // s2.Get()->Print("v");
  // s3.Get()->Print("v");

  // LLR test
  /*
  TF1 *f2 = new TF1("f2","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)",1.8,5.0);
  f2->SetParNames("a","b","c","d","e","f");
  f2->SetParameters(1,1,1,0,0,0);

  for (int i=5;i>=0;--i) {
    TFitResultPtr s1 = h1->Fit(f2,"RMEIQSN");
    TFitResultPtr s2 = h1->Fit(f2,"MEIQSN","",2.0,4.5);
    TFitResultPtr s3 = h1->Fit(f2,"MEIQSN","",2.2,4.2);
    if (!s1->IsValid())cout << s1->Status() << endl;
    if (!s2->IsValid())cout << s2->Status() << endl;
    if (!s3->IsValid())cout << s3->Status() << endl;
    cout << i <<": " << s1->Chi2() << "\t" << s2->Chi2() << "\t" << s3->Chi2() << endl;
    cout << i <<": " << s1->Ndf() << "\t" << s2->Ndf() << "\t" << s3->Ndf() << endl;
    f2->FixParameter(i,0);
  }
  */
  return;
}
