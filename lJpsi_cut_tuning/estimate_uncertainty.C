#include <iostream>
#include "TRandom3.h"
#include "TH1F.h"
#include "TCanvas.h"

void estimate_uncertainty(int N=1000)
{
  double a,b,c,d;
  TRandom3 rand;
  double mean = 0.0;
  double sigma = 0.25;
  double sigma1 = 0.22;

  TH1F *hA = new TH1F("hA","numerator 1;A;counts", 200,0,2);
  TH1F *hB = new TH1F("hB","denominator 1;B;counts", 200,0,2);
  TH1F *hC = new TH1F("hC","numerator 2;C;counts", 200,0,2);
  TH1F *hD = new TH1F("hD","denominator 2;D;counts", 200,0,2);
  TH1F *hR1 = new TH1F("hR1","ratio 1;single ratio;counts",200,0,2);
  TH1F *hR2 = new TH1F("hR2","ratio 2;single raito;counts",200,0,2);
  TH1F *hDR = new TH1F("hDR","double ratio;double ratio;counts",200,0,2);

  hA->Sumw2();
  hB->Sumw2();
  hC->Sumw2();
  hD->Sumw2();
  hR1->Sumw2();
  hR2->Sumw2();
  hDR->Sumw2();

  hA->SetMarkerColor(kBlack);
  hB->SetMarkerColor(kRed);
  hC->SetMarkerColor(kBlue);
  hD->SetMarkerColor(kGreen+2);

  hB->SetMarkerStyle(21);
  hC->SetMarkerStyle(24);
  hD->SetMarkerStyle(25);

  hR2->SetMarkerColor(kRed);
  hR2->SetMarkerStyle(24);

  hA->GetXaxis()->CenterTitle(1);
  hR1->GetXaxis()->CenterTitle(1);
  hDR->GetXaxis()->CenterTitle(1);

  rand.SetSeed(0);

  for (int i=0; i<N; ++i) {
    // a = rand.Gaus(mean, sigma1);
    // b = a;
    a = rand.Gaus(mean, sigma1);
    a *= rand.Gaus(mean, sigma);
    b = rand.Gaus(mean, sigma1);
    b *= rand.Gaus(mean, sigma);

    // c = rand.Gaus(mean, sigma1);
    // d = c;
    c = rand.Gaus(mean, sigma1);
    c *= rand.Gaus(mean, sigma);
    d = rand.Gaus(mean, sigma1);
    d *= rand.Gaus(mean, sigma);

    a += 1;
    b += 1;
    c += 1;
    d += 1;

    hA->Fill(a);
    hB->Fill(b);
    hC->Fill(c);
    hD->Fill(d);

    hR1->Fill( (a/b) );
    hR2->Fill( (c/d) );

    hDR->Fill( (a/b) / (c/d) );
  }

  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Divide(3,1);
  c1->cd(1);
  hA->Draw();
  hB->Draw("same");
  hC->Draw("same");
  hD->Draw("same");
  c1->cd(2);
  hR1->Draw();
  hR2->Draw("same");
  c1->cd(3);
  hDR->Draw();

  cout << "Mean = " << hDR->GetMean() << " RMS = " << hDR->GetRMS() << endl;

  return;
}
