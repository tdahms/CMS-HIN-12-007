#include <iostream>
#include <fstream>
#include "TH1F.h"
#include "TCanvas.h"

using namespace std;

void read_ratio(string ifname)
{
  ifstream ifs;
  ifs.open(ifname.c_str(),ifstream::in);

  if (!ifs.good()) {
    cout << "Failed to open " << ifname << endl;
    return;
  }
  
  double fracP;
  double fracPerror;
  TH1F *histo = new TH1F("histo","histo",400,-0.1,0.3);
  while (ifs.peek() != EOF) {
    ifs >> fracP >> fracPerror;
    histo->Fill(fracP);
  }

  histo->GetXaxis()->SetTitle("R_{#psi(2S)}");
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->cd();
  histo->Draw();

  int minBin = histo->FindFirstBinAbove(0);
  int maxBin = histo->FindLastBinAbove(0);

  cout << "min R = " << histo->GetBinLowEdge(minBin) << endl;
  cout << "max R = " << histo->GetBinLowEdge(maxBin) << endl;

  return;
}
