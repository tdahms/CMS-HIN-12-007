#include <iostream>
#include <fstream>

using namespace std;

void plot_pulls(string ifname, int nlines) {
  ifstream ifs;
  ifs.open(ifname.c_str(),ifstream::in);

  if (!ifs.good()) {
    cout << "Failed to open " << ifname << endl;
    return;
  }

  string binName;
  double binValue;
  int i = 0;
  int nbins = nlines - 1; // don't count last empty line
  TH1F *histo = new TH1F("histo","histo",nbins,0,nbins);
  while (ifs.peek() != EOF) {
    if (i>nbins) break;
    ifs >> binName >> binValue;
    cout << binName << ": " << binValue << endl;
    histo->GetXaxis()->SetBinLabel(i,binName.c_str());
    histo->SetBinContent(i++, binValue);
  }

  histo->GetYaxis()->SetTitle(ifname.c_str());
  TCanvas *c1 = new TCanvas("c1","c1");
  histo->Draw();

  return;
}

/* to produce input list do:
 * grep 'Mean' 20140129_M2242_FixCBtail_NoSumw2_DblMu0_NoMinos/ppFracLogCBG_*0.txt | awk -F / '{print $2}' | awk -F: '{print $1 " " $2}' | awk -F 'CBG_' '{print $2}'| awk -F ' ' '{print $1 " " $3 }'| awk -F '.txt' '{print $1 " " $2}' >> pp_CBG_Mean.txt
 */
