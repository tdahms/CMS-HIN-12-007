#include <iostream>


void plot_toyMC(string fname="fracP_HI_pp.txt", double val=1.0)
{
  TH1F *hToyMC_HI = new TH1F("hToyMC_HI","R_{#psi(2S)}(PbPb) from ToyMC;R_{#psi(2S)}(PbPb);Number of pseudo-experiments",100,-1,1);
  TH1F *hToyMC_pp = new TH1F("hToyMC_pp","R_{#psi(2S)}(pp) from ToyMC;R_{#psi(2S)}(pp);Number of pseudo-experiments",100,-1,1);
  TH1F *hToyMC = new TH1F("hToyMC","Double ratio from ToyMC;R_{#psi(2S)}(PbPb)/R_{#psi(2S)}(pp);Number of pseudo-experiments",20000,-10,10);
  double fracP_HI=0.0;
  double fracP_pp=0.0;

  ifstream fin(fname.c_str(),ifstream::in);

  while (fin.good()) {
    fin >> fracP_HI >> fracP_pp;
    hToyMC_HI->Fill(fracP_HI);
    hToyMC_pp->Fill(fracP_pp);
    hToyMC->Fill(fracP_HI/fracP_pp);
  }

  fin.close();

  int max = hToyMC->GetXaxis()->FindBin(val-0.0001);
  cout << "cutting between " <<  hToyMC->GetXaxis()->GetBinLowEdge(max)
       << " and " << hToyMC->GetXaxis()->GetBinLowEdge(max+1) << endl;

  cout << "signficance of being larger than observed value" << endl;
  cout << "p-value = " << hToyMC->Integral(max,101)/hToyMC->GetEntries() << endl;
  cout << "sigma = " << TMath::ErfInverse(1.0-hToyMC->Integral(max,101)/hToyMC->GetEntries())*sqrt(2) << endl;


  cout << "signficance of being smaller than observed value" << endl;
  cout << "p-value = " << hToyMC->Integral(-1,max)/hToyMC->GetEntries() << endl;
  cout << "sigma = " << TMath::ErfInverse(1.0-hToyMC->Integral(-1,max)/hToyMC->GetEntries())*sqrt(2) << endl;

  hToyMC->Rebin(200);

  TH1F* hToyMC_above = (TH1F*) hToyMC->Clone("hToyMC_above");
  hToyMC_above->SetFillStyle(3013);
  hToyMC_above->SetFillColor(kBlue);
  hToyMC_above->SetLineColor(kBlue);
  if (val>1)
    hToyMC_above->GetXaxis()->SetRangeUser(val,10);
  else
    hToyMC_above->GetXaxis()->SetRangeUser(-10,val);

  double ymax = hToyMC->GetMaximum();

  hToyMC_HI->GetXaxis()->CenterTitle(true);
  hToyMC_pp->GetXaxis()->CenterTitle(true);
  hToyMC->GetXaxis()->CenterTitle(true);

  hToyMC_HI->GetXaxis()->SetRangeUser(-0.5,0.5);
  hToyMC_pp->GetXaxis()->SetRangeUser(-0.2,0.2);


  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  c1->Divide(2,1);
  c1->cd(1);
  hToyMC_HI->Draw();
  c1->cd(2);
  hToyMC_pp->Draw();


  TCanvas *c2 = new TCanvas("c2","c2");
  hToyMC->Draw();
  hToyMC_above->Draw("same");

  TLatex *cms = new TLatex(val+0.5, 0.9*ymax, "CMS");
  cms->SetTextFont(42);
  cms->SetTextSize(0.06);
  cms->SetTextColor(kRed);

  if (val<1)
    cms->SetX(val-3.5);

  cms->Draw();

  TLine *l = new TLine(val,0,val,ymax);
  l->SetVertical(1);
  l->SetLineColor(kRed);
  l->SetLineWidth(2);
  l->Draw();

  return;
}
