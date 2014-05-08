#include <iostream>

//#ifndef __CINT__
#include "RooGlobalFunc.h"
//#endif

#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooConstVar.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooFormulaVar.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooDataSet.h"
#include "RooFitResult.h"

#include "RooGenericPdf.h"
#include "RooMinuit.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"

using namespace RooFit;

void check_same_sign_fit_bias(const int N=1, string infname="fwd_pbpb_cent0-20_workspace.root")
{
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  TFile *inf = new TFile(infname.c_str(),"READ");
  RooWorkspace *ws = (RooWorkspace*) inf->Get("workspace");
  ws->var("Ndummy")->setConstant(true);

  RooRealVar *Jpsi_Mass = ws->var("Jpsi_Mass");
  RooAbsPdf *model = ws->pdf("model_A");

  double Nevents =  ws->var("Nbkg_A")->getVal();

  TH1F *h0 = new TH1F("h0","h0",(int) 1000,0,0.1);
  TH1F *h1 = new TH1F("h1","h1",(int) 1000,0,0.1);
  TH1F *h2 = new TH1F("h2","h2",(int) 1000,0,0.1);
  h0->Sumw2();
  h1->Sumw2();
  h1->SetMarkerColor(2);
  h1->SetLineColor(2);
  h2->SetMarkerColor(4);
  h2->SetLineColor(4);

  //  Jpsi_Mass->setRange("signal",3.6,3.76);
  Jpsi_Mass->setRange("M2045",2.0,4.5);
  Jpsi_Mass->setRange("M2242",2.2,4.2);
  cout << Jpsi_Mass->getBinning("signal") << endl;
  //  cout << "signal: " << Jpsi_Mass->min("signal") << "--" << Jpsi_Mass->max("signal") << endl;
  RooFitResult *fitMall;
  RooFitResult *fitM;

  // RooFit cannot normalize Chebychev polynomials to a subrange (no analytic integral?)
  // using normal polynomials instead
  RooRealVar a("a","a",0.5);a.setConstant(false);
  RooRealVar b("b","b",0.01);b.setConstant(false);
  RooRealVar c("c","c",-0.005);c.setConstant(false);
  RooPolynomial bkg("bkg","bkg",*Jpsi_Mass,RooArgSet(a,b,c));
  RooPolynomial bkg_M2242("bkg_M2242","bkg_M2242",*Jpsi_Mass,RooArgSet(a,b,c));
  ws->import(bkg);
  ws->import(bkg_M2242);
  //  ws->factory("SUM::sigMassPDF_M2242(Ndummy*dummy,Nbkg_A*bkg)");

  RooDataSet *data[N];
  RooDataSet *redData;

  RooPlot *frame = Jpsi_Mass->frame();
  double yield=0;
  RooAbsReal *yield_bkg = ws->pdf("bkg")->createIntegral(*Jpsi_Mass,NormSet(*Jpsi_Mass),Range("signal"));
  RooAbsReal *yield_bkg_M2242 = ws->pdf("bkg_M2242")->createIntegral(*Jpsi_Mass,NormSet(*Jpsi_Mass),Range("signal"));

  for (int i=0;i<N;++i) {
    cout << "Generating event " << i << "/" << N << endl;
    data[i] = model->generate(*Jpsi_Mass,Nevents);
    cout << data[i]->sumEntries() << endl;
    yield = data[i]->sumEntries("Jpsi_Mass>3.6&&Jpsi_Mass<3.76")/data[i]->sumEntries();
    cout << yield << endl;
    h2->Fill(yield);
  }

  for (int i=0;i<N;++i) {
    cout << "Fitting event " << i << "/" << N << endl;
    fitMall = ws->pdf("bkg")->fitTo(*data[i],Extended(0),Hesse(1),Save(1),NumCPU(8),PrintEvalErrors(-1),Verbose(0),PrintLevel(-1));
    if (fitMall->statusCodeHistory(fitMall->numStatusHistory()-1) != 0) {i--; continue;}

    cout << yield_bkg->getVal()  << endl;
    h0->Fill(yield_bkg->getVal());

    redData = (RooDataSet*)data[i]->reduce("Jpsi_Mass>2.2&&Jpsi_Mass<4.2");

    fitM = ws->pdf("bkg_M2242")->fitTo(*redData,Extended(0),Hesse(1),Save(1),NumCPU(8),PrintEvalErrors(-1),Verbose(0),PrintLevel(-1),Range("M2242"));
    if (fitM->statusCodeHistory(fitM->numStatusHistory()-1) != 0) {i--; continue;}

    cout << yield_bkg_M2242->getVal() << endl;
    h1->Fill(yield_bkg_M2242->getVal());
  }
  
  // data[N-1]->plotOn(frame);
  // ws->pdf("bkg")->plotOn(frame);
  // redData->plotOn(frame,MarkerStyle(24),MarkerColor(kRed));
  // ws->pdf("bkg_M2242")->plotOn(frame,LineColor(kRed),LineStyle(7));

  TCanvas *c1 = new TCanvas("c1","c1");
  c1->cd();
  frame->Draw();

  TCanvas *c2 = new TCanvas("c2","c2");
  c2->cd();
  h0->Draw();
  h1->Draw("histsame");
  h2->Draw("histsame");
  cout << h0->GetMean() << "\t" << h0->GetRMS() << endl;
  cout << h1->GetMean() << "\t" << h1->GetRMS() << endl;

  c2->SaveAs(Form("toy_fits_N%i.pdf",N));

  TFile *outf = new TFile(Form("toy_same_sign_fits_N%i.root",N),"RECREATE");
  h0->Write();
  h1->Write();
  outf->Close();


  return;
}
