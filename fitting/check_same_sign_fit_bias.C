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

  TH1F *h0 = new TH1F("h0","h0",500,0,1000);
  TH1F *h1 = new TH1F("h1","h1",500,0,1000);
  h0->Sumw2();
  h1->Sumw2();
  h1->SetMarkerColor(2);
  h1->SetLineColor(2);

  //  Jpsi_Mass->setRange("signal",3.6,3.76);
  Jpsi_Mass->setRange("M2045",2.0,4.5);
  Jpsi_Mass->setRange("M2242",2.2,4.2);

  RooFitResult *fitMall;
  RooFitResult *fitM;

  // RooFit cannot normalize Chebychev polynomials to a subrange (no analytic integral?)
  // using normal polynomials instead
  RooRealVar a("a","a",0.5);a.setConstant(false);
  RooRealVar b("b","b",0.01);b.setConstant(false);
  RooRealVar c("c","c",-0.005);c.setConstant(false);
  RooPolynomial bkg("bkg","bkg",*Jpsi_Mass,RooArgSet(a,b,c));
  ws->import(bkg);
  //  ws->factory("SUM::sigMassPDF_M2242(Ndummy*dummy,Nbkg_A*bkg)");

  RooDataSet *data[N];
  RooDataSet *redData;

  for (int i=0;i<N;++i) {
    cout << "Generating event " << i << "/" << N << endl;
    // cout << "Generate N = " << Nevents << " events with R_{psi} = " << ws->function("fracP_HI020")->getVal() << endl;
    data[i] = model->generate(*Jpsi_Mass,Nevents);
  }
    // cout << data->sumEntries() << endl;

  for (int i=0;i<N;++i) {
    cout << "Fitting event " << i << "/" << N << endl;
    fitMall = ws->pdf("bkg")->fitTo(*data[i],Extended(0),Hesse(1),Save(1),NumCPU(8),PrintEvalErrors(-1),Verbose(0),PrintLevel(-1));
    if (fitMall->statusCodeHistory(fitMall->numStatusHistory()-1) != 0) {i--; continue;}

    //  cout << "Fitted R_{psi} = " << ws->function("fracP_HI020")->getVal() << " +/- " << ws->function("fracP_HI020")->getPropagatedError(*fitMall) << endl;
    double yield = ws->pdf("bkg")->analyticalIntegral(1,"signal");
    cout << yield << endl;
    h0->Fill(yield);


    redData = (RooDataSet*)data[i]->reduce("Jpsi_Mass>2.2&&Jpsi_Mass<4.2");
    //  cout << redData->sumEntries() << endl;

    fitM = ws->pdf("bkg")->fitTo(*redData,Extended(0),Hesse(1),Save(1),NumCPU(8),PrintEvalErrors(-1),Verbose(0),PrintLevel(-1),Range("M2242"));
    if (fitM->statusCodeHistory(fitM->numStatusHistory()-1) != 0) {i--; continue;}

    // cout << "Fit over M2242: R_{psi} = " << ws->function("fracP_HI020")->getVal() << " +/- " << ws->function("fracP_HI020")->getPropagatedError(*fitM) << endl;
    // cout << "double ratio = " << ws->var("doubleRatio_HI020")->getVal() << " +/- " << ws->var("doubleRatio_HI020")->getError() << endl;

    // RooPlot *frame = Jpsi_Mass->frame(Range("M2242"));
    // redData->plotOn(frame);
    // ws->pdf("sigMassPDF_M2242")->plotOn(frame,Range("M2242"));
    //  redData->plotOn(frame,MarkerStyle(24),MarkerColor(kRed),Range("M2242"));

    yield = ws->pdf("bkg")->analyticalIntegral(1,"signal");
    cout << yield << endl;
    h1->Fill(yield);
  }

  TCanvas *c1 = new TCanvas("c1","c1");
  c1->cd();
  h0->Draw();
  h1->Draw("same");
  //  frame->Draw();
  cout << h0->GetMean() << "\t" << h0->GetRMS() << endl;
  cout << h1->GetMean() << "\t" << h1->GetRMS() << endl;

  c1->SaveAs(Form("toy_fits_N%i.pdf",N));

  TFile *outf = new TFile(Form("toy_same_sign_fits_N%i.root",N),"RECREATE");
  h0->Write();
  h1->Write();
  outf->Close();


  return;
}
