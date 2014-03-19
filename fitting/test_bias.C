#include <iostream>
#include <string>

#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"

#include "RooFit.h"
#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooGenericPdf.h"
#include "RooRandom.h"
#include "RooPlot.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooArgset.h"
#include "RooFitResult.h"

#include "RooStats/ModelConfig.h"

using namespace RooFit;
using namespace std;

void test_bias(string fname="20140311_M2242_DblMu0/ppFracLogCBG_pol1_rap16-24_pT3-30_cent0-100",
	       string truthPdfName="pol1",
	       string altPdfName="expPol1",
	       unsigned int Ntoys = 1000,
	       unsigned int seed = 0)
{
  TFile *infW = new TFile((fname+"_Workspace.root").c_str(),"READ");

  RooWorkspace * w = dynamic_cast<RooWorkspace*>( infW->Get("workspace") );
  //  RooStats::ModelConfig* bModel = (RooStats::ModelConfig*) w->obj("B_only_model");
  RooAbsPdf *truthBkg = (RooAbsPdf*) w->obj(truthPdfName.c_str());
  RooAbsPdf *altBkg = (RooAbsPdf*) w->obj(altPdfName.c_str());
  // RooRealVar *nTruthBkg = new RooRealVar("nTrughBkg","nTruthBkg",100);
  // RooRealVar *nAltBkg = new RooRealVar("nAltBkg","nAltBkg",100);
  // nAltBkg.setConstant(0);
  //  RooGenericPdf altBkgModel("altBkgModel","@0*@1",RooArgList(nAltBkg,*altBkg));

  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  RooRandom::randomGenerator()->SetSeed(seed);


  RooArgSet * pars = w->pdf("sigMassPDF")->getParameters(w->data("data"));
  pars->readFromFile((fname+"_allVars.txt").c_str());

  // std::cout << "\nstarting parameters from fit\n";
  // pars->Print("v");

  RooFitResult * nll = 0;
  RooArgSet * newPars = 0;
  TFile *nllf = new TFile((fname+"_fitResult.root").c_str());
  //  nllf->GetObject("sigMassPDFSim", nll);
  nll= dynamic_cast<RooFitResult*>( nllf->Get("fitresult_sigMassPDF_data") );
  //  assert(nll);
  newPars = new RooArgSet(nll->randomizePars());
  //  newPars->Print("v");
  nllf->Close();
  

  //  w->Print();
  // truthBkg->Print("v");
  // altBkg->Print("v");
  TH1F *hBias = new TH1F("hBias","hBias;N_{true} - N_{fit};events",100000,-1,1);
  hBias->Sumw2();

  for (unsigned int i=0; i<Ntoys; ++i) {
    //    cout << "Poisson( N = " << w->var("NBkg")->getVal() << ")" << endl;
    double Nevt = RooRandom::randomGenerator()->Poisson(w->var("NBkg")->getVal());
    //    cout << "Generate N = " << int(Nevt+0.5) << " background events." << endl;
  
    pars->readFromFile((fname+"_allVars.txt").c_str());
    nll->randomizePars();
    pars->assignValueOnly(*newPars);
    //    pars->Print("v");
    
    RooDataSet *pseudoData = truthBkg->generate(RooArgSet(*w->var("Jpsi_Mass")),NumEvents(int(Nevt+0.5)));
    RooPlot *frame = w->var("Jpsi_Mass")->frame();
    pseudoData->plotOn(frame);
    truthBkg->plotOn(frame);

    //  RooFitResult * fitM;
    altBkg->fitTo(*pseudoData,Extended(0),Hesse(1),Minos(0),Save(1),Verbose(0),PrintEvalErrors(-1),PrintLevel(-1));
    //  fitM->printMultiline(cout,1,1,"");
    altBkg->plotOn(frame,LineColor(kRed));

    RooRealVar *resolution = new RooRealVar("resolution","resolution",0.3);
    resolution->setVal( sqrt( w->var("coeffGaus")->getVal()*pow(w->function("sigmaSig2")->getVal(),2) +
			      (1 - w->var("coeffGaus")->getVal())*pow(w->function("sigmaSig1")->getVal(),2) ) );

    // +/- 1.5 sigma (+/- 1.177 sigma = FWHM)
    w->var("Jpsi_Mass")->setRange("signal",w->function("meanSig1P")->getVal() - 1.5*resolution->getVal(),w->function("meanSig1P")->getVal() + 1.5*resolution->getVal());

    RooAbsReal *nTruthBkg = truthBkg->createIntegral(*w->var("Jpsi_Mass"),NormSet(*w->var("Jpsi_Mass")),Range("signal")) ;
    RooAbsReal *nAltBkg = altBkg->createIntegral(*w->var("Jpsi_Mass"),NormSet(*w->var("Jpsi_Mass")),Range("signal")) ;

    double bias = (nTruthBkg->getVal()-nAltBkg->getVal())/nTruthBkg->getVal();
    // cout << "Truth: " << nTruthBkg->getVal() << endl;
    // cout << "Fit: " << nAltBkg->getVal() << endl;
    // cout << "Bias: " << bias/nTruthBkg->getVal() << endl;

    hBias->Fill(bias);
  }

  double NJpsi = w->var("NJpsi")->getVal();
  double fracP = w->var("fracP")->getVal();
  double errNJpsi = w->var("NJpsi")->getError();
  double errFracP = w->var("fracP")->getError();
  TFile *infF = new TFile((fname+"_fitResult.root").c_str(),"READ");
  RooFitResult *origFit= dynamic_cast<RooFitResult*>( infF->Get("fitresult_sigMassPDF_data") );
  //  w->import(origFit);
  double corr = origFit->correlation( *w->var("NJpsi") , *w->var("fracP") );
  //  std::cout << "Correlation between NJpsi and fracP: " << corr << std::endl;
  double Npsi2S = w->function("NPsiP")->getVal();
  double errNpsi2S = sqrt( pow(errNJpsi/NJpsi,2) + pow(errFracP/fracP,2) + 2*errNJpsi*errFracP*corr/Npsi2S )*Npsi2S;

  // cout << "psi(2S) yield: " << Npsi2S << "+/-" << errNpsi2S << endl;
  // cout << "relative uncertainty: " << errNpsi2S/Npsi2S*100 << "%" << endl;

  TCanvas *c1 = new TCanvas("c1","c1");
  c1->cd();
  //  frame->Draw();
  hBias->Draw();
  double median[1] = {0.0};
  double probSum[1] = {0.5};
  hBias->GetQuantiles(1,median,probSum);
  cout << "|Median (bias)| = " << median[0] << endl;


  cout << "|Median (bias)| / uncertainty: " << median[0]/(errNpsi2S/Npsi2S) << endl;

  //  infF->Close();
  //  infW->Close();
  return;
}
