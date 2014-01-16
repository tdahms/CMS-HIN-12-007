#include "RooSimWSTool.h"
#include "RooSimPdfBuilder.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooCategory.h"
#include "RooWorkspace.h"
#include "RooFormulaVar.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooExtendPdf.h"
#include "RooKeysPdf.h"
#include "RooSimultaneous.h"
#include "RooGenericPdf.h"

#include "TTree.h"
#include "TFile.h"

#include "RooRandom.h"
#include "RooSimultaneous.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "TCanvas.h"
#include "RooMsgService.h"

//#include "fit1DData_pbpb.cpp"

#include <iostream>

bool readData(RooWorkspace &ws, TString HIfilename, TString ppFilename);
double computeRatio(RooRealVar& x, RooRealVar& y);
double computeRatioError(RooRealVar& x, RooRealVar& y, double correlation = 0.);
void defineMassBkg(RooWorkspace &ws);
void defineMassSig(RooWorkspace &ws);
void buildSimPdf(RooWorkspace& ws);

RooDataSet * genOppositeSignBackground(RooWorkspace& ws, int& Nhi,
				       int& Npp, bool doPoisson = true) {
  RooRealVar * Jpsi_mass = ws.var("Jpsi_Mass");
  TIter types(ws.cat("sample")->typeIterator());
  RooCatType * type;
  types.Reset();
  RooDataSet * osBkgData = 
    (RooDataSet *)ws.data("data")->emptyClone("toy_Data");
  RooDataSet * protoData;
  RooAbsPdf * bkg;
  int N;
  if (doPoisson) {
    Nhi = RooRandom::randomGenerator()->Poisson(Nhi);
    Npp = RooRandom::randomGenerator()->Poisson(Npp);
  }
  while ((type=(RooCatType*)types.Next())) {
    protoData = 
      (RooDataSet*)ws.data(TString::Format("data_%s", type->GetName()));
    bkg = ws.pdf(TString::Format("expFunct_%s", type->GetName()));
    if ((*type) == "HI")
      N = Nhi;
    else
      N = Npp;
    RooDataSet * tmpData = bkg->generate(RooArgSet(*Jpsi_mass), N,
					 RooFit::ProtoData(*protoData));
    osBkgData->append(*tmpData);
    delete tmpData;
  }
  return osBkgData;
}


RooDataSet * genOppositeSignSignal(RooWorkspace& ws, int Nhi, int Npp,double x2 = 1.0) {
  
  RooRealVar * Jpsi_mass = ws.var("Jpsi_Mass");
  RooDataSet * osSigData = 
    (RooDataSet *)ws.data("data")->emptyClone("toy_Data");
  RooAbsPdf * sig = ws.pdf("sigMassPDF_pp");

  //  RooAbsReal * f2_pp = ws.function("fracP_pp");
  RooRealVar * fracP_pp = ws.var("fracP_pp");

  RooDataSet * protoData;
  int N;
  static TString types[2] = {TString("pp"), TString("HI")};
  for (int i = 0; i<2; ++i) {
    protoData = 
      (RooDataSet*)ws.data(TString::Format("data_%s", types[i].Data()));
    if (types[i] == "HI") {
      N = Nhi;
      sig = ws.pdf("sigMassPDF_HI");
      ws.var("fracP_HI")->setVal(fracP_pp->getVal()*x2);
    } else {
      N = Npp;
      sig = ws.pdf("sigMassPDF_pp");
    }
    RooDataSet * tmpData = sig->generate(RooArgSet(*Jpsi_mass), N,
					 RooFit::ProtoData(*protoData));
    osSigData->append(*tmpData);
    delete tmpData;
  }
  return osSigData;
}

void toyMC(int Ntoys, TString parVals, TString outfname, 
	    TString randfname = "", TString systfname1 = "", TString systfname2 = "", TString systfname3 = "",
	    double ppMult = 1.0, unsigned int seed = 0,
	   double x2 = 1.0, double systErr = 0.0)
{
  RooWorkspace w;
  const double Jpsi_MassMin=2.6;
  const double Jpsi_MassMax=4.2;
  RooRealVar * Jpsi_Mass = new RooRealVar("Jpsi_Mass", "J/#psi mass", Jpsi_MassMin, Jpsi_MassMax, "GeV/c^{2}");
  Jpsi_Mass->setBins(40);
  w.import(*Jpsi_Mass);

  TString hidatafile("20120421/Data2011_DblMu0_cent0-100.root");
  TString ppdatafile("20120421/ppData2011_DblMu0_cent0-100.root");

  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  RooRandom::randomGenerator()->SetSeed(seed);

  readData(w, hidatafile, ppdatafile);
  defineMassSig(w);
  defineMassBkg(w);
  //  buildPdf(ws, true, useKeys);
  //  buildPdf(ws, false, useKeys);

  RooCategory* sample = w.cat("sample");//new RooCategory("sample","sample") ;
  // sample->defineType("HI") ;
  // sample->defineType("pp") ;
  //  w.factory("sample[HI,pp]");  //Index for simultaneous fit


  buildSimPdf(w);
  RooArgSet * pars = w.pdf("sigMassPDFSim")->getParameters(w.data("data"));
  pars->readFromFile(parVals);

  std::cout << "\nstarting parameters from fit\n";
  pars->Print("v");

  RooFitResult * nll = 0;
  RooArgSet * newPars = 0;
  if (randfname.Length() > 0) {
    TFile nllf(randfname);
    nllf.GetObject("sigMassPDFSim", nll);
    assert(nll);
    newPars = new RooArgSet(nll->randomizePars());
    newPars->Print("v");
    nllf.Close();
  }

  RooFitResult * systRes1 = 0;
  RooArgSet * systPars1 = 0;
  if (systfname1.Length() > 0) {
    TFile systf(systfname1);
    systf.GetObject("sigMassPDFSim", systRes1);
    assert(systRes1);
    systPars1 = new RooArgSet(systRes1->randomizePars());
    std::cout << "\nsystematic parameters: alpha and n\n";
    systPars1->Print("v");
    systf.Close();
  }

  RooFitResult * systRes2 = 0;
  RooArgSet * systPars2 = 0;
  if (systfname2.Length() > 0) {
    TFile systf(systfname2);
    systf.GetObject("sigMassPDFSim", systRes2);
    assert(systRes2);
    systPars2 = new RooArgSet(systRes2->randomizePars());
    std::cout << "\nsystematic parameters: cutx\n";
    systPars2->Print("v");
    systf.Close();
  }

  RooFitResult * systRes3 = 0;
  RooArgSet * systPars3 = 0;
  if (systfname1.Length() > 0) {
    TFile systf(systfname3);
    systf.GetObject("jpsiMassPDFSim", systRes3);
    assert(systRes3);
    systPars3 = new RooArgSet(systRes3->randomizePars());
    std::cout << "\nsystematic parameters: coeffGauss and sigmaSig1\n";
    systPars3->Print("v");
    systf.Close();
  }

  int Npp_tot(w.var("NJpsi_pp")->getVal() +
	      w.function("NPsiP_pp")->getVal() +
	      w.var("NBkg_pp")->getVal() + 0.5);
  int Nhi_tot(w.var("NJpsi_HI")->getVal() +
	      w.function("NPsiP_HI")->getVal() +
	      w.var("NBkg_HI")->getVal() + 0.5);

  Npp_tot *= ppMult;
  std::cout << "Npp_tot: " << Npp_tot << '\n'
	    << "Nhi_tot: " << Nhi_tot << '\n';


  TIter dataCats(sample->typeIterator());
  RooCatType * dataStr;
  dataCats.Reset();
  TString theCut;
  RooDataSet * tmpData;
  while ((dataStr=(RooCatType*)dataCats.Next())) {
    theCut = TString::Format("(sample == sample::%s)", dataStr->GetName());

    tmpData = (RooDataSet *)
      w.data("data")->reduce(RooFit::SelectVars(RooArgSet(*sample)), 
			     RooFit::Name(TString::Format("data_%s", 
							  dataStr->GetName())),
			     RooFit::Cut(theCut));
    w.import(*tmpData);
    delete tmpData;
  }

  ofstream fout(outfname.Data(), ios_base::out|ios_base::trunc); 

  int Npp_bkg(w.var("NBkg_pp")->getVal() + 0.5);
  int Nhi_bkg(w.var("NBkg_HI")->getVal() + 0.5);
  double x2_smeared = x2;
  RooDataSet * toyData;
  RooCategory * toyCat;
  RooSimultaneous * simPdfToy;
  RooArgSet * toyPars;
  RooRealVar * par;
  RooFitResult * fr;
  for (int i = 0; i < Ntoys; ++i) {
    pars->readFromFile(parVals);
    if (nll) {
      nll->randomizePars();
      pars->assignValueOnly(*newPars);
      //      pars->Print("v");
    }

    pars->setRealValue("NJpsi_pp", pars->getRealValue("NJpsi_pp")*ppMult);
    pars->setRealValue("NBkg_pp", pars->getRealValue("NBkg_pp")*ppMult);
    if (systRes1) {
      systRes1->randomizePars();
      pars->setRealValue("alpha", systPars1->getRealValue("alpha", 1.6));
      pars->setRealValue("enneW", systPars1->getRealValue("enneW", 2.3));
      pars->setRealValue("alpha_pp", systPars1->getRealValue("alpha_pp", 1.6));
      pars->setRealValue("enneW_pp", systPars1->getRealValue("enneW_pp", 2.3));
    }

    // if (systRes2) {
    //   //      systRes2->randomizePars();
    //   pars->setRealValue("cutx", systPars2->getRealValue("cutx", 1.6));
    // }

    if (systRes3) {
      systRes3->randomizePars();
      pars->setRealValue("sigmaSig1", systPars3->getRealValue("sigmaSig1", 0.092));
      pars->setRealValue("coeffGaus", systPars3->getRealValue("coeffGaus1", 0.092));
      pars->setRealValue("sigmaSig1_pp", systPars3->getRealValue("sigmaSig1_pp", 0.092));
      pars->setRealValue("coeffGaus_pp", systPars3->getRealValue("coeffGaus_pp", 0.092));
    }

    // std::cout << "Generating parameters\n";
    // pars->Print("v");
    Npp_bkg = int(w.var("NBkg_pp")->getVal() + 0.5);
    Nhi_bkg = int(w.var("NBkg_HI")->getVal() + 0.5);
    toyData = (RooDataSet *)w.data("data")->emptyClone();


    if (systErr>0.0)
      x2_smeared = RooRandom::randomGenerator()->Gaus(x2,systErr);
    else
      x2_smeared = x2;


    cout << "double ratio: " << x2 << " smeared: " << x2_smeared << endl;

    tmpData = genOppositeSignBackground(w, Nhi_bkg, Npp_bkg);
    toyData->append(*tmpData);
    delete tmpData;
    tmpData = genOppositeSignSignal(w, Nhi_tot-Nhi_bkg, Npp_tot-Npp_bkg,
     				    x2_smeared);
    toyData->append(*tmpData);
    delete tmpData;

    RooWorkspace wsToy("wsToy", "wsToy");
    wsToy.import(*toyData);
    delete toyData;
    //    wsToy.Print("v");

    defineMassSig(wsToy);
    defineMassBkg(wsToy);

    toyCat = wsToy.cat("sample");
    //    toyCat->Print("v");
    buildSimPdf(wsToy);
    simPdfToy = dynamic_cast<RooSimultaneous *> (wsToy.pdf("sigMassPDFSim"));
    toyPars = simPdfToy->getParameters(wsToy.data("data"));
    toyPars->readFromFile(parVals);
    // std::cout << "Toy parameters\n";
    // toyPars->Print("v");

    wsToy.var("alpha")->setConstant(kTRUE);
    wsToy.var("enneW")->setConstant(kTRUE);
    wsToy.var("alpha_pp")->setConstant(kTRUE);
    wsToy.var("enneW_pp")->setConstant(kTRUE);
    //    wsToy.var("cutx")->setConstant(kTRUE);
    wsToy.var("sigmaSig1")->setConstant(kTRUE);
    wsToy.var("sigmaSig1_pp")->setConstant(kTRUE);
    wsToy.var("coeffGaus")->setConstant(kTRUE);
    wsToy.var("coeffGaus_pp")->setConstant(kTRUE);

    toyPars->setRealValue("NJpsi_pp", 
     			  toyPars->getRealValue("NJpsi_pp")*ppMult);
    toyPars->setRealValue("NBkg_pp", toyPars->getRealValue("NBkg_pp")*ppMult);
    toyPars->setRealValue("fracP_HI", wsToy.function("fracP_pp")->getVal()*x2_smeared);
    //    toyPars->setRealValue("DoubleRatio", x2);

    // std::cout << "Modified toy parameters\n";
    // toyPars->Print("v");

    RooArgSet truePars;
    toyPars->snapshot(truePars, true);


    // wsToy.var("NJpsi_HI")->setConstant(kTRUE);
    // wsToy.var("NJpsi_pp")->setConstant(kTRUE);
    // wsToy.var("fracP_HI")->setConstant(kTRUE);


    toyData = (RooDataSet *)wsToy.data("data");

    fr = wsToy.pdf("sigMassPDFSim")->fitTo(*toyData, RooFit::Extended(1),
					   RooFit::Minos(false),
					   RooFit::PrintLevel((Ntoys==1) ? 1 : -1),
					   RooFit::SumW2Error(kTRUE),
					   RooFit::Save(true),RooFit::NumCPU(4));
    fout << "nll " << fr->minNll() << ' '
	 << "edm " << fr->edm() << ' '
	 << "covQual " << fr->covQual() << ' ';

    RooDataSet *toyData_HI = (RooDataSet *)toyData->reduce("sample==sample::HI");
    RooDataSet *toyData_pp = (RooDataSet *)toyData->reduce("sample==sample::pp");
    RooPlot *mframe_HI = wsToy.var("Jpsi_Mass")->frame();
    RooPlot *mframe_pp = wsToy.var("Jpsi_Mass")->frame();
    toyData_HI->plotOn(mframe_HI,RooFit::DataError(RooAbsData::SumW2),RooFit::XErrorSize(0),RooFit::MarkerSize(1));
    toyData_pp->plotOn(mframe_pp,RooFit::DataError(RooAbsData::SumW2),RooFit::XErrorSize(0),RooFit::MarkerSize(1));
    
    wsToy.pdf("sigMassPDF_HI")->plotOn(mframe_HI,RooFit::DrawOption("L"),RooFit::LineColor(kBlue));
    wsToy.pdf("sigMassPDF_pp")->plotOn(mframe_pp,RooFit::DrawOption("L"),RooFit::LineColor(kBlue));
    TCanvas *c1 = new TCanvas("c1","c1");
    c1->Divide(2,1);
    c1->cd(1);
    mframe_HI->Draw();
    c1->cd(2);
    mframe_pp->Draw();


    TIter finalPar(fr->floatParsFinal().createIterator());
    while ((par = (RooRealVar *)finalPar())) {
      double trueVal = truePars.getRealValue(par->GetName(), 0.);
      fout << par->GetName() << ' '
    	   << par->getVal() << ' '
    	   << par->getError() << ' '
    	   << trueVal << ' ';
    }
    
    fout << "x2 "
	 << computeRatio(*(wsToy.var("fracP_HI")), *(wsToy.var("fracP_pp"))) << ' '
	 << computeRatioError(*(wsToy.var("fracP_HI")), *(wsToy.var("fracP_pp")),
			      fr->correlation("fracP_HI", "fracP_pp"))
	 << ' ' << x2_smeared << ' ';
    fout << '\n';
    //    delete toyData;
    delete toyPars;
    delete fr;
  }
  delete pars;
  return;
}


bool readData(RooWorkspace &ws, TString HIfilename, TString ppFilename) {

  RooRealVar *Jpsi_Mass = ws.var("Jpsi_Mass");

  RooArgSet cols(*Jpsi_Mass);


  TFile *HIfile = new TFile(HIfilename);
  TFile *ppFile = new TFile(ppFilename);
  RooDataSet *hidata = (RooDataSet*)HIfile->Get("dataJpsi");
  RooDataSet *ppdata = (RooDataSet*)ppFile->Get("dataJpsi");
  delete HIfile;
  delete ppFile;

  RooCategory* sample = new RooCategory("sample","sample") ;
  sample->defineType("HI") ;
  sample->defineType("pp") ;
  //  ws.factory("sample[HI,pp]");  //Index for simultaneous fit

  RooDataSet data("data","data",cols,RooFit::Index(*(sample)),RooFit::Import("HI",*hidata),RooFit::Import("pp",*ppdata));

  return ws.import(data);  
}

double computeRatio(RooRealVar& x, RooRealVar& y) {
  assert(y.getVal() != 0);
  return x.getVal()/y.getVal();
}

double computeRatioError(RooRealVar& x, RooRealVar& y, double correlation) {
  double err2 = x.getError()*x.getError()/x.getVal()/x.getVal() +
    y.getError()*y.getError()/y.getVal()/y.getVal() - 
    2.*x.getError()*y.getError()/x.getVal()/y.getVal()*correlation;

  return fabs(computeRatio(x,y))*sqrt(err2);
}


void defineMassBkg(RooWorkspace &ws) {
  // 1st order polynomial
  ws.factory("Chebychev::pol1_HI(Jpsi_Mass,{coeffPol1[-0.8,-150.,150.]})");
  ws.factory("Chebychev::pol1_pp(Jpsi_Mass,{coeffPol1_pp[-0.8,-150.,150.]})");
  // 2nd order polynomial
  ws.factory("Chebychev::pol2_HI(Jpsi_Mass,{coeffPol1, coeffPol2[0.00,-15.,15.]})");
  ws.factory("Chebychev::pol2_pp(Jpsi_Mass,{coeffPol1_pp, coeffPol2_pp[0.00,-15.,15.]})");
  // Exponential
  ws.factory("Exponential::expFunct_HI(Jpsi_Mass,coeffExp[-0.1,-3.,1.])");
  ws.factory("Exponential::expFunct_pp(Jpsi_Mass,coeffExp_pp[-0.1,-3.,1.])");

  // 2nd Exponential
  ws.factory("Exponential::expFunct2_HI(Jpsi_Mass,coeffExp2[-0.2.,-3.,1.])");
  ws.factory("Exponential::expFunct2_pp(Jpsi_Mass,coeffExp2_pp[-0.2.,-3.,1.])");


  // Sum of two exponentials
  RooRealVar cutx("cutx","cutx",3.4,3.2,3.5); cutx.setConstant(false); ws.import(cutx); //Region below than this cut will use cTau1, above than it will use cTau2. determined by free fit.
  RooRealVar cutx_pp("cutx_pp","cutx_pp",3.4,3.2,3.5); cutx_pp.setConstant(false); ws.import(cutx_pp); //Region below than this cut will use cTau1, above than it will use cTau2. determined by free fit.

  RooGenericPdf twoExpFunct_HI("twoExpFunct_HI","twoExpFunct_HI","1.0/abs(@1)*(exp(@0*@1)*(@0<@3)+exp(@0*@2 - @3*(@2-@1) )*(@0>=@3))", RooArgList(*(ws.var("Jpsi_Mass")), *(ws.var("coeffExp")), *(ws.var("coeffExp2")), *(ws.var("cutx"))));
  ws.import(twoExpFunct_HI);

  RooGenericPdf twoExpFunct_pp("twoExpFunct_pp","twoExpFunct_pp","1.0/abs(@1)*(exp(@0*@1)*(@0<@3)+exp(@0*@2 - @3*(@2-@1) )*(@0>=@3))", RooArgList(*(ws.var("Jpsi_Mass")), *(ws.var("coeffExp_pp")), *(ws.var("coeffExp2_pp")), *(ws.var("cutx_pp"))));
  ws.import(twoExpFunct_pp);

  return;
}

void defineMassSig(RooWorkspace &ws) {
  //////// Candidates for signal
  // Normal gaussians
  ws.factory("Gaussian::signalG1(Jpsi_Mass,meanSig1[3.0969,3.05,3.15],sigmaSig1[0.08,0.06,0.2])");
  ws.factory("Gaussian::signalG1_pp(Jpsi_Mass,meanSig1_pp[3.0969,3.05,3.15],sigmaSig1_pp[0.08,0.06,0.2])");

  // Crystall Ball
  ws.factory("CBShape::signalCB2WN(Jpsi_Mass,meanSig1,sigmaSig2[0.03,0.01,0.06],alpha[1.,0.,3.],enneW[5.,1.,50.])");
  ws.factory("CBShape::signalCB2WN_pp(Jpsi_Mass,meanSig1_pp,sigmaSig2_pp[0.03,0.01,0.06],alpha_pp[1.,0.,3.],enneW_pp[5.,1.,50.])");

  // Fix Jpsi-psi' mass difference
  //  RooFormulaVar meanSig1P("meanSig1P","@0+0.58917",RooArgList(*(ws.var("meanSig1")))); ws.import(meanSig1P);
  // Fix Jpsi-psi' mass ratio
  RooFormulaVar meanSig1P("meanSig1P","@0*1.1902",RooArgList(*(ws.var("meanSig1")))); ws.import(meanSig1P);
  // Fix resolution scale: sigma_MJpsi/MJpsi = sigma_Mpsi'/Mpsi'
  RooFormulaVar sigmaSig1P("sigmaSig1P","@0*1.1902",RooArgList(*(ws.var("sigmaSig1")))); ws.import(sigmaSig1P);
  RooFormulaVar sigmaSig2P("sigmaSig2P","@0*1.1902",RooArgList(*(ws.var("sigmaSig2")))); ws.import(sigmaSig2P);
  ws.factory("Gaussian::signalG1P(Jpsi_Mass,meanSig1P,sigmaSig1P)");
  ws.factory("CBShape::signalCB2WNP(Jpsi_Mass,meanSig1P,sigmaSig2P,alpha,enneW)");

  RooFormulaVar meanSig1P_pp("meanSig1P_pp","@0*1.1902",RooArgList(*(ws.var("meanSig1_pp")))); ws.import(meanSig1P_pp);
  RooFormulaVar sigmaSig1P_pp("sigmaSig1P_pp","@0*1.1902",RooArgList(*(ws.var("sigmaSig1_pp")))); ws.import(sigmaSig1P_pp);
  RooFormulaVar sigmaSig2P_pp("sigmaSig2P_pp","@0*1.1902",RooArgList(*(ws.var("sigmaSig2_pp")))); ws.import(sigmaSig2P_pp);
  ws.factory("Gaussian::signalG1P_pp(Jpsi_Mass,meanSig1P_pp,sigmaSig1P_pp)");
  ws.factory("CBShape::signalCB2WNP_pp(Jpsi_Mass,meanSig1P_pp,sigmaSig2P_pp,alpha_pp,enneW_pp)");

  //////// Sum of signal functions
  // Sum of gaussian 1 and crystall ball 2 with wide n
  ws.factory("SUM::sigCB2WNG1(coeffGaus[0.1,0.05,0.93]*signalG1,signalCB2WN)");
  ws.factory("SUM::sigCB2WNG1P(coeffGaus*signalG1P,signalCB2WNP)");

  ws.factory("SUM::sigCB2WNG1_pp(coeffGaus_pp[0.1,0.05,0.93]*signalG1_pp,signalCB2WN_pp)");
  ws.factory("SUM::sigCB2WNG1P_pp(coeffGaus_pp*signalG1P_pp,signalCB2WNP_pp)");

  // ws.factory("SUM::sigCBGJpsiPsiP(fracP*sigCB2WNG1P,sigCB2WNG1)");
  // ws.factory("SUM::sigCBJpsiPsiP(fracP*signalCB2WNP,signalCB2WN)");

  // ws.factory("SUM::sigCBGJpsiPsiP_pp(fracP_pp*sigCB2WNG1P_pp,sigCB2WNG1_pp)");
  // ws.factory("SUM::sigCBJpsiPsiP_pp(fracP_pp*signalCB2WNP_pp,signalCB2WN_pp)");

  return;
}


void buildSimPdf(RooWorkspace& ws)
{
  RooRealVar *NBkg_pp  = new RooRealVar("NBkg_pp","pp background yield", 2000.0,1.0,500000.0);   ws.import(*NBkg_pp);
  RooRealVar *NBkg_HI  = new RooRealVar("NBkg_HI","PbPb background yield", 2000.0,1.0,500000.0);   ws.import(*NBkg_HI);


  RooRealVar *NJpsi_pp  = new RooRealVar("NJpsi_pp","pp J/psi yield", 4000.0, 1.0, 50000.0);  ws.import(*NJpsi_pp);
  RooRealVar *NJpsi_HI  = new RooRealVar("NJpsi_HI","PbPb J/psi yield", 4000.0, 1.0, 50000.0);  ws.import(*NJpsi_HI);
  RooRealVar *fracP_pp  = new RooRealVar("fracP_pp","pp psi(2S) fraction" ,0.1);fracP_pp->setConstant(false); ws.import(*fracP_pp);
  RooRealVar *fracP_HI  = new RooRealVar("fracP_HI","PbPb psi(2S) fraction" ,0.1);fracP_HI->setConstant(false); ws.import(*fracP_HI);

  RooFormulaVar *NPsiP_HI  = new RooFormulaVar("NPsiP_HI","@0*@1",RooArgList(*(ws.var("NJpsi_HI")),*(ws.var("fracP_HI"))));   ws.import(*NPsiP_HI);
  RooFormulaVar *NPsiP_pp = new RooFormulaVar("NPsiP_pp", "@0*@1", RooArgList(*(ws.var("NJpsi_pp")),*(ws.var("fracP_pp"))));  ws.import(*NPsiP_pp);

  RooFormulaVar *DoubleRatio  = new RooFormulaVar("DoubleRatio","@0/@1",RooArgList(*(ws.var("fracP_HI")),*(ws.var("fracP_pp"))));   ws.import(*DoubleRatio);

  RooFormulaVar *NPsiP_mix  = new RooFormulaVar("NPsiP_mix","@0*@1",RooArgList(*(ws.var("NJpsi_pp")),*(ws.function("fracP_HI"))));   ws.import(*NPsiP_mix);

  ws.factory("SUM::sigMassPDF_pp(NJpsi_pp*sigCB2WNG1_pp,NPsiP_pp*sigCB2WNG1P_pp,NBkg_pp*expFunct_pp)");
  ws.factory("SUM::sigMassPDF_HI(NJpsi_HI*sigCB2WNG1,NPsiP_HI*sigCB2WNG1P,NBkg_HI*twoExpFunct_HI)");
  ws.factory("SUM::sigMassPDF_mix(NJpsi_pp*sigCB2WNG1_pp,NPsiP_mix*sigCB2WNG1P_pp,NBkg_pp*expFunct_pp)");
  
  ws.factory("SIMUL::sigMassPDFSim(sample,HI=sigMassPDF_HI,pp=sigMassPDF_pp)");

  return;
}
