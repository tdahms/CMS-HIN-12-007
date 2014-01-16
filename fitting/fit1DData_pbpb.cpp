#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <TROOT.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TRandom.h>

#include <RooFit.h>
#include <RooGlobalFunc.h>
#include <RooCategory.h>
#include "RooGenericPdf.h"

#include "RooAbsPdf.h"
//#include "RooRealProxy.h"
#include "RooDataHist.h"

//#include <RooHistPdfConv.h>
//#include <RooExp2.h>
#include <RooWorkspace.h>
#include <RooBinning.h>
#include <RooHistPdf.h>
#include <RooProdPdf.h>
#include <RooAddPdf.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooHist.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooChebychev.h>

using namespace RooFit;
using namespace std;

bool superImpose = true;
bool analyticBlifetime = true;
bool narrowSideband = false;
bool oneGaussianResol = false;

void getOptRange(string &ran,float *min,float *max);
void setWSRange(RooWorkspace *ws);
//RooBinning setCtBinning(float lmin,float lmax);
void defineMassBkg(RooWorkspace *ws);
void defineMassSig(RooWorkspace *ws);


int main(int argc, char* argv[]) {
  gROOT->Macro("./rootlogon.C");

  string FileName, FileNameMC1, FileNameMC2;
  string mBkgFunct, mJpsiFunct, mPsiPFunct;
  bool prefitMass = false;
  bool prefitSignalCTau = false;
  bool prefitBkg = false;
  int  isGG = 0;
  bool fracfix = true;
  int isMB = false;
  int isMBCtau = false;
  bool isPT = false;
  string prange, lrange, yrange, crange, phirange, errrange;
  string dirPre;
  string rpmethod = "etHF";
  bool fitRatio=false;
  bool logScale=false;
  bool fitSubRange=false;
  bool isPbPb=true;
  bool zoom = false;
  bool isPaper = false;
  bool fixCBtoMC=true;
  bool cutNonPrompt=false;

  // *** Check options
  for (int i=1; i<argc; i++) {
    char *tmpargv = argv[i];
    switch (tmpargv[0]) {
    case '-':{
      switch (tmpargv[1]) {
      case 'f':
	FileName = argv[i+1];
	cout << "Fitted data file: " << FileName << endl;
	break;
      case 'v':
	mJpsiFunct = argv[i+1];
	mPsiPFunct = argv[i+2];
	mBkgFunct = argv[i+3];
	cout << "Mass J/psi function: " << mJpsiFunct << endl;
	cout << "Mass psi(2S) function: " << mPsiPFunct << endl;
	cout << "Mass background function: " << mBkgFunct << endl;
	break;
      case 'd':
	dirPre = argv[i+1];
	cout << "Prefix for all result files: " << dirPre << endl;
	break;
      case 'u':
	prefitMass = true;
	cout << "Turn on: signal(=data, depends on muon pair type) mass pre-fitting" << endl;
	break;
      case 'p':
	prange = argv[i+1];
	cout << "pT range: " << prange << " GeV/c" << endl;
	break;
      case 'y':
	yrange = argv[i+1];
	cout << "Rapidity range: " << yrange << "" << endl;
	break;
      case 't':
	crange = argv[i+1];
	cout << "Centrality range: " << crange << " %" << endl;
	break;
      case 'x':
	isMB = atoi(argv[i+1]);
	cout << "Inclusive fit option: " << isMB << endl;
	break;
      case 'r':
	fitRatio = atoi(argv[i+1]);
	cout << "Fit ratio (instead of 2S yield): " << fitRatio << endl;
	break;
      case 'l':
	logScale = atoi(argv[i+1]);
	cout << "plot on log scale: " << logScale << endl;
	break;
      case 's':
	fitSubRange = atoi(argv[i+1]);
	cout << "fit sub ranges: " << fitSubRange << endl;
	break;
      case 'n':
	cutNonPrompt = atoi(argv[i+1]);
	cout << "remove non-prompt J/psi: " << cutNonPrompt << endl;
	break;
      case 'z':
	isPaper = atoi(argv[i+1]);
	cout << "paper plots: " << isPaper << endl;
	break;
      }
    }
    }
  }// End check options
 
  string mBkgFunctP = mBkgFunct + "P";

  float pmin=0, pmax=0, ymin=0, ymax=0, cmin=0, cmax=0, psmax=0, psmin=0;
  getOptRange(prange,&pmin,&pmax);
  getOptRange(crange,&cmin,&cmax);
  getOptRange(yrange,&ymin,&ymax);
  getOptRange(phirange,&psmin,&psmax);


  // *** TFile for saving fitting results
  string resultFN;
  resultFN = dirPre + "_" + mBkgFunct + "_rap" + yrange + "_pT" + prange + "_cent" + crange + "_dPhi" + phirange + "_fitResult.root";
  TFile resultF(resultFN.c_str(),"RECREATE");

  // *** Read Data files

  if (FileName.find("pp")!=string::npos)
    isPbPb = false;
  TFile fInData(FileName.c_str());
  cout << FileName.c_str() << endl;
  if (fInData.IsZombie()) { cout << "CANNOT open data root file\n"; return 1; }
  fInData.cd();
  RooDataSet *data = (RooDataSet*)fInData.Get("dataJpsi");
  data->SetName("data");

  // Create workspace to play with
  RooWorkspace *ws = new RooWorkspace("workspace");

  // Reduce "dataMC" with given ranges/cuts
  char reduceDS[300];
  if (cutNonPrompt && !isPbPb) {
    if (prange=="6.5-30.0" && yrange=="0.0-2.4")
      sprintf(reduceDS,"Jpsi_Pt>%.2f && Jpsi_Pt<%.2f && abs(Jpsi_Y)>%.2f && abs(Jpsi_Y)<%.2f && Jpsi_Ct<0.04",pmin,pmax,ymin,ymax);
    else if (prange=="6.5-30.0" && yrange=="0.0-1.6")
      sprintf(reduceDS,"Jpsi_Pt>%.2f && Jpsi_Pt<%.2f && abs(Jpsi_Y)>%.2f && abs(Jpsi_Y)<%.2f && Jpsi_Ct<0.04",pmin,pmax,ymin,ymax);
    else if (prange=="3.0-30.0" && yrange=="1.6-2.4")
      sprintf(reduceDS,"Jpsi_Pt>%.2f && Jpsi_Pt<%.2f && abs(Jpsi_Y)>%.2f && abs(Jpsi_Y)<%.2f && Jpsi_Ct<0.11",pmin,pmax,ymin,ymax);
    else if (prange=="3.0-6.5" && yrange=="1.6-2.4")
      sprintf(reduceDS,"Jpsi_Pt>%.2f && Jpsi_Pt<%.2f && abs(Jpsi_Y)>%.2f && abs(Jpsi_Y)<%.2f && Jpsi_Ct<0.11",pmin,pmax,ymin,ymax);
    else if (prange=="6.5-30.0" && yrange=="1.6-2.4")
      sprintf(reduceDS,"Jpsi_Pt>%.2f && Jpsi_Pt<%.2f && abs(Jpsi_Y)>%.2f && abs(Jpsi_Y)<%.2f && Jpsi_Ct<0.06",pmin,pmax,ymin,ymax);
  }
  else if (cutNonPrompt && isPbPb) {
    if (prange=="6.5-30.0" && yrange=="0.0-2.4")
      sprintf(reduceDS,"Jpsi_Pt>%.2f && Jpsi_Pt<%.2f && abs(Jpsi_Y)>%.2f && abs(Jpsi_Y)<%.2f && Jpsi_Ct<0.04",pmin,pmax,ymin,ymax);
    else if (prange=="6.5-30.0" && yrange=="0.0-1.6")
      sprintf(reduceDS,"Jpsi_Pt>%.2f && Jpsi_Pt<%.2f && abs(Jpsi_Y)>%.2f && abs(Jpsi_Y)<%.2f && Jpsi_Ct<0.04",pmin,pmax,ymin,ymax);
    else if (prange=="3.0-30.0" && yrange=="1.6-2.4")
      sprintf(reduceDS,"Jpsi_Pt>%.2f && Jpsi_Pt<%.2f && abs(Jpsi_Y)>%.2f && abs(Jpsi_Y)<%.2f && Jpsi_Ct<0.08",pmin,pmax,ymin,ymax);
    else if (prange=="3.0-6.5" && yrange=="1.6-2.4")
      sprintf(reduceDS,"Jpsi_Pt>%.2f && Jpsi_Pt<%.2f && abs(Jpsi_Y)>%.2f && abs(Jpsi_Y)<%.2f && Jpsi_Ct<0.10",pmin,pmax,ymin,ymax);
    else if (prange=="6.5-30.0" && yrange=="1.6-2.4")
      sprintf(reduceDS,"Jpsi_Pt>%.2f && Jpsi_Pt<%.2f && abs(Jpsi_Y)>%.2f && abs(Jpsi_Y)<%.2f && Jpsi_Ct<0.06",pmin,pmax,ymin,ymax);
  }
  else
    sprintf(reduceDS,"Jpsi_Pt>%.2f && Jpsi_Pt<%.2f && abs(Jpsi_Y)>%.2f && abs(Jpsi_Y)<%.2f",pmin,pmax,ymin,ymax);

  cout << "reduceDS for MC and data: " << reduceDS << endl;

  RooDataSet *redData = (RooDataSet*)data->reduce(reduceDS);
  ws->import(*redData);
 
  setWSRange(ws);

  // Draw data
  string titlestr;

  int nbins = 125;
  // Define binning for true lifetime
  //  ws->var("Jpsi_Mass")->setBins(80);
  ws->var("Jpsi_Mass")->setBins(nbins);
  //  ws->var("Jpsi_Mass")->setBins(50);
  //  ws->var("Jpsi_Mass")->setBins(55);

  // Additional cuts on data and get sub-datasets/histograms
  RooDataSet *redDataCut;
  string reduceDSstr;
  if (isGG == 0) {
    redDataCut = (RooDataSet*)redData->reduce("Jpsi_Type == Jpsi_Type::GG");
  } else if (isGG == 1) {
    redDataCut = (RooDataSet*)redData->reduce("Jpsi_Type == Jpsi_Type::GT || Jpsi_Type == Jpsi_Type::GG");
  } else {
    redDataCut = (RooDataSet*)redData->reduce("Jpsi_Ct < 600000.");
  }

  RooDataHist *binData = new RooDataHist("binData","binData",RooArgSet( *(ws->var("Jpsi_Mass")) ), *redDataCut);
  cout << "DATA :: N events to fit: " << binData->sumEntries() << endl;

  // SYSTEMATICS 1 (very sidebands)
  RooDataSet *redDataSB;
  if (narrowSideband) {
   redDataSB = (RooDataSet*) redDataCut->reduce("Jpsi_Mass<2.8 || Jpsi_Mass>3.4");
  } else {
   redDataSB = (RooDataSet*) redDataCut->reduce("Jpsi_Mass<2.9 || Jpsi_Mass>3.3");
  }
  RooDataHist *binDataSB = new RooDataHist("binDataSB","Data distribution for background",RooArgSet( *(ws->var("Jpsi_Mass")) ),*redDataSB);
  RooDataSet *redDataSIG = (RooDataSet*)redDataCut->reduce("Jpsi_Mass > 2.9 && Jpsi_Mass < 3.3");

  //  ws->factory("Chebychev::test(Jpsi_Mass,{coeffPol[-0.8,-1.0,1.0]})");

  // *** Define PDFs with parameters (mass and ctau)
  // J/psi mass parameterization

  RooRealVar aa("aa","aa",0.5,-1,1);
  RooRealVar ab("ab","ab",-0.5,-1,1);
  RooChebychev tmpPol("tmpPol","tmpPol",*(ws->var("Jpsi_Mass")),RooArgSet(aa,ab));

  defineMassBkg(ws);
  defineMassSig(ws);

  char funct[125];
  string partTit, partFile;
  if (isGG == 0) { partTit = "glb-glb"; partFile = "GG"; }
  else if (isGG == 1) { partTit = "glb-trk"; partFile = "GT"; }
  else { partTit = "all"; partFile = "ALL"; }

  // Binning for invariant mass distribution

  RooBinning rbm(2.0,4.5);
  //  rbm.addUniform(80,2.6,4.2);
  rbm.addUniform(nbins,2.0,4.5);
  //  rbm.addUniform(50,2.6,4.6);
  //  rbm.addUniform(55,2.5,4.7);

  // Global TLatex, TH1, TGraph objects for drawing
  TLatex *t = new TLatex();
  t->SetNDC(); t->SetTextAlign(12);

  Double_t fx[2], fy[2], fex[2], fey[2];
  TGraphErrors *gfake1 = new TGraphErrors(2,fx,fy,fex,fey);
  gfake1->SetMarkerStyle(20); gfake1->SetMarkerSize(0.9);
  TH1F hfake11 = TH1F("hfake1","hfake1",100,200,300);
  hfake11.SetLineColor(kBlue); hfake11.SetLineWidth(4); hfake11.SetLineStyle(7); hfake11.SetFillColor(kAzure-9); hfake11.SetFillStyle(1001);
  TH1F hfake21 = TH1F("hfake21","hfake21",100,200,300);
  hfake21.SetLineColor(kBlack); hfake21.SetLineWidth(4); hfake21.SetFillColor(kBlack); hfake21.SetFillStyle(3354);
  TH1F hfake22 = TH1F("hfake22","hfake22",100,200,300);
  hfake22.SetLineColor(kRed); hfake22.SetLineWidth(4); hfake22.SetLineStyle(2);
  TH1F hfake23 = TH1F("hfake23","hfake23",100,200,300);
  hfake23.SetLineColor(kGreen+2); hfake23.SetLineWidth(4); hfake23.SetFillColor(kGreen); hfake23.SetFillStyle(3354);
  TH1F hfake31 = TH1F("hfake3","hfake3",100,200,300);
  hfake31.SetLineColor(kRed); hfake31.SetMarkerStyle(kCircle); hfake31.SetLineWidth(4); hfake31.SetMarkerColor(kRed); hfake31.SetLineStyle(9); hfake31.SetFillColor(kRed-7); hfake31.SetFillStyle(3444);


  RooFitResult *fitM;
  RooFitResult *fitM_pre;
  RooFitResult *fitP;
  struct PARAM {
    double fracP; double fracPErr;
    double coeffExp;  double coeffExpErr;
    double coeffExp2;  double coeffExp2Err;
    double cutx;  double cutxErr;
    double bkgMean; double bkgMeanErr;
    double bkgSigma; double bkgSigmaErr;
    double coeffPol1; double coeffPol1Err;
    double coeffPol2; double coeffPol2Err;
    double coeffPol3; double coeffPol3Err;
    double coeffPol4; double coeffPol4Err;
    double coeffPol5; double coeffPol5Err;
    double coeffGaus; double coeffGausErr;
    double meanSig1;  double meanSig1Err;
    double sigmaSig1; double sigmaSig1Err;
    double sigmaSig2; double sigmaSig2Err;
    double wideFactor; double wideFactorErr;
    double alpha;     double alphaErr;
    double enne;      double enneErr;
    //    double enneW;     double enneWErr;
  };

  PARAM inputCB;
  PARAM inputCBG;

  bool centConst = false;  //False: fit w/o any constrained parameters (centrality dep.)
  bool dPhiConst = false;  //False: fit w/o any constrained parameters (dPhi dep.)
  double inputN[4] = {0};  //Number of Jpsi, psiP, fracP, and background in the 0-1.571 rad bin
  if (isMB != 0) {
    if (isMB != 4 && (mJpsiFunct.compare("sigCB2WNG1") || mBkgFunct.compare("expFunct"))) {
      cout << "For fit systematics:\n";
      cout << "1) Signal shape check: should use default runOpt (4)\n";
      cout << "2) Background shape check: should use defatul runOpt (4)\n";
      cout << "3) Constrain fit: should use sigCB2WNG1 for signal shape, expFunct for bkg shape with runOpt (3)\n";
      return -1;
    }
  } //End of isMB != 0

  // RooRealVar *NJpsi  = new RooRealVar("NJpsi","J/psi yield", 4000.0, 0.0, 50000.0);  ws->import(*NJpsi);
  // RooRealVar *NBkg  = new RooRealVar("NBkg","Brackground yield", 2000.0,0.0,500000.0);   ws->import(*NBkg);
  // RooRealVar *NBkg2  = new RooRealVar("NBkg2","Second brackground yield", 2000.0,1.0,500000.0);   ws->import(*NBkg2);
  RooRealVar *NJpsi  = new RooRealVar("NJpsi","J/psi yield",0.5*binData->sumEntries(),0.0,2.0*binData->sumEntries());  ws->import(*NJpsi);
  RooRealVar *NBkg  = new RooRealVar("NBkg","Brackground yield", 0.5*binData->sumEntries(),0.0,2.0*binData->sumEntries());   ws->import(*NBkg);
  RooRealVar *NBkg2  = new RooRealVar("NBkg2","Second brackground yield",0.5*binData->sumEntries(),0.0,2.0*binData->sumEntries());   ws->import(*NBkg2);

  if (fitRatio) {
    RooRealVar *fracP  = new RooRealVar("fracP","psi(2S) fraction" ,0.01);
    //    fracP->setVal(0.0);fracP->setConstant(true);ws->import(*fracP);
    fracP->setConstant(false); ws->import(*fracP);
    RooFormulaVar *NPsiP = new RooFormulaVar("NPsiP", "@0*@1", RooArgList(*(ws->var("NJpsi")),*(ws->var("fracP"))));  ws->import(*NPsiP);
  }
  else {
    //    RooRealVar *NPsiP  = new RooRealVar("NPsiP","psi(2S) yield", 400.0, 1.0, 50000.0);  ws->import(*NPsiP);
    RooRealVar *NPsiP  = new RooRealVar("NPsiP","psi(2S) yield", 0.01*binData->sumEntries(),0.0,2.0*binData->sumEntries());  ws->import(*NPsiP);
  }
  

  if (fitSubRange) {
    sprintf(funct,"SUM::jpsiMassPDF(NJpsi*%s,NBkg*%s)",mJpsiFunct.c_str(),mBkgFunct.c_str());
    ws->factory(funct);

    sprintf(funct,"SUM::psipMassPDF(NPsiP*%s,NBkg2*%s)",mPsiPFunct.c_str(),mBkgFunctP.c_str());
    ws->factory(funct);
  }
  else {
    if (prefitMass) {
      sprintf(funct,"SUM::jpsiMassPDF(NJpsi*%s,NBkg*%s)",mJpsiFunct.c_str(),mBkgFunct.c_str());
      ws->factory(funct);
    }

    sprintf(funct,"SUM::sigMassPDF(NJpsi*%s,NPsiP*%s,NBkg*%s)",mJpsiFunct.c_str(),mPsiPFunct.c_str(),mBkgFunct.c_str());
    ws->factory(funct);
    if (!isPbPb) {
      RooRealVar *fracP_HI  = new RooRealVar("fracP_HI","psi(2S) fraction" ,0.1,0.0,1.0);   ws->import(*fracP_HI);
    
      // update these numbers later for the textbook cartoon
      if (prange=="6.5-30.0" && yrange=="0.0-2.4")
	ws->var("fracP_HI")->setVal(0.031);
      else if (prange=="6.5-30.0" && yrange=="0.0-1.6")
	ws->var("fracP_HI")->setVal(0.024);
      else if (prange=="3.0-30.0" && yrange=="1.6-2.4")
	ws->var("fracP_HI")->setVal(0.105);
      else if (prange=="3.0-6.5" && yrange=="1.6-2.4")
	ws->var("fracP_HI")->setVal(0.170);
      else if (prange=="6.5-30.0" && yrange=="1.6-2.4")
	ws->var("fracP_HI")->setVal(0.0516);
      else
	ws->var("fracP_HI")->setVal(0.0);

      RooFormulaVar *NPsiP_HI = new RooFormulaVar("NPsiP_HI", "@0*@1", RooArgList(*(ws->var("NJpsi")),*(ws->var("fracP_HI"))));  ws->import(*NPsiP_HI);
      sprintf(funct,"SUM::sigMassPDF_mix(NJpsi*%s,NPsiP_HI*%s,NBkg*%s)",mJpsiFunct.c_str(),mPsiPFunct.c_str(),mBkgFunct.c_str());
      ws->factory(funct);
    }
  }
  // if (!yrange.compare("0.0-1.2")) {
  //   sprintf(funct,"SUM::sigMassPDF(NJpsi[4000.0,1.0,50000.0]*%s,NBkg[2000.0,1.0,500000.0]*%s)",mSigFunct.c_str(),mBkgFunct.c_str());
  // } else if (!yrange.compare("1.6-2.4")) {
  //   sprintf(funct,"SUM::sigMassPDF(NJpsi[1000.0,1.0,50000.0]*%s,NBkg[25000.0,1.0,500000.0]*%s)",mSigFunct.c_str(),mBkgFunct.c_str());
  // } else {
  //   sprintf(funct,"SUM::sigMassPDF(NJpsi[5000.0,1.0,50000.0]*%s,NBkg[6000.0,1.0,500000.0]*%s)",mSigFunct.c_str(),mBkgFunct.c_str());
  // }


    
  string str("CBG");
  string dirPre2 = dirPre;
  size_t found = dirPre2.find(str);
  if (found!=string::npos)
    dirPre2.replace(found,str.length(),"CB");
  
  // if ( found==string::npos )
  //   ws->var("sigmaSig2")->setMax(0.060);

  // fix CB parameters to MC
  if (fixCBtoMC) {
    if (isPbPb) {
      if (prange=="6.5-30.0" && yrange=="0.0-2.4") {
	ws->var("alpha")->setVal(1.967);
	ws->var("enne")->setVal(1.380);
      }
      else if (prange=="6.5-30.0" && yrange=="0.0-1.6") {
	ws->var("alpha")->setVal(1.796);
	ws->var("enne")->setVal(1.479);
      }
      else if (prange=="3.0-30.0" && yrange=="1.6-2.4") {
	ws->var("alpha")->setVal(2.196);
	ws->var("enne")->setVal(1.280);
      }
      else if (prange=="3.0-6.5" && yrange=="1.6-2.4") {
	ws->var("alpha")->setVal(2.151);
	ws->var("enne")->setVal(1.380);
      }
      else if (prange=="6.5-30.0" && yrange=="1.6-2.4") {
	ws->var("alpha")->setVal(2.257);
	ws->var("enne")->setVal(1.170);
      }
      else {
	ws->var("alpha")->setVal(2.0);
	ws->var("enne")->setVal(1.4);
      }
    }
    else {
      if (prange=="6.5-30.0" && yrange=="0.0-2.4") {
	ws->var("alpha")->setVal(1.967);
	ws->var("enne")->setVal(1.380);
      }
      else if (prange=="6.5-30.0" && yrange=="0.0-1.6") {
	ws->var("alpha")->setVal(1.796);
	ws->var("enne")->setVal(1.479);
      }
      else if (prange=="3.0-30.0" && yrange=="1.6-2.4") {
	ws->var("alpha")->setVal(2.196);
	ws->var("enne")->setVal(1.280);
      }
      else if (prange=="3.0-6.5" && yrange=="1.6-2.4") {
	ws->var("alpha")->setVal(2.151);
	ws->var("enne")->setVal(1.380);
      }
      else if (prange=="6.5-30.0" && yrange=="1.6-2.4") {
	ws->var("alpha")->setVal(2.257);
	ws->var("enne")->setVal(1.170);
      }
      else {
	ws->var("alpha")->setVal(2.0);
	ws->var("enne")->setVal(1.4);
      }
    }
    // ws->var("alpha")->setConstant(kTRUE);
    // ws->var("enne")->setConstant(kTRUE);
  }

  /* 20131218
  // PbPb centrality dependence: CB + CBG
  if ( (crange.compare("0-100") || found!=string::npos) ) {
    // fix CB parameters within centrality bins to MB (use CB also for CBG)
    string inputFNcb;

    inputFNcb =  dirPre2 + "_rap" + yrange + "_pT" + prange + "_cent0-100_dPhi.txt";

    ifstream input;
    input.open(inputFNcb.c_str());
    if (!input.good()) { cout << "Failed to open: " <<inputFNcb << endl; return 1; }
    string tmp;
    double inputTmp[2] = {0};
    PARAM inputCB;
    input >> tmp >> inputTmp[0] >> inputTmp[1]; //NJpsi
    inputN[0] = inputTmp[0];  //NJpsi
    if (fitRatio) {
      input >> tmp >> inputTmp[0]; //NPsiP
      inputN[1] = inputTmp[0];  //NPsiP
      input >> tmp >> inputTmp[0] >> inputTmp[1]; //fracP
    }
    else {
      input >> tmp >> inputTmp[0] >> inputTmp[1]; //NPsiP
      inputN[1] = inputTmp[0];  //NPsiP
      input >> tmp >> inputTmp[0]; //fracP
    }

    input >> tmp >> inputTmp[0] >> inputTmp[1]; //Resolution
    input >> tmp >> inputTmp[0] >> inputTmp[1]; //NBkg
    inputN[2] = inputTmp[0];  //NBkg

    for (int p=0; p<11; p++) {   //Mass signal parameters
      input >> tmp >> inputTmp[0] >> inputTmp[1];
      cout << tmp << " " << inputTmp[0] << endl;
      if (!tmp.compare("coeffGaus")) inputCB.coeffGaus = inputTmp[0];
      else if (!tmp.compare("coeffExp")) inputCB.coeffExp = inputTmp[0];
      else if (!tmp.compare("coeffExp2")) inputCB.coeffExp2 = inputTmp[0];
      else if (!tmp.compare("cutx")) inputCB.cutx = inputTmp[0];
      else if (!tmp.compare("meanSig1")) inputCB.meanSig1 = inputTmp[0];
      else if (!tmp.compare("sigmaSig1")) inputCB.sigmaSig1 = inputTmp[0];
      else if (!tmp.compare("sigmaSig2")) inputCB.sigmaSig2 = inputTmp[0];
      else if (!tmp.compare("alpha")) inputCB.alpha = inputTmp[0];
      else if (!tmp.compare("enne")) inputCB.enne = inputTmp[0];
      //      else if (!tmp.compare("enneW")) inputCB.enneW = inputTmp[0];
    }


    ws->var("alpha")->setVal(inputCB.alpha);
    ws->var("enne")->setVal(inputCB.enne);
    //    ws->var("enneW")->setVal(inputCB.enneW);

    // ws->var("coeffExp")->setVal(inputCB.coeffExp);
    // ws->var("coeffExp2")->setVal(inputCB.coeffExp2);
    // ws->var("cutx")->setVal(inputCB.cutx);

    ws->var("alpha")->setConstant(kTRUE);
    ws->var("enne")->setConstant(kTRUE);
    //    ws->var("enneW")->setConstant(kTRUE);

    // ws->var("coeffExp")->setConstant(kTRUE);
    // ws->var("coeffExp2")->setConstant(kTRUE);
    // ws->var("cutx")->setConstant(kTRUE);

  }
  */

  /* 20131218
  if ( found!=string::npos ) {
    // fix cutx parameters within centrality bins to MB (use CB for CBG)
    string inputFNcb;

    inputFNcb =  dirPre2 + "_rap" + yrange + "_pT" + prange + "_cent" + crange + "_dPhi.txt";

    ifstream input;
    input.open(inputFNcb.c_str());
    if (!input.good()) { cout << "Failed to open: " <<inputFNcb << endl; return 1; }
    string tmp;
    double inputTmp[2] = {0};
    PARAM inputCB;
    input >> tmp >> inputTmp[0] >> inputTmp[1]; //NJpsi
    inputN[0] = inputTmp[0];  //NJpsi
    if (fitRatio) {
      input >> tmp >> inputTmp[0]; //NPsiP
      inputN[1] = inputTmp[0];  //NPsiP
      input >> tmp >> inputTmp[0] >> inputTmp[1]; //fracP
    }
    else {
      input >> tmp >> inputTmp[0] >> inputTmp[1]; //NPsiP
      inputN[1] = inputTmp[0];  //NPsiP
      input >> tmp >> inputTmp[0]; //fracP
    }

    input >> tmp >> inputTmp[0] >> inputTmp[1]; //Resolution
    input >> tmp >> inputTmp[0] >> inputTmp[1]; //NBkg
    inputN[2] = inputTmp[0];  //NBkg

    for (int p=0; p<11; p++) {   //Mass signal parameters
      input >> tmp >> inputTmp[0] >> inputTmp[1];
      cout << tmp << " " << inputTmp[0] << endl;
      if (!tmp.compare("cutx")) inputCB.cutx = inputTmp[0];
    }
    // ws->var("cutx")->setVal(inputCB.cutx);
    // ws->var("cutx")->setConstant(kTRUE);  
  }
  */
  
  // PbPb centrality dependence: CBG only
  if ( (crange.compare("0-100") && found!=string::npos) ) {
    // fix Gaussian parameters within centrality bins to MB (only needed for CBG)
    string inputFNcb;
    inputFNcb =  dirPre + "_" + mBkgFunct + "_rap" + yrange + "_pT" + prange + "_cent0-100_dPhi.txt";

    ifstream input;
    input.open(inputFNcb.c_str());
    if (!input.good()) { cout << "Failed to open: " <<inputFNcb << endl; return 1; }
    string tmp;
    double inputTmp[2] = {0};

    input >> tmp >> inputTmp[0] >> inputTmp[1]; //NJpsi
    inputN[0] = inputTmp[0];  //NJpsi
    if (fitRatio) {
      input >> tmp >> inputTmp[0]; //NPsiP
      inputN[1] = inputTmp[0];  //NPsiP
      input >> tmp >> inputTmp[0] >> inputTmp[1]; //fracP
    }
    else {
      input >> tmp >> inputTmp[0] >> inputTmp[1]; //NPsiP
      inputN[1] = inputTmp[0];  //NPsiP
      input >> tmp >> inputTmp[0]; //fracP
    }

    input >> tmp >> inputTmp[0] >> inputTmp[1]; //Resolution
    input >> tmp >> inputTmp[0] >> inputTmp[1]; //NBkg
    inputN[2] = inputTmp[0];  //NBkg

    for (int p=0; p<11; p++) {   //Mass signal parameters
      input >> tmp >> inputTmp[0] >> inputTmp[1];
      cout << tmp << " " << inputTmp[0] << endl;
      if (!tmp.compare("coeffGaus")) {
	inputCBG.coeffGaus = inputTmp[0];
	inputCBG.coeffGausErr = inputTmp[1];
      }
      else if (!tmp.compare("sigmaSig1")) {
	inputCBG.sigmaSig1 = inputTmp[0];
	inputCBG.sigmaSig1Err = inputTmp[1];
      }
      else if (!tmp.compare("sigmaSig2")) {
	inputCBG.sigmaSig2 = inputTmp[0];
	inputCBG.sigmaSig2Err = inputTmp[1];
      }
      else if (!tmp.compare("wideFactor")) {
	inputCBG.wideFactor = inputTmp[0];
	inputCBG.wideFactorErr = inputTmp[1];
      }
      else if (!tmp.compare("alpha")) inputCBG.alpha = inputTmp[0];
      else if (!tmp.compare("enne")) inputCBG.enne = inputTmp[0];
      //      else if (!tmp.compare("enneW")) inputCBG.enneW = inputTmp[0];
      else if (!tmp.compare("coeffExp")) inputCBG.coeffExp = inputTmp[0];
      else if (!tmp.compare("coeffExp2")) inputCBG.coeffExp2 = inputTmp[0];
      else if (!tmp.compare("cutx")) inputCBG.cutx = inputTmp[0];
      else if (!tmp.compare("bkgMean")) inputCB.bkgMean = inputTmp[0];
      else if (!tmp.compare("bkgSigma")) inputCB.bkgSigma = inputTmp[0];
      else if (!tmp.compare("coeffPol1")) inputCB.coeffPol1 = inputTmp[0];
      else if (!tmp.compare("coeffPol2")) inputCB.coeffPol2 = inputTmp[0];
      else if (!tmp.compare("coeffPol3")) inputCB.coeffPol3 = inputTmp[0];
      else if (!tmp.compare("coeffPol4")) inputCB.coeffPol4 = inputTmp[0];
      else if (!tmp.compare("coeffPol5")) inputCB.coeffPol5 = inputTmp[0];
    }

    ws->var("sigmaSig1")->setVal(inputCBG.sigmaSig1);
    //    ws->var("sigmaSig2")->setVal(inputCBG.sigmaSig2);
    ws->var("wideFactor")->setVal(inputCBG.wideFactor);
    ws->var("coeffGaus")->setVal(inputCBG.coeffGaus);
    // ws->var("alpha")->setVal(inputCBG.alpha);
    // ws->var("enne")->setVal(inputCBG.enne);
    // ws->var("enneW")->setVal(inputCBG.enneW);

    //    ws->var("coeffExp")->setVal(inputCBG.coeffExp);
    //    ws->var("coeffExp2")->setVal(inputCBG.coeffExp2);
    //    ws->var("cutx")->setVal(inputCBG.cutx);

    // ws->var("bkgMean")->setVal(inputCBG.bkgMean);
    // ws->var("bkgSigma")->setVal(inputCBG.bkgSigma);

    ws->var("coeffPol1")->setVal(inputCBG.coeffPol1);
    ws->var("coeffPol2")->setVal(inputCBG.coeffPol2);
    ws->var("coeffPol3")->setVal(inputCBG.coeffPol3);
    ws->var("coeffPol4")->setVal(inputCBG.coeffPol4);
    //    ws->var("coeffPol5")->setVal(inputCBG.coeffPol5);

    //    ws->var("wideFactor")->setConstant(kTRUE);
    //    ws->var("coeffGaus")->setConstant(kTRUE);
    // ws->var("alpha")->setConstant(kTRUE);
    // ws->var("enne")->setConstant(kTRUE);
    // ws->var("enneW")->setConstant(kTRUE);

    // ws->var("coeffExp")->setConstant(kTRUE);
    // ws->var("coeffExp2")->setConstant(kTRUE);
    // ws->var("cutx")->setConstant(kTRUE);

    // ws->var("bkgMean")->setConstant(kTRUE);
    // ws->var("bkgSigma")->setConstant(kTRUE);

  }

  // ws->var("Jpsi_Mass")->setRange("jpsi",2.6,3.55);
  // ws->var("Jpsi_Mass")->setRange("psip",3.4,4.2);

  /* 20131218
  if (fitSubRange) {
    // ws->var("fracP")->setVal(0.0);
    // ws->var("fracP")->setConstant(kTRUE);
    fitM = ws->pdf("jpsiMassPDF")->fitTo(*redDataCut,Extended(1),Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(8),Range("jpsi"));

    cout << "Fitting psi(2S)" << endl;
    ws->var("alpha")->setConstant(kTRUE);
    ws->var("enne")->setConstant(kTRUE);
    //    ws->var("enneW")->setConstant(kTRUE);
    ws->var("sigmaSig2")->setConstant(kTRUE);
    ws->var("sigmaSig1")->setConstant(kTRUE);
    ws->var("meanSig1")->setConstant(kTRUE);
    ws->var("coeffGaus")->setConstant(kTRUE);
    ws->var("NJpsi")->setConstant(kTRUE);

    fitP = ws->pdf("psipMassPDF")->fitTo(*redDataCut,Extended(1),Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(8),Range("psip"));
    resultF.cd();
    cout << "fitM->Write(jpsiMassPDF):" << fitM->Write("jpsiMassPDF") << endl;
    cout << "fitP->Write(psipMassPDF):" << fitP->Write("psipMassPDF") << endl;
  }
  else {
    if (prefitMass) {
      // ws->var("alpha")->setConstant(false);
      // ws->var("enne")->setConstant(false);
      // ws->var("enneW")->setConstant(false);
      // ws->var("fracP")->setVal(0.0);
      // ws->var("fracP")->setConstant(kTRUE);
      // ws->var("NBkg")->setVal(0.0);
      // ws->var("NBkg")->setConstant(true);
      fitM_pre = ws->pdf("jpsiMassPDF")->fitTo(*redDataCut,Extended(1),Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(8),Range("jpsi"));
      fitM_pre->Print("v");
      string resultFN_pre = dirPre + "_rap" + yrange + "_pT" + prange + "_cent" + crange + "_dPhi" + phirange + "_preFitResult.root";
      TFile resultPreF(resultFN_pre.c_str(),"RECREATE");
      resultPreF.cd(); cout << "fitM_pre->Write(jpsiMassPDF):" << fitM_pre->Write("jpsiMassPDF") << endl;resultPreF.Close();
      // ws->var("coeffExp")->setVal(ws->var("coeffExp")->getVal());
      // ws->var("enne")->setConstant(true);
      // ws->var("enneW")->setConstant(true);
      // ws->var("alpha")->setConstant(true);
      // ws->var("meanSig1")->setConstant(kTRUE);
      // ws->var("meanSig2")->setConstant(kTRUE);

      // if ( (!crange.compare("0-100") && found==string::npos) )
      // 	ws->var("enne")->setConstant(true);
      // 	ws->var("enneW")->setConstant(true);

	 if (!crange.compare("0-100")){
	 inputCBG.sigmaSig1Err = ws->var("sigmaSig1")->getError();
	 inputCBG.coeffGausErr = ws->var("coeffGaus")->getError();

	 ws->var("sigmaSig1")->setConstant(kTRUE);
	 // ws->var("sigmaSig2")->setConstant(kTRUE);
	 ws->var("coeffGaus")->setConstant(kTRUE);
	 }
    }
    */
    //    ws->var("NJpsi")->setVal(0.0);ws->var("NJpsi")->setConstant(true);
  fitM = ws->pdf("sigMassPDF")->fitTo(*redDataCut,Extended(1),Minos(1),Save(1),SumW2Error(kTRUE),NumCPU(8));
  resultF.cd(); cout << "fitM->Write(sigMassPDF):" << fitM->Write("sigMassPDF") << endl;
    //  }

  // *** Draw mass plot before do ctau fit
  RooPlot *mframe_wob = ws->var("Jpsi_Mass")->frame();
  redDataCut->plotOn(mframe_wob,DataError(RooAbsData::SumW2),XErrorSize(0),MarkerSize(0.9),Binning(rbm));
  titlestr = "2D fit for" + partTit + "muons (mass projection), p_{T} = " + prange  + "_dPhi" + phirange + " GeV/c and |y| = " + yrange;
  mframe_wob->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  mframe_wob->GetXaxis()->CenterTitle(1);
  double max = mframe_wob->GetMaximum() * 1.3;
  double min = 0.0;
  mframe_wob->SetMaximum(max);
  min = ws->var("NBkg")->getVal()/(double(nbins)) * 0.7;
  RooHist *hpull_jpsi;
  RooHist *hpull_psip;
  RooHist *hpull; 
  if (fitSubRange) {
    ws->pdf("jpsiMassPDF")->plotOn(mframe_wob,DrawOption("F"),FillColor(kBlack),FillStyle(3354));

    //    ws->pdf("jpsiMassPDF")->plotOn(mframe_wob,Components(mBkgFunct.c_str()),DrawOption("F"),FillColor(kAzure-9),FillStyle(1001));
    ws->pdf("jpsiMassPDF")->plotOn(mframe_wob,Components(mBkgFunct.c_str()),LineColor(kBlue),LineStyle(7),LineWidth(5));
    ws->pdf("jpsiMassPDF")->plotOn(mframe_wob,LineColor(kBlack),LineWidth(2));
    hpull_jpsi = mframe_wob->pullHist(0,0,true); hpull_jpsi->SetName("hpullhist");

    ws->pdf("psipMassPDF")->plotOn(mframe_wob,DrawOption("F"),FillColor(kGreen+2),FillStyle(3354));
    
    //    ws->pdf("psipMassPDF")->plotOn(mframe_wob,Components(mBkgFunctP.c_str()),DrawOption("F"),FillColor(kRed-9),FillStyle(1001));
    
    ws->pdf("psipMassPDF")->plotOn(mframe_wob,Components(mBkgFunctP.c_str()),LineColor(kRed),LineStyle(7),LineWidth(5));
    ws->pdf("psipMassPDF")->plotOn(mframe_wob,LineColor(kGreen+2),LineWidth(2));
    hpull_psip = mframe_wob->pullHist(0,0,true); hpull_psip->SetName("hpullhistP");
    hpull_psip->SetMarkerColor(kGreen+2);
  }
  else {
    ws->pdf("sigMassPDF")->plotOn(mframe_wob,DrawOption("F"),FillColor(kBlack),FillStyle(3354),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));

    ws->pdf("sigMassPDF")->plotOn(mframe_wob,Components(mBkgFunct.c_str()),DrawOption("F"),FillColor(kAzure-9),FillStyle(1001),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
    
    ws->pdf("sigMassPDF")->plotOn(mframe_wob,Components(mBkgFunct.c_str()),LineColor(kBlue),LineStyle(7),LineWidth(5),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
    if (!isPaper && found!=string::npos) { 
      ws->pdf("sigMassPDF")->plotOn(mframe_wob,Components((mBkgFunct+",signalG2,signalG2P").c_str()),LineColor(kGreen+2),LineStyle(kDashed),LineWidth(3));
      //      ws->pdf("sigMassPDF")->plotOn(mframe_wob,Components((mBkgFunct+",signalG2P").c_str()),LineColor(kOrange+2),LineStyle(kDashed),LineWidth(3));
    }

    ws->pdf("sigMassPDF")->plotOn(mframe_wob,LineColor(kBlack),LineWidth(2),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
    
    hpull = mframe_wob->pullHist(0,0,true); hpull->SetName("hpullhist");
    
    // if (prefitMass) {
    //   ws->pdf("jpsiMassPDF")->plotOn(mframe_wob,DrawOption("F"),FillColor(kBlack),FillStyle(3354));

    //   ws->pdf("jpsiMassPDF")->plotOn(mframe_wob,Components(mBkgFunct.c_str()),DrawOption("F"),FillColor(kAzure-9),FillStyle(1001));
    
    //   ws->pdf("jpsiMassPDF")->plotOn(mframe_wob,Components(mBkgFunct.c_str()),LineColor(kBlue),LineStyle(7),LineWidth(5));
    //   ws->pdf("jpsiMassPDF")->plotOn(mframe_wob,LineColor(kRed),LineWidth(2));
    // }

  }
  redDataCut->plotOn(mframe_wob,DataError(RooAbsData::SumW2),XErrorSize(0),MarkerSize(0.9),Binning(rbm));

  TH1 *hdata = redDataCut->createHistogram("hdata",*ws->var("Jpsi_Mass"),Binning(rbm));
  // *** Calculate chi2/nDof for mass fitting

  unsigned int nBins = hdata->GetNbinsX();
  double Chi2 = 0;
  int nFullBinsPull = 0;
  double *ypull;
  double *ypull_psip;
  if (fitSubRange) {
    ypull = hpull_jpsi->GetY();
    ypull_psip = hpull_psip->GetY();
  }
  else
    ypull = hpull->GetY();
  
  double Chi2P = 0;
  int nFullBinsPullP = 0;

  if (fitSubRange) {
    int minBins = hdata->GetXaxis()->FindBin(3.4);
    for (unsigned int i=minBins; i < nBins; i++) {
      if (hdata->GetBinContent(i+1) == 0) continue;
      nFullBinsPullP++;
      Chi2P = Chi2P + pow(ypull_psip[i-minBins],2);
    }

    nBins = hdata->GetXaxis()->FindBin(3.5);

    for (unsigned int i=0; i < nBins; i++) {
      if (hdata->GetBinContent(i+1) == 0) continue;
      nFullBinsPull++;
      Chi2 = Chi2 + pow(ypull[i],2);
    }
  }
  else {
    for (unsigned int i=0; i < nBins; i++) {
      if (hdata->GetBinContent(i+1) == 0) continue;
      nFullBinsPull++;
      Chi2 = Chi2 + pow(ypull[i],2);
    }
  }

  double UnNormChi2 = Chi2;
  int nFitParam = fitM->floatParsFinal().getSize();
  int Dof = nFullBinsPull - nFitParam;
  if (Dof!=0)
    Chi2 /= Dof;

  double theNLL=0;
  theNLL = fitM->minNll();

  double chi2FromRoo = mframe_wob->chiSquare(nFitParam);

  double UnNormChi2P = Chi2P;
  int nFitParamP = 0;
  int DofP = 1;
  if (fitSubRange) {
    nFitParamP = fitP->floatParsFinal().getSize();
    DofP = nFullBinsPullP - nFitParamP;
    if (DofP!=0)
      Chi2P /= DofP;
  }

  RooPlot* mframepull =  ws->var("Jpsi_Mass")->frame(Title("Pull")) ;
  mframepull->GetYaxis()->SetTitle("Pull");
  mframepull->SetLabelSize(0.08,"XYZ");
  mframepull->SetTitleSize(0.1,"XYZ");
  mframepull->SetTitleOffset(0.55,"Y");
  if (fitSubRange) {
    mframepull->addPlotable(hpull_jpsi,"PX") ;
    mframepull->addPlotable(hpull_psip,"PX") ;
  }
  else 
    mframepull->addPlotable(hpull,"PX") ;

  if ( fabs(mframepull->GetMaximum()) > fabs(mframepull->GetMinimum()) )
    mframepull->SetMinimum(-(mframepull->GetMaximum())); 
  else
    mframepull->SetMaximum(-(mframepull->GetMinimum())); 
  mframepull->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c^{2})");
  mframepull->GetXaxis()->CenterTitle(1);

  // *** Check in narrower signal region NSig
  //  ws->var("Jpsi_Mass")->setRange("sigpeak",2.9,3.3);

  //  RooAbsReal *inteAll = ws->pdf(mJpsiFunct.c_str())->createIntegral(RooArgSet(*ws->var("Jpsi_Mass")),NormSet(RooArgSet(*ws->var("Jpsi_Mass"))));
  //  RooAbsReal *inteSig = ws->pdf(mSigFunct.c_str())->createIntegral(RooArgSet(*ws->var("Jpsi_Mass")),Range("sigpeak"),NormSet(RooArgSet(*ws->var("Jpsi_Mass"))),Range("sigpeak"));

  // mframe_wob->SetMinimum(0.03*max);
  // mframe_wob->SetMaximum(5.0*max);
  TCanvas c00; c00.cd();

  TPad *pad1 = new TPad("pad1","This is pad1",0.05,0.35,0.95,0.97);
  pad1->SetBottomMargin(0); 
  TPad *pad2 = new TPad("pad2","This is pad2",0.05,0.03,0.95,0.35);
  pad2->SetTopMargin(0);  pad2->SetBottomMargin(0.24); 
  if (!isPaper) {
    pad1->Draw();
    pad2->Draw();
    pad2->cd(); mframepull->Draw();
    pad1->cd();
  }

  // if (!isPbPb)
  //   ws->pdf("sigMassPDF_mix")->plotOn(mframe_wob,LineColor(kRed),LineStyle(2),LineWidth(4),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));

  if (logScale) {
    if (isPaper)
      c00.SetLogy(1);
    else
      pad1->SetLogy(1);

    mframe_wob->SetMinimum(0.2*min);
    if (yrange == "0.0-1.6" && isPbPb) {
      mframe_wob->SetMaximum(2.5*max);
      mframe_wob->SetMinimum(0.8*min);
    }
    else if (yrange == "1.6-2.4" && isPbPb) {
      mframe_wob->SetMaximum(1.5*max);
      mframe_wob->SetMinimum(0.9*min);
    }
    else 
      mframe_wob->SetMaximum(25*max);
  }
  else if (zoom) {
    if (min>0)
      mframe_wob->SetMinimum(0.5*min);

    if (yrange == "1.6-2.4") {
      if (prange == "3.0-6.5" && crange == "0-20")
       	mframe_wob->SetMaximum(max);
      else if (prange == "3.0-30.0" && crange == "20-40")
       	mframe_wob->SetMaximum(max*0.3);
      else if (prange == "3.0-30.0" && crange=="40-100")
       	mframe_wob->SetMaximum(max*0.2);
      else 
	mframe_wob->SetMaximum(max*0.5);
    }
    else
      mframe_wob->SetMaximum(max*0.15);
  }

  mframe_wob->Draw();
  t->SetTextSize(0.05);
  //  t->DrawLatex(0.17,0.90,"CMS Preliminary");
  if (isPbPb)
    t->DrawLatex(0.17,0.90,"CMS PbPb #sqrt{s_{NN}} = 2.76 TeV");
  else 
    t->DrawLatex(0.17,0.90,"CMS pp #sqrt{s} = 2.76 TeV");
  t->SetTextSize(0.04);
  if (isPbPb)
    t->DrawLatex(0.17,0.83,"L_{int} = 150 #mub^{-1}");
  else
    t->DrawLatex(0.17,0.83,"L_{int} = 5.3 pb^{-1}");
  t->SetTextSize(0.035);
  if (isPbPb) { 
    if (ymin==0.0)
      sprintf(reduceDS,"%.0f-%.0f%%, |y| < %.1f",cmin,cmax,ymax);
    else
      sprintf(reduceDS,"%.0f-%.0f%%, %0.1f < |y| < %.1f",cmin,cmax,ymin,ymax);
  }
  else {
    if (ymin==0.0)
      sprintf(reduceDS,"|y| < %.1f",ymax);
    else
      sprintf(reduceDS,"%0.1f < |y| < %.1f",ymin,ymax);
  }
  // if (zoom)
  //   t->DrawLatex(0.56,0.70,reduceDS);
  // else
  //   t->DrawLatex(0.56,0.50,reduceDS);
    t->DrawLatex(0.17,0.78,reduceDS);
  if (pmin==6.5)
    sprintf(reduceDS,"%.1f < p_{T} < %.0f GeV/c",pmin,pmax);
  else if (pmax==6.5)
    sprintf(reduceDS,"%.0f < p_{T} < %.1f GeV/c",pmin,pmax);
  else
    sprintf(reduceDS,"%.0f < p_{T} < %.0f GeV/c",pmin,pmax);
  
  // if (zoom)
  //   t->DrawLatex(0.56,0.65,reduceDS);
  // else
  //   t->DrawLatex(0.56,0.45,reduceDS);
  t->DrawLatex(0.17,0.73,reduceDS);
  // sprintf(reduceDS,"%.2f < |#phi_{J/#psi}-#Psi_{RP}| < %.2f",psmin,psmax);
  // t->DrawLatex(0.17,0.59,reduceDS);
  // sprintf(reduceDS,"EP: %s",rpmethod.c_str());
  // t->DrawLatex(0.17,0.54,reduceDS);
  //  t->SetTextSize(0.04);

  if (!isPaper) {
    sprintf(reduceDS,"Min. NLL = %0.1f",theNLL);
    t->DrawLatex(0.17,0.68,reduceDS);
  }

  if (fitSubRange) {
    sprintf(reduceDS,"#chi^{2}/dof (J/#psi) = %0.1f/%d   #chi^{2}/dof (#psi(2S)) = %0.1f/%d",UnNormChi2,Dof,UnNormChi2P,DofP);
    t->DrawLatex(0.45,0.90,reduceDS);
  }
  else {
    sprintf(reduceDS,"#chi^{2}/dof = %0.1f/%d (Test: %0.1f)",UnNormChi2,Dof,chi2FromRoo);
    if (!isPaper)
      t->DrawLatex(0.62,0.90,reduceDS);
    sprintf(reduceDS,"p-value = %0.4f (%0.4f)",TMath::Prob(UnNormChi2,Dof),TMath::Prob(chi2FromRoo,Dof));
    if (!isPaper)
      t->DrawLatex(0.62,0.85,reduceDS);
  }

  if (ws->var("NJpsi")->hasAsymError() && abs(-1.0*ws->var("NJpsi")->getErrorLo()/ws->var("NJpsi")->getErrorHi() - 1)>0.1)
    sprintf(reduceDS,"N_{J/#psi} = %0.0f^{+%0.0f}_{%0.0f}",ws->var("NJpsi")->getVal(),ws->var("NJpsi")->getErrorHi(),ws->var("NJpsi")->getErrorLo());
  else
    sprintf(reduceDS,"N_{J/#psi} = %0.0f #pm %0.0f",ws->var("NJpsi")->getVal(),ws->var("NJpsi")->getError());
  t->DrawLatex(0.62,0.80,reduceDS);
  // sprintf(reduceDS,"N_{#psi(2S)} = %0.0f #pm %0.0f",ws->var("NPsiP")->getVal(),ws->var("NPsiP")->getError());
  // t->DrawLatex(0.62,0.80,reduceDS);

  if (fitRatio) {
    if (ws->var("fracP")->hasAsymError() && abs(-1.0*ws->var("fracP")->getErrorLo()/ws->var("fracP")->getErrorHi() - 1)>0.1)
      sprintf(reduceDS,"R_{#psi(2S)} = %0.3f^{+%0.3f}_{%0.3f}",ws->var("fracP")->getVal(),ws->var("fracP")->getErrorHi(),ws->var("fracP")->getErrorLo());
    else
      sprintf(reduceDS,"R_{#psi(2S)} = %0.3f #pm %0.3f",ws->var("fracP")->getVal(),ws->var("fracP")->getError());
    t->DrawLatex(0.62,0.75,reduceDS);
    if (ws->var("fracP")->hasAsymError() && abs(-1.0*ws->var("fracP")->getErrorLo()/ws->var("fracP")->getErrorHi() - 1)>0.1)
      sprintf(reduceDS,"N_{#psi(2S)} = %0.1f^{+%0.1f}_{%0.1f}",ws->function("NPsiP")->getVal(),ws->var("fracP")->getErrorHi()*ws->var("NJpsi")->getVal(),ws->var("fracP")->getErrorLo()*ws->var("NJpsi")->getVal());
    else
      sprintf(reduceDS,"N_{#psi(2S)} = %0.1f #pm %0.1f",ws->function("NPsiP")->getVal(),ws->var("fracP")->getError()*ws->var("NJpsi")->getVal());
    t->DrawLatex(0.62,0.70,reduceDS);
  }
  else {
    if (ws->var("NPsiP")->hasAsymError() && abs(-1.0*ws->var("NPsiP")->getErrorLo()/ws->var("NPsiP")->getErrorHi() - 1)>0.1)
      sprintf(reduceDS,"N_{#psi(2S)} = %0.0f^{+%0.0f}_{%0.0f}",ws->var("NPsiP")->getVal(),ws->var("NPsiP")->getErrorHi(),ws->var("NPsiP")->getErrorLo());
    else
      sprintf(reduceDS,"N_{#psi(2S)} = %0.1f #pm %0.1f",ws->var("NPsiP")->getVal(),ws->var("NPsiP")->getError());
    t->DrawLatex(0.62,0.75,reduceDS);
  }

  double coeffGaus = ws->var("coeffGaus")->getVal();
  double sigmaSig1 = ws->var("sigmaSig1")->getVal();
  double sigmaSig2 = ws->function("sigmaSig2")->getVal();
  double ErrcoeffGaus = ws->var("coeffGaus")->getError();
  double ErrsigmaSig1 = ws->var("sigmaSig1")->getError();
  double ErrsigmaSig2 = sigmaSig1*ws->var("wideFactor")->getError();

  if (ErrcoeffGaus == 0.0) ErrcoeffGaus = inputCBG.coeffGausErr;
  if (ErrsigmaSig1 == 0.0) ErrsigmaSig1 = inputCBG.sigmaSig1Err;
  if (ErrsigmaSig2 == 0.0) ErrsigmaSig2 = inputCBG.sigmaSig2Err;

  if (!isPaper) {
    if (ErrsigmaSig1<0.0005) 
      sprintf(reduceDS,"#sigma_{CB} = (%0.0f #pm %0.0f) MeV/c^{2}",sigmaSig1*1000.0,1.0);
    else
      sprintf(reduceDS,"#sigma_{CB} = (%0.0f #pm %0.0f) MeV/c^{2}",sigmaSig1*1000.0,ErrsigmaSig1*1000.0);
    t->DrawLatex(0.62,0.65,reduceDS);
    // if (ErrsigmaSig2<0.0005) 
    //   sprintf(reduceDS,"#sigma_{G} = (%0.0f #pm %0.0f) MeV/c^{2}",sigmaSig2*1000.0,1.0);
    // else
    if  (found!=string::npos) {
      sprintf(reduceDS,"#sigma_{G} = %0.0f MeV/c^{2}",sigmaSig2*1000.0);
      t->DrawLatex(0.62,0.60,reduceDS);
      if (ws->var("wideFactor")->hasAsymError() && abs(-1.0*ws->var("wideFactor")->getErrorLo()/ws->var("wideFactor")->getErrorHi() - 1)>0.1)
	sprintf(reduceDS,"n_{G} = %0.2f^{+%0.2f}_{%0.2f}",ws->var("wideFactor")->getVal(),ws->var("wideFactor")->getErrorHi(),ws->var("wideFactor")->getErrorLo());
      else
	sprintf(reduceDS,"n_{G} = %0.2f #pm %0.2f",ws->var("wideFactor")->getVal(),ws->var("wideFactor")->getError());
      t->DrawLatex(0.62,0.55,reduceDS);
    }
  }

  double resol = sigmaSig1;
  double Errresol = ErrsigmaSig1;

  if  (found!=string::npos) {
    resol = sqrt( coeffGaus*pow(sigmaSig2,2) + (1-coeffGaus)*pow(sigmaSig1,2) );
    Errresol = (0.5/resol) * sqrt( pow(sigmaSig2*coeffGaus*ErrsigmaSig2,2) +
				   pow(sigmaSig1*(1-coeffGaus)*ErrsigmaSig1,2) +
				   pow(0.5*(pow(sigmaSig2,2)-pow(sigmaSig1,2))*ErrcoeffGaus,2) );
  }


  if (!isPaper) {
    if  (found!=string::npos) {
      if (Errresol<0.0005) 
	sprintf(reduceDS,"#sigma_{CB+G} = (%0.0f #pm %0.0f) MeV/c^{2}",resol*1000.0,1.0);
      else
	sprintf(reduceDS,"#sigma_{CB+G} = (%0.0f #pm %0.0f) MeV/c^{2}",resol*1000.0,Errresol*1000.0);

      t->DrawLatex(0.62,0.50,reduceDS);
      if (ws->var("alpha")->isConstant())
	sprintf(reduceDS,"#alpha = %0.2f (fixed)",ws->var("alpha")->getVal());
      else if (ws->var("alpha")->hasAsymError() && abs(-1.0*ws->var("alpha")->getErrorLo()/ws->var("alpha")->getErrorHi() - 1)>0.1)
	sprintf(reduceDS,"#alpha = %0.2f^{+%0.2f}_{%0.2f}",ws->var("alpha")->getVal(),ws->var("alpha")->getErrorHi(),ws->var("alpha")->getErrorLo());
      else
	sprintf(reduceDS,"#alpha = %0.2f #pm %0.2f",ws->var("alpha")->getVal(),ws->var("alpha")->getError());
      t->DrawLatex(0.62,0.45,reduceDS);
      if (ws->var("enne")->isConstant())
	sprintf(reduceDS,"n = %0.2f (fixed)",ws->var("enne")->getVal());
      else if (ws->var("enne")->hasAsymError() && abs(-1.0*ws->var("enne")->getErrorLo()/ws->var("enne")->getErrorHi() - 1)>0.1)
	sprintf(reduceDS,"n = %0.2f^{+%0.2f}_{%0.2f}",ws->var("enne")->getVal(),ws->var("enne")->getErrorHi(),ws->var("enne")->getErrorLo());
      else
	sprintf(reduceDS,"n = %0.2f #pm %0.2f",ws->var("enne")->getVal(),ws->var("enne")->getError());
      t->DrawLatex(0.62,0.40,reduceDS);
      if (ws->var("coeffGaus")->hasAsymError() && abs(-1.0*ws->var("coeffGaus")->getErrorLo()/ws->var("coeffGaus")->getErrorHi() - 1)>0.1)
	sprintf(reduceDS,"f_{G} = %0.3f^{+%0.3f}_{%0.3f}",ws->var("coeffGaus")->getVal(),ws->var("coeffGaus")->getErrorHi(),ws->var("coeffGaus")->getErrorLo());
      else
	sprintf(reduceDS,"f_{G} = %0.2f #pm %0.2f",ws->var("coeffGaus")->getVal(),ws->var("coeffGaus")->getError());
      t->DrawLatex(0.62,0.35,reduceDS);
    }
    else {
      if (ws->var("alpha")->isConstant())
	sprintf(reduceDS,"#alpha = %0.2f (fixed)",ws->var("alpha")->getVal());
      else
	sprintf(reduceDS,"#alpha = %0.2f #pm %0.2f",ws->var("alpha")->getVal(),ws->var("alpha")->getError());
      t->DrawLatex(0.62,0.60,reduceDS);
      if (ws->var("enne")->isConstant())
	sprintf(reduceDS,"n = %0.2f (fixed)",ws->var("enne")->getVal());
      else
	sprintf(reduceDS,"n = %0.2f #pm %0.2f",ws->var("enne")->getVal(),ws->var("enne")->getError());
      t->DrawLatex(0.62,0.55,reduceDS);
    }
  }
  //  sprintf(reduceDS,"#sigma = (%0.0f #pm %0.0f) MeV/c^{2}",ws->var("sigmaSig2")->getVal()*1000,ws->var("sigmaSig2")->getError()*1000);
  // t->DrawLatex(0.62,0.55,reduceDS);

  TLegend *leg1;
  //  if (isPbPb) {
  if (fitSubRange)
    leg1 = new TLegend(0.63,0.53,0.92,0.72,NULL,"brNDC");
  else
    leg1 = new TLegend(0.63,0.53,0.92,0.68,NULL,"brNDC");
  //  }
  // else
  //   leg1 = new TLegend(0.63,0.53,0.92,0.72,NULL,"brNDC");
  leg1->SetFillStyle(0); leg1->SetBorderSize(0); leg1->SetShadowColor(0); leg1->SetMargin(0.2);
  leg1->AddEntry(gfake1,"data","p");
  if (fitSubRange) {
    leg1->AddEntry(&hfake21,"total fit J/#psi","lf");
    leg1->AddEntry(&hfake11,"background J/#psi","l");
    leg1->AddEntry(&hfake23,"total fit #psi(2S)","lf");
    leg1->AddEntry(&hfake22,"background #psi(2S)","l");
  }
  else {
    leg1->AddEntry(&hfake21,"total fit","lf");
    leg1->AddEntry(&hfake11,"background","lf");
  }
  // if (!isPbPb)
  //   leg1->AddEntry(&hfake22,"with R_{#psi(2S)}^{0-20%}(PbPb)","lf");

  if (!zoom && isPaper)
    leg1->Draw("same");

  string yrange_str, prange_str;
  if (yrange == "0.0-2.4")
    yrange_str = "0-24";
  else if (yrange == "0.0-1.2")
    yrange_str = "0-12";
  else if (yrange == "0.0-1.6")
    yrange_str = "0-16";
  else if (yrange == "1.2-1.6")
    yrange_str = "12-16";
  else if (yrange == "1.6-2.4")
    yrange_str = "16-24";
  else
    yrange_str = yrange;

  if (prange == "0.0-40.0")
    prange_str = "0-40";
  else if (prange == "0.0-30.0")
    prange_str = "0-30";
  else if (prange == "3.0-40.0")
    prange_str = "3-40";
  else if (prange == "3.0-30.0")
    prange_str = "3-30";
  else if (prange == "6.5-40.0")
    prange_str = "65-40";
  else if (prange == "6.5-30.0")
    prange_str = "65-30";
  else if (prange == "3.0-6.5")
    prange_str = "3-65";
  else
    prange_str = prange;

  titlestr = dirPre + "_" + mBkgFunct + "_rap" + yrange_str  + "_pT" + prange_str + "_cent" + crange + "_massfit_wob.pdf";
  c00.SaveAs(titlestr.c_str());


  Double_t NJpsi_fin = ws->var("NJpsi")->getVal();
  Double_t ErrNJpsi_fin = ws->var("NJpsi")->getError();
  Double_t NPsiP_fin = 0.0;
  Double_t ErrNPsiP_fin = 0.0;
  Double_t fracP_fin = 0.0;
  Double_t ErrFracP_fin = 0.0;
  if (fitRatio) {
    fracP_fin = ws->var("fracP")->getVal();
    ErrFracP_fin = ws->var("fracP")->getError();
  }
  else {
    NPsiP_fin = ws->var("NPsiP")->getVal();
    ErrNPsiP_fin = ws->var("NPsiP")->getError();
  }
  Double_t NBkg_fin = ws->var("NBkg")->getVal();
  Double_t ErrNBkg_fin = ws->var("NBkg")->getError();

  // To check values of fit parameters
  cout << endl << "J/psi yields:" << endl;
  cout << "NJpsi:       Fit: "  << NJpsi_fin << " +/- " << ErrNJpsi_fin << endl;
  if (fitRatio) {
    cout << "R_psi(2S)    Fit: "  << fracP_fin << " +/- " << ErrFracP_fin << endl;
    cout << "NPsiP:       Fit: "  << fracP_fin*NJpsi_fin << " +/- " << ErrFracP_fin*NJpsi_fin << endl;
  }
  else
    cout << "NPsiP:       Fit: "  << NPsiP_fin << " +/- " << ErrNPsiP_fin << endl;

  titlestr = dirPre + "_" + mBkgFunct + "_rap" + yrange + "_pT" + prange + "_cent" + crange + "_dPhi" + phirange + ".txt";

  ofstream outputFile(titlestr.c_str());
  if (!outputFile.good()) {cout << "Fail to open result file." << endl; return 1;}
  outputFile
    << "NJpsi "        << NJpsi_fin                       << " " << ErrNJpsi_fin << "\n";
  if (fitRatio) {
    outputFile
      << "NPsiP "        << NJpsi_fin*fracP_fin << "\n"
      << "R_psi(2S) "    << fracP_fin            << " " << ErrFracP_fin << "\n";
  }
  else {
    outputFile
      << "NPsiP "        << NPsiP_fin            << " " << ErrNPsiP_fin << "\n"
      << "R_psi(2S) "    << NPsiP_fin/NJpsi_fin << "\n";
  }
  outputFile
    << "Resolution "   << resol*1000.0                      << " " << Errresol*1000.0 << "\n"
    << "NBkg "         << NBkg_fin                          << " " << ErrNBkg_fin << "\n";
  if (!mBkgFunct.compare("expFunct")) {
    outputFile
      << "coeffExp "      << ws->var("coeffExp")->getVal()      << " " << ws->var("coeffExp")->getError() << "\n";
  }
  if (!mBkgFunct.compare("twoExpFunct")) {
    outputFile
      << "coeffExp "      << ws->var("coeffExp")->getVal()      << " " << ws->var("coeffExp")->getError() << "\n"
      << "coeffExp2 "     << ws->var("coeffExp2")->getVal()     << " " << ws->var("coeffExp2")->getError() << "\n"
      << "cutx "          << ws->var("cutx")->getVal()          << " " << ws->var("cutx")->getError() << "\n";
  }
  if (!mBkgFunct.compare("bkgGaus")) {
    outputFile
      << "meanBkg "      << ws->var("meanBkg")->getVal()      << " " << ws->var("meanBkg")->getError() << "\n"
      << "sigmaBkg "     << ws->var("sigmaBkg")->getVal()     << " " << ws->var("sigmaBkg")->getError() << "\n";
  }
  if (!mBkgFunct.compare("pol1")) {
    outputFile
      << "coeffPol1 "     << ws->var("coeffPol1")->getVal()     << " " << ws->var("coeffPol1")->getError() << "\n";
  }
  if (!mBkgFunct.compare("pol2")) {
    outputFile
      << "coeffPol1 "     << ws->var("coeffPol1")->getVal()     << " " << ws->var("coeffPol1")->getError() << "\n"
      << "coeffPol2 "     << ws->var("coeffPol2")->getVal()     << " " << ws->var("coeffPol2")->getError() << "\n";
  }
  if (!mBkgFunct.compare("pol3")) {
    outputFile
      << "coeffPol1 "     << ws->var("coeffPol1")->getVal()     << " " << ws->var("coeffPol1")->getError() << "\n"
      << "coeffPol2 "     << ws->var("coeffPol2")->getVal()     << " " << ws->var("coeffPol2")->getError() << "\n"
      << "coeffPol3 "     << ws->var("coeffPol3")->getVal()     << " " << ws->var("coeffPol3")->getError() << "\n";
  }
  if (!mBkgFunct.compare("pol4")) {
    outputFile
      << "coeffPol1 "     << ws->var("coeffPol1")->getVal()     << " " << ws->var("coeffPol1")->getError() << "\n"
      << "coeffPol2 "     << ws->var("coeffPol2")->getVal()     << " " << ws->var("coeffPol2")->getError() << "\n"
      << "coeffPol3 "     << ws->var("coeffPol3")->getVal()     << " " << ws->var("coeffPol3")->getError() << "\n"
      << "coeffPol4 "     << ws->var("coeffPol4")->getVal()     << " " << ws->var("coeffPol4")->getError() << "\n";
  }
  if (!mBkgFunct.compare("pol5")) {
    outputFile
      << "coeffPol1 "     << ws->var("coeffPol1")->getVal()     << " " << ws->var("coeffPol1")->getError() << "\n"
      << "coeffPol2 "     << ws->var("coeffPol2")->getVal()     << " " << ws->var("coeffPol2")->getError() << "\n"
      << "coeffPol3 "     << ws->var("coeffPol3")->getVal()     << " " << ws->var("coeffPol3")->getError() << "\n"
      << "coeffPol4 "     << ws->var("coeffPol4")->getVal()     << " " << ws->var("coeffPol4")->getError() << "\n"
      << "coeffPol5 "     << ws->var("coeffPol5")->getVal()     << " " << ws->var("coeffPol5")->getError() << "\n";
  }
  outputFile
    << "coeffGaus "    << ws->var("coeffGaus")->getVal()    << " " << ws->var("coeffGaus")->getError() << "\n"
    << "meanSig1 "     << ws->var("meanSig1")->getVal()     << " " << ws->var("meanSig1")->getError() << "\n"
    //    << "meanSig1P "    << ws->var("meanSig1P")->getVal()     << " " << ws->var("meanSig1P")->getError() << "\n"
    << "sigmaSig1 "    << ws->var("sigmaSig1")->getVal()    << " " << ws->var("sigmaSig1")->getError() << "\n"
    << "wideFactor "   << ws->var("wideFactor")->getVal()   << " " << ws->var("wideFactor")->getError() << "\n"
    << "sigmaSig2 "    << ws->function("sigmaSig2")->getVal() << " " << ws->function("sigmaSig2")->getVal()*ws->var("wideFactor")->getError()/ws->var("wideFactor")->getVal() << "\n"
    << "alpha "        << ws->var("alpha")->getVal()        << " " << ws->var("alpha")->getError() << "\n"
    << "enne "         << ws->var("enne")->getVal()         << " " << ws->var("enne")->getError() << "\n"
    //    << "enneW "        << ws->var("enneW")->getVal()        << " " << ws->var("enneW")->getError() << "\n"
    << "chi2 "         << UnNormChi2                        << "\n"
    << "DOF "          << Dof                               << "\n"
    << "NLL "          << theNLL                            << endl;


  fInData.Close();
  resultF.Close();
  outputFile.close();

  return 0;
}




/////////////////////////////////////////////////////////
//////////////////// Sub-routines ///////////////////////
/////////////////////////////////////////////////////////
void getOptRange(string &ran, float *min, float *max) {
  if (sscanf(ran.c_str(), "%f-%f", min, max) == 0) {
    cout << ran.c_str() << ": not valid!" << endl;
    assert(0);
  }
  return ;
}

void setWSRange(RooWorkspace *ws){
  ws->cat("Jpsi_Type")->setRange("glbglb","GG");
  return;
}


void defineMassBkg(RooWorkspace *ws) {
  // 1st order polynomial
  ws->factory("Chebychev::pol1(Jpsi_Mass,{coeffPol1[-0.8,-1.0,1.0]})");
  // 2nd order polynomial
  ws->factory("Chebychev::pol2(Jpsi_Mass,{coeffPol1, coeffPol2[0.0,-1.0,1.0]})");
  // 3rd order polynomial
  ws->factory("Chebychev::pol3(Jpsi_Mass,{coeffPol1, coeffPol2, coeffPol3[0.0,-1.0,1.0]})");
  // 4th order polynomial
  ws->factory("Chebychev::pol4(Jpsi_Mass,{coeffPol1, coeffPol2, coeffPol3, coeffPol4[0.0,-1.0,1.0]})");
  // 5th order polynomial
  ws->factory("Chebychev::pol5(Jpsi_Mass,{coeffPol1, coeffPol2, coeffPol3, coeffPol4, coeffPol5[0.0,-1.0,1.0]})");
  // 6th order polynomial
  ws->factory("Chebychev::pol6(Jpsi_Mass,{coeffPol1, coeffPol2, coeffPol3, coeffPol4, coeffPol5, coeffPol6[0.0,-1.0,1.0]})");

  // expo
  ws->factory("Exponential::expFunct(Jpsi_Mass,coeffExp[-0.1,-3.0,1.0])");

  // gauss
  ws->factory("Gaussian::bkgGaus(Jpsi_Mass,meanBkg[0.0,0.0,10.0],sigmaBkg[1.0,0.5,5.0])");

  return;
}

/* old
void defineMassBkg(RooWorkspace *ws) {
  // 1st order polynomial
  ws->factory("Chebychev::pol1(Jpsi_Mass,{coeffPol1[0.8,-1.,1.]})");
  //  ws->factory("Chebychev::pol1(Jpsi_Mass,{coeffPol1[-0.8]})");
  ws->var("coeffPol1")->setVal(-0.8);
  ws->var("coeffPol1")->setConstant(false);

  // 2nd order polynomial
  ws->factory("Chebychev::pol2(Jpsi_Mass,{coeffPol1, coeffPol2[0.00,-1.,1.]})");
  ws->factory("Chebychev::pol3(Jpsi_Mass,{coeffPol1, coeffPol2, coeffPol3[0.00,-1.,1.]})");
  //  ws->factory("Chebychev::pol2(Jpsi_Mass,{coeffPol1, coeffPol2[0.00]})");
  //  ws->var("coeffPol2")->setConstant(false);
  //  ws->factory("Chebychev::pol5(Jpsi_Mass,{coeffPol1, coeffPol2, coeffPol3[0.00,-15.,15.], coeffPol4[0.00,-15.,15.], coeffPol5[0.00,-15.,15.]})");
  //  RooFormulaVar coeffPol2("coeffPol2","@0/(-2.0*3.4)",RooArgList(*(ws->var("coeffPol1")))); ws->import(coeffPol2);
  //  ws->factory("Polynomial::pol2(Jpsi_Mass,{coeffPol1, coeffPol2})");
  // Exponential
  //  RooRealVar coeffExp("coeffExp","coeffExp",-0.1,-3.0,1.0);coeffExp.setConstant(false);ws->import("coeffExp");
  ws->factory("Exponential::expFunct(Jpsi_Mass,coeffExp[-0.1,-3.0,1.0])");
  ws->factory("Exponential::expFunctP(Jpsi_Mass,coeffExpP[-1.,-3.,0.])");
  // 2nd Exponential
  //  RooRealVar coeffExp2("coeffExp2","coeffExp2",-0.2,-3.0,1.0); coeffExp2.setConstant(false);ws->import("coeffExp");
  ws->factory("Exponential::expFunct2(Jpsi_Mass,coeffExp2[-0.2,-3.0,1.0])");
  // Sum of two exponentials
  //  ws->factory("SUM::twoExpFunct(coeffBkg[0.9,0.1,1.0]*expFunct,expFunct2)");


  // RooRealVar turnOn("turnOn","turnOn", 6.);
  // turnOn.setConstant(false);
  // RooRealVar width("width","width", 1., 0., 20.);
  // RooRealVar decay("decay","decay", 7.);
  // decay.setConstant(false);
  // RooGenericPdf bkgErfExp("bkgErfExp","bkg","exp(-@0/@3)*(TMath::Erf((@0-@1)/@2)+1)", RooArgList(*(ws->var("Jpsi_Mass")), turnOn, width, decay));
  // ws->import(bkgErfExp);
  // RooGenericPdf bkgErf("bkgErf","bkgErf","(TMath::Erf((@0-@1)/@2)+1)", RooArgList(*(ws->var("Jpsi_Mass")), turnOn, width));
  // ws->import(bkgErf);

  RooRealVar cutx("cutx","cutx",3.4,3.2,3.5); cutx.setConstant(false); ws->import(cutx); //Region below than this cut will use cTau1, above than it will use cTau2. determined by free fit.

  RooGenericPdf twoExpFunct("twoExpFunct","twoExpFunct","1.0/abs(@1)*(exp(@0*@1)*(@0<@3)+exp(@0*@2 - @3*(@2-@1) )*(@0>=@3))", RooArgList(*(ws->var("Jpsi_Mass")), *(ws->var("coeffExp")), *(ws->var("coeffExp2")), *(ws->var("cutx"))));
  ws->import(twoExpFunct);
  
  
  RooRealVar coeffExpPol1("coeffExpPol1","coeffExpPol1",-0.1);coeffExpPol1.setConstant(false);
  RooRealVar coeffExpPol2("coeffExpPol2","coeffExpPol2",-0.01);coeffExpPol2.setConstant(false); 

  RooGenericPdf expPol2Funct("expPol2Funct","expPol2Funct","exp(@0*@1+pow(@0*@2,2))*(3.5>@0)", RooArgList(*(ws->var("Jpsi_Mass")), coeffExpPol1, coeffExpPol2));
  ws->import(expPol2Funct);


  // RooRealVar cutx("cutx","cutx",4.2,3.3,4.2); cutx.setConstant(true); ws->import(cutx); //Region below than this cut will use cTau1, above than it will use cTau2. determined by free fit.
  // RooRealVar coeffExp("coeffExp","coeffExp",-0.1,-3.0,1.0); coeffExp.setConstant(true);ws->import(coeffExp);
  // RooRealVar coeffExp2("coeffExp2","coeffExp2",-0.2,-3.0,1.0); coeffExp2.setConstant(true);ws->import(coeffExp2);
  // RooExp2 twoExpFunct("twoExpFunct","twoExpFunct",*(ws->var("Jpsi_Mass")),cutx,coeffExp,coeffExp2); ws->import(twoExpFunct);


  ws->factory("Gaussian::bkgGaus(Jpsi_Mass,meanBkg[0.0,0.0,10.0],sigmaBkg[1.0,0.5,5.0])");
  

  // RooRealVar *m_step = new RooRealVar("m_step","m_step",3.4,3.3,3.5); 
  // //  m_step->setConstant(true);
  // ws->import(*m_step);
  // RooStats::Heaviside *step = new RooStats::Heaviside("step","split two expoenentials",*(ws->var("Jpsi_Mass")),*(ws->var("m_step")));  
  // RooFormulaVar *coeffBkg(
  // ws->import(*coeffBkg);
  
  return;
}
*/

void defineMassSig(RooWorkspace *ws) {
  // narrow Gauss
  ws->factory("Gaussian::signalG1(Jpsi_Mass,meanSig1[3.0969,3.05,3.15],sigmaSig1[0.03,0.005,0.080])");
  // narrow CB
  ws->factory("CBShape::signalCB1(Jpsi_Mass,meanSig1,sigmaSig1,alpha[1.0,0.0,3.0],enne[5.0,1.0,50.0])");

  RooRealVar wideFactor("wideFactor","wideFactor",1.2,1.0,5.0);ws->import(wideFactor);
  RooFormulaVar sigmaSig2("sigmaSig2","@0*@1",RooArgList(*(ws->var("sigmaSig1")),wideFactor));ws->import(sigmaSig2);
  // wide Gauss
  ws->factory("Gaussian::signalG2(Jpsi_Mass,meanSig1,sigmaSig2)");

  // CB(narrow) + Gauss(wide)
  ws->factory("SUM::sigCB1G2(coeffGaus[0.1,0.0,1.0]*signalG2,signalCB1)");

  // Fix Jpsi-psi' mass difference
  //  RooFormulaVar meanSig1P("meanSig1P","@0+0.58919",RooArgList(*(ws->var("meanSig1")))); ws->import(meanSig1P);
  // Fix Jpsi/psi' mass ratio
  RooFormulaVar meanSig1P("meanSig1P","@0*1.19025",RooArgList(*(ws->var("meanSig1")))); ws->import(meanSig1P);
  // Fix resolution scale: sigma_MJpsi/MJpsi = sigma_Mpsi'/Mpsi'
  RooFormulaVar sigmaSig1P("sigmaSig1P","@0*1.19025",RooArgList(*(ws->var("sigmaSig1")))); ws->import(sigmaSig1P);
  RooFormulaVar sigmaSig2P("sigmaSig2P","@0*1.19025",RooArgList(*(ws->function("sigmaSig2")))); ws->import(sigmaSig2P);
  // narrow Gauss
  ws->factory("Gaussian::signalG1P(Jpsi_Mass,meanSig1P,sigmaSig1P)");
  // wide Gauss
  ws->factory("Gaussian::signalG2P(Jpsi_Mass,meanSig1P,sigmaSig2P)");
  // narrow CB
  ws->factory("CBShape::signalCB1P(Jpsi_Mass,meanSig1P,sigmaSig1P,alpha,enne)");

  // CB(narrow) + Gauss(wide)
  ws->factory("SUM::sigCB1G2P(coeffGaus*signalG2P,signalCB1P)");

  return;
}

/* old
void defineMassSig(RooWorkspace *ws) {
  //////// Candidates for signal
  // Crystall Ball
  ws->factory("CBShape::signalCB2(Jpsi_Mass,meanSig1[3.0969,3.05,3.15],sigmaSig2[0.03,0.01,0.060],alpha[1.,0.,3.],enne[5.,1.,30.])");
  ws->factory("CBShape::signalCB2WN(Jpsi_Mass,meanSig1,sigmaSig2,alpha,enneW[5.,1.,50.])");

  // Normal gaussians
  //  RooRealVar wideFactor("wideFactor","wideFactor",1.2,1.0,5.0);ws->import(wideFactor);
  //  RooFormulaVar sigmaSig1("sigmaSig1","@0*@1",RooArgList(*(ws->var("sigmaSig2")),wideFactor));ws->import(sigmaSig1);

  //  ws->factory("Gaussian::signalG1(Jpsi_Mass,meanSig1,sigmaSig1)");
  ws->factory("Gaussian::signalG1(Jpsi_Mass,meanSig1,sigmaSig1[0.08,0.06,0.2])");
//  ws->factory("Gaussian::signalG1(Jpsi_Mass,meanSig1[3.0975,3.05,3.15],sigmaSig1[0.1,0.06,0.2])");
//  ws->factory("Gaussian::signalG2(Jpsi_Mass,meanSig2[3.0975,3.05,3.15],sigmaSig2[0.03,0.008,0.2])");
//  ws->factory("Gaussian::signalGm1s2(Jpsi_Mass,meanSig1,sigmaSig2)");

  ws->factory("CBShape::signalCB(Jpsi_Mass,meanSig1,sigmaSig1,alpha,enne)");
  ws->factory("CBShape::signalCBWN(Jpsi_Mass,meanSig1,sigmaSig1,alpha,enneW)");


  // Fix Jpsi-psi' mass difference
  //  RooFormulaVar meanSig1P("meanSig1P","@0+0.58919",RooArgList(*(ws->var("meanSig1")))); ws->import(meanSig1P);
  // Fix Jpsi-psi' mass ratio
  RooFormulaVar meanSig1P("meanSig1P","@0*1.19025",RooArgList(*(ws->var("meanSig1")))); ws->import(meanSig1P);
  // Fix resolution scale: sigma_MJpsi/MJpsi = sigma_Mpsi'/Mpsi'
  RooFormulaVar sigmaSig1P("sigmaSig1P","@0*1.19025",RooArgList(*(ws->function("sigmaSig1")))); ws->import(sigmaSig1P);
  RooFormulaVar sigmaSig2P("sigmaSig2P","@0*1.19025",RooArgList(*(ws->var("sigmaSig2")))); ws->import(sigmaSig2P);
  ws->factory("Gaussian::signalG1P(Jpsi_Mass,meanSig1P,sigmaSig1P)");
  //  ws->factory("Gaussian::signalG1P(Jpsi_Mass,meanSig1P[3.6861,3.4,4.0],sigmaSig1P)");
  ws->factory("CBShape::signalCB2WNP(Jpsi_Mass,meanSig1P,sigmaSig2P,alpha,enneW)");

  //////// Sum of signal functions
  // Sum of gaussian 1 and a crystall ball
  ws->factory("SUM::sigCBG1(coeffGaus[0.1,0.05,0.93]*signalG1,signalCB)");
  // Sum of gaussian 1 and crystall ball 2
  ws->factory("SUM::sigCB2G1(coeffGaus*signalG1,signalCB2)");
  // Sum of gaussian 1 and crystall ball with wide n
  ws->factory("SUM::sigCBWNG1(coeffGaus*signalG1,signalCBWN)");
  // Sum of gaussian 1 and crystall ball 2 with wide n
  ws->factory("SUM::sigCB2WNG1(coeffGaus*signalG1,signalCB2WN)");
  ws->factory("SUM::sigCB2WNG1P(coeffGaus*signalG1P,signalCB2WNP)");

  //  ws->factory("SUM::sigCBGJpsiPsiP(fracP[0.1,0.0,1.0]*sigCB2WNG1P,sigCB2WNG1)");
  //  ws->factory("SUM::sigCBJpsiPsiP(fracP*signalCB2WNP,signalCB2WN)");

  return;
}

*/
