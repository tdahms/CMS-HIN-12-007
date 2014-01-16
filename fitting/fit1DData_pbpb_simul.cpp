#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

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
#include <RooHistPdfConv.h>
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

using namespace RooFit;
using namespace std;

bool superImpose = true;
bool analyticBlifetime = true;
bool narrowSideband = false;
bool oneGaussianResol = false;

void getOptRange(string &ran,float *min,float *max);
void setWSRange(RooWorkspace *ws);
RooBinning setCtBinning(float lmin,float lmax);
void defineMassBkg(RooWorkspace *ws);
void defineMassSig(RooWorkspace *ws);
double computeRatio(RooRealVar& x, RooRealVar& y);
double computeRatioError(RooRealVar& x, RooRealVar& y, double correlation = 0.);

int main (int argc, char* argv[]) {
  gROOT->Macro("~/logon.C");

  string FileNameMC1, FileNameMC2;
  vector<string> FileName;
  vector<string> mBkgFunct;
  string mJpsiFunct, mPsiPFunct;
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
  bool logScale=false;

  // *** Check options
  for (int i=1; i<argc; i++) {
    char *tmpargv = argv[i];
    switch (tmpargv[0]) {
    case '-':{
      switch (tmpargv[1]) {
      case 'f':
	FileName.push_back(argv[i+1]);
	FileName.push_back(argv[i+2]);
	cout << "Fitted PbPb data file: " << FileName[0] << endl;
	cout << "Fitted pp data file: " << FileName[1] << endl;
	break;
      case 'v':
	mJpsiFunct = argv[i+1];
	mPsiPFunct = argv[i+2];
	mBkgFunct.push_back(argv[i+3]);
	mBkgFunct.push_back(argv[i+4]);
	cout << "Mass J/psi function: " << mJpsiFunct << endl;
	cout << "Mass psi(2S) function: " << mPsiPFunct << endl;
	cout << "Mass PbPb background function: " << mBkgFunct[0] << endl;
	cout << "Mass pp background function: " << mBkgFunct[1] << endl;
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
      case 'l':
	logScale = atoi(argv[i+1]);
	cout << "plot on log scale: " << logScale << endl;
	break;
        }
      }
    }
  }// End check options
 
  float pmin=0, pmax=0, ymin=0, ymax=0, cmin=0, cmax=0, psmax=0, psmin=0;
  getOptRange(prange,&pmin,&pmax);
  getOptRange(crange,&cmin,&cmax);
  getOptRange(yrange,&ymin,&ymax);
  getOptRange(phirange,&psmin,&psmax);


  // *** TFile for saving fitting results
  string resultFN;
  resultFN = dirPre + "_rap" + yrange + "_pT" + prange + "_cent" + crange + "_dPhi" + phirange + "_fitResult.root";
  TFile resultF(resultFN.c_str(),"RECREATE");

  // *** Read Data files
  const int nFiles = FileName.size();
  TFile *fInData[nFiles];
  RooDataSet *data[nFiles];
  for (int i=0; i<nFiles; i++) {
    char name[100] = {0};
    fInData[i] = new TFile(FileName[i].c_str());
    cout << FileName[i].c_str() << endl;
    if (fInData[i]->IsZombie()) { cout << "CANNOT open data root file: " << FileName[i] << "\n"; return 1; }
    fInData[i]->cd();
    data[i] = (RooDataSet*)fInData[i]->Get("dataJpsi");
    sprintf(name,"data%d",i+1);
    data[i]->SetName(name);
  }

  // Create workspace to play with
  RooWorkspace *ws = new RooWorkspace("workspace");

  // Reduce "dataMC" with given ranges/cuts
  char reduceDS[300];
  sprintf(reduceDS,"Jpsi_Pt>%.2f && Jpsi_Pt<%.2f && abs(Jpsi_Y)>%.2f && abs(Jpsi_Y)<%.2f",pmin,pmax,ymin,ymax);
  cout << "reduceDS for MC and data: " << reduceDS << endl;

  RooDataSet *redData[nFiles];
  for (int i=0; i<FileName.size(); i++) {
    redData[i] = (RooDataSet*)data[i]->reduce(reduceDS);
    ws->import(*redData[i]);
  }
 
  setWSRange(ws);

  // Draw data
  string titlestr;

  // Define binning for true lifetime
  ws->var("Jpsi_Mass")->setBins(40);

  // Additional cuts on data and get sub-datasets/histograms
  RooDataSet *redDataCut[nFiles];
  string reduceDSstr;
  if (isGG == 0) {
    for (int i=0; i<nFiles; i++) {
      redDataCut[i] = (RooDataSet*)redData[i]->reduce("Jpsi_Type == Jpsi_Type::GG");
    }
  } else if (isGG == 1) {
    for (int i=0; i<nFiles; i++) {
      redDataCut[i] = (RooDataSet*)redData[i]->reduce("Jpsi_Type == Jpsi_Type::GT || Jpsi_Type == Jpsi_Type::GG");
    }
  } else {
    for (int i=0; i<nFiles; i++) {
      redDataCut[i] = (RooDataSet*)redData[i]->reduce("Jpsi_Ct < 600000.");
    }
  }

  RooCategory* sample = new RooCategory("sample","sample") ;
  sample->defineType("HI") ;
  sample->defineType("pp") ;
  ws->factory("sample[HI,pp]");  //Index for simultaneous fit

  RooDataSet *redDataCutSim = new RooDataSet("redDataCutSim","redDataCutSim",RooArgSet(*(ws->var("Jpsi_Mass"))),Index(*(ws->cat("sample"))),Import("HI",*redDataCut[0]),Import("pp",*redDataCut[1]));

  RooDataHist *binData[nFiles];
  for (int i=0; i<nFiles; i++) {
    char name[100] = {0};
    sprintf(name,"binData%d",i+1);
    binData[i]= new RooDataHist(name,name,RooArgSet( *(ws->var("Jpsi_Mass")) ), *redDataCut[i]);
    cout << "DATA" << i+1 << " :: N events to fit: " << binData[i]->sumEntries() << endl;
  }

  // SYSTEMATICS 1 (very sidebands)
  RooDataSet *redDataSB[nFiles];
  RooDataHist *binDataSB[nFiles];
  RooDataSet *redDataSIG[nFiles];
  for (int i=0; i<nFiles; i++) {
    if (narrowSideband) {
     redDataSB[i] = (RooDataSet*) redDataCut[i]->reduce("Jpsi_Mass<2.8 || Jpsi_Mass>3.4");
    } else {
     redDataSB[i] = (RooDataSet*) redDataCut[i]->reduce("Jpsi_Mass<2.9 || Jpsi_Mass>3.3");
    }
    char name[100] = {0};
    sprintf(name,"binDataSB%d",i+1);
    binDataSB[i] = new RooDataHist(name,"Data distribution for background",RooArgSet( *(ws->var("Jpsi_Mass")) ),*redDataSB[i]);
    redDataSIG[i] = (RooDataSet*)redDataCut[i]->reduce("Jpsi_Mass > 2.9 && Jpsi_Mass < 3.3");
  }

  // *** Define PDFs with parameters (mass and ctau)
  // J/psi mass parameterization
  defineMassBkg(ws);
  defineMassSig(ws);

  char funct[100];
  string partTit, partFile;
  if (isGG == 0) { partTit = "glb-glb"; partFile = "GG"; }
  else if (isGG = 1) { partTit = "glb-trk"; partFile = "GT"; }
  else { partTit = "all"; partFile = "ALL"; }

  // Binning for invariant mass distribution
  RooBinning rbm(2.6,4.2);
  rbm.addUniform(40,2.6,4.2);

  // Global TLatex, TH1, TGraph objects for drawing
  TLatex *t = new TLatex();
  t->SetNDC(); t->SetTextAlign(12);

  Double_t fx[2], fy[2], fex[2], fey[2];
  TGraphErrors *gfake1 = new TGraphErrors(2,fx,fy,fex,fey);
  gfake1->SetMarkerStyle(20); gfake1->SetMarkerSize(1);
  TH1F hfake11 = TH1F("hfake1","hfake1",100,200,300);
  hfake11.SetLineColor(kBlue); hfake11.SetLineWidth(4); hfake11.SetLineStyle(7); hfake11.SetFillColor(kAzure-9); hfake11.SetFillStyle(1001);
  TH1F hfake21 = TH1F("hfake2","hfake2",100,200,300);
  hfake21.SetLineColor(kBlack); hfake21.SetLineWidth(4); hfake21.SetFillColor(kBlack); hfake21.SetFillStyle(3354);
  TH1F hfake22 = TH1F("hfake22","hfake22",100,200,300);
  hfake22.SetLineColor(kRed); hfake22.SetLineWidth(5); hfake22.SetLineStyle(2);
  TH1F hfake31 = TH1F("hfake3","hfake3",100,200,300);
  hfake31.SetLineColor(kRed); hfake31.SetMarkerStyle(kCircle); hfake31.SetLineWidth(4); hfake31.SetMarkerColor(kRed); hfake31.SetLineStyle(9); hfake31.SetFillColor(kRed-7); hfake31.SetFillStyle(3444);


  RooFitResult *fitM;
  RooFitResult *fitM_pre;
  struct PARAM {
    double fracP; double fracPErr;
    double fracP_pp; double fracPErr_pp;
    double coeffExp;  double coeffExpErr;
    double coeffExp_pp;  double coeffExpErr_pp;
    double coeffExp2;  double coeffExp2Err;
    double coeffExp2_pp;  double coeffExp2Err_pp;
    double cutx;  double cutxErr;
    double cutx_pp;  double cutxErr_pp;
    double coeffGaus; double coeffGausErr;
    double coeffGaus_pp; double coeffGausErr_pp;
    double meanSig1;  double meanSig1Err;
    double sigmaSig1; double sigmaSig1Err;
    double sigmaSig2; double sigmaSig2Err;
    double alpha;     double alphaErr;
    double enne;      double enneErr;
    double enneW;     double enneWErr;
    double meanSig1_pp;  double meanSig1Err_pp;
    double sigmaSig1_pp; double sigmaSig1Err_pp;
    double sigmaSig2_pp; double sigmaSig2Err_pp;
    double alpha_pp;     double alphaErr_pp;
    double enne_pp;      double enneErr_pp;
    double enneW_pp;     double enneWErr_pp;
  };

  PARAM inputCB;
  PARAM inputCBG;


  bool centConst = false;  //False: fit w/o any constrained parameters (centrality dep.)
  bool dPhiConst = false;  //False: fit w/o any constrained parameters (dPhi dep.)
  double inputN[6] = {0};  //Number of Jpsi, psiP, fracP, double ratio, and background in the 0-1.571 rad bin
  if (isMB != 0) {
    if (isMB != 4 && (mJpsiFunct.compare("sigCB2WNG1") || mBkgFunct[0].compare("expFunct"))) {
      cout << "For fit systematics:\n";
      cout << "1) Signal shape check: should use default runOpt (4)\n";
      cout << "2) Background shape check: should use defatul runOpt (4)\n";
      cout << "3) Constrain fit: should use sigCB2WNG1 for signal shape, expFunct for bkg shape with runOpt (3)\n";
      return -1;
    }
  } //End of isMB != 0
    
  RooRealVar *NBkg_pp  = new RooRealVar("NBkg_pp","pp background yield", 2000.0,1.0,500000.0);   ws->import(*NBkg_pp);
  RooRealVar *NBkg_HI  = new RooRealVar("NBkg_HI","PbPb background yield", 200.0,1.0,500000.0);   ws->import(*NBkg_HI);

  //  RooRealVar *DoubleRatio  = new RooRealVar("DoubleRatio","psi(2S) RAA / J/psi RAA" ,1.0,0.0,15.0);   ws->import(*DoubleRatio);

  RooRealVar *NJpsi_pp  = new RooRealVar("NJpsi_pp","pp J/psi yield", 1000.0, 1.0, 50000.0);  ws->import(*NJpsi_pp);
  RooRealVar *NJpsi_HI  = new RooRealVar("NJpsi_HI","PbPb J/psi yield", 4000.0, 1.0, 50000.0);  ws->import(*NJpsi_HI);
  RooRealVar *fracP_pp  = new RooRealVar("fracP_pp","pp psi(2S) fraction" ,0.1);fracP_pp->setConstant(false); ws->import(*fracP_pp);
  RooRealVar *fracP_HI  = new RooRealVar("fracP_HI","PbPb psi(2S) fraction" ,0.1);fracP_HI->setConstant(false); ws->import(*fracP_HI);

  //  RooRealVar *fracP_HI  = new RooRealVar("fracP_HI","PbPb psi(2S) fraction" ,1.0,0.0,1.0);   ws->import(*fracP_pp);

  //  RooRealVar *fracP_pp  = new RooRealVar("fracP_pp","pp psi(2S) fraction" ,0.1,0.0,1.0);   ws->import(*fracP_pp);
  //  RooFormulaVar *NPsiP_pp = new RooFormulaVar("NPsiP_pp", "@0*@1", RooArgList(*(ws->var("NJpsi_pp")),*(ws->var("fracP_pp"))));  ws->import(*NPsiP_pp);

  //  RooFormulaVar *fracP_HI  = new RooFormulaVar("fracP_HI","@0*@1",RooArgList(*(ws->var("DoubleRatio")),*(ws->var("fracP_pp"))));   ws->import(*fracP_HI);

  //  RooRealVar *fracP_HI  = new RooRealVar("fracP_HI","PbPb psi(2S) fraction" ,0.1,0.0,1.0);   ws->import(*fracP_HI);
  //RooFormulaVar *fracP_HI  = new RooFormulaVar("fracP_HI","@0*@1",RooArgList(*(ws->var("DoubleRatio")),*(ws->var("fracP_pp"))));   ws->import(*fracP_HI);

  RooFormulaVar *NPsiP_HI  = new RooFormulaVar("NPsiP_HI","@0*@1",RooArgList(*(ws->var("NJpsi_HI")),*(ws->var("fracP_HI"))));   ws->import(*NPsiP_HI);
  //  RooFormulaVar *NPsiP_HI  = new RooFormulaVar("NPsiP_HI","@0*@1*@2",RooArgList(*(ws->var("NJpsi_HI")),*(ws->var("DoubleRatio")),*(ws->var("fracP_pp"))));   ws->import(*NPsiP_HI);

  /* fracP_HI = free paramter
     RooFormulaVar *fracP_pp  = new RooFormulaVar("fracP_pp","@1/@0",RooArgList(*(ws->var("DoubleRatio")),*(ws->var("fracP_HI"))));   ws->import(*fracP_pp);
     RooFormulaVar *NPsiP_pp = new RooFormulaVar("NPsiP_pp", "@0*@2/@1", RooArgList(*(ws->var("NJpsi_pp")),*(ws->var("DoubleRatio")),*(ws->var("fracP_HI"))));  ws->import(*NPsiP_pp);
  */
  RooFormulaVar *NPsiP_pp = new RooFormulaVar("NPsiP_pp", "@0*@1", RooArgList(*(ws->var("NJpsi_pp")),*(ws->var("fracP_pp"))));  ws->import(*NPsiP_pp);

  RooFormulaVar *DoubleRatio  = new RooFormulaVar("DoubleRatio","@0/@1",RooArgList(*(ws->var("fracP_HI")),*(ws->var("fracP_pp"))));   ws->import(*DoubleRatio);
  RooFormulaVar *NPsiP_mix  = new RooFormulaVar("NPsiP_mix","@0*@1",RooArgList(*(ws->var("NJpsi_pp")),*(ws->function("fracP_HI"))));   ws->import(*NPsiP_mix);

  // double ratio = RAA(psiP)/RAA(Jpsi)
  //  RooFormulaVar *RatioJpsi  = new RooFormulaVar("RatioJpsi","@0/@1",RooArgList(*(ws->var("RatioPsiP")),*(ws->var("DoubleRatio"))));   ws->import(*RatioJpsi);

  //  RooFormulaVar *NJpsi_HI  = new RooFormulaVar("NJpsi_HI","@0*@1/@2",RooArgList(*(ws->var("RatioPsiP")),*(ws->var("NJpsi_pp")),*(ws->var("DoubleRatio"))));   ws->import(*NJpsi_HI);




  // RooRealVar *NJpsi_HI  = new RooRealVar("NJpsi_HI","PbPb J/psi yield", 4000.0, 1.0, 50000.0);  ws->import(*NJpsi);
  // RooRealVar *NBkg_HI  = new RooRealVar("NBkg_HI","PbPb Brackground yield", 2000.0,1.0,500000.0);   ws->import(*NBkg);
  // RooRealVar *RatioJpsi  = new RooRealVar("RatioJpsi","J/psi RAA" ,1.0,0.0,10.0);   ws->import(*RatioJpsi);


  if (prefitMass) {
    sprintf(funct,"SUM::jpsiMassPDF_HI(NJpsi_HI*%s,NBkg_HI*%s)",mJpsiFunct.c_str(),mBkgFunct[0].c_str());
    ws->factory(funct);
    sprintf(funct,"SUM::jpsiMassPDF_pp(NJpsi_pp*%s,NBkg_pp*%s)",(mJpsiFunct+"_pp").c_str(),mBkgFunct[1].c_str());
    ws->factory(funct);
    ws->factory("SIMUL::jpsiMassPDFSim(sample,HI=jpsiMassPDF_HI,pp=jpsiMassPDF_pp)");
  }


  sprintf(funct,"SUM::sigMassPDF_pp(NJpsi_pp*%s,NPsiP_pp*%s,NBkg_pp*%s)",(mJpsiFunct+"_pp").c_str(),(mPsiPFunct+"_pp").c_str(),mBkgFunct[1].c_str());
  ws->factory(funct);
  sprintf(funct,"SUM::sigMassPDF_HI(NJpsi_HI*%s,NPsiP_HI*%s,NBkg_HI*%s)",mJpsiFunct.c_str(),mPsiPFunct.c_str(),mBkgFunct[0].c_str());
  ws->factory(funct);
  sprintf(funct,"SUM::sigMassPDF_mix(NJpsi_pp*%s,NPsiP_mix*%s,NBkg_pp*%s)",(mJpsiFunct+"_pp").c_str(),(mPsiPFunct+"_pp").c_str(),mBkgFunct[1].c_str());
  ws->factory(funct);
  
  ws->factory("SIMUL::sigMassPDFSim(sample,HI=sigMassPDF_HI,pp=sigMassPDF_pp)");



  string str("CBG");
  string dirPre2 = dirPre;
  size_t found = dirPre2.find(str);
  if (found!=string::npos)
    dirPre2.replace(found,str.length(),"CB");
  
  if ( (crange.compare("0-100") || found!=string::npos) ) {
    // fix CB parameters within centrality bins to MB (use CB also for CBG)
    string inputFNcb;

    inputFNcb =  dirPre2 + "_rap" + yrange + "_pT" + prange + "_cent0-100_dPhi.txt";
    ifstream input;
    input.open(inputFNcb.c_str());
    if (!input.good()) { cout << "Failed to open: " <<inputFNcb << endl; return 1; }
    string tmp;
    double inputTmp[2] = {0};
    input >> tmp >> inputTmp[0] >> inputTmp[1]; //NJpsi1
    inputN[0] = inputTmp[0];  //NJpsi1
    input >> tmp >> inputTmp[0] >> inputTmp[1]; //NBkg1
    inputN[1] = inputTmp[0];  //NBkg1
    input >> tmp >> inputTmp[0] >> inputTmp[1]; //double ratio
    inputN[2] = inputTmp[0];  //double ratio
    input >> tmp >> inputTmp[0] >> inputTmp[1]; //R_psi2S1
    inputN[3] = inputTmp[0];  //R_psi2S1
    input >> tmp >> inputTmp[0] >> inputTmp[1]; //NJpsi2
    inputN[4] = inputTmp[0];  //NJpsi2
    input >> tmp >> inputTmp[0] >> inputTmp[1]; //NBkg2
    inputN[5] = inputTmp[0];  //NBkg2
    input >> tmp >> inputTmp[0] >> inputTmp[1]; //R_psi2S2
    inputN[6] = inputTmp[0];  //R_psi2S2


    for (int p=0; p<16; p++) {   //Mass signal parameters
      input >> tmp >> inputTmp[0] >> inputTmp[1];
      cout << tmp << " " << inputTmp[0] << endl;
      if (!tmp.compare("coeffGaus")) inputCB.coeffGaus = inputTmp[0];
      else if (!tmp.compare("coeffGauss_pp")) inputCB.coeffGaus_pp = inputTmp[0];
      else if (!tmp.compare("coeffExp")) inputCB.coeffExp = inputTmp[0];
      else if (!tmp.compare("coeffExp_pp")) inputCB.coeffExp_pp = inputTmp[0];
      else if (!tmp.compare("coeffExp2")) inputCB.coeffExp2 = inputTmp[0];
      else if (!tmp.compare("coeffExp2_pp")) inputCB.coeffExp2_pp = inputTmp[0];
      else if (!tmp.compare("meanSig1")) inputCB.meanSig1 = inputTmp[0];
      else if (!tmp.compare("sigmaSig1")) inputCB.sigmaSig1 = inputTmp[0];
      else if (!tmp.compare("sigmaSig2")) inputCB.sigmaSig2 = inputTmp[0];
      else if (!tmp.compare("alpha")) inputCB.alpha = inputTmp[0];
      else if (!tmp.compare("enneW")) inputCB.enneW = inputTmp[0];
      else if (!tmp.compare("meanSig1_pp")) inputCB.meanSig1_pp = inputTmp[0];
      else if (!tmp.compare("sigmaSig1_pp")) inputCB.sigmaSig1_pp = inputTmp[0];
      else if (!tmp.compare("sigmaSig2_pp")) inputCB.sigmaSig2_pp = inputTmp[0];
      else if (!tmp.compare("alpha_pp")) inputCB.alpha_pp = inputTmp[0];
      else if (!tmp.compare("enneW_pp")) inputCB.enneW_pp = inputTmp[0];
    }


    ws->var("alpha")->setVal(inputCB.alpha);
    ws->var("enneW")->setVal(inputCB.enneW);

    ws->var("alpha_pp")->setVal(inputCB.alpha_pp);
    ws->var("enneW_pp")->setVal(inputCB.enneW_pp);

    ws->var("alpha")->setConstant(kTRUE);
    ws->var("enneW")->setConstant(kTRUE);

    ws->var("alpha_pp")->setConstant(kTRUE);
    ws->var("enneW_pp")->setConstant(kTRUE);
  }


  if ( found!=string::npos ) {
    // fix CB parameters within centrality bins to MB (use CB also for CBG)
    string inputFNcb;

    inputFNcb =  dirPre2 + "_rap" + yrange + "_pT" + prange + "_cent" + crange + "_dPhi.txt";
    ifstream input;
    input.open(inputFNcb.c_str());
    if (!input.good()) { cout << "Failed to open: " <<inputFNcb << endl; return 1; }
    string tmp;
    double inputTmp[2] = {0};
    PARAM inputCB;
    input >> tmp >> inputTmp[0] >> inputTmp[1]; //NJpsi1
    inputN[0] = inputTmp[0];  //NJpsi1
    input >> tmp >> inputTmp[0] >> inputTmp[1]; //NBkg1
    inputN[1] = inputTmp[0];  //NBkg1
    input >> tmp >> inputTmp[0] >> inputTmp[1]; //double ratio
    inputN[2] = inputTmp[0];  //double ratio
    input >> tmp >> inputTmp[0] >> inputTmp[1]; //R_psi2S1
    inputN[3] = inputTmp[0];  //R_psi2S1
    input >> tmp >> inputTmp[0] >> inputTmp[1]; //NJpsi2
    inputN[4] = inputTmp[0];  //NJpsi2
    input >> tmp >> inputTmp[0] >> inputTmp[1]; //NBkg2
    inputN[5] = inputTmp[0];  //NBkg2
    input >> tmp >> inputTmp[0] >> inputTmp[1]; //R_psi2S2
    inputN[6] = inputTmp[0];  //R_psi2S2


    for (int p=0; p<16; p++) {   //Mass signal parameters
      input >> tmp >> inputTmp[0] >> inputTmp[1];
      cout << tmp << " " << inputTmp[0] << endl;
      if (!tmp.compare("cutx")) inputCB.cutx = inputTmp[0];
      else if (!tmp.compare("cutx_pp")) inputCB.cutx_pp = inputTmp[0];
    }


    ws->var("cutx")->setVal(inputCB.cutx);
    ws->var("cutx_pp")->setVal(inputCB.cutx_pp);

    ws->var("cutx")->setConstant(kTRUE);
    ws->var("cutx_pp")->setConstant(kTRUE);
  }


  if ( (crange.compare("0-100") && found!=string::npos) ) {
    // fix CB parameters within centrality bins to MB (use CB also for CBG)
    string inputFNcb;

    inputFNcb =  dirPre + "_rap" + yrange + "_pT" + prange + "_cent0-100_dPhi.txt";
    ifstream input;
    input.open(inputFNcb.c_str());
    if (!input.good()) { cout << "Failed to open: " <<inputFNcb << endl; return 1; }
    string tmp;
    double inputTmp[2] = {0};
    PARAM inputCBG;
    input >> tmp >> inputTmp[0] >> inputTmp[1]; //NJpsi1
    inputN[0] = inputTmp[0];  //NJpsi1
    input >> tmp >> inputTmp[0] >> inputTmp[1]; //NBkg1
    inputN[1] = inputTmp[0];  //NBkg1
    input >> tmp >> inputTmp[0] >> inputTmp[1]; //double ratio
    inputN[2] = inputTmp[0];  //double ratio
    input >> tmp >> inputTmp[0] >> inputTmp[1]; //R_psi2S1
    inputN[3] = inputTmp[0];  //R_psi2S1
    input >> tmp >> inputTmp[0] >> inputTmp[1]; //NJpsi2
    inputN[4] = inputTmp[0];  //NJpsi2
    input >> tmp >> inputTmp[0] >> inputTmp[1]; //NBkg2
    inputN[5] = inputTmp[0];  //NBkg2
    input >> tmp >> inputTmp[0] >> inputTmp[1]; //R_psi2S2
    inputN[6] = inputTmp[0];  //R_psi2S2


    for (int p=0; p<16; p++) {   //Mass signal parameters
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
      else if (!tmp.compare("coeffGaus_pp")) {
	inputCBG.coeffGaus_pp = inputTmp[0];
	inputCBG.coeffGausErr_pp = inputTmp[1];
      }
      else if (!tmp.compare("sigmaSig1_pp")) {
	inputCBG.sigmaSig1_pp = inputTmp[0];
	inputCBG.sigmaSig1Err_pp = inputTmp[1];
      }
      else if (!tmp.compare("sigmaSig2_pp")) {
	inputCBG.sigmaSig2_pp = inputTmp[0];
	inputCBG.sigmaSig2Err_pp = inputTmp[1];
      }
      else if (!tmp.compare("alpha")) inputCBG.alpha = inputTmp[0];
      else if (!tmp.compare("enneW")) inputCBG.enneW = inputTmp[0];
      else if (!tmp.compare("coeffExp")) inputCBG.coeffExp = inputTmp[0];
      else if (!tmp.compare("coeffExp2")) inputCBG.coeffExp2 = inputTmp[0];
      else if (!tmp.compare("alpha_pp")) inputCBG.alpha_pp = inputTmp[0];
      else if (!tmp.compare("enneW_pp")) inputCBG.enneW_pp = inputTmp[0];
      else if (!tmp.compare("coeffExp_pp")) inputCBG.coeffExp_pp = inputTmp[0];
      else if (!tmp.compare("coeffExp2_pp")) inputCBG.coeffExp2_pp = inputTmp[0];
    }

    ws->var("sigmaSig1")->setVal(inputCBG.sigmaSig1);
    ws->var("coeffGaus")->setVal(inputCBG.coeffGaus);
    ws->var("sigmaSig1_pp")->setVal(inputCBG.sigmaSig1_pp);
    ws->var("coeffGaus_pp")->setVal(inputCBG.coeffGaus_pp);

    ws->var("sigmaSig1")->setConstant(kTRUE);
    ws->var("coeffGaus")->setConstant(kTRUE);
    ws->var("sigmaSig1_pp")->setConstant(kTRUE);
    ws->var("coeffGaus_pp")->setConstant(kTRUE);

    ws->var("sigmaSig2")->setVal(inputCBG.sigmaSig2);
    ws->var("coeffExp")->setVal(inputCBG.coeffExp);
    ws->var("coeffExp2")->setVal(inputCBG.coeffExp2);

    ws->var("sigmaSig2_pp")->setVal(inputCBG.sigmaSig2_pp);
    ws->var("coeffExp_pp")->setVal(inputCBG.coeffExp_pp);
    ws->var("coeffExp2_pp")->setVal(inputCBG.coeffExp2_pp);
  }

  ws->var("Jpsi_Mass")->setRange("jpsi",2.6,3.55);
  ws->var("Jpsi_Mass")->setRange("psip",3.4,4.2);

  if (prefitMass) {

    if (!yrange.compare("0.0-1.6")) {
	ws->var("coeffGaus")->setVal(0.0);
	ws->var("sigmaSig1")->setVal(0.0);
	ws->var("coeffGaus")->setError(0.0);
	ws->var("sigmaSig1")->setError(0.0);
	ws->var("coeffGaus")->setConstant(true);
	ws->var("sigmaSig1")->setConstant(true);
      }
    fitM_pre = ws->pdf("jpsiMassPDFSim")->fitTo(*redDataCutSim,Extended(1),Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(6),Range("jpsi"));
    fitM_pre->Print("v");
    string resultFN_pre = dirPre + "_rap" + yrange + "_pT" + prange + "_cent" + crange + "_dPhi" + phirange + "_preFitResult.root";
    TFile resultPreF(resultFN_pre.c_str(),"RECREATE");
    resultPreF.cd(); cout << "fitM_pre->Write(jpsiMassPDFSim):" << fitM_pre->Write("jpsiMassPDFSim") << endl;resultPreF.Close();
    
    if (!crange.compare("0-100")){
      inputCBG.sigmaSig1Err = ws->var("sigmaSig1")->getError();
      inputCBG.coeffGausErr = ws->var("coeffGaus")->getError();
      
      ws->var("sigmaSig1")->setConstant(kTRUE);
      ws->var("coeffGaus")->setConstant(kTRUE);
    }

    inputCBG.sigmaSig1Err_pp = ws->var("sigmaSig1_pp")->getError();
    inputCBG.coeffGausErr_pp = ws->var("coeffGaus_pp")->getError();
    
    ws->var("sigmaSig1_pp")->setConstant(kTRUE);
    ws->var("coeffGaus_pp")->setConstant(kTRUE);
  }


  fitM = ws->pdf("sigMassPDFSim")->fitTo(*redDataCutSim,Extended(1),Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(6));
  resultF.cd(); cout << "fitM->Write(sigMassPDFSim):" << fitM->Write("sigMassPDFSim") << endl;
  string fitFileStr = dirPre + "_rap" + yrange + "_pT" + prange + "_cent" + crange + "_FitResult.txt";
  RooArgSet * pars = ws->pdf("sigMassPDFSim")->getParameters(ws->data("data"));
  pars->writeToFile(fitFileStr.c_str());

  // *** Draw mass plot before do ctau fit
  for (int i=0; i<nFiles; i++) {
    char name[100] = {0};
    RooPlot *mframe_wob = ws->var("Jpsi_Mass")->frame();
    redDataCut[i]->plotOn(mframe_wob,DataError(RooAbsData::SumW2),XErrorSize(0),MarkerSize(1),Binning(rbm));
    titlestr = "2D fit for" + partTit + "muons (mass projection), p_{T} = " + prange  + "_dPhi" + phirange + " GeV/c and |y| = " + yrange;
    mframe_wob->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c^{2})");
    mframe_wob->GetXaxis()->CenterTitle(1);
    double max = mframe_wob->GetMaximum() * 1.3;
    double min = 0.0;
    if (i==0)
      sprintf(name,"NBkg_HI");
    else
      sprintf(name,"NBkg_pp");
    min = ws->var(name)->getVal()/80.0 * 0.7;
    cout << "MINIMUM: " << min << endl;
    mframe_wob->SetMaximum(max);
    if (i==0)
      sprintf(name,"sigMassPDF_HI");
    else
      sprintf(name,"sigMassPDF_pp");

    ws->pdf(name)->plotOn(mframe_wob,DrawOption("F"),FillColor(kBlack),FillStyle(3354),Normalization(redDataCut[i]->sumEntries(),RooAbsReal::NumEvent));
    
    ws->pdf(name)->plotOn(mframe_wob,Components(mBkgFunct[i].c_str()),DrawOption("F"),FillColor(kAzure-9),FillStyle(1001),Normalization(redDataCut[i]->sumEntries(),RooAbsReal::NumEvent));
    
    ws->pdf(name)->plotOn(mframe_wob,Components(mBkgFunct[i].c_str()),LineColor(kBlue),LineStyle(7),LineWidth(5),Normalization(redDataCut[i]->sumEntries(),RooAbsReal::NumEvent));

    ws->pdf(name)->plotOn(mframe_wob,LineColor(kBlack),LineWidth(2),Normalization(redDataCut[i]->sumEntries(),RooAbsReal::NumEvent));
    redDataCut[i]->plotOn(mframe_wob,DataError(RooAbsData::SumW2),XErrorSize(0),MarkerSize(1),Binning(rbm));
    
    TH1 *hdata = redDataCut[i]->createHistogram("hdata",*ws->var("Jpsi_Mass"),Binning(rbm));

    if (i==1)
      ws->pdf("sigMassPDF_mix")->plotOn(mframe_wob,LineColor(kRed),LineStyle(2),LineWidth(5),Normalization(1.06*redDataCut[i]->sumEntries(),RooAbsReal::NumEvent));

    // *** Calculate chi2/nDof for mass fitting
    int nBins = hdata->GetNbinsX();
    RooHist *hpull; hpull = mframe_wob->pullHist(); hpull->SetName("hpullhist");
    double Chi2 = 0;
    int nFullBinsPull = 0;
    double *ypull = hpull->GetY();
    for (unsigned int j=0; j < nBins; j++) {
      if (hdata->GetBinContent(j) == 0) continue;
      nFullBinsPull++;
      Chi2 = Chi2 + pow(ypull[j],2);
    }
    double UnNormChi2 = Chi2;
    int nFitParam = fitM->floatParsFinal().getSize();
    int Dof = nFullBinsPull - nFitParam;
    Chi2 /= (nFullBinsPull - nFitParam);

    // *** Check in narrower signal region NSig
    ws->var("Jpsi_Mass")->setRange("sigpeak",2.9,3.3);
    RooAbsReal *inteAll = ws->pdf(mJpsiFunct.c_str())->createIntegral(RooArgSet(*ws->var("Jpsi_Mass")),NormSet(RooArgSet(*ws->var("Jpsi_Mass"))));
    //  RooAbsReal *inteSig = ws->pdf(mSigFunct.c_str())->createIntegral(RooArgSet(*ws->var("Jpsi_Mass")),Range("sigpeak"),NormSet(RooArgSet(*ws->var("Jpsi_Mass"))),Range("sigpeak"));

    if (i==0) sprintf(name,"PbPb");
    else if (i==1) sprintf(name,"pp");

    TCanvas c00; c00.cd();

    if (logScale) {
      c00.SetLogy(1);
      mframe_wob->SetMinimum(0.5*min);
      mframe_wob->SetMaximum(8.0*max);
    }

    mframe_wob->Draw();
    //    cout << "Set MINIMUM: " << mframe_wob->GetMinimum() << endl;
    t->SetTextSize(0.05);
    t->DrawLatex(0.15,0.90,"CMS Preliminary");
    if (i==0)
      t->DrawLatex(0.15,0.82,"PbPb  #sqrt{s_{NN}} = 2.76 TeV");
    else
      t->DrawLatex(0.15,0.82,"pp  #sqrt{s} = 2.76 TeV");
    t->SetTextSize(0.04);
    if (i==0)
      t->DrawLatex(0.15,0.75,"L_{int} = 150 #mub^{-1}");
    else
      t->DrawLatex(0.15,0.75,"L_{int} = 231 nb^{-1}");
    t->SetTextSize(0.035);
    if (i==0) {
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
    if (i==0) t->DrawLatex(0.56,0.50,reduceDS);
    if (pmin==6.5)
      sprintf(reduceDS,"%.1f < p_{T} < %.0f GeV/c",pmin,pmax);
    else if (pmax==6.5)
      sprintf(reduceDS,"%.0f < p_{T} < %.1f GeV/c",pmin,pmax);
    else
      sprintf(reduceDS,"%.0f < p_{T} < %.0f GeV/c",pmin,pmax);
    if (i==0) t->DrawLatex(0.56,0.45,reduceDS);
    // sprintf(reduceDS,"%.2f < |#phi_{J/#psi}-#Psi_{RP}| < %.2f",psmin,psmax);
    // t->DrawLatex(0.17,0.59,reduceDS);
    // sprintf(reduceDS,"EP: %s",rpmethod.c_str());
    // t->DrawLatex(0.17,0.54,reduceDS);
    //  t->SetTextSize(0.04);
    sprintf(reduceDS,"#chi^{2}/dof = %0.1f/%d",UnNormChi2,Dof);
    //    t->DrawLatex(0.62,0.90,reduceDS);

    double DoubleRatio_tmp = computeRatio(*(ws->var("fracP_HI")),*(ws->var("fracP_pp")));
    double ErrDoubleRatio_tmp = computeRatioError(*(ws->var("fracP_HI")),*(ws->var("fracP_pp")),fitM->correlation("fracP_HI","fracP_pp"));

    if (i==0) {
      sprintf(name,"NJpsi_HI");
      sprintf(reduceDS,"N_{J/#psi}: %0.0f #pm %0.0f",ws->var(name)->getVal(),ws->var(name)->getError());
      t->DrawLatex(0.62,0.85,reduceDS);
      sprintf(reduceDS,"#chi_{#psi(2S)}: %0.3f #pm %0.3f",DoubleRatio_tmp, ErrDoubleRatio_tmp);
      t->DrawLatex(0.62,0.80,reduceDS);
      sprintf(name,"fracP_HI");
      sprintf(reduceDS,"R_{#psi(2S)}^{PbPb}: %0.3f #pm %0.3f",ws->var(name)->getVal(),ws->var(name)->getError());
      t->DrawLatex(0.62,0.75,reduceDS);

      double coeffGaus = ws->var("coeffGaus")->getVal();
      double sigmaSig1 = ws->var("sigmaSig1")->getVal();
      double sigmaSig2 = ws->var("sigmaSig2")->getVal();
      double ErrcoeffGaus = ws->var("coeffGaus")->getError();
      double ErrsigmaSig1 = ws->var("sigmaSig1")->getError();
      double ErrsigmaSig2 = ws->var("sigmaSig2")->getError();

      if (ErrcoeffGaus == 0.0) ErrcoeffGaus = inputCBG.coeffGausErr;
      if (ErrsigmaSig1 == 0.0) ErrsigmaSig1 = inputCBG.sigmaSig1Err;
      if (ErrsigmaSig2 == 0.0) ErrsigmaSig2 = inputCBG.sigmaSig2Err;
      
      double resol = sigmaSig2;
      double Errresol = ErrsigmaSig2;

      if  (found!=string::npos) {
	resol = sqrt( coeffGaus*pow(sigmaSig1,2) + (1-coeffGaus)*pow(sigmaSig2,2) );
	Errresol = (0.5/resol) * sqrt( pow(sigmaSig1*coeffGaus*ErrsigmaSig1,2) +
				       pow(sigmaSig2*(1-coeffGaus)*ErrsigmaSig2,2) +
				       pow(0.5*(pow(sigmaSig1,2)-pow(sigmaSig2,2))*ErrcoeffGaus,2) );
      }

      if (Errresol<0.5) 
	sprintf(reduceDS,"#sigma = (%0.0f #pm %0.0f) MeV/c^{2}",resol*1000.0,1.0);
      else
	sprintf(reduceDS,"#sigma = (%0.0f #pm %0.0f) MeV/c^{2}",resol*1000.0,Errresol*1000.0);

      //      sprintf(reduceDS,"#sigma = (%0.0f #pm %0.0f) MeV/c^{2}",ws->var("sigmaSig2")->getVal()*1000,ws->var("sigmaSig2")->getError()*1000);
    }
    else {
      sprintf(name,"NJpsi_pp");
      sprintf(reduceDS,"N_{J/#psi}: %0.0f #pm %0.0f",ws->var(name)->getVal(),ws->var(name)->getError());
      t->DrawLatex(0.62,0.85,reduceDS);
      sprintf(reduceDS,"#chi_{#psi(2S)}: %0.3f #pm %0.3f",DoubleRatio_tmp, ErrDoubleRatio_tmp);
      t->DrawLatex(0.62,0.80,reduceDS);
      sprintf(name,"fracP_pp");
      sprintf(reduceDS,"R_{#psi(2S)}^{pp}: %0.3f #pm %0.3f",ws->var(name)->getVal(),ws->var(name)->getError());
      t->DrawLatex(0.62,0.75,reduceDS);

      double coeffGaus_pp = ws->var("coeffGaus_pp")->getVal();
      double sigmaSig1_pp = ws->var("sigmaSig1_pp")->getVal();
      double sigmaSig2_pp = ws->var("sigmaSig2_pp")->getVal();
      double ErrcoeffGaus_pp = ws->var("coeffGaus_pp")->getError();
      double ErrsigmaSig1_pp = ws->var("sigmaSig1_pp")->getError();
      double ErrsigmaSig2_pp = ws->var("sigmaSig2_pp")->getError();

      if (ErrcoeffGaus_pp == 0.0) ErrcoeffGaus_pp = inputCBG.coeffGausErr_pp;
      if (ErrsigmaSig1_pp == 0.0) ErrsigmaSig1_pp = inputCBG.sigmaSig1Err_pp;
      if (ErrsigmaSig2_pp == 0.0) ErrsigmaSig2_pp = inputCBG.sigmaSig2Err_pp;
      
      double resol_pp = sigmaSig2_pp;
      double Errresol_pp = ErrsigmaSig2_pp;

      if  (found!=string::npos) {
	resol_pp = sqrt( coeffGaus_pp*pow(sigmaSig1_pp,2) + (1-coeffGaus_pp)*pow(sigmaSig2_pp,2) );
	Errresol_pp = (0.5/resol_pp) * sqrt( pow(sigmaSig1_pp*coeffGaus_pp*ErrsigmaSig1_pp,2) +
				       pow(sigmaSig2_pp*(1-coeffGaus_pp)*ErrsigmaSig2_pp,2) +
				       pow(0.5*(pow(sigmaSig1_pp,2)-pow(sigmaSig2_pp,2))*ErrcoeffGaus_pp,2) );
      }

      if (Errresol_pp<0.5) 
	sprintf(reduceDS,"#sigma = (%0.0f #pm %0.0f) MeV/c^{2}",resol_pp*1000.0,1.0);
      else
	sprintf(reduceDS,"#sigma = (%0.0f #pm %0.0f) MeV/c^{2}",resol_pp*1000.0,Errresol_pp*1000.0);

      //      sprintf(reduceDS,"#sigma = (%0.0f #pm %0.0f) MeV/c^{2}",ws->var("sigmaSig2_pp")->getVal()*1000,ws->var("sigmaSig2_pp")->getError()*1000);
    }
    t->DrawLatex(0.62,0.70,reduceDS);

    TLegend *leg1;
    if (i==0)
      leg1 = new TLegend(0.63,0.53,0.92,0.68,NULL,"brNDC");
    else
      leg1 = new TLegend(0.63,0.48,0.92,0.68,NULL,"brNDC");
    leg1->SetFillStyle(0); leg1->SetBorderSize(0); leg1->SetShadowColor(0); leg1->SetMargin(0.2);
    leg1->AddEntry(gfake1,"data","p");
    leg1->AddEntry(&hfake21,"total fit","lf");
    leg1->AddEntry(&hfake11,"background","lf");
    if(i==1)
      leg1->AddEntry(&hfake22,"PbPb shape","l");
    
    leg1->Draw("same");
    if (i==0)
      sprintf(name,"PbPb");
    else
      sprintf(name,"pp");

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

    titlestr = dirPre + "_rap" + yrange_str  + "_pT" + prange_str + "_cent" + crange + "_massfit_wob" + name + ".pdf";
    c00.SaveAs(titlestr.c_str());
    delete hpull;
    delete hdata;
  }

  Double_t NJpsi_fin[2] = {0,0};
  Double_t ErrNJpsi_fin[2] = {0,0};
  Double_t NBkg_fin[2] = {0,0};
  Double_t ErrNBkg_fin[2] = {0,0};
  Double_t fracP_fin[2] = {0,0};
  Double_t ErrFracP_fin[2] = {0,0};
  Double_t DoubleRatio_fin = 0.0;
  Double_t ErrDoubleRatio_fin = 0.0;

  for (int i=0; i<nFiles; i++) {
    char name[100] = {0};

    if (i==0) {
      sprintf(name,"NJpsi_HI");
      NJpsi_fin[i]= ws->var(name)->getVal();
      ErrNJpsi_fin[i] = ws->var(name)->getError();
      sprintf(name,"NBkg_HI");
      NBkg_fin[i] = ws->var(name)->getVal();
      ErrNBkg_fin[i] = ws->var(name)->getError();
      DoubleRatio_fin = computeRatio(*(ws->var("fracP_HI")),*(ws->var("fracP_pp")));
      ErrDoubleRatio_fin = computeRatioError(*(ws->var("fracP_HI")),*(ws->var("fracP_pp")),fitM->correlation("fracP_HI","fracP_pp"));
      sprintf(name,"DoubleRatio");
      sprintf(name,"fracP_HI");
      fracP_fin[i]= ws->function(name)->getVal();
      ErrFracP_fin[i] = ws->var(name)->getError();
    }
    else {
      sprintf(name,"NJpsi_pp");
      NJpsi_fin[i]= ws->var(name)->getVal();
      ErrNJpsi_fin[i] = ws->var(name)->getError();
      sprintf(name,"NBkg_pp");
      NBkg_fin[i] = ws->var(name)->getVal();
      ErrNBkg_fin[i] = ws->var(name)->getError();
      sprintf(name,"fracP_pp");
      fracP_fin[i]= ws->var(name)->getVal();
      ErrFracP_fin[i] = ws->var(name)->getError();
    }
  }

  //  Double_t NJpsi_fin = ws->var("NJpsi")->getVal();
  //  Double_t ErrNJpsi_fin = ws->var("NJpsi")->getError();
  //  Double_t NBkg_fin = ws->var("NBkg")->getVal();
  //  Double_t ErrNBkg_fin = ws->var("NBkg")->getError();


  // To check values of fit parameters
  cout << endl << "J/psi yields:" << endl;
  for (int i=0; i<nFiles; i++) {
    cout << "NJpsi" << i <<" :       Fit :"  << NJpsi_fin[i] << " +/- " << ErrNJpsi_fin[i] << endl;
  }
  cout << "Double Ratio: " << DoubleRatio_fin << " +/- " << ErrDoubleRatio_fin << endl;

  titlestr = dirPre + "_rap" + yrange + "_pT" + prange + "_cent" + crange + "_dPhi" + phirange + ".txt";

  double theNLL=0;
  theNLL = fitM->minNll();

  ofstream outputFile(titlestr.c_str());
  if (!outputFile.good()) {cout << "Fail to open result file." << endl; return 1;}

  for (int i=0; i<nFiles; i++) {
    outputFile << "NJpsi" << i+1 << " " << NJpsi_fin[i]                << " " << ErrNJpsi_fin[i] << "\n"
	       << "NBkg" << i+1 << " " << NBkg_fin[i]                << " " << ErrNBkg_fin[i] << "\n";
    if (i==0)
      outputFile << "DoubleRatio " << DoubleRatio_fin          << " " << ErrDoubleRatio_fin << "\n";
    outputFile << "Rpsi2S" << i+1 << " " << fracP_fin[i]          << " " << ErrFracP_fin[i] << "\n";
  }
  if (!mBkgFunct[0].compare("expFunct_HI")) {
    outputFile
      << "coeffExp "   << ws->var("coeffExp")->getVal()      << " " << ws->var("coeffExp")->getError() << "\n"
      << "coeffExp_pp "   << ws->var("coeffExp_pp")->getVal()      << " " << ws->var("coeffExp_pp")->getError() << "\n";
  }
  if (!mBkgFunct[0].compare("twoExpFunct_HI")) {
    outputFile
      << "coeffExp "      << ws->var("coeffExp")->getVal()      << " " << ws->var("coeffExp")->getError() << "\n"
      << "coeffExp2 "     << ws->var("coeffExp2")->getVal()     << " " << ws->var("coeffExp2")->getError() << "\n"
      << "cutx "          << ws->var("cutx")->getVal()          << " " << ws->var("cutx")->getError() << "\n"
      << "coeffExp_pp "      << ws->var("coeffExp_pp")->getVal()      << " " << ws->var("coeffExp_pp")->getError() << "\n";
      // << "coeffExp2_pp "     << ws->var("coeffExp2_pp")->getVal()     << " " << ws->var("coeffExp2_pp")->getError() << "\n"
      // << "cutx_pp "          << ws->var("cutx_pp")->getVal()          << " " << ws->var("cutx_pp")->getError() << "\n";
  }
  if (!mBkgFunct[0].compare("pol1_HI")) {
    outputFile
      << "coeffPol1 "  << ws->var("coeffPol1_HI")->getVal()     << " " << ws->var("coeffPol1_HI")->getError() << "\n"
      << "coeffPol1_pp "  << ws->var("coeffPol1_pp")->getVal()     << " " << ws->var("coeffPol1_pp")->getError() << "\n";
  }
  if (!mBkgFunct[0].compare("pol2_HI")) {
    outputFile
      << "coeffPol1 "  << ws->var("coeffPol1")->getVal()     << " " << ws->var("coeffPol1")->getError() << "\n"
      << "coeffPol1_pp "  << ws->var("coeffPol1_pp")->getVal()     << " " << ws->var("coeffPol1_pp")->getError() << "\n"
      << "coeffPol2 "  << ws->var("coeffPol2")->getVal()     << " " << ws->var("coeffPol2")->getError() << "\n"
      << "coeffPol2_pp "  << ws->var("coeffPol2_pp")->getVal()     << " " << ws->var("coeffPol2_pp")->getError() << "\n";
  }
  outputFile
    << "coeffGaus " << ws->var("coeffGaus")->getVal()    << " " << ws->var("coeffGaus")->getError() << "\n"
    << "meanSig1 "  << ws->var("meanSig1")->getVal()     << " " << ws->var("meanSig1")->getError() << "\n"
    //    << "meanSig1P "    << ws->var("meanSig1P")->getVal()     << " " << ws->var("meanSig1P")->getError() << "\n"
    << "sigmaSig1 " << ws->var("sigmaSig1")->getVal()    << " " << ws->var("sigmaSig1")->getError() << "\n"
    << "sigmaSig2 " << ws->var("sigmaSig2")->getVal()    << " " << ws->var("sigmaSig2")->getError()<< "\n"
    << "alpha "     << ws->var("alpha")->getVal()        << " " << ws->var("alpha")->getError() << "\n"
    << "enneW "     << ws->var("enneW")->getVal()        << " " << ws->var("enneW")->getError() << "\n"
    << "coeffGaus_pp " << ws->var("coeffGaus_pp")->getVal()   << " " << ws->var("coeffGaus_pp")->getError() << "\n"
    << "meanSig1_pp "  << ws->var("meanSig1_pp")->getVal()    << " " << ws->var("meanSig1_pp")->getError() << "\n"
    //    << "meanSig1P_pp "    << ws->var("meanSig1P")->getVal()     << " " << ws->var("meanSig1P")->getError() << "\n"
    << "sigmaSig1_pp " << ws->var("sigmaSig1_pp")->getVal()   << " " << ws->var("sigmaSig1_pp")->getError() << "\n"
    << "sigmaSig2_pp " << ws->var("sigmaSig2_pp")->getVal()   << " " << ws->var("sigmaSig2_pp")->getError()<< "\n"
    << "alpha_pp "     << ws->var("alpha_pp")->getVal()       << " " << ws->var("alpha_pp")->getError() << "\n"
    << "enneW_pp "     << ws->var("enneW_pp")->getVal()       << " " << ws->var("enneW_pp")->getError() << "\n"
    << "NLL "       << theNLL                              << endl;

  for (int i=0; i<FileName.size(); i++) {
    fInData[i]->Close();
  }
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
  ws->factory("Chebychev::pol1_HI(Jpsi_Mass,{coeffPol1[-0.8,-150.,150.]})");
  ws->factory("Chebychev::pol1_pp(Jpsi_Mass,{coeffPol1_pp[-0.8,-150.,150.]})");
  // 2nd order polynomial
  ws->factory("Chebychev::pol2_HI(Jpsi_Mass,{coeffPol1, coeffPol2[0.00,-15.,15.]})");
  ws->factory("Chebychev::pol2_pp(Jpsi_Mass,{coeffPol1_pp, coeffPol2_pp[0.00,-15.,15.]})");
  // Exponential
  ws->factory("Exponential::expFunct_HI(Jpsi_Mass,coeffExp[-0.1,-3.,1.])");
  ws->factory("Exponential::expFunct_pp(Jpsi_Mass,coeffExp_pp[-0.1,-3.,1.])");

  // 2nd Exponential
  ws->factory("Exponential::expFunct2_HI(Jpsi_Mass,coeffExp2[-0.2.,-3.,1.])");
  ws->factory("Exponential::expFunct2_pp(Jpsi_Mass,coeffExp2_pp[-0.2.,-3.,1.])");


  // Sum of two exponentials
  RooRealVar cutx("cutx","cutx",3.4,3.2,3.5); cutx.setConstant(false); ws->import(cutx); //Region below than this cut will use cTau1, above than it will use cTau2. determined by free fit.
  RooRealVar cutx_pp("cutx_pp","cutx_pp",3.4,3.2,3.5); cutx_pp.setConstant(false); ws->import(cutx_pp); //Region below than this cut will use cTau1, above than it will use cTau2. determined by free fit.

  RooGenericPdf twoExpFunct_HI("twoExpFunct_HI","twoExpFunct_HI","1.0/abs(@1)*(exp(@0*@1)*(@0<@3)+exp(@0*@2 - @3*(@2-@1) )*(@0>=@3))", RooArgList(*(ws->var("Jpsi_Mass")), *(ws->var("coeffExp")), *(ws->var("coeffExp2")), *(ws->var("cutx"))));
  ws->import(twoExpFunct_HI);

  RooGenericPdf twoExpFunct_pp("twoExpFunct_pp","twoExpFunct_pp","1.0/abs(@1)*(exp(@0*@1)*(@0<@3)+exp(@0*@2 - @3*(@2-@1) )*(@0>=@3))", RooArgList(*(ws->var("Jpsi_Mass")), *(ws->var("coeffExp_pp")), *(ws->var("coeffExp2_pp")), *(ws->var("cutx_pp"))));
  ws->import(twoExpFunct_pp);



  // ws->factory("SUM::twoExpFunct_HI(coeffBkg_HI[0.9,0.1,1.0]*expFunct_HI,expFunct2_HI)");
  // ws->factory("SUM::twoExpFunct_pp(coeffBkg_pp[0.9,0.1,1.0]*expFunct_pp,expFunct2_pp)");


  return;
}

void defineMassSig(RooWorkspace *ws) {
  //////// Candidates for signal
  // Normal gaussians
  ws->factory("Gaussian::signalG1(Jpsi_Mass,meanSig1[3.0969,3.05,3.15],sigmaSig1[0.08,0.06,0.2])");
  ws->factory("Gaussian::signalG1_pp(Jpsi_Mass,meanSig1_pp[3.0969,3.05,3.15],sigmaSig1_pp[0.08,0.06,0.2])");

  // Crystall Ball
  ws->factory("CBShape::signalCB2WN(Jpsi_Mass,meanSig1,sigmaSig2[0.03,0.01,0.060],alpha[1.,0.,3.],enneW[5.,1.,50.])");
  ws->factory("CBShape::signalCB2WN_pp(Jpsi_Mass,meanSig1_pp,sigmaSig2_pp[0.03,0.01,0.06],alpha_pp[1.,0.,3.],enneW_pp[5.,1.,50.])");

  // Fix Jpsi-psi' mass difference
  //  RooFormulaVar meanSig1P("meanSig1P","@0+0.58917",RooArgList(*(ws->var("meanSig1")))); ws->import(meanSig1P);
  // Fix Jpsi-psi' mass ratio
  RooFormulaVar meanSig1P("meanSig1P","@0*1.1902",RooArgList(*(ws->var("meanSig1")))); ws->import(meanSig1P);
  // Fix resolution scale: sigma_MJpsi/MJpsi = sigma_Mpsi'/Mpsi'
  RooFormulaVar sigmaSig1P("sigmaSig1P","@0*1.1902",RooArgList(*(ws->var("sigmaSig1")))); ws->import(sigmaSig1P);
  RooFormulaVar sigmaSig2P("sigmaSig2P","@0*1.1902",RooArgList(*(ws->var("sigmaSig2")))); ws->import(sigmaSig2P);
  ws->factory("Gaussian::signalG1P(Jpsi_Mass,meanSig1P,sigmaSig1P)");
  ws->factory("CBShape::signalCB2WNP(Jpsi_Mass,meanSig1P,sigmaSig2P,alpha,enneW)");

  RooFormulaVar meanSig1P_pp("meanSig1P_pp","@0*1.1902",RooArgList(*(ws->var("meanSig1_pp")))); ws->import(meanSig1P_pp);
  RooFormulaVar sigmaSig1P_pp("sigmaSig1P_pp","@0*1.1902",RooArgList(*(ws->var("sigmaSig1_pp")))); ws->import(sigmaSig1P_pp);
  RooFormulaVar sigmaSig2P_pp("sigmaSig2P_pp","@0*1.1902",RooArgList(*(ws->var("sigmaSig2_pp")))); ws->import(sigmaSig2P_pp);
  ws->factory("Gaussian::signalG1P_pp(Jpsi_Mass,meanSig1P_pp,sigmaSig1P_pp)");
  ws->factory("CBShape::signalCB2WNP_pp(Jpsi_Mass,meanSig1P_pp,sigmaSig2P_pp,alpha_pp,enneW_pp)");

  //////// Sum of signal functions
  // Sum of gaussian 1 and crystall ball 2 with wide n
  ws->factory("SUM::sigCB2WNG1(coeffGaus[0.1,0.05,0.93]*signalG1,signalCB2WN)");
  ws->factory("SUM::sigCB2WNG1P(coeffGaus*signalG1P,signalCB2WNP)");

  ws->factory("SUM::sigCB2WNG1_pp(coeffGaus_pp[0.1,0.05,0.93]*signalG1_pp,signalCB2WN_pp)");
  ws->factory("SUM::sigCB2WNG1P_pp(coeffGaus_pp*signalG1P_pp,signalCB2WNP_pp)");

  // ws->factory("SUM::sigCBGJpsiPsiP(fracP[0.1,0.0,1.0]*sigCB2WNG1P,sigCB2WNG1)");
  // ws->factory("SUM::sigCBJpsiPsiP(fracP*signalCB2WNP,signalCB2WN)");

  // ws->factory("SUM::sigCBGJpsiPsiP_pp(fracP_pp[0.1,0.0,1.0]*sigCB2WNG1P_pp,sigCB2WNG1_pp)");
  // ws->factory("SUM::sigCBJpsiPsiP_pp(fracP_pp*signalCB2WNP_pp,signalCB2WN_pp)");

  return;
}

double computeRatio(RooRealVar& x, RooRealVar& y) {
  assert(y.getVal() != 0);
  return x.getVal()/y.getVal();
}

double computeRatioError(RooRealVar& x, RooRealVar& y, 
			 double correlation) {
  double err2 = x.getError()*x.getError()/x.getVal()/x.getVal() +
    y.getError()*y.getError()/y.getVal()/y.getVal() - 
    2.*x.getError()*y.getError()/x.getVal()/y.getVal()*correlation;

  return fabs(computeRatio(x,y))*sqrt(err2);
}
