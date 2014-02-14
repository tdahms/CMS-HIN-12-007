#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <TROOT.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TPaveStats.h>
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
#include "RooDataHist.h"

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

#include <RooStats/ModelConfig.h>

using namespace RooFit;
using namespace std;

bool superImpose = true;
bool analyticBlifetime = true;
bool narrowSideband = false;
bool oneGaussianResol = false;

void getOptRange(string &ran,float *min,float *max);
void setWSRange(RooWorkspace *ws);
void defineMassBkg(RooWorkspace *ws);
void defineMassSig(RooWorkspace *ws);


int main(int argc, char* argv[]) {
  gROOT->Macro("./rootlogon.C");

  string FileName, FileNameMC1, FileNameMC2;
  string mBkgFunct, mJpsiFunct, mPsiPFunct;
  bool prefitMass = false;
  int  isGG = 0;
  string prange, lrange, yrange, crange;
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
  bool fixAlpha = true;
  bool fixN = true;
  bool fixGwidth = true; // only used in PbPb

  // *** Check options
  for (int i=1; i<argc; ++i) {
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
      case 'b':
	cutNonPrompt = atoi(argv[i+1]);
	cout << "remove non-prompt J/psi: " << cutNonPrompt << endl;
	break;
      case 'a':
	fixAlpha = false;
	cout << "fix CB alpha: " << fixAlpha << endl;
	break;
      case 'n':
	fixN = false;
	cout << "fix CB n: " << fixN << endl;
	break;
      case 'g':
	fixGwidth = false;
	cout << "fix wideFactor: " << fixGwidth << endl;
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

  float pmin=0, pmax=0, ymin=0, ymax=0, cmin=0, cmax=0;
  getOptRange(prange,&pmin,&pmax);
  getOptRange(crange,&cmin,&cmax);
  getOptRange(yrange,&ymin,&ymax);

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

  string fix_str="";
  if (!fixAlpha)
    fix_str+="_freeAlpha";
  if (!fixN)
    fix_str+="_freeN";
  if (!fixGwidth)
    fix_str+="_freeGwidth";
  if (fix_str != "")
    fix_str+="_";

  // *** TFile for saving fitting results
  string resultFN;
  resultFN = dirPre + "_" + mBkgFunct + "_rap" + yrange_str + "_pT" + prange_str + "_cent" + crange + fix_str + "fitResult.root";
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
      sprintf(reduceDS,"Jpsi_Pt>%.2f && Jpsi_Pt<%.2f && abs(Jpsi_Y)>%.2f && abs(Jpsi_Y)<%.2f && Jpsi_Ct<0.05",pmin,pmax,ymin,ymax);
    else if (prange=="6.5-30.0" && yrange=="0.0-1.6")
      sprintf(reduceDS,"Jpsi_Pt>%.2f && Jpsi_Pt<%.2f && abs(Jpsi_Y)>%.2f && abs(Jpsi_Y)<%.2f && Jpsi_Ct<0.04",pmin,pmax,ymin,ymax);
    else if (prange=="3.0-30.0" && yrange=="1.6-2.4")
      sprintf(reduceDS,"Jpsi_Pt>%.2f && Jpsi_Pt<%.2f && abs(Jpsi_Y)>%.2f && abs(Jpsi_Y)<%.2f && Jpsi_Ct<0.10",pmin,pmax,ymin,ymax);
    else if (prange=="3.0-6.5" && yrange=="1.6-2.4")
      sprintf(reduceDS,"Jpsi_Pt>%.2f && Jpsi_Pt<%.2f && abs(Jpsi_Y)>%.2f && abs(Jpsi_Y)<%.2f && Jpsi_Ct<0.11",pmin,pmax,ymin,ymax);
    else if (prange=="6.5-30.0" && yrange=="1.6-2.4")
      sprintf(reduceDS,"Jpsi_Pt>%.2f && Jpsi_Pt<%.2f && abs(Jpsi_Y)>%.2f && abs(Jpsi_Y)<%.2f && Jpsi_Ct<0.06",pmin,pmax,ymin,ymax);
  }
  else if (cutNonPrompt && isPbPb) {
    if (prange=="6.5-30.0" && yrange=="0.0-2.4")
      sprintf(reduceDS,"Jpsi_Pt>%.2f && Jpsi_Pt<%.2f && abs(Jpsi_Y)>%.2f && abs(Jpsi_Y)<%.2f && Jpsi_Ct<0.05",pmin,pmax,ymin,ymax);
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

  // Binning for invariant mass distribution
  RooBinning rbm(2.2,4.2);
  int nbins = 100;
  bool highStats=true;
  if (highStats)
    rbm.addUniform(nbins,2.2,4.2);
  else {
    // rbm.addUniform(14,2.2,2.760);
    // rbm.addUniform(24,2.760,3.240);
    // rbm.addUniform(24,3.240,4.2);
    rbm.addUniform(11,2.2,2.860); // 60 MeV
    rbm.addUniform(19,2.860,3.240); // 20 MeV
    rbm.addUniform(6,3.240,3.600); // 60 MeV
    rbm.addUniform(4,3.600,3.720); // 30 MeV
    rbm.addUniform(8,3.720,4.2); // 60 MeV
  }

  ws->var("Jpsi_Mass")->setBinning(rbm);

  RooDataHist *binData = new RooDataHist("binData","binData",RooArgSet( *(ws->var("Jpsi_Mass")) ), *redData);
  cout << "DATA :: N events to fit: " << binData->sumEntries() << endl;

  // SYSTEMATICS 1 (very sidebands)
  RooDataSet *redDataSB;
  if (narrowSideband) {
   redDataSB = (RooDataSet*) redData->reduce("Jpsi_Mass<2.8 || Jpsi_Mass>3.4");
  } else {
   redDataSB = (RooDataSet*) redData->reduce("Jpsi_Mass<2.9 || Jpsi_Mass>3.3");
  }
  //  RooDataHist *binDataSB = new RooDataHist("binDataSB","Data distribution for background",RooArgSet( *(ws->var("Jpsi_Mass")) ),*redDataSB);
  //  RooDataSet *redDataSIG = (RooDataSet*)redData->reduce("Jpsi_Mass > 2.9 && Jpsi_Mass < 3.3");

  // *** Define PDFs with parameters (mass and ctau)
  // Just so RooFit does not crash on Ubuntu
  RooRealVar aa("aa","aa",0.5,-1,1);
  RooRealVar ab("ab","ab",-0.5,-1,1);
  RooChebychev tmpPol("tmpPol","tmpPol",*(ws->var("Jpsi_Mass")),RooArgSet(aa,ab));

  // J/psi mass parameterization
  defineMassBkg(ws);
  defineMassSig(ws);

  char funct[125];
  string partTit, partFile;
  if (isGG == 0) { partTit = "glb-glb"; partFile = "GG"; }
  else if (isGG == 1) { partTit = "glb-trk"; partFile = "GT"; }
  else { partTit = "all"; partFile = "ALL"; }

  // Global TLatex, TH1, TGraph objects for drawing
  TLatex *lCMS = new TLatex();
  lCMS->SetNDC(); lCMS->SetTextAlign(12);
  TLatex *lLumi = new TLatex();
  lLumi->SetNDC(); lLumi->SetTextAlign(12);
  TLatex *lRap = new TLatex();
  lRap->SetNDC(); lRap->SetTextAlign(12);
  TLatex *lPt = new TLatex();
  lPt->SetNDC(); lPt->SetTextAlign(12);
  TLatex *lNLL = new TLatex();
  lNLL->SetNDC(); lNLL->SetTextAlign(12);
  TLatex *lChi = new TLatex();
  lChi->SetNDC(); lChi->SetTextAlign(12);
  TLatex *lPval = new TLatex();
  lPval->SetNDC(); lPval->SetTextAlign(12);
  TLatex *lNJpsi = new TLatex();
  lNJpsi->SetNDC(); lNJpsi->SetTextAlign(12);
  TLatex *lRpsi = new TLatex();
  lRpsi->SetNDC(); lRpsi->SetTextAlign(12);
  TLatex *lNpsiP = new TLatex();
  lNpsiP->SetNDC(); lNpsiP->SetTextAlign(12);
  TLatex *lSigCB = new TLatex();
  lSigCB->SetNDC(); lSigCB->SetTextAlign(12);
  TLatex *lSigG = new TLatex();
  lSigG->SetNDC(); lSigG->SetTextAlign(12);
  TLatex *lNG = new TLatex();
  lNG->SetNDC(); lNG->SetTextAlign(12);
  TLatex *lSigma = new TLatex();
  lSigma->SetNDC(); lSigma->SetTextAlign(12);
  TLatex *lAlpha = new TLatex();
  lAlpha->SetNDC(); lAlpha->SetTextAlign(12);
  TLatex *lN = new TLatex();
  lN->SetNDC(); lN->SetTextAlign(12);
  TLatex *lFG = new TLatex();
  lFG->SetNDC(); lFG->SetTextAlign(12);

  Double_t fx[2], fy[2], fex[2], fey[2];
  TGraphErrors *gfake1 = new TGraphErrors(2,fx,fy,fex,fey);
  gfake1->SetMarkerStyle(20); gfake1->SetMarkerSize(0.8);
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
  RooFitResult *fitP = NULL;

  RooRealVar *NJpsi  = new RooRealVar("NJpsi","J/psi yield",0.5*binData->sumEntries(),0.0,2.0*binData->sumEntries());  ws->import(*NJpsi);
  RooRealVar *NBkg  = new RooRealVar("NBkg","Brackground yield", 0.5*binData->sumEntries(),0.0,2.0*binData->sumEntries());   ws->import(*NBkg);
  RooRealVar *NBkg2  = new RooRealVar("NBkg2","Second brackground yield",0.5*binData->sumEntries(),0.0,2.0*binData->sumEntries());   ws->import(*NBkg2);

  if (fitRatio) {
    RooRealVar *fracP  = new RooRealVar("fracP","psi(2S) fraction" ,0.01);
    fracP->setConstant(false); ws->import(*fracP);
    RooFormulaVar *NPsiP = new RooFormulaVar("NPsiP", "@0*@1", RooArgList(*(ws->var("NJpsi")),*(ws->var("fracP"))));  ws->import(*NPsiP);
  }
  else {
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
	ws->var("alpha")->setVal(1.489);
	ws->var("enne")->setVal(1.730);
	ws->var("wideFactor")->setVal(1.901);
      }
      else if (prange=="6.5-30.0" && yrange=="0.0-1.6") {
	ws->var("alpha")->setVal(1.516);
	ws->var("wideFactor")->setVal(1.699);
	ws->var("enne")->setVal(1.683);
      }
      else if (prange=="3.0-30.0" && yrange=="1.6-2.4") {
	ws->var("alpha")->setVal(1.880);
	ws->var("enne")->setVal(1.430);
	ws->var("wideFactor")->setVal(1.503);
      }
      else if (prange=="3.0-6.5" && yrange=="1.6-2.4") {
	ws->var("alpha")->setVal(1.860);
	ws->var("enne")->setVal(1.330);
	ws->var("wideFactor")->setVal(1.540);
      }
      else if (prange=="6.5-30.0" && yrange=="1.6-2.4") {
	ws->var("alpha")->setVal(2.187);
	ws->var("enne")->setVal(1.295);
	ws->var("wideFactor")->setVal(1.830);
      }
      else {
	ws->var("alpha")->setVal(2.0);
	ws->var("enne")->setVal(1.4);
	ws->var("wideFactor")->setVal(2.0);
      }
    }
    else {
      if (prange=="6.5-30.0" && yrange=="0.0-2.4") {
	ws->var("alpha")->setVal(1.693);
	ws->var("enne")->setVal(1.709);
      }
      else if (prange=="6.5-30.0" && yrange=="0.0-1.6") {
	ws->var("alpha")->setVal(1.710);
	ws->var("enne")->setVal(1.646);
      }
      else if (prange=="3.0-30.0" && yrange=="1.6-2.4") {
	ws->var("alpha")->setVal(2.079);
	ws->var("enne")->setVal(1.390);
      }
      else if (prange=="3.0-6.5" && yrange=="1.6-2.4") {
	ws->var("alpha")->setVal(2.000);
	ws->var("enne")->setVal(1.350);
      }
      else if (prange=="6.5-30.0" && yrange=="1.6-2.4") {
	ws->var("alpha")->setVal(2.161);
	ws->var("enne")->setVal(1.370);
      }
      else {
	ws->var("alpha")->setVal(2.0);
	ws->var("enne")->setVal(1.4);
      }
    }
    if (fixAlpha)
      ws->var("alpha")->setConstant(kTRUE);
    if (fixN)
      ws->var("enne")->setConstant(kTRUE);
    if (isPbPb && fixGwidth)
      ws->var("wideFactor")->setConstant(kTRUE);
  }

  // 20140128: seed with CB fit parameters
  if ( found!=string::npos ) {
    string inputFNcb;
    inputFNcb =  dirPre2 + "_" + mBkgFunct + "_rap" + yrange_str + "_pT" + prange_str + "_cent" + crange + fix_str + "allVars.txt";
    RooArgSet *set = ws->pdf("sigMassPDF")->getParameters(*(ws->var("Jpsi_Mass")));
    //    set->Print("v");
    set->readFromFile(inputFNcb.c_str());
    //    cout << set->getRealValue("sigmaSig1",0,1) << endl;
    cout << "Import variable values from CB fit: " << inputFNcb << endl;
    ws->import(*set);
    // set->Print("v");
    // ws->Print("v");

    /* set it to MC instead
    // try to set wideFactor to the result from the pp fit.
      if (isPbPb) {
      string dirPrePp = dirPre;
      string strPbPb = "frac";
      size_t found = dirPre.find(strPbPb);
      if (found!=string::npos) {
      dirPrePp.replace(found,strPbPb.length(),"ppFrac");
      //	cout << "pp string: " << dirPrePp << endl; 
      inputFNcb =  dirPrePp + "_" + mBkgFunct + "_rap" + yrange_str + "_pT" + prange_str + "_cent0-100_allVars.txt";
      cout << "pp file name: " << inputFNcb << endl;
      RooArgSet *ppSet = ws->pdf("sigMassPDF")->getParameters(*(ws->var("Jpsi_Mass")));
      //	ppSet->Print("v");
      //	cout << "read from file: " << ppSet->readFromFile(inputFNcb.c_str(),"","",1) << endl;
      ppSet->readFromFile(inputFNcb.c_str());
      //	ppSet->Print("v");
      cout << "setting wideFactor to " << ppSet->getRealValue("wideFactor",0.0,kFALSE) << " from pp."<< endl;
      ws->var("wideFactor")->setVal(ppSet->getRealValue("wideFactor",0.0,kFALSE));
      //	ws->var("wideFactor")->setConstant(kTRUE);
      }
      else
      cout << "String " << strPbPb << "not found in " << dirPre << endl;
      }
    */
  }

  fitM = ws->pdf("sigMassPDF")->fitTo(*redData,Extended(1),Hesse(1),Minos(0),Save(1),SumW2Error(0),NumCPU(8),PrintEvalErrors(0),Verbose(0));
  fitM->printMultiline(cout,1,1,"");

  int edmStatus = 1;
  if (fitM->edm()<1.0e-03)
    edmStatus = 0;

  int covStatus = 1;
  if (fitM->covQual()==3)
    covStatus = 0;

  int fitStatus = fitM->status();
  if (found!=string::npos)
    cout << "CBG_" << mBkgFunct << "_rap" << yrange_str << "_pT" << prange_str << "_cent" << crange << fix_str << endl;
  else
    cout << "CB_" << mBkgFunct << "_rap" << yrange_str << "_pT" << prange_str << "_cent" << crange << fix_str << endl;

  cout << "---FIT result summary: " << edmStatus << covStatus << fitStatus;
  for (unsigned int i=0; i<fitM->numStatusHistory(); ++i) {
    cout << fitM->statusCodeHistory(i);
  }
  cout << " ---" << endl;

  int nFitParam = fitM->floatParsFinal().getSize();
  resultF.cd();
  fitM->Write();
  cout << "fitM->Write() into: " << resultFN << endl;

  // *** Draw mass plot
  RooPlot *mframe = ws->var("Jpsi_Mass")->frame();
  redData->plotOn(mframe,DataError(RooAbsData::SumW2),XErrorSize(0),MarkerSize(0.8),Binning(rbm));
  titlestr = "2D fit for" + partTit + "muons (mass projection), p_{T} = " + prange + " GeV/c and |y| = " + yrange;
  mframe->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  mframe->GetXaxis()->CenterTitle(1);
  double max = mframe->GetMaximum() * 1.3;
  double min = 0.0;
  mframe->SetMaximum(max);
  min = ws->var("NBkg")->getVal()/(double(nbins)) * 0.7;
  RooHist *hpull_jpsi;
  RooHist *hpull_psip;
  RooHist *hpull; 

  if (fitSubRange) {
    ws->pdf("jpsiMassPDF")->plotOn(mframe,DrawOption("F"),FillColor(kBlack),FillStyle(3354));
    ws->pdf("jpsiMassPDF")->plotOn(mframe,Components(mBkgFunct.c_str()),LineColor(kBlue),LineStyle(7),LineWidth(5));
    ws->pdf("jpsiMassPDF")->plotOn(mframe,LineColor(kBlack),LineWidth(2));
    hpull_jpsi = mframe->pullHist(0,0,true); hpull_jpsi->SetName("hpullhist");

    ws->pdf("psipMassPDF")->plotOn(mframe,DrawOption("F"),FillColor(kGreen+2),FillStyle(3354));
    
    ws->pdf("psipMassPDF")->plotOn(mframe,Components(mBkgFunctP.c_str()),LineColor(kRed),LineStyle(7),LineWidth(5));
    ws->pdf("psipMassPDF")->plotOn(mframe,LineColor(kGreen+2),LineWidth(2));
    hpull_psip = mframe->pullHist(0,0,true); hpull_psip->SetName("hpullhistP");
    hpull_psip->SetMarkerColor(kGreen+2);
  }
  else {
    ws->pdf("sigMassPDF")->plotOn(mframe,VisualizeError(*fitM,1,kFALSE),FillColor(kGray),Normalization(redData->sumEntries(),RooAbsReal::NumEvent));
    ws->pdf("sigMassPDF")->plotOn(mframe,LineColor(kBlack),LineWidth(2),Normalization(redData->sumEntries(),RooAbsReal::NumEvent));

    hpull = mframe->pullHist(0,0,true); hpull->SetName("hpullhist");
    
    ws->pdf("sigMassPDF")->plotOn(mframe,Components(mBkgFunct.c_str()),VisualizeError(*fitM,1,kFALSE),FillColor(kCyan),Normalization(redData->sumEntries(),RooAbsReal::NumEvent));
    ws->pdf("sigMassPDF")->plotOn(mframe,Components(mBkgFunct.c_str()),LineColor(kBlue),LineStyle(kDashed),LineWidth(2),Normalization(redData->sumEntries(),RooAbsReal::NumEvent));

    if (!isPaper && found!=string::npos) { 
      ws->pdf("sigMassPDF")->plotOn(mframe,VisualizeError(*fitM,1,kFALSE),FillColor(kRed-9),Components((mBkgFunct+",signalCB1,signalCB1P").c_str()));
      ws->pdf("sigMassPDF")->plotOn(mframe,Components((mBkgFunct+",signalCB1,signalCB1P").c_str()),LineColor(kRed),LineStyle(kDashed),LineWidth(2));
      ws->pdf("sigMassPDF")->plotOn(mframe,VisualizeError(*fitM,1,kFALSE),FillColor(kGreen-9),Components((mBkgFunct+",signalG2,signalG2P").c_str()));
      ws->pdf("sigMassPDF")->plotOn(mframe,Components((mBkgFunct+",signalG2,signalG2P").c_str()),LineColor(kGreen+2),LineStyle(kDashed),LineWidth(2));
    }
  }
  redData->plotOn(mframe,DataError(RooAbsData::SumW2),XErrorSize(0),MarkerSize(0.8),Binning(rbm));

  TH1 *hdata = redData->createHistogram("hdata",*ws->var("Jpsi_Mass"),Binning(rbm));

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
  
  // Make a histogram with the pull distribution
  TH1F *hpull_proj = new TH1F("hpull_proj","hpull_proj;Pull;Events",100,-10,10);
  for (unsigned int i=0; i < nBins; ++i) {
    hpull_proj->Fill(ypull[i]);
  }

  double Chi2P = 0;
  int nFullBinsPullP = 0;

  if (fitSubRange) {
    int minBins = hdata->GetXaxis()->FindBin(3.4);
    for (unsigned int i=minBins; i < nBins; ++i) {
      if (hdata->GetBinContent(i+1) == 0) continue;
      nFullBinsPullP++;
      Chi2P = Chi2P + pow(ypull_psip[i-minBins],2);
    }

    nBins = hdata->GetXaxis()->FindBin(3.5);

    for (unsigned int i=0; i < nBins; ++i) {
      if (hdata->GetBinContent(i+1) == 0) continue;
      nFullBinsPull++;
      Chi2 = Chi2 + pow(ypull[i],2);
    }
  }
  else {
    for (unsigned int i=0; i < nBins; ++i) {
      if (hdata->GetBinContent(i+1) == 0) continue;
      nFullBinsPull++;
      Chi2 = Chi2 + pow(ypull[i],2);
    }
  }

  double UnNormChi2 = Chi2;
  int Dof = nFullBinsPull - nFitParam;
  if (Dof!=0)
    Chi2 /= Dof;

  double theNLL=0;
  theNLL = fitM->minNll();

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

  TF1 *f0 = new TF1("f0","0",2.2,4.2);
  f0->SetLineWidth(2);
  f0->SetLineStyle(1);
  f0->SetLineColor(kBlack);

  TF1 *f1 = new TF1("f1","2",2.2,4.2);
  f1->SetLineWidth(1);
  f1->SetLineStyle(2);
  f1->SetLineColor(kBlack);

  TF1 *f2 = new TF1("f2","-2",2.2,4.2);
  f2->SetLineWidth(1);
  f2->SetLineStyle(2);
  f2->SetLineColor(kBlack);

  mframepull->addObject(f0,"same");
  mframepull->addObject(f1,"same");
  mframepull->addObject(f2,"same");

  // *** Check in narrower signal region NSig
  //  ws->var("Jpsi_Mass")->setRange("sigpeak",2.9,3.3);

  //  RooAbsReal *inteAll = ws->pdf(mJpsiFunct.c_str())->createIntegral(RooArgSet(*ws->var("Jpsi_Mass")),NormSet(RooArgSet(*ws->var("Jpsi_Mass"))));
  //  RooAbsReal *inteSig = ws->pdf(mSigFunct.c_str())->createIntegral(RooArgSet(*ws->var("Jpsi_Mass")),Range("sigpeak"),NormSet(RooArgSet(*ws->var("Jpsi_Mass"))),Range("sigpeak"));

  // mframe->SetMinimum(0.03*max);
  // mframe->SetMaximum(5.0*max);
  TCanvas c00; c00.cd();

  TPad *pad1 = new TPad("pad1","This is pad1",0.0,0.35,1.00,1.00);
  pad1->SetBottomMargin(0); 
  TPad *pad2 = new TPad("pad2","This is pad2",0.0,0.00,1.00,0.35);
  pad2->SetTopMargin(0);  pad2->SetBottomMargin(0.24); 
  TPad *pad3 = new TPad("pad3","This is pad3",0.16,0.52,0.40,0.75);
  pad3->SetLeftMargin(0.2);
  pad3->SetBottomMargin(0.2);
  if (!isPaper) {
    pad1->Draw();
    pad2->Draw();
    pad3->Draw();
    pad2->cd(); mframepull->Draw();
    pad1->cd();
  }

  // if (!isPbPb)
  //   ws->pdf("sigMassPDF_mix")->plotOn(mframe,LineColor(kRed),LineStyle(2),LineWidth(4),Normalization(redData->sumEntries(),RooAbsReal::NumEvent));

  if (logScale) {
    if (isPaper)
      c00.SetLogy(1);
    else
      pad1->SetLogy(1);

    mframe->SetMinimum(0.2*min);
    if ( (yrange == "0.0-1.6" || yrange == "0.0-2.4") && isPbPb) {
      if (crange == "20-40") {
	mframe->SetMaximum(3.5*max);
	mframe->SetMinimum(0.5*min);
      }
      else if ( crange == "40-100" ) {
	mframe->SetMaximum(60*max);
	mframe->SetMinimum(0.1*min);
      }
      else{
	mframe->SetMaximum(2.5*max);
	mframe->SetMinimum(0.8*min);
      }
    }
    else if (yrange == "1.6-2.4" && isPbPb) {
      if ( crange == "40-100" ) {
	mframe->SetMaximum(10.0*max);
	mframe->SetMinimum(0.5*min);
      }
      else if ( crange == "20-40" ) {
	mframe->SetMaximum(2.0*max);
	mframe->SetMinimum(0.9*min);
      }
      else {
	mframe->SetMaximum(1.5*max);
	mframe->SetMinimum(0.95*min);
      }
    }
    else {
      mframe->SetMaximum(20*max);
      mframe->SetMinimum(0.5*min);
    }
  }
  else if (zoom) {
    if (min>0)
      mframe->SetMinimum(0.5*min);

    if (yrange == "1.6-2.4") {
      if (prange == "3.0-6.5" && crange == "0-20")
       	mframe->SetMaximum(max);
      else if (prange == "3.0-30.0" && crange == "20-40")
       	mframe->SetMaximum(max*0.3);
      else if (prange == "3.0-30.0" && crange=="40-100")
       	mframe->SetMaximum(max*0.2);
      else 
	mframe->SetMaximum(max*0.5);
    }
    else
      mframe->SetMaximum(max*0.15);
  }

  if (!highStats) {
    max = mframe->GetMaximum();
    min = mframe->GetMinimum();
    mframe->SetMaximum(2.0*max);
    mframe->SetMinimum(2.0*min);
  }

  lCMS->SetTextSize(0.05);
  if (isPbPb)
    lCMS->SetText(0.17,0.90,"CMS PbPb #sqrt{s_{NN}} = 2.76 TeV");
  else
    lCMS->SetText(0.17,0.90,"CMS pp #sqrt{s} = 2.76 TeV");
  mframe->addObject(lCMS,"");

  lLumi->SetTextSize(0.035);
  if (isPbPb)
    lLumi->SetText(0.17,0.83,"L_{int} = 150 #mub^{-1}");
  else
    lLumi->SetText(0.17,0.83,"L_{int} = 5.3 pb^{-1}");
  mframe->addObject(lLumi,"");

  lRap->SetTextSize(0.035);
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
  lRap->SetText(0.17,0.78,reduceDS);
  mframe->addObject(lRap,"");

  lPt->SetTextSize(0.035);
  if (pmin==6.5)
    sprintf(reduceDS,"%.1f < p_{T} < %.0f GeV/c",pmin,pmax);
  else if (pmax==6.5)
    sprintf(reduceDS,"%.0f < p_{T} < %.1f GeV/c",pmin,pmax);
  else
    sprintf(reduceDS,"%.0f < p_{T} < %.0f GeV/c",pmin,pmax);
  
  lPt->SetText(0.17,0.73,reduceDS);
  mframe->addObject(lPt,"");

  lNLL->SetTextSize(0.035);
  if (!isPaper) {
    sprintf(reduceDS,"Min. NLL = %0.2f",theNLL);
    lNLL->SetText(0.17,0.68,reduceDS);
    mframe->addObject(lNLL,"");
  }

  if (fitSubRange) {
    sprintf(reduceDS,"#chi^{2}/dof (J/#psi) = %0.1f/%d   #chi^{2}/dof (#psi(2S)) = %0.1f/%d",UnNormChi2,Dof,UnNormChi2P,DofP);
    lChi->SetText(0.45,0.90,reduceDS);
  }
  else {
    if (!isPaper) {
      sprintf(reduceDS,"#chi^{2}/dof = %0.1f/%d",UnNormChi2,Dof);
      lChi->SetText(0.62,0.90,reduceDS);
      mframe->addObject(lChi,"");
      sprintf(reduceDS,"p-value = %0.4f",TMath::Prob(UnNormChi2,Dof));
      lPval->SetText(0.62,0.85,reduceDS);
      mframe->addObject(lPval,"");
    }
  }

  if (ws->var("NJpsi")->hasAsymError() && abs(-1.0*ws->var("NJpsi")->getErrorLo()/ws->var("NJpsi")->getErrorHi() - 1)>0.1)
    sprintf(reduceDS,"N_{J/#psi} = %0.0f^{+%0.0f}_{%0.0f}",ws->var("NJpsi")->getVal(),ws->var("NJpsi")->getErrorHi(),ws->var("NJpsi")->getErrorLo());
  else
    sprintf(reduceDS,"N_{J/#psi} = %0.0f #pm %0.0f",ws->var("NJpsi")->getVal(),ws->var("NJpsi")->getError());
  lNJpsi->SetText(0.62,0.80,reduceDS);
  mframe->addObject(lNJpsi,"");

  if (fitRatio) {
    if (ws->var("fracP")->hasAsymError() && abs(-1.0*ws->var("fracP")->getErrorLo()/ws->var("fracP")->getErrorHi() - 1)>0.1)
      sprintf(reduceDS,"R_{#psi(2S)} = %0.3f^{+%0.3f}_{%0.3f}",ws->var("fracP")->getVal(),ws->var("fracP")->getErrorHi(),ws->var("fracP")->getErrorLo());
    else
      sprintf(reduceDS,"R_{#psi(2S)} = %0.3f #pm %0.3f",ws->var("fracP")->getVal(),ws->var("fracP")->getError());
    lRpsi->SetText(0.62,0.75,reduceDS);
    mframe->addObject(lRpsi,"");

    double NJpsi = ws->var("NJpsi")->getVal();
    double fracP = ws->var("fracP")->getVal();
    double errNJpsi = ws->var("NJpsi")->getError();
    double errFracP = ws->var("fracP")->getError();
    double corr = fitM->correlation( *ws->var("NJpsi") , *ws->var("fracP") );
    std::cout << "Correlation between NJpsi and fracP: " << corr << std::endl;
    double Npsi2S = ws->function("NPsiP")->getVal();
    double errNpsi2S = sqrt( pow(errNJpsi/NJpsi,2) + pow(errFracP/fracP,2) + 2*errNJpsi*errFracP*corr/Npsi2S )*Npsi2S;

    // if (ws->var("fracP")->hasAsymError() && abs(-1.0*ws->var("fracP")->getErrorLo()/ws->var("fracP")->getErrorHi() - 1)>0.1)
    //   sprintf(reduceDS,"N_{#psi(2S)} = %0.1f^{+%0.1f}_{%0.1f}",ws->function("NPsiP")->getVal(),ws->var("fracP")->getErrorHi()*ws->var("NJpsi")->getVal(),ws->var("fracP")->getErrorLo()*ws->var("NJpsi")->getVal());
    // else
    sprintf(reduceDS,"N_{#psi(2S)} = %0.1f #pm %0.1f",Npsi2S,errNpsi2S);
    lNpsiP->SetText(0.62,0.70,reduceDS);
  }
  else {
    if (ws->var("NPsiP")->hasAsymError() && abs(-1.0*ws->var("NPsiP")->getErrorLo()/ws->var("NPsiP")->getErrorHi() - 1)>0.1)
      sprintf(reduceDS,"N_{#psi(2S)} = %0.0f^{+%0.0f}_{%0.0f}",ws->var("NPsiP")->getVal(),ws->var("NPsiP")->getErrorHi(),ws->var("NPsiP")->getErrorLo());
    else
      sprintf(reduceDS,"N_{#psi(2S)} = %0.1f #pm %0.1f",ws->var("NPsiP")->getVal(),ws->var("NPsiP")->getError());
    lNpsiP->SetText(0.62,0.75,reduceDS);
  }
  mframe->addObject(lNpsiP,"");

  double coeffGaus = ws->var("coeffGaus")->getVal();
  double sigmaSig1 = ws->var("sigmaSig1")->getVal();
  double sigmaSig2 = ws->function("sigmaSig2")->getVal();
  double ErrcoeffGaus = ws->var("coeffGaus")->getError();
  double ErrsigmaSig1 = ws->var("sigmaSig1")->getError();
  double ErrsigmaSig2 = sigmaSig1*ws->var("wideFactor")->getError();

  // if (ErrcoeffGaus == 0.0) ErrcoeffGaus = inputCBG.coeffGausErr;
  // if (ErrsigmaSig1 == 0.0) ErrsigmaSig1 = inputCBG.sigmaSig1Err;
  // if (ErrsigmaSig2 == 0.0) ErrsigmaSig2 = inputCBG.sigmaSig2Err;

  if (!isPaper) {
    if (ErrsigmaSig1<0.0005) 
      sprintf(reduceDS,"#sigma_{CB} = (%0.0f #pm %0.0f) MeV/c^{2}",sigmaSig1*1000.0,1.0);
    else
      sprintf(reduceDS,"#sigma_{CB} = (%0.0f #pm %0.0f) MeV/c^{2}",sigmaSig1*1000.0,ErrsigmaSig1*1000.0);
    lSigCB->SetText(0.62,0.65,reduceDS);

    if  (found!=string::npos) {
      sprintf(reduceDS,"#sigma_{G} = %0.0f MeV/c^{2}",sigmaSig2*1000.0);
      lSigG->SetText(0.62,0.60,reduceDS);
      if (ws->var("wideFactor")->isConstant())
	sprintf(reduceDS,"n_{G} = %0.2f (fixed)",ws->var("wideFactor")->getVal());
      else if (ws->var("wideFactor")->hasAsymError() && abs(-1.0*ws->var("wideFactor")->getErrorLo()/ws->var("wideFactor")->getErrorHi() - 1)>0.1)
	sprintf(reduceDS,"n_{G} = %0.2f^{+%0.2f}_{%0.2f}",ws->var("wideFactor")->getVal(),ws->var("wideFactor")->getErrorHi(),ws->var("wideFactor")->getErrorLo());
      else
	sprintf(reduceDS,"n_{G} = %0.2f #pm %0.2f",ws->var("wideFactor")->getVal(),ws->var("wideFactor")->getError());
      lNG->SetText(0.62,0.55,reduceDS);
    }
    mframe->addObject(lSigCB,"");
    mframe->addObject(lSigG,"");
    mframe->addObject(lNG,"");
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

      lSigma->SetText(0.62,0.50,reduceDS);
      mframe->addObject(lSigma,"");

      if (ws->var("alpha")->isConstant())
	sprintf(reduceDS,"#alpha = %0.2f (fixed)",ws->var("alpha")->getVal());
      else if (ws->var("alpha")->hasAsymError() && abs(-1.0*ws->var("alpha")->getErrorLo()/ws->var("alpha")->getErrorHi() - 1)>0.1)
	sprintf(reduceDS,"#alpha = %0.2f^{+%0.2f}_{%0.2f}",ws->var("alpha")->getVal(),ws->var("alpha")->getErrorHi(),ws->var("alpha")->getErrorLo());
      else
	sprintf(reduceDS,"#alpha = %0.2f #pm %0.2f",ws->var("alpha")->getVal(),ws->var("alpha")->getError());
      lAlpha->SetText(0.62,0.45,reduceDS);
      if (ws->var("enne")->isConstant())
	sprintf(reduceDS,"n = %0.2f (fixed)",ws->var("enne")->getVal());
      else if (ws->var("enne")->hasAsymError() && abs(-1.0*ws->var("enne")->getErrorLo()/ws->var("enne")->getErrorHi() - 1)>0.1)
	sprintf(reduceDS,"n = %0.2f^{+%0.2f}_{%0.2f}",ws->var("enne")->getVal(),ws->var("enne")->getErrorHi(),ws->var("enne")->getErrorLo());
      else
	sprintf(reduceDS,"n = %0.2f #pm %0.2f",ws->var("enne")->getVal(),ws->var("enne")->getError());
      lN->SetText(0.62,0.40,reduceDS);
      if (ws->var("coeffGaus")->hasAsymError() && abs(-1.0*ws->var("coeffGaus")->getErrorLo()/ws->var("coeffGaus")->getErrorHi() - 1)>0.1)
	sprintf(reduceDS,"f_{G} = %0.3f^{+%0.3f}_{%0.3f}",ws->var("coeffGaus")->getVal(),ws->var("coeffGaus")->getErrorHi(),ws->var("coeffGaus")->getErrorLo());
      else
	sprintf(reduceDS,"f_{G} = %0.2f #pm %0.2f",ws->var("coeffGaus")->getVal(),ws->var("coeffGaus")->getError());
      lFG->SetText(0.62,0.35,reduceDS);
      mframe->addObject(lFG,"");
    }
    else {
      if (ws->var("alpha")->isConstant())
	sprintf(reduceDS,"#alpha = %0.2f (fixed)",ws->var("alpha")->getVal());
      else
	sprintf(reduceDS,"#alpha = %0.2f #pm %0.2f",ws->var("alpha")->getVal(),ws->var("alpha")->getError());
      lAlpha->SetText(0.62,0.60,reduceDS);
      if (ws->var("enne")->isConstant())
	sprintf(reduceDS,"n = %0.2f (fixed)",ws->var("enne")->getVal());
      else
	sprintf(reduceDS,"n = %0.2f #pm %0.2f",ws->var("enne")->getVal(),ws->var("enne")->getError());
      lN->SetText(0.62,0.55,reduceDS);
    }

    mframe->addObject(lAlpha,"");
    mframe->addObject(lN,"");
  }

  TLegend *leg1;
  if (fitSubRange)
    leg1 = new TLegend(0.63,0.53,0.92,0.72,NULL,"brNDC");
  else
    leg1 = new TLegend(0.63,0.53,0.92,0.68,NULL,"brNDC");

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
    mframe->addObject(leg1,"sa,e");

  mframe->Draw();
  if (!isPaper) {
    pad3->cd();
    pad3->SetFillStyle(0);
    gStyle->SetOptStat("mr");
    hpull_proj->Draw();
    hpull_proj->GetXaxis()->CenterTitle(1);
    hpull_proj->GetXaxis()->SetRangeUser(-5,5);
    hpull_proj->GetXaxis()->SetLabelSize(0.07);
    hpull_proj->GetYaxis()->SetLabelSize(0.07);
    hpull_proj->GetXaxis()->SetTitleSize(0.07);
    hpull_proj->GetYaxis()->SetTitleSize(0.07);
    pad3->Update();
    TPaveStats *st = (TPaveStats*) hpull_proj->FindObject("stats");
    st->SetX1NDC(0.64);
    st->SetX2NDC(0.99);
    st->SetY1NDC(0.75);
    st->SetY2NDC(0.99);
    st->SetTextSize(0.065);
    st->SetTextFont(42);
  }

  titlestr = dirPre + "_" + mBkgFunct + "_rap" + yrange_str  + "_pT" + prange_str + "_cent" + crange + fix_str +  "massfit.pdf";
  c00.SaveAs(titlestr.c_str());
  resultF.cd();
  mframe->Write();
  mframepull->Write();
  hpull_proj->Write();

  // prepare models for CLs calculations
  RooStats::ModelConfig *model = new RooStats::ModelConfig("model");
  model->SetWorkSpace(ws);
  model->SetPdf("sigMassPDF");
  model->SetObservables(RooArgSet("Jpsi_Mass"));
  RooArgSet *pars = ws->pdf("sigMassPDF")->getParameters(redData);
  //  pars.assignFast(pars);
  pars->Print("v");
  RooArgSet *nuisances = RooARgSet(pars);
  nuisances->remove(ws->var("fracP"));
  model->SetNuisanceParameters(nuisances);
  RooArgSet poi = RooArgSet(ws->var("fracP"));
  model->SetParametersOfInterest(poi);
  model->SetSnapshot(pars);

  // prepare models for CLs calculations
  RooStats::ModelConfig *model_null = (RooStats::ModelConfig*)model->Clone();
  //  model_null->SetWorkSpace(ws);
  //  model_null->SetPdf("sigMassPDF");
  //  model_null->SetObservables(RooArgSet("Jpsi_Mass"));
  //  RooArgSet *pars = ws->pdf("sigMassPDF")->getParameters(redData);
  //  pars.assignFast(pars);
  //  pars->Print("v");
  //  RooArgSet *nuisances = RooARgSet(pars);
  //  nuisances->remove(ws->var("fracP"));
  //  model->SetNuisanceParameters(nuisances);
  RooArgSet poiNew = RooArgSet(ws->var("fracP"));
  double oldVal = ws->var("fracP")->getVal();
  ws->var("fracP")->setVal(0.0);
  model_null->SetSnapshot(poiNew);
  ws->var("fracP")->setVal(oldVal);
  // to be continued...

  string fname = dirPre + "_" + mBkgFunct + "_rap" + yrange_str  + "_pT" + prange_str + "_cent" + crange + fix_str + "Workspace.root";
  ws->writeToFile(fname.c_str(),kTRUE);
  fname = dirPre + "_" + mBkgFunct + "_rap" + yrange_str  + "_pT" + prange_str + "_cent" + crange + fix_str + "allVars.txt";
  ws->allVars().writeToFile(fname.c_str());

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

  titlestr = dirPre + "_" + mBkgFunct + "_rap" + yrange_str + "_pT" + prange_str + "_cent" + crange + fix_str + "_fitResult.txt";

  ofstream outputFile(titlestr.c_str());
  if (!outputFile.good()) {cout << "Fail to open result file." << endl; return 1;}
  outputFile.precision(10);
  //  outputFile.setf( std::ios::fixed, std::ios::floatfield)
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
  if (!mBkgFunct.compare("gausBkg")) {
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
    << "sigmaSig1 "    << ws->var("sigmaSig1")->getVal()    << " " << ws->var("sigmaSig1")->getError() << "\n"
    << "wideFactor "   << ws->var("wideFactor")->getVal()   << " " << ws->var("wideFactor")->getError() << "\n"
    << "sigmaSig2 "    << ws->function("sigmaSig2")->getVal() << " " << ws->function("sigmaSig2")->getVal()*ws->var("wideFactor")->getError()/ws->var("wideFactor")->getVal() << "\n"
    << "alpha "        << ws->var("alpha")->getVal()        << " " << ws->var("alpha")->getError() << "\n"
    << "enne "         << ws->var("enne")->getVal()         << " " << ws->var("enne")->getError() << "\n"
    << "Mean(Pull) "   << hpull_proj->GetMean()             << "\n"
    << "RMS(Pull) "    << hpull_proj->GetRMS()              << "\n"
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
  ws->factory("Gaussian::gausBkg(Jpsi_Mass,meanBkg[0.0,0.0,10.0],sigmaBkg[1.0,0.5,5.0])");

  return;
}

void defineMassSig(RooWorkspace *ws) {
  // narrow Gauss
  ws->factory("Gaussian::signalG1(Jpsi_Mass,meanSig1[3.0969,3.05,3.15],sigmaSig1[0.03,0.005,0.080])");
  // narrow CB
  ws->factory("CBShape::signalCB1(Jpsi_Mass,meanSig1,sigmaSig1,alpha[1.0,0.0,3.0],enne[5.0,1.0,50.0])");

  RooRealVar wideFactor("wideFactor","wideFactor",2.0,1.0,6.0);ws->import(wideFactor);
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
