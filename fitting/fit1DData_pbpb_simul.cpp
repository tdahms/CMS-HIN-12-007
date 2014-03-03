#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

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
void defineMassBkgHI(RooWorkspace *ws);
void defineMassSigHI(RooWorkspace *ws);


int main(int argc, char* argv[]) {
  gROOT->Macro("./rootlogon.C");

  vector<string> FileName;
  vector<string> mBkgFunct;
  string mJpsiFunct, mPsiPFunct;
  int  isGG = 0;
  string prange, lrange, yrange, crange;
  string dirPre;
  string rpmethod = "etHF";
  bool logScale=false;
  bool isPbPb=true;
  bool isPaper = false;
  bool fixCBtoMC=true;
  bool cutNonPrompt=false;
  bool fixAlpha = true;
  bool fixN = true;
  bool fixGwidth = true; // only used in PbPb
  bool shareShape = false;

  // *** Check options
  for (int i=1; i<argc; ++i) {
    char *tmpargv = argv[i];
    switch (tmpargv[0]) {
    case '-':{
      switch (tmpargv[1]) {
      case 'f':
	FileName.push_back(argv[i+1]);
	FileName.push_back(argv[i+2]);
	if ( argv[i+3][0] != '-' ) {
	  FileName.push_back(argv[i+3]);
	  if ( argv[i+4][0] != '-' ) 
	    FileName.push_back(argv[i+4]);
	}
	cout << "Fitted PbPb data file: " << FileName.front() << endl;
	cout << "Fitted pp data file: " << FileName.back() << endl;
	break;
      case 'v':
	mJpsiFunct = argv[i+1];
	mPsiPFunct = argv[i+2];
	if ( argv[i+3][0] != '-' ) {
	  mBkgFunct.push_back(argv[i+3]);
	  if ( argv[i+4][0] != '-' ) {
	    mBkgFunct.push_back(argv[i+4]);
	    if ( argv[i+5][0] != '-' ) {
	      mBkgFunct.push_back(argv[i+5]);
	      if ( argv[i+6][0] != '-' ) 
		mBkgFunct.push_back(argv[i+6]);
	    }
	  }
	}
       	cout << "Mass J/psi function: " << mJpsiFunct << endl;
	cout << "Mass psi(2S) function: " << mPsiPFunct << endl;
	cout << "Mass PbPb background function: " << mBkgFunct.front();
	for (unsigned int j=1;j<mBkgFunct.size()-1; ++j) {cout << " " << mBkgFunct.at(j);}
	cout << endl;
	cout << "Mass pp background function: " << mBkgFunct.back() << endl;
	break;
      case 'd':
	dirPre = argv[i+1];
	cout << "Prefix for all result files: " << dirPre << endl;
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
      case 'l':
	logScale = atoi(argv[i+1]);
	cout << "plot on log scale: " << logScale << endl;
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
      case 's':
	shareShape = atoi(argv[i+1]);
	cout << "Use same shape for pp and PbPb: " << shareShape << endl;
	break;
      case 'z':
	isPaper = atoi(argv[i+1]);
	cout << "paper plots: " << isPaper << endl;
	break;
      }
    }
    }
  }// End check options
 
  float pmin=0, pmax=0, ymin=0, ymax=0, cmin=0, cmax=0;
  getOptRange(prange,&pmin,&pmax);
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

  string fix_str="_";
  if (!fixAlpha)
    fix_str+="freeAlpha_";
  if (!fixN)
    fix_str+="freeN_";
  if (!fixGwidth)
    fix_str+="freeGwidth_";
  
  const unsigned int nFiles = FileName.size();
  if (nFiles<2) {
    cout << "Need at least two files for simultaneous fit. Only " << nFiles << " provided." << endl;
    return 1;
  }

  // *** TFile for saving fitting results
  string resultFN;
  resultFN = dirPre + "_PbPb" + mBkgFunct.front();
  for (unsigned int i=1; i<nFiles-1;++i) {
    resultFN +=  "_" + mBkgFunct.at(i);
  }
  resultFN += "_pp" + mBkgFunct.back() + "_rap" + yrange_str + "_pT" + prange_str + "_cent" + crange + fix_str + "fitResult.root";
  TFile resultF(resultFN.c_str(),"RECREATE");


  vector<string> varSuffix;
  if (nFiles==2) {
    varSuffix.push_back("HI");
    varSuffix.push_back("pp");
  }
  else if(nFiles==4) {
    varSuffix.push_back("HI020");
    varSuffix.push_back("HI2040");
    varSuffix.push_back("HI40100");
    varSuffix.push_back("pp");
  }
  else {
    char substr[5];
    for (unsigned int i=0; i<nFiles; ++i) {
      sprintf(substr,"%u",i);
      varSuffix.push_back(substr);
    }
  }

  // *** Read Data files
  TFile *fInData[nFiles];
  RooDataSet *data[nFiles];
  for (unsigned int i=0; i<nFiles; i++) {
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
  RooDataSet *redData[nFiles];
  char reduceDS[300];
  if (cutNonPrompt) {
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

  cout << "reduceDS for PbPb data: " << reduceDS << endl;
  for (unsigned int i=0; i<nFiles-1; ++i) {
    redData[i] = (RooDataSet*)data[i]->reduce(reduceDS);
    ws->import(*redData[i]);
  }

  if (cutNonPrompt) {
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
  else
    sprintf(reduceDS,"Jpsi_Pt>%.2f && Jpsi_Pt<%.2f && abs(Jpsi_Y)>%.2f && abs(Jpsi_Y)<%.2f",pmin,pmax,ymin,ymax);

  cout << "reduceDS for pp data: " << reduceDS << endl;

  redData[nFiles-1] = (RooDataSet*)data[nFiles-1]->reduce(reduceDS);
  ws->import(*redData[nFiles-1]);

  setWSRange(ws);

  RooCategory* sample = new RooCategory("sample","sample") ;
  string index;
  if (nFiles==2) {
    sample->defineType("HI");
    sample->defineType("pp");
    index = "sample[HI,pp]";
  }
  else if (nFiles==4) {
    sample->defineType("HI020");
    sample->defineType("HI2040");
    sample->defineType("HI40100");
    sample->defineType("pp");
    index = "sample[HI020,HI2040,HI40100,pp]";
  }
  else {
    index = "sample[";
    char substr[5];
    for (unsigned int i=0;i<nFiles;++i) {
      sprintf(substr,"%u",i);
      sample->defineType(substr);
      if (i>0)
	index += ",";

      index += substr;
    }
    index += "]";
  }
  ws->factory(index.c_str());  //Index for simultaneous fit

  RooDataSet *redDataSim;
  if (nFiles==2) {
    redDataSim = new RooDataSet("redDataSim","redDataSim",RooArgSet(*(ws->var("Jpsi_Mass"))),Index(*(ws->cat("sample"))),Import("HI",*redData[0]),Import("pp",*redData[nFiles-1]));
  }
  else if (nFiles==4) {
    redDataSim = new RooDataSet("redDataSim","redDataSim",RooArgSet(*(ws->var("Jpsi_Mass"))),Index(*(ws->cat("sample"))),Import("HI020",*redData[0]),Import("HI2040",*redData[1]),Import("HI40100",*redData[2]),Import("pp",*redData[3]));
  }
  else {
    redDataSim = NULL;
    cout << nFiles << " files not yet supported." << endl;
    return 1;
  }

  ws->import(*redDataSim);

  // Draw data
  string titlestr;

  // Binning for invariant mass distribution
  RooBinning rbm(2.2,4.2);
  int nbins = 100;
  // RooBinning rbm(2.0,4.5);
  // int nbins = 150;
  bool highStats=true;
  if (highStats)
    rbm.addUniform(nbins,2.2,4.2);
  //    rbm.addUniform(nbins,2.0,4.5);
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


  RooDataHist *binData[nFiles];
  for (unsigned int i=0; i<nFiles; i++) {
    char name[100] = {0};
    sprintf(name,"binData%d",i+1);
    binData[i]= new RooDataHist(name,name,RooArgSet( *(ws->var("Jpsi_Mass")) ), *redData[i]);
    cout << "DATA" << i+1 << " :: N events to fit: " << binData[i]->sumEntries() << endl;
  }
  
  // *** Define PDFs with parameters (mass and ctau)
  // Just so RooFit does not crash on Ubuntu
  RooRealVar aa("aa","aa",0.5,-1,1);
  RooRealVar ab("ab","ab",-0.5,-1,1);
  RooChebychev tmpPol("tmpPol","tmpPol",*(ws->var("Jpsi_Mass")),RooArgSet(aa,ab));

  // J/psi mass parameterization
  defineMassBkg(ws);
  defineMassSig(ws);
  defineMassBkgHI(ws);
  defineMassSigHI(ws);

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
  
  RooRealVar *NJpsi[nFiles];
  RooRealVar *NBkg[nFiles];
  RooRealVar *fracP_pp;
  RooFormulaVar *NPsiP[nFiles];
  RooRealVar *doubleRatio[nFiles-1];
  RooFormulaVar *fracP_HI[nFiles-1];
  RooFormulaVar *NPsiP_mix[nFiles-1];

  for (int i=nFiles-1; i>=0; --i) {
    NJpsi[i]  = new RooRealVar(("NJpsi_"+varSuffix.at(i)).c_str(),("J/psi yield in "+varSuffix.at(i)).c_str(),0.5*binData[i]->sumEntries(),0.0,2.0*binData[i]->sumEntries()); ws->import(*NJpsi[i]);
    
    NBkg[i]  = new RooRealVar(("NBkg_"+varSuffix.at(i)).c_str(),("Brackground yield in "+varSuffix.at(i)).c_str(), 0.5*binData[i]->sumEntries(),0.0,2.0*binData[i]->sumEntries()); ws->import(*NBkg[i]);

    if (i == (int)nFiles-1) {
      fracP_pp = new RooRealVar(("fracP_"+varSuffix.at(i)).c_str(),("psi(2S) fraction in "+varSuffix.at(i)).c_str(),0.01);
      fracP_pp->setConstant(false); ws->import(*fracP_pp);
      NPsiP[i] = new RooFormulaVar(("NPsiP_"+varSuffix.at(i)).c_str(), "@0*@1", RooArgList(*(ws->var(("NJpsi_"+varSuffix.at(i)).c_str())),*(ws->var(("fracP_"+varSuffix.at(i)).c_str()))));  ws->import(*NPsiP[i]);
    }
    else {
      doubleRatio[i] = new RooRealVar(("doubleRatio_"+varSuffix.at(i)).c_str(),("psi(2S) double ratio in "+varSuffix.at(i)).c_str(),1.0);
      doubleRatio[i]->setConstant(false); ws->import(*doubleRatio[i]);

      fracP_HI[i] = new RooFormulaVar(("fracP_"+varSuffix.at(i)).c_str(), "@0*@1", RooArgList(*(ws->var(("doubleRatio_"+varSuffix.at(i)).c_str())),*(ws->var("fracP_pp")))); ws->import(*fracP_HI[i]);

      NPsiP[i] = new RooFormulaVar(("NPsiP_"+varSuffix.at(i)).c_str(), "@0*@1", RooArgList(*(ws->var(("NJpsi_"+varSuffix.at(i)).c_str())),*(ws->function(("fracP_"+varSuffix.at(i)).c_str()))));  ws->import(*NPsiP[i]);

      NPsiP_mix[i]  = new RooFormulaVar(("NPsiP_mix_"+varSuffix.at(i)).c_str(),"@0*@1",RooArgList(*(ws->var(("NJpsi_"+varSuffix.at(nFiles-1)).c_str())),*(ws->function(("fracP_"+varSuffix.at(i)).c_str()))));   ws->import(*NPsiP_mix[i]);
    }

    
    sprintf(funct,"SUM::sigMassPDF_%s(NJpsi_%s*%s,NPsiP_%s*%s,NBkg_%s*%s)",varSuffix.at(i).c_str(),varSuffix.at(i).c_str(),mJpsiFunct.c_str(),varSuffix.at(i).c_str(),mPsiPFunct.c_str(),varSuffix.at(i).c_str(),mBkgFunct[i].c_str());

    cout << funct << endl;
    ws->factory(funct);
    
    if (i < (int)nFiles-1) {
      sprintf(funct,"SUM::sigMassPDF_mix_%s(NJpsi_pp*%s,NPsiP_mix_%s*%s,NBkg_pp*%s)",varSuffix.at(i).c_str(),mJpsiFunct.c_str(),varSuffix.at(i).c_str(),mPsiPFunct.c_str(),mBkgFunct[nFiles-1].c_str());
      ws->factory(funct);
      cout << funct << endl;
    }
  }
  
  if (nFiles==2) {
    ws->factory("SIMUL::sigMassPDFSim(sample,HI=sigMassPDF_HI,pp=sigMassPDF_pp)");
  }
  else if (nFiles==4) {
    ws->factory("SIMUL::sigMassPDFSim(sample,HI020=sigMassPDF_HI020,HI2040=sigMassPDF_HI2040,HI40100=sigMassPDF_HI40100,pp=sigMassPDF_pp)");
  }
  else {
    index = "SIMUL::sigMassPDFSim(sample";
    char substr[5];
    for (unsigned int i=0;i<nFiles;++i) {
      sprintf(substr,",%u=sigMassPDF_%u",i,i);
      index += substr;
    }
    index += ")";
    ws->factory(index.c_str());
  }

  string str("CBG");
  string dirPre2 = dirPre;
  size_t found = dirPre2.find(str);
  if (found!=string::npos)
    dirPre2.replace(found,str.length(),"CB");

  // fix CB parameters to MC
  if (fixCBtoMC) {
    if (prange=="6.5-30.0" && yrange=="0.0-2.4") {
      ws->var("alpha_HI")->setVal(1.489);
      ws->var("enne_HI")->setVal(1.730);
      ws->var("wideFactor_HI")->setVal(1.901);
      ws->var("alpha")->setVal(1.693);
      ws->var("enne")->setVal(1.709);
    }
    else if (prange=="6.5-30.0" && yrange=="0.0-1.6") {
      ws->var("alpha_HI")->setVal(1.516);
      ws->var("wideFactor_HI")->setVal(1.699);
      ws->var("enne_HI")->setVal(1.683);
      ws->var("alpha")->setVal(1.710);
      ws->var("enne")->setVal(1.646);
    }
    else if (prange=="3.0-30.0" && yrange=="1.6-2.4") {
      ws->var("alpha_HI")->setVal(1.880);
      ws->var("enne_HI")->setVal(1.430);
      ws->var("wideFactor_HI")->setVal(1.503);
      ws->var("alpha")->setVal(2.079);
      ws->var("enne")->setVal(1.390);
    }
    else if (prange=="3.0-6.5" && yrange=="1.6-2.4") {
      ws->var("alpha_HI")->setVal(1.860);
      ws->var("enne_HI")->setVal(1.330);
      ws->var("wideFactor_HI")->setVal(1.540);
      ws->var("alpha")->setVal(2.000);
      ws->var("enne")->setVal(1.350);
    }
    else if (prange=="6.5-30.0" && yrange=="1.6-2.4") {
      ws->var("alpha_HI")->setVal(2.187);
      ws->var("enne_HI")->setVal(1.295);
      ws->var("wideFactor_HI")->setVal(1.830);
      ws->var("alpha")->setVal(2.161);
      ws->var("enne")->setVal(1.370);
    }
    else {
      ws->var("alpha_HI")->setVal(2.0);
      ws->var("enne_HI")->setVal(1.4);
      ws->var("wideFactor_HI")->setVal(2.0);
      ws->var("alpha")->setVal(2.0);
      ws->var("enne")->setVal(1.4);
    }
    
    if (fixAlpha) {
      ws->var("alpha")->setConstant(kTRUE);
      ws->var("alpha_HI")->setConstant(kTRUE);
    }
    else {
      ws->var("alpha")->setConstant(kFALSE);
      ws->var("alpha_HI")->setConstant(kFALSE);
    }
    if (fixN) {
      ws->var("enne")->setConstant(kTRUE);
      ws->var("enne_HI")->setConstant(kTRUE);
    }
    else {
      ws->var("alpha")->setConstant(kFALSE);
      ws->var("alpha_HI")->setConstant(kFALSE);
    }
    if (fixGwidth)
      ws->var("wideFactor_HI")->setConstant(kTRUE);
    else
      ws->var("wideFactor_HI")->setConstant(kFALSE);
  }

  // 20140128: seed with CB fit parameters
  if ( found!=string::npos ) {
    string inputFNcb;
    if ( !fixAlpha || !fixN || !fixGwidth) {// read for free alpha, n, or wideFactor from the default CBG fit
      inputFNcb =  dirPre + "_PbPb" + mBkgFunct[0];
      for (unsigned int i=1; i<nFiles-1; ++i) {
	inputFNcb += "_" + mBkgFunct[i];
      }
      inputFNcb += "_pp" + mBkgFunct[nFiles-1] + "_rap" + yrange_str + "_pT" + prange_str + "_cent" + crange + "_allVars.txt";
    }
    else { // read results from CB fit
      inputFNcb =  dirPre2 + "_PbPb" + mBkgFunct[0];
      for (unsigned int i=1; i<nFiles-1; ++i) {
	inputFNcb += "_" + mBkgFunct[i];
      }
      inputFNcb += "_pp" + mBkgFunct[nFiles-1] + "_rap" + yrange_str + "_pT" + prange_str + "_cent" + crange + "_allVars.txt";
    }

    RooArgSet *set = ws->pdf("sigMassPDFSim")->getParameters(*(ws->var("Jpsi_Mass")));
    set->readFromFile(inputFNcb.c_str());
    cout << "Import variable values from: " << inputFNcb << endl;
    ws->import(*set);
  }

  fitM = ws->pdf("sigMassPDFSim")->fitTo(*redDataSim,Extended(1),Hesse(1),Minos(0),Save(1),SumW2Error(0),NumCPU(8),PrintEvalErrors(0),Verbose(0));
  fitM->printMultiline(cout,1,1,"");

  int edmStatus = 1;
  if (fitM->edm()<1.0e-03)
    edmStatus = 0;

  int covStatus = 1;
  if (fitM->covQual()==3)
    covStatus = 0;

  int fitStatus = fitM->status();
  if (found!=string::npos)
    cout << "CBG_PbPb";
  else
    cout << "CB_PbPb";

  for (unsigned int i=0; i<nFiles-1; ++i) {
    cout << "_" << mBkgFunct[i];
  }
  cout << "_pp" << mBkgFunct[nFiles-1] << "_rap" << yrange_str << "_pT" << prange_str << "_cent" << crange << fix_str << endl;

  cout << "---FIT result summary: " << edmStatus << covStatus << fitStatus;
  for (unsigned int j=0; j<fitM->numStatusHistory(); ++j) {
    cout << fitM->statusCodeHistory(j);
  }
  cout << " ---" << endl;

  int nFitParam = fitM->floatParsFinal().getSize();
  resultF.cd();
  fitM->Write();
  cout << "fitM->Write() into: " << resultFN << endl;
  resultF.Close();

  // *** Draw mass plot
  RooPlot *mframe[nFiles];
  RooPlot *mframepull[nFiles];
  RooHist *hpull[nFiles]; 
  TH1 *hdata[nFiles];
  TH1F *hpull_proj[nFiles];

  // *** TFile for saving plots
  string plotFN;
  plotFN = dirPre + "_PbPb" + mBkgFunct[0];
  for (unsigned int i=1; i<nFiles-1; ++i) {
    plotFN += "_" + mBkgFunct[i];
  }
  plotFN += "_pp" + mBkgFunct[nFiles-1] + "_rap" + yrange_str + "_pT" + prange_str + "_cent" + crange + fix_str + "plots.root";
  TFile plotF(plotFN.c_str(),"RECREATE");
  //    plotF.cd();
  TDirectory *dir[nFiles];

  for (unsigned int i=0; i<nFiles; ++i) {
    if (i<nFiles-1)
      isPbPb = true;
    else
      isPbPb = false;

    string crange_tmp = crange;
    if (nFiles==2) {
      switch (i) {
      case 0:
	crange = crange_tmp;
	break;
      case 1:
	crange = "pp";
	break;
      default:
	crange = "0-100";
	break;
      }
    }
    else if (nFiles==4) {
      switch (i) {
      case 0:
	crange = "0-20";
	break;
      case 1:
	crange = "20-40";
	break;
      case 2:
	crange = "40-100";
	break;
      case 3:
	crange = "pp";
	break;
      default:
	crange = "0-100";
	break;
      }
    }
    else
      crange = crange_tmp;

    if (isPbPb)
      getOptRange(crange,&cmin,&cmax);

    mframe[i] = ws->var("Jpsi_Mass")->frame();
    mframe[i]->SetName(("frame_"+varSuffix.at(i)).c_str());
    mframe[i]->SetTitle(("A RooPlot of \"J/#psi mass\" in "+varSuffix.at(i)).c_str());

    redData[i]->plotOn(mframe[i],DataError(RooAbsData::SumW2),XErrorSize(0),MarkerSize(0.8),Binning(rbm));
    titlestr = "2D fit for" + partTit + "muons (mass projection), p_{T} = " + prange + " GeV/c and |y| = " + yrange;
    mframe[i]->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
    mframe[i]->GetXaxis()->CenterTitle(1);
    double max = mframe[i]->GetMaximum() * 1.3;
    double min = 0.0;
    mframe[i]->SetMaximum(max);

    string name = "NBkg_"+varSuffix.at(i);
    min = ws->var(name.c_str())->getVal()/(double(nbins)) * 0.7;

    ws->pdf(("sigMassPDF_"+varSuffix.at(i)).c_str())->plotOn(mframe[i],VisualizeError(*fitM,1,kFALSE),FillColor(kGray),LineWidth(2),Normalization(redData[i]->sumEntries(),RooAbsReal::NumEvent));
    ws->pdf(("sigMassPDF_"+varSuffix.at(i)).c_str())->plotOn(mframe[i],LineColor(kBlack),LineWidth(2),Normalization(redData[i]->sumEntries(),RooAbsReal::NumEvent));

    hpull[i] = mframe[i]->pullHist(0,0,true);
    hpull[i]->SetName(("hpullhist_"+varSuffix.at(i)).c_str());
    hpull[i]->SetTitle(("hpullhist_"+varSuffix.at(i)).c_str());
    
    ws->pdf(("sigMassPDF_"+varSuffix.at(i)).c_str())->plotOn(mframe[i],Components(mBkgFunct[i].c_str()),VisualizeError(*fitM,1,kFALSE),FillColor(kCyan),Normalization(redData[i]->sumEntries(),RooAbsReal::NumEvent));
    ws->pdf(("sigMassPDF_"+varSuffix.at(i)).c_str())->plotOn(mframe[i],Components(mBkgFunct[i].c_str()),LineColor(kBlue),LineStyle(kDashed),LineWidth(2),Normalization(redData[i]->sumEntries(),RooAbsReal::NumEvent));

    if (!isPaper && found!=string::npos) { 
      string signal;
      if (isPbPb && !shareShape)
       	signal = ",signalCB1_HI,signalCB1P_HI";
      else
	signal = ",signalCB1,signalCB1P";
      ws->pdf(("sigMassPDF_"+varSuffix.at(i)).c_str())->plotOn(mframe[i],VisualizeError(*fitM,1,kFALSE),FillColor(kRed-9),Components((mBkgFunct[i]+signal).c_str()));
      ws->pdf(("sigMassPDF_"+varSuffix.at(i)).c_str())->plotOn(mframe[i],Components((mBkgFunct[i]+signal).c_str()),LineColor(kRed),LineStyle(kDashed),LineWidth(2));
      if (isPbPb && !shareShape)
       	signal = ",signalG2_HI,signalG2P_HI";
      else
	signal = ",signalG2,signalG2P";
      ws->pdf(("sigMassPDF_"+varSuffix.at(i)).c_str())->plotOn(mframe[i],VisualizeError(*fitM,1,kFALSE),FillColor(kGreen-9),Components((mBkgFunct[i]+signal).c_str()));
      ws->pdf(("sigMassPDF_"+varSuffix.at(i)).c_str())->plotOn(mframe[i],Components((mBkgFunct[i]+signal).c_str()),LineColor(kGreen+2),LineStyle(kDashed),LineWidth(2));
    }

    if (i==nFiles-1) {
      ws->pdf(("sigMassPDF_mix_"+varSuffix.at(0)).c_str())->plotOn(mframe[i],VisualizeError(*fitM,1,kFALSE),FillColor(kOrange-9));
      redData[i]->plotOn(mframe[i],DataError(RooAbsData::SumW2),XErrorSize(0),MarkerStyle(20),MarkerSize(0.8),Binning(rbm));
      ws->pdf(("sigMassPDF_mix_"+varSuffix.at(0)).c_str())->plotOn(mframe[i],LineColor(kOrange+2),LineWidth(2));
    }


    hdata[i] = redData[i]->createHistogram(("hdata"+varSuffix.at(i)).c_str(),*ws->var("Jpsi_Mass"),Binning(rbm));


    // *** Calculate chi2/nDof for mass fitting
    unsigned int nBins = hdata[i]->GetNbinsX();
    double Chi2 = 0;
    int nFullBinsPull = 0;
    double *ypull;
    ypull = hpull[i]->GetY();
  
    // Make a histogram with the pull distribution
    hpull_proj[i] = new TH1F("hpull_proj","hpull_proj;Pull;Events",100,-10,10);
    hpull_proj[i]->SetName(("hpull_proj_"+varSuffix.at(i)).c_str());
    hpull_proj[i]->SetTitle(("hpull_proj_"+varSuffix.at(i)).c_str());

    for (unsigned int j=0; j < nBins; ++j) {
      hpull_proj[i]->Fill(ypull[j]);
      if (hdata[i]->GetBinContent(j+1) == 0) continue;
      nFullBinsPull++;
      Chi2 = Chi2 + pow(ypull[j],2);
    }

    double UnNormChi2 = Chi2;
    cout << "nFullBinsPull = " << nFullBinsPull << endl;
    // only consider the free parameters relevant for particular dataset
    nFitParam = ws->pdf(("sigMassPDF_"+varSuffix.at(i)).c_str())->getParameters(redData[i])->selectByAttrib("Constant",kFALSE)->getSize();
    if (i!=nFiles-1) nFitParam--; // do not count fracP_pp as free parameter in PbPb fits
    
    int Dof = nFullBinsPull - nFitParam;
    if (Dof!=0)
      Chi2 /= Dof;

    double theNLL=0;
    theNLL = fitM->minNll();
    
    mframepull[i] =  ws->var("Jpsi_Mass")->frame(Title(("Pull of "+varSuffix.at(i)).c_str()));
    mframepull[i]->SetName(("pull_frame_"+varSuffix.at(i)).c_str());
    mframepull[i]->GetYaxis()->SetTitle("Pull");
    mframepull[i]->SetLabelSize(0.08,"XYZ");
    mframepull[i]->SetTitleSize(0.1,"XYZ");
    mframepull[i]->SetTitleOffset(0.55,"Y");
    mframepull[i]->addPlotable(hpull[i],"PX") ;

    if ( fabs(mframepull[i]->GetMaximum()) > fabs(mframepull[i]->GetMinimum()) )
      mframepull[i]->SetMinimum(-(mframepull[i]->GetMaximum())); 
    else
      mframepull[i]->SetMaximum(-(mframepull[i]->GetMinimum())); 
    mframepull[i]->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
    mframepull[i]->GetXaxis()->CenterTitle(1);

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

    mframepull[i]->addObject(f0,"same");
    mframepull[i]->addObject(f1,"same");
    mframepull[i]->addObject(f2,"same");

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
      pad2->cd(); mframepull[i]->Draw();
      pad1->cd();
    }
    

    if (logScale) {
      if (isPaper)
	c00.SetLogy(1);
      else
	pad1->SetLogy(1);

      mframe[i]->SetMinimum(0.2*min);
      if ( (yrange == "0.0-1.6" || yrange == "0.0-2.4") && isPbPb) {
	if (crange == "20-40") {
	  mframe[i]->SetMaximum(3.5*max);
	  mframe[i]->SetMinimum(0.5*min);
	}
	else if ( crange == "40-100" ) {
	  mframe[i]->SetMaximum(60*max);
	  mframe[i]->SetMinimum(0.1*min);
	}
	else{
	  mframe[i]->SetMaximum(2.5*max);
	  mframe[i]->SetMinimum(0.8*min);
	}
      }
      else if (yrange == "1.6-2.4" && isPbPb) {
	if ( crange == "40-100" ) {
	  mframe[i]->SetMaximum(10.0*max);
	  mframe[i]->SetMinimum(0.5*min);
	}
	else if ( crange == "20-40" ) {
	  mframe[i]->SetMaximum(2.0*max);
	  mframe[i]->SetMinimum(0.9*min);
	}
	else {
	  mframe[i]->SetMaximum(1.5*max);
	  mframe[i]->SetMinimum(0.95*min);
	}
      }
      else {
	mframe[i]->SetMaximum(20*max);
	mframe[i]->SetMinimum(0.5*min);
      }
    }

    if (!highStats) {
      max = mframe[i]->GetMaximum();
      min = mframe[i]->GetMinimum();
      mframe[i]->SetMaximum(2.0*max);
      mframe[i]->SetMinimum(2.0*min);
    }

    lCMS->SetTextSize(0.05);
    if (isPbPb)
      lCMS->SetText(0.17,0.90,"CMS PbPb #sqrt{s_{NN}} = 2.76 TeV");
    else
      lCMS->SetText(0.17,0.90,"CMS pp #sqrt{s} = 2.76 TeV");
    mframe[i]->addObject(lCMS,"");

    lLumi->SetTextSize(0.035);
    if (isPbPb)
      lLumi->SetText(0.17,0.83,"L_{int} = 150 #mub^{-1}");
    else
      lLumi->SetText(0.17,0.83,"L_{int} = 5.3 pb^{-1}");
    mframe[i]->addObject(lLumi,"");

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
    mframe[i]->addObject(lRap,"");

    lPt->SetTextSize(0.035);
    if (pmin==6.5)
      sprintf(reduceDS,"%.1f < p_{T} < %.0f GeV/c",pmin,pmax);
    else if (pmax==6.5)
      sprintf(reduceDS,"%.0f < p_{T} < %.1f GeV/c",pmin,pmax);
    else
      sprintf(reduceDS,"%.0f < p_{T} < %.0f GeV/c",pmin,pmax);
  
    lPt->SetText(0.17,0.73,reduceDS);
    mframe[i]->addObject(lPt,"");

    lNLL->SetTextSize(0.035);
    if (!isPaper) {
      sprintf(reduceDS,"Min. NLL = %0.2f",theNLL);
      lNLL->SetText(0.17,0.68,reduceDS);
      mframe[i]->addObject(lNLL,"");
    }

    if (!isPaper) {
      sprintf(reduceDS,"#chi^{2}/dof = %0.1f/%d",UnNormChi2,Dof);
      lChi->SetText(0.62,0.90,reduceDS);
      mframe[i]->addObject(lChi,"");
      sprintf(reduceDS,"p-value = %0.4f",TMath::Prob(UnNormChi2,Dof));
      lPval->SetText(0.62,0.85,reduceDS);
      mframe[i]->addObject(lPval,"");
    }

    name = "NJpsi_"+varSuffix.at(i);
    
    if (ws->var(name.c_str())->hasAsymError() && abs(-1.0*ws->var(name.c_str())->getErrorLo()/ws->var(name.c_str())->getErrorHi() - 1)>0.1)
      sprintf(reduceDS,"N_{J/#psi} = %0.0f^{+%0.0f}_{%0.0f}",ws->var(name.c_str())->getVal(),ws->var(name.c_str())->getErrorHi(),ws->var(name.c_str())->getErrorLo());
    else
      sprintf(reduceDS,"N_{J/#psi} = %0.0f #pm %0.0f",ws->var(name.c_str())->getVal(),ws->var(name.c_str())->getError());
    lNJpsi->SetText(0.62,0.80,reduceDS);
    mframe[i]->addObject(lNJpsi,"");

    if (isPbPb) {
      name = "doubleRatio_"+varSuffix.at(i);
      
      double DoubleRatio = ws->var(name.c_str())->getVal();
      double fracP_pp = ws->var("fracP_pp")->getVal();
      double errDoubleRatio = ws->var(name.c_str())->getError();
      double errFracP_pp = ws->var("fracP_pp")->getError();
      double corr = fitM->correlation( *ws->var(name.c_str()) , *ws->var("fracP_pp") );
      std::cout << "Correlation between double ratio and fracP_pp: " << corr << std::endl;
      double fracP_HI = ws->function( ("fracP_"+varSuffix.at(i)).c_str())->getVal();
      double errFracP_HI = sqrt( pow(errDoubleRatio/DoubleRatio,2) + pow(errFracP_pp/fracP_pp,2) + 2*errDoubleRatio*errFracP_pp*corr/fracP_HI )*fracP_HI;

      sprintf(reduceDS,"R_{#psi(2S)} = %0.3f #pm %0.3f",fracP_HI,errFracP_HI);
    }
    else {
      name = "fracP_pp";
      
      if (ws->var(name.c_str())->hasAsymError() && abs(-1.0*ws->var(name.c_str())->getErrorLo()/ws->var(name.c_str())->getErrorHi() - 1)>0.1)
	sprintf(reduceDS,"R_{#psi(2S)} = %0.3f^{+%0.3f}_{%0.3f}",ws->var(name.c_str())->getVal(),ws->var(name.c_str())->getErrorHi(),ws->var(name.c_str())->getErrorLo());
      else
	sprintf(reduceDS,"R_{#psi(2S)} = %0.3f #pm %0.3f",ws->var(name.c_str())->getVal(),ws->var(name.c_str())->getError());
    }
    lRpsi->SetText(0.62,0.75,reduceDS);
    mframe[i]->addObject(lRpsi,"");

    if (!isPbPb) {
      double NJpsi = ws->var("NJpsi_pp")->getVal();
      double fracP = ws->var("fracP_pp")->getVal();
      double errNJpsi = ws->var("NJpsi_pp")->getError();
      double errFracP = ws->var("fracP_pp")->getError();
      double corr = fitM->correlation( *ws->var("NJpsi_pp") , *ws->var("fracP_pp") );
      std::cout << "Correlation between NJpsi and fracP_pp: " << corr << std::endl;
      double Npsi2S = ws->function("NPsiP_pp")->getVal();
      double errNpsi2S = sqrt( pow(errNJpsi/NJpsi,2) + pow(errFracP/fracP,2) + 2*errNJpsi*errFracP*corr/Npsi2S )*Npsi2S;
      // if (ws->var("fracP")->hasAsymError() && abs(-1.0*ws->var("fracP")->getErrorLo()/ws->var("fracP")->getErrorHi() - 1)>0.1)
      //   sprintf(reduceDS,"N_{#psi(2S)} = %0.1f^{+%0.1f}_{%0.1f}",ws->function("NPsiP")->getVal(),ws->var("fracP")->getErrorHi()*ws->var("NJpsi")->getVal(),ws->var("fracP")->getErrorLo()*ws->var("NJpsi")->getVal());
      // else
      sprintf(reduceDS,"N_{#psi(2S)} = %0.1f #pm %0.1f",Npsi2S,errNpsi2S);
    }
    else {
      sprintf(reduceDS,"#chi_{#psi(2S)} = %0.3f #pm %0.3f",ws->var(("doubleRatio_"+varSuffix.at(i)).c_str())->getVal(),ws->var(("doubleRatio_"+varSuffix.at(i)).c_str())->getError());
    }
    lNpsiP->SetText(0.62,0.70,reduceDS);
    mframe[i]->addObject(lNpsiP,"");

    double coeffGaus;
    double sigmaSig1;
    double sigmaSig2;
    double ErrcoeffGaus;
    double ErrsigmaSig1;
    double ErrsigmaSig2;

    if (isPbPb && !shareShape) {
      coeffGaus = ws->var("coeffGaus_HI")->getVal();
      sigmaSig1 = ws->var("sigmaSig1_HI")->getVal();
      sigmaSig2 = ws->function("sigmaSig2_HI")->getVal();
      ErrcoeffGaus = ws->var("coeffGaus_HI")->getError();
      ErrsigmaSig1 = ws->var("sigmaSig1_HI")->getError();
      ErrsigmaSig2 = sigmaSig1*ws->var("wideFactor_HI")->getError();
    }
    else {
      coeffGaus = ws->var("coeffGaus")->getVal();
      sigmaSig1 = ws->var("sigmaSig1")->getVal();
      sigmaSig2 = ws->function("sigmaSig2")->getVal();
      ErrcoeffGaus = ws->var("coeffGaus")->getError();
      ErrsigmaSig1 = ws->var("sigmaSig1")->getError();
      ErrsigmaSig2 = sigmaSig1*ws->var("wideFactor")->getError();
    }

    if (!isPaper) {
      if (ErrsigmaSig1<0.0005) 
	sprintf(reduceDS,"#sigma_{CB} = (%0.0f #pm %0.0f) MeV/c^{2}",sigmaSig1*1000.0,1.0);
      else
	sprintf(reduceDS,"#sigma_{CB} = (%0.0f #pm %0.0f) MeV/c^{2}",sigmaSig1*1000.0,ErrsigmaSig1*1000.0);
      lSigCB->SetText(0.62,0.65,reduceDS);

      if  (found!=string::npos) {
	sprintf(reduceDS,"#sigma_{G} = %0.0f MeV/c^{2}",sigmaSig2*1000.0);
	lSigG->SetText(0.62,0.60,reduceDS);
	if (isPbPb && !shareShape)
	  name = "wideFactor_HI";
	else
	  name = "wideFactor";

	if (ws->var(name.c_str())->isConstant())
	  sprintf(reduceDS,"n_{G} = %0.2f (fixed)",ws->var(name.c_str())->getVal());
	else if (ws->var(name.c_str())->hasAsymError() && abs(-1.0*ws->var(name.c_str())->getErrorLo()/ws->var(name.c_str())->getErrorHi() - 1)>0.1)
	  sprintf(reduceDS,"n_{G} = %0.2f^{+%0.2f}_{%0.2f}",ws->var(name.c_str())->getVal(),ws->var(name.c_str())->getErrorHi(),ws->var(name.c_str())->getErrorLo());
	else
	  sprintf(reduceDS,"n_{G} = %0.2f #pm %0.2f",ws->var(name.c_str())->getVal(),ws->var(name.c_str())->getError());
	lNG->SetText(0.62,0.55,reduceDS);
      }
      mframe[i]->addObject(lSigCB,"");
      mframe[i]->addObject(lSigG,"");
      mframe[i]->addObject(lNG,"");
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
	mframe[i]->addObject(lSigma,"");

	if (isPbPb && !shareShape)
	  name = "alpha_HI";
	else
	  name = "alpha";
	if (ws->var(name.c_str())->isConstant())
	  sprintf(reduceDS,"#alpha = %0.2f (fixed)",ws->var(name.c_str())->getVal());
	else if (ws->var(name.c_str())->hasAsymError() && abs(-1.0*ws->var(name.c_str())->getErrorLo()/ws->var(name.c_str())->getErrorHi() - 1)>0.1)
	  sprintf(reduceDS,"#alpha = %0.2f^{+%0.2f}_{%0.2f}",ws->var(name.c_str())->getVal(),ws->var(name.c_str())->getErrorHi(),ws->var(name.c_str())->getErrorLo());
	else
	  sprintf(reduceDS,"#alpha = %0.2f #pm %0.2f",ws->var(name.c_str())->getVal(),ws->var(name.c_str())->getError());
	lAlpha->SetText(0.62,0.45,reduceDS);

	if (isPbPb && !shareShape)
	  name = "enne_HI";
	else
	  name = "enne";
	if (ws->var(name.c_str())->isConstant())
	  sprintf(reduceDS,"n = %0.2f (fixed)",ws->var(name.c_str())->getVal());
	else if (ws->var(name.c_str())->hasAsymError() && abs(-1.0*ws->var(name.c_str())->getErrorLo()/ws->var(name.c_str())->getErrorHi() - 1)>0.1)
	  sprintf(reduceDS,"n = %0.2f^{+%0.2f}_{%0.2f}",ws->var(name.c_str())->getVal(),ws->var(name.c_str())->getErrorHi(),ws->var(name.c_str())->getErrorLo());
	else
	  sprintf(reduceDS,"n = %0.2f #pm %0.2f",ws->var(name.c_str())->getVal(),ws->var(name.c_str())->getError());
	lN->SetText(0.62,0.40,reduceDS);

	if (isPbPb && !shareShape)
	  name = "coeffGaus_HI";
	else
	  name = "coeffGaus";
	if (ws->var(name.c_str())->hasAsymError() && abs(-1.0*ws->var(name.c_str())->getErrorLo()/ws->var(name.c_str())->getErrorHi() - 1)>0.1)
	  sprintf(reduceDS,"f_{G} = %0.3f^{+%0.3f}_{%0.3f}",ws->var(name.c_str())->getVal(),ws->var(name.c_str())->getErrorHi(),ws->var(name.c_str())->getErrorLo());
	else
	  sprintf(reduceDS,"f_{G} = %0.2f #pm %0.2f",ws->var(name.c_str())->getVal(),ws->var(name.c_str())->getError());
	lFG->SetText(0.62,0.35,reduceDS);
	mframe[i]->addObject(lFG,"");
      }
      else {
	if (isPbPb && !shareShape)
	  name = "alpha_HI";
	else
	  name = "alpha";
	if (ws->var(name.c_str())->isConstant())
	  sprintf(reduceDS,"#alpha = %0.2f (fixed)",ws->var(name.c_str())->getVal());
	else
	  sprintf(reduceDS,"#alpha = %0.2f #pm %0.2f",ws->var(name.c_str())->getVal(),ws->var(name.c_str())->getError());
	lAlpha->SetText(0.62,0.60,reduceDS);

	if (isPbPb && !shareShape)
	  name = "enne_HI";
	else
	  name = "enne";
	if (ws->var(name.c_str())->isConstant())
	  sprintf(reduceDS,"n = %0.2f (fixed)",ws->var(name.c_str())->getVal());
	else
	  sprintf(reduceDS,"n = %0.2f #pm %0.2f",ws->var(name.c_str())->getVal(),ws->var(name.c_str())->getError());
	lN->SetText(0.62,0.55,reduceDS);
      }

      mframe[i]->addObject(lAlpha,"");
      mframe[i]->addObject(lN,"");
    }

    TLegend *leg1 = new TLegend(0.63,0.53,0.92,0.68,NULL,"brNDC");

    leg1->SetFillStyle(0); leg1->SetBorderSize(0); leg1->SetShadowColor(0); leg1->SetMargin(0.2);
    leg1->AddEntry(gfake1,"data","p");
    leg1->AddEntry(&hfake21,"total fit","lf");
    leg1->AddEntry(&hfake11,"background","lf");
    if (!isPbPb)
      leg1->AddEntry(&hfake22,"with R_{#psi(2S)}(PbPb)","lf");

    if (isPaper)
      mframe[i]->addObject(leg1,"same");

    mframe[i]->Draw();
    if (!isPaper) {
      pad3->cd();
      pad3->SetFillStyle(0);
      gStyle->SetOptStat("mr");
      hpull_proj[i]->Draw();
      hpull_proj[i]->GetXaxis()->CenterTitle(1);
      hpull_proj[i]->GetXaxis()->SetRangeUser(-5,5);
      hpull_proj[i]->GetXaxis()->SetLabelSize(0.07);
      hpull_proj[i]->GetYaxis()->SetLabelSize(0.07);
      hpull_proj[i]->GetXaxis()->SetTitleSize(0.07);
      hpull_proj[i]->GetYaxis()->SetTitleSize(0.07);
      pad3->Update();
      TPaveStats *st = (TPaveStats*) hpull_proj[i]->FindObject("stats");
      st->SetX1NDC(0.64);
      st->SetX2NDC(0.99);
      st->SetY1NDC(0.75);
      st->SetY2NDC(0.99);
      st->SetTextSize(0.065);
      st->SetTextFont(42);
    }
 
    titlestr = dirPre + "_" + mBkgFunct[i] + "_rap" + yrange_str  + "_pT" + prange_str + "_cent" + crange + fix_str +  "massfit.pdf";
    c00.SaveAs(titlestr.c_str());


    ws->import(*mframe[i]);
    ws->import(*mframepull[i]);
    ws->import(*hpull[i]);
    ws->import(*hpull_proj[i]);

    dir[i] = plotF.mkdir(varSuffix.at(i).c_str());
    dir[i]->cd();

    mframe[i]->Write();
    mframepull[i]->Write();
    hpull[i]->Write();
    hpull_proj[i]->Write();
    plotF.cd();

    crange = crange_tmp; // set back to orignal crange
  }
  plotF.Close();

  // prepare models for CLs calculations
  RooStats::ModelConfig *model[nFiles-1];
  RooStats::ModelConfig *bmodel[nFiles-1];
  RooArgSet *nuisances[nFiles-1];

  for (unsigned int i=0;i<nFiles-1;++i) {
    model[i] = new RooStats::ModelConfig(("model_"+varSuffix.at(i)).c_str(),ws);
    model[i]->SetPdf(*ws->pdf("sigMassPDFSim"));
    model[i]->SetParametersOfInterest(*ws->var(("doubleRatio_"+varSuffix.at(i)).c_str()));
    model[i]->SetObservables(RooArgSet(*ws->var("Jpsi_Mass"),*ws->cat("sample")));

    RooArgSet *pars = ws->pdf("sigMassPDFSim")->getParameters(redData[nFiles-1]);
    pars->add(*ws->pdf("sigMassPDFSim")->getParameters(redData[i]));
    cout << "All Parameters:" << endl;
    pars->Print("v");
    nuisances[i] = new RooArgSet(*pars);
    nuisances[i]->remove(*ws->var(("doubleRatio_"+varSuffix.at(i)).c_str()));
    nuisances[i]->remove(*ws->cat("sample"));
    cout << "Nuisance Parameters:" << endl;
    nuisances[i]->Print("v");
    ws->defineSet(("nuisParameters_"+varSuffix.at(i)).c_str(),*nuisances[i]);
    model[i]->SetNuisanceParameters(*ws->set(("nuisParameters_"+varSuffix.at(i)).c_str()));
    RooRealVar* poi = (RooRealVar*) model[i]->GetParametersOfInterest()->first();
    model[i]->SetSnapshot(*poi);

    ws->import(*model[i]);


    bmodel[i] = (RooStats::ModelConfig*)model[i]->Clone(("bmodel_"+varSuffix.at(i)).c_str());
    bmodel[i]->SetName(("B_only_model_"+varSuffix.at(i)).c_str());
    poi->setVal(1.0);
    bmodel[i]->SetSnapshot(*poi);
    ws->import(*bmodel[i]);
  }
  
  string fname = dirPre + "_PbPb" + mBkgFunct[0];
  for (unsigned int i=1; i<nFiles-1;++i) {
    fname += "_" + mBkgFunct[i];
  }
  fname += "_pp" + mBkgFunct[nFiles-1] + "_rap" + yrange_str  + "_pT" + prange_str + "_cent" + crange + fix_str + "Workspace.root";
  ws->writeToFile(fname.c_str(),kTRUE);

  fname = dirPre + "_PbPb" + mBkgFunct[0];
  for (unsigned int i=1; i<nFiles-1;++i) {
    fname += "_" + mBkgFunct[i];
  }
  fname += "_pp" + mBkgFunct[nFiles-1] + "_rap" + yrange_str  + "_pT" + prange_str + "_cent" + crange + fix_str + "allVars.txt";
  ws->allVars().writeToFile(fname.c_str());

  double NJpsi_fin[2] = {0,0};
  double ErrNJpsi_fin[2] = {0,0};
  // double NPsiP_fin[2] = {0,0};
  // double ErrNPsiP_fin[2] = {0,0};
  double fracP_fin[2] = {0,0};
  double ErrFracP_fin[2] = {0,0};
  double NBkg_fin[2] = {0,0};
  double ErrNBkg_fin[2] = {0,0};
  double DoubleRatio_fin = 0.0;
  double ErrDoubleRatio_fin = 0.0;

  for (unsigned int i=0; i<nFiles; i++) {
    string name;

    name = "NJpsi_"+varSuffix.at(i);
    NJpsi_fin[i]= ws->var(name.c_str())->getVal();
    ErrNJpsi_fin[i] = ws->var(name.c_str())->getError();
    name = "NBkg_"+varSuffix.at(i);
    NBkg_fin[i] = ws->var(name.c_str())->getVal();
    ErrNBkg_fin[i] = ws->var(name.c_str())->getError();

    if (i<nFiles-1) {
      name = "doubleRatio_"+varSuffix.at(i);
      DoubleRatio_fin = ws->var(name.c_str())->getVal();
      ErrDoubleRatio_fin = ws->var(name.c_str())->getError();
      name = "fracP_"+varSuffix.at(i);

      double fracP_pp = ws->var("fracP_pp")->getVal();
      double errFracP_pp = ws->var("fracP_pp")->getError();
      double corr = fitM->correlation( *ws->var(("doubleRatio_"+varSuffix.at(i)).c_str()) , *ws->var("fracP_pp") );
      fracP_fin[i] = ws->function(("fracP_"+varSuffix.at(i)).c_str())->getVal();
      ErrFracP_fin[i] = sqrt( pow(ErrDoubleRatio_fin/DoubleRatio_fin,2) + pow(errFracP_pp/fracP_pp,2) + 2*ErrDoubleRatio_fin*errFracP_pp*corr/fracP_fin[i] )*fracP_fin[i];
    }
    else {
      name = "fracP_pp";
      fracP_fin[i]= ws->var(name.c_str())->getVal();
      ErrFracP_fin[i] = ws->var(name.c_str())->getError();
    }
  }


  // To check values of fit parameters
  cout << endl << "J/psi yields:" << endl;
  for (unsigned int i=0; i<nFiles; i++) {
    cout << "NJpsi"     << i << " :    Fit :"  << NJpsi_fin[i] << " +/- " << ErrNJpsi_fin[i] << endl;
    cout << "R_psi(2S)" << i << " :    Fit: "  << fracP_fin[i] << " +/- " << ErrFracP_fin[i] << endl;
  }
  cout << "Double Ratio: " << DoubleRatio_fin << " +/- " << ErrDoubleRatio_fin << endl;
 
  /*
  titlestr = dirPre + "_PbPb" + mBkgFunct[0] + "_" + mBkgFunct[1] + "_" + mBkgFunct[2] + "_pp" + mBkgFunct[3] + "_rap" + yrange_str + "_pT" + prange_str + "_cent" + crange + fix_str + "fitResult.txt";

  ofstream outputFile(titlestr.c_str());
  if (!outputFile.good()) {cout << "Fail to open result file." << endl; return 1;}
  outputFile.precision(10);
  //  outputFile.setf( std::ios::fixed, std::ios::floatfield)

  for (int i=0; i<nFiles; i++) {
  outputFile
    << "NJpsi " << i+1 << NJpsi_fin[i]                       << " " << ErrNJpsi_fin[i] << "\n"
    << "NPsiP " << i+1 << NJpsi_fin[i]*fracP_fin[i] << "\n"
    << "R_psi(2S) " << i+1 << fracP_fin[i]            << " " << ErrFracP_fin[i] << "\n"
    << "NBkg "  << i+1 << NBkg_fin[i]                          << " " << ErrNBkg_fin[i] << "\n";
  }
  if (!mBkgFunct[0].compare("expFunct_HI")) {
    outputFile
      << "coeffExp_HI "      << ws->var("coeffExp_HI")->getVal()      << " " << ws->var("coeffExp_HI")->getError() << "\n";
  }
  if (!mBkgFunct[0].compare("twoExpFunct_HI")) {
    outputFile
      << "coeffExp_HI "      << ws->var("coeffExp_HI")->getVal()      << " " << ws->var("coeffExp_HI")->getError() << "\n"
      << "coeffExp2_HI "     << ws->var("coeffExp2_HI")->getVal()     << " " << ws->var("coeffExp2_HI")->getError() << "\n"
      << "cutx_HI "          << ws->var("cutx_HI")->getVal()          << " " << ws->var("cutx_HI")->getError() << "\n";
  }
  if (!mBkgFunct[0].compare("gausBkg_HI")) {
    outputFile
      << "meanBkg_HI "      << ws->var("meanBkg_HI")->getVal()      << " " << ws->var("meanBkg_HI")->getError() << "\n"
      << "sigmaBkg_HI "     << ws->var("sigmaBkg_HI")->getVal()     << " " << ws->var("sigmaBkg_HI")->getError() << "\n";
  }
  if (!mBkgFunct[0].compare("pol1_HI")) {
    outputFile
      << "coeffPol1_HI "     << ws->var("coeffPol1_HI")->getVal()     << " " << ws->var("coeffPol1_HI")->getError() << "\n";
  }
  if (!mBkgFunct[0].compare("pol2_HI")) {
    outputFile
      << "coeffPol1_HI "     << ws->var("coeffPol1_HI")->getVal()     << " " << ws->var("coeffPol1_HI")->getError() << "\n"
      << "coeffPol2_HI "     << ws->var("coeffPol2_HI")->getVal()     << " " << ws->var("coeffPol2_HI")->getError() << "\n";
  }
  if (!mBkgFunct[0].compare("pol3_HI")) {
    outputFile
      << "coeffPol1_HI "     << ws->var("coeffPol1_HI")->getVal()     << " " << ws->var("coeffPol1_HI")->getError() << "\n"
      << "coeffPol2_HI "     << ws->var("coeffPol2_HI")->getVal()     << " " << ws->var("coeffPol2_HI")->getError() << "\n"
      << "coeffPol3_HI "     << ws->var("coeffPol3_HI")->getVal()     << " " << ws->var("coeffPol3_HI")->getError() << "\n";
  }
  if (!mBkgFunct[0].compare("pol4_HI")) {
    outputFile
      << "coeffPol1_HI "     << ws->var("coeffPol1_HI")->getVal()     << " " << ws->var("coeffPol1_HI")->getError() << "\n"
      << "coeffPol2_HI "     << ws->var("coeffPol2_HI")->getVal()     << " " << ws->var("coeffPol2_HI")->getError() << "\n"
      << "coeffPol3_HI "     << ws->var("coeffPol3_HI")->getVal()     << " " << ws->var("coeffPol3_HI")->getError() << "\n"
      << "coeffPol4_HI "     << ws->var("coeffPol4_HI")->getVal()     << " " << ws->var("coeffPol4_HI")->getError() << "\n";
  }
  if (!mBkgFunct[0].compare("pol5_HI")) {
    outputFile
      << "coeffPol1_HI "     << ws->var("coeffPol1_HI")->getVal()     << " " << ws->var("coeffPol1_HI")->getError() << "\n"
      << "coeffPol2_HI "     << ws->var("coeffPol2_HI")->getVal()     << " " << ws->var("coeffPol2_HI")->getError() << "\n"
      << "coeffPol3_HI "     << ws->var("coeffPol3_HI")->getVal()     << " " << ws->var("coeffPol3_HI")->getError() << "\n"
      << "coeffPol4_HI "     << ws->var("coeffPol4_HI")->getVal()     << " " << ws->var("coeffPol4_HI")->getError() << "\n"
      << "coeffPol5_HI "     << ws->var("coeffPol5_HI")->getVal()     << " " << ws->var("coeffPol5_HI")->getError() << "\n";
  }
  if (!mBkgFunct[1].compare("expFunct")) {
    outputFile
      << "coeffExp "      << ws->var("coeffExp")->getVal()      << " " << ws->var("coeffExp")->getError() << "\n";
  }
  if (!mBkgFunct[1].compare("twoExpFunct")) {
    outputFile
      << "coeffExp "      << ws->var("coeffExp")->getVal()      << " " << ws->var("coeffExp")->getError() << "\n"
      << "coeffExp2 "     << ws->var("coeffExp2")->getVal()     << " " << ws->var("coeffExp2")->getError() << "\n"
      << "cutx "          << ws->var("cutx")->getVal()          << " " << ws->var("cutx")->getError() << "\n";
  }
  if (!mBkgFunct[1].compare("gausBkg")) {
    outputFile
      << "meanBkg "      << ws->var("meanBkg")->getVal()      << " " << ws->var("meanBkg")->getError() << "\n"
      << "sigmaBkg "     << ws->var("sigmaBkg")->getVal()     << " " << ws->var("sigmaBkg")->getError() << "\n";
  }
  if (!mBkgFunct[1].compare("pol1")) {
    outputFile
      << "coeffPol1 "     << ws->var("coeffPol1")->getVal()     << " " << ws->var("coeffPol1")->getError() << "\n";
  }
  if (!mBkgFunct[1].compare("pol2")) {
    outputFile
      << "coeffPol1 "     << ws->var("coeffPol1")->getVal()     << " " << ws->var("coeffPol1")->getError() << "\n"
      << "coeffPol2 "     << ws->var("coeffPol2")->getVal()     << " " << ws->var("coeffPol2")->getError() << "\n";
  }
  if (!mBkgFunct[1].compare("pol3")) {
    outputFile
      << "coeffPol1 "     << ws->var("coeffPol1")->getVal()     << " " << ws->var("coeffPol1")->getError() << "\n"
      << "coeffPol2 "     << ws->var("coeffPol2")->getVal()     << " " << ws->var("coeffPol2")->getError() << "\n"
      << "coeffPol3 "     << ws->var("coeffPol3")->getVal()     << " " << ws->var("coeffPol3")->getError() << "\n";
  }
  if (!mBkgFunct[1].compare("pol4")) {
    outputFile
      << "coeffPol1 "     << ws->var("coeffPol1")->getVal()     << " " << ws->var("coeffPol1")->getError() << "\n"
      << "coeffPol2 "     << ws->var("coeffPol2")->getVal()     << " " << ws->var("coeffPol2")->getError() << "\n"
      << "coeffPol3 "     << ws->var("coeffPol3")->getVal()     << " " << ws->var("coeffPol3")->getError() << "\n"
      << "coeffPol4 "     << ws->var("coeffPol4")->getVal()     << " " << ws->var("coeffPol4")->getError() << "\n";
  }
  if (!mBkgFunct[1].compare("pol5")) {
    outputFile
      << "coeffPol1 "     << ws->var("coeffPol1")->getVal()     << " " << ws->var("coeffPol1")->getError() << "\n"
      << "coeffPol2 "     << ws->var("coeffPol2")->getVal()     << " " << ws->var("coeffPol2")->getError() << "\n"
      << "coeffPol3 "     << ws->var("coeffPol3")->getVal()     << " " << ws->var("coeffPol3")->getError() << "\n"
      << "coeffPol4 "     << ws->var("coeffPol4")->getVal()     << " " << ws->var("coeffPol4")->getError() << "\n"
      << "coeffPol5 "     << ws->var("coeffPol5")->getVal()     << " " << ws->var("coeffPol5")->getError() << "\n";
  }


  outputFile
    << "coeffGaus_HI "    << ws->var("coeffGaus_HI")->getVal()    << " " << ws->var("coeffGaus_HI")->getError() << "\n"
    << "meanSig1_HI "     << ws->var("meanSig1_HI")->getVal()     << " " << ws->var("meanSig1_HI")->getError() << "\n"
    << "sigmaSig1_HI "    << ws->var("sigmaSig1_HI")->getVal()    << " " << ws->var("sigmaSig1_HI")->getError() << "\n"
    << "wideFactor_HI "   << ws->var("wideFactor_HI")->getVal()   << " " << ws->var("wideFactor_HI")->getError() << "\n"
    << "sigmaSig2_HI "    << ws->function("sigmaSig2_HI")->getVal() << " " << ws->function("sigmaSig2_HI")->getVal()*ws->var("wideFactor_HI")->getError()/ws->var("wideFactor_HI")->getVal() << "\n"
    << "alpha_HI "        << ws->var("alpha_HI")->getVal()        << " " << ws->var("alpha_HI")->getError() << "\n"
    << "enne_HI "         << ws->var("enne_HI")->getVal()         << " " << ws->var("enne_HI")->getError() << "\n"
    << "coeffGaus "    << ws->var("coeffGaus")->getVal()    << " " << ws->var("coeffGaus")->getError() << "\n"
    << "meanSig1 "     << ws->var("meanSig1")->getVal()     << " " << ws->var("meanSig1")->getError() << "\n"
    << "sigmaSig1 "    << ws->var("sigmaSig1")->getVal()    << " " << ws->var("sigmaSig1")->getError() << "\n"
    << "wideFactor "   << ws->var("wideFactor")->getVal()   << " " << ws->var("wideFactor")->getError() << "\n"
    << "sigmaSig2 "    << ws->function("sigmaSig2")->getVal() << " " << ws->function("sigmaSig2")->getVal()*ws->var("wideFactor")->getError()/ws->var("wideFactor")->getVal() << "\n"
    << "alpha "        << ws->var("alpha")->getVal()        << " " << ws->var("alpha")->getError() << "\n"
    << "enne "         << ws->var("enne")->getVal()         << " " << ws->var("enne")->getError() << endl;
  */


  for (unsigned int i=0; i<nFiles; i++) {
    fInData[i]->Close();
  }

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


void defineMassBkgHI(RooWorkspace *ws) {
  // 1st order polynomial
  ws->factory("Chebychev::pol1_HI(Jpsi_Mass,{coeffPol1_HI[-0.8,-1.0,1.0]})");
  // 2nd order polynomial
  ws->factory("Chebychev::pol2_HI(Jpsi_Mass,{coeffPol1_HI, coeffPol2_HI[0.0,-1.0,1.0]})");
  // 3rd order polynomial
  ws->factory("Chebychev::pol3_HI(Jpsi_Mass,{coeffPol1_HI, coeffPol2_HI, coeffPol3_HI[0.0,-1.0,1.0]})");
  // 4th order polynomial
  ws->factory("Chebychev::pol4_HI(Jpsi_Mass,{coeffPol1_HI, coeffPol2_HI, coeffPol3_HI, coeffPol4_HI[0.0,-1.0,1.0]})");
  // 5th order polynomial
  ws->factory("Chebychev::pol5_HI(Jpsi_Mass,{coeffPol1_HI, coeffPol2_HI, coeffPol3_HI, coeffPol4_HI, coeffPol5_HI[0.0,-1.0,1.0]})");
  // 6th order polynomial
  ws->factory("Chebychev::pol6_HI(Jpsi_Mass,{coeffPol1_HI, coeffPol2_HI, coeffPol3_HI, coeffPol4_HI, coeffPol5_HI, coeffPol6_HI[0.0,-1.0,1.0]})");

  // expo
  ws->factory("Exponential::expFunct_HI(Jpsi_Mass,coeffExp_HI[-0.1,-3.0,1.0])");

  // gauss
  ws->factory("Gaussian::gausBkg_HI(Jpsi_Mass,meanBkg_HI[0.0,0.0,10.0],sigmaBkg_HI[1.0,0.5,5.0])");

  // 0-20%
  // 1st order polynomial
  ws->factory("Chebychev::pol1_HI020(Jpsi_Mass,{coeffPol1_HI020[-0.8,-1.0,1.0]})");
  // 2nd order polynomial
  ws->factory("Chebychev::pol2_HI020(Jpsi_Mass,{coeffPol1_HI020, coeffPol2_HI020[0.0,-1.0,1.0]})");
  // 3rd order polynomial
  ws->factory("Chebychev::pol3_HI020(Jpsi_Mass,{coeffPol1_HI020, coeffPol2_HI020, coeffPol3_HI020[0.0,-1.0,1.0]})");
  // 4th order polynomial
  ws->factory("Chebychev::pol4_HI020(Jpsi_Mass,{coeffPol1_HI020, coeffPol2_HI020, coeffPol3_HI020, coeffPol4_HI020[0.0,-1.0,1.0]})");
  // 5th order polynomial
  ws->factory("Chebychev::pol5_HI020(Jpsi_Mass,{coeffPol1_HI020, coeffPol2_HI020, coeffPol3_HI020, coeffPol4_HI020, coeffPol5_HI020[0.0,-1.0,1.0]})");
  // 6th order polynomial
  ws->factory("Chebychev::pol6_HI020(Jpsi_Mass,{coeffPol1_HI020, coeffPol2_HI020, coeffPol3_HI020, coeffPol4_HI020, coeffPol5_HI020, coeffPol6_HI020[0.0,-1.0,1.0]})");

  // expo
  ws->factory("Exponential::expFunct_HI020(Jpsi_Mass,coeffExp_HI020[-0.1,-3.0,1.0])");

  // gauss
  ws->factory("Gaussian::gausBkg_HI020(Jpsi_Mass,meanBkg_HI020[0.0,0.0,10.0],sigmaBkg_HI020[1.0,0.5,5.0])");


  // 20-40%
  // 1st order polynomial
  ws->factory("Chebychev::pol1_HI2040(Jpsi_Mass,{coeffPol1_HI2040[-0.8,-1.0,1.0]})");
  // 2nd order polynomial
  ws->factory("Chebychev::pol2_HI2040(Jpsi_Mass,{coeffPol1_HI2040, coeffPol2_HI2040[0.0,-1.0,1.0]})");
  // 3rd order polynomial
  ws->factory("Chebychev::pol3_HI2040(Jpsi_Mass,{coeffPol1_HI2040, coeffPol2_HI2040, coeffPol3_HI2040[0.0,-1.0,1.0]})");
  // 4th order polynomial
  ws->factory("Chebychev::pol4_HI2040(Jpsi_Mass,{coeffPol1_HI2040, coeffPol2_HI2040, coeffPol3_HI2040, coeffPol4_HI2040[0.0,-1.0,1.0]})");
  // 5th order polynomial
  ws->factory("Chebychev::pol5_HI2040(Jpsi_Mass,{coeffPol1_HI2040, coeffPol2_HI2040, coeffPol3_HI2040, coeffPol4_HI2040, coeffPol5_HI2040[0.0,-1.0,1.0]})");
  // 6th order polynomial
  ws->factory("Chebychev::pol6_HI2040(Jpsi_Mass,{coeffPol1_HI2040, coeffPol2_HI2040, coeffPol3_HI2040, coeffPol4_HI2040, coeffPol5_HI2040, coeffPol6_HI2040[0.0,-1.0,1.0]})");

  // expo
  ws->factory("Exponential::expFunct_HI2040(Jpsi_Mass,coeffExp_HI2040[-0.1,-3.0,1.0])");

  // gauss
  ws->factory("Gaussian::gausBkg_HI2040(Jpsi_Mass,meanBkg_HI2040[0.0,0.0,10.0],sigmaBkg_HI2040[1.0,0.5,5.0])");


  // 40-100%
  // 1st order polynomial
  ws->factory("Chebychev::pol1_HI40100(Jpsi_Mass,{coeffPol1_HI40100[-0.8,-1.0,1.0]})");
  // 2nd order polynomial
  ws->factory("Chebychev::pol2_HI40100(Jpsi_Mass,{coeffPol1_HI40100, coeffPol2_HI40100[0.0,-1.0,1.0]})");
  // 3rd order polynomial
  ws->factory("Chebychev::pol3_HI40100(Jpsi_Mass,{coeffPol1_HI40100, coeffPol2_HI40100, coeffPol3_HI40100[0.0,-1.0,1.0]})");
  // 4th order polynomial
  ws->factory("Chebychev::pol4_HI40100(Jpsi_Mass,{coeffPol1_HI40100, coeffPol2_HI40100, coeffPol3_HI40100, coeffPol4_HI40100[0.0,-1.0,1.0]})");
  // 5th order polynomial
  ws->factory("Chebychev::pol5_HI40100(Jpsi_Mass,{coeffPol1_HI40100, coeffPol2_HI40100, coeffPol3_HI40100, coeffPol4_HI40100, coeffPol5_HI40100[0.0,-1.0,1.0]})");
  // 6th order polynomial
  ws->factory("Chebychev::pol6_HI40100(Jpsi_Mass,{coeffPol1_HI40100, coeffPol2_HI40100, coeffPol3_HI40100, coeffPol4_HI40100, coeffPol5_HI40100, coeffPol6_HI40100[0.0,-1.0,1.0]})");

  // expo
  ws->factory("Exponential::expFunct_HI40100(Jpsi_Mass,coeffExp_HI40100[-0.1,-3.0,1.0])");

  // gauss
  ws->factory("Gaussian::gausBkg_HI40100(Jpsi_Mass,meanBkg_HI40100[0.0,0.0,10.0],sigmaBkg_HI40100[1.0,0.5,5.0])");

  return;
}

void defineMassSigHI(RooWorkspace *ws) {
  // narrow Gauss
  ws->factory("Gaussian::signalG1_HI(Jpsi_Mass,meanSig1_HI[3.0969,3.05,3.15],sigmaSig1_HI[0.03,0.005,0.080])");
  // narrow CB
  ws->factory("CBShape::signalCB1_HI(Jpsi_Mass,meanSig1_HI,sigmaSig1_HI,alpha_HI[1.0,0.0,3.0],enne_HI[5.0,1.0,50.0])");

  RooRealVar wideFactor_HI("wideFactor_HI","wideFactor_HI",2.0,1.0,6.0);ws->import(wideFactor_HI);
  RooFormulaVar sigmaSig2_HI("sigmaSig2_HI","@0*@1",RooArgList(*(ws->var("sigmaSig1_HI")),wideFactor_HI));ws->import(sigmaSig2_HI);
  // wide Gauss
  ws->factory("Gaussian::signalG2_HI(Jpsi_Mass,meanSig1_HI,sigmaSig2_HI)");

  // CB(narrow) + Gauss(wide)
  ws->factory("SUM::sigCB1G2_HI(coeffGaus_HI[0.1,0.0,1.0]*signalG2_HI,signalCB1_HI)");

  // Fix Jpsi/psi' mass ratio
  RooFormulaVar meanSig1P_HI("meanSig1P_HI","@0*1.19025",RooArgList(*(ws->var("meanSig1_HI")))); ws->import(meanSig1P_HI);
  // Fix resolution scale: sigma_MJpsi/MJpsi = sigma_Mpsi'/Mpsi'
  RooFormulaVar sigmaSig1P_HI("sigmaSig1P_HI","@0*1.19025",RooArgList(*(ws->var("sigmaSig1_HI")))); ws->import(sigmaSig1P_HI);
  RooFormulaVar sigmaSig2P_HI("sigmaSig2P_HI","@0*1.19025",RooArgList(*(ws->function("sigmaSig2_HI")))); ws->import(sigmaSig2P_HI);
  // narrow Gauss
  ws->factory("Gaussian::signalG1P_HI(Jpsi_Mass,meanSig1P_HI,sigmaSig1P_HI)");
  // wide Gauss
  ws->factory("Gaussian::signalG2P_HI(Jpsi_Mass,meanSig1P_HI,sigmaSig2P_HI)");
  // narrow CB
  ws->factory("CBShape::signalCB1P_HI(Jpsi_Mass,meanSig1P_HI,sigmaSig1P_HI,alpha_HI,enne_HI)");

  // CB(narrow) + Gauss(wide)
  ws->factory("SUM::sigCB1G2P_HI(coeffGaus_HI*signalG2P_HI,signalCB1P_HI)");

  return;
}
