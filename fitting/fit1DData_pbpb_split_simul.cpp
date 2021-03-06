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

#include "RooAbsReal.h"
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
#include <RooPolynomial.h>

#include <RooStats/ModelConfig.h>

#include "RooMsgService.h"

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
  gStyle->SetLineStyleString(11,"8 12"); // thick dotted
  gStyle->SetLineStyleString(12,"20 12"); // thick dashed

  //  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  vector<string> FileName;
  vector<string> mBkgFunct;
  vector<string> mBkgPFunct;
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
  bool fitCentIntegrated = false;
  bool useSystematics = false;
  bool twoCB=false;
  bool overlay=false;
  bool showInsert=false;

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
      case 'u':
	useSystematics = atoi(argv[i+1]);
	cout << "Set systematic uncertainties as nuisance parameters: " << useSystematics << endl;
	break;
      case 'x':
	fitCentIntegrated = atoi(argv[i+1]);
	cout << "Fit with double ratio in 0-100% centrality instead of 0-20%: " << fitCentIntegrated << endl;
	break;
      case 'z':
	isPaper = atoi(argv[i+1]);
	cout << "paper plots: " << isPaper << endl;
	break;
      }
    }
    }
  }// End check options

  if (mJpsiFunct.find("CB2")!=string::npos)
    twoCB=true;

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

  for (unsigned int i=0; i<nFiles; ++i) {
    mBkgPFunct.push_back(mBkgFunct.at(i) + "P");
  }

  // *** TFile for saving fitting results
  string resultFN;
  resultFN = dirPre + "_PbPb" + mBkgFunct.front() + "_" + mBkgPFunct.front();
  for (unsigned int i=1; i<nFiles-1;++i) {
    resultFN +=  "_" + mBkgFunct.at(i);
    resultFN +=  "_" + mBkgPFunct.at(i);
  }
  resultFN += "_pp" + mBkgFunct.back() + "_" + mBkgPFunct.back() + "_rap" + yrange_str + "_pT" + prange_str + "_cent" + crange + fix_str + "fitResult.root";
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
    fInData[i] = new TFile(FileName.at(i).c_str());
    cout << FileName.at(i) << endl;
    if (fInData[i]->IsZombie()) { cout << "CANNOT open data root file: " << FileName.at(i) << "\n"; return 1; }
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
      sprintf(reduceDS,"Jpsi_Pt>%.2f && Jpsi_Pt<%.2f && abs(Jpsi_Y)>%.2f && abs(Jpsi_Y)<%.2f && Jpsi_Ct<0.04",pmin,pmax,ymin,ymax);
    else if (prange=="6.5-30.0" && yrange=="0.0-1.6")
      sprintf(reduceDS,"Jpsi_Pt>%.2f && Jpsi_Pt<%.2f && abs(Jpsi_Y)>%.2f && abs(Jpsi_Y)<%.2f && Jpsi_Ct<0.04",pmin,pmax,ymin,ymax);
    else if (prange=="3.0-30.0" && yrange=="1.6-2.4")
      sprintf(reduceDS,"Jpsi_Pt>%.2f && Jpsi_Pt<%.2f && abs(Jpsi_Y)>%.2f && abs(Jpsi_Y)<%.2f && Jpsi_Ct<0.08",pmin,pmax,ymin,ymax);
    else if (prange=="3.0-6.5" && yrange=="1.6-2.4")
      sprintf(reduceDS,"Jpsi_Pt>%.2f && Jpsi_Pt<%.2f && abs(Jpsi_Y)>%.2f && abs(Jpsi_Y)<%.2f && Jpsi_Ct<0.09",pmin,pmax,ymin,ymax);
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
      sprintf(reduceDS,"Jpsi_Pt>%.2f && Jpsi_Pt<%.2f && abs(Jpsi_Y)>%.2f && abs(Jpsi_Y)<%.2f && Jpsi_Ct<0.04",pmin,pmax,ymin,ymax);
    else if (prange=="6.5-30.0" && yrange=="0.0-1.6")
      sprintf(reduceDS,"Jpsi_Pt>%.2f && Jpsi_Pt<%.2f && abs(Jpsi_Y)>%.2f && abs(Jpsi_Y)<%.2f && Jpsi_Ct<0.04",pmin,pmax,ymin,ymax);
    else if (prange=="3.0-30.0" && yrange=="1.6-2.4")
      sprintf(reduceDS,"Jpsi_Pt>%.2f && Jpsi_Pt<%.2f && abs(Jpsi_Y)>%.2f && abs(Jpsi_Y)<%.2f && Jpsi_Ct<0.09",pmin,pmax,ymin,ymax);
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
  RooBinning *rbm[nFiles];
  int nbins = 100;
  bool highStats=false;
  for (unsigned int i=0;i<nFiles;++i) {
    rbm[i] = new RooBinning(2.2,4.2);

    if (highStats)
      rbm[i]->addUniform(nbins,2.2,4.2);
    else if (logScale && false) {
      switch (i) {
      case 0: // 0-20%
	rbm[i]->addUniform(22,2.2,2.860);   // 30 MeV
	rbm[i]->addUniform(19,2.860,3.240); // 20 MeV
	rbm[i]->addUniform(12,3.240,3.600); // 30 MeV
	rbm[i]->addUniform(3,3.600,3.720);  // 40 MeV
	rbm[i]->addUniform(16,3.720,4.2);   // 30 MeV
	break;
      case 1: // 20-40%
	rbm[i]->addUniform(22,2.2,2.860);   // 30 MeV
	rbm[i]->addUniform(19,2.860,3.240); // 20 MeV
	rbm[i]->addUniform(12,3.240,3.600); // 30 MeV
	rbm[i]->addUniform(4,3.600,3.720);  // 30 MeV
	rbm[i]->addUniform(16,3.720,4.2);   // 30 MeV
	break;
      case 2: // 40-100%
	rbm[i]->addUniform(11,2.2,2.860);   // 60 MeV
	rbm[i]->addUniform(19,2.860,3.240); // 20 MeV
	rbm[i]->addUniform(6,3.240,3.600);  // 60 MeV
	rbm[i]->addUniform(3,3.600,3.720);  // 40 MeV
	rbm[i]->addUniform(8,3.720,4.2);    // 60 MeV
	break;
      case 3: // pp
	rbm[i]->addUniform(22,2.2,2.860);   // 30 MeV
	rbm[i]->addUniform(19,2.860,3.240); // 20 MeV
	rbm[i]->addUniform(12,3.240,3.600); // 30 MeV
	rbm[i]->addUniform(6,3.600,3.720);  // 20 MeV
	rbm[i]->addUniform(16,3.720,4.2);   // 30 MeV
	break;
      default:
	rbm[i]->addUniform(nbins,2.2,4.2);
	break;
      }
    }
    else {
      rbm[i]->addUniform(22,2.2,2.860);   // 30 MeV
      rbm[i]->addUniform(19,2.860,3.240); // 20 MeV
      rbm[i]->addUniform(12,3.240,3.600); // 30 MeV
      rbm[i]->addUniform(5,3.600,3.720);  // 24 MeV
      rbm[i]->addUniform(16,3.720,4.2);   // 30 MeV
    }
  }
  
  RooBinning rbmZoom(2.2,4.2);
  rbmZoom.addUniform(22,2.2,2.860); // 30 MeV
  rbmZoom.addUniform(19,2.860,3.240); // 20 MeV
  rbmZoom.addUniform(2,3.240,3.300); // 30 MeV
  rbmZoom.addUniform(5,3.300,3.600); // 60 MeV
  rbmZoom.addUniform(3,3.600,3.720); // 40 MeV
  rbmZoom.addUniform(8,3.720,4.2); // 60 MeV

  ws->var("Jpsi_Mass")->SetTitle("m_{#mu^{+}#mu^{-}}");
  ws->var("Jpsi_Mass")->setBinning(*rbm[nFiles-1]);


  ws->var("Jpsi_Mass")->setRange("jpsi",2.2,3.5);
  ws->var("Jpsi_Mass")->setRange("psip",3.3,4.2);


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
  RooPolynomial tmpPol("tmpPol","tmpPol",*(ws->var("Jpsi_Mass")),RooArgSet(aa,ab));

  // J/psi mass parameterization
  defineMassBkg(ws);
  defineMassSig(ws);
  defineMassBkgHI(ws);
  defineMassSigHI(ws);

  char funct[250];
  char functP[250];
  string partTit, partFile;
  if (isGG == 0) { partTit = "glb-glb"; partFile = "GG"; }
  else if (isGG == 1) { partTit = "glb-trk"; partFile = "GT"; }
  else { partTit = "all"; partFile = "ALL"; }

  // Global TLatex, TH1, TGraph objects for drawing
  TLatex *lPre = new TLatex();
  lPre->SetNDC(); lPre->SetTextAlign(11);
  TLatex *lCMS = new TLatex();
  lCMS->SetNDC(); lCMS->SetTextAlign(11);
  TLatex *lLumi = new TLatex();
  lLumi->SetNDC(); lLumi->SetTextAlign(11);
  TLatex *lRap = new TLatex();
  lRap->SetNDC(); lRap->SetTextAlign(11);
  TLatex *lPt = new TLatex();
  lPt->SetNDC(); lPt->SetTextAlign(11);
  TLatex *lNLL = new TLatex();
  lNLL->SetNDC(); lNLL->SetTextAlign(11);
  TLatex *lChi = new TLatex();
  lChi->SetNDC(); lChi->SetTextAlign(11);
  TLatex *lPval = new TLatex();
  lPval->SetNDC(); lPval->SetTextAlign(11);
  TLatex *lNJpsi = new TLatex();
  lNJpsi->SetNDC(); lNJpsi->SetTextAlign(11);
  TLatex *lRpsi = new TLatex();
  lRpsi->SetNDC(); lRpsi->SetTextAlign(11);
  TLatex *lNPsiP = new TLatex();
  lNPsiP->SetNDC(); lNPsiP->SetTextAlign(11);
  TLatex *lXpsi = new TLatex();
  lXpsi->SetNDC(); lXpsi->SetTextAlign(11);
  TLatex *lSigCB = new TLatex();
  lSigCB->SetNDC(); lSigCB->SetTextAlign(11);
  TLatex *lSigG = new TLatex();
  lSigG->SetNDC(); lSigG->SetTextAlign(11);
  TLatex *lNG = new TLatex();
  lNG->SetNDC(); lNG->SetTextAlign(11);
  TLatex *lSigma = new TLatex();
  lSigma->SetNDC(); lSigma->SetTextAlign(11);
  TLatex *lAlpha = new TLatex();
  lAlpha->SetNDC(); lAlpha->SetTextAlign(11);
  TLatex *lN = new TLatex();
  lN->SetNDC(); lN->SetTextAlign(11);
  TLatex *lFG = new TLatex();
  lFG->SetNDC(); lFG->SetTextAlign(11);

  Double_t fx[2], fy[2], fex[2], fey[2];
  TGraphErrors *gDataLegend = new TGraphErrors(2,fx,fy,fex,fey);
  gDataLegend->SetName("gDataLegend");
  gDataLegend->SetTitle("gDataLegend");
  gDataLegend->SetMarkerStyle(20); gDataLegend->SetMarkerSize(0.8);
  TH1F *hBkgLegend = new TH1F("hBkgLegend","hBkgLegend",100,0,1);
  hBkgLegend->SetLineColor(kBlack); hBkgLegend->SetLineWidth(2); hBkgLegend->SetLineStyle(12); hBkgLegend->SetFillColor(kGray);
  TH1F *hTotalLegend = new TH1F("hTotalLegend","hTotalLegend",100,0,1);
  hTotalLegend->SetLineColor(kRed); hTotalLegend->SetLineWidth(2); hTotalLegend->SetFillColor(kYellow);
  TH1F *hMixLegend = new TH1F("hMixLegend","hMixLegend",100,0,1);
  hMixLegend->SetLineColor(kBlue); hMixLegend->SetLineWidth(2); hMixLegend->SetLineStyle(11);hMixLegend->SetFillColor(kAzure-9);
  TH1F *hTotalLegend_pp = new TH1F("hTotalLegend_pp","hTotalLegend_pp",100,0,1);
  hTotalLegend_pp->SetLineColor(kBlue); hTotalLegend_pp->SetLineWidth(2); hTotalLegend_pp->SetFillColor(kAzure-9);
  TH1F *hMixLegend_pp = new TH1F("hMixLegend_pp","hMixLegend_pp",100,0,1);
  hMixLegend_pp->SetLineColor(kRed); hMixLegend_pp->SetLineWidth(2); hMixLegend_pp->SetLineStyle(11);hMixLegend_pp->SetFillColor(kYellow);

  RooFitResult *fitM;
  RooFitResult *fitMP;
  
  RooRealVar *NJpsi[nFiles];
  RooRealVar *NBkg[nFiles];
  RooRealVar *NPsiP[nFiles];
  RooRealVar *NBkgP[nFiles];
  RooFormulaVar *fracP_pp;
  RooFormulaVar *doubleRatio[nFiles-1];
  RooFormulaVar *fracP_HI[nFiles-1];

  for (int i=nFiles-1; i>=0; --i) {
    NJpsi[i]  = new RooRealVar(("NJpsi_"+varSuffix.at(i)).c_str(),("J/psi yield in "+varSuffix.at(i)).c_str(),0.5*binData[i]->sumEntries(),0.0,2.0*binData[i]->sumEntries()); ws->import(*NJpsi[i]);
    
    NBkg[i]  = new RooRealVar(("NBkg_"+varSuffix.at(i)).c_str(),("Background yield in "+varSuffix.at(i)).c_str(), 0.5*binData[i]->sumEntries(),0.0,2.0*binData[i]->sumEntries()); ws->import(*NBkg[i]);

    NPsiP[i]  = new RooRealVar(("NPsiP_"+varSuffix.at(i)).c_str(),("psi(2S) yield in "+varSuffix.at(i)).c_str(),0.01*binData[i]->sumEntries(),0.0,2.0*binData[i]->sumEntries()); ws->import(*NPsiP[i]);

    NBkgP[i]  = new RooRealVar(("NBkgP_"+varSuffix.at(i)).c_str(),("Background yield in "+varSuffix.at(i)).c_str(), 0.5*binData[i]->sumEntries(),0.0,2.0*binData[i]->sumEntries()); ws->import(*NBkgP[i]);

    if (i == (int)nFiles-1) {
      fracP_pp = new RooFormulaVar(("fracP_"+varSuffix.at(i)).c_str(),("psi(2S) fraction in "+varSuffix.at(i)).c_str(),
				   "@0/@1",
				   RooArgList(*(ws->var(("NPsiP_"+varSuffix.at(i)).c_str())),								  *(ws->var(("NJpsi_"+varSuffix.at(i)).c_str()))));
      ws->import(*fracP_pp);
    }
    else {
      fracP_HI[i] = new RooFormulaVar(("fracP_"+varSuffix.at(i)).c_str(),("psi(2S) fraction in "+varSuffix.at(i)).c_str(),
				      "@0/@1",
				      RooArgList(*(ws->var(("NPsiP_"+varSuffix.at(i)).c_str())),
						 *(ws->var(("NJpsi_"+varSuffix.at(i)).c_str()))));
      ws->import(*fracP_HI[i]);

      doubleRatio[i] = new RooFormulaVar(("doubleRatio_"+varSuffix.at(i)).c_str(),
					 ("psi(2S) double ratio in "+varSuffix.at(i)).c_str(),"@0/@1",
					 RooArgList(*(ws->function(("fracP_"+varSuffix.at(i)).c_str())),
						    *(ws->function("fracP_pp"))));
      ws->import(*doubleRatio[i]);
    }

    if ( i == (int)nFiles-1 ) {
      sprintf(funct,"SUM::jpsiMassPDF_%s(NJpsi_%s*%s,NBkg_%s*%s)",
	      varSuffix.at(i).c_str(),
	      varSuffix.at(i).c_str(),
	      mJpsiFunct.c_str(),
	      varSuffix.at(i).c_str(),
	      mBkgFunct.at(i).c_str());

      sprintf(functP,"SUM::psipMassPDF_%s(NPsiP_%s*%s,NBkgP_%s*%s)",
	      varSuffix.at(i).c_str(),
	      varSuffix.at(i).c_str(),
	      mPsiPFunct.c_str(),
	      varSuffix.at(i).c_str(),
	      mBkgPFunct.at(i).c_str());
    }
    else {
      sprintf(funct,"SUM::jpsiMassPDF_%s(NJpsi_%s*%s_%s,NBkg_%s*%s)",
	      varSuffix.at(i).c_str(),
	      varSuffix.at(i).c_str(),
	      mJpsiFunct.c_str(),"HI",
	      varSuffix.at(i).c_str(),
	      mBkgFunct.at(i).c_str());

      sprintf(functP,"SUM::psipMassPDF_%s(NPsiP_%s*%s_%s,NBkgP_%s*%s)",
	      varSuffix.at(i).c_str(),
	      varSuffix.at(i).c_str(),
	      mPsiPFunct.c_str(),"HI",
	      varSuffix.at(i).c_str(),
	      mBkgPFunct.at(i).c_str());
    }
    cout << funct << endl;
    ws->factory(funct);

    cout << functP << endl;
    ws->factory(functP);
  }

  // Define quantities for 0-100%
  RooFormulaVar *NJpsi_HI0100 = new RooFormulaVar("NJpsi_HI0100","(@0+@1+@2)",
						  RooArgList(*(ws->var(("NJpsi_"+varSuffix.at(0)).c_str())),
							     *(ws->var(("NJpsi_"+varSuffix.at(1)).c_str())),
							     *(ws->var(("NJpsi_"+varSuffix.at(2)).c_str()))));
  ws->import(*NJpsi_HI0100);

  RooFormulaVar *NPsiP_HI0100 = new RooFormulaVar("NPsiP_HI0100", "@0+@1+@2",
						  RooArgList(*(ws->var(("NPsiP_"+varSuffix.at(0)).c_str())),
							     *(ws->var(("NPsiP_"+varSuffix.at(1)).c_str())),
							     *(ws->var(("NPsiP_"+varSuffix.at(2)).c_str()))));
  ws->import(*NPsiP_HI0100);

  RooFormulaVar *fracP_HI0100 = new RooFormulaVar("fracP_HI0100", "@0/@1",
						  RooArgList(*(ws->function("NPsiP_HI0100")),
							     *(ws->function("NJpsi_HI0100"))));
  ws->import(*fracP_HI0100);
  RooFormulaVar *doubleRatio_HI0100 = new RooFormulaVar("doubleRatio_HI0100","@0/@1",
							RooArgList(*(ws->function("fracP_HI0100")),
								   *(ws->function("fracP_pp"))));
  ws->import(*doubleRatio_HI0100);

  if (nFiles==2) {
    ws->factory("SIMUL::jpsiMassPDFSim(sample,HI=jpsiMassPDF_HI,pp=jpsiMassPDF_pp)");
    ws->factory("SIMUL::psipMassPDFSim(sample,HI=psipMassPDF_HI,pp=psipMassPDF_pp)");
  }
  else if (nFiles==4) {
    ws->factory("SIMUL::jpsiMassPDFSim(sample,HI020=jpsiMassPDF_HI020,HI2040=jpsiMassPDF_HI2040,HI40100=jpsiMassPDF_HI40100,pp=jpsiMassPDF_pp)");
    ws->factory("SIMUL::psipMassPDFSim(sample,HI020=psipMassPDF_HI020,HI2040=psipMassPDF_HI2040,HI40100=psipMassPDF_HI40100,pp=psipMassPDF_pp)");
  }
  else {
    index = "SIMUL::jpsiMassPDFSim(sample";
    char substr[5];
    for (unsigned int i=0;i<nFiles;++i) {
      sprintf(substr,",%u=jpsiMassPDF_%u",i,i);
      index += substr;
    }
    index += ")";
    ws->factory(index.c_str());

    index = "SIMUL::psipMassPDFSim(sample";
    for (unsigned int i=0;i<nFiles;++i) {
      sprintf(substr,",%u=psipMassPDF_%u",i,i);
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
    if (twoCB) {
      if (prange=="6.5-30.0" && yrange=="0.0-1.6") {
	ws->var("alpha_HI")->setVal(1.518);
	ws->var("enne_HI")->setVal(1.683);
	ws->var("wideFactor_HI")->setVal(1.698);
	ws->var("alpha")->setVal(1.935);
	ws->var("enne")->setVal(1.689);
	ws->var("wideFactor")->setVal(1.75);
      }
      else if (prange=="3.0-30.0" && yrange=="1.6-2.4") {
	ws->var("alpha_HI")->setVal(1.890);
	ws->var("enne_HI")->setVal(1.440);
	ws->var("wideFactor_HI")->setVal(1.506);
	ws->var("alpha")->setVal(2.152);
	ws->var("enne")->setVal(1.510);
	ws->var("wideFactor")->setVal(1.553);
      }
      else {
	ws->var("alpha_HI")->setVal(2.0);
	ws->var("enne_HI")->setVal(1.4);
	ws->var("wideFactor_HI")->setVal(2.0);
	ws->var("alpha")->setVal(2.0);
	ws->var("enne")->setVal(1.4);
	ws->var("wideFactor")->setVal(2.0);
      }
    }
    else {
      if (prange=="6.5-30.0" && yrange=="0.0-2.4") {
	ws->var("alpha_HI")->setVal(1.492);
	ws->var("enne_HI")->setVal(1.733);
	ws->var("wideFactor_HI")->setVal(1.895);
	ws->var("alpha")->setVal(1.696);
	ws->var("enne")->setVal(1.718);
	ws->var("wideFactor")->setVal(1.880);
      }
      else if (prange=="6.5-30.0" && yrange=="0.0-1.6") {
	ws->var("alpha_HI")->setVal(1.518);
	ws->var("enne_HI")->setVal(1.683);
	ws->var("wideFactor_HI")->setVal(1.698);
	ws->var("alpha")->setVal(1.722);
	ws->var("enne")->setVal(1.640);
	ws->var("wideFactor")->setVal(1.719);
      }
      else if (prange=="3.0-30.0" && yrange=="1.6-2.4") {
	ws->var("alpha_HI")->setVal(1.890);
	ws->var("enne_HI")->setVal(1.440);
	ws->var("wideFactor_HI")->setVal(1.506);
	ws->var("alpha")->setVal(2.059);
	ws->var("enne")->setVal(1.410);
	ws->var("wideFactor")->setVal(1.553);
      }
      else if (prange=="3.0-6.5" && yrange=="1.6-2.4") {
	ws->var("alpha_HI")->setVal(1.860);
	ws->var("enne_HI")->setVal(1.380);
	ws->var("wideFactor_HI")->setVal(1.550);
	ws->var("alpha")->setVal(2.040);
	ws->var("enne")->setVal(1.310);
	ws->var("wideFactor")->setVal(1.508);
      }
      else if (prange=="6.5-30.0" && yrange=="1.6-2.4") {
	ws->var("alpha_HI")->setVal(2.199);
	ws->var("enne_HI")->setVal(1.276);
	ws->var("wideFactor_HI")->setVal(1.830);
	ws->var("alpha")->setVal(2.143);
	ws->var("enne")->setVal(1.410);
	ws->var("wideFactor")->setVal(1.970);
      }
      else {
	ws->var("alpha_HI")->setVal(2.0);
	ws->var("enne_HI")->setVal(1.4);
	ws->var("wideFactor_HI")->setVal(2.0);
	ws->var("alpha")->setVal(2.0);
	ws->var("enne")->setVal(1.4);
	ws->var("wideFactor")->setVal(2.0);
      }
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
    // wideFactor is always free
    ws->var("wideFactor")->setConstant(kFALSE);
  }

  // 20140128: seed with CB fit parameters
  if ( found!=string::npos ) {
    string inputFNcb;
    if ( !fixAlpha || !fixN || !fixGwidth) {// read for free alpha, n, or wideFactor from the default CBG fit
      inputFNcb =  dirPre + "_PbPb" + mBkgFunct.front() + "_" + mBkgPFunct.front();
      for (unsigned int i=1; i<nFiles-1; ++i) {
	inputFNcb += "_" + mBkgFunct.at(i);
	inputFNcb += "_" + mBkgPFunct.at(i);
      }
      inputFNcb += "_pp" + mBkgFunct.back() + "_" + mBkgPFunct.back() + "_rap" + yrange_str + "_pT" + prange_str + "_cent" + crange + "_allVars.txt";
    }
    else { // read results from CB fit
      inputFNcb =  dirPre2 + "_PbPb" + mBkgFunct.front() + "_" + mBkgPFunct.front();
      for (unsigned int i=1; i<nFiles-1; ++i) {
	inputFNcb += "_" + mBkgFunct.at(i);
	inputFNcb += "_" + mBkgPFunct.at(i);
      }
      inputFNcb += "_pp" + mBkgFunct.back() + "_" + mBkgPFunct.back() + "_rap" + yrange_str + "_pT" + prange_str + "_cent" + crange + "_allVars.txt";
    }

    RooArgSet *set = ws->pdf("jpsiMassPDFSim")->getParameters(*(ws->var("Jpsi_Mass")));
    set->readFromFile(inputFNcb.c_str());
    // remove CB tail parameters, keep seed from MC instead
    //    set->remove(RooArgList("alpha,alpha_HI,enne,enne_HI"));
    cout << "Import variable values from: " << inputFNcb << endl;
    ws->import(*set);

    RooArgSet *setP = ws->pdf("psipMassPDFSim")->getParameters(*(ws->var("Jpsi_Mass")));
    setP->readFromFile(inputFNcb.c_str());
    // remove CB tail parameters, keep seed from MC instead
    //    set->remove(RooArgList("alpha,alpha_HI,enne,enne_HI"));
    ws->import(*setP);

    if (!fixAlpha) {
      ws->var("alpha")->setConstant(kFALSE);
      ws->var("alpha_HI")->setConstant(kFALSE);
    }
    if (!fixN) {
      ws->var("enne")->setConstant(kFALSE);
      ws->var("enne_HI")->setConstant(kFALSE);
    }
    if (!fixGwidth) {
      ws->var("wideFactor")->setConstant(kFALSE);
      ws->var("wideFactor_HI")->setConstant(kFALSE);
    }
  }


  ws->pdf("jpsiMassPDFSim")->Print("v");
  fitM = ws->pdf("jpsiMassPDFSim")->fitTo(*redDataSim,Extended(1),Hesse(1),Save(1),SumW2Error(0),NumCPU(8),Range("jpsi"),SplitRange(0));
  // PrintEvalErrors(0),Verbose(0),
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
    cout << "_" << mBkgFunct.at(i);
  }
  cout << "_pp" << mBkgFunct.back() << "_rap" << yrange_str << "_pT" << prange_str << "_cent" << crange << fix_str << endl;

  cout << "---FIT result summary: " << edmStatus << covStatus << fitStatus;
  for (unsigned int j=0; j<fitM->numStatusHistory(); ++j) {
    cout << fitM->statusCodeHistory(j);
  }
  cout << " ---" << endl;

  int nFitParam = fitM->floatParsFinal().getSize();
 
  /*
  ws->var("meanSig1")->setConstant(true);
  ws->var("sigmaSig1")->setConstant(true);
  ws->var("wideFactor")->setConstant(true);
  fitMP=0;

  fitMP = ws->pdf("psipMassPDFSim")->fitTo(*redDataSim,Extended(1),Hesse(1),Save(1),SumW2Error(0),NumCPU(8),PrintEvalErrors(0),Verbose(0),Range("psip"));
  fitMP->printMultiline(cout,1,1,"");

  ws->var("meanSig1")->setConstant(false);
  ws->var("sigmaSig1")->setConstant(false);
  ws->var("wideFactor")->setConstant(false);

  int edmStatusP = 1;
  if (fitMP->edm()<1.0e-03)
    edmStatusP = 0;

  int covStatusP = 1;
  if (fitMP->covQual()==3)
    covStatusP = 0;

  int fitStatusP = fitMP->status();
  if (found!=string::npos)
    cout << "CBG_PbPb";
  else
    cout << "CB_PbPb";

  for (unsigned int i=0; i<nFiles-1; ++i) {
    cout << "_" << mBkgPFunct.at(i);
  }
  cout << "_pp" << mBkgPFunct.back() << "_rap" << yrange_str << "_pT" << prange_str << "_cent" << crange << fix_str << endl;

  cout << "---FIT result summary: " << edmStatusP << covStatusP << fitStatusP;
  for (unsigned int j=0; j<fitMP->numStatusHistory(); ++j) {
    cout << fitMP->statusCodeHistory(j);
  }
  cout << " ---" << endl;

  //  int nFitParamP = fitMP->floatParsFinal().getSize();
  */
  resultF.cd();
  fitM->Write();
  //  fitMP->Write();
  cout << "fitM->Write() into: " << resultFN << endl;
  cout << "fitMP->Write() into: " << resultFN << endl;
  resultF.Close();

  // *** Draw mass plot
  RooPlot *mframe[nFiles];
  RooPlot *mframezoom[nFiles];
  RooPlot *mframepull[nFiles];
  RooHist *hpull[nFiles]; 
  TH1 *hdata[nFiles];
  TH1F *hpull_proj[nFiles];

  // *** TFile for saving plots
  string plotFN;
  plotFN = dirPre + "_PbPb" + mBkgFunct.front() + "_" + mBkgPFunct.front();
  for (unsigned int i=1; i<nFiles-1; ++i) { 
    plotFN += "_" + mBkgFunct.at(i);
    plotFN += "_" + mBkgPFunct.at(i);
 }
  plotFN += "_pp" + mBkgFunct.back() + "_" + mBkgPFunct.back() + "_rap" + yrange_str + "_pT" + prange_str + "_cent" + crange + fix_str + "plots.root";
  TFile plotF(plotFN.c_str(),"RECREATE");
  //    plotF.cd();
  TDirectory *dir[nFiles];
  string name;

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
    //    mframe[i]->GetYaxis()->SetTitle("Events");
    mframe[i]->SetTitleOffset(1.10,"X");
    mframe[i]->SetTitleOffset(1.05,"Y");
    mframe[i]->GetYaxis()->SetTitle("Events / (0.027 GeV/c^{2})");


    mframezoom[i] = ws->var("Jpsi_Mass")->frame(3.42,3.960);//3.3,4.1
    mframezoom[i]->SetName(("zoom_frame_"+varSuffix.at(i)).c_str());
    mframezoom[i]->SetTitle(("A RooPlot of \"J/#psi mass\" in "+varSuffix.at(i)).c_str());
    mframezoom[i]->GetYaxis()->SetTitle(mframe[i]->GetYaxis()->GetTitle());

    if (yrange=="1.6-2.4" && i==1) {
      redData[i]->plotOn(mframe[i],DataError(RooAbsData::SumW2),XErrorSize(0),MarkerSize(0.8),Binning(*rbm[i+1]));
      redData[i]->plotOn(mframezoom[i],DataError(RooAbsData::SumW2),XErrorSize(0),MarkerSize(0.8),Binning(*rbm[i+1]));
    }
    else {
      redData[i]->plotOn(mframe[i],DataError(RooAbsData::SumW2),XErrorSize(0),MarkerSize(0.8),Binning(*rbm[i]));
      redData[i]->plotOn(mframezoom[i],DataError(RooAbsData::SumW2),XErrorSize(0),MarkerSize(0.8),Binning(*rbm[i]));
    }
    titlestr = "2D fit for" + partTit + "muons (mass projection), p_{T} = " + prange + " GeV/c and |y| = " + yrange;
    mframe[i]->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
    mframe[i]->GetXaxis()->CenterTitle(1);
    mframezoom[i]->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
    mframezoom[i]->GetXaxis()->CenterTitle(1);
    double max = mframe[i]->GetMaximum() * 1.3;
    double min = 0.0;
    mframe[i]->SetMaximum(max);

    name = "NBkg_"+varSuffix.at(i);
    min = ws->var(name.c_str())->getVal()/(double(nbins)) * 0.7;

    string signal, signalP;
    // plot lines
    if (isPbPb) {
      ws->pdf(("jpsiMassPDF_"+varSuffix.at(i)).c_str())->plotOn(mframe[i],LineColor(kRed),LineWidth(2),Precision(1e-4),Range("jpsi"));
      //      ws->pdf(("jpsiMassPDF_"+varSuffix.at(i)).c_str())->plotOn(mframezoom[i],LineColor(kRed),LineWidth(2),Precision(1e-4),Range("jpsi"));

      ws->pdf(("psipMassPDF_"+varSuffix.at(i)).c_str())->plotOn(mframe[i],LineColor(kRed),LineWidth(2),Precision(1e-4),Range("psip"));
      ws->pdf(("psipMassPDF_"+varSuffix.at(i)).c_str())->plotOn(mframezoom[i],LineColor(kRed),LineWidth(2),Precision(1e-4),Range("psip"));
    }
    else {
      ws->pdf(("jpsiMassPDF_"+varSuffix.at(i)).c_str())->plotOn(mframe[i],LineColor(kBlue),LineWidth(2),Precision(1e-4),Range("jpsi"));
      //      ws->pdf(("jpsiMassPDF_"+varSuffix.at(i)).c_str())->plotOn(mframezoom[i],LineColor(kBlue),LineWidth(2),Precision(1e-4),Range("jpsi"));


      ws->pdf(("psipMassPDF_"+varSuffix.at(i)).c_str())->plotOn(mframe[i],LineColor(kBlue),LineWidth(2),Precision(1e-4),Range("psip"));
      ws->pdf(("psipMassPDF_"+varSuffix.at(i)).c_str())->plotOn(mframezoom[i],LineColor(kBlue),LineWidth(2),Precision(1e-4),Range("psip"));
    }

    hpull[i] = mframe[i]->pullHist(0,0,true);
    hpull[i]->SetName(("hpullhist_"+varSuffix.at(i)).c_str());
    hpull[i]->SetTitle(("hpullhist_"+varSuffix.at(i)).c_str());
    
    if (!isPaper && found!=string::npos) { 
      if (isPbPb && !shareShape) {
       	signal = ",signalCB1_HI";
	signalP = "signalCB1P_HI";
      }
      else {
	signal = ",signalCB1";
	signalP = ",signalCB1P";
      }
      ws->pdf(("jpsiMassPDF_"+varSuffix.at(i)).c_str())->plotOn(mframe[i],Components((mBkgFunct.at(i)+signal).c_str()),LineColor(kOrange+2),LineStyle(kDashed),LineWidth(2),Precision(1e-4),Range("jpsi"));
      //      ws->pdf(("jpsiMassPDF_"+varSuffix.at(i)).c_str())->plotOn(mframezoom[i],Components((mBkgFunct.at(i)+signal).c_str()),LineColor(kOrange+2),LineStyle(kDashed),LineWidth(2),Precision(1e-4),Range("jpsi"));

      ws->pdf(("psipMassPDF_"+varSuffix.at(i)).c_str())->plotOn(mframe[i],Components((mBkgPFunct.at(i)+signalP).c_str()),LineColor(kOrange+2),LineStyle(kDashed),LineWidth(2),Precision(1e-4),Range("psip"));
      ws->pdf(("psipMassPDF_"+varSuffix.at(i)).c_str())->plotOn(mframezoom[i],Components((mBkgPFunct.at(i)+signalP).c_str()),LineColor(kOrange+2),LineStyle(kDashed),LineWidth(2),Precision(1e-4),Range("psip"));
      if (isPbPb && !shareShape) {
	if (twoCB) {
	  signal = ",signalCB2_HI";
	  signalP = ",signalCB2P_HI";
	}
	else {
	  signal = ",signalG2_HI";
	  signalP = ",signalG2P_HI";
	}
      }
      else {
	if (twoCB) {
	  signal = ",signalCB2";
	  signalP = ",signalCB2P";
	}
	else {
	  signal = ",signalG2";
	  signalP = ",signalG2P";
	}
      }
      ws->pdf(("jpsiMassPDF_"+varSuffix.at(i)).c_str())->plotOn(mframe[i],Components((mBkgFunct.at(i)+signal).c_str()),LineColor(kGreen+2),LineStyle(kDashed),LineWidth(2),Precision(1e-4),Range("jpsi"));
      //      ws->pdf(("jpsiMassPDF_"+varSuffix.at(i)).c_str())->plotOn(mframezoom[i],Components((mBkgFunct.at(i)+signal).c_str()),LineColor(kGreen+2),LineStyle(kDashed),LineWidth(2),Precision(1e-4),Range("jpsi"));

      ws->pdf(("psipMassPDF_"+varSuffix.at(i)).c_str())->plotOn(mframe[i],Components((mBkgPFunct.at(i)+signalP).c_str()),LineColor(kGreen+2),LineStyle(kDashed),LineWidth(2),Precision(1e-4),Range("psip"));
      ws->pdf(("psipMassPDF_"+varSuffix.at(i)).c_str())->plotOn(mframezoom[i],Components((mBkgPFunct.at(i)+signalP).c_str()),LineColor(kGreen+2),LineStyle(kDashed),LineWidth(2),Precision(1e-4),Range("psip"));
    }

    ws->pdf(("jpsiMassPDF_"+varSuffix.at(i)).c_str())->plotOn(mframe[i],Components(mBkgFunct.at(i).c_str()),LineColor(kBlack),LineStyle(12),LineWidth(2),Precision(1e-4),Range("jpsi"));
    //    ws->pdf(("jpsiMassPDF_"+varSuffix.at(i)).c_str())->plotOn(mframezoom[i],Components(mBkgFunct.at(i).c_str()),LineColor(kBlack),LineStyle(12),LineWidth(2),Precision(1e-4),Range("jpsi"));


    ws->pdf(("psipMassPDF_"+varSuffix.at(i)).c_str())->plotOn(mframe[i],Components(mBkgPFunct.at(i).c_str()),LineColor(kBlack),LineStyle(12),LineWidth(2),Precision(1e-4),Range("psip"));
    ws->pdf(("psipMassPDF_"+varSuffix.at(i)).c_str())->plotOn(mframezoom[i],Components(mBkgPFunct.at(i).c_str()),LineColor(kBlack),LineStyle(12),LineWidth(2),Precision(1e-4),Range("psip"));

    // redraw total line
    if (isPbPb) {
      ws->pdf(("jpsiMassPDF_"+varSuffix.at(i)).c_str())->plotOn(mframe[i],LineColor(kRed),LineWidth(2),Precision(1e-4),Range("jpsi"));
      //      ws->pdf(("jpsiMassPDF_"+varSuffix.at(i)).c_str())->plotOn(mframezoom[i],LineColor(kRed),LineWidth(2),Precision(1e-4),Range("jpsi"));

      ws->pdf(("psipMassPDF_"+varSuffix.at(i)).c_str())->plotOn(mframe[i],LineColor(kRed),LineWidth(2),Precision(1e-4),Range("psip"));
      ws->pdf(("psipMassPDF_"+varSuffix.at(i)).c_str())->plotOn(mframezoom[i],LineColor(kRed),LineWidth(2),Precision(1e-4),Range("psip"));
    }
    else {
      ws->pdf(("jpsiMassPDF_"+varSuffix.at(i)).c_str())->plotOn(mframe[i],LineColor(kBlue),LineWidth(2),Precision(1e-4),Range("jpsi"));
      //      ws->pdf(("jpsiMassPDF_"+varSuffix.at(i)).c_str())->plotOn(mframezoom[i],LineColor(kBlue),LineWidth(2),Precision(1e-4),Range("jpsi"));

      ws->pdf(("psipMassPDF_"+varSuffix.at(i)).c_str())->plotOn(mframe[i],LineColor(kBlue),LineWidth(2),Precision(1e-4),Range("psip"));
      ws->pdf(("psipMassPDF_"+varSuffix.at(i)).c_str())->plotOn(mframezoom[i],LineColor(kBlue),LineWidth(2),Precision(1e-4),Range("psip"));
    }

    hdata[i] = redData[i]->createHistogram(("hdata"+varSuffix.at(i)).c_str(),*ws->var("Jpsi_Mass"),Binning(*rbm[i]));

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
    nFitParam = ws->pdf(("jpsiMassPDF_"+varSuffix.at(i)).c_str())->getParameters(redData[i])->selectByAttrib("Constant",kFALSE)->getSize();
    if (i!=nFiles-1) nFitParam--; // do not count fracP_pp as free parameter in PbPb fits
    
    int Dof = nFullBinsPull - nFitParam;
    if (Dof!=0)
      Chi2 /= Dof;

    double theNLL=0;
    theNLL = fitM->minNll();
    
    mframepull[i] =  ws->var("Jpsi_Mass")->frame(Title(("Pull of "+varSuffix.at(i)).c_str()));
    mframepull[i]->SetName(("pull_frame_"+varSuffix.at(i)).c_str());
    mframepull[i]->GetYaxis()->SetTitle("Pull ");
    mframepull[i]->SetLabelSize(0.08,"XYZ");
    mframepull[i]->SetTitleSize(0.1,"XYZ");
    mframepull[i]->SetTitleOffset(1.05,"X");
    mframepull[i]->SetTitleOffset(0.51,"Y");
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

    TCanvas c00("c00","c00"); c00.cd();

    TPad *pad1 = new TPad("pad1","This is pad1",0.0,0.35,1.00,1.00);
    pad1->SetBottomMargin(0); 
    TPad *pad2 = new TPad("pad2","This is pad2",0.0,0.00,1.00,0.35);
    pad2->SetTopMargin(0);  pad2->SetBottomMargin(0.24); 
    TPad *pad3 = new TPad("pad3","This is pad3",0.16,0.52,0.40,0.75);
    pad3->SetLeftMargin(0.2);
    pad3->SetBottomMargin(0.2);
    //    TPad *pad4 = new TPad("pad4","This is pad4",0.52,0.33,0.90,0.71);
    TPad *pad4 = new TPad("pad4","This is pad4",0.52,0.37,0.92,0.77);
    pad4->SetFillStyle(0);
    pad4->SetLeftMargin(0.2);
    pad4->SetBottomMargin(0.2);
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
      
      if (highStats) {
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
	    mframe[i]->GetYaxis()->SetRangeUser(10,2000);
	  }
	  else if ( crange == "20-40" ) {
	    mframe[i]->GetYaxis()->SetRangeUser(100,3000);
	  }
	  else {
	    mframe[i]->GetYaxis()->SetRangeUser(300,3000);
	  }
	}
	else {
	  mframe[i]->SetMaximum(20*max);
	  mframe[i]->SetMinimum(0.5*min);
	}
      }
      else {
	switch (i) {
	case 0:
	  if (yrange == "0.0-1.6")
	    mframe[i]->GetYaxis()->SetRangeUser(40,5000);
	  else if (yrange == "1.6-2.4")
	    mframe[i]->GetYaxis()->SetRangeUser(300,3000);
	  break;
	case 1:
	  if (yrange == "0.0-1.6")
	    mframe[i]->GetYaxis()->SetRangeUser(8,5000);
	  else if (yrange == "1.6-2.4")
	    // paper
	    //	    mframe[i]->GetYaxis()->SetRangeUser(120,2000);
	    // preliminary
	    mframe[i]->GetYaxis()->SetRangeUser(120,3000);
	  break;
	case 2:
	  if (yrange == "0.0-1.6")
	    mframe[i]->GetYaxis()->SetRangeUser(3,5000);
	  else if (yrange == "1.6-2.4")
	    mframe[i]->GetYaxis()->SetRangeUser(20,2000);
	  break;
	case 3:
	  if (yrange == "0.0-1.6")
	    mframe[i]->GetYaxis()->SetRangeUser(10,20000);
	  else if (yrange == "1.6-2.4")
	    mframe[i]->GetYaxis()->SetRangeUser(20,10000);
	  break;
	default:
	  mframe[i]->GetYaxis()->SetRangeUser(1,10000);
	  break;
	}
      }
    }
    else {
      c00.SetLeftMargin(0.14);
      c00.SetRightMargin(0.045);
      mframe[i]->GetYaxis()->SetTitleOffset(1.35);
      switch (i) {
      case 0:
	if (yrange == "0.0-1.6")
	  mframe[i]->GetYaxis()->SetRangeUser(0,1600);
	else if (yrange == "1.6-2.4")
	  mframe[i]->GetYaxis()->SetRangeUser(300,1800);
	break;
      case 1:
	if (yrange == "0.0-1.6")
	  mframe[i]->GetYaxis()->SetRangeUser(0,1000);
	else if (yrange == "1.6-2.4")
	  mframe[i]->GetYaxis()->SetRangeUser(100,800);
	break;
      case 2:
	if (yrange == "0.0-1.6")
	  mframe[i]->GetYaxis()->SetRangeUser(0,500);
	else if (yrange == "1.6-2.4")
	  mframe[i]->GetYaxis()->SetRangeUser(0,400);
	break;
      case 3:
	if (yrange == "0.0-1.6")
	  mframe[i]->GetYaxis()->SetRangeUser(0,6000);
	else if (yrange == "1.6-2.4")
	  mframe[i]->GetYaxis()->SetRangeUser(0,4500);
	break;
      default:
	mframe[i]->GetYaxis()->SetRangeUser(1,10000);
	break;
      }
    }

    double minX = 0.18;
    double maxX = 0.48;

    double minX2 = 0.64;
    double maxX2 = 0.94;

    if (logScale) {
      minX = 0.16;
      maxX = 0.46;
      minX2 = 0.62;
      maxX2 = 0.92;
    }

    if (isPaper) {
      // paper
      minX2 = 0.60;
      if (!logScale)
	  minX2 = 0.62;
      // preliminary
      // minX2 = 0.65;
    }

    double maxY = 0.9;
    double step = 0.05;
    double stepLarge = 0.07;

    lCMS->SetTextSize(0.05);
    if (true) { // paper
      if (isPbPb)
	lCMS->SetText(minX,maxY-0.01,"CMS PbPb #sqrt{s_{NN}} = 2.76 TeV");
      else
	lCMS->SetText(minX,maxY-0.01,"CMS pp #sqrt{s} = 2.76 TeV");

      mframe[i]->addObject(lCMS,"");

      maxY-=stepLarge;
    }
    else { // prelimimary
      lCMS->SetTextSize(0.045);
      lPre->SetTextSize(0.045);
      lPre->SetText(0.60,maxY-0.01,"CMS Preliminary");
      mframe[i]->addObject(lPre,"");
      //      maxY-=0.05;
      if (isPbPb)
	lCMS->SetText(minX,maxY-0.01,"PbPb #sqrt{s_{NN}} = 2.76 TeV");
      else
	lCMS->SetText(minX,maxY-0.01,"pp #sqrt{s} = 2.76 TeV");
      
      mframe[i]->addObject(lCMS,"");
      maxY-=stepLarge;
    }

    lLumi->SetTextSize(0.035);
    if (isPbPb)
      lLumi->SetText(minX,maxY,"L_{int} = 150 #mub^{-1}");
    else
      lLumi->SetText(minX,maxY,"L_{int} = 5.4 pb^{-1}");
    mframe[i]->addObject(lLumi,"");

    maxY-=step;
    lRap->SetTextSize(0.035);
    if (isPbPb) { 
      if (ymin==0.0)
	sprintf(reduceDS,"Cent. %.0f-%.0f%%, |y| < %.1f",cmin,cmax,ymax);
      else
	sprintf(reduceDS,"Cent. %.0f-%.0f%%, %0.1f < |y| < %.1f",cmin,cmax,ymin,ymax);
    }
    else {
      if (ymin==0.0)
	sprintf(reduceDS,"|y| < %.1f",ymax);
      else
	sprintf(reduceDS,"%0.1f < |y| < %.1f",ymin,ymax);
    }
    lRap->SetText(minX,maxY,reduceDS);
    mframe[i]->addObject(lRap,"");

    maxY-=step;
    lPt->SetTextSize(0.035);
    if (pmin==6.5)
      sprintf(reduceDS,"%.1f < p_{T} < %.0f GeV/c",pmin,pmax);
    else if (pmax==6.5)
      sprintf(reduceDS,"%.0f < p_{T} < %.1f GeV/c",pmin,pmax);
    else
      sprintf(reduceDS,"%.0f < p_{T} < %.0f GeV/c",pmin,pmax);
  
    lPt->SetText(minX,maxY,reduceDS);
    mframe[i]->addObject(lPt,"");

    lNLL->SetTextSize(0.035);
    if (!isPaper) {
      sprintf(reduceDS,"Min. NLL = %0.2f",theNLL);
      lNLL->SetText(minX,0.68,reduceDS);
      mframe[i]->addObject(lNLL,"");
    }

    maxY=0.9;
    if (!isPaper) {
      sprintf(reduceDS,"#chi^{2}/dof = %0.1f/%d",UnNormChi2,Dof);
      lChi->SetText(minX2,maxY,reduceDS);
      maxY-=step;//0.85
      mframe[i]->addObject(lChi,"");
      sprintf(reduceDS,"p-value = %0.4f",TMath::Prob(UnNormChi2,Dof));
      lPval->SetText(minX2,maxY,reduceDS);
      maxY-=step;//0.8
      mframe[i]->addObject(lPval,"");

      name = "NJpsi_"+varSuffix.at(i);
      if (ws->var(name.c_str())->hasAsymError() && abs(-1.0*ws->var(name.c_str())->getErrorLo()/ws->var(name.c_str())->getErrorHi() - 1)>0.1)
	sprintf(reduceDS,"N_{J/#psi} = %0.0f^{+%0.0f}_{%0.0f}",ws->var(name.c_str())->getVal(),ws->var(name.c_str())->getErrorHi(),ws->var(name.c_str())->getErrorLo());
      else
	sprintf(reduceDS,"N_{J/#psi} = %0.0f #pm %0.0f",ws->var(name.c_str())->getVal(),ws->var(name.c_str())->getError());
      lNJpsi->SetText(minX2,maxY,reduceDS);
      maxY-=step;//0.75
      mframe[i]->addObject(lNJpsi,"");

      name = "NPsiP_"+varSuffix.at(i);
      if (ws->var(name.c_str())->hasAsymError() && abs(-1.0*ws->var(name.c_str())->getErrorLo()/ws->var(name.c_str())->getErrorHi() - 1)>0.1)
	sprintf(reduceDS,"N_{#psi(2S)} = %0.0f^{+%0.0f}_{%0.0f}",ws->var(name.c_str())->getVal(),ws->var(name.c_str())->getErrorHi(),ws->var(name.c_str())->getErrorLo());
      else
	sprintf(reduceDS,"N_{#psi(2S)} = %0.0f #pm %0.0f",ws->var(name.c_str())->getVal(),ws->var(name.c_str())->getError());
      lNPsiP->SetText(minX2,maxY,reduceDS);
      maxY-=step;//0.7
      mframe[i]->addObject(lNPsiP,"");
    }

    if (!isPaper) {
      name = "fracP_"+varSuffix.at(i);
      //      sprintf(reduceDS,"R_{#psi} = %0.3f #pm %0.3f",ws->function(name.c_str())->getVal(),ws->function(name.c_str())->getPropagatedError(*fitMP));
      lRpsi->SetText(minX2,maxY,reduceDS);
      maxY-=step;//0.65
      mframe[i]->addObject(lRpsi,"");

      if (isPbPb) {
      name = "doubleRatio_"+varSuffix.at(i);
      //      sprintf(reduceDS,"#chi_{#psi} = %0.3f #pm %0.3f",ws->function(name.c_str())->getVal(),ws->function(name.c_str())->getPropagatedError(*fitMP));
	lXpsi->SetText(minX2,maxY,reduceDS);
	maxY-=step;//0.6
	mframe[i]->addObject(lXpsi,"");
      }
    }

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
      lSigCB->SetText(minX2,maxY,reduceDS);
      maxY-=step;//0.55

      if  (found!=string::npos) {
	if (twoCB)
	  sprintf(reduceDS,"#sigma_{CB2} = %0.0f MeV/c^{2}",sigmaSig2*1000.0);
	else
	  sprintf(reduceDS,"#sigma_{G} = %0.0f MeV/c^{2}",sigmaSig2*1000.0);
	lSigG->SetText(minX2,maxY,reduceDS);
	maxY-=step;//0.5
	if (isPbPb && !shareShape)
	  name = "wideFactor_HI";
	else
	  name = "wideFactor";

	if (twoCB) {
	  if (ws->var(name.c_str())->isConstant())
	    sprintf(reduceDS,"n_{CB2} = %0.2f (fixed)",ws->var(name.c_str())->getVal());
	  else if (ws->var(name.c_str())->hasAsymError() && abs(-1.0*ws->var(name.c_str())->getErrorLo()/ws->var(name.c_str())->getErrorHi() - 1)>0.1)
	    sprintf(reduceDS,"n_{CB2} = %0.2f^{+%0.2f}_{%0.2f}",ws->var(name.c_str())->getVal(),ws->var(name.c_str())->getErrorHi(),ws->var(name.c_str())->getErrorLo());
	  else
	    sprintf(reduceDS,"n_{CB2} = %0.2f #pm %0.2f",ws->var(name.c_str())->getVal(),ws->var(name.c_str())->getError());
	}
	else {
	  if (ws->var(name.c_str())->isConstant())
	    sprintf(reduceDS,"n_{G} = %0.2f (fixed)",ws->var(name.c_str())->getVal());
	  else if (ws->var(name.c_str())->hasAsymError() && abs(-1.0*ws->var(name.c_str())->getErrorLo()/ws->var(name.c_str())->getErrorHi() - 1)>0.1)
	    sprintf(reduceDS,"n_{G} = %0.2f^{+%0.2f}_{%0.2f}",ws->var(name.c_str())->getVal(),ws->var(name.c_str())->getErrorHi(),ws->var(name.c_str())->getErrorLo());
	  else
	    sprintf(reduceDS,"n_{G} = %0.2f #pm %0.2f",ws->var(name.c_str())->getVal(),ws->var(name.c_str())->getError());
	}
	lNG->SetText(minX2,maxY,reduceDS);
	maxY-=step;//0.45
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
	if (twoCB) {
	  if (Errresol<0.0005) 
	    sprintf(reduceDS,"#sigma_{CB+CB2} = (%0.0f #pm %0.0f) MeV/c^{2}",resol*1000.0,1.0);
	  else
	    sprintf(reduceDS,"#sigma_{CB+CB2} = (%0.0f #pm %0.0f) MeV/c^{2}",resol*1000.0,Errresol*1000.0);
	}
	else {
	  if (Errresol<0.0005) 
	    sprintf(reduceDS,"#sigma_{CB+G} = (%0.0f #pm %0.0f) MeV/c^{2}",resol*1000.0,1.0);
	  else
	    sprintf(reduceDS,"#sigma_{CB+G} = (%0.0f #pm %0.0f) MeV/c^{2}",resol*1000.0,Errresol*1000.0);
	}
	lSigma->SetText(minX2,maxY,reduceDS);
	maxY-=step;//0.4
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
	lAlpha->SetText(minX2,maxY,reduceDS);
	maxY-=step;//0.35

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
	lN->SetText(minX2,maxY,reduceDS);
	maxY-=step;//0.3

	if (isPbPb && !shareShape)
	  name = "coeffGaus_HI";
	else
	  name = "coeffGaus";
	if (twoCB) {
	  if (ws->var(name.c_str())->hasAsymError() && abs(-1.0*ws->var(name.c_str())->getErrorLo()/ws->var(name.c_str())->getErrorHi() - 1)>0.1)
	    sprintf(reduceDS,"f_{CB2} = %0.3f^{+%0.3f}_{%0.3f}",ws->var(name.c_str())->getVal(),ws->var(name.c_str())->getErrorHi(),ws->var(name.c_str())->getErrorLo());
	  else
	    sprintf(reduceDS,"f_{CB2} = %0.2f #pm %0.2f",ws->var(name.c_str())->getVal(),ws->var(name.c_str())->getError());
	}
	else {
	  if (ws->var(name.c_str())->hasAsymError() && abs(-1.0*ws->var(name.c_str())->getErrorLo()/ws->var(name.c_str())->getErrorHi() - 1)>0.1)
	    sprintf(reduceDS,"f_{G} = %0.3f^{+%0.3f}_{%0.3f}",ws->var(name.c_str())->getVal(),ws->var(name.c_str())->getErrorHi(),ws->var(name.c_str())->getErrorLo());
	  else
	    sprintf(reduceDS,"f_{G} = %0.2f #pm %0.2f",ws->var(name.c_str())->getVal(),ws->var(name.c_str())->getError());
	}
	lFG->SetText(minX2,maxY,reduceDS);
	maxY-=step;//0.25

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
	lAlpha->SetText(minX2,maxY,reduceDS);
	maxY-=step;//0.5

	if (isPbPb && !shareShape)
	  name = "enne_HI";
	else
	  name = "enne";
	if (ws->var(name.c_str())->isConstant())
	  sprintf(reduceDS,"n = %0.2f (fixed)",ws->var(name.c_str())->getVal());
	else
	  sprintf(reduceDS,"n = %0.2f #pm %0.2f",ws->var(name.c_str())->getVal(),ws->var(name.c_str())->getError());
	lN->SetText(minX2,maxY,reduceDS);
	maxY-=step;//0.45
      }

      mframe[i]->addObject(lAlpha,"");
      mframe[i]->addObject(lN,"");
    }

    TLegend *leg1;
    if (showInsert)
      {
	if (overlay)  // paper+insert+overlay
	  leg1 = new TLegend(minX,0.6895-4*step,maxX,0.680,NULL,"brNDC");
	else // paper+insert
	  leg1 = new TLegend(minX,0.6895-3*step,maxX,0.680,NULL,"brNDC");
      }
    else if (overlay)
      // paper + overlay
      //      leg1 = new TLegend(minX,0.8625-4*step,0.92,0.853,NULL,"brNDC");
      leg1 = new TLegend(minX2,0.868-4*step,maxX2,0.868,NULL,"brNDC");
    else
      // paper
      leg1 = new TLegend(minX2,0.868-3*step,maxX2,0.868,NULL,"brNDC");

    leg1->SetFillStyle(0);
    leg1->SetFillColor(0);
    leg1->SetBorderSize(0);
    leg1->SetMargin(0.14);
    leg1->SetTextSize(0.035);
    leg1->AddEntry(gDataLegend,"data","P");
    if (isPbPb)
      leg1->AddEntry(hTotalLegend,"total fit","L");
    else
      leg1->AddEntry(hTotalLegend_pp,"total fit","L");
    leg1->AddEntry(hBkgLegend,"background","L");
    if (!isPbPb && overlay) {
      //      leg1->AddEntry(hMixLegend,"total with","L");
      leg1->AddEntry(hMixLegend_pp,"R_{#psi}^{PbPb 0-20%}","L");
    }
    else if (overlay)
      leg1->AddEntry(hMixLegend,"R_{#psi}^{pp}","L");

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
    else {
      if (showInsert && isPbPb) {
	switch (i) {
	case 0:
	  if (yrange == "0.0-1.6")
	    mframezoom[i]->GetYaxis()->SetRangeUser(50,130);
	  else if (yrange == "1.6-2.4")
	    mframezoom[i]->GetYaxis()->SetRangeUser(300,600);
	  break;
	case 1:
	  if (yrange == "0.0-1.6")
	    mframezoom[i]->GetYaxis()->SetRangeUser(0,60);
	  else if (yrange == "1.6-2.4")
	    mframezoom[i]->GetYaxis()->SetRangeUser(150,300);
	  break;
	case 2:
	  if (yrange == "0.0-1.6")
	    mframezoom[i]->GetYaxis()->SetRangeUser(0,30);
	  else if (yrange == "1.6-2.4")
	    mframezoom[i]->GetYaxis()->SetRangeUser(20,70);
	  break;
	default:
	  break;
	}
	pad4->Draw();
	pad4->cd();mframezoom[i]->Draw();
	pad4->RedrawAxis();
	c00.cd();
      }
    }
    gPad->RedrawAxis(); 
    titlestr = dirPre + "_" + mBkgFunct.at(i) + "_" + mBkgPFunct.at(i) + "_rap" + yrange_str  + "_pT" + prange_str + "_cent" + crange + fix_str +  "massfit.pdf";
    c00.SaveAs(titlestr.c_str());


    ws->import(*mframe[i]);
    ws->import(*mframezoom[i]);
    ws->import(*mframepull[i]);
    ws->import(*hpull[i]);
    ws->import(*hpull_proj[i]);

    dir[i] = plotF.mkdir(varSuffix.at(i).c_str());
    dir[i]->cd();

    mframe[i]->Write();
    mframezoom[i]->Write();
    mframepull[i]->Write();
    hpull[i]->Write();
    hpull_proj[i]->Write();
    plotF.cd();

    crange = crange_tmp; // set back to orignal crange
  }
  plotF.Close();


  string fname = dirPre + "_PbPb" + mBkgFunct.front() + "_" + mBkgPFunct.front();
  for (unsigned int i=1; i<nFiles-1;++i) {
    fname += "_" + mBkgFunct.at(i);
    fname += "_" + mBkgPFunct.at(i);
  }
  fname += "_pp" + mBkgFunct.back() + "_" + mBkgPFunct.back() + "_rap" + yrange_str  + "_pT" + prange_str + "_cent" + crange + fix_str + "Workspace.root";
  ws->writeToFile(fname.c_str(),kTRUE);

  fname = dirPre + "_PbPb" + mBkgFunct.front() + "_" + mBkgPFunct.front();
  for (unsigned int i=1; i<nFiles-1;++i) {
    fname += "_" + mBkgFunct.at(i);
    fname += "_" + mBkgPFunct.at(i);
  }
  fname += "_pp" + mBkgFunct.back() + "_" + mBkgPFunct.back() + "_rap" + yrange_str  + "_pT" + prange_str + "_cent" + crange + fix_str + "allVars.txt";
  ws->allVars().writeToFile(fname.c_str());

  double NJpsi_fin[nFiles];
  double ErrNJpsi_fin[nFiles];
  double NPsiP_fin[nFiles];
  double ErrNPsiP_fin[nFiles];

  double fracP_fin[nFiles];
  double ErrFracP_fin[nFiles];

  double NBkg_fin[nFiles];
  double ErrNBkg_fin[nFiles];

  double NBkgP_fin[nFiles];
  double ErrNBkgP_fin[nFiles];

  double DoubleRatio_fin[nFiles-1];
  double ErrDoubleRatio_fin[nFiles-1];

  for (unsigned int i=0; i<nFiles; i++) {
    name = "NJpsi_"+varSuffix.at(i);
    NJpsi_fin[i]= ws->var(name.c_str())->getVal();
    ErrNJpsi_fin[i] = ws->var(name.c_str())->getError();
    name = "NBkg_"+varSuffix.at(i);
    NBkg_fin[i] = ws->var(name.c_str())->getVal();
    ErrNBkg_fin[i] = ws->var(name.c_str())->getError();


    name = "NPsiP_"+varSuffix.at(i);
    NPsiP_fin[i]= ws->var(name.c_str())->getVal();
    ErrNPsiP_fin[i] = ws->var(name.c_str())->getError();
    name = "NBkgP_"+varSuffix.at(i);
    NBkgP_fin[i] = ws->var(name.c_str())->getVal();
    ErrNBkgP_fin[i] = ws->var(name.c_str())->getError();

    if (i<nFiles-1) {
      name = "doubleRatio_"+varSuffix.at(i);
      DoubleRatio_fin[i] = ws->function(name.c_str())->getVal();
      ErrDoubleRatio_fin[i] = 0;//ws->function(name.c_str())->getPropagatedError(*fitMP);
    }

    name = "fracP_"+varSuffix.at(i);
    fracP_fin[i]= ws->function(name.c_str())->getVal();
    ErrFracP_fin[i] = 0;//ws->function(name.c_str())->getPropagatedError(*fitMP);
  }
  
  // To check values of fit parameters
  cout << endl << "J/psi yields:\n" << endl;
  if (nFiles>2) {
    cout << "HI0100:" << endl;
    cout << "---------------------" << endl;
    cout << "NJpsi: " << ws->function("NJpsi_HI0100")->getVal() << " +/- " << 0/*ws->function("NJpsi_HI0100")->getPropagatedError(*fitMP)*/ << endl;
    cout << "R_psi(2S): " << ws->function("fracP_HI0100")->getVal() << " +/- " << 0/*ws->function("fracP_HI0100")->getPropagatedError(*fitMP)*/ << endl;
    cout << "Double Ratio: " << ws->function("doubleRatio_HI0100")->getVal() << " +/- " << 0/*ws->function("doubleRatio_HI0100")->getPropagatedError(*fitMP)*/ << endl;
    cout << endl;
  }

  for (unsigned int i=0; i<nFiles; i++) {
    cout << varSuffix.at(i) << ":"<< endl;
    cout << "---------------------" << endl;
    cout << "NBkg: "  << NBkg_fin[i] << " +/- " << ErrNBkg_fin[i] << endl;
    cout << "NBkgP: "  << NBkgP_fin[i] << " +/- " << ErrNBkgP_fin[i] << endl;
    cout << "NJpsi: " << NJpsi_fin[i] << " +/- " << ErrNJpsi_fin[i] << endl;
    cout << "NPsiP: " << NPsiP_fin[i] << " +/- " << ErrNPsiP_fin[i] << endl;
    cout << "R_psi(2S) Fit: "  << fracP_fin[i] << " +/- " << ErrFracP_fin[i] << endl;
    if (i<nFiles-1) {
      cout << "Double Ratio: " << DoubleRatio_fin[i] << " +/- " << ErrDoubleRatio_fin[i] << endl;
    }
    cout << endl;
  }
  
  for (unsigned int i=0; i<nFiles; ++i) {
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
  // 0th order polynomial
  ws->factory("Polynomial::pol0(Jpsi_Mass,{coeffPol0[0.0]})");
  // 1st order polynomial
  ws->factory("Polynomial::pol1(Jpsi_Mass,{coeffPol1[0.0]})");
  // 2nd order polynomial
  ws->factory("Polynomial::pol2(Jpsi_Mass,{coeffPol1, coeffPol2[0.0]})");
  // 3rd order polynomial
  ws->factory("Polynomial::pol3(Jpsi_Mass,{coeffPol1, coeffPol2, coeffPol3[0.0]})");
  // 4th order polynomial
  ws->factory("Polynomial::pol4(Jpsi_Mass,{coeffPol1, coeffPol2, coeffPol3, coeffPol4[0.0]})");
  // 5th order polynomial
  ws->factory("Polynomial::pol5(Jpsi_Mass,{coeffPol1, coeffPol2, coeffPol3, coeffPol4, coeffPol5[0.0]})");
  // 6th order polynomial
  ws->factory("Polynomial::pol6(Jpsi_Mass,{coeffPol1, coeffPol2, coeffPol3, coeffPol4, coeffPol5, coeffPol6[0.0]})");
  // 7th order polynomial
  ws->factory("Polynomial::pol7(Jpsi_Mass,{coeffPol1, coeffPol2, coeffPol3, coeffPol4, coeffPol5, coeffPol6, coeffPol7[0.0]})");
  ws->var("coeffPol0")->setConstant(true);
  ws->var("coeffPol1")->setConstant(false);
  ws->var("coeffPol2")->setConstant(false);
  ws->var("coeffPol3")->setConstant(false);
  ws->var("coeffPol4")->setConstant(false);
  ws->var("coeffPol5")->setConstant(false);
  ws->var("coeffPol6")->setConstant(false);
  ws->var("coeffPol7")->setConstant(false);

  // psi(2S) part
  ws->factory("Polynomial::pol0P(Jpsi_Mass,{coeffPol0P[0.0]})");
  ws->factory("Polynomial::pol1P(Jpsi_Mass,{coeffPol1P[0.0]})");
  ws->factory("Polynomial::pol2P(Jpsi_Mass,{coeffPol1P, coeffPol2P[0.0]})");
  ws->factory("Polynomial::pol3P(Jpsi_Mass,{coeffPol1P, coeffPol2P, coeffPol3P[0.0]})");
  ws->factory("Polynomial::pol4P(Jpsi_Mass,{coeffPol1P, coeffPol2P, coeffPol3P, coeffPol4P[0.0]})");
  ws->factory("Polynomial::pol5P(Jpsi_Mass,{coeffPol1P, coeffPol2P, coeffPol3P, coeffPol4P, coeffPol5P[0.0]})");
  ws->factory("Polynomial::pol6P(Jpsi_Mass,{coeffPol1P, coeffPol2P, coeffPol3P, coeffPol4P, coeffPol5P, coeffPol6P[0.0]})");
  ws->factory("Polynomial::pol7P(Jpsi_Mass,{coeffPol1P, coeffPol2P, coeffPol3P, coeffPol4P, coeffPol5P, coeffPol6P, coeffPol7P[0.0]})");
  ws->var("coeffPol0P")->setConstant(true);
  ws->var("coeffPol1P")->setConstant(false);
  ws->var("coeffPol2P")->setConstant(false);
  ws->var("coeffPol3P")->setConstant(false);
  ws->var("coeffPol4P")->setConstant(false);
  ws->var("coeffPol5P")->setConstant(false);
  ws->var("coeffPol6P")->setConstant(false);
  ws->var("coeffPol7P")->setConstant(false);
  
  // 1st order polynomial in exponential function
  ws->factory("Polynomial::expPol1Arg(Jpsi_Mass,{expCoeffPol1[-0.1]})");
  ws->factory("Exponential::expPol1(expPol1Arg,expCoeffPol0[1.0])");
  // 2nd order polynomial in exponential function
  ws->factory("Polynomial::expPol2Arg(Jpsi_Mass,{expCoeffPol1, expCoeffPol2[0.0]})");
  ws->factory("Exponential::expPol2(expPol2Arg,expCoeffPol0)");
  // 3rd order polynomial in exponential function
  ws->factory("Polynomial::expPol3Arg(Jpsi_Mass,{expCoeffPol1, expCoeffPol2, expCoeffPol3[0.0]})");
  ws->factory("Exponential::expPol3(expPol3Arg,expCoeffPol0)");
  // 4th order polynomial in exponential function
  ws->factory("Polynomial::expPol4Arg(Jpsi_Mass,{expCoeffPol1, expCoeffPol2, expCoeffPol3, expCoeffPol4[0.0]})");
  ws->factory("Exponential::expPol4(expPol4Arg,expCoeffPol0)");
  // 5th order polynomial in exponential function
  ws->factory("Polynomial::expPol5Arg(Jpsi_Mass,{expCoeffPol1, expCoeffPol2, expCoeffPol3, expCoeffPol4, expCoeffPol5[0.0]})");
  ws->factory("Exponential::expPol5(expPol5Arg,expCoeffPol0)");
  // 6th order polynomial in exponential function
  ws->factory("Polynomial::expPol6Arg(Jpsi_Mass,{expCoeffPol1, expCoeffPol2, expCoeffPol3, expCoeffPol4, expCoeffPol5, expCoeffPol6[0.0]})");
  ws->factory("Exponential::expPol6(expPol6Arg,expCoeffPol0)");
  // 7th order polynomial in exponential function
  ws->factory("Polynomial::expPol7Arg(Jpsi_Mass,{expCoeffPol1, expCoeffPol2, expCoeffPol3, expCoeffPol4, expCoeffPol5, expCoeffPol6, expCoeffPol7[0.0]})");
  ws->factory("Exponential::expPol7(expPol7Arg,expCoeffPol0)");

  ws->var("expCoeffPol0")->setConstant(true);
  ws->var("expCoeffPol1")->setConstant(false);
  ws->var("expCoeffPol2")->setConstant(false);
  ws->var("expCoeffPol3")->setConstant(false);
  ws->var("expCoeffPol4")->setConstant(false);
  ws->var("expCoeffPol5")->setConstant(false);
  ws->var("expCoeffPol6")->setConstant(false);
  ws->var("expCoeffPol7")->setConstant(false);


  // 1st order polynomial in exponential function
  ws->factory("Polynomial::expPol1ArgP(Jpsi_Mass,{expCoeffPol1P[-0.1]})");
  ws->factory("Exponential::expPol1P(expPol1ArgP,expCoeffPol0P[1.0])");
  ws->factory("Polynomial::expPol2ArgP(Jpsi_Mass,{expCoeffPol1P, expCoeffPol2P[0.0]})");
  ws->factory("Exponential::expPol2P(expPol2ArgP,expCoeffPol0P)");
  ws->factory("Polynomial::expPol3ArgP(Jpsi_Mass,{expCoeffPol1P, expCoeffPol2P, expCoeffPol3P[0.0]})");
  ws->factory("Exponential::expPol3P(expPol3ArgP,expCoeffPol0P)");
  ws->factory("Polynomial::expPol4ArgP(Jpsi_Mass,{expCoeffPol1P, expCoeffPol2P, expCoeffPol3P, expCoeffPol4P[0.0]})");
  ws->factory("Exponential::expPol4P(expPol4ArgP,expCoeffPol0P)");
  ws->factory("Polynomial::expPol5ArgP(Jpsi_Mass,{expCoeffPol1P, expCoeffPol2P, expCoeffPol3P, expCoeffPol4P, expCoeffPol5P[0.0]})");
  ws->factory("Exponential::expPol5P(expPol5ArgP,expCoeffPol0P)");
  ws->factory("Polynomial::expPol6ArgP(Jpsi_Mass,{expCoeffPol1P, expCoeffPol2P, expCoeffPol3P, expCoeffPol4P, expCoeffPol5P, expCoeffPol6P[0.0]})");
  ws->factory("Exponential::expPol6P(expPol6ArgP,expCoeffPol0P)");
  ws->factory("Polynomial::expPol7ArgP(Jpsi_Mass,{expCoeffPol1P, expCoeffPol2P, expCoeffPol3P, expCoeffPol4P, expCoeffPol5P, expCoeffPol6P, expCoeffPol7P[0.0]})");
  ws->factory("Exponential::expPol7P(expPol7ArgP,expCoeffPol0P)");

  ws->var("expCoeffPol0P")->setConstant(true);
  ws->var("expCoeffPol1P")->setConstant(false);
  ws->var("expCoeffPol2P")->setConstant(false);
  ws->var("expCoeffPol3P")->setConstant(false);
  ws->var("expCoeffPol4P")->setConstant(false);
  ws->var("expCoeffPol5P")->setConstant(false);
  ws->var("expCoeffPol6P")->setConstant(false);
  ws->var("expCoeffPol7P")->setConstant(false);

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
  // wide CB
  ws->factory("CBShape::signalCB2(Jpsi_Mass,meanSig1,sigmaSig2,alpha,enne)");

  // CB(narrow) + Gauss(wide)
  ws->factory("SUM::sigCB1G2(coeffGaus[0.1,0.0,1.0]*signalG2,signalCB1)");

  // CB(narrow) + CB(wide)
  ws->factory("SUM::sigCB1CB2(coeffGaus[0.1,0.0,1.0]*signalCB2,signalCB1)");

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
  // wide CB
  ws->factory("CBShape::signalCB2P(Jpsi_Mass,meanSig1P,sigmaSig2P,alpha,enne)");

  // CB(narrow) + Gauss(wide)
  ws->factory("SUM::sigCB1G2P(coeffGaus*signalG2P,signalCB1P)");
  // CB(narrow) + CB(wide)
  ws->factory("SUM::sigCB1CB2P(coeffGaus*signalCB2P,signalCB1P)");

  return;
}


void defineMassBkgHI(RooWorkspace *ws) {
  // 0th order polynomial
  ws->factory("Polynomial::pol0_HI(Jpsi_Mass,{coeffPol0_HI[0.0]})");
  // 1st order polynomial
  ws->factory("Polynomial::pol1_HI(Jpsi_Mass,{coeffPol1_HI[0.0]})");
  // 2nd order polynomial
  ws->factory("Polynomial::pol2_HI(Jpsi_Mass,{coeffPol1_HI, coeffPol2_HI[0.0]})");
  // 3rd order polynomial
  ws->factory("Polynomial::pol3_HI(Jpsi_Mass,{coeffPol1_HI, coeffPol2_HI, coeffPol3_HI[0.0]})");
  // 4th order polynomial
  ws->factory("Polynomial::pol4_HI(Jpsi_Mass,{coeffPol1_HI, coeffPol2_HI, coeffPol3_HI, coeffPol4_HI[0.0]})");
  // 5th order polynomial
  ws->factory("Polynomial::pol5_HI(Jpsi_Mass,{coeffPol1_HI, coeffPol2_HI, coeffPol3_HI, coeffPol4_HI, coeffPol5_HI[0.0]})");
  // 6th order polynomial
  ws->factory("Polynomial::pol6_HI(Jpsi_Mass,{coeffPol1_HI, coeffPol2_HI, coeffPol3_HI, coeffPol4_HI, coeffPol5_HI, coeffPol6_HI[0.0]})");
  // 7th order polynomial
  ws->factory("Polynomial::pol7_HI(Jpsi_Mass,{coeffPol1_HI, coeffPol2_HI, coeffPol3_HI, coeffPol4_HI, coeffPol5_HI, coeffPol6_HI, coeffPol7_HI[0.0]})");
  ws->var("coeffPol0_HI")->setConstant(true);
  ws->var("coeffPol1_HI")->setConstant(false);
  ws->var("coeffPol2_HI")->setConstant(false);
  ws->var("coeffPol3_HI")->setConstant(false);
  ws->var("coeffPol4_HI")->setConstant(false);
  ws->var("coeffPol5_HI")->setConstant(false);
  ws->var("coeffPol6_HI")->setConstant(false);
  ws->var("coeffPol7_HI")->setConstant(false);

  ws->factory("Polynomial::pol0_HIP(Jpsi_Mass,{coeffPol0_HIP[0.0]})");
  ws->factory("Polynomial::pol1_HIP(Jpsi_Mass,{coeffPol1_HIP[0.0]})");
  ws->factory("Polynomial::pol2_HIP(Jpsi_Mass,{coeffPol1_HIP, coeffPol2_HIP[0.0]})");
  ws->factory("Polynomial::pol3_HIP(Jpsi_Mass,{coeffPol1_HIP, coeffPol2_HIP, coeffPol3_HIP[0.0]})");
  ws->factory("Polynomial::pol4_HIP(Jpsi_Mass,{coeffPol1_HIP, coeffPol2_HIP, coeffPol3_HIP, coeffPol4_HIP[0.0]})");
  ws->factory("Polynomial::pol5_HIP(Jpsi_Mass,{coeffPol1_HIP, coeffPol2_HIP, coeffPol3_HIP, coeffPol4_HIP, coeffPol5_HIP[0.0]})");
  ws->factory("Polynomial::pol6_HIP(Jpsi_Mass,{coeffPol1_HIP, coeffPol2_HIP, coeffPol3_HIP, coeffPol4_HIP, coeffPol5_HIP, coeffPol6_HIP[0.0]})");
  ws->factory("Polynomial::pol7_HIP(Jpsi_Mass,{coeffPol1_HIP, coeffPol2_HIP, coeffPol3_HIP, coeffPol4_HIP, coeffPol5_HIP, coeffPol6_HIP, coeffPol7_HIP[0.0]})");
  ws->var("coeffPol0_HIP")->setConstant(true);
  ws->var("coeffPol1_HIP")->setConstant(false);
  ws->var("coeffPol2_HIP")->setConstant(false);
  ws->var("coeffPol3_HIP")->setConstant(false);
  ws->var("coeffPol4_HIP")->setConstant(false);
  ws->var("coeffPol5_HIP")->setConstant(false);
  ws->var("coeffPol6_HIP")->setConstant(false);
  ws->var("coeffPol7_HIP")->setConstant(false);

  // 1st order polynomial in exponential function
  ws->factory("Polynomial::expPol1Arg_HI(Jpsi_Mass,{expCoeffPol1_HI[-0.1]})");
  ws->factory("Exponential::expPol1_HI(expPol1Arg_HI,expCoeffPol0_HI[1.0])");
  // 2nd order polynomial in exponential function
  ws->factory("Polynomial::expPol2Arg_HI(Jpsi_Mass,{expCoeffPol1_HI, expCoeffPol2_HI[0.0]})");
  ws->factory("Exponential::expPol2_HI(expPol2Arg_HI,expCoeffPol0_HI)");
  // 3rd order polynomial in exponential function
  ws->factory("Polynomial::expPol3Arg_HI(Jpsi_Mass,{expCoeffPol1_HI, expCoeffPol2_HI, expCoeffPol3_HI[0.0]})");
  ws->factory("Exponential::expPol3_HI(expPol3Arg_HI,expCoeffPol0_HI)");
  // 4th order polynomial in exponential function
  ws->factory("Polynomial::expPol4Arg_HI(Jpsi_Mass,{expCoeffPol1_HI, expCoeffPol2_HI, expCoeffPol3_HI, expCoeffPol4_HI[0.0]})");
  ws->factory("Exponential::expPol4_HI(expPol4Arg_HI,expCoeffPol0_HI)");
  // 5th order polynomial in exponential function
  ws->factory("Polynomial::expPol5Arg_HI(Jpsi_Mass,{expCoeffPol1_HI, expCoeffPol2_HI, expCoeffPol3_HI, expCoeffPol4_HI, expCoeffPol5_HI[0.0]})");
  ws->factory("Exponential::expPol5_HI(expPol5Arg_HI,expCoeffPol0_HI)");
  // 6th order polynomial in exponential function
  ws->factory("Polynomial::expPol6Arg_HI(Jpsi_Mass,{expCoeffPol1_HI, expCoeffPol2_HI, expCoeffPol3_HI, expCoeffPol4_HI, expCoeffPol5_HI, expCoeffPol6_HI[0.0]})");
  ws->factory("Exponential::expPol6_HI(expPol6Arg_HI,expCoeffPol0_HI)");
  // 7th order polynomial in exponential function
  ws->factory("Polynomial::expPol7Arg_HI(Jpsi_Mass,{expCoeffPol1_HI, expCoeffPol2_HI, expCoeffPol3_HI, expCoeffPol4_HI, expCoeffPol5_HI, expCoeffPol6_HI, expCoeffPol7_HI[0.0]})");
  ws->factory("Exponential::expPol7_HI(expPol7Arg_HI,expCoeffPol0_HI)");

  ws->var("expCoeffPol0_HI")->setConstant(true);
  ws->var("expCoeffPol1_HI")->setConstant(false);
  ws->var("expCoeffPol2_HI")->setConstant(false);
  ws->var("expCoeffPol3_HI")->setConstant(false);
  ws->var("expCoeffPol4_HI")->setConstant(false);
  ws->var("expCoeffPol5_HI")->setConstant(false);
  ws->var("expCoeffPol6_HI")->setConstant(false);
  ws->var("expCoeffPol7_HI")->setConstant(false);


  // 0-20%
  // 0th order polynomial
  ws->factory("Polynomial::pol0_HI020(Jpsi_Mass,{coeffPol0_HI020[0.0]})");
  // 1st order polynomial
  ws->factory("Polynomial::pol1_HI020(Jpsi_Mass,{coeffPol1_HI020[0.0]})");
  // 2nd order polynomial
  ws->factory("Polynomial::pol2_HI020(Jpsi_Mass,{coeffPol1_HI020, coeffPol2_HI020[0.0]})");
  // 3rd order polynomial
  ws->factory("Polynomial::pol3_HI020(Jpsi_Mass,{coeffPol1_HI020, coeffPol2_HI020, coeffPol3_HI020[0.0]})");
  // 4th order polynomial
  ws->factory("Polynomial::pol4_HI020(Jpsi_Mass,{coeffPol1_HI020, coeffPol2_HI020, coeffPol3_HI020, coeffPol4_HI020[0.0]})");
  // 5th order polynomial
  ws->factory("Polynomial::pol5_HI020(Jpsi_Mass,{coeffPol1_HI020, coeffPol2_HI020, coeffPol3_HI020, coeffPol4_HI020, coeffPol5_HI020[0.0]})");
  // 6th order polynomial
  ws->factory("Polynomial::pol6_HI020(Jpsi_Mass,{coeffPol1_HI020, coeffPol2_HI020, coeffPol3_HI020, coeffPol4_HI020, coeffPol5_HI020, coeffPol6_HI020[0.0]})");
  // 7th order polynomial
  ws->factory("Polynomial::pol7_HI020(Jpsi_Mass,{coeffPol1_HI020, coeffPol2_HI020, coeffPol3_HI020, coeffPol4_HI020, coeffPol5_HI020, coeffPol6_HI020, coeffPol7_HI020[0.0]})");
  ws->var("coeffPol0_HI020")->setConstant(true);
  ws->var("coeffPol1_HI020")->setConstant(false);
  ws->var("coeffPol2_HI020")->setConstant(false);
  ws->var("coeffPol3_HI020")->setConstant(false);
  ws->var("coeffPol4_HI020")->setConstant(false);
  ws->var("coeffPol5_HI020")->setConstant(false);
  ws->var("coeffPol6_HI020")->setConstant(false);
  ws->var("coeffPol7_HI020")->setConstant(false);

  ws->factory("Polynomial::pol0_HI020P(Jpsi_Mass,{coeffPol0_HI020P[0.0]})");
  ws->factory("Polynomial::pol1_HI020P(Jpsi_Mass,{coeffPol1_HI020P[0.0]})");
  ws->factory("Polynomial::pol2_HI020P(Jpsi_Mass,{coeffPol1_HI020P, coeffPol2_HI020P[0.0]})");
  ws->factory("Polynomial::pol3_HI020P(Jpsi_Mass,{coeffPol1_HI020P, coeffPol2_HI020P, coeffPol3_HI020P[0.0]})");
  ws->factory("Polynomial::pol4_HI020P(Jpsi_Mass,{coeffPol1_HI020P, coeffPol2_HI020P, coeffPol3_HI020P, coeffPol4_HI020P[0.0]})");
  ws->factory("Polynomial::pol5_HI020P(Jpsi_Mass,{coeffPol1_HI020P, coeffPol2_HI020P, coeffPol3_HI020P, coeffPol4_HI020P, coeffPol5_HI020P[0.0]})");
  ws->factory("Polynomial::pol6_HI020P(Jpsi_Mass,{coeffPol1_HI020P, coeffPol2_HI020P, coeffPol3_HI020P, coeffPol4_HI020P, coeffPol5_HI020P, coeffPol6_HI020P[0.0]})");
  ws->factory("Polynomial::pol7_HI020P(Jpsi_Mass,{coeffPol1_HI020P, coeffPol2_HI020P, coeffPol3_HI020P, coeffPol4_HI020P, coeffPol5_HI020P, coeffPol6_HI020P, coeffPol7_HI020P[0.0]})");
  ws->var("coeffPol0_HI020P")->setConstant(true);
  ws->var("coeffPol1_HI020P")->setConstant(false);
  ws->var("coeffPol2_HI020P")->setConstant(false);
  ws->var("coeffPol3_HI020P")->setConstant(false);
  ws->var("coeffPol4_HI020P")->setConstant(false);
  ws->var("coeffPol5_HI020P")->setConstant(false);
  ws->var("coeffPol6_HI020P")->setConstant(false);
  ws->var("coeffPol7_HI020P")->setConstant(false);


  // 1st order polynomial in exponential function
  ws->factory("Polynomial::expPol1Arg_HI020(Jpsi_Mass,{expCoeffPol1_HI020[-0.1]})");
  ws->factory("Exponential::expPol1_HI020(expPol1Arg_HI020,expCoeffPol0_HI020[1.0])");
  // 2nd order polynomial in exponential function
  ws->factory("Polynomial::expPol2Arg_HI020(Jpsi_Mass,{expCoeffPol1_HI020, expCoeffPol2_HI020[0.0]})");
  ws->factory("Exponential::expPol2_HI020(expPol2Arg_HI020,expCoeffPol0_HI020)");
  // 3rd order polynomial in exponential function
  ws->factory("Polynomial::expPol3Arg_HI020(Jpsi_Mass,{expCoeffPol1_HI020, expCoeffPol2_HI020, expCoeffPol3_HI020[0.0]})");
  ws->factory("Exponential::expPol3_HI020(expPol3Arg_HI020,expCoeffPol0_HI020)");
  // 4th order polynomial in exponential function
  ws->factory("Polynomial::expPol4Arg_HI020(Jpsi_Mass,{expCoeffPol1_HI020, expCoeffPol2_HI020, expCoeffPol3_HI020, expCoeffPol4_HI020[0.0]})");
  ws->factory("Exponential::expPol4_HI020(expPol4Arg_HI020,expCoeffPol0_HI020)");
  // 5th order polynomial in exponential function
  ws->factory("Polynomial::expPol5Arg_HI020(Jpsi_Mass,{expCoeffPol1_HI020, expCoeffPol2_HI020, expCoeffPol3_HI020, expCoeffPol4_HI020, expCoeffPol5_HI020[0.0]})");
  ws->factory("Exponential::expPol5_HI020(expPol5Arg_HI020,expCoeffPol0_HI020)");
  // 6th order polynomial in exponential function
  ws->factory("Polynomial::expPol6Arg_HI020(Jpsi_Mass,{expCoeffPol1_HI020, expCoeffPol2_HI020, expCoeffPol3_HI020, expCoeffPol4_HI020, expCoeffPol5_HI020, expCoeffPol6_HI020[0.0]})");
  ws->factory("Exponential::expPol6_HI020(expPol6Arg_HI020,expCoeffPol0_HI020)");
  // 7th order polynomial in exponential function
  ws->factory("Polynomial::expPol7Arg_HI020(Jpsi_Mass,{expCoeffPol1_HI020, expCoeffPol2_HI020, expCoeffPol3_HI020, expCoeffPol4_HI020, expCoeffPol5_HI020, expCoeffPol6_HI020, expCoeffPol7_HI020[0.0]})");
  ws->factory("Exponential::expPol7_HI020(expPol7Arg_HI020,expCoeffPol0_HI020)");

  ws->var("expCoeffPol0_HI020")->setConstant(true);
  ws->var("expCoeffPol1_HI020")->setConstant(false);
  ws->var("expCoeffPol2_HI020")->setConstant(false);
  ws->var("expCoeffPol3_HI020")->setConstant(false);
  ws->var("expCoeffPol4_HI020")->setConstant(false);
  ws->var("expCoeffPol5_HI020")->setConstant(false);
  ws->var("expCoeffPol6_HI020")->setConstant(false);
  ws->var("expCoeffPol7_HI020")->setConstant(false);

  // 20-40%
  // 0th order polynomial
  ws->factory("Polynomial::pol0_HI2040(Jpsi_Mass,{coeffPol0_HI2040[0.0]})");
  // 1st order polynomial
  ws->factory("Polynomial::pol1_HI2040(Jpsi_Mass,{coeffPol1_HI2040[0.0]})");
  // 2nd order polynomial
  ws->factory("Polynomial::pol2_HI2040(Jpsi_Mass,{coeffPol1_HI2040, coeffPol2_HI2040[0.0]})");
  // 3rd order polynomial
  ws->factory("Polynomial::pol3_HI2040(Jpsi_Mass,{coeffPol1_HI2040, coeffPol2_HI2040, coeffPol3_HI2040[0.0]})");
  // 4th order polynomial
  ws->factory("Polynomial::pol4_HI2040(Jpsi_Mass,{coeffPol1_HI2040, coeffPol2_HI2040, coeffPol3_HI2040, coeffPol4_HI2040[0.0]})");
  // 5th order polynomial
  ws->factory("Polynomial::pol5_HI2040(Jpsi_Mass,{coeffPol1_HI2040, coeffPol2_HI2040, coeffPol3_HI2040, coeffPol4_HI2040, coeffPol5_HI2040[0.0]})");
  // 6th order polynomial
  ws->factory("Polynomial::pol6_HI2040(Jpsi_Mass,{coeffPol1_HI2040, coeffPol2_HI2040, coeffPol3_HI2040, coeffPol4_HI2040, coeffPol5_HI2040, coeffPol6_HI2040[0.0]})");
  // 7th order polynomial
  ws->factory("Polynomial::pol7_HI2040(Jpsi_Mass,{coeffPol1_HI2040, coeffPol2_HI2040, coeffPol3_HI2040, coeffPol4_HI2040, coeffPol5_HI2040, coeffPol6_HI2040, coeffPol7_HI2040[0.0]})");
  ws->var("coeffPol0_HI2040")->setConstant(true);
  ws->var("coeffPol1_HI2040")->setConstant(false);
  ws->var("coeffPol2_HI2040")->setConstant(false);
  ws->var("coeffPol3_HI2040")->setConstant(false);
  ws->var("coeffPol4_HI2040")->setConstant(false);
  ws->var("coeffPol5_HI2040")->setConstant(false);
  ws->var("coeffPol6_HI2040")->setConstant(false);
  ws->var("coeffPol7_HI2040")->setConstant(false);

  ws->factory("Polynomial::pol0_HI2040P(Jpsi_Mass,{coeffPol0_HI2040P[0.0]})");
  ws->factory("Polynomial::pol1_HI2040P(Jpsi_Mass,{coeffPol1_HI2040P[0.0]})");
  ws->factory("Polynomial::pol2_HI2040P(Jpsi_Mass,{coeffPol1_HI2040P, coeffPol2_HI2040P[0.0]})");
  ws->factory("Polynomial::pol3_HI2040P(Jpsi_Mass,{coeffPol1_HI2040P, coeffPol2_HI2040P, coeffPol3_HI2040P[0.0]})");
  ws->factory("Polynomial::pol4_HI2040P(Jpsi_Mass,{coeffPol1_HI2040P, coeffPol2_HI2040P, coeffPol3_HI2040P, coeffPol4_HI2040P[0.0]})");
  ws->factory("Polynomial::pol5_HI2040P(Jpsi_Mass,{coeffPol1_HI2040P, coeffPol2_HI2040P, coeffPol3_HI2040P, coeffPol4_HI2040P, coeffPol5_HI2040P[0.0]})");
  ws->factory("Polynomial::pol6_HI2040P(Jpsi_Mass,{coeffPol1_HI2040P, coeffPol2_HI2040P, coeffPol3_HI2040P, coeffPol4_HI2040P, coeffPol5_HI2040P, coeffPol6_HI2040P[0.0]})");
  ws->factory("Polynomial::pol7_HI2040P(Jpsi_Mass,{coeffPol1_HI2040P, coeffPol2_HI2040P, coeffPol3_HI2040P, coeffPol4_HI2040P, coeffPol5_HI2040P, coeffPol6_HI2040P, coeffPol7_HI2040P[0.0]})");
  ws->var("coeffPol0_HI2040P")->setConstant(true);
  ws->var("coeffPol1_HI2040P")->setConstant(false);
  ws->var("coeffPol2_HI2040P")->setConstant(false);
  ws->var("coeffPol3_HI2040P")->setConstant(false);
  ws->var("coeffPol4_HI2040P")->setConstant(false);
  ws->var("coeffPol5_HI2040P")->setConstant(false);
  ws->var("coeffPol6_HI2040P")->setConstant(false);
  ws->var("coeffPol7_HI2040P")->setConstant(false);

  // 1st order polynomial in exponential function
  ws->factory("Polynomial::expPol1Arg_HI2040(Jpsi_Mass,{expCoeffPol1_HI2040[-0.1]})");
  ws->factory("Exponential::expPol1_HI2040(expPol1Arg_HI2040,expCoeffPol0_HI2040[1.0])");
  // 2nd order polynomial in exponential function
  ws->factory("Polynomial::expPol2Arg_HI2040(Jpsi_Mass,{expCoeffPol1_HI2040, expCoeffPol2_HI2040[0.0]})");
  ws->factory("Exponential::expPol2_HI2040(expPol2Arg_HI2040,expCoeffPol0_HI2040)");
  // 3rd order polynomial in exponential function
  ws->factory("Polynomial::expPol3Arg_HI2040(Jpsi_Mass,{expCoeffPol1_HI2040, expCoeffPol2_HI2040, expCoeffPol3_HI2040[0.0]})");
  ws->factory("Exponential::expPol3_HI2040(expPol3Arg_HI2040,expCoeffPol0_HI2040)");
  // 4th order polynomial in exponential function
  ws->factory("Polynomial::expPol4Arg_HI2040(Jpsi_Mass,{expCoeffPol1_HI2040, expCoeffPol2_HI2040, expCoeffPol3_HI2040, expCoeffPol4_HI2040[0.0]})");
  ws->factory("Exponential::expPol4_HI2040(expPol4Arg_HI2040,expCoeffPol0_HI2040)");
  // 5th order polynomial in exponential function
  ws->factory("Polynomial::expPol5Arg_HI2040(Jpsi_Mass,{expCoeffPol1_HI2040, expCoeffPol2_HI2040, expCoeffPol3_HI2040, expCoeffPol4_HI2040, expCoeffPol5_HI2040[0.0]})");
  ws->factory("Exponential::expPol5_HI2040(expPol5Arg_HI2040,expCoeffPol0_HI2040)");
  // 6th order polynomial in exponential function
  ws->factory("Polynomial::expPol6Arg_HI2040(Jpsi_Mass,{expCoeffPol1_HI2040, expCoeffPol2_HI2040, expCoeffPol3_HI2040, expCoeffPol4_HI2040, expCoeffPol5_HI2040, expCoeffPol6_HI2040[0.0]})");
  ws->factory("Exponential::expPol6_HI2040(expPol6Arg_HI2040,expCoeffPol0_HI2040)");
  // 7th order polynomial in exponential function
  ws->factory("Polynomial::expPol7Arg_HI2040(Jpsi_Mass,{expCoeffPol1_HI2040, expCoeffPol2_HI2040, expCoeffPol3_HI2040, expCoeffPol4_HI2040, expCoeffPol5_HI2040, expCoeffPol6_HI2040, expCoeffPol7_HI2040[0.0]})");
  ws->factory("Exponential::expPol7_HI2040(expPol7Arg_HI2040,expCoeffPol0_HI2040)");

  ws->var("expCoeffPol0_HI2040")->setConstant(true);
  ws->var("expCoeffPol1_HI2040")->setConstant(false);
  ws->var("expCoeffPol2_HI2040")->setConstant(false);
  ws->var("expCoeffPol3_HI2040")->setConstant(false);
  ws->var("expCoeffPol4_HI2040")->setConstant(false);
  ws->var("expCoeffPol5_HI2040")->setConstant(false);
  ws->var("expCoeffPol6_HI2040")->setConstant(false);
  ws->var("expCoeffPol7_HI2040")->setConstant(false);


  // 40-100%
  // 0th order polynomial
  ws->factory("Polynomial::pol0_HI40100(Jpsi_Mass,{coeffPol0_HI40100[0.0]})");
  // 1st order polynomial
  ws->factory("Polynomial::pol1_HI40100(Jpsi_Mass,{coeffPol1_HI40100[0.0]})");
  // 2nd order polynomial
  ws->factory("Polynomial::pol2_HI40100(Jpsi_Mass,{coeffPol1_HI40100, coeffPol2_HI40100[0.0]})");
  // 3rd order polynomial
  ws->factory("Polynomial::pol3_HI40100(Jpsi_Mass,{coeffPol1_HI40100, coeffPol2_HI40100, coeffPol3_HI40100[0.0]})");
  // 4th order polynomial
  ws->factory("Polynomial::pol4_HI40100(Jpsi_Mass,{coeffPol1_HI40100, coeffPol2_HI40100, coeffPol3_HI40100, coeffPol4_HI40100[0.0]})");
  // 5th order polynomial
  ws->factory("Polynomial::pol5_HI40100(Jpsi_Mass,{coeffPol1_HI40100, coeffPol2_HI40100, coeffPol3_HI40100, coeffPol4_HI40100, coeffPol5_HI40100[0.0]})");
  // 6th order polynomial
  ws->factory("Polynomial::pol6_HI40100(Jpsi_Mass,{coeffPol1_HI40100, coeffPol2_HI40100, coeffPol3_HI40100, coeffPol4_HI40100, coeffPol5_HI40100, coeffPol6_HI40100[0.0]})");
  // 7th order polynomial
  ws->factory("Polynomial::pol7_HI40100(Jpsi_Mass,{coeffPol1_HI40100, coeffPol2_HI40100, coeffPol3_HI40100, coeffPol4_HI40100, coeffPol5_HI40100, coeffPol6_HI40100, coeffPol7_HI40100[0.0]})");
  ws->var("coeffPol0_HI40100")->setConstant(true);
  ws->var("coeffPol1_HI40100")->setConstant(false);
  ws->var("coeffPol2_HI40100")->setConstant(false);
  ws->var("coeffPol3_HI40100")->setConstant(false);
  ws->var("coeffPol4_HI40100")->setConstant(false);
  ws->var("coeffPol5_HI40100")->setConstant(false);
  ws->var("coeffPol6_HI40100")->setConstant(false);
  ws->var("coeffPol7_HI40100")->setConstant(false);

  ws->factory("Polynomial::pol0_HI40100P(Jpsi_Mass,{coeffPol0_HI40100P[0.0]})");
  ws->factory("Polynomial::pol1_HI40100P(Jpsi_Mass,{coeffPol1_HI40100P[0.0]})");
  ws->factory("Polynomial::pol2_HI40100P(Jpsi_Mass,{coeffPol1_HI40100P, coeffPol2_HI40100P[0.0]})");
  ws->factory("Polynomial::pol3_HI40100P(Jpsi_Mass,{coeffPol1_HI40100P, coeffPol2_HI40100P, coeffPol3_HI40100P[0.0]})");
  ws->factory("Polynomial::pol4_HI40100P(Jpsi_Mass,{coeffPol1_HI40100P, coeffPol2_HI40100P, coeffPol3_HI40100P, coeffPol4_HI40100P[0.0]})");
  ws->factory("Polynomial::pol5_HI40100P(Jpsi_Mass,{coeffPol1_HI40100P, coeffPol2_HI40100P, coeffPol3_HI40100P, coeffPol4_HI40100P, coeffPol5_HI40100P[0.0]})");
  ws->factory("Polynomial::pol6_HI40100P(Jpsi_Mass,{coeffPol1_HI40100P, coeffPol2_HI40100P, coeffPol3_HI40100P, coeffPol4_HI40100P, coeffPol5_HI40100P, coeffPol6_HI40100P[0.0]})");
  ws->factory("Polynomial::pol7_HI40100P(Jpsi_Mass,{coeffPol1_HI40100P, coeffPol2_HI40100P, coeffPol3_HI40100P, coeffPol4_HI40100P, coeffPol5_HI40100P, coeffPol6_HI40100P, coeffPol7_HI40100P[0.0]})");
  ws->var("coeffPol0_HI40100P")->setConstant(true);
  ws->var("coeffPol1_HI40100P")->setConstant(false);
  ws->var("coeffPol2_HI40100P")->setConstant(false);
  ws->var("coeffPol3_HI40100P")->setConstant(false);
  ws->var("coeffPol4_HI40100P")->setConstant(false);
  ws->var("coeffPol5_HI40100P")->setConstant(false);
  ws->var("coeffPol6_HI40100P")->setConstant(false);
  ws->var("coeffPol7_HI40100P")->setConstant(false);

  // 1st order polynomial in exponential function
  ws->factory("Polynomial::expPol1Arg_HI40100(Jpsi_Mass,{expCoeffPol1_HI40100[-0.1]})");
  ws->factory("Exponential::expPol1_HI40100(expPol1Arg_HI40100,expCoeffPol0_HI40100[1.0])");
  // 2nd order polynomial in exponential function
  ws->factory("Polynomial::expPol2Arg_HI40100(Jpsi_Mass,{expCoeffPol1_HI40100, expCoeffPol2_HI40100[0.0]})");
  ws->factory("Exponential::expPol2_HI40100(expPol2Arg_HI40100,expCoeffPol0_HI40100)");
  // 3rd order polynomial in exponential function
  ws->factory("Polynomial::expPol3Arg_HI40100(Jpsi_Mass,{expCoeffPol1_HI40100, expCoeffPol2_HI40100, expCoeffPol3_HI40100[0.0]})");
  ws->factory("Exponential::expPol3_HI40100(expPol3Arg_HI40100,expCoeffPol0_HI40100)");
  // 4th order polynomial in exponential function
  ws->factory("Polynomial::expPol4Arg_HI40100(Jpsi_Mass,{expCoeffPol1_HI40100, expCoeffPol2_HI40100, expCoeffPol3_HI40100, expCoeffPol4_HI40100[0.0]})");
  ws->factory("Exponential::expPol4_HI40100(expPol4Arg_HI40100,expCoeffPol0_HI40100)");
  // 5th order polynomial in exponential function
  ws->factory("Polynomial::expPol5Arg_HI40100(Jpsi_Mass,{expCoeffPol1_HI40100, expCoeffPol2_HI40100, expCoeffPol3_HI40100, expCoeffPol4_HI40100, expCoeffPol5_HI40100[0.0]})");
  ws->factory("Exponential::expPol5_HI40100(expPol5Arg_HI40100,expCoeffPol0_HI40100)");
  // 6th order polynomial in exponential function
  ws->factory("Polynomial::expPol6Arg_HI40100(Jpsi_Mass,{expCoeffPol1_HI40100, expCoeffPol2_HI40100, expCoeffPol3_HI40100, expCoeffPol4_HI40100, expCoeffPol5_HI40100, expCoeffPol6_HI40100[0.0]})");
  ws->factory("Exponential::expPol6_HI40100(expPol6Arg_HI40100,expCoeffPol0_HI40100)");
  // 7th order polynomial in exponential function
  ws->factory("Polynomial::expPol7Arg_HI40100(Jpsi_Mass,{expCoeffPol1_HI40100, expCoeffPol2_HI40100, expCoeffPol3_HI40100, expCoeffPol4_HI40100, expCoeffPol5_HI40100, expCoeffPol6_HI40100, expCoeffPol7_HI40100[0.0]})");
  ws->factory("Exponential::expPol7_HI40100(expPol7Arg_HI40100,expCoeffPol0_HI40100)");

  ws->var("expCoeffPol0_HI40100")->setConstant(true);
  ws->var("expCoeffPol1_HI40100")->setConstant(false);
  ws->var("expCoeffPol2_HI40100")->setConstant(false);
  ws->var("expCoeffPol3_HI40100")->setConstant(false);
  ws->var("expCoeffPol4_HI40100")->setConstant(false);
  ws->var("expCoeffPol5_HI40100")->setConstant(false);
  ws->var("expCoeffPol6_HI40100")->setConstant(false);
  ws->var("expCoeffPol7_HI40100")->setConstant(false);

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
  // wide CB
  ws->factory("CBShape::signalCB2_HI(Jpsi_Mass,meanSig1_HI,sigmaSig2_HI,alpha_HI,enne_HI)");

  // CB(narrow) + Gauss(wide)
  ws->factory("SUM::sigCB1G2_HI(coeffGaus_HI[0.1,0.0,1.0]*signalG2_HI,signalCB1_HI)");
  // CB(narrow) + CB(wide)
  ws->factory("SUM::sigCB1CB2_HI(coeffGaus_HI[0.1,0.0,1.0]*signalCB2_HI,signalCB1_HI)");

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
  // narrow CB
  ws->factory("CBShape::signalCB2P_HI(Jpsi_Mass,meanSig1P_HI,sigmaSig2P_HI,alpha_HI,enne_HI)");

  // CB(narrow) + Gauss(wide)
  ws->factory("SUM::sigCB1G2P_HI(coeffGaus_HI*signalG2P_HI,signalCB1P_HI)");
  // CB(narrow) + CB(wide)
  ws->factory("SUM::sigCB1CB2P_HI(coeffGaus_HI*signalCB2P_HI,signalCB1P_HI)");

  return;
}
