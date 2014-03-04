#include <iostream>

//#ifndef __CINT__
#include "RooGlobalFunc.h"
//#endif

#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooConstVar.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooFormulaVar.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooBinning.h"

#include "RooGenericPdf.h"
#include "RooChi2Var.h"
#include "RooMinuit.h"

#include "TCut.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLatex.h"
#include "TFile.h"
#include "TAxis.h"
#include "TH1.h"
#include "TF1.h"

using namespace RooFit;

void fitSideBands(bool isHI=false, double ptmin=0.0, double ptmax=30.0, double ymin=0.0, double ymax=2.4, bool absRapidity=true, bool savePlot=false)
{
  TCut ptCut = Form("Jpsi_Pt>%3.1f&&Jpsi_Pt<%3.1f",ptmin,ptmax);
  TCut rapCut;
  if (absRapidity)
    rapCut = Form("abs(Jpsi_Y)>%3.1f&&abs(Jpsi_Y)<%3.1f",ymin,ymax);
  else
    rapCut = Form("Jpsi_Y>%3.1f&&Jpsi_Y<%3.1f",ymin,ymax);

  RooRealVar Jpsi_Mass("Jpsi_Mass","m_{#mu#mu}",2.2,3.4,"GeV/c^{2}");
  RooRealVar Jpsi_Pt("Jpsi_Pt","J/#psi pt",0,30,"GeV/c");
  RooRealVar Jpsi_Y("Jpsi_Y","J/#psi y",-2.4,2.4);
  RooRealVar Jpsi_Ct("Jpsi_Ct","J/#psi c#tau",-3.0,3.5,"mm");
  RooRealVar Jpsi_CtErr("Jpsi_CtErr","J/#psi c#tau error",-1.,1.,"mm");

  Jpsi_Mass.setRange("SBLow",2.2,2.8);
  Jpsi_Mass.setRange("SBMid",3.3,3.5);
  Jpsi_Mass.setRange("SBHigh",3.85,4.2);

  //  Jpsi_Mass.setBins(100); for 2.2--4.2 GeV
  RooBinning rbm(2.2,4.2);
  int nbins = 100;
  bool highStats=true;
  if (highStats)
    rbm.addUniform(nbins,2.2,4.2);
  else {
    rbm.addUniform(11,2.2,2.860); // 60 MeV
    rbm.addUniform(19,2.860,3.240); // 20 MeV
    rbm.addUniform(6,3.240,3.600); // 60 MeV
    rbm.addUniform(4,3.600,3.720); // 30 MeV
    rbm.addUniform(8,3.720,4.2); // 60 MeV
  }

  Jpsi_Mass.setBinning(rbm);
  // Jpsi_Pt.setBins(60);
  // Jpsi_Y.setBins(1);
  // Jpsi_Ct.setBins(1);
  // Jpsi_CtErr.setBins(1);

  RooRealVar a0("a0","a0",0.0,-1.0,1.0);
  RooRealVar a1("a1","a1",0.0,-1.0,1.0);
  RooRealVar a2("a2","a2",0.0,-1.0,1.0);
  RooRealVar a3("a3","a3",0.0,-1.0,1.0);
  RooRealVar a4("a4","a4",0.0,-1.0,1.0);
  RooRealVar a5("a5","a5",0.0,-1.0,1.0);
  RooChebychev pol1("pol1","Background (pol1)",Jpsi_Mass,RooArgSet(a0));
  RooChebychev pol2("pol2","Background (pol2)",Jpsi_Mass,RooArgSet(a0,a1));
  RooChebychev pol3("pol3","Background (pol3)",Jpsi_Mass,RooArgSet(a0,a1,a2));
  RooChebychev pol4("pol4","Background (pol4)",Jpsi_Mass,RooArgSet(a0,a1,a2,a3));
  RooChebychev pol5("pol5","Background (pol5)",Jpsi_Mass,RooArgSet(a0,a1,a2,a3,a4));
  RooChebychev pol6("pol6","Background (pol6)",Jpsi_Mass,RooArgSet(a0,a1,a2,a3,a4,a5));
  RooRealVar nbkg("nbkg","nbkg",1000,0,100000);
  RooAddPdf model("model","model",RooArgList(pol6),RooArgList(nbkg));

  std::string fname;
  if (isHI)
    fname = "../root_files/PbPbData2011_DblMu0_cent0-100_M22-42.root";
  else
    fname = "../root_files/ppData2013_DblMu0_cent0-100_M22-42.root";

  TFile *inf = new TFile(fname.c_str());

  RooDataSet *data = (RooDataSet*)inf->Get("dataJpsi");
  data->SetName("data");

  RooDataSet *redData;
  RooArgList list(Jpsi_Mass);
  redData = (RooDataSet*)data->reduce(list, (ptCut+rapCut).GetTitle());
  redData->Print("v");

  RooFitResult *fitM1;
  RooFitResult *fitM2;
  RooFitResult *fitM3;
  RooFitResult *fitM4;
  RooFitResult *fitM5;
  RooFitResult *fitM6;
  //  fitM1 = pol1.fitTo(*redData,Extended(0),Minos(0),Save(1),SumW2Error(kFALSE),NumCPU(4),Range("SBLow,SBMid,SBHigh"));
  //  fitM2 = pol2.fitTo(*redData,Extended(0),Minos(0),Save(1),SumW2Error(kFALSE),NumCPU(4),Range("SBLow,SBMid,SBHigh"));
  fitM3 = pol3.fitTo(*redData,Extended(0),Minos(0),Save(1),SumW2Error(kFALSE),NumCPU(4),Range("SBLow"));
  // fitM4 = pol4.fitTo(*redData,Extended(0),Minos(0),Save(1),SumW2Error(kFALSE),NumCPU(4),Range("SBLow,SBMid,SBHigh"));
  // fitM5 = pol5.fitTo(*redData,Extended(0),Minos(0),Save(1),SumW2Error(kFALSE),NumCPU(4),Range("SBLow,SBMid,SBHigh"));
  //  fitM6 = pol6.fitTo(*redData,Extended(0),Minos(0),Save(1),SumW2Error(kFALSE),NumCPU(4),Range("SBLow"));
  //  fitM6->printMultiline(std::cout,1,1,"");


  // RooDataHist* binnedData = redData->binnedClone("binnedData","binnedData");
  // RooChi2Var chi2("chi2","chi2",pol3,*binnedData,Range("SBLow,SBMid,SBHigh"));
  // RooMinuit m(chi2);
  // m.migrad();
  // m.hesse();
  // m.hesse();
  // //    m.minos();
  // fitM3 = m.save();


  RooPlot* xframe = Jpsi_Mass.frame(Title("data"));
  redData->plotOn(xframe,RooFit::Binning(rbm));

  //  model_g.plotOn(xframe,LineWidth(1),LineColor(kBlack));
  //  model_cb.plotOn(xframe,LineWidth(1), LineColor(kBlue));

  // pol1.plotOn(xframe, LineColor(kRed));
  pol3.plotOn(xframe, LineColor(kRed),Range("SBLow"));
  pol3.plotOn(xframe, LineColor(kBlue),LineStyle(kDashed),Range("Full"));
  // pol3.plotOn(xframe, LineColor(kGreen+2), LineWidth(2));
  // pol4.plotOn(xframe, LineColor(kOrange+2), LineWidth(2));
  // pol5.plotOn(xframe, LineColor(kMagenta-1), LineStyle(kDashed));
  //pol6.plotOn(xframe, LineColor(kCyan+2), LineStyle(kDotted),Range("Full"));

  //  redData->plotOn(xframe,RooFit::Binning(rbm));

  xframe->SetMinimum(90.0);
  xframe->SetMaximum(1e4);


  TCanvas *c1 = new TCanvas("c1","c1");
  c1->SetLogy();
  xframe->Draw();

  TLatex *lcoll;
  if (isHI)
    lcoll = new TLatex(0.15,0.9,"CMS PbPb #sqrt{s_{NN}} = 2.76 TeV");
  else
    lcoll = new TLatex(0.15,0.9,"CMS pp #sqrt{s} = 2.76 TeV");

  lcoll->SetNDC(kTRUE);
  lcoll->Draw();

  TLatex *lpt;
  if (ptmin==0.0)
    lpt = new TLatex(0.15,0.84,Form("p_{T} < %3.1f GeV/c",ptmax));
  else
    lpt = new TLatex(0.15,0.84,Form("%3.1f < p_{T} < %3.1f GeV/c",ptmin,ptmax));
  TLatex *lrap;
  if (absRapidity){
    if (ymin==0.0)
      lrap = new TLatex(0.15,0.78,Form("|y| < %3.1f",ymax));
    else
      lrap = new TLatex(0.15,0.78,Form("%3.1f < |y| < %3.1f",ymin,ymax));
  }
  else {
    if (ymin==0.0)
      lrap = new TLatex(0.15,0.78,Form("y < %3.1f",ymax));
    else
      lrap = new TLatex(0.15,0.78,Form("%3.1f < y < %3.1f",ymin,ymax));
  }

  lpt->SetNDC(kTRUE);
  lrap->SetNDC(kTRUE);
  lpt->Draw();
  lrap->Draw();


  TString outfname;
  if (isHI)
    outfname = Form("PbPb_Mass_Sidebands_Rap_%3.1f-%3.1f_Pt_%3.1f-%3.1f.pdf",ymin,ymax,ptmin,ptmax);
  else
    outfname = Form("pp_Mass_SideBands_Rap_%3.1f-%3.1f_Pt_%3.1f-%3.1f.pdf",ymin,ymax,ptmin,ptmax);

  if (savePlot) {
    std::cout << outfname << std::endl;
    c1->SaveAs(outfname);
  }

  return;
}
