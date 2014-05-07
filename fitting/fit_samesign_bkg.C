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
#include "RooPolynomial.h"
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

#include <RooWorkspace.h>

#include "TCut.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLatex.h"
#include "TFile.h"
#include "TAxis.h"
#include "TH1.h"
#include "TF1.h"

using namespace std;
using namespace RooFit;

void getOptRange(string &ran, float *min, float *max);

void fit_samesign_bkg(bool isHI=true, string prange="6.5-30.0", string yrange="0.0-1.6")
{
  string infname;
  if (isHI)
    infname = "../root_files/PbPbData2011_DblMu0_cent0-20_M18-50.root";
  else
    infname = "../root_files/ppData2013_DblMu0_cent0-100_M18-50.root";
  TFile *inf = new TFile(infname.c_str(),"READ");

  float pmin=0, pmax=0, ymin=0, ymax=0;
  getOptRange(prange,&pmin,&pmax);
  getOptRange(yrange,&ymin,&ymax);

  // Create workspace to play with
  RooWorkspace *ws = new RooWorkspace("workspace");

  RooRealVar Jpsi_Mass("Jpsi_Mass","J/#psi mass",1.8,5.0,"GeV/c^{2}");
  RooRealVar Jpsi_Pt("Jpsi_Pt","J/#psi pt",0,30,"GeV/c");
  RooRealVar Jpsi_Y("Jpsi_Y","J/#psi y",-2.4,2.4);
  RooRealVar Jpsi_Ct("Jpsi_Ct","J/#psi c#tau",-3.0,3.5,"mm");

  Jpsi_Mass.setBins(160);
  Jpsi_Pt.setBins(60);
  Jpsi_Y.setBins(1);
  Jpsi_Ct.setBins(1);

  Jpsi_Mass.setRange("signal",3.6,3.76);
  // Jpsi_Mass.setRange("M2045",2.0,4.5);
  // Jpsi_Mass.setRange("M2242",2.2,4.2);
  // cout << "Min(full) = " << Jpsi_Mass.min() << endl;
  // cout << "Max(full) = " << Jpsi_Mass.max() << endl;
  // cout << "Min(M2045) = " << Jpsi_Mass.min("M2045") << endl;
  // cout << "Max(M2045) = " << Jpsi_Mass.max("M2045") << endl;

  RooRealVar a0("a0","a0",0.0);a0.setConstant(false);
  RooRealVar a1("a1","a1",0.0);a1.setConstant(false);
  RooRealVar a2("a2","a2",0.0);a2.setConstant(false);
  RooRealVar b0("b0","b0",0.0);b0.setConstant(false);
  RooRealVar b1("b1","b1",0.0);b1.setConstant(false);
  RooRealVar b2("b2","b2",0.0);b2.setConstant(false);
  RooRealVar c0("c0","c0",0.0);c0.setConstant(false);
  RooRealVar c1("c1","c1",0.0);c1.setConstant(false);
  RooRealVar c2("c2","c2",0.0);c2.setConstant(false);
  RooChebychev bkg1("bkg1","Background pol3",Jpsi_Mass,RooArgSet(a0,a1,a2));
  RooChebychev bkg2("bkg2","Background pol3",Jpsi_Mass,RooArgSet(b0,b1,b2));
  RooChebychev bkg3("bkg3","Background pol3",Jpsi_Mass,RooArgSet(c0,c1,c2));
  // RooChebychev bkg2("bkg2","Background pol2",Jpsi_Mass,RooArgSet(a0,a1));
  // RooChebychev bkg3("bkg3","Background pol3",Jpsi_Mass,RooArgSet(a0,a1,a2));
  // RooChebychev bkg4("bkg4","Background pol4",Jpsi_Mass,RooArgSet(a0,a1,a2,a3));
  // RooChebychev bkg5("bkg5","Background pol5",Jpsi_Mass,RooArgSet(a0,a1,a2,a3,a4));

  RooRealVar zero_A("zero_A","zero_A",0.0);zero_A.setConstant(true);
  RooRealVar zero_B("zero_B","zero_B",0.0);zero_B.setConstant(true);
  RooRealVar zero_C("zero_C","zero_C",0.0);zero_C.setConstant(true);
  RooPolynomial dummy("dummy","dummy pdf",Jpsi_Mass,RooArgSet(zero_A));

  RooRealVar Ndummy("Ndummy","Ndummy",0);  Ndummy.setConstant(true);
  RooRealVar Nbkg_A("Nbkg_A","Nbkg_A",100,0,1e6);
  RooRealVar Nbkg_B("Nbkg_B","Nbkg_B",100,0,1e6);
  RooRealVar Nbkg_C("Nbkg_C","Nbkg_C",100,0,1e6);

  RooAddPdf model_A("model_A","model_A",RooArgList(bkg1,dummy),RooArgList(Nbkg_A,Ndummy));
  RooAddPdf model_B("model_B","model_B",RooArgList(bkg2,zero_B),RooArgList(Nbkg_B,Ndummy));
  RooAddPdf model_C("model_C","model_C",RooArgList(bkg3,zero_C),RooArgList(Nbkg_C,Ndummy));

  ws->import(model_A);
  ws->Print("v");

  RooDataSet *data = (RooDataSet*)inf->Get("dataJpsiSame");
  char reduceDS[300];
  if (isHI) {
    if (prange=="6.5-30.0" && yrange=="0.0-1.6")
      sprintf(reduceDS,"Jpsi_Pt>%.2f && Jpsi_Pt<%.2f && abs(Jpsi_Y)>%.2f && abs(Jpsi_Y)<%.2f && Jpsi_Ct<0.04",pmin,pmax,ymin,ymax);
    else if (prange=="3.0-30.0" && yrange=="1.6-2.4")
      sprintf(reduceDS,"Jpsi_Pt>%.2f && Jpsi_Pt<%.2f && abs(Jpsi_Y)>%.2f && abs(Jpsi_Y)<%.2f && Jpsi_Ct<0.09",pmin,pmax,ymin,ymax);
  }
  else {
    if (prange=="6.5-30.0" && yrange=="0.0-1.6")
      sprintf(reduceDS,"Jpsi_Pt>%.2f && Jpsi_Pt<%.2f && abs(Jpsi_Y)>%.2f && abs(Jpsi_Y)<%.2f && Jpsi_Ct<0.04",pmin,pmax,ymin,ymax);
    else if (prange=="3.0-30.0" && yrange=="1.6-2.4")
      sprintf(reduceDS,"Jpsi_Pt>%.2f && Jpsi_Pt<%.2f && abs(Jpsi_Y)>%.2f && abs(Jpsi_Y)<%.2f && Jpsi_Ct<0.09",pmin,pmax,ymin,ymax);
  }
  RooDataSet *tmp = (RooDataSet*)data->reduce(reduceDS);

  RooArgList list(Jpsi_Mass);
  RooDataSet *redData = (RooDataSet*) tmp->reduce(list);
  redData->SetName("redData");
  redData->SetTitle("redData");

  ws->import(*redData);

  // RooDataSet *redDataM2045 = (RooDataSet*) tmp->reduce("Jpsi_Mass>2.0&&Jpsi_Mass<4.5");
  // redDataM2045->SetName("redDataM2045");
  // redDataM2045->SetTitle("redDataM2045");

  // RooDataSet *redDataM2242 = (RooDataSet*) tmp->reduce("Jpsi_Mass>2.2&&Jpsi_Mass<4.2");
  // redDataM2242->SetName("redDataM2242");
  // redDataM2242->SetTitle("redDataM2242");

  
  RooFitResult *fitM;
  fitM = ws->pdf("model_A")->fitTo(*redData,Extended(1),Minos(1),Save(1),SumW2Error(kFALSE));//,NumCPU(4),Range("full"));
  double relYield_all = ws->pdf("model_A")->analyticalIntegral(1,"signal")/ws->pdf("model_A")->analyticalIntegral(1,"full");
  fitM->Print();


  /*
  Jpsi_Mass.setMin(2.0);
  Jpsi_Mass.setMax(4.5);
  b0.setVal(a0.getVal());
  b1.setVal(a1.getVal());
  b2.setVal(a2.getVal());
  fitM = model_B.fitTo(*redDataM2045,Extended(1),Minos(1),Save(1),SumW2Error(kFALSE));//,NumCPU(4),Range("full"));
  double relYield_M2045 = model_B.analyticalIntegral(1,"signal")/model_B.analyticalIntegral(1,"full");


  Jpsi_Mass.setMin(2.2);
  Jpsi_Mass.setMax(4.2);
  c0.setVal(b0.getVal());
  c1.setVal(b1.getVal());
  c2.setVal(b2.getVal());
  fitM = model_C.fitTo(*redDataM2242,Extended(1),Minos(1),Save(1),SumW2Error(kFALSE),NumCPU(4));//,Range("full"));
  
  double relYield_M2242 = model_C.analyticalIntegral(1,"signal")/model_C.analyticalIntegral(1,"full");


  Jpsi_Mass.setMin(1.8);
  Jpsi_Mass.setMax(5);
  */
  RooPlot* mframe = Jpsi_Mass.frame(Title("data"));
  ws->data("redData")->plotOn(mframe);
  ws->pdf("model_A")->plotOn(mframe, LineColor(kBlue));//,Range("full"),NormRange("full"));
  // model_B.plotOn(mframe, LineColor(kRed),Range("M2045"),NormRange("M2045"));
  // model_C.plotOn(mframe, LineColor(kGreen+2),Range("M2242"),NormRange("M2242"));

  cout << "Mass range 1.8-5.0 GeV: " << relYield_all << endl;
  // cout << "Mass range 2.0-4.5 GeV: " << relYield_M2045 << endl;
  // cout << "Mass range 2.2-4.2 GeV: " << relYield_M2242 << endl;
  
  char *cut = NULL;
  cout << "Actual fraction in data" << endl;
  cout << ws->data("redData")->sumEntries(cut,"signal")/ws->data("redData")->sumEntries(cut,"full") <<endl;
  // cout << redData->sumEntries(cut,"signal")/redData->sumEntries(cut,"M2045") <<endl;
  // cout << redData->sumEntries(cut,"signal")/redData->sumEntries(cut,"M2242") <<endl;

  TCanvas *canv = new TCanvas("canv","canv");
  canv->SetLogy(0);
  mframe->Draw();

  // TH1 *h1 = (TH1*) redData->createHistogram("Jpsi_Mass")->Clone("h1");
  // TFile *outf = new TFile("fwd_pbpb_cent0-20.root","RECREATE");
  // h1->Write();
  // outf->Close();

  ws->writeToFile("fwd_pbpb_cent0-20_workspace.root",kTRUE);

  return;
}

void getOptRange(string &ran, float *min, float *max) {
  if (sscanf(ran.c_str(), "%f-%f", min, max) == 0) {
    cout << ran.c_str() << ": not valid!" << endl;
    assert(0);
  }
  return ;
}
