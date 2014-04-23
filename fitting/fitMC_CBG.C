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

void fitMC_CBG(bool isHI=false, double ptmin=0.0, double ptmax=30.0, double ymin=0.0, double ymax=2.4, bool absRapidity=true, bool fitChi2=false, bool savePlot=false)
{
  TCut ptCut = Form("Jpsi_Pt>%3.1f&&Jpsi_Pt<%3.1f",ptmin,ptmax);
  TCut rapCut;
  if (absRapidity)
    rapCut = Form("abs(Jpsi_Y)>%3.1f&&abs(Jpsi_Y)<%3.1f",ymin,ymax);
  else
    rapCut = Form("Jpsi_Y>%3.1f&&Jpsi_Y<%3.1f",ymin,ymax);

  RooRealVar Jpsi_Mass("Jpsi_Mass","J/#psi mass",2.2,3.4,"GeV/c^{2}");
  RooRealVar Jpsi_Pt("Jpsi_Pt","J/#psi pt",0,30,"GeV/c");
  RooRealVar Jpsi_Y("Jpsi_Y","J/#psi y",-2.4,2.4);
  RooRealVar Jpsi_Ct("Jpsi_Ct","J/#psi c#tau",-3.0,3.5,"mm");
  RooRealVar Jpsi_CtErr("Jpsi_CtErr","J/#psi c#tau error",-1.,1.,"mm");
  RooRealVar Jpsi_CtTrue("Jpsi_CtTrue","J/#psi c#tau true",-3.0,3.5,"mm");
  RooRealVar Gen_Pt("Gen_Pt","Generated J/#psi pt",0,30,"GeV/c");

  //  Jpsi_Mass.setBins(100); for 2.2--4.2 GeV
  RooBinning bins(2.2,3.4);
  //  bins.addUniform(14,2.2,2.9);
  bins.addUniform(12,2.2,2.8);
  //  bins.addUniform(20,2.9,3.3);
  bins.addUniform(25,2.8,3.3);
  bins.addUniform(1,3.3,3.4);
  Jpsi_Mass.setBinning(bins);
  //  Jpsi_Mass.setBins(35);
  Jpsi_Pt.setBins(60);
  Jpsi_Y.setBins(1);
  Jpsi_Ct.setBins(1);
  Jpsi_CtErr.setBins(1);
  Jpsi_CtTrue.setBins(1);
  Gen_Pt.setBins(60);

  RooRealVar mean("mean","mean", 3.0969,3.05,3.15,"GeV/c^{2}");
  RooRealVar sigma1("#sigma_{CB}","#sigma_{CB}", 0.03,0.005,0.080,"GeV/c^{2}");
  RooGaussian gaussS("gaussS","gaussS",Jpsi_Mass,mean,sigma1);

  RooRealVar alpha("#alpha","#alpha", 1.0,0.0,3.0);
  RooRealVar n("n","n", 5,1.0,50.0);
  RooCBShape cballS("cballS","cballS",Jpsi_Mass,mean,sigma1,alpha,n);
  RooCBShape cball("cball","cball",Jpsi_Mass,mean,sigma1,alpha,n);

  RooRealVar wideFactor("n_{G}","n_{G}",2.0,1.0,6.0);
  RooFormulaVar sigma2("#sigma_{G}","@0*@1",RooArgList(sigma1,wideFactor));
  RooGaussian gauss2("gauss","gauss",Jpsi_Mass,mean,sigma2);

  RooRealVar coeffGauss("f_{G}","f_{G}",0.1,0.0,1.0);
  RooAddPdf cbg("cbg","cbg",RooArgList(gauss2,cball),coeffGauss);

  RooRealVar a0("a0","a0",0.0,-1.0,1.0);
  RooRealVar a1("a1","a1",0.0,-1.0,1.0);
  RooChebychev bkg("bkg","Background",Jpsi_Mass,RooArgSet(a0));

  RooRealVar Nsignal("Nsignal","Nsginal",10000);
  RooRealVar Nbackground("Nbackground","Nbackground",100,0,1e6);
  Nsignal.setConstant(0);
  //  if (!isHI) {
  // Nbackground.setVal(0);
  // Nbackground.setConstant(1);
  //  }

  RooRealVar bkgfrac("bkgfrac","fraction of background",0.001,0.0,1.0);
  // RooAddPdf model_g("model_g","cb+bkg",RooArgList(bkg,gaussS),bkgfrac);
  // RooAddPdf model_cb("model_cb","cb+bkg",RooArgList(bkg,cballS),bkgfrac);
  // RooAddPdf model_cbg("model_cbg","cb+bkg",RooArgList(bkg,cbg),bkgfrac);
  // Extended definitions
  RooAddPdf model_g("model_g","cb+bkg",RooArgList(bkg,gaussS),RooArgList(Nbackground,Nsignal));
  RooAddPdf model_cb("model_cb","cb+bkg",RooArgList(bkg,cballS),RooArgList(Nbackground,Nsignal));
  RooAddPdf model_cbg("model_cbg","cb+bkg",RooArgList(bkg,cbg),RooArgList(Nbackground,Nsignal));
  RooAddPdf model_cbg2("model_cbg2","cb+bkg2",RooArgList(bkg,cbg),RooArgList(Nbackground,Nsignal));
  a0.setVal(0);
  //  if (!isHI)
  a0.setConstant(1);
  bkgfrac.setVal(0);
  //  bkgfrac.setConstant(1);

  if (isHI && false) {
    std::cout << ptmin << " " << ptmax << " " << ymin << " " << ymax << std::endl;
    if (fabs(ptmin-6.5)/ptmin<1e-5 && fabs(ptmax-30)/ptmax<1e-5 && fabs(ymin-0)<1e-5 && fabs(ymax-2.4)/ymax<1e-5) {
      std::cout << "GOOD" << std::endl;
      alpha.setVal(1.818);
      n.setVal(1.588);
      wideFactor.setVal(1.951);
    }
    else if (fabs(ptmin-6.5)/ptmin<1e-5 && fabs(ptmax-30)/ptmax<1e-5 && fabs(ymin-0)<1e-5 && fabs(ymax-1.6)/ymax<1e-5) {
      std::cout << "GOOD" << std::endl;
      alpha.setVal(1.735);
      n.setVal(1.628);
      wideFactor.setVal(1.734);
    }
    else if (fabs(ptmin-3)/ptmin<1e-5 && fabs(ptmax-30)/ptmax<1e-5 && fabs(ymin-1.6)/ymin<1e-5 && fabs(ymax-2.4)/ymax<1e-5) {
      std::cout << "GOOD" << std::endl;
      alpha.setVal(2.136);
      n.setVal(1.320);
      //      alpha.setVal(2.120);
      //      n.setVal(1.347);
      wideFactor.setVal(1.681);
    }
    else if (fabs(ptmin-3)/ptmin<1e-5 && fabs(ptmax-6.5)/ptmax<1e-5 && fabs(ymin-1.6)/ymin<1e-5 && fabs(ymax-2.4)/ymax<1e-5) {
      std::cout << "GOOD" << std::endl;
      alpha.setVal(2.091);
      n.setVal(1.270);
      wideFactor.setVal(1.546);
    }
    else if (fabs(ptmin-6.5)/ptmin<1e-5 && fabs(ptmax-30)/ptmax<1e-5 && fabs(ymin-1.6)/ymin<1e-5 && fabs(ymax-2.4)/ymax<1e-5) {
      std::cout << "GOOD" << std::endl;
      alpha.setVal(2.161);
      n.setVal(1.370);
      wideFactor.setVal(1.970);
    }
    else {
      std::cout << "BAD" << std::endl;
      alpha.setVal(2.0);
      n.setVal(1.4);
    }
    // alpha.setConstant(kTRUE);
    // n.setConstant(kTRUE);
    // wideFactor.setConstant(kTRUE);
  }

  std::string fname;
  if (isHI)
    fname = "../root_files/PbPbPromptJpsiMC_DblMu0_cent0-100_M2234.root";
  else
    fname = "../root_files/ppPromptJpsiMC_DblMu0_cent0-100_M2234.root";

  TFile *inf = new TFile(fname.c_str());

  RooDataSet *data = (RooDataSet*)inf->Get("dataJpsi");
  data->SetName("data");
  RooDataSet *redData;
  RooDataSet *tmp;
  if (isHI) {
    RooFormulaVar wFunc("w","event weight","Gen_Pt>30?0:(Gen_Pt>15.0?(0.00579/166960.0):(Gen_Pt>12?(0.00979/143944.0):(Gen_Pt>9?(0.0379/167336.0):(Gen_Pt>6?(0.178/172176.0):(Gen_Pt>3?(0.915/124802.0):(1.0/117512.0))))))",Gen_Pt);
    RooRealVar *w = (RooRealVar*) data->addColumn(wFunc);
    tmp = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),(ptCut+rapCut).GetTitle(),w->GetName());
    RooArgList list(Jpsi_Mass);
    redData = (RooDataSet*) tmp->reduce(list);
    redData->SetName("redData");
    redData->SetTitle("redData");
  }
  else {
    RooArgList list(Jpsi_Mass);
    redData = (RooDataSet*)data->reduce(list, (ptCut+rapCut).GetTitle());
  }

  redData->Print("v");
  RooDataHist* binnedData = redData->binnedClone("binnedData","binnedData");
  binnedData->Print("v");

  // if (isHI) {
  //   model_g.fitTo(*binnedData,Extended(1),Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(4));
  //   //   model_cb.fitTo(*binnedData,Extended(1),Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(4));
  // }
  // else {
  //   model_g.fitTo(*redData,Extended(1),Minos(0),Save(1),SumW2Error(kFALSE),NumCPU(4));
  //   //   model_cb.fitTo(*redData,Extended(1),Minos(0),Save(1),SumW2Error(kFALSE),NumCPU(4));
  // }

  TString *CBalpha = new TString();
  CBalpha = alpha.format(1,"LNEA");
  TString *CBn = new TString();
  CBn = n.format(1,"LNEA");
  TString *CBsigma = new TString();
  CBsigma = sigma1.format(1,"LNEAU");
  RooFitResult *fitM;

  if (fitChi2) {
    // chi2 fit
    RooChi2Var chi2("chi2","chi2",model_cbg,*binnedData);
    RooMinuit m(chi2);
    m.migrad();
    m.hesse();
    //    m.hesse();
    //    m.minos();
    fitM = m.save();
  }
  else {
    if (isHI)
      fitM = model_cbg.fitTo(*binnedData,Extended(1),Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(4));
    else
      fitM = model_cbg.fitTo(*redData,Extended(1),Minos(0),Save(1),SumW2Error(kFALSE),NumCPU(4));
  }

  fitM->printMultiline(std::cout,1,1,"");
  TString *CBGalpha = new TString();
  CBGalpha = alpha.format(2,"LNEA");
  TString *CBGn = new TString();
  CBGn = n.format(2,"LNEA");
  TString *CBGsigma = new TString();
  CBGsigma = sigma1.format(1,"LNEAU");
  TString *CBGwideFactor = new TString();
  CBGwideFactor = wideFactor.format(1,"LNEA");
  TString *CBGcoeffGauss = new TString();
  CBGcoeffGauss = coeffGauss.format(2,"LNEA");

  RooPlot* xframe = Jpsi_Mass.frame(Title("data"));
  if (isHI)
    binnedData->plotOn(xframe,RooFit::Binning(bins));
  else
    redData->plotOn(xframe,RooFit::Binning(bins));

  //  model_g.plotOn(xframe,LineWidth(1),LineColor(kBlack));
  //  model_cb.plotOn(xframe,LineWidth(1), LineColor(kBlue));

  model_cbg.plotOn(xframe, LineColor(kBlack));
  model_cbg.plotOn(xframe, Components(RooArgSet(gauss2)), LineColor(kGreen+2), LineStyle(kDashed));
  model_cbg.plotOn(xframe, Components(RooArgSet(cball)), LineColor(kRed), LineStyle(kDashed));
  model_cbg.plotOn(xframe, Components(RooArgSet(bkg)), LineColor(kBlue), LineStyle(kDashed));


  if (isHI)
    binnedData->plotOn(xframe,RooFit::Binning(bins));
  else
    redData->plotOn(xframe,RooFit::Binning(bins));

  if (isHI)
    xframe->SetMinimum(1e-6);
  else
    xframe->SetMinimum(0.1);

  xframe->SetMaximum(20.0*xframe->GetMaximum());

  //  cbg.Print("t");
  std::cout << "CB only:\t" << *CBalpha << "\t" << *CBn << std::endl;
  std::cout << "CB + Gauss:\t" << *CBGalpha << "\t" << *CBGn << std::endl;

  TCanvas *c1 = new TCanvas("c1","c1");
  //  if (!isHI)
  c1->SetLogy();
  xframe->Draw();

  TLatex *lcoll;
  if (isHI)
    lcoll = new TLatex(0.15,0.9,"PYTHIA+HYDJET: PbPb #sqrt{s_{NN}} = 2.76 TeV");
  else
    lcoll = new TLatex(0.15,0.9,"PYTHIA: pp #sqrt{s} = 2.76 TeV");

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
      lrap = new TLatex(0.15,0.79,Form("|y| < %3.1f",ymax));
    else
      lrap = new TLatex(0.15,0.79,Form("%3.1f < |y| < %3.1f",ymin,ymax));
  }
  else {
    if (ymin==0.0)
      lrap = new TLatex(0.15,0.79,Form("y < %3.1f",ymax));
    else
      lrap = new TLatex(0.15,0.79,Form("%3.1f < y < %3.1f",ymin,ymax));
  }

  lpt->SetNDC(kTRUE);
  lrap->SetNDC(kTRUE);
  lpt->Draw();
  lrap->Draw();

  TLatex *lCB = new TLatex(0.15,0.72,"CB:");
  TLatex *lCBalpha = new TLatex(0.15,0.68,CBalpha->Data());
  TLatex *lCBn = new TLatex(0.15,0.64,CBn->Data());
  TLatex *lCBsigma = new TLatex(0.15,0.60,CBsigma->Data());
  lCB->SetTextSize(0.035);
  lCBalpha->SetTextSize(0.035);
  lCBn->SetTextSize(0.035);
  lCBsigma->SetTextSize(0.035);

  // TLatex *lCBG = new TLatex(0.15,0.54,"CB + Gauss:");
  // TLatex *lCBGalpha = new TLatex(0.15,0.50,CBGalpha->Data());
  // TLatex *lCBGn = new TLatex(0.15,0.46,CBGn->Data());
  // TLatex *lCBGsigma = new TLatex(0.15,0.42,CBGsigma->Data());
  // TLatex *lCBGwideFactor = new TLatex(0.15,0.38,CBGwideFactor->Data());
  TLatex *lCBG = new TLatex(0.15,0.72,"CB + Gauss:");
  TLatex *lCBGalpha = new TLatex(0.15,0.68,CBGalpha->Data());
  TLatex *lCBGn = new TLatex(0.15,0.64,CBGn->Data());
  TLatex *lCBGsigma = new TLatex(0.15,0.60,CBGsigma->Data());
  TLatex *lCBGwideFactor = new TLatex(0.15,0.56,CBGwideFactor->Data());
  TLatex *lCBGcoeffGauss = new TLatex(0.15,0.52,CBGcoeffGauss->Data());
  lCBG->SetTextSize(0.035);
  lCBGalpha->SetTextSize(0.035);
  lCBGn->SetTextSize(0.035);
  lCBGsigma->SetTextSize(0.035);
  lCBGwideFactor->SetTextSize(0.035);
  lCBGcoeffGauss->SetTextSize(0.035);

  lCB->SetNDC(kTRUE);
  lCBalpha->SetNDC(kTRUE);
  lCBn->SetNDC(kTRUE);
  lCBsigma->SetNDC(kTRUE);
  lCBG->SetNDC(kTRUE);
  lCBGalpha->SetNDC(kTRUE);
  lCBGn->SetNDC(kTRUE);
  lCBGsigma->SetNDC(kTRUE);
  lCBGwideFactor->SetNDC(kTRUE);
  lCBGcoeffGauss->SetNDC(kTRUE);

  // lCB->Draw();
  // lCBalpha->Draw();
  // lCBn->Draw();
  // lCBsigma->Draw();
  lCBG->Draw();
  lCBGalpha->Draw();
  lCBGn->Draw();
  lCBGsigma->Draw();
  lCBGwideFactor->Draw();
  lCBGcoeffGauss->Draw();

  // TCanvas *c2 = new TCanvas("c2","c2");
  // c2->SetLogy();
  // RooArgList parameters(mean, sigma1, alpha, n, wideFactor, coeffGauss, a0, Nsignal, Nbackground);
  // TH1 *h1 = (TH1*) redData->createHistogram("Jpsi_Mass")->Clone("h1");
  // TF1 *f1 = (TF1*) model_cbg.asTF(RooArgList(Jpsi_Mass), parameters)->Clone("f1");
  // h1->Scale(50.0/h1->GetEntries());
  // h1->Draw();
  // //f1->SetParameter(7,f1->GetParameter(7)*1e4);
  // h1->Fit(f1,"RMEI");
  // //  f1->SetParameters(3.0968,2.9504e-02,1.7104,1.3338e+01,2.1972,1.0422e-01,0,9.9117e-03);
  // f1->Draw("same");

  TString outfname;
  if (fitChi2) {
    if (isHI)
      outfname = Form("20140326_MC_plots/Jpsi_PbPb_MCshape_Rap_%3.1f-%3.1f_Pt_%3.1f-%3.1f_M2234_chi2.pdf",ymin,ymax,ptmin,ptmax);
    else
      outfname = Form("20140326_MC_plots/Jpsi_pp_MCshape_Rap_%3.1f-%3.1f_Pt_%3.1f-%3.1f_M2234_chi2.pdf",ymin,ymax,ptmin,ptmax);
  }
  else {
    if (isHI)
      outfname = Form("20140326_MC_plots/Jpsi_PbPb_MCshape_Rap_%3.1f-%3.1f_Pt_%3.1f-%3.1f_M2234.pdf",ymin,ymax,ptmin,ptmax);
    else
      outfname = Form("20140326_MC_plots/Jpsi_pp_MCshape_Rap_%3.1f-%3.1f_Pt_%3.1f-%3.1f_M2234.pdf",ymin,ymax,ptmin,ptmax);
  }

  if (savePlot) {
    std::cout << outfname << std::endl;
    c1->SaveAs(outfname);
  }


  /*
    TCanvas *c2 = new TCanvas("c2","c2");
    c2->SetLogy();
    RooPlot* recPtframe = Jpsi_Pt.frame(Title("data"));
    if (isHI)
    tmp->plotOn(recPtframe);
    else
    redData->plotOn(recPtframe);
    recPtframe->Draw();
  
    TCanvas *c3 = new TCanvas("c3","c3");
    //  c3->SetLogy();
    RooPlot* genPtframe = Gen_Pt.frame(Title("data"));
    if (isHI)
    tmp->plotOn(genPtframe);
    else
    redData->plotOn(genPtframe);
    genPtframe->Draw();
  */
  return;
}
