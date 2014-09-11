#include <iostream>

//#ifndef __CINT__
#include "RooGlobalFunc.h"
//#endif

#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooExtendPdf.h"
#include "RooConstVar.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooPolynomial.h"
#include "RooChebychev.h"
#include "RooFormulaVar.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHist.h"
#include "RooFitResult.h"
#include "RooBinning.h"

#include "RooGenericPdf.h"
#include "RooChi2Var.h"
#include "RooMinuit.h"

#include "RooMsgService.h"

#include "TCut.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLatex.h"
#include "TFile.h"
#include "TAxis.h"
#include "TH1.h"
#include "TF1.h"

using namespace RooFit;

void fitSideBands(bool isHI=false, double ptmin=0.0, double ptmax=30.0, double ymin=0.0, double ymax=2.4, bool absRapidity=true, bool savePlot=false, bool fitAll=false, unsigned int n=2, double N=1.5)
{
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  TCut ptCut = Form("Jpsi_Pt>%3.1f&&Jpsi_Pt<%3.1f",ptmin,ptmax);
  TCut rapCut;
  if (absRapidity)
    rapCut = Form("abs(Jpsi_Y)>%3.1f&&abs(Jpsi_Y)<%3.1f",ymin,ymax);
  else
    rapCut = Form("Jpsi_Y>%3.1f&&Jpsi_Y<%3.1f",ymin,ymax);

  double ctmax = 1;
  if (isHI) {
    if (ptmin==6.5 && ptmax==30.0 && ymin==0.0 && ymax==2.4)
      ctmax = 0.04;
    else if (ptmin==6.5 && ptmax==30.0 && ymin==0.0 && ymax==1.6)
      ctmax = 0.04;
    else if (ptmin==3.0 && ptmax==30.0 && ymin==1.6 && ymax==2.4)
      ctmax = 0.08;
    else if (ptmin==3.0 && ptmax==6.5 && ymin==1.6 && ymax==2.4)
      ctmax = 0.09;
    else if (ptmin==6.5 && ptmax==30.0 && ymin==1.6 && ymax==2.4)
      ctmax = 0.06;
  }
  else {
    if (ptmin==6.5 && ptmax==30.0 && ymin==0.0 && ymax==2.4)
      ctmax = 0.04;
    else if (ptmin==6.5 && ptmax==30.0 && ymin==0.0 && ymax==1.6)
      ctmax = 0.04;
    else if (ptmin==3.0 && ptmax==30.0 && ymin==1.6 && ymax==2.4)
      ctmax = 0.09;
    else if (ptmin==3.0 && ptmax==6.5 && ymin==1.6 && ymax==2.4)
      ctmax = 0.11;
    else if (ptmin==6.5 && ptmax==30.0 && ymin==1.6 && ymax==2.4)
      ctmax = 0.06;
  }

  TCut ctauCut = Form("Jpsi_Ct<%3.3f",ctmax);

  RooRealVar Jpsi_Mass("Jpsi_Mass","m_{#mu^{+}#mu^{-}}",2.2,3.4,"GeV/c^{2}");
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
    rbm.addUniform(5,3.600,3.720); // 24 MeV
    rbm.addUniform(8,3.720,4.2); // 60 MeV
  }

  Jpsi_Mass.setBinning(rbm);
  // Jpsi_Pt.setBins(60);
  // Jpsi_Y.setBins(1);
  // Jpsi_Ct.setBins(1);
  // Jpsi_CtErr.setBins(1);

  RooRealVar A0("A0","A0",0.0);A0.setConstant(true);
  RooRealVar a0("a0","a0",0.0);a0.setConstant(false);

  RooRealVar b0("b0","b0",0.0);b0.setConstant(false);
  RooRealVar b1("b1","b1",0.0);b1.setConstant(false);

  RooRealVar c0("c0","c0",0.0);c0.setConstant(false);
  RooRealVar c1("c1","c1",0.0);c1.setConstant(false);
  RooRealVar c2("c2","c2",0.0);c2.setConstant(false);

  RooRealVar d0("d0","d0",0.0);d0.setConstant(false);
  RooRealVar d1("d1","d1",0.0);d1.setConstant(false);
  RooRealVar d2("d2","d2",0.0);d2.setConstant(false);
  RooRealVar d3("d3","d3",0.0);d3.setConstant(false);

  RooRealVar e0("e0","e0",0.0);e0.setConstant(false);
  RooRealVar e1("e1","e1",0.0);e1.setConstant(false);
  RooRealVar e2("e2","e2",0.0);e2.setConstant(false);
  RooRealVar e3("e3","e3",0.0);e3.setConstant(false);
  RooRealVar e4("e4","e4",0.0);e4.setConstant(false);

  RooRealVar f0("f0","f0",0.0);f0.setConstant(false);
  RooRealVar f1("f1","f1",0.0);f1.setConstant(false);
  RooRealVar f2("f2","f2",0.0);f2.setConstant(false);
  RooRealVar f3("f3","f3",0.0);f3.setConstant(false);
  RooRealVar f4("f4","f4",0.0);f4.setConstant(false);
  RooRealVar f5("f5","f5",0.0);f5.setConstant(false);

  RooPolynomial pol0("pol0","Background (pol0)",Jpsi_Mass,RooArgSet(A0));

  RooPolynomial pol1("pol1","Background (pol1)",Jpsi_Mass,RooArgSet(a0));
  RooPolynomial pol2("pol2","Background (pol2)",Jpsi_Mass,RooArgSet(b0,b1));

  RooPolynomial pol3("pol3","Background (pol3)",Jpsi_Mass,RooArgSet(c0,c1,c2));
  //  RooPolynomial pol3a("pol3a","Background (pol3a)",Jpsi_Mass,RooArgSet(d0,d1,d2));

  RooPolynomial pol4("pol4","Background (pol4)",Jpsi_Mass,RooArgSet(d0,d1,d2,d3));
  RooPolynomial pol5("pol5","Background (pol5)",Jpsi_Mass,RooArgSet(e0,e1,e2,e3,e4));
  RooPolynomial pol6("pol6","Background (pol6)",Jpsi_Mass,RooArgSet(f0,f1,f2,f3,f4,f5));
  RooRealVar nbkg("nbkg","nbkg",1000,0,100000);
  RooRealVar nbkga("nbkga","nbkga",1000,0,100000);
  RooRealVar nzero("nzero","nzero",0);
  RooRealVar F0("F0","F0",0.0);F0.setConstant(true);
  RooPolynomial fake("fake","fake background",Jpsi_Mass,RooArgSet(F0));
  RooAddPdf model0("model0","model0",RooArgList(pol0,fake),RooArgList(nbkg,nzero));
  RooAddPdf model1("model1","model1",RooArgList(pol1,fake),RooArgList(nbkg,nzero));
  RooAddPdf model2("model2","model2",RooArgList(pol2,fake),RooArgList(nbkg,nzero));
  RooAddPdf model3("model3","model3",RooArgList(pol3,fake),RooArgList(nbkg,nzero));
  RooAddPdf model4("model4","model4",RooArgList(pol4,fake),RooArgList(nbkg,nzero));
  RooAddPdf model5("model5","model5",RooArgList(pol5,fake),RooArgList(nbkg,nzero));
  RooAddPdf model6("model6","model6",RooArgList(pol6,fake),RooArgList(nbkg,nzero));

  vector<RooAddPdf> models;
  models.push_back(model0);
  models.push_back(model1);
  models.push_back(model2);
  models.push_back(model3);
  models.push_back(model4);
  models.push_back(model5);
  models.push_back(model6);
  vector<double> nll;

  //  RooExtendPdf model2("model2","model2",pol3a,nbkga);

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
  redData = (RooDataSet*)data->reduce(list, (ptCut+rapCut+ctauCut).GetTitle());
  redData->Print("v");

  cout << (ptCut+rapCut+ctauCut).GetTitle() << endl;

  // RooFitResult *fitM1;
  // RooFitResult *fitM2;
  // RooFitResult *fitM3;
  // RooFitResult *fitM4;
  // RooFitResult *fitM5;
  // RooFitResult *fitM6;
  vector<RooFitResult *>fitM;
  //  RooFitResult *fitMa;
  // fitM1 = pol1.fitTo(*redData,Extended(0),Minos(0),Save(1),SumW2Error(kFALSE),NumCPU(4),Range("SBLow,SBMid"));
  // fitM2 = pol2.fitTo(*redData,Extended(0),Minos(0),Save(1),SumW2Error(kFALSE),NumCPU(4),Range("SBLow,SBMid"));
  // fitM3 = pol3.fitTo(*redData,Extended(0),Minos(0),Save(1),SumW2Error(kFALSE),NumCPU(4),Range("SBLow,SBMid"));
  // fitM4 = pol4.fitTo(*redData,Extended(0),Minos(0),Save(1),SumW2Error(kFALSE),NumCPU(4),Range("SBLow,SBMid"));
  // fitM5 = pol5.fitTo(*redData,Extended(0),Minos(0),Save(1),SumW2Error(kFALSE),NumCPU(4),Range("SBLow,SBMid"));
  // fitM6 = pol6.fitTo(*redData,Extended(0),Minos(0),Save(1),SumW2Error(kFALSE),NumCPU(4),Range("SBLow,SBMid"));
  for (unsigned int i=0; i<models.size(); ++i) {
    double theNLL=0;
    if (fitAll || i==n) {
      fitM.push_back( models.at(i).fitTo(*redData,Extended(1),Minos(0),Save(1),SumW2Error(kFALSE),NumCPU(4),Range("SBLow,SBMid,SBHigh"),PrintEvalErrors(0),Verbose(0)));
    //  fitMa = model2.fitTo(*redData,Extended(1),Minos(0),Save(1),SumW2Error(kFALSE),NumCPU(4),Range("SBLow,SBHigh"));
      theNLL = fitM.at(i)->minNll();
    }
    else {
      theNLL = 0;
      fitM.push_back(NULL);
    }
    nll.push_back(theNLL);

    //      fitM.at(i)->printMultiline(std::cout,1,1,"");
    //  fitMa->printMultiline(std::cout,1,1,"");

  }
  cout <<  models.at(n).GetName() << endl;
  fitM.at(n)->printMultiline(std::cout,1,1,"");

  for (unsigned int i=0;i<nll.size();++i) {
    cout << "pol" << i << ": " << nll.at(i) << endl;
  }

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

  models.at(n).plotOn(xframe, LineColor(kRed));
  models.at(n).plotOn(xframe, LineColor(kBlue),LineStyle(kDashed),Range("Full"));


  // model2.plotOn(xframe, LineColor(kRed));
  // model2.plotOn(xframe, LineColor(kRed),LineStyle(kDashed),Range("Full"));
  // pol1.plotOn(xframe, LineColor(kBlue));
  // pol2.plotOn(xframe, LineColor(kRed));
  // //  pol3.plotOn(xframe, LineColor(kRed),Range("SBLow,SBMid,SBHigh"));
  // //  pol3.plotOn(xframe, LineColor(kBlue),LineStyle(kDashed),Range("Full"));
  // pol3.plotOn(xframe, LineColor(kGreen+2), LineWidth(2));
  // pol4.plotOn(xframe, LineColor(kOrange+2), LineWidth(2));
  // pol5.plotOn(xframe, LineColor(kMagenta-1), LineStyle(kDashed));
  // pol6.plotOn(xframe, LineColor(kCyan+2), LineStyle(kDotted));

  //  redData->plotOn(xframe,RooFit::Binning(rbm));

  // xframe->SetMinimum(90.0);
  // xframe->SetMaximum(1e4);

  RooHist *h =  xframe->residHist();
  h->GetXaxis()->SetRangeUser(2.200001,4.199999);
  h->GetXaxis()->CenterTitle(true);
  h->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  h->GetYaxis()->SetTitle("Events / (0.02 GeV/c^{2})");
  h->GetYaxis()->SetTitleOffset(1.4);

  TCanvas *cc1 = new TCanvas("cc1","cc1",1200,600);
  cc1->Divide(2,1);
  cc1->GetPad(1)->SetLeftMargin(0.14);
  cc1->GetPad(1)->SetRightMargin(0.045);
  cc1->GetPad(2)->SetLeftMargin(0.14);
  cc1->GetPad(2)->SetRightMargin(0.045);

  cc1->cd(1);
  //  cc1->SetLogy();
  xframe->GetXaxis()->CenterTitle(true);
  xframe->GetYaxis()->SetTitleOffset(1.4);
  double max = xframe->GetMaximum();
  xframe->SetMaximum(1.2*max);
  xframe->Draw();

  //  TCanvas *cc2 = new TCanvas("cc2","cc2");
  cc1->cd(2);
  h->Draw("AP");

  cc1->cd();

  double sigma = 0.030;
  if (ymin==1.6 && ymax==2.4)
    sigma = 0.051;

  double meanJpsi = 3.097;
  double meanPsiP = 3.686;

  double NJpsi = h->getFitRangeNEvt(meanJpsi-N*sigma,meanJpsi+N*sigma);
  double NpsiP = h->getFitRangeNEvt(meanPsiP-N*1.19025*sigma,meanPsiP+N*1.19025*sigma);

  double errJpsi = 0.0;
  double errPsiP = 0.0;

  for (int i=0; i<h->GetN(); ++i) {
    double x,y;
    h->GetPoint(i,x,y);
    if (x>=meanJpsi-1.5*sigma && x<=meanJpsi+1.5*sigma)
      errJpsi+=pow(h->GetErrorY(i),2);
    else if (x>=meanPsiP-1.5*1.19025*sigma && x<=meanPsiP+1.5*1.19025*sigma)
      errPsiP+=pow(h->GetErrorY(i),2);
  }
  errJpsi = sqrt(errJpsi);
  errPsiP = sqrt(errPsiP);

  cout << "Integrated J/psi yield in " << meanJpsi-1.5*sigma << "--" << meanJpsi+1.5*sigma << endl;
  cout << "Integrated psi(2S) yield in " << meanPsiP-1.5*1.19025*sigma << "--" << meanPsiP+1.5*1.19025*sigma << endl;

  cout << "NJpsi: " << NJpsi << " +/- " << errJpsi << "\tNpsi(2S): " << NpsiP << " +/- " << errPsiP  << endl;
  cout << "Rpsi: " << NpsiP/NJpsi << " +/- " << sqrt(pow(errJpsi/NJpsi,2)+pow(errPsiP/NpsiP,2))*NpsiP/NJpsi  << endl;

  TLatex *lcoll;
  if (isHI)
    lcoll = new TLatex(0.1,0.9,"CMS PbPb #sqrt{s_{NN}} = 2.76 TeV");
  else
    lcoll = new TLatex(0.1,0.9,"CMS pp #sqrt{s} = 2.76 TeV");

  lcoll->SetNDC(kTRUE);
  lcoll->Draw();

  TLatex *lpt;
  if (ptmin==0.0)
    lpt = new TLatex(0.1,0.84,Form("p_{T} < %3.1f GeV/c",ptmax));
  else
    lpt = new TLatex(0.1,0.84,Form("%3.1f < p_{T} < %3.1f GeV/c",ptmin,ptmax));
  TLatex *lrap;
  if (absRapidity){
    if (ymin==0.0)
      lrap = new TLatex(0.1,0.78,Form("|y| < %3.1f",ymax));
    else
      lrap = new TLatex(0.1,0.78,Form("%3.1f < |y| < %3.1f",ymin,ymax));
  }
  else {
    if (ymin==0.0)
      lrap = new TLatex(0.1,0.78,Form("y < %3.1f",ymax));
    else
      lrap = new TLatex(0.1,0.78,Form("%3.1f < y < %3.1f",ymin,ymax));
  }

  lpt->SetNDC(kTRUE);
  lrap->SetNDC(kTRUE);
  lpt->Draw();
  lrap->Draw();


  TString outfname;
  if (isHI)
    outfname = Form("PbPb_Mass_Sidebands_pol%d_cent0-100_Rap_%3.1f-%3.1f_Pt_%3.1f-%3.1f.pdf",n,ymin,ymax,ptmin,ptmax);
  else
    outfname = Form("pp_Mass_SideBands_pol%d_Rap_%3.1f-%3.1f_Pt_%3.1f-%3.1f.pdf",n,ymin,ymax,ptmin,ptmax);

  if (savePlot) {
    std::cout << outfname << std::endl;
    cc1->SaveAs(outfname);
  }

  return;
}
