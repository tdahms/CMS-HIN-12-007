#include <iostream>

//#ifndef __CINT__
#include "RooGlobalFunc.h"
//#endif

#include "RooWorkspace.h"
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
#include "RooFitResult.h"

#include "RooGenericPdf.h"
#include "RooMinuit.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"

using namespace RooFit;

void check_fit_bias_bkgmodel(const int N=1, string infname="20140514_SimFits_M2242_DblMu0_AllCent_WithSyst_final_forPreliminary/fracLogCBG_PbPbexpPol3_HI020_expPol2_HI2040_expPol1_HI40100_ppexpPol3_rap16-24_pT3-30_centMult_Workspace.root")
{
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  TFile *inf = new TFile(infname.c_str(),"READ");
  RooWorkspace *ws = (RooWorkspace*) inf->Get("workspace");
  ws->var("fracP_pp")->setConstant(true);

  RooRealVar *Jpsi_Mass = ws->var("Jpsi_Mass");
  RooAbsPdf *model = ws->pdf("sigMassPDF_HI020");

  double Nevents =  ws->var("NJpsi_HI020")->getVal()+ws->function("NPsiP_HI020")->getVal()+ws->var("NBkg_HI020")->getVal();

  TH1F *h0 = new TH1F("h0","h0",300,0.0,0.3);
  TH1F *h1 = new TH1F("h1","h1",300,0.0,0.3);
  TH2F *h2 = new TH2F("h2","h2;R_{#psi} (exp);R_{#psi} (lin)",300,0.0,0.3,300,0.0,0.3);
  h0->Sumw2();
  h1->Sumw2();
  h2->Sumw2();
  h1->SetMarkerColor(2);
  h1->SetLineColor(2);

  Jpsi_Mass->setRange("signal",3.6,3.76);
  //  Jpsi_Mass->setRange("M2045",2.0,4.5);
  //  Jpsi_Mass->setRange("M2242",2.2,4.2);

  RooFitResult *fitMall;
  RooFitResult *fitM;

  // RooFit cannot normalize Chebychev polynomials to a subrange (no analytic integral?)
  // using normal polynomials instead
  RooRealVar a("a","a",0.5);a.setConstant(false);
  RooRealVar b("b","b",0.01);b.setConstant(false);
  RooRealVar c("c","c",-0.005);c.setConstant(false);
  RooPolynomial bkg("bkg","bkg",*Jpsi_Mass,RooArgSet(a,b,c));
  ws->import(bkg);
  ws->factory("SUM::sigMassPDF_lin(NJpsi_HI020*sigCB1G2_HI,NPsiP_HI020*sigCB1G2P_HI,NBkg_HI020*bkg)");

  RooDataSet *data[N];
  //  RooDataSet *redData;

  for (int i=0;i<N;++i) {
    cout << "Generating event " << i << "/" << N << endl;
    // cout << "Generate N = " << Nevents << " events with R_{psi} = " << ws->function("fracP_HI020")->getVal() << endl;
    data[i] = model->generate(*Jpsi_Mass,Nevents);
  }
    // cout << data->sumEntries() << endl;

  for (int i=0;i<N;++i) {
    cout << "Fitting event " << i << "/" << N << endl;
    fitMall = ws->pdf("sigMassPDF_HI020")->fitTo(*data[i],Extended(1),Hesse(1),Save(1),NumCPU(8),PrintEvalErrors(-1),Verbose(0),PrintLevel(-1));
    if (fitMall->statusCodeHistory(fitMall->numStatusHistory()-1) != 0) {i--; continue;}

    //  cout << "Fitted R_{psi} = " << ws->function("fracP_HI020")->getVal() << " +/- " << ws->function("fracP_HI020")->getPropagatedError(*fitMall) << endl;
    double Rexp = ws->function("fracP_HI020")->getVal();
    h0->Fill(Rexp);


    //    redData = (RooDataSet*)data[i]->reduce("Jpsi_Mass>2.2&&Jpsi_Mass<4.2");
    //  cout << redData->sumEntries() << endl;

    fitM = ws->pdf("sigMassPDF_lin")->fitTo(*data[i],Extended(1),Hesse(1),Save(1),NumCPU(8),PrintEvalErrors(-1),Verbose(0),PrintLevel(-1));
    if (fitM->statusCodeHistory(fitM->numStatusHistory()-1) != 0) {i--; continue;}

    double Rlin = ws->function("fracP_HI020")->getVal();
    h1->Fill(Rlin);
    h2->Fill(Rexp,Rlin);
  }
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->cd();
  h0->Draw();
  h1->Draw("same");
  //  frame->Draw();
  cout << h0->GetMean() << "\t" << h0->GetRMS() << endl;
  cout << h1->GetMean() << "\t" << h1->GetRMS() << endl;

  c1->SaveAs(Form("toy_fits_bkgModel_N%i.pdf",N));

  TCanvas *c2 = new TCanvas("c2","c2");
  c2->cd();
  h2->Draw("colz");

  TF1 *f3 = new TF1("f3","x",0,0.3);
  f3->SetLineWidth(1);
  f3->Draw("same");

  TFile *outf = new TFile(Form("toy_fits_bkgModel_new_N%i.root",N),"RECREATE");
  h0->Write();
  h1->Write();
  h2->Write();
  outf->Close();


  return;
}
