#include <iostream>

//#ifndef __CINT__
#include "RooGlobalFunc.h"
//#endif

#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"
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

void check_fit_bias_sim(const int N=1, string infname="20140409_SimFits_M1850_DblMu0_AllCent/fracLogCBG_PbPbpol3_HI020_pol2_HI2040_pol2_HI40100_pppol2_rap16-24_pT3-30_centMult_Workspace.root")
{
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  TFile *inf = new TFile(infname.c_str(),"READ");
  RooWorkspace *ws = (RooWorkspace*) inf->Get("workspace");
  // ws->var("doubleRatio_HI020")->setConstant(true);
  // ws->var("doubleRatio_HI2040")->setConstant(true);
  // ws->var("doubleRatio_HI40100")->setConstant(true);

  RooRealVar *Jpsi_Mass = ws->var("Jpsi_Mass");
  RooCategory *sample = ws->cat("sample");

  // TIter types(ws->cat("sample")->typeIterator());
  // RooCatType * type;
  // types.Reset();
  RooDataSet * protoData = (RooDataSet*) ws->data("redDataSim");
  RooDataSet * protoData_pp = (RooDataSet*) ws->data("data4");
  RooDataSet * protoData_HI020 = (RooDataSet*) ws->data("data1");
  RooDataSet * protoData_HI2040 = (RooDataSet*) ws->data("data2");
  RooDataSet * protoData_HI40100 = (RooDataSet*) ws->data("data3");

  //  RooSimultaneous *model = (RooSimultaneous*) ws->pdf("sigMassPDFSim");
  RooAbsPdf *model_pp = (RooAbsPdf*) ws->pdf("sigMassPDF_pp");
  RooAbsPdf *model_HI020 = (RooAbsPdf*) ws->pdf("pdf_HI020");
  RooAbsPdf *model_HI2040 = (RooAbsPdf*) ws->pdf("pdf_HI2040");
  RooAbsPdf *model_HI40100 = (RooAbsPdf*) ws->pdf("pdf_HI40100");

  // ws->var("sigma_pol")->setVal(0.0);
  // ws->var("sigma_b")->setVal(0.0);
  // ws->var("sigma_pol")->setConstant(true);
  // ws->var("sigma_b")->setConstant(true);

  // ws->var("sigma_fit_HI020")->setVal(0.0);
  // ws->var("sigma_eff_HI020")->setVal(0.0);
  // ws->var("sigma_fit_HI020")->setConstant(true);
  // ws->var("sigma_eff_HI020")->setConstant(true);

  // ws->var("sigma_fit_HI2040")->setVal(0.0);
  // ws->var("sigma_eff_HI2040")->setVal(0.0);
  // ws->var("sigma_fit_HI2040")->setConstant(true);
  // ws->var("sigma_eff_HI2040")->setConstant(true);

  // ws->var("sigma_fit_HI40100")->setVal(0.0);
  // ws->var("sigma_eff_HI40100")->setVal(0.0);
  // ws->var("sigma_fit_HI40100")->setConstant(true);
  // ws->var("sigma_eff_HI40100")->setConstant(true);

  RooPlot *frame = Jpsi_Mass->frame(Range("M2242"));
  protoData->plotOn(frame, Cut("sample==sample::pp"));
  protoData->plotOn(frame, Cut("sample==sample::HI020"), MarkerColor(kRed));
  protoData->plotOn(frame, Cut("sample==sample::HI2040"), MarkerStyle(25), MarkerColor(kBlue));
  protoData->plotOn(frame, Cut("sample==sample::HI40100"), MarkerStyle(24), MarkerColor(kGreen+2));

  //  model->plotOn(frame, Slice(*sample,"pp"), ProjWData(*sample,*protoData));
  TCanvas *c0 = new TCanvas("c0","c0");
  c0->cd();
  frame->Draw();
  //  return;
  //  ws->pdf("sigMassPDF_M2242")->plotOn(frame,Range("M2242"));
  //  redData->plotOn(frame,MarkerStyle(24),MarkerColor(kRed),Range("M2242"));


  // while ((type=(RooCatType*)types.Next())) {
  //   protoData = 
  //     (RooDataSet*)ws->data(TString::Format("data_%s", type->GetName()));
  //   bkg = ws.pdf(TString::Format("expFunct_%s", type->GetName()));
  //   if ((*type) == "HI")
  //     N = Nhi;
  //   else
  //     N = Npp;
  //   RooDataSet * tmpData = bkg->generate(RooArgSet(*Jpsi_mass), N,
  // 					 RooFit::ProtoData(*protoData));
  //   osBkgData->append(*tmpData);
  //   delete tmpData;
  // }



  //  double Nevents =  ws->var("NJpsi_pp")->getVal()+ws->function("NPsiP_pp")->getVal()+ws->var("NBkg_pp")->getVal();

  TH1F *h0_pp = new TH1F("h0_pp","h0_pp;R_{#psi}^{pp};Events",1200,-0.3,0.3);
  TH1F *h1_pp = new TH1F("h1_pp","h1_pp;R_{#psi}^{pp};Events",1200,-0.3,0.3);
  TH2F *h2_pp = new TH2F("h2_pp","h2_pp;R_{#psi}^{pp} (M1850);R_{#psi}^{pp} (M2242)",600,-0.3,0.3,600,-0.3,0.3);
  h0_pp->Sumw2();
  h1_pp->Sumw2();
  h2_pp->Sumw2();
  h1_pp->SetMarkerColor(2);
  h1_pp->SetLineColor(2);
  h0_pp->GetXaxis()->CenterTitle(true);
  h1_pp->GetXaxis()->CenterTitle(true);
  h2_pp->GetXaxis()->CenterTitle(true);

  TH1F *h0_HI020 = new TH1F("h0_HI020","h0_HI020;#chi_{#psi}^{HI020};Events",5000,-1.0,4.0);
  TH1F *h1_HI020 = new TH1F("h1_HI020","h1_HI020;#chi_{#psi}^{HI020};Events",5000,-1.0,4.0);
  TH2F *h2_HI020 = new TH2F("h2_HI020","h2_HI020;#chi_{#psi}^{HI020} (M1850);#chi_{#psi}^{HI020} (M2242)",500,-1,4,500,-1,4);
  h0_HI020->Sumw2();
  h1_HI020->Sumw2();
  h2_HI020->Sumw2();
  h1_HI020->SetMarkerColor(2);
  h1_HI020->SetLineColor(2);
  h0_HI020->GetXaxis()->CenterTitle(true);
  h1_HI020->GetXaxis()->CenterTitle(true);
  h2_HI020->GetXaxis()->CenterTitle(true);

  TH1F *h0_HI2040 = new TH1F("h0_HI2040","h0_HI2040;#chi_{#psi}^{HI2040};Events",5000,-1.0,4.0);
  TH1F *h1_HI2040 = new TH1F("h1_HI2040","h1_HI2040;#chi_{#psi}^{HI2040};Events",5000,-1.0,4.0);
  TH2F *h2_HI2040 = new TH2F("h2_HI2040","h2_HI2040;#chi_{#psi}^{HI2040} (M1850);#chi_{#psi}^{HI2040} (M2242)",500,-1,4,500,-1,4);
  h0_HI2040->Sumw2();
  h1_HI2040->Sumw2();
  h2_HI2040->Sumw2();
  h1_HI2040->SetMarkerColor(2);
  h1_HI2040->SetLineColor(2);
  h0_HI2040->GetXaxis()->CenterTitle(true);
  h1_HI2040->GetXaxis()->CenterTitle(true);
  h2_HI2040->GetXaxis()->CenterTitle(true);

  TH1F *h0_HI40100 = new TH1F("h0_HI40100","h0_HI40100;#chi_{#psi}^{HI40100};Events",5000,-1.0,4.0);
  TH1F *h1_HI40100 = new TH1F("h1_HI40100","h1_HI40100;#chi_{#psi}^{HI40100};Events",5000,-1.0,4.0);
  TH2F *h2_HI40100 = new TH2F("h2_HI40100","h2_HI40100;#chi_{#psi}^{HI40100} (M1850);#chi_{#psi}^{HI40100} (M2242)",500,-1,4,500,-1,4);
  h0_HI40100->Sumw2();
  h1_HI40100->Sumw2();
  h2_HI40100->Sumw2();
  h1_HI40100->SetMarkerColor(2);
  h1_HI40100->SetLineColor(2);
  h0_HI40100->GetXaxis()->CenterTitle(true);
  h1_HI40100->GetXaxis()->CenterTitle(true);
  h2_HI40100->GetXaxis()->CenterTitle(true);


  Jpsi_Mass->setRange("signal",3.6,3.76);
  Jpsi_Mass->setRange("M2045",2.0,4.5);
  Jpsi_Mass->setRange("M2242",2.2,4.2);

  RooFitResult *fitMall;
  RooFitResult *fitM;

  // RooFit cannot normalize Chebychev polynomials to a subrange (no analytic integral?)
  // using normal polynomials instead
  RooRealVar a("a","a",0.0);a.setConstant(false);
  RooRealVar b("b","b",0.01);b.setConstant(false);
  RooRealVar c("c","c",-0.005);c.setConstant(false);
  // mid
  //  RooPolynomial bkg_pp("bkg_pp","bkg_pp",*Jpsi_Mass,RooArgSet(a,b,c));
  RooPolynomial bkg_pp("bkg_pp","bkg_pp",*Jpsi_Mass,RooArgSet(a));//,b,c));

  RooRealVar a_HI020("a_HI020","a_HI020",0.0);a_HI020.setConstant(false);
  RooRealVar b_HI020("b_HI020","b_HI020",0.01);b_HI020.setConstant(false);
  RooRealVar c_HI020("c_HI020","c_HI020",-0.005);c_HI020.setConstant(false);
  // mid
  //  RooPolynomial bkg_HI020("bkg_HI020","bkg_HI020",*Jpsi_Mass,RooArgSet(a_HI020));//,b_HI020,c_HI020));
  // fwd
  RooPolynomial bkg_HI020("bkg_HI020","bkg_HI020",*Jpsi_Mass,RooArgSet(a_HI020,b_HI020,c_HI020));

  RooRealVar a_HI2040("a_HI2040","a_HI2040",0.0);a_HI2040.setConstant(false);
  RooRealVar b_HI2040("b_HI2040","b_HI2040",0.01);b_HI2040.setConstant(false);
  RooRealVar c_HI2040("c_HI2040","c_HI2040",-0.005);c_HI2040.setConstant(false);
  // mid
  //  RooPolynomial bkg_HI2040("bkg_HI2040","bkg_HI2040",*Jpsi_Mass,RooArgSet(a_HI2040));//,b_HI2040,c_HI2040));
  // fwd
  RooPolynomial bkg_HI2040("bkg_HI2040","bkg_HI2040",*Jpsi_Mass,RooArgSet(a_HI2040,b_HI2040));//,c_HI2040));

  RooRealVar a_HI40100("a_HI40100","a_HI40100",0.0);a_HI40100.setConstant(false);
  RooRealVar b_HI40100("b_HI40100","b_HI40100",0.01);b_HI40100.setConstant(false);
  RooRealVar c_HI40100("c_HI40100","c_HI40100",-0.005);c_HI40100.setConstant(false);
  // mid
  // a_HI40100.setConstant(true);
  //  RooPolynomial bkg_HI40100("bkg_HI40100","bkg_HI40100",*Jpsi_Mass,RooArgSet(a_HI40100));//,b_HI40100,c_HI40100));
  // fwd
  RooPolynomial bkg_HI40100("bkg_HI40100","bkg_HI40100",*Jpsi_Mass,RooArgSet(a_HI40100));//,b_HI40100,c_HI40100));

  ws->import(bkg_pp);
  ws->import(bkg_HI020);
  ws->import(bkg_HI2040);
  ws->import(bkg_HI40100);

  ws->factory("SUM::sigMassPDF_pp_M2242(NJpsi_pp*sigCB1G2_HI,NPsiP_pp*sigCB1G2P_HI,NBkg_pp*bkg_pp)");
  ws->factory("SUM::sigMassPDF_HI020_M2242(NJpsi_HI020*sigCB1G2_HI,NPsiP_HI020*sigCB1G2P_HI,NBkg_HI020*bkg_HI020)");
  ws->factory("SUM::sigMassPDF_HI2040_M2242(NJpsi_HI2040*sigCB1G2_HI,NPsiP_HI2040*sigCB1G2P_HI,NBkg_HI2040*bkg_HI2040)");
  ws->factory("SUM::sigMassPDF_HI40100_M2242(NJpsi_HI40100*sigCB1G2_HI,NPsiP_HI40100*sigCB1G2P_HI,NBkg_HI40100*bkg_HI40100)");

  ws->factory("SIMUL::sigMassPDFSim_M2242(sample,HI020=sigMassPDF_HI020_M2242,HI2040=sigMassPDF_HI2040_M2242,HI40100=sigMassPDF_HI40100_M2242,pp=sigMassPDF_pp_M2242)");
    
  RooDataSet *data[N];
  RooDataSet *data_pp[N];
  RooDataSet *data_HI020[N];
  RooDataSet *data_HI2040[N];
  RooDataSet *data_HI40100[N];
  RooDataSet *redData;

  for (int i=0;i<N;++i) {
    cout << "Generating event " << i << "/" << N << endl;
    // cout << "Generate N = " << Nevents << " events with R_{psi} = " << ws->var("fracP_pp")->getVal() << endl;
    //    data[i] = model->generateSimGlobal(*Jpsi_Mass,ProtoData(*protoData), Verbose(1));
    data_pp[i] = model_pp->generate(*Jpsi_Mass, ProtoData(*protoData_pp),Verbose(0));
    data_HI020[i] = model_HI020->generate(*Jpsi_Mass, ProtoData(*protoData_HI020),Verbose(0));
    data_HI2040[i] = model_HI2040->generate(*Jpsi_Mass, ProtoData(*protoData_HI2040),Verbose(0));
    data_HI40100[i] = model_HI40100->generate(*Jpsi_Mass, ProtoData(*protoData_HI40100),Verbose(0));

    data[i] = new RooDataSet("data","data",RooArgSet(*Jpsi_Mass),Index(*sample),Import("HI020",*data_HI020[i]),Import("HI2040",*data_HI2040[i]),Import("HI40100",*data_HI40100[i]),Import("pp",*data_pp[i]));
  }
    
  RooPlot *frame2 = Jpsi_Mass->frame();
  data_pp[0]->plotOn(frame2);
  data[N-1]->plotOn(frame2, Cut("sample==sample::pp"),MarkerColor(kRed), MarkerStyle(24));
  data[N-1]->plotOn(frame2, Cut("sample==sample::HI020"), MarkerColor(kRed));
  data[N-1]->plotOn(frame2, Cut("sample==sample::HI2040"), MarkerStyle(25), MarkerColor(kBlue));
  data[N-1]->plotOn(frame2, Cut("sample==sample::HI40100"), MarkerStyle(24), MarkerColor(kGreen+2));
  // data[N-1]->plotOn(frame2, Cut("sample!=sample::pp"), MarkerColor(kRed));
  //  model->plotOn(frame, Slice(*sample,"pp"), ProjWData(*sample,*protoData));
  TCanvas *c0a = new TCanvas("c0a","c0a");
  c0a->cd();
  frame2->Draw();
  //  return;
  // cout << data->sumEntries() << endl;

  for (int i=0;i<N;++i) {
    cout << "Fitting event " << i << "/" << N << endl;
    fitMall = ws->pdf("sigMassPDFSim")->fitTo(*data[i],Extended(1),Hesse(1),Save(1),NumCPU(8),PrintEvalErrors(-1),Verbose(0),PrintLevel(-1));
    if (fitMall->statusCodeHistory(fitMall->numStatusHistory()-1) != 0) {i--; continue;}
    //    fitMall->Print("v");

    //    cout << "Fitted R_{psi} = " << ws->var("fracP_pp")->getVal() << " +/- " << ws->var("fracP_pp")->getPropagatedError(*fitMall) << endl;
    double R1850_pp = ws->var("fracP_pp")->getVal();
    double R1850_HI020 = ws->var("doubleRatio_HI020")->getVal();
    double R1850_HI2040 = ws->var("doubleRatio_HI2040")->getVal();
    double R1850_HI40100 = ws->var("doubleRatio_HI40100")->getVal();
    h0_pp->Fill(R1850_pp);
    h0_HI020->Fill(R1850_HI020);
    h0_HI2040->Fill(R1850_HI2040);
    h0_HI40100->Fill(R1850_HI40100);


    redData = (RooDataSet*)data[i]->reduce("Jpsi_Mass>2.2&&Jpsi_Mass<4.2");
    //  cout << redData->sumEntries() << endl;

    fitM = ws->pdf("sigMassPDFSim_M2242")->fitTo(*redData,Extended(1),Hesse(1),Save(1),NumCPU(8),PrintEvalErrors(-1),Verbose(0),PrintLevel(-1),Range("M2242"));
    if (fitM->statusCodeHistory(fitM->numStatusHistory()-1) != 0) {i--; continue;}

    //    fitM->Print("v");
    // cout << "Fit over M2242: R_{psi} = " << ws->var("fracP_pp")->getVal() << " +/- " << ws->var("fracP_pp")->getPropagatedError(*fitM) << endl;

    // RooPlot *frame = Jpsi_Mass->frame(Range("M2242"));
    // redData->plotOn(frame);
    // ws->pdf("sigMassPDF_M2242")->plotOn(frame,Range("M2242"));
    //  redData->plotOn(frame,MarkerStyle(24),MarkerColor(kRed),Range("M2242"));

    double R2242_pp = ws->var("fracP_pp")->getVal();
    double R2242_HI020 = ws->var("doubleRatio_HI020")->getVal();
    double R2242_HI2040 = ws->var("doubleRatio_HI2040")->getVal();
    double R2242_HI40100 = ws->var("doubleRatio_HI40100")->getVal();
    h1_pp->Fill(R2242_pp);
    h1_HI020->Fill(R2242_HI020);
    h1_HI2040->Fill(R2242_HI2040);
    h1_HI40100->Fill(R2242_HI40100);

    h2_pp->Fill(R1850_pp,R2242_pp);
    h2_HI020->Fill(R1850_HI020,R2242_HI020);
    h2_HI2040->Fill(R1850_HI2040,R2242_HI2040);
    h2_HI40100->Fill(R1850_HI40100,R2242_HI40100);
  }
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->Divide(2,2);
  c1->cd(1);
  h0_pp->Draw();
  h1_pp->Draw("same");
  c1->cd(2);
  h0_HI020->Draw();
  h1_HI020->Draw("same");
  c1->cd(3);
  h0_HI2040->Draw();
  h1_HI2040->Draw("same");
  c1->cd(4);
  h0_HI40100->Draw();
  h1_HI40100->Draw("same");
  //  frame->Draw();
  cout << "pp" << endl;
  cout << h0_pp->GetMean() << "\t" << h0_pp->GetRMS() << endl;
  cout << h1_pp->GetMean() << "\t" << h1_pp->GetRMS() << endl;
  cout << 1-(h0_pp->GetMean()/h1_pp->GetMean()) << endl;

  cout << "HI020" << endl;
  cout << h0_HI020->GetMean() << "\t" << h0_HI020->GetRMS() << endl;
  cout << h1_HI020->GetMean() << "\t" << h1_HI020->GetRMS() << endl;
  cout << 1-(h0_HI020->GetMean()/h1_HI020->GetMean()) << endl;

  cout << "HI2040" << endl;
  cout << h0_HI2040->GetMean() << "\t" << h0_HI2040->GetRMS() << endl;
  cout << h1_HI2040->GetMean() << "\t" << h1_HI2040->GetRMS() << endl;
  cout << 1-(h0_HI2040->GetMean()/h1_HI2040->GetMean()) << endl;

  cout << "HI40100" << endl;
  cout << h0_HI40100->GetMean() << "\t" << h0_HI40100->GetRMS() << endl;
  cout << h1_HI40100->GetMean() << "\t" << h1_HI40100->GetRMS() << endl;
  cout << 1-(h0_HI40100->GetMean()/h1_HI40100->GetMean()) << endl;

  c1->SaveAs(Form("toy_fits_dblRatio_fwd_M1850_N%i.pdf",N));

  TCanvas *c2 = new TCanvas("c2","c2");
  c2->Divide(2,2);
  c2->cd(1);
  h2_pp->Draw("colz");
  c2->cd(2);
  h2_HI020->Draw("colz");
  c2->cd(3);
  h2_HI2040->Draw("colz");
  c2->cd(4);
  h2_HI40100->Draw("colz");

  TF1 *f3 = new TF1("f3","x",-0.5,3.0);
  f3->SetLineWidth(1);
  f3->Draw("same");

  TFile *outf = new TFile(Form("toy_fits_dblRatio_fwd_M1850_N%i.root",N),"RECREATE");
  h0_pp->Write();
  h1_pp->Write();
  h2_pp->Write();
  h0_HI020->Write();
  h1_HI020->Write();
  h2_HI020->Write();
  h0_HI2040->Write();
  h1_HI2040->Write();
  h2_HI2040->Write();
  h0_HI40100->Write();
  h1_HI40100->Write();
  h2_HI40100->Write();
  outf->Close();


  return;
}
