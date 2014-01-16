
using namespace RooFit;

void fitMC_CBG(bool isHI=false, float ptmin=0.0, float ptmax=30.0, float ymin=0.0, float ymax=2.4, bool absRapidity=true, bool savePlot=false)
{
  TCut ptCut = Form("Jpsi_Pt>%3.1f&&Jpsi_Pt<%3.1f",ptmin,ptmax);
  TCut rapCut;
  if (absRapidity)
    rapCut = Form("abs(Jpsi_Y)>%3.1f&&abs(Jpsi_Y)<%3.1f",ymin,ymax);
  else
    rapCut = Form("Jpsi_Y>%3.1f&&Jpsi_Y<%3.1f",ymin,ymax);

  RooRealVar Jpsi_Mass("Jpsi_Mass","J/#psi mass",2.5,4.2,"GeV/c^{2}");
  RooRealVar Jpsi_Pt("Jpsi_Pt","J/#psi pt",0,30,"GeV/c");
  RooRealVar Jpsi_Y("Jpsi_Y","J/#psi y",-2.4,2.4);
  RooRealVar Jpsi_Ct("Jpsi_Ct","J/#psi c#tau",-3.0,3.5,"mm");
  RooRealVar Jpsi_CtErr("Jpsi_CtErr","J/#psi c#tau error",-1.,1.,"mm");
  RooRealVar Jpsi_CtTrue("Jpsi_CtTrue","J/#psi c#tau true",-3.0,3.5,"mm");
  RooRealVar Gen_Pt("Gen_Pt","Generated J/#psi pt",0,30,"GeV/c");

  Jpsi_Mass.setBins(85);
  Jpsi_Pt.setBins(60);
  Jpsi_Y.setBins(1);
  Jpsi_Ct.setBins(1);
  Jpsi_CtErr.setBins(1);
  Jpsi_CtTrue.setBins(1);
  Gen_Pt.setBins(60);

  RooRealVar mean("mean","mean", 3.0969,3.05,3.15);
  RooRealVar sigma1("sigma1","sigma1", 0.03,0.01,0.060);
  RooGaussian gaussS("gaussS","gaussS",Jpsi_Mass,mean,sigma1);

  RooRealVar alpha("alpha","alpha", 1.0,0.0,5.0);
  RooRealVar n("n","n", 2,1.0,100.0);
  RooCBShape cballS("cballS","cballS",Jpsi_Mass,mean,sigma1,alpha,n);
  RooCBShape cball("cball","cball",Jpsi_Mass,mean,sigma1,alpha,n);

  RooRealVar wideFactor("wideFactor","wideFactor",1.2,1.0,5.0);
  RooFormulaVar sigma2("sigma2","@0*@1",RooArgList(sigma1,wideFactor));
  RooGaussian gauss2("gauss","gauss",Jpsi_Mass,mean,sigma2);

  RooRealVar coeffGauss("coeffGauss","coeffGauss",0.1,0.0,1.0);
  RooAddPdf cbg("cbg","cbg",RooArgList(gauss2,cball),coeffGauss);

  RooRealVar a0("a0","a0",0.0,-1.0,1.0);
  RooRealVar a1("a1","a1",0.0,-1.0,1.0);
  RooChebychev bkg("bkg","Background",Jpsi_Mass,RooArgSet(a0));

  RooRealVar bkgfrac("bkgfrac","fraction of background",0.001,0.0,1.0);
  RooAddPdf model_g("model_g","cb+bkg",RooArgList(bkg,gaussS),bkgfrac);
  RooAddPdf model_cb("model_cb","cb+bkg",RooArgList(bkg,cballS),bkgfrac);
  RooAddPdf model_cbg("model_cbg","cb+bkg",RooArgList(bkg,cbg),bkgfrac);
  a0.setVal(0);
  a0.setConstant(1);
  bkgfrac.setVal(0);
  bkgfrac.setConstant(1);

  string fname;
  if (isHI)
    fname = "/home/tdahms/CMS/HIN-12-007/analysis/20131213_DblMu0_MC/PbPbPromptJpsiMC_DblMu0_cent0-100.root";
  else
    fname = "/home/tdahms/CMS/HIN-12-007/analysis/20131213_DblMu0_MC/ppPromptJpsiMC_DblMu0_cent0-100.root";

  TFile *inf = new TFile(fname.c_str());

  RooDataSet *data = (RooDataSet*)inf->Get("dataJpsi");
  data->SetName("data");
  RooDataSet *redData;
  RooDataSet *tmp;
  if (isHI) {
    RooFormulaVar wFunc("w","event weight","Gen_Pt>30?0:(Gen_Pt>15.0?(0.00579/166960.0):(Gen_Pt>12?(0.00979/143944.0):(Gen_Pt>9?(0.0379/167336.0):(Gen_Pt>6?(0.178/172176.0):(Gen_Pt>3?(0.915/124802.0):(1.0/117512.0))))))",Gen_Pt);
    RooRealVar *w = (RooRealVar*) data->addColumn(wFunc);
    data->Print(); 
    //    tmp = (RooDataSet*)data->reduce((ptCut+rapCut).GetTitle());
    tmp = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),(ptCut+rapCut).GetTitle(),w->GetName());
    RooArgList list(Jpsi_Mass,Gen_Pt,*w);
    RooDataSet wdata(tmp->GetName(),tmp->GetTitle(),tmp,list,0,w->GetName());
    redData = (RooDataSet*)wdata;
  }
  else {
    data->Print(); 
    redData = (RooDataSet*)data->reduce((ptCut+rapCut).GetTitle());
  }
  redData->Print(); 
  //  redData->binnedClone()->Print("v");
  if (isHI)
    RooDataHist* binnedData = redData->binnedClone();
  else
    RooDataSet *binnedData = (RooDataSet*)redData;
  binnedData->Print("v");

  model_g.fitTo(*binnedData,Extended(0),Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(4));
  model_cb.fitTo(*binnedData,Extended(0),Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(4));

  TString *CBalpha = new TString();
  CBalpha = alpha.format(2,"NEA");
  TString *CBn = new TString();
  CBn = n.format(2,"NEA");
  model_cbg.fitTo(*binnedData,Extended(0),Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(4));
  TString *CBGalpha = new TString();
  CBGalpha = alpha.format(2,"NEA");
  TString *CBGn = new TString();
  CBGn = n.format(2,"NEA");

  RooPlot* xframe = Jpsi_Mass.frame(Title("data"));
  binnedData->plotOn(xframe);
  model_g.plotOn(xframe,LineWidth(1),LineColor(kBlack));
  model_cb.plotOn(xframe,LineWidth(1), LineColor(kBlue));
  model_cbg.plotOn(xframe, LineColor(kRed));
  model_cbg.plotOn(xframe, Components(RooArgSet(gauss2)), LineColor(kBlue), LineStyle(kDashed));
  model_cbg.plotOn(xframe, Components(RooArgSet(cball)), LineColor(kGreen+2), LineStyle(kDashed));
  binnedData->plotOn(xframe);

  //  cbg.Print("t");
  cout << "CB only:\t" << *CBalpha << "\t" << *CBn << endl;
  cout << "CB + Gauss:\t" << *CBGalpha << "\t" << *CBGn << endl;

  TCanvas *c1 = new TCanvas("c1","c1");
  //  c1->SetLogy();
  xframe->Draw();

  TLatex *lcoll = new TLatex(0.5,0.9,"PYTHIA: pp #sqrt{s} = 2.76 TeV");
  lcoll->SetNDC(kTRUE);
  lcoll->Draw();

  TLatex *lpt;
  if (ptmin==0.0)
    lpt = new TLatex(0.5,0.84,Form("p_{T} < %3.1f GeV/c",ptmax));
  else
    lpt = new TLatex(0.5,0.84,Form("%3.1f < p_{T} < %3.1f GeV/c",ptmin,ptmax));
  TLatex *lrap;
  if (absRapidity){
    if (ymin==0.0)
      lrap = new TLatex(0.5,0.79,Form("|y| < %3.1f",ymax));
    else
      lrap = new TLatex(0.5,0.79,Form("%3.1f < |y| < %3.1f",ymin,ymax));
  }
  else {
    if (ymin==0.0)
      lrap = new TLatex(0.5,0.79,Form("y < %3.1f",ymax));
    else
      lrap = new TLatex(0.5,0.79,Form("%3.1f < y < %3.1f",ymin,ymax));
  }

  lpt->SetNDC(kTRUE);
  lrap->SetNDC(kTRUE);
  lpt->Draw();
  lrap->Draw();

  TLatex *lCB = new TLatex(0.5,0.72,"CB:");
  TLatex *lCBalpha = new TLatex(0.5,0.67,CBalpha->Data());
  TLatex *lCBn = new TLatex(0.5,0.62,CBn->Data());
  TLatex *lCBG = new TLatex(0.5,0.55,"CB + Gauss:");
  TLatex *lCBGalpha = new TLatex(0.5,0.50,CBGalpha->Data());
  TLatex *lCBGn = new TLatex(0.5,0.45,CBGn->Data());

  lCB->SetNDC(kTRUE);
  lCBalpha->SetNDC(kTRUE);
  lCBn->SetNDC(kTRUE);
  lCBG->SetNDC(kTRUE);
  lCBGalpha->SetNDC(kTRUE);
  lCBGn->SetNDC(kTRUE);

  lCB->Draw();
  lCBalpha->Draw();
  lCBn->Draw();
  lCBG->Draw();
  lCBGalpha->Draw();
  lCBGn->Draw();

  TString outfname;
  outfname = Form("Jpsi_pp_MCshape_Rap_%3.1f-%3.1f_Pt_%3.1f-%3.1f.pdf",ymin,ymax,ptmin,ptmax);

  std::cout << outfname << std::endl;
  if (savePlot)
    c1->SaveAs(outfname);

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
