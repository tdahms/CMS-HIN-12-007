#include <iostream>
#include <sstream>
#include <string>

#include <TROOT.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
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
#include <RooHistPdfConv.h>
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

using namespace RooFit;
using namespace std;

bool superImpose = true;
bool analyticBlifetime = true;
bool narrowSideband = false;
bool oneGaussianResol = false;

void getOptRange(string &ran,float *min,float *max);
void setWSRange(RooWorkspace *ws, float lmin, float lmax, float errmin, float errmax);
RooBinning setCtBinning(float lmin,float lmax);
void defineMassBkg(RooWorkspace *ws);
void defineMassSig(RooWorkspace *ws);
void getMCTrueLifetime(RooWorkspace *ws, RooDataSet *redMCCutNP, float *bgmcVal, float *bctauVal, string titlestr);
void defineCTResol(RooWorkspace *ws);
void defineCTBkg(RooWorkspace *ws);
void defineCTSig(RooWorkspace *ws, RooDataSet *redMCCutNP, string titlestr);
RooDataHist* subtractSidebands(RooWorkspace* ws, RooDataHist* subtrData, RooDataHist* all, RooDataHist* side, float scalefactor, string varName);


int main (int argc, char* argv[]) {

  gROOT->Macro("/afs/cern.ch/user/m/miheejo/public/JpsiV2/JpsiStyle.C");
  string FileName, FileNameMC1, FileNameMC2;
  string mBkgFunct, mSigFunct;
  bool prefitMass = false;
  bool prefitSignalCTau = false;
  bool prefitBkg = false;
  int  isGG = 0;
  bool fracfix = true;
  int isMB = false;
  int isMBCtau = false;
  bool isPT = false;
  string prange, lrange, yrange, crange, phirange, errrange;
  string dirPre;
  string rpmethod = "etHF";

  // *** Check options
  for (int i=1; i<argc; i++) {
    char *tmpargv = argv[i];
    switch (tmpargv[0]) {
      case '-':{
        switch (tmpargv[1]) {
          case 'f':
            FileName = argv[i+1];
            cout << "Fitted data file: " << FileName << endl;
            break;
          case 'm':
            FileNameMC1 = argv[i+1];
            FileNameMC2 = argv[i+2];
            cout << "MC data file 1: " << FileNameMC1 << endl;
            cout << "MC data file 2: " << FileNameMC2 << endl;
            break;
          case 'v':
            mSigFunct = argv[i+1];
            mBkgFunct = argv[i+2];
            cout << "Mass signal function: " << mSigFunct << endl;
            cout << "Mass background function: " << mBkgFunct << endl;
            break;
          case 'r':
            rpmethod = argv[i+1];
            cout << "Reaction plane: " << rpmethod << endl;
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
            cout << "Rapidity range: " << yrange << " rad" << endl;
            break;
          case 't':
            crange = argv[i+1];
            cout << "Centrality range: " << crange << " %" << endl;
            break;
          case 'l':
            lrange = argv[i+1];
            cout << "l(J/psi) range: " << lrange << " mm" << endl; 
            break;
          case 'e':
            errrange = argv[i+1];
            cout << "Range for sigma_l(J/psi) is " << errrange << " mm" << endl;
            break;    
          case 's':
            phirange = argv[i+1];
            cout << "dPhi(J/psi) range: " << phirange << " rad" << endl;
            break;
          case 'a':
            if (atoi(argv[i+1]) == 0) {
              analyticBlifetime = false;
              cout << "Turn off: RooHistPdf from MC template of J/psi Ctau lifetime will be used" << endl;
            } else {
              analyticBlifetime = true;
              cout << "Turn on: Analytical MC J/psi Ctau lifetime PDF will be used" << endl;
            }
            break;
          case 'c':
            prefitSignalCTau = true;
            cout << "Turn on: prompt signal ctau distribution pre-fitting" << endl;
            break;
          case 'b':
            prefitBkg = true;
            cout << "Turn on: Background(=sideband data) ctau distribution pre-fitting" << endl;
            break;
          case 'x':
            isMB = atoi(argv[i+1]);
            cout << "Inclusive fit option: " << isMB << endl;
            break;
          case 'z':
            if (0 == atoi(argv[i+1])) {
              fracfix = true;
              cout << "Fix frac2,3 BEFORE ctau bkg fitting" << endl;
            } else {
              fracfix = false;
              cout << "Fix frac2,3 AFTER ctau bkg fitting" << endl;
            }
            break;
        }
      }
    }
  }// End check options
 
  float pmin=0, pmax=0, ymin=0, ymax=0, lmin=0, lmax=0, cmin=0, cmax=0, psmax=0, psmin=0, errmin=0, errmax=0;
  getOptRange(prange,&pmin,&pmax);
  getOptRange(lrange,&lmin,&lmax);
  getOptRange(errrange,&errmin,&errmax);
  getOptRange(crange,&cmin,&cmax);
  getOptRange(yrange,&ymin,&ymax);
  getOptRange(phirange,&psmin,&psmax);


  // *** TFile for saving fitting results
  string resultFN;
  resultFN = dirPre + "_rap" + yrange + "_pT" + prange + "_cent" + crange + "_dPhi" + phirange + "_fitResult.root";
  TFile resultF(resultFN.c_str(),"RECREATE");

  // *** Read MC and Data files
  TFile fInMC(FileNameMC1.c_str());   //Non-prompt J/psi MC
  cout << FileNameMC1.c_str() << endl;
  if (fInMC.IsZombie()) { cout << "CANNOT open MC1 root file\n"; return 1; }
  fInMC.cd();
  RooDataSet *dataMC = (RooDataSet*)fInMC.Get("dataJpsi");
  dataMC->SetName("dataMC");

  TFile fInMC2(FileNameMC2.c_str());  //Prompt J/psi MC
  cout << FileNameMC2.c_str() << endl;
  if (fInMC2.IsZombie()) { cout << "CANNOT open MC2 root file\n"; return 1; }
  fInMC2.cd();
  RooDataSet *dataMC2 = (RooDataSet*)fInMC2.Get("dataJpsi");
  dataMC2->SetName("dataMC2");

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
  sprintf(reduceDS,"Jpsi_Pt>%.2f && Jpsi_Pt<%.2f && abs(Jpsi_Y)>%.2f && abs(Jpsi_Y)<%.2f && Jpsi_CtErr > %.2f && Jpsi_CtErr < %.2f",pmin,pmax,ymin,ymax,errmin,errmax);
  cout << "reduceDS for MC and data: " << reduceDS << endl;

  RooDataSet *redMC = (RooDataSet*)dataMC->reduce(reduceDS);
  ws->import(*redMC);
  RooDataSet *redMC2 = (RooDataSet*)dataMC2->reduce(reduceDS);
  ws->import(*redMC2);
  RooDataSet *redData = (RooDataSet*)data->reduce(reduceDS);
  ws->import(*redData);
 
  setWSRange(ws, lmin, lmax, errmin, errmax);

  // Draw data
  ws->var("Jpsi_Ct")->SetTitle("#font[12]{l}_{J/#psi}");
//  ws->var("Jpsi_CtTrue")->setBins(2000);

  // Test true lifetimes
  RooPlot *trueFrame = ws->var("Jpsi_CtTrue")->frame();
  ws->data("dataMC")->plotOn(trueFrame,DataError(RooAbsData::SumW2),Cut("MCType==MCType::NP"));

  string titlestr;
  TCanvas c0;
  c0.cd(); trueFrame->Draw();
  titlestr = dirPre + "_rap" + yrange + "_pT" + prange + "_cent" + crange + "_dPhi" + phirange + "_testTrueLife_Lin.pdf";
  c0.SaveAs(titlestr.c_str());

  // Define binning for true lifetime
  ws->var("Jpsi_Mass")->setBins(60);
//  ws->var("Jpsi_Mass")->setBins(133);
  ws->var("Jpsi_CtErr")->setBins(25);
  if (pmin > 40.) ws->var("Jpsi_CtErr")->setBins(8);

  // define binning for true lifetime
  RooBinning rb(-0.1,4.0);
  rb.addUniform(5,-0.1,0.0);
  rb.addUniform(100,0.0,0.5);
  rb.addUniform(15,0.5,1.0);
  rb.addUniform(20,1.0,2.5);
  rb.addUniform(5,2.5,4.0);
  if (analyticBlifetime) {
    ws->var("Jpsi_CtTrue")->setBins(200);
  } else {
    ws->var("Jpsi_CtTrue")->setBinning(rb);
  }

  // Define binning for lifetime
  RooBinning rb2 = setCtBinning(lmin,lmax);
  ws->var("Jpsi_Ct")->setBinning(rb2);
  
  // Define ctau binning for plotting (coarser bin)
  RooBinning rb3(-lmin,lmax);
  rb3.addBoundary(-1.0);
  rb3.addBoundary(-0.7);
  rb3.addBoundary(-0.6);
  rb3.addBoundary(-0.5);
  rb3.addUniform(5,-0.5,-0.2);
  rb3.addUniform(15,-0.2,0.2);
  rb3.addUniform(5,0.2,0.5);
  rb3.addUniform(5,0.5,1.0);

  // Additional cuts on data and get sub-datasets/histograms
  RooDataSet *redDataCut;
  string reduceDSstr;
  if (isGG == 0) {
    reduceDSstr = "Jpsi_Type == Jpsi_Type::GG &&\
                  (MCType != MCType::NP || Jpsi_CtTrue>0.0001) &&\
                  (MCType == MCType::PR || MCType == MCType::NP)";
    redDataCut = (RooDataSet*)redData->reduce("Jpsi_Type == Jpsi_Type::GG");
  } else if (isGG == 1) {
    reduceDSstr = "(Jpsi_Type == Jpsi_Type::GT || Jpsi_Type == Jpsi_Type::GG) &&\
                  (MCType != MCType::NP || Jpsi_CtTrue > 0.0001) &&\
                  (MCType == MCType::PR || MCType == MCType::NP)";
    redDataCut = (RooDataSet*)redData->reduce("Jpsi_Type == Jpsi_Type::GT || Jpsi_Type == Jpsi_Type::GG");
  } else {
    reduceDSstr = "(MCType != MCType::NP || Jpsi_CtTrue>0.0001) &&\
                   (MCType == MCType::PR || MCType == MCType::NP)";
    redDataCut = (RooDataSet*)redData->reduce("Jpsi_Ct < 600000.");
  }

  RooDataHist *binData = new RooDataHist("binData","binData",RooArgSet( *(ws->var("Jpsi_Mass")),*(ws->var("Jpsi_Ct")),*(ws->var("Jpsi_CtErr")) ), *redDataCut);
  RooDataHist *binDataCtErr = new RooDataHist("binDataCtErr","binDataCtErr",RooArgSet(*(ws->var("Jpsi_CtErr"))),*redDataCut);
  cout << "DATA :: N events to fit: " << binData->sumEntries() << endl;

  // *** Get MC sub-datasets and its histograms corresponds to data
  RooDataSet *redMCCut = (RooDataSet*) redMC->reduce(reduceDSstr.c_str());
  RooDataSet *redMCCutNP = (RooDataSet*) redMCCut->reduce(RooArgSet(*(ws->var("Jpsi_CtTrue"))),"MCType == MCType::NP");
  RooDataSet *redMCCut2 = (RooDataSet*) redMC2->reduce(reduceDSstr.c_str());
  RooDataSet *redMCCutPR = (RooDataSet*) redMCCut2->reduce("MCType == MCType::PR");

  // SYSTEMATICS 1 (very sidebands)
  RooDataSet *redDataSB;
  if (narrowSideband) {
   redDataSB = (RooDataSet*) redDataCut->reduce("Jpsi_Mass<2.8 || Jpsi_Mass>3.4");
  } else {
   redDataSB = (RooDataSet*) redDataCut->reduce("Jpsi_Mass<2.9 || Jpsi_Mass>3.3");
  }
  RooDataHist *binDataSB = new RooDataHist("binDataSB","Data distribution for background",RooArgSet( *(ws->var("Jpsi_Mass")),*(ws->var("Jpsi_Ct")) ),*redDataSB);
  RooDataSet *redDataSIG = (RooDataSet*)redDataCut->reduce("Jpsi_Mass > 2.9 && Jpsi_Mass < 3.3");

  RooDataHist *binMCPR = new RooDataHist("binMCPR","MC distribution for PR signal",RooArgSet( *(ws->var("Jpsi_Mass")),*(ws->var("Jpsi_Ct")) ),*redMCCutPR);
  cout << "PRMC :: N events to fit: " << binMCPR->sumEntries() << endl;
  RooDataHist *binMCNP = new RooDataHist("binMCNP","MC distribution for NP signal",RooArgSet(*(ws->var("Jpsi_CtTrue"))),*redMCCutNP);
  cout << "NPMC :: N events to fit: " << binMCNP->sumEntries() << endl;

  RooDataHist *binDataCtErrSB = new RooDataHist("binDataCtErrSB","Data ct error distribution for bkg",RooArgSet(*(ws->var("Jpsi_CtErr"))),*redDataSB);
  RooDataHist *binDataCtErrSIG = new RooDataHist("binDataCtErrSIG","Data ct error distribution for sig",RooArgSet(*(ws->var("Jpsi_CtErr"))),*redDataSIG);
  RooHistPdf errPdfBkg("errPdfBkg","Error PDF bkg",RooArgSet(*(ws->var("Jpsi_CtErr"))),*binDataCtErrSB);  ws->import(errPdfBkg);

  // ** Test Ct error distribution on the sideband region
  RooPlot *errframe2 = ws->var("Jpsi_CtErr")->frame();
  binDataCtErrSB->plotOn(errframe2,DataError(RooAbsData::SumW2));
  ws->pdf("errPdfBkg")->plotOn(errframe2,LineColor(kBlue),Normalization(binDataCtErrSB->sumEntries(),RooAbsReal::NumEvent));

  c0.Clear(); c0.cd(); c0.SetLogy(0); errframe2->Draw();
  titlestr = dirPre + "_rap" + yrange + "_pT" + prange + "_cent" + crange  + "_dPhi" + phirange + "_testErrPdfBkg_Lin.pdf";
  c0.SaveAs(titlestr.c_str());
  c0.SetLogy(1); errframe2->Draw();
  titlestr = dirPre + "_rap" + yrange + "_pT" + prange + "_cent" + crange + "_dPhi" + phirange + "_testErrPdfBkg_Log.pdf";
  c0.SaveAs(titlestr.c_str());


  // *** Define PDFs with parameters (mass and ctau)
  // J/psi mass parameterization
  defineMassBkg(ws);
  defineMassSig(ws);
  // J/psi CTau parameterization
  defineCTResol(ws);              // R(l) : resolution function
  defineCTBkg(ws);                // theta(l') convolution R(l')
  titlestr = dirPre + "_rap" + yrange + "_pT" + prange + "_cent" + crange + "_dPhi" + phirange + "_testTrueLifeFit";
  defineCTSig(ws,redMCCutNP,titlestr); // F_B(l) : R(l') convolution X_mc(l')

  RooProdPdf bkgCtauTOT_PEE("bkgCtauTOT_PEE","PDF with PEE", *(ws->pdf("errPdfBkg")),
                            Conditional(*(ws->pdf("bkgCtTot")),RooArgList(*(ws->var("Jpsi_Ct"))))
                           );  ws->import(bkgCtauTOT_PEE);  
  char funct[100];
  sprintf(funct,"PROD::totBKG(%s,bkgCtTot)",mBkgFunct.c_str());
  ws->factory(funct); // F_bkg(l) : exp*theta(l')

  string partTit, partFile;
  if (isGG == 0) { partTit = "glb-glb"; partFile = "GG"; }
  else if (isGG = 1) { partTit = "glb-trk"; partFile = "GT"; }
  else { partTit = "all"; partFile = "ALL"; }

  // Binning for invariant mass distribution
  RooBinning rbm(2.6,3.5);
  rbm.addUniform(45,2.6,3.5);

  // Global TLatex, TH1, TGraph objects for drawing
  TLatex *t = new TLatex();
  t->SetNDC(); t->SetTextAlign(12);

  Double_t fx[2], fy[2], fex[2], fey[2];
  TGraphErrors *gfake1 = new TGraphErrors(2,fx,fy,fex,fey);
  gfake1->SetMarkerStyle(20); gfake1->SetMarkerSize(1);
  TH1F hfake11 = TH1F("hfake1","hfake1",100,200,300);
  hfake11.SetLineColor(kBlue); hfake11.SetLineWidth(4); hfake11.SetLineStyle(7); hfake11.SetFillColor(kAzure-9); hfake11.SetFillStyle(1001);
  TH1F hfake21 = TH1F("hfake2","hfake2",100,200,300);
  hfake21.SetLineColor(kBlack); hfake21.SetLineWidth(4); hfake21.SetFillColor(kBlack); hfake21.SetFillStyle(3354);
  TH1F hfake31 = TH1F("hfake3","hfake3",100,200,300);
  hfake31.SetLineColor(kRed); hfake31.SetMarkerStyle(kCircle); hfake31.SetLineWidth(4); hfake31.SetMarkerColor(kRed); hfake31.SetLineStyle(9); hfake31.SetFillColor(kRed-7); hfake31.SetFillStyle(3444);


  // Set some fitting variables to constant. It depends on the prefitting options.
  RooFitResult *fitM;
  if (prefitMass) {
    struct PARAM {
      double coeffGaus; double coeffGausErr;
      double meanSig1;  double meanSig1Err;
      double sigmaSig1; double sigmaSig1Err;
      double sigmaSig2; double sigmaSig2Err;
      double alpha;     double alphaErr;
      double enne;      double enneErr;
      double enneW;     double enneWErr;
    };

    bool centConst = false;  //False: fit w/o any constrained parameters (centrality dep.)
    bool dPhiConst = false;  //False: fit w/o any constrained parameters (dPhi dep.)
    double inputN[2] = {0};  //Number of yield/background in the 0-1.571 rad bin
    if (isMB != 0) {
      if (isMB != 4 && (mSigFunct.compare("sigCB2WNG1") || mBkgFunct.compare("expFunct"))) {
        cout << "For fit systematics:\n";
        cout << "1) Signal shape check: should use default runOpt (4)\n";
        cout << "2) Background shape check: should use defatul runOpt (4)\n";
        cout << "3) Constrain fit: should use sigCB2WNG1 for signal shape, expFunct for bkg shape with runOpt (3)\n";
        return -1;
      }

      string inputFN, inputFNcb;
      bool centest = true;  //False: fit w/o any fixed parameters (centrality dep.)
      bool dPhitest = true;  //False: fit w/o any fixed parameters (dPhi dep.)
      if (isMB == 1) {
        // Except the *_dPhi0.00-1.57 bins, all other delta Phi bins read/set signal parameters of this bin
        centest = false;
        centConst = false;
        dPhiConst = false;
        if (!phirange.compare("0.000-1.571")) {
          dPhitest = false;
        } else {
          dPhitest = true;
          inputFN =  dirPre + "_rap" + yrange + "_pT" + prange + "_cent" + crange + "_dPhi0.000-1.571.txt";
        }

      } else if (isMB == 2) {
        // Except the *_cent0-100_dPhi0.000-1.571 bin, all other bins read/set signal parameters of this bin
        centest = false;
        centConst = false;
        dPhiConst = false;
        if (!crange.compare("0-100") && !phirange.compare("0.000-1.571")) {
          dPhitest = false;
        } else {
          dPhitest = true;
          inputFN =  dirPre + "_rap" + yrange + "_pT" + prange + "_cent0-100_dPhi0.000-1.571.txt";
        }

      } else if (isMB == 3) {
        // Constrain signal parameters that were fixed in the default method.
        // Constrained mean values are coming from *_cent0-100_dPhi0.000-1.571 bin result.
        centest = false;
        dPhitest = false;
        if (!crange.compare("0-100")) {
          centConst = false;
        } else {
          centConst = true;
          inputFNcb =  dirPre + "_rap" + yrange + "_pT" + prange + "_cent0-100_dPhi0.000-1.571.txt";
        }
        if (!phirange.compare("0.000-1.571")) {
          dPhiConst = false;
        } else {
          dPhiConst = true;
          inputFN =  dirPre + "_rap" + yrange + "_pT" + prange + "_cent" + crange + "_dPhi0.000-1.571.txt";
        }

      } else if (isMB == 4) {
        centConst = false;
        dPhiConst = false;
        if (!crange.compare("0-100")) {
          centest = false;
        } else {
          centest = true;
          inputFNcb =  dirPre + "_rap" + yrange + "_pT" + prange + "_cent0-100_dPhi0.000-1.571.txt";
        }
        if (!phirange.compare("0.000-1.571")) {
          dPhitest = false;
        } else {
          dPhitest = true;
          inputFN =  dirPre + "_rap" + yrange + "_pT" + prange + "_cent" + crange + "_dPhi0.000-1.571.txt";
        }
      }

      if (centConst || dPhiConst) {
        ifstream input;
        if (centConst && !dPhiConst) input.open(inputFNcb.c_str());
        else if (dPhiConst) input.open(inputFN.c_str());
        if (!input.good()) { cout << "Failed to open: " <<inputFNcb << endl; return 1; }
        string tmp;
        double inputTmp[2] = {0};
        PARAM inputP;
        input >> tmp >> inputTmp[0] >> inputTmp[1]; //NSig
        inputN[0] = inputTmp[0];  //NSig
        input >> tmp >> inputTmp[0] >> inputTmp[1]; //NBkg
        inputN[1] = inputTmp[0];  //NBkg
        for (int p=0; p<9; p++) {   //Mass signal parameters
          input >> tmp >> inputTmp[0] >> inputTmp[1];
          cout << tmp << " " << inputTmp[0] << " " << inputTmp[1]<< endl;
          if (!tmp.compare("coeffGaus")) { 
            inputP.coeffGaus = inputTmp[0]; inputP.coeffGausErr = inputTmp[1];
          } else if (!tmp.compare("meanSig1")) {
            inputP.meanSig1 = inputTmp[0];  inputP.meanSig1Err = inputTmp[1];
          } else if (!tmp.compare("sigmaSig1")) {
            inputP.sigmaSig1 = inputTmp[0]; inputP.sigmaSig1Err = inputTmp[1];
          } else if (!tmp.compare("sigmaSig2")) {
            inputP.sigmaSig2 = inputTmp[0]; inputP.sigmaSig2Err = inputTmp[1];
          } else if (!tmp.compare("alpha")) {
            inputP.alpha = inputTmp[0];     inputP.alphaErr = inputTmp[1];
          } else if (!tmp.compare("enne")) {
            inputP.enne = inputTmp[0];      inputP.enneErr = inputTmp[1];
          } else if (!tmp.compare("enneW")) {
            inputP.enneW = inputTmp[0];     inputP.enneWErr = inputTmp[1];
          }
        }

        char confunct[1000]={0};
        if (dPhiConst) {  //Constrain CB width only for dPhi bins
          sprintf(confunct,"Gaussian::sigmaSig2Con(sigmaSig2,%f,%f)",inputP.sigmaSig2,inputP.sigmaSig2Err);
          ws->factory(confunct);
        }
        sprintf(confunct,"Gaussian::sigmaSig1Con(sigmaSig1,%f,%f)",inputP.sigmaSig1,inputP.sigmaSig1Err);
        ws->factory(confunct);
        sprintf(confunct,"Gaussian::meanSig1Con(meanSig1,%f,%f)",inputP.meanSig1,inputP.meanSig1Err);
        ws->factory(confunct);
        sprintf(confunct,"Gaussian::coeffGausCon(coeffGaus,%f,%f)",inputP.coeffGaus,inputP.coeffGausErr);
        ws->factory(confunct);
        sprintf(confunct,"Gaussian::alphaCon(alpha,%f,%f)",inputP.alpha,inputP.alphaErr);
        ws->factory(confunct);
        sprintf(confunct,"Gaussian::enneCon(enne,%f,%f)",inputP.enne,inputP.enneErr);
        ws->factory(confunct);
        sprintf(confunct,"Gaussian::enneWCon(enneW,%f,%f)",inputP.enneW,inputP.enneWErr);
        ws->factory(confunct);

        input.close();
      }

      if (centest) {
        ifstream input(inputFNcb.c_str());
        if (!input.good()) { cout << "Failed to open: " <<inputFNcb << endl; return 1; }
        string tmp;
        double inputNS[2] = {0};
        double inputTmp[2] = {0};
        PARAM inputP;
        input >> tmp >> inputNS[0] >> inputNS[1]; //NSig
        input >> tmp >> inputTmp[0] >> inputTmp[1]; //NBkg
        for (int p=0; p<9; p++) {   //Mass signal parameters
          input >> tmp >> inputTmp[0] >> inputTmp[1];
          cout << tmp << " " << inputTmp[0] << endl;
          if (!tmp.compare("coeffGaus")) inputP.coeffGaus = inputTmp[0];
          else if (!tmp.compare("meanSig1")) inputP.meanSig1 = inputTmp[0];
          else if (!tmp.compare("sigmaSig1")) inputP.sigmaSig1 = inputTmp[0];
          else if (!tmp.compare("sigmaSig2")) inputP.sigmaSig2 = inputTmp[0];
          else if (!tmp.compare("alpha")) inputP.alpha = inputTmp[0];
          else if (!tmp.compare("enne")) inputP.enne = inputTmp[0];
          else if (!tmp.compare("enneW")) inputP.enneW = inputTmp[0];
        }

//        ws->var("sigmaSig2")->setVal(inputP.sigmaSig2);
        ws->var("coeffGaus")->setVal(inputP.coeffGaus);
        ws->var("meanSig1")->setVal(inputP.meanSig1);
        ws->var("sigmaSig1")->setVal(inputP.sigmaSig1);
        ws->var("alpha")->setVal(inputP.alpha);
        ws->var("enne")->setVal(inputP.enne);
        ws->var("enneW")->setVal(inputP.enneW);

//        ws->var("sigmaSig2")->setConstant(kTRUE);
        ws->var("coeffGaus")->setConstant(kTRUE);
        ws->var("meanSig1")->setConstant(kTRUE);
        ws->var("sigmaSig1")->setConstant(kTRUE);
        ws->var("alpha")->setConstant(kTRUE);
        ws->var("enne")->setConstant(kTRUE);
        ws->var("enneW")->setConstant(kTRUE);

        input.close();
      }

      if (dPhitest) {
        ifstream input(inputFN.c_str());
        if (!input.good()) { cout <<"Failed to open: " <<  inputFN << endl; return 1; }
        string tmp;
        PARAM inputP;
        double inputNS[2] = {0};
        double inputTmp[2]={0};
        input >> tmp >> inputNS[0] >> inputNS[1]; //NSig
        input >> tmp >> inputTmp[0] >> inputTmp[1]; //NBkg
        for (int p=0; p<9; p++) {   //Mass signal parameters
          input >> tmp >> inputTmp[0] >> inputTmp[1];
          cout << tmp << " " << inputTmp[0] << endl;
          if (!tmp.compare("coeffGaus")) inputP.coeffGaus = inputTmp[0];
          else if (!tmp.compare("meanSig1")) inputP.meanSig1 = inputTmp[0];
          else if (!tmp.compare("sigmaSig1")) inputP.sigmaSig1 = inputTmp[0];
          else if (!tmp.compare("sigmaSig2")) inputP.sigmaSig2 = inputTmp[0];
          else if (!tmp.compare("alpha")) inputP.alpha = inputTmp[0];
          else if (!tmp.compare("enne")) inputP.enne = inputTmp[0];
          else if (!tmp.compare("enneW")) inputP.enneW = inputTmp[0];
        }

        ws->var("sigmaSig2")->setVal(inputP.sigmaSig2);
        ws->var("coeffGaus")->setVal(inputP.coeffGaus);
        ws->var("meanSig1")->setVal(inputP.meanSig1);
        ws->var("sigmaSig1")->setVal(inputP.sigmaSig1);
        ws->var("alpha")->setVal(inputP.alpha);
        ws->var("enne")->setVal(inputP.enne);
        ws->var("enneW")->setVal(inputP.enneW);

        ws->var("sigmaSig2")->setConstant(kTRUE);
        ws->var("coeffGaus")->setConstant(kTRUE);
        ws->var("meanSig1")->setConstant(kTRUE);
        ws->var("sigmaSig1")->setConstant(kTRUE);
        ws->var("alpha")->setConstant(kTRUE);
        ws->var("enne")->setConstant(kTRUE);
        ws->var("enneW")->setConstant(kTRUE);

        input.close();
      }

    } //End of isMB != 0
    
    if (!yrange.compare("0.0-1.2")) {
      sprintf(funct,"SUM::sigMassPDF(NSig[4000.0,1.0,50000.0]*%s,NBkg[2000.0,1.0,500000.0]*%s)",mSigFunct.c_str(),mBkgFunct.c_str());
    } else if (!yrange.compare("1.6-2.4")) {
      sprintf(funct,"SUM::sigMassPDF(NSig[1000.0,1.0,50000.0]*%s,NBkg[25000.0,1.0,500000.0]*%s)",mSigFunct.c_str(),mBkgFunct.c_str());
    } else {
      sprintf(funct,"SUM::sigMassPDF(NSig[5000.0,1.0,50000.0]*%s,NBkg[6000.0,1.0,500000.0]*%s)",mSigFunct.c_str(),mBkgFunct.c_str());
    }

/*    if (!yrange.compare("0.0-1.2")) {
        sprintf(funct,"SUM::sigMassPDF(NSig[3000.0,1.0,50000.0]*%s,NBkg[1500.0,1.0,500000.0]*%s)",mSigFunct.c_str(),mBkgFunct.c_str());
    } else if (!yrange.compare("1.6-2.4")) {
      if (pmin > 6.5 )
        sprintf(funct,"SUM::sigMassPDF(NSig[1000.0,1.0,50000.0]*%s,NBkg[1500.0,1.0,500000.0]*%s)",mSigFunct.c_str(),mBkgFunct.c_str());
      else
        sprintf(funct,"SUM::sigMassPDF(NSig[1000.0,1.0,50000.0]*%s,NBkg[25000.0,1.0,500000.0]*%s)",mSigFunct.c_str(),mBkgFunct.c_str());
    } else if (!yrange.compare("0.0-2.4")) {
      if (pmin < 10.0)
        sprintf(funct,"SUM::sigMassPDF(NSig[2000.0,1.0,50000.0]*%s,NBkg[1000.0,1.0,500000.0]*%s)",mSigFunct.c_str(),mBkgFunct.c_str());
      else if (pmin == 10.0 && pmax == 13.0)
        sprintf(funct,"SUM::sigMassPDF(NSig[900.0,1.0,50000.0]*%s,NBkg[360.0,1.0,500000.0]*%s)",mSigFunct.c_str(),mBkgFunct.c_str());
      else
        sprintf(funct,"SUM::sigMassPDF(NSig[1500.0,1.0,50000.0]*%s,NBkg[1000.0,1.0,500000.0]*%s)",mSigFunct.c_str(),mBkgFunct.c_str());
    }
*/
    ws->factory(funct);

    if (dPhiConst) { //sigmaSig2 will be constrained too!
      fitM = ws->pdf("sigMassPDF")->fitTo(*redDataCut,ExternalConstraints(RooArgSet(*(ws->pdf("sigmaSig2Con")),*(ws->pdf("sigmaSig1Con")),*(ws->pdf("meanSig1Con")),*(ws->pdf("coeffGausCon")),*(ws->pdf("alphaCon")),*(ws->pdf("enneWCon")))),Extended(1),Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(6));
    } else if (centConst && !dPhiConst) { //sigmaSig2 will be NOT constrained!
      fitM = ws->pdf("sigMassPDF")->fitTo(*redDataCut,ExternalConstraints(RooArgSet(*(ws->pdf("sigmaSig1Con")),*(ws->pdf("meanSig1Con")),*(ws->pdf("coeffGausCon")),*(ws->pdf("alphaCon")),*(ws->pdf("enneWCon")))),Extended(1),Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(6));
    } else { // all free fit bin
      fitM = ws->pdf("sigMassPDF")->fitTo(*redDataCut,Extended(1),Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(6));
    }
    resultF.cd(); cout << "fitM->Write(sigMassPDF):" << fitM->Write("sigMassPDF") << endl;
 

    // *** Draw mass plot before do ctau fit
    RooPlot *mframe_wob = ws->var("Jpsi_Mass")->frame();
    redDataCut->plotOn(mframe_wob,DataError(RooAbsData::SumW2),XErrorSize(0),MarkerSize(1),Binning(rbm));
    titlestr = "2D fit for" + partTit + "muons (mass projection), p_{T} = " + prange  + "_dPhi" + phirange + " GeV/c and |y| = " + yrange;
    mframe_wob->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c^{2})");
    mframe_wob->GetXaxis()->CenterTitle(1);
    const double max = mframe_wob->GetMaximum() * 1.3;
    mframe_wob->SetMaximum(max);
    ws->pdf("sigMassPDF")->plotOn(mframe_wob,DrawOption("F"),FillColor(kBlack),FillStyle(3354),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));

    ws->pdf("sigMassPDF")->plotOn(mframe_wob,Components(mBkgFunct.c_str()),DrawOption("F"),FillColor(kAzure-9),FillStyle(1001),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
    
    ws->pdf("sigMassPDF")->plotOn(mframe_wob,Components(mBkgFunct.c_str()),LineColor(kBlue),LineStyle(7),LineWidth(5),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
    ws->pdf("sigMassPDF")->plotOn(mframe_wob,LineColor(kBlack),LineWidth(2),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
    redDataCut->plotOn(mframe_wob,DataError(RooAbsData::SumW2),XErrorSize(0),MarkerSize(1),Binning(rbm));

    TH1 *hdata = redDataCut->createHistogram("hdata",*ws->var("Jpsi_Mass"),Binning(rbm));
    // *** Calculate chi2/nDof for mass fitting
    int nBins = hdata->GetNbinsX();
    RooHist *hpull; hpull = mframe_wob->pullHist(); hpull->SetName("hpullhist");
    double Chi2 = 0;
    int nFullBinsPull = 0;
    double *ypull = hpull->GetY();
    for (unsigned int i=0; i < nBins; i++) {
      if (hdata->GetBinContent(i) == 0) continue;
      nFullBinsPull++;
      Chi2 = Chi2 + pow(ypull[i],2);
    }
    double UnNormChi2 = Chi2;
    int nFitParam = fitM->floatParsFinal().getSize();
    int Dof = nFullBinsPull - nFitParam;
    Chi2 /= (nFullBinsPull - nFitParam);

    // *** Check in narrower signal region NSig
    ws->var("Jpsi_Mass")->setRange("sigpeak",2.9,3.3);
    RooAbsReal *inteAll = ws->pdf(mSigFunct.c_str())->createIntegral(RooArgSet(*ws->var("Jpsi_Mass")),NormSet(RooArgSet(*ws->var("Jpsi_Mass"))));
    RooAbsReal *inteSig = ws->pdf(mSigFunct.c_str())->createIntegral(RooArgSet(*ws->var("Jpsi_Mass")),Range("sigpeak"),NormSet(RooArgSet(*ws->var("Jpsi_Mass"))),Range("sigpeak"));

    TCanvas c00; c00.cd(); mframe_wob->Draw();
    t->SetTextSize(0.05);
    t->DrawLatex(0.17,0.90,"CMS Preliminary");
    t->DrawLatex(0.17,0.82,"PbPb  #sqrt{s_{NN}} = 2.76 TeV");
    t->SetTextSize(0.04);
    t->DrawLatex(0.17,0.75,"L_{int} =  150 #mub^{-1}");
    t->SetTextSize(0.035);
    //    sprintf(reduceDS,"%.0f-%.0f%, %0.1f < |y| < %.1f",cmin,cmax,ymin,ymax);
    t->DrawLatex(0.17,0.69,reduceDS);
    //    sprintf(reduceDS,"%.1f < p_{T} < %.1f GeV/c",pmin,pmax);
    t->DrawLatex(0.17,0.64,reduceDS);
    //    sprintf(reduceDS,"%.2f < |#phi_{J/#psi}-#Psi_{RP}| < %.2f",psmin,psmax);
    t->DrawLatex(0.17,0.59,reduceDS);
    //    sprintf(reduceDS,"EP: %s",rpmethod.c_str());
    t->DrawLatex(0.17,0.54,reduceDS);
    t->SetTextSize(0.04);
    //    sprintf(reduceDS,"#chi^{2}/dof = %0.1f/%d",UnNormChi2,Dof);
    t->DrawLatex(0.65,0.85,reduceDS);
    //    sprintf(reduceDS,"yield: %0.0f",ws->var("NSig")->getVal());
    t->DrawLatex(0.65,0.80,reduceDS);
    //    sprintf(reduceDS,"#sigma = %0.0f MeV/c^{2}",ws->var("sigmaSig2")->getVal()*1000);
    t->DrawLatex(0.65,0.75,reduceDS);

    TLegend leg1(0.63,0.6,0.92,0.73,NULL,"brNDC");
    leg1.SetFillStyle(0); leg1.SetBorderSize(0); leg1.SetShadowColor(0); leg1.SetMargin(0.2);
    leg1.AddEntry(gfake1,"data","p");
    leg1.AddEntry(&hfake21,"total fit","lf");
    leg1.AddEntry(&hfake11,"background","lf");
    leg1.Draw("same");
    titlestr = dirPre + "_rap" + yrange  + "_pT" + prange + "_cent" + crange + "_dPhi" + phirange + "_massfit_wob.pdf";
    c00.SaveAs(titlestr.c_str());

  } else {
    RooRealVar NSig("NSig","dummy total signal events",0.);
    ws->import(NSig);
  }

  Double_t NSig_fin = ws->var("NSig")->getVal();
  Double_t ErrNSig_fin = ws->var("NSig")->getError();
  Double_t NBkg_fin = ws->var("NBkg")->getVal();
  Double_t ErrNBkg_fin = ws->var("NBkg")->getError();

//  if (false) {  // skip ctau fitting
  // *** scaleF is defined to scale down ct err dist in 2.9-3.3 GeV/c2 signal region
  float bc;
  if (!mBkgFunct.compare("expFunct")) bc = ws->var("coefExp")->getVal();
  else if (!mBkgFunct.compare("polFunct")) bc = ws->var("coefPol1")->getVal();
  float scaleF = (exp(2.9*bc)-exp(3.3*bc))/(exp(2.6*bc)-exp(2.9*bc)+exp(3.3*bc)-exp(3.5*bc));
  RooDataHist* subtrData = new RooDataHist("subtrData","Subtracted data",RooArgSet(*(ws->var("Jpsi_CtErr")))); 
  RooDataHist* scaledBkg = subtractSidebands(ws,subtrData,binDataCtErrSIG,binDataCtErrSB,scaleF,"Jpsi_CtErr");
  subtrData->SetName("subtrData");
  RooHistPdf errPdfSig("errPdfSig","Error PDF signal",RooArgSet(*(ws->var("Jpsi_CtErr"))),*subtrData);  ws->import(errPdfSig);

  if (prefitMass) {
    ws->var("alpha")->setConstant(kTRUE);
    ws->var("enne")->setConstant(kTRUE);
    ws->var("enneW")->setConstant(kTRUE);
    ws->var("coeffGaus")->setConstant(kTRUE);
    ws->var("sigmaSig1")->setConstant(kTRUE);
    ws->var("sigmaSig2")->setConstant(kTRUE);
    ws->var("meanSig1")->setConstant(kTRUE);
    ws->var("coefExp")->setConstant(kTRUE);   //mh
    ws->var("NSig")->setConstant(kTRUE);
    ws->var("NBkg")->setConstant(kTRUE);

    RooFormulaVar fBkg("fBkg","@0/(@0+@1)",RooArgList(*(ws->var("NBkg")),*(ws->var("NSig"))));    ws->import(fBkg);
    sprintf(funct,"PROD::totSIGPR(%s,sigPR)",mSigFunct.c_str()); ws->factory(funct);
    sprintf(funct,"PROD::totSIGNP(%s,sigNP)",mSigFunct.c_str()); ws->factory(funct);
    RooProdPdf totSIGPR_PEE("totSIGPR_PEE","PDF with PEE", *(ws->pdf("errPdfSig")),
                           Conditional( *(ws->pdf("totSIGPR")), RooArgList(*(ws->var("Jpsi_Ct")),*(ws->var("Jpsi_Mass"))) )
                           );  ws->import(totSIGPR_PEE);
    RooProdPdf totSIGNP_PEE("totSIGNP_PEE","PDF with PEE", *(ws->pdf("errPdfSig")),
                           Conditional( *(ws->pdf("totSIGNP")), RooArgList(*(ws->var("Jpsi_Ct")),*(ws->var("Jpsi_Mass"))))
                           );  ws->import(totSIGNP_PEE);    
    RooProdPdf totBKG_PEE("totBKG_PEE","PDF with PEE", *(ws->pdf("errPdfBkg")),
                         Conditional( *(ws->pdf("totBKG")), RooArgList(*(ws->var("Jpsi_Ct")),*(ws->var("Jpsi_Mass"))) )
                         );  ws->import(totBKG_PEE);
    
    ws->factory("RSUM::totPDF_PEE(fBkg*totBKG_PEE,Bfrac[0.25,0.,1.]*totSIGNP_PEE,totSIGPR_PEE)");

    // ** Check scaled ct err dist curves and data points in signal region are matched
    RooPlot *errframe3 = ws->var("Jpsi_CtErr")->frame();
    subtrData->plotOn(errframe3,DataError(RooAbsData::SumW2));
    ws->pdf("errPdfSig")->plotOn(errframe3,LineColor(kBlue),Normalization(subtrData->sumEntries(),RooAbsReal::NumEvent));
    binDataCtErrSIG->plotOn(errframe3,DataError(RooAbsData::SumW2),LineColor(kRed),MarkerColor(kRed));  //Not subtracted D_sig
    
    TCanvas ctest3;
    ctest3.cd(); errframe3->Draw();
    titlestr = dirPre + "_rap" + yrange + "_pT" + prange + "_cent" + crange + "_dPhi" + phirange + "_testErrPdfSig_Lin.pdf";
    ctest3.SaveAs(titlestr.c_str());
    ctest3.Clear(); ctest3.SetLogy(1); errframe3->Draw();
    titlestr = dirPre + "_rap" + yrange + "_pT" + prange + "_cent" + crange + "_dPhi" + phirange + "_testErrPdfSig_Log.pdf";
    ctest3.SaveAs(titlestr.c_str());

    RooPlot *errframe = ws->var("Jpsi_CtErr")->frame();
    binDataCtErrSIG->plotOn(errframe,DataError(RooAbsData::SumW2),MarkerColor(kRed),LineColor(kRed));
    binDataCtErrSB->plotOn(errframe,DataError(RooAbsData::SumW2),MarkerColor(kGreen+2),LineColor(kGreen+2),MarkerStyle(24));
    scaledBkg->plotOn(errframe,DataError(RooAbsData::SumW2),MarkerColor(kBlue),MarkerStyle(24),LineColor(kBlue));
    subtrData->plotOn(errframe,DataError(RooAbsData::SumW2));

    ctest3.Clear(); ctest3.SetLogy(1); errframe->Draw();
    titlestr = dirPre + "_rap" + yrange + "_pT" + prange + "_cent" + crange + "_dPhi" + phirange + "_testErrPdfSigDetail_Log.pdf";
    ctest3.SaveAs(titlestr.c_str());
   
  } else {
    cout << "Should use prefit mass option. exit.\n";
    return -1;

    ws->factory("PROD::totSigPR(sigCBG1,sigPR)");
    ws->factory("PROD::totSigNP(sigCBG1,sigNP)");
    ws->factory("SUM::totPDF(NSigPR[4000.0,10.0,1000000.0]*totSigPR,NSigNP[900.0,10.,1000000.]*totSigNP,NBkg[1400.,10.,1000000.]*totBKG)");   //Final F(l,m)
  }

  // *** Start prefit on the signal ctau function
  RooFitResult *fitPR, *fitSB;
  double RSS = 0;
  unsigned int nFullBinsResid = 0;
  if (prefitSignalCTau) {
    RooProdPdf sigPR_PEE("sigPR_PEE","PDF with PEE", *(ws->pdf("errPdfSig")),
                        Conditional(*(ws->pdf("sigPR")), RooArgList(*(ws->var("Jpsi_Ct"))))
                        );  ws->import(sigPR_PEE);

    fitPR = ws->pdf("sigPR_PEE")->fitTo(*redMCCutPR,Range("promptfit"),SumW2Error(kTRUE),ConditionalObservables(RooArgSet(*(ws->var("Jpsi_CtErr")))),Save(1),NumCPU(6));
//    resultF.cd();    cout <<"fitPR->Write(): " << fitPR->Write("sigPR") << endl;
    if (ws->var("sigmaResSigW")) ws->var("sigmaResSigW")->setConstant(kTRUE);
    ws->var("meanResSigW")->setConstant(kTRUE);

    // Check prompt fit is fine with per event error fit. CtWeighted means l/err l
    RooRealVar* CtWeighted = new RooRealVar("CtWeighted","#font[12]{l}_{J/#psi} / #sigma( #font[12]{l}_{J/#psi} )",-5.,5.);
    ws->import(*CtWeighted);
    const RooArgSet* thisRow = (RooArgSet*)redMCCutPR->get(0); 
    RooArgSet* newRow = new RooArgSet(*CtWeighted);
    RooDataSet* tempJpsi = new RooDataSet("tempJpsi","new data",*newRow);
    for (Int_t iSamp = 0; iSamp < redMCCutPR->numEntries(); iSamp++) {
      thisRow = (RooArgSet*)redMCCutPR->get(iSamp);
      RooRealVar* myct = (RooRealVar*)thisRow->find("Jpsi_Ct");
      RooRealVar* mycterr = (RooRealVar*)thisRow->find("Jpsi_CtErr");
      CtWeighted->setVal(myct->getVal()/mycterr->getVal());
      RooArgSet* tempRow = new RooArgSet(*CtWeighted);
      tempJpsi->add(*tempRow);
    }

    if (oneGaussianResol) {
      ws->factory("Gaussian::tempsigPR(CtWeighted,meanResSigW,sigmaResSigN)");
    } else {
      ws->factory("Gaussian::tempresGW(CtWeighted,meanResSigW,sigmaResSigW)");
      ws->factory("Gaussian::tempresGN(CtWeighted,meanResSigW,sigmaResSigN)");
      ws->factory("SUM::tempsigPR(fracRes*tempresGW,tempresGN)");
    }  

    RooPlot *tframePR = ws->var("CtWeighted")->frame();
    tempJpsi->plotOn(tframePR,DataError(RooAbsData::SumW2));
    ws->pdf("tempsigPR")->plotOn(tframePR,LineColor(kBlue),Normalization(tempJpsi->sumEntries(),RooAbsReal::NumEvent));
    TCanvas c00; c00.cd();tframePR->Draw();
    sprintf(reduceDS,"Entry: %0.0f",tempJpsi->sumEntries());
    t->DrawLatex(0.5,0.5,reduceDS);
    titlestr = dirPre + "_rap" + yrange + "_pT" + prange + "_cent" + crange + "_dPhi" + phirange + "_resolInit_Lin.pdf";
    c00.SaveAs(titlestr.c_str());
    c00.Clear(); c00.SetLogy(1); tframePR->Draw();
    sprintf(reduceDS,"Entry: %0.0f",tempJpsi->sumEntries());
    t->DrawLatex(0.5,0.5,reduceDS);
    titlestr = dirPre + "_rap" + yrange + "_pT" + prange + "_cent" + crange + "_dPhi" + phirange + "_resolInit_Log.pdf";
    c00.SaveAs(titlestr.c_str());

    if (ws->var("sigmaResSigW")) ws->var("sigmaResSigW")->setConstant(kFALSE);
    ws->var("meanResSigW")->setConstant(kFALSE);

    // Get sum of squared residual to fisher's F-test
/*    RooHist *hresid;
    hresid = tframePR->residHist();
    hresid->SetName("residual");
    double *yresid = hresid->GetY();
    unsigned int nBins = ws->var("Jpsi_Ct")->getBinning().numBins();
    for (unsigned int i=0; i < nBins; i++) {
      RSS = RSS + pow(yresid[i],2);
      if (fabs(yresid[i]) > 0.0001) nFullBinsResid++;
    }

    RooPlot* tframeres =  ws->var("Jpsi_Ct")->frame(Title("Residual Distribution")) ;
    tframeres->SetLabelSize(0.08,"XYZ");
    tframeres->SetTitleSize(0.1,"XYZ");
    tframeres->SetTitleOffset(0.55,"Y");
    tframeres->addPlotable(hresid,"P") ;
    tframeres->SetMaximum(-(tframeres->GetMinimum())); 
    tframeres->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi}/#font[12]{#sigma({l}_{J/#psi})} (mm)");
    tframeres->GetXaxis()->CenterTitle(1);

    c00.cd();
    TPad padpr1("padpr1","This is pad1",0.05,0.35,0.95,0.97);
    padpr1.SetBottomMargin(0);
    padpr1.Draw();
    TPad padpr2("padpr2","This is pad2",0.05,0.03,0.95,0.35);
    padpr2.SetTopMargin(0);
    padpr2.SetBottomMargin(0.24);
    padpr2.Draw();

    padpr1.cd(); tframePR->Draw();
    padpr2.cd();
    tframeres->SetMinimum(0.5); 
    tframeres->Draw();

    TLatex tpr;
    tpr.SetNDC(); tpr.SetTextAlign(22);
    tpr.SetTextSize(0.07);
    sprintf(reduceDS,"RSS = %.0f, nFullBinsResid = %d",RSS,nFullBinsResid);
    tpr.DrawLatex(0.7,0.90,reduceDS);

    titlestr = dirPre + "_rap" + yrange + "_pT" + prange + "_cent" + crange + "_dPhi" + phirange + "_2D_" + partFile + "resolfit_Lin.pdf";
    c00.SaveAs(titlestr.c_str());
    
    padpr1.cd(); padpr1.SetLogy(1); tframePR->Draw();
    padpr2.cd(); padpr2.SetLogy(1); tframeres->Draw();
    tpr.DrawLatex(0.7,0.90,reduceDS);
    titlestr = dirPre + "_rap" + yrange + "_pT" + prange + "_cent" + crange + "_dPhi" + phirange + "_2D_" + partFile + "resolfit_Log.pdf";
    c00.SaveAs(titlestr.c_str());
*/
  }

  double bfraction[2] = {0};
  if (prefitBkg) {
    cout << "DATA :: N events to fit on the sidebands: " << binDataSB->sumEntries() << endl;
//    ws->var("fpm")->setConstant(kTRUE); //mh

    if (prefitSignalCTau) {
      if (ws->var("fracRes")) ws->var("fracRes")->setConstant(kTRUE);
      ws->var("meanResSigW")->setConstant(kTRUE);
      if (ws->var("sigmaResBkgW")) ws->var("sigmaResBkgW")->setVal(ws->var("sigmaResSigW")->getVal());
      if (ws->var("sigmaResBkgN")) ws->var("sigmaResBkgN")->setVal(ws->var("sigmaResSigN")->getVal());
    }


    fitSB = ws->pdf("bkgCtauTOT_PEE")->fitTo(*redDataSB,SumW2Error(kTRUE),Minos(0),NumCPU(6),Save(1),ConditionalObservables(RooArgSet(*(ws->var("Jpsi_CtErr")))));
//    resultF.cd(); cout << fitSB->Write("bkgCtTot") << endl;
   
    ws->var("fpm")->setConstant(kTRUE);
    ws->var("fLiving")->setConstant(kTRUE);
    ws->var("fbkgCtTot")->setConstant(kTRUE);
    if (ws->var("sigmaResBkgN")) ws->var("sigmaResBkgN")->setConstant(kTRUE);
    if (ws->var("sigmaResBkgW")) ws->var("sigmaResBkgW")->setConstant(kTRUE);
    if (ws->var("meanResBkgW")) ws->var("meanResBkgW")->setConstant(kTRUE);
    ws->var("lambdap")->setConstant(kTRUE);
    ws->var("lambdam")->setConstant(kTRUE);
    ws->var("lambdasym")->setConstant(kTRUE);
    
    if (prefitSignalCTau) {
      if (ws->var("fracRes")) ws->var("fracRes")->setConstant(kFALSE);
      ws->var("meanResSigW")->setConstant(kFALSE);
    }

    RooPlot *tframe1 = ws->var("Jpsi_Ct")->frame();
    titlestr = "2D fit for" + partTit + "muons (J/ #psi c  #tau projection, mass sidebands), p_{T} = " + prange + " GeV/c and |y| = " + yrange;
    tframe1->SetTitle(titlestr.c_str());
    redDataSB->plotOn(tframe1,DataError(RooAbsData::SumW2),Binning(rb2));
    ws->pdf("bkgCtauTOT_PEE")->plotOn(tframe1,ProjWData(RooArgList(*(ws->var("Jpsi_CtErr"))),*binDataCtErr,kTRUE),NumCPU(6),Normalization(redDataSB->sumEntries(),RooAbsReal::NumEvent));

    RooHist *hpull = tframe1->pullHist(); hpull->SetName("hpull");
    int nFitPar = fitSB->floatParsFinal().getSize();

    TCanvas *c3 = new TCanvas("c3","The Canvas",200,10,600,880);
    c3->cd();
    TPad *pad1 = new TPad("pad1","This is pad1",0.05,0.35,0.95,0.97);
    pad1->SetBottomMargin(0);  pad1->Draw();
    TPad *pad2 = new TPad("pad2","This is pad2",0.05,0.03,0.95,0.35);
    pad2->SetTopMargin(0);  pad2->SetBottomMargin(0.24);  pad2->Draw();

    pad1->cd(); tframe1->Draw();
    t->SetTextSize(0.05);
    t->DrawLatex(0.44,0.90,"CMS Preliminary");
    t->DrawLatex(0.44,0.83,"PbPb #sqrt{s_{NN}} = 2.76 TeV");
    t->SetTextSize(0.04);
    t->DrawLatex(0.47,0.76,"L_{int} =  150 #mub^{-1}");
    t->SetTextSize(0.030);
    //    sprintf(reduceDS,"%.0f-%.0f%, %0.1f < |y| < %.1f",cmin,cmax,ymin,ymax);
    t->DrawLatex(0.47,0.70,reduceDS);
    //    sprintf(reduceDS,"%.1f < p_{T} < %.1f GeV/c",pmin,pmax);
    t->DrawLatex(0.47,0.66,reduceDS);
    //    sprintf(reduceDS,"%.2f < |#phi_{J/#psi}-#Psi_{RP}| < %.2f",psmin,psmax);
    t->DrawLatex(0.47,0.62,reduceDS);
    //    sprintf(reduceDS,"EP: %s",rpmethod.c_str());
    t->DrawLatex(0.47,0.58,reduceDS);

    double chi2 = 0, unNormChi2 = 0;
    int dof = 0;
    double *ypulls = hpull->GetY();
    unsigned int nBins = ws->var("Jpsi_Ct")->getBinning().numBins();
    unsigned int nFullBins = 0;
    for (unsigned int i = 0; i < nBins; i++) {
      chi2 += ypulls[i]*ypulls[i];
      if (fabs(ypulls[i]) > 0.0001) nFullBins++;
    }
    unNormChi2 = chi2;
    dof = nFullBins - nFitPar;
    chi2 /= (nFullBins - nFitPar);
    for (unsigned int i = 0; i < nBins; i++) {
      if (fabs(ypulls[i]) < 0.0001) ypulls[i] = 999.; 
    } 
    int nDOF = ws->var("Jpsi_Ct")->getBinning().numBins() - nFitPar;

    RooPlot* tframepull =  ws->var("Jpsi_Ct")->frame(Title("Pull")) ;
    tframepull->GetYaxis()->SetTitle("Pull");
    tframepull->SetLabelSize(0.08,"XYZ");
    tframepull->SetTitleSize(0.1,"XYZ");
    tframepull->SetTitleOffset(0.55,"Y");
    tframepull->addPlotable(hpull,"P") ;
    tframepull->SetMaximum(-(tframepull->GetMinimum())); 
    tframepull->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
    tframepull->GetXaxis()->CenterTitle(1);

    pad2->cd(); tframepull->Draw();

    TLatex *t2 = new TLatex();
    t2->SetNDC(); t2->SetTextAlign(22); t2->SetTextSize(0.07);
    sprintf(reduceDS,"#chi^{2}/dof = %.2f/%d",unNormChi2,dof);
    t2->DrawLatex(0.76,0.90,reduceDS);

    titlestr = dirPre + "_rap" + yrange + "_pT" + prange + "_cent" + crange + "_dPhi" + phirange + "_timeside_Lin.pdf";
    c3->SaveAs(titlestr.c_str());

    TCanvas* c3a = new TCanvas("c3a","The Canvas",200,10,600,880);
    c3a->cd();
    TPad *pad1a = new TPad("pad1a","This is pad1",0.05,0.35,0.95,0.97);
    pad1a->SetBottomMargin(0);
    pad1a->Draw();
    TPad *pad2a = new TPad("pad2a","This is pad2",0.05,0.03,0.95,0.35);
    pad2a->SetTopMargin(0);
    pad2a->SetBottomMargin(0.24);
    pad2a->Draw();

    pad1a->cd(); pad1a->SetLogy(1);
    tframe1->SetMaximum(tframe1->GetMaximum()*9); 
    tframe1->SetMinimum(0.5); 
    tframe1->Draw();
   
    t->SetTextSize(0.05);
    t->DrawLatex(0.17,0.90,"CMS Preliminary");
    t->DrawLatex(0.17,0.85,"PbPb  #sqrt{s_{NN}} = 2.76 TeV");
    t->SetTextSize(0.04);
    t->DrawLatex(0.17,0.79,"L_{int} =  150 #mub^{-1}");
    t->SetTextSize(0.035);
    //    sprintf(reduceDS,"%.0f-%.0f%, %0.1f < |y| < %.1f",cmin,cmax,ymin,ymax);
    t->DrawLatex(0.52,0.80,reduceDS);
    //    sprintf(reduceDS,"%.1f < p_{T} < %.1f GeV/c",pmin,pmax);
    t->DrawLatex(0.52,0.75,reduceDS);
    //    sprintf(reduceDS,"%.2f < |#phi_{J/#psi}-#Psi_{RP}| < %.2f",psmin,psmax);
    t->DrawLatex(0.52,0.70,reduceDS);
    //    sprintf(reduceDS,"EP: %s",rpmethod.c_str());
    t->DrawLatex(0.52,0.65,reduceDS);

    pad2a->cd(); tframepull->Draw();
    sprintf(reduceDS,"#chi^{2}/dof = %.2f/%d",unNormChi2,dof);
    t2->DrawLatex(0.76,0.90,reduceDS);

    titlestr = dirPre + "_rap" + yrange + "_pT" + prange + "_cent" + crange + "_dPhi" + phirange + "_timeside_Log.pdf";
    c3a->SaveAs(titlestr.c_str());

    delete pad1;
    delete pad2;
    delete c3;
    delete pad1a;
    delete pad2a;
    delete c3a;
  }

  // Fix below bkg variables in any case
  ws->var("fpm")->setConstant(kTRUE);
  ws->var("fLiving")->setConstant(kTRUE);

  // *** Get NSig, NBkg, Bfraction and their errors
  Double_t NSigPR_fin, ErrNSigPR_fin;
  Double_t NSigNP_fin, ErrNSigNP_fin;
  Double_t Bfrac_fin, ErrBfrac_fin;
  int nFitPar;
  Double_t theNLL;

  if (prefitMass) {
    RooFitResult *rfr;
    if (redDataCut->sumEntries() < 5000) {
      rfr = ws->pdf("totPDF_PEE")->fitTo(*redDataCut,Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(6),ConditionalObservables(RooArgSet(*(ws->var("Jpsi_CtErr")))));
    } else {
      rfr = ws->pdf("totPDF_PEE")->fitTo(*binData,Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(6),ConditionalObservables(RooArgSet(*(ws->var("Jpsi_CtErr")))));
    }
//    resultF.cd(); cout << rfr->Write("totPDF") << endl;
    nFitPar = rfr->floatParsFinal().getSize();
    theNLL = rfr->minNll();
    Bfrac_fin = ws->var("Bfrac")->getVal();
    ErrBfrac_fin = ws->var("Bfrac")->getError();
    if (ws->var("fracRes")) ws->var("fracRes")->setConstant(kTRUE);
    ws->var("meanResSigW")->setConstant(kTRUE);
    NSigNP_fin = NSig_fin * Bfrac_fin;
    NSigPR_fin = NSig_fin * (1-Bfrac_fin);
    ErrNSigNP_fin = NSigNP_fin * sqrt( pow(ErrNSig_fin/NSig_fin,2)+pow(ErrBfrac_fin/Bfrac_fin,2) );
    ErrNSigPR_fin = NSigPR_fin * sqrt ( pow(ErrNSig_fin/NSig_fin,2)+pow(ErrBfrac_fin/(1.0-Bfrac_fin),2) );
  } else {
    RooFitResult *rfr = ws->pdf("totPDF")->fitTo(*redDataCut,Extended(1),Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(6));
    nFitPar = rfr->floatParsFinal().getSize();
    NSigNP_fin = ws->var("NSigNP")->getVal();
    NSigPR_fin = ws->var("NSigPR")->getVal();
    ErrNSigNP_fin = ws->var("NSigNP")->getError();
    ErrNSigPR_fin = ws->var("NSigPR")->getError();
    Bfrac_fin = NSigNP_fin/(NSigNP_fin+NSigPR_fin);
    ErrBfrac_fin = sqrt( pow(NSigNP_fin*ErrNSigPR_fin,2) + pow(NSigPR_fin*ErrNSigNP_fin,2) ) / pow(NSigNP_fin+NSigPR_fin,2);
  }

  const double coeffGaus = ws->var("fracRes")->getVal();
  const double sigmaSig1 = ws->var("sigmaResSigW")->getVal();
  const double sigmaSig2 = ws->var("sigmaResSigN")->getVal();
  const double ErrcoeffGaus = 0;
  const double ErrsigmaSig1 = ws->var("sigmaResSigW")->getError();
  const double ErrsigmaSig2 = ws->var("sigmaResSigN")->getError();

  double resol = sqrt( coeffGaus*pow(sigmaSig1,2) + (1-coeffGaus)*pow(sigmaSig2,2) );
  double Errresol = (0.5/resol) *
                sqrt( pow(sigmaSig1*coeffGaus*ErrsigmaSig1,2) +
                      pow(sigmaSig2*(1-coeffGaus)*ErrsigmaSig2,2) +
                      pow(0.5*(pow(sigmaSig1,2)-pow(sigmaSig2,2))*ErrcoeffGaus,2) );

  // *** Plot various fit results and data points
  // Temporary variables for plotting
  RooRealVar tmpVar1("tmpVar1","tmpVar1",NSigNP_fin);
  RooRealVar tmpVar2("tmpVar2","tmpVar2",NBkg_fin);

  // Mass plot
  RooPlot *mframe = ws->var("Jpsi_Mass")->frame();
  redDataCut->plotOn(mframe,DataError(RooAbsData::SumW2),XErrorSize(0),MarkerSize(1),Binning(rbm));
  titlestr = "2D fit for" + partTit + "muons (mass projection), p_{T} = " + prange  + "_dPhi" + phirange + " GeV/c and |y| = " + yrange;
  mframe->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c^{2})");
  mframe->GetXaxis()->CenterTitle(1);
  const double max = mframe->GetMaximum() * 1.3;
  mframe->SetMaximum(max);

  if (prefitMass) {
    // Fill color
    ws->pdf("totPDF_PEE")->plotOn(mframe,DrawOption("F"),FillColor(kBlack),FillStyle(3354),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
    RooAddPdf tmpPDF("tmpPDF","tmpPDF",RooArgList(*(ws->pdf(mSigFunct.c_str())),*(ws->pdf(mBkgFunct.c_str()))),RooArgList(tmpVar1,tmpVar2));
    tmpPDF.plotOn(mframe,LineColor(kRed),DrawOption("F"),FillColor(kWhite),FillStyle(1001),Normalization(NSigNP_fin+NBkg_fin,RooAbsReal::NumEvent));
    tmpPDF.plotOn(mframe,LineColor(kRed),DrawOption("F"),FillColor(kRed),FillStyle(3444),Normalization(NSigNP_fin+NBkg_fin,RooAbsReal::NumEvent));
    gStyle->SetHatchesLineWidth(2);
    ws->pdf("totPDF_PEE")->plotOn(mframe,Components(mBkgFunct.c_str()),DrawOption("F"),FillColor(kAzure-9),FillStyle(1001),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
    //Line color
    ws->pdf("totPDF_PEE")->plotOn(mframe,Components(mBkgFunct.c_str()),LineColor(kBlue),LineStyle(7),LineWidth(5),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
    tmpPDF.plotOn(mframe,LineColor(kRed),LineStyle(9),LineWidth(5),Normalization(NSigNP_fin+NBkg_fin,RooAbsReal::NumEvent));
    ws->pdf("totPDF_PEE")->plotOn(mframe,LineColor(kBlack),LineWidth(2),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
  } else {
    ws->pdf("totPDF")->plotOn(mframe,Components("totSigNP,totBKG"),LineColor(kRed),Normalization(1.0,RooAbsReal::RelativeExpected));
    ws->pdf("totPDF")->plotOn(mframe,Components("totBKG"),LineColor(kBlue),LineColor(7),Normalization(1.0,RooAbsReal::RelativeExpected));
    ws->pdf("totPDF")->plotOn(mframe,LineColor(kBlack),Normalization(1.0,RooAbsReal::RelativeExpected));
  }
  redDataCut->plotOn(mframe,DataError(RooAbsData::SumW2),XErrorSize(0),MarkerSize(1),Binning(rbm));

  TCanvas c1; c1.cd(); mframe->Draw();
  t->SetTextSize(0.05);
  t->DrawLatex(0.17,0.90,"CMS Preliminary");
  t->DrawLatex(0.17,0.82,"PbPb  #sqrt{s_{NN}} = 2.76 TeV");
  t->SetTextSize(0.04);
  t->DrawLatex(0.17,0.75,"L_{int} =  150 #mub^{-1}");
  t->SetTextSize(0.03);
  //  sprintf(reduceDS,"%.0f-%.0f%, %0.1f < |y| < %.1f",cmin,cmax,ymin,ymax);
  t->DrawLatex(0.17,0.69,reduceDS);
  //  sprintf(reduceDS,"%.1f < p_{T} < %.1f GeV/c",pmin,pmax);
  t->DrawLatex(0.17,0.64,reduceDS);
  //  sprintf(reduceDS,"%.2f < |#phi_{J/#psi}-#Psi_{RP}| < %.2f",psmin,psmax);
  t->DrawLatex(0.17,0.59,reduceDS);
  //  sprintf(reduceDS,"EP: %s",rpmethod.c_str());
  t->DrawLatex(0.17,0.54,reduceDS);
  t->SetTextSize(0.04);
  //  sprintf(reduceDS,"yield: %0.0f",ws->var("NSig")->getVal());
  t->DrawLatex(0.65,0.80,reduceDS);
  //  sprintf(reduceDS,"#sigma = %0.0f MeV/c^{2}",ws->var("sigmaSig2")->getVal()*1000);
  t->DrawLatex(0.65,0.75,reduceDS);


  TLegend * leg11 = new TLegend(0.63,0.55,0.92,0.73,NULL,"brNDC");
  leg11->SetFillStyle(0); leg11->SetBorderSize(0); leg11->SetShadowColor(0);
  leg11->SetMargin(0.2);
  leg11->AddEntry(gfake1,"data","p");
  leg11->AddEntry(&hfake21,"total fit","lf");
  leg11->AddEntry(&hfake31,"bkgd + non-prompt","lf"); 
  leg11->AddEntry(&hfake11,"background","lf");
  leg11->Draw("same");
  titlestr = dirPre + "_rap" + yrange + "_pT" + prange + "_cent" + crange + "_dPhi" + phirange + "_massfit.pdf";
  c1.SaveAs(titlestr.c_str());

  // Ctau plots with/without pull distribution
  // Check prompt fit is fine with per event error fit. CtWeighted means l/err l
  RooRealVar *CtWeighted = (RooRealVar*)ws->var("CtWeighted");
  const RooArgSet *thisRow = (RooArgSet*)redDataCut->get(0); 
  RooArgSet *newRow = new RooArgSet(*CtWeighted);
  RooDataSet *tempJpsiD = new RooDataSet("tempJpsiD","new data",*newRow);
  for (Int_t iSamp = 0; iSamp < redDataCut->numEntries(); iSamp++) {
    thisRow = (RooArgSet*)redDataCut->get(iSamp);
    RooRealVar *myct = (RooRealVar*)thisRow->find("Jpsi_Ct");
    RooRealVar *mycterr = (RooRealVar*)thisRow->find("Jpsi_CtErr");
    CtWeighted->setVal(myct->getVal()/mycterr->getVal());
    RooArgSet* tempRowD = new RooArgSet(*CtWeighted);
    tempJpsiD->add(*tempRowD);
  }

  RooPlot *tframePR = ws->var("CtWeighted")->frame();
  tempJpsiD->plotOn(tframePR,DataError(RooAbsData::SumW2));
  ws->pdf("tempsigPR")->plotOn(tframePR,LineColor(kBlue),Normalization(tempJpsiD->sumEntries(),RooAbsReal::NumEvent));
  TCanvas c00; c00.cd();tframePR->Draw();
  sprintf(reduceDS,"Entry: %.0f",tempJpsiD->sumEntries());
  t->DrawLatex(0.5,0.5,reduceDS);
  titlestr = dirPre + "_rap" + yrange + "_pT" + prange + "_cent" + crange + "_dPhi" + phirange + "_resolfit_Lin.pdf";
  c00.SaveAs(titlestr.c_str());
  c00.Clear(); c00.SetLogy(1); tframePR->Draw();
  sprintf(reduceDS,"Entry: %.0f",tempJpsiD->sumEntries());
  t->DrawLatex(0.5,0.5,reduceDS);
  titlestr = dirPre + "_rap" + yrange + "_pT" + prange + "_cent" + crange + "_dPhi" + phirange + "_resolfit_Log.pdf";
  c00.SaveAs(titlestr.c_str());

  ws->var("Jpsi_Ct")->setBinning(rb3);
  RooPlot *tframe = ws->var("Jpsi_Ct")->frame();
  titlestr = "2D fit for" + partTit + "muons (c#tau projection), p_{T} = " + prange  + "_dPhi" + phirange + " GeV/c and |y| = " + yrange;
  tframe->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
  tframe->GetXaxis()->CenterTitle(1);
  tframe->GetYaxis()->SetTitle("Events / (0.088 mm)");

  // Ctau total distributions
  RooHist *hpull;
//  redDataCut->plotOn(tframe,DataError(RooAbsData::SumW2),Binning(rb2),MarkerSize(1));
  redDataCut->plotOn(tframe,DataError(RooAbsData::SumW2),Binning(rb3),MarkerSize(1));

  if (prefitMass) {
    ws->pdf("totPDF_PEE")->plotOn(tframe,LineColor(kBlack),LineWidth(2),ProjWData(RooArgList(*(ws->var("Jpsi_CtErr"))),*binDataCtErr,kTRUE),NumCPU(6),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
    hpull = tframe->pullHist(); hpull->SetName("hpull");
    ws->pdf("totPDF_PEE")->plotOn(tframe,Components("totBKG"),LineColor(kBlue),LineWidth(5),ProjWData(RooArgList(*(ws->var("Jpsi_CtErr"))),*binDataCtErr,kTRUE),NumCPU(6),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent),LineStyle(7));
    if (superImpose) {
      RooAddPdf tmpPDF2("tmpPDF2","tmpPDF2",RooArgList(*(ws->pdf("totSIGNP")),*(ws->pdf("totBKG"))),RooArgList(tmpVar1,tmpVar2));
      tmpPDF2.plotOn(tframe,LineColor(kRed),ProjWData(RooArgList(*(ws->var("Jpsi_CtErr"))),*binDataCtErr,kTRUE),NumCPU(6),Normalization(NSigNP_fin+NBkg_fin,RooAbsReal::NumEvent),LineWidth(5),LineStyle(9));
    } else {
      ws->pdf("totPDF_PEE")->plotOn(tframe,Components("totSIGNP"),LineColor(kRed),ProjWData(RooArgList(*(ws->var("Jpsi_CtErr"))),*binDataCtErr,kTRUE),NumCPU(6),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent),LineStyle(kDashed));
      ws->pdf("totPDF_PEE")->plotOn(tframe,Components("totSIGPR"),LineColor(kGreen),ProjWData(RooArgList(*(ws->var("Jpsi_CtErr"))),*binDataCtErr,kTRUE),NumCPU(6),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent),LineStyle(kDashDotted));
    }
    ws->pdf("totPDF_PEE")->plotOn(tframe,LineColor(kBlack),LineWidth(2),ProjWData(RooArgList(*(ws->var("Jpsi_CtErr"))),*binDataCtErr,kTRUE),NumCPU(6),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
  } else {
    ws->pdf("totPDF")->plotOn(tframe,LineColor(kBlack),Normalization(1.0,RooAbsReal::RelativeExpected));
    hpull = tframe->pullHist(); hpull->SetName("hpull");
    ws->pdf("totPDF")->plotOn(tframe,Components("totSigNP,totBKG"),LineColor(kRed),Normalization(1.0,RooAbsReal::RelativeExpected),LineStyle(9));
    ws->pdf("totPDF")->plotOn(tframe,Components("totBKG"),LineColor(kBlue),Normalization(1.0,RooAbsReal::RelativeExpected),LineStyle(7));
    ws->pdf("totPDF")->plotOn(tframe,LineColor(kBlack),Normalization(1.0,RooAbsReal::RelativeExpected));
  }

  double chi2 = 0, unNormChi2 = 0;
  int dof = 0;
  double *ypulls = hpull->GetY();
  unsigned int nBins = ws->var("Jpsi_Ct")->getBinning().numBins();
  unsigned int nFullBins = 0;
  for (unsigned int i = 0; i < nBins; i++) {
    chi2 += ypulls[i]*ypulls[i];
    if (fabs(ypulls[i]) > 0.0001) nFullBins++;
  }
  unNormChi2 = chi2;
  dof = nFullBins - nFitPar;
  chi2 /= (nFullBins - nFitPar);
  for (unsigned int i = 0; i < nBins; i++) {
    if (fabs(ypulls[i]) < 0.0001) ypulls[i] = 999.; 
  } 

  // WITH RESIDUALS
  TCanvas* c2 = new TCanvas("c2","The Canvas",200,10,600,880);
  c2->cd();
  TPad *pad1 = new TPad("pad1","This is pad1",0.05,0.35,0.95,0.97);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  TPad *pad2 = new TPad("pad2","This is pad2",0.05,0.03,0.95,0.35);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.24);
  pad2->Draw();

  pad1->cd(); tframe->Draw();

  t->SetTextSize(0.05);
  t->DrawLatex(0.44,0.90,"CMS Preliminary");
  t->DrawLatex(0.44,0.83,"PbPb #sqrt{s_{NN}} = 2.76 TeV");
  t->SetTextSize(0.04);
  t->DrawLatex(0.47,0.76,"L_{int} =  150 #mub^{-1}");
  t->SetTextSize(0.030);
  //  sprintf(reduceDS,"%.0f-%.0f%, %0.1f < |y| < %.1f",cmin,cmax,ymin,ymax);
  t->DrawLatex(0.47,0.70,reduceDS);
  sprintf(reduceDS,"%.1f < p_{T} < %.1f GeV/c",pmin,pmax);
  t->DrawLatex(0.47,0.66,reduceDS);
  sprintf(reduceDS,"%.2f < |#phi_{J/#psi}-#Psi_{RP}| < %.2f",psmin,psmax);
  t->DrawLatex(0.47,0.62,reduceDS);
  sprintf(reduceDS,"EP: %s",rpmethod.c_str());
  t->DrawLatex(0.47,0.58,reduceDS);

  TLegend * leg = new TLegend(0.455,0.35,0.85,0.54,NULL,"brNDC");
  leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetShadowColor(0);
  leg->SetMargin(0.2);
  leg->AddEntry(gfake1,"data","p");
  leg->AddEntry(&hfake21,"total fit","lf");
  leg->AddEntry(&hfake31,"bkgd + non-prompt","lf"); 
  leg->AddEntry(&hfake11,"background","lf");
  leg->Draw("same"); 

  RooPlot* tframepull =  ws->var("Jpsi_Ct")->frame(Title("Pull")) ;
  tframepull->GetYaxis()->SetTitle("Pull");
  tframepull->SetLabelSize(0.08,"XYZ");
  tframepull->SetTitleSize(0.1,"XYZ");
  tframepull->SetTitleOffset(0.55,"Y");
//  tframepull->SetTitleOffset(0.95,"X");
  tframepull->addPlotable(hpull,"P") ;
  tframepull->SetMaximum(-(tframepull->GetMinimum())); 
  tframepull->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
  tframepull->GetXaxis()->CenterTitle(1);

  pad2->cd(); tframepull->Draw();

  int nDOF = ws->var("Jpsi_Ct")->getBinning().numBins() - nFitPar;

  TLatex *t2 = new TLatex();
  t2->SetNDC(); t2->SetTextAlign(22); t2->SetTextSize(0.07);
  sprintf(reduceDS,"#chi^{2}/dof = %.2f/%d",unNormChi2,dof);
  t2->DrawLatex(0.76,0.90,reduceDS);
  
  c2->Update();
  titlestr = dirPre + "_rap" + yrange + "_pT" + prange + "_cent" + crange + "_dPhi" + phirange + "_timefit_Lin.pdf";
  c2->SaveAs(titlestr.c_str());

  TCanvas* c2a = new TCanvas("c2a","The Canvas",200,10,600,880);
  c2a->cd();
  TPad *pad1a = new TPad("pad1a","This is pad1",0.05,0.35,0.95,0.97);
  pad1a->SetBottomMargin(0);
  pad1a->Draw();
  TPad *pad2a = new TPad("pad2a","This is pad2",0.05,0.03,0.95,0.35);
  pad2a->SetTopMargin(0);
  pad2a->SetBottomMargin(0.24);
  pad2a->Draw();

  pad1a->cd(); pad1a->SetLogy(1);
  tframe->SetMaximum(tframe->GetMaximum()*9); 
  tframe->SetMinimum(0.5); 
  tframe->Draw();
 
  t->SetTextSize(0.05);
  t->DrawLatex(0.17,0.90,"CMS Preliminary");
  t->DrawLatex(0.17,0.85,"PbPb  #sqrt{s_{NN}} = 2.76 TeV");
  t->SetTextSize(0.04);
  t->DrawLatex(0.17,0.79,"L_{int} =  150 #mub^{-1}");
  t->SetTextSize(0.035);
  //  sprintf(reduceDS,"%.0f-%.0f%, %0.1f < |y| < %.1f",cmin,cmax,ymin,ymax);
  t->DrawLatex(0.52,0.80,reduceDS);
  sprintf(reduceDS,"%.1f < p_{T} < %.1f GeV/c",pmin,pmax);
  t->DrawLatex(0.52,0.75,reduceDS);
  sprintf(reduceDS,"%.2f < |#phi_{J/#psi}-#Psi_{RP}| < %.2f",psmin,psmax);
  t->DrawLatex(0.52,0.70,reduceDS);
  sprintf(reduceDS,"EP: %s",rpmethod.c_str());
  t->DrawLatex(0.52,0.65,reduceDS);

  leg->SetX1NDC(0.52);
  leg->SetY1NDC(0.50);
  leg->SetX2NDC(0.92);
  leg->SetY2NDC(0.62);
  leg->Draw("same");

  pad2a->cd(); tframepull->Draw();
  sprintf(reduceDS,"#chi^{2}/dof = %.2f/%d",unNormChi2,dof);
  t2->DrawLatex(0.76,0.90,reduceDS);
  
  c2a->Update();
  titlestr = dirPre + "_rap" + yrange + "_pT" + prange +  "_cent" + crange + "_dPhi" + phirange + "_timefit_Log.pdf";
  c2a->SaveAs(titlestr.c_str());

  TCanvas* c2b = new TCanvas("c2b","The Canvas",200,10,540,546);
  c2b->cd(); c2b->Draw(); c2b->SetLogy(1);

  RooPlot *tframefill = ws->var("Jpsi_Ct")->frame();
  redDataCut->plotOn(tframefill,DataError(RooAbsData::SumW2),Binning(rb3),MarkerSize(1));
  if (prefitMass) {
    ws->pdf("totPDF_PEE")->plotOn(tframefill,DrawOption("F"),FillColor(kBlack),FillStyle(3354),ProjWData(RooArgList(*(ws->var("Jpsi_CtErr"))),*binDataCtErr,kTRUE),NumCPU(6),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
    ws->pdf("totPDF_PEE")->plotOn(tframefill,LineColor(kBlack),LineWidth(2),ProjWData(RooArgList(*(ws->var("Jpsi_CtErr"))),*binDataCtErr,kTRUE),NumCPU(6),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
    RooAddPdf tmpPDF2("tmpPDF2","tmpPDF2",RooArgList(*(ws->pdf("totSIGNP")),*(ws->pdf("totBKG"))),RooArgList(tmpVar1,tmpVar2));
    tmpPDF2.plotOn(tframefill,DrawOption("F"),FillColor(kWhite),FillStyle(1001),ProjWData(RooArgList(*(ws->var("Jpsi_CtErr"))),*binDataCtErr,kTRUE),NumCPU(6),Normalization(NSigNP_fin+NBkg_fin,RooAbsReal::NumEvent));
    tmpPDF2.plotOn(tframefill,DrawOption("F"),FillColor(kRed),FillStyle(3444),ProjWData(RooArgList(*(ws->var("Jpsi_CtErr"))),*binDataCtErr,kTRUE),NumCPU(6),Normalization(NSigNP_fin+NBkg_fin,RooAbsReal::NumEvent));
    gStyle->SetHatchesLineWidth(2);
    ws->pdf("totPDF_PEE")->plotOn(tframefill,Components("totBKG"),DrawOption("F"),FillColor(kAzure-9),FillStyle(1001),ProjWData(RooArgList(*(ws->var("Jpsi_CtErr"))),*binDataCtErr,kTRUE),NumCPU(6),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
    ws->pdf("totPDF_PEE")->plotOn(tframefill,Components("totBKG"),LineColor(kBlue),LineWidth(5),ProjWData(RooArgList(*(ws->var("Jpsi_CtErr"))),*binDataCtErr,kTRUE),NumCPU(6),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent),LineStyle(7));
    tmpPDF2.plotOn(tframefill,LineColor(kRed),ProjWData(RooArgList(*(ws->var("Jpsi_CtErr"))),*binDataCtErr,kTRUE),NumCPU(6),Normalization(NSigNP_fin+NBkg_fin,RooAbsReal::NumEvent),LineWidth(5),LineStyle(9));
  } else {
    ws->pdf("totPDF")->plotOn(tframefill,LineColor(kBlack),Normalization(1.0,RooAbsReal::RelativeExpected));
    ws->pdf("totPDF")->plotOn(tframefill,Components("totSigNP,totBKG"),LineColor(kRed),Normalization(1.0,RooAbsReal::RelativeExpected),LineStyle(9));
    ws->pdf("totPDF")->plotOn(tframefill,Components("totBKG"),LineColor(kBlue),Normalization(1.0,RooAbsReal::RelativeExpected),LineStyle(7));
    ws->pdf("totPDF")->plotOn(tframefill,LineColor(kBlack),Normalization(1.0,RooAbsReal::RelativeExpected));
  }
  redDataCut->plotOn(tframefill,DataError(RooAbsData::SumW2),Binning(rb3),MarkerSize(1));

  tframefill->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
  tframefill->GetXaxis()->CenterTitle(1);
  tframefill->GetYaxis()->SetTitle("Events / (0.088 mm)");
  tframefill->SetMaximum(tframefill->GetMaximum()*9); 
  tframefill->SetMinimum(0.5); 
  tframefill->Draw();

  t->SetTextSize(0.05);
  t->DrawLatex(0.17,0.9,"CMS Preliminary");
  t->DrawLatex(0.17,0.83,"PbPb  #sqrt{s_{NN}} = 2.76 TeV"); 
  t->SetTextSize(0.04);
  t->DrawLatex(0.17,0.76,"L_{int} =  150 #mub^{-1}"); 
  t->SetTextSize(0.035);
  //  sprintf(reduceDS,"%.0f-%.0f%, %.1f < |y| < %.1f",cmin,cmax,ymin,ymax);
  t->DrawLatex(0.55,0.76,reduceDS);
  sprintf(reduceDS,"%.1f < p_{T} < %.1f GeV/c",pmin,pmax);
  t->DrawLatex(0.55,0.71,reduceDS);
  sprintf(reduceDS,"%.2f < |#phi_{J/#psi}-#Psi_{RP}| < %.2f",psmin,psmax);
  t->DrawLatex(0.55,0.66,reduceDS);
  sprintf(reduceDS,"EP: %s",rpmethod.c_str());
  t->DrawLatex(0.55,0.62,reduceDS);

  leg->SetX1NDC(0.57);
  leg->SetY1NDC(0.45);
  leg->SetX2NDC(0.92);
  leg->SetY2NDC(0.59);
  leg->Draw("same");

  titlestr = dirPre + "_rap" + yrange + "_pT" + prange + "_cent" + crange + "_dPhi" + phirange + "_timefit_Log_wopull.pdf";
  c2b->SaveAs(titlestr.c_str());
//  }// end of skip ctau fitting

/*  RooRealVar jpsiM = ws->var("Jpsi_Mass");
  RooRealVar jpsiCt = ws->var("Jpsi_Ct");
  RooDataHist hist2D("2DtotPDF","2D totPDF with mass and ctau",RooArgSet(jpsiM,jpsiCt),redDataCut);
  TH2F *htotPDF = (TH2F*)hist2D.createHistogram();
  TCanvas* c2c = new TCanvas("c2c","The Canvas",200,10,540,546);
  c2c->cd();
  htotPDF->Draw("COLZ");
  titlestr = dirPre + "_rap" + yrange + "_pT" + prange + "_cent" + crange + "_dPhi" + phirange + "_2D_" + partFile + "2D_totFit.pdf";
  c2c->SaveAs(titlestr.c_str());
*/

  // To check values of fit parameters
  cout << endl << "J/psi yields:" << endl;
  cout << "NSig :       Fit :"  << NSig_fin << " +/- " << ErrNSig_fin << endl;
  cout << "PROMPT :     Fit : " << NSigPR_fin << " +/- " << ErrNSigPR_fin << endl;
  cout << "NON-PROMPT : Fit : " << NSigNP_fin << " +/- " << ErrNSigNP_fin << endl;
  cout << "Bfraction : Fit : " << Bfrac_fin << " +/- " << ErrBfrac_fin << endl;
  cout << "Resolution : Fit : " << resol*1000. << " +/- " << Errresol*1000. << " mum" << endl;

  titlestr = dirPre + "_rap" + yrange + "_pT" + prange + "_cent" + crange + "_dPhi" + phirange + ".txt";

  ofstream outputFile(titlestr.c_str());
  if (!outputFile.good()) {cout << "Fail to open result file." << endl; return 1;}
  outputFile
  << "NSig "         << NSig_fin                          << " " << ErrNSig_fin << "\n"
  << "NBkg "         << NBkg_fin                          << " " << ErrNBkg_fin << "\n"
  << "coefExp "      << ws->var("coefExp")->getVal()      << " " << ws->var("coefExp")->getError() << "\n"
  << "coefPol1 "     << ws->var("coefPol1")->getVal()     << " " << ws->var("coefPol1")->getError() << "\n"
  << "coeffGaus "    << ws->var("coeffGaus")->getVal()    << " " << ws->var("coeffGaus")->getError() << "\n"
  << "meanSig1 "     << ws->var("meanSig1")->getVal()     << " " << ws->var("meanSig1")->getError() << "\n"
  << "sigmaSig1 "    << ws->var("sigmaSig1")->getVal()    << " " << ws->var("sigmaSig1")->getError() << "\n"
  << "sigmaSig2 "    << ws->var("sigmaSig2")->getVal()    << " " << ws->var("sigmaSig2")->getError()<< "\n"
  << "alpha "        << ws->var("alpha")->getVal()        << " " << ws->var("alpha")->getError() << "\n"
  << "enne "         << ws->var("enne")->getVal()         << " " << ws->var("enne")->getError() << "\n"
  << "enneW "        << ws->var("enneW")->getVal()        << " " << ws->var("enneW")->getError() << "\n";
  if (analyticBlifetime) {
    outputFile
    << "Gmc "        << ws->var("Gmc")->getVal()          << " " << ws->var("Gmc")->getError() <<  "\n"
    << "bTau "       << ws->var("bTau")->getVal()         << " " << ws->var("bTau")->getError() <<  "\n";
  } else {
    outputFile
    << "Gmc "        << "0"                               << " " << "0" <<  "\n"
    << "bTau "       << "0"                               << " " << "0" <<  "\n";
  }
  outputFile
  << "fracRes "      << ws->var("fracRes")->getVal()      << " " << ws->var("fracRes")->getError() <<  "\n"
  << "meanResSigW "  << ws->var("meanResSigW")->getVal()  << " " << ws->var("meanResSigW")->getError() << "\n"
  << "sigmaResSigW " << ws->var("sigmaResSigW")->getVal() << " " << ws->var("sigmaResSigW")->getError() <<  "\n"
  << "sigmaResSigN " << ws->var("sigmaResSigN")->getVal() << " " << ws->var("sigmaResSigN")->getError() << "\n"
  << "fLiving "      << ws->var("fLiving")->getVal()      << " " << ws->var("fLiving")->getError() << "\n"
  << "fpm "          << ws->var("fpm")->getVal()          << " " << ws->var("fpm")->getError() << "\n"
  << "fbkgCtTot "    << ws->var("fbkgCtTot")->getVal()    << " " << ws->var("fbkgCtTot")->getError() << "\n"
  << "lambdam "      << ws->var("lambdam")->getVal()      << " " << ws->var("lambdam")->getError() << "\n"
  << "lambdap "      << ws->var("lambdap")->getVal()      << " " << ws->var("lambdap")->getError() << "\n"
  << "lambdasym "    << ws->var("lambdasym")->getVal()    << " " << ws->var("lambdasym")->getError() << "\n"
//  << "PROMPT "       <<"0 0" << endl
//  << "NON-PROMPT "   <<"0 0" << endl;
  << "NLL "          << theNLL                            << endl
  << "PROMPT "       << NSigPR_fin                        << " " << ErrNSigPR_fin << endl
  << "NON-PROMPT "   << NSigNP_fin                        << " " << ErrNSigNP_fin << endl
  << "Bfraction "    << Bfrac_fin                         << " " << ErrBfrac_fin << endl
  << "Resolution "   << resol*1000.                       << " " << Errresol*1000. << endl;

/*  << "UnNormChi2 "   << UnNormChi2                        << " " << Dof << endl
  << "All "          << inteAll->getVal()                 << endl
  << "Signal "       << inteSig->getVal()                 << endl;
//  << "nFullBinsResid "<< nFullBinsResid                   << endl
//  << "RSS "          << RSS                               << endl;
*/
  
/*  titlestr = dirPre + "_rap" + yrange + "_pT" + prange + "_cent" + crange + "_dPhi" + phirange + "_workspace.root";
  if (ws->importClassCode("RooHistPdfConv")) cout << "importClassCode() succeed" << endl;
  else cout << "importClassCode() failed" << endl;
  if (ws->writeToFile(titlestr.c_str(),kTRUE)) cout << "Workspace saved" << endl;
  else cout << "Workspace NOT saved" << endl;*/

  fInMC.Close();
  fInData.Close();
  resultF.Close();
  outputFile.close();

/*  string resultRoot;
  resultRoot = dirPre + "_rap" + yrange + "_pT" + prange + "_cent" + crange + "_dPhi" + phirange + "_workspace.root";
  Bool_t ok = ws->writeToFile(resultRoot.c_str());
  if (!ok) { cout << "CANNOT write on workspace.root file\n"; }*/

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

void setWSRange(RooWorkspace *ws, float lmin, float lmax, float errmin, float errmax){
  float minRangeForPF = -4*errmax;
  if (minRangeForPF < -lmin) minRangeForPF = -lmin;

  ws->var("Jpsi_CtTrue")->setRange(-0.1,3.5);
  ws->var("Jpsi_Ct")->setRange("promptfit",minRangeForPF,4*errmax);
  ws->var("Jpsi_Ct")->setRange(-lmin,lmax);
  ws->var("Jpsi_CtErr")->setRange(errmin,errmax);

  ws->cat("Jpsi_Type")->setRange("glbglb","GG");
  ws->cat("MCType")->setRange("prompt","PR");
  ws->cat("MCType")->setRange("nonprompt","NP");
  ws->cat("MCType")->setRange("bkg","BK");

  return;
}

RooBinning setCtBinning(float lmin,float lmax) {
  RooBinning rb2(-lmin,lmax);

  if (lmax+lmin>4.9) {
    rb2.addBoundary(-1.5);
    rb2.addBoundary(-1.0);
    rb2.addBoundary(-0.8);
    rb2.addBoundary(-0.6);
    rb2.addBoundary(-0.5);
    rb2.addUniform(6,-0.5,-0.2);
    rb2.addUniform(18,-0.2,0.2);
    rb2.addUniform(9,0.2,0.5);
    rb2.addUniform(5,0.5,1.0);
    rb2.addUniform(15,1.0,lmax);
  } else if (lmax+lmin > 4.4) {
    rb2.addBoundary(-1.5);
    rb2.addBoundary(-1.0);
    rb2.addBoundary(-0.8);
    rb2.addBoundary(-0.6);
    rb2.addBoundary(-0.5);
    rb2.addUniform(9,-0.5,-0.2);
    rb2.addUniform(36,-0.2,0.2);
    rb2.addUniform(12,0.2,0.5);
    rb2.addUniform(5,0.5,1.0);
    rb2.addUniform(5,1.0,lmax);
  } else {
    rb2.addBoundary(-1.0);
    rb2.addBoundary(-0.7);
    rb2.addBoundary(-0.6);
    rb2.addBoundary(-0.5);
    rb2.addUniform(9,-0.5,-0.2);
    rb2.addUniform(40,-0.2,0.2);
    rb2.addUniform(12,0.2,0.5);
    rb2.addUniform(10,0.5,1.0);
    rb2.addUniform(4,1.0,lmax);
  }
  return rb2;
}

void defineMassBkg(RooWorkspace *ws) {
  // 1st order polynomial
  ws->factory("Polynomial::polFunct(Jpsi_Mass,{coefPol1[-0.05,-150.,150.]})");
  // Exponential
  ws->factory("Exponential::expFunct(Jpsi_Mass,coefExp[-1.,-3.,1.])");

  return;
}

void defineMassSig(RooWorkspace *ws) {
  //////// Candidates for signal
  // Normal gaussians
  ws->factory("Gaussian::signalG1(Jpsi_Mass,meanSig1[3.0975,3.05,3.15],sigmaSig1[0.1,0.04,0.2])");
//  ws->factory("Gaussian::signalG1(Jpsi_Mass,meanSig1[3.0975,3.05,3.15],sigmaSig1[0.1,0.06,0.2])");
//  ws->factory("Gaussian::signalG2(Jpsi_Mass,meanSig2[3.0975,3.05,3.15],sigmaSig2[0.03,0.008,0.2])");
//  ws->factory("Gaussian::signalGm1s2(Jpsi_Mass,meanSig1,sigmaSig2)");

  // Crystall Ball
  ws->factory("CBShape::signalCB(Jpsi_Mass,meanSig1,sigmaSig1,alpha[1.,0.,3.],enne[5.,1.,30.])");
  ws->factory("CBShape::signalCB2(Jpsi_Mass,meanSig1,sigmaSig2[0.03,0.01,0.04],alpha,enne)");
  ws->factory("CBShape::signalCBWN(Jpsi_Mass,meanSig1,sigmaSig1,alpha,enneW[5.,1.,50.])");
  ws->factory("CBShape::signalCB2WN(Jpsi_Mass,meanSig1,sigmaSig2,alpha,enneW)");

  //////// Sum of signal functions
  // Sum of gaussian 1 and a crystall ball
  ws->factory("SUM::sigCBG1(coeffGaus[0.1,0.0,0.93]*signalG1,signalCB)");
  // Sum of gaussian 1 and crystall ball 2
  ws->factory("SUM::sigCB2G1(coeffGaus*signalG1,signalCB2)");
  // Sum of gaussian 1 and crystall ball with wide n
  ws->factory("SUM::sigCBWNG1(coeffGaus*signalG1,signalCBWN)");
  // Sum of gaussian 1 and crystall ball 2 with wide n
  ws->factory("SUM::sigCB2WNG1(coeffGaus*signalG1,signalCB2WN)");
  return;
}

RooDataHist* subtractSidebands(RooWorkspace* ws, RooDataHist* subtrData, RooDataHist* all, RooDataHist* side, float scalefactor, string varName = "Jpsi_CtErr") {
  const RooArgSet* aRow;
  const RooArgSet* aRowS;
 
  if (all->numEntries() != side->numEntries()) {
    cout << "ERROR subtractSidebands : different binning!" << endl;
    return 0;
  }

  RooDataHist* scaledBkg = new RooDataHist("scaledBkg","Scale applied sideband data",RooArgSet(*(ws->var(varName.c_str())))); 

  for (Int_t i=0; i<all->numEntries(); i++) {
    aRow = all->get(i);
    aRowS = side->get(i);
    RooRealVar* thisVar = (RooRealVar*)aRow->find(varName.c_str());
    ws->var(varName.c_str())->setVal(thisVar->getVal());
    float sfBkg = scalefactor*side->weight(*aRowS,0,false);
    float newWeight = all->weight(*aRow,0,false) - scalefactor*side->weight(*aRowS,0,false);
    if (newWeight <= 2.0) newWeight = 2.0;
    subtrData->add(RooArgSet(*(ws->var(varName.c_str()))),newWeight);
    scaledBkg->add(RooArgSet(*(ws->var(varName.c_str()))),sfBkg);
  }

  return scaledBkg;
}


void getMCTrueLifetime(RooWorkspace *ws, RooDataSet *redMCCutNP, float *bgmcVal, float *bctauVal, string titlestr) {
  ws->pdf("bMCTrue")->fitTo(*redMCCutNP,Minos(0),SumW2Error(kTRUE),NumCPU(6));
  *bgmcVal = ws->var("Gmc")->getVal();
  *bctauVal = ws->var("bTau")->getVal();

  // *** Draw MC true Lifetime fit
  RooPlot *trueframef = ws->var("Jpsi_CtTrue")->frame();
  redMCCutNP->plotOn(trueframef);
  ws->pdf("bMCTrue")->plotOn(trueframef,LineColor(kBlue),Normalization(redMCCutNP->sumEntries(),RooAbsReal::NumEvent));

  TCanvas c0f; c0f.cd(); c0f.SetLogy(0); trueframef->Draw();
  string titlestr_lin = titlestr + "_Lin.pdf";
  c0f.SaveAs(titlestr_lin.c_str());
  c0f.cd(); c0f.SetLogy(1); trueframef->Draw();
  titlestr_lin = titlestr + "_Log.pdf";
  c0f.SaveAs(titlestr_lin.c_str());

  return ;
}


void defineCTResol(RooWorkspace *ws) {
  if (oneGaussianResol) {
    ws->factory("GaussModel::sigPR(Jpsi_Ct,meanResSigW[0.,-0.01,0.01],sigmaResSigN[0.8,0.6,2.0],one[1.0],Jpsi_CtErr)");
  } else {
    ws->factory("GaussModel::resGW(Jpsi_Ct,meanResSigW[0.,-0.01,0.01],sigmaResSigW[2.3,1.3,3.5],one[1.0],Jpsi_CtErr)");
    ws->factory("GaussModel::resGN(Jpsi_Ct,meanResSigW,sigmaResSigN[0.8,0.6,1.1],one,Jpsi_CtErr)");
    ws->factory("AddModel::sigPR({resGW,resGN},{fracRes[0.05,0.001,0.3]})");
  }

  return;
}

void defineCTBkg(RooWorkspace *ws) {
  ws->factory("Decay::bkg2(Jpsi_Ct,lambdap[0.42,0.05,1.5],sigPR,RooDecay::SingleSided)");
  ws->factory("Decay::bkg3(Jpsi_Ct,lambdam[0.79,0.02,1.5],sigPR,RooDecay::Flipped)");
  ws->factory("Decay::bkg4(Jpsi_Ct,lambdasym[0.69,0.02,5.0],sigPR,RooDecay::DoubleSided)");

  ws->factory("SUM::bkgPart1(fpm[1.0,0.0,1.0]*bkg2,bkg3)");
  ws->factory("SUM::bkgPart2(fLiving[0.9,0.0,1.0]*bkgPart1,bkg4)");
  ws->factory("SUM::bkgCtTot(fbkgCtTot[0.29,0.0,1.0]*sigPR,bkgPart2)");

  return;
}

void defineCTSig(RooWorkspace *ws, RooDataSet *redMCCutNP, string titlestr) {
  if (analyticBlifetime) {

    float GmcVal, bTauVal;
    ws->factory("GaussModel::bresGTrue(Jpsi_CtTrue,mean[0.0],Gmc[0.002,0.00001,0.02])");
    ws->factory("Decay::bMCTrue(Jpsi_CtTrue,bTau[0.04,0.01,1.0],bresGTrue,RooDecay::SingleSided)");

    getMCTrueLifetime(ws, redMCCutNP, &GmcVal, &bTauVal, titlestr);
    RooRealVar gmc("gmc","Sigma of MC Gaussian",GmcVal);
    RooRealVar btauFix1("btauFix1","Slope of MC exponential",bTauVal);   ws->import(btauFix1);
    RooFormulaVar bResSigN("bResSigN", "sqrt((@0*@1)**2+(@2)**2)", RooArgList(*(ws->var("sigmaResSigN")), *(ws->var("Jpsi_CtErr")),gmc));  ws->import(bResSigN);
    if (oneGaussianResol) {
      ws->factory("GaussModel::bresG(Jpsi_Ct,meanResSigW,bResSigN)");
    } else {
      RooFormulaVar bResSigW("bResSigW", "sqrt((@0*@1)**2+(@2)**2)", RooArgList(*(ws->var("sigmaResSigW")), *(ws->var("Jpsi_CtErr")),gmc));  ws->import(bResSigW);
      
      ws->factory("GaussModel::bresGN(Jpsi_Ct,meanResSigW,bResSigN)");
      ws->factory("GaussModel::bresGW(Jpsi_Ct,meanResSigW,bResSigW)");
      ws->factory("AddModel::bresG({bresGW,bresGN},{fracRes})");
    }
    // fix tau_B
    // ws->factory("Decay::sigNP1(Jpsi_Ct,btauFix1,bresG,RooDecay::SingleSided)");
    // ws->factory("Decay::sigNP2(Jpsi_Ct,btauFix2,bresG,RooDecay::SingleSided)");
    // float tau_B
    ws->factory("Decay::sigNP(Jpsi_Ct,bTau,bresG,RooDecay::SingleSided)");
    
  } else {

    RooDataHist* binMCCutNP = new RooDataHist("binMCCutNP","MC distribution for NP signal",RooArgSet(*(ws->var("Jpsi_CtTrue"))),*redMCCutNP);
    if (oneGaussianResol) {
      RooHistPdfConv sigNP("sigNP","Non-prompt signal with narrow gaussian",*(ws->var("Jpsi_Ct")),*(ws->var("meanResSigW")),*(ws->var("sigmaResSigN")),*(ws->var("one")),*(ws->var("Jpsi_CtErr")),*binMCCutNP);  ws->import(sigNP);
    } else {
      RooHistPdfConv sigNPW("sigNPW","Non-prompt signal with wide gaussian",*(ws->var("Jpsi_Ct")),*(ws->var("meanResSigW")),*(ws->var("sigmaResSigW")),*(ws->var("one")),*(ws->var("Jpsi_CtErr")),*binMCCutNP);  ws->import(sigNPW);
      RooHistPdfConv sigNPN("sigNPN","Non-prompt signal with narrow gaussian",*(ws->var("Jpsi_Ct")),*(ws->var("meanResSigW")),*(ws->var("sigmaResSigN")),*(ws->var("one")),*(ws->var("Jpsi_CtErr")),*binMCCutNP);  ws->import(sigNPN);
      RooAddPdf sigNP("sigNP","Non-prompt signal",RooArgSet(*(ws->pdf("sigNPW")),*(ws->pdf("sigNPN"))),RooArgSet(*(ws->var("fracRes"))));  ws->import(sigNP); 
    }

    // *** Draw MC true Lifetime fit
    RooPlot *trueframef = ws->var("Jpsi_CtTrue")->frame();
    redMCCutNP->plotOn(trueframef);
    ws->pdf("sigNP")->plotOn(trueframef,LineColor(kBlue),Normalization(redMCCutNP->sumEntries(),RooAbsReal::NumEvent));

    TCanvas c0f; c0f.cd(); c0f.SetLogy(0); trueframef->Draw();
    string titlestr_lin = titlestr + "_Lin.pdf";
    c0f.SaveAs(titlestr_lin.c_str());
    c0f.cd(); c0f.SetLogy(1); trueframef->Draw();
    titlestr_lin = titlestr + "_Log.pdf";
    c0f.SaveAs(titlestr_lin.c_str());

  }
  
  return;
}

