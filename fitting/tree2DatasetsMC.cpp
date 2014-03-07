#include <iostream>
#include <fstream>

#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include "TH2.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TChain.h"

#include "RooFit.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooCategory.h"

#include <TCanvas.h>
#include "TH2.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include <TVector3.h>
#include <math.h>
#include <TPaveStats.h>

using namespace RooFit;
using namespace std;


bool isAccept(const TLorentzVector* aMuon) {
  // *USE* muon kinematical cuts (eta dependent momentum / pT cuts )
  /*   return (fabs(aMuon->Eta()) < 2.4 &&
       ((fabs(aMuon->Eta()) < 1.3 && aMuon->Pt() > 3.3) ||
       (fabs(aMuon->Eta()) > 1.3 && fabs(aMuon->Eta()) < 2.2 && aMuon->P() > 2.9) ||
       (fabs(aMuon->Eta()) > 2.2 && aMuon->Pt() > 0.8)));
  */
  // *REMOVE* muon kinematical cuts (eta dependent momentum / pT cuts )
  // by just returning TRUE
  return true;
}

double CorrectMass(const TLorentzVector* mu1,const TLorentzVector* mu2, int mode){  
  double CMass=0;
  const double mumass=0.105658;
  double k1,k2;
  double pt1=mu1->Pt();
  double pt2=mu2->Pt();
  double eta1=mu1->Eta();
  double eta2=mu2->Eta();
  if (mode==1){
    k1=1.0009;//constant scale correction
    k2=1.0009;
  }
  if (mode==2){
    k1=1.0019-0.0004*pt1;
    k2=1.0019-0.0004*pt2; // pt dependent correction
  }
  if (mode==3){
    double a0=0.00038; //3.8 * pow(10,-4);
    double a1=0.0;
    double a2=0.0003; //3.0 * pow(10,-4);
    double a3=0.0;

    k1=1+a0+a1*fabs(eta1)+a2*eta1*eta1+a3*pt1;
    k2=1+a0+a1*fabs(eta2)+a2*eta2*eta2+a3*pt2;// pt and eta dependent
  }

  if (mode == 4){
    double a0=1.002;
    double a1=-0.002;
    double a2=0.001;
    double a3=-0.0001;

    k1=a0+a1*fabs(eta1)+a2*eta1*eta1+a3*pt1;
    k2=a0+a1*fabs(eta2)+a2*eta2*eta2+a3*pt2;// pt and eta dependent
  }

  TVector3 mom1=mu1->Vect();
  TVector3 mom2=mu2->Vect();
  mom1=k1*mom1; 
  mom2=k2*mom2;
  double E1=sqrt(mom1.Mag2()+(mumass*mumass));
  double E2=sqrt(mom2.Mag2()+(mumass*mumass));
  TVector3 momtot=mom1+mom2;
  CMass=sqrt((E1+E2)*(E1+E2)-momtot.Mag2());
  return CMass;
}



int main(int argc, char* argv[]) {

  const double Jpsi_MassMin=2.2;
  const double Jpsi_MassMax=3.4;
  const double Jpsi_PtMin=0;
  const double Jpsi_PtMax=50;
  const double Jpsi_YMin=0;
  const double Jpsi_YMax=2.4;
  const double Jpsi_CtMin = -3.0;
  const double Jpsi_CtMax = 3.5;

  bool isHI = true;

  string fileName;
  if ( argc != 5 && argc != 3 ){
    char msg[300];
    sprintf(msg,"Usage1: %s [input file] [output file directory] [start event #] [end event # (-1 for total)]",argv[0]);
    cout << msg << endl; 
    sprintf(msg,"Usage2: %s [input file] [output file directory]",argv[0]);
    cout << msg << endl; 
    return 1;
  }

  fileName = argv[1];

  if (fileName.find("pp")!=string::npos)
    isPbPb = false;


  TFile *file= TFile::Open(fileName);
  TTree * Tree=(TTree*)file->Get("myTree");

  UInt_t          eventNb;
  //  Bool_t          goodEvent;
  Int_t           Centrality;
  Int_t           Reco_QQ_size;
  Int_t           Reco_QQ_type[20];   //[Reco_QQ_size]
  Int_t           Reco_QQ_sign[20];   //[Reco_QQ_size]
  Int_t           Reco_QQ_trig[20];   //[Reco_QQ_size]
  Int_t           Gen_QQ_size;
  Int_t           Gen_QQ_type[20];
  TClonesArray    *Gen_QQ_4mom;
  TClonesArray    *Reco_QQ_4mom;
  TClonesArray    *Reco_QQ_mupl_4mom;
  TClonesArray    *Reco_QQ_mumi_4mom;
  // Bool_t           Reco_QQ_mupl_acc[20];   //[Reco_QQ_size]
  // Bool_t           Reco_QQ_mumi_acc[20];   //[Reco_QQ_size]
  // Bool_t           Reco_QQ_mupl_cuts[20];   //[Reco_QQ_size]
  // Bool_t           Reco_QQ_mumi_cuts[20];   //[Reco_QQ_size]
  Float_t         Reco_QQ_ctau[20];   //[Reco_QQ_size]
  Float_t         Reco_QQ_ctauErr[20];   //[Reco_QQ_size]
  Float_t         Reco_QQ_ctauTrue[20];   //[Reco_QQ_size]
  Float_t         Reco_QQ_VtxProb[20];   //[Reco_QQ_size]
  //  Float_t         rpAng[38];   //[Reco_QQ_size]

  TBranch        *b_eventNb;
  //  TBranch        *b_goodEvent;   //!
  TBranch        *b_Centrality;   //!
  TBranch        *b_Reco_QQ_size;   //!
  TBranch        *b_Reco_QQ_type;   //!
  TBranch        *b_Reco_QQ_sign;   //!
  TBranch        *b_HLTriggers;   //!
  TBranch        *b_Reco_QQ_trig;   //!
  TBranch        *b_Gen_QQ_size;   //!
  TBranch        *b_Gen_QQ_type;
  TBranch        *b_Gen_QQ_4mom;   //!
  TBranch        *b_Reco_QQ_4mom;   //!
  TBranch        *b_Reco_QQ_mupl_4mom;   //!
  TBranch        *b_Reco_QQ_mumi_4mom;   //!
  // TBranch        *b_Reco_QQ_mupl_acc;   //!
  // TBranch        *b_Reco_QQ_mumi_acc;   //!
  // TBranch        *b_Reco_QQ_mupl_cuts;   //!
  // TBranch        *b_Reco_QQ_mumi_cuts;   //!
  TBranch        *b_Reco_QQ_ctau;   //!
  TBranch        *b_Reco_QQ_ctauErr;   //!
  TBranch        *b_Reco_QQ_ctauTrue;   //!
  TBranch        *b_Reco_QQ_VtxProb;   //!
  //  TBranch        *b_rpAng;   //!

  TLorentzVector* JP= new TLorentzVector;
  TLorentzVector* m1P= new TLorentzVector;
  TLorentzVector* m2P= new TLorentzVector;

  double vprob, theCt, theCtErr, theCtTrue;
  int HLTriggers,theCat,Jq,genType;

  static const unsigned int centRegions = 3;
  //  int centLimits[centRegions+1] = {0, 8, 16, 24, 40}; // 0 20 40 60 100
  int centLimits[centRegions+1] = {0, 8, 16, 40}; // 0 20 40 100
  //  int centLimits[centRegions+1] = {0, 4, 8, 20, 40}; // 0 10 20 50 100
  //  int centLimits[centRegions+1] = {0, 8, 20, 40}; // 0 20 50 100
  //  int centLimits[centRegions+1] = {0, 4, 12, 24, 40};  //0 10 30 60 100
  //  int centLimits[centRegions+1] = {0, 2, 4, 12, 24, 40};  //0 5 10 30 60 100
  //  int centLimits[centRegions+1] = {0, 4, 40};  //0 10
  //  int centLimits[centRegions+1] = {4, 24};  //10 60
  //  int centLimits[centRegions+1] = {0, 2, 4};  //0 5 10

  //  static const unsigned int rapRegions = 3;
  //  float rapLimits[rapRegions+1] = {Jpsi_YMin, Jpsi_YMax};
  //  float rapLimits[rapRegions+1] = {Jpsi_YMin,1.2,1.6,Jpsi_YMax};
  //   static const unsigned int rapRegions = 1;
  //   float rapLimits[rapRegions+1] = {Jpsi_YMin,Jpsi_YMax};

  TH1D *hSig[centRegions+1];
  TH1D *hBkg[centRegions+1];
  RooDataSet* dataJpsi[centRegions+1];
  RooDataSet* dataJpsiSame[centRegions+1];
  RooDataSet* dataPsip[centRegions+1];
  RooRealVar* Jpsi_Mass;
  RooRealVar* Psip_Mass;      
  RooRealVar* Jpsi_Pt;
  RooRealVar* Jpsi_Ct;
  RooRealVar* Jpsi_CtErr;
  RooRealVar* Jpsi_CtTrue;
  RooRealVar* Jpsi_Y;
  RooCategory* Jpsi_Type;
  RooCategory* Jpsi_Sign;
  RooCategory* MCType;
  RooRealVar* Gen_Pt;

  Jpsi_Mass = new RooRealVar("Jpsi_Mass","J/#psi mass",Jpsi_MassMin,Jpsi_MassMax,"GeV/c^{2}");
  Psip_Mass = new RooRealVar("Psip_Mass","#psi' mass",3.3,Jpsi_MassMax,"GeV/c^{2}");
  Jpsi_Pt = new RooRealVar("Jpsi_Pt","J/#psi pt",Jpsi_PtMin,Jpsi_PtMax,"GeV/c");
  Jpsi_Y = new RooRealVar("Jpsi_Y","J/#psi y",-Jpsi_YMax,Jpsi_YMax);
  Jpsi_Type = new RooCategory("Jpsi_Type","Category of Jpsi_");
  Jpsi_Sign = new RooCategory("Jpsi_Sign","Charge combination of Jpsi_");
  MCType = new RooCategory("MCType","Type of generated Jpsi_");
  Jpsi_Ct = new RooRealVar("Jpsi_Ct","J/#psi c#tau",Jpsi_CtMin,Jpsi_CtMax,"mm");
  Jpsi_CtErr = new RooRealVar("Jpsi_CtErr","J/#psi c#tau error",-1.,1.,"mm");
  Jpsi_CtTrue = new RooRealVar("Jpsi_CtTrue","J/#psi c#tau true",Jpsi_CtMin,Jpsi_CtMax,"mm");
  Gen_Pt = new RooRealVar("Gen_Pt","Generated J/#psi pt",Jpsi_PtMin,Jpsi_PtMax,"GeV/c");

  Jpsi_Type->defineType("GG",0);
  Jpsi_Type->defineType("GT",1);
  Jpsi_Type->defineType("TT",2);

  Jpsi_Sign->defineType("OS",0);
  Jpsi_Sign->defineType("PP",1);
  Jpsi_Sign->defineType("MM",2);

  MCType->defineType("PR",0);
  MCType->defineType("NP",1);

  Reco_QQ_4mom = 0;
  Reco_QQ_mupl_4mom = 0;
  Reco_QQ_mumi_4mom = 0;

  Tree->SetBranchAddress("eventNb", &eventNb, &b_eventNb);
  //  Tree->SetBranchAddress("goodEvent", &goodEvent, &b_goodEvent);
  Tree->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
  //  Tree->SetBranchAddress("NfRpAng", rpAng, &b_rpAng);
  //  Tree->SetBranchAddress("rpAng", rpAng, &b_rpAng);   //Flatten reaction plane
  Tree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
  Tree->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
  Tree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
  Tree->SetBranchAddress("Reco_QQ_type", Reco_QQ_type, &b_Reco_QQ_type);
  Tree->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
  Tree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
  Tree->SetBranchAddress("Gen_QQ_type", Gen_QQ_type, &b_Gen_QQ_type);
  Tree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
  Tree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
  Tree->SetBranchAddress("Reco_QQ_mupl_4mom", &Reco_QQ_mupl_4mom, &b_Reco_QQ_mupl_4mom);
  Tree->SetBranchAddress("Reco_QQ_mumi_4mom", &Reco_QQ_mumi_4mom, &b_Reco_QQ_mumi_4mom);
  // Tree->SetBranchAddress("Reco_QQ_mupl_acc", Reco_QQ_mupl_acc, &b_Reco_QQ_mupl_acc);
  // Tree->SetBranchAddress("Reco_QQ_mumi_acc", Reco_QQ_mumi_acc, &b_Reco_QQ_mumi_acc);
  // Tree->SetBranchAddress("Reco_QQ_mupl_cuts", Reco_QQ_mupl_cuts, &b_Reco_QQ_mupl_cuts);
  // Tree->SetBranchAddress("Reco_QQ_mumi_cuts", Reco_QQ_mumi_cuts, &b_Reco_QQ_mumi_cuts);


  Tree->SetBranchAddress("Reco_QQ_ctau", Reco_QQ_ctau, &b_Reco_QQ_ctau);
  Tree->SetBranchAddress("Reco_QQ_ctauErr", Reco_QQ_ctauErr, &b_Reco_QQ_ctauErr);
  Tree->SetBranchAddress("Reco_QQ_ctauTrue", Reco_QQ_ctauTrue, &b_Reco_QQ_ctauTrue);
  Tree->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);
  /*
    RooArgList varlist(*Jpsi_Mass,*Jpsi_Pt,*Jpsi_Y,*Jpsi_Type,*Jpsi_Sign,*Jpsi_Ct,*Jpsi_CtErr);
    RooArgList varlistSame(*Jpsi_Mass,*Jpsi_Pt,*Jpsi_Y,*Jpsi_Type,*Jpsi_Sign,*Jpsi_Ct,*Jpsi_CtErr);
    RooArgList varlist2(*Psip_Mass,*Jpsi_Pt,*Jpsi_Y,*Jpsi_Type,*Jpsi_Sign,*Jpsi_Ct,*Jpsi_CtErr);
  */
  RooArgList varlist(*Jpsi_Mass,*Jpsi_Pt,*Jpsi_Y,*Jpsi_Type,*Jpsi_Sign,*MCType,*Jpsi_Ct,*Jpsi_CtErr,*Jpsi_CtTrue);
  RooArgList varlistSame(*Jpsi_Mass,*Jpsi_Pt,*Jpsi_Y,*Jpsi_Type,*Jpsi_Sign,*MCType,*Jpsi_Ct,*Jpsi_CtErr,*Jpsi_CtTrue);
  RooArgList varlist2(*Psip_Mass,*Jpsi_Pt,*Jpsi_Y,*Jpsi_Type,*Jpsi_Sign,*MCType,*Jpsi_Ct,*Jpsi_CtErr,*Jpsi_CtTrue);
  
  varlist.add(*Gen_Pt);
  varlistSame.add(*Gen_Pt);
  varlist2.add(*Gen_Pt);

  if (isHI) {
    for (unsigned int j = 0; j <= centRegions; j++) {
      char namefile[200] = {0};
      if (j==centRegions) {
	sprintf(namefile,"cent%d-%d",
		int(centLimits[0]*2.5),int(centLimits[j]*2.5));
      }
      else {
	sprintf(namefile,"cent%d-%d",
		int(centLimits[j]*2.5),int(centLimits[j+1]*2.5));
      }
      hSig[j] = new TH1D(namefile,namefile,50,-0.5,2);
      hSig[j]->GetXaxis()->SetTitle("#sigma(#font[12]{l}_{J/#psi})");
      strcat(namefile,"_nbkg");
      hBkg[j] = new TH1D(namefile,namefile,50,-0.5,2);
      hBkg[j]->GetXaxis()->SetTitle("#sigma(#font[12]{l}_{J/#psi})");
      dataJpsi[j] = new RooDataSet("dataJpsi","A sample",varlist);
      dataJpsiSame[j] = new RooDataSet("dataJpsiSame","A sample",varlistSame);
      dataPsip[j] = new RooDataSet("dataPsip","A sample",varlist2);
    }
  }
  else {
    unsigned int j=0;
    char namefile[200] = {0};
    sprintf(namefile,"cent%d-%d",0,100);
    hSig[j] = new TH1D(namefile,namefile,50,-0.5,2);
    hSig[j]->GetXaxis()->SetTitle("#sigma(#font[12]{l}_{J/#psi})");
    strcat(namefile,"_nbkg");
    hBkg[j] = new TH1D(namefile,namefile,50,-0.5,2);
    hBkg[j]->GetXaxis()->SetTitle("#sigma(#font[12]{l}_{J/#psi})");
    dataJpsi[j] = new RooDataSet("dataJpsi","A sample",varlist);
    dataJpsiSame[j] = new RooDataSet("dataJpsiSame","A sample",varlistSame);
    dataPsip[j] = new RooDataSet("dataPsip","A sample",varlist2);
  }

  int initev = 0;
  int nevt =  Tree->GetEntries();
  if (argc == 5) {
    initev = atoi(argv[3]);
    if (atoi(argv[4]) != -1)
      nevt = atoi(argv[4]);
  }

  TH1D *JpsiPt = new TH1D("JpsiPt","JpsiPt;p_{T} (J/#psi)",35,0,140);

  int PassingEvent[centRegions+1] = {0};
  for (int ev=initev; ev<nevt; ++ev) {
    if (ev%50000==0) cout << ">>>>> EVENT " << ev << " / " << Tree->GetEntries() <<  endl;

    Tree->GetEntry(ev);

    int theCentrality=Centrality;
    //    float theRPAng = rpAng[RPNUM];

    // MC
    //    if (!goodEvent) continue;

    //    if (Reco_QQ_size != Gen_QQ_size) continue;
    //    cout << Reco_QQ_size << endl;
    for (int i=0; i<Reco_QQ_size; ++i) {
      JP = (TLorentzVector*) Reco_QQ_4mom->At(i);
      m1P = (TLorentzVector*) Reco_QQ_mupl_4mom->At(i);
      m2P = (TLorentzVector*) Reco_QQ_mumi_4mom->At(i);
      vprob = Reco_QQ_VtxProb[i];
      theCat = Reco_QQ_type[i];
      Jq = Reco_QQ_sign[i];
      genType = Gen_QQ_type[i];
      theCt = Reco_QQ_ctau[i];
      theCtErr = Reco_QQ_ctauErr[i];
      theCtTrue = Reco_QQ_ctauTrue[i];

      double theMass =JP->M();
      double theRapidity=JP->Rapidity();
      double thePt=JP->Pt();
      double theGenPt;
      if (Gen_QQ_size>0)
	theGenPt = ((TLorentzVector*)Gen_QQ_4mom->At(0))->Pt();
      else
	theGenPt = 40.0;

      bool ok1=isAccept(m1P);
      bool ok2=isAccept(m2P);
     
      if (theMass>3.350) {
	theCt*=3.686109/3.096916; // use proper mass for l_psi(2S)
	theCtErr*=3.686109/3.096916; // use proper mass for error on l_psi(2S)
      }

      bool isTriggered = false;
      if (isPbPb)
	isTriggered = ((HLTriggers&1)==1 && (Reco_QQ_trig[i]&1)==1); // HLT_HIL1DoubleMu0_HighQ_v*
      else
	isTriggered = ((HLTriggers&2)==2 && (Reco_QQ_trig[i]&2)==2); // HLT_PAL1DoubleMu0_HighQ_v1

      if (theMass > Jpsi_MassMin && theMass < Jpsi_MassMax && 
	  //          theCt > Jpsi_CtMin && theCt < Jpsi_CtMax && 
	  //          thePt > Jpsi_PtMin && thePt < Jpsi_PtMax && 
	  //          fabs(theRapidity) > Jpsi_YMin && fabs(theRapidity) < Jpsi_YMax &&
	  //          theRPAng != -10 &&
	  ok2 && ok1 &&
	  isTriggered
	  //MC
	  //	  && (Reco_QQ_mupl_acc && Reco_QQ_mumi_acc && Reco_QQ_mupl_cuts && Reco_QQ_mumi_cuts)
          ) {

        JpsiPt->Fill(thePt);
        Jpsi_Pt->setVal(thePt); 
        Jpsi_Y->setVal(theRapidity); 
        Jpsi_Mass->setVal(theMass);
        Psip_Mass->setVal(theMass);
        Jpsi_Ct->setVal(theCt);
        Jpsi_CtErr->setVal(theCtErr);
	Jpsi_CtTrue->setVal(theCtTrue);
        Jpsi_Type->setIndex(theCat,kTRUE);
        Gen_Pt->setVal(theGenPt); 
        if (Jq == 0){ Jpsi_Sign->setIndex(Jq,kTRUE); }
        else { Jpsi_Sign->setIndex(Jq,kTRUE); }
	MCType->setIndex(0,kTRUE);//genType

	RooArgList varlist_tmp(*Jpsi_Mass,*Jpsi_Pt,*Jpsi_Y,*Jpsi_Type,*Jpsi_Sign,*MCType,*Jpsi_Ct,*Jpsi_CtErr,*Jpsi_CtTrue);
	RooArgList varlist2_tmp(*Psip_Mass,*Jpsi_Pt,*Jpsi_Y,*Jpsi_Type,*Jpsi_Sign,*MCType,*Jpsi_Ct,*Jpsi_CtErr,*Jpsi_CtTrue);
        // RooArgList varlist_tmp(*Jpsi_Mass,*Jpsi_Pt,*Jpsi_Y,*Jpsi_Type,*Jpsi_Sign,*Jpsi_Ct,*Jpsi_CtErr);
        // RooArgList varlist2_tmp(*Psip_Mass,*Jpsi_Pt,*Jpsi_Y,*Jpsi_Type,*Jpsi_Sign,*Jpsi_Ct,*Jpsi_CtErr);
	varlist_tmp.add(*Gen_Pt);
	varlist2_tmp.add(*Gen_Pt);

	if (isHI) {       
	  for (unsigned int j = 0; j <= centRegions; j++) {
	    if ( (j==centRegions && theCentrality < centLimits[j] && theCentrality >= centLimits[0]) ||
		 (theCentrality < centLimits[j+1] && theCentrality >= centLimits[j]) ) {
	      if (Jq == 0) {
		if (theMass < 4.2) {
		  dataJpsi[j]->add(varlist_tmp);
		  PassingEvent[j] = PassingEvent[j] + 1;
		}
		if (theMass > 3.35) {
		  dataPsip[j]->add(varlist2_tmp);
		}
		if (theMass > 2.9 && theMass < 3.35)
		  hSig[j]->Fill(theCtErr);
		else if (theMass < 2.9 || theMass > 3.35)
		  hBkg[j]->Fill(theCtErr);
	      } else {
		if (theMass < 4.2) {
		  dataJpsiSame[j]->add(varlist_tmp);
		}
	      }
	    }
	  } // End of j, i loop
	}
	else {
	  unsigned int j=0;
	  if (Jq == 0) {
	    if (theMass < 4.2) {
	      dataJpsi[j]->add(varlist_tmp);
	      PassingEvent[j] = PassingEvent[j] + 1;
	    }
	    if (theMass > 3.35) {
	      dataPsip[j]->add(varlist2_tmp);
	    }
	    if (theMass > 2.9 && theMass < 3.35)
	      hSig[j]->Fill(theCtErr);
	    else if (theMass < 2.9 || theMass > 3.35)
	      hBkg[j]->Fill(theCtErr);
	  } else {
	    if (theMass < 4.2) {
	      dataJpsiSame[j]->add(varlist_tmp);
	    }
	  }
	}
      } // End of if() statement for cuts
    } // End of Reco_QQ_size loop
  }


  char namefile[200];
  TCanvas *canv = new TCanvas("canv","canv",800,600);
  canv->cd();
  canv->SetLogy();
  JpsiPt->Draw("text");
  sprintf(namefile,"%s/JpsiPt.pdf",argv[2]);
  canv->SaveAs(namefile);

  delete canv;

  canv = new TCanvas("canv","canv",4000,3000);
  if (isHI) 
    canv->Divide(1,centRegions+1);
  canv->Draw();

  /// *** Fill TFiles with RooDataSet
  TFile* Out[centRegions+1];
  int padnum = 1;
  if (isHI) {
    for (unsigned int j = 0; j <= centRegions; j++) {
      if (j==centRegions) {
	sprintf(namefile,"%s/PbPbPromptJpsiMC_DblMu0_cent%d-%d_M2234.root",argv[2],
		int(centLimits[0]*2.5),int(centLimits[j]*2.5));
      }
      else {
	sprintf(namefile,"%s/PbPbPromptJpsiMC_DblMu0_cent%d-%d_M2234.root",argv[2],
		int(centLimits[j]*2.5),int(centLimits[j+1]*2.5));
      }
      Out[j] = new TFile(namefile,"RECREATE");
      Out[j]->cd();
      dataJpsi[j]->Write();
      dataJpsiSame[j]->Write();
      dataPsip[j]->Write();
      Out[j]->Close();
      cout << namefile << endl;
      cout << "PassingEvent[" << j << "]: " << PassingEvent[j] << endl;
      
      canv->cd(padnum);
      gPad->SetLogy(1);
      hSig[j]->SetLineColor(kBlack);
      hSig[j]->SetMaximum(hSig[j]->GetMaximum()*5);
      hSig[j]->Draw();
      hBkg[j]->SetMarkerColor(kRed);
      hBkg[j]->SetLineColor(kRed);
      hBkg[j]->SetMarkerSize(1.5);
      hBkg[j]->SetMarkerStyle(20);
      hBkg[j]->Draw("sames");

      gPad->Update();
      TPaveStats *s = (TPaveStats*)hSig[j]->FindObject("stats");
      s->SetX1NDC(s->GetX1NDC()-0.18);
      s->SetY1NDC(s->GetY2NDC()-0.3);
      TPaveStats *sr = (TPaveStats*)hBkg[j]->FindObject("stats");
      sr->SetX1NDC(sr->GetX1NDC()-0.18);
      sr->SetY2NDC(s->GetY1NDC()-0.03);
      sr->SetY1NDC(sr->GetY2NDC() - (s->GetY2NDC()-s->GetY1NDC()) );
      sr->SetLineColor(kRed);
      sr->SetTextColor(kRed);
      TText txt1;
      txt1.DrawTextNDC(0.2,0.85,"Signal");
      TText txt2;
      txt2.SetTextColor(kRed);
      txt2.DrawTextNDC(0.2,0.8,"Sideband");
      gPad->Update();
      padnum++;
    }
  }
  else {
    unsigned int j=0;
    sprintf(namefile,"%s/ppPromptJpsiMC_DblMu0_cent%d-%d_M2234.root",argv[2],0,100);
    Out[j] = new TFile(namefile,"RECREATE");
    Out[j]->cd();
    dataJpsi[j]->Write();
    dataJpsiSame[j]->Write();
    dataPsip[j]->Write();
    Out[j]->Close();
    cout << namefile << endl;
    cout << "PassingEvent[" << j << "]: " << PassingEvent[j] << endl;
      
    canv->cd(padnum);
    gPad->SetLogy(1);
    hSig[j]->SetLineColor(kBlack);
    hSig[j]->SetMaximum(hSig[j]->GetMaximum()*5);
    hSig[j]->Draw();
    hBkg[j]->SetMarkerColor(kRed);
    hBkg[j]->SetLineColor(kRed);
    hBkg[j]->SetMarkerSize(1.5);
    hBkg[j]->SetMarkerStyle(20);
    hBkg[j]->Draw("sames");

    gPad->Update();
    TPaveStats *s = (TPaveStats*)hSig[j]->FindObject("stats");
    s->SetX1NDC(s->GetX1NDC()-0.18);
    s->SetY1NDC(s->GetY2NDC()-0.3);
    TPaveStats *sr = (TPaveStats*)hBkg[j]->FindObject("stats");
    sr->SetX1NDC(sr->GetX1NDC()-0.18);
    sr->SetY2NDC(s->GetY1NDC()-0.03);
    sr->SetY1NDC(sr->GetY2NDC() - (s->GetY2NDC()-s->GetY1NDC()) );
    sr->SetLineColor(kRed);
    sr->SetTextColor(kRed);
    TText txt1;
    txt1.DrawTextNDC(0.2,0.85,"Signal");
    TText txt2;
    txt2.SetTextColor(kRed);
    txt2.DrawTextNDC(0.2,0.8,"Sideband");
    gPad->Update();
  }
  
  if (isHI)
    sprintf(namefile,"%s/PbPbCtauErr_centBins%d.pdf",argv[2],centRegions);
  else
    sprintf(namefile,"%s/ppCtauErr.pdf",argv[2]);

  canv->SaveAs(namefile);

  return 0;
}


