#include <iostream>
#include <string>

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1.h"
#include "TMath.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"

TF1 *f1;
TF1 *f2;
TF1 *f3;
TF1 *f4;
TF1 *f5;
TF1 *f5a;
TF1 *f6;
TF1 *f7;
TF1 *f8;

TFitResultPtr res_00100;
TFitResultPtr res_0010;
TFitResultPtr res_1020;
TFitResultPtr res_2030;
TFitResultPtr res_3040;
TFitResultPtr res_4050;
TFitResultPtr res_5060;
TFitResultPtr res_60100;
TFitResultPtr res_0020;
TFitResultPtr res_2040;
TFitResultPtr res_40100;

TFitResultPtr fitAndDrawHisto(TH1 *h, bool fitCBGsum);


const double binWidth = 0.020;
bool fitRatio = false;
bool fitBgPol2 = false;
bool fitBgPol1 = false;
bool fit2Expo = false;

double crystal_ball_expo(double *x, double *par)
{
  double N = par[0];

  double mean = par[1];
  double sigma = par[2];
  double alpha = par[3];
  double n = par[4];

  double A = pow(n/fabs(alpha),n)*exp(-pow(alpha,2)/2);
  double B = n/fabs(alpha) - fabs(alpha);

  double NBkg = par[5]*fabs(par[6]);

  if ((x[0]-mean)/sigma>-alpha)
    return binWidth*(N*TMath::Gaus(x[0],mean,sigma,1)+NBkg*TMath::Exp(x[0]*par[6]));
  else
    return binWidth*(N/(sqrt(TMath::TwoPi())*sigma)*A*pow(B-(x[0]-mean)/sigma,-n)+NBkg*TMath::Exp(x[0]*par[6]));
}

double two_crystal_ball_expo(double *x, double *par)
{
  double N1 = par[0];

  double mean = par[1];
  double sigma = par[2];
  double alpha = par[3];
  double n = par[4];

  double A = pow(n/fabs(alpha),n)*exp(-pow(alpha,2)/2);
  double B = n/fabs(alpha) - fabs(alpha);

  double massRatio = 3.686/3.0969;

  double N2;
  if (fitRatio)
    N2 = par[7]*N1;
  else
    N2 = par[7];

  double NBkg = par[5]*fabs(par[6]);

  if ((x[0]-mean)/sigma>-alpha)
    return binWidth*( N1*TMath::Gaus(x[0],mean,sigma,1) + 
		      N2*TMath::Gaus(x[0],mean*massRatio,sigma*massRatio,1) + 
		      NBkg*TMath::Exp(x[0]*par[6]) );
  else {
    return binWidth*( N1/(sqrt(TMath::TwoPi())*sigma)*A*pow(B-(x[0]-mean)/sigma,-n) + 
		      N2/(sqrt(TMath::TwoPi())*sigma*massRatio)*A*pow(B-(x[0]-mean*massRatio)/(sigma*massRatio),-n) + 
		      NBkg*TMath::Exp(x[0]*par[6]) );
  }
}

double two_crystal_ball_two_expo(double *x, double *par)
{
  double N1 = par[0];

  double mean = par[1];
  double sigma = par[2];
  double alpha = par[3];
  double n = par[4];

  double A = pow(n/fabs(alpha),n)*exp(-pow(alpha,2)/2);
  double B = n/fabs(alpha) - fabs(alpha);

  double massRatio = 3.686/3.0969;

  double N2;
  if (fitRatio)
    N2 = par[7]*N1;
  else
    N2 = par[7];

  double NBkg = fabs(par[5]*par[6]);
  double NBkg2 = fabs(par[8]*par[9]);


  double xJoin = log(NBkg/NBkg2)/(par[6]-par[9]);

  double value;
  if ((x[0]-mean)/sigma>-alpha)
    if (x[0]<xJoin) {
      value = binWidth*( N1*TMath::Gaus(x[0],mean,sigma,1) + 
			 N2*TMath::Gaus(x[0],mean*massRatio,sigma*massRatio,1) + 
			 NBkg*TMath::Exp(x[0]*par[6]) );
    }
    else {
      value = binWidth*( N1*TMath::Gaus(x[0],mean,sigma,1) + 
			 N2*TMath::Gaus(x[0],mean*massRatio,sigma*massRatio,1) + 
			 NBkg2*TMath::Exp(x[0]*par[9]) );
    }
  else {
    if (x[0]<xJoin) {
      value = binWidth*( N1/(sqrt(TMath::TwoPi())*sigma)*A*pow(B-(x[0]-mean)/sigma,-n) + 
			 N2/(sqrt(TMath::TwoPi())*sigma*massRatio)*A*pow(B-(x[0]-mean*massRatio)/(sigma*massRatio),-n) + 
			 NBkg*TMath::Exp(x[0]*par[6]) );
    }
    else {
      value = binWidth*( N1/(sqrt(TMath::TwoPi())*sigma)*A*pow(B-(x[0]-mean)/sigma,-n) + 
			 N2/(sqrt(TMath::TwoPi())*sigma*massRatio)*A*pow(B-(x[0]-mean*massRatio)/(sigma*massRatio),-n) + 
			 NBkg2*TMath::Exp(x[0]*par[9]) );
    }
  }

  return value;
}


double two_crystal_ball_gauss_expo(double *x, double *par)
{
  double N1 = par[0];

  double mean = par[1];
  double sigmaCB = par[2];
  double alpha = par[3];
  double n = par[4];

  double A = pow(n/fabs(alpha),n)*exp(-pow(alpha,2)/2);
  double B = n/fabs(alpha) - fabs(alpha);

  double massRatio = 3.686/3.0969;

  double N2;
  if (fitRatio)
    N2 = par[7]*N1;
  else
    N2 = par[7];

  double fracG = par[8];
  double sigmaG = par[9];

  double NBkg = par[5]*fabs(par[6]);

  if ((x[0]-mean)/sigmaCB>-alpha) {
    return binWidth*( N1*( (1-fracG)*TMath::Gaus(x[0],mean,sigmaCB,1) + 
			   fracG*TMath::Gaus(x[0],mean,sigmaG,1) ) +
		      N2*( (1-fracG)*TMath::Gaus(x[0],mean*massRatio,sigmaCB*massRatio,1) + 
			   fracG*TMath::Gaus(x[0],mean*massRatio,sigmaG*massRatio,1) ) + 
		      NBkg*TMath::Exp(x[0]*par[6]) );
  }
  else {
    return binWidth*( N1*( (1-fracG)/(sqrt(TMath::TwoPi())*sigmaCB)*A*pow(B-(x[0]-mean)/sigmaCB,-n) + 
			   fracG*TMath::Gaus(x[0],mean,sigmaG,1) ) + 
		      N2*( (1-fracG)/(sqrt(TMath::TwoPi())*sigmaCB)*A*pow(B-(x[0]-mean*massRatio)/(sigmaCB*massRatio),-n) + 
			   fracG*TMath::Gaus(x[0],mean*massRatio,sigmaG*massRatio,1) ) + 
		      NBkg*TMath::Exp(x[0]*par[6]) );
  }
}

double crystal_ball_pol2(double *x, double *par)
{
  double N = par[0];

  double mean = par[1];
  double sigma = par[2];
  double alpha = par[3];
  double n = par[4];

  double A = pow(n/fabs(alpha),n)*exp(-pow(alpha,2)/2);
  double B = n/fabs(alpha) - fabs(alpha);

  double NBkg = par[5];///(x[0]+0.5*x[0]*x[0]*par[6]+(1.0/3.0)*x[0]*x[0]*x[0]*par[7]);

  if ((x[0]-mean)/sigma>-alpha)
    return binWidth*(N*TMath::Gaus(x[0],mean,sigma,1)+NBkg*(1+x[0]*par[5]+x[0]*x[0]*par[6]));
  else
    return binWidth*(N/(sqrt(TMath::TwoPi())*sigma)*A*pow(B-(x[0]-mean)/sigma,-n)+NBkg*(1+x[0]*par[5]+x[0]*x[0]*par[6]));
}

double two_crystal_ball_pol2(double *x, double *par)
{
  double N1 = par[0];

  double mean = par[1];
  double sigma = par[2];
  double alpha = par[3];
  double n = par[4];

  double A = pow(n/fabs(alpha),n)*exp(-pow(alpha,2)/2);
  double B = n/fabs(alpha) - fabs(alpha);

  double massRatio = 3.686/3.0969;

  double N2;
  if (fitRatio)
    N2 = par[8]*N1;
  else
    N2 = par[8];

  double NBkg = par[5];///(x[0]+0.5*x[0]*x[0]*par[6]+(1.0/3.0)*x[0]*x[0]*x[0]*par[7]);

  if ((x[0]-mean)/sigma>-alpha)
    return binWidth*( N1*TMath::Gaus(x[0],mean,sigma,1) + 
		      N2*TMath::Gaus(x[0],mean*massRatio,sigma*massRatio,1) + 
		      NBkg*(1+x[0]*par[6]+x[0]*x[0]*par[7]) );
  else {
    return binWidth*( N1/(sqrt(TMath::TwoPi())*sigma)*A*pow(B-(x[0]-mean)/sigma,-n) + 
		      N2/(sqrt(TMath::TwoPi())*sigma*massRatio)*A*pow(B-(x[0]-mean*massRatio)/(sigma*massRatio),-n) + 
		      NBkg*(1+x[0]*par[6]+x[0]*x[0]*par[7]) );
  }
}


double two_crystal_ball_gauss_pol2(double *x, double *par)
{
  double N1 = par[0];

  double mean = par[1];
  double sigmaCB = par[2];
  double alpha = par[3];
  double n = par[4];

  double A = pow(n/fabs(alpha),n)*exp(-pow(alpha,2)/2);
  double B = n/fabs(alpha) - fabs(alpha);

  double massRatio = 3.686/3.0969;

  double N2;
  if (fitRatio)
    N2 = par[8]*N1;
  else
    N2 = par[8];

  double fracG = par[9];
  double sigmaG = par[10];

  double NBkg = par[5];///(x[0]+0.5*x[0]*x[0]*par[6]+(1.0/3.0)*x[0]*x[0]*x[0]*par[7]);

  if ((x[0]-mean)/sigmaCB>-alpha) {
    return binWidth*( N1*( (1-fracG)*TMath::Gaus(x[0],mean,sigmaCB,1) + 
			   fracG*TMath::Gaus(x[0],mean,sigmaG,1) ) +
		      N2*( (1-fracG)*TMath::Gaus(x[0],mean*massRatio,sigmaCB*massRatio,1) + 
			   fracG*TMath::Gaus(x[0],mean*massRatio,sigmaG*massRatio,1) ) + 
		      NBkg*(1+x[0]*par[6]+x[0]*x[0]*par[7]) );
  }
  else {
    return binWidth*( N1*( (1-fracG)/(sqrt(TMath::TwoPi())*sigmaCB)*A*pow(B-(x[0]-mean)/sigmaCB,-n) + 
			   fracG*TMath::Gaus(x[0],mean,sigmaG,1) ) + 
		      N2*( (1-fracG)/(sqrt(TMath::TwoPi())*sigmaCB)*A*pow(B-(x[0]-mean*massRatio)/(sigmaCB*massRatio),-n) + 
			   fracG*TMath::Gaus(x[0],mean*massRatio,sigmaG*massRatio,1) ) + 
		      NBkg*(1+x[0]*par[6]+x[0]*x[0]*par[7]) );
  }
}



void plot_Jpsi_psi2s_fit(std::string fname="../Jpsi_Histos_150ub_mass2-5.root", bool isHI=true, bool ratioFit=false, bool twoExpoFit=false, bool bgPol2Fit=false, bool bgPol1Fit=false)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);

  fitRatio = ratioFit;
  fit2Expo = twoExpoFit;
  fitBgPol2 = (!fit2Expo && bgPol2Fit) ? 1 : 0;
  fitBgPol1 = (bgPol1Fit && fitBgPol2) ? 1 : 0;

  f1 = new TF1("f1","0.020*[0]*TMath::Gaus(x,[1],[2],1)",2.6,3.5);
  f1->SetParNames("Yield","Mean","Sigma"); 
  f1->SetNpx(10000);
  f1->SetLineWidth(2);

  if (fitBgPol2) {
    f2 = new TF1("f2","0.020*([0]*TMath::Gaus(x,[1],[2],1)+[3]*(1+x*[4]+x*x*[5]))",2.6,3.5);
    f2->SetParNames("Yield","Mean","Sigma","N_{Bkg}","a","b");
  }
  else {
    f2 = new TF1("f2","0.020*([0]*TMath::Gaus(x,[1],[2],1)+[3]*abs([4])*TMath::Exp(x*[4]))",2.6,3.5);
    f2->SetParNames("Yield","Mean","Sigma","N_{Bkg}","a");
  }
  f2->SetNpx(10000);
  f2->SetLineWidth(2);
  f2->SetLineColor(4);
  if (fitBgPol2) {
    f4 = new TF1("f4","0.020*([0]*TMath::Gaus(x,[1],[2],1)+[3]*(1+x*[4]+x*x*[5]))",3.3,4.2);
    f4->SetParNames("Yield (#psi')","Mean","Sigma","N_{Bkg}","a","b");
  }
  else {
    f4 = new TF1("f4","0.020*([0]*TMath::Gaus(x,[1],[2],1)+[3]*abs([4])*TMath::Exp(x*[4]))",3.3,4.2);
    f4->SetParNames("Yield (#psi')","Mean","Sigma","N_{Bkg}","a");
  }
  f4->SetNpx(10000);
  f4->SetLineWidth(2);
  f4->SetLineColor(4);


  // f5 = new TF1("f5","0.020*([0]*TMath::Gaus(x,[1],[2],1)+[3]*TMath::Gaus(x,[4],[5],1))+TMath::Exp([6]+x*[7])",2.6,3.5);
  // f5->SetParNames("Yield (J/#psi)","Mean (J/#psi)","Sigma (J/#psi)","Yield (#psi')","Mean (#psi')","Sigma (#psi')","a","b");
  // f5->SetNpx(10000);
  // f5->SetLineWidth(2);
  f5 = new TF1("f5","0.020*([0]*TMath::Gaus(x,[1],[2],1)+[5]*[0]*TMath::Gaus(x,[1]*3.686/3.0969,[2]*3.686/3.0969,1)+[3]*TMath::Exp(x*[4]))",2.6,3.5);
  if (fitRatio)
    f5->SetParNames("Yield (J/#psi)","Mean (J/#psi)","Sigma (J/#psi)","N_{Bkg}","a","R_{#psi'}");
  else
    f5->SetParNames("Yield (J/#psi)","Mean (J/#psi)","Sigma (J/#psi)","N_{Bkg}","a","Yield (#psi')");
  f5->SetNpx(10000);
  f5->SetLineWidth(2);
  f5->SetLineColor(4);

  if (fitBgPol2) {
    f5a = new TF1("f5",two_crystal_ball_pol2,2.6,4.2,9);
    if (fitRatio)
      f5a->SetParNames("Yield (J/#psi)","Mean (J/#psi)","Sigma (J/#psi)","#alpha","n","N_{Bkg}","a","b","R_{#psi'}");
    else
      f5a->SetParNames("Yield (J/#psi)","Mean (J/#psi)","Sigma (J/#psi)","#alpha","n","N_{Bkg}","a","b","Yield (#psi')");
  }
  else if (fit2Expo){
    f5a = new TF1("f5",two_crystal_ball_two_expo,2.6,4.2,10);
    if (fitRatio)
      f5a->SetParNames("Yield (J/#psi)","Mean (J/#psi)","Sigma (J/#psi)","#alpha","n","N_{Bkg}","a","R_{#psi'}","N'_{Bkg}","a'");
    else
      f5a->SetParNames("Yield (J/#psi)","Mean (J/#psi)","Sigma (J/#psi)","#alpha","n","N_{Bkg}","a","Yield (#psi')","N'_{Bkg}","a'");
  }
  else {
    f5a = new TF1("f5",two_crystal_ball_expo,2.6,4.2,8);
    if (fitRatio)
      f5a->SetParNames("Yield (J/#psi)","Mean (J/#psi)","Sigma (J/#psi)","#alpha","n","N_{Bkg}","a","R_{#psi'}");
    else
      f5a->SetParNames("Yield (J/#psi)","Mean (J/#psi)","Sigma (J/#psi)","#alpha","n","N_{Bkg}","a","Yield (#psi')");
  }
  f5a->SetNpx(10000);
  f5a->SetLineWidth(2);
  f5a->SetLineColor(4);


  if (fitBgPol2) {
    f6 = new TF1("f6",crystal_ball_pol2,2.6,3.5,8);
    f6->SetParNames("Yield (J/#psi)","Mean","Sigma","#alpha","n","N_{Bkg}","a","b");
  }
  else {
    f6 = new TF1("f6",crystal_ball_expo,2.6,3.5,7);
    f6->SetParNames("Yield (J/#psi)","Mean","Sigma","#alpha","n","N_{Bkg}","a");
  }
  f6->SetNpx(10000);
  f6->SetLineWidth(2);
  f6->SetLineColor(4);

  if (fitBgPol2) {
    f7 = new TF1("f7",crystal_ball_pol2,3.3,4.2,8);
    f7->SetParNames("Yield (#psi')","Mean","Sigma","#alpha","n","N_{Bkg}","a","b");
  }
  else {
    f7 = new TF1("f7",crystal_ball_expo,3.3,4.2,7);
    f7->SetParNames("Yield (#psi')","Mean","Sigma","#alpha","n","N_{Bkg}","a");
  }
  f7->SetNpx(10000);
  f7->SetLineWidth(2);
  f7->SetLineColor(4);

  if (fitBgPol2) {
    f8 = new TF1("f8",two_crystal_ball_gauss_pol2,2.6,4.2,11);
    if (fitRatio)
      f8->SetParNames("Yield (J/#psi)","Mean (J/#psi)","Sigma_{CB} (J/#psi)","#alpha","n","N_{Bkg}","a","b","R_{#psi'}","f_{Gauss}","Sigma_{Gauss}");
    else
      f8->SetParNames("Yield (J/#psi)","Mean (J/#psi)","Sigma_{CB} (J/#psi)","#alpha","n","N_{Bkg}","a","b","Yield (#psi')","f_{Gauss}","Sigma_{Gauss}");
  }
  else {
    f8 = new TF1("f8",two_crystal_ball_gauss_expo,2.6,4.2,10);
    if (fitRatio)
      f8->SetParNames("Yield (J/#psi)","Mean (J/#psi)","Sigma_{CB} (J/#psi)","#alpha","n","N_{Bkg}","a","R_{#psi'}","f_{Gauss}","Sigma_{Gauss}");
    else
      f8->SetParNames("Yield (J/#psi)","Mean (J/#psi)","Sigma_{CB} (J/#psi)","#alpha","n","N_{Bkg}","a","Yield (#psi')","f_{Gauss}","Sigma_{Gauss}");
  }
  f8->SetNpx(10000);
  f8->SetLineWidth(2);
  f8->SetLineColor(4);

  if (fitBgPol2) {
    f3 = new TF1("f3","0.020*[0]*(1+[1]*x+[2]*x*x)",2.6,4.2);
    f3->SetParNames("N_{Bkg}","a","b");
  }
  else {
    if (fit2Expo) {
      f3 = new TF1("f3","0.020*(abs([0]*[1])*TMath::Exp(x*[1])*(x<log([0]/[2])/([1]-[3]))+abs([2]*[3])*TMath::Exp(x*[3])*(x>log([0]/[2])/([1]-[3])))",2.6,4.2);
      f3->SetParNames("N_{Bkg}","a","N'_{Bkg}","a'");
    }
    else {
      f3 = new TF1("f3","0.020*[0]*abs([1])*TMath::Exp(x*[1])",2.6,4.2);
      f3->SetParNames("N_{Bkg}","a");
    }
  }
  f3->SetNpx(10000);
  f3->SetLineWidth(2);
  f3->SetLineColor(4);
  
  TFile *inf;
  inf = new TFile(fname.c_str(),"READ");
  TCanvas *c1 = new TCanvas("c1","c1");
  if (isHI)
    c1->Divide(3,2);//4,2
  else
    c1->SetLogy();

  c1->cd(1);
  int nbins = 3.0/binWidth;
  TH1F *h1 = new TH1F("h1","Same sign pairs (0-100%);m_{inv} (#mu^{#pm}#mu^{#pm}) [GeV/c^{2}];counts",nbins,2.0,5.0);
  TH1F *h2 = new TH1F("h2","Opposite sign pairs (0-100%);m_{inv} (#mu^{+}#mu^{-}) [GeV/c^{2}];counts",nbins,2.0,5.0);

  TH1F *h3 = new TH1F("h3","h3;m_{inv} (#mu^{#pm}#mu^{#pm}) [GeV/c^{2}];counts",nbins,2.0,5.0);
  TH1F *h4 = new TH1F("h4","h4;m_{inv} (#mu^{+}#mu^{-}) [GeV/c^{2}];counts",nbins,2.0,5.0);
  TH1F *h5 = new TH1F("h5","h5;m_{inv} (#mu^{#pm}#mu^{#pm}) [GeV/c^{2}];counts",nbins,2.0,5.0);
  TH1F *h6 = new TH1F("h6","h6;m_{inv} (#mu^{+}#mu^{-}) [GeV/c^{2}];counts",nbins,2.0,5.0);
  TH1F *h7 = new TH1F("h7","h7;m_{inv} (#mu^{#pm}#mu^{#pm}) [GeV/c^{2}];counts",nbins,2.0,5.0);
  TH1F *h8 = new TH1F("h8","h8;m_{inv} (#mu^{+}#mu^{-}) [GeV/c^{2}];counts",nbins,2.0,5.0);

  TH1F *hPt = new TH1F("hPt","hPt;p_{T} (#mu^{+}#mu^{-}) [GeV/c];counts",300,0.0,30.0);
  TH1F *hPt2 = new TH1F("hPt2","hPt2;p_{T} (#mu^{#pm}#mu^{#pm}) [GeV/c];counts",300,0.0,30.0);

  TH1F *hRap = new TH1F("hRap","hRap;|y| (#mu^{+}#mu^{-});counts",24,0.0,2.4);
  TH1F *hRap2 = new TH1F("hRap2","hRap2;|y| (#mu^{#pm}#mu^{#pm});counts",24,0.0,2.4);

  h1->Sumw2();
  h1->SetMarkerSize(0.8);
  h1->SetMarkerStyle(24);
  h1->SetMarkerColor(kRed);
  h1->SetLineColor(kRed);
  h2->Sumw2();
  h2->SetMarkerSize(1.0);
  h2->SetMarkerStyle(20);
  h2->SetMarkerColor(kBlack);
  h2->SetLineColor(kBlack);

  h2->GetXaxis()->SetRangeUser(2.6,4.19);

  TH1F *h0010os = (TH1F*) h2->Clone("h0010os");
  TH1F *h1020os = (TH1F*) h2->Clone("h1020os");
  TH1F *h2030os = (TH1F*) h2->Clone("h2030os");
  TH1F *h3040os = (TH1F*) h2->Clone("h3040os");
  TH1F *h4050os = (TH1F*) h2->Clone("h4050os");
  TH1F *h5060os = (TH1F*) h2->Clone("h5060os");
  TH1F *h60100os = (TH1F*) h2->Clone("h60100os");
  TH1F *h0020os = (TH1F*) h2->Clone("h0020os");
  TH1F *h2040os = (TH1F*) h2->Clone("h2040os");
  TH1F *h40100os = (TH1F*) h2->Clone("h40100os");

  TH1F *h0010ss = (TH1F*) h1->Clone("h0010ss");
  TH1F *h1020ss = (TH1F*) h1->Clone("h1020ss");
  TH1F *h2030ss = (TH1F*) h1->Clone("h2030ss");
  TH1F *h3040ss = (TH1F*) h1->Clone("h3040ss");
  TH1F *h4050ss = (TH1F*) h1->Clone("h4050ss");
  TH1F *h5060ss = (TH1F*) h1->Clone("h5060ss");
  TH1F *h60100ss = (TH1F*) h1->Clone("h60100ss");
  TH1F *h0020ss = (TH1F*) h1->Clone("h0020ss");
  TH1F *h2040ss = (TH1F*) h1->Clone("h2040ss");
  TH1F *h40100ss = (TH1F*) h1->Clone("h40100ss");

  h0010os->SetTitle("Opposite sign pairs (0-10%)");
  h1020os->SetTitle("Opposite sign pairs (10-20%)");
  h2030os->SetTitle("Opposite sign pairs (20-30%)");
  h3040os->SetTitle("Opposite sign pairs (30-40%)");
  h4050os->SetTitle("Opposite sign pairs (40-50%)");
  h5060os->SetTitle("Opposite sign pairs (50-60%)");
  h60100os->SetTitle("Opposite sign pairs (60-100%)");
  h0020os->SetTitle("Opposite sign pairs (0-20%)");
  h2040os->SetTitle("Opposite sign pairs (20-40%)");
  h40100os->SetTitle("Opposite sign pairs (40-100%)");

  h0010ss->SetTitle("Same sign pairs (0-10%)");
  h1020ss->SetTitle("Same sign pairs (10-20%)");
  h2030ss->SetTitle("Same sign pairs (20-30%)");
  h3040ss->SetTitle("Same sign pairs (30-40%)");
  h4050ss->SetTitle("Same sign pairs (40-50%)");
  h5060ss->SetTitle("Same sign pairs (50-60%)");
  h60100ss->SetTitle("Same sign pairs (60-100%)");
  h0020ss->SetTitle("Same sign pairs (0-20%)");
  h2040ss->SetTitle("Same sign pairs (20-40%)");
  h40100ss->SetTitle("Same sign pairs (40-100%)");

  h3->Sumw2();
  h3->SetMarkerStyle(24);
  h3->SetMarkerColor(kRed);
  h3->SetLineColor(kRed);
  h4->Sumw2();
  h4->SetMarkerStyle(20);
  h4->SetMarkerColor(kBlack);
  h4->SetLineColor(kBlack);
  h5->Sumw2();
  h5->SetMarkerStyle(24);
  h5->SetMarkerColor(kRed);
  h5->SetLineColor(kRed);
  h6->Sumw2();
  h6->SetMarkerStyle(20);
  h6->SetMarkerColor(kBlack);
  h6->SetLineColor(kBlack);
  h7->Sumw2();
  h7->SetMarkerStyle(24);
  h7->SetMarkerColor(kRed);
  h7->SetLineColor(kRed);
  h8->Sumw2();
  h8->SetMarkerStyle(20);
  h8->SetMarkerColor(kBlack);
  h8->SetLineColor(kBlack);
  // h1->Rebin();
  // h2->Rebin();
  // h3->Rebin();
  hPt->Sumw2();
  hPt2->Sumw2();
  hRap->Sumw2();
  hRap2->Sumw2();

  // default
  // 38305
  ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>h2","Reco_QQ_sign==0&&Reco_QQ_type==0&&Reco_QQ_VtxProb>0.01&&Reco_QQ_4mom.Pt()>3.0&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4","e");
  ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>h1","Reco_QQ_sign!=0&&Reco_QQ_type==0&&Reco_QQ_VtxProb>0.01&&Reco_QQ_4mom.Pt()>3.0&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4","esame");

  // ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>h2","Reco_QQ_sign==0&&Reco_QQ_type==0&&Reco_QQ_VtxProb>0.01&&Reco_QQ_4mom.Pt()>6.5&&abs(Reco_QQ_4mom.Rapidity())<1.6","e");
  // ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>h1","Reco_QQ_sign!=0&&Reco_QQ_type==0&&Reco_QQ_VtxProb>0.01&&Reco_QQ_4mom.Pt()>6.5&&abs(Reco_QQ_4mom.Rapidity())<1.6","esame");

  h3->Add(h2,h1,1,-1);
  //  h2->Add(h1,-1);

  /*
  // 0-10%
  c1->cd(2);
  ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>h0010os","Reco_QQ_sign==0&&Reco_QQ_type==0&&Reco_QQ_VtxProb>0.05&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Centrality<4","e");
  ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>h0010ss","Reco_QQ_sign!=0&&Reco_QQ_type==0&&Reco_QQ_VtxProb>0.05&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Centrality<4","esame");
  //  h0010os->Add(h0010ss,-1);
  h0010ss->Rebin(5);
  h0010ss->Scale(0.2);
  // 10-20%
  c1->cd(3);
  ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>h1020os","Reco_QQ_sign==0&&Reco_QQ_type==0&&Reco_QQ_VtxProb>0.05&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Centrality>=4&&Centrality<8","e");
  ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>h1020ss","Reco_QQ_sign!=0&&Reco_QQ_type==0&&Reco_QQ_VtxProb>0.05&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Centrality>=4&&Centrality<8","esame");
  //  h1020os->Add(h1020ss,-1);
  h1020ss->Rebin(5);
  h1020ss->Scale(0.2);
  // 20-30%
  c1->cd(4);
  ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>h2030os","Reco_QQ_sign==0&&Reco_QQ_type==0&&Reco_QQ_VtxProb>0.05&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Centrality>=8&&Centrality<12","e");
  ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>h2030ss","Reco_QQ_sign!=0&&Reco_QQ_type==0&&Reco_QQ_VtxProb>0.05&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Centrality>=8&&Centrality<12","esame");
  //  h2030os->Add(h2030ss,-1);
  h2030ss->Rebin(5);
  h2030ss->Scale(0.2);
  // 30-40%
  c1->cd(4);
  ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>h3040os","Reco_QQ_sign==0&&Reco_QQ_type==0&&Reco_QQ_VtxProb>0.05&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Centrality>=12&&Centrality<16","e");
  ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>h3040ss","Reco_QQ_sign!=0&&Reco_QQ_type==0&&Reco_QQ_VtxProb>0.05&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Centrality>=12&&Centrality<16","esame");
  //  h3040os->Add(h3040ss,-1);
  h3040ss->Rebin(5);
  h3040ss->Scale(0.2);
  // 40-50%
  c1->cd(4);
  ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>h4050os","Reco_QQ_sign==0&&Reco_QQ_type==0&&Reco_QQ_VtxProb>0.05&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Centrality>=16&&Centrality<20","e");
  ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>h4050ss","Reco_QQ_sign!=0&&Reco_QQ_type==0&&Reco_QQ_VtxProb>0.05&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Centrality>=16&&Centrality<20","esame");  // 50-60%
  //  h4050os->Add(h4050ss,-1);
  h4050ss->Rebin(5);
  h4050ss->Scale(0.2);
  // 50-60%
  c1->cd(4);
  ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>h5060os","Reco_QQ_sign==0&&Reco_QQ_type==0&&Reco_QQ_VtxProb>0.05&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Centrality>=20&&Centrality<24","e");
  ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>h5060ss","Reco_QQ_sign!=0&&Reco_QQ_type==0&&Reco_QQ_VtxProb>0.05&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Centrality>=20&&Centrality<24","esame");  // 60-100%
  //  h5060os->Add(h5060ss,-1);
  h5060ss->Rebin(5);
  h5060ss->Scale(0.2);
  // 60-100%
  //  c1->cd(8);
  c1->cd(4);
  ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>h60100os","Reco_QQ_sign==0&&Reco_QQ_type==0&&Reco_QQ_VtxProb>0.05&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Centrality>=24","e");
  ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>h60100ss","Reco_QQ_sign!=0&&Reco_QQ_type==0&&Reco_QQ_VtxProb>0.05&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Centrality>=24","esame");
  //  h60100os->Add(h60100ss,-1);
  h60100ss->Rebin(5);
  h60100ss->Scale(0.2);
*/

  if (isHI) {
    // 0-20%
    c1->cd(2);
    ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>h0020os","Reco_QQ_sign==0&&Reco_QQ_type==0&&Reco_QQ_VtxProb>0.05&&Reco_QQ_4mom.Pt()>3.0&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Centrality<8","e");
    ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>h0020ss","Reco_QQ_sign!=0&&Reco_QQ_type==0&&Reco_QQ_VtxProb>0.05&&Reco_QQ_4mom.Pt()>3.0&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Centrality<8","esame");
    //  h0020os->Add(h0020ss,-1);
    h0020ss->Rebin(5);
    h0020ss->Scale(0.2);
    // 20-40%
    c1->cd(3);
    ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>h2040os","Reco_QQ_sign==0&&Reco_QQ_type==0&&Reco_QQ_VtxProb>0.05&&Reco_QQ_4mom.Pt()>3.0&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Centrality>=8&&Centrality<16","e");
    ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>h2040ss","Reco_QQ_sign!=0&&Reco_QQ_type==0&&Reco_QQ_VtxProb>0.05&&Reco_QQ_4mom.Pt()>3.0&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Centrality>=8&&Centrality<16","esame");
    //  h2040os->Add(h2040ss,-1);
    h2040ss->Rebin(5);
    h2040ss->Scale(0.2);
    c1->cd(4);
    // 40-100%
    ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>h40100os","Reco_QQ_sign==0&&Reco_QQ_type==0&&Reco_QQ_VtxProb>0.05&&Reco_QQ_4mom.Pt()>3.0&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Centrality>=16","e");
    ((TTree*) gROOT->FindObject("myTree"))->Draw("Reco_QQ_4mom.M()>>h40100ss","Reco_QQ_sign!=0&&Reco_QQ_type==0&&Reco_QQ_VtxProb>0.05&&Reco_QQ_4mom.Pt()>3.0&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4&&Centrality>=16","esame");
    //  h40100os->Add(h40100ss,-1);
    h40100ss->Rebin(5);
    h40100ss->Scale(0.2);
  }

  TCanvas *c2 = new TCanvas("c2","c2");
  if (isHI)
    c2->Divide(2,2);
  else
    c2->SetLogy();

  c2->cd(1);
  res_00100 = fitAndDrawHisto(h2,0);
  if (isHI) {
    c2->cd(2);
    res_0020 = fitAndDrawHisto(h0020os,0);
    //  res_0010 = fitAndDrawHisto(h0010os);
    c2->cd(3);
    //  res_1020 = fitAndDrawHisto(h1020os);
    //  c2->cd(4);
    //  res_2030 = fitAndDrawHisto(h2030os); // failure in last fit
    res_2040 = fitAndDrawHisto(h2040os,0);
    //  c2->cd(5);
    c2->cd(4);
    //  res_3040 = fitAndDrawHisto(h3040os); // failure in last fit
    res_40100 = fitAndDrawHisto(h40100os,0);
    // c2->cd(6);
    // res_4050 = fitAndDrawHisto(h4050os); // failure in last fit
    // c2->cd(7);
    // res_5060 = fitAndDrawHisto(h5060os); // failure at combined fit with CB
    // c2->cd(8);
    // res_60100 = fitAndDrawHisto(h60100os); // failure at Gaus fit
  }
  return;

  


  double signal = f2->Integral(f2->GetParameter(1)-2*f2->GetParameter(3), f2->GetParameter(1)+2*f2->GetParameter(3));
  double bg = f3->Integral(f2->GetParameter(1)-2*f2->GetParameter(3), f2->GetParameter(1)+2*f2->GetParameter(3))/0.050;

  std::cout << "S/B = " << signal << " / " << bg << " = " << signal/bg << std::endl;


  TCanvas *c3 = new TCanvas("c3","c3");
  c3->SetLogy();
  hPt->Add(hPt2,-1.0);
  hPt->Draw();

  // All+Barrel
  hPt->GetXaxis()->SetRangeUser(6.5,30.0);
  std::cout << "<pT> = "<< hPt->GetMean() << " +/- " << hPt->GetRMS() << " (RMS)" << std::endl;
  hPt->GetXaxis()->SetRangeUser(6.5,10.0);
  std::cout << "<pT> = "<< hPt->GetMean() << " +/- " << hPt->GetRMS() << " (RMS)" << std::endl;
  hPt->GetXaxis()->SetRangeUser(10.0,30.0);
  std::cout << "<pT> = "<< hPt->GetMean() << " +/- " << hPt->GetRMS() << " (RMS)" << std::endl;
  hPt->GetXaxis()->SetRangeUser(0.0,30.0);

  // Overlap
  // hPt->GetXaxis()->SetRangeUser(6.5,30.0);
  // std::cout << "<pT> = "<< hPt->GetMean() << " +/- " << hPt->GetRMS() << " (RMS)" << std::endl;
  // hPt->GetXaxis()->SetRangeUser(5.5,30.0);
  // std::cout << "<pT> = "<< hPt->GetMean() << " +/- " << hPt->GetRMS() << " (RMS)" << std::endl;

  // EndCap
  // hPt->GetXaxis()->SetRangeUser(6.5,30.0);
  // std::cout << "<pT> = "<< hPt->GetMean() << " +/- " << hPt->GetRMS() << " (RMS)" << std::endl;
  // hPt->GetXaxis()->SetRangeUser(3.0,30.0);
  // std::cout << "<pT> = "<< hPt->GetMean() << " +/- " << hPt->GetRMS() << " (RMS)" << std::endl;

  double error;
  double errorLS;

  std::cout << "Yields" << std::endl;
  hPt2->IntegralAndError(66,300,errorLS);
  hPt->IntegralAndError(66,300,error);
  std::cout << "6.5 < pT < 30.0: "<< hPt->IntegralAndError(66,300,error) << " +/- " << error << " +/- " << errorLS << std::endl;
  hPt2->IntegralAndError(66,100,errorLS);
  std::cout << "6.5 <= pT < 10.0: "<< hPt->IntegralAndError(66,100,error) << " +/- " << error << " +/- " << errorLS << std::endl;
  hPt2->IntegralAndError(101,300,errorLS);
  std::cout << "10.0 <= pT < 30.0: "<< hPt->IntegralAndError(101,300,error) << " +/- " << error << " +/- " << errorLS << std::endl;

  hPt2->IntegralAndError(1,65,errorLS);
  std::cout << "pT < 6.5: "<< hPt->IntegralAndError(1,65,error) << " +/- " << error << " +/- " << errorLS << std::endl;
  hPt2->IntegralAndError(31,300,errorLS);
  std::cout << "3.0 <= pT < 30.0: "<< hPt->IntegralAndError(31,300,error) << " +/- " << error << " +/- " << errorLS << std::endl;
  hPt2->IntegralAndError(56,300,errorLS);
  std::cout << "5.5 <= pT < 30.0: "<< hPt->IntegralAndError(56,300,error) << " +/- " << error << " +/- " << errorLS << std::endl;


  // cross check
  std::cout << "Yields" << std::endl;
  std::cout << "3-3.2: "<< h3->IntegralAndError(21,24,error) << " +/- " << error << std::endl;



  TCanvas *c4 = new TCanvas("c4","c4");
  c4->SetLogy();
  c4->cd();
  hRap->Add(hRap2,-1.0);
  hRap->Draw();


  hRap->GetXaxis()->SetRangeUser(0.0,2.4);
  std::cout << "<|y|> = "<< hRap->GetMean() << " +/- " << hRap->GetRMS() << " (RMS)" << std::endl;
  hRap->GetXaxis()->SetRangeUser(0.0,1.2);
  std::cout << "<|y|> = "<< hRap->GetMean() << " +/- " << hRap->GetRMS() << " (RMS)" << std::endl;
  hRap->GetXaxis()->SetRangeUser(1.2,1.6);
  std::cout << "<|y|> = "<< hRap->GetMean() << " +/- " << hRap->GetRMS() << " (RMS)" << std::endl;
  hRap->GetXaxis()->SetRangeUser(1.6,2.4);
  std::cout << "<|y|> = "<< hRap->GetMean() << " +/- " << hRap->GetRMS() << " (RMS)" << std::endl;
  hRap->GetXaxis()->SetRangeUser(0.0,2.4);


  return;
}

TFitResultPtr fitAndDrawHisto(TH1 *h, bool fitCBGsum)
{
  TFitResultPtr result;

  double epsilon = 0.0001;
  int massMin = h->GetXaxis()->FindBin(2.9+epsilon);
  int massMax = h->GetXaxis()->FindBin(3.2-epsilon);

  f1->SetParameters(h->Integral(massMin, massMax),3.0969,0.050);

  // fit of simple Gauss to J/psi
  h->Fit(f1,"QMEN","",2.9,3.2);
  f2->SetParameter(0, f1->GetParameter(0));
  f2->SetParameter(1, f1->GetParameter(1));
  f2->SetParameter(2, f1->GetParameter(2));

  if (fitBgPol2) {
    f2->SetParameter(3, 10000.0);
    f2->SetParameter(4, -0.15);
    if (fitBgPol1)
      f2->FixParameter(5, 0.0);
    else
      f2->SetParameter(5, 1.0);
  }
  else {
    f2->SetParameter(3, 1000.0);
    f2->SetParameter(4, -0.1);
  }
  // fit of Gauss + expo to J/psi
  h->Fit(f2,"QMEN","",2.6,3.5);

  massMin = h->GetXaxis()->FindBin(3.6+epsilon);
  massMax = h->GetXaxis()->FindBin(3.8-epsilon);
  
  f4->SetParameter(0, h->Integral(massMin, massMax));
  f4->SetParameter(1, 3.686);
  f4->SetParameter(2, f2->GetParameter(2)*3.686/3.0969);
  f4->SetParameter(3, f2->GetParameter(3));
  f4->SetParameter(4, f2->GetParameter(4));
  if (fitBgPol2) {
    if (fitBgPol1)
      f4->FixParameter(5, 0.0);
    else
      f4->SetParameter(5, f2->GetParameter(5));
  }

  // fit of Gaus + expo to psi(2S)
  h->Fit(f4,"QMEN","",3.4,4.2);

  f6->SetParameter(0, f2->GetParameter(0));
  f6->SetParameter(1, f2->GetParameter(1));
  f6->SetParameter(2, f2->GetParameter(2));
  f6->SetParameter(3, 1.5);
  f6->SetParameter(4, 2.0);
  f6->SetParameter(5, f2->GetParameter(3));
  f6->SetParameter(6, f2->GetParameter(4));
  if (fitBgPol2) {
    if (fitBgPol1)
      f6->FixParameter(7, 0.0);
    else
      f6->SetParameter(7, f2->GetParameter(5));
  }

  // fit of CB + expo to J/psi
  h->Fit(f6,"QRMEN","");

  f7->SetParameter(0, f4->GetParameter(0)); // yield
  f7->SetParameter(1, f4->GetParameter(1)); // mean
  f7->SetParameter(2, f4->GetParameter(2)); // sigma
  f7->SetParameter(3, f6->GetParameter(3)); // alpha
  f7->SetParameter(4, f6->GetParameter(4)); // n
  f7->SetParameter(5, f4->GetParameter(3)); // Nbkg
  f7->SetParameter(6, f4->GetParameter(4)); // a
  if (fitBgPol2) {
    if (fitBgPol1)
      f7->FixParameter(7, 0.0);
    else
      f7->SetParameter(7, f4->GetParameter(5));
  }

  // fit of CB + expo to psi(2S)
  h->Fit(f7,"QRMEN","sames");
    
  f5a->SetParameter(0, f6->GetParameter(0));
  f5a->SetParameter(1, f6->GetParameter(1));
  f5a->SetParameter(2, f6->GetParameter(2));
  f5a->SetParameter(3, f6->GetParameter(3));
  f5a->SetParameter(4, f6->GetParameter(4));
  f5a->SetParameter(5, f6->GetParameter(5));
  f5a->SetParameter(6, f6->GetParameter(6));
  if (fitBgPol2) {
    if (fitBgPol1)
      f5a->FixParameter(7, 0.0);
    else
      f5a->SetParameter(7, f6->GetParameter(7));
    if (fitRatio)
      f5a->SetParameter(8, f4->GetParameter(0)/f6->GetParameter(0));
    else
      f5a->SetParameter(8, f4->GetParameter(0));
  }
  else if (fit2Expo){
    if (fitRatio)
      f5a->SetParameter(7, f4->GetParameter(0)/f6->GetParameter(0));
    else
      f5a->SetParameter(7, f4->GetParameter(0));

    f5a->SetParameter(8, f6->GetParameter(5));
    f5a->SetParameter(9, f6->GetParameter(6));
  }
  else {
    if (fitRatio)
      f5a->SetParameter(7, f4->GetParameter(0)/f6->GetParameter(0));
    else
      f5a->SetParameter(7, f4->GetParameter(0));
  }

  // fit of CB + expo to J/psi and psi(2S)
  if (fitCBGsum)
    h->Fit(f5a,"QMEN","",2.6,4.2);
  else {
    result = h->Fit(f5a,"LMES","",2.6,4.2);
  }

  if (fitCBGsum) {
    f8->SetParameter(0, f5a->GetParameter(0)); // J/psi yield
    f8->SetParameter(1, f5a->GetParameter(1)); // J/psi mean
    f8->SetParameter(2, f5a->GetParameter(2)); // CB width
    f8->SetParameter(3, f5a->GetParameter(3)); // alpha
    f8->SetParameter(4, f5a->GetParameter(4)); // n
    f8->SetParameter(5, f5a->GetParameter(5)); // NBkg 
    f8->SetParameter(6, f5a->GetParameter(6)); // a

    if (fitBgPol2) {
      if (fitBgPol1)
	f8->FixParameter(7, 0.0);
      else
	f8->SetParameter(7, f5a->GetParameter(7)); // b
      f8->SetParameter(8, f5a->GetParameter(8)); // psi' yield
    }
    else {
      f8->SetParameter(7, f5a->GetParameter(7)); // psi' yield
    }
    //  if ( strcmp(h->GetName(),"h2") == 0) {
    f8->SetParameter(8, 0.1); // f_Gauss
    f8->SetParameter(9, 0.1); // Sigma_Gauss
    // }
    // else {
    //   f8->FixParameter(3, res_00100->Parameter(3)); // alpha
    //   f8->FixParameter(4, res_00100->Parameter(4)); // n
    //   f8->FixParameter(8, res_00100->Parameter(8)); // f_Gauss
    //   f8->FixParameter(9, res_00100->Parameter(9)); // Sigma_Gauss
    // }


    // fit of CB + expo to J/psi and psi(2S)
    result = h->Fit(f8,"MES","",2.6,4.2);

    // plot CB + bkg J/psi
    f6->SetParameter(0, (1.0-f8->GetParameter(8))*f8->GetParameter(0));
    f6->SetParameter(1, f8->GetParameter(1));
    f6->SetParameter(2, f8->GetParameter(2));
    f6->SetParameter(3, f8->GetParameter(3));
    f6->SetParameter(4, f8->GetParameter(4));    
    f6->SetParameter(5, f8->GetParameter(5));
    f6->SetParameter(6, f8->GetParameter(6));
    if (fitBgPol2) {
      if (fitBgPol1)
	f6->FixParameter(7, 0.0);
      else
	f6->SetParameter(7, f8->GetParameter(7));
    }

    // plot CB + bkg psi'
    if (fitBgPol2) {
      if (fitRatio)
	f7->SetParameter(0, (1.0-f8->GetParameter(9))*f8->GetParameter(8)*f8->GetParameter(0));
      else
	f7->SetParameter(0, (1.0-f8->GetParameter(9))*f8->GetParameter(8));
    }
    else {
      if (fitRatio)
	f7->SetParameter(0, (1.0-f8->GetParameter(8))*f8->GetParameter(7)*f8->GetParameter(0));
      else
	f7->SetParameter(0, (1.0-f8->GetParameter(8))*f8->GetParameter(7));
    }
    f7->SetParameter(1, f8->GetParameter(1)*3.686/3.0969);
    f7->SetParameter(2, f8->GetParameter(2)*3.686/3.0969);
    f7->SetParameter(3, f8->GetParameter(3));
    f7->SetParameter(4, f8->GetParameter(4));    
    f7->SetParameter(5, f8->GetParameter(5));
    f7->SetParameter(6, f8->GetParameter(6));
    if (fitBgPol2)  {
      if (fitBgPol1)
	f7->FixParameter(7, 0.0);
      else
	f7->SetParameter(7, f8->GetParameter(7));
    }
    // plot Gauss + bkg J/psi
    if (fitBgPol2)
      f2->SetParameter(0, f8->GetParameter(9)*f8->GetParameter(0));
    else 
      f2->SetParameter(0, f8->GetParameter(8)*f8->GetParameter(0));

    f2->SetParameter(1, f8->GetParameter(1));
    f2->SetParameter(2, f8->GetParameter(9));
    f2->SetParameter(3, f8->GetParameter(5));
    f2->SetParameter(4, f8->GetParameter(6));    

    // plot Gauss + bkg psi'
    if (fitBgPol2) {
      if (fitRatio)
	f4->SetParameter(0, f8->GetParameter(9)*f8->GetParameter(8)*f8->GetParameter(0));
      else
	f4->SetParameter(0, f8->GetParameter(8)*f8->GetParameter(7));
    }
    else {
      if (fitRatio)
	f4->SetParameter(0, f8->GetParameter(8)*f8->GetParameter(7)*f8->GetParameter(0));
      else
	f4->SetParameter(0, f8->GetParameter(8)*f8->GetParameter(7));
    }
    f4->SetParameter(1, f8->GetParameter(1)*3.686/3.0969);
    f4->SetParameter(2, f8->GetParameter(9)*3.686/3.0969);
    f4->SetParameter(3, f8->GetParameter(5));
    f4->SetParameter(4, f8->GetParameter(6));    
    if (fitBgPol2) {
      if (fitBgPol1)
	f4->FixParameter(5, 0.0);
      else
	f4->SetParameter(5, f8->GetParameter(7));
    }
    f3->SetParameter(0, f8->GetParameter(5));
    f3->SetParameter(1, f8->GetParameter(6));
    if (fitBgPol2) {
      if (fitBgPol1)
	f3->FixParameter(2, 0.0);
      else
	f3->SetParameter(2, f8->GetParameter(7));
    }
  }
  else {
    f6->SetParameter(0, f5a->GetParameter(0));
    f6->SetParameter(1, f5a->GetParameter(1));
    f6->SetParameter(2, f5a->GetParameter(2));
    f6->SetParameter(3, f5a->GetParameter(3));
    f6->SetParameter(4, f5a->GetParameter(4));    
    f6->SetParameter(5, f5a->GetParameter(5));
    f6->SetParameter(6, f5a->GetParameter(6));
    if (fitBgPol2) {
      if (fitBgPol1)
	f6->FixParameter(7, 0.0);
      else
	f6->SetParameter(7, f5a->GetParameter(7));
    }
    if (fitBgPol2) {
      if (fitRatio)
	f7->SetParameter(0, f8->GetParameter(8)*f8->GetParameter(0));
      else
	f7->SetParameter(0, f5a->GetParameter(8));
    }
    else {
      if (fitRatio)
	f7->SetParameter(0, f8->GetParameter(7)*f8->GetParameter(0));
      else
	f7->SetParameter(0, f5a->GetParameter(7));
    }
    f7->SetParameter(1, f5a->GetParameter(1)*3.686/3.0969);
    f7->SetParameter(2, f5a->GetParameter(2)*3.686/3.0969);
    f7->SetParameter(3, f5a->GetParameter(3));
    f7->SetParameter(4, f5a->GetParameter(4));    
    f7->SetParameter(5, f5a->GetParameter(5));
    f7->SetParameter(6, f5a->GetParameter(6));
    if (fitBgPol2) {
      if (fitBgPol1)
	f7->FixParameter(7, 0.0);
      else
	f7->SetParameter(7, f5a->GetParameter(7));
    }
    f3->SetParameter(0, f5a->GetParameter(5));
    f3->SetParameter(1, f5a->GetParameter(6));
    if (fitBgPol2) {
      if (fitBgPol1)
	f3->FixParameter(2, 0.0);
      else
	f3->SetParameter(2, f5a->GetParameter(7));
    }
    else if (fit2Expo) {
	f3->SetParameter(2, f5a->GetParameter(8));
	f3->SetParameter(3, f5a->GetParameter(9));
    }
  }

  f3->SetLineStyle(2);

  f6->SetLineColor(kGreen+2);
  f7->SetLineColor(kGreen+2);
  f2->SetLineColor(2);
  f4->SetLineColor(2);

  if (fitCBGsum) {
    f6->DrawCopy("same");
    f7->DrawCopy("same");
    f2->DrawCopy("same");
    f4->DrawCopy("same");
  }
  f3->DrawCopy("same");

  return result;
}
