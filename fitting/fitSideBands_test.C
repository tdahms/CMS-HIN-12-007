const int n=3;

double poly_func(double *x, double *par)
{
  double result;
  if ( (x[0] > 2.8 && x[0] < 3.3) ||
       (x[0] > 3.5 && x[0] < 3.85)) {
    TF1::RejectPoint();
    return 0;
  }

  switch (n) {
  case 0:
    result = par[0];
    break;
  case 1:
    result = par[0]*(1+par[1]*x[0]);
    break;
  case 2:
    result = par[0]*(1+par[1]*x[0]+par[2]*pow(x[0],2));
    break;
  case 3:
    result = par[0]*(1+par[1]*x[0]+par[2]*pow(x[0],2)+par[3]*pow(x[0],3));
    break;
  case 4:
    result = par[0]*(1+par[1]*x[0]+par[2]*pow(x[0],2)+par[3]*pow(x[0],3)+par[4]*pow(x[0],4));
    break;
  case 5:
    result = par[0]*(1+par[1]*x[0]+par[2]*pow(x[0],2)+par[3]*pow(x[0],3)+par[4]*pow(x[0],4)+par[5]*pow(x[0],5));
    break;
  default:
    result = par[0]*(1+par[1]*x[0]+par[2]*x[0]*x[0]);
  }

  return result;
}

void fitSideBands_test(bool isHI=true)
{
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->SetLogy();
  TFile *inf;
  if (isHI) {
    inf = new TFile("/Users/tdahms/CMS/HIN-12-007/Jpsi_Histos_PbPb_RegIT_glbglb.root","READ");
    myTree->Draw("Reco_QQ_4mom.M()>>h1(100,2.2,4.2)","(Reco_QQ_trig&1)==1&&Reco_QQ_sign==0&&Reco_QQ_4mom.Pt()>3.0&&Reco_QQ_4mom.Pt()<6.5&&abs(Reco_QQ_4mom.Rapidity())>1.6&&abs(Reco_QQ_4mom.Rapidity())<2.4","e");
  }
  else {
    inf = new TFile("/Users/tdahms/CMS/HIN-12-007/All_v2.24_Histos_Runs_211739-211831_GlbGlb_woPileUpRej_muLessPV.root","READ");
    myTree->Draw("Reco_QQ_4mom.M()>>h1(100,2.2,4.2)","(Reco_QQ_trig&2)==2&&Reco_QQ_sign==0&&Reco_QQ_4mom.Pt()>6.5&&Reco_QQ_4mom.Pt()<30.0&&abs(Reco_QQ_4mom.Rapidity())>0.0&&abs(Reco_QQ_4mom.Rapidity())<1.6","e");
  }


  TF1 *f1 = new TF1("f1",poly_func,2.2,4.2,n);
  f1->SetLineWidth(2);
  f1->SetLineColor(kRed);
  h1->Fit(f1,"IMER");

  h1->Draw();

  return;
}
