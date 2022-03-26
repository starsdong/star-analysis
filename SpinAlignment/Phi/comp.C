#include "style.C+"
#include "draw.C+"
#include "/Users/starsdong/work/work/C/plot6.C"
#include "/Users/starsdong/work/work/fitfun/rho00.C"

void comp()
{
  //  gROOT->ProcessLine(".L ~/work/work/C/plot6.C");

  style();
  TFile *f1 = new TFile("accept_1.root");
  TH2D *hPtCosThetaRc = (TH2D *)f1->Get("hPtCosThetaRc");
  TH1D *hCosThetaRc = (TH1D *)hPtCosThetaRc->ProjectionY("cosTheta_rc",51,60);
  hCosThetaRc->Rebin(5);
  TF1 *fRho00 = new TF1("rho00",CosThetaStar,-1,1,2);
  double par[2] = {hCosThetaRc->GetBinContent(1), 0.333};  
  plot(hCosThetaRc, -1., 1., hCosThetaRc->GetMaximum()*0.8, hCosThetaRc->GetMaximum()/0.8, "cos(#theta^{*})","Counts","test",0,fRho00,-1.,1.,par);

}
