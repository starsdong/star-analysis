#include "style.C+"
#include "draw.C+"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TMath.h"

//Double_t InvMass(Double_t p_P, Double_t p_Pi)

void test()
{
  style();
  //  const Double_t M_P = 0.93827;
  //  const Double_t M_Pi = 0.13957;

  const Int_t N = 10;
  double xp[N] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
  double yp[N] = {1.36, 1.73, 2.94, 4.15, 5.33, 5.40, 7.03, 7.91, 9.09, 10.39};
  double ype[N] = {0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3};

  TCanvas *c1 = new TCanvas("c1", "c1", 0, 0, 600, 600);
  c1->Draw();
  c1->cd();

  TH1D *h0 = new TH1D("h0", "", 1, 0, 1.2);
  h0->SetMinimum(0);
  h0->SetMaximum(12);
  h0->Draw();

  TGraphErrors *gr = new TGraphErrors(10, xp, yp, 0, ype);
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(1.5);
  gr->SetLineWidth(2);
  gr->Draw("p");

  gr->Fit("pol1","R");
  TF1 *polf = (TF1*)gr->GetListOfFunctions()->FindObject("pol1");
  double par0 = polf->GetParameter(0);
  double par1 = polf->GetParameter(1);
  double err0 = polf->GetParError(0);
  double err1 = polf->GetParError(1);


  TH2D *hab = new TH2D("hab","",200,5,15,200,-2.5,2.5);
  
  // Markov Chain Monte Carlo - Metropolis-Hastings Random Walk algorithm with Gibbs sampling method
  // function:  y = a*x + b
  TF1 *func = new TF1("func","[0]*x+[1]",0,12);
  TRandom3 *gRandom = new TRandom3();
  
  const Int_t NM = 2000000;
  double a_old = 10;
  double b_old = 0;
  for(int i=0;i<NM;i++) {
    double det = 1.;
    for(int j=0;j<N;j++) {
      det *= ype[j]*ype[j];
    }

    // try a new "a"
    double a = 10 + (gRandom->Rndm() - 0.5)*10.; // a range:  (5, 15);
    //    double a = gRandom->Gaus(a_old, 3.0);
    func->SetParameters(a_old, b_old);
    double logp_old = 0;
    for(int j=0;j<N;j++) {
      logp_old += -0.5*TMath::Power((func->Eval(xp[j]) - yp[j])/ype[j], 2.0);
    }
    //    logp_old += -0.5*TMath::Log(det);
    double prob_old = TMath::Exp(logp_old);

    func->SetParameters(a, b_old);
    double logp = 0;   // log probability
    for(int j=0;j<N;j++) {
      logp += -0.5*TMath::Power((func->Eval(xp[j]) - yp[j])/ype[j], 2.0);
    }
    //    logp += -0.5*TMath::Log(det);
    double prob = TMath::Exp(logp);

    if(gRandom->Rndm()<prob/prob_old) {
      a_old = a;
      hab->Fill(a, b_old);
    }

    // try a new "b"
    double b = 0 + (gRandom->Rndm() - 0.5)*5;     // b range:  (-2.5, 2.5);
    //    double b = gRandom->Gaus(b_old, 2.0);
    func->SetParameters(a_old, b_old);
    logp_old = 0.;
    for(int j=0;j<N;j++) {
      logp_old += -0.5*TMath::Power((func->Eval(xp[j]) - yp[j])/ype[j], 2.0);
    }
    prob_old = TMath::Exp(logp_old);

    func->SetParameters(a_old, b);
    logp = 0.;   // log probability
    for(int j=0;j<N;j++) {
      logp += -0.5*TMath::Power((func->Eval(xp[j]) - yp[j])/ype[j], 2.0);
    }
    //    logp += -0.5*TMath::Log(det);
    prob = TMath::Exp(logp);

    if(gRandom->Rndm()<prob/prob_old) {
      b_old = b;
      hab->Fill(a_old, b);
    }

  }

  TCanvas *c2 = new TCanvas("c2", "c2", 0, 0, 600, 600);
  c2->Draw();
  c2->cd();


  hab->Draw("colz");

  TH1D *ha = (TH1D *)hab->ProjectionX("ha");
  ha->Fit("gaus");
  TF1 *ga = (TF1*)ha->GetListOfFunctions()->FindObject("gaus");
  double a_fit = ga->GetParameter(1);
  double a_err_fit = ga->GetParameter(2);
  
  TH1D *hb = (TH1D *)hab->ProjectionY("hb");
  hb->Fit("gaus");
  TF1 *gb = (TF1*)hb->GetListOfFunctions()->FindObject("gaus");
  double b_fit = gb->GetParameter(1);
  double b_err_fit = gb->GetParameter(2);

  cout << " ========================================= " << endl;
  cout << "    Results from direct fit " << endl;
  cout << par0 << "+/-" << err0 << "\t" << par1 << "+/-" << err1 << endl;
  cout << " ========================================= " << endl;
  
  cout << " ========================================= " << endl;
  cout << "    Results from MCMC fit " << endl;
  cout << b_fit << "+/-" << b_err_fit << "\t" << a_fit << "+/-" << a_err_fit << endl;
  cout << " ========================================= " << endl;

  TFile *fout = new TFile("test.root","recreate");
  hab->Write();
  fout->Close();
  
  
}
