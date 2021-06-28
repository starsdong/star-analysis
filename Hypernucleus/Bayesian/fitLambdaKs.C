#include "style.C"
#include "draw.C"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TLine.h"

void generateDecay(double p, double M0, double M1, double M2, TLorentzVector *dau1, TLorentzVector *dau2)
{
  // const Double_t M1 = 0.938272;
  // const Double_t M2 = 0.139570;
  //  const Double_t M0 = 1.115683;
  //  cout << " ++ Input ++ " << p << " " << M0 << " " << M1 << " " << M2 << endl;
  const Double_t y0 = -1.05; // beam rapidity
  double E = TMath::Sqrt(p*p + M0*M0);
  double y, mT;
  do {
    y = gRandom->Rndm()*fabs(y0) + y0;  // first try - flat from (y0, 0);
    mT = E/TMath::CosH(y);
  } while (mT<M0||y<y0||y>0);  
  double pT = TMath::Sqrt(mT*mT - M0*M0);
  double phi = gRandom->Rndm()*TMath::Pi()*2.0;
  double pz = mT*TMath::SinH(y);
  double eta = 0.5*TMath::Log((p+pz)/(p-pz));

  TLorentzVector La(0.,0.,0.,0.);
  La.SetPtEtaPhiE(pT, eta, phi, E);
  //  cout << " ++ Parent (pT, eta, phi) " << pT << " " << eta << " " << phi << endl;

  // decay in parent rest frame
  double p_pri = TMath::Sqrt((M0+M1+M2)*(M0+M1-M2)*(M0-M1+M2)*(M0-M1-M2))/2./M0;
  double phi_pri = gRandom->Rndm()*TMath::Pi()*2.0;
  double costheta_pri = gRandom->Rndm()*2.-1;
  double pz_pri = p_pri*costheta_pri;
  double pT_pri = TMath::Sqrt(p_pri*p_pri-pz_pri*pz_pri);
  double eta_pri = 0.5*TMath::Log((p_pri+pz_pri)/(p_pri-pz_pri));
  
  dau1->SetPtEtaPhiM(pT_pri, eta_pri, phi_pri, M1);
  dau2->SetPtEtaPhiM(pT_pri, -eta_pri, phi_pri+TMath::Pi(), M2);
  
  //  transform to the lab frame
  dau1->Boost(La.BoostVector());
  dau2->Boost(La.BoostVector());

  //  cout << " ++ dau1 = " << dau1->X() << " " << dau1->Y() << " " << dau1->Z() << endl;
  //  cout << " ++ dau2 = " << dau2->X() << " " << dau2->Y() << " " << dau2->Z() << endl;

  return;
}

Bool_t isAccept(TLorentzVector dau1, TLorentzVector dau2)
{
  return dau1.Pt()>0.15 && dau2.Pt()>0.15 && dau1.Eta()>-1.8 && dau1.Eta()<0. && dau1.Eta()>-1.8 && dau2.Eta()<0.;
}

// pion correction function
TLorentzVector corrPion(TLorentzVector pi, TF1 *f_corr, Double_t p0, Double_t p1, Double_t p2)
{
  const Double_t M_Pi = 0.139570;
  double pi_corr_pt =  ( pi.Pt() - ( p0*f_corr->Eval(pi.P()/M_Pi) + p1 ) ) * p2;
  double s = pi_corr_pt / pi.Pt();
  TLorentzVector pion_corr(0.,0.,0.,0.);
  pion_corr.SetXYZM(pi.X()*s, pi.Y()*s, pi.Z()*s, M_Pi);
  return pion_corr;
}

// proton correction function
TLorentzVector corrProton(TLorentzVector pro, TF1 *f_corr, TGraph *gr_corr, Double_t p0, Double_t p1, Double_t p2, Double_t p3)
{
  const Double_t M_Pi = 0.139570;
  const Double_t M_P = 0.938272;
  double pro_corr_pt = ( pro.Pt() - ( p0*f_corr->Eval(pro.P()/M_Pi) + p1 + p3*gr_corr->Eval(pro.P()/M_P) ) ) * p2;
  double s = pro_corr_pt / pro.Pt();
  TLorentzVector proton_corr(0.,0.,0.,0.);
  proton_corr.SetXYZM(pro.X()*s, pro.Y()*s, pro.Z()*s, M_P);
  return proton_corr;
}

void fitLambdaKs(const Int_t NM=10000)
{
  gSystem->Load("style_C.so");
  gSystem->Load("draw_C.so");
  
  style();
  const Double_t M_P = 0.938272;
  const Double_t M_Pi = 0.139570;
  const Double_t M_La = 1.115683;
  const Double_t M_Ks = 0.497611;
  const Double_t s_sig = 10.; //10.;   // scale up the data uncertainties
  const Double_t p_binW = 0.25;
  
  TFile *fin = new TFile("dE_corr.root");
  TF1 *f_corr_p = (TF1 *)fin->Get("Corr_p");
  TF1 *f_corr_pi = (TF1 *)fin->Get("Corr_pi");
  TGraph *gr_p_dcorr = (TGraph *)fin->Get("dCorr_p");
  // correction functions - see dEFit.C for details
  // for pions: use f_corr_pi directly
  // for protons:   use f_corr_pi(p/m_pi) + gr_p_dcorr(p/m_p)
  fin->Close();

  fin = new TFile("Mass_Lamdba_Ks_Uncorr.root");
  TGraphErrors *gr_La = (TGraphErrors *)fin->Get("M_Uncorr_Lambda");
  gr_La->RemovePoint(0);
  //  gr_La->RemovePoint(0);
  //  gr_La->RemovePoint(0);
  TGraphErrors *gr_Ks = (TGraphErrors *)fin->Get("M_Uncorr_KS");
  gr_Ks->RemovePoint(0);
  //  gr_Ks->RemovePoint(0);
  //  gr_Ks->RemovePoint(0);
  const Int_t NSkip = 0; // skip the last X data points
  
  // Markov Chain Monte Carlo - Metropolis-Hastings Random Walk algorithm with Gibbs sampling method
  // function:  p_corr = ( p - [ a*f(p/m_pi) + b] ) * c    --- f(p/m_pi): energy loss correction function from embedding for pions
  // for pions:  p_corr = ( p - [ a*f(p/m_pi) + b] ) * c
  // for protons:  p_corr = ( p - [ a*f(p/m_pi) + b + d*gr_p_dcorr(p/m_p) ] ) * c  - version 20210608+, additional coefficient d need to be applied to gr_p_dcorr
  // a: default - 1, range - [0., 2.0]
  // b: default - 0, range - [-0.005, 0.005]
  // c: default - 0.98, range - [0.97, 1.03]
  // d: default - 1, range - [0., 2.0]
  const Int_t NP = 4;   // number of parameters
  const Int_t Nbin = 100;
  double par[NP] = {1.0, 0.0, 1.0, 1.0};  // initial values
  double par_p[NP];
  double err[NP] = {0.01, 0.0001, 0.0001, 0.01};  // dummy, not used
  double min[NP] = {-1.0, -0.002, 0.99, -1.0};
  double max[NP] = {2.0, 0.002, 1.01, 2.0};
  TH2D *h_corr[NP][NP]; // correlation matrices
  for(int i=0;i<NP;i++) {
    for(int j=0;j<NP;j++) {
      h_corr[i][j] = new TH2D(Form("Corr_%d_%d",i,j), "", Nbin, min[j], max[j], Nbin, min[i], max[i]);
    }
  }
  
  
  
  TRandom3 *gRandom = new TRandom3();

  for(int i=0;i<NM;i++) {
    if(i%100==0)
      cout << " Processing " << i << "-th event " << endl;

    for(int j=0;j<NP;j++) {
      ////////////////////////
      double logProb = 0;
      double logProb_p = 0;
      ////////////////////////
      // try a new par[j]
      ////////////////////////
      par_p[j] = gRandom->Rndm() * (max[j] - min[j]) + min[j];  // Random walk to a new position
      for(int k=0;k<gr_La->GetN()-NSkip;k++) {
	// Generate a Lambda decay
	TLorentzVector pion, proton;
	do {
	  generateDecay(gr_La->GetX()[k], gr_La->GetY()[k], M_P, M_Pi, &proton, &pion);
	  //	  generateDecay(gr_La->GetX()[k] + (gRandom->Rndm()-0.5)*p_binW, gr_La->GetY()[k], M_P, M_Pi, &proton, &pion);
	} while (!isAccept(pion, proton));
	
	// calculate the Prob at the current par[j] position
	TLorentzVector pion_corr = corrPion(pion, f_corr_pi, par[0], par[1], par[2]);
	TLorentzVector proton_corr = corrProton(proton, f_corr_pi, gr_p_dcorr, par[0], par[1], par[2], par[3]);
	TLorentzVector Lambda_corr = proton_corr + pion_corr;
	double La_Mass_corr = Lambda_corr.M();
	double sigma = gr_La->GetEY()[k]*s_sig;      
	logProb += -0.5*TMath::Power((La_Mass_corr - M_La)/sigma, 2.0);
	
	// calculate the Prob at the next par_p[j] position
	switch(j) {
	case 0:
	  pion_corr = corrPion(pion, f_corr_pi, par_p[0], par[1], par[2]);
	  proton_corr = corrProton(proton, f_corr_pi, gr_p_dcorr, par_p[0], par[1], par[2], par[3]);
	  break;
	case 1:
	  pion_corr = corrPion(pion, f_corr_pi, par[0], par_p[1], par[2]);
	  proton_corr = corrProton(proton, f_corr_pi, gr_p_dcorr, par[0], par_p[1], par[2], par[3]);
	  break;
	case 2:
	  pion_corr = corrPion(pion, f_corr_pi, par[0], par[1], par_p[2]);
	  proton_corr = corrProton(proton, f_corr_pi, gr_p_dcorr, par[0], par[1], par_p[2], par[3]);
	case 3:
	  pion_corr = corrPion(pion, f_corr_pi, par[0], par[1], par_p[2]);
	  proton_corr = corrProton(proton, f_corr_pi, gr_p_dcorr, par[0], par[1], par[2], par_p[3]);	  
	  break;
	default:
	  break;
	};
	
	Lambda_corr = proton_corr + pion_corr;
	La_Mass_corr = Lambda_corr.M();      
	logProb_p += -0.5*TMath::Power((La_Mass_corr - M_La)/sigma, 2.0);
      } // end loop Lambda
      
      for(int k=0;k<gr_Ks->GetN()-NSkip;k++) {
	// Generate a Ks decay
	TLorentzVector pip, pim;
	do {
	  generateDecay(gr_Ks->GetX()[k], gr_Ks->GetY()[k], M_Pi, M_Pi, &pip, &pim);
	  //	  generateDecay(gr_Ks->GetX()[k]+(gRandom->Rndm()-0.5)*p_binW, gr_Ks->GetY()[k], M_Pi, M_Pi, &pip, &pim);
	} while (!isAccept(pip, pim));
      
	// calculate the Prob at the current par[j] position	
	TLorentzVector pip_corr = corrPion(pip, f_corr_pi, par[0], par[1], par[2]);
	TLorentzVector pim_corr = corrPion(pim, f_corr_pi, par[0], par[1], par[2]);
	TLorentzVector Ks_corr = pip_corr + pim_corr;
	double Ks_Mass_corr = Ks_corr.M();
	double sigma = gr_Ks->GetEY()[k]*s_sig;      
	logProb += -0.5*TMath::Power((Ks_Mass_corr - M_Ks)/sigma, 2.0);

	// calculate the Prob at the next par_p[j] position
	switch(j) {
	case 0:
	  pip_corr = corrPion(pip, f_corr_pi, par_p[0], par[1], par[2]);
	  pim_corr = corrPion(pim, f_corr_pi, par_p[0], par[1], par[2]);
	  break;
	case 1:
	  pip_corr = corrPion(pip, f_corr_pi, par[0], par_p[1], par[2]);
	  pim_corr = corrPion(pim, f_corr_pi, par[0], par_p[1], par[2]);
	  break;
	case 2:
	  pip_corr = corrPion(pip, f_corr_pi, par[0], par[1], par_p[2]);
	  pim_corr = corrPion(pim, f_corr_pi, par[0], par[1], par_p[2]);
	  break;
	default:
	  break;
	};
	
	Ks_corr = pip_corr + pim_corr;
	Ks_Mass_corr = Ks_corr.M();      
	logProb_p += -0.5*TMath::Power((Ks_Mass_corr - M_Ks)/sigma, 2.0);
      } // end loop Ks
      
      if(gRandom->Rndm()<TMath::Exp(logProb_p - logProb)) {
	for(int l=0;l<NP;l++) {
	  h_corr[j][l]->Fill(par[l], par_p[j]);
	}
	par[j] = par_p[j];
      }
    } // end int j - number of parameters
  } // end int i - number of events

  TCanvas *c2 = new TCanvas("c2", "c2", 0, 0, 800, 800);
  c2->Draw();
  c2->Divide(NP, NP);
  for(int i=0;i<NP;i++) {
    for(int j=0;j<NP;j++) {
      c2->cd(i*NP+j+1);
      h_corr[i][j]->Draw("cont");
    }
  }
  c2->Update();
  c2->SaveAs("LambdaKsFitCorrMatrix.pdf");
  c2->SaveAs("LambdaKsFitCorrMatrix.png");

  double par_fit[NP], err_fit[NP];
  for(int i=0;i<NP;i++) {
    par_fit[i] = h_corr[i][i]->ProjectionY()->GetMean();
    err_fit[i] = h_corr[i][i]->ProjectionY()->GetRMS();
  }
  cout << " ++++++++++ Fit Results ++++++++++ " << endl;
  for (int i=0;i<NP;i++) {
    cout << " Par[" << i << "] = " << par_fit[i] << "+/-" << err_fit[i] << endl;
  }
  cout << " +++++++++++++++++++++++++++++++++ " << endl;
  
  double la_p[100], la_m[100], la_m_err[100], la_m_corr[100], la_m_corr_err[100];
  double ks_p[100], ks_m[100], ks_m_err[100], ks_m_corr[100], ks_m_corr_err[100];
  for(int k=0;k<gr_La->GetN();k++) {
    // Generate a Lambda decay
    TLorentzVector pion, proton;
    do {
      generateDecay(gr_La->GetX()[k], gr_La->GetY()[k], M_P, M_Pi, &proton, &pion);
      //      generateDecay(gr_La->GetX()[k] + (gRandom->Rndm()-0.5)*p_binW, gr_La->GetY()[k], M_P, M_Pi, &proton, &pion);
    } while (!isAccept(pion, proton));
    
    TLorentzVector pion_corr = corrPion(pion, f_corr_pi, par_fit[0], par_fit[1], par_fit[2]);
    TLorentzVector proton_corr = corrProton(proton, f_corr_pi, gr_p_dcorr, par_fit[0], par_fit[1], par_fit[2], par_fit[3]);
    TLorentzVector Lambda_corr = proton_corr + pion_corr;
    double La_Mass_corr = Lambda_corr.M();

    la_p[k] = gr_La->GetX()[k];
    la_m[k] = gr_La->GetY()[k];
    la_m_err[k] = gr_La->GetEY()[k]*s_sig;

    la_m_corr[k] = Lambda_corr.M();
    la_m_corr_err[k] = gr_La->GetEY()[k]*s_sig;
  }  

  for(int k=0;k<gr_Ks->GetN();k++) {
    // Generate a Ks decay
    TLorentzVector pip, pim;
    do {
      generateDecay(gr_Ks->GetX()[k], gr_Ks->GetY()[k], M_Pi, M_Pi, &pip, &pim);
      //      generateDecay(gr_Ks->GetX()[k] + (gRandom->Rndm()-0.5)*p_binW, gr_Ks->GetY()[k], M_Pi, M_Pi, &pip, &pim);
    } while (!isAccept(pip, pim));
    
    TLorentzVector pip_corr = corrPion(pip, f_corr_pi, par_fit[0], par_fit[1], par_fit[2]);
    TLorentzVector pim_corr = corrPion(pim, f_corr_pi, par_fit[0], par_fit[1], par_fit[2]);
    TLorentzVector Ks_corr = pip_corr + pim_corr;
    double Ks_Mass_corr = Ks_corr.M();

    ks_p[k] = gr_Ks->GetX()[k];
    ks_m[k] = gr_Ks->GetY()[k];
    ks_m_err[k] = gr_Ks->GetEY()[k]*s_sig;

    ks_m_corr[k] = Ks_corr.M();
    ks_m_corr_err[k] = gr_Ks->GetEY()[k]*s_sig;    
  }  
  TGraphErrors *gr_la_m = new TGraphErrors(gr_La->GetN(), la_p, la_m, 0, la_m_err);
  gr_la_m->Print();
  TGraphErrors *gr_la_m_corr = new TGraphErrors(gr_La->GetN(), la_p, la_m_corr, 0, la_m_corr_err);
  gr_la_m_corr->Print();
  
  TGraphErrors *gr_ks_m = new TGraphErrors(gr_Ks->GetN(), ks_p, ks_m, 0, ks_m_err);
  gr_ks_m->Print();
  TGraphErrors *gr_ks_m_corr = new TGraphErrors(gr_Ks->GetN(), ks_p, ks_m_corr, 0, ks_m_corr_err);
  gr_ks_m_corr->Print();
  
  TCanvas *c1 = new TCanvas("c1", "c1", 600, 600, 600, 800);
  c1->Divide(1,2);
  c1->Draw();

  c1->cd(1);
  TH1D *h0 = new TH1D("h0", "", 1, 0, 5.0);
  h0->SetMinimum(1.115);
  h0->SetMaximum(1.117);
  h0->Draw("c");

  TLine *l0 = new TLine(0., M_La, 5.0, M_La);
  l0->SetLineWidth(2);
  l0->SetLineStyle(4);
  l0->Draw("same");

  gr_la_m->SetMarkerStyle(24);
  gr_la_m->SetMarkerSize(1.5);
  gr_la_m->SetLineWidth(2);
  gr_la_m->Draw("p");

  gr_la_m_corr->SetMarkerStyle(20);
  gr_la_m_corr->SetMarkerSize(1.5);
  gr_la_m_corr->SetLineWidth(2);
  gr_la_m_corr->Draw("p");

  //
  c1->cd(2);

  TH1D *h1 = new TH1D("h1", "", 1, 0, 5.0);
  h1->SetMinimum(0.495);
  h1->SetMaximum(0.500);
  h1->Draw("c");
  
  TLine *l1 = new TLine(0., M_Ks, 5.0, M_Ks);
  l1->SetLineWidth(2);
  l1->SetLineStyle(4);
  l1->Draw("same");

  gr_ks_m->SetMarkerStyle(24);
  gr_ks_m->SetMarkerSize(1.5);
  gr_ks_m->SetLineWidth(2);
  gr_ks_m->Draw("p");

  gr_ks_m_corr->SetMarkerStyle(20);
  gr_ks_m_corr->SetMarkerSize(1.5);
  gr_ks_m_corr->SetLineWidth(2);
  gr_ks_m_corr->Draw("p");
  
  c1->Update();
  c1->SaveAs("LambdaKsFit.pdf");
  c1->SaveAs("LambdaKsFit.png");
  
  
  
}
