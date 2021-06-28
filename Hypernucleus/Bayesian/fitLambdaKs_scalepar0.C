//#include "style.C+"
//#include "draw.C+"
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

void generateLambdaDecay(double p, double M0, TLorentzVector *dau1, TLorentzVector *dau2)
{
  const Double_t M1 = 0.938272;
  const Double_t M2 = 0.139570;
  //  const Double_t M0 = 1.115683;
  const Double_t y0 = -1.05; // beam rapidity
  double E = TMath::Sqrt(p*p + M0*M0);
  double y, mT;
  do {
    y = gRandom->Rndm()*fabs(y0) + y0;  // first try - flat from (y0, 0);
    mT = E/TMath::CosH(y);
  } while (mT<M0);  
  double pT = TMath::Sqrt(mT*mT - M0*M0);
  double phi = gRandom->Rndm()*TMath::Pi()*2.0;
  double pz = mT*TMath::SinH(y);
  double eta = 0.5*TMath::Log((p+pz)/(p-pz));

  TLorentzVector La(0.,0.,0.,0.);
  La.SetPtEtaPhiE(pT, eta, phi, E);
  //  cout << " ++ Lambda (pT, eta, phi) " << pT << " " << eta << " " << phi << endl;

  // decay in Lambda rest frame
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

  // cout << " ++ dau1 = " << dau1->X() << " " << dau1->Y() << " " << dau1->Z() << endl;
  // cout << " ++ dau2 = " << dau2->X() << " " << dau2->Y() << " " << dau2->Z() << endl;

  return;
}

void fitLambdaKs()
{
  // gSystem->Load("style_C.so");
  // gSystem->Load("draw_C.so");
  
  //  style();
  const Double_t M_P = 0.938272;
  const Double_t M_Pi = 0.139570;
  const Double_t M_La = 1.115683;
  const Double_t M_Ks = 0.497611;
  const Double_t s_sig = 10.;
  
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
  TGraphErrors *gr_Ks = (TGraphErrors *)fin->Get("M_Uncorr_KS");
  
  // Markov Chain Monte Carlo - Metropolis-Hastings Random Walk algorithm with Gibbs sampling method
  // function:  p_corr = ( p - [ a*f(p/m_pi) + b] ) * c    --- f(p/m_pi): energy loss correction function from embedding for pions
  // for pions:  p_corr = ( p - [ a*f(p/m_pi) + b] ) * c
  // for protons:  p_corr = ( p - [ a*f(p/m_pi) + b + gr_p_dcorr(p/m_p) ] ) * c
  // a: default - 1, range - [0.5, 1.5]
  // b: default - 0, range - [-0.005, 0.005]
  // c: default - 0.98, range - [0.94, 1.02]
  const Int_t NP = 3;
  const Int_t Nbin = 100;
  double par[NP] = {1.0, 0.0, 1.0};  // initial values
  double par_p[NP];
  double err[NP] = {0.01, 0.01, 0.001};  // dummy, not used
  double min[NP] = {0.1, -0.005, 0.97};
  double max[NP] = {2.0, 0.005, 1.03};
  TH2D *h_corr[NP*(NP+1)/2]; // correlation matrices
  int index = 0;
  for(int i=0;i<NP;i++) {
    for(int j=0;j<=i;j++) {
      h_corr[index++] = new TH2D(Form("Corr_%d_%d",i,j), "", Nbin, min[j], max[j], Nbin, min[i], max[i]);
      // 0:  par[0]:par[0]
      // 1:  par[1]:par[0]    2:  par[1]:par[1]
      // 3:  par[2]:par[0]    4:  par[2]:par[1]    5:  par[2]:par[2] 
    }
  }
  
  
  
  const Int_t NM = 100000;
  TRandom3 *gRandom = new TRandom3();

  
  for(int i=0;i<NM;i++) {
    if(i%1000==0) cout << " Processing " << i << "-th event " << endl;

    ////////////////////////
    double logProb = 0;
    double logProb_p = 0;
    ////////////////////////
    // try a new par[0]
    ////////////////////////
    //    cout << " Test ==== par[0-2] = " << par[0] << "\t" << par[1] << "\t" << par[2] << endl;
    par_p[0] = gRandom->Rndm() * (max[0] - min[0]) + min[0];  // Random walk to a new position
    //    cout << " Try  ==== par[0-2] = " << par_p[0] << "\t" << par[1] << "\t" << par[2] << endl;
    for(int k=0;k<gr_La->GetN();k++) {
      // Generate a Lambda decay
      TLorentzVector daughters[2];
      do {
	generateLambdaDecay(gr_La->GetX()[k], gr_La->GetY()[k], &daughters[0], &daughters[1]);
	TLorentzVector pion = daughters[1];
	TLorentzVector proton = daughters[0];
	//	cout << " pion " << pion.X() << " " << pion.Y() << " " << pion.Z() << " " << pion.E() << "\t proton " << proton.X() << " " << proton.Y() << " " << proton.Z() << " " << proton.E() << endl;
      } while (daughters[0].Pt()<0.1 || daughters[1].Pt()<0.1 || daughters[0].Eta()<-2 || daughters[0].Eta()>0 || daughters[1].Eta()<-2 || daughters[1].Eta()>0);

      TLorentzVector pion = daughters[1];
      TLorentzVector proton = daughters[0];
      //      cout << " pion " << pion.X() << " " << pion.Y() << " " << pion.Z() << " " << pion.E() << "\t proton " << proton.X() << " " << proton.Y() << " " << proton.Z() << " " << proton.E() << endl;

      // calculate the Prob at the current par[0] position
      double dp_pion = par[0]*f_corr_pi->Eval(pion.P()/M_Pi) + par[1];
      double pion_p_corr = ( pion.P() - ( f_corr_pi->Eval(pion.P()/M_Pi*par[0]) + par[1] ) ) * par[2];
      double s_pion = pion_p_corr/pion.P();
      TLorentzVector pion_corr(0.,0.,0.,0.);
      pion_corr.SetXYZM(pion.X()*s_pion, pion.Y()*s_pion, pion.Z()*s_pion, M_Pi);
      
      double proton_p_corr = ( proton.P() - ( f_corr_pi->Eval(proton.P()/M_Pi*par[0]) + par[1] + gr_p_dcorr->Eval(proton.P()/M_P*par[0]) ) ) * par[2];
      double s_proton = proton_p_corr/proton.P();
      TLorentzVector proton_corr(0.,0.,0.,0.);
      proton_corr.SetXYZM(proton.X()*s_proton, proton.Y()*s_proton, proton.Z()*s_proton, M_P);
      
      TLorentzVector Lambda_corr = proton_corr + pion_corr;
      double La_Mass_corr = Lambda_corr.M();
      // cout << " At default position: momentum " << pion.P() << " " << dp_pion << endl;
      // cout << " At default position: scaling factor " << s_pion << " " << s_proton << endl;
      //      cout << " Lambda mom = " << gr_La->GetX()[k] << " Mass = " << gr_La->GetY()[k] << " Mass_Corr = " << La_Mass_corr << endl;
      double sigma = gr_La->GetEY()[k]*s_sig;
      
      logProb += -0.5*TMath::Power((La_Mass_corr - M_La)/sigma, 2.0);

      // calculate the Prob at the next par_p[0] position
      dp_pion = par_p[0]*f_corr_pi->Eval(pion.P()/M_Pi) + par[1];
      pion_p_corr = ( pion.P() - ( f_corr_pi->Eval(pion.P()/M_Pi*par_p[0]) + par[1] ) ) * par[2];
      s_pion = pion_p_corr/pion.P();
      pion_corr.SetXYZM(pion.X()*s_pion, pion.Y()*s_pion, pion.Z()*s_pion, M_Pi);
      
      proton_p_corr = ( proton.P() - ( f_corr_pi->Eval(proton.P()/M_Pi*par_p[0]) + par[1] + gr_p_dcorr->Eval(proton.P()/M_P*par[0]) ) ) * par[2];
      s_proton = proton_p_corr/proton.P();
      proton_corr.SetXYZM(proton.X()*s_proton, proton.Y()*s_proton, proton.Z()*s_proton, M_P);
      
      Lambda_corr = proton_corr + pion_corr;
      La_Mass_corr = Lambda_corr.M();
      // cout << " At try position: momentum " << pion.P() << " " << dp_pion << endl;
      // cout << " At try position: scaling factor " << s_pion << " " << s_proton << endl;
      //      cout << " Lambda mom = " << gr_La->GetX()[k] << " Mass = " << gr_La->GetY()[k] << " Mass_Corr_Try = " << La_Mass_corr << endl;
      
      logProb_p += -0.5*TMath::Power((La_Mass_corr - M_La)/sigma, 2.0);
    }
    //    cout << "LogProb = " << logProb << "\t" << logProb_p << endl;
    double prob = TMath::Exp(logProb);
    double prob_p = TMath::Exp(logProb_p);
    if(gRandom->Rndm()<prob_p/prob) {
      h_corr[0]->Fill(par[0], par_p[0]);
      par[0] = par_p[0];
      h_corr[1]->Fill(par[0], par[1]);
      h_corr[3]->Fill(par[0], par[2]);
      //      cout << " Filling .... " << endl;
    }

    ////////////////////////
    logProb = 0.;
    logProb_p = 0.;
    ////////////////////////
    // try a new par[1]
    ////////////////////////
    //    cout << " Test ==== par[0-2] = " << par[0] << "\t" << par[1] << "\t" << par[2] << endl;
    par_p[1] = gRandom->Rndm() * (max[1] - min[1]) + min[1];  // Random walk to a new position
    //    cout << " Try  ==== par[0-2] = " << par[0] << "\t" << par_p[1] << "\t" << par[2] << endl;
    for(int k=0;k<gr_La->GetN();k++) {
      // Generate a Lambda decay
      TLorentzVector daughters[2];
      do {
	generateLambdaDecay(gr_La->GetX()[k], gr_La->GetY()[k], &daughters[0], &daughters[1]);
	TLorentzVector pion = daughters[1];
	TLorentzVector proton = daughters[0];
	//	cout << " pion " << pion.X() << " " << pion.Y() << " " << pion.Z() << " " << pion.E() << "\t proton " << proton.X() << " " << proton.Y() << " " << proton.Z() << " " << proton.E() << endl;
      } while (daughters[0].Pt()<0.1 || daughters[1].Pt()<0.1 || daughters[0].Eta()<-2 || daughters[0].Eta()>0 || daughters[1].Eta()<-2 || daughters[1].Eta()>0);

      TLorentzVector pion = daughters[1];
      TLorentzVector proton = daughters[0];

      // calculate the Prob at the current par[0] position
      double pion_p_corr = ( pion.P() - ( f_corr_pi->Eval(pion.P()/M_Pi*par[0]) + par[1] ) ) * par[2];
      double s_pion = pion_p_corr/pion.P();
      TLorentzVector pion_corr(0.,0.,0.,0.);
      pion_corr.SetXYZM(pion.X()*s_pion, pion.Y()*s_pion, pion.Z()*s_pion, M_Pi);
      
      double proton_p_corr = ( proton.P() - ( f_corr_pi->Eval(proton.P()/M_Pi*par[0]) + par[1] + gr_p_dcorr->Eval(proton.P()/M_P*par[0]) ) ) * par[2];
      double s_proton = proton_p_corr/proton.P();
      TLorentzVector proton_corr(0.,0.,0.,0.);
      proton_corr.SetXYZM(proton.X()*s_proton, proton.Y()*s_proton, proton.Z()*s_proton, M_P);
      
      TLorentzVector Lambda_corr = proton_corr + pion_corr;
      double La_Mass_corr = Lambda_corr.M();
      double sigma = gr_La->GetEY()[k]*s_sig;
      
      logProb += -0.5*TMath::Power((La_Mass_corr - M_La)/sigma, 2.0);

      // calculate the Prob at the next par_p[1] position
      pion_p_corr = ( pion.P() - ( f_corr_pi->Eval(pion.P()/M_Pi*par[0]) + par_p[1] ) ) * par[2];
      s_pion = pion_p_corr/pion.P();
      pion_corr.SetXYZM(pion.X()*s_pion, pion.Y()*s_pion, pion.Z()*s_pion, M_Pi);
      
      proton_p_corr = ( proton.P() - ( f_corr_pi->Eval(proton.P()/M_Pi*par[0]) + par_p[1] + gr_p_dcorr->Eval(proton.P()/M_P*par[0]) ) ) * par[2];
      s_proton = proton_p_corr/proton.P();
      proton_corr.SetXYZM(proton.X()*s_proton, proton.Y()*s_proton, proton.Z()*s_proton, M_P);
      
      Lambda_corr = proton_corr + pion_corr;
      La_Mass_corr = Lambda_corr.M();
      
      logProb_p += -0.5*TMath::Power((La_Mass_corr - M_La)/sigma, 2.0);
    }
    //    cout << "LogProb = " << logProb << "\t" << logProb_p << endl;
    prob = TMath::Exp(logProb);
    prob_p = TMath::Exp(logProb_p);
    if(gRandom->Rndm()<prob_p/prob) {
      h_corr[2]->Fill(par[1], par_p[1]);
      par[1] = par_p[1];
      h_corr[1]->Fill(par[0], par[1]);
      h_corr[4]->Fill(par[1], par[2]);
    }
    
    ////////////////////////
    logProb = 0.;
    logProb_p = 0.;
    ////////////////////////
    // try a new par[2]
    ////////////////////////
    //    cout << " Test ==== par[0-2] = " << par[0] << "\t" << par[1] << "\t" << par[2] << endl;
    par_p[2] = gRandom->Rndm() * (max[2] - min[2]) + min[2];  // Random walk to a new position
    //    cout << " Try  ==== par[0-2] = " << par[0] << "\t" << par[1] << "\t" << par_p[2] << endl;
    for(int k=0;k<gr_La->GetN();k++) {
      // Generate a Lambda decay
      TLorentzVector daughters[2];
      do {
	generateLambdaDecay(gr_La->GetX()[k], gr_La->GetY()[k], &daughters[0], &daughters[1]);
	TLorentzVector pion = daughters[1];
	TLorentzVector proton = daughters[0];
	//	cout << " pion " << pion.X() << " " << pion.Y() << " " << pion.Z() << " " << pion.E() << "\t proton " << proton.X() << " " << proton.Y() << " " << proton.Z() << " " << proton.E() << endl;
      } while (daughters[0].Pt()<0.1 || daughters[1].Pt()<0.1 || daughters[0].Eta()<-2 || daughters[0].Eta()>0 || daughters[1].Eta()<-2 || daughters[1].Eta()>0);

      TLorentzVector pion = daughters[1];
      TLorentzVector proton = daughters[0];

      // calculate the Prob at the current par[0] position
      double pion_p_corr = ( pion.P() - ( f_corr_pi->Eval(pion.P()/M_Pi*par[0]) + par[1] ) ) * par[2];
      double s_pion = pion_p_corr/pion.P();
      TLorentzVector pion_corr(0.,0.,0.,0.);
      pion_corr.SetXYZM(pion.X()*s_pion, pion.Y()*s_pion, pion.Z()*s_pion, M_Pi);
      
      double proton_p_corr = ( proton.P() - ( f_corr_pi->Eval(proton.P()/M_Pi*par[0]) + par[1] + gr_p_dcorr->Eval(proton.P()/M_P*par[0]) ) ) * par[2];
      double s_proton = proton_p_corr/proton.P();
      TLorentzVector proton_corr(0.,0.,0.,0.);
      proton_corr.SetXYZM(proton.X()*s_proton, proton.Y()*s_proton, proton.Z()*s_proton, M_P);
      
      TLorentzVector Lambda_corr = proton_corr + pion_corr;
      double La_Mass_corr = Lambda_corr.M();
      double sigma = gr_La->GetEY()[k]*s_sig;
      
      logProb += -0.5*TMath::Power((La_Mass_corr - M_La)/sigma, 2.0);

      // calculate the Prob at the next par_p[2] position
      pion_p_corr = ( pion.P() - ( f_corr_pi->Eval(pion.P()/M_Pi*par[0]) + par[1] ) ) * par_p[2];
      s_pion = pion_p_corr/pion.P();
      pion_corr.SetXYZM(pion.X()*s_pion, pion.Y()*s_pion, pion.Z()*s_pion, M_Pi);
      
      proton_p_corr = ( proton.P() - ( f_corr_pi->Eval(proton.P()/M_Pi*par[0]) + par[1] + gr_p_dcorr->Eval(proton.P()/M_P*par[0]) ) ) * par_p[2];
      s_proton = proton_p_corr/proton.P();
      proton_corr.SetXYZM(proton.X()*s_proton, proton.Y()*s_proton, proton.Z()*s_proton, M_P);
      
      Lambda_corr = proton_corr + pion_corr;
      La_Mass_corr = Lambda_corr.M();
      
      logProb_p += -0.5*TMath::Power((La_Mass_corr - M_La)/sigma, 2.0);
    }
    //    cout << "LogProb = " << logProb << "\t" << logProb_p << endl;
    prob = TMath::Exp(logProb);
    prob_p = TMath::Exp(logProb_p);
    if(gRandom->Rndm()<prob_p/prob) {
      h_corr[5]->Fill(par[2], par_p[2]);
      par[2] = par_p[2];
      h_corr[3]->Fill(par[0], par[2]);
      h_corr[4]->Fill(par[1], par[2]);
    }

  }

  TCanvas *c2 = new TCanvas("c2", "c2", 0, 0, 800, 800);
  c2->Draw();
  c2->Divide(3,3);
  c2->cd(1);
  h_corr[0]->Draw("colz");
  c2->cd(4);
  h_corr[1]->Draw("colz");
  c2->cd(5);
  h_corr[2]->Draw("colz");
  c2->cd(7);
  h_corr[3]->Draw("colz");
  c2->cd(8);
  h_corr[4]->Draw("colz");
  c2->cd(9);
  h_corr[5]->Draw("colz");
  c2->Update();

  double par_fit[3], err_fit[3];
  par_fit[0] = h_corr[0]->ProjectionY()->GetMean();
  err_fit[0] = h_corr[0]->ProjectionY()->GetRMS();
  par_fit[1] = h_corr[2]->ProjectionY()->GetMean();
  err_fit[1] = h_corr[2]->ProjectionY()->GetRMS();
  par_fit[2] = h_corr[5]->ProjectionY()->GetMean();
  err_fit[2] = h_corr[5]->ProjectionY()->GetRMS();
  cout << " ++++++++++ Fit Results ++++++++++ " << endl;
  for (int i=0;i<NP;i++) {
    cout << " Par[" << i << "] = " << par_fit[i] << "+/-" << err_fit[i] << endl;
  }
  cout << " +++++++++++++++++++++++++++++++++ " << endl;
  
  double la_p_corr[100], la_m[100], la_m_err[100], la_m_corr[100], la_m_corr_err[100];
  for(int k=0;k<gr_La->GetN();k++) {
    // Generate a Lambda decay
    TLorentzVector daughters[2];
    do {
      generateLambdaDecay(gr_La->GetX()[k], gr_La->GetY()[k], &daughters[0], &daughters[1]);
      TLorentzVector pion = daughters[1];
      TLorentzVector proton = daughters[0];
      //	cout << " pion " << pion.X() << " " << pion.Y() << " " << pion.Z() << " " << pion.E() << "\t proton " << proton.X() << " " << proton.Y() << " " << proton.Z() << " " << proton.E() << endl;
    } while (daughters[0].Pt()<0.1 || daughters[1].Pt()<0.1 || daughters[0].Eta()<-2 || daughters[0].Eta()>0 || daughters[1].Eta()<-2 || daughters[1].Eta()>0);
    
    TLorentzVector pion = daughters[1];
    TLorentzVector proton = daughters[0];
    
    // calculate the Prob at the current par[0] position
    double pion_p_corr = ( pion.P() - ( par_fit[0]*f_corr_pi->Eval(pion.P()/M_Pi) + par_fit[1] ) ) * par_fit[2];
    double s_pion = pion_p_corr/pion.P();
    TLorentzVector pion_corr(0.,0.,0.,0.);
    pion_corr.SetXYZM(pion.X()*s_pion, pion.Y()*s_pion, pion.Z()*s_pion, M_Pi);
    
    double proton_p_corr = ( proton.P() - ( par_fit[0]*f_corr_pi->Eval(proton.P()/M_Pi) + par_fit[1] + gr_p_dcorr->Eval(proton.P()/M_P) ) ) * par_fit[2];
    double s_proton = proton_p_corr/proton.P();
    TLorentzVector proton_corr(0.,0.,0.,0.);
    proton_corr.SetXYZM(proton.X()*s_proton, proton.Y()*s_proton, proton.Z()*s_proton, M_P);
    
    TLorentzVector Lambda_corr = proton_corr + pion_corr;

    la_p_corr[k] = gr_La->GetX()[k];
    la_m[k] = gr_La->GetY()[k];
    la_m_err[k] = gr_La->GetEY()[k]*s_sig;

    la_m_corr[k] = Lambda_corr.M();
    la_m_corr_err[k] = gr_La->GetEY()[k]*s_sig;    
  }  

  TGraphErrors *gr_la_m = new TGraphErrors(gr_La->GetN(), la_p_corr, la_m, 0, la_m_err);
  gr_la_m->Print();
  TGraphErrors *gr_la_m_corr = new TGraphErrors(gr_La->GetN(), la_p_corr, la_m_corr, 0, la_m_corr_err);
  gr_la_m_corr->Print();
  
  
  TCanvas *c1 = new TCanvas("c1", "c1", 600, 600, 600, 600);
  c1->Draw();

  c1->cd();
  /*
  TH1D *h0 = new TH1D("h0", "", 1, 0, 5.0);
  h0->SetMinimum(0.49);
  h0->SetMaximum(0.51);
  h0->Draw("c");

  gr_Ks->SetMarkerStyle(20);
  gr_Ks->SetMarkerSize(1.5);
  gr_Ks->SetLineWidth(2);
  gr_Ks->Draw("p");
  */
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
  
  
}
