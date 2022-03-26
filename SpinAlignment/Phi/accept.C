Double_t powerlaw(Double_t *x, Double_t *par)
{
  double n = par[2];
  double p0 = par[1]*(n-3)/2.;
  double A = 2.*(n-1)*(n-2)/TMath::Pi()/(n-3)/(n-3)/par[1]/par[1]*par[0];
  if(x[0]<5.0) return A*x[0]*pow(1+x[0]/p0,-n);
  else return A*100.*x[0]*pow(1+x[0]/p0,-n);
}

Double_t tpcEffEtaFunc(Double_t *x, Double_t *par)
{
  Double_t y = (TMath::Abs(x[0]-par[0])-par[2])/par[3];
  return par[1]/(1.+TMath::Exp(y))+par[4];
}

Double_t tofEffEtaFunc(Double_t *x, Double_t *par)
{
  Double_t y = (TMath::Abs(x[0]-par[0])-par[2])/par[3];
  Double_t z = (TMath::Abs(x[0]+par[0])-par[2])/par[3];
  return par[1]/(1.+TMath::Exp(y))+par[1]/(1.+TMath::Exp(z))+par[4];
} // to account for the eta gap in the middle

void accept(const Int_t Nevt = 10000000)
{

  const Double_t MassPhi = 1.01946;
  const Double_t SigmaPhi = 0.00182; // 1.82 MeV ~ Gamma = 4.2 MeV
  const Double_t MassPi = 0.13957;
  const Double_t MassK = 0.49368;

  const Double_t tpcAveEff = 0.8; // tpc ave eff at plateu now assuming 80%
  TF1 *fTpcEffEta = new TF1("tpcEffEtaFun",tpcEffEtaFunc,-1.,1.,5.);
  fTpcEffEta->SetParameters(0., tpcAveEff, 1.05, 0.05, 0.); // initial
  // tof eff vs eta is evaluated from Run10 data
  TF1 *fTofEffEta = new TF1("tofEffEtaFun",tofEffEtaFunc,-1.,1.,5.);
  fTofEffEta->SetParameters(0.49, 0.67, 0.48, 0.044, 0.); // initial
  double tofEffAve = fTofEffEta->Integral(-1.,1.)/2.;
  
  // TFile *fin = new TFile("input/Eff_pik.root");
  // TH1D *hEff_Pi = (TH1D *)fin->Get("hpiplus1");
  // TH1D *hEff_K = (TH1D *)fin->Get("hkplus1");
  // TGraph *fTpcEff_Pi = new TGraph(hEff_Pi);
  // TGraph *fTpcEff_K = new TGraph(hEff_K);
      

  TF1 *funMom = new TF1("funMom",powerlaw,0.,10.,3);
  funMom->SetParameters(1, 1.0, 11);
  funMom->SetNpx(1000);

  TF1 *funCosTheta = new TF1("funCosTheta","(1-[0])+(3*[0]-1)*x",-1,1);
  funCosTheta->SetParameter(0, 1./3.);
  
  TRandom3 *gRandom = new TRandom3();

  //-------------------
  // Book histograms
  //-------------------
  TH2D *hPtY = new TH2D("hPtY","",200,0.,2.,200,-2.,2.);
  TH2D *hPtEta = new TH2D("hPtEta","",200,0.,2.,200,-2.,2.);
  TH2D *hPtCosTheta = new TH2D("hPtCosTheta","",200,0.,2.,200,-1.,1.);

  TH2D *hCosThetaDiff = new TH2D("hCosThetaDiff","",200,-1.,1.,200,-0.01,0.01);
  
  //|eta|<1
  TH2D *hPtYRc = new TH2D("hPtYRc","",200,0.,2.0,200,-2.,2.);
  TH2D *hPtEtaRc = new TH2D("hPtEtaRc","",200,0.,2.0,200,-2.,2.);
  TH2D *hPtCosThetaRc = new TH2D("hPtCosThetaRc","",200,0.,2.,200,-1.,1.);

  //|eta|<1 pT>0.2
  TH2D *hPtYRc1 = new TH2D("hPtYRc1","",200,0.,2.0,200,-2.,2.);
  TH2D *hPtEtaRc1 = new TH2D("hPtEtaRc1","",200,0.,2.0,200,-2.,2.);
  TH2D *hPtCosThetaRc1 = new TH2D("hPtCosThetaRc1","",200,0.,2.,200,-1.,1.);

  TH2D *hPtPhiK1 = new TH2D("hPtPhiK1","",200,0.,2.0,200,0.,2.0);
  TH2D *hPtPhiK2 = new TH2D("hPtPhiK2","",200,0.,2.0,200,0.,2.0);

  //|eta|<1 pT>0.2, Tof Efficiency vs eta
  TH2D *hPtYRc2 = new TH2D("hPtYRc2","",200,0.,2.0,200,-2.,2.);
  TH2D *hPtEtaRc2 = new TH2D("hPtEtaRc2","",200,0.,2.0,200,-2.,2.);
  TH2D *hPtCosThetaRc2 = new TH2D("hPtCosThetaRc2","",200,0.,2.,200,-1.,1.);

  //|eta|<1 pT>0.2, Tof Efficiency constant vs eta
  TH2D *hPtYRc3 = new TH2D("hPtYRc3","",200,0.,2.0,200,-2.,2.);
  TH2D *hPtEtaRc3 = new TH2D("hPtEtaRc3","",200,0.,2.0,200,-2.,2.);
  TH2D *hPtCosThetaRc3 = new TH2D("hPtCosThetaRc3","",200,0.,2.,200,-1.,1.);
  
  //|eta|<1 pT>0.2, Tof Efficiency vs eta, Tpc efficiency vs eta
  TH2D *hPtYRc4 = new TH2D("hPtYRc4","",200,0.,2.0,200,-2.,2.);
  TH2D *hPtEtaRc4 = new TH2D("hPtEtaRc4","",200,0.,2.0,200,-2.,2.);
  TH2D *hPtCosThetaRc4 = new TH2D("hPtCosThetaRc4","",200,0.,2.,200,-1.,1.);

  TLorentzVector Phimom(0.,0.,0.,0.);
  TLorentzVector K1mom(0.,0.,0.,0.);   // Kaon momentum in the lab frame
  TLorentzVector K2mom(0.,0.,0.,0.);   // Kaon momentum in the lab frame
  TLorentzVector PhimomRc(0.,0.,0.,0.);   // Phi momentum Rc in the lab frame
  TLorentzVector K1momStar(0.,0.,0.,0.);   // Kaon momentum Rc in the Phi rest frame

  for(int ievt = 0; ievt<Nevt; ievt++) {
    if(ievt%100000==0) cout << " processing " << ievt << "-th event ..." << endl;
    double MassMom = gRandom->Gaus(MassPhi, SigmaPhi);
    double pT = gRandom->Rndm()*2.; //funMom->GetRandom();
    //    double y = gRandom->Rndm()*2. - 1. ; // rapidity
    double y;
    do {
      y = gRandom->Gaus(0, 1.0);
    } while (fabs(y)>1.0);
    double phi = gRandom->Rndm()*TMath::Pi()*2.;
    double mT = TMath::Sqrt(pT*pT+MassMom*MassMom);
    double pz = mT * TMath::SinH(y);
    double p = TMath::Sqrt(pT*pT+pz*pz);
    double eta = 0.5*TMath::Log((p+pz)/(p-pz));
    Phimom.SetPtEtaPhiM(pT, eta, phi, MassMom);

    hPtY->Fill(pT, y);
    hPtEta->Fill(pT, eta);

    double p_pri = TMath::Sqrt((MassMom+2*MassK)*(MassMom-2*MassK))/2.;
    double phi_pri = gRandom->Rndm()*TMath::Pi()*2.0;   // phi x-z
    double costheta_pri = gRandom->Rndm()*2. - 1.;
    double py_pri = p_pri*costheta_pri;    // theta - angle between Y-axis (pol axis) and momentum
    double pxz_pri = TMath::Sqrt(p_pri*p_pri-py_pri*py_pri);
    double px_pri = pxz_pri*TMath::Cos(phi_pri);
    double pz_pri = pxz_pri*TMath::Sin(phi_pri);
    K1mom.SetXYZM(px_pri, py_pri, pz_pri, MassK);
    K2mom.SetXYZM(-px_pri, -py_pri, -pz_pri, MassK);
    hPtCosTheta->Fill(pT, costheta_pri);

    // transform to the lab frame
    K1mom.Boost(Phimom.BoostVector());
    K2mom.Boost(Phimom.BoostVector());

    PhimomRc = K1mom + K2mom;
    K1momStar = K1mom;
    K1momStar.Boost(-PhimomRc.BoostVector());
    double costheta_rc = K1momStar.Y()/K1momStar.Vect().Mag();
    hCosThetaDiff->Fill(costheta_pri, costheta_rc-costheta_pri);

    // Check the momentum
    if(0) {
      cout << " Phi momentum = " << Phimom.Px() << " " << Phimom.Py() << " " << Phimom.Pz() << " " << Phimom.E() << endl;
      cout << " K1 momentum = " << K1mom.Px() << " " << K1mom.Py() << " " << K1mom.Pz() << " " << K1mom.E() << endl;
      cout << " K2 momentum = " << K2mom.Px() << " " << K2mom.Py() << " " << K2mom.Pz() << " " << K2mom.E() << endl;
    }

    if(fabs(K1mom.Eta())<1. && fabs(K2mom.Eta())<1.) {
      hPtYRc->Fill(pT, y);
      hPtEtaRc->Fill(pT, eta);
      hPtCosThetaRc->Fill(pT, costheta_rc);

      if(K1mom.Pt()>0.2 && K2mom.Pt()>0.2) {
	hPtYRc1->Fill(pT, y);
	hPtEtaRc1->Fill(pT, eta);
	hPtCosThetaRc1->Fill(pT, costheta_rc);

	hPtPhiK1->Fill(pT, K1mom.Pt());
	hPtPhiK2->Fill(pT, K2mom.Pt());

	// if(gRandom->Rndm()<fTofEffEta->Eval(K1mom.Eta()) &&
	//    gRandom->Rndm()<fTofEffEta->Eval(K2mom.Eta()) ) {
	  hPtYRc2->Fill(pT, y);
	  hPtEtaRc2->Fill(pT, eta);
	  hPtCosThetaRc2->Fill(pT, costheta_rc);

	  // if(gRandom->Rndm()<fTpcEff_K->Eval(K1mom.Pt()) &&
	  //    gRandom->Rndm()<fTpcEff_Pi->Eval(K2mom.Pt()) ) {
	    hPtYRc4->Fill(pT, y);
	    hPtEtaRc4->Fill(pT, eta);
	    hPtCosThetaRc4->Fill(pT, costheta_rc);
	//	  }
	// }

	// if(gRandom->Rndm()<tofEffAve && gRandom->Rndm()<tofEffAve) {
	//   hPtYRc3->Fill(pT, y);
	//   hPtEtaRc3->Fill(pT, eta);
	// }

      }
    }

  }

  TFile *fout = new TFile("accept.root","recreate");
  hPtY->Write();
  hPtEta->Write();
  hPtCosTheta->Write();
  hPtYRc->Write();
  hPtEtaRc->Write();
  hPtCosThetaRc->Write();
  hCosThetaDiff->Write();
  hPtYRc1->Write();
  hPtEtaRc1->Write();
  hPtCosThetaRc1->Write();
  hPtYRc2->Write();
  hPtEtaRc2->Write();
  hPtCosThetaRc2->Write();
  hPtYRc3->Write();
  hPtEtaRc3->Write();
  hPtCosThetaRc3->Write();
  hPtYRc4->Write();
  hPtEtaRc4->Write();
  hPtCosThetaRc4->Write();

  hPtPhiK1->Write();
  hPtPhiK2->Write();

  fout->Close();

}
