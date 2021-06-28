#include "style.C+"
#include "draw.C+"

void dEFit()
{
  style();
  const Double_t MassPi = 0.13957;
  const Double_t MassP = 0.93827;
  
  // (RcP-McP)*sin(theta)/q^2 vs RcP/m
  TF1 *f_corr_p = new TF1("corr_proton", "[0]*exp([1]+[2]*x)+[3]", 0.1, 100);
  f_corr_p->SetParameters(-4.593658e-01, -3.564137e-01, -6.922880e+00, -1.614707e-03);
  
  TF1 *f_corr_pi = new TF1("corr_pion", "[0]*exp([1]+[2]*x)+[3]", 0.1, 100);
  f_corr_pi->SetParameters(-2.797414e-02, -2.041152e+00, -8.372173e-01, 1.336872e-04);

  const Int_t N = 1000;
  double poverm[N], dcorr[N];
  for(int i=0;i<N;i++) {
    poverm[i] = (i+0.5)*0.1;   // proton p/m
    double povermpi = poverm[i]*MassP/MassPi;   // converted to pion p/m_pi
    dcorr[i] = f_corr_p->Eval(poverm[i]) - f_corr_pi->Eval(povermpi);  //  
  }

  TGraph *gr_p_dcorr = new TGraph(N, poverm, dcorr);

  TCanvas *c1 = new TCanvas("c1", "c1", 0, 0, 600, 600);
  c1->Draw();
  c1->SetLogx();
  c1->cd();

  TH1D *h0 = new TH1D("h0", "", 1, 0.5, 20.);
  h0->SetMinimum(-0.006);
  h0->SetMaximum(0.002);
  h0->Draw("c");

  f_corr_p->SetLineColor(4);
  f_corr_p->SetName("Corr_p");
  f_corr_p->Draw("same");
  
  f_corr_pi->SetLineColor(2);
  f_corr_pi->SetName("Corr_pi");
  f_corr_pi->Draw("same");

  gr_p_dcorr->SetName("dCorr_p");
  gr_p_dcorr->SetLineWidth(2);
  gr_p_dcorr->SetLineColor(1);
  gr_p_dcorr->Draw("c");

  TFile *fout = new TFile("dE_corr.root","recreate");
  f_corr_p->Write();
  f_corr_pi->Write();
  gr_p_dcorr->Write();
  fout->Close();
}
