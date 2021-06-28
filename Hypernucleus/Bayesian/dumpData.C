void dumpData()
{
 // Lambda vs p
  TGraphErrors *gre = new TGraphErrors(20);
   gre->SetName("M_Uncorr_Lambda");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);
   gre->SetMarkerStyle(3);
   gre->SetPoint(0,0.125,1.1231);
   gre->SetPointError(0,0,4.27129e-05);
   gre->SetPoint(1,0.375,1.11618);
   gre->SetPointError(1,0,0.000138468);
   gre->SetPoint(2,0.625,1.11601);
   gre->SetPointError(2,0,1.08535e-05);
   gre->SetPoint(3,0.875,1.11589);
   gre->SetPointError(3,0,3.77992e-06);
   gre->SetPoint(4,1.125,1.11587);
   gre->SetPointError(4,0,2.68026e-06);
   gre->SetPoint(5,1.375,1.11586);
   gre->SetPointError(5,0,2.54277e-06);
   gre->SetPoint(6,1.625,1.11586);
   gre->SetPointError(6,0,2.77833e-06);
   gre->SetPoint(7,1.875,1.11585);
   gre->SetPointError(7,0,3.23459e-06);
   gre->SetPoint(8,2.125,1.11584);
   gre->SetPointError(8,0,3.98581e-06);
   gre->SetPoint(9,2.375,1.11584);
   gre->SetPointError(9,0,5.00734e-06);
   gre->SetPoint(10,2.625,1.11582);
   gre->SetPointError(10,0,6.43905e-06);
   gre->SetPoint(11,2.875,1.11583);
   gre->SetPointError(11,0,8.4499e-06);
   gre->SetPoint(12,3.125,1.11584);
   gre->SetPointError(12,0,1.14694e-05);
   gre->SetPoint(13,3.375,1.11585);
   gre->SetPointError(13,0,1.58178e-05);
   gre->SetPoint(14,3.625,1.11589);
   gre->SetPointError(14,0,2.30086e-05);
   gre->SetPoint(15,3.875,1.11585);
   gre->SetPointError(15,0,3.32718e-05);
   gre->SetPoint(16,4.125,1.11588);
   gre->SetPointError(16,0,5.34829e-05);
   gre->SetPoint(17,4.375,1.11591);
   gre->SetPointError(17,0,7.43407e-05);
   gre->SetPoint(18,4.625,1.11579);
   gre->SetPointError(18,0,0.000165387);
   gre->SetPoint(19,4.875,1.11594);
   gre->SetPointError(19,0,0.000125276);


// K_S0 vs p
   TGraphErrors *gre1 = new TGraphErrors(12);
   gre1->SetName("M_Uncorr_KS");
   gre1->SetTitle("Graph");
   gre1->SetFillColor(1);
   gre1->SetMarkerStyle(3);
   gre1->SetPoint(0,0.125,0.497053);
   gre1->SetPointError(0,0,0.00017417);
   gre1->SetPoint(1,0.375,0.497654);
   gre1->SetPointError(1,0,1.82465e-05);
   gre1->SetPoint(2,0.625,0.498007);
   gre1->SetPointError(2,0,7.82082e-06);
   gre1->SetPoint(3,0.875,0.498235);
   gre1->SetPointError(3,0,7.05824e-06);
   gre1->SetPoint(4,1.125,0.498323);
   gre1->SetPointError(4,0,7.5182e-06);
   gre1->SetPoint(5,1.375,0.498275);
   gre1->SetPointError(5,0,1.00961e-05);
   gre1->SetPoint(6,1.625,0.498117);
   gre1->SetPointError(6,0,1.5416e-05);
   gre1->SetPoint(7,1.875,0.498004);
   gre1->SetPointError(7,0,2.42703e-05);
   gre1->SetPoint(8,2.125,0.497875);
   gre1->SetPointError(8,0,3.86575e-05);
   gre1->SetPoint(9,2.375,0.497607);
   gre1->SetPointError(9,0,6.51931e-05);
   gre1->SetPoint(10,2.625,0.497864);
   gre1->SetPointError(10,0,0.000100459);
   gre1->SetPoint(11,2.875,0.497652);
   gre1->SetPointError(11,0,0.000179768);

   TFile *fout = new TFile("Mass_Lamdba_Ks_Uncorr1.root","recreate");
   gre->Write();
   gre1->Write();
   fout->Close();
}
