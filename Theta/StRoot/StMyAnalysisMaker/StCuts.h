#ifndef CUTS_H
#define CUTS_H

namespace cuts
{
  //constant 
  const Float_t onePi         = 3.14159265;
  const Float_t twoPi         = 6.28318530;
  const Float_t piMass        = 0.13957;
  const Float_t kMass         = 0.493677;
  const Float_t TofResolution = 0.012;
  //event
  const Float_t vz            = 30.0; //30 - for NPE analysis, 60 - for Dmeson analysis. 
  const Float_t vzVpdVz       = 6.0;

  //global tracking 
  const Int_t   nHitsFit         = 20; //20
  const Float_t nHitsFitnHitsMax = 0.52;
  const Float_t dca              = 2.0; //global Dca Cut.
  const Float_t eta              = 1.0;
  const UShort_t adc             = 18; //HT1 threshold 

  //primary electron and partner electron track cuts 
  const Float_t gelectronPt = 0.2;
  const Float_t gnSigmaElectronHigh = 5.0;
  const Float_t gnSigmaElectronLow = -5.0;
  const Float_t gelectronsBeta = 0.035; // for tracks with pT< 1.0
  const Float_t partnerElecNSigmaElectronLow = -2.0;
  const Float_t partnerElecNSigmaElectronHig = 3.0;
  const Float_t partnerElecGlobalPt = 0.15;

  //tof 
  const Float_t tofYlocal = 1.8;

  //bemc
  const Float_t ePtCut   = 2.0;
  const Float_t PEMin    = 0.3;
  const Float_t PEMax    = 1.5;
  const Int_t   NEta     = 1;
  const Int_t   NPhi     = 1;
  const Float_t enEtaMin = 2.0;
  const Float_t enPhiMin = 2.0;
  const Float_t ZDist    = 3;
  const Float_t PhiDist  = 0.015;
  const Float_t EEtaCut  = 2.6;

  // electrons cuts
  const Float_t eDca       = 1.5;
  const Int_t   eHitsFit   = 20;
  const Int_t   eHitsDedx  = 15;
  const Float_t eFirstR    = 73;//first point R 
  const Float_t ePt        = 0.2;
  const Float_t eEta       = 0.7; //1.0
  const Float_t nSigmaEHig = 3.0;
  const Float_t nSigmaELow = -1.0;
  const Float_t eBeta      = 0.03; //cut for "1/beta-1"
  const Float_t poelow     = 0.3;
  const Float_t poehig     = 1.5;

  //electron pair cuts 
  const Float_t pairmass   = 0.3; //0.3
  const Float_t pairdca    = 1.0;
  const Float_t photonmass = 0.20;//0.15 ; 0.10

  //pion pid
  const Float_t nSigmaPion  = 3.0;//2.0
  const Float_t pionPt      = 0.2; 
  const Float_t pionEta     = 1.0; 
  const Float_t nsigtofpion = 3.; 

  //kaon pid
  const Float_t kaonPt        = 0.25; //0.2GeV
  const Float_t kaonEta       = 1.0;
  const Float_t nSigmaKaon    = 2.0;//3.0
  const Float_t tofkaonptcuts = 2.07; //look at David's analysis note .

  //Dmeson track cuts 
  const Int_t   pikHitsFit = 20;
  const Float_t pidca      = 2.0;
  const Float_t spidca     = 2.0;
  const Float_t kdca       = 2.0;

  //D0 cuts 
  const Float_t d0Masslow        = 1.2;//1.6
  const Float_t d0Masshig        = 2.4;//2.2
  const Float_t d0deltamlow      = 0.144;
  const Float_t d0deltamhig      = 0.147; 
  const Float_t d0costhetaStarMB = 0.85;//MB Trigger 0.85 
  const Float_t d0costhetaStarHT = 0.6; //HT Trigger 0.6 
  const Float_t d0pionnsigpi   	 = 2;
  const Float_t d0kaonnsigk    	 = 2.0;

  const Float_t D0MassCut1   = 1.83;//1.82
  const Float_t D0MassCut2   = 1.90;//1.91
  const Float_t D0MassSBCut1 = 1.72;
  const Float_t D0MassSBCut2 = 1.80;
  const Float_t D0MassSBCut3 = 1.92;
  const Float_t D0MassSBCut4 = 2.00;

  //Dstar cuts
  const Float_t dstarMasslow = 1.6;
  const Float_t dstarMasshig = 2.4;
  const Float_t deltamlow    = 0.138;
  const Float_t deltamhig    = 0.178; 
  const Float_t rapidity     = 1.0;
  const Float_t dstarcosthetaStarHT = 0.6;//HT Trigger 0.6
  const Float_t dstarcosthetaStarMB = 0.77;//MB Trigger 0.85. reduce conbination background .
  //soft pion cuts
  const Float_t softpionpt        = 0.2; //0.15 
  const Float_t softpionnsigpi    = 2;
  const Float_t softpionnsigtofpi = 3;

  //topological cuts
  const Float_t d0ptOfspionptlow = 7;//10
  const Float_t d0ptOfspionpthig = 20;//20

  //eD* correlation D* mass window  
  //means and sigma are coming from David's fit results.
  //means = 145.4 MeV, sigma = 0.37 MeV(fit results)
  const Float_t DeltaMassLow = 0.1446; //0.1443(mean-3*sigma) , 0.1446(mean-2*sigma), 0.138(mean-5*sigma)--good
  const Float_t DeltaMassHig = 0.1461; //0.1465(mean+3*sigma) , 0.1461(mean+2*sigma), 0.152(means+5*sigma)--good

  //eD0 correlation D0 mass window
  //Range mean+/- 2*Sigma
  //Mean:1.8647 (from my fit result)
  //Sigma:1.573e-2 (from my fit result)
  const Float_t D0MassWinLow = 1.833;//1.800(mean-3*sigma);1.833 (mean-2*sigma)
  const Float_t D0MassWinHig = 1.895;//1.940(mena+3*sigma);1.895(mean+2*sigma)
  //Default value;
  const Float_t Default = -999;
}
#endif
