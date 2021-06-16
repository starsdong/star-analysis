#include "LambdaTrack.h"
ClassImp(LambdaTrack)

LambdaTrack::LambdaTrack()
{
   m = 0.;
   pt = 0.;
   pz = 0.;
   ist=-2;
   pxl=-2;
   decaylen1=0.;
   decaylen2=0.;
   decaylen3=0.;
   decaylenavg=0.;
   dcabetdau01=0.;
   dcabetdau02=0.;
   dcabetdau03=0.;
   dcabetdau12=0.;
   dcabetdau31=0.;
   dcabetdau32=0.;
   poiangle=0;
   dca2vtx=0;
   dca1 = 0.;
   dca2 = 0.;
   dca3 = 0.;
   Qxqx_pos=0.;
   Qyqy_pos=0.;
   Qxqx_neg=0.;
   Qyqy_neg=0.;
   Qxqx=0.;
   Qyqy=0.;
   phimass = 0.;

}

LambdaTrack::LambdaTrack(Float_t mm, 
			 Float_t mpt, 
			 Float_t mpz, 
			 Int_t mist,
			 Int_t mpxl,
			 Float_t mdecaylen1,
			 Float_t mdecaylen2,
			 Float_t mdecaylen3,
			 Float_t mdecaylenavg,
			 Float_t mdcabetdau01,
			 Float_t mdcabetdau02,
			 Float_t mdcabetdau03,
			 Float_t mdcabetdau12,
			 Float_t mdcabetdau31,
			 Float_t mdcabetdau32,
			 Float_t mpoiangle,
			 Float_t mdca2vtx,
			 Float_t mdca1,
			 Float_t mdca2,
			 Float_t mdca3,
			 Float_t mQxqx_pos, 
			 Float_t mQyqy_pos, 
			 Float_t mQxqx_neg, 
			 Float_t mQyqy_neg, 
			 Float_t mQxqx, 
			 Float_t mQyqy, 
			 Float_t mphimass
) 
{
  m = mm;
  pt = mpt;
  pz = mpz;
  ist=mist;
  pxl=mpxl;
  decaylen1=mdecaylen1;
  decaylen2=mdecaylen2;
  decaylen3=mdecaylen3;
  decaylenavg=mdecaylenavg;
  
  dcabetdau01=mdcabetdau01;
  dcabetdau02=mdcabetdau02;
  dcabetdau03=mdcabetdau03;
  dcabetdau12=mdcabetdau12;
  dcabetdau31=mdcabetdau31;
  dcabetdau32=mdcabetdau32;

  poiangle=mpoiangle;
  dca2vtx= mdca2vtx;
  dca1=mdca1;		
  dca2=mdca2;		
  dca3=mdca3;
  Qxqx_pos=mQxqx_pos;
  Qyqy_pos= mQyqy_pos;	
  Qxqx_neg=mQxqx_neg;
  Qyqy_neg= mQyqy_neg;	
  Qxqx=mQxqx;
  Qyqy= mQyqy;
  phimass=mphimass;

}
