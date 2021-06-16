#ifndef ROOT_LambdaTrack
#define ROOT_LambdaTrack

#include "TObject.h"

class LambdaTrack : public TObject
{
  private:
    Float_t m;
    Float_t pt;
    Float_t pz;
    Int_t ist;
    Int_t pxl;
    Float_t decaylen1;
    Float_t decaylen2;
    Float_t decaylen3;
    Float_t decaylenavg;
    Float_t dcabetdau01;
    Float_t dcabetdau02;
    Float_t dcabetdau03;
    Float_t dcabetdau12;
    Float_t dcabetdau31;
    Float_t dcabetdau32;
    Float_t poiangle;
    Float_t dca2vtx;
    Float_t dca1;
    Float_t dca2;
    Float_t dca3;
    Float_t Qxqx_pos; 
    Float_t Qyqy_pos; 
    Float_t Qxqx_neg; 
    Float_t Qyqy_neg;
    Float_t Qxqx; 
    Float_t Qyqy; 
    Float_t phimass;
    
  public:
    LambdaTrack();
    LambdaTrack(Float_t mm, 
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


	  );

    Float_t M() { return m; }
    Float_t Pt() { return pt; }
    Float_t Pz() { return pz; }
    Int_t Ist() {return ist;}
    Int_t Pxl() {return pxl;}
    Float_t Decaylen1() {return decaylen1;}
    Float_t Decaylen2() {return decaylen2;}
    Float_t Decaylen3() {return decaylen3;}
    Float_t DecaylenAvg() {return decaylenavg;}

    Float_t Dcabetdau01() {return dcabetdau01;}
    Float_t Dcabetdau02() {return dcabetdau02;}
    Float_t Dcabetdau03() {return dcabetdau03;}
    Float_t Dcabetdau12() {return dcabetdau12;}
    Float_t Dcabetdau31() {return dcabetdau31;}
    Float_t Dcabetdau32() {return dcabetdau32;}


    Float_t Poiangle() {  return poiangle;}
	
    Float_t Dca2vtx() {return dca2vtx;}

    Float_t Dca1() {return dca1;}
    Float_t Dca2() {return dca2;}
    Float_t Dca3() {return dca3;}
    Float_t QX_pos() {return Qxqx_pos ;} 
    Float_t QY_pos() {return Qyqy_pos ;} 
    Float_t QX_neg() {return Qxqx_neg ;} 
    Float_t QY_neg() {return Qyqy_neg ;} 
    Float_t QX() {return Qxqx ;} 
    Float_t QY() {return Qyqy ;} 

		
    Float_t Phimass() {return phimass;}
   


    ClassDef(LambdaTrack,2) //Lambda track class
};

#endif
