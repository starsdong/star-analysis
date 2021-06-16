#ifndef ROOT_LambdaEvent
#define ROOT_LambdaEvent

#include "TObject.h"

class LambdaEvent : public TObject
{
  private:
    Float_t vtx_x;
    Float_t vtx_y;
    Float_t vtx_z;

    UInt_t Nprim;
    UInt_t trig;
    UInt_t Nrun;

    Float_t zdcx;
    Float_t bbcx;
    Int_t tofmult;

    Float_t Qx;
    Float_t Qy;

    Float_t Qx_west;
    Float_t Qy_west;

    Float_t Qx_east;
    Float_t Qy_east;

    Int_t centrality;
    Float_t reweight;
  public:
    LambdaEvent();
    LambdaEvent(Float_t mvtx_x, Float_t mvtx_y, Float_t mvtx_z,
	UInt_t mNprim, UInt_t mtrig, UInt_t mNrun, Float_t mzdcx, Float_t mbbcx, Int_t mtofmult,
	Float_t mQx, Float_t mQy,
	Float_t mQx_west, Float_t mQy_west,
	Float_t mQx_east, Float_t mQy_east,	         
	Int_t mcentrality, Float_t mreweight
	);

    void SetLambdaData(Float_t mvtx_x, Float_t mvtx_y, Float_t mvtx_z,
	UInt_t mNprim, UInt_t mtrig, UInt_t mNrun, Float_t mzdcx, Float_t mbbcx, Int_t mtofmult,
	Float_t mQx, Float_t mQy,
	Float_t mQx_west, Float_t mQy_west,
	Float_t mQx_east, Float_t mQy_east,		
	Int_t mcentrality, Float_t mreweight
	);

    Float_t VtxX() { return vtx_x; }
    Float_t VtxY() { return vtx_y; }
    Float_t VtxZ() { return vtx_z; }

    UInt_t NPrimaries() { return Nprim; }
    UInt_t Trigger() { return trig; }
    UInt_t NRun() { return Nrun; }
    Float_t ZDCx() { return zdcx; }
    Float_t BBCx() { return bbcx; }
    Int_t TOfMult() {return tofmult;}

    Float_t Q_x() { return Qx; }
    Float_t Q_y() { return Qy; }
    Float_t Q_xwest() { return Qx_west; }
    Float_t Q_ywest() { return Qy_west; }
    
    Float_t Q_xeast() { return Qx_east; }
    Float_t Q_yeast() { return Qy_east; }
    
    Int_t Centrality() { return centrality; }
    Float_t Reweight() { return reweight; }


    ClassDef(LambdaEvent,1) //Lambda event class
};

#endif
