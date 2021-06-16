#include "LambdaEvent.h"
ClassImp(LambdaEvent)

LambdaEvent::LambdaEvent()
{
  vtx_x = 0.;
  vtx_y = 0.;
  vtx_z = 0.;

  Nprim = 0;
  trig = 0;
  Nrun = 0;
  zdcx=0.;
  bbcx=0.;
  tofmult=0;
  Qx = 0.;
  Qy = 0.;
  Qx_west = 0.;
  Qy_west = 0.;
  Qx_east = 0.;
  Qy_east = 0.;
  centrality = -2;
  reweight = 0.;
}

LambdaEvent::LambdaEvent(Float_t mvtx_x, Float_t mvtx_y, Float_t mvtx_z,
			 UInt_t mNprim, UInt_t mtrig, UInt_t mNrun, Float_t mzdcx, Float_t mbbcx, Int_t mtofmult,
			 Float_t mQx, Float_t mQy,
			 Float_t mQx_west, Float_t mQy_west,
			 Float_t mQx_east, Float_t mQy_east,
			 Int_t mcentrality, Float_t mreweight
    )
{
  vtx_x = mvtx_x;
  vtx_y = mvtx_y;
  vtx_z = mvtx_z;

  Nprim = mNprim;
  trig = mtrig;
  Nrun = mNrun;
  zdcx=mzdcx;
  bbcx=mbbcx;
  tofmult=mtofmult;
  Qx = mQx;
  Qy = mQy;
  
  Qx_west = mQx_west;
  Qy_west = mQy_west;
  
  Qx_east = mQx_east;
  Qy_east = mQy_east;
  
  centrality = mcentrality;
  reweight = mreweight;

}

void LambdaEvent::SetLambdaData(Float_t mvtx_x, Float_t mvtx_y, Float_t mvtx_z,
				UInt_t mNprim, UInt_t mtrig, UInt_t mNrun, Float_t mzdcx, Float_t mbbcx, Int_t mtofmult,
				Float_t mQx, Float_t mQy,
				Float_t mQx_west, Float_t mQy_west,
				Float_t mQx_east, Float_t mQy_east,
				Int_t mcentrality, Float_t mreweight

    )
{
  vtx_x = mvtx_x;
  vtx_y = mvtx_y;
  vtx_z = mvtx_z;

  Nprim = mNprim;
  trig = mtrig;
  Nrun = mNrun;
  zdcx=mzdcx;
  bbcx=mbbcx;
  tofmult=mtofmult;
  Qx = mQx;
  Qy = mQy;
  

  Qx_west = mQx_west;
  Qy_west = mQy_west;
  

  Qx_east = mQx_east;
  Qy_east = mQy_east;
  
  centrality = mcentrality;

  reweight = mreweight;

}
