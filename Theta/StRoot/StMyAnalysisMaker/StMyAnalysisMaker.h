#ifndef StMyAnalysisMaker_h
#define StMyAnalysisMaker_h

#include "StMaker.h"
#include "StThreeVectorF.hh"
class StPicoDst;
class StPicoDstMaker;
class TString;
class TH1F;
class TH2F;
class TTree;

class StMyAnalysisMaker : public StMaker {
  public:
    StMyAnalysisMaker(const char *name, StPicoDstMaker *picoMaker, const char *outName);
    virtual ~StMyAnalysisMaker();
    
    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();
    
    void    DeclareHistograms();
    void    WriteHistograms();
    
    Bool_t  makePair(StPicoDst *picoDst, int i1, int i2, StThreeVectorF vtx, Double_t B, int flag, vector<Int_t> proton, int iV0);
    
  private:
    StPicoDstMaker *mPicoDstMaker;
    StPicoDst      *mPicoDst;
    
    TString    mOutName;
    TTree  *mTree;
    TH2F   *mVzNBTofMatch;
    TH2F   *mVyVx;
    TH2F   *mNBTofHitvsNT0;
    TH2F   *mNBTofHitvsNBTofMatch;
    TH2F   *mNBTofHitvsNBTofPidTraits;
    TH2F   *mNBTofHitvsFXTMult;
    TH2F   *mNT0vsFXTMult;
    TH2F   *mDedxvsp;
    TH2F   *mLogDedxvsp;
    TH2F   *mNSigPivsp;
    TH2F   *mNSigPvsp;
    TH2F   *mInvBetavsp;
    TH2F   *mInvBetaDiffPivsNT0;
    TH2F   *mInvBetaDiffPvsNT0;
    TH2F   *mLogInvBetavsp;
    TH2F   *mNHitsvsEta;
    TH2F   *mDedxvsp_T0old;
    TH2F   *mDedxvsp_T0Pi;
    TH2F   *mDedxvsp_T0P;
    TH2F   *mKsMassPt;
    TH2F   *mPKsMassPt;
    
  struct   V0TREE{
    Int_t   mRunId;
    Int_t   mEvtId;
    Float_t mVx;
    Float_t mVy;
    Float_t mVz;
    Short_t mFxtMult;

    Int_t   mNV0;
    UChar_t mFlag[10000];
    Float_t mPt1[10000];
    Float_t mPt2[10000];
    Float_t mEta1[10000];
    Float_t mEta2[10000];
    Float_t mPhi1[10000];
    Float_t mPhi2[10000];
    Float_t mDca1[10000];
    Float_t mDca2[10000];
    Float_t mDca12[10000];
    Float_t mDecayL[10000];
    Float_t mPt[10000];
    Float_t mEta[10000];
    Float_t mPhi[10000];
    Float_t mMass[10000];
    Float_t mDca2Vtx[10000];
    Float_t mTheta[10000];
  };
  V0TREE v0Tree;

  /////////////////////////////////////
  // Flag of different combinations
  // 0 - pi+ pi-   (Ks decay)
  // 1 - p pi-     (Lambda decay)
  /////////////////////////////////////    
  enum{
    pos = 0,
    neg = 1
  };

  static const Int_t NFLAG = 2; 
  static const Int_t NBODY = 2;
  static const Int_t PDGID[NFLAG][NBODY];
  static const Double_t ParMass[NFLAG][NBODY];

    ClassDef(StMyAnalysisMaker, 1)
};

#endif
