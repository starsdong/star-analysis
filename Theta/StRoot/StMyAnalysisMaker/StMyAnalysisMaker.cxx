#include "StMyAnalysisMaker.h"
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "PhysicalConstants.h"
#include "phys_constants.h"
#include "StLorentzVector.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TRVector.h"
#include "TLorentzVector.h"
#include "StPhysicalHelixD.hh"
#include "StDcaService.h"

ClassImp(StMyAnalysisMaker)

const Int_t StMyAnalysisMaker::PDGID[NFLAG][NBODY] = {{211, -211},{2212, -211}};
const Double_t StMyAnalysisMaker::ParMass[NFLAG][NBODY] = {{M_PION_PLUS, M_PION_MINUS},   // Ks
                                                           {M_PROTON, M_PION_MINUS}};  // Lambda

//-----------------------------------------------------------------------------
StMyAnalysisMaker::StMyAnalysisMaker(const char* name, StPicoDstMaker *picoMaker, const char* outName)
  : StMaker(name)
{
  mPicoDstMaker = picoMaker;
  mPicoDst = 0;
  mOutName = outName;
}

//----------------------------------------------------------------------------- 
StMyAnalysisMaker::~StMyAnalysisMaker()
{ /*  */ }

//----------------------------------------------------------------------------- 
Int_t StMyAnalysisMaker::Init() {
  DeclareHistograms();
  return kStOK;
}

//----------------------------------------------------------------------------- 
Int_t StMyAnalysisMaker::Finish() {
  if(mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(),"RECREATE");
    fout->cd();
    WriteHistograms();
    fout->Close();
  }
  return kStOK;
}

//-----------------------------------------------------------------------------
void StMyAnalysisMaker::DeclareHistograms() {
  mVzNBTofMatch = new TH2F("VzNBTofMatch","",200,0,200,200,180,220);
  mVyVx = new TH2F("VyVx","",200,-5.,5.,200,-5.,5.);
  mNBTofHitvsNT0 = new TH2F("NBTofHitvsNT0","",200,0,200,500,0,500);
  mNBTofHitvsNBTofMatch = new TH2F("NBTofHitvsNBTofMatch","",200,0,200,500,0,500);
  mNBTofHitvsNBTofPidTraits = new TH2F("NBTofHitvsNBTofPidTraits","",200,0,500,200,0,500);
  mNBTofHitvsFXTMult = new TH2F("NBTofHitvsFXTMult","",500,0,500,500,0,500);
  mNT0vsFXTMult = new TH2F("NT0vsFXTMult","",500,0,500,200,0,200);
  mDedxvsp = new TH2F("Dedxvsp","",500,-5.,5.,500,0.,10.);
  mLogDedxvsp = new TH2F("LogDedxvsp","",500,-5.,5.,500,0.,2.5);
  mNSigPivsp = new TH2F("NSigPivsp","",500,-5.,5.,500,-10.,10.);
  mNSigPvsp = new TH2F("NSigPvsp","",500,-5.,5.,500,-10.,10.);
  mInvBetavsp = new TH2F("InvBetavsp","",500,-5.,5.,500,0.,5.0);
  mInvBetaDiffPivsNT0 = new TH2F("InvBetaDiffPivsNT0","",100,0,100,500,-0.25,0.25);
  mInvBetaDiffPvsNT0 = new TH2F("InvBetaDiffPvsNT0","",100,0,100,500,-0.25,0.25);
  mLogInvBetavsp = new TH2F("LogInvBetavsp","",500,-5.,5.,500,-0.5,1.5);
  mDedxvsp_T0old = new TH2F("Dedxvsp_T0old","",500,-5.,5.,500,0.,10.);
  mDedxvsp_T0Pi = new TH2F("Dedxvsp_T0Pi","",500,-5.,5.,500,0.,10.);
  mDedxvsp_T0P = new TH2F("Dedxvsp_T0P","",500,-5.,5.,500,0.,10.);
  mKsMassPt = new TH2F("KsMassPt","",400,0.,4.0,400,0.3,0.8);
  mPKsMassPt = new TH2F("PKsMassPt","",400,0.,4.0,400,1.4,1.8);
  
  mTree = new TTree("mTree","V0 Tree");
  mTree->Branch("mRunId",&v0Tree.mRunId,"mRunId/I");
  mTree->Branch("mEvtId",&v0Tree.mEvtId,"mEvtId/I");
  mTree->Branch("mVx",&v0Tree.mVx,"mVx/F");
  mTree->Branch("mVy",&v0Tree.mVy,"mVy/F");
  mTree->Branch("mVz",&v0Tree.mVz,"mVz/F");
  mTree->Branch("mFxtMult",&v0Tree.mFxtMult,"mFxtMult/S");
  mTree->Branch("mNV0",&v0Tree.mNV0,"mNV0/I");
  mTree->Branch("mFlag",v0Tree.mFlag,"mFlag[mNV0]/b");
  mTree->Branch("mPt1",v0Tree.mPt1,"mPt1[mNV0]/F");
  mTree->Branch("mPt2",v0Tree.mPt2,"mPt2[mNV0]/F");
  mTree->Branch("mEta1",v0Tree.mEta1,"mEta1[mNV0]/F");
  mTree->Branch("mEta2",v0Tree.mEta2,"mEta2[mNV0]/F");
  mTree->Branch("mPhi1",v0Tree.mPhi1,"mPhi1[mNV0]/F");
  mTree->Branch("mPhi2",v0Tree.mPhi2,"mPhi2[mNV0]/F");  
  mTree->Branch("mDca1",v0Tree.mDca1,"mDca1[mNV0]/F");
  mTree->Branch("mDca2",v0Tree.mDca2,"mDca2[mNV0]/F");
  mTree->Branch("mDca12",v0Tree.mDca12,"mDca12[mNV0]/F");
  mTree->Branch("mDecayL",v0Tree.mDecayL,"mDecayL[mNV0]/F");
  mTree->Branch("mPt",v0Tree.mPt,"mPt[mNV0]/F");
  mTree->Branch("mEta",v0Tree.mEta,"mEta[mNV0]/F");
  mTree->Branch("mPhi",v0Tree.mPhi,"mPhi[mNV0]/F");
  mTree->Branch("mMass",v0Tree.mMass,"mMass[mNV0]/F");
  mTree->Branch("mDca2Vtx",v0Tree.mDca2Vtx,"mDca2Vtx[mNV0]/F");
  mTree->Branch("mTheta",v0Tree.mTheta,"mTheta[mNV0]/F");

  
}

//-----------------------------------------------------------------------------
void StMyAnalysisMaker::WriteHistograms() {
  mVzNBTofMatch->Write();
  mVyVx->Write();
  mNBTofHitvsNT0->Write();
  mNBTofHitvsNBTofMatch->Write();
  mNBTofHitvsNBTofPidTraits->Write();
  mNBTofHitvsFXTMult->Write();
  mNT0vsFXTMult->Write();
  mDedxvsp->Write();
  mDedxvsp_T0old->Write();
  mDedxvsp_T0Pi->Write();
  mDedxvsp_T0P->Write();
  mLogDedxvsp->Write();
  mNSigPivsp->Write();
  mNSigPvsp->Write();
  mInvBetavsp->Write();
  mInvBetaDiffPivsNT0->Write();
  mInvBetaDiffPvsNT0->Write();
  mLogInvBetavsp->Write();
  mKsMassPt->Write();
  mPKsMassPt->Write();
  mTree->Write();
}  

//----------------------------------------------------------------------------- 
void StMyAnalysisMaker::Clear(Option_t *opt) {
}

//----------------------------------------------------------------------------- 
Int_t StMyAnalysisMaker::Make() {
  if(!mPicoDstMaker) {
    LOG_WARN << " No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  mPicoDst = mPicoDstMaker->picoDst();
  if(!mPicoDst) {
    LOG_WARN << " No PicoDst! Skip! " << endm;
    return kStWarn;
  }

  StPicoEvent *mPicoEvent = mPicoDst->event();
  if(!mPicoEvent) {
    LOG_WARN << " No PicoEvent! Skip! " << endm;
    return kStWarn;
  }

  if(!mPicoEvent->isTrigger(620052)) return kStOK;
  
  int nTofT0 = mPicoEvent->nTofT0();
  TVector3 pVtx = mPicoEvent->primaryVertex();
  int nBTofMatch = mPicoEvent->nBTOFMatch();    
  
  mVzNBTofMatch->Fill(nBTofMatch, pVtx[2]);
  mVyVx->Fill(pVtx[0], pVtx[1]);
  mNBTofHitvsNT0->Fill(nTofT0, mPicoDst->numberOfBTofHits());
  mNBTofHitvsNBTofMatch->Fill(nBTofMatch, mPicoDst->numberOfBTofHits());
  mNBTofHitvsNBTofPidTraits->Fill(mPicoDst->numberOfBTofPidTraits(), mPicoDst->numberOfBTofHits());  
  
  if(pVtx[2]<198 || pVtx[2]>202) return kStOK;
  if(sqrt(pVtx[0]*pVtx[0]+(pVtx[1]+2.)*(pVtx[1]+2.))>1.5) return kStOK;
  
  int fxtmult = 0;
  for(int i=0;i<(int)mPicoDst->numberOfTracks();i++) {
    StPicoTrack *t = (StPicoTrack *)mPicoDst->track(i);
    if(!t) continue;
    float p = t->pMom().Mag();
    if(p<1.e-5) continue;
    int q = t->charge();
    if(q==0) continue;
    int nHits = t->nHitsFit();
    int nHitsMax = t->nHitsMax();
    float dca = t->gDCA(pVtx).Mag();
    if(1.0*nHits/(1.0*nHitsMax)<0.51) continue;
    if(dca>3.0) continue;
    fxtmult++;
    
    float dedx = t->dEdx();
    float nsig_pi = t->nSigmaPion();
    float nsig_p = t->nSigmaProton();

    mDedxvsp->Fill(p*q, dedx);
    if(dedx>1.e-6) mLogDedxvsp->Fill(p*q,log10(dedx));
    mNSigPivsp->Fill(p*q, nsig_pi);
    mNSigPvsp->Fill(p*q, nsig_p);
    
    //
    if(p>0.2 && p<0.6 && fabs(nsig_pi)<2.0) {
      mDedxvsp_T0old->Fill(p*q, dedx);
      mDedxvsp_T0Pi->Fill(p*q, dedx);
    }
    else if( ((q<0 && p>0.6) || (q>0 && p>0.6 && p<1.0)) && fabs(nsig_pi)<2.0) mDedxvsp_T0Pi->Fill(p*q, dedx);
    else if( q>0 && fabs(nsig_p)<2.0) mDedxvsp_T0P->Fill(p*q, dedx);
    
    if(!t->isTofTrack()) continue;
    int index2btof = t->bTofPidTraitsIndex();
    StPicoBTofPidTraits *pid = mPicoDst->btofPidTraits(index2btof);
    if(!pid) continue;
    float beta = pid->btofBeta();
    if(beta<0) continue;
    mInvBetavsp->Fill(p*q, 1./beta);
    mLogInvBetavsp->Fill(p*q, log10(1./beta));
    
    float invbeta_pi = sqrt(0.13957*0.13957/p/p+1.);
    float invbeta_p = sqrt(0.93827*0.93827/p/p+1.);
    if( p>0.2 && p<1.0 && fabs(nsig_pi)<2.0) {
      mInvBetaDiffPivsNT0->Fill(nTofT0, 1./beta-invbeta_pi);
    }  
    if( q>0 && fabs(nsig_p)<2.0) {
      mInvBetaDiffPvsNT0->Fill(nTofT0, 1./beta-invbeta_p);    
    }
  }
  
  mNBTofHitvsFXTMult->Fill(fxtmult, mPicoDst->numberOfBTofHits());
  mNT0vsFXTMult->Fill(fxtmult, nTofT0);

  // find candidate daughters
  Float_t B = mPicoEvent->bField();
  vector<Int_t> pipos;
  vector<Int_t> pineg;
  vector<Int_t> proton;
  vector<Int_t> primp;
  
  pipos.clear();
  pineg.clear();
  proton.clear();
  primp.clear();
  
  for(int i=0;i<(int)mPicoDst->numberOfTracks();i++) {
    StPicoTrack *t = (StPicoTrack *)mPicoDst->track(i);
    if(!t) continue;
    int q = t->charge();
    if(q==0) continue;
    int nHits = t->nHitsFit();
    int nHitsMax = t->nHitsMax();
    float dca = t->gDCA(pVtx).Mag();
//    if(1.0*nHits/(1.0*nHitsMax)<0.51) continue;

    float dedx = t->dEdx();
    float nsig_pi = t->nSigmaPion();
    float nsig_p = t->nSigmaProton();

    if(fabs(nsig_pi)<2.5 && dca>0.7) {
      if(q==+1) {
        pipos.push_back(i);
      } else {
        pineg.push_back(i);
      }
    }
    if(fabs(nsig_p)<2.5 && q==+1 && dca>0.7) {
      proton.push_back(i);
    }
    if(fabs(nsig_p)<2.5 && q==+1 && dca<2.0) {
      primp.push_back(i);
    }
  }
  
  //FILL TREE WITH EVENT INFORMATION
  v0Tree.mRunId = mPicoEvent->runId();
  v0Tree.mEvtId = mPicoEvent->eventId();
  v0Tree.mVx = pVtx[0];
  v0Tree.mVy = pVtx[1];
  v0Tree.mVz = pVtx[2];
  v0Tree.mFxtMult = fxtmult;
  
//  cout<<"# of pi+ "<<pipos.size()<<" and pi- "<<pineg.size()<<" # of 2nd proton "<<proton.size()<<" primary proton" << primp.size() << endl;
  StThreeVectorF vtx(pVtx[0], pVtx[1], pVtx[2]);

  int nV0 = 0;
  for(size_t i=0;i<pipos.size();i++) {
    for(size_t j=0;j<pineg.size();j++) {
      if(makePair(mPicoDst, pipos[i], pineg[j], vtx, B, 0, primp, nV0)) nV0++;
    }
  }
  
  v0Tree.mNV0 = nV0;
//  mTree->Fill();
  
  return kStOK;
}

Bool_t StMyAnalysisMaker::makePair(StPicoDst *picoDst, int i_pos, int i_neg, StThreeVectorF vtx, Double_t B, int flag, vector<Int_t> proton, int iV0)
{
  if(!picoDst) return false;
  StPicoTrack *t_pos = (StPicoTrack*)picoDst->track(i_pos);
  StPicoTrack *t_neg = (StPicoTrack*)picoDst->track(i_neg);
  if(!t_pos || !t_neg) return false;
  if((t_pos->charge() + t_neg->charge())!=0) return false;
  
  StThreeVectorF mom_pos(t_pos->gMom().X(), t_pos->gMom().Y(), t_pos->gMom().Z());
  StThreeVectorF ori_pos(t_pos->origin().X(), t_pos->origin().Y(), t_pos->origin().Z());
  StPhysicalHelixD h_pos(mom_pos, ori_pos, B*kilogauss, +1);
  
  StThreeVectorF mom_neg(t_neg->gMom().X(), t_neg->gMom().Y(), t_neg->gMom().Z());
  StThreeVectorF ori_neg(t_neg->origin().X(), t_neg->origin().Y(), t_neg->origin().Z());
  StPhysicalHelixD h_neg(mom_neg, ori_neg, B*kilogauss, -1);

  StThreeVectorF avg_v0, moms2_pos, moms2_neg;
  float dca12_a, dca12_b;
/*  
//  if(!mDcaAlgoLong) {
    pair<double,double> s2 = h_pos.pathLengths(h_neg);
    StThreeVectorF dcas2_pos = h_pos.at(s2.first);
    StThreeVectorF dcas2_neg = h_neg.at(s2.second);
    dca12_a = (dcas2_pos - dcas2_neg).mag();

    avg_v0 = (dcas2_pos+dcas2_neg)*0.5;
    moms2_pos = h_pos.momentumAt(s2.first, B*kilogauss);
    moms2_neg = h_neg.momentumAt(s2.second, B*kilogauss);
//  }
*/
//  if(mDcaAlgoLong) 
     dca12_b = closestDistance(h_pos, h_neg, B, vtx, avg_v0, moms2_pos, moms2_neg);
  
//  cout << " DCA12 = " << dca12_a << " " << dca12_b << endl;
  float dca12 = dca12_b;
  
  if(dca12 > 2.0) return kFALSE;

  //decay lenght
  float avg_decaylen = (avg_v0 - vtx).mag();

  //calculate momentum                                                                                                                      
  double temp1 = h_pos.pathLength(avg_v0);
  StThreeVectorF avg_p1 = h_pos.momentumAt(temp1,B*kilogauss);

  double temp2 = h_neg.pathLength(avg_v0);
  StThreeVectorF avg_p2 = h_neg.momentumAt(temp2,B*kilogauss);

  StThreeVectorF avg_v0mom = avg_p1 + avg_p2; // momentum  vector of V0

  float Angle = (avg_v0-vtx).angle(avg_v0mom);
  float Dca2vtx = (avg_v0-vtx).mag()*TMath::Sin(Angle);
  if(Dca2vtx>2.0) return kFALSE;


  float e_pos = TMath::Sqrt(avg_p1.mag2() + ParMass[flag][0]*ParMass[flag][0]);
  float e_neg = TMath::Sqrt(avg_p2.mag2() + ParMass[flag][1]*ParMass[flag][1]);

  StLorentzVector<Float_t> fourMom_pos(avg_p1, e_pos);
  StLorentzVector<Float_t> fourMom_neg(avg_p2, e_neg);
  StLorentzVector<Float_t> fourMom_v0 = fourMom_pos + fourMom_neg; // four momentum of V0

  v0Tree.mFlag[iV0] = flag;
  v0Tree.mDca1[iV0] = fabs(h_pos.geometricSignedDistance(vtx));
  v0Tree.mDca2[iV0] = fabs(h_neg.geometricSignedDistance(vtx)); 
  v0Tree.mPt1[iV0] = avg_p1.perp();
  v0Tree.mPt2[iV0] = avg_p2.perp();
  v0Tree.mEta1[iV0] = avg_p1.pseudoRapidity();
  v0Tree.mEta2[iV0] = avg_p2.pseudoRapidity();
  v0Tree.mPhi1[iV0] = avg_p1.phi();
  v0Tree.mPhi2[iV0] = avg_p2.phi();
  v0Tree.mDca12[iV0]= dca12;
  v0Tree.mDecayL[iV0] = avg_decaylen;
  v0Tree.mPt[iV0] = fourMom_v0.perp();
  v0Tree.mEta[iV0] = fourMom_v0.pseudoRapidity();
  v0Tree.mPhi[iV0] = fourMom_v0.phi();
  v0Tree.mMass[iV0] = fourMom_v0.m();
  v0Tree.mDca2Vtx[iV0] = Dca2vtx;
  v0Tree.mTheta[iV0] = Angle;
  
  if(v0Tree.mDca1[iV0]>1.2 && v0Tree.mDca2[iV0]>1.2 && dca12<1.5 && Dca2vtx<1.0 && avg_p1.mag()<1.0) {
    mKsMassPt->Fill(fourMom_v0.perp(), fourMom_v0.m());
    if(fourMom_v0.m()>0.49 && fourMom_v0.m()<0.508) {
      for(size_t k=0;k<proton.size();k++) {
        StPicoTrack *tp = (StPicoTrack*)picoDst->track(proton[k]);
        if(!tp) continue;
        float p = tp->pMom().Mag();
        if(p<1.e-5) continue;
        StThreeVectorF mom_p(tp->pMom().X(),tp->pMom().Y(),tp->pMom().Z());
        float e_p = TMath::Sqrt(mom_p.mag2() + M_PROTON*M_PROTON);
        StLorentzVector<Float_t> fourMom_p(mom_p, e_p);
        StLorentzVector<Float_t> fourMom_pKs = fourMom_p + fourMom_v0;
        mPKsMassPt->Fill(fourMom_pKs.perp(), fourMom_pKs.m());
      }
    }
  }

  
  return kTRUE;

}