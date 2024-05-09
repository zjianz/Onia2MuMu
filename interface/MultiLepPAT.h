// -*- C++ -*-
//
// Package:    MultiLepPAT
// Class:      MultiLepPAT
// 
/**\class MultiLepPAT MultiLepPAT.cc myAnalyzers/MultiLepPAT/src/MultiLepPAT.cc

 Description: <one line class summary>
Make rootTuple for JPsiKKK reconstruction

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  
//
//

#ifndef _MultiLepPAT_h
#define _MultiLepPAT_h

// system include files
#include <memory>
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // xining MINIAODtest

// user include files
#include "../interface/VertexReProducer.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/Error.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Vector3D.h"


#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/V0Candidate/interface/V0Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
//#include "RecoVertex/V0Producer/interface/V0Producer.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"

#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"
#include "PhysicsTools/RecoUtils/interface/CheckHitPattern.h"

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

#include <vector>
#include <utility>

#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include "DataFormats/MuonReco/interface/MuonSelectors.h"

//
// class decleration
//

using std::vector;
using namespace edm;
using namespace reco;
using namespace std;

class MultiLepPAT : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit MultiLepPAT(const ParameterSet&);
  ~MultiLepPAT();

  
private:
  virtual void beginJob() ;
  virtual void beginRun(Run const & iRun, EventSetup const& iSetup);
  virtual void analyze(const Event&, const EventSetup&);
  virtual void endJob() ;

 
//add token here
  edm::EDGetTokenT<L1GlobalTriggerReadoutRecord> gtRecordToken_;
  edm::EDGetTokenT<BeamSpot> gtbeamspotToken_;
  edm::EDGetTokenT<VertexCollection> gtprimaryVtxToken_;
  edm::EDGetTokenT<edm::View<pat::Muon> > gtpatmuonToken_; // MINIAOD
  edm::EDGetTokenT<edm::TriggerResults> gttriggerToken_;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate> > trackToken_; // MINIAOD
  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;

  // ----------member data ---------------------------

 virtual int GetHitsBefore(const EventSetup& setup, const vector<reco::Track>& trks, RefCountedKinematicTree& fitTr)
        {
        return -1;
        }

virtual int GetMissesAfter(const EventSetup& setup, const vector<reco::Track>& trks, RefCountedKinematicTree& VrtxTree)
        {
        return -1;
        }

  //get ctau from beamspot
  virtual double GetcTau(RefCountedKinematicVertex& decayVrtx, RefCountedKinematicParticle& kinePart, Vertex& bs)
  {	TVector3 vtx;
    TVector3 pvtx;
    vtx.SetXYZ((*decayVrtx).position().x(), (*decayVrtx).position().y(), 0);
    pvtx.SetXYZ(bs.position().x(), bs.position().y(), 0);
    VertexDistanceXY vdistXY;
    TVector3 pperp(kinePart->currentState().globalMomentum().x(),
		   kinePart->currentState().globalMomentum().y(), 0);
    
    TVector3 vdiff = vtx - pvtx;
    double cosAlpha = vdiff.Dot(pperp) / (vdiff.Perp() * pperp.Perp());
    Measurement1D distXY = vdistXY.distance(Vertex(*decayVrtx), Vertex(bs));
    double ctauPV = distXY.value() * cosAlpha * kinePart->currentState().mass() / pperp.Perp();
    return ctauPV;    
  }
  
  virtual double GetcTauErr(RefCountedKinematicVertex& decayVrtx, RefCountedKinematicParticle& kinePart, Vertex& bs)
  {       
    TVector3 pperp(kinePart->currentState().globalMomentum().x(),
		   kinePart->currentState().globalMomentum().y(), 0);
    AlgebraicVector vpperp(3);
    vpperp[0] = pperp.x();
    vpperp[1] = pperp.y();
    vpperp[2] = 0.;
    
    GlobalError v1e = (Vertex(*decayVrtx)).error();
    GlobalError v2e = bs.error();
    AlgebraicSymMatrix vXYe = asHepMatrix(v1e.matrix()) + asHepMatrix(v2e.matrix());
    double ctauErrPV = sqrt(vXYe.similarity(vpperp)) * kinePart->currentState().mass() / (pperp.Perp2());
    
    return ctauErrPV;    
  }
  


  double deltaR(double eta1, double phi1, double eta2, double phi2) {
    double deta = eta1 - eta2;
    double dphi = phi1 - phi2;
    while (dphi >   M_PI) dphi -= 2*M_PI;
    while (dphi <= -M_PI) dphi += 2*M_PI;
    return sqrt(deta*deta + dphi*dphi);
  }
  


  // ----------member data ---------------------------
  string proccessName_;
  HLTConfigProvider hltConfig_;

  InputTag hlTriggerResults_;
  InputTag inputGEN_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord>  magneticFieldToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord>  theTTBuilderToken_;
  string vtxSample;
  bool doMC;
  int MCParticle;
  bool doJPsiMassCost;
  int MuPixHits_c;
  int MuSiHits_c;
  double MuNormChi_c;
  double MuD0_c;

  double JMaxM_c;
  double JMinM_c;
  int PiSiHits_c;
  double MuPt_c;
  double JPiPiDR_c;
  double XPiPiDR_c;
  bool UseXDr_c;	
  double JPiPiMax_c;
  double JPiPiMin_c;
	
  bool resolveAmbiguity_; 
  bool addXlessPrimaryVertex_;
  vector<string>      TriggersForJpsi_;
  vector<string>      FiltersForJpsi_;
  vector<string>      TriggersForUpsilon_;
  vector<string>      FiltersForUpsilon_;

  int JpsiMatchTrig[50], UpsilonMatchTrig[50];

  vector<string>      TriggersForMatching_;
  vector<string>      FiltersForMatching_;
  int  MatchingTriggerResult[50];
  bool Debug_;
  double Chi_Track_;

  TTree* X_One_Tree_;

  unsigned int        runNum, evtNum, lumiNum;
  unsigned int        nGoodPrimVtx;
  
  vector<unsigned int>* trigRes;
  vector<std::string>* trigNames;
  vector<unsigned int>* L1TT;
  vector<std::string>* MatchTriggerNames;

  float               priVtxX, priVtxY, priVtxZ, priVtxXE, priVtxYE, priVtxZE, priVtxChiNorm, priVtxChi, priVtxCL;
  vector<float>       *PriVtxXCorrX, *PriVtxXCorrY, *PriVtxXCorrZ;
  vector<double>      *PriVtxXCorrEX, *PriVtxXCorrEY, *PriVtxXCorrEZ;
  vector<float>	      *PriVtxXCorrC2, *PriVtxXCorrCL;

  unsigned int         nMu;
  vector<float>       *muPx, *muPy, *muPz, *muD0, *muD0E, *muDz, *muChi2, *muGlChi2,  *mufHits;   
  vector<bool>        *muFirstBarrel, *muFirstEndCap;
  vector<float>       *muDzVtx, *muDxyVtx;   
  vector<int>         *muNDF, *muGlNDF, *muPhits, *muShits, *muGlMuHits, *muType, *muQual;
  vector<int>         *muTrack;
  vector<float>       *muCharge;
  vector<float>       *muIsoratio;
  vector<int>         *muIsGoodLooseMuon, *muIsGoodLooseMuonNew, *muIsGoodSoftMuonNewIlse,*muIsGoodSoftMuonNewIlseMod, *muIsGoodTightMuon,*muIsJpsiTrigMatch, *muIsUpsTrigMatch, *munMatchedSeg;
  vector<int>         *muIsPatLooseMuon, *muIsPatTightMuon, *muIsPatSoftMuon, *muIsPatMediumMuon;

  //for Maksat trigger match
  vector<int> *muUpsVrtxMatch, *muL3TriggerMatch;

  //added by zhenhu for MuonID
  vector<float>       *muMVAMuonID, *musegmentCompatibility; 
  //for Stoyan slope pull
  vector<float>  *mupulldXdZ_pos_noArb, *mupulldYdZ_pos_noArb;
  //addition 3,4,5
  vector<float>  *mupulldXdZ_pos_ArbDef, *mupulldYdZ_pos_ArbDef;
  vector<float>  *mupulldXdZ_pos_ArbST, *mupulldYdZ_pos_ArbST;
  vector<float>  *mupulldXdZ_pos_noArb_any, *mupulldYdZ_pos_noArb_any;
 



 //xining
  unsigned int         Jpsi1_nmumuonly;
  vector<float> *Jpsi1_mumuonlyMass,*Jpsi1_mumuonlyMassErr, 
    * Psi2S_JPiPiMass, *Psi2S_JPiPiMassErr, 
    * Psi2S_cs_JPiPiMass,* Psi2S_cs_JPiPiVtxProb,* Psi2S_cs_JPiPiFit_Chi2 ,*Psi2S_cs_JPiPiFit_ndof,
    *Psi2S_px,*Psi2S_py,*Psi2S_pz	,*Psi2S_cs_px,*Psi2S_cs_py,*Psi2S_cs_pz  ,
    *Psi2S_Pi1_pt, *Psi2S_Pi2_pt,*Psi2S_Pi1_deltaR, *Psi2S_Pi2_deltaR ,
    *Psi2S_mu1Idx, *Psi2S_mu2Idx, 
    *Jpsi2_mumuonlyMass34, *Jpsi2_mumuonlyMassErr34,
    *Jpsi2_mumuonly34VtxProb, *Jpsi2_mumuonly34Fit_Chi2, *Jpsi2_mumuonly34Fit_ndof, 
    *Jpsi2_mumuonly34_px, *Jpsi2_mumuonly34_py, *Jpsi2_mumuonly34_pz, 
    *Jpsi2_cs_mumuonlyMass34, *Jpsi2_cs_mumuonly34VtxProb, *Jpsi2_cs_mumuonly34Fit_Chi2, *Jpsi2_cs_mumuonly34Fit_ndof,
    *Jpsi2_cs_mumuonly34_px, *Jpsi2_cs_mumuonly34_py, *Jpsi2_cs_mumuonly34_pz, 
    *XMass, *XMassErr,  
    *Jpsi1_mumuonlyVtxProb, *Jpsi1_mumuonlyFit_Chi2 , *Jpsi1_mumuonlyFit_ndof, 
    *Psi2S_JPiPiVtxProb, *Psi2S_JPiPiFit_Chi2 ,*Psi2S_JPiPiFit_ndof, *Jpsi1_mumuonlyPx, *Jpsi1_mumuonlyPy, *Jpsi1_mumuonlyPz,  
    *Jpsi1_mumuonlymu1Idx, *Jpsi1_mumuonlymu2Idx, *Jpsi2_mumuonlymu3Idx, *Jpsi2_mumuonlymu4Idx,
    *X_mu1Idx, *X_mu2Idx,*X_mu3Idx, *X_mu4Idx,
    *X_mass, *X_Fit_VtxProb, *X_Fit_Chi2, *X_Fit_ndof,
    *X_JPiPi_mass, *X_JPiPi_VtxProb, *X_JPiPi_px, *X_JPiPi_py, *X_JPiPi_pz,
    *X_JPiPi_massErr, 
    *X_px, *X_py,*X_pz, *X_JPiPimass, *X_Jpsi1mass, *X_Jpsi1prob ,
    *X_Jpsi2mass,*X_Jpsi2prob, *X_Jpsi2px,  *X_Jpsi2py,  *X_Jpsi2pz,*X_Jpsi2massErr, 
    *X_4mumass, *X_4muprob,*X_4mupx, *X_4mupy,*X_4mupz,*X_4mumassErr, 
    *X_Jpsi2_mumuonly34VtxProb, *X_Psi2S_px, *X_Psi2S_py, *X_Psi2S_pz,
    *X_JPiPi_Pi1px, *X_JPiPi_Pi1py,*X_JPiPi_Pi1pz,
    *X_JPiPi_Pi2px, *X_JPiPi_Pi2py,*X_JPiPi_Pi2pz, *X_JPiPi_Pi1DeltaR, *X_JPiPi_Pi2DeltaR, *X_JPiPi_Pi1pt, *X_JPiPi_Pi2pt, 
    *X_Jpsi1px, *X_Jpsi1py, *X_Jpsi1pz, *X_Jpsi1massErr, 
    *X_mu12cs_px, *X_mu12cs_py, *X_mu12cs_pz, 
    *X_mu1px, *X_mu1py, *X_mu1pz,
    *X_mu2px, *X_mu2py, *X_mu2pz,
    *X_JPiPics_px, *X_JPiPics_py, *X_JPiPics_pz,
    *X_mu34cs_px, *X_mu34cs_py, *X_mu34cs_pz,
    *Jpsi1_mumuonlyctau, *Jpsi1_mumuonlyctauerr;  

  vector<int> *mumuonlymuoverlapped, *Jpsi1_mumuonlyChg;
  
//xining  

//added for mass constrain on 1208
 vector<float> *cs_X_mass, *cs_X_Fit_VtxProb, *cs_X_Fit_Chi2, *cs_X_Fit_ndof,
 *cs_X_mu1Idx, *cs_X_mu2Idx,*cs_X_mu3Idx, *cs_X_mu4Idx,
 *cs_X_Jpsi1mass, *cs_X_Jpsi1prob, *cs_X_Jpsi1px, *cs_X_Jpsi1py, *cs_X_Jpsi1pz, *cs_X_Jpsi1massErr,
 *cs_X_Jpsi2mass,*cs_X_Jpsi2prob, *cs_X_Jpsi2px,  *cs_X_Jpsi2py,  *cs_X_Jpsi2pz,*cs_X_Jpsi2massErr,
 *cs_X_JPiPi_mass, *cs_X_JPiPi_VtxProb, *cs_X_JPiPi_px, *cs_X_JPiPi_py, *cs_X_JPiPi_pz,*cs_X_JPiPi_massErr,
 *cs_X_JPiPinoMC_mass, *cs_X_JPiPinoMC_massErr,
 *cs_X_JPiPi_Pi1px, *cs_X_JPiPi_Pi1py,*cs_X_JPiPi_Pi1pz,
 *cs_X_JPiPi_Pi2px, *cs_X_JPiPi_Pi2py,*cs_X_JPiPi_Pi2pz;

//doMC
vector<float> 
*MC_X_px,
*MC_X_py,
*MC_X_pz,
*MC_X_mass,
*MC_Dau_jpsipx,
*MC_Dau_jpsipy,
*MC_Dau_jpsipz,
*MC_Dau_jpsimass,
*MC_Dau_psi2spx,
*MC_Dau_psi2spy,
*MC_Dau_psi2spz,
*MC_Dau_psi2smass,
*MC_Granddau_mu1px,
*MC_Granddau_mu1py,
*MC_Granddau_mu1pz,
*MC_Granddau_mu2px,
*MC_Granddau_mu2py,
*MC_Granddau_mu2pz,
*MC_Granddau_jpsipx,
*MC_Granddau_jpsipy,
*MC_Granddau_jpsipz,
*MC_Granddau_jpsimass,
*MC_Granddau_pi1px,
*MC_Granddau_pi1py,
*MC_Granddau_pi1pz,
*MC_Granddau_pi2px,
*MC_Granddau_pi2py,
*MC_Granddau_pi2pz,
*MC_Grandgranddau_mu3px,
*MC_Grandgranddau_mu3py,
*MC_Grandgranddau_mu3pz,
*MC_Grandgranddau_mu4px,
*MC_Grandgranddau_mu4py,
*MC_Grandgranddau_mu4pz;

vector<int>
*MC_X_chg,
*MC_Dau_jpsipdgId,
*MC_Dau_psi2spdgId,
*MC_Granddau_mu1pdgId,
*MC_Granddau_mu2pdgId,
*MC_Granddau_jpsipdgId,
*MC_Granddau_pi1pdgId,
*MC_Granddau_pi2pdgId,
*MC_Grandgranddau_mu3pdgId,
*MC_Grandgranddau_mu4pdgId;

vector<float> 
*Match_mu1px,
*Match_mu1py,
*Match_mu1pz,
*Match_mu2px,
*Match_mu2py,
*Match_mu2pz,
*Match_mu3px,
*Match_mu3py,
*Match_mu3pz,
*Match_mu4px,
*Match_mu4py,
*Match_mu4pz,

*Match_pi1px,
*Match_pi1py,
*Match_pi1pz,
*Match_pi2px,
*Match_pi2py,
*Match_pi2pz; 

  float mybxlumicorr,myrawbxlumi;
  ////////////////////////  

  



  
};

#endif
