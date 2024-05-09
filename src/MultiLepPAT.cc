// system include files
#include <memory>
#include "TLorentzVector.h"
// user include files
#include "../interface/MultiLepPAT.h"
#include "../interface/VertexReProducer.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtFdlWord.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Math/interface/Error.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/DeepCopyPointer.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/CLHEP/interface/AlgebraicObjects.h"
#include "DataFormats/CLHEP/interface/Migration.h"

#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "MagneticField/Engine/interface/MagneticField.h"

//////////This is necessary for lumicalc///////
#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "RecoLuminosity/LumiProducer/interface/LumiCorrectionParam.h"
#include "RecoLuminosity/LumiProducer/interface/LumiCorrectionParamRcd.h"

#include "PhysicsTools/RecoUtils/interface/CheckHitPattern.h"

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include "TrackingTools/IPTools/interface/IPTools.h"

#include "../data/TMVAClassification_BDT.class.C"

// about photon
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include <boost/foreach.hpp>

#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // MINIAOD

typedef math::Error<3>::type CovarianceMatrix;
typedef ROOT::Math::SVector<double, 3> SVector3;
typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>> SMatrixSym3D;

// constructors and destructor
MultiLepPAT::MultiLepPAT(const edm::ParameterSet &iConfig)
	: hlTriggerResults_(iConfig.getUntrackedParameter<edm::InputTag>("HLTriggerResults", edm::InputTag("TriggerResults::HLT"))),
	  inputGEN_(iConfig.getUntrackedParameter<edm::InputTag>("inputGEN", edm::InputTag("genParticles"))),
	  magneticFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()),
	  theTTBuilderToken_(esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"))),
	  vtxSample(iConfig.getUntrackedParameter<std::string>("VtxSample", std::string("offlinePrimaryVertices"))),
	  doMC(iConfig.getUntrackedParameter<bool>("DoMonteCarloTree", false)),
	  MCParticle(iConfig.getUntrackedParameter<int>("MonteCarloParticleId", 20443)), // 20443 X, 100443 Psi(2S), 9120443  // X from B
	  doJPsiMassCost(iConfig.getUntrackedParameter<bool>("DoJPsiMassConstraint")),
	  MuPixHits_c(iConfig.getUntrackedParameter<int>("MinNumMuPixHits", 0)),
	  MuSiHits_c(iConfig.getUntrackedParameter<int>("MinNumMuSiHits", 0)),
	  MuNormChi_c(iConfig.getUntrackedParameter<double>("MaxMuNormChi2", 1000)),
	  MuD0_c(iConfig.getUntrackedParameter<double>("MaxMuD0", 1000)),
	  JMaxM_c(iConfig.getUntrackedParameter<double>("MaxJPsiMass", 4)),
	  JMinM_c(iConfig.getUntrackedParameter<double>("MinJPsiMass", 2.2)),
	  PiSiHits_c(iConfig.getUntrackedParameter<int>("MinNumTrSiHits", 0)),
	  MuPt_c(iConfig.getUntrackedParameter<double>("MinMuPt", 0)),
	  JPiPiDR_c(iConfig.getUntrackedParameter<double>("JPsiKKKMaxDR", 1)),
	  XPiPiDR_c(iConfig.getUntrackedParameter<double>("XCandPiPiMaxDR", 1.1)),
	  UseXDr_c(iConfig.getUntrackedParameter<bool>("UseXDr", false)),
	  JPiPiMax_c(iConfig.getUntrackedParameter<double>("JPsiKKKMaxMass", 50)),
	  JPiPiMin_c(iConfig.getUntrackedParameter<double>("JPsiKKKMinMass", 0)),
	  resolveAmbiguity_(iConfig.getUntrackedParameter<bool>("resolvePileUpAmbiguity", true)),
	  addXlessPrimaryVertex_(iConfig.getUntrackedParameter<bool>("addXlessPrimaryVertex", true)),
	  TriggersForJpsi_(iConfig.getUntrackedParameter<std::vector<std::string>>("TriggersForJpsi")),
	  FiltersForJpsi_(iConfig.getUntrackedParameter<std::vector<std::string>>("FiltersForJpsi")),
	  TriggersForUpsilon_(iConfig.getUntrackedParameter<std::vector<std::string>>("TriggersForUpsilon")),
	  FiltersForUpsilon_(iConfig.getUntrackedParameter<std::vector<std::string>>("FiltersForUpsilon")),
	  Debug_(iConfig.getUntrackedParameter<bool>("Debug_Output", false)),
	  Chi_Track_(iConfig.getUntrackedParameter<double>("Chi2NDF_Track", 10)),
	  X_One_Tree_(0),

	  runNum(0), evtNum(0), lumiNum(0), nGoodPrimVtx(0),
	  trigRes(0), trigNames(0), L1TT(0), MatchTriggerNames(0),

	  priVtxX(0), priVtxY(0), priVtxZ(0), priVtxXE(0), priVtxYE(0), priVtxZE(0), priVtxChiNorm(0), priVtxChi(0), priVtxCL(0),
	  PriVtxXCorrX(0), PriVtxXCorrY(0), PriVtxXCorrZ(0),
	  PriVtxXCorrEX(0), PriVtxXCorrEY(0), PriVtxXCorrEZ(0), PriVtxXCorrC2(0), PriVtxXCorrCL(0),

	  nMu(0),
	  muPx(0), muPy(0), muPz(0), muD0(0), muD0E(0), muDz(0), muChi2(0), muGlChi2(0), mufHits(0),
	  muFirstBarrel(0), muFirstEndCap(0), muDzVtx(0), muDxyVtx(0),
	  muNDF(0), muGlNDF(0), muPhits(0), muShits(0), muGlMuHits(0), muType(0), muQual(0),
	  muTrack(0), muCharge(0), muIsoratio(0), muIsGoodLooseMuon(0), muIsGoodLooseMuonNew(0),
	  muIsGoodSoftMuonNewIlse(0), muIsGoodSoftMuonNewIlseMod(0), muIsGoodTightMuon(0), muIsJpsiTrigMatch(0), muIsUpsTrigMatch(0), munMatchedSeg(0),

	  muIsPatLooseMuon(0), muIsPatTightMuon(0), muIsPatSoftMuon(0), muIsPatMediumMuon(0),
	  muUpsVrtxMatch(0), muL3TriggerMatch(0),

	  muMVAMuonID(0), musegmentCompatibility(0),
	  mupulldXdZ_pos_noArb(0), mupulldYdZ_pos_noArb(0),
	  mupulldXdZ_pos_ArbDef(0), mupulldYdZ_pos_ArbDef(0),
	  mupulldXdZ_pos_ArbST(0), mupulldYdZ_pos_ArbST(0),
	  mupulldXdZ_pos_noArb_any(0), mupulldYdZ_pos_noArb_any(0),
	  Jpsi1_nmumuonly(0),
	  Jpsi1_mumuonlyMass(0), Jpsi1_mumuonlyMassErr(0),

	  Jpsi1_mumuonlyVtxProb(0), Jpsi1_mumuonlyFit_Chi2(0), Jpsi1_mumuonlyFit_ndof(0),
	  Jpsi1_mumuonlyPx(0), Jpsi1_mumuonlyPy(0), Jpsi1_mumuonlyPz(0),
	  Jpsi1_mumuonlymu1Idx(0), Jpsi1_mumuonlymu2Idx(0),
	  X_mu1Idx(0), X_mu2Idx(0), X_mu3Idx(0), X_mu4Idx(0),
	  X_mass(0), X_Fit_VtxProb(0), X_Fit_Chi2(0), X_Fit_ndof(0),
	  X_JPiPi_mass(0), X_JPiPi_VtxProb(0), X_JPiPi_px(0), X_JPiPi_py(0), X_JPiPi_pz(0), X_JPiPi_massErr(0),
	  X_px(0), X_py(0), X_pz(0), X_JPiPimass(0),
	  X_Jpsi1mass(0), X_Jpsi1prob(0),
	  X_Jpsi2mass(0), X_Jpsi2prob(0), X_Jpsi2px(0), X_Jpsi2py(0), X_Jpsi2pz(0), X_Jpsi2massErr(0),
	  X_4mumass(0), X_4muprob(0), X_4mupx(0), X_4mupy(0), X_4mupz(0), X_4mumassErr(0),
	  X_JPiPi_Pi1px(0), X_JPiPi_Pi1py(0), X_JPiPi_Pi1pz(0),
	  X_JPiPi_Pi2px(0), X_JPiPi_Pi2py(0), X_JPiPi_Pi2pz(0), X_JPiPi_Pi1DeltaR(0), X_JPiPi_Pi2DeltaR(0),
	  X_JPiPi_Pi1pt(0), X_JPiPi_Pi2pt(0),
	  //
	  X_Jpsi1px(0), X_Jpsi1py(0), X_Jpsi1pz(0), X_Jpsi1massErr(0),
	  X_mu12cs_px(0), X_mu12cs_py(0), X_mu12cs_pz(0),
	  X_mu1px(0), X_mu1py(0), X_mu1pz(0),
	  X_mu2px(0), X_mu2py(0), X_mu2pz(0),
	  X_JPiPics_px(0), X_JPiPics_py(0), X_JPiPics_pz(0),
	  X_mu34cs_px(0), X_mu34cs_py(0), X_mu34cs_pz(0),
	  Jpsi1_mumuonlyctau(0), Jpsi1_mumuonlyctauerr(0),
	  mumuonlymuoverlapped(0),
	  Jpsi1_mumuonlyChg(0),
	  // mass constrain variables on 1208
	  cs_X_mass(0), cs_X_Fit_VtxProb(0), cs_X_Fit_Chi2(0), cs_X_Fit_ndof(0),
	  cs_X_mu1Idx(0), cs_X_mu2Idx(0), cs_X_mu3Idx(0), cs_X_mu4Idx(0),
	  cs_X_Jpsi1mass(0), cs_X_Jpsi1prob(0), cs_X_Jpsi1px(0), cs_X_Jpsi1py(0), cs_X_Jpsi1pz(0), cs_X_Jpsi1massErr(0),
	  cs_X_Jpsi2mass(0), cs_X_Jpsi2prob(0), cs_X_Jpsi2px(0), cs_X_Jpsi2py(0), cs_X_Jpsi2pz(0), cs_X_Jpsi2massErr(0),
	  cs_X_JPiPi_mass(0), cs_X_JPiPi_VtxProb(0), cs_X_JPiPi_px(0), cs_X_JPiPi_py(0), cs_X_JPiPi_pz(0), cs_X_JPiPi_massErr(0),
	  cs_X_JPiPinoMC_mass(0), cs_X_JPiPinoMC_massErr(0),
	  cs_X_JPiPi_Pi1px(0), cs_X_JPiPi_Pi1py(0), cs_X_JPiPi_Pi1pz(0),
	  cs_X_JPiPi_Pi2px(0), cs_X_JPiPi_Pi2py(0), cs_X_JPiPi_Pi2pz(0),

	  // doMC
	  MC_X_px(0),
	  MC_X_py(0),
	  MC_X_pz(0),
	  MC_X_mass(0),
	  MC_Dau_jpsipx(0),
	  MC_Dau_jpsipy(0),
	  MC_Dau_jpsipz(0),
	  MC_Dau_jpsimass(0),
	  MC_Dau_psi2spx(0),
	  MC_Dau_psi2spy(0),
	  MC_Dau_psi2spz(0),
	  MC_Dau_psi2smass(0),
	  MC_Granddau_mu1px(0),
	  MC_Granddau_mu1py(0),
	  MC_Granddau_mu1pz(0),
	  MC_Granddau_mu2px(0),
	  MC_Granddau_mu2py(0),
	  MC_Granddau_mu2pz(0),
	  MC_Granddau_jpsipx(0),
	  MC_Granddau_jpsipy(0),
	  MC_Granddau_jpsipz(0),
	  MC_Granddau_jpsimass(0),
	  MC_Granddau_pi1px(0),
	  MC_Granddau_pi1py(0),
	  MC_Granddau_pi1pz(0),
	  MC_Granddau_pi2px(0),
	  MC_Granddau_pi2py(0),
	  MC_Granddau_pi2pz(0),
	  MC_Grandgranddau_mu3px(0),
	  MC_Grandgranddau_mu3py(0),
	  MC_Grandgranddau_mu3pz(0),
	  MC_Grandgranddau_mu4px(0),
	  MC_Grandgranddau_mu4py(0),
	  MC_Grandgranddau_mu4pz(0),

	  MC_X_chg(0),
	  MC_Dau_jpsipdgId(0),
	  MC_Dau_psi2spdgId(0),
	  MC_Granddau_mu1pdgId(0),
	  MC_Granddau_mu2pdgId(0),
	  MC_Granddau_jpsipdgId(0),
	  MC_Granddau_pi1pdgId(0),
	  MC_Granddau_pi2pdgId(0),
	  MC_Grandgranddau_mu3pdgId(0),
	  MC_Grandgranddau_mu4pdgId(0),

	  Match_mu1px(0),
	  Match_mu1py(0),
	  Match_mu1pz(0),
	  Match_mu2px(0),
	  Match_mu2py(0),
	  Match_mu2pz(0),
	  Match_mu3px(0),
	  Match_mu3py(0),
	  Match_mu3pz(0),
	  Match_mu4px(0),
	  Match_mu4py(0),
	  Match_mu4pz(0),

	  Match_pi1px(0),
	  Match_pi1py(0),
	  Match_pi1pz(0),
	  Match_pi2px(0),
	  Match_pi2py(0),
	  Match_pi2pz(0),

	  mybxlumicorr(0), myrawbxlumi(0)
{
	// get token here for four-muon;
	gtRecordToken_ = consumes<L1GlobalTriggerReadoutRecord>(edm::InputTag("gtDigis"));
	gtbeamspotToken_ = consumes<BeamSpot>(edm::InputTag("offlineBeamSpot"));
	gtprimaryVtxToken_ = consumes<VertexCollection>(edm::InputTag("offlineSlimmedPrimaryVertices")); //  MINIAOD
	gtpatmuonToken_ = consumes<edm::View<pat::Muon>>(edm::InputTag("slimmedMuons"));				 //  MINIAOD
	gttriggerToken_ = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults::HLT"));
	trackToken_ = consumes<edm::View<pat::PackedCandidate>>(edm::InputTag("packedPFCandidates")); //  MINIAOD
	genParticlesToken_ = consumes<reco::GenParticleCollection>(edm::InputTag("genParticles"));
}

MultiLepPAT::~MultiLepPAT()
{
	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)
}
// member functions

//    ofstream myfile("comparison.txt");
// ------------ method called to for each event  ------------
void MultiLepPAT::analyze(const edm::Event &iEvent,
						  const edm::EventSetup &iSetup)
{
	// PDG 2020
	const double myJmass = 3.0969, myJmasserr = 0.00004;
	const double mypsi2smass = 3.686097, mypsi2smasserr = 0.00003;
	const double myMumass = 0.1056583745;
	const double myMumasserr = myMumass * 1e-6;
	const double myPimass = 0.13957039;
	// try
	const double myPimasserr = myPimass * 1e-6;

	TLorentzVector MC_mu1p4, MC_mu2p4, MC_mu3p4, MC_mu4p4, MC_pi1p4, MC_pi2p4;
	if (doMC)
	{
		edm::Handle<reco::GenParticleCollection> genParticles;
		iEvent.getByToken(genParticlesToken_, genParticles);

		for (const auto &particle : *(genParticles.product()))
		{
			if (std::abs(particle.pdgId()) == 35 && particle.numberOfDaughters() >= 2)
			{
				MC_X_px->push_back(particle.px());
				MC_X_py->push_back(particle.py());
				MC_X_pz->push_back(particle.pz());
				MC_X_mass->push_back(particle.mass());
				MC_X_chg->push_back(particle.charge());
				// particle.daughter(0)->pdgId() == 443 jpsi
				MC_Dau_jpsipdgId->push_back(particle.daughter(0)->pdgId());
				MC_Dau_jpsipx->push_back(particle.daughter(0)->px());
				MC_Dau_jpsipy->push_back(particle.daughter(0)->py());
				MC_Dau_jpsipz->push_back(particle.daughter(0)->pz());
				MC_Dau_jpsimass->push_back(particle.daughter(0)->mass());
				MC_Dau_psi2spdgId->push_back(particle.daughter(1)->pdgId());
				MC_Dau_psi2spx->push_back(particle.daughter(1)->px());
				MC_Dau_psi2spy->push_back(particle.daughter(1)->py());
				MC_Dau_psi2spz->push_back(particle.daughter(1)->pz());
				MC_Dau_psi2smass->push_back(particle.daughter(1)->mass());
				// particle.daughter(0)->daughter(1)->pdgId() == -13 mu+
				MC_Granddau_mu1pdgId->push_back(particle.daughter(0)->daughter(1)->pdgId());
				MC_Granddau_mu1px->push_back(particle.daughter(0)->daughter(1)->px());
				MC_Granddau_mu1py->push_back(particle.daughter(0)->daughter(1)->py());
				MC_Granddau_mu1pz->push_back(particle.daughter(0)->daughter(1)->pz());
				MC_Granddau_mu2pdgId->push_back(particle.daughter(0)->daughter(0)->pdgId());
				MC_Granddau_mu2px->push_back(particle.daughter(0)->daughter(0)->px());
				MC_Granddau_mu2py->push_back(particle.daughter(0)->daughter(0)->py());
				MC_Granddau_mu2pz->push_back(particle.daughter(0)->daughter(0)->pz());
				// particle.daughter(1)->daughter(0)->pdgId() == 443 jpsi from psi2s
				MC_Granddau_jpsipdgId->push_back(particle.daughter(1)->daughter(0)->pdgId());
				MC_Granddau_jpsipx->push_back(particle.daughter(1)->daughter(0)->px());
				MC_Granddau_jpsipy->push_back(particle.daughter(1)->daughter(0)->py());
				MC_Granddau_jpsipz->push_back(particle.daughter(1)->daughter(0)->pz());
				MC_Granddau_jpsimass->push_back(particle.daughter(1)->daughter(0)->mass());
				// particle.daughter(1)->daughter(1)->pdgId() == 211 pi+ from psi2s
				MC_Granddau_pi1pdgId->push_back(particle.daughter(1)->daughter(1)->pdgId());
				MC_Granddau_pi1px->push_back(particle.daughter(1)->daughter(1)->px());
				MC_Granddau_pi1py->push_back(particle.daughter(1)->daughter(1)->py());
				MC_Granddau_pi1pz->push_back(particle.daughter(1)->daughter(1)->pz());
				MC_Granddau_pi2pdgId->push_back(particle.daughter(1)->daughter(2)->pdgId());
				MC_Granddau_pi2px->push_back(particle.daughter(1)->daughter(2)->px());
				MC_Granddau_pi2py->push_back(particle.daughter(1)->daughter(2)->py());
				MC_Granddau_pi2pz->push_back(particle.daughter(1)->daughter(2)->pz());
				// particle.daughter(1)->daughter(0)->daughter(1)->pdgId() == -13 mu+ from jpsi from psi2s
				MC_Grandgranddau_mu3pdgId->push_back(particle.daughter(1)->daughter(0)->daughter(1)->pdgId());
				MC_Grandgranddau_mu3px->push_back(particle.daughter(1)->daughter(0)->daughter(1)->px());
				MC_Grandgranddau_mu3py->push_back(particle.daughter(1)->daughter(0)->daughter(1)->py());
				MC_Grandgranddau_mu3pz->push_back(particle.daughter(1)->daughter(0)->daughter(1)->pz());
				MC_Grandgranddau_mu4pdgId->push_back(particle.daughter(1)->daughter(0)->daughter(0)->pdgId());
				MC_Grandgranddau_mu4px->push_back(particle.daughter(1)->daughter(0)->daughter(0)->px());
				MC_Grandgranddau_mu4py->push_back(particle.daughter(1)->daughter(0)->daughter(0)->py());
				MC_Grandgranddau_mu4pz->push_back(particle.daughter(1)->daughter(0)->daughter(0)->pz());

				MC_mu1p4.SetXYZM((*MC_Granddau_mu1px)[0], (*MC_Granddau_mu1py)[0], (*MC_Granddau_mu1pz)[0], myMumass);
				MC_mu2p4.SetXYZM((*MC_Granddau_mu2px)[0], (*MC_Granddau_mu2py)[0], (*MC_Granddau_mu2pz)[0], myMumass);
				MC_mu3p4.SetXYZM((*MC_Grandgranddau_mu3px)[0], (*MC_Grandgranddau_mu3py)[0], (*MC_Grandgranddau_mu3pz)[0], myMumass);
				MC_mu4p4.SetXYZM((*MC_Grandgranddau_mu4px)[0], (*MC_Grandgranddau_mu4py)[0], (*MC_Grandgranddau_mu4pz)[0], myMumass);
				MC_pi1p4.SetXYZM((*MC_Granddau_pi1px)[0], (*MC_Granddau_pi1py)[0], (*MC_Granddau_pi1pz)[0], myPimass);
				MC_pi2p4.SetXYZM((*MC_Granddau_pi2px)[0], (*MC_Granddau_pi2py)[0], (*MC_Granddau_pi2pz)[0], myPimass);
			} // for (const auto& particle: *(genParticles.product()))
		}	  // if ( std::abs(particle.pdgId())  == 35 && particle.numberOfDaughters() ==2 )
	}		  // doMC

	// get event content information

	using std::vector;
	using namespace edm;
	using namespace reco;
	using namespace std;

	runNum = iEvent.id().run();
	evtNum = iEvent.id().event();
	lumiNum = iEvent.id().luminosityBlock();

	const MagneticField &bFieldHandle = iSetup.getData(magneticFieldToken_);

	edm::Handle<edm::TriggerResults> hltresults;
	try
	{
		iEvent.getByToken(gttriggerToken_, hltresults);
	}
	catch (...)
	{
		cout << "Couldn't get handle on HLT Trigger!" << endl;
	}
	if (!hltresults.isValid())
	{
		cout << "No Trigger Results!" << endl;
	}
	else
	{
		int ntrigs = hltresults->size();
		if (ntrigs == 0)
		{
			cout << "No trigger name given in TriggerResults of the input " << endl;
		}

		edm::TriggerNames triggerNames_;
		triggerNames_ = iEvent.triggerNames(*hltresults);

		int nUpstrigger = TriggersForUpsilon_.size();
		int nJpsitrigger = TriggersForJpsi_.size();

		for (int JpsiTrig = 0; JpsiTrig < nJpsitrigger; JpsiTrig++)
		{
			JpsiMatchTrig[JpsiTrig] = 0;
		} // jpsi trigger

		for (int UpsTrig = 0; UpsTrig < nUpstrigger; UpsTrig++)
		{
			UpsilonMatchTrig[UpsTrig] = 0;
		} // upsilon trig

		for (int itrig = 0; itrig < ntrigs; itrig++)
		{
			string trigName = triggerNames_.triggerName(itrig);
			int hltflag = (*hltresults)[itrig].accept();
			trigRes->push_back(hltflag);
			trigNames->push_back(trigName);

			for (unsigned int JpsiTrig = 0; JpsiTrig < TriggersForJpsi_.size(); JpsiTrig++)
			{
				if (TriggersForJpsi_[JpsiTrig] == triggerNames_.triggerName(itrig))
				{
					JpsiMatchTrig[JpsiTrig] = hltflag;
					break;
				}

			} // Jpsi Trigger

			for (unsigned int UpsTrig = 0; UpsTrig < TriggersForUpsilon_.size(); UpsTrig++)
			{
				if (TriggersForUpsilon_[UpsTrig] == triggerNames_.triggerName(itrig))
				{
					UpsilonMatchTrig[UpsTrig] = hltflag;
					break;
				}
			} // Upsilon Trigger
		}

		for (int MatchTrig = 0; MatchTrig < nJpsitrigger; MatchTrig++)
		{
			MatchTriggerNames->push_back(TriggersForJpsi_[MatchTrig]);
		}

	} // end of HLT trigger info

	std::string vrtxFilter("hltVertexmumuFilterUpsilonMuon");
	std::string L3Filter("hltTripleMuL3PreFiltered0");

	edm::ESHandle<TransientTrackBuilder> theB = iSetup.getHandle(theTTBuilderToken_);

	edm::Handle<L1GlobalTriggerReadoutRecord> gtRecord;
	iEvent.getByToken(gtRecordToken_, gtRecord);
	const DecisionWord dWord = gtRecord->decisionWord();

	const TechnicalTriggerWord ttWord = gtRecord->technicalTriggerWord();
	for (unsigned int l1i = 0; l1i != ttWord.size(); ++l1i)
	{
		L1TT->push_back(ttWord.at(l1i));
	}

	Vertex thePrimaryV;
	Vertex theRecoVtx;
	Vertex theBeamSpotV;
	BeamSpot beamSpot;
	math::XYZPoint RefVtx;

	// get BeamSplot
	edm::Handle<reco::BeamSpot> beamSpotHandle;
	iEvent.getByToken(gtbeamspotToken_, beamSpotHandle);
	if (beamSpotHandle.isValid())
	{
		beamSpot = *beamSpotHandle;
		theBeamSpotV = Vertex(beamSpot.position(), beamSpot.covariance3D());
	}
	else
		cout << "No beam spot available from EventSetup" << endl;

	Handle<VertexCollection> recVtxs;
	iEvent.getByToken(gtprimaryVtxToken_, recVtxs);
	unsigned int nVtxTrks = 0;

	///////////////////////////////////////////////////////////////////////
	////////////////Check Lines below for Primary Vertex///////////////////
	///////////////////////////////////////////////////////////////////////

	int mynGoodPrimVtx = 0;
	for (unsigned myi = 0; myi < recVtxs->size(); myi++)
	{
		if ((*recVtxs)[myi].ndof() >= 5 && fabs((*recVtxs)[myi].z()) <= 24 && fabs((*recVtxs)[myi].position().rho()) <= 2.0)
		{
			mynGoodPrimVtx++;
		}
	}
	nGoodPrimVtx = mynGoodPrimVtx;

	if (recVtxs->begin() != recVtxs->end())
	{
		if (addXlessPrimaryVertex_ || resolveAmbiguity_)
		{
			thePrimaryV = Vertex(*(recVtxs->begin()));
		}
		else
		{
			for (reco::VertexCollection::const_iterator vtx = recVtxs->begin(); vtx != recVtxs->end(); ++vtx)
			{
				if (nVtxTrks < vtx->tracksSize())
				{
					nVtxTrks = vtx->tracksSize();
					thePrimaryV = Vertex(*vtx);
				}
			}
		}
	}
	else
	{
		thePrimaryV = Vertex(beamSpot.position(), beamSpot.covariance3D());
	}

	RefVtx = thePrimaryV.position();
	priVtxX = (thePrimaryV.position().x());
	priVtxY = (thePrimaryV.position().y());
	priVtxZ = (thePrimaryV.position().z());
	priVtxXE = (thePrimaryV.xError());
	priVtxYE = (thePrimaryV.yError());
	priVtxZE = (thePrimaryV.zError());
	priVtxChiNorm = (thePrimaryV.normalizedChi2());
	priVtxChi = thePrimaryV.chi2();
	priVtxCL = ChiSquaredProbability((double)(thePrimaryV.chi2()), (double)(thePrimaryV.ndof()));

	edm::Handle<edm::View<pat::Muon>> thePATMuonHandle; //  MINIAOD
	iEvent.getByToken(gtpatmuonToken_, thePATMuonHandle);
	edm::Handle<edm::View<pat::PackedCandidate>> theTrackHandle; //  MINIAOD
	iEvent.getByToken(trackToken_, theTrackHandle);				 //  MINIAOD

	if (thePATMuonHandle->size() >= 2)
	{
		vector<std::string> theInputVariables;
		theInputVariables.push_back("validFrac");
		theInputVariables.push_back("globalChi2");
		theInputVariables.push_back("pt");
		theInputVariables.push_back("eta");
		theInputVariables.push_back("segComp");
		theInputVariables.push_back("chi2LocMom");
		theInputVariables.push_back("chi2LocPos");
		theInputVariables.push_back("glbTrackProb");
		theInputVariables.push_back("NTrkVHits");
		theInputVariables.push_back("NTrkEHitsOut");
		ReadBDT muonID(theInputVariables);
		vector<double> inputValues;
		inputValues.resize(10, 0.);
		// fill muon track block
		for (edm::View<pat::Muon>::const_iterator iMuonP = thePATMuonHandle->begin(); //  MINIAOD
			 iMuonP != thePATMuonHandle->end(); ++iMuonP)
		{
			// push back all muon information
			++nMu;
			muIsPatLooseMuon->push_back(iMuonP->isLooseMuon());
			muIsPatTightMuon->push_back(iMuonP->isTightMuon(thePrimaryV));
			muIsPatSoftMuon->push_back(iMuonP->isSoftMuon(thePrimaryV));
			muIsPatMediumMuon->push_back(iMuonP->isMediumMuon());

			muPx->push_back(iMuonP->px());
			muPy->push_back(iMuonP->py());
			muPz->push_back(iMuonP->pz());
			muCharge->push_back(iMuonP->charge());

			int goodSoftMuonNewIlseMod = 0;
			int goodSoftMuonNewIlse = 0;
			int goodLooseMuonNew = 0;
			int goodLooseMuon = 0;
			int goodTightMuon = 0;

			// float mymuMVABs = -1;

			bool JpsiTrigger = false;

			for (unsigned int JpsiTrig = 0; JpsiTrig < TriggersForJpsi_.size(); JpsiTrig++)
			{
				if (JpsiMatchTrig[JpsiTrig] != 0)
				{
					const pat::TriggerObjectStandAloneCollection muJpsiHLTMatches =
						iMuonP->triggerObjectMatchesByFilter(FiltersForJpsi_[JpsiTrig]);
					bool pass1 = muJpsiHLTMatches.size() > 0;
					if (pass1)
						JpsiTrigger = true;
				}
			}

			muIsJpsiTrigMatch->push_back(JpsiTrigger);

			bool UpsTrigger = false;

			for (unsigned int UpsTrig = 0; UpsTrig < TriggersForUpsilon_.size(); UpsTrig++)
			{
				if (UpsilonMatchTrig[UpsTrig] != 0)
				{
					const pat::TriggerObjectStandAloneCollection muUpsHLTMatches =
						iMuonP->triggerObjectMatchesByFilter(FiltersForUpsilon_[UpsTrig]);
					bool pass1 = muUpsHLTMatches.size() > 0;
					if (pass1)
						UpsTrigger = true;
				}
			}

			muIsUpsTrigMatch->push_back(UpsTrigger);

			munMatchedSeg->push_back(-1); // MINIOAOD

			int muL3TriMuonVrtxFilter = 0, muSingleMuL3Filter = 0;

			// Checking Single Trigger
			for (unsigned int UpsTrig = 0; UpsTrig < TriggersForUpsilon_.size(); UpsTrig++)
			{
				if (UpsilonMatchTrig[UpsTrig] != 0)
				{
					const pat::TriggerObjectStandAloneCollection muMatchVrxtFilter = iMuonP->triggerObjectMatchesByFilter(vrtxFilter);
					const pat::TriggerObjectStandAloneCollection muMatchL3Filter = iMuonP->triggerObjectMatchesByFilter(L3Filter);

					muL3TriMuonVrtxFilter = muMatchVrxtFilter.size() > 0;
					muSingleMuL3Filter = muMatchL3Filter.size() > 0;
				}
			}
			muUpsVrtxMatch->push_back(muL3TriMuonVrtxFilter); //  MINIAOD
			muL3TriggerMatch->push_back(muSingleMuL3Filter);  //  MINIAOD
		}
	} // if two muons

	if (doMC)
	{
		// pion loop
		for (View<pat::PackedCandidate>::const_iterator iTrack = theTrackHandle->begin(); // MINIAOD
			 iTrack != theTrackHandle->end(); ++iTrack)
		{
			TLorentzVector RECO_pip4;
			RECO_pip4.SetXYZM(iTrack->px(), iTrack->py(), iTrack->pz(), myPimass);
			if (fabs(MC_pi1p4.Pt() - RECO_pip4.Pt()) < 0.08 * MC_pi1p4.Pt() && MC_pi1p4.DeltaR(RECO_pip4) < 0.1)
			{
				Match_pi1px->push_back(RECO_pip4.Px());
				Match_pi1py->push_back(RECO_pip4.Py());
				Match_pi1pz->push_back(RECO_pip4.Pz());
			}
			if ((fabs(MC_pi2p4.Pt() - RECO_pip4.Pt()) < 0.08 * MC_pi2p4.Pt() && MC_pi2p4.DeltaR(RECO_pip4) < 0.1))
			{
				Match_pi2px->push_back(RECO_pip4.Px());
				Match_pi2py->push_back(RECO_pip4.Py());
				Match_pi2pz->push_back(RECO_pip4.Pz());
			}
		}
	} // if(doMC)

	// Data
	// It takes a lot less memory out of the loop
	KinematicConstraint *jpsi_cs = new MassKinematicConstraint(myJmass, myJmasserr);
	KinematicConstraint *JPiPi_cs = new MassKinematicConstraint(mypsi2smass, mypsi2smasserr);
	KinematicConstraint *jpsi_cs34 = new MassKinematicConstraint(myJmass, myJmasserr);
	double pionptcut = 0.25;
	double pionDRcut = 0.7;
	double vtxprobprecut = 1.0e-7;

	// get the four moun fit, but first also fit dimuon
	if (thePATMuonHandle->size() < 4)
	{
		return;
	}

	//  get X and MyFourMuon cands
	for (edm::View<pat::Muon>::const_iterator iMuon1 = thePATMuonHandle->begin(); // MINIAOD
		 iMuon1 != thePATMuonHandle->end(); ++iMuon1)
	{
		TrackRef muTrack1 = iMuon1->track();
		if (muTrack1.isNull())
		{
			continue;
		}

		reco::Track recoMu1 = *iMuon1->track();

		// next check for mu2
		for (edm::View<pat::Muon>::const_iterator iMuon2 = iMuon1 + 1; // MINIAOD
			 iMuon2 != thePATMuonHandle->end(); ++iMuon2)
		{
			TrackRef muTrack2 = iMuon2->track();
			if (muTrack2.isNull())
			{
				continue;
			}

			reco::Track recoMu2 = *iMuon2->track();

			if (!(1. < (iMuon1->p4() + iMuon2->p4()).mass() && (iMuon1->p4() + iMuon2->p4()).mass() < 4.))
			{
				continue;
			}
			if ((iMuon1->charge() + iMuon2->charge()) != 0)
			{
				continue;
			}

			TransientTrack muon1TT(muTrack1, &(bFieldHandle)); // MINIAOD
			TransientTrack muon2TT(muTrack2, &(bFieldHandle)); // MINIAOD
			KinematicParticleFactoryFromTransientTrack pmumuFactory;
			// The mass of a muon and the insignificant mass sigma
			// to avoid singularities in the covariance matrix.
			ParticleMass muon_mass = myMumass; // pdg mass
			float muon_sigma = myMumasserr;
			// initial chi2 and ndf before kinematic fits.
			float chi = 0.;
			float ndf = 0.;
			vector<RefCountedKinematicParticle> muonParticles;
			muonParticles.push_back(pmumuFactory.particle(muon1TT, muon_mass, chi, ndf, muon_sigma));
			muonParticles.push_back(pmumuFactory.particle(muon2TT, muon_mass, chi, ndf, muon_sigma));
			KinematicParticleVertexFitter fitter;
			RefCountedKinematicTree jpsi1VertexFitTree;
			jpsi1VertexFitTree = fitter.fit(muonParticles);

			if (!jpsi1VertexFitTree->isValid())
			{
				continue;
			}
			jpsi1VertexFitTree->movePointerToTheTop();
			RefCountedKinematicParticle jpsi1_vFit_noMC = jpsi1VertexFitTree->currentParticle();
			RefCountedKinematicVertex jpsi1_vFit_vertex_noMC = jpsi1VertexFitTree->currentDecayVertex();
			KinematicParameters mymumupara = jpsi1_vFit_noMC->currentState().kinematicParameters();
			double jpsi1_vtxprob = ChiSquaredProbability((double)(jpsi1_vFit_vertex_noMC->chiSquared()), (double)(jpsi1_vFit_vertex_noMC->degreesOfFreedom()));
			if (jpsi1_vFit_noMC->currentState().mass() > 3.4 || jpsi1_vFit_noMC->currentState().mass() < 2.8)
			{
				continue;
			}
			if (jpsi1_vtxprob < vtxprobprecut)
			{
				continue;
			}

			// mass constrain for jpsi1 from psi2s:
			KinematicParticleFitter csfitter;
			jpsi1VertexFitTree = csfitter.fit(jpsi_cs, jpsi1VertexFitTree);
			jpsi1VertexFitTree->movePointerToTheTop();
			RefCountedKinematicParticle jpsi1_vFit_cs = jpsi1VertexFitTree->currentParticle();
			RefCountedKinematicVertex jpsi1_vFit_vertex_cs = jpsi1VertexFitTree->currentDecayVertex();
			KinematicParameters mymumupara_cs = jpsi1_vFit_cs->currentState().kinematicParameters();
			double jpsi1_vtxprob_cs = ChiSquaredProbability((double)(jpsi1_vFit_vertex_cs->chiSquared()), (double)(jpsi1_vFit_vertex_cs->degreesOfFreedom()));

			// float myJpsi1_mumuonlyctau = GetcTau(jpsi1_vFit_vertex_noMC, jpsi1_vFit_noMC, theBeamSpotV);
			// float myJpsi1_mumuonlyctauerr = GetcTauErr(jpsi1_vFit_vertex_noMC, jpsi1_vFit_noMC, theBeamSpotV);
			Jpsi1_mumuonlyChg->push_back((iMuon1->charge() + iMuon2->charge()));

			// we need
			Jpsi1_mumuonlyVtxProb->push_back(jpsi1_vtxprob);
			Jpsi1_mumuonlyFit_Chi2->push_back(jpsi1_vFit_vertex_noMC->chiSquared());
			Jpsi1_mumuonlyFit_ndof->push_back(jpsi1_vFit_vertex_noMC->degreesOfFreedom());
			Jpsi1_mumuonlyMass->push_back(jpsi1_vFit_noMC->currentState().mass());
			if (jpsi1_vFit_noMC->currentState().kinematicParametersError().matrix()(6, 6) > 0)
			{
				Jpsi1_mumuonlyMassErr->push_back(sqrt(jpsi1_vFit_noMC->currentState().kinematicParametersError().matrix()(6, 6)));
			}
			else
			{
				Jpsi1_mumuonlyMassErr->push_back(-9);
			}
			Jpsi1_mumuonlyPx->push_back(mymumupara.momentum().x());
			Jpsi1_mumuonlyPy->push_back(mymumupara.momentum().y());
			Jpsi1_mumuonlyPz->push_back(mymumupara.momentum().z());
			Jpsi1_mumuonlymu1Idx->push_back(std::distance(thePATMuonHandle->begin(), iMuon1));
			Jpsi1_mumuonlymu2Idx->push_back(std::distance(thePATMuonHandle->begin(), iMuon2));
			++Jpsi1_nmumuonly;

			TLorentzVector P4_mu1;
			P4_mu1.SetPtEtaPhiM(iMuon1->track()->pt(), iMuon1->track()->eta(), iMuon1->track()->phi(), myMumass);
			TLorentzVector P4_mu2;
			P4_mu2.SetPtEtaPhiM(iMuon2->track()->pt(), iMuon2->track()->eta(), iMuon2->track()->phi(), myMumass);

			// mu3mu4(X6900->jpsi)
			for (edm::View<pat::Muon>::const_iterator iMuon3 = thePATMuonHandle->begin(); // MINIAOD
				 iMuon3 != thePATMuonHandle->end(); ++iMuon3)
			{
				TrackRef muTrack3 = iMuon3->track();
				if (muTrack3.isNull())
				{
					continue;
				}


				if (iMuon3 == iMuon1 || iMuon3 == iMuon2)
				{
					continue;
				}

				reco::Track recoMu3 = *iMuon3->track();
				for (edm::View<pat::Muon>::const_iterator iMuon4 = iMuon3 + 1; // MINIAOD
					 iMuon4 != thePATMuonHandle->end(); ++iMuon4)
				{
					TrackRef muTrack4 = iMuon4->track();
					if (muTrack4.isNull())
					{
						continue;
					}

					if (iMuon4 == iMuon1 || iMuon4 == iMuon2)
					{
						continue;
					}

					reco::Track recoMu4 = *iMuon4->track(); // MINIAOD

					if (!(1. < (iMuon3->p4() + iMuon4->p4()).mass() && (iMuon3->p4() + iMuon4->p4()).mass() < 4.5))
					{
						continue;
					}
					if ((iMuon3->charge() + iMuon4->charge()) != 0)
					{
						continue;
					}

					TransientTrack muon3TT(muTrack3, &(bFieldHandle)); // MINIAOD
					TransientTrack muon4TT(muTrack4, &(bFieldHandle)); // MINIAOD
					KinematicParticleFactoryFromTransientTrack pmumuFactory34;
					ParticleMass muon_mass = myMumass; // pdg mass
					float muon_sigma = myMumasserr;
					float chi = 0.;
					float ndf = 0.;
					vector<RefCountedKinematicParticle> muonParticles34;
					muonParticles34.push_back(pmumuFactory34.particle(muon3TT, muon_mass, chi, ndf, muon_sigma));
					muonParticles34.push_back(pmumuFactory34.particle(muon4TT, muon_mass, chi, ndf, muon_sigma));
					KinematicParticleVertexFitter fitter34;
					RefCountedKinematicTree jpsi2VertexFitTree;
					jpsi2VertexFitTree = fitter34.fit(muonParticles34);
					if (!jpsi2VertexFitTree->isValid())
					{
						continue;
					}
					jpsi2VertexFitTree->movePointerToTheTop();
					RefCountedKinematicParticle jpsi2_vFit_noMC = jpsi2VertexFitTree->currentParticle();
					RefCountedKinematicVertex jpsi2_vFit_vertex_noMC = jpsi2VertexFitTree->currentDecayVertex();
					KinematicParameters mymumupara2 = jpsi2_vFit_noMC->currentState().kinematicParameters();
					double jpsi2_vtxprob = ChiSquaredProbability((double)(jpsi2_vFit_vertex_noMC->chiSquared()), (double)(jpsi2_vFit_vertex_noMC->degreesOfFreedom()));
					if (jpsi2_vFit_noMC->currentState().mass() > 3.4 || jpsi2_vFit_noMC->currentState().mass() < 2.8)
					{
						continue;
					} // change 4.0 to 3.4 in 1208
					if (jpsi2_vtxprob < vtxprobprecut)
					{
						continue;
					}

					// mass constrain for jpsi from X6900:
					KinematicParticleFitter csfitter34;
					jpsi2VertexFitTree = csfitter34.fit(jpsi_cs34, jpsi2VertexFitTree);
					jpsi2VertexFitTree->movePointerToTheTop();
					RefCountedKinematicParticle jpsi2_vFit_cs = jpsi2VertexFitTree->currentParticle();
					RefCountedKinematicVertex jpsi2_vFit_vertex_cs = jpsi2VertexFitTree->currentDecayVertex();
					KinematicParameters mymumupara2_cs = jpsi2_vFit_cs->currentState().kinematicParameters();
					double jpsi2_vtxprob_cs = ChiSquaredProbability((double)(jpsi2_vFit_vertex_cs->chiSquared()), (double)(jpsi2_vFit_vertex_cs->degreesOfFreedom()));

					// actually the code looped the same 4muon particles twice due to different combinations. Let it be
					vector<RefCountedKinematicParticle> X_4muParticles;
					X_4muParticles.push_back(pmumuFactory.particle(muon1TT, muon_mass, chi, ndf, muon_sigma));
					X_4muParticles.push_back(pmumuFactory.particle(muon2TT, muon_mass, chi, ndf, muon_sigma));
					X_4muParticles.push_back(pmumuFactory34.particle(muon3TT, muon_mass, chi, ndf, muon_sigma));
					X_4muParticles.push_back(pmumuFactory34.particle(muon4TT, muon_mass, chi, ndf, muon_sigma));
					KinematicParticleVertexFitter fitterX_4mu;
					RefCountedKinematicTree X_4muVertexFitTree;
					X_4muVertexFitTree = fitterX_4mu.fit(X_4muParticles);
					if (!X_4muVertexFitTree->isValid())
					{
						continue;
					}
					X_4muVertexFitTree->movePointerToTheTop();
					RefCountedKinematicParticle X_4mu_vFit_noMC = X_4muVertexFitTree->currentParticle();
					RefCountedKinematicVertex X_4mu_vFit_vertex_noMC = X_4muVertexFitTree->currentDecayVertex();
					double X_4mu_vtxprob = ChiSquaredProbability((double)(X_4mu_vFit_vertex_noMC->chiSquared()), (double)(X_4mu_vFit_vertex_noMC->degreesOfFreedom()));
					if (X_4mu_vtxprob < vtxprobprecut)
					{
						continue;
					}

					// psi2s
					for (View<pat::PackedCandidate>::const_iterator iTrack1 = theTrackHandle->begin(); // MINIAOD
						 iTrack1 != theTrackHandle->end(); ++iTrack1)
					{
						if (((iTrack1->px() == iMuon2->track()->px() && iTrack1->py() == iMuon2->track()->py() && iTrack1->pz() == iMuon2->track()->pz()) || (iTrack1->px() == iMuon1->track()->px() && iTrack1->py() == iMuon1->track()->py() && iTrack1->pz() == iMuon1->track()->pz())))
						{
							continue;
						}
						if (((iTrack1->px() == iMuon3->track()->px() && iTrack1->py() == iMuon3->track()->py() && iTrack1->pz() == iMuon3->track()->pz()) || (iTrack1->px() == iMuon4->track()->px() && iTrack1->py() == iMuon4->track()->py() && iTrack1->pz() == iMuon4->track()->pz())))
						{
							continue;
						}
						if (iTrack1->pt() < pionptcut)
						{
							continue;
						}

						for (View<pat::PackedCandidate>::const_iterator iTrack2 = iTrack1 + 1; // MINIAOD
							 iTrack2 != theTrackHandle->end(); ++iTrack2)
						{
							if (((iTrack2->px() == iMuon2->track()->px() && iTrack2->py() == iMuon2->track()->py() && iTrack2->pz() == iMuon2->track()->pz()) || (iTrack2->px() == iMuon1->track()->px() && iTrack2->py() == iMuon1->track()->py() && iTrack2->pz() == iMuon1->track()->pz())))
							{
								continue;
							}
							if (((iTrack2->px() == iMuon3->track()->px() && iTrack2->py() == iMuon3->track()->py() && iTrack2->pz() == iMuon3->track()->pz()) || (iTrack2->px() == iMuon4->track()->px() && iTrack2->py() == iMuon4->track()->py() && iTrack2->pz() == iMuon4->track()->pz())))
							{
								continue;
							}
							if (iTrack2->pt() < pionptcut)
							{
								continue;
							}
							if ((iTrack1->charge() + iTrack2->charge()) != 0)
								continue;

							// MINIAOD begin
							if (!iTrack1->hasTrackDetails() || iTrack1->charge() == 0)
							{
								// cout << "iTrack1->hasTrackDetails() = " << iTrack1->hasTrackDetails() << endl;
								// cout << "iTrack1->charge() = " << iTrack1->charge() << endl;
								continue;
							}
							if (!iTrack2->hasTrackDetails() || iTrack2->charge() == 0)
							{
								// cout << "iTrack2->hasTrackDetails() = " << iTrack2->hasTrackDetails() << endl;
								// cout << "iTrack2->charge() = " << iTrack2->charge() << endl;
								continue;
							}
							// MINIAOD end

							TLorentzVector P4_Track1, P4_Track2, P4_jpsipipi;
							P4_Track1.SetPtEtaPhiM(iTrack1->pt(), iTrack1->eta(), iTrack1->phi(), myPimass);
							P4_Track2.SetPtEtaPhiM(iTrack2->pt(), iTrack2->eta(), iTrack2->phi(), myPimass);
							P4_jpsipipi = P4_mu1 + P4_mu2 + P4_Track1 + P4_Track2;

							if (P4_Track1.DeltaR(P4_jpsipipi) > pionDRcut)
							{
								continue;
							}
							if (P4_Track2.DeltaR(P4_jpsipipi) > pionDRcut)
							{
								continue;
							}

							// TransientTrack trackTT1(*iTrack1, &(bFieldHandle));
							// TransientTrack trackTT2(*iTrack2, &(bFieldHandle));
							TransientTrack trackTT1(*(iTrack1->bestTrack()), &(bFieldHandle)); // MINIAOD
							TransientTrack trackTT2(*(iTrack2->bestTrack()), &(bFieldHandle)); // MINIAOD
							KinematicParticleFactoryFromTransientTrack JPiPiFactory;
							// The mass of a muon and the insignificant mass sigma
							// to avoid singularities in the covariance matrix.
							ParticleMass pion_mass = myPimass; // pdg mass
							float pion_sigma = myPimasserr;
							// initial chi2 and ndf before kinematic fits.
							float chi = 0.;
							float ndf = 0.;

							// mass constrain for psi2s from X6900, now noMC
							// first fit to have a psi2s
							vector<RefCountedKinematicParticle> JPiPiParticles;
							JPiPiParticles.push_back(JPiPiFactory.particle(trackTT1, pion_mass, chi, ndf, pion_sigma));
							JPiPiParticles.push_back(JPiPiFactory.particle(trackTT2, pion_mass, chi, ndf, pion_sigma));
							JPiPiParticles.push_back(pmumuFactory.particle(muon1TT, muon_mass, chi, ndf, muon_sigma));
							JPiPiParticles.push_back(pmumuFactory.particle(muon2TT, muon_mass, chi, ndf, muon_sigma));

							KinematicParticleVertexFitter fitter2;
							RefCountedKinematicTree JPiPiVertexFitTree;
							JPiPiVertexFitTree = fitter2.fit(JPiPiParticles);
							if (!(JPiPiVertexFitTree->isValid()))
							{
								continue;
							}
							JPiPiVertexFitTree->movePointerToTheTop();
							RefCountedKinematicParticle JPiPi_vFit_noMC = JPiPiVertexFitTree->currentParticle();
							RefCountedKinematicVertex JPiPi_vFit_vertex_noMC = JPiPiVertexFitTree->currentDecayVertex();

							double Jpsipipivtxprob = ChiSquaredProbability((double)(JPiPi_vFit_vertex_noMC->chiSquared()), (double)(JPiPi_vFit_vertex_noMC->degreesOfFreedom()));
							if (JPiPi_vFit_noMC->currentState().mass() > 4.5)
							{
								continue;
							}
							if (Jpsipipivtxprob < vtxprobprecut)
							{
								continue;
							}

							// fit the 6 track together to a vertex
							vector<RefCountedKinematicParticle> JPiPiJParticles;
							JPiPiJParticles.push_back(JPiPiFactory.particle(trackTT1, pion_mass, chi, ndf, pion_sigma));
							JPiPiJParticles.push_back(JPiPiFactory.particle(trackTT2, pion_mass, chi, ndf, pion_sigma));
							JPiPiJParticles.push_back(pmumuFactory.particle(muon1TT, muon_mass, chi, ndf, muon_sigma));
							JPiPiJParticles.push_back(pmumuFactory.particle(muon2TT, muon_mass, chi, ndf, muon_sigma));
							JPiPiJParticles.push_back(pmumuFactory.particle(muon3TT, muon_mass, chi, ndf, muon_sigma));
							JPiPiJParticles.push_back(pmumuFactory.particle(muon4TT, muon_mass, chi, ndf, muon_sigma));

							KinematicParticleVertexFitter fitter3;
							RefCountedKinematicTree JPiPiJVertexFitTree;
							JPiPiJVertexFitTree = fitter3.fit(JPiPiJParticles);
							if (!(JPiPiJVertexFitTree->isValid()))
							{
								continue;
							}
							JPiPiJVertexFitTree->movePointerToTheTop();
							RefCountedKinematicParticle JPiPiJ_vFit_noMC = JPiPiJVertexFitTree->currentParticle();
							RefCountedKinematicVertex JPiPiJ_vFit_vertex_noMC = JPiPiJVertexFitTree->currentDecayVertex();

							double JpsipipiJpsivtxprob = ChiSquaredProbability((double)(JPiPiJ_vFit_vertex_noMC->chiSquared()), (double)(JPiPiJ_vFit_vertex_noMC->degreesOfFreedom()));

							JPiPiJVertexFitTree->movePointerToTheFirstChild();
							RefCountedKinematicParticle JPiPiJ_pi1 = JPiPiJVertexFitTree->currentParticle();
							JPiPiJVertexFitTree->movePointerToTheNextChild();
							RefCountedKinematicParticle JPiPiJ_pi2 = JPiPiJVertexFitTree->currentParticle();

							KinematicParameters JPiPiJ_pi1_KP = JPiPiJ_pi1->currentState().kinematicParameters();
							KinematicParameters JPiPiJ_pi2_KP = JPiPiJ_pi2->currentState().kinematicParameters();

							X_mass->push_back(JPiPiJ_vFit_noMC->currentState().mass());
							X_Fit_VtxProb->push_back(JpsipipiJpsivtxprob);
							X_Fit_Chi2->push_back((double)(JPiPiJ_vFit_vertex_noMC->chiSquared()));
							X_Fit_ndof->push_back((double)(JPiPiJ_vFit_vertex_noMC->degreesOfFreedom()));

							X_mu1Idx->push_back(std::distance(thePATMuonHandle->begin(), iMuon1));
							X_mu2Idx->push_back(std::distance(thePATMuonHandle->begin(), iMuon2));
							X_mu3Idx->push_back(std::distance(thePATMuonHandle->begin(), iMuon3));
							X_mu4Idx->push_back(std::distance(thePATMuonHandle->begin(), iMuon4));

							X_Jpsi1mass->push_back(jpsi1_vFit_noMC->currentState().mass());
							X_Jpsi1prob->push_back(jpsi1_vtxprob);
							X_Jpsi1px->push_back(mymumupara.momentum().x());
							X_Jpsi1py->push_back(mymumupara.momentum().y());
							X_Jpsi1pz->push_back(mymumupara.momentum().z());
							if (jpsi1_vFit_noMC->currentState().kinematicParametersError().matrix()(6, 6) > 0)
							{
								X_Jpsi1massErr->push_back(sqrt(jpsi1_vFit_noMC->currentState().kinematicParametersError().matrix()(6, 6)));
							}
							else
							{
								X_Jpsi1massErr->push_back(-9);
							}

							X_Jpsi2mass->push_back(jpsi2_vFit_noMC->currentState().mass());
							X_Jpsi2prob->push_back(jpsi2_vtxprob);
							X_Jpsi2px->push_back(mymumupara2.momentum().x());
							X_Jpsi2py->push_back(mymumupara2.momentum().y());
							X_Jpsi2pz->push_back(mymumupara2.momentum().z());
							if (jpsi2_vFit_noMC->currentState().kinematicParametersError().matrix()(6, 6) > 0)
							{
								X_Jpsi2massErr->push_back(sqrt(jpsi2_vFit_noMC->currentState().kinematicParametersError().matrix()(6, 6)));
							}
							else
							{
								X_Jpsi2massErr->push_back(-9);
							}

							X_4mumass->push_back(X_4mu_vFit_noMC->currentState().mass());
							X_4muprob->push_back(X_4mu_vtxprob);
							X_4mupx->push_back(X_4mu_vFit_noMC->currentState().kinematicParameters().momentum().x());
							X_4mupy->push_back(X_4mu_vFit_noMC->currentState().kinematicParameters().momentum().y());
							X_4mupz->push_back(X_4mu_vFit_noMC->currentState().kinematicParameters().momentum().z());
							if (X_4mu_vFit_noMC->currentState().kinematicParametersError().matrix()(6, 6) > 0)
							{
								X_4mumassErr->push_back(sqrt(X_4mu_vFit_noMC->currentState().kinematicParametersError().matrix()(6, 6)));
							}
							else
							{
								X_4mumassErr->push_back(-9);
							}

							X_JPiPi_mass->push_back(JPiPi_vFit_noMC->currentState().mass());
							X_JPiPi_VtxProb->push_back(Jpsipipivtxprob);
							X_JPiPi_px->push_back(JPiPi_vFit_noMC->currentState().kinematicParameters().momentum().x());
							X_JPiPi_py->push_back(JPiPi_vFit_noMC->currentState().kinematicParameters().momentum().y());
							X_JPiPi_pz->push_back(JPiPi_vFit_noMC->currentState().kinematicParameters().momentum().z());
							if (JPiPi_vFit_noMC->currentState().kinematicParametersError().matrix()(6, 6) > 0)
							{
								X_JPiPi_massErr->push_back(sqrt(JPiPi_vFit_noMC->currentState().kinematicParametersError().matrix()(6, 6)));
							}
							else
							{
								X_JPiPi_massErr->push_back(-9);
							}

							X_JPiPi_Pi1px->push_back(JPiPiJ_pi1_KP.momentum().x()); // iTrack1->px()
							X_JPiPi_Pi1py->push_back(JPiPiJ_pi1_KP.momentum().y());
							X_JPiPi_Pi1pz->push_back(JPiPiJ_pi1_KP.momentum().z());

							X_JPiPi_Pi2px->push_back(JPiPiJ_pi2_KP.momentum().x());
							X_JPiPi_Pi2py->push_back(JPiPiJ_pi2_KP.momentum().y());
							X_JPiPi_Pi2pz->push_back(JPiPiJ_pi2_KP.momentum().z());

							// mass constrain for psi2s from X6900 on 1208:
							int flag = 1;
							vector<RefCountedKinematicParticle> cs_JPiPiParticles;
							cs_JPiPiParticles.push_back(JPiPiFactory.particle(trackTT1, pion_mass, chi, ndf, pion_sigma));
							cs_JPiPiParticles.push_back(JPiPiFactory.particle(trackTT2, pion_mass, chi, ndf, pion_sigma));
							cs_JPiPiParticles.push_back(jpsi1_vFit_cs);

							KinematicParticleVertexFitter cs_fitter2;
							RefCountedKinematicTree cs_JPiPiVertexFitTree;
							cs_JPiPiVertexFitTree = cs_fitter2.fit(cs_JPiPiParticles);
							if (!(cs_JPiPiVertexFitTree->isValid()))
							{
								flag = 0;
							}
							RefCountedKinematicParticle cs_JPiPi_vFit_noMC;
							RefCountedKinematicVertex cs_JPiPi_vFit_vertex_noMC;
							double cs_Jpsipipivtxprob;
							if (flag == 1)
							{
								cs_JPiPiVertexFitTree->movePointerToTheTop();
								cs_JPiPi_vFit_noMC = cs_JPiPiVertexFitTree->currentParticle();
								cs_JPiPi_vFit_vertex_noMC = cs_JPiPiVertexFitTree->currentDecayVertex();
								cs_Jpsipipivtxprob = ChiSquaredProbability((double)(cs_JPiPi_vFit_vertex_noMC->chiSquared()), (double)(cs_JPiPi_vFit_vertex_noMC->degreesOfFreedom()));
								if (cs_JPiPi_vFit_noMC->currentState().mass() > 4.5)
								{
									flag = 0;
								}
								if (cs_Jpsipipivtxprob < vtxprobprecut)
								{
									flag = 0;
								}
							} // if(flag==1)
							KinematicParticleFitter csfitterJPiPi;
							RefCountedKinematicParticle JPiPi_vFit_cs;
							RefCountedKinematicVertex JPiPi_vFit_vertex_cs;
							double cs_Jpsipipivtxprob_cs;
							// X block
							vector<RefCountedKinematicParticle> XParticles;
							KinematicParticleVertexFitter fitterX;
							RefCountedKinematicTree XVertexFitTree;
							RefCountedKinematicParticle JPiPi_pi1, JPiPi_pi2;
							KinematicParameters JPiPi_pi1_KP, JPiPi_pi2_KP;
							if (flag == 1)
							{
								cs_JPiPiVertexFitTree = csfitterJPiPi.fit(JPiPi_cs, cs_JPiPiVertexFitTree);
								cs_JPiPiVertexFitTree->movePointerToTheTop();
								JPiPi_vFit_cs = cs_JPiPiVertexFitTree->currentParticle();
								JPiPi_vFit_vertex_cs = cs_JPiPiVertexFitTree->currentDecayVertex();
								cs_Jpsipipivtxprob_cs = ChiSquaredProbability((double)(JPiPi_vFit_vertex_cs->chiSquared()), (double)(JPiPi_vFit_vertex_cs->degreesOfFreedom()));
								// 1222
								cs_JPiPiVertexFitTree->movePointerToTheFirstChild();
								JPiPi_pi1 = cs_JPiPiVertexFitTree->currentParticle();
								cs_JPiPiVertexFitTree->movePointerToTheNextChild();
								JPiPi_pi2 = cs_JPiPiVertexFitTree->currentParticle();
								JPiPi_pi1_KP = JPiPi_pi1->currentState().kinematicParameters();
								JPiPi_pi2_KP = JPiPi_pi2->currentState().kinematicParameters();

								// X block
								XParticles.push_back(jpsi2_vFit_cs);
								XParticles.push_back(JPiPi_vFit_cs);
								XVertexFitTree = fitterX.fit(XParticles);
								if (!(XVertexFitTree->isValid()))
								{
									flag = 0;
								}
							} // if(flag==1)
							if (flag == 1)
							{
								XVertexFitTree->movePointerToTheTop();
								RefCountedKinematicParticle X_vFit_noMC = XVertexFitTree->currentParticle();
								RefCountedKinematicVertex X_vFit_vertex_noMC = XVertexFitTree->currentDecayVertex();
								// cout << "Xmass: " << X_vFit_noMC->currentState().mass() << endl;

								double Xvtxprob = ChiSquaredProbability((double)(X_vFit_vertex_noMC->chiSquared()), (double)(X_vFit_vertex_noMC->degreesOfFreedom()));
								cs_X_mass->push_back(X_vFit_noMC->currentState().mass());
								cs_X_Fit_VtxProb->push_back(Xvtxprob);
								cs_X_Fit_Chi2->push_back((double)(X_vFit_vertex_noMC->chiSquared()));
								cs_X_Fit_ndof->push_back((double)(X_vFit_vertex_noMC->degreesOfFreedom()));

								cs_X_Jpsi1mass->push_back(jpsi1_vFit_cs->currentState().mass());
								cs_X_Jpsi1prob->push_back(jpsi1_vtxprob_cs);
								cs_X_Jpsi1px->push_back(mymumupara_cs.momentum().x());
								cs_X_Jpsi1py->push_back(mymumupara_cs.momentum().y());
								cs_X_Jpsi1pz->push_back(mymumupara_cs.momentum().z());
								if (jpsi1_vFit_cs->currentState().kinematicParametersError().matrix()(6, 6) > 0)
								{
									cs_X_Jpsi1massErr->push_back(sqrt(jpsi1_vFit_cs->currentState().kinematicParametersError().matrix()(6, 6)));
								}
								else
								{
									cs_X_Jpsi1massErr->push_back(-9);
								}

								cs_X_Jpsi2mass->push_back(jpsi2_vFit_cs->currentState().mass());
								cs_X_Jpsi2prob->push_back(jpsi2_vtxprob_cs);
								cs_X_Jpsi2px->push_back(mymumupara2_cs.momentum().x());
								cs_X_Jpsi2py->push_back(mymumupara2_cs.momentum().y());
								cs_X_Jpsi2pz->push_back(mymumupara2_cs.momentum().z());
								if (jpsi2_vFit_cs->currentState().kinematicParametersError().matrix()(6, 6) > 0)
								{
									cs_X_Jpsi2massErr->push_back(sqrt(jpsi2_vFit_cs->currentState().kinematicParametersError().matrix()(6, 6)));
								}
								else
								{
									cs_X_Jpsi2massErr->push_back(-9);
								}

								cs_X_JPiPi_mass->push_back(JPiPi_vFit_cs->currentState().mass());
								cs_X_JPiPi_VtxProb->push_back(cs_Jpsipipivtxprob_cs);
								cs_X_JPiPi_px->push_back(JPiPi_vFit_cs->currentState().kinematicParameters().momentum().x());
								cs_X_JPiPi_py->push_back(JPiPi_vFit_cs->currentState().kinematicParameters().momentum().y());
								cs_X_JPiPi_pz->push_back(JPiPi_vFit_cs->currentState().kinematicParameters().momentum().z());
								if (JPiPi_vFit_cs->currentState().kinematicParametersError().matrix()(6, 6) > 0)
								{
									cs_X_JPiPi_massErr->push_back(sqrt(JPiPi_vFit_cs->currentState().kinematicParametersError().matrix()(6, 6)));
								}
								else
								{
									cs_X_JPiPi_massErr->push_back(-9);
								} // 1222

								cs_X_JPiPinoMC_mass->push_back(cs_JPiPi_vFit_noMC->currentState().mass());
								if (cs_JPiPi_vFit_noMC->currentState().kinematicParametersError().matrix()(6, 6) > 0)
								{
									cs_X_JPiPinoMC_massErr->push_back(sqrt(cs_JPiPi_vFit_noMC->currentState().kinematicParametersError().matrix()(6, 6)));
								}
								else
								{
									cs_X_JPiPinoMC_massErr->push_back(-9);
								}

								cs_X_JPiPi_Pi1px->push_back(JPiPi_pi1_KP.momentum().x());
								cs_X_JPiPi_Pi1py->push_back(JPiPi_pi1_KP.momentum().y());
								cs_X_JPiPi_Pi1pz->push_back(JPiPi_pi1_KP.momentum().z());

								cs_X_JPiPi_Pi2px->push_back(JPiPi_pi2_KP.momentum().x());
								cs_X_JPiPi_Pi2py->push_back(JPiPi_pi2_KP.momentum().y());
								cs_X_JPiPi_Pi2pz->push_back(JPiPi_pi2_KP.momentum().z());
							} // if(flag==1)
							if (flag == 0)
							{
								// cout << "-1" << endl;
								cs_X_mass->push_back(-9);
								cs_X_Fit_VtxProb->push_back(-9);
								cs_X_Fit_Chi2->push_back(-9);
								cs_X_Fit_ndof->push_back(-9);
								// cs_X_mu1Idx->push_back(-9);//no use, the same with no mass constrain, and no -9
								// cs_X_mu2Idx->push_back(-9);
								// cs_X_mu3Idx->push_back(-9);
								// cs_X_mu4Idx->push_back(-9);
								cs_X_Jpsi1mass->push_back(-9);
								cs_X_Jpsi1prob->push_back(-9);
								cs_X_Jpsi1px->push_back(-9999);
								cs_X_Jpsi1py->push_back(-9999);
								cs_X_Jpsi1pz->push_back(-9999);
								cs_X_Jpsi1massErr->push_back(-9999);
								cs_X_Jpsi2mass->push_back(-9);
								cs_X_Jpsi2prob->push_back(-9);
								cs_X_Jpsi2px->push_back(-9999);
								cs_X_Jpsi2py->push_back(-9999);
								cs_X_Jpsi2pz->push_back(-9999);
								cs_X_Jpsi2massErr->push_back(-9999);
								cs_X_JPiPi_mass->push_back(-9);
								cs_X_JPiPi_VtxProb->push_back(-9);
								cs_X_JPiPi_px->push_back(-9999);
								cs_X_JPiPi_py->push_back(-9999);
								cs_X_JPiPi_pz->push_back(-9999);
								cs_X_JPiPi_massErr->push_back(-9999);
								cs_X_JPiPinoMC_mass->push_back(-9999);
								cs_X_JPiPinoMC_massErr->push_back(-9999);
								cs_X_JPiPi_Pi1px->push_back(-9999);
								cs_X_JPiPi_Pi1py->push_back(-9999);
								cs_X_JPiPi_Pi1pz->push_back(-9999);
								cs_X_JPiPi_Pi2px->push_back(-9999);
								cs_X_JPiPi_Pi2py->push_back(-9999);
								cs_X_JPiPi_Pi2pz->push_back(-9999);
							} // if(flag==0)
							  // 1208over
						}	  // itrack2
					}		  // itrack1
				}			  // mu4_loop
			}				  // mu3_loop
		}					  // for (std::vector < pat::Muon >::const_iterator iMuon2 = iMuon1 + 1;
	}						  // for (std::vector < pat::Muon >::const_iterator iMuon1 = thePATMuonHandle->begin();

	if (X_Fit_VtxProb->size() > 0 || doMC)
	{
		X_One_Tree_->Fill();
	}

	if (Debug_)
	{
	}

	if (doMC)
	{
		MC_X_px->clear();
		MC_X_py->clear();
		MC_X_pz->clear();
		MC_X_mass->clear();
		MC_X_chg->clear();
		MC_Dau_jpsipdgId->clear();
		MC_Dau_jpsipx->clear();
		MC_Dau_jpsipy->clear();
		MC_Dau_jpsipz->clear();
		MC_Dau_jpsimass->clear();
		MC_Dau_psi2spdgId->clear();
		MC_Dau_psi2spx->clear();
		MC_Dau_psi2spy->clear();
		MC_Dau_psi2spz->clear();
		MC_Dau_psi2smass->clear();
		MC_Granddau_mu1pdgId->clear();
		MC_Granddau_mu1px->clear();
		MC_Granddau_mu1py->clear();
		MC_Granddau_mu1pz->clear();
		MC_Granddau_mu2pdgId->clear();
		MC_Granddau_mu2px->clear();
		MC_Granddau_mu2py->clear();
		MC_Granddau_mu2pz->clear();
		MC_Granddau_jpsipdgId->clear();
		MC_Granddau_jpsipx->clear();
		MC_Granddau_jpsipy->clear();
		MC_Granddau_jpsipz->clear();
		MC_Granddau_jpsimass->clear();
		MC_Granddau_pi1pdgId->clear();
		MC_Granddau_pi1px->clear();
		MC_Granddau_pi1py->clear();
		MC_Granddau_pi1pz->clear();
		MC_Granddau_pi2pdgId->clear();
		MC_Granddau_pi2px->clear();
		MC_Granddau_pi2py->clear();
		MC_Granddau_pi2pz->clear();
		MC_Grandgranddau_mu3pdgId->clear();
		MC_Grandgranddau_mu3px->clear();
		MC_Grandgranddau_mu3py->clear();
		MC_Grandgranddau_mu3pz->clear();
		MC_Grandgranddau_mu4pdgId->clear();
		MC_Grandgranddau_mu4px->clear();
		MC_Grandgranddau_mu4py->clear();
		MC_Grandgranddau_mu4pz->clear();

		Match_mu1px->clear();
		Match_mu1py->clear();
		Match_mu1pz->clear();
		Match_mu2px->clear();
		Match_mu2py->clear();
		Match_mu2pz->clear();
		Match_mu3px->clear();
		Match_mu3py->clear();
		Match_mu3pz->clear();
		Match_mu4px->clear();
		Match_mu4py->clear();
		Match_mu4pz->clear();

		Match_pi1px->clear();
		Match_pi1py->clear();
		Match_pi1pz->clear();
		Match_pi2px->clear();
		Match_pi2py->clear();
		Match_pi2pz->clear();
	}

	trigRes->clear();
	trigNames->clear();
	L1TT->clear();
	MatchTriggerNames->clear();
	muIsJpsiTrigMatch->clear();
	muIsUpsTrigMatch->clear();
	runNum = 0;
	evtNum = 0;
	lumiNum = 0;
	nGoodPrimVtx = 0;
	priVtxX = 0;
	priVtxY = 0;
	priVtxZ = 0;
	priVtxXE = 0;
	priVtxYE = 0;
	priVtxZE = 0;
	priVtxChiNorm = 0;
	priVtxChi = 0;
	priVtxCL = 0;
	mybxlumicorr = 0;
	myrawbxlumi = 0;
	PriVtxXCorrX->clear();
	PriVtxXCorrY->clear();
	PriVtxXCorrZ->clear();
	PriVtxXCorrEX->clear();
	PriVtxXCorrEY->clear();
	PriVtxXCorrEZ->clear();
	PriVtxXCorrC2->clear();
	PriVtxXCorrCL->clear();

	nMu = 0;
	muPx->clear();
	muPy->clear();
	muPz->clear();
	muD0->clear();
	muD0E->clear();
	muDz->clear();
	muChi2->clear();
	muGlChi2->clear();
	mufHits->clear();
	muFirstBarrel->clear();
	muFirstEndCap->clear();
	muDzVtx->clear();
	muDxyVtx->clear();
	muNDF->clear();
	muGlNDF->clear();
	muPhits->clear();
	muShits->clear();
	muGlMuHits->clear();
	muType->clear();
	muQual->clear();
	muTrack->clear();
	muCharge->clear();
	muIsoratio->clear();
	muIsGoodLooseMuon->clear();
	muIsGoodLooseMuonNew->clear();
	muIsGoodSoftMuonNewIlse->clear();
	muIsGoodSoftMuonNewIlseMod->clear();
	muIsGoodTightMuon->clear();
	munMatchedSeg->clear();
	muMVAMuonID->clear();
	musegmentCompatibility->clear();

	mupulldXdZ_pos_noArb->clear();
	mupulldYdZ_pos_noArb->clear();
	mupulldXdZ_pos_ArbDef->clear();
	mupulldYdZ_pos_ArbDef->clear();
	mupulldXdZ_pos_ArbST->clear();
	mupulldYdZ_pos_ArbST->clear();
	mupulldXdZ_pos_noArb_any->clear();
	mupulldYdZ_pos_noArb_any->clear();

	muIsPatLooseMuon->clear();
	muIsPatTightMuon->clear();
	muIsPatSoftMuon->clear();
	muIsPatMediumMuon->clear();
	muUpsVrtxMatch->clear();
	muL3TriggerMatch->clear();

	Jpsi1_nmumuonly = 0;
	Jpsi1_mumuonlyMass->clear();
	Jpsi1_mumuonlyMassErr->clear();
	Jpsi1_mumuonlyVtxProb->clear();
	Jpsi1_mumuonlyFit_Chi2->clear();
	Jpsi1_mumuonlyFit_ndof->clear();

	Jpsi1_mumuonlyPx->clear();
	Jpsi1_mumuonlyPy->clear();
	Jpsi1_mumuonlyPz->clear();
	Jpsi1_mumuonlymu1Idx->clear();
	Jpsi1_mumuonlymu2Idx->clear();

	X_mu1Idx->clear();
	X_mu2Idx->clear();
	X_mu3Idx->clear();
	X_mu4Idx->clear();
	X_mass->clear();
	X_Fit_VtxProb->clear();
	X_Fit_Chi2->clear();
	X_Fit_ndof->clear();
	X_JPiPi_mass->clear();
	X_JPiPi_VtxProb->clear();
	X_JPiPi_px->clear();
	X_JPiPi_py->clear();
	X_JPiPi_pz->clear();
	X_JPiPi_massErr->clear();
	X_px->clear();
	X_py->clear();
	X_pz->clear();
	X_JPiPimass->clear();
	X_Jpsi1mass->clear();
	X_Jpsi1prob->clear();
	X_Jpsi2mass->clear();
	X_Jpsi2prob->clear();
	X_Jpsi2px->clear();
	X_Jpsi2py->clear();
	X_Jpsi2pz->clear();
	X_Jpsi2massErr->clear();
	X_4mumass->clear();
	X_4muprob->clear();
	X_4mupx->clear();
	X_4mupy->clear();
	X_4mupz->clear();
	X_4mumassErr->clear();
	X_JPiPi_Pi1px->clear();
	X_JPiPi_Pi1py->clear();
	X_JPiPi_Pi1pz->clear();
	X_JPiPi_Pi2px->clear();
	X_JPiPi_Pi2py->clear();
	X_JPiPi_Pi2pz->clear();
	X_JPiPi_Pi1DeltaR->clear();
	X_JPiPi_Pi2DeltaR->clear();
	X_JPiPi_Pi1pt->clear();
	X_JPiPi_Pi2pt->clear();
	X_Jpsi1px->clear();
	X_Jpsi1py->clear();
	X_Jpsi1pz->clear();
	X_Jpsi1massErr->clear();
	X_mu12cs_px->clear();
	X_mu12cs_py->clear();
	X_mu12cs_pz->clear();
	X_mu1px->clear();
	X_mu1py->clear();
	X_mu1pz->clear();
	X_mu2px->clear();
	X_mu2py->clear();
	X_mu2pz->clear();
	X_JPiPics_px->clear();
	X_JPiPics_py->clear();
	X_JPiPics_pz->clear();
	X_mu34cs_px->clear();
	X_mu34cs_py->clear();
	X_mu34cs_pz->clear();
	Jpsi1_mumuonlyctau->clear();
	Jpsi1_mumuonlyctauerr->clear();
	mumuonlymuoverlapped->clear();
	Jpsi1_mumuonlyChg->clear();

	// mass constrain variables
	cs_X_mass->clear();
	cs_X_Fit_VtxProb->clear();
	cs_X_Fit_Chi2->clear();
	cs_X_Fit_ndof->clear();
	cs_X_mu1Idx->clear();
	cs_X_mu2Idx->clear();
	cs_X_mu3Idx->clear();
	cs_X_mu4Idx->clear();
	cs_X_Jpsi1mass->clear();
	cs_X_Jpsi1prob->clear();
	cs_X_Jpsi1px->clear();
	cs_X_Jpsi1py->clear();
	cs_X_Jpsi1pz->clear();
	cs_X_Jpsi1massErr->clear();
	cs_X_Jpsi2mass->clear();
	cs_X_Jpsi2prob->clear();
	cs_X_Jpsi2px->clear();
	cs_X_Jpsi2py->clear();
	cs_X_Jpsi2pz->clear();
	cs_X_Jpsi2massErr->clear();
	cs_X_JPiPi_mass->clear();
	cs_X_JPiPi_VtxProb->clear();
	cs_X_JPiPi_px->clear();
	cs_X_JPiPi_py->clear();
	cs_X_JPiPi_pz->clear();
	cs_X_JPiPi_massErr->clear();
	cs_X_JPiPinoMC_mass->clear();
	cs_X_JPiPinoMC_massErr->clear();
	cs_X_JPiPi_Pi1px->clear();
	cs_X_JPiPi_Pi1py->clear();
	cs_X_JPiPi_Pi1pz->clear();
	cs_X_JPiPi_Pi2px->clear();
	cs_X_JPiPi_Pi2py->clear();
	cs_X_JPiPi_Pi2pz->clear();

} // analyze

// ------------ method called once each job just before starting event loop  ------------
void MultiLepPAT::beginRun(edm::Run const &iRun, edm::EventSetup const &iSetup)
{
	// bool changed = true;
	// proccessName_="HLT";
	// hltConfig_.init(iRun,iSetup,proccessName_,changed);
}

void MultiLepPAT::beginJob()
{
	edm::Service<TFileService> fs;

	// estree_ = fs->make<TTree>("eventSummary", "General Event Summary");
	X_One_Tree_ = fs->make<TTree>("X_data", "X(3872) Data");

	X_One_Tree_->Branch("TrigRes", &trigRes);
	X_One_Tree_->Branch("TrigNames", &trigNames);
	X_One_Tree_->Branch("MatchTriggerNames", &MatchTriggerNames);
	X_One_Tree_->Branch("L1TrigRes", &L1TT);

	X_One_Tree_->Branch("evtNum", &evtNum, "evtNum/i");
	X_One_Tree_->Branch("runNum", &runNum, "runNum/i");
	X_One_Tree_->Branch("lumiNum", &lumiNum, "lumiNum/i");
	X_One_Tree_->Branch("nGoodPrimVtx", &nGoodPrimVtx, "nGoodPrimVtx/i");

	// inst. lumi is here
	X_One_Tree_->Branch("mybxlumicorr", &mybxlumicorr, "mybxlumicorr/f");
	X_One_Tree_->Branch("myrawbxlumi", &myrawbxlumi, "myrawbxlumi/f");

	X_One_Tree_->Branch("priVtxX", &priVtxX, "priVtxX/f");
	X_One_Tree_->Branch("priVtxY", &priVtxY, "priVtxY/f");
	X_One_Tree_->Branch("priVtxZ", &priVtxZ, "priVtxZ/f");
	X_One_Tree_->Branch("priVtxXE", &priVtxXE, "priVtxXE/f");
	X_One_Tree_->Branch("priVtxYE", &priVtxYE, "priVtxYE/f");
	X_One_Tree_->Branch("priVtxZE", &priVtxZE, "priVtxZE/f");
	X_One_Tree_->Branch("priVtxChiNorm", &priVtxChiNorm, "priVtxChiNorm/f");
	X_One_Tree_->Branch("priVtxChi", &priVtxChi, "priVtxChi/f");
	X_One_Tree_->Branch("priVtxCL", &priVtxCL, "priVtxCL/f");

	X_One_Tree_->Branch("PriVtxXCorrX", &PriVtxXCorrX);
	X_One_Tree_->Branch("PriVtxXCorrY", &PriVtxXCorrY);
	X_One_Tree_->Branch("PriVtxXCorrZ", &PriVtxXCorrZ);
	X_One_Tree_->Branch("PriVtxXCorrEX", &PriVtxXCorrEX);
	X_One_Tree_->Branch("PriVtxXCorrEY", &PriVtxXCorrEY);
	X_One_Tree_->Branch("PriVtxXCorrEZ", &PriVtxXCorrEZ);
	X_One_Tree_->Branch("PriVtxXCorrC2", &PriVtxXCorrC2);
	X_One_Tree_->Branch("PriVtxXCorrCL", &PriVtxXCorrCL);

	X_One_Tree_->Branch("nMu", &nMu, "nMu/i");
	X_One_Tree_->Branch("muPx", &muPx);
	X_One_Tree_->Branch("muPy", &muPy);
	X_One_Tree_->Branch("muPz", &muPz);
	X_One_Tree_->Branch("muD0", &muD0);
	X_One_Tree_->Branch("muD0E", &muD0E);
	X_One_Tree_->Branch("muDz", &muDz);
	X_One_Tree_->Branch("muChi2", &muChi2);
	X_One_Tree_->Branch("muGlChi2", &muGlChi2);
	X_One_Tree_->Branch("mufHits", &mufHits);
	X_One_Tree_->Branch("muFirstBarrel", &muFirstBarrel);
	X_One_Tree_->Branch("muFirstEndCap", &muFirstEndCap);
	X_One_Tree_->Branch("muDzVtx", &muDzVtx);
	X_One_Tree_->Branch("muDxyVtx", &muDxyVtx);
	X_One_Tree_->Branch("muNDF", &muNDF);
	X_One_Tree_->Branch("muGlNDF", &muGlNDF);
	X_One_Tree_->Branch("muPhits", &muPhits);
	X_One_Tree_->Branch("muShits", &muShits);
	X_One_Tree_->Branch("muGlMuHits", &muGlMuHits);
	X_One_Tree_->Branch("muType", &muType);
	X_One_Tree_->Branch("muQual", &muQual);
	X_One_Tree_->Branch("muTrack", &muTrack);
	X_One_Tree_->Branch("muCharge", &muCharge);
	X_One_Tree_->Branch("muIsoratio", &muIsoratio);
	X_One_Tree_->Branch("munMatchedSeg", &munMatchedSeg);
	X_One_Tree_->Branch("muIsGoodSoftMuonNewIlseMod", &muIsGoodSoftMuonNewIlseMod);
	X_One_Tree_->Branch("muIsGoodSoftMuonNewIlse", &muIsGoodSoftMuonNewIlse);
	X_One_Tree_->Branch("muIsGoodLooseMuonNew", &muIsGoodLooseMuonNew);
	X_One_Tree_->Branch("muIsGoodLooseMuon", &muIsGoodLooseMuon);
	X_One_Tree_->Branch("muIsGoodTightMuon", &muIsGoodTightMuon);

	X_One_Tree_->Branch("muIsPatLooseMuon", &muIsPatLooseMuon);
	X_One_Tree_->Branch("muIsPatTightMuon", &muIsPatTightMuon);
	X_One_Tree_->Branch("muIsPatSoftMuon", &muIsPatSoftMuon);
	X_One_Tree_->Branch("muIsPatMediumMuon", &muIsPatMediumMuon);

	X_One_Tree_->Branch("muIsJpsiTrigMatch", &muIsJpsiTrigMatch);
	X_One_Tree_->Branch("muIsUpsTrigMatch", &muIsUpsTrigMatch);
	X_One_Tree_->Branch("muMVAMuonID", &muMVAMuonID);
	X_One_Tree_->Branch("musegmentCompatibility", &musegmentCompatibility);

	X_One_Tree_->Branch("mupulldXdZ_pos_noArb", &mupulldXdZ_pos_noArb);
	X_One_Tree_->Branch("mupulldYdZ_pos_noArb", &mupulldYdZ_pos_noArb);
	X_One_Tree_->Branch("mupulldXdZ_pos_ArbDef", &mupulldXdZ_pos_ArbDef);
	X_One_Tree_->Branch("mupulldYdZ_pos_ArbDef", &mupulldYdZ_pos_ArbDef);
	X_One_Tree_->Branch("mupulldXdZ_pos_ArbST", &mupulldXdZ_pos_ArbST);
	X_One_Tree_->Branch("mupulldYdZ_pos_ArbST", &mupulldYdZ_pos_ArbST);
	X_One_Tree_->Branch("mupulldXdZ_pos_noArb_any", &mupulldXdZ_pos_noArb_any);
	X_One_Tree_->Branch("mupulldYdZ_pos_noArb_any", &mupulldYdZ_pos_noArb_any);

	X_One_Tree_->Branch("muUpsVrtxMatch", &muUpsVrtxMatch);
	X_One_Tree_->Branch("muL3TriggerMatch", &muL3TriggerMatch);

	X_One_Tree_->Branch("Jpsi1_nmumuonly", &Jpsi1_nmumuonly);
	X_One_Tree_->Branch("Jpsi1_mumuonlyMass", &Jpsi1_mumuonlyMass);
	X_One_Tree_->Branch("Jpsi1_mumuonlyMassErr", &Jpsi1_mumuonlyMassErr);
	X_One_Tree_->Branch("Jpsi1_mumuonlyVtxProb", &Jpsi1_mumuonlyVtxProb);
	X_One_Tree_->Branch("Jpsi1_mumuonlyFit_Chi2", &Jpsi1_mumuonlyFit_Chi2);
	X_One_Tree_->Branch("Jpsi1_mumuonlyFit_ndof", &Jpsi1_mumuonlyFit_ndof);
	X_One_Tree_->Branch("Jpsi1_mumuonlyPx", &Jpsi1_mumuonlyPx);
	X_One_Tree_->Branch("Jpsi1_mumuonlyPy", &Jpsi1_mumuonlyPy);
	X_One_Tree_->Branch("Jpsi1_mumuonlyPz", &Jpsi1_mumuonlyPz);
	X_One_Tree_->Branch("Jpsi1_mumuonlymu1Idx", &Jpsi1_mumuonlymu1Idx);
	X_One_Tree_->Branch("Jpsi1_mumuonlymu2Idx", &Jpsi1_mumuonlymu2Idx);
	X_One_Tree_->Branch("X_mu1Idx", &X_mu1Idx);
	X_One_Tree_->Branch("X_mu2Idx", &X_mu2Idx);
	X_One_Tree_->Branch("X_mu3Idx", &X_mu3Idx);
	X_One_Tree_->Branch("X_mu4Idx", &X_mu4Idx);
	X_One_Tree_->Branch("X_mass", &X_mass);
	X_One_Tree_->Branch("X_Fit_VtxProb", &X_Fit_VtxProb);
	X_One_Tree_->Branch("X_Fit_Chi2", &X_Fit_Chi2);
	X_One_Tree_->Branch("X_Fit_ndof", &X_Fit_ndof);
	X_One_Tree_->Branch("X_JPiPi_mass", &X_JPiPi_mass);
	X_One_Tree_->Branch("X_JPiPi_VtxProb", &X_JPiPi_VtxProb);
	X_One_Tree_->Branch("X_JPiPi_px", &X_JPiPi_px);
	X_One_Tree_->Branch("X_JPiPi_py", &X_JPiPi_py);
	X_One_Tree_->Branch("X_JPiPi_pz", &X_JPiPi_pz);
	X_One_Tree_->Branch("X_JPiPi_massErr", &X_JPiPi_massErr);
	X_One_Tree_->Branch("X_px", &X_px);
	X_One_Tree_->Branch("X_py", &X_py);
	X_One_Tree_->Branch("X_pz", &X_pz);
	X_One_Tree_->Branch("X_Jpsi1mass", &X_Jpsi1mass);
	X_One_Tree_->Branch("X_Jpsi1prob", &X_Jpsi1prob);
	X_One_Tree_->Branch("X_JPiPimass", &X_JPiPimass);
	X_One_Tree_->Branch("X_Jpsi2mass", &X_Jpsi2mass);
	X_One_Tree_->Branch("X_Jpsi2prob", &X_Jpsi2prob);
	X_One_Tree_->Branch("X_Jpsi2px", &X_Jpsi2px);
	X_One_Tree_->Branch("X_Jpsi2py", &X_Jpsi2py);
	X_One_Tree_->Branch("X_Jpsi2pz", &X_Jpsi2pz);
	X_One_Tree_->Branch("X_Jpsi2massErr", &X_Jpsi2massErr);
	X_One_Tree_->Branch("X_4mumass", &X_4mumass);
	X_One_Tree_->Branch("X_4muprob", &X_4muprob);
	X_One_Tree_->Branch("X_4mupx", &X_4mupx);
	X_One_Tree_->Branch("X_4mupy", &X_4mupy);
	X_One_Tree_->Branch("X_4mupz", &X_4mupz);
	X_One_Tree_->Branch("X_4mumassErr", &X_4mumassErr);
	X_One_Tree_->Branch("X_JPiPi_Pi1px", &X_JPiPi_Pi1px);
	X_One_Tree_->Branch("X_JPiPi_Pi1py", &X_JPiPi_Pi1py);
	X_One_Tree_->Branch("X_JPiPi_Pi1pz", &X_JPiPi_Pi1pz);
	X_One_Tree_->Branch("X_JPiPi_Pi2px", &X_JPiPi_Pi2px);
	X_One_Tree_->Branch("X_JPiPi_Pi2py", &X_JPiPi_Pi2py);
	X_One_Tree_->Branch("X_JPiPi_Pi2pz", &X_JPiPi_Pi2pz);
	X_One_Tree_->Branch("X_JPiPi_Pi1DeltaR", &X_JPiPi_Pi1DeltaR);
	X_One_Tree_->Branch("X_JPiPi_Pi2DeltaR", &X_JPiPi_Pi2DeltaR);
	X_One_Tree_->Branch("X_JPiPi_Pi1pt", &X_JPiPi_Pi1pt);
	X_One_Tree_->Branch("X_JPiPi_Pi2pt", &X_JPiPi_Pi2pt);
	X_One_Tree_->Branch("X_Jpsi1px", &X_Jpsi1px);
	X_One_Tree_->Branch("X_Jpsi1py", &X_Jpsi1py);
	X_One_Tree_->Branch("X_Jpsi1pz", &X_Jpsi1pz);
	X_One_Tree_->Branch("X_Jpsi1massErr", &X_Jpsi1massErr);
	X_One_Tree_->Branch("X_mu12cs_px", &X_mu12cs_px);
	X_One_Tree_->Branch("X_mu12cs_py", &X_mu12cs_py);
	X_One_Tree_->Branch("X_mu12cs_pz", &X_mu12cs_pz);
	X_One_Tree_->Branch("X_mu1px", &X_mu1px);
	X_One_Tree_->Branch("X_mu1py", &X_mu1py);
	X_One_Tree_->Branch("X_mu1pz", &X_mu1pz);
	X_One_Tree_->Branch("X_mu2px", &X_mu2px);
	X_One_Tree_->Branch("X_mu2py", &X_mu2py);
	X_One_Tree_->Branch("X_mu2pz", &X_mu2pz);
	X_One_Tree_->Branch("X_JPiPics_px", &X_JPiPics_px);
	X_One_Tree_->Branch("X_JPiPics_py", &X_JPiPics_py);
	X_One_Tree_->Branch("X_JPiPics_pz", &X_JPiPics_pz);
	X_One_Tree_->Branch("X_mu34cs_px", &X_mu34cs_px);
	X_One_Tree_->Branch("X_mu34cs_py", &X_mu34cs_py);
	X_One_Tree_->Branch("X_mu34cs_pz", &X_mu34cs_pz);
	X_One_Tree_->Branch("Jpsi1_mumuonlyctau", &Jpsi1_mumuonlyctau);
	X_One_Tree_->Branch("Jpsi1_mumuonlyctauerr", &Jpsi1_mumuonlyctauerr);
	X_One_Tree_->Branch("mumuonlymuoverlapped", &mumuonlymuoverlapped);
	X_One_Tree_->Branch("Jpsi1_mumuonlyChg", &Jpsi1_mumuonlyChg);

	// mass constrain variables
	X_One_Tree_->Branch("cs_X_mass", &cs_X_mass);
	X_One_Tree_->Branch("cs_X_Fit_VtxProb", &cs_X_Fit_VtxProb);
	X_One_Tree_->Branch("cs_X_Fit_Chi2", &cs_X_Fit_Chi2);
	X_One_Tree_->Branch("cs_X_Fit_ndof", &cs_X_Fit_ndof);
	X_One_Tree_->Branch("cs_X_mu1Idx", &cs_X_mu1Idx);
	X_One_Tree_->Branch("cs_X_mu2Idx", &cs_X_mu2Idx);
	X_One_Tree_->Branch("cs_X_mu3Idx", &cs_X_mu3Idx);
	X_One_Tree_->Branch("cs_X_mu4Idx", &cs_X_mu4Idx);
	X_One_Tree_->Branch("cs_X_Jpsi1mass", &cs_X_Jpsi1mass);
	X_One_Tree_->Branch("cs_X_Jpsi1prob", &cs_X_Jpsi1prob);
	X_One_Tree_->Branch("cs_X_Jpsi1px", &cs_X_Jpsi1px);
	X_One_Tree_->Branch("cs_X_Jpsi1py", &cs_X_Jpsi1py);
	X_One_Tree_->Branch("cs_X_Jpsi1pz", &cs_X_Jpsi1pz);
	X_One_Tree_->Branch("cs_X_Jpsi1massErr", &cs_X_Jpsi1massErr);
	X_One_Tree_->Branch("cs_X_Jpsi2mass", &cs_X_Jpsi2mass);
	X_One_Tree_->Branch("cs_X_Jpsi2prob", &cs_X_Jpsi2prob);
	X_One_Tree_->Branch("cs_X_Jpsi2px", &cs_X_Jpsi2px);
	X_One_Tree_->Branch("cs_X_Jpsi2py", &cs_X_Jpsi2py);
	X_One_Tree_->Branch("cs_X_Jpsi2pz", &cs_X_Jpsi2pz);
	X_One_Tree_->Branch("cs_X_Jpsi2massErr", &cs_X_Jpsi2massErr);
	X_One_Tree_->Branch("cs_X_JPiPi_mass", &cs_X_JPiPi_mass);
	X_One_Tree_->Branch("cs_X_JPiPi_VtxProb", &cs_X_JPiPi_VtxProb);
	X_One_Tree_->Branch("cs_X_JPiPi_px", &cs_X_JPiPi_px);
	X_One_Tree_->Branch("cs_X_JPiPi_py", &cs_X_JPiPi_py);
	X_One_Tree_->Branch("cs_X_JPiPi_pz", &cs_X_JPiPi_pz);
	X_One_Tree_->Branch("cs_X_JPiPi_massErr", &cs_X_JPiPi_massErr);
	X_One_Tree_->Branch("cs_X_JPiPinoMC_mass", &cs_X_JPiPinoMC_mass);
	X_One_Tree_->Branch("cs_X_JPiPinoMC_massErr", &cs_X_JPiPinoMC_massErr);
	X_One_Tree_->Branch("cs_X_JPiPi_Pi1px", &cs_X_JPiPi_Pi1px);
	X_One_Tree_->Branch("cs_X_JPiPi_Pi1py", &cs_X_JPiPi_Pi1py);
	X_One_Tree_->Branch("cs_X_JPiPi_Pi1pz", &cs_X_JPiPi_Pi1pz);
	X_One_Tree_->Branch("cs_X_JPiPi_Pi2px", &cs_X_JPiPi_Pi2px);
	X_One_Tree_->Branch("cs_X_JPiPi_Pi2py", &cs_X_JPiPi_Pi2py);
	X_One_Tree_->Branch("cs_X_JPiPi_Pi2pz", &cs_X_JPiPi_Pi2pz);

	if (doMC)
	{
		X_One_Tree_->Branch("MC_X_px", &MC_X_px);
		X_One_Tree_->Branch("MC_X_py", &MC_X_py);
		X_One_Tree_->Branch("MC_X_pz", &MC_X_pz);
		X_One_Tree_->Branch("MC_X_mass", &MC_X_mass);
		X_One_Tree_->Branch("MC_X_chg", &MC_X_chg);
		X_One_Tree_->Branch("MC_Dau_jpsipdgId", &MC_Dau_jpsipdgId);
		X_One_Tree_->Branch("MC_Dau_jpsipx", &MC_Dau_jpsipx);
		X_One_Tree_->Branch("MC_Dau_jpsipy", &MC_Dau_jpsipy);
		X_One_Tree_->Branch("MC_Dau_jpsipz", &MC_Dau_jpsipz);
		X_One_Tree_->Branch("MC_Dau_jpsimass", &MC_Dau_jpsimass);
		X_One_Tree_->Branch("MC_Dau_psi2spdgId", &MC_Dau_psi2spdgId);
		X_One_Tree_->Branch("MC_Dau_psi2spx", &MC_Dau_psi2spx);
		X_One_Tree_->Branch("MC_Dau_psi2spy", &MC_Dau_psi2spy);
		X_One_Tree_->Branch("MC_Dau_psi2spz", &MC_Dau_psi2spz);
		X_One_Tree_->Branch("MC_Dau_psi2smass", &MC_Dau_psi2smass);
		X_One_Tree_->Branch("MC_Granddau_mu1pdgId", &MC_Granddau_mu1pdgId);
		X_One_Tree_->Branch("MC_Granddau_mu1px", &MC_Granddau_mu1px);
		X_One_Tree_->Branch("MC_Granddau_mu1py", &MC_Granddau_mu1py);
		X_One_Tree_->Branch("MC_Granddau_mu1pz", &MC_Granddau_mu1pz);
		X_One_Tree_->Branch("MC_Granddau_mu2pdgId", &MC_Granddau_mu2pdgId);
		X_One_Tree_->Branch("MC_Granddau_mu2px", &MC_Granddau_mu2px);
		X_One_Tree_->Branch("MC_Granddau_mu2py", &MC_Granddau_mu2py);
		X_One_Tree_->Branch("MC_Granddau_mu2pz", &MC_Granddau_mu2pz);
		X_One_Tree_->Branch("MC_Granddau_jpsipdgId", &MC_Granddau_jpsipdgId);
		X_One_Tree_->Branch("MC_Granddau_jpsipx", &MC_Granddau_jpsipx);
		X_One_Tree_->Branch("MC_Granddau_jpsipy", &MC_Granddau_jpsipy);
		X_One_Tree_->Branch("MC_Granddau_jpsipz", &MC_Granddau_jpsipz);
		X_One_Tree_->Branch("MC_Granddau_jpsimass", &MC_Granddau_jpsimass);
		X_One_Tree_->Branch("MC_Granddau_pi1pdgId", &MC_Granddau_pi1pdgId);
		X_One_Tree_->Branch("MC_Granddau_pi1px", &MC_Granddau_pi1px);
		X_One_Tree_->Branch("MC_Granddau_pi1py", &MC_Granddau_pi1py);
		X_One_Tree_->Branch("MC_Granddau_pi1pz", &MC_Granddau_pi1pz);
		X_One_Tree_->Branch("MC_Granddau_pi2pdgId", &MC_Granddau_pi2pdgId);
		X_One_Tree_->Branch("MC_Granddau_pi2px", &MC_Granddau_pi2px);
		X_One_Tree_->Branch("MC_Granddau_pi2py", &MC_Granddau_pi2py);
		X_One_Tree_->Branch("MC_Granddau_pi2pz", &MC_Granddau_pi2pz);
		X_One_Tree_->Branch("MC_Grandgranddau_mu3pdgId", &MC_Grandgranddau_mu3pdgId);
		X_One_Tree_->Branch("MC_Grandgranddau_mu3px", &MC_Grandgranddau_mu3px);
		X_One_Tree_->Branch("MC_Grandgranddau_mu3py", &MC_Grandgranddau_mu3py);
		X_One_Tree_->Branch("MC_Grandgranddau_mu3pz", &MC_Grandgranddau_mu3pz);
		X_One_Tree_->Branch("MC_Grandgranddau_mu4pdgId", &MC_Grandgranddau_mu4pdgId);
		X_One_Tree_->Branch("MC_Grandgranddau_mu4px", &MC_Grandgranddau_mu4px);
		X_One_Tree_->Branch("MC_Grandgranddau_mu4py", &MC_Grandgranddau_mu4py);
		X_One_Tree_->Branch("MC_Grandgranddau_mu4pz", &MC_Grandgranddau_mu4pz);
		X_One_Tree_->Branch("Match_mu1px", &Match_mu1px);
		X_One_Tree_->Branch("Match_mu1py", &Match_mu1py);
		X_One_Tree_->Branch("Match_mu1pz", &Match_mu1pz);
		X_One_Tree_->Branch("Match_mu2px", &Match_mu2px);
		X_One_Tree_->Branch("Match_mu2py", &Match_mu2py);
		X_One_Tree_->Branch("Match_mu2pz", &Match_mu2pz);
		X_One_Tree_->Branch("Match_mu3px", &Match_mu3px);
		X_One_Tree_->Branch("Match_mu3py", &Match_mu3py);
		X_One_Tree_->Branch("Match_mu3pz", &Match_mu3pz);
		X_One_Tree_->Branch("Match_mu4px", &Match_mu4px);
		X_One_Tree_->Branch("Match_mu4py", &Match_mu4py);
		X_One_Tree_->Branch("Match_mu4pz", &Match_mu4pz);

		X_One_Tree_->Branch("Match_pi1px", &Match_pi1px);
		X_One_Tree_->Branch("Match_pi1py", &Match_pi1py);
		X_One_Tree_->Branch("Match_pi1pz", &Match_pi1pz);
		X_One_Tree_->Branch("Match_pi2px", &Match_pi2px);
		X_One_Tree_->Branch("Match_pi2py", &Match_pi2py);
		X_One_Tree_->Branch("Match_pi2pz", &Match_pi2pz);
	} // if(doMC)

} // begin Job

// ------------ method called once each job just after ending the event loop  ------------
void MultiLepPAT::endJob()
{
	X_One_Tree_->GetDirectory()->cd();
	X_One_Tree_->Write();
}

// define this as a plug-in
DEFINE_FWK_MODULE(MultiLepPAT);
