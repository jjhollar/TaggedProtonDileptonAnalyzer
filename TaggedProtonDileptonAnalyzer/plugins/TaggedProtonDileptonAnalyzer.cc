// -*- C++ -*-
//
// Package:    TaggedProtonDileptonAnalyzer/TaggedProtonDileptonAnalyzer
// Class:      TaggedProtonDileptonAnalyzer
//
/**\class TaggedProtonDileptonAnalyzer TaggedProtonDileptonAnalyzer.cc TaggedProtonDileptonAnalyzer/TaggedProtonDileptonAnalyzer/plugins/TaggedProtonDileptonAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Simple ntuple maker for PPS timing detector studies. Runs on AOD, and writes flat ntuples with arrays of RecHit, local track, and proton object (if available) variables. 
// Also saves information on central CMS vertices for low-pileup studies
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/Registry.h"

#include "FWCore/Framework/interface/EventSetup.h"

// HLT information                                   
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

// Central                                                                                                                                                                                   
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Vertexing.h"

// RP                                                                                                                                                                                        
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"

// Timing                                                                                                                                                                                    
#include "DataFormats/CTPPSDetId/interface/CTPPSDiamondDetId.h"
#include "DataFormats/CTPPSReco/interface/CTPPSDiamondRecHit.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

// Protons                                                                                                                                                                                   
#include "DataFormats/ProtonReco/interface/ForwardProton.h"
#include "DataFormats/ProtonReco/interface/ForwardProtonFwd.h"

#include "CondFormats/RunInfo/interface/LHCInfo.h"
#include "CondFormats/DataRecord/interface/LHCInfoRcd.h"

#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <TLorentzVector.h>

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class TaggedProtonDileptonAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit TaggedProtonDileptonAnalyzer(const edm::ParameterSet&);
  ~TaggedProtonDileptonAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ----------member data ---------------------------
  // Input tags                                                                                                                                                                               

  std::string outputFile_;

  edm::EDGetTokenT<std::vector<CTPPSLocalTrackLite> > pps_tracklite_token_;

  edm::EDGetTokenT< edm::View<reco::Vertex> > verticesToken_;
  edm::EDGetTokenT< edm::View<pat::Muon> > muonsToken_;

  edm::EDGetTokenT<std::vector<reco::ForwardProton> > recoProtonsSingleRPToken_;
  edm::EDGetTokenT<std::vector<reco::ForwardProton> > recoProtonsMultiRPToken_;

  edm::ESGetToken<LHCInfo, LHCInfoRcd> lhcInfoToken_;

  std::string lhcInfoLabel_;

  // Tree contents                                                                                                                                                           

  // Run+event quantities                                                                                                                                                                   
  Int_t BX, Run, LumiSection, EventNum;
  Float_t CrossingAngle;

  Int_t nVertices, nArmsTiming, nTracksTiming;

  Int_t nLiteTracks;
  Float_t TrackLiteX[100], TrackLiteY[100], TrackLiteZ[100], TrackLiteTime[100], TrackLiteTimeUnc[100];
  Int_t TrackLiteRPID[100], TrackLiteArm[100];

  Double_t PrimVertexZ[100], PrimVertexY[100], PrimVertexX[100];
  Int_t PrimVertexIsBS[100], PrimVertexNtracks[100], PrimVertexNdof[100];

  Int_t nProtons;
  Float_t ProtonXi[100];
  Float_t ProtonThY[100];
  Float_t ProtonThX[100];
  Float_t Protont[100];
  Int_t ProtonIsMultiRP[100];
  Int_t ProtonRPID[100];
  Int_t ProtonArm[100];
  Float_t ProtonTime[100];

  Int_t nMuons = 0;
  Float_t MuPlusPt = 0;
  Float_t MuMinusPt = 0;
  Float_t MuPlusEta = 0;
  Float_t MuMinusEta = 0;
  Float_t MuPlusPhi = 0;
  Float_t MuMinusPhi = 0;
  Float_t MuPlusM = 0;
  Float_t MuMinusM = 0;

  Float_t MuMuXi = 0;
  Float_t MuMuY = 0;
  Float_t MuMuM = 0;
  Float_t MuMuAcop = 0;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//
TFile* file;
TTree* tree;

//
// constructors and destructor
//
TaggedProtonDileptonAnalyzer::TaggedProtonDileptonAnalyzer(const edm::ParameterSet& iConfig)
  : pps_tracklite_token_ ( consumes<std::vector<CTPPSLocalTrackLite>>(iConfig.getParameter<edm::InputTag>("tagTrackLites") ) ),
    verticesToken_    ( consumes< edm::View<reco::Vertex> >( iConfig.getParameter<edm::InputTag>( "verticesTag" ) ) ),
    muonsToken_    ( consumes< edm::View<pat::Muon> >( iConfig.getParameter<edm::InputTag>( "muonsTag" ) ) ),
    recoProtonsSingleRPToken_   ( consumes<std::vector<reco::ForwardProton> >      ( iConfig.getParameter<edm::InputTag>( "ppsRecoProtonSingleRPTag" ) ) ),
    recoProtonsMultiRPToken_   ( consumes<std::vector<reco::ForwardProton> >      ( iConfig.getParameter<edm::InputTag>( "ppsRecoProtonMultiRPTag" ) ) ) {

  //now do what ever initialization is needed
  outputFile_ = iConfig.getUntrackedParameter<std::string>("outfilename", "output.root");
  file = new TFile(outputFile_.c_str(), "recreate");
  file->cd();
  // tree definition                                                                                                                                                         
  tree = new TTree("ntp1", "ntp1");

}

TaggedProtonDileptonAnalyzer::~TaggedProtonDileptonAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  file->Write();
  file->Close();

}

//
// member functions
//

// ------------ method called for each event  ------------
void TaggedProtonDileptonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  // Vertices in the central CMS tracker                
  edm::Handle< edm::View<reco::Vertex> > vertices;
  iEvent.getByToken( verticesToken_, vertices );

  edm::Handle< edm::View<pat::Muon> > muons;
  iEvent.getByToken( muonsToken_, muons );



  // Run and BX information                                                                                                                          
  BX = iEvent.bunchCrossing();
  Run = iEvent.id().run();
  LumiSection = iEvent.luminosityBlock();
  EventNum = iEvent.id().event();

  // Initializing counters                                                                                                                           
  nProtons = 0;
  nTracksTiming = 0;
  nLiteTracks = 0;
  nVertices = 0;

  nMuons = 0;

  // Initializing arrays                                                                                                                             
  for(int i = 0; i < 100; i++)
    {
      PrimVertexZ[i] = -999.;
      PrimVertexX[i] = -999.;
      PrimVertexY[i] = -999.;
      PrimVertexIsBS[i] = -999;
      PrimVertexNtracks[i] = -999;
      PrimVertexNdof[i] = -999;
      
      ProtonXi[i] = -999;
      ProtonThY[i] = -999;
      ProtonThX[i] = -999;
      Protont[i] = -999;
      ProtonIsMultiRP[i] = -999;
      ProtonRPID[i] = -999;
      ProtonArm[i] = -999;
      ProtonTime[i] = -999.;
      TrackLiteX[i] = 0;
      TrackLiteY[i] = 0;
      TrackLiteZ[i] = 0;
      TrackLiteArm[i] = 0;
      TrackLiteTime[i] = 0;
      TrackLiteTimeUnc[i] = 0;
      TrackLiteRPID[i] = 0;
    }

  /* Primary vertices in central CMS tracker*/
  for ( const auto& vtx : *vertices ) {
    PrimVertexZ[nVertices] = vtx.z();
    PrimVertexX[nVertices] = vtx.x();
    PrimVertexY[nVertices] = vtx.y();
    PrimVertexNtracks[nVertices] = vtx.nTracks(0);
    PrimVertexNdof[nVertices] = vtx.ndof();

    if(vtx.isFake() == 1)
      PrimVertexIsBS[nVertices] = 1;
    else
      PrimVertexIsBS[nVertices] = 0;
    nVertices++;
  }


  MuPlusPt = 0;
  MuMinusPt = 0;
  MuPlusEta = 0;
  MuMinusEta = 0;
  MuPlusPhi = 0;
  MuMinusPhi = 0;
  MuPlusM = 0;
  MuMinusM = 0;

  MuMuXi = 0;
  MuMuY = 0;
  MuMuM = 0;

  TLorentzVector mup, mum, mumu;

  /* Muons */
  for( const auto& mu1 : *muons )
    {
      if(mu1.pt() > 25 && mu1.charge()==1 && mu1.isTightMuon(vertices->at(0)))
	{
	  if(mu1.pt() > MuPlusPt)
	    {
	      MuPlusPt = mu1.pt();
	      MuPlusEta = mu1.eta();
              MuPlusPhi = mu1.phi();
              MuPlusM = mu1.mass();
	    }
	}
      if(mu1.pt() > 25 && mu1.charge()==-1 && mu1.isTightMuon(vertices->at(0)))
	{
	  if(mu1.pt() > MuMinusPt)
	    {
	      MuMinusPt = mu1.pt();
	      MuMinusEta = mu1.eta();
	      MuMinusPhi = mu1.phi();
	      MuMinusM = mu1.mass();
	    }
	}
    }
  
  if(MuPlusPt>0 && MuMinusPt>0)
    {
      mup.SetPtEtaPhiM(MuPlusPt,MuPlusEta,MuPlusPhi,MuPlusM);
      mum.SetPtEtaPhiM(MuMinusPt,MuMinusEta,MuMinusPhi,MuMinusM);
      mumu = mup+mum;
      MuMuM = mumu.M();
      MuMuY = mumu.Rapidity();

      if(MuMuY > 0)
	MuMuXi = (1.0/13600.0)*((MuPlusPt*TMath::Exp(MuPlusEta)) + (MuMinusPt*TMath::Exp(MuMinusEta)));
      if(MuMuY < 0)
        MuMuXi = (1.0/13600.0)*((MuPlusPt*TMath::Exp(-1.0*MuPlusEta)) + (MuMinusPt*TMath::Exp(-1.0*MuMinusEta)));

      float dphi = fabs(MuPlusPhi-MuMinusPhi);
      if(dphi > 3.14159)
	dphi = (2.0*3.14159) - dphi;
      float acop = 1.0 - (dphi/3.14159);
      MuMuAcop = acop;

    }

  /* Proton lite tracks */
  edm::Handle<std::vector<CTPPSLocalTrackLite> > ppsTracksLite;
  iEvent.getByToken( pps_tracklite_token_, ppsTracksLite );

  for ( const auto& trklite : *ppsTracksLite )
    {
      const CTPPSDetId detid( trklite.rpId() );

      // transform the raw, 32-bit unsigned integer detId into the TOTEM "decimal" notation                                                            
      const unsigned short raw_id = 100*detid.arm()+10*detid.station()+detid.rp();

      TrackLiteX[nLiteTracks] = trklite.x();
      TrackLiteY[nLiteTracks] = trklite.y();
      TrackLiteTime[nLiteTracks] = trklite.time();
      TrackLiteTimeUnc[nLiteTracks] = trklite.timeUnc();
      TrackLiteRPID[nLiteTracks] = raw_id;

      if(raw_id == 3)
	{
	  TrackLiteZ[nLiteTracks] = 210;
	  TrackLiteArm[nLiteTracks] = 0;
	}
      if(raw_id == 23)
        {
          TrackLiteZ[nLiteTracks] = 220;
          TrackLiteArm[nLiteTracks] = 0;
	}
      if(raw_id == 16)
        {
          TrackLiteZ[nLiteTracks] = 216;
          TrackLiteArm[nLiteTracks] = 0;
        }
      if(raw_id == 103)
        {
          TrackLiteZ[nLiteTracks] = 210;
          TrackLiteArm[nLiteTracks] = 1;
        }
      if(raw_id == 123)
        {
          TrackLiteZ[nLiteTracks] = 220;
          TrackLiteArm[nLiteTracks] = 1;
        }
      if(raw_id == 116)
        {
          TrackLiteZ[nLiteTracks] = 216;
          TrackLiteArm[nLiteTracks] = 1;
        }


      nLiteTracks++;
    }

  /* Full Reco protons */
  edm::Handle<std::vector<reco::ForwardProton>> recoMultiRPProtons;
  iEvent.getByToken(recoProtonsMultiRPToken_, recoMultiRPProtons);
  edm::Handle<std::vector<reco::ForwardProton>> recoSingleRPProtons;
  iEvent.getByToken(recoProtonsSingleRPToken_, recoSingleRPProtons);

  int ismultirp = -999;
  unsigned int decRPId = -999;
  unsigned int armId = -999;
  float th_y = -999;
  float th_x = -999;
  float t = -999;
  float xi = -999.;
  float protontime = -999.;

  // Proton object with Single-RP algorithm                                                                                                          
  for (const auto & proton : *recoSingleRPProtons)
    {
      if (proton.validFit())
        {
          th_y = proton.thetaY();
          th_x = proton.thetaX();
          xi = proton.xi();
          //      t = proton.t();                                                                                                                      
          protontime = proton.time();

          CTPPSDetId rpId((*proton.contributingLocalTracks().begin())->rpId());
          decRPId = rpId.arm()*100 + rpId.station()*10 + rpId.rp();

          ismultirp = 0;

          ProtonXi[nProtons] = xi;
          ProtonThX[nProtons] = th_x;
          ProtonThY[nProtons] = th_y;
          Protont[nProtons] = t;
          ProtonIsMultiRP[nProtons] = ismultirp;
          ProtonRPID[nProtons] = decRPId;
          ProtonArm[nProtons] = armId;
          ProtonTime[nProtons] = protontime;
          nProtons++;
        }
    }

  // Proton object with Multi-RP algorithm                                                                                                          
  for (const auto & proton : *recoMultiRPProtons)
    {
      if (proton.validFit())
        {
          th_y = proton.thetaY();
          th_x = proton.thetaX();
          xi = proton.xi();
          //      t = proton.t();                                                                                                                    
          protontime = proton.time();

          CTPPSDetId rpId((*proton.contributingLocalTracks().begin())->rpId());
          armId = rpId.arm();

          ismultirp = 1;

          ProtonXi[nProtons] = xi;
          ProtonThX[nProtons] = th_x;
          ProtonThY[nProtons] = th_y;
          Protont[nProtons] = t;
          ProtonIsMultiRP[nProtons] = ismultirp;
          ProtonRPID[nProtons] = decRPId;
          ProtonArm[nProtons] = armId;
          ProtonTime[nProtons] = protontime;
          nProtons++;
        }
    }

  if(MuMuM > 0)
    tree->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void TaggedProtonDileptonAnalyzer::beginJob() {
  // Booking the ntuple                                                                                                                                

  tree->Branch("Run", &Run, "Run/I");
  tree->Branch("LumiSection", &LumiSection, "LumiSection/I");
  tree->Branch("BX", &BX, "BX/I");
  tree->Branch("EventNum", &EventNum, "EventNum/I");
  tree->Branch("CrossingAngle", &CrossingAngle, "CrossingAngle/F");

  tree->Branch("nLiteTracks", &nLiteTracks, "nLiteTracks/I");
  tree->Branch("TrackLiteX", &TrackLiteX, "TrackLiteX[nLiteTracks]/F");
  tree->Branch("TrackLiteY", &TrackLiteY, "TrackLiteY[nLiteTracks]/F");
  tree->Branch("TrackLiteZ", &TrackLiteZ, "TrackLiteZ[nLiteTracks]/F");
  tree->Branch("TrackLiteTime", &TrackLiteTime, "TrackLiteTime[nLiteTracks]/F");
  tree->Branch("TrackLiteTimeUnc", &TrackLiteTimeUnc, "TrackLiteTimeUnc[nLiteTracks]/F");
  tree->Branch("TrackLiteRPID", &TrackLiteRPID, "TrackLiteRPID[nLiteTracks]/I");
  tree->Branch("TrackLiteArm", &TrackLiteArm, "TrackLiteArm[nLiteTracks]/I");

  tree->Branch("nVertices", &nVertices, "nVertices/I");
  tree->Branch("PrimVertexZ", &PrimVertexZ, "PrimVertexZ[nVertices]/D");
  tree->Branch("PrimVertexX", &PrimVertexX, "PrimVertexX[nVertices]/D");
  tree->Branch("PrimVertexY", &PrimVertexY, "PrimVertexY[nVertices]/D");
  tree->Branch("PrimVertexNtracks", &PrimVertexNtracks, "PrimVertexNtracks[nVertices]/I");
  tree->Branch("PrimVertexIsBS", &PrimVertexIsBS, "PrimVertexIsBS[nVertices]/I");
  tree->Branch("PrimVertexNdof", &PrimVertexNdof, "PrimVertexNdof[nVertices]/D");

  tree->Branch("nProtons", &nProtons, "nProtons/I");
  tree->Branch("ProtonXi", &ProtonXi, "ProtonXi[nProtons]/F");
  tree->Branch("ProtonThX", &ProtonThX, "ProtonThX[nProtons]/F");
  tree->Branch("ProtonThY", &ProtonThY, "ProtonThY[nProtons]/F");
  tree->Branch("Protont", &Protont, "Protont[nProtons]/F");
  tree->Branch("ProtonIsMultiRP", &ProtonIsMultiRP, "ProtonIsMultiRP[nProtons]/I");
  tree->Branch("ProtonRPID", &ProtonRPID, "ProtonRPID[nProtons]/I");
  tree->Branch("ProtonArm", &ProtonArm, "ProtonArm[nProtons]/I");
  tree->Branch("ProtonTime", &ProtonTime, "ProtonTime[nProtons]/F");

  tree->Branch("MuMuM",&MuMuM, "MuMuM/F");
  tree->Branch("MuMuY",&MuMuY, "MuMuY/F");
  tree->Branch("MuMuXi",&MuMuXi,"MuMuXi/F");
  tree->Branch("MuMuAcop",&MuMuAcop,"MuMuAcop/F");

}

// ------------ method called once each job just after ending the event loop  ------------
void TaggedProtonDileptonAnalyzer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TaggedProtonDileptonAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TaggedProtonDileptonAnalyzer);
