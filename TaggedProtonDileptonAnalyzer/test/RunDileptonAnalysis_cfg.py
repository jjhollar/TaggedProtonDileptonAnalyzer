
import FWCore.ParameterSet.Config as cms
import copy

process = cms.Process('ProtonsDimuons')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring(
                                 '/store/data/Run2022D/Muon/MINIAOD/PromptReco-v2/000/357/802/00000/80450074-ef3d-48d1-a5e1-42cd71acaa20.root'
                        )
)


process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
from Configuration.StandardSequences.FrontierConditions_GlobalTag_cff import GlobalTag
process.GlobalTag.globaltag = '124X_dataRun3_Prompt_frozen_v4'
process.GlobalTag.toGet = cms.VPSet()

process.mydimuons = cms.EDAnalyzer(
    'TaggedProtonDileptonAnalyzer',
    #    lhcInfoLabel = cms.string(''),
    verticesTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muonsTag = cms.InputTag("slimmedMuons"),
    tagTrackLites = cms.InputTag( "ctppsLocalTrackLiteProducer"),
    ppsRecoProtonSingleRPTag = cms.InputTag("ctppsProtons", "singleRP"),
    ppsRecoProtonMultiRPTag = cms.InputTag("ctppsProtons", "multiRP"),
    outfilename = cms.untracked.string( "output_Dimuon.root" )
)

# Trigger                                                                                                                  
from HLTrigger.HLTfilters.hltHighLevel_cfi import *
process.hltFilter = copy.deepcopy(hltHighLevel)
process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilter.HLTPaths = ['HLT_IsoMu24_v*']

process.ALL = cms.Path(
    process.hltFilter * 
    # Uncomment these lines, to re-run the PPS local+proton timing reconstruction starting from AOD
    process.mydimuons
                       )

process.schedule = cms.Schedule(process.ALL)

#print(process.dumpPython())
