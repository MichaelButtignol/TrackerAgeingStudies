import FWCore.ParameterSet.Config as cms

process = cms.Process("Tree")

process.load("FWCore.MessageService.MessageLogger_cfi")


 # Conditions (Global Tag is used here):
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
# CMSSW523
process.GlobalTag.globaltag = 'GR_P_V41::All'

#Geometry and field
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("TrackingTools.RecoGeometry.RecoGeometries_cff")
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")


# reduce verbosity
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    '/store/data/Run2012D/SingleElectron/RECO/PromptReco-v1/000/203/835/02692CE4-760B-E211-B9A1-BCAEC518FF8A.root'
    )
)


#name of the output file containing the tree
process.TFileService = cms.Service("TFileService", fileName = cms.string("clustersTree.root") )


process.demo = cms.EDAnalyzer('SignalBiasScan',
    trackLabel       = cms.InputTag('generalTracks'),
    tkTraj           = cms.InputTag('TrackRefitter'),
    labelTrajToTrack = cms.InputTag('TrackRefitter'),
    siStripClusters  = cms.InputTag('TrackRefitter'),
    primaryVertexColl= cms.InputTag('offlinePrimaryVertices'),
    fullHitInfo_TIB       = cms.bool(False),
    fullHitInfo_TOB       = cms.bool(False),
    fullHitInfo_TID       = cms.bool(False),
    fullHitInfo_TEC       = cms.bool(False),
    isSmallBiasScan       = cms.untracked.bool(False)
)


process.p = cms.Path(process.TrackRefitter*process.demo)

