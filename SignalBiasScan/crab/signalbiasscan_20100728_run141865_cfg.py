import FWCore.ParameterSet.Config as cms

process = cms.Process("Tree")

process.load("FWCore.MessageService.MessageLogger_cfi")


 # Conditions (Global Tag is used here):
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = 'GR_R_42_V14B::All'
process.GlobalTag.globaltag = 'FT_R_44_V9::All'

#Geometry and field
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("TrackingTools.RecoGeometry.RecoGeometries_cff")
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")


# reduce verbosity
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    '/store/data/Run2011A/MinimumBias/RECO/May10ReReco-v2/0000/E6E18C6D-718A-E011-9B93-002618943867.root'
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
    fullHitInfo_TEC       = cms.bool(False)
)


process.p = cms.Path(process.TrackRefitter*process.demo)

