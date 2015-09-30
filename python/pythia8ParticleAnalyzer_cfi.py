import FWCore.ParameterSet.Config as cms

pythia8ParticleAnalyzer = cms.EDAnalyzer("Pythia8ParticleAnalyzer",
    src = cms.InputTag("genParticles"),
    numberShowInfo = cms.untracked.int32(1),
    numberShowProcess = cms.untracked.int32(1),
    numberShowEvent = cms.untracked.int32(1),
    useMessageLogger = cms.untracked.bool(False),
    histFile = cms.untracked.string("histograms.root"),
    makePlots = cms.untracked.bool(True),
)
