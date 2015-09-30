import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("testParticle")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring("file:TQToHCQ_HToMuMu_M-30_8TeV_pythia8175_step0.root"),
    skipEvents = cms.untracked.uint32(0),
)

dirname = "/eos/uscms/store/user/lpchzg/FCNH_MC/jiafulow/TQToHCQ_HToMuMu_M-30_8TeV_pythia8175_20150903/TQToHCQ_HToMuMu_M-30_8TeV_pythia8175_20150903/c3f6222a1159f05a0c6ed04f258385d9/"
fileNames = []
for fname in os.listdir(dirname):
    fname1 = os.path.join(dirname[10:], fname)  # strips away "/eos/uscms" (10 chars)
    fileNames.append(fname1)
#for fname in fileNames: print fname

process.source.fileNames  = fileNames
#process.maxEvents.input   = 1000
#process.source.skipEvents = 700

process.printTree1 = cms.EDAnalyzer("ParticleListDrawer",
    src = cms.InputTag("genParticles"),
    maxEventsToPrint  = cms.untracked.int32(1),
    printVertex = cms.untracked.bool(False),
    printOnlyHardInteraction = cms.untracked.bool(False), # Print only status=3 particles. This will not work for Pythia8, which does not have any such particles.
    useMessageLogger = cms.untracked.bool(False),
)

process.printTree2 = cms.EDAnalyzer("ParticleTreeDrawer",
    src = cms.InputTag("genParticles"),
    printP4 = cms.untracked.bool(False),
    printPtEtaPhi = cms.untracked.bool(False),
    printVertex = cms.untracked.bool(False),
    printStatus = cms.untracked.bool(False),
    printIndex  = cms.untracked.bool(False),
    status = cms.untracked.vint32(),  # Print only status=3 particles. This will not work for Pythia8, which does not have any such particles.
)

process.load("NWU.GenAnalyzers.pythia8ParticleAnalyzer_cfi")
process.pythia8ParticleAnalyzer.numberShowInfo    = 10
process.pythia8ParticleAnalyzer.numberShowProcess = 10
process.pythia8ParticleAnalyzer.numberShowEvent   = 10

process.printEventNumber = cms.OutputModule("AsciiOutputModule")

#process.p = cms.Path(process.printTree1*process.printTree2)
process.p = cms.Path(process.printTree1*process.pythia8ParticleAnalyzer)
#process.outpath = cms.EndPath(process.printEventNumber)
#process.MessageLogger.destinations = cms.untracked.vstring('cout','cerr')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

