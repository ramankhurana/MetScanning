import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(
        'root://eoscms.cern.ch//store/group/phys_jetmet/schoef/private0TSkim_v3/ZeroBias8/crab_ZeroBias8_Run2015A-PromptReco-v1_RECO/150610_094947/0000/skim_98.root'
        )
                            )

process.demo = cms.EDAnalyzer('MetAnalyzer'
                              )


process.p = cms.Path(process.demo)
