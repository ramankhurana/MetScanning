import FWCore.ParameterSet.Config as cms

process = cms.Process("METTailAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:skim0922test.root'
        #'file:/tmp/khurana/MET.root'
        #'file:../skim/python/skim_247072_132_132378357_Jet.root'
        #'file:/afs/cern.ch/work/k/khurana/METScanners/WithBugFix/CMSSW_7_4_4_patch2/src/MetScanning/skim/python/MET.root'
        #'file:/afs/cern.ch/work/k/khurana/METScanners/WithBugFix/CMSSW_7_4_4_patch2/src/MetScanning/skim/python/MET.root'
        #"root://eoscms.cern.ch//store/group/phys_jetmet/schoef/private0TSkim_v3/ZeroBias8/crab_ZeroBias8_Run2015A-PromptReco-v1_RECO/150610_094947/0000/skim_1.root",
        #'file:/afs/cern.ch/work/c/chchen/public/skim.root'
        #'file:/afs/cern.ch/work/k/khurana/METScanners/CMSSW_7_4_4_patch2/src/MetScanning/skim/python/skim_247072_132_132378357_Jet.root'
        )
)

process.Metanalyzer = cms.EDAnalyzer('MetAnalyzer'
)


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("METScanning0727.root")
                                   )
process.p = cms.Path(process.Metanalyzer)
