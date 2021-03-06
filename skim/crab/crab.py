from CRABClient.UserUtilities import config
#submit with 'python crab.py'
#Don't write to my directory (schoef), though

config = config()
config.General.requestName = 'ZeroBias1_Run2015A-PromptReco-v1_RECO'
config.General.workArea = 'private0TSkim_v3'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../python/skim.py'

config.Data.inputDataset = '/ZeroBias1/Run2015A-PromptReco-v1/RECO'
config.Data.inputDBS = 'global'
config.Data.lumiMask = 'json_DCSONLY_150710_skim_v3.txt' 
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 10

config.Data.publication = False
#config.Data.outLFNDirBase = '' 
#config.Data.publishDataName = ''

config.Data.outLFNDirBase = '/store/group/phys_jetmet/schoef/private0TSkim_v3/'
config.Site.storageSite = 'T2_CH_CERN'

datasets=[
'/EGamma/Run2015A-PromptReco-v1/RECO',
'/Jet/Run2015A-PromptReco-v1/RECO',
'/ZeroBias1/Run2015A-PromptReco-v1/RECO',
'/ZeroBias2/Run2015A-PromptReco-v1/RECO',
'/ZeroBias3/Run2015A-PromptReco-v1/RECO',
'/ZeroBias4/Run2015A-PromptReco-v1/RECO',
'/ZeroBias5/Run2015A-PromptReco-v1/RECO',
'/ZeroBias6/Run2015A-PromptReco-v1/RECO',
'/ZeroBias7/Run2015A-PromptReco-v1/RECO',
'/ZeroBias8/Run2015A-PromptReco-v1/RECO',
'/SingleMu/Run2015A-PromptReco-v1/RECO',
]

if __name__ == '__main__':
    from CRABAPI.RawCommand import crabCommand
    for dataset in datasets:
        config.Data.inputDataset = dataset
        config.General.requestName = dataset.rstrip('/').lstrip('/').replace('/','_')
#        print config.General.requestName
        crabCommand('submit', config = config)
