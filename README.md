#MetAnalyzer 
## Install 
```
cmsrel CMSSW_7_4_4_patch2
cd CMSSW_7_4_4_patch2/src
cmsenv
git clone git@github.com:ramankhurana/MetScanning.git
scram b -j9
cd MetScanning/MetAnalyzer
cmsRun ConfFile_cfg.py
```


# MetScanning
## Install
```
  cmsrel CMSSW_7_4_4_patch2
  cd CMSSW_7_4_4_patch2/src
  cmsenv
  git clone https://github.com/cms-met/MetScanning
  scram b -j9
  ```
## Run on local file
```
  cmsRun MetScanning/skim/python/skim.py
```
## Run with crab
In ``MetScanning/skim/crab/`` edit crab.py and adjust samples, JSON, and the EOS directory. 
Then do
```
  cd MetScanning/skim/crab/
  python crab.py
```
