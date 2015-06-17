// -*- C++ -*-
//
// Package:    MetScanning/MetAnalyzer
// Class:      MetAnalyzer
// 
/**\class MetAnalyzer MetAnalyzer.cc MetScanning/MetAnalyzer/plugins/MetAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Raman Khurana
//         Created:  Wed, 10 Jun 2015 18:44:03 GMT
//
//

#define nfilter  4
// system include files
#include <memory>
#include <iostream>
// user include files
#include "TProfile2D.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/PFClusterMET.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "FWCore/Common/interface/TriggerNames.h" 

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TGraphAsymmErrors.h"
//
// class declaration
//

class MetAnalyzer : public edm::EDAnalyzer {
   public:
      explicit MetAnalyzer(const edm::ParameterSet&);
      ~MetAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
  //TFile* f;
  edm::Service<TFileService> fs;
  bool hbhet;
  bool csct;
  TTree* metTree;

  Double_t  calomet_        ;
  Double_t  pfclustermet_   ;
  Double_t  pfcalomet_      ; 
  
  Double_t  calometphi_     ; 
  Double_t  pfclustermetphi_;
  Double_t  pfcalometphi_   ; 
  
  Bool_t   hbhet_          ; 
  Bool_t   csct_           ; 
  
  ULong64_t  run      ; 
  ULong64_t  lumi     ; 
  ULong64_t  event    ; 
  
  
  TH1F* caloMET[nfilter];
  TH1F* PFclusterMET[nfilter];
  TH1F* PFcaloMET[nfilter];
  
  TH1F* caloMETPhi[nfilter];
  TH1F* PFclusterMETPhi[nfilter];
  TH1F* PFcaloMETPhi[nfilter];
  
  TH2F* caloMET_vs_PFclusterMET[nfilter];
  TH2F* caloMET_vs_PFcaloMET[nfilter];
  TH2F* PFclusterMET_vs_PFcaloMET[nfilter];
  
  TH2F* caloMET_vs_caloMETPhi[nfilter];
  TProfile2D* pfClECALMap;
  TProfile2D* pfClHCALMap;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MetAnalyzer::MetAnalyzer(const edm::ParameterSet& iConfig)

{
  
  metTree = fs->make<TTree>("metTree","metTree");
  //f = new TFile("MetScaning_1.root","RECREATE");
  std::vector<TString> postfix;
  postfix.clear();
  postfix.push_back("NoFilter");
  postfix.push_back("hbhet");
  postfix.push_back("csct");
  postfix.push_back("all");
  
  for(size_t i=0; i<postfix.size(); i++){
    caloMET[i] = fs->make<TH1F>("caloMET"+postfix[i],"caloMET;caloMET;# of events (normalized to 1)",4000,0,8000);
    PFclusterMET[i] = fs->make<TH1F>("PFclusterMET"+postfix[i],"PFclusterMET;caloMET;# of events (normalized to 1)",4000,0,8000);
    PFcaloMET[i] = fs->make<TH1F>("PFcaloMET"+postfix[i],"PFcaloMET;caloMET;# of events (normalized to 1)",4000,0,8000);

    caloMETPhi[i] = fs->make<TH1F>("caloMETPhi"+postfix[i],"caloMET;caloMETPhi;# of events (normalized to 1)",70,-3.5,3.5);
    PFclusterMETPhi[i] = fs->make<TH1F>("PFclusterMETPhi"+postfix[i],"PFclusterMETPhi;caloMET;# of events (normalized to 1)",70,-3.5,3.5);
    PFcaloMETPhi[i] = fs->make<TH1F>("PFcaloMETPhi"+postfix[i],"PFcaloMET;caloMETPhi;# of events (normalized to 1)",70,-3.5,3.5);
    
    caloMET_vs_PFclusterMET[i]      = fs->make<TH2F>("caloMET_vs_PFclusterMET"+postfix[i],"caloMET_vs_PFclusterMET;caloMET;PFclusterMET",4000,0,8000,4000,0,8000);
    caloMET_vs_PFcaloMET[i]         = fs->make<TH2F>("caloMET_vs_PFcaloMET"+postfix[i],"caloMET_vs_PFcaloMET;caloMET;PFcaloMET",4000,0,8000,4000,0,8000);
    PFclusterMET_vs_PFcaloMET[i]    = fs->make<TH2F>("PFclusterMET_vs_PFcaloMET"+postfix[i],"PFclusterMET_vs_PFcaloMET;PFclusterMET;PFcaloMET",4000,0,8000,4000,0,8000);
    
    caloMET_vs_caloMETPhi[i]        = fs->make<TH2F>("caloMET_vs_caloMETPhi"+postfix[i],"caloMET_vs_caloMETPhi;caloMET #phi;caloMET",70,-3.5,3.5,4000,0,8000);
    //    if(i>0) {
    //MetEff[i] = fs->make<TGraphAsymmErrors>("");
    //    }
  }
  
  pfClECALMap   =  fs->make<TProfile2D>("pfClECALMap","pfClECALMap;#eta;#phi;<p_{T}> (in GeV)",140,-3.14,3.14,140,-3.14,3.14);
  pfClHCALMap   =  fs->make<TProfile2D>("pfClHCALMap","pfClHCALMap;#eta;#phi;<p_{T}> (in GeV)",140,-3.14,3.14,140,-3.14,3.14);
  //now do what ever initialization is needed

}


MetAnalyzer::~MetAnalyzer()
{
  
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  std::cout<<" working "<<std::endl;
  hbhet = false;
  csct  = false;
  
  edm::Handle<vector<reco::CaloMET> > caloMETH;
  iEvent.getByLabel("caloMet",caloMETH);
  std::vector<reco::CaloMET>::const_iterator cmet = caloMETH.product()->begin();
  //std::cout<<" MET = "<<cmet->et()<<std::endl;
  
  
  edm::Handle<vector<reco::PFClusterMET> > PFClusterMETH;
  iEvent.getByLabel("pfClusterMet",PFClusterMETH);
  std::vector<reco::PFClusterMET>::const_iterator pfclustermet = PFClusterMETH.product()->begin();
  //std::cout<<" PFClusterMETH = "<<pfclustermet->et()<<std::endl;
  
  edm::Handle<vector<reco::PFMET> > PFCaloMETH;
  iEvent.getByLabel("pfCaloMet",PFCaloMETH);
  std::vector<reco::PFMET>::const_iterator pfcalomet = PFCaloMETH.product()->begin();
  //std::cout<<" PFCaloMET = "<<pfcalomet->et()<<std::endl;
  
  // HBHE Tight Filter
  edm::Handle<bool> HBHET;
  edm::InputTag  hbhetag("HBHENoiseFilterResultProducer","HBHENoiseFilterResultRun2Tight","SKIM");
  iEvent.getByLabel(hbhetag,HBHET);
  hbhet = (*HBHET.product());
  if(false) std::cout<<" HBHE = "<<(*HBHET.product())<<std::endl;


  edm::Handle<bool> CSCT;
  edm::InputTag  csctag("CSCTightHaloFilter");
  iEvent.getByLabel(csctag,CSCT);
  csct = (*CSCT.product()) ;
  if(false) std::cout<<" CSCT = "<<(*CSCT.product())<<std::endl;
  
  
  // particleFlowClusterECAL
  edm::Handle<vector<reco::PFCluster> > H_particleFlowClusterECAL;
  iEvent.getByLabel("particleFlowClusterECAL",H_particleFlowClusterECAL);
  std::vector<reco::PFCluster>::const_iterator pfClECAL;
  
  // particleFlowClusterHCAL
  edm::Handle<vector<reco::PFCluster> > H_particleFlowClusterHCAL;
  iEvent.getByLabel("particleFlowClusterHCAL",H_particleFlowClusterHCAL);
  std::vector<reco::PFCluster>::const_iterator pfClHCAL;
  
  // PF ECAL 
  for(pfClECAL = H_particleFlowClusterECAL->begin(); pfClECAL != H_particleFlowClusterECAL->end(); ++pfClECAL){
    if(false) std::cout<<" pfClECAL = "<<pfClECAL->eta()<<"  "<<pfClECAL->phi()<<"  "<<pfClECAL->pt()<<std::endl;
    pfClECALMap->Fill(pfClECAL->eta(),pfClECAL->phi(),pfClECAL->pt());
  }

  // PF HCAL
  for(pfClHCAL = H_particleFlowClusterHCAL->begin(); pfClHCAL != H_particleFlowClusterHCAL->end(); ++pfClHCAL){
    if(false) std::cout<<" pfClHCAL = "<<pfClHCAL->eta()<<"  "<<pfClHCAL->phi()<<"  "<<pfClHCAL->pt()<<std::endl;
    pfClHCALMap->Fill(pfClHCAL->eta(),pfClHCAL->phi(),pfClHCAL->pt());
  }
  
  calomet_         = (Double_t) cmet->et();
  std::cout<<" calomet_ = "<<calomet_<<std::endl;
  pfclustermet_    = (Double_t) pfclustermet->et();
  pfcalomet_       = (Double_t) pfcalomet->et();

  calometphi_      = (Double_t) cmet->phi();
  pfclustermetphi_ = (Double_t) pfclustermet->phi();
  pfcalometphi_    = (Double_t) pfcalomet->phi();

  hbhet_            = (Bool_t) hbhet;
  csct_             = (Bool_t) csct ;
  
  run               = (ULong64_t) iEvent.id().run();
  lumi              = (ULong64_t) iEvent.id().luminosityBlock();
  event             = (ULong64_t) iEvent.id().event();
  
  std::vector<bool> filtervec;
  filtervec.clear();
  filtervec.push_back(true);// This is without any filter.
  filtervec.push_back(hbhet);// For HBHE Tight
  filtervec.push_back(csct);//For CSC Tight)
  filtervec.push_back(hbhet && csct);//hbhet && csct
  
  bool cat1 = (cmet->et()>200 ) || (pfclustermet->et() > 200) || (pfcalomet->et() > 200); //  Any MET is > 250 GeV
  // pfclustermet vs calomet
  bool cat2 = (((pfclustermet->et()     > 1.5*cmet->et()) && (cmet->et()         > 30.0) )); //
  bool cat3 =  ((1.5*pfclustermet->et() < cmet->et())     && (pfclustermet->et() > 30.0) ); //
  //bool cat2 = ((pfclustermet->et()   > 3*cmet->et() ) && (cmet->et()         < 30.0 )) || ((pfclustermet->et()     > 1.5*cmet->et()) && (cmet->et()         > 30.0) ); //
  //bool cat3 = ((3*pfclustermet->et() < cmet->et()   ) && (pfclustermet->et() < 30.0 )) || ((1.5*pfclustermet->et() < cmet->et())     && (pfclustermet->et() > 30.0) ); //
  // pfcalomet vs calomet
  bool cat4 =  ((pfcalomet->et()     > 2.5*cmet->et()) && (cmet->et()         > 30.0) ); //
  bool cat5 =  ((2.5*pfcalomet->et() < cmet->et())     && (pfcalomet->et()    > 30.0) ); //
  //bool cat4 = ((pfcalomet->et()   > 4*cmet->et() ) && (cmet->et()         < 30.0 )) || ((pfcalomet->et()     > 2.5*cmet->et()) && (cmet->et()         > 30.0) ); //
  //bool cat5 = ((4*pfcalomet->et() < cmet->et()   ) && (pfcalomet->et() < 30.0 ))    || ((2.5*pfcalomet->et() < cmet->et())     && (pfcalomet->et()    > 30.0) ); //
  //pfcalomet vs pfclustermet
  bool cat6 =  ((pfcalomet->et()     > 2.5*pfclustermet->et()) && (pfclustermet->et()  > 30.0) ); //
  bool cat7 = ((1.5*pfcalomet->et() < pfclustermet->et())     && (pfcalomet->et()     > 30.0) ); //
  //bool cat6 = ((pfcalomet->et()   > 4*pfclustermet->et() ) && (pfclustermet->et() < 30.0 )) || ((pfcalomet->et()     > 2.5*pfclustermet->et()) && (pfclustermet->et()  > 30.0) ); //
  //bool cat7 = ((4*pfcalomet->et() < pfclustermet->et()   ) && (pfcalomet->et()    < 30.0 )) || ((1.5*pfcalomet->et() < pfclustermet->et())     && (pfcalomet->et()     > 30.0) ); //
  
  std::vector<bool> categories;
  categories.clear();
  categories.push_back(cat1);
  categories.push_back(cat2);
  categories.push_back(cat3);
  categories.push_back(cat4);
  categories.push_back(cat5);
  categories.push_back(cat6);
  categories.push_back(cat7);
  for (size_t icat=0; icat<categories.size(); icat++){
    if(categories[icat]) {
      std::cout<<" category "<<icat+1<<"  "<<cmet->et()
	       <<" : "<<pfclustermet->et()
	       <<" : "<<pfcalomet->et()
	       <<" : "<<hbhet
	       <<" : "<<csct
	       <<" : "<< iEvent.id().run()<<":"
	       << iEvent.id().luminosityBlock()<<":"
	       << iEvent.id().event()
	       <<std::endl;
      break;
    }
  }
for(size_t ifilter=0; ifilter<filtervec.size(); ifilter++){
  if(!filtervec[ifilter]) continue;
  // 1D Histograms
    caloMET[ifilter]->Fill(cmet->et());
    PFclusterMET[ifilter]->Fill(pfclustermet->et());
    PFcaloMET[ifilter]->Fill(pfcalomet->et());

    caloMETPhi[ifilter]->Fill(cmet->phi());
    PFclusterMETPhi[ifilter]->Fill(pfclustermet->phi());
    PFcaloMETPhi[ifilter]->Fill(pfcalomet->phi());
    
    // 2D Histograms
    caloMET_vs_PFclusterMET[ifilter]->Fill(cmet->et(),pfclustermet->et());
    caloMET_vs_PFcaloMET[ifilter]->Fill(cmet->et(),pfcalomet->et());
    PFclusterMET_vs_PFcaloMET[ifilter]->Fill(pfclustermet->et(),pfcalomet->et());

    caloMET_vs_caloMETPhi[ifilter]->Fill(cmet->phi(),cmet->et());
  }
  
  
  /*
  edm::Handle<edm::TriggerResults> trigResults;
  edm::InputTag trigTag("TriggerResults::SKIM");
  iEvent.getByLabel(trigTag, trigResults);
  const edm::TriggerNames & trigNames = iEvent.triggerNames(*trigResults);
  for (unsigned int i=0; i<trigResults->size(); i++)    {
    if(false)  std::cout<<" trig name = "<<trigNames.triggerName(i)
			<<"    "<<trigResults->accept(i)
			<<std::endl;
  }
  */

 metTree->Fill();
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
MetAnalyzer::beginJob()
{

  metTree->Branch("calomet_",&calomet_,"calomet_/D");
  metTree->Branch("pfclustermet_",&pfclustermet_,"pfclustermet_/D");
  metTree->Branch("pfcalomet_",&pfcalomet_,"pfcalomet_/D");
  metTree->Branch("calometphi_",&calometphi_,"calometphi_/D");
  metTree->Branch("pfclustermetphi_",&pfclustermetphi_,"pfclustermetphi_/D");
  metTree->Branch("pfcalometphi_",&pfcalometphi_,"pfcalometphi_/D");
  
  metTree->Branch("hbhet_",&hbhet_,"hbhet_/O");
  metTree->Branch("csct_",&csct_,"csct_/O");
  
  metTree->Branch("run",&run,"run/L");
  metTree->Branch("lumi",&lumi,"lumi/L");
  metTree->Branch("event",&event,"event/L");

}

// ------------ method called once each job just after ending the event loop  ------------
void 
MetAnalyzer::endJob() 
  
{
  
  //metTree->Write();
  //f->cd();
  //for(int i=0;i<4;i++){
  //  caloMET[i]->Write();
  //  PFclusterMET[i]->Write();
  //  PFcaloMET[i]->Write();
  //  caloMET_vs_PFclusterMET[i]->Write();
  //  caloMET_vs_PFcaloMET[i]->Write();
  //  PFclusterMET_vs_PFcaloMET[i]->Write();
  //}
  //f->Close();
}

// ------------ method called when starting to processes a run  ------------
/*
void 
MetAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
MetAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
MetAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
MetAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MetAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MetAnalyzer);
