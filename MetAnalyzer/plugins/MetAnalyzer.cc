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
#include "utils.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/METReco/interface/PFClusterMET.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "FWCore/Common/interface/TriggerNames.h" 
#include "DataFormats/METReco/interface/HcalNoiseSummary.h"
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
  void Setetaphi();
  bool FindEvent(float eta, float phi);
  std::vector<float> etaVec;
  std::vector<float> phiVec;
  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
  // ----------member data ---------------------------
  //TFile* f;
  edm::Service<TFileService> fs;
  bool hbhet;
  bool csct;
  bool hbhetr1;
  

  TTree* metTree;

  // Add variables to be stored in the branch
  Double_t  calomet_        ;
  Double_t  pfclustermet_   ;
  Double_t  pfcalomet_      ; 
  Double_t  pfcalometPhi_          ; // 
  
  Double_t  calometphi_     ; 
  Double_t  pfclustermetphi_;
  Double_t  pfcalometphi_   ; 

  Double_t caloSumEt_       ;
  Double_t pfclusterSumEt_  ;
  Double_t pfcalocaloSumEt_ ;
  
  Bool_t   hbhet_          ; 
  Bool_t   csct_           ; 
  Bool_t hbhetR1_          ;
  
  Bool_t ecaldeadcellTPFilter_;
  Bool_t eebadscf_;
  Float_t NIsolatedNoiseChannels_;
  Float_t IsolatedNoiseSumE_;
  Float_t IsolatedNoiseSumEt_;
  Bool_t isolationfilterstatus_;
  ULong64_t  run      ; 
  ULong64_t  lumi     ; 
  ULong64_t  event    ; 
  
  
  TH1F* caloMET[nfilter];
  TH1F* PFclusterMET[nfilter];
  TH1F* PFcaloMET[nfilter];
  
  TH1F* caloSumET[nfilter];
  TH1F* PFclusterSumET[nfilter];
  TH1F* PFcaloSumET[nfilter];

  TH1F* caloMET_over_SumET[nfilter];
  TH1F* PFclusterMET_over_SumET[nfilter];
  TH1F* PFcaloMET_over_SumET[nfilter];
  
  TH1F* caloMETPhi[nfilter];
  TH1F* PFclusterMETPhi[nfilter];
  TH1F* PFcaloMETPhi[nfilter];
  
  TH2F* caloMET_vs_PFclusterMET[nfilter];
  TH2F* caloMET_vs_PFcaloMET[nfilter];
  TH2F* PFclusterMET_vs_PFcaloMET[nfilter];
  
  TH2F* caloMET_vs_caloMETPhi[nfilter];
  TProfile2D* pfClECALMap;
  TProfile2D* pfClHCALMap;
  TH2F* calojet_map;
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
    PFclusterMET[i] = fs->make<TH1F>("PFclusterMET"+postfix[i],"PFclusterMET;PFClusterMET;# of events (normalized to 1)",4000,0,8000);
    PFcaloMET[i] = fs->make<TH1F>("PFcaloMET"+postfix[i],"PFcaloMET;PFcaloMET;# of events (normalized to 1)",4000,0,8000);


    caloSumET[i] = fs->make<TH1F>("caloSumET"+postfix[i],"caloSumET;caloSumET;# of events (normalized to 1)",4000,0,8000);
    PFclusterSumET[i] = fs->make<TH1F>("PFclusterSumET"+postfix[i],"PFclusterSumET;caloSumET;# of events (normalized to 1)",4000,0,8000);
    PFcaloSumET[i] = fs->make<TH1F>("PFcaloSumET"+postfix[i],"PFcaloSumET;PFcaloSumET;# of events (normalized to 1)",4000,0,8000);

    caloMET_over_SumET[i] = fs->make<TH1F>("caloMET_over_SumET"+postfix[i],"caloMET_over_SumET;caloSumET;# of events (normalized to 1)",400,0,1.0);
    PFclusterMET_over_SumET[i] = fs->make<TH1F>("PFclusterMET_over_SumET"+postfix[i],"PFclusterMET_over_SumET;caloSumET;# of events (normalized to 1)",400,0,1.0);
    PFcaloMET_over_SumET[i] = fs->make<TH1F>("PFcaloMET_over_SumET"+postfix[i],"PFcaloMET_over_SumET;PFcaloSumET;# of events (normalized to 1)",400,0,1.0);

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

  calojet_map   = fs->make<TH2F>("calojet_map","calojet_map;#eta;#phi",70,-3.5,3.5, 70,-3.5,3.5);
  //now do what ever initialization is needed
  Setetaphi();
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
  //std::cout<<" working "<<std::endl;
  hbhet = false;
  csct  = false;
  hbhetr1=false;
  
  ecaldeadcellTPFilter_=false;
  eebadscf_ = false;
  NIsolatedNoiseChannels_=-999.;
  IsolatedNoiseSumE_=-999.0;
  IsolatedNoiseSumEt_=-999.0;
  isolationfilterstatus_=true;;

  
  edm::Handle<vector<reco::CaloMET> > caloMETH;
  iEvent.getByLabel("caloMet",caloMETH);
  std::vector<reco::CaloMET>::const_iterator cmet = caloMETH.product()->begin();
  //std::cout<<" MET = "<<cmet->et()<<std::endl;


  edm::Handle<vector<reco::PFMET> > pfMETH;
  iEvent.getByLabel("pfMet",pfMETH);
  std::vector<reco::PFMET>::const_iterator pfmet = pfMETH.product()->begin();
  std::cout<<" PFMET ========== "<<pfmet->et()<<std::endl;
  
  
  edm::Handle<vector<reco::PFClusterMET> > PFClusterMETH;
  iEvent.getByLabel("pfClusterMet",PFClusterMETH);
  std::vector<reco::PFClusterMET>::const_iterator pfclustermet = PFClusterMETH.product()->begin();
  //std::cout<<" PFClusterMETH = "<<pfclustermet->et()<<std::endl;
  
  edm::Handle<vector<reco::PFMET> > PFCaloMETH;
  iEvent.getByLabel("pfCaloMet",PFCaloMETH);
  std::vector<reco::PFMET>::const_iterator pfcalomet = PFCaloMETH.product()->begin();
  //std::cout<<" PFCaloMET = "<<pfcalomet->et()<<std::endl;
  
  
  // CaloJets 
  edm::Handle<std::vector<reco::CaloJet> > CaloJetH;
  iEvent.getByLabel("ak4CaloJets",CaloJetH);
  std::vector<reco::CaloJet> jetcoll(*(CaloJetH.product()));
  std::sort(jetcoll.begin(), jetcoll.end(), PtGreater());
  std::vector<reco::CaloJet>::const_iterator calojet = jetcoll.begin();
  


  // HBHE Tight Filter
  edm::Handle<bool> HBHET;
  edm::InputTag  hbhetag("HBHENoiseFilterResultProducer","HBHENoiseFilterResultRun2Tight","SKIM");
  iEvent.getByLabel(hbhetag,HBHET);
  hbhet = (*HBHET.product());
  if(false) std::cout<<" HBHE = "<<(*HBHET.product())<<std::endl;
  
  // HBHET Run 1 Filter 
  edm::Handle<bool> HBHETR1;
  edm::InputTag  hbhetagR1("HBHENoiseFilterResultProducer","HBHENoiseFilterResultRun1","SKIM");
  iEvent.getByLabel(hbhetagR1,HBHETR1);
  hbhetr1 = (*HBHET.product());
  
//edm::Handle<HcalNoiseSummary> hSummary;
//iEvent.getByLabel("hcalnoise", hSummary);
//
//std::cout<<" num of iso channels = "<<hSummary->numIsolatedNoiseChannels()
//	   <<" iso noise sum E = "<<hSummary->isolatedNoiseSumE()
//	   <<" iso noise sum Et = "<<hSummary->isolatedNoiseSumEt()
//	   <<std::endl;
  //if( hSummary->numIsolatedNoiseChannels() >=10 ) return false;
  //if( hSummary->isolatedNoiseSumE() >=50        ) return false;
  //if( hSummary->isolatedNoiseSumEt() >=25       ) return false;
  


  edm::Handle<bool> CSCT;
  edm::InputTag  csctag("CSCTightHaloFilter");
  iEvent.getByLabel(csctag,CSCT);
  csct = (*CSCT.product()) ;
  if(false) std::cout<<" CSCT = "<<(*CSCT.product())<<std::endl;
  
  
  //ECAL Dead Cell Filter
  // Add branch for this.
  edm::Handle<bool>  ECALDeadCellTP;
  edm::InputTag ecaldeadcelltp("EcalDeadCellTriggerPrimitiveFilter");
  iEvent.getByLabel(ecaldeadcelltp,ECALDeadCellTP);
  ecaldeadcellTPFilter_ = (*ECALDeadCellTP.product()); // branch  Bool_t

  // eeBadScFilter
  // Add branch for this.
  edm::Handle<bool>  EEBADSCFILTER;
  edm::InputTag eebasscfilter("eeBadScFilter");
  iEvent.getByLabel(eebasscfilter,EEBADSCFILTER);
  eebadscf_   = (*EEBADSCFILTER.product()); // branch Bool_t
  
  // Isolation Filter Vars
  edm::Handle<HcalNoiseSummary> hSummary;
  iEvent.getByLabel("hcalnoise", hSummary);
  NIsolatedNoiseChannels_ = hSummary->numIsolatedNoiseChannels() ; // branch Float_t
  IsolatedNoiseSumE_      = hSummary->isolatedNoiseSumE(); // branch  Float_t
  IsolatedNoiseSumEt_     = hSummary->isolatedNoiseSumEt(); // branch  Float_t
  
  isolationfilterstatus_ = true; // branch Bool_t
  if( hSummary->numIsolatedNoiseChannels() >=10 ) isolationfilterstatus_  = false;
  if( hSummary->isolatedNoiseSumE() >=50        ) isolationfilterstatus_  = false;
  if( hSummary->isolatedNoiseSumEt() >=25       ) isolationfilterstatus_  = false;

  
  
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
  
  
  
  for(calojet = jetcoll.begin(); calojet != jetcoll.end(); ++calojet){
    if(false) std::cout<<" jet pt = "<<calojet->pt()
	     <<", "<<calojet->eta()
	     <<", "<<calojet->phi()
	     <<std::endl;
    if(cmet->et()<90) continue;
    if(!hbhet) continue;
    if(!csct) continue;
    if(false) std::cout<<" found event = "<<FindEvent(calojet->eta(),calojet->phi())<<std::endl;
    calojet_map->Fill(calojet->eta(),calojet->phi());
    break;
  }
  calomet_         = (Double_t) cmet->et();
  pfclustermet_    = (Double_t) pfclustermet->et();
  pfcalomet_       = (Double_t) pfcalomet->et();
  
  caloSumEt_       =  (Double_t)   cmet->sumEt(); 
  pfclusterSumEt_  =  (Double_t)   pfclustermet->sumEt(); 
  pfcalocaloSumEt_ =  (Double_t)   pfcalomet->sumEt(); 
  
  calometphi_      = (Double_t) cmet->phi();
  pfclustermetphi_ = (Double_t) pfclustermet->phi();
  pfcalometphi_    = (Double_t) pfcalomet->phi();
  
  pfcalometPhi_    = (Double_t) pfcalomet->phi();
  
  hbhet_            = (Bool_t) hbhet;
  csct_             = (Bool_t) csct ;
  hbhetR1_          = (Bool_t) hbhetr1;
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
  if(false) std::cout<<" category "<<"  "<<cmet->et()
	   <<" : "<<pfclustermet->et()
	   <<" : "<<pfcalomet->et()
	   <<" : "<<hbhet
	   <<" : "<<csct
	   <<" : "<<hbhetr1
	   <<" : "<< iEvent.id().run()<<":"
	   << iEvent.id().luminosityBlock()<<":"
	   << iEvent.id().event()
	   <<std::endl;
  for (size_t icat=0; icat<categories.size(); icat++){
    if(categories[icat]) {
      if(false) std::cout<<" category "<<icat+1<<"  "<<cmet->et()
	       <<" : "<<pfclustermet->et()
	       <<" : "<<pfcalomet->et()
	       <<" : "<<hbhet
	       <<" : "<<csct
	       <<" : "<<hbhetr1
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

    caloSumET[ifilter]->Fill(cmet->sumEt());
    PFclusterSumET[ifilter]->Fill(pfclustermet->sumEt());
    PFcaloSumET[ifilter]->Fill(pfcalomet->sumEt());

    caloMET_over_SumET[ifilter]->Fill(cmet->et()/cmet->sumEt());
    PFclusterMET_over_SumET[ifilter]->Fill(pfclustermet->et()/pfclustermet->sumEt());
    PFcaloMET_over_SumET[ifilter]->Fill(pfcalomet->et()/pfcalomet->sumEt());

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
  
  metTree->Branch("pfcalometPhi_",&pfcalometPhi_,"pfcalometPhi_/D");

  metTree->Branch("caloSumEt_",&caloSumEt_,"caloSumEt_/D");
  metTree->Branch("pfclusterSumEt_",&pfclusterSumEt_,"pfclusterSumEt_/D");
  metTree->Branch("pfcalocaloSumEt_",&pfcalocaloSumEt_,"pfcalocaloSumEt_/D");

  metTree->Branch("hbhet_",&hbhet_,"hbhet_/O");
  metTree->Branch("csct_",&csct_,"csct_/O");
  metTree->Branch("hbhetR1_",&hbhetR1_,"hbhetR1_/O");
  
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




bool MetAnalyzer::FindEvent(float eta1, float phi1){
  float eta;
  float phi;
  float deta;
  float dphi;
  bool result = false;
  //std::cout<<" inside findevent func "<<std::endl;
  for (int i =0; i<(int)etaVec.size();i++){
    //std::cout<<" inside loop "<<std::endl;
    eta = etaVec[i];  phi = phiVec[i];
    deta = eta-eta1;
    dphi = phi-phi1;
    float dr = TMath::Sqrt(deta*deta + dphi*dphi);
    if(dr<0.10) {
      std::cout<<" eta, phi,dr = "<<eta1<<", "<<phi1
	       <<", "<<eta<<", "<<phi
	       <<", "<<dr<<std::endl;
      result = true;
      break;
    }
  }  
  return result;
}





void MetAnalyzer::Setetaphi(){
etaVec.clear(); phiVec.clear();
etaVec.push_back(-1.470);  phiVec.push_back( +1.841);
etaVec.push_back(-1.470);  phiVec.push_back( +2.190);
etaVec.push_back(-1.470);  phiVec.push_back( +2.295);
etaVec.push_back(-1.470);  phiVec.push_back( +2.679);
etaVec.push_back(-1.470);  phiVec.push_back( +3.081);
etaVec.push_back(-1.470);  phiVec.push_back( +3.133);
etaVec.push_back(-1.470);  phiVec.push_back( -3.133);
etaVec.push_back(-1.470);  phiVec.push_back( -3.115);
etaVec.push_back(-1.470);  phiVec.push_back( -3.098);
etaVec.push_back(-1.470);  phiVec.push_back( -3.081);
etaVec.push_back(-1.470);  phiVec.push_back( -3.063);
etaVec.push_back(-1.470);  phiVec.push_back( -1.649);
etaVec.push_back(-1.470);  phiVec.push_back( -1.632);
etaVec.push_back(-1.470);  phiVec.push_back( -1.614);
etaVec.push_back(-1.453);  phiVec.push_back( +0.777);
etaVec.push_back(-1.453);  phiVec.push_back( +1.300);
etaVec.push_back(-1.453);  phiVec.push_back( -3.133);
etaVec.push_back(-1.453);  phiVec.push_back( -3.115);
etaVec.push_back(-1.453);  phiVec.push_back( -3.098);
etaVec.push_back(-1.453);  phiVec.push_back( -3.081);
etaVec.push_back(-1.453);  phiVec.push_back( -3.063);
etaVec.push_back(-1.436);  phiVec.push_back( -3.133);
etaVec.push_back(-1.436);  phiVec.push_back( -3.115);
etaVec.push_back(-1.436);  phiVec.push_back( -3.098);
etaVec.push_back(-1.436);  phiVec.push_back( -3.081);
etaVec.push_back(-1.436);  phiVec.push_back( -3.063);
etaVec.push_back(-1.436);  phiVec.push_back( -1.580);
etaVec.push_back(-1.418);  phiVec.push_back( -3.133);
etaVec.push_back(-1.418);  phiVec.push_back( -3.115);
etaVec.push_back(-1.418);  phiVec.push_back( -3.098);
etaVec.push_back(-1.418);  phiVec.push_back( -3.081);
etaVec.push_back(-1.418);  phiVec.push_back( -3.063);
etaVec.push_back(-1.418);  phiVec.push_back( -1.597);
etaVec.push_back(-1.401);  phiVec.push_back( -3.133);
etaVec.push_back(-1.401);  phiVec.push_back( -3.115);
etaVec.push_back(-1.401);  phiVec.push_back( -3.098);
etaVec.push_back(-1.401);  phiVec.push_back( -3.081);
etaVec.push_back(-1.401);  phiVec.push_back( -3.063);
etaVec.push_back(-1.401);  phiVec.push_back( -1.580);
etaVec.push_back(-1.383);  phiVec.push_back( +1.510);
etaVec.push_back(-1.314);  phiVec.push_back( -0.166);
etaVec.push_back(-1.314);  phiVec.push_back( -0.148);
etaVec.push_back(-1.314);  phiVec.push_back( -0.131);
etaVec.push_back(-1.314);  phiVec.push_back( -0.113);
etaVec.push_back(-1.314);  phiVec.push_back( -0.096);
etaVec.push_back(-1.296);  phiVec.push_back( -2.697);
etaVec.push_back(-1.296);  phiVec.push_back( -2.679);
etaVec.push_back(-1.296);  phiVec.push_back( -2.662);
etaVec.push_back(-1.296);  phiVec.push_back( -2.644);
etaVec.push_back(-1.296);  phiVec.push_back( -2.627);
etaVec.push_back(-1.279);  phiVec.push_back( -2.697);
etaVec.push_back(-1.279);  phiVec.push_back( -2.679);
etaVec.push_back(-1.279);  phiVec.push_back( -2.662);
etaVec.push_back(-1.279);  phiVec.push_back( -2.644);
etaVec.push_back(-1.279);  phiVec.push_back( -2.627);
etaVec.push_back(-1.279);  phiVec.push_back( -1.161);
etaVec.push_back(-1.262);  phiVec.push_back( -2.697);
etaVec.push_back(-1.262);  phiVec.push_back( -2.679);
etaVec.push_back(-1.262);  phiVec.push_back( -2.662);
etaVec.push_back(-1.262);  phiVec.push_back( -2.644);
etaVec.push_back(-1.262);  phiVec.push_back( -2.627);
etaVec.push_back(-1.209);  phiVec.push_back( +2.923);
etaVec.push_back(-1.209);  phiVec.push_back( -1.213);
etaVec.push_back(-1.209);  phiVec.push_back( -1.196);
etaVec.push_back(-1.209);  phiVec.push_back( -1.178);
etaVec.push_back(-1.209);  phiVec.push_back( -1.161);
etaVec.push_back(-1.209);  phiVec.push_back( -1.143);
etaVec.push_back(-1.192);  phiVec.push_back( -1.213);
etaVec.push_back(-1.192);  phiVec.push_back( -1.196);
etaVec.push_back(-1.192);  phiVec.push_back( -1.178);
etaVec.push_back(-1.192);  phiVec.push_back( -1.161);
etaVec.push_back(-1.192);  phiVec.push_back( -1.143);
etaVec.push_back(-1.174);  phiVec.push_back( -1.213);
etaVec.push_back(-1.174);  phiVec.push_back( -1.196);
etaVec.push_back(-1.174);  phiVec.push_back( -1.178);
etaVec.push_back(-1.174);  phiVec.push_back( -1.161);
etaVec.push_back(-1.174);  phiVec.push_back( -1.143);
etaVec.push_back(-1.157);  phiVec.push_back( -1.213);
etaVec.push_back(-1.157);  phiVec.push_back( -1.196);
etaVec.push_back(-1.157);  phiVec.push_back( -1.178);
etaVec.push_back(-1.157);  phiVec.push_back( -1.161);
etaVec.push_back(-1.157);  phiVec.push_back( -1.143);
etaVec.push_back(-1.140);  phiVec.push_back( +2.138);
etaVec.push_back(-1.140);  phiVec.push_back( -1.213);
etaVec.push_back(-1.140);  phiVec.push_back( -1.196);
etaVec.push_back(-1.140);  phiVec.push_back( -1.178);
etaVec.push_back(-1.140);  phiVec.push_back( -1.161);
etaVec.push_back(-1.140);  phiVec.push_back( -1.143);
etaVec.push_back(-1.087);  phiVec.push_back( +3.115);
etaVec.push_back(-1.070);  phiVec.push_back( -1.161);
etaVec.push_back(-1.000);  phiVec.push_back( -0.340);
etaVec.push_back(-0.983);  phiVec.push_back( +2.295);
etaVec.push_back(-0.966);  phiVec.push_back( +1.562);
etaVec.push_back(-0.948);  phiVec.push_back( +1.614);
etaVec.push_back(-0.948);  phiVec.push_back( +2.016);
etaVec.push_back(-0.948);  phiVec.push_back( +2.033);
etaVec.push_back(-0.948);  phiVec.push_back( +2.051);
etaVec.push_back(-0.948);  phiVec.push_back( +2.068);
etaVec.push_back(-0.948);  phiVec.push_back( +2.086);
etaVec.push_back(-0.948);  phiVec.push_back( +2.976);
etaVec.push_back(-0.948);  phiVec.push_back( +2.993);
etaVec.push_back(-0.948);  phiVec.push_back( +3.011);
etaVec.push_back(-0.948);  phiVec.push_back( +3.028);
etaVec.push_back(-0.948);  phiVec.push_back( +3.046);
etaVec.push_back(-0.931);  phiVec.push_back( +2.016);
etaVec.push_back(-0.931);  phiVec.push_back( +2.033);
etaVec.push_back(-0.931);  phiVec.push_back( +2.051);
etaVec.push_back(-0.931);  phiVec.push_back( +2.068);
etaVec.push_back(-0.931);  phiVec.push_back( +2.086);
etaVec.push_back(-0.931);  phiVec.push_back( +2.976);
etaVec.push_back(-0.931);  phiVec.push_back( +2.993);
etaVec.push_back(-0.931);  phiVec.push_back( +3.011);
etaVec.push_back(-0.931);  phiVec.push_back( +3.028);
etaVec.push_back(-0.931);  phiVec.push_back( +3.046);
etaVec.push_back(-0.914);  phiVec.push_back( +2.016);
etaVec.push_back(-0.914);  phiVec.push_back( +2.033);
etaVec.push_back(-0.914);  phiVec.push_back( +2.051);
etaVec.push_back(-0.914);  phiVec.push_back( +2.068);
etaVec.push_back(-0.914);  phiVec.push_back( +2.086);
etaVec.push_back(-0.914);  phiVec.push_back( +2.976);
etaVec.push_back(-0.914);  phiVec.push_back( +2.993);
etaVec.push_back(-0.914);  phiVec.push_back( +3.011);
etaVec.push_back(-0.914);  phiVec.push_back( +3.028);
etaVec.push_back(-0.914);  phiVec.push_back( +3.046);
etaVec.push_back(-0.896);  phiVec.push_back( +2.016);
etaVec.push_back(-0.896);  phiVec.push_back( +2.033);
etaVec.push_back(-0.896);  phiVec.push_back( +2.051);
etaVec.push_back(-0.896);  phiVec.push_back( +2.068);
etaVec.push_back(-0.896);  phiVec.push_back( +2.086);
etaVec.push_back(-0.896);  phiVec.push_back( +2.976);
etaVec.push_back(-0.896);  phiVec.push_back( +2.993);
etaVec.push_back(-0.896);  phiVec.push_back( +3.011);
etaVec.push_back(-0.896);  phiVec.push_back( +3.028);
etaVec.push_back(-0.896);  phiVec.push_back( +3.046);
etaVec.push_back(-0.879);  phiVec.push_back( +0.113);
etaVec.push_back(-0.879);  phiVec.push_back( +2.016);
etaVec.push_back(-0.879);  phiVec.push_back( +2.033);
etaVec.push_back(-0.879);  phiVec.push_back( +2.051);
etaVec.push_back(-0.879);  phiVec.push_back( +2.068);
etaVec.push_back(-0.879);  phiVec.push_back( +2.086);
etaVec.push_back(-0.879);  phiVec.push_back( +2.976);
etaVec.push_back(-0.879);  phiVec.push_back( +2.993);
etaVec.push_back(-0.879);  phiVec.push_back( +3.011);
etaVec.push_back(-0.879);  phiVec.push_back( +3.028);
etaVec.push_back(-0.879);  phiVec.push_back( +3.046);
etaVec.push_back(-0.844);  phiVec.push_back( +0.951);
etaVec.push_back(-0.826);  phiVec.push_back( +2.714);
etaVec.push_back(-0.774);  phiVec.push_back( +2.068);
etaVec.push_back(-0.774);  phiVec.push_back( -2.278);
etaVec.push_back(-0.774);  phiVec.push_back( -0.777);
etaVec.push_back(-0.774);  phiVec.push_back( -0.759);
etaVec.push_back(-0.774);  phiVec.push_back( -0.742);
etaVec.push_back(-0.774);  phiVec.push_back( -0.730);
etaVec.push_back(-0.774);  phiVec.push_back( -0.707);
etaVec.push_back(-0.757);  phiVec.push_back( -0.777);
etaVec.push_back(-0.757);  phiVec.push_back( -0.759);
etaVec.push_back(-0.757);  phiVec.push_back( -0.742);
etaVec.push_back(-0.757);  phiVec.push_back( -0.724);
etaVec.push_back(-0.757);  phiVec.push_back( -0.707);
etaVec.push_back(-0.739);  phiVec.push_back( +1.108);
etaVec.push_back(-0.739);  phiVec.push_back( -0.777);
etaVec.push_back(-0.739);  phiVec.push_back( -0.759);
etaVec.push_back(-0.739);  phiVec.push_back( -0.742);
etaVec.push_back(-0.739);  phiVec.push_back( -0.724);
etaVec.push_back(-0.739);  phiVec.push_back( -0.707);
etaVec.push_back(-0.722);  phiVec.push_back( -0.777);
etaVec.push_back(-0.722);  phiVec.push_back( -0.759);
etaVec.push_back(-0.722);  phiVec.push_back( -0.742);
etaVec.push_back(-0.722);  phiVec.push_back( -0.724);
etaVec.push_back(-0.722);  phiVec.push_back( -0.707);
etaVec.push_back(-0.705);  phiVec.push_back( +2.662);
etaVec.push_back(-0.705);  phiVec.push_back( -0.777);
etaVec.push_back(-0.705);  phiVec.push_back( -0.759);
etaVec.push_back(-0.705);  phiVec.push_back( -0.742);
etaVec.push_back(-0.705);  phiVec.push_back( -0.724);
etaVec.push_back(-0.705);  phiVec.push_back( -0.707);
etaVec.push_back(-0.652);  phiVec.push_back( +3.098);
etaVec.push_back(-0.652);  phiVec.push_back( -0.602);
etaVec.push_back(-0.583);  phiVec.push_back( +0.951);
etaVec.push_back(-0.583);  phiVec.push_back( +1.929);
etaVec.push_back(-0.513);  phiVec.push_back( +2.714);
etaVec.push_back(-0.513);  phiVec.push_back( +2.731);
etaVec.push_back(-0.513);  phiVec.push_back( +2.749);
etaVec.push_back(-0.513);  phiVec.push_back( +2.766);
etaVec.push_back(-0.513);  phiVec.push_back( +2.784);
etaVec.push_back(-0.513);  phiVec.push_back( -2.260);
etaVec.push_back(-0.513);  phiVec.push_back( -2.190);
etaVec.push_back(-0.513);  phiVec.push_back( -1.038);
etaVec.push_back(-0.496);  phiVec.push_back( +2.714);
etaVec.push_back(-0.496);  phiVec.push_back( +2.731);
etaVec.push_back(-0.496);  phiVec.push_back( +2.749);
etaVec.push_back(-0.496);  phiVec.push_back( +2.766);
etaVec.push_back(-0.496);  phiVec.push_back( +2.784);
etaVec.push_back(-0.479);  phiVec.push_back( +1.318);
etaVec.push_back(-0.479);  phiVec.push_back( +2.714);
etaVec.push_back(-0.479);  phiVec.push_back( +2.731);
etaVec.push_back(-0.479);  phiVec.push_back( +2.749);
etaVec.push_back(-0.479);  phiVec.push_back( +2.766);
etaVec.push_back(-0.479);  phiVec.push_back( +2.784);
etaVec.push_back(-0.461);  phiVec.push_back( +2.714);
etaVec.push_back(-0.461);  phiVec.push_back( +2.731);
etaVec.push_back(-0.461);  phiVec.push_back( +2.749);
etaVec.push_back(-0.461);  phiVec.push_back( +2.766);
etaVec.push_back(-0.461);  phiVec.push_back( +2.784);
etaVec.push_back(-0.444);  phiVec.push_back( +2.714);
etaVec.push_back(-0.444);  phiVec.push_back( +2.731);
etaVec.push_back(-0.444);  phiVec.push_back( +2.749);
etaVec.push_back(-0.444);  phiVec.push_back( +2.766);
etaVec.push_back(-0.444);  phiVec.push_back( +2.784);
etaVec.push_back(-0.426);  phiVec.push_back( +2.278);
etaVec.push_back(-0.426);  phiVec.push_back( -1.946);
etaVec.push_back(-0.409);  phiVec.push_back( +2.033);
etaVec.push_back(-0.409);  phiVec.push_back( +2.958);
etaVec.push_back(-0.391);  phiVec.push_back( +1.580);
etaVec.push_back(-0.391);  phiVec.push_back( +2.278);
etaVec.push_back(-0.357);  phiVec.push_back( +3.081);
etaVec.push_back(-0.357);  phiVec.push_back( -2.155);
etaVec.push_back(-0.339);  phiVec.push_back( +0.183);
etaVec.push_back(-0.339);  phiVec.push_back( +0.201);
etaVec.push_back(-0.339);  phiVec.push_back( +0.218);
etaVec.push_back(-0.339);  phiVec.push_back( +0.236);
etaVec.push_back(-0.339);  phiVec.push_back( +0.253);
etaVec.push_back(-0.322);  phiVec.push_back( +0.183);
etaVec.push_back(-0.322);  phiVec.push_back( +0.201);
etaVec.push_back(-0.322);  phiVec.push_back( +0.218);
etaVec.push_back(-0.322);  phiVec.push_back( +0.236);
etaVec.push_back(-0.322);  phiVec.push_back( +0.253);
etaVec.push_back(-0.322);  phiVec.push_back( -1.213);
etaVec.push_back(-0.305);  phiVec.push_back( +0.183);
etaVec.push_back(-0.305);  phiVec.push_back( +0.201);
etaVec.push_back(-0.305);  phiVec.push_back( +0.218);
etaVec.push_back(-0.305);  phiVec.push_back( +0.236);
etaVec.push_back(-0.305);  phiVec.push_back( +0.253);
etaVec.push_back(-0.287);  phiVec.push_back( +0.183);
etaVec.push_back(-0.287);  phiVec.push_back( +0.201);
etaVec.push_back(-0.287);  phiVec.push_back( +0.218);
etaVec.push_back(-0.287);  phiVec.push_back( +0.236);
etaVec.push_back(-0.287);  phiVec.push_back( +0.253);
etaVec.push_back(-0.287);  phiVec.push_back( -1.422);
etaVec.push_back(-0.270);  phiVec.push_back( +0.183);
etaVec.push_back(-0.270);  phiVec.push_back( +0.201);
etaVec.push_back(-0.270);  phiVec.push_back( +0.218);
etaVec.push_back(-0.270);  phiVec.push_back( +0.236);
etaVec.push_back(-0.270);  phiVec.push_back( +0.253);
etaVec.push_back(-0.252);  phiVec.push_back( +0.096);
etaVec.push_back(-0.252);  phiVec.push_back( +0.113);
etaVec.push_back(-0.252);  phiVec.push_back( +0.131);
etaVec.push_back(-0.252);  phiVec.push_back( +0.148);
etaVec.push_back(-0.252);  phiVec.push_back( +0.166);
etaVec.push_back(-0.235);  phiVec.push_back( +0.096);
etaVec.push_back(-0.235);  phiVec.push_back( +0.113);
etaVec.push_back(-0.235);  phiVec.push_back( +0.131);
etaVec.push_back(-0.235);  phiVec.push_back( +0.148);
etaVec.push_back(-0.235);  phiVec.push_back( +0.166);
etaVec.push_back(-0.235);  phiVec.push_back( -1.876);
etaVec.push_back(-0.218);  phiVec.push_back( +0.096);
etaVec.push_back(-0.218);  phiVec.push_back( +0.113);
etaVec.push_back(-0.218);  phiVec.push_back( +0.131);
etaVec.push_back(-0.218);  phiVec.push_back( +0.148);
etaVec.push_back(-0.218);  phiVec.push_back( +0.166);
etaVec.push_back(-0.218);  phiVec.push_back( +0.864);
etaVec.push_back(-0.200);  phiVec.push_back( +0.096);
etaVec.push_back(-0.200);  phiVec.push_back( +0.113);
etaVec.push_back(-0.200);  phiVec.push_back( +0.131);
etaVec.push_back(-0.200);  phiVec.push_back( +0.148);
etaVec.push_back(-0.200);  phiVec.push_back( +0.166);
etaVec.push_back(-0.183);  phiVec.push_back( +0.096);
etaVec.push_back(-0.183);  phiVec.push_back( +0.113);
etaVec.push_back(-0.183);  phiVec.push_back( +0.131);
etaVec.push_back(-0.183);  phiVec.push_back( +0.148);
etaVec.push_back(-0.183);  phiVec.push_back( +0.166);
etaVec.push_back(-0.183);  phiVec.push_back( -3.098);
etaVec.push_back(-0.183);  phiVec.push_back( -1.806);
etaVec.push_back(-0.165);  phiVec.push_back( +0.497);
etaVec.push_back(-0.165);  phiVec.push_back( +0.707);
etaVec.push_back(-0.165);  phiVec.push_back( +0.724);
etaVec.push_back(-0.165);  phiVec.push_back( +0.742);
etaVec.push_back(-0.165);  phiVec.push_back( +0.759);
etaVec.push_back(-0.165);  phiVec.push_back( +0.777);
etaVec.push_back(-0.165);  phiVec.push_back( -2.609);
etaVec.push_back(-0.165);  phiVec.push_back( -2.592);
etaVec.push_back(-0.165);  phiVec.push_back( -2.574);
etaVec.push_back(-0.165);  phiVec.push_back( -2.557);
etaVec.push_back(-0.165);  phiVec.push_back( -2.539);
etaVec.push_back(-0.148);  phiVec.push_back( +0.707);
etaVec.push_back(-0.148);  phiVec.push_back( +0.724);
etaVec.push_back(-0.148);  phiVec.push_back( +0.742);
etaVec.push_back(-0.148);  phiVec.push_back( +0.759);
etaVec.push_back(-0.148);  phiVec.push_back( +0.777);
etaVec.push_back(-0.148);  phiVec.push_back( -2.609);
etaVec.push_back(-0.148);  phiVec.push_back( -2.592);
etaVec.push_back(-0.148);  phiVec.push_back( -2.574);
etaVec.push_back(-0.148);  phiVec.push_back( -2.557);
etaVec.push_back(-0.148);  phiVec.push_back( -2.539);
etaVec.push_back(-0.131);  phiVec.push_back( +0.707);
etaVec.push_back(-0.131);  phiVec.push_back( +0.724);
etaVec.push_back(-0.131);  phiVec.push_back( +0.742);
etaVec.push_back(-0.131);  phiVec.push_back( +0.759);
etaVec.push_back(-0.131);  phiVec.push_back( +0.777);
etaVec.push_back(-0.131);  phiVec.push_back( -2.609);
etaVec.push_back(-0.131);  phiVec.push_back( -2.592);
etaVec.push_back(-0.131);  phiVec.push_back( -2.574);
etaVec.push_back(-0.131);  phiVec.push_back( -2.557);
etaVec.push_back(-0.131);  phiVec.push_back( -2.539);
etaVec.push_back(-0.113);  phiVec.push_back( +0.707);
etaVec.push_back(-0.113);  phiVec.push_back( +0.724);
etaVec.push_back(-0.113);  phiVec.push_back( +0.742);
etaVec.push_back(-0.113);  phiVec.push_back( +0.759);
etaVec.push_back(-0.113);  phiVec.push_back( +0.777);
etaVec.push_back(-0.113);  phiVec.push_back( -2.609);
etaVec.push_back(-0.113);  phiVec.push_back( -2.592);
etaVec.push_back(-0.113);  phiVec.push_back( -2.574);
etaVec.push_back(-0.113);  phiVec.push_back( -2.557);
etaVec.push_back(-0.113);  phiVec.push_back( -2.539);
etaVec.push_back(-0.113);  phiVec.push_back( -0.689);
etaVec.push_back(-0.096);  phiVec.push_back( +0.393);
etaVec.push_back(-0.096);  phiVec.push_back( +0.707);
etaVec.push_back(-0.096);  phiVec.push_back( +0.724);
etaVec.push_back(-0.096);  phiVec.push_back( +0.742);
etaVec.push_back(-0.096);  phiVec.push_back( +0.759);
etaVec.push_back(-0.096);  phiVec.push_back( +0.777);
etaVec.push_back(-0.096);  phiVec.push_back( +3.063);
etaVec.push_back(-0.096);  phiVec.push_back( -2.609);
etaVec.push_back(-0.096);  phiVec.push_back( -2.592);
etaVec.push_back(-0.096);  phiVec.push_back( -2.574);
etaVec.push_back(-0.096);  phiVec.push_back( -2.557);
etaVec.push_back(-0.096);  phiVec.push_back( -2.539);
etaVec.push_back(-0.078);  phiVec.push_back( +0.881);
etaVec.push_back(-0.078);  phiVec.push_back( +0.899);
etaVec.push_back(-0.078);  phiVec.push_back( +0.916);
etaVec.push_back(-0.078);  phiVec.push_back( +0.934);
etaVec.push_back(-0.078);  phiVec.push_back( +0.951);
etaVec.push_back(-0.061);  phiVec.push_back( +0.881);
etaVec.push_back(-0.061);  phiVec.push_back( +0.899);
etaVec.push_back(-0.061);  phiVec.push_back( +0.916);
etaVec.push_back(-0.061);  phiVec.push_back( +0.934);
etaVec.push_back(-0.061);  phiVec.push_back( +0.951);
etaVec.push_back(-0.043);  phiVec.push_back( +0.881);
etaVec.push_back(-0.043);  phiVec.push_back( +0.899);
etaVec.push_back(-0.043);  phiVec.push_back( +0.916);
etaVec.push_back(-0.043);  phiVec.push_back( +0.934);
etaVec.push_back(-0.043);  phiVec.push_back( +0.951);
etaVec.push_back(-0.026);  phiVec.push_back( +0.881);
etaVec.push_back(-0.026);  phiVec.push_back( +0.899);
etaVec.push_back(-0.026);  phiVec.push_back( +0.916);
etaVec.push_back(-0.026);  phiVec.push_back( +0.934);
etaVec.push_back(-0.026);  phiVec.push_back( +0.951);
etaVec.push_back(-0.026);  phiVec.push_back( -2.976);
etaVec.push_back(-0.009);  phiVec.push_back( +0.881);
etaVec.push_back(-0.009);  phiVec.push_back( +0.899);
etaVec.push_back(-0.009);  phiVec.push_back( +0.916);
etaVec.push_back(-0.009);  phiVec.push_back( +0.934);
etaVec.push_back(-0.009);  phiVec.push_back( +0.951);
etaVec.push_back(-0.009);  phiVec.push_back( +3.063);
etaVec.push_back(-0.009);  phiVec.push_back( -1.894);
etaVec.push_back(+0.061);  phiVec.push_back( -0.759);
etaVec.push_back(+0.096);  phiVec.push_back( +1.754);
etaVec.push_back(+0.096);  phiVec.push_back( +1.772);
etaVec.push_back(+0.096);  phiVec.push_back( +1.789);
etaVec.push_back(+0.096);  phiVec.push_back( +1.806);
etaVec.push_back(+0.096);  phiVec.push_back( +1.824);
etaVec.push_back(+0.096);  phiVec.push_back( +2.697);
etaVec.push_back(+0.096);  phiVec.push_back( -1.388);
etaVec.push_back(+0.096);  phiVec.push_back( -1.126);
etaVec.push_back(+0.113);  phiVec.push_back( +1.754);
etaVec.push_back(+0.113);  phiVec.push_back( +1.772);
etaVec.push_back(+0.113);  phiVec.push_back( +1.789);
etaVec.push_back(+0.113);  phiVec.push_back( +1.806);
etaVec.push_back(+0.113);  phiVec.push_back( +1.824);
etaVec.push_back(+0.113);  phiVec.push_back( -0.497);
etaVec.push_back(+0.131);  phiVec.push_back( +1.754);
etaVec.push_back(+0.131);  phiVec.push_back( +1.772);
etaVec.push_back(+0.131);  phiVec.push_back( +1.789);
etaVec.push_back(+0.131);  phiVec.push_back( +1.806);
etaVec.push_back(+0.131);  phiVec.push_back( +1.824);
etaVec.push_back(+0.131);  phiVec.push_back( +2.627);
etaVec.push_back(+0.131);  phiVec.push_back( -0.916);
etaVec.push_back(+0.148);  phiVec.push_back( +1.754);
etaVec.push_back(+0.148);  phiVec.push_back( +1.772);
etaVec.push_back(+0.148);  phiVec.push_back( +1.789);
etaVec.push_back(+0.148);  phiVec.push_back( +1.806);
etaVec.push_back(+0.148);  phiVec.push_back( +1.824);
etaVec.push_back(+0.165);  phiVec.push_back( +1.754);
etaVec.push_back(+0.165);  phiVec.push_back( +1.772);
etaVec.push_back(+0.165);  phiVec.push_back( +1.789);
etaVec.push_back(+0.165);  phiVec.push_back( +1.806);
etaVec.push_back(+0.165);  phiVec.push_back( +1.824);
etaVec.push_back(+0.287);  phiVec.push_back( -1.056);
etaVec.push_back(+0.305);  phiVec.push_back( -0.759);
etaVec.push_back(+0.339);  phiVec.push_back( +2.801);
etaVec.push_back(+0.444);  phiVec.push_back( +1.876);
etaVec.push_back(+0.444);  phiVec.push_back( -0.393);
etaVec.push_back(+0.531);  phiVec.push_back( +0.742);
etaVec.push_back(+0.566);  phiVec.push_back( -0.079);
etaVec.push_back(+0.583);  phiVec.push_back( -0.812);
etaVec.push_back(+0.583);  phiVec.push_back( -0.777);
etaVec.push_back(+0.618);  phiVec.push_back( +0.655);
etaVec.push_back(+0.635);  phiVec.push_back( +0.079);
etaVec.push_back(+0.652);  phiVec.push_back( +0.881);
etaVec.push_back(+0.687);  phiVec.push_back( -2.173);
etaVec.push_back(+0.722);  phiVec.push_back( -1.422);
etaVec.push_back(+0.739);  phiVec.push_back( -2.208);
etaVec.push_back(+0.774);  phiVec.push_back( +0.113);
etaVec.push_back(+0.774);  phiVec.push_back( +0.497);
etaVec.push_back(+0.774);  phiVec.push_back( -2.697);
etaVec.push_back(+0.774);  phiVec.push_back( -1.702);
etaVec.push_back(+0.774);  phiVec.push_back( -0.480);
etaVec.push_back(+0.792);  phiVec.push_back( -0.567);
etaVec.push_back(+0.792);  phiVec.push_back( -0.183);
etaVec.push_back(+0.879);  phiVec.push_back( +0.358);
etaVec.push_back(+0.879);  phiVec.push_back( +0.375);
etaVec.push_back(+0.879);  phiVec.push_back( +0.393);
etaVec.push_back(+0.879);  phiVec.push_back( +0.410);
etaVec.push_back(+0.879);  phiVec.push_back( +0.428);
etaVec.push_back(+0.879);  phiVec.push_back( +1.667);
etaVec.push_back(+0.879);  phiVec.push_back( +1.684);
etaVec.push_back(+0.879);  phiVec.push_back( +1.702);
etaVec.push_back(+0.879);  phiVec.push_back( +1.719);
etaVec.push_back(+0.879);  phiVec.push_back( +1.737);
etaVec.push_back(+0.879);  phiVec.push_back( +2.801);
etaVec.push_back(+0.879);  phiVec.push_back( +2.819);
etaVec.push_back(+0.879);  phiVec.push_back( +2.836);
etaVec.push_back(+0.879);  phiVec.push_back( +2.854);
etaVec.push_back(+0.879);  phiVec.push_back( +2.871);
etaVec.push_back(+0.896);  phiVec.push_back( +1.667);
etaVec.push_back(+0.896);  phiVec.push_back( +1.684);
etaVec.push_back(+0.896);  phiVec.push_back( +1.702);
etaVec.push_back(+0.896);  phiVec.push_back( +1.719);
etaVec.push_back(+0.896);  phiVec.push_back( +1.737);
etaVec.push_back(+0.896);  phiVec.push_back( +2.801);
etaVec.push_back(+0.896);  phiVec.push_back( +2.819);
etaVec.push_back(+0.896);  phiVec.push_back( +2.836);
etaVec.push_back(+0.896);  phiVec.push_back( +2.854);
etaVec.push_back(+0.896);  phiVec.push_back( +2.871);
etaVec.push_back(+0.896);  phiVec.push_back( -1.963);
etaVec.push_back(+0.914);  phiVec.push_back( +1.091);
etaVec.push_back(+0.914);  phiVec.push_back( +1.667);
etaVec.push_back(+0.914);  phiVec.push_back( +1.684);
etaVec.push_back(+0.914);  phiVec.push_back( +1.702);
etaVec.push_back(+0.914);  phiVec.push_back( +1.719);
etaVec.push_back(+0.914);  phiVec.push_back( +1.737);
etaVec.push_back(+0.914);  phiVec.push_back( +2.801);
etaVec.push_back(+0.914);  phiVec.push_back( +2.819);
etaVec.push_back(+0.914);  phiVec.push_back( +2.836);
etaVec.push_back(+0.914);  phiVec.push_back( +2.854);
etaVec.push_back(+0.914);  phiVec.push_back( +2.871);
etaVec.push_back(+0.931);  phiVec.push_back( +1.667);
etaVec.push_back(+0.931);  phiVec.push_back( +1.684);
etaVec.push_back(+0.931);  phiVec.push_back( +1.702);
etaVec.push_back(+0.931);  phiVec.push_back( +1.719);
etaVec.push_back(+0.931);  phiVec.push_back( +1.737);
etaVec.push_back(+0.931);  phiVec.push_back( +2.801);
etaVec.push_back(+0.931);  phiVec.push_back( +2.819);
etaVec.push_back(+0.931);  phiVec.push_back( +2.836);
etaVec.push_back(+0.931);  phiVec.push_back( +2.854);
etaVec.push_back(+0.931);  phiVec.push_back( +2.871);
etaVec.push_back(+0.931);  phiVec.push_back( +3.098);
etaVec.push_back(+0.948);  phiVec.push_back( +1.667);
etaVec.push_back(+0.948);  phiVec.push_back( +1.684);
etaVec.push_back(+0.948);  phiVec.push_back( +1.702);
etaVec.push_back(+0.948);  phiVec.push_back( +1.719);
etaVec.push_back(+0.948);  phiVec.push_back( +1.737);
etaVec.push_back(+0.948);  phiVec.push_back( +2.801);
etaVec.push_back(+0.948);  phiVec.push_back( +2.819);
etaVec.push_back(+0.948);  phiVec.push_back( +2.836);
etaVec.push_back(+0.948);  phiVec.push_back( +2.854);
etaVec.push_back(+0.948);  phiVec.push_back( +2.871);
etaVec.push_back(+0.966);  phiVec.push_back( +0.620);
etaVec.push_back(+0.966);  phiVec.push_back( +0.637);
etaVec.push_back(+0.966);  phiVec.push_back( +0.655);
etaVec.push_back(+0.966);  phiVec.push_back( +0.672);
etaVec.push_back(+0.966);  phiVec.push_back( +0.689);
etaVec.push_back(+0.983);  phiVec.push_back( +0.620);
etaVec.push_back(+0.983);  phiVec.push_back( +0.637);
etaVec.push_back(+0.983);  phiVec.push_back( +0.655);
etaVec.push_back(+0.983);  phiVec.push_back( +0.672);
etaVec.push_back(+0.983);  phiVec.push_back( +0.689);
etaVec.push_back(+1.000);  phiVec.push_back( +0.620);
etaVec.push_back(+1.000);  phiVec.push_back( +0.637);
etaVec.push_back(+1.000);  phiVec.push_back( +0.655);
etaVec.push_back(+1.000);  phiVec.push_back( +0.672);
etaVec.push_back(+1.000);  phiVec.push_back( +0.689);
etaVec.push_back(+1.018);  phiVec.push_back( +0.620);
etaVec.push_back(+1.018);  phiVec.push_back( +0.637);
etaVec.push_back(+1.018);  phiVec.push_back( +0.655);
etaVec.push_back(+1.018);  phiVec.push_back( +0.672);
etaVec.push_back(+1.018);  phiVec.push_back( +0.689);
etaVec.push_back(+1.018);  phiVec.push_back( +2.417);
etaVec.push_back(+1.035);  phiVec.push_back( +0.620);
etaVec.push_back(+1.035);  phiVec.push_back( +0.637);
etaVec.push_back(+1.035);  phiVec.push_back( +0.655);
etaVec.push_back(+1.035);  phiVec.push_back( +0.672);
etaVec.push_back(+1.035);  phiVec.push_back( +0.689);
etaVec.push_back(+1.053);  phiVec.push_back( -3.133);
etaVec.push_back(+1.053);  phiVec.push_back( -3.115);
etaVec.push_back(+1.053);  phiVec.push_back( -3.098);
etaVec.push_back(+1.053);  phiVec.push_back( -3.081);
etaVec.push_back(+1.053);  phiVec.push_back( -3.063);
etaVec.push_back(+1.070);  phiVec.push_back( -3.133);
etaVec.push_back(+1.070);  phiVec.push_back( -3.115);
etaVec.push_back(+1.070);  phiVec.push_back( -3.098);
etaVec.push_back(+1.070);  phiVec.push_back( -3.081);
etaVec.push_back(+1.070);  phiVec.push_back( -3.063);
etaVec.push_back(+1.087);  phiVec.push_back( -3.133);
etaVec.push_back(+1.087);  phiVec.push_back( -3.115);
etaVec.push_back(+1.087);  phiVec.push_back( -3.098);
etaVec.push_back(+1.087);  phiVec.push_back( -3.081);
etaVec.push_back(+1.087);  phiVec.push_back( -3.063);
etaVec.push_back(+1.105);  phiVec.push_back( -3.133);
etaVec.push_back(+1.105);  phiVec.push_back( -3.115);
etaVec.push_back(+1.105);  phiVec.push_back( -3.098);
etaVec.push_back(+1.105);  phiVec.push_back( -3.081);
etaVec.push_back(+1.105);  phiVec.push_back( -3.063);
etaVec.push_back(+1.122);  phiVec.push_back( +3.063);
etaVec.push_back(+1.122);  phiVec.push_back( -3.133);
etaVec.push_back(+1.122);  phiVec.push_back( -3.115);
etaVec.push_back(+1.122);  phiVec.push_back( -3.098);
etaVec.push_back(+1.122);  phiVec.push_back( -3.081);
etaVec.push_back(+1.122);  phiVec.push_back( -3.063);
etaVec.push_back(+1.122);  phiVec.push_back( -2.871);
etaVec.push_back(+1.140);  phiVec.push_back( +0.026);
etaVec.push_back(+1.140);  phiVec.push_back( +2.365);
etaVec.push_back(+1.140);  phiVec.push_back( -0.846);
etaVec.push_back(+1.157);  phiVec.push_back( -1.230);
etaVec.push_back(+1.174);  phiVec.push_back( +1.667);
etaVec.push_back(+1.174);  phiVec.push_back( -1.894);
etaVec.push_back(+1.192);  phiVec.push_back( +1.091);
etaVec.push_back(+1.227);  phiVec.push_back( +2.435);
etaVec.push_back(+1.227);  phiVec.push_back( -0.515);
etaVec.push_back(+1.227);  phiVec.push_back( -0.497);
etaVec.push_back(+1.227);  phiVec.push_back( -0.480);
etaVec.push_back(+1.227);  phiVec.push_back( -0.463);
etaVec.push_back(+1.227);  phiVec.push_back( -0.445);
etaVec.push_back(+1.244);  phiVec.push_back( +1.649);
etaVec.push_back(+1.244);  phiVec.push_back( -0.986);
etaVec.push_back(+1.244);  phiVec.push_back( -0.515);
etaVec.push_back(+1.244);  phiVec.push_back( -0.497);
etaVec.push_back(+1.244);  phiVec.push_back( -0.480);
etaVec.push_back(+1.244);  phiVec.push_back( -0.463);
etaVec.push_back(+1.244);  phiVec.push_back( -0.445);
etaVec.push_back(+1.262);  phiVec.push_back( -0.515);
etaVec.push_back(+1.262);  phiVec.push_back( -0.497);
etaVec.push_back(+1.262);  phiVec.push_back( -0.480);
etaVec.push_back(+1.262);  phiVec.push_back( -0.463);
etaVec.push_back(+1.262);  phiVec.push_back( -0.445);
etaVec.push_back(+1.279);  phiVec.push_back( -0.515);
etaVec.push_back(+1.279);  phiVec.push_back( -0.497);
etaVec.push_back(+1.279);  phiVec.push_back( -0.480);
etaVec.push_back(+1.279);  phiVec.push_back( -0.463);
etaVec.push_back(+1.279);  phiVec.push_back( -0.445);
etaVec.push_back(+1.296);  phiVec.push_back( -0.515);
etaVec.push_back(+1.296);  phiVec.push_back( -0.497);
etaVec.push_back(+1.296);  phiVec.push_back( -0.480);
etaVec.push_back(+1.296);  phiVec.push_back( -0.463);
etaVec.push_back(+1.296);  phiVec.push_back( -0.445);
etaVec.push_back(+1.296);  phiVec.push_back( -0.288);
etaVec.push_back(+1.314);  phiVec.push_back( -0.201);
etaVec.push_back(+1.401);  phiVec.push_back( +2.714);
etaVec.push_back(+1.401);  phiVec.push_back( +2.731);
etaVec.push_back(+1.401);  phiVec.push_back( +2.749);
etaVec.push_back(+1.401);  phiVec.push_back( +2.766);
etaVec.push_back(+1.401);  phiVec.push_back( +2.784);
etaVec.push_back(+1.401);  phiVec.push_back( -2.609);
etaVec.push_back(+1.401);  phiVec.push_back( -2.592);
etaVec.push_back(+1.401);  phiVec.push_back( -2.574);
etaVec.push_back(+1.401);  phiVec.push_back( -2.557);
etaVec.push_back(+1.401);  phiVec.push_back( -2.539);
etaVec.push_back(+1.401);  phiVec.push_back( -0.253);
etaVec.push_back(+1.401);  phiVec.push_back( -0.236);
etaVec.push_back(+1.401);  phiVec.push_back( -0.218);
etaVec.push_back(+1.401);  phiVec.push_back( -0.201);
etaVec.push_back(+1.401);  phiVec.push_back( -0.183);
etaVec.push_back(+1.418);  phiVec.push_back( +0.375);
etaVec.push_back(+1.418);  phiVec.push_back( +2.714);
etaVec.push_back(+1.418);  phiVec.push_back( +2.731);
etaVec.push_back(+1.418);  phiVec.push_back( +2.749);
etaVec.push_back(+1.418);  phiVec.push_back( +2.766);
etaVec.push_back(+1.418);  phiVec.push_back( +2.784);
etaVec.push_back(+1.418);  phiVec.push_back( -2.609);
etaVec.push_back(+1.418);  phiVec.push_back( -2.592);
etaVec.push_back(+1.418);  phiVec.push_back( -2.574);
etaVec.push_back(+1.418);  phiVec.push_back( -2.557);
etaVec.push_back(+1.418);  phiVec.push_back( -2.539);
etaVec.push_back(+1.418);  phiVec.push_back( -0.253);
etaVec.push_back(+1.418);  phiVec.push_back( -0.236);
etaVec.push_back(+1.418);  phiVec.push_back( -0.218);
etaVec.push_back(+1.418);  phiVec.push_back( -0.201);
etaVec.push_back(+1.418);  phiVec.push_back( -0.183);
etaVec.push_back(+1.436);  phiVec.push_back( +2.714);
etaVec.push_back(+1.436);  phiVec.push_back( +2.731);
etaVec.push_back(+1.436);  phiVec.push_back( +2.749);
etaVec.push_back(+1.436);  phiVec.push_back( +2.766);
etaVec.push_back(+1.436);  phiVec.push_back( +2.784);
etaVec.push_back(+1.436);  phiVec.push_back( -2.609);
etaVec.push_back(+1.436);  phiVec.push_back( -2.592);
etaVec.push_back(+1.436);  phiVec.push_back( -2.574);
etaVec.push_back(+1.436);  phiVec.push_back( -2.557);
etaVec.push_back(+1.436);  phiVec.push_back( -2.539);
etaVec.push_back(+1.436);  phiVec.push_back( -0.253);
etaVec.push_back(+1.436);  phiVec.push_back( -0.236);
etaVec.push_back(+1.436);  phiVec.push_back( -0.218);
etaVec.push_back(+1.436);  phiVec.push_back( -0.201);
etaVec.push_back(+1.436);  phiVec.push_back( -0.183);
etaVec.push_back(+1.453);  phiVec.push_back( +2.714);
etaVec.push_back(+1.453);  phiVec.push_back( +2.731);
etaVec.push_back(+1.453);  phiVec.push_back( +2.749);
etaVec.push_back(+1.453);  phiVec.push_back( +2.766);
etaVec.push_back(+1.453);  phiVec.push_back( +2.784);
etaVec.push_back(+1.453);  phiVec.push_back( -0.253);
etaVec.push_back(+1.453);  phiVec.push_back( -0.236);
etaVec.push_back(+1.453);  phiVec.push_back( -0.218);
etaVec.push_back(+1.453);  phiVec.push_back( -0.201);
etaVec.push_back(+1.453);  phiVec.push_back( -0.183);
etaVec.push_back(+1.470);  phiVec.push_back( +2.714);
etaVec.push_back(+1.470);  phiVec.push_back( +2.731);
etaVec.push_back(+1.470);  phiVec.push_back( +2.749);
etaVec.push_back(+1.470);  phiVec.push_back( +2.766);
etaVec.push_back(+1.470);  phiVec.push_back( +2.784);
etaVec.push_back(+1.470);  phiVec.push_back( +2.801);
etaVec.push_back(+1.470);  phiVec.push_back( -0.253);
etaVec.push_back(+1.470);  phiVec.push_back( -0.236);
etaVec.push_back(+1.470);  phiVec.push_back( -0.218);
etaVec.push_back(+1.470);  phiVec.push_back( -0.201);
etaVec.push_back(+1.470);  phiVec.push_back( -0.183);
etaVec.push_back(+1.507);  phiVec.push_back( -3.092);
etaVec.push_back(+1.526);  phiVec.push_back( -3.091);
etaVec.push_back(+1.546);  phiVec.push_back( -3.090);
etaVec.push_back(+1.566);  phiVec.push_back( -3.089);
etaVec.push_back(+1.587);  phiVec.push_back( -3.088);
etaVec.push_back(+1.574);  phiVec.push_back( -2.871);
etaVec.push_back(+1.595);  phiVec.push_back( -2.865);
etaVec.push_back(+1.615);  phiVec.push_back( -2.860);
etaVec.push_back(+1.637);  phiVec.push_back( -2.854);
etaVec.push_back(+1.658);  phiVec.push_back( -2.847);
etaVec.push_back(+1.597);  phiVec.push_back( +2.679);
etaVec.push_back(+1.588);  phiVec.push_back( +2.660);
etaVec.push_back(+1.899);  phiVec.push_back( +2.757);
etaVec.push_back(+1.996);  phiVec.push_back( +2.748);
etaVec.push_back(+1.555);  phiVec.push_back( -2.162);
etaVec.push_back(+1.674);  phiVec.push_back( +2.256);
etaVec.push_back(+1.584);  phiVec.push_back( -2.157);
etaVec.push_back(+1.854);  phiVec.push_back( +2.401);
etaVec.push_back(+1.689);  phiVec.push_back( +2.238);
etaVec.push_back(+2.307);  phiVec.push_back( -3.027);
etaVec.push_back(+1.563);  phiVec.push_back( +2.020);
etaVec.push_back(+2.113);  phiVec.push_back( +2.198);
etaVec.push_back(+2.107);  phiVec.push_back( +1.978);
etaVec.push_back(+2.835);  phiVec.push_back( +2.413);
etaVec.push_back(+2.782);  phiVec.push_back( +2.357);
etaVec.push_back(+1.522);  phiVec.push_back( +1.703);
etaVec.push_back(+1.525);  phiVec.push_back( +1.662);
etaVec.push_back(+1.609);  phiVec.push_back( -1.605);
etaVec.push_back(+1.695);  phiVec.push_back( +1.462);
etaVec.push_back(+1.505);  phiVec.push_back( -1.461);
etaVec.push_back(+1.843);  phiVec.push_back( -1.412);
etaVec.push_back(+2.164);  phiVec.push_back( -1.348);
etaVec.push_back(+2.884);  phiVec.push_back( -0.900);
etaVec.push_back(+2.835);  phiVec.push_back( -0.838);
etaVec.push_back(+2.891);  phiVec.push_back( -0.784);
etaVec.push_back(+2.782);  phiVec.push_back( -0.785);
etaVec.push_back(+2.835);  phiVec.push_back( -0.729);
etaVec.push_back(+2.052);  phiVec.push_back( +1.229);
etaVec.push_back(+1.820);  phiVec.push_back( +0.961);
etaVec.push_back(+1.567);  phiVec.push_back( -0.996);
etaVec.push_back(+1.817);  phiVec.push_back( +0.782);
etaVec.push_back(+1.798);  phiVec.push_back( +0.801);
etaVec.push_back(+1.780);  phiVec.push_back( +0.820);
etaVec.push_back(+1.761);  phiVec.push_back( +0.837);
etaVec.push_back(+1.743);  phiVec.push_back( +0.855);
etaVec.push_back(+1.551);  phiVec.push_back( +1.009);
etaVec.push_back(+1.798);  phiVec.push_back( +0.764);
etaVec.push_back(+1.780);  phiVec.push_back( +0.783);
etaVec.push_back(+1.762);  phiVec.push_back( +0.801);
etaVec.push_back(+1.744);  phiVec.push_back( +0.819);
etaVec.push_back(+1.726);  phiVec.push_back( +0.836);
etaVec.push_back(+1.779);  phiVec.push_back( +0.746);
etaVec.push_back(+1.761);  phiVec.push_back( +0.765);
etaVec.push_back(+1.744);  phiVec.push_back( +0.784);
etaVec.push_back(+1.727);  phiVec.push_back( +0.801);
etaVec.push_back(+1.709);  phiVec.push_back( +0.818);
etaVec.push_back(+1.760);  phiVec.push_back( +0.729);
etaVec.push_back(+1.743);  phiVec.push_back( +0.748);
etaVec.push_back(+1.726);  phiVec.push_back( +0.766);
etaVec.push_back(+1.709);  phiVec.push_back( +0.784);
etaVec.push_back(+1.693);  phiVec.push_back( +0.801);
etaVec.push_back(+1.740);  phiVec.push_back( +0.713);
etaVec.push_back(+1.724);  phiVec.push_back( +0.731);
etaVec.push_back(+1.708);  phiVec.push_back( +0.750);
etaVec.push_back(+1.692);  phiVec.push_back( +0.767);
etaVec.push_back(+1.676);  phiVec.push_back( +0.785);
etaVec.push_back(+1.663);  phiVec.push_back( -0.764);
etaVec.push_back(+1.678);  phiVec.push_back( -0.747);
etaVec.push_back(+1.694);  phiVec.push_back( -0.729);
etaVec.push_back(+1.709);  phiVec.push_back( -0.711);
etaVec.push_back(+1.724);  phiVec.push_back( -0.692);
etaVec.push_back(+1.739);  phiVec.push_back( -0.671);
etaVec.push_back(+1.754);  phiVec.push_back( -0.651);
etaVec.push_back(+1.769);  phiVec.push_back( -0.631);
etaVec.push_back(+1.783);  phiVec.push_back( -0.609);
etaVec.push_back(+1.798);  phiVec.push_back( -0.587);
etaVec.push_back(+1.646);  phiVec.push_back( -0.749);
etaVec.push_back(+1.661);  phiVec.push_back( -0.732);
etaVec.push_back(+1.676);  phiVec.push_back( -0.714);
etaVec.push_back(+1.691);  phiVec.push_back( -0.696);
etaVec.push_back(+1.706);  phiVec.push_back( -0.677);
etaVec.push_back(+1.719);  phiVec.push_back( -0.656);
etaVec.push_back(+1.734);  phiVec.push_back( -0.636);
etaVec.push_back(+1.748);  phiVec.push_back( -0.616);
etaVec.push_back(+1.762);  phiVec.push_back( -0.595);
etaVec.push_back(+1.776);  phiVec.push_back( -0.573);
etaVec.push_back(+1.630);  phiVec.push_back( -0.734);
etaVec.push_back(+1.644);  phiVec.push_back( -0.717);
etaVec.push_back(+1.658);  phiVec.push_back( -0.699);
etaVec.push_back(+1.673);  phiVec.push_back( -0.681);
etaVec.push_back(+1.687);  phiVec.push_back( -0.663);
etaVec.push_back(+1.700);  phiVec.push_back( -0.642);
etaVec.push_back(+1.714);  phiVec.push_back( -0.622);
etaVec.push_back(+1.727);  phiVec.push_back( -0.602);
etaVec.push_back(+1.741);  phiVec.push_back( -0.581);
etaVec.push_back(+1.754);  phiVec.push_back( -0.560);
etaVec.push_back(+1.613);  phiVec.push_back( -0.719);
etaVec.push_back(+1.627);  phiVec.push_back( -0.702);
etaVec.push_back(+1.641);  phiVec.push_back( -0.685);
etaVec.push_back(+1.655);  phiVec.push_back( -0.667);
etaVec.push_back(+1.668);  phiVec.push_back( -0.649);
etaVec.push_back(+1.681);  phiVec.push_back( -0.628);
etaVec.push_back(+1.694);  phiVec.push_back( -0.609);
etaVec.push_back(+1.707);  phiVec.push_back( -0.589);
etaVec.push_back(+1.720);  phiVec.push_back( -0.568);
etaVec.push_back(+1.733);  phiVec.push_back( -0.547);
etaVec.push_back(+1.596);  phiVec.push_back( -0.705);
etaVec.push_back(+1.610);  phiVec.push_back( -0.689);
etaVec.push_back(+1.623);  phiVec.push_back( -0.671);
etaVec.push_back(+1.637);  phiVec.push_back( -0.654);
etaVec.push_back(+1.650);  phiVec.push_back( -0.635);
etaVec.push_back(+1.662);  phiVec.push_back( -0.615);
etaVec.push_back(+1.674);  phiVec.push_back( -0.595);
etaVec.push_back(+1.687);  phiVec.push_back( -0.576);
etaVec.push_back(+1.700);  phiVec.push_back( -0.555);
etaVec.push_back(+1.712);  phiVec.push_back( -0.535);
etaVec.push_back(+1.849);  phiVec.push_back( -0.072);
etaVec.push_back(+1.520);  phiVec.push_back( -0.766);
etaVec.push_back(+1.533);  phiVec.push_back( -0.751);
etaVec.push_back(+1.546);  phiVec.push_back( -0.736);
etaVec.push_back(+1.559);  phiVec.push_back( -0.721);
etaVec.push_back(+1.572);  phiVec.push_back( -0.705);
etaVec.push_back(+1.584);  phiVec.push_back( -0.687);
etaVec.push_back(+1.597);  phiVec.push_back( -0.671);
etaVec.push_back(+1.610);  phiVec.push_back( -0.653);
etaVec.push_back(+1.622);  phiVec.push_back( -0.636);
etaVec.push_back(+1.635);  phiVec.push_back( -0.618);
etaVec.push_back(+1.506);  phiVec.push_back( -0.753);
etaVec.push_back(+1.519);  phiVec.push_back( -0.738);
etaVec.push_back(+1.531);  phiVec.push_back( -0.723);
etaVec.push_back(+1.544);  phiVec.push_back( -0.708);
etaVec.push_back(+1.557);  phiVec.push_back( -0.692);
etaVec.push_back(+1.568);  phiVec.push_back( -0.675);
etaVec.push_back(+1.580);  phiVec.push_back( -0.658);
etaVec.push_back(+1.592);  phiVec.push_back( -0.641);
etaVec.push_back(+1.605);  phiVec.push_back( -0.623);
etaVec.push_back(+1.617);  phiVec.push_back( -0.605);
etaVec.push_back(+1.551);  phiVec.push_back( -0.662);
etaVec.push_back(+1.563);  phiVec.push_back( -0.646);
etaVec.push_back(+1.575);  phiVec.push_back( -0.629);
etaVec.push_back(+1.587);  phiVec.push_back( -0.611);
etaVec.push_back(+1.599);  phiVec.push_back( -0.593);
etaVec.push_back(+1.535);  phiVec.push_back( -0.650);
etaVec.push_back(+1.547);  phiVec.push_back( -0.634);
etaVec.push_back(+1.558);  phiVec.push_back( -0.617);
etaVec.push_back(+1.570);  phiVec.push_back( -0.600);
etaVec.push_back(+1.581);  phiVec.push_back( -0.582);
etaVec.push_back(+1.519);  phiVec.push_back( -0.638);
etaVec.push_back(+1.530);  phiVec.push_back( -0.622);
etaVec.push_back(+1.542);  phiVec.push_back( -0.605);
etaVec.push_back(+1.553);  phiVec.push_back( -0.588);
etaVec.push_back(+1.564);  phiVec.push_back( -0.571);
etaVec.push_back(+1.508);  phiVec.push_back( -0.622);
etaVec.push_back(+1.518);  phiVec.push_back( -0.606);
etaVec.push_back(+1.529);  phiVec.push_back( -0.590);
etaVec.push_back(+1.539);  phiVec.push_back( -0.573);
etaVec.push_back(+1.550);  phiVec.push_back( -0.555);
etaVec.push_back(+1.605);  phiVec.push_back( -0.441);
etaVec.push_back(+1.492);  phiVec.push_back( -0.611);
etaVec.push_back(+1.502);  phiVec.push_back( -0.595);
etaVec.push_back(+1.512);  phiVec.push_back( -0.579);
etaVec.push_back(+1.523);  phiVec.push_back( -0.562);
etaVec.push_back(+1.533);  phiVec.push_back( -0.545);
etaVec.push_back(-1.495);  phiVec.push_back( -2.975);
etaVec.push_back(-1.508);  phiVec.push_back( -3.131);
etaVec.push_back(-1.530);  phiVec.push_back( -2.948);
etaVec.push_back(-1.563);  phiVec.push_back( +2.764);
etaVec.push_back(-1.825);  phiVec.push_back( -3.127);
etaVec.push_back(-1.957);  phiVec.push_back( -2.965);
etaVec.push_back(-1.744);  phiVec.push_back( -2.323);
etaVec.push_back(-1.743);  phiVec.push_back( -2.287);
etaVec.push_back(-1.584);  phiVec.push_back( +2.157);
etaVec.push_back(-1.567);  phiVec.push_back( +2.145);
etaVec.push_back(-1.950);  phiVec.push_back( -2.471);
etaVec.push_back(-1.969);  phiVec.push_back( -2.496);
etaVec.push_back(-1.987);  phiVec.push_back( -2.522);
etaVec.push_back(-2.006);  phiVec.push_back( -2.549);
etaVec.push_back(-2.024);  phiVec.push_back( -2.578);
etaVec.push_back(-1.974);  phiVec.push_back( -2.451);
etaVec.push_back(-1.994);  phiVec.push_back( -2.477);
etaVec.push_back(-2.014);  phiVec.push_back( -2.503);
etaVec.push_back(-2.034);  phiVec.push_back( -2.531);
etaVec.push_back(-2.053);  phiVec.push_back( -2.559);
etaVec.push_back(-1.999);  phiVec.push_back( -2.430);
etaVec.push_back(-2.020);  phiVec.push_back( -2.456);
etaVec.push_back(-2.041);  phiVec.push_back( -2.483);
etaVec.push_back(-2.062);  phiVec.push_back( -2.511);
etaVec.push_back(-2.082);  phiVec.push_back( -2.540);
etaVec.push_back(-2.024);  phiVec.push_back( -2.408);
etaVec.push_back(-2.047);  phiVec.push_back( -2.434);
etaVec.push_back(-2.069);  phiVec.push_back( -2.461);
etaVec.push_back(-2.090);  phiVec.push_back( -2.489);
etaVec.push_back(-2.112);  phiVec.push_back( -2.519);
etaVec.push_back(-2.050);  phiVec.push_back( -2.385);
etaVec.push_back(-2.073);  phiVec.push_back( -2.411);
etaVec.push_back(-2.096);  phiVec.push_back( -2.439);
etaVec.push_back(-2.119);  phiVec.push_back( -2.467);
etaVec.push_back(-2.142);  phiVec.push_back( -2.497);
etaVec.push_back(-1.820);  phiVec.push_back( +2.180);
etaVec.push_back(-2.063);  phiVec.push_back( -2.205);
etaVec.push_back(-1.894);  phiVec.push_back( -1.893);
etaVec.push_back(-1.513);  phiVec.push_back( -1.761);
etaVec.push_back(-1.505);  phiVec.push_back( -1.681);
etaVec.push_back(-2.061);  phiVec.push_back( -1.732);
etaVec.push_back(-2.108);  phiVec.push_back( -1.628);
etaVec.push_back(-1.853);  phiVec.push_back( +1.587);
etaVec.push_back(-2.187);  phiVec.push_back( -1.549);
etaVec.push_back(-2.228);  phiVec.push_back( -1.548);
etaVec.push_back(-2.270);  phiVec.push_back( -1.547);
etaVec.push_back(-2.315);  phiVec.push_back( -1.545);
etaVec.push_back(-2.361);  phiVec.push_back( -1.544);
etaVec.push_back(-1.722);  phiVec.push_back( -1.532);
etaVec.push_back(-2.185);  phiVec.push_back( -1.508);
etaVec.push_back(-2.226);  phiVec.push_back( -1.505);
etaVec.push_back(-2.268);  phiVec.push_back( -1.503);
etaVec.push_back(-2.312);  phiVec.push_back( -1.500);
etaVec.push_back(-2.359);  phiVec.push_back( -1.497);
etaVec.push_back(-2.182);  phiVec.push_back( -1.469);
etaVec.push_back(-2.222);  phiVec.push_back( -1.464);
etaVec.push_back(-2.264);  phiVec.push_back( -1.460);
etaVec.push_back(-2.308);  phiVec.push_back( -1.455);
etaVec.push_back(-2.354);  phiVec.push_back( -1.450);
etaVec.push_back(-2.177);  phiVec.push_back( -1.429);
etaVec.push_back(-2.217);  phiVec.push_back( -1.423);
etaVec.push_back(-2.259);  phiVec.push_back( -1.417);
etaVec.push_back(-2.302);  phiVec.push_back( -1.410);
etaVec.push_back(-2.348);  phiVec.push_back( -1.403);
etaVec.push_back(-2.171);  phiVec.push_back( -1.390);
etaVec.push_back(-2.210);  phiVec.push_back( -1.383);
etaVec.push_back(-2.251);  phiVec.push_back( -1.375);
etaVec.push_back(-2.294);  phiVec.push_back( -1.366);
etaVec.push_back(-2.339);  phiVec.push_back( -1.357);
etaVec.push_back(-2.542);  phiVec.push_back( -0.906);
etaVec.push_back(-2.506);  phiVec.push_back( +0.706);
etaVec.push_back(-1.743);  phiVec.push_back( -0.855);
etaVec.push_back(-1.761);  phiVec.push_back( -0.837);
etaVec.push_back(-1.551);  phiVec.push_back( +1.009);
etaVec.push_back(-1.572);  phiVec.push_back( +0.967);
etaVec.push_back(-1.555);  phiVec.push_back( +0.979);
etaVec.push_back(-2.072);  phiVec.push_back( +0.019);
etaVec.push_back(-1.962);  phiVec.push_back( +0.144);
etaVec.push_back(-1.635);  phiVec.push_back( +0.618);
etaVec.push_back(-1.622);  phiVec.push_back( +0.636);
etaVec.push_back(-1.610);  phiVec.push_back( +0.653);
etaVec.push_back(-1.597);  phiVec.push_back( +0.671);
etaVec.push_back(-1.584);  phiVec.push_back( +0.687);
etaVec.push_back(-1.617);  phiVec.push_back( +0.605);
etaVec.push_back(-1.605);  phiVec.push_back( +0.623);
etaVec.push_back(-1.592);  phiVec.push_back( +0.641);
etaVec.push_back(-1.580);  phiVec.push_back( +0.658);
etaVec.push_back(-1.568);  phiVec.push_back( +0.675);
etaVec.push_back(-1.599);  phiVec.push_back( +0.593);
etaVec.push_back(-1.587);  phiVec.push_back( +0.611);
etaVec.push_back(-1.575);  phiVec.push_back( +0.629);
etaVec.push_back(-1.563);  phiVec.push_back( +0.646);
etaVec.push_back(-1.551);  phiVec.push_back( +0.662);
etaVec.push_back(-1.581);  phiVec.push_back( +0.582);
etaVec.push_back(-1.570);  phiVec.push_back( +0.600);
etaVec.push_back(-1.558);  phiVec.push_back( +0.617);
etaVec.push_back(-1.547);  phiVec.push_back( +0.634);
etaVec.push_back(-1.535);  phiVec.push_back( +0.650);
etaVec.push_back(-1.564);  phiVec.push_back( +0.571);
etaVec.push_back(-1.553);  phiVec.push_back( +0.588);
etaVec.push_back(-1.542);  phiVec.push_back( +0.605);
etaVec.push_back(-1.530);  phiVec.push_back( +0.622);
etaVec.push_back(-1.519);  phiVec.push_back( +0.638);
etaVec.push_back(-1.522);  phiVec.push_back( +0.112);
etaVec.push_back(-1.504);  phiVec.push_back( -0.089);
}



//define this as a plug-in
DEFINE_FWK_MODULE(MetAnalyzer);
