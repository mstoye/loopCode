//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Feb 27 16:22:54 2014 by ROOT version 5.32/00
// from TTree SUSYLooperHistsSoftBase/SUSYLooperHistsSoftBase
// found on file: /afs/cern.ch/work/a/agapitos/public/CMG/CMSSW_5_3_12_patch3/src/CMGTools/TTHAnalysis/cfg/RPV_1200_500/AtoniousSUSY/SUSYLooperHistsSoftBase/SUSYLooperHistsSoftBase_tree.root
//////////////////////////////////////////////////////////

#ifndef A7LoopBase_h
#define A7LoopBase_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "SampleInfo.h" 
#include <iostream>
#include <TLorentzVector.h>


// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class A7LoopBase {

 public :
  TString outputFileName;
  float weight;
  unsigned int maxEvents;
 

public :

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types

 Int_t           run;
   Int_t           lumi;
   Int_t           evt;
   Int_t           isData;  
   Int_t           HLT_SingleMu;
   Float_t         rho;
   Int_t           nVert;
   Int_t           nJet25;
   Int_t           nBJetLoose25;
   Int_t           nBJetMedium25;
   Int_t           nBJetTight25;
   Int_t           nJet40;
   Int_t           nJet40a;
   Int_t           nBJetLoose40;
   Int_t           nBJetMedium40;
   Int_t           nBJetTight40;
   Int_t           nLepGood20;
   Int_t           nLepGood15;
   Int_t           nLepGood10;
   Int_t           checkEcalDead;
   Int_t           checkhcalLaser;
   Int_t           checktrackingFailure;
   Int_t           checkprimaryVertex;
   Int_t           checknoscraping;
   Int_t           checktrackIsolationMaker;
   Int_t           checkmetNoiseCleaning;
   Int_t           checkeeBadSc;
   Int_t           checkecalLaser;
   Int_t           checktotalKinematics;
   Int_t           checkCSCTightHalo;
   Int_t           checkHBHENoise;
   Float_t         htJet25;
   Float_t         mhtJet25;
   Float_t         htJet40j;
   Float_t         htJet40ja;
   Float_t         htJet40;
   Float_t         htJet40a;
   Float_t         mhtJet40;
   Float_t         mhtJet40a;
   Float_t         mZ1;
   Float_t         mZ1SFSS;
   Float_t         minMllSFOS;
   Float_t         maxMllSFOS;
   Float_t         minMllAFOS;
   Float_t         maxMllAFOS;
   Float_t         minMllAFSS;
   Float_t         maxMllAFSS;
   Float_t         minMllAFAS;
   Float_t         maxMllAFAS;
   Float_t         m2l;
   Int_t           minMWjj;
   Int_t           minMWjjPt;
   Int_t           bestMWjj;
   Int_t           bestMWjjPt;
   Int_t           bestMTopHad;
   Int_t           bestMTopHadPt;
   Float_t         met_pt;
   Float_t         met_eta;
   Float_t         met_phi;
   Float_t         met_mass;
   Float_t         met_sumEt;
   Float_t         metNoPU_pt;
   Float_t         metNoPU_eta;
   Float_t         metNoPU_phi;
   Float_t         metNoPU_mass;
   Int_t           nLepOther;
   Float_t         LepOther_pt[6];   //[nLepOther]
   Float_t         LepOther_eta[6];   //[nLepOther]
   Float_t         LepOther_phi[6];   //[nLepOther]
   Float_t         LepOther_mass[6];   //[nLepOther]
   Int_t           LepOther_pdgId[6];   //[nLepOther]
   Int_t           LepOther_charge[6];   //[nLepOther]
   Float_t         LepOther_dxy[6];   //[nLepOther]
   Float_t         LepOther_dz[6];   //[nLepOther]
   Float_t         LepOther_edxy[6];   //[nLepOther]
   Float_t         LepOther_edz[6];   //[nLepOther]
   Float_t         LepOther_ip3d[6];   //[nLepOther]
   Float_t         LepOther_sip3d[6];   //[nLepOther]
   Int_t           LepOther_tightId[6];   //[nLepOther]
   Int_t           LepOther_convVeto[6];   //[nLepOther]
   Int_t           LepOther_lostHits[6];   //[nLepOther]
   Int_t           LepOther_looseIdSusy[6];   //[nLepOther]
   Float_t         LepOther_relIso03[6];   //[nLepOther]
   Float_t         LepOther_relIso04[6];   //[nLepOther]
   Float_t         LepOther_chargedHadRelIso03[6];   //[nLepOther]
   Float_t         LepOther_chargedHadRelIso04[6];   //[nLepOther]
   Int_t           LepOther_convVetoFull[6];   //[nLepOther]
   Int_t           LepOther_eleCutId[6];   //[nLepOther]
   Int_t           LepOther_eleMVAId[6];   //[nLepOther]
   Int_t           LepOther_tightCharge[6];   //[nLepOther]
   Float_t         LepOther_mvaId[6];   //[nLepOther]
   Float_t         LepOther_mvaIdTrig[6];   //[nLepOther]
   Float_t         LepOther_nStations[6];   //[nLepOther]
   Float_t         LepOther_trkKink[6];   //[nLepOther]
   Float_t         LepOther_caloCompatibility[6];   //[nLepOther]
   Float_t         LepOther_globalTrackChi2[6];   //[nLepOther]
   Int_t           LepOther_trackerLayers[6];   //[nLepOther]
   Int_t           LepOther_pixelLayers[6];   //[nLepOther]
   Float_t         LepOther_mvaTTH[6];   //[nLepOther]
   Float_t         LepOther_jetPtRatio[6];   //[nLepOther]
   Float_t         LepOther_jetBTagCSV[6];   //[nLepOther]
   Float_t         LepOther_jetDR[6];   //[nLepOther]
   Int_t           LepOther_softMuID[6];   //[nLepOther]
   Int_t           nTauGood;
   Float_t         TauGood_pt[3];   //[nTauGood]
   Float_t         TauGood_eta[3];   //[nTauGood]
   Float_t         TauGood_phi[3];   //[nTauGood]
   Float_t         TauGood_mass[3];   //[nTauGood]
   Int_t           TauGood_pdgId[3];   //[nTauGood]
   Int_t           TauGood_charge[3];   //[nTauGood]
   Float_t         TauGood_dxy[3];   //[nTauGood]
   Float_t         TauGood_dz[3];   //[nTauGood]
   Int_t           TauGood_idMVA2[3];   //[nTauGood]
   Int_t           TauGood_idCI3hit[3];   //[nTauGood]
   Float_t         TauGood_isoMVA2[3];   //[nTauGood]
   Int_t           nJet;
   Float_t         Jet_pt[8];   //[nJet]
   Float_t         Jet_eta[8];   //[nJet]
   Float_t         Jet_phi[8];   //[nJet]
   Float_t         Jet_mass[8];   //[nJet]
   Float_t         Jet_btagCSV[8];   //[nJet]
   Float_t         Jet_rawPt[8];   //[nJet]
   Float_t         Jet_quarkGluonID[8];   //[nJet]
   Int_t           Jet_puId[8];   //[nJet]
   Float_t         Jet_area[8];   //[nJet]
   Int_t           Jet_id[8];   //[nJet]
   Float_t         Jet_CHEF[8];   //[nJet]
   Float_t         Jet_NHEF[8];   //[nJet]
   Float_t         Jet_PHEF[8];   //[nJet]
   Float_t         Jet_MUEF[8];   //[nJet]
   Float_t         Jet_ELEF[8];   //[nJet]
   Int_t           nLepGood;
   Float_t         LepGood_pt[5];   //[nLepGood]
   Float_t         LepGood_eta[5];   //[nLepGood]
   Float_t         LepGood_phi[5];   //[nLepGood]
   Float_t         LepGood_mass[5];   //[nLepGood]
   Int_t           LepGood_pdgId[5];   //[nLepGood]
   Int_t           LepGood_charge[5];   //[nLepGood]
   Float_t         LepGood_dxy[5];   //[nLepGood]
   Float_t         LepGood_dz[5];   //[nLepGood]
   Float_t         LepGood_edxy[5];   //[nLepGood]
   Float_t         LepGood_edz[5];   //[nLepGood]
   Float_t         LepGood_ip3d[5];   //[nLepGood]
   Float_t         LepGood_sip3d[5];   //[nLepGood]
   Int_t           LepGood_tightId[5];   //[nLepGood]
   Int_t           LepGood_convVeto[5];   //[nLepGood]
   Int_t           LepGood_lostHits[5];   //[nLepGood]
   Int_t           LepGood_looseIdSusy[5];   //[nLepGood]
   Float_t         LepGood_relIso03[5];   //[nLepGood]
   Float_t         LepGood_relIso04[5];   //[nLepGood]
   Float_t         LepGood_chargedHadRelIso03[5];   //[nLepGood]
   Float_t         LepGood_chargedHadRelIso04[5];   //[nLepGood]
   Int_t           LepGood_convVetoFull[5];   //[nLepGood]
   Int_t           LepGood_eleCutId[5];   //[nLepGood]
   Int_t           LepGood_eleMVAId[5];   //[nLepGood]
   Int_t           LepGood_tightCharge[5];   //[nLepGood]
   Float_t         LepGood_mvaId[5];   //[nLepGood]
   Float_t         LepGood_mvaIdTrig[5];   //[nLepGood]
   Float_t         LepGood_nStations[5];   //[nLepGood]
   Float_t         LepGood_trkKink[5];   //[nLepGood]
   Float_t         LepGood_caloCompatibility[5];   //[nLepGood]
   Float_t         LepGood_globalTrackChi2[5];   //[nLepGood]
   Int_t           LepGood_trackerLayers[5];   //[nLepGood]
   Int_t           LepGood_pixelLayers[5];   //[nLepGood]
   Float_t         LepGood_mvaTTH[5];   //[nLepGood]
   Float_t         LepGood_jetPtRatio[5];   //[nLepGood]
   Float_t         LepGood_jetBTagCSV[5];   //[nLepGood]
   Float_t         LepGood_jetDR[5];   //[nLepGood]
   Int_t           LepGood_softMuID[5];   //[nLepGood]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_isData;   //!
   TBranch        *b_HLT_SingleMu;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_nVert;   //!
   TBranch        *b_nJet25;   //!
   TBranch        *b_nBJetLoose25;   //!
   TBranch        *b_nBJetMedium25;   //!
   TBranch        *b_nBJetTight25;   //!
   TBranch        *b_nJet40;   //!
   TBranch        *b_nJet40a;   //!
   TBranch        *b_nBJetLoose40;   //!
   TBranch        *b_nBJetMedium40;   //!
   TBranch        *b_nBJetTight40;   //!
   TBranch        *b_nLepGood20;   //!
   TBranch        *b_nLepGood15;   //!
   TBranch        *b_nLepGood10;   //!
   TBranch        *b_checkEcalDead;   //!
   TBranch        *b_checkhcalLaser;   //!
   TBranch        *b_checktrackingFailure;   //!
   TBranch        *b_checkprimaryVertex;   //!
   TBranch        *b_checknoscraping;   //!
   TBranch        *b_checktrackIsolationMaker;   //!
   TBranch        *b_checkmetNoiseCleaning;   //!
   TBranch        *b_checkeeBadSc;   //!
   TBranch        *b_checkecalLaser;   //!
   TBranch        *b_checktotalKinematics;   //!
   TBranch        *b_checkCSCTightHalo;   //!
   TBranch        *b_checkHBHENoise;   //!
   TBranch        *b_htJet25;   //!
   TBranch        *b_mhtJet25;   //!
   TBranch        *b_htJet40j;   //!
   TBranch        *b_htJet40ja;   //!
   TBranch        *b_htJet40;   //!
   TBranch        *b_htJet40a;   //!
   TBranch        *b_mhtJet40;   //!
   TBranch        *b_mhtJet40a;   //!
   TBranch        *b_mZ1;   //!
   TBranch        *b_mZ1SFSS;   //!
   TBranch        *b_minMllSFOS;   //!
   TBranch        *b_maxMllSFOS;   //!
   TBranch        *b_minMllAFOS;   //!
   TBranch        *b_maxMllAFOS;   //!
   TBranch        *b_minMllAFSS;   //!
   TBranch        *b_maxMllAFSS;   //!
   TBranch        *b_minMllAFAS;   //!
   TBranch        *b_maxMllAFAS;   //!
   TBranch        *b_m2l;   //!
   TBranch        *b_minMWjj;   //!
   TBranch        *b_minMWjjPt;   //!
   TBranch        *b_bestMWjj;   //!
   TBranch        *b_bestMWjjPt;   //!
   TBranch        *b_bestMTopHad;   //!
   TBranch        *b_bestMTopHadPt;   //!
   TBranch        *b_met_pt;   //!
   TBranch        *b_met_eta;   //!
   TBranch        *b_met_phi;   //!
   TBranch        *b_met_mass;   //!
   TBranch        *b_met_sumEt;   //!
   TBranch        *b_metNoPU_pt;   //!
   TBranch        *b_metNoPU_eta;   //!
   TBranch        *b_metNoPU_phi;   //!
   TBranch        *b_metNoPU_mass;   //!
   TBranch        *b_nLepOther;   //!
   TBranch        *b_LepOther_pt;   //!
   TBranch        *b_LepOther_eta;   //!
   TBranch        *b_LepOther_phi;   //!
   TBranch        *b_LepOther_mass;   //!
   TBranch        *b_LepOther_pdgId;   //!
   TBranch        *b_LepOther_charge;   //!
   TBranch        *b_LepOther_dxy;   //!
   TBranch        *b_LepOther_dz;   //!
   TBranch        *b_LepOther_edxy;   //!
   TBranch        *b_LepOther_edz;   //!
   TBranch        *b_LepOther_ip3d;   //!
   TBranch        *b_LepOther_sip3d;   //!
   TBranch        *b_LepOther_tightId;   //!
   TBranch        *b_LepOther_convVeto;   //!
   TBranch        *b_LepOther_lostHits;   //!
   TBranch        *b_LepOther_looseIdSusy;   //!
   TBranch        *b_LepOther_relIso03;   //!
   TBranch        *b_LepOther_relIso04;   //!
   TBranch        *b_LepOther_chargedHadRelIso03;   //!
   TBranch        *b_LepOther_chargedHadRelIso04;   //!
   TBranch        *b_LepOther_convVetoFull;   //!
   TBranch        *b_LepOther_eleCutId;   //!
   TBranch        *b_LepOther_eleMVAId;   //!
   TBranch        *b_LepOther_tightCharge;   //!
   TBranch        *b_LepOther_mvaId;   //!
   TBranch        *b_LepOther_mvaIdTrig;   //!
   TBranch        *b_LepOther_nStations;   //!
   TBranch        *b_LepOther_trkKink;   //!
   TBranch        *b_LepOther_caloCompatibility;   //!
   TBranch        *b_LepOther_globalTrackChi2;   //!
   TBranch        *b_LepOther_trackerLayers;   //!
   TBranch        *b_LepOther_pixelLayers;   //!
   TBranch        *b_LepOther_mvaTTH;   //!
   TBranch        *b_LepOther_jetPtRatio;   //!
   TBranch        *b_LepOther_jetBTagCSV;   //!
   TBranch        *b_LepOther_jetDR;   //!
   TBranch        *b_LepOther_softMuID;   //!
   TBranch        *b_nTauGood;   //!
   TBranch        *b_TauGood_pt;   //!
   TBranch        *b_TauGood_eta;   //!
   TBranch        *b_TauGood_phi;   //!
   TBranch        *b_TauGood_mass;   //!
   TBranch        *b_TauGood_pdgId;   //!
   TBranch        *b_TauGood_charge;   //!
   TBranch        *b_TauGood_dxy;   //!
   TBranch        *b_TauGood_dz;   //!
   TBranch        *b_TauGood_idMVA2;   //!
   TBranch        *b_TauGood_idCI3hit;   //!
   TBranch        *b_TauGood_isoMVA2;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_mass;   //!
   TBranch        *b_Jet_btagCSV;   //!
   TBranch        *b_Jet_rawPt;   //!
   TBranch        *b_Jet_quarkGluonID;   //!
   TBranch        *b_Jet_puId;   //!
   TBranch        *b_Jet_area;   //!
   TBranch        *b_Jet_id;   //!
   TBranch        *b_Jet_CHEF;   //!
   TBranch        *b_Jet_NHEF;   //!
   TBranch        *b_Jet_PHEF;   //!
   TBranch        *b_Jet_MUEF;   //!
   TBranch        *b_Jet_ELEF;   //!
   TBranch        *b_nLepGood;   //!
   TBranch        *b_LepGood_pt;   //!
   TBranch        *b_LepGood_eta;   //!
   TBranch        *b_LepGood_phi;   //!
   TBranch        *b_LepGood_mass;   //!
   TBranch        *b_LepGood_pdgId;   //!
   TBranch        *b_LepGood_charge;   //!
   TBranch        *b_LepGood_dxy;   //!
   TBranch        *b_LepGood_dz;   //!
   TBranch        *b_LepGood_edxy;   //!
   TBranch        *b_LepGood_edz;   //!
   TBranch        *b_LepGood_ip3d;   //!
   TBranch        *b_LepGood_sip3d;   //!
   TBranch        *b_LepGood_tightId;   //!
   TBranch        *b_LepGood_convVeto;   //!
   TBranch        *b_LepGood_lostHits;   //!
   TBranch        *b_LepGood_looseIdSusy;   //!
   TBranch        *b_LepGood_relIso03;   //!
   TBranch        *b_LepGood_relIso04;   //!
   TBranch        *b_LepGood_chargedHadRelIso03;   //!
   TBranch        *b_LepGood_chargedHadRelIso04;   //!
   TBranch        *b_LepGood_convVetoFull;   //!
   TBranch        *b_LepGood_eleCutId;   //!
   TBranch        *b_LepGood_eleMVAId;   //!
   TBranch        *b_LepGood_tightCharge;   //!
   TBranch        *b_LepGood_mvaId;   //!
   TBranch        *b_LepGood_mvaIdTrig;   //!
   TBranch        *b_LepGood_nStations;   //!
   TBranch        *b_LepGood_trkKink;   //!
   TBranch        *b_LepGood_caloCompatibility;   //!
   TBranch        *b_LepGood_globalTrackChi2;   //!
   TBranch        *b_LepGood_trackerLayers;   //!
   TBranch        *b_LepGood_pixelLayers;   //!
   TBranch        *b_LepGood_mvaTTH;   //!
   TBranch        *b_LepGood_jetPtRatio;   //!
   TBranch        *b_LepGood_jetBTagCSV;   //!
   TBranch        *b_LepGood_jetDR;   //!
   TBranch        *b_LepGood_softMuID;   //!

   A7LoopBase(SampleInfo mySample);
   void setTree(SampleInfo mySample);
   A7LoopBase();
   virtual ~A7LoopBase();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    CutLepGenAcc(Long64_t entry);
   virtual Int_t    CutRECO();

   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual TFile*   Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
  
};
 A7LoopBase::A7LoopBase(SampleInfo mySample)
{
  setTree(mySample );
};
A7LoopBase::A7LoopBase( )
{
  cout << "Empty contructor: should not be used"<<endl;
};


#endif
