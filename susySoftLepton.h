//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jul  9 13:57:58 2014 by ROOT version 5.32/00
// from TTree treeProducerSusySoftlepton/treeProducerSusySoftlepton
// found on file: /afs/cern.ch/work/m/mstoye/CMGTools_SUSY_git/CMSSW_5_3_14/src/CMGTools/TTHAnalysis/cfg/ttFully/TTJetsLep/treeProducerSusySoftlepton/treeProducerSusySoftlepton_tree.root
//////////////////////////////////////////////////////////

#ifndef susySoftLepton_h
#define susySoftLepton_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "SampleInfo.h" 
#include <TString.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class susySoftLepton {

 private:
  TString outputFileName;
  float weight;

public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           run;
   Int_t           lumi;
   Int_t           evt;
   Int_t           isData;
   Int_t           HLT_MetTrigger;
   Int_t           HLT_MET15;
   Int_t           HLT_HT650;
   Int_t           HLT_triggers_HTMET;
   Int_t           HLT_MetTriggerParked;
   Int_t           HLT_SingleMu;
   Int_t           HLT_MuEG;
   Int_t           HLT_TripleEl;
   Int_t           HLT_DoubleEl;
   Int_t           HLT_DoubleMu;
   Float_t         puWeight;
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
   Int_t           GenHeaviestQCDFlavour;
   Float_t         LepEff_1lep;
   Float_t         LepEff_2lep;
   Int_t           GenSusyMScan1;
   Int_t           GenSusyMScan2;
   Int_t           GenSusyMScan3;
   Int_t           GenSusyMScan4;
   Int_t           GenSusyMGluino;
   Int_t           GenSusyMGravitino;
   Int_t           GenSusyMStop;
   Int_t           GenSusyMSbottom;
   Int_t           GenSusyMStop2;
   Int_t           GenSusyMSbottom2;
   Int_t           GenSusyMSquark;
   Int_t           GenSusyMNeutralino;
   Int_t           GenSusyMNeutralino2;
   Int_t           GenSusyMNeutralino3;
   Int_t           GenSusyMNeutralino4;
   Int_t           GenSusyMChargino;
   Int_t           GenSusyMChargino2;
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
   Float_t         met_genPt;
   Float_t         met_genPhi;
   Float_t         met_genEta;
   Float_t         metNoPU_pt;
   Float_t         metNoPU_eta;
   Float_t         metNoPU_phi;
   Float_t         metNoPU_mass;
   Int_t           nGenLep;
   Float_t         GenLep_pt[2];   //[nGenLep]
   Float_t         GenLep_eta[2];   //[nGenLep]
   Float_t         GenLep_phi[2];   //[nGenLep]
   Float_t         GenLep_mass[2];   //[nGenLep]
   Int_t           GenLep_pdgId[2];   //[nGenLep]
   Float_t         GenLep_charge[2];   //[nGenLep]
   Int_t           GenLep_sourceId[2];   //[nGenLep]
   Int_t           nGenQuark;
   Float_t         GenQuark_pt[1];   //[nGenQuark]
   Float_t         GenQuark_eta[1];   //[nGenQuark]
   Float_t         GenQuark_phi[1];   //[nGenQuark]
   Float_t         GenQuark_mass[1];   //[nGenQuark]
   Int_t           GenQuark_pdgId[1];   //[nGenQuark]
   Float_t         GenQuark_charge[1];   //[nGenQuark]
   Int_t           GenQuark_sourceId[1];   //[nGenQuark]
   Int_t           nGenLepFromTau;
   Float_t         GenLepFromTau_pt[2];   //[nGenLepFromTau]
   Float_t         GenLepFromTau_eta[2];   //[nGenLepFromTau]
   Float_t         GenLepFromTau_phi[2];   //[nGenLepFromTau]
   Float_t         GenLepFromTau_mass[2];   //[nGenLepFromTau]
   Int_t           GenLepFromTau_pdgId[2];   //[nGenLepFromTau]
   Float_t         GenLepFromTau_charge[2];   //[nGenLepFromTau]
   Int_t           GenLepFromTau_sourceId[2];   //[nGenLepFromTau]
   Int_t           nLepOther;
   Float_t         LepOther_pt[8];   //[nLepOther]
   Float_t         LepOther_eta[8];   //[nLepOther]
   Float_t         LepOther_phi[8];   //[nLepOther]
   Float_t         LepOther_mass[8];   //[nLepOther]
   Int_t           LepOther_pdgId[8];   //[nLepOther]
   Int_t           LepOther_charge[8];   //[nLepOther]
   Float_t         LepOther_dxy[8];   //[nLepOther]
   Float_t         LepOther_dz[8];   //[nLepOther]
   Float_t         LepOther_edxy[8];   //[nLepOther]
   Float_t         LepOther_edz[8];   //[nLepOther]
   Float_t         LepOther_ip3d[8];   //[nLepOther]
   Float_t         LepOther_sip3d[8];   //[nLepOther]
   Int_t           LepOther_tightId[8];   //[nLepOther]
   Int_t           LepOther_convVeto[8];   //[nLepOther]
   Int_t           LepOther_lostHits[8];   //[nLepOther]
   Int_t           LepOther_looseIdSusy[8];   //[nLepOther]
   Float_t         LepOther_relIso03[8];   //[nLepOther]
   Float_t         LepOther_relIso04[8];   //[nLepOther]
   Float_t         LepOther_chargedHadRelIso03[8];   //[nLepOther]
   Float_t         LepOther_chargedHadRelIso04[8];   //[nLepOther]
   Int_t           LepOther_convVetoFull[8];   //[nLepOther]
   Int_t           LepOther_eleCutId[8];   //[nLepOther]
   Int_t           LepOther_eleMVAId[8];   //[nLepOther]
   Int_t           LepOther_tightCharge[8];   //[nLepOther]
   Float_t         LepOther_mvaId[8];   //[nLepOther]
   Float_t         LepOther_mvaIdTrig[8];   //[nLepOther]
   Float_t         LepOther_nStations[8];   //[nLepOther]
   Float_t         LepOther_trkKink[8];   //[nLepOther]
   Int_t           LepOther_trackerLayers[8];   //[nLepOther]
   Int_t           LepOther_pixelLayers[8];   //[nLepOther]
   Float_t         LepOther_mvaTTH[8];   //[nLepOther]
   Float_t         LepOther_jetPtRatio[8];   //[nLepOther]
   Float_t         LepOther_jetBTagCSV[8];   //[nLepOther]
   Float_t         LepOther_jetDR[8];   //[nLepOther]
   Int_t           LepOther_mcMatchId[8];   //[nLepOther]
   Int_t           LepOther_mcMatchAny[8];   //[nLepOther]
   Int_t           LepOther_mcMatchTau[8];   //[nLepOther]
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
   Int_t           TauGood_mcMatchId[3];   //[nTauGood]
   Int_t           nGenBQuark;
   Float_t         GenBQuark_pt[2];   //[nGenBQuark]
   Float_t         GenBQuark_eta[2];   //[nGenBQuark]
   Float_t         GenBQuark_phi[2];   //[nGenBQuark]
   Float_t         GenBQuark_mass[2];   //[nGenBQuark]
   Int_t           GenBQuark_pdgId[2];   //[nGenBQuark]
   Float_t         GenBQuark_charge[2];   //[nGenBQuark]
   Int_t           nJet;
   Float_t         Jet_pt[8];   //[nJet]
   Float_t         Jet_eta[8];   //[nJet]
   Float_t         Jet_phi[8];   //[nJet]
   Float_t         Jet_mass[8];   //[nJet]
   Float_t         Jet_btagCSV[8];   //[nJet]
   Float_t         Jet_rawPt[8];   //[nJet]
   Float_t         Jet_mcPt[8];   //[nJet]
   Int_t           Jet_mcFlavour[8];   //[nJet]
   Float_t         Jet_quarkGluonID[8];   //[nJet]
   Int_t           Jet_mcMatchId[8];   //[nJet]
   Int_t           Jet_mcMatchFlav[8];   //[nJet]
   Int_t           Jet_PuId_full[8];   //[nJet]
   Int_t           Jet_PuId_simple[8];   //[nJet]
   Int_t           Jet_PuId_cut_based[8];   //[nJet]
   Int_t           Jet_Id[8];   //[nJet]
   Int_t           nGenTop;
   Float_t         GenTop_pt[2];   //[nGenTop]
   Float_t         GenTop_eta[2];   //[nGenTop]
   Float_t         GenTop_phi[2];   //[nGenTop]
   Float_t         GenTop_mass[2];   //[nGenTop]
   Int_t           GenTop_pdgId[2];   //[nGenTop]
   Float_t         GenTop_charge[2];   //[nGenTop]
   Int_t           nGenP6StatusThree;
   Float_t         GenP6StatusThree_pt[19];   //[nGenP6StatusThree]
   Float_t         GenP6StatusThree_eta[19];   //[nGenP6StatusThree]
   Float_t         GenP6StatusThree_phi[19];   //[nGenP6StatusThree]
   Float_t         GenP6StatusThree_mass[19];   //[nGenP6StatusThree]
   Int_t           GenP6StatusThree_pdgId[19];   //[nGenP6StatusThree]
   Float_t         GenP6StatusThree_charge[19];   //[nGenP6StatusThree]
   Int_t           GenP6StatusThree_motherId[19];   //[nGenP6StatusThree]
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
   Int_t           LepGood_trackerLayers[5];   //[nLepGood]
   Int_t           LepGood_pixelLayers[5];   //[nLepGood]
   Float_t         LepGood_mvaTTH[5];   //[nLepGood]
   Float_t         LepGood_jetPtRatio[5];   //[nLepGood]
   Float_t         LepGood_jetBTagCSV[5];   //[nLepGood]
   Float_t         LepGood_jetDR[5];   //[nLepGood]
   Int_t           LepGood_mcMatchId[5];   //[nLepGood]
   Int_t           LepGood_mcMatchAny[5];   //[nLepGood]
   Int_t           LepGood_mcMatchTau[5];   //[nLepGood]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_isData;   //!
   TBranch        *b_HLT_MetTrigger;   //!
   TBranch        *b_HLT_MET15;   //!
   TBranch        *b_HLT_HT650;   //!
   TBranch        *b_HLT_triggers_HTMET;   //!
   TBranch        *b_HLT_MetTriggerParked;   //!
   TBranch        *b_HLT_SingleMu;   //!
   TBranch        *b_HLT_MuEG;   //!
   TBranch        *b_HLT_TripleEl;   //!
   TBranch        *b_HLT_DoubleEl;   //!
   TBranch        *b_HLT_DoubleMu;   //!
   TBranch        *b_puWeight;   //!
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
   TBranch        *b_GenHeaviestQCDFlavour;   //!
   TBranch        *b_LepEff_1lep;   //!
   TBranch        *b_LepEff_2lep;   //!
   TBranch        *b_GenSusyMScan1;   //!
   TBranch        *b_GenSusyMScan2;   //!
   TBranch        *b_GenSusyMScan3;   //!
   TBranch        *b_GenSusyMScan4;   //!
   TBranch        *b_GenSusyMGluino;   //!
   TBranch        *b_GenSusyMGravitino;   //!
   TBranch        *b_GenSusyMStop;   //!
   TBranch        *b_GenSusyMSbottom;   //!
   TBranch        *b_GenSusyMStop2;   //!
   TBranch        *b_GenSusyMSbottom2;   //!
   TBranch        *b_GenSusyMSquark;   //!
   TBranch        *b_GenSusyMNeutralino;   //!
   TBranch        *b_GenSusyMNeutralino2;   //!
   TBranch        *b_GenSusyMNeutralino3;   //!
   TBranch        *b_GenSusyMNeutralino4;   //!
   TBranch        *b_GenSusyMChargino;   //!
   TBranch        *b_GenSusyMChargino2;   //!
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
   TBranch        *b_met_genPt;   //!
   TBranch        *b_met_genPhi;   //!
   TBranch        *b_met_genEta;   //!
   TBranch        *b_metNoPU_pt;   //!
   TBranch        *b_metNoPU_eta;   //!
   TBranch        *b_metNoPU_phi;   //!
   TBranch        *b_metNoPU_mass;   //!
   TBranch        *b_nGenLep;   //!
   TBranch        *b_GenLep_pt;   //!
   TBranch        *b_GenLep_eta;   //!
   TBranch        *b_GenLep_phi;   //!
   TBranch        *b_GenLep_mass;   //!
   TBranch        *b_GenLep_pdgId;   //!
   TBranch        *b_GenLep_charge;   //!
   TBranch        *b_GenLep_sourceId;   //!
   TBranch        *b_nGenQuark;   //!
   TBranch        *b_GenQuark_pt;   //!
   TBranch        *b_GenQuark_eta;   //!
   TBranch        *b_GenQuark_phi;   //!
   TBranch        *b_GenQuark_mass;   //!
   TBranch        *b_GenQuark_pdgId;   //!
   TBranch        *b_GenQuark_charge;   //!
   TBranch        *b_GenQuark_sourceId;   //!
   TBranch        *b_nGenLepFromTau;   //!
   TBranch        *b_GenLepFromTau_pt;   //!
   TBranch        *b_GenLepFromTau_eta;   //!
   TBranch        *b_GenLepFromTau_phi;   //!
   TBranch        *b_GenLepFromTau_mass;   //!
   TBranch        *b_GenLepFromTau_pdgId;   //!
   TBranch        *b_GenLepFromTau_charge;   //!
   TBranch        *b_GenLepFromTau_sourceId;   //!
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
   TBranch        *b_LepOther_trackerLayers;   //!
   TBranch        *b_LepOther_pixelLayers;   //!
   TBranch        *b_LepOther_mvaTTH;   //!
   TBranch        *b_LepOther_jetPtRatio;   //!
   TBranch        *b_LepOther_jetBTagCSV;   //!
   TBranch        *b_LepOther_jetDR;   //!
   TBranch        *b_LepOther_mcMatchId;   //!
   TBranch        *b_LepOther_mcMatchAny;   //!
   TBranch        *b_LepOther_mcMatchTau;   //!
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
   TBranch        *b_TauGood_mcMatchId;   //!
   TBranch        *b_nGenBQuark;   //!
   TBranch        *b_GenBQuark_pt;   //!
   TBranch        *b_GenBQuark_eta;   //!
   TBranch        *b_GenBQuark_phi;   //!
   TBranch        *b_GenBQuark_mass;   //!
   TBranch        *b_GenBQuark_pdgId;   //!
   TBranch        *b_GenBQuark_charge;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_mass;   //!
   TBranch        *b_Jet_btagCSV;   //!
   TBranch        *b_Jet_rawPt;   //!
   TBranch        *b_Jet_mcPt;   //!
   TBranch        *b_Jet_mcFlavour;   //!
   TBranch        *b_Jet_quarkGluonID;   //!
   TBranch        *b_Jet_mcMatchId;   //!
   TBranch        *b_Jet_mcMatchFlav;   //!
   TBranch        *b_Jet_PuId_full;   //!
   TBranch        *b_Jet_PuId_simple;   //!
   TBranch        *b_Jet_PuId_cut_based;   //!
   TBranch        *b_Jet_Id;   //!
   TBranch        *b_nGenTop;   //!
   TBranch        *b_GenTop_pt;   //!
   TBranch        *b_GenTop_eta;   //!
   TBranch        *b_GenTop_phi;   //!
   TBranch        *b_GenTop_mass;   //!
   TBranch        *b_GenTop_pdgId;   //!
   TBranch        *b_GenTop_charge;   //!
   TBranch        *b_nGenP6StatusThree;   //!
   TBranch        *b_GenP6StatusThree_pt;   //!
   TBranch        *b_GenP6StatusThree_eta;   //!
   TBranch        *b_GenP6StatusThree_phi;   //!
   TBranch        *b_GenP6StatusThree_mass;   //!
   TBranch        *b_GenP6StatusThree_pdgId;   //!
   TBranch        *b_GenP6StatusThree_charge;   //!
   TBranch        *b_GenP6StatusThree_motherId;   //!
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
   TBranch        *b_LepGood_trackerLayers;   //!
   TBranch        *b_LepGood_pixelLayers;   //!
   TBranch        *b_LepGood_mvaTTH;   //!
   TBranch        *b_LepGood_jetPtRatio;   //!
   TBranch        *b_LepGood_jetBTagCSV;   //!
   TBranch        *b_LepGood_jetDR;   //!
   TBranch        *b_LepGood_mcMatchId;   //!
   TBranch        *b_LepGood_mcMatchAny;   //!
   TBranch        *b_LepGood_mcMatchTau;   //!

   susySoftLepton(SampleInfo mySample);
   virtual ~susySoftLepton();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef susySoftLepton_cxx
susySoftLepton::susySoftLepton(SampleInfo mySample)
{
  TFile *f = TFile::Open(mySample.FilePath);
  TTree* tree = (TTree*) f->Get("markusTreeProducer"); 
  outputFileName = mySample.OutputFileNameTag+".root";
  weight = mySample.weight();
  Init(tree);
}

susySoftLepton::~susySoftLepton()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t susySoftLepton::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t susySoftLepton::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void susySoftLepton::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("isData", &isData, &b_isData);
   fChain->SetBranchAddress("HLT_MetTrigger", &HLT_MetTrigger, &b_HLT_MetTrigger);
   fChain->SetBranchAddress("HLT_MET15", &HLT_MET15, &b_HLT_MET15);
   fChain->SetBranchAddress("HLT_HT650", &HLT_HT650, &b_HLT_HT650);
   fChain->SetBranchAddress("HLT_triggers_HTMET", &HLT_triggers_HTMET, &b_HLT_triggers_HTMET);
   fChain->SetBranchAddress("HLT_MetTriggerParked", &HLT_MetTriggerParked, &b_HLT_MetTriggerParked);
   fChain->SetBranchAddress("HLT_SingleMu", &HLT_SingleMu, &b_HLT_SingleMu);
   fChain->SetBranchAddress("HLT_MuEG", &HLT_MuEG, &b_HLT_MuEG);
   fChain->SetBranchAddress("HLT_TripleEl", &HLT_TripleEl, &b_HLT_TripleEl);
   fChain->SetBranchAddress("HLT_DoubleEl", &HLT_DoubleEl, &b_HLT_DoubleEl);
   fChain->SetBranchAddress("HLT_DoubleMu", &HLT_DoubleMu, &b_HLT_DoubleMu);
   fChain->SetBranchAddress("puWeight", &puWeight, &b_puWeight);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("nVert", &nVert, &b_nVert);
   fChain->SetBranchAddress("nJet25", &nJet25, &b_nJet25);
   fChain->SetBranchAddress("nBJetLoose25", &nBJetLoose25, &b_nBJetLoose25);
   fChain->SetBranchAddress("nBJetMedium25", &nBJetMedium25, &b_nBJetMedium25);
   fChain->SetBranchAddress("nBJetTight25", &nBJetTight25, &b_nBJetTight25);
   fChain->SetBranchAddress("nJet40", &nJet40, &b_nJet40);
   fChain->SetBranchAddress("nJet40a", &nJet40a, &b_nJet40a);
   fChain->SetBranchAddress("nBJetLoose40", &nBJetLoose40, &b_nBJetLoose40);
   fChain->SetBranchAddress("nBJetMedium40", &nBJetMedium40, &b_nBJetMedium40);
   fChain->SetBranchAddress("nBJetTight40", &nBJetTight40, &b_nBJetTight40);
   fChain->SetBranchAddress("nLepGood20", &nLepGood20, &b_nLepGood20);
   fChain->SetBranchAddress("nLepGood15", &nLepGood15, &b_nLepGood15);
   fChain->SetBranchAddress("nLepGood10", &nLepGood10, &b_nLepGood10);
   fChain->SetBranchAddress("GenHeaviestQCDFlavour", &GenHeaviestQCDFlavour, &b_GenHeaviestQCDFlavour);
   fChain->SetBranchAddress("LepEff_1lep", &LepEff_1lep, &b_LepEff_1lep);
   fChain->SetBranchAddress("LepEff_2lep", &LepEff_2lep, &b_LepEff_2lep);
   fChain->SetBranchAddress("GenSusyMScan1", &GenSusyMScan1, &b_GenSusyMScan1);
   fChain->SetBranchAddress("GenSusyMScan2", &GenSusyMScan2, &b_GenSusyMScan2);
   fChain->SetBranchAddress("GenSusyMScan3", &GenSusyMScan3, &b_GenSusyMScan3);
   fChain->SetBranchAddress("GenSusyMScan4", &GenSusyMScan4, &b_GenSusyMScan4);
   fChain->SetBranchAddress("GenSusyMGluino", &GenSusyMGluino, &b_GenSusyMGluino);
   fChain->SetBranchAddress("GenSusyMGravitino", &GenSusyMGravitino, &b_GenSusyMGravitino);
   fChain->SetBranchAddress("GenSusyMStop", &GenSusyMStop, &b_GenSusyMStop);
   fChain->SetBranchAddress("GenSusyMSbottom", &GenSusyMSbottom, &b_GenSusyMSbottom);
   fChain->SetBranchAddress("GenSusyMStop2", &GenSusyMStop2, &b_GenSusyMStop2);
   fChain->SetBranchAddress("GenSusyMSbottom2", &GenSusyMSbottom2, &b_GenSusyMSbottom2);
   fChain->SetBranchAddress("GenSusyMSquark", &GenSusyMSquark, &b_GenSusyMSquark);
   fChain->SetBranchAddress("GenSusyMNeutralino", &GenSusyMNeutralino, &b_GenSusyMNeutralino);
   fChain->SetBranchAddress("GenSusyMNeutralino2", &GenSusyMNeutralino2, &b_GenSusyMNeutralino2);
   fChain->SetBranchAddress("GenSusyMNeutralino3", &GenSusyMNeutralino3, &b_GenSusyMNeutralino3);
   fChain->SetBranchAddress("GenSusyMNeutralino4", &GenSusyMNeutralino4, &b_GenSusyMNeutralino4);
   fChain->SetBranchAddress("GenSusyMChargino", &GenSusyMChargino, &b_GenSusyMChargino);
   fChain->SetBranchAddress("GenSusyMChargino2", &GenSusyMChargino2, &b_GenSusyMChargino2);
   fChain->SetBranchAddress("htJet25", &htJet25, &b_htJet25);
   fChain->SetBranchAddress("mhtJet25", &mhtJet25, &b_mhtJet25);
   fChain->SetBranchAddress("htJet40j", &htJet40j, &b_htJet40j);
   fChain->SetBranchAddress("htJet40ja", &htJet40ja, &b_htJet40ja);
   fChain->SetBranchAddress("htJet40", &htJet40, &b_htJet40);
   fChain->SetBranchAddress("htJet40a", &htJet40a, &b_htJet40a);
   fChain->SetBranchAddress("mhtJet40", &mhtJet40, &b_mhtJet40);
   fChain->SetBranchAddress("mhtJet40a", &mhtJet40a, &b_mhtJet40a);
   fChain->SetBranchAddress("mZ1", &mZ1, &b_mZ1);
   fChain->SetBranchAddress("mZ1SFSS", &mZ1SFSS, &b_mZ1SFSS);
   fChain->SetBranchAddress("minMllSFOS", &minMllSFOS, &b_minMllSFOS);
   fChain->SetBranchAddress("maxMllSFOS", &maxMllSFOS, &b_maxMllSFOS);
   fChain->SetBranchAddress("minMllAFOS", &minMllAFOS, &b_minMllAFOS);
   fChain->SetBranchAddress("maxMllAFOS", &maxMllAFOS, &b_maxMllAFOS);
   fChain->SetBranchAddress("minMllAFSS", &minMllAFSS, &b_minMllAFSS);
   fChain->SetBranchAddress("maxMllAFSS", &maxMllAFSS, &b_maxMllAFSS);
   fChain->SetBranchAddress("minMllAFAS", &minMllAFAS, &b_minMllAFAS);
   fChain->SetBranchAddress("maxMllAFAS", &maxMllAFAS, &b_maxMllAFAS);
   fChain->SetBranchAddress("m2l", &m2l, &b_m2l);
   fChain->SetBranchAddress("minMWjj", &minMWjj, &b_minMWjj);
   fChain->SetBranchAddress("minMWjjPt", &minMWjjPt, &b_minMWjjPt);
   fChain->SetBranchAddress("bestMWjj", &bestMWjj, &b_bestMWjj);
   fChain->SetBranchAddress("bestMWjjPt", &bestMWjjPt, &b_bestMWjjPt);
   fChain->SetBranchAddress("bestMTopHad", &bestMTopHad, &b_bestMTopHad);
   fChain->SetBranchAddress("bestMTopHadPt", &bestMTopHadPt, &b_bestMTopHadPt);
   fChain->SetBranchAddress("met_pt", &met_pt, &b_met_pt);
   fChain->SetBranchAddress("met_eta", &met_eta, &b_met_eta);
   fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
   fChain->SetBranchAddress("met_mass", &met_mass, &b_met_mass);
   fChain->SetBranchAddress("met_sumEt", &met_sumEt, &b_met_sumEt);
   fChain->SetBranchAddress("met_genPt", &met_genPt, &b_met_genPt);
   fChain->SetBranchAddress("met_genPhi", &met_genPhi, &b_met_genPhi);
   fChain->SetBranchAddress("met_genEta", &met_genEta, &b_met_genEta);
   fChain->SetBranchAddress("metNoPU_pt", &metNoPU_pt, &b_metNoPU_pt);
   fChain->SetBranchAddress("metNoPU_eta", &metNoPU_eta, &b_metNoPU_eta);
   fChain->SetBranchAddress("metNoPU_phi", &metNoPU_phi, &b_metNoPU_phi);
   fChain->SetBranchAddress("metNoPU_mass", &metNoPU_mass, &b_metNoPU_mass);
   fChain->SetBranchAddress("nGenLep", &nGenLep, &b_nGenLep);
   fChain->SetBranchAddress("GenLep_pt", GenLep_pt, &b_GenLep_pt);
   fChain->SetBranchAddress("GenLep_eta", GenLep_eta, &b_GenLep_eta);
   fChain->SetBranchAddress("GenLep_phi", GenLep_phi, &b_GenLep_phi);
   fChain->SetBranchAddress("GenLep_mass", GenLep_mass, &b_GenLep_mass);
   fChain->SetBranchAddress("GenLep_pdgId", GenLep_pdgId, &b_GenLep_pdgId);
   fChain->SetBranchAddress("GenLep_charge", GenLep_charge, &b_GenLep_charge);
   fChain->SetBranchAddress("GenLep_sourceId", GenLep_sourceId, &b_GenLep_sourceId);
   fChain->SetBranchAddress("nGenQuark", &nGenQuark, &b_nGenQuark);
   fChain->SetBranchAddress("GenQuark_pt", &GenQuark_pt, &b_GenQuark_pt);
   fChain->SetBranchAddress("GenQuark_eta", &GenQuark_eta, &b_GenQuark_eta);
   fChain->SetBranchAddress("GenQuark_phi", &GenQuark_phi, &b_GenQuark_phi);
   fChain->SetBranchAddress("GenQuark_mass", &GenQuark_mass, &b_GenQuark_mass);
   fChain->SetBranchAddress("GenQuark_pdgId", &GenQuark_pdgId, &b_GenQuark_pdgId);
   fChain->SetBranchAddress("GenQuark_charge", &GenQuark_charge, &b_GenQuark_charge);
   fChain->SetBranchAddress("GenQuark_sourceId", &GenQuark_sourceId, &b_GenQuark_sourceId);
   fChain->SetBranchAddress("nGenLepFromTau", &nGenLepFromTau, &b_nGenLepFromTau);
   fChain->SetBranchAddress("GenLepFromTau_pt", GenLepFromTau_pt, &b_GenLepFromTau_pt);
   fChain->SetBranchAddress("GenLepFromTau_eta", GenLepFromTau_eta, &b_GenLepFromTau_eta);
   fChain->SetBranchAddress("GenLepFromTau_phi", GenLepFromTau_phi, &b_GenLepFromTau_phi);
   fChain->SetBranchAddress("GenLepFromTau_mass", GenLepFromTau_mass, &b_GenLepFromTau_mass);
   fChain->SetBranchAddress("GenLepFromTau_pdgId", GenLepFromTau_pdgId, &b_GenLepFromTau_pdgId);
   fChain->SetBranchAddress("GenLepFromTau_charge", GenLepFromTau_charge, &b_GenLepFromTau_charge);
   fChain->SetBranchAddress("GenLepFromTau_sourceId", GenLepFromTau_sourceId, &b_GenLepFromTau_sourceId);
   fChain->SetBranchAddress("nLepOther", &nLepOther, &b_nLepOther);
   fChain->SetBranchAddress("LepOther_pt", LepOther_pt, &b_LepOther_pt);
   fChain->SetBranchAddress("LepOther_eta", LepOther_eta, &b_LepOther_eta);
   fChain->SetBranchAddress("LepOther_phi", LepOther_phi, &b_LepOther_phi);
   fChain->SetBranchAddress("LepOther_mass", LepOther_mass, &b_LepOther_mass);
   fChain->SetBranchAddress("LepOther_pdgId", LepOther_pdgId, &b_LepOther_pdgId);
   fChain->SetBranchAddress("LepOther_charge", LepOther_charge, &b_LepOther_charge);
   fChain->SetBranchAddress("LepOther_dxy", LepOther_dxy, &b_LepOther_dxy);
   fChain->SetBranchAddress("LepOther_dz", LepOther_dz, &b_LepOther_dz);
   fChain->SetBranchAddress("LepOther_edxy", LepOther_edxy, &b_LepOther_edxy);
   fChain->SetBranchAddress("LepOther_edz", LepOther_edz, &b_LepOther_edz);
   fChain->SetBranchAddress("LepOther_ip3d", LepOther_ip3d, &b_LepOther_ip3d);
   fChain->SetBranchAddress("LepOther_sip3d", LepOther_sip3d, &b_LepOther_sip3d);
   fChain->SetBranchAddress("LepOther_tightId", LepOther_tightId, &b_LepOther_tightId);
   fChain->SetBranchAddress("LepOther_convVeto", LepOther_convVeto, &b_LepOther_convVeto);
   fChain->SetBranchAddress("LepOther_lostHits", LepOther_lostHits, &b_LepOther_lostHits);
   fChain->SetBranchAddress("LepOther_looseIdSusy", LepOther_looseIdSusy, &b_LepOther_looseIdSusy);
   fChain->SetBranchAddress("LepOther_relIso03", LepOther_relIso03, &b_LepOther_relIso03);
   fChain->SetBranchAddress("LepOther_relIso04", LepOther_relIso04, &b_LepOther_relIso04);
   fChain->SetBranchAddress("LepOther_chargedHadRelIso03", LepOther_chargedHadRelIso03, &b_LepOther_chargedHadRelIso03);
   fChain->SetBranchAddress("LepOther_chargedHadRelIso04", LepOther_chargedHadRelIso04, &b_LepOther_chargedHadRelIso04);
   fChain->SetBranchAddress("LepOther_convVetoFull", LepOther_convVetoFull, &b_LepOther_convVetoFull);
   fChain->SetBranchAddress("LepOther_eleCutId", LepOther_eleCutId, &b_LepOther_eleCutId);
   fChain->SetBranchAddress("LepOther_eleMVAId", LepOther_eleMVAId, &b_LepOther_eleMVAId);
   fChain->SetBranchAddress("LepOther_tightCharge", LepOther_tightCharge, &b_LepOther_tightCharge);
   fChain->SetBranchAddress("LepOther_mvaId", LepOther_mvaId, &b_LepOther_mvaId);
   fChain->SetBranchAddress("LepOther_mvaIdTrig", LepOther_mvaIdTrig, &b_LepOther_mvaIdTrig);
   fChain->SetBranchAddress("LepOther_nStations", LepOther_nStations, &b_LepOther_nStations);
   fChain->SetBranchAddress("LepOther_trkKink", LepOther_trkKink, &b_LepOther_trkKink);
   fChain->SetBranchAddress("LepOther_trackerLayers", LepOther_trackerLayers, &b_LepOther_trackerLayers);
   fChain->SetBranchAddress("LepOther_pixelLayers", LepOther_pixelLayers, &b_LepOther_pixelLayers);
   fChain->SetBranchAddress("LepOther_mvaTTH", LepOther_mvaTTH, &b_LepOther_mvaTTH);
   fChain->SetBranchAddress("LepOther_jetPtRatio", LepOther_jetPtRatio, &b_LepOther_jetPtRatio);
   fChain->SetBranchAddress("LepOther_jetBTagCSV", LepOther_jetBTagCSV, &b_LepOther_jetBTagCSV);
   fChain->SetBranchAddress("LepOther_jetDR", LepOther_jetDR, &b_LepOther_jetDR);
   fChain->SetBranchAddress("LepOther_mcMatchId", LepOther_mcMatchId, &b_LepOther_mcMatchId);
   fChain->SetBranchAddress("LepOther_mcMatchAny", LepOther_mcMatchAny, &b_LepOther_mcMatchAny);
   fChain->SetBranchAddress("LepOther_mcMatchTau", LepOther_mcMatchTau, &b_LepOther_mcMatchTau);
   fChain->SetBranchAddress("nTauGood", &nTauGood, &b_nTauGood);
   fChain->SetBranchAddress("TauGood_pt", TauGood_pt, &b_TauGood_pt);
   fChain->SetBranchAddress("TauGood_eta", TauGood_eta, &b_TauGood_eta);
   fChain->SetBranchAddress("TauGood_phi", TauGood_phi, &b_TauGood_phi);
   fChain->SetBranchAddress("TauGood_mass", TauGood_mass, &b_TauGood_mass);
   fChain->SetBranchAddress("TauGood_pdgId", TauGood_pdgId, &b_TauGood_pdgId);
   fChain->SetBranchAddress("TauGood_charge", TauGood_charge, &b_TauGood_charge);
   fChain->SetBranchAddress("TauGood_dxy", TauGood_dxy, &b_TauGood_dxy);
   fChain->SetBranchAddress("TauGood_dz", TauGood_dz, &b_TauGood_dz);
   fChain->SetBranchAddress("TauGood_idMVA2", TauGood_idMVA2, &b_TauGood_idMVA2);
   fChain->SetBranchAddress("TauGood_idCI3hit", TauGood_idCI3hit, &b_TauGood_idCI3hit);
   fChain->SetBranchAddress("TauGood_isoMVA2", TauGood_isoMVA2, &b_TauGood_isoMVA2);
   fChain->SetBranchAddress("TauGood_mcMatchId", TauGood_mcMatchId, &b_TauGood_mcMatchId);
   fChain->SetBranchAddress("nGenBQuark", &nGenBQuark, &b_nGenBQuark);
   fChain->SetBranchAddress("GenBQuark_pt", GenBQuark_pt, &b_GenBQuark_pt);
   fChain->SetBranchAddress("GenBQuark_eta", GenBQuark_eta, &b_GenBQuark_eta);
   fChain->SetBranchAddress("GenBQuark_phi", GenBQuark_phi, &b_GenBQuark_phi);
   fChain->SetBranchAddress("GenBQuark_mass", GenBQuark_mass, &b_GenBQuark_mass);
   fChain->SetBranchAddress("GenBQuark_pdgId", GenBQuark_pdgId, &b_GenBQuark_pdgId);
   fChain->SetBranchAddress("GenBQuark_charge", GenBQuark_charge, &b_GenBQuark_charge);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet_eta", Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_phi", Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_mass", Jet_mass, &b_Jet_mass);
   fChain->SetBranchAddress("Jet_btagCSV", Jet_btagCSV, &b_Jet_btagCSV);
   fChain->SetBranchAddress("Jet_rawPt", Jet_rawPt, &b_Jet_rawPt);
   fChain->SetBranchAddress("Jet_mcPt", Jet_mcPt, &b_Jet_mcPt);
   fChain->SetBranchAddress("Jet_mcFlavour", Jet_mcFlavour, &b_Jet_mcFlavour);
   fChain->SetBranchAddress("Jet_quarkGluonID", Jet_quarkGluonID, &b_Jet_quarkGluonID);
   fChain->SetBranchAddress("Jet_mcMatchId", Jet_mcMatchId, &b_Jet_mcMatchId);
   fChain->SetBranchAddress("Jet_mcMatchFlav", Jet_mcMatchFlav, &b_Jet_mcMatchFlav);
   fChain->SetBranchAddress("Jet_PuId_full", Jet_PuId_full, &b_Jet_PuId_full);
   fChain->SetBranchAddress("Jet_PuId_simple", Jet_PuId_simple, &b_Jet_PuId_simple);
   fChain->SetBranchAddress("Jet_PuId_cut_based", Jet_PuId_cut_based, &b_Jet_PuId_cut_based);
   fChain->SetBranchAddress("Jet_Id", Jet_Id, &b_Jet_Id);
   fChain->SetBranchAddress("nGenTop", &nGenTop, &b_nGenTop);
   fChain->SetBranchAddress("GenTop_pt", GenTop_pt, &b_GenTop_pt);
   fChain->SetBranchAddress("GenTop_eta", GenTop_eta, &b_GenTop_eta);
   fChain->SetBranchAddress("GenTop_phi", GenTop_phi, &b_GenTop_phi);
   fChain->SetBranchAddress("GenTop_mass", GenTop_mass, &b_GenTop_mass);
   fChain->SetBranchAddress("GenTop_pdgId", GenTop_pdgId, &b_GenTop_pdgId);
   fChain->SetBranchAddress("GenTop_charge", GenTop_charge, &b_GenTop_charge);
   fChain->SetBranchAddress("nGenP6StatusThree", &nGenP6StatusThree, &b_nGenP6StatusThree);
   fChain->SetBranchAddress("GenP6StatusThree_pt", GenP6StatusThree_pt, &b_GenP6StatusThree_pt);
   fChain->SetBranchAddress("GenP6StatusThree_eta", GenP6StatusThree_eta, &b_GenP6StatusThree_eta);
   fChain->SetBranchAddress("GenP6StatusThree_phi", GenP6StatusThree_phi, &b_GenP6StatusThree_phi);
   fChain->SetBranchAddress("GenP6StatusThree_mass", GenP6StatusThree_mass, &b_GenP6StatusThree_mass);
   fChain->SetBranchAddress("GenP6StatusThree_pdgId", GenP6StatusThree_pdgId, &b_GenP6StatusThree_pdgId);
   fChain->SetBranchAddress("GenP6StatusThree_charge", GenP6StatusThree_charge, &b_GenP6StatusThree_charge);
   fChain->SetBranchAddress("GenP6StatusThree_motherId", GenP6StatusThree_motherId, &b_GenP6StatusThree_motherId);
   fChain->SetBranchAddress("nLepGood", &nLepGood, &b_nLepGood);
   fChain->SetBranchAddress("LepGood_pt", LepGood_pt, &b_LepGood_pt);
   fChain->SetBranchAddress("LepGood_eta", LepGood_eta, &b_LepGood_eta);
   fChain->SetBranchAddress("LepGood_phi", LepGood_phi, &b_LepGood_phi);
   fChain->SetBranchAddress("LepGood_mass", LepGood_mass, &b_LepGood_mass);
   fChain->SetBranchAddress("LepGood_pdgId", LepGood_pdgId, &b_LepGood_pdgId);
   fChain->SetBranchAddress("LepGood_charge", LepGood_charge, &b_LepGood_charge);
   fChain->SetBranchAddress("LepGood_dxy", LepGood_dxy, &b_LepGood_dxy);
   fChain->SetBranchAddress("LepGood_dz", LepGood_dz, &b_LepGood_dz);
   fChain->SetBranchAddress("LepGood_edxy", LepGood_edxy, &b_LepGood_edxy);
   fChain->SetBranchAddress("LepGood_edz", LepGood_edz, &b_LepGood_edz);
   fChain->SetBranchAddress("LepGood_ip3d", LepGood_ip3d, &b_LepGood_ip3d);
   fChain->SetBranchAddress("LepGood_sip3d", LepGood_sip3d, &b_LepGood_sip3d);
   fChain->SetBranchAddress("LepGood_tightId", LepGood_tightId, &b_LepGood_tightId);
   fChain->SetBranchAddress("LepGood_convVeto", LepGood_convVeto, &b_LepGood_convVeto);
   fChain->SetBranchAddress("LepGood_lostHits", LepGood_lostHits, &b_LepGood_lostHits);
   fChain->SetBranchAddress("LepGood_looseIdSusy", LepGood_looseIdSusy, &b_LepGood_looseIdSusy);
   fChain->SetBranchAddress("LepGood_relIso03", LepGood_relIso03, &b_LepGood_relIso03);
   fChain->SetBranchAddress("LepGood_relIso04", LepGood_relIso04, &b_LepGood_relIso04);
   fChain->SetBranchAddress("LepGood_chargedHadRelIso03", LepGood_chargedHadRelIso03, &b_LepGood_chargedHadRelIso03);
   fChain->SetBranchAddress("LepGood_chargedHadRelIso04", LepGood_chargedHadRelIso04, &b_LepGood_chargedHadRelIso04);
   fChain->SetBranchAddress("LepGood_convVetoFull", LepGood_convVetoFull, &b_LepGood_convVetoFull);
   fChain->SetBranchAddress("LepGood_eleCutId", LepGood_eleCutId, &b_LepGood_eleCutId);
   fChain->SetBranchAddress("LepGood_eleMVAId", LepGood_eleMVAId, &b_LepGood_eleMVAId);
   fChain->SetBranchAddress("LepGood_tightCharge", LepGood_tightCharge, &b_LepGood_tightCharge);
   fChain->SetBranchAddress("LepGood_mvaId", LepGood_mvaId, &b_LepGood_mvaId);
   fChain->SetBranchAddress("LepGood_mvaIdTrig", LepGood_mvaIdTrig, &b_LepGood_mvaIdTrig);
   fChain->SetBranchAddress("LepGood_nStations", LepGood_nStations, &b_LepGood_nStations);
   fChain->SetBranchAddress("LepGood_trkKink", LepGood_trkKink, &b_LepGood_trkKink);
   fChain->SetBranchAddress("LepGood_trackerLayers", LepGood_trackerLayers, &b_LepGood_trackerLayers);
   fChain->SetBranchAddress("LepGood_pixelLayers", LepGood_pixelLayers, &b_LepGood_pixelLayers);
   fChain->SetBranchAddress("LepGood_mvaTTH", LepGood_mvaTTH, &b_LepGood_mvaTTH);
   fChain->SetBranchAddress("LepGood_jetPtRatio", LepGood_jetPtRatio, &b_LepGood_jetPtRatio);
   fChain->SetBranchAddress("LepGood_jetBTagCSV", LepGood_jetBTagCSV, &b_LepGood_jetBTagCSV);
   fChain->SetBranchAddress("LepGood_jetDR", LepGood_jetDR, &b_LepGood_jetDR);
   fChain->SetBranchAddress("LepGood_mcMatchId", LepGood_mcMatchId, &b_LepGood_mcMatchId);
   fChain->SetBranchAddress("LepGood_mcMatchAny", LepGood_mcMatchAny, &b_LepGood_mcMatchAny);
   fChain->SetBranchAddress("LepGood_mcMatchTau", LepGood_mcMatchTau, &b_LepGood_mcMatchTau);
   Notify();
}

Bool_t susySoftLepton::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void susySoftLepton::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t susySoftLepton::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef susySoftLepton_cxx
