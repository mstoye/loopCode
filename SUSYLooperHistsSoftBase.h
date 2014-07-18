//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Feb 27 16:22:54 2014 by ROOT version 5.32/00
// from TTree SUSYLooperHistsSoftBase/SUSYLooperHistsSoftBase
// found on file: /afs/cern.ch/work/a/agapitos/public/CMG/CMSSW_5_3_12_patch3/src/CMGTools/TTHAnalysis/cfg/RPV_1200_500/AtoniousSUSY/SUSYLooperHistsSoftBase/SUSYLooperHistsSoftBase_tree.root
//////////////////////////////////////////////////////////

#ifndef SUSYLooperHistsSoftBase_h
#define SUSYLooperHistsSoftBase_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "SampleInfo.h" 
#include <iostream>
#include <TLorentzVector.h>


// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class SUSYLooperHistsSoftBase {

 public :
  TString outputFileName;
  float weight;
  unsigned int maxEvents;
  float  Mstop;
  float  MLSP;

public :

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           run;
   Int_t           lumi;
   Int_t           evt;
   Int_t           HLT_MetTrigger;
   Int_t           HLT_SingleMuNonIso;
   Int_t           HLT_MetTriggerParked;
   Int_t           HLT_SingleMu;
   Int_t           HLT_MuEG;
   Int_t           HLT_TripleEl;
   Int_t           HLT_DoubleEl;
   Int_t           HLT_DoubleMu;
   Float_t         puWeight;
   Int_t           nVert;
   Int_t           nJet30;
   Int_t           nBJetLoose30;
   Int_t           nBJetMedium30;
   Int_t           nBJetTight30;
   Float_t         LepEff_1lep;
   Float_t         LepEff_2lep;
   Float_t         LepEff_3lep;
   Float_t         LepEff_4lep;
   Int_t           GenSusyMStop;
   Int_t           GenSusyMNeutralino;
   Float_t         met_pt;
   Float_t         met_eta;
   Float_t         met_phi;
   Float_t         met_mass;
   Float_t         met_sumEt;
   Float_t         met_genPt;
   Float_t         met_genPhi;
   Float_t         met_genEta;
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
   Float_t         LepGood_relIso[5];   //[nLepGood]
   Float_t         LepGood_chargedRelIso[5];   //[nLepGood]
   Int_t           LepGood_tightId[5];   //[nLepGood]
   Int_t           LepGood_convVeto[5];   //[nLepGood]
   Int_t           LepGood_lostHits[5];   //[nLepGood]
   Int_t           nGenLepFromTau;
   Float_t         GenLepFromTau_pt[2];   //[nGenLepFromTau]
   Float_t         GenLepFromTau_eta[2];   //[nGenLepFromTau]
   Float_t         GenLepFromTau_phi[2];   //[nGenLepFromTau]
   Float_t         GenLepFromTau_mass[2];   //[nGenLepFromTau]
   Int_t           GenLepFromTau_pdgId[2];   //[nGenLepFromTau]
   Float_t         GenLepFromTau_charge[2];   //[nGenLepFromTau]
   Int_t           GenLepFromTau_sourceId[2];   //[nGenLepFromTau]
   Int_t           nGenStatusThree;
   Float_t         GenStatusThree_pt[20];   //[nGenStatusThree]
   Float_t         GenStatusThree_eta[20];   //[nGenStatusThree]
   Float_t         GenStatusThree_phi[20];   //[nGenStatusThree]
   Float_t         GenStatusThree_mass[20];   //[nGenStatusThree]
   Int_t           GenStatusThree_pdgId[20];   //[nGenStatusThree]
   Float_t         GenStatusThree_charge[20];   //[nGenStatusThree]
   Int_t           GenStatusThree_sourceId[20];   //[nGenStatusThree]
   Int_t           nGenLep;
   Float_t         GenLep_pt[2];   //[nGenLep]
   Float_t         GenLep_eta[2];   //[nGenLep]
   Float_t         GenLep_phi[2];   //[nGenLep]
   Float_t         GenLep_mass[2];   //[nGenLep]
   Int_t           GenLep_pdgId[2];   //[nGenLep]
   Float_t         GenLep_charge[2];   //[nGenLep]
   Int_t           GenLep_sourceId[2];   //[nGenLep]
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
   Int_t           nJet;
   Float_t         Jet_pt[12];   //[nJet]
   Float_t         Jet_eta[12];   //[nJet]
   Float_t         Jet_phi[12];   //[nJet]
   Float_t         Jet_mass[12];   //[nJet]
   Float_t         Jet_btagCSV[12];   //[nJet]
   Float_t         Jet_rawPt[12];   //[nJet]
   Float_t         Jet_mcPt[12];   //[nJet]
   Int_t           Jet_mcFlavour[12];   //[nJet]
   Float_t         Jet_quarkGluonID[12];   //[nJet]
   Int_t           Jet_PuJetId_full[12];   //[nJet]
   Int_t           Jet_PuJetId_simple[12];   //[nJet]
   Int_t           Jet_PuJetId_cut_based[12];   //[nJet]
   Int_t           Jet_mcMatchId[12];   //[nJet]
   Int_t           Jet_mcMatchFlav[12];   //[nJet]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_HLT_MetTrigger;   //!
   TBranch        *b_HLT_SingleMuNonIso;   //!
   TBranch        *b_HLT_MetTriggerParked;   //!
   TBranch        *b_HLT_SingleMu;   //!
   TBranch        *b_HLT_MuEG;   //!
   TBranch        *b_HLT_TripleEl;   //!
   TBranch        *b_HLT_DoubleEl;   //!
   TBranch        *b_HLT_DoubleMu;   //!
   TBranch        *b_puWeight;   //!
   TBranch        *b_nVert;   //!
   TBranch        *b_nJet30;   //!
   TBranch        *b_nBJetLoose30;   //!
   TBranch        *b_nBJetMedium30;   //!
   TBranch        *b_nBJetTight30;   //!
   TBranch        *b_LepEff_1lep;   //!
   TBranch        *b_LepEff_2lep;   //!
   TBranch        *b_LepEff_3lep;   //!
   TBranch        *b_LepEff_4lep;   //!
   TBranch        *b_GenSusyMStop;   //!
   TBranch        *b_GenSusyMNeutralino;   //!
   TBranch        *b_met_pt;   //!
   TBranch        *b_met_eta;   //!
   TBranch        *b_met_phi;   //!
   TBranch        *b_met_mass;   //!
   TBranch        *b_met_sumEt;   //!
   TBranch        *b_met_genPt;   //!
   TBranch        *b_met_genPhi;   //!
   TBranch        *b_met_genEta;   //!
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
   TBranch        *b_LepGood_relIso;   //!
   TBranch        *b_LepGood_chargedRelIso;   //!
   TBranch        *b_LepGood_tightId;   //!
   TBranch        *b_LepGood_convVeto;   //!
   TBranch        *b_LepGood_lostHits;   //!
   TBranch        *b_nGenLepFromTau;   //!
   TBranch        *b_GenLepFromTau_pt;   //!
   TBranch        *b_GenLepFromTau_eta;   //!
   TBranch        *b_GenLepFromTau_phi;   //!
   TBranch        *b_GenLepFromTau_mass;   //!
   TBranch        *b_GenLepFromTau_pdgId;   //!
   TBranch        *b_GenLepFromTau_charge;   //!
   TBranch        *b_GenLepFromTau_sourceId;   //!
   TBranch        *b_nGenStatusThree;   //!
   TBranch        *b_GenStatusThree_pt;   //!
   TBranch        *b_GenStatusThree_eta;   //!
   TBranch        *b_GenStatusThree_phi;   //!
   TBranch        *b_GenStatusThree_mass;   //!
   TBranch        *b_GenStatusThree_pdgId;   //!
   TBranch        *b_GenStatusThree_charge;   //!
   TBranch        *b_GenStatusThree_sourceId;   //!
   TBranch        *b_nGenLep;   //!
   TBranch        *b_GenLep_pt;   //!
   TBranch        *b_GenLep_eta;   //!
   TBranch        *b_GenLep_phi;   //!
   TBranch        *b_GenLep_mass;   //!
   TBranch        *b_GenLep_pdgId;   //!
   TBranch        *b_GenLep_charge;   //!
   TBranch        *b_GenLep_sourceId;   //!
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
   TBranch        *b_Jet_PuJetId_full;   //!
   TBranch        *b_Jet_PuJetId_simple;   //!
   TBranch        *b_Jet_PuJetId_cut_based;   //!
   TBranch        *b_Jet_mcMatchId;   //!
   TBranch        *b_Jet_mcMatchFlav;   //!
   SUSYLooperHistsSoftBase(SampleInfo mySample, float Mstop_ = 0, float MLSP_ = 0);
   void setTree(SampleInfo mySample, float Mstop_ = 0, float MLSP_ = 0);
   SUSYLooperHistsSoftBase();
   virtual ~SUSYLooperHistsSoftBase();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    JetCuts();
   virtual Int_t    SecondSoftLeptonCuts();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual TFile*   Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   double DiTau_InvMass( TLorentzVector Met, TLorentzVector L1, TLorentzVector L2, float Al );

};
 
#endif

#ifdef SUSYLooperHistsSoftBase_cxx


//SUSYLooperHistsSoftBase::SUSYLooperHistsSoftBase(TTree *tree) : fChain(0) 
SUSYLooperHistsSoftBase::SUSYLooperHistsSoftBase(SampleInfo mySample, float Mstop_ , float MLSP_ )
{
  setTree(mySample, Mstop_ , MLSP_ );
}

void SUSYLooperHistsSoftBase::setTree(SampleInfo mySample, float Mstop_ , float MLSP_ )
{

  TFile *f = TFile::Open(mySample.FilePath);
  TTree* tree = (TTree*) f->Get("markusTreeProducer"); 
  outputFileName = mySample.OutputFileNameTag+".root";
  weight = mySample.weight();
  Mstop= Mstop_;
  MLSP=MLSP_;
  Init(tree);
  cout << "looping over: \""<< mySample.FilePath<< "\" to put plots into: \""<< outputFileName<<"\""<<endl;;

}

SUSYLooperHistsSoftBase::SUSYLooperHistsSoftBase( )
{

  cout << "Empty contructor: should not be used"<<endl;

}


SUSYLooperHistsSoftBase::~SUSYLooperHistsSoftBase()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t SUSYLooperHistsSoftBase::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t SUSYLooperHistsSoftBase::LoadTree(Long64_t entry)
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

void SUSYLooperHistsSoftBase::Init(TTree *tree)
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
   fChain->SetBranchAddress("HLT_MetTrigger", &HLT_MetTrigger, &b_HLT_MetTrigger);
   fChain->SetBranchAddress("HLT_SingleMuNonIso", &HLT_SingleMuNonIso, &b_HLT_SingleMuNonIso);
   fChain->SetBranchAddress("HLT_MetTriggerParked", &HLT_MetTriggerParked, &b_HLT_MetTriggerParked);
   fChain->SetBranchAddress("HLT_SingleMu", &HLT_SingleMu, &b_HLT_SingleMu);
   fChain->SetBranchAddress("HLT_MuEG", &HLT_MuEG, &b_HLT_MuEG);
   fChain->SetBranchAddress("HLT_TripleEl", &HLT_TripleEl, &b_HLT_TripleEl);
   fChain->SetBranchAddress("HLT_DoubleEl", &HLT_DoubleEl, &b_HLT_DoubleEl);
   fChain->SetBranchAddress("HLT_DoubleMu", &HLT_DoubleMu, &b_HLT_DoubleMu);
   fChain->SetBranchAddress("puWeight", &puWeight, &b_puWeight);
   fChain->SetBranchAddress("nVert", &nVert, &b_nVert);
   fChain->SetBranchAddress("nJet30", &nJet30, &b_nJet30);
   fChain->SetBranchAddress("nBJetLoose30", &nBJetLoose30, &b_nBJetLoose30);
   fChain->SetBranchAddress("nBJetMedium30", &nBJetMedium30, &b_nBJetMedium30);
   fChain->SetBranchAddress("nBJetTight30", &nBJetTight30, &b_nBJetTight30);
   fChain->SetBranchAddress("LepEff_1lep", &LepEff_1lep, &b_LepEff_1lep);
   fChain->SetBranchAddress("LepEff_2lep", &LepEff_2lep, &b_LepEff_2lep);
   fChain->SetBranchAddress("LepEff_3lep", &LepEff_3lep, &b_LepEff_3lep);
   fChain->SetBranchAddress("LepEff_4lep", &LepEff_4lep, &b_LepEff_4lep);
   fChain->SetBranchAddress("GenSusyMStop", &GenSusyMStop, &b_GenSusyMStop);
   fChain->SetBranchAddress("GenSusyMNeutralino", &GenSusyMNeutralino, &b_GenSusyMNeutralino);
   fChain->SetBranchAddress("met_pt", &met_pt, &b_met_pt);
   fChain->SetBranchAddress("met_eta", &met_eta, &b_met_eta);
   fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
   fChain->SetBranchAddress("met_mass", &met_mass, &b_met_mass);
   fChain->SetBranchAddress("met_sumEt", &met_sumEt, &b_met_sumEt);
   fChain->SetBranchAddress("met_genPt", &met_genPt, &b_met_genPt);
   fChain->SetBranchAddress("met_genPhi", &met_genPhi, &b_met_genPhi);
   fChain->SetBranchAddress("met_genEta", &met_genEta, &b_met_genEta);
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
   fChain->SetBranchAddress("LepGood_relIso", LepGood_relIso, &b_LepGood_relIso);
   fChain->SetBranchAddress("LepGood_chargedRelIso", LepGood_chargedRelIso, &b_LepGood_chargedRelIso);
   fChain->SetBranchAddress("LepGood_tightId", LepGood_tightId, &b_LepGood_tightId);
   fChain->SetBranchAddress("LepGood_convVeto", LepGood_convVeto, &b_LepGood_convVeto);
   fChain->SetBranchAddress("LepGood_lostHits", LepGood_lostHits, &b_LepGood_lostHits);
   fChain->SetBranchAddress("nGenLepFromTau", &nGenLepFromTau, &b_nGenLepFromTau);
   fChain->SetBranchAddress("GenLepFromTau_pt", GenLepFromTau_pt, &b_GenLepFromTau_pt);
   fChain->SetBranchAddress("GenLepFromTau_eta", GenLepFromTau_eta, &b_GenLepFromTau_eta);
   fChain->SetBranchAddress("GenLepFromTau_phi", GenLepFromTau_phi, &b_GenLepFromTau_phi);
   fChain->SetBranchAddress("GenLepFromTau_mass", GenLepFromTau_mass, &b_GenLepFromTau_mass);
   fChain->SetBranchAddress("GenLepFromTau_pdgId", GenLepFromTau_pdgId, &b_GenLepFromTau_pdgId);
   fChain->SetBranchAddress("GenLepFromTau_charge", GenLepFromTau_charge, &b_GenLepFromTau_charge);
   fChain->SetBranchAddress("GenLepFromTau_sourceId", GenLepFromTau_sourceId, &b_GenLepFromTau_sourceId);
   fChain->SetBranchAddress("nGenStatusThree", &nGenStatusThree, &b_nGenStatusThree);
   fChain->SetBranchAddress("GenStatusThree_pt", GenStatusThree_pt, &b_GenStatusThree_pt);
   fChain->SetBranchAddress("GenStatusThree_eta", GenStatusThree_eta, &b_GenStatusThree_eta);
   fChain->SetBranchAddress("GenStatusThree_phi", GenStatusThree_phi, &b_GenStatusThree_phi);
   fChain->SetBranchAddress("GenStatusThree_mass", GenStatusThree_mass, &b_GenStatusThree_mass);
   fChain->SetBranchAddress("GenStatusThree_pdgId", GenStatusThree_pdgId, &b_GenStatusThree_pdgId);
   fChain->SetBranchAddress("GenStatusThree_charge", GenStatusThree_charge, &b_GenStatusThree_charge);
   fChain->SetBranchAddress("GenStatusThree_sourceId", GenStatusThree_sourceId, &b_GenStatusThree_sourceId);
   fChain->SetBranchAddress("nGenLep", &nGenLep, &b_nGenLep);
   fChain->SetBranchAddress("GenLep_pt", GenLep_pt, &b_GenLep_pt);
   fChain->SetBranchAddress("GenLep_eta", GenLep_eta, &b_GenLep_eta);
   fChain->SetBranchAddress("GenLep_phi", GenLep_phi, &b_GenLep_phi);
   fChain->SetBranchAddress("GenLep_mass", GenLep_mass, &b_GenLep_mass);
   fChain->SetBranchAddress("GenLep_pdgId", GenLep_pdgId, &b_GenLep_pdgId);
   fChain->SetBranchAddress("GenLep_charge", GenLep_charge, &b_GenLep_charge);
   fChain->SetBranchAddress("GenLep_sourceId", GenLep_sourceId, &b_GenLep_sourceId);
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
   fChain->SetBranchAddress("Jet_PuJetId_full", Jet_PuJetId_full, &b_Jet_PuJetId_full);
   fChain->SetBranchAddress("Jet_PuJetId_simple", Jet_PuJetId_simple, &b_Jet_PuJetId_simple);
   fChain->SetBranchAddress("Jet_PuJetId_cut_based", Jet_PuJetId_cut_based, &b_Jet_PuJetId_cut_based);
   fChain->SetBranchAddress("Jet_mcMatchId", Jet_mcMatchId, &b_Jet_mcMatchId);
   fChain->SetBranchAddress("Jet_mcMatchFlav", Jet_mcMatchFlav, &b_Jet_mcMatchFlav);


   Notify();
}

Bool_t SUSYLooperHistsSoftBase::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void SUSYLooperHistsSoftBase::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t SUSYLooperHistsSoftBase::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
  cout << " can be used to put preselection cuts for " << entry<<endl;
// returns -1 otherwise.
   return 1;
}

Int_t SUSYLooperHistsSoftBase::JetCuts()
{
// This function may be called from Loop.
 if(nJet==0) return 0;
  if(nJet>0){
    if (Jet_pt[0] < 110) return 0;
    if(fabs(Jet_eta[0]) > 2.4) return 0;
  }
  if(nJet>2){
    if (Jet_pt[2] > 60) return 0;
  }
  return 1;
}

Int_t SUSYLooperHistsSoftBase::SecondSoftLeptonCuts()
{
// This function may  be called from Loop.
  if (nLepGood<2) return 0;
  if(fabs(LepGood_dz[1]) > 0.02||fabs(LepGood_dxy[1])> 0.02) return 0;
  if(LepGood_relIso[1]*LepGood_pt[1]>5.)  return 0;
  if(LepGood_relIso[1]>.5)  return 0;
  if(LepGood_pt[1]<3||fabs(LepGood_eta[1])>1.5) return 0;
  if(fabs(LepGood_pt[1])==11&&LepGood_pt[1]<5) return 0;
  if(LepGood_pt[1]>15) return 0;
  return 1;
}




#endif // #ifdef SUSYLooperHistsSoftBase_cxx
