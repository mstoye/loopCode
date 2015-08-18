#ifndef A7LoopBase_cxx 
#define A7LoopBase_cxx
#include "A7LoopBase.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include  <vector>

#include <TF1.h>
#include <TMath.h>
#include <TMatrixF.h>
#include <TVectorF.h>

TFile* A7LoopBase::Loop()
{
  TFile* myFile = new TFile(outputFileName,"recreate");
   if (fChain == 0) return myFile;
   Long64_t nentries = fChain->GetEntriesFast();
  
 
  
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

    
   }


  
   return (myFile);
      // if (Cut(ientry) < 0) continue;
}




void A7LoopBase::setTree(SampleInfo mySample )
{

  TFile *f = TFile::Open(mySample.FilePath);
  f->ls();
  TTree* tree = (TTree*) f->Get("treeProducerA7W"); 
  outputFileName = mySample.OutputFileNameTag+".root";
  weight = mySample.weight();
  Init(tree);
  cout << "looping over: \""<< mySample.FilePath<< "\" to put plots into: \""<< outputFileName<<"\""<<endl;;

}




A7LoopBase::~A7LoopBase()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t A7LoopBase::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t A7LoopBase::LoadTree(Long64_t entry)
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

void A7LoopBase::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
  
  if (!tree) {
    cout << "A7LoopBase::Init: error: tree name likely wrong!" <<endl;
    return;
  } 
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   //   fChain->SetBranchAddress("run", &run, &b_run);
   //  fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   //  fChain->SetBranchAddress("evt", &evt, &b_evt);
   //   fChain->SetBranchAddress("isData", &isData, &b_isData);
   fChain->SetBranchAddress("HLT_SingleMu", &HLT_SingleMu, &b_HLT_SingleMu);
   //   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("nVert", &nVert, &b_nVert);
   fChain->SetBranchAddress("nJet25", &nJet25, &b_nJet25);
   fChain->SetBranchAddress("nBJetLoose25", &nBJetLoose25, &b_nBJetLoose25);
   fChain->SetBranchAddress("nBJetMedium25", &nBJetMedium25, &b_nBJetMedium25);
   // fChain->SetBranchAddress("nBJetTight25", &nBJetTight25, &b_nBJetTight25);
   fChain->SetBranchAddress("nJet40", &nJet40, &b_nJet40);
   // fChain->SetBranchAddress("nJet40a", &nJet40a, &b_nJet40a);
   // fChain->SetBranchAddress("nBJetLoose40", &nBJetLoose40, &b_nBJetLoose40);
   //fChain->SetBranchAddress("nBJetMedium40", &nBJetMedium40, &b_nBJetMedium40);
   //  fChain->SetBranchAddress("nBJetTight40", &nBJetTight40, &b_nBJetTight40);
   // fChain->SetBranchAddress("nLepGood20", &nLepGood20, &b_nLepGood20);
   // fChain->SetBranchAddress("nLepGood15", &nLepGood15, &b_nLepGood15);
   //fChain->SetBranchAddress("nLepGood10", &nLepGood10, &b_nLepGood10);
   // fChain->SetBranchAddress("checkEcalDead", &checkEcalDead, &b_checkEcalDead);
   // fChain->SetBranchAddress("checkhcalLaser", &checkhcalLaser, &b_checkhcalLaser);
   // fChain->SetBranchAddress("checktrackingFailure", &checktrackingFailure, &b_checktrackingFailure);
   //  fChain->SetBranchAddress("checkprimaryVertex", &checkprimaryVertex, &b_checkprimaryVertex);
   //  fChain->SetBranchAddress("checknoscraping", &checknoscraping, &b_checknoscraping);
   //  fChain->SetBranchAddress("checktrackIsolationMaker", &checktrackIsolationMaker, &b_checktrackIsolationMaker);
   // fChain->SetBranchAddress("checkmetNoiseCleaning", &checkmetNoiseCleaning, &b_checkmetNoiseCleaning);
   // fChain->SetBranchAddress("checkeeBadSc", &checkeeBadSc, &b_checkeeBadSc);
   //  fChain->SetBranchAddress("checkecalLaser", &checkecalLaser, &b_checkecalLaser);
   //  fChain->SetBranchAddress("checktotalKinematics", &checktotalKinematics, &b_checktotalKinematics);
   //  fChain->SetBranchAddress("checkCSCTightHalo", &checkCSCTightHalo, &b_checkCSCTightHalo);
   // fChain->SetBranchAddress("checkHBHENoise", &checkHBHENoise, &b_checkHBHENoise);
   //  fChain->SetBranchAddress("htJet25", &htJet25, &b_htJet25);
   // fChain->SetBranchAddress("mhtJet25", &mhtJet25, &b_mhtJet25);
   // fChain->SetBranchAddress("htJet40j", &htJet40j, &b_htJet40j);
   //  fChain->SetBranchAddress("htJet40ja", &htJet40ja, &b_htJet40ja);
   // fChain->SetBranchAddress("htJet40", &htJet40, &b_htJet40);
   //  fChain->SetBranchAddress("htJet40a", &htJet40a, &b_htJet40a);
   //  fChain->SetBranchAddress("mhtJet40", &mhtJet40, &b_mhtJet40);
   //  fChain->SetBranchAddress("mhtJet40a", &mhtJet40a, &b_mhtJet40a);
   //  fChain->SetBranchAddress("mZ1", &mZ1, &b_mZ1);
   //   fChain->SetBranchAddress("mZ1SFSS", &mZ1SFSS, &b_mZ1SFSS);
   //   fChain->SetBranchAddress("minMllSFOS", &minMllSFOS, &b_minMllSFOS);
   //   fChain->SetBranchAddress("maxMllSFOS", &maxMllSFOS, &b_maxMllSFOS);
   //  fChain->SetBranchAddress("minMllAFOS", &minMllAFOS, &b_minMllAFOS);
   //   fChain->SetBranchAddress("maxMllAFOS", &maxMllAFOS, &b_maxMllAFOS);
   //   fChain->SetBranchAddress("minMllAFSS", &minMllAFSS, &b_minMllAFSS);
   //   fChain->SetBranchAddress("maxMllAFSS", &maxMllAFSS, &b_maxMllAFSS);
   //   fChain->SetBranchAddress("minMllAFAS", &minMllAFAS, &b_minMllAFAS);
   //   fChain->SetBranchAddress("maxMllAFAS", &maxMllAFAS, &b_maxMllAFAS);
   //   fChain->SetBranchAddress("m2l", &m2l, &b_m2l);
   //   fChain->SetBranchAddress("minMWjj", &minMWjj, &b_minMWjj);
   //   //   fChain->SetBranchAddress("minMWjjPt", &minMWjjPt, &b_minMWjjPt);
   //   fChain->SetBranchAddress("bestMWjj", &bestMWjj, &b_bestMWjj);
   //   fChain->SetBranchAddress("bestMWjjPt", &bestMWjjPt, &b_bestMWjjPt);
   //   fChain->SetBranchAddress("bestMTopHad", &bestMTopHad, &b_bestMTopHad);
   //   fChain->SetBranchAddress("bestMTopHadPt", &bestMTopHadPt, &b_bestMTopHadPt);
   fChain->SetBranchAddress("met_pt", &met_pt, &b_met_pt);
   fChain->SetBranchAddress("met_eta", &met_eta, &b_met_eta);
   fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
   fChain->SetBranchAddress("met_mass", &met_mass, &b_met_mass);
   //fChain->SetBranchAddress("met_sumEt", &met_sumEt, &b_met_sumEt);
   /*   fChain->SetBranchAddress("metNoPU_pt", &metNoPU_pt, &b_metNoPU_pt);
   fChain->SetBranchAddress("metNoPU_eta", &metNoPU_eta, &b_metNoPU_eta);
   fChain->SetBranchAddress("metNoPU_phi", &metNoPU_phi, &b_metNoPU_phi);
   fChain->SetBranchAddress("metNoPU_mass", &metNoPU_mass, &b_metNoPU_mass);
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
   fChain->SetBranchAddress("LepOther_caloCompatibility", LepOther_caloCompatibility, &b_LepOther_caloCompatibility);
   fChain->SetBranchAddress("LepOther_globalTrackChi2", LepOther_globalTrackChi2, &b_LepOther_globalTrackChi2);
   fChain->SetBranchAddress("LepOther_trackerLayers", LepOther_trackerLayers, &b_LepOther_trackerLayers);
   fChain->SetBranchAddress("LepOther_pixelLayers", LepOther_pixelLayers, &b_LepOther_pixelLayers);
   fChain->SetBranchAddress("LepOther_mvaTTH", LepOther_mvaTTH, &b_LepOther_mvaTTH);
   fChain->SetBranchAddress("LepOther_jetPtRatio", LepOther_jetPtRatio, &b_LepOther_jetPtRatio);
   fChain->SetBranchAddress("LepOther_jetBTagCSV", LepOther_jetBTagCSV, &b_LepOther_jetBTagCSV);
   fChain->SetBranchAddress("LepOther_jetDR", LepOther_jetDR, &b_LepOther_jetDR);
   fChain->SetBranchAddress("LepOther_softMuID", LepOther_softMuID, &b_LepOther_softMuID);
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
   */
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet_eta", Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_phi", Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_mass", Jet_mass, &b_Jet_mass);
   fChain->SetBranchAddress("Jet_btagCSV", Jet_btagCSV, &b_Jet_btagCSV);
   fChain->SetBranchAddress("Jet_rawPt", Jet_rawPt, &b_Jet_rawPt);
   //  fChain->SetBranchAddress("Jet_quarkGluonID", Jet_quarkGluonID, &b_Jet_quarkGluonID);
   fChain->SetBranchAddress("Jet_puId", Jet_puId, &b_Jet_puId);
   fChain->SetBranchAddress("Jet_area", Jet_area, &b_Jet_area);
   fChain->SetBranchAddress("Jet_id", Jet_id, &b_Jet_id);
   fChain->SetBranchAddress("Jet_CHEF", Jet_CHEF, &b_Jet_CHEF);
   fChain->SetBranchAddress("Jet_NHEF", Jet_NHEF, &b_Jet_NHEF);
   fChain->SetBranchAddress("Jet_PHEF", Jet_PHEF, &b_Jet_PHEF);
   fChain->SetBranchAddress("Jet_MUEF", Jet_MUEF, &b_Jet_MUEF);
   fChain->SetBranchAddress("Jet_ELEF", Jet_ELEF, &b_Jet_ELEF);
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
   /*   fChain->SetBranchAddress("LepGood_chargedHadRelIso03", LepGood_chargedHadRelIso03, &b_LepGood_chargedHadRelIso03);
   fChain->SetBranchAddress("LepGood_chargedHadRelIso04", LepGood_chargedHadRelIso04, &b_LepGood_chargedHadRelIso04);
   fChain->SetBranchAddress("LepGood_convVetoFull", LepGood_convVetoFull, &b_LepGood_convVetoFull);
   fChain->SetBranchAddress("LepGood_eleCutId", LepGood_eleCutId, &b_LepGood_eleCutId);
   fChain->SetBranchAddress("LepGood_eleMVAId", LepGood_eleMVAId, &b_LepGood_eleMVAId);
   fChain->SetBranchAddress("LepGood_tightCharge", LepGood_tightCharge, &b_LepGood_tightCharge);
   fChain->SetBranchAddress("LepGood_mvaId", LepGood_mvaId, &b_LepGood_mvaId);
   fChain->SetBranchAddress("LepGood_mvaIdTrig", LepGood_mvaIdTrig, &b_LepGood_mvaIdTrig);
   fChain->SetBranchAddress("LepGood_nStations", LepGood_nStations, &b_LepGood_nStations);
   fChain->SetBranchAddress("LepGood_trkKink", LepGood_trkKink, &b_LepGood_trkKink);
   fChain->SetBranchAddress("LepGood_caloCompatibility", LepGood_caloCompatibility, &b_LepGood_caloCompatibility);
   fChain->SetBranchAddress("LepGood_globalTrackChi2", LepGood_globalTrackChi2, &b_LepGood_globalTrackChi2);
   fChain->SetBranchAddress("LepGood_trackerLayers", LepGood_trackerLayers, &b_LepGood_trackerLayers);
   fChain->SetBranchAddress("LepGood_pixelLayers", LepGood_pixelLayers, &b_LepGood_pixelLayers);
   fChain->SetBranchAddress("LepGood_mvaTTH", LepGood_mvaTTH, &b_LepGood_mvaTTH);
   fChain->SetBranchAddress("LepGood_jetPtRatio", LepGood_jetPtRatio, &b_LepGood_jetPtRatio);
   fChain->SetBranchAddress("LepGood_jetBTagCSV", LepGood_jetBTagCSV, &b_LepGood_jetBTagCSV);
   fChain->SetBranchAddress("LepGood_jetDR", LepGood_jetDR, &b_LepGood_jetDR);
   fChain->SetBranchAddress("LepGood_softMuID", LepGood_softMuID, &b_LepGood_softMuID);*/


   Notify();
}

Bool_t A7LoopBase::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void A7LoopBase::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t A7LoopBase::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
  cout << " can be used to put preselection cuts for " << entry<<endl;
// returns -1 otherwise.
   return 1;
}

Int_t A7LoopBase::CutLepGenAcc(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
  cout << " can be used to put preselection cuts for " << entry<<endl;
// returns -1 otherwise.
   return 1;
}
 
Int_t A7LoopBase::CutRECO()
{
// This function may be called from Loop.
  int keep = 1;

  if(nLepGood==0)  keep = -1;
  //if(!HLT_SingleMu) continue;
  if(LepGood_pt[0]<25) keep = -1;
  if(fabs(LepGood_pdgId[0])!=13) keep = -1;
  if(LepGood_relIso04[0]>0.12) keep = -1;
  if( fabs(LepGood_eta[0]) >2.1) keep = -1;
  if(!LepGood_tightId[0]) keep = -1;
 
  if(nLepGood> 1) {
    for (int i =1 ;i<nLepGood;i++)
      {
	if(LepGood_pt[i]>10&& fabs(LepGood_eta[i])<2.4 && LepGood_relIso04[i]<0.2)  keep = -1;
      }
  }
    if(nJet==0) keep = -1;
    if(Jet_pt[0]<30) keep = -1;
    if(fabs(Jet_eta[0])>2.5) keep = -1;
    if( nJet40 > 2 ) keep = -1;
    if(nBJetLoose25>0) keep = -1;

// returns -1 otherwise.
   return keep;
  }
  


#endif

