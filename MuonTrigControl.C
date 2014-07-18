#define MuonTrigControl_cxx
#include "MuonTrigControl.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include  <vector>
#include <TF1.h>
#include <TMath.h>


TFile* MuonTrigControl::Loop()
{
  TFile* myFile = new TFile(outputFileName,"recreate");
  if (fChain == 0) return myFile;
  Long64_t nentries = fChain->GetEntriesFast();

  TH1D *h_CutFlow  =new TH1D("h_CutFlow","h_CutFlow", 25, -0.5, 24.5 ); 
  h_CutFlow->Sumw2();
 
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    float HT30 =0.;
    int Nj30=0;
    int Nb=0;
   
    for(int i=0;i<nJet;i++){//Loop Over Selected Jets, measure Nj, Nb, HT;
      if(Jet_pt[i]>30.&& fabs(Jet_eta[i])<4.5 ) {
	Nj30++;
	HT30=HT30+Jet_pt[i]; 
	if(Jet_btagCSV[i]>0.679/*0.244*/)Nb++;  //CSV[M] b-tagging
      } 
    }
    if( Nb>0                )continue;  h_CutFlow->Fill(11);


    
     if( MLSP!=0 ){//-----(Cuts for scans to select point)---------- 
      
      if( (GenSusyMStop<=(Mstop-12)) || (GenSusyMStop>=(Mstop+12)) )     continue; //<---[M_stop, M_LSP: set by by the last 2 arguments of the SampleInfo objects]
      if( (GenSusyMNeutralino<=(MLSP - 4)) || (GenSusyMNeutralino>=(MLSP  +4)) )     continue; 
    }//-------[ CUTS for Non-Scan samples ]-------------------------
    if( nJet==0 &&  nLepGood<2                          )continue;  h_CutFlow->Fill(1);
    if( fabs(LepGood_pdgId[0]*LepGood_pdgId[1])!=121  )continue;  h_CutFlow->Fill(2);
    if( LepGood_pt[0]<5 || fabs(LepGood_eta[0])>1.5   )continue;  h_CutFlow->Fill(3);
    if( LepGood_pt[1]<3 || fabs(LepGood_eta[1])>1.5   )continue;  h_CutFlow->Fill(4);
    if( fabs(LepGood_pdgId[1])==11 && LepGood_pt[1]<5 )continue;  h_CutFlow->Fill(5); //Threashold for electrons.
    if( fabs(LepGood_pdgId[0])==11 && LepGood_pt[0]<7 )continue;  h_CutFlow->Fill(6); //Threashold for electrons.
    if( LepGood_relIso[1]>.5 || LepGood_relIso[1]*LepGood_pt[1]>5 )continue;  h_CutFlow->Fill(7); //2nd lep: *Universal* iso requirements!
    if( fabs(LepGood_dz[0])>0.02 || fabs(LepGood_dxy[0])>0.02 )    continue; 
    if( fabs(LepGood_dz[1])>0.02 || fabs(LepGood_dxy[1])>0.02 )    continue;  h_CutFlow->Fill(8);
    if( Jet_pt[0]<150||fabs(Jet_eta[0])>2.4           )continue;  h_CutFlow->Fill(9); //1st Jet requirements.
    if( nJet>2 && Jet_pt[2]>60                          )continue;  h_CutFlow->Fill(10);//2nd Jet requirements.
    for(int i=0;i<nJet;i++){//Loop Over Selected Jets, measure Nj, Nb, HT;
      if(Jet_pt[i]>30.&& fabs(Jet_eta[i])<4.5 ) {
	Nj30++;
	HT30=HT30+Jet_pt[i]; 
	if(Jet_btagCSV[i]>0.679/*0.244*/)Nb++;  //CSV[M] b-tagging
      } 
    }
    if( Nb>0                )continue;  h_CutFlow->Fill(11);
    
  }
  
  myFile->Write();
  myFile->Close();
  
  return (myFile);
      // if (Cut(ientry) < 0) continue;
}
