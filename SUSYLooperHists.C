#define SUSYLooperHists_cxx
#include "SUSYLooperHists.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include  <vector>

#include <TF1.h>
#include <TMath.h>
#include <TMatrixF.h>
#include <TVectorF.h>

TFile* SUSYLooperHists::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L SUSYLooperHists.C
//      Root > SUSYLooperHists t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch


   if (fChain == 0) return NULL;
   Long64_t nentries = fChain->GetEntriesFast();
  
   TFile* myFile = new TFile(outputFileName,"recreate");
   // checks the Mtt mass ++ simple plot example ++
   TH1D* LepTwoZmass= new TH1D("LepTwoZmass",";approximate Z(#tau#tau) mass [GeV];",100,0,2000);
   LepTwoZmass->Sumw2();

   // baseline plots +++++++++++++++++++++
   // allow to fill cutflow to sychronize to others
   TH1D* CutFlow = new TH1D("CutFlow",";CutFlow [unweighted];",20,-0.5,19.5);
   CutFlow->Sumw2();

 // a scan should be filled without weights before any cut to get the efficiency. The reason is that for a scan each points have different x-section, which are not accessable during the loop
   TH2D* scan = new TH2D("scan","scan",8,112.5,312.5,40,112.5,312.5);
   // fill after some cuts the scans to get efficiency, i.e. ->Divide(scan) after the loop
   TH2D* scanA = new TH2D("scanA","scan",8,112.5,312.5,40,112.5,312.5);
   TH2D* scanB = new TH2D("scanB","scan",8,112.5,312.5,40,112.5,312.5);
   TH2D* scanC = new TH2D("scanC","scan",8,112.5,312.5,40,112.5,312.5);


   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if(jentry%100000==0) cout << "Event: "<<jentry <<endl;
      scan->Fill(GenSusyMStop,GenSusyMNeutralino);
      // select lepton and put it into TLorentzvector 

      TLorentzVector lep;
      TLorentzVector met;
      vector <TLorentzVector> lepVec;
      for(int i =0;i<nLepGood;i++)
	{
	  TLorentzVector aLep;
	  aLep.SetPtEtaPhiM(LepGood_pt[i],LepGood_eta[i],LepGood_phi[i],0.1);
	  lepVec.push_back(aLep);
	}
      if(nLepGood>0) lep.SetPtEtaPhiM(LepGood_pt[0],LepGood_eta[0],LepGood_phi[0],0.1);
      else lep.SetPtEtaPhiM(0,0,0,0.1);
      met.SetPtEtaPhiM(met_pt,0.,met_phi,0);
      float MT =  sqrt(lep.Pt()*2*met.Pt()*(1-cos(lep.DeltaPhi(met))));
      TLorentzVector jetSum;
      vector <TLorentzVector> myVecJets;
      jetSum.SetPtEtaPhiM(0,0,0,0);
      unsigned int nb40=0;
      unsigned int nb30L=0;

      float HT30 = 0;
      for (int i=0;i<nJet;i++)
	{
	  TLorentzVector myjet;
	  //cout <<" pt: " << Jet_pt[i] <<  " "<< Jet_eta[i]<< " "<<  Jet_phi[i]<< " QG "<< Jet_quarkGluonID[i] <<endl;
	  myjet.SetPtEtaPhiM(Jet_pt[i],Jet_eta[i],Jet_phi[i],Jet_mass[i]);
	  jetSum=jetSum+myjet;
	  myVecJets.push_back(myjet);
	  if(Jet_btagCSV[i]>0.679) {
	   
	                      nb40++;
	  }

	  if(Jet_btagCSV[i]>0.244) {
	  
	    nb30L++;
	  }

	  if(fabs(Jet_eta[i])<4.5&&Jet_pt[i]>30) HT30=HT30+Jet_pt[i];

	  //	   QG->Fill(Jet_quarkGluonID[i]);
	}

      // preselection syncronizd to Antonis ++++++++++++++++++++
      if( nLepGood == 0  ) continue;
      CutFlow->Fill(0.);
      if(met.Pt()<200) continue;
      if(nLepGood<2)  continue;
      if(nJet==0)  continue;
      if(nTauGood!=0)  continue;
      CutFlow->Fill(1.);
      if(fabs(LepGood_pdgId[0]* LepGood_pdgId[1])<122)   continue;
      CutFlow->Fill(2.);
      if(lep.Pt()<5||fabs(lep.Eta())>1.5||lep.Pt()>60) continue;
      CutFlow->Fill(3.);
      if(fabs(LepGood_pdgId[0])==11&&lep.Pt()<7) continue;
      CutFlow->Fill(4.);
      // check that a cut was ahead in nLepton !!!!!!!!
      TLorentzVector secondLep;
      secondLep.SetPtEtaPhiM(LepGood_pt[1], LepGood_eta[1], LepGood_phi[1],0.1);
      if(secondLep.Pt()<3||fabs(secondLep.Eta())>1.5||secondLep.Pt()>60) continue;
      CutFlow->Fill(5.);
      // lepton cuts
      if(LepGood_relIso[0]*lep.Pt()>10.) continue;
      if(LepGood_relIso[1]*secondLep.Pt()>5.) continue;
      if(LepGood_relIso[1]>.5) continue;
      CutFlow->Fill(6.);
      // JET CUTS
      if(nJet>0){
	if (Jet_pt[0] < 150) continue;
	if(fabs(Jet_eta[0]) > 2.4) continue;
      }
      if(nJet>2){
	if (Jet_pt[2] > 60) continue;
      }
      CutFlow->Fill(7.);
      if(met.Pt()/HT30<2./3.) continue;
      if(nb40!=0) continue;    
      CutFlow->Fill(9.);
    // preselection syncronizd to  ++++++++++++++++++++

      if(MT>60&&MT<100) continue;
 

      float pairmass = DiTau_InvMass(met,lep,secondLep,0);     
      if(lep.Pt()>25||secondLep.Pt()>15) continue;
 
      // define some over and underflow 
      if (pairmass>2000) pairmass=1999.;
      if (pairmass<0) pairmass=0.;

      LepTwoZmass->Fill(pairmass,weight*puWeight);


	if((LepGood_pdgId[0]* LepGood_pdgId[1]<-121))
	  {
	 
	    if(pairmass>160||pairmass<20){
	      CutFlow->Fill(10., weight*puWeight);
	      scanB->Fill(GenSusyMStop,GenSusyMNeutralino);
	    
	      if(nb30L==0)
		{	
		  scanC->Fill(GenSusyMStop,GenSusyMNeutralino);
		}
	      for (int i=0;i<nJet;i++)
		{
		  if(Jet_btagCSV[i]>0.244 ) {
		    // LooseBPt->Fill(Jet_pt[i],weight );		    
		  }
		  //	   QG->Fill(Jet_quarkGluonID[i]);
		}

	  }


	else // same sign leptons
	  {	  
	   
	  }
      }
          
   }


   scanA->Divide(scan);
   scanB->Divide(scan);
   scanC->Divide(scan);

   CutFlow->DrawCopy();

   myFile->Write();
   return (myFile);
      // if (Cut(ientry) < 0) continue;
}

double SUSYLooperHists::DiTau_InvMass( TLorentzVector Met, TLorentzVector L1, TLorentzVector L2, float Al ){   
  TLorentzVector T1,T2;
  double DiTauMass; 
  TMatrixF A(2,2); 
  TVectorF C(2),X(2); 
  A(0,0)=L1.Px();
  A(0,1)=L2.Px();
  A(1,0)=L1.Py();
  A(1,1)=L2.Py();
  A=A.Invert();
  C(0)=(Met+L1+L2).Px();
  C(1)=(Met+L1+L2).Py();
  X=A*C;// double X0i=X(0), X1i=X(1);
  //---------------[ MET ReAlignement subsection ]------------------------------
  if(X(0)<0||X(1)<0){ 
    if     ( fabs(L1.DeltaPhi(Met))>Al && fabs(L2.DeltaPhi(Met))>Al                                                  ) {/*DO NOTHING just normaly a non-Z event!*/}
    else if( fabs(L1.DeltaPhi(Met))<Al && fabs(L2.DeltaPhi(Met))>Al                                                  ) Met.SetPtEtaPhiM(Met.Pt(),0,L1.Phi(),0);
    else if( fabs(L2.DeltaPhi(Met))<Al && fabs(L1.DeltaPhi(Met))>Al                                                  ) Met.SetPtEtaPhiM(Met.Pt(),0,L2.Phi(),0);
    else if( fabs(L1.DeltaPhi(Met))<Al && fabs(L2.DeltaPhi(Met))<Al && fabs(L1.DeltaPhi(Met))<fabs(L2.DeltaPhi(Met)) ) Met.SetPtEtaPhiM(Met.Pt(),0,L1.Phi(),0);
    else if( fabs(L1.DeltaPhi(Met))<Al && fabs(L2.DeltaPhi(Met))<Al && fabs(L1.DeltaPhi(Met))>fabs(L2.DeltaPhi(Met)) ) Met.SetPtEtaPhiM(Met.Pt(),0,L2.Phi(),0);
  }//---------------------------------------------------------------------------
  C(0)=(Met+L1+L2).Px(); C(1)=(Met+L1+L2).Py(); X=A*C;
  T1.SetPxPyPzE( L1.Px()*X(0), L1.Py()*X(0), L1.Pz()*X(0), sqrt( 3.1571 +L1.P()*L1.P()*X(0)*X(0) ) );
  T2.SetPxPyPzE( L2.Px()*X(1), L2.Py()*X(1), L2.Pz()*X(1), sqrt( 3.1571 +L2.P()*L2.P()*X(1)*X(1) ) );
  if( X(0)>0 && X(1)>0 ) DiTauMass=(T1+T2).M();  else DiTauMass=-(T1+T2).M();  return DiTauMass;
  //if((X(0)!=X0i||X(1)!=X1i))std::cout<<X(0)<<" "<<X(1)<<" <--"<<X0i<<" "<<X1i<<" RMETal.phi="<<(T1+T2-L1-L2).Phi()<<" RMETal.eta"<<(T1+T2-L1-L2).Eta()<<" MZ="<<DiTauMass<<endl; 
}

 
