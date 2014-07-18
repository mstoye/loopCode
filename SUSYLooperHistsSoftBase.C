#define SUSYLooperHistsSoftBase_cxx
#include "SUSYLooperHistsSoftBase.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include  <vector>

#include <TF1.h>
#include <TMath.h>
#include <TMatrixF.h>
#include <TVectorF.h>

TFile* SUSYLooperHistsSoftBase::Loop()
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


double SUSYLooperHistsSoftBase::DiTau_InvMass( TLorentzVector Met, TLorentzVector L1, TLorentzVector L2, float Al ){   
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

 

