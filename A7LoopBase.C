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




