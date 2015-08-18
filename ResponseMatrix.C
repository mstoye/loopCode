#define ResponseMatrix_cxx
#include "ResponseMatrix.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include  <vector>
#include <TF1.h>
#include <TMath.h>
#include <iostream>


TFile* ResponseMatrix::Loop()
{
  TFile* myFile = new TFile(outputFileName,"recreate");
  if (fChain == 0) return myFile;
  Long64_t nentries = fChain->GetEntriesFast();

 
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
 
    //  int isGen = CutLepGenAcc();    
    //  int isRECO =  CutRECO();
    cout << " check "<<endl;

  }
  
  myFile->Write();
  myFile->Close();
  
  return (myFile);
      // if (Cut(ientry) < 0) continue;
}
