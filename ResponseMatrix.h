#ifndef ResponseMatrix_h
#define ResponseMatrix_h

#include "SUSYLooperHistsSoftBase.h"
#include "SampleInfo.h" 

class ResponseMatrix : public SUSYLooperHistsSoftBase
{
   public:
   ResponseMatrix(SampleInfo mySample) {setTree(mySample);};
  TFile*   Loop();
};
#endif 
