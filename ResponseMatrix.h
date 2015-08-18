#include "A7LoopBase.h"
#include "SampleInfo.h" 

#ifndef ResponseMatrix_h
#define ResponseMatrix_h



class ResponseMatrix : public A7LoopBase
{
   public:
   ResponseMatrix(SampleInfo mySample) {setTree(mySample);};
  TFile*   Loop();
};
#endif 
