#ifndef MuonTrigControl_h
#define MuonTrigControl_h

#include "SUSYLooperHistsSoftBase.h"
#include "SampleInfo.h" 

class MuonTrigControl : public SUSYLooperHistsSoftBase
{
   public:
  MuonTrigControl(SampleInfo mySample, float Mstop_ = 0, float MLSP_ = 0) {setTree(mySample, Mstop_, MLSP_);};
  TFile*   Loop();
};
#endif 
