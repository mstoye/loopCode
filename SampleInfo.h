#ifndef SampleInfo_h
#define SampleInfo_h

#include "TString.h"

class SampleInfo
{

 public:
  SampleInfo( float Lumi_,float NEvents_,float CrossSection_,TString FilePath_, TString LegendLabel_,TString OutputFileNameTag_,int color_): Lumi(Lumi_),NEvents(NEvents_), CrossSection(CrossSection_),FilePath(FilePath_),LegendLabel(LegendLabel_),OutputFileNameTag(OutputFileNameTag_), color(color_))
{;}
    ;
  float weight();


 private:
  float Lumi;
  float NEvents; // number of events before any skimmming in MC
  float CrossSection; // put -1 for data
  TString FilePath; 
  TString LegendLabel;
  TString OutputFileNameTag;
  int color;
 

};


float SampleInfo::weight()
{

  if(CrossSection==-2) return 1.; // for scans
  if(CrossSection==-1) return 1.; // for data
  else
    {
      return (CrossSection/NEvents*Lumi*1000);
    }
};
#endif 
