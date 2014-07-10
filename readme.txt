myplots is an example of making plots from ntupels.


SampleInfo.h is a class that stores all information you need for a sample. It has the path to the rootfile, but also the name of the sample, x-section .... This can generally be used and extended as desired.

SUSYLooperHists loops over a ntupel (TTree) to make plots and gets all the information it needs from SampleInfo for each sample. From any Ttree one can make automated code for looping using makeClass(). SUSYLooperHists is such a code with minor modification and some examples of plots. You can make your oen makeClass code and than apply the modification neede are described below. I.e. if you have a different TTree just add these modifications to the automated code.

myplots uses the above two classes to make loops over some samples (ntupels) and overlays the plots. 

When running the code make sure your root version understand eos and or you use only local samples. The example pathes to rootfiles (.root) need of course to be adopted.


_______________________________________________________________________________________________________
changes to makeClass code, e.g. XYZ.C 


// defaults  from makeClass() that need to be changed

class XYZ {

 public:
...
  void      Loop(); // return nothing
...

}

XYZ::XYZ(TTree *tree) : fChain(0) // need to give it a tree
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("XYZ.root"); 
      if (!f || !f->IsOpen()) {
         f = new TFile("data.root");
      }
      f->GetObject("XYZ",tree);

   }
   Init(tree);
}


Here the changes for SUSYLooperHists!!

//*** return a TFile pointer for the loop, allows to handle the  TFile withplots afterwards easily***
// add this +++++ to header
#include "SampleInfo.h" 
#include <TString.h>

class SUSYLooperHists {

// add this +++++ (members)
 private:
  TString outputFileName;
  float weight;
//  ++++++++

...
  TFile*      Loop();// returns Tfile with plots: replace void  Loop();
  susySoftLepton(SampleInfo mySample);// constructor that takes sample info as input: replace default contructor susySoftLepton(TTree *tree=0);
...

}

// *** get all the information like x-section weigt, rootfile path from class SampleInfo ***
SUSYLooperHists::SUSYLooperHists(SampleInfo mySample)// need to give it a SampleInfo
{

  TFile *f = TFile::Open(mySample.FilePath);
  TTree* tree = (TTree*) f->Get(YOURTREENAME); 
  outputFileName = mySample.OutputFileNameTag+".root";
  weight = mySample.weight();
  Init(tree);
  cout << "looping over: \""<< mySample.FilePath<< "\" to put plots into: \""<< outputFileName<<"\""<<endl;;

}
