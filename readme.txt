myplots is an example of making plots from ntupels.


SampleInfo.h is a class that stores all information you need for a sample. It has the path to the rootfile, but also the name of the sample, x-section .... This can generally be used and extended as desired.

SUSYLooperHists loops over a ntupel (TTree) to make plots and gets all the information it needs from SampleInfo for each sample. From any Ttree one can make automated code for looping using makeClass(). SUSYLooperHists is such a code with minor modification and some examples of plots. You can make your oen makeClass code and than apply the modification neede are described below. I.e. if you have a different TTree just add these modifications to the automated code.

myplots uses the above two classes to make loops over some samples (ntupels) and overlays the plots. 

When running the code make sure your root version understand eos and or you use only local samples. The example pathes to rootfiles (.root) need of course to be adopted.


_______________________________________________________________________________________________________
changes to makeClass code:


// default from makeClass():

class markusTreeProducer {
...
  void      Loop(); // return nothing
...

}

markusTreeProducer::markusTreeProducer(TTree *tree) : fChain(0) // need to give it a tree
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("data.root"); 
      if (!f || !f->IsOpen()) {
         f = new TFile("data.root");
      }
      f->GetObject("markusTreeProducer",tree);

   }
   Init(tree);
}


Here the changes for SUSYLooperHists!!

//*** return a TFile pointer for the loop, allows to handle the  TFile withplots afterwards easily***

class SUSYLooperHists {
...
  TFile*      Loop();// returns Tfile with plots
...

}

// *** get all the information like x-section weigt, rootfile path from class SampleInfo ***
SUSYLooperHists::SUSYLooperHists(SampleInfo mySample)// need to give it a SampleInfo
{

  TFile *f = TFile::Open(mySample.FilePath);
  TTree* tree = (TTree*) f->Get("markusTreeProducer"); 
  outputFileName = mySample.OutputFileNameTag+".root";
  weight = mySample.weight();
  Init(tree);
  cout << "looping over: \""<< mySample.FilePath<< "\" to put plots into: \""<< outputFileName<<"\""<<endl;;

}
