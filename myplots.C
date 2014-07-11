#include "SampleInfo.h" 
#include "SUSYLooperHists.C" 
#include "TLegend.h" 
#include <vector>


// just to overlay some plots, should be more automated, e.g. how to check which histo is data, axis ranges, overflows, ratio plots
void plotHist (  vector <TFile*> Files,  vector <SampleInfo> Names, TString Histname)
  {
   TCanvas* canvas = new  TCanvas(Histname+"can",Histname+"can");
   TLegend* alegend = new TLegend(0.6,0.6,0.99,0.99);
   alegend->SetFillColor(0);

   TH1D* histMC = NULL; // histogramm for all MCs added
   int first=0;
   for (unsigned int i =0 ;i< Names.size(); i++)
     {
       
       TH1D* hist = (TH1D*) Files[i]->Get(Histname);
       hist->SetLineColor(Names[i].color);
       cout << Names[i].CrossSection <<endl;
      //  hist->Scale(Names[i].weight());
       alegend->AddEntry(hist, Names[i].LegendLabel,"lp");
       if(i==0) {
	 // can be done more automated, i.e. check CrossSection<0 to use "data" style
	 hist->SetMarkerStyle(20);
	 hist->SetMarkerSize(1.5);
	 hist->Draw("P");
       }   
       else   hist->Draw("HIST same");
       if(Names[i].CrossSection>0) { // for data the x-section is set to -1
	 if(first==0){
	   histMC = ( TH1D*)  hist->Clone();
	   first++;
	 }
	 else{
	   histMC->Add(hist,1);
	 }
       }
     }
   if(histMC==NULL) 
     {
       cout <<"WARNING: The plotHist aim to overlay MC and data, no MC has been givem (no .0 X-section found)"<<endl;
     }
   else{ 
     histMC->SetLineColor(kRed);
     alegend->AddEntry(histMC,"all MC","lp");
     histMC->SetLineWidth(3);
     histMC->Draw("same");
   }
   alegend->Draw("same");
}

void plotHist2D (  vector <TFile*> Files,  vector <SampleInfo> Names, TString Histname)
  {
   TCanvas* canvas = new  TCanvas(Histname+"can",Histname+"can");
   canvas->Divide(2);
  
   TH2D* histMC = NULL; // histogramm for all MCs added
   int first=0;
   for (unsigned int i =0 ;i< Names.size(); i++)
     {
       
       TH2D* hist = (TH2D*) Files[i]->Get(Histname);
       if(Names[i].CrossSection>0) { // for data the x-section is set to -1
	 if(first==0){
	   histMC = ( TH2D*)  hist->Clone();
	   first++;
	 }
	 else{
	   histMC->Add(hist,1);
	 }
       }
       else{ 
	 canvas->cd(1);
	 hist->Draw("colz");
       }
     }
   if(histMC==NULL) 
     {
       cout <<"WARNING: The plotHist aim to overlay MC and data, no MC has been givem (no .0 X-section found)"<<endl;
     }
   else{ 
     canvas->cd(2);
     histMC->Draw("colz");
    
   }
   
}


void myplots()
{

  // initialising the samples to run on 

  // all files are at eos, so start with these
  SampleInfo eosZll(19.7,2655795,40.5,TString("root://eoscms.cern.ch//eos/cms/store/cmst3/group/susy/markus/preliminary_samples_11thApril_2014/DYJetsToLLPtZ100/markusTreeProducer/markusTreeProducer_tree.root"),TString("Z_ll"),TString("Z(ll)"),2);

  // You can copy them over to your local directory to gain speed
  SampleInfo Zll(19.7,2655795,40.5,TString("../data_test_antonios/DYJetsToLLPtZ100/markusTreeProducer/markusTreeProducer_tree.root"),TString("Z(ll)"),TString("Zll"),2);

  // pathes below to local directories 
  SampleInfo wjet(19.7 ,12742383,282,TString("../data_test_antonios/WPT100JETS/markusTreeProducer/markusTreeProducer_tree.root"),TString("W_PT_100"),TString("W(l#nu)"),2);
  SampleInfo ttsemi(19.7,320./963.*24895259,106.15,TString("../data_test_antonios/myhaddedFileTsemi.root"),TString("TTbarSemiLep"),TString("tt(l)"),7);
  SampleInfo ttlep(19.7,11947326*0.467,25.8,TString("../data_test_antonios/myhaddedFileTlep.root"),TString("TTbarFullyLep"),TString("tt(ll)"),8);
  //  SampleInfo Znunu50(19.7,4040980,435.71,TString("../data_test_antonios/ZJetsToNuNu100/markusTreeProducer/markusTreeProducer_tree.root"),TString("Znunu50"),TString("Z(#nu#nu)"),5);
  //  SampleInfo Znunu100(19.7,4416646,186.84,TString("../data_test_antonios/ZJetsToNuNu100/markusTreeProducer/markusTreeProducer_tree.root"),TString("Znunu100"),TString("Z(#nu#nu)"),5);
  //  SampleInfo Znunu200(19.7,5055885,45.60,TString("../data_test_antonios/ZJetsToNuNu100/markusTreeProducer/markusTreeProducer_tree.root"),TString("Znunu200"),TString("Z(#nu#nu)"),5);  SampleInfo Znunu400(19.7,1006928,6.25,TString("../data_test_antonios/ZJetsToNuNu100/markusTreeProducer/markusTreeProducer_tree.root"),TString("Znunu400"),TString("Z(#nu#nu)"),5);
  //  SampleInfo ZZ4mu(19.7,1006928,6.25,TString("/afs/cern.ch/work/m/mstoye/CMGTools_SUSY/CMSSW_5_3_14/src/CMGTools/TTHAnalysis/cfg/Zz4mu/ZZTo4mu/markusTreeProducer/markusTreeProducer_tree.root"),TString("ZZ4m"),TString("Z(#nu#nu)"),5);
  SampleInfo data(19.7,947236  ,-1,TString("data.root"),TString("datapl"),TString("datapl"),1);
  SampleInfo ZZll(19.7,947236  ,0.28,TString("/afs/cern.ch/work/m/mstoye/CMGTools_SUSY/CMSSW_5_3_14/src/CMGTools/TTHAnalysis/cfg/ZZllnunumad/ZZJetsTo2L2Nu/markusTreeProducer/markusTreeProducer_tree.root"),TString("ZZ(ll#nu#nu)"),TString("ZZll"),3);
  SampleInfo WZll(19.7, 2007116 ,0.8674,TString("/afs/cern.ch/work/m/mstoye/CMGTools_SUSY/CMSSW_5_3_14/src/CMGTools/TTHAnalysis/cfg/WZmad/WZJets/markusTreeProducer/markusTreeProducer_tree.root"),TString("WZ(ll)"),TString("WZll"),4);

  SampleInfo WW(19.7, 9976539,54.80,TString("/disk2/WWtree.root"),TString("WW"),TString("WW"),7);
  SampleInfo singleTop(19.7,991118 ,22,TString("/disk2/singleTop.root"),TString("t"),TString("t"),12);
  // signal 

  // for a quick plot a Tfiles are stored in vectors, can be done in many ways
  vector <TFile*> Files;
  vector <SampleInfo> Names;

  //  SUSYLooperHists myZll(eosZll);
  //  TFile* myZllfile = myZll.Loop();
  //  Files.push_back(myZllfile);
  //  Names.push_back(eosZll);
  //  cout << endl<< "see the difference local vs. eos" <<endl <<endl;

  SUSYLooperHists mydata(data);
  TFile* mydataplots = mydata.Loop();
  Files.push_back(mydataplots);
  Names.push_back(data);

  SUSYLooperHists myLocalZll(Zll);
  TFile* myLocalZllfile = myLocalZll.Loop();
  Files.push_back(myLocalZllfile);
  Names.push_back(Zll);

  SUSYLooperHists myWZll(WZll);
  TFile* myWZllfile = myWZll.Loop();
  Files.push_back(myWZllfile);
  Names.push_back(WZll);
  
  SUSYLooperHists myZZll(ZZll);
  TFile* myZZllfile = myZZll.Loop();
  Files.push_back(myZZllfile);
  Names.push_back(ZZll);

  SUSYLooperHists mytt2l(ttlep);
  TFile* myttlep = mytt2l.Loop();
  Files.push_back(myttlep);
  Names.push_back(ttlep);
 
  SUSYLooperHists mytt1l( ttsemi);
  TFile* mytt1lep = mytt1l.Loop();
  Files.push_back(mytt1lep);
  Names.push_back(ttsemi);

  SUSYLooperHists myWs(wjet);
  TFile* myWplots = myWs.Loop();
  Files.push_back(myWplots);
  Names.push_back(wjet);

  SUSYLooperHists myWW(WW);
  TFile* myWWplots = myWW.Loop();
  Files.push_back(myWWplots);
  Names.push_back(WW);

SUSYLooperHists mysingleTop(singleTop);
  TFile* mysingleTopplots = mysingleTop.Loop();
  Files.push_back(mysingleTopplots);
  Names.push_back(singleTop);




  /*
    SampleInfo T23(1.99608,156000,19.7,TString("../data_test_antonios/T2degTest/markusTreeProducer/markusTreeProducer_tree.root"),TString("PrelimT2_out"),TString("PrelimT2_out"),4);
    SUSYLooperHists T2loop3(T23);
    TFile* T2loopFile3 = T2loop3.Loop();
    Files.push_back(T2loopFile3);
    Names.push_back(T23);
  */

  // this is a plot to test the loops, for regular/frequent plotting with overlaying histograms code should be written in a automated manner starting from stored Tfiles
  plotHist(Files,Names,TString("Sig"));
  plotHist(Files,Names,TString("IP3D"));

 plotHist(Files,Names,TString("Dz"));
  plotHist(Files,Names,TString("Dxy"));
plotHist2D(Files,Names,TString("DxyDz"));


}


