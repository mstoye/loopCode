#include "SampleInfo.h" 
#include "SUSYLooperHists.C" 
#include <vector>

void myplots()
{


  // 
  SampleInfo wjet(282,12742383,19.7 ,TString("../data_test_antonios/WPT100JETS/markusTreeProducer/markusTreeProducer_tree.root"),TString("W_PT_100"),TString("W PT 100"),2);
  SampleInfo ttsemi(106.15,320./963.*24895259,19.7,TString("../data_test_antonios/myhaddedFileTsemi.root"),TString("TTbarSemiLep"),TString("tt(l)"),4);
  SampleInfo ttlep(25.8,11947326*0.467,19.7,TString("../data_test_antonios/myhaddedFileTlep.root"),TString("TTbarFullyLep"),TString("tt(l)"),5);
  SampleInfo Zll(40.5,2655795,19.7,TString("../data_test_antonios/DYJetsToLLPtZ100/markusTreeProducer/markusTreeProducer_tree.root"),TString("Zll"),TString("Z(ll)"),5);
  SampleInfo Znunu50(435.71,4040980,19.7,TString("../data_test_antonios/ZJetsToNuNu100/markusTreeProducer/markusTreeProducer_tree.root"),TString("Znunu50"),TString("Z(#nu#nu)"),5);
  SampleInfo Znunu100(186.84,4416646,19.7,TString("../data_test_antonios/ZJetsToNuNu100/markusTreeProducer/markusTreeProducer_tree.root"),TString("Znunu100"),TString("Z(#nu#nu)"),5);
  SampleInfo Znunu200(45.60,5055885,19.7,TString("../data_test_antonios/ZJetsToNuNu100/markusTreeProducer/markusTreeProducer_tree.root"),TString("Znunu200"),TString("Z(#nu#nu)"),5);
  SampleInfo Znunu400(6.25,1006928,19.7,TString("../data_test_antonios/ZJetsToNuNu100/markusTreeProducer/markusTreeProducer_tree.root"),TString("Znunu400"),TString("Z(#nu#nu)"),5);
  SampleInfo ZZ4mu(6.25,1006928,19.7,TString("/afs/cern.ch/work/m/mstoye/CMGTools_SUSY/CMSSW_5_3_14/src/CMGTools/TTHAnalysis/cfg/Zz4mu/ZZTo4mu/markusTreeProducer/markusTreeProducer_tree.root"),TString("ZZ4m"),TString("Z(#nu#nu)"),5);
  SampleInfo ZZ(6.25,1006928,19.7,TString("ZZpythia.root"),TString("ZZ4m"),TString("Z(#nu#nu)"),5);


  vector <TFile*> Files;
  vector <SampleInfo> Names;
  SUSYLooperHists myZll(Zll);
  TFile* myZllfile = myZll.Loop();
  Files.push_back(myZll);
  Names.push_back(Zll);
  SUSYLooperHists mytt2l(ttlep);
  TFile* myttlep = mytt2l.Loop();
  Files.push_back(myttlep);
  Names.push_back(ttlep);
  SUSYLooperHists mytt1l(  ttsemi);
  TFile* mytt1lep = mytt1l.Loop();
  Files.push_back(mytt1lep);
  Names.push_back(ttsemi);
  
  /*
   SUSYLooperHists myZnunu50(Znunu50);
   TFile* myZnunu50file = myZnunu50.Loop();
   Files.push_back(myZnunu50file);
   Names.push_back(Znunu50);
   SUSYLooperHists myZnunu100(Znunu100);
   TFile* myZnunu100file = myZnunu100.Loop();
   Files.push_back(myZnunu100file);
   Names.push_back(Znunu100);
   SUSYLooperHists myZnunu200(Znunu200);
   TFile* myZnunu200file = myZnunu200.Loop();
   Files.push_back(myZnunu200file);
   Names.push_back(Znunu200);
   SUSYLooperHists myZnunu400(Znunu400);
   TFile* myZnunu400file = myZnunu400.Loop();
   Files.push_back(myZnunu400file);
   Names.push_back(Znunu400);
 
   SUSYLooperHists myZll(Zll);
   TFile* myZllfile = myZll.Loop();
   Files.push_back(myZllfile);
   Names.push_back(Zll);
  

   SampleInfo T23(1.99608,156000,TString("../data_test_antonios/T2degTest/markusTreeProducer/markusTreeProducer_tree.root"),TString("PrelimT2_out"),TString("PrelimT2_out"),4);
    SUSYLooperHists T2loop3(T23);
    TFile* T2loopFile3 = T2loop3.Loop();
    Files.push_back(T2loopFile3);
    Names.push_back(T23);
  */
    /*
   SUSYLooperHists myZll(Zll);
   TFile* myZllfile = myZll.Loop();
   Files.push_back(myZllfile);
   Names.push_back(Zll);
   
    SUSYLooperHists myWs(wjet);
     TFile* myWplots = myWs.Loop();
    Files.push_back(myWplots);
    Names.push_back(wjet);
   */
  

    /*
    
  
    
   TCanvas* canvas = new  TCanvas("LepTwoZmass","LepTwoZmass");
   


   for (int i =0 ;i< Names.size(); i++)
     {
       
      TH1D* hist = (TH1D*) Files[i]->Get("Iso1D");
      hist->SetLineColor(Names[i].SampleColor);
      
      //  hist->Scale(Names[i].weight());
      if(i==0) hist->Draw();
      else   hist->Draw("same");
    }
  
  TCanvas* canvaos = new  TCanvas("LepPtW","LeoptW");
  
  for (int i =0 ;i< Names.size(); i++)
    {
      
      TH1D* hist = (TH1D*) Files[i]->Get("MHTHT");
      hist->SetLineColor(Names[i].SampleColor);
      
      //   hist->Scale(Names[i].weight());
      if(i==0) hist->Draw();
      else   hist->Draw("same");
    }
  TCanvas* can76vaos = new  TCanvas("LetW","LtW");
  
  for (int i =0 ;i< Names.size(); i++)
    {
      
      TH1D* hist = (TH1D*) Files[i]->Get("LepPtTwo");
      hist->SetLineColor(Names[i].SampleColor);
      
      //    hist->Scale(Names[i].weight());
      if(i==0) hist->Draw();
      else   hist->Draw("same");
    }


 TCanvas* SBcan = new  TCanvas("SB","SB");


 TH2D* hist = (TH2D*) Files[0]->Get("STLepM");
 hist->Scale(Names[0].weight());
 TH2D* histStop = (TH2D*) Files[1]->Get("STLepM");
 histStop->Scale(Names[1].weight());
 SBcan->Divide(2);
 SBcan->cd(1);
 histStop->DrawCopy("colz");
 histStop->Divide(hist);
 SBcan->cd(2);
 histStop->DrawCopy("colz");
   

    */
}

