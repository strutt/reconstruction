// -*- C++ -*-.
/***********************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description:
             Reconstruct decimated data set.
********************************************************************************************************* */

#include "TFile.h"
#include "TChain.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TProfile2D.h"
#include "THnSparse.h"

#include "RawAnitaHeader.h"
#include "UsefulAdu5Pat.h"
#include "UsefulAnitaEvent.h"
#include "CalibratedAnitaEvent.h"
#include "AnitaEventCalibrator.h"

#include "ProgressBar.h"
#include "CrossCorrelator.h"
#include "OutputConvention.h"
#include "AnitaEventSummary.h"
#include "FFTtools.h"


int main(int argc, char *argv[]){

  // 135 06:39:00 23940
  // 136 03:43:11 13391
  // 137 07:15:59 26159
  // 138 08:32:09 30729
  // 139 07:56:04 28564
  // 140 08:22:25 30145
  // 145 00:56:12 3372
  // 150 07:11:23 25883
  // 151 07:35:38 27338
  // 153 04:13:13 15193
  // 155 06:28:13 23293
  // 156 03:31:45 12705
  // 161 02:06:50 7610
  // 162 06:41:55 24115
  // 180 04:51:36 17496
  // 181 02:03:16 7396
  // 182 00:30:41 1841
  // 183 01:06:33 3993
  // 190 09:45:35 35135
  // 197 04:30:17 16217
  // 198 00:46:34 2794
  // 199 06:47:32 24452
  // 205 01:05:20 3920
  // 208 07:08:48 25728

  std::map<Int_t, Double_t> runTimes;
  runTimes[135] = 23940;
  runTimes[136] = 13391;
  runTimes[137] = 26159;
  runTimes[138] = 30729;
  runTimes[139] = 28564;
  runTimes[140] = 30145;
  runTimes[145] = 3372;
  runTimes[150] = 25883;
  runTimes[151] = 27338;
  runTimes[153] = 15193;
  runTimes[155] = 23293;
  runTimes[156] = 12705;
  runTimes[161] = 7610;
  runTimes[162] = 24115;
  runTimes[180] = 17496;
  runTimes[181] = 7396;
  runTimes[182] = 1841;
  runTimes[183] = 3993;
  runTimes[190] = 35135;
  runTimes[197] = 16217;
  runTimes[198] = 2794;
  runTimes[199] = 24452;
  runTimes[205] = 3920;
  runTimes[208] = 25728;
  std::map<Int_t ,Double_t>::iterator theEnd = runTimes.end();
  
  if(!(argc==3 || argc==2)){
    std::cerr << "Usage 1: " << argv[0] << " [run]" << std::endl;
    std::cerr << "Usage 2: " << argv[0] << " [firstRun] [lastRun]" << std::endl;
    return 1;
  }
  
  std::cout << argv[0] << "\t" << argv[1];
  if(argc==3){std::cout << "\t" << argv[2];}
  std::cout << std::endl;
  const Int_t firstRun = atoi(argv[1]);
  const Int_t lastRun = argc==3 ? atoi(argv[2]) : firstRun;

  Double_t maxTime = 0;
  Long64_t totalEvents  = 0;
  Long64_t doneEvents  = 0;

  TGraph* gr = new TGraph();  

  std::map<Int_t, Long64_t> numEntriesRun;
  
  for(Int_t run=firstRun; run<=lastRun; run++){
    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    TFile* headFile = TFile::Open(fileName);
    if(headFile){
      TTree* headTree = (TTree*) headFile->Get("headTree");
      const Long64_t numEntries = headTree->GetEntries();
      std::cout << run << "\t" << numEntries << std::endl;
      headFile->Close();

      totalEvents += numEntries;

      std::map<Int_t ,Double_t>::iterator it = runTimes.find(run);
      if(it!=theEnd){
	// Double_t timePerEvent = it->second / numEntries;
	maxTime += it->second;
	doneEvents += numEntries;
	gr->SetPoint(gr->GetN(), numEntries, it->second/3600);
	std::cerr << "numEntries = " << numEntries << "\thours = " << it->second/3600 << std::endl;
      }
      else{
	numEntriesRun[run] = numEntries;
      }
    }
  }  

  Double_t meanRate = doneEvents / (maxTime/3600);

  OutputConvention oc(argc, argv);
  TString outFileName = oc.getOutputFileName();
  TFile* outFile = new TFile(outFileName, "recreate");
  if(outFile->IsZombie()){
    std::cerr << "Error! Unable to open output file " << outFileName.Data() << std::endl;
    return 1;
  }

  
  
  TF1* fLine = new TF1("fLine", "[0]*x+[1]");
  gr->Fit(fLine);
  gr->SetName("grRunTimes");
  gr->SetTitle("Run times vs. Num Entries; Number of entries; Run time (hours)");  
  // std::cout << "I reckon you'll need " << (maxTime / doneEvents)* (totalEvents - doneEvents)/(3600) << " hours of CPU time" << std::endl;
  gr->Write();

  TGraph* gr2 = new TGraph();
  Double_t remainingTimeEst = 0;
  for(auto& it : numEntriesRun){
    // Double_t timeEst = it.second*fLine->GetParameter(0) + it.second*fLine->GetParameter(1);
    Double_t timeEst = it.second*meanRate;
    gr2->SetPoint(gr2->GetN(), it.second, timeEst);
    remainingTimeEst += timeEst;
  }

  gr2->SetName("grEstimate");
  gr2->Write();

  std::cout << "I reckon there's " << remainingTimeEst << " remaining" << std::endl;

  outFile->Write();
  outFile->Close();  
  
  return 0;
}






