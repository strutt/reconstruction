// -*- C++ -*-.
/***********************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description:
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
#include "TGraph2D.h"
#include "TRandom3.h"

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
#include "AntarcticaMapPlotter.h"

#include "AnitaClusterer.h"
#include "AnalysisCuts.h"


int main(int argc, char *argv[]){

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

  TChain* eventSummaryChain = new TChain("cutTree");
  // TChain* headChain = new TChain("headTree");
  // TChain* indexedHeadChain = new TChain("headTree");
  // TChain* gpsChain = new TChain("adu5PatTree");

  for(Int_t run=firstRun; run<=lastRun; run++){
    if(run>=257 && run<=263){
      continue;
    }
    if(run == 198 || run == 287){
      continue;
    }

    TString fileName = TString::Format("projectionPlots/projectionPlots_%d_*", run);
    eventSummaryChain->Add(fileName);

  }

  if(eventSummaryChain->GetEntries()==0){
    std::cerr << "Unable to find eventSummary files!" << std::endl;
    return 1;
  }


  TFile* fRes = TFile::Open("resolutionPlots.root");
  if(fRes==NULL){
    std::cerr << "Warning! Unable to find resolution plots!" << std::endl;
  }
  TH1D* hResTheta[AnitaPol::kNotAPol];
  TH1D* hResPhi[AnitaPol::kNotAPol];
  hResTheta[AnitaPol::kHorizontal] = (TH1D*) fRes->Get("hDeltaThetaWais2_2__4__4");
  hResPhi[AnitaPol::kHorizontal] = (TH1D*) fRes->Get("hDeltaPhiWais2_2__3__3");
  hResTheta[AnitaPol::kVertical] = (TH1D*) fRes->Get("hDeltaThetaLdb2_2__2__2");
  hResPhi[AnitaPol::kVertical] = (TH1D*) fRes->Get("hDeltaPhiLdb2_2__1__1");


  TF1* fThetaFit = (TF1*) (hResTheta[AnitaPol::kHorizontal]->GetListOfFunctions()->At(0));
  TF1* fPhiFit = (TF1*) (hResPhi[AnitaPol::kVertical]->GetListOfFunctions()->At(0));


  AnitaEventSummary* eventSummary = NULL;
  eventSummaryChain->SetBranchAddress("eventSummary", &eventSummary);
  AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal;
  eventSummaryChain->SetBranchAddress("pol", &pol);
  int sunCut=0;
  eventSummaryChain->SetBranchAddress("sunCut", &sunCut);



  // AnalysisCuts::Status_t l3TriggerCut;
  // eventSummaryChain->SetBranchAddress("l3TriggerCut", (int*)&l3TriggerCut);
  // AnalysisCuts::Status_t surfSaturation;
  // eventSummaryChain->SetBranchAddress("surfSaturation", (int*)&surfSaturation);
  // AnalysisCuts::Status_t selfTriggeredBlastCut;
  // eventSummaryChain->SetBranchAddress("selfTriggeredBlastCut", (int*)&selfTriggeredBlastCut);
  // AnalysisCuts::Status_t thermalCut;
  // eventSummaryChain->SetBranchAddress("thermalCut", (int*)&thermalCut);
  // AnalysisCuts::Status_t thetaAngleCut;
  // eventSummaryChain->SetBranchAddress("thetaAngleCut", (int*)&thetaAngleCut);
  // AnalysisCuts::Status_t peakRatioCut;
  // eventSummaryChain->SetBranchAddress("peakRatioCut", (int*)&peakRatioCut);


  Adu5Pat* pat = NULL;
  // gpsChain->SetBranchAddress("pat", &pat);
  eventSummaryChain->SetBranchAddress("pat", &pat);
  // UInt_t eventNumberPat = 0;
  // gpsChain->SetBranchAddress("eventNumber", &eventNumberPat);

  OutputConvention oc(argc, argv);
  TString outFileName = oc.getOutputFileName();
  TFile* outFile = new TFile(outFileName, "recreate");
  if(outFile->IsZombie()){
    std::cerr << "Error! Unable to open output file " << outFileName.Data() << std::endl;
    return 1;
  }

  Long64_t nEntries = eventSummaryChain->GetEntries();
  Long64_t maxEntry = 0; //2500;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  const int K = 50;
  const int numIterations = 10;
  AnitaClusterer clusterer(K, numIterations, nEntries);

  for(Long64_t entry = startEntry; entry < maxEntry; entry++){
    eventSummaryChain->GetEntry(entry);

    // if(l3TriggerCut==0 && surfSaturation==0 && selfTriggeredBlastCut==0 && thermalCut==0 && thetaAngleCut == 0 && peakRatioCut==0){

      // Long64_t entry2 = indexedHeadChain->GetEntryNumberWithIndex(eventSummary->eventNumber);
      // gpsChain->GetEntry(entry);

      // if(eventSummary->eventNumber!=eventNumberPat){
      //   std::cerr << "fuck" << std::endl;
      // }
    // if(sunCut!=0){
    //   continue;
    // }
    const int peakInd = 0;

    Double_t sourceLat = eventSummary->peak[pol][peakInd].latitude;
    Double_t sourceLon = eventSummary->peak[pol][peakInd].longitude;
    Double_t sourceAlt = eventSummary->peak[pol][peakInd].altitude;
    if(sourceLat > -999 && sourceLon > -999){
      if(eventSummary->peak[pol][peakInd].distanceToSource > 1 && eventSummary->peak[pol][peakInd].distanceToSource < 1e6){

	// std::cout << sourceLat << "\t" << sourceLon << "\t" << sourceAlt << std::endl;

	// Double_t snr = eventSummary->peak[pol][peakInd].snr;
	// // Double_t snr = eventSummary->coherent[pol][peakInd].snr;
	// Double_t sigmaTheta = fThetaFit->Eval(snr);
	// Double_t sigmaPhi = fPhiFit->Eval(snr);
	Double_t sigmaTheta = fThetaFit->Eval(0);
	Double_t sigmaPhi = fPhiFit->Eval(0);


	// std::cout << snr << "\t" << sigmaTheta << "\t" << sigmaPhi << std::endl;

	clusterer.addPoint(pat, sourceLat,sourceLon,sourceAlt, eventSummary->run, eventSummary->eventNumber, sigmaTheta, sigmaPhi, pol);
      }
    }
    // }

    p++;
  }

  clusterer.initializeBaseList();


  // clusterer.mergeClusters();

  Int_t sumNumEntries2 = 0;
  const int minEnergy= 19;
  const int maxEnergy= 21;
  const int numHalves = 2;
  // TChain* mcEventSummaryChain = new TChain("cutTree");
  // TChain* mcEventSummaryChain = new TChain("eventSummaryTree");
  for(int energyInd=minEnergy; energyInd <= maxEnergy; energyInd++){
  // for(int energyInd=18; energyInd <= 20; energyInd++){
    for(int halfInd=numHalves-1; halfInd >= 0; halfInd--){
      // for(int mcRun=1; mcRun < 30; mcRun++){
      double energy = double(energyInd) - 0.5*halfInd;

      TString fileName = TString::Format("projectionMonteCarloPlots_%d_%d_*.root", energyInd, halfInd);
      TFile* f = OutputConvention::getFile(fileName);
      //TFile::Open(fileName);
      if(f==nullptr){
	continue;
      }
      TTree* mcEventSummaryChain = (TTree*) f->Get("cutTree");

      const Long64_t numEntries2 = mcEventSummaryChain->GetEntries();
      sumNumEntries2 += numEntries2;
      std::cerr << fileName << "\tnumEntries2 = " << numEntries2 << std::endl;
      ProgressBar prog2(numEntries2);
      AnitaEventSummary* eventSummary2 = NULL;
      Adu5Pat* pat2 = NULL;
      Double_t weight = 0;
      mcEventSummaryChain->SetBranchAddress("weight", &weight);
      mcEventSummaryChain->SetBranchAddress("pat", &pat2);
      mcEventSummaryChain->SetBranchAddress("eventSummary", &eventSummary2);
      for(Long64_t entry2=0; entry2 < mcEventSummaryChain->GetEntries(); entry2++){
	mcEventSummaryChain->GetEntry(entry2);

	for(Int_t polInd=0; polInd < NUM_POL; polInd++){
	  for(int peakInd=0; peakInd < 5; peakInd++){
	    Double_t sourceLat = eventSummary2->peak[polInd][peakInd].latitude;
	    Double_t sourceLon = eventSummary2->peak[polInd][peakInd].longitude;
	    Double_t sourceAlt = eventSummary2->peak[polInd][peakInd].altitude;
	    if(sourceLat > -999 && sourceLon > -999){
	      // std::cout << (eventSummary->peak[polInd][peakInd].distanceToSource < 1e6 )<< std::endl;
	      // if(eventSummary->peak[polInd][peakInd].theta < 0 && eventSummary->peak[polInd][peakInd].distanceToSource < 1e6){
	      if(eventSummary2->peak[polInd][peakInd].distanceToSource < 1e6){

		// Double_t snr = eventSummary2->peak[polInd][peakInd].snr;

		// Double_t snr = eventSummary2->coherent[polInd][peakInd].snr;
		// Double_t sigmaTheta = hResTheta[polInd]->GetBinContent(hResPhi[polInd]->FindBin(snr));
		// Double_t sigmaPhi = hResPhi[polInd]->GetBinContent(hResPhi[polInd]->FindBin(snr));
		Double_t sigmaTheta = fThetaFit->Eval(0);
		Double_t sigmaPhi = fPhiFit->Eval(0);
		// Double_t sigmaTheta = fThetaFit->Eval(snr);
		// Double_t sigmaPhi = fPhiFit->Eval(snr);


		// std::cout << polInd << "\t" << snr << "\t" << sigmaTheta << "\t" << sigmaPhi << std::endl;
		// std::cout << sourceLat << "\t" << sourceLon << "\t" << sourceAlt << std::endl;

		clusterer.addMCPoint(pat2, sourceLat, sourceLon, sourceAlt,
				     eventSummary2->run, eventSummary2->eventNumber,
				     sigmaTheta, sigmaPhi, (AnitaPol::AnitaPol_t)polInd, weight, (double)energy);
		// std::cerr << pat2->altitude << "\t" << pat2->latitude << "\t" << pat->longitude << std::endl;
		// Int_t n = clusterer.addPoint(sourceLat,sourceLon,sourceAlt);
		// std::cout << n << std::endl;
	      }
	    }
	  }
	}
	prog2++;
      }
      f->Close();
    }
  }
  outFile->cd();
  if(sumNumEntries2==0){
    std::cerr << "Error: read in " << sumNumEntries2 << " MC events" << std::endl;
    return 1;
  }


  // TGraph* grNumPoints = new TGraph();
  // TGraph* grNumMcSinglets = new TGraph();
  std::vector<TGraph*> grNumIsolateds(clusterer.maxRetestClusterSize, NULL);
  std::vector<TGraph*> grNumBaseIsolateds(clusterer.maxRetestClusterSize, NULL);
  for(int i=0; i < clusterer.maxRetestClusterSize; i++){
    grNumIsolateds.at(i) = new TGraph();
    grNumBaseIsolateds.at(i) = new TGraph();
  }

  const int numLLs = 1;
  // const double deltaLL = 1;
  double llCuts[numLLs] = {100};

  // // double llCuts[numLLs] = {0};//, 10, 20, 40, 80, 160, 320, 640, 1280, 2560};
  // for(int j=0; j < numLLs; j++){
  //   llCuts[j] = (j+1)*deltaLL;
  // }

  // Double_t sumMcWeights = clusterer.getSumOfMcWeights();
  for(int i=0; i < numLLs; i++){
    clusterer.llCut = llCuts[i];
    clusterer.resetClusters();
    clusterer.recursivelyAddClusters(0);

    std::cout << "********************************" << std::endl;
    std::cout << "llCut = " << llCuts[i] << std::endl;
    std::cout << "********************************" << std::endl;
    clusterer.findClosestPointToClustersOfSizeOne();

    clusterer.assignMCPointsToClusters();

    for(int j=0; j < clusterer.maxRetestClusterSize; j++){
      grNumIsolateds.at(j)->SetPoint(grNumIsolateds.at(j)->GetN(),
				     clusterer.llCut, clusterer.numIsolatedSmallClusters.at(j));
      grNumBaseIsolateds.at(j)->SetPoint(grNumBaseIsolateds.at(j)->GetN(),
					 clusterer.llCut, clusterer.numIsolatedSmallBaseClusters.at(j));
    }
    // grNumMcSinglets->SetPoint(grNumMcSinglets->GetN(), llCuts[i], clusterer.numMcIsolatedSinglets/sumMcWeights);
  }

  for(int i=0; i < clusterer.maxRetestClusterSize; i++){
    TString grName = TString::Format("grNumBaseIsolateds_%d", i+1);
    grNumBaseIsolateds.at(i)->SetName(grName);
    grNumBaseIsolateds.at(i)->Write();

    if(i > 0){
      grName = TString::Format("grNumIsolateds_%d", i+1);
      grNumIsolateds.at(i)->SetName(grName);
      grNumIsolateds.at(i)->Write();
    }
  }
  // grNumMcSinglets->SetName("grNumMcSinglets");
  // grNumMcSinglets->Write();

  outFile->cd();

  // for(int clusterInd=0; clusterInd < clusterer.getNumClusters(); clusterInd++){
  //   TGraph* gr = clusterer.makeClusterSummaryTGraph(clusterInd);

  //   if(gr){
  //     gr->Write();
  //     delete gr;
  //   }
  //   else{
  //     std::cerr << "??????? " << clusterInd << std::endl;
  //   }
  // }

  outFileName.ReplaceAll("Plots", "HiddenPlots");
  TFile* signalBox = new TFile(outFileName, "recreate");
  // no need to write?
  TTree* clusterTree = clusterer.makeClusterSummaryTree(outFile, signalBox);
  // clusterTree->BuildIndex("eventNumber");
  std::cout << "Made clusterTree with " << clusterTree->GetEntries() << " entries." << std::endl;

  ClusteredAnitaEvent* clusteredEvent = 0;
  clusterTree->SetBranchAddress("clusteredEvent", &clusteredEvent);
  // eventSummaryChain->BuildIndex("eventNumber");



  for(Long64_t entry = startEntry; entry < maxEntry; entry++){
    clusterTree->GetEntry(entry);

    if(clusteredEvent->inCluster < 0){
      eventSummaryChain->GetEntryWithIndex(clusteredEvent->eventNumber);
      if(clusteredEvent->eventNumber == eventSummary->eventNumber){

	std::cout << clusteredEvent->eventNumber  << "\t" << eventSummary->eventNumber << std::endl;
      }
    }
  }

  outFile->Write();
  outFile->Close();

  return 0;
}
