// -*- C++ -*-.
/***********************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description:
             Reconstruct entire data set.
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

    const Int_t nDivisions = 10; // divide job up into this many subjobs

  // if(!(argc==3 || argc==2)){
  if(argc!=3){
    std::cerr << "Usage: " << argv[0] << TString::Format(" [firstRun] [subDivisionOfEvents (0-%d)]", nDivisions)
	      << std::endl;
    return 1;
  }

  std::cout << argv[0] << "\t" << argv[1];
  if(argc==3){std::cout << "\t" << argv[2];}
  std::cout << std::endl;
  const Int_t firstRun = atoi(argv[1]);
  const Int_t lastRun = firstRun; //argc==3 ? atoi(argv[2]) : firstRun;

  const Int_t whichDivisionOfEvents = atoi(argv[2]); // should be a number from 0-9 inclusive
  if(whichDivisionOfEvents >= nDivisions){
    return 1;
  }

  CrossCorrelator* cc = new CrossCorrelator();

  TChain* headChain = new TChain("headTree");
  // TChain* indexedHeadChain = new TChain("headTree");
  TChain* gpsChain = new TChain("adu5PatTree");
  TChain* calEventChain = new TChain("eventTree");

  for(Int_t run=firstRun; run<=lastRun; run++){
    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    // TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/decimatedHeadFile%d.root", run, run);
    headChain->Add(fileName);

    // fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/indexedBlindHeadFile%d.root", run, run);
    // indexedHeadChain->Add(fileName);

    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsEvent%d.root", run, run);
    gpsChain->Add(fileName);

    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/calEventFile%d.root", run, run);
    calEventChain->Add(fileName);
  }

  if(headChain->GetEntries()==0){
    std::cerr << "Unable to find header files!" << std::endl;
    return 1;
  }
  if(gpsChain->GetEntries()==0){
    std::cerr << "Unable to find gps files!" << std::endl;
    return 1;
  }
  if(calEventChain->GetEntries()==0){
    std::cerr << "Unable to find calEvent files!" << std::endl;
    return 1;
  }

  // indexedHeadChain->BuildIndex("eventNumber");
  // calEventChain->BuildIndex("eventNumber");
  // gpsChain->BuildIndex("eventNumber");

  RawAnitaHeader* header = NULL;
  headChain->SetBranchAddress("header", &header);
  Adu5Pat* pat = NULL;
  gpsChain->SetBranchAddress("pat", &pat);
  CalibratedAnitaEvent* calEvent = NULL;
  calEventChain->SetBranchAddress("event", &calEvent);

  OutputConvention oc(argc, argv);
  TString outFileName = oc.getOutputFileName();
  TFile* outFile = new TFile(outFileName, "recreate");
  if(outFile->IsZombie()){
    std::cerr << "Error! Unable to open output file " << outFileName.Data() << std::endl;
    return 1;
  }

  CrossCorrelator::SimpleNotch notch260("n260Notch", "260MHz Satellite And 200MHz Notch Notch",
  					260-26, 260+26);
  CrossCorrelator::SimpleNotch notch370("n370Notch", "370MHz Satellite Notch",
					370-26, 370+26);
  CrossCorrelator::SimpleNotch notch400("n400Notch", "400 MHz Satellite Notch",
					 400-10, 410);
  CrossCorrelator::SimpleNotch notch762("n762Notch", "762MHz Satellite Notch (one bin wide)",
					762-8, 762+8);
  CrossCorrelator::SimpleNotch notch200("n200Notch", "200 MHz high pass band",
					0, 200);
  CrossCorrelator::SimpleNotch notch1200("n1200Notch", "1200 MHz low pass band",
					 1200, 9999);

  cc->addNotch(notch260);
  cc->addNotch(notch370);
  cc->addNotch(notch400);
  cc->addNotch(notch762);
  cc->addNotch(notch200);
  cc->addNotch(notch1200);

  const Int_t myNumPeaksCoarse = 2; //5;
  const Int_t myNumPeaksFine = 2; //5;
  const Int_t coherentDeltaPhi = 0;

  TNamed* comments = new TNamed("comments", "Applied simple, static notch at 260#pm26 MHz, 370#pm26 MHz, 762#pm8 MHz");
  comments->Write();
  delete comments;

  TNamed* comments2 = new TNamed("comments2",
				 TString::Format("%d coarse peaks, %d fine peaks",
						 myNumPeaksCoarse, myNumPeaksFine).Data());
  comments2->Write();
  delete comments2;

  TTree* eventSummaryTree = new TTree("eventSummaryTree", "eventSummaryTree");
  // AnitaEventSummary* eventSummary = new AnitaEventSummary();
  AnitaEventSummary* eventSummary = NULL; //new AnitaEventSummary();
  eventSummaryTree->Branch("eventSummary", &eventSummary);

  Long64_t nEntries = headChain->GetEntries();
  Long64_t startEntry = whichDivisionOfEvents*floor(Double_t(nEntries)/nDivisions);
  Long64_t maxEntry = (whichDivisionOfEvents+1)*floor(Double_t(nEntries)/nDivisions);
  maxEntry = whichDivisionOfEvents==nDivisions-1 ? nEntries : maxEntry;

  std::cout << "Processing " << maxEntry-startEntry << " of " << nEntries << " entries." << std::endl;
  std::cout << "Starting at entry " << startEntry << " up to  " << maxEntry << " entry." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  for(Long64_t entry = startEntry; entry < maxEntry; entry++){

    headChain->GetEntry(entry);

    // Int_t isMinBias = RootTools::isMinBiasSampleEvent(header);
    // if(isMinBias == 0){ // don't want to do min bias events twice
    {

      // Int_t entry2 = indexedHeadChain->GetEntryNumberWithIndex(header->eventNumber);
      // gpsChain->GetEntry(entry2);
      // calEventChain->GetEntry(entry2);

      gpsChain->GetEntry(entry);
      calEventChain->GetEntry(entry);

      UsefulAnitaEvent* usefulEvent = new UsefulAnitaEvent(calEvent);

      UsefulAdu5Pat usefulPat(pat);

      cc->reconstructEvent(usefulEvent, myNumPeaksCoarse, myNumPeaksFine);
      // AnitaPol::AnitaPol_t peakPol = cc->reconstructEventPeakPol(usefulEvent, myNumPeaksCoarse, myNumPeaksFine);


      eventSummary = new AnitaEventSummary(header, &usefulPat);
      // std::cout << eventSummary->sun.theta << "\t" << eventSummary->sun.phi << std::endl;

      Double_t minY = 0;


      for(Int_t polInd=0; polInd < AnitaPol::kNotAPol; polInd++){

	AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;

	// if(pol!=peakPol){
	for(Int_t peakInd=0; peakInd < myNumPeaksCoarse; peakInd++){

	  const int thisPeakInd = peakInd+myNumPeaksFine;

	  cc->getCoarsePeakInfo(pol, peakInd,
				eventSummary->peak[pol][thisPeakInd].value,
				eventSummary->peak[pol][thisPeakInd].phi,
				eventSummary->peak[pol][thisPeakInd].theta);

	  TGraph* grGlobal0 = cc->makeCoherentlySummedWaveform(pol,
							       eventSummary->peak[pol][thisPeakInd].phi,
							       eventSummary->peak[pol][thisPeakInd].theta,
							       coherentDeltaPhi,
							       eventSummary->peak[pol][thisPeakInd].snr);


	  if(grGlobal0!=NULL){

	    TGraph* grGlobal0Hilbert = FFTtools::getHilbertEnvelope(grGlobal0);

	    RootTools::getMaxMin(grGlobal0Hilbert, eventSummary->coherent[pol][thisPeakInd].peakHilbert, minY);

	    delete grGlobal0;
	    delete grGlobal0Hilbert;
	  }
	}
	// }
	// else{
	for(Int_t peakInd=0; peakInd < myNumPeaksFine; peakInd++){
	  cc->getFinePeakInfo(pol, peakInd,
			      eventSummary->peak[pol][peakInd].value,
			      eventSummary->peak[pol][peakInd].phi,
			      eventSummary->peak[pol][peakInd].theta);

	  // std::cout << eventSummary->peak[pol][peakInd].value << "\t"
	  // 	    << eventSummary->peak[pol][peakInd].phi << "\t"
	  // 	    << eventSummary->peak[pol][peakInd].theta << std::endl;


	  TGraph* grZ0 = cc->makeUpsampledCoherentlySummedWaveform(pol,
								   eventSummary->peak[pol][peakInd].phi,
								   eventSummary->peak[pol][peakInd].theta,
								   coherentDeltaPhi,
								   // eventSummary->peak[pol][peakInd].snr);
								   eventSummary->coherent[pol][peakInd].snr);

	  // std::cout << grZ0 << "\t" << eventSummary->peak[pol][peakInd].snr << std::endl;


	  if(grZ0!=NULL){
	    TGraph* grZ0Hilbert = FFTtools::getHilbertEnvelope(grZ0);

	    RootTools::getMaxMin(grZ0Hilbert, eventSummary->coherent[pol][peakInd].peakHilbert, minY);

	    delete grZ0;
	    delete grZ0Hilbert;
	  }
	  // }
	}
      }
      // Flags

      eventSummary->flags.isGood = 1;
      eventSummary->flags.isPayloadBlast = 0; //!< To be determined.
      eventSummary->flags.nadirFlag = 0; //!< Not sure I will use this.
      eventSummary->flags.strongCWFlag = 0; //!< Not sure I will use this.
      eventSummary->flags.isVarner = 0; //!< Not sure I will use this.
      eventSummary->flags.isVarner2 = 0; //!< Not sure I will use this.
      eventSummary->flags.pulser = AnitaEventSummary::EventFlags::NONE; //!< Not yet.

      delete usefulEvent;

      eventSummaryTree->Fill();
      // delete eventSummary;
    }
    p.inc(entry, nEntries);
  }

  // saves time later
  eventSummaryTree->BuildIndex("eventNumber");

  outFile->Write();
  outFile->Close();

  return 0;
}
