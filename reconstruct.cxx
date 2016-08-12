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

  if(!(argc==3 || argc==2)){
    std::cerr << "Usage 1: " << argv[0] << " [run]" << std::endl;
    std::cerr << "Usage 2: " << argv[0] << " [firstRun] [lastRun]" << std::endl;
    return 1;
  }
  
  std::cout << argv[0] << "\t" << argv[1];
  if(argc==3){std::cout << "\t" << argv[2];}
  std::cout << std::endl;
  const Int_t firstRun = atoi(argv[1]);
  const Int_t lastRun = firstRun; //argc==3 ? atoi(argv[2]) : firstRun;

  CrossCorrelator* cc = new CrossCorrelator();

  TChain* headChain = new TChain("headTree");
  TChain* gpsChain = new TChain("adu5PatTree");
  TChain* calEventChain = new TChain("eventTree");
  
  for(Int_t run=firstRun; run<=lastRun; run++){
    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    // TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/decimatedHeadFile%d.root", run, run);
    headChain->Add(fileName);
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
  cc->addNotch(notch260);
  cc->addNotch(notch370);

  TNamed* comments = new TNamed("comments", "Applied simple, static notch at 260#pm26 MHz and 370#pm26");
  comments->Write();
  delete comments;
    
  TTree* eventSummaryTree = new TTree("eventSummaryTree", "eventSummaryTree");
  // AnitaEventSummary* eventSummary = new AnitaEventSummary();
  AnitaEventSummary* eventSummary = NULL; //new AnitaEventSummary();
  eventSummaryTree->Branch("eventSummary", &eventSummary);
  
  Long64_t nEntries = headChain->GetEntries();
  Long64_t maxEntry = 0; //2513; //33;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  for(Long64_t entry = startEntry; entry < maxEntry; entry++){

    headChain->GetEntry(entry);
    gpsChain->GetEntry(entry);
    calEventChain->GetEntry(entry);
    // gpsChain->GetEntryWithIndex(header->eventNumber);
    // calEventChain->GetEntryWithIndex(header->eventNumber);

    UsefulAnitaEvent* usefulEvent = new UsefulAnitaEvent(calEvent);

    UsefulAdu5Pat usefulPat(pat);
    // cc->correlateEvent(usefulEvent);

    const Int_t myNumPeaks = 2;
    cc->reconstructEvent(usefulEvent, myNumPeaks, myNumPeaks);

    
    eventSummary = new AnitaEventSummary(header, &usefulPat);
    // std::cout << eventSummary->sun.theta << "\t" << eventSummary->sun.phi << std::endl;

    for(Int_t polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
      AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;

      Int_t pointInd=0;
      for(Int_t peakInd=0; peakInd < myNumPeaks; peakInd++){
	
	Double_t minY = 0;

	cc->getCoarsePeakInfo(pol, peakInd,
			      eventSummary->peak[pol][pointInd].value,
			      eventSummary->peak[pol][pointInd].phi,
			      eventSummary->peak[pol][pointInd].theta);

	
	const Int_t coherentDeltaPhi = 0;
	TGraph* grGlobal0 = cc->makeCoherentlySummedWaveform(pol,
							     eventSummary->peak[pol][pointInd].phi,
							     eventSummary->peak[pol][pointInd].theta,
							     coherentDeltaPhi,
							     eventSummary->peak[pol][pointInd].snr);
	
	TGraph* grGlobal0Hilbert = FFTtools::getHilbertEnvelope(grGlobal0);
      
	RootTools::getMaxMin(grGlobal0Hilbert, eventSummary->coherent[pol][peakInd].peakHilbert, minY);
      
	delete grGlobal0;
	delete grGlobal0Hilbert;

	pointInd++;

	cc->getFinePeakInfo(pol, peakInd, 
			    eventSummary->peak[pol][pointInd].value,
			    eventSummary->peak[pol][pointInd].phi,
			    eventSummary->peak[pol][pointInd].theta);

	TGraph* grZ0 = cc->makeUpsampledCoherentlySummedWaveform(pol,
								 eventSummary->peak[pol][pointInd].phi,
								 eventSummary->peak[pol][pointInd].theta,
								 coherentDeltaPhi,
								 eventSummary->peak[pol][pointInd].snr);

	TGraph* grZ0Hilbert = FFTtools::getHilbertEnvelope(grZ0);

	RootTools::getMaxMin(grZ0Hilbert, eventSummary->coherent[pol][pointInd].peakHilbert, minY);

	delete grZ0;
	delete grZ0Hilbert;

	pointInd++;
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
    delete eventSummary;
    p.inc(entry, nEntries);
  }
  
  outFile->Write();
  outFile->Close();

  return 0;
}
