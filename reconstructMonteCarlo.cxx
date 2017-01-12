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

  if(argc!=4){
    std::cerr << "Usage: " << argv[0] << " [run] [energyExponent] [half]" << std::endl;
    std::cerr << "CBA to read in floats... set [half] to non-zero to do half decade above [energyExponent]" << std::endl;
    return 1;
  }

  std::cout << argv[0] << "\t" << argv[1];
  if(argc==3){std::cout << "\t" << argv[2];}
  std::cout << std::endl;
  const Int_t firstRun = atoi(argv[1]);
  const Int_t lastRun = firstRun; //argc==3 ? atoi(argv[2]) : firstRun;
  const Int_t energy = atoi(argv[2]);; //argc==3 ? atoi(argv[2]) : firstRun;
  const Int_t andAHalf = atoi(argv[3]);; //argc==3 ? atoi(argv[2]) : firstRun;

  TString mcDir = (energy < 18 || energy > 22) ? "JustNoise" : TString::Format("Energy_E%d", energy);
  if(energy >= 18 && andAHalf > 0){
    mcDir += "half";
  }
  std::cout << mcDir.Data() << std::endl;
  // return 0;

  CrossCorrelator* cc = new CrossCorrelator();

  TChain* headChain = new TChain("headTree");
  TChain* gpsChain = new TChain("adu5PatTree");
  // TChain* calEventChain = new TChain("eventTree");
  TChain* usefulChain = new TChain("eventTree");
  // TChain* iceChain = new TChain("passing_events");
  TChain* iceChain = new TChain("truthAnitaTree");

  for(Int_t run=firstRun; run<=lastRun; run++){
    // TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    TString fileName = TString::Format("~/UCL/ANITA/newMonteCarlo/%s/run%d/SimulatedAnitaHeadFile%d.root", mcDir.Data(), run, run);
    headChain->Add(fileName);

    fileName = TString::Format("~/UCL/ANITA/newMonteCarlo/%s/run%d/SimulatedAnitaGpsFile%d.root", mcDir.Data(), run, run);
    gpsChain->Add(fileName);

    fileName = TString::Format("~/UCL/ANITA/newMonteCarlo/%s/run%d/SimulatedAnitaEventFile%d.root", mcDir.Data(), run, run);
    usefulChain->Add(fileName);

    // fileName = TString::Format("~/UCL/ANITA/newMonteCarlo/%s/run%d/icefinal%d.root", mcDir.Data(), run, run);
    fileName = TString::Format("~/UCL/ANITA/newMonteCarlo/%s/run%d/SimulatedAnitaTruthFile%d.root", mcDir.Data(), run, run);
    iceChain->Add(fileName);
  }

  if(headChain->GetEntries()==0){
    std::cerr << "Unable to find header files!" << std::endl;
    return 1;
  }
  if(gpsChain->GetEntries()==0){
    std::cerr << "Unable to find gps files!" << std::endl;
    return 1;
  }
  if(usefulChain->GetEntries()==0){
    std::cerr << "Unable to find calEvent files!" << std::endl;
    return 1;
  }
  if(iceChain->GetEntries()==0){
    std::cerr << "Unable to find icefinal files!" << std::endl;
    return 1;
  }

  RawAnitaHeader* header = NULL;
  headChain->SetBranchAddress("header", &header);
  double weight = 0;
  headChain->SetBranchAddress("weight",       &weight           );

  Adu5Pat* pat = NULL;
  gpsChain->SetBranchAddress("pat", &pat);
  UsefulAnitaEvent* usefulEvent = NULL;
  usefulChain->SetBranchAddress("event", &usefulEvent);

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

  const Int_t myNumPeaksCoarse = 2;
  const Int_t myNumPeaksFine = 2;
  const Int_t coherentDeltaPhi = 0;

  TNamed* comments = new TNamed("comments", "Applied simple, static notch at 260#pm26 MHz, 370#pm26 MHz, 762#pm8 MHz");
  comments->Write();
  delete comments;

  TNamed* comments2 = new TNamed("comments2",
				 TString::Format("%d coarse peaks, %d fine peaks",
						 myNumPeaksCoarse, myNumPeaksFine).Data());
  comments2->Write();
  delete comments2;

  double sourceLat, sourceLon, sourceAlt, sourceMag;
  double thetaWave,phiWave;
  double weight2;
  // double posnu[3];
  // double rfexit[5][3];
  // double r_bn[3], r_enterice[3], r_in[3];
  double e_component;
  // int inu;

  // iceChain->SetBranchAddress("inu",          &inu              );
  iceChain->SetBranchAddress("sourceLon",    &sourceLon        );
  iceChain->SetBranchAddress("sourceLat",    &sourceLat        );
  iceChain->SetBranchAddress("sourceAlt",    &sourceAlt        );
  // iceChain->SetBranchAddress("sourceMag",    &sourceMag        );
  // iceChain->SetBranchAddress("posnu",        &posnu            );
  // iceChain->SetBranchAddress("rfexit",       &rfexit           );
  // iceChain->SetBranchAddress("r_bn",         &r_bn             );
  // iceChain->SetBranchAddress("r_enterice",   &r_enterice       );
  // iceChain->SetBranchAddress("r_in",         &r_in             );
  iceChain->SetBranchAddress("e_component",  &e_component      );

  TTree* eventSummaryTree = new TTree("eventSummaryTree", "eventSummaryTree");
  // AnitaEventSummary* eventSummary = new AnitaEventSummary();
  AnitaEventSummary* eventSummary = NULL; //new AnitaEventSummary();
  eventSummaryTree->Branch("eventSummary"	, &eventSummary);
  eventSummaryTree->Branch("header"		, &header);
  eventSummaryTree->Branch("pat"		, &pat    );
  eventSummaryTree->Branch("thetaExpectedDeg"	, &thetaWave );
  eventSummaryTree->Branch("phiExpectedDeg"	, &phiWave   );
  eventSummaryTree->Branch("weight"		, &weight2    );
  eventSummaryTree->Branch("e_component"	, &e_component);

  Long64_t nEntries = headChain->GetEntries();
  Long64_t maxEntry = 0; //2513; //33;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  for(Long64_t entry = startEntry; entry < maxEntry; entry++){

    headChain->GetEntry(entry);

    gpsChain->GetEntry(entry);
    usefulChain->GetEntry(entry);
    iceChain->GetEntry(entry);
    weight2 = weight;

    UsefulAdu5Pat usefulPat(pat);

    // cc->reconstructEvent(usefulEvent, myNumPeaksCoarse, myNumPeaksFine);
    AnitaPol::AnitaPol_t peakPol = cc->reconstructEventPeakPol(usefulEvent, myNumPeaksCoarse, myNumPeaksFine);

    eventSummary = new AnitaEventSummary(header, &usefulPat);
    // std::cout << eventSummary->sun.theta << "\t" << eventSummary->sun.phi << std::endl;

    Double_t minY = 0;

    for(Int_t polInd=0; polInd < AnitaPol::kNotAPol; polInd++){

      AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;

      if(pol!=peakPol){
	for(Int_t peakInd=0; peakInd < myNumPeaksCoarse; peakInd++){


	  cc->getCoarsePeakInfo(pol, peakInd,
				eventSummary->peak[pol][peakInd].value,
				eventSummary->peak[pol][peakInd].phi,
				eventSummary->peak[pol][peakInd].theta);

	  TGraph* grGlobal0 = cc->makeCoherentlySummedWaveform(pol,
							       eventSummary->peak[pol][peakInd].phi,
							       eventSummary->peak[pol][peakInd].theta,
							       coherentDeltaPhi,
							       eventSummary->peak[pol][peakInd].snr);



	  if(grGlobal0!=NULL){

	    TGraph* grGlobal0Hilbert = FFTtools::getHilbertEnvelope(grGlobal0);

	    RootTools::getMaxMin(grGlobal0Hilbert, eventSummary->coherent[pol][peakInd].peakHilbert, minY);

	    delete grGlobal0;
	    delete grGlobal0Hilbert;
	  }
	}
      }
      else{
	for(Int_t peakInd=0; peakInd < myNumPeaksFine; peakInd++){
	  cc->getFinePeakInfo(pol, peakInd,
			      eventSummary->peak[pol][peakInd].value,
			      eventSummary->peak[pol][peakInd].phi,
			      eventSummary->peak[pol][peakInd].theta);

	  TGraph* grZ0 = cc->makeUpsampledCoherentlySummedWaveform(pol,
								   eventSummary->peak[pol][peakInd].phi,
								   eventSummary->peak[pol][peakInd].theta,
								   coherentDeltaPhi,
								   eventSummary->peak[pol][peakInd].snr);

	  // std::cout << grZ0 << "\t" << eventSummary->peak[pol][peakInd].snr << std::endl;

	  if(grZ0!=NULL){
	    TGraph* grZ0Hilbert = FFTtools::getHilbertEnvelope(grZ0);

	    RootTools::getMaxMin(grZ0Hilbert, eventSummary->coherent[pol][peakInd].peakHilbert, minY);

	    delete grZ0;
	    delete grZ0Hilbert;
	  }
	}
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


    usefulPat.getThetaAndPhiWave(sourceLon, sourceLat, sourceAlt, thetaWave,phiWave);
    thetaWave*=TMath::RadToDeg();
    phiWave*=TMath::RadToDeg();


    // delete usefulEvent;

    eventSummaryTree->Fill();
    delete eventSummary;
    p.inc(entry, nEntries);

  }

  outFile->Write();
  outFile->Close();

  return 0;
}
