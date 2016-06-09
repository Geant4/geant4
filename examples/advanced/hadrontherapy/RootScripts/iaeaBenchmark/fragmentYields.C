#include "Riostream.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TApplication.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "THStack.h"
#include "TCut.h"
#include "TString.h"
#include "TMath.h"

/***************************
 * Root script that produces the yield calculations for fragments
 * below the angle of 10 degrees.
 * This is the same numbers E.Haettner has calculated
 * However, there is room for much improvement to be made here.
 * 
 * @author Gillis Danielsen
 * **************************/


void fragmentYields() {


//   gROOT->SetStyle("clearRetro");
 //this will be used as base for pulling the experimental data
   TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
   dir.ReplaceAll("fragmentAngularDistribution.C","");
   dir.ReplaceAll("/./","/");
   ifstream in;

   TString pDepth, fragment, Znum, normToOneAtZeroAngle;
   cout << "Enter phantom depth (eg. 27.9, see experimentalData directory for choices): ";
   cin >> pDepth;

   TString simulationDataPath = "IAEA_" + pDepth + ".root";

   //Let's pull in the simulation-data
   //TFile *MCData = TFile::Open("IAEA_" + pDepth + ".root");
   TFile *MCData = TFile::Open(simulationDataPath);
   TNtuple *fragments = (TNtuple*) MCData->Get("fragmentNtuple");

   //Block bellow pulls out the simulation's metadata from the metadata ntuple.
   TNtuple *metadata = (TNtuple*) MCData->Get("metaData");
   Float_t events, detectorDistance,waterThickness,beamEnergy,energyError,phantomCenterDistance;
   metadata->SetBranchAddress("events",&events);
   metadata->SetBranchAddress("waterThickness",&waterThickness);
   metadata->SetBranchAddress("detectorDistance",&detectorDistance);
   metadata->SetBranchAddress("beamEnergy",&beamEnergy);
   metadata->SetBranchAddress("energyError",&energyError);
   metadata->SetBranchAddress("phantomCenterDistance",&phantomCenterDistance);
   metadata->GetEntry(0); //there is just one row to consider.
   
	//good to keep for ref. G4 might give weird units due to change.
	metadata->Scan();
	
	//ALL UNITS ARE cm!
	Double_t scatteringDistance = detectorDistance - phantomCenterDistance; //temporarily hard-coded, should be distance from target-center to detector

	Double_t r;
	Double_t degrees = 10.0;
	Double_t rMin;
	Double_t rMax;
	TString rMinString;
	TString rMaxString;

	rMinString = "0.00";
	rMaxString = Form("%f", scatteringDistance*TMath::ATan(degrees*TMath::DegToRad()));
	std::cout << "Computes comparable figures to E.haettners yield calculations for zero to ten degree angles." << endl;
	std::cout << "Firsta value is yield of Z=1 per C12 incident paritcle second is for Z=2 etc." << endl;
	std::cout << "verboseness has been limited to make pasting easyer." << endl;
	for(int i = 1; i <= 6; i++){
		int numEntries = fragments->GetEntries("(Z == " + TString::Format("%i",i) + "  && sqrt(posY^2 + posZ^2) < " + rMaxString + "&& sqrt(posY*posY + posZ*posZ) > " + rMinString + ")");
		std::cout << numEntries / (events * .999) << endl; //< for calculation of FWHM
		}
}
