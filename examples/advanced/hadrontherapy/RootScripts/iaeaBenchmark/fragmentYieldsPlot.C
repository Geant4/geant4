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
 * Then plots the data oin comparison to Haettner
 * 
 * There is no fit here yet.
 * 
 * @author Gillis Danielsen
 * **************************/


void fragmentYieldsPlot() {
   gStyle->SetOptStat(0000000000); //remove the for this graphs totally redundant statbox
   int ZNumGiven;
   cout << "Enter fragment Z-number (eg. 1): ";
   cin >> ZNumGiven;

   TCanvas *c1 = new TCanvas("fragmentYieldsPlot", "Total yield of fragments zero to ten degrees as function of depth");
   
   TString fragmentNameChoices[6] = {"H","He","Li","Be","B","C"};
   
   TString fragmentName = fragmentNameChoices[ZNumGiven-1];  
   
   std::cout << fragmentName << endl;
   
   TH1F* dummyHisto = new TH1F("dummyHisto", fragmentName + " yields 0-10 degrees" ,100, 0.0,40); //Dummyhisto fix for missing TNtuple methods.
   dummyHisto->SetXTitle("Depth (cm)");
   dummyHisto->SetYTitle("N/N0");

   ifstream in;
   TString experimentalDataPath = "experimentalData/iaeaBenchmark/yields/TDK" + fragmentName + ".dat";
   
   ifstream in;

   //Pull in ascii/exfor-style data
   in.open(experimentalDataPath);

   Float_t f1,f2;
   Int_t nlines = 0;
   TFile *f = new TFile("fragmentAngularDistribution.root","RECREATE");
   TNtuple *ntuple = new TNtuple("ntuple","Data from ascii file","x:y");
	  
   Char_t DATAFLAG[4];
   Int_t NDATA;
   Char_t n1[15], n2[15];
   in >> DATAFLAG >> NDATA ; // Read EXFOR line: 'DATA 6'
   in >> n1 >> n2; // Read  column titles: 'Energy He B [...]'

   cout <<n1<<"   "<<n2<<"\n";
   while (1) {
      in >> f1 >> f2;
      if (!in.good()) break;
      if (nlines < 500 ) printf("%f %f\n",f1,f2);
      ntuple->Fill(f1,f2);
      nlines++;
   }
   std::cout << "Imported " << nlines << " lines from data-file" << endl;




   TNtuple *simData = new TNtuple("ntuple","Data from ascii file","depth:H:He:Li:Be:B:C");
   
//   gROOT->SetStyle("clearRetro");
 //this will be used as base for pulling the experimental data
   TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
   dir.ReplaceAll("fragmentAngularDistribution.C","");
   dir.ReplaceAll("/./","/");
   ifstream in;
   simData->Fill(0.0,0.0,0.0,0.0,0.0,0.0,1.0);
for(int j = 1; j <= 40;j=j=j+1){
   TString pDepth, fragment, Znum, normToOneAtZeroAngle;
   pDepth = Form("%i",j);
/*
   cout << "Enter phantom depth (eg. 27.9, see experimentalData directory for choices): ";
   cin >> pDepth;
*/
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
	//ALL UNITS ARE cm!
	Double_t scatteringDistance = detectorDistance - phantomCenterDistance; //temporarily hard-coded, should be distance from target-center to detector

	Double_t degrees = 10.0;
	Double_t r, rMin, rMax, graphMaximum = 0.0;
	Double_t norming = events*.999;
	TString rMinString;
	TString rMaxString;

	rMinString = "0.00";
	rMaxString = Form("%f", scatteringDistance*TMath::ATan(degrees*TMath::DegToRad()));

		Double_t H = ((Double_t*) fragments->GetEntries("(Z == " + TString::Format("%i",1) + "  && sqrt(posY^2 + posZ^2) < " + rMaxString + "&& sqrt(posY*posY + posZ*posZ) > " + rMinString + ")")) / norming;
		Double_t He = ((Double_t*) fragments->GetEntries("(Z == " + TString::Format("%i",2) + "  && sqrt(posY^2 + posZ^2) < " + rMaxString + "&& sqrt(posY*posY + posZ*posZ) > " + rMinString + ")")) / norming;
		Double_t Li = ((Double_t*) fragments->GetEntries("(Z == " + TString::Format("%i",3) + "  && sqrt(posY^2 + posZ^2) < " + rMaxString + "&& sqrt(posY*posY + posZ*posZ) > " + rMinString + ")")) / norming;
		Double_t Be = ((Double_t*) fragments->GetEntries("(Z == " + TString::Format("%i",4) + "  && sqrt(posY^2 + posZ^2) < " + rMaxString + "&& sqrt(posY*posY + posZ*posZ) > " + rMinString + ")")) / norming;
		Double_t B = ((Double_t*) fragments->GetEntries("(Z == " + TString::Format("%i",5) + "  && sqrt(posY^2 + posZ^2) < " + rMaxString + "&& sqrt(posY*posY + posZ*posZ) > " + rMinString + ")")) / norming;
		Double_t C = ((Double_t*) fragments->GetEntries("(Z == " + TString::Format("%i",6) + "  && sqrt(posY^2 + posZ^2) < " + rMaxString + "&& sqrt(posY*posY + posZ*posZ) > " + rMinString + ")")) / norming;
		simData->Fill(waterThickness,H,He,Li,Be,B,C);

	}
	simData->Scan();
	simData->SetMarkerStyle(2); //filled dot
	simData->SetMarkerColor(kBlue);
	graphMaximum = TMath::Max(graphMaximum, simData->GetMaximum(fragmentName));
	graphMaximum = TMath::Max(graphMaximum, ntuple->GetMaximum("y"));
	dummyHisto->SetMaximum(graphMaximum + .05*graphMaximum);
	dummyHisto->Draw();
	
	simData->Draw(fragmentName + ":depth","","p,same");

	ntuple->SetMarkerStyle(22); //triangle
    ntuple->SetMarkerColor(kRed);
	ntuple->Draw("y:x","","p,same");
	c1->SaveAs("fragmentYieldsFor" + fragmentName + ".png");
}
