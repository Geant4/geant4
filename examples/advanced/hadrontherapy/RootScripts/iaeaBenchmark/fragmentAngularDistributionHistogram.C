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
 * Root script that produces a histogram of the angular distribution of
 * a certain type of charged fragment (default Z=1).
 * 
 * Notice that this produces a histogram with equally spaced bins for angles.
 * 
 * Comparison to E.Haettner will not yield good results, but shows 
 * alternative approach for analyzing the same data.
 * 
 * @author Gillis Danielsen
 * **************************/

void fragmentAngularDistributionHistogram() {
gStyle->SetOptStat(0000000000); //remove the redundant statbox
   gROOT->SetStyle("clearRetro");
 //this will be used as base for pulling the experimental data
   TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
   dir.ReplaceAll("basic.C","");
   dir.ReplaceAll("/./","/");
   ifstream in;
   
   TString pDepth, fragment, Znum;
   cout << "Enter phantom depth (eg. 27.9, see experimentalData directory for choices): ";
   cin >> pDepth;
   cout << "Enter fragment Z-number (eg. 1): ";
   cin >> Znum;
   cout << "Enter fragment name (Znum 1 -> H,Znum 2->He...): ";
   cin >> fragment;
   
   TString experimentalDataPath = "experimentalData/iaeaBenchmark/angularDistributions/" + pDepth + "/" + fragment + "" + pDepth +".dat";
   TString simulationDataPath = "IAEA_" + pDepth + ".root";
   
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
	   printf(" found %d points\n",nlines);
   //Let's pull in the simulation-data
   //TCanvas *mc = new TCanvas("mc", "Simulation");
   TFile *simulation = TFile::Open(simulationDataPath);
   TNtuple *fragments = (TNtuple*) simulation->Get("fragmentNtuple");

   //Block bellow pulls out the simulation's metadata from the metadata ntuple.
   TNtuple *metadata = (TNtuple*) simulation->Get("metaData");
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


	TString Znum = Form("%i", 1); //hydrogen seems to be the only interesting fragment due to largest amounts
	/* 
	 * A logical bin amount could be calculated as follows, E.Haettner used a detector with 4cm side length, 
	 * the scattering distance is 295.85, thus arctan(4 / 295.85) = 0.774612657 degrees. As 14 degrees is shown
	 * in Haettners graphs and is chosen here as well that would mean approx 14/.77 ~ 18.
	 * This should mean moving the detector one detectorlength. so that the bins "do not overlap". This is ab it of an approx though
	 * because the change in degrees is not linear.
	 * So far these results are though very much inconsistent with Haettners.
	 * 
	 * Although here one should note that the degrees are not linear.
	 * 
	 * The negativeHist is a cludge-fix to get the negative degrees (will be exactly the same)
	*/
	Double_t maxDegrees = 14; //highest plotted degree amount
	int binAmount = 25; //amount of bins in plotted histogram
    TH1F *hist1 = new TH1F("hist1", "Fragment angular distribution", binAmount, 0, maxDegrees); //The histogram for the angular distribution at a set length
	//This histogram needs to be hist1 with a symmetric negative side
	TH1F *symmetricHist = new TH1F("symmetricHist", "Fragment angular distr.", 2*binAmount, -1*maxDegrees, maxDegrees); //bin amount must be even

	//Things needed for the analysis (based on the metadata pulled earlier)
	//ALL UNITS ARE FOR NOW cm
	Double_t detectorSideLength = 4; //40mm
	Double_t scatteringDistance = detectorDistance - phantomCenterDistance;
	TString sdstring = Form("%.4f", scatteringDistance);
	
	TCanvas *c1 = new TCanvas("histograms", "Angular distribution at certain phantom thickness, sd" + sdstring); //This is where we will plot
	
////////////////////////////////////////
//////        Analysis         /////////
////////////////////////////////////////
TString middleY;
//Projection from ntuple to histogram, so that the angle is trigonometrically pulled as the binned variable
	fragments->Project("hist1","57.29577*atan((posZ^2+posY^2)/" + sdstring + ")", "(Z == " + Znum + ")");

//Secondly all bins need to be scaled according to the solid angle of that phi+-deltaPhi circlepart. Before this the curve is "bragg-curvish", what I want is "concaveish".
//this could be done inside the reading of the tuple into the histo, however doing it afterwards improves
	Double_t value, width, deltaPhi, degrees;
	Double_t binNormalization = 1;
	std::cout << "bin-number" << "\t" << "value" << "\t" << "solid angle segment" << endl;
for(int bin = 0; bin <= hist1->GetNbinsX(); bin++){
		value = hist1->GetBinContent(bin); //the incident-particle normalized amount of hits
		width = hist1->GetBinWidth(bin); //so this is degrees/radians
		degrees = hist1->GetBinCenter(bin);
		deltaPhi = width/2;
		binNormalization = 2*TMath::Pi()*(TMath::Cos(TMath::DegToRad()*(degrees-deltaPhi)) - TMath::Cos(TMath::DegToRad()*(degrees+deltaPhi))); //Gunzer-marx uses this , which is a tad of an approximation
		symmetricHist->SetBinContent(bin+hist1->GetNbinsX(), value/(binNormalization*events)); //Solid angle and amount of events
		symmetricHist->SetBinContent(hist1->GetNbinsX()-bin+1, value/(binNormalization*events)); //Solid angle and amount of events
	}
	symmetricHist->SetAxisRange(-2,14);
	symmetricHist->SetXTitle("Angle (degrees)");
	symmetricHist->SetYTitle("(N/N0) / [sr]");
	
	///fragments->Scan("posY:posZ:atan((posZ^2+posY^2)/" + sdstring + ")");
	TF1* fitgaus = new TF1("fitgaus","gaus");
	TF2* fitexpo = new TF1("fitexpo","expo");
	fitgaus->SetLineColor(2);
	fitexpo->SetLineColor(4);	
	symmetricHist->Fit(fitgaus,""); // data should be reshaped a bit for the gaussian (as root seems to have some problems here)
	symmetricHist->Fit(fitexpo, "+"); //aparently two fits on the same histo seem to much for root
	symmetricHist->Draw();
	//symmetricHist->SetMaximum(35); //needed to get E.Haettners data visible
	/*
	//Plots of experimental data are kind of meaningless with a constant-angle-binned histogram.
	ntuple->SetMarkerStyle(22);
    ntuple->SetMarkerColor(kRed);
	ntuple->Draw("y:x","","p,same");
	*/
	
   c1->SaveAs("angularDistributionHistogramWithFits.png");

}
