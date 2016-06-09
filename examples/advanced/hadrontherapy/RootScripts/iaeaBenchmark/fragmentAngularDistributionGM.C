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
 * Root script that produces a graph of the angular distribution of
 * a certain type of charged fragment (default Z=1).
 * 
 * Results are not stored in a histogram.
 * 
 * This can be compared to measurements made with a square 
 * detector that is being moved around. Such as that of E.Haettner[1].
 * 
 * Results are normalized to the 0-angle because documentation on E.Haettner's 
 * normalization is not found.
 * 
 * @author Gillis Danielsen
 * **************************/


void fragmentAngularDistributionGM() {

	TCanvas *c1 = new TCanvas("AngularDistribution", "Angular distribution with discrete measurement annuluses");
	
  gROOT->SetStyle("clearRetro");
 //this will be used as base for pulling the experimental data
   TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
   ifstream in;

   //Settings for analysis
   TString pDepth = "5.9"; //set ehre what depth to analyse, requries suitable root file and data file
   TString fragment = "H"; //string of what fragment is being looked at H, He, Li ...
   TString Anum = "1"; //set here what will be put in the selection for z, does not automatically change imoprted data.
   //TString experimentalDataFile = "experimentalData/iaeaBenchmark/angularDistributions/" + pDepth + "/" + fragment + pDepth + ".dat";
   //in.open(Form("experimentalData/iaeaBenchmark/angularDistributions/" + pDepth + "/" + fragment + pDepth + ".dat",dir.Data()));
   in.open(Form("experimentalData/iaeaBenchmark/alternativeGMBeamAD.dat",dir.Data())); //fixme: redundant dir.data()
   std::cout << "experimentalData/iaeaBenchmark/angularDistributions/" + pDepth + "/" + fragment + pDepth + ".dat" << endl;
   std::cout << "didata" << dir.Data() << endl;
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

   //Let's pull in the simulation-data
   //TFile *MCData = TFile::Open("IAEA_" + pDepth + ".root");
   TFile *MCData = TFile::Open("IAEA_G-M.root");
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
	Double_t detectorSideLength = 4; //40mm, as e.haettner H1 detector, G-M uses 1.6cm sidelength (which is however to small for < 1 million events files
	Double_t scatteringDistance = detectorDistance - phantomCenterDistance; //temporarily hard-coded, should be distance from target-center to detector

	Double_t r;
	Double_t degrees;
	Double_t rMin;
	Double_t rMax;
	TString rMinString;
	TString rMaxString;
	Double_t maxValue = 0.0; //When normalizing to 1 this will allways end up being 1)
	
	int i = 0; //so that the degree steps can be varied to unevenly spaced values separate counter is used
	TNtuple* distrib = new TNtuple("angularDistrib","FragmentAngularDistrib","angle:particleAmount:normalized");
	std::cout << "Fragments comparison to the graphs in appendices of E.Haettner\n";
	std::cout << "Scattering distance: " << scatteringDistance << " cm" << endl ;
	std::cout << "(scattering distance may vary with data-files too, see haettner A.1." << endl << endl;
	//This will norm it to the zero degree entry to get rid of Emma's weird normalization
	//fixme detectorsidelengthstring
	rMinString = "0.00";
	rMaxString = Form("%f", detectorSideLength/2);
	//SA of a square Double_t zeroSAsquare = 4 * TMath::ASin(pow(detectorSideLength,2.0) / (4*pow(scatteringDistance,2) + pow(detectorSideLength,2)) );
	//SA with square approx squareApprox = (4*4) / (4*3.14*scatteringDistance*scatteringDistance);
	
	//First calculates the normalization from the zero position
	//Normalization by events becomes redundant but is left in place for future needs

	Double_t deltaPhi = TMath::ATan((detectorSideLength/2)/scatteringDistance); //Angle where side of detector is found
	/*
	//Alternative normalization, here zero position is also done with annulus where rMin=0, the actual detector is a square though
	//Difference with this approach and the other is very small
	Double_t normEntries = fragments->GetEntries("(Z == " + Anum + " && energy > 0 && sqrt(posY^2 + posZ^2) < " + rMaxString + "&& sqrt(posY*posY + posZ*posZ) > " + rMinString + ")");
	Double_t zeroSA = 2*TMath::Pi()*(TMath::Cos(0) - TMath::Cos(deltaPhi));
	*/
	//fragments->Scan();
	//Results are normalized by a square detector mimicing H1 with center at 0 degrees.
	Double_t normEntries = fragments->GetEntries("(A == " + Anum + " && posY < " + rMaxString + " && posY > -" + rMaxString + " &&  posZ > -" + rMaxString + " && posZ < " + rMaxString + ")");
	Double_t zeroSA = 4 * TMath::ASin(pow(detectorSideLength,2.0) / (4*pow(scatteringDistance,2) + pow(detectorSideLength,2)) );
	Double_t zeroNorm = normEntries / (events * zeroSA);
	distrib->Fill(0,normEntries,zeroNorm); //< degrees, entyamount, normalized result for graph
	//fragments->Scan();
	std::cout << "Norming events: " << normEntries << endl;
	//Loop through all other wanted angles, too large angles will fall outside reach of phantom window.
	for(Double_t j = deltaPhi*TMath::RadToDeg(); j <= 15.0; j=j+.05){
		i++;
		degrees = j * TMath::DegToRad();
		//Distance from straight beam at the requested angle
		r = scatteringDistance * TMath::Tan(degrees);
		//now the "detector is rotated around all possible perpendicularlynangle values to beamline".
		//This forms an annulus with rMin and Rmax as outer and inner radiuses
		//this will give a bit of approximation at small angles and at 0 degrees this gives a completely round sensor.
		Double_t deltaPhi = TMath::ATan((TMath::Cos(degrees)*detectorSideLength)/(2*scatteringDistance));
		rMin = TMath::Max(0.0,r - (detectorSideLength/(2*TMath::Cos(degrees))));
		rMax = rMin + ((detectorSideLength*TMath::Sin(degrees))/TMath::Tan((TMath::Pi()/2) - degrees - deltaPhi)) + (detectorSideLength*TMath::Cos(degrees));
		rMinString = Form("%f", rMin);
		rMaxString = Form("%f", rMax);		/*
		* From Gunzert-marx. Solid angle of annulus with rmin trmax, 
		* a bit of an aproximation especially at small phi.
		*/
		Double_t deltaOmega = 2*TMath::Pi()*(TMath::Cos(TMath::Max(0.0,degrees-deltaPhi)) - TMath::Cos(degrees+deltaPhi));
		int numEntries = fragments->GetEntries("(A == " + Anum + " && sqrt(posY^2 + posZ^2) < " + rMaxString + "&& sqrt(posY*posY + posZ*posZ) > " + rMinString + ")");
		//distrib->Fill(j,numEntries,numEntries/(deltaOmega * events * zeroNorm)); //< degrees, entyamount, normalized result for graph
		distrib->Fill(j,numEntries,numEntries/(deltaOmega * events));
		//distrib->Fill(-j,numEntries,numEntries/(deltaOmega * events * zeroNorm)); //< To get gaussian shape better visible
		distrib->Fill(-j,numEntries,numEntries/(deltaOmega * events));
		maxValue = TMath::Max(maxValue, numEntries/(deltaOmega * events)); //< for calculation of FWHM
		}
	distrib->SetMarkerStyle(2); //filled dot
	distrib->SetMarkerColor(kBlue);
	distrib->Draw("normalized:angle","angle > -3 && angle < 14","p"); //similar axises to e.haettner
	ntuple->SetMarkerStyle(22); //triangle
    ntuple->SetMarkerColor(kRed);
	Float_t zeroPosData; //This is where we store what we norm the experimental data with
	Float_t zeroPosAngle; //okay, so this should be zero, but regrettably is not allways that
	ntuple->SetBranchAddress("y",&zeroPosData);
	ntuple->SetBranchAddress("x",&zeroPosAngle);
	int row = 0;	
	ntuple->GetEntry(row); //Pull the first row, usually is the right one
	while(zeroPosAngle*zeroPosAngle > .01){
		row++;
		ntuple->GetEntry(row);
		if(row == ntuple->GetEntries()){
			std::cerr << "Could not find zero angle data in imported experimental data. Change normalization or relax exactness of this check." << endl;
			exit();
			}
		}
	std::cout << "For zero-position of experimental data using angle " << zeroPosAngle << " with amount " << zeroPosData << " on row " << row << endl;
	TString experimentalNorm = Form("(1/%f)*", zeroPosData);
	ntuple->Draw(experimentalNorm + "y:x","","p,same");

	//Calculate FWHM, ineffective solution. (but works also without normalizing to 1)
	Float_t fwhm = 0.0, middle = 0.0, currentX, currentY;
	distrib->SetBranchAddress("normalized",&currentY);
	distrib->SetBranchAddress("angle",&currentX);
	for(int i = 0; i < distrib->GetEntries(); i++){
		distrib->GetEntry(i);
		if(pow(maxValue/2 - middle, 2.0) > pow(maxValue/2 - currentY, 2.0)){
				fwhm = 2*currentX;
				middle = currentY;
			}else{
				}
		}
		std::cout << "Calculated (closest point) FWHM of Monte-Carlo simulation to be: " << fwhm << " degrees" << endl; 
	
	/*
	 * This code is left here because it allows to calcualte the values without using annuluses
	 * this of course is a bit more like the experimental data but statistically less precise.
	for(Double_t p = 0.0;p < 14.0; p = p + 1.0){
	rMinString = "0.00";
	rMaxString = "2.00";
	TString plusY = Form("(posY - %f)", p);
	TString plusZ = Form("(posZ - %f)", p);
		Double_t deltaPhi = TMath::ATan((detectorSideLength/2)/scatteringDistance);
	Double_t normEntries = fragments->GetEntries("(Z == " + Anum + " && energy > 0 && sqrt(" + plusY + "^2 + " + plusZ + "^2) < " + rMaxString + "&& sqrt("+plusY + "^2 + " + plusZ + "^2) > " + rMinString + ")");
	//std::cout << "(Z == " + Anum + " && energy > 0 && sqrt(" + plusY + "^2 + " + plusZ + "^2) < " + rMaxString + "&& sqrt("+plusY + "^2 + " + plusZ + "^2) > " + rMinString + ")";
	Double_t zeroSA = 2*TMath::Pi()*(TMath::Cos(0) - TMath::Cos(deltaPhi));
	std::cout << "with " << p << "cm the amount is " << normEntries << " / " << normEntries /(events*zeroSA) << endl;
		}
		*/
	
	
	c1->SaveAs("AD_for_Z_" + Anum + "_ComparedToGM.png");
	in.close();
	f->Write();
}
