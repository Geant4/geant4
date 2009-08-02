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


void fragmentAngularDistribution() {

//   gROOT->SetStyle("clearRetro");
 //this will be used as base for pulling the experimental data
   TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
   dir.ReplaceAll("basic.C","");
   dir.ReplaceAll("/./","/");
   ifstream in;
   in.open(Form("experimentalData/iaeaBenchmark/H27.9.dat",dir.Data()));
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
   TFile *simulation = TFile::Open("IAEA.root");
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
	*/
    TH1F *hist1 = new TH1F("hist1", "Fragment angular distr.", 18, 0.0, 14); //The histogram for the angular distribution at a set length
	//Things needed for the analysis (based on the metadata pulled earlier)
	//ALL UNITS ARE FOR NOW cm
	Double_t detectorSideLength = 4; //40mm
	Double_t scatteringDistance = detectorDistance - phantomCenterDistance; //temporarily hard-coded, should be distance from target-center to detector
	TString sdstring = Form("%.4f", scatteringDistance);
	
	TCanvas *c1 = new TCanvas("histograms", "Angular distribtuion at certain phantom thickness, sd" + sdstring); //This is where we will plot
	
////////////////////////////////////////
//////        Analysis         /////////
////////////////////////////////////////
TString middleY;
Double_t sinuse;
std::cout << "Check, single detector, different degrees, amounts of hits:\n\n\n";
for(float p = 0.0; p < 12; p = p + 1.0){
middleY = Form("%f",scatteringDistance * TMath::Sin(p/57.3));
std::cout << "\n" << p << " degrees: " << fragments->GetEntries("posY > (" + middleY + " -2.00) && posY < (" + middleY + "+2.00) && posZ > -2 && posZ < 2 && Z == 1");
}
std::cout << "\n";
std::cout << "\n";
//Projection from ntuple to histogram, so that the angle is trigonometrically pulled as the binend variable
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
		std::cout << bin << "\t" << value << "\t" << binNormalization << endl;
		hist1->SetBinContent(bin, value/(binNormalization*events)); //Solid angle and amount of events
		//hist1->SetBinContent(bin, value/(binNormalization*events*width)); //normalized to solid angle, bin width and event count
	}

	///fragments->Scan("posY:posZ:atan((posZ^2+posY^2)/" + sdstring + ")");
	TF1* fitgaus = new TF1("fitgaus","gaus");
	TF2* fitexpo = new TF1("fitexpo","expo");
	//fitgaus->SetLineColor(2);
	fitexpo->SetLineColor(2);	
	hist1->Fit(fitgaus,""); // data should be reshaped a bit for the gaussian (as root seems to have some problems here)
	hist1->Fit(fitexpo, "+"); //aparently two fits on the same histo seem to much for root
	hist1->Draw();
	hist1->SetMaximum(35);
	ntuple->SetMarkerStyle(22);
    ntuple->SetMarkerColor(kRed);
	ntuple->Draw("y:x","","p,same");


   c1->SaveAs("angulardistribution.png");

}
