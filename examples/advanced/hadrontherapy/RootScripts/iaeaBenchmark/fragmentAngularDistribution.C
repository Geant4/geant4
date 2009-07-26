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
/* //this will be used as base for pulling the experimental data
   TString±±±±± dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
   dir.ReplaceAll("basic.C","");
   dir.ReplaceAll("/./","/");
   ifstream in;
   in.open(Form("experimentalData/iaeaBenchmark/fragmentEnergySpctra279mmWater0deg.dat",dir.Data()));
   Float_t f1,f2,f3, f4,f5,f6;
   Int_t nlines = 0;
   TFile *f = new TFile("fragmentEnergyWithAngularDistribution.root","RECREATE");
   TNtuple *ntuple = new TNtuple("ntuple","Data from ascii file","Energy:He:B:H:Li:Be");
	  
   Char_t DATAFLAG[4];
   Int_t NDATA;
   Char_t n1[6], n2[2], n3[2], n4[2], n5[2], n6[2];
   in >> DATAFLAG >> NDATA ; // Read EXFOR line: 'DATA 6'
   in >> n1 >> n2 >> n3 >> n4 >> n5 >> n6; // Read  column titles: 'Energy He B [...]'

   cout <<n1<<" "<<n2<<" "<<n3<<"    "<<n4<<"    "<<n5<<"   "<<n6<<"\n";
   while (1) {
      in >> f1 >> f2 >> f3 >>f4 >> f5 >> f6;
      if (!in.good()) break;
      if (nlines < 500 ) printf("%f  %0.2f %0.2f %0.2f %0.2f %0.2f \n",f1,f2,f3,f4,f5,f6);
      ntuple->Fill(f1,f2,f3,f4,f5,f6);
      nlines++;
   }
*/
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
	 * So far these results are though very much inconsistent with Haettners.
	*/
    TH1F *hist1 = new TH1F("hist1", "Fragment angular distr.", 18, 0.0, 12); //The histogram for the angular distribution at a set length
	//Things needed for the analysis (based on the metadata pulled earlier)
	//ALL UNITS ARE FOR NOW cm
	Double_t detectorSideLength = 4; //40mm
	Double_t scatteringDistance = detectorDistance - phantomCenterDistance; //temporarily hard-coded, should be distance from target-center to detector
	TString sdstring = Form("%.4f", scatteringDistance);
	
	TCanvas *c1 = new TCanvas("histograms", "Angular distribtuion at certain phantom thickness"); //This is where we will plot
	
////////////////////////////////////////
//////        Analysis         /////////
////////////////////////////////////////

//Projection from ntuple to histogram, so that the angle is trigonometrically pulled as the binend variable
	fragments->Project("hist1","57.29577*atan((posZ^2+posY^2)/" + sdstring + ")", "(Z == " + Znum + ")");

//Secondly all bins need to be scaled according to the steradian of that phi+-deltaPhi circlepart. Before this the curve is "bragg-curvish", what I want is "concaveish".
//this could be done inside the reading of the tuple into the histo, however doing it afterwards improves
	Double_t value, width, deltaPhi, degrees;
	Double_t binNormalization = 1;
	std::cout << scatteringDistance << "bin-number" << "\t" << "value" << "\t" << "solid angle segment" << endl;
for(int bin = 0; bin <= hist1->GetNbinsX(); bin++){
		value = hist1->GetBinContent(bin); //the incident-particle normalized amount of hits
		width = hist1->GetBinWidth(bin); //so this is degrees/radians
		degrees = hist1->GetBinCenter(bin);
		deltaPhi = width/2;
		binNormalization = 2*TMath::Pi()*(TMath::Cos(TMath::DegToRad()*(degrees-deltaPhi)) - TMath::Cos(TMath::DegToRad()*(degrees+deltaPhi))); //Gunzer-marx uses this , which is a tad of an approximation
		std::cout << bin << "\t" << value << "\t" << binNormalization << endl;
		hist1->SetBinContent(bin, value/(binNormalization*events));
	}

	///fragments->Scan("posY:posZ:atan((posZ^2+posY^2)/" + sdstring + ")");
	// hist1.Fit("gaus"); // hmm we should reform the data a bit for this
	hist1.Fit("expo");
	hist1->Draw();


   c1->SaveAs("angulardistribution.png");

}
