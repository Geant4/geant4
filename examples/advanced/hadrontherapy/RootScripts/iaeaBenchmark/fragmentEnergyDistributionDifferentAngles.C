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

void fragmentEnergyDistributionDifferentAngles() {

//   gROOT->SetStyle("clearRetro");

   TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
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
   

   //Let's pull in the simulation-data
   TFile *MCData = TFile::Open("IAEA.root");
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
	std::cout << "Recieved metadata-row: " << events << " " <<  detectorDistance << " " << waterThickness << " " << beamEnergy << " " << energyError << " " << phantomCenterDistance;


//A lot of hardcoded histograms, ugly
   Double_t binAmount = 50.0; //casting from int failed somehow, so in float temporarily, fixme
   Double_t maxEnergy = 450.0;
   Double_t binWidth = maxEnergy / binAmount;
   TH1F *hist1 = new TH1F("hist1", "1 degrees", binAmount, 0.0, maxEnergy);
   TH1F *hist2 = new TH1F("hist2", "2 degrees", binAmount, 0.0, maxEnergy);   
   TH1F *hist3 = new TH1F("hist3", "3 degrees", binAmount, 0.0, maxEnergy);
   TH1F *hist4 = new TH1F("hist4", "4 degrees", binAmount, 0.0, maxEnergy);
   TH1F *hist5 = new TH1F("hist5", "5 degrees", binAmount, 0.0, maxEnergy);
   TH1F *hist6 = new TH1F("hist6", "6 degrees", binAmount, 0.0, maxEnergy);
   TH1F *hist7 = new TH1F("hist7", "7 degrees", binAmount, 0.0, maxEnergy);
   TH1F *hist8 = new TH1F("hist8", "8 degrees", binAmount, 0.0, maxEnergy);
   TH1F *hist9 = new TH1F("hist9", "9 degrees", binAmount, 0.0, maxEnergy);
for(int k = 1; k <= 6; k++){
		TString Znum = Form("%i", k);
	hist1->SetTitle("Z=" + Znum);


	//ALL UNITS ARE cm!
	Double_t detectorSideLength = 4; //40mm, as e.haettner H1 detector
	Double_t scatteringDistance = detectorDistance - phantomCenterDistance; //temporarily hard-coded, should be distance from target-center to detector
	Double_t degrees = 1.0 * TMath::DegToRad(); //< actually radians
		
	//Distance from straight beam at the requested angle
	Double_t r = scatteringDistance * TMath::Tan(degrees);
	//now the "detector is rotated around all possible perpendicularlynangle values to beamline"
	Double_t rMin = TMath::Max(0.0,r - (detectorSideLength/2));
	Double_t rMax = r + (detectorSideLength/2);
	
	TString rMinString(Form("%f", rMin));
	TString rMaxString(Form("%f", rMax));
	
	//"(Z == 1 && energy > 45 && sqrt(posY^2 + posZ^2) < " + rMaxString + "&& sqrt(posY*posY + posZ*posZ) > " + rMinString + ")"
	TString same = "";
	TString histName;

	TCanvas *c3 = new TCanvas("histograms", "Distribution (at different angles)");
	THStack hs("hs","test stacked histograms");
	int i = 0; //so that the degree steps can be varied to unevenly spaced values separate counter is used
	for(Double_t j = 1.0; j <= 9.0; j=j+1.0){
		i++;
		degrees = j * TMath::DegToRad();
		std::cout << "plotting for Z = " << Znum << " at " << j << " degrees\n";
		//Distance from straight beam at the requested angle
		r = scatteringDistance * TMath::Tan(degrees);
		//now the "detector is rotated around all possible perpendicularlynangle values to beamline".
		//This forms an annulus with rMin and RMax as otuer and inner radiuses
		//Notice this will give a bit of approximation at small angles where at 0 degrees this gives a round sensor.
		rMin = TMath::Max(0.0,r - (detectorSideLength/2));
		rMax = r + (detectorSideLength/2);
		rMinString = Form("%f", rMin);
		rMaxString = Form("%f", rMax);
		
		//normalization of the bins.
		deltaPhi = degrees - TMath::ATan(TMath::Tan(degrees) - detectorSideLength/scatteringDistance); // this should be around arctan(detectorsidelength/sd)
		deltaOmega = 2*TMath::Pi()*(TMath::Cos(degrees-deltaPhi) - TMath::Cos(degrees+deltaPhi));
		normFloat = deltaOmega * events * binWidth;
		normString = Form("/%f", normFloat);

	// The following is veryvery ugly relies on a bunch of hardcoded histograms	because other solutions did not work
		histName = Form("hist%i", i);
		fragments->Project(histName,"energy", "(Z == " + Znum + " && energy > 0 && sqrt(posY^2 + posZ^2) < " + rMaxString + "&& sqrt(posY*posY + posZ*posZ) > " + rMinString + ")" + normString);
		}

//the ugly hardcoded histograms being plotted
		hist1->Draw();
		hist2->Draw("same");
		hist3->Draw("same");
		hist4->Draw("same");
		hist5->Draw("same");
		hist6->Draw("same");
		hist7->Draw("same");
		hist8->Draw("same");
		hist9->Draw("same");

   c3->SaveAs("distrib" + Znum + ".png");
}

   in.close();

   f->Write();
}
