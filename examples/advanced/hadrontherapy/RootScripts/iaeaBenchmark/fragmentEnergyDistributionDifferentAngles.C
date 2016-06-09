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

gStyle->SetOptStat(0000000000); //remove the for this graphs totally redundant statbox

//   gROOT->SetStyle("clearRetro");

   TString pDepth, fragment, Znum, normToOneAtZeroAngle;
   cout << "Enter phantom depth (eg. 27.9, see experimentalData directory for choices): ";
   cin >> pDepth;  
   TString simulationDataPath = "IAEA_" + pDepth + ".root";
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
   TFile *MCData = TFile::Open("IAEA_200000.root");
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
   TH1F *hist1 = new TH1F("hist1", "", binAmount, 0.0, maxEnergy);
   TH1F *hist2 = new TH1F("hist2", "", binAmount, 0.0, maxEnergy);   
   TH1F *hist3 = new TH1F("hist3", "", binAmount, 0.0, maxEnergy);
   TH1F *hist4 = new TH1F("hist4", "", binAmount, 0.0, maxEnergy);
   TH1F *hist5 = new TH1F("hist5", "", binAmount, 0.0, maxEnergy);
   TH1F *hist6 = new TH1F("hist6", "", binAmount, 0.0, maxEnergy);
   TH1F *hist7 = new TH1F("hist7", "", binAmount, 0.0, maxEnergy);
   TH1F *hist8 = new TH1F("hist8", "", binAmount, 0.0, maxEnergy);
   TH1F *hist9 = new TH1F("hist9", "", binAmount, 0.0, maxEnergy);
for(int k = 1; k <= 6; k++){
		TString Znum = Form("%i", k);
	hist1->SetTitle("Z=" + Znum);


	//ALL UNITS ARE cm!
	Double_t detectorSideLength = 4; //40mm, as e.haettner H1 detector
	Double_t scatteringDistance = detectorDistance - phantomCenterDistance; //temporarily hard-coded, should be distance from target-center to detector
	Double_t degrees; //< actually radians
	
	Double_t r, rMin, rMax, deltaOmega, normFloat;
	TString rMinString, rMaxString, normString;
	TString same = "";
	TString histName;
	TCanvas *c3 = new TCanvas("histograms", "Distribution (at different angles)");
	int i = 0; //so that the degree steps can be varied to unevenly spaced values separate counter is used

	std::cout << "The following numbers also make it possible to make number of fragments comparison to the graph in A1 of E.Haettner\n";
	for(Double_t j = 0.0; j <= 8.0; j=j+1.0){
		i++;
		degrees = j * TMath::DegToRad();
		//std::cout << "plotting for Z = " << Znum << " at " << j << " degrees\n";
		//Distance from straight beam at the requested angle
		r = scatteringDistance * TMath::Tan(degrees);
		//now the "detector is rotated around all possible perpendicularlynangle values to beamline".
		//This forms an annulus with rMin and RMax as otuer and inner radiuses
		//Notice this will give a bit of approximation at small angles where at 0 degrees this gives a round sensor.
		Double_t deltaPhi = TMath::ATan((TMath::Cos(degrees)*detectorSideLength)/(2*scatteringDistance));
		rMin = TMath::Max(0.0,r - (detectorSideLength/(2*TMath::Cos(degrees))));
		rMax = rMin + ((detectorSideLength*TMath::Sin(degrees))/TMath::Tan((TMath::Pi()/2) - degrees - deltaPhi)) + (detectorSideLength*TMath::Cos(degrees));
		rMinString = Form("%f", rMin);
		rMaxString = Form("%f", rMax);
		//normalization of the bins.
		deltaPhi = degrees - TMath::ATan(TMath::Tan(degrees) - detectorSideLength/(2*scatteringDistance)); // this should be around arctan(detectorsidelength/sd)
		if(j != 0.0){
		deltaOmega = 2*TMath::Pi()*(TMath::Cos(TMath::Max(0.0,degrees-deltaPhi)) - TMath::Cos(degrees+deltaPhi));
		}else{
		deltaOmega = 4 * TMath::ASin(pow(detectorSideLength,2.0) / (4*pow(scatteringDistance,2) + pow(detectorSideLength,2)) );
		}
		normFloat = deltaOmega * events * binWidth;
		normString = Form("/%f", normFloat);

	// The following is veryvery ugly relies on a bunch of hardcoded histograms	because other solutions did not work
		histName = Form("hist%i", i);
		if(j != 0.0){
		fragments->Project(histName,"energy", "(Z == " + Znum + " && energy > 0 && sqrt(posY^2 + posZ^2) < " + rMaxString + "&& sqrt(posY*posY + posZ*posZ) > " + rMinString + ")" + normString);
		}else{
		fragments->Project(histName,"energy", "(Z == " + Znum + " && energy > 0 && posZ < " + rMaxString + "&& posY < " + rMaxString + " && posY > 0 && posZ > 0)" + normString);
		}
		int numEntries = fragments->GetEntries("(Z == " + Znum + " && energy > 0 && sqrt(posY^2 + posZ^2) < " + rMaxString + "&& sqrt(posY*posY + posZ*posZ) > " + rMinString + ")");
		std::cout << "\nj: "<< numEntries/(deltaOmega * events) << " entries for " << j;
		
		}

//the ugly hardcoded histograms being plotted
		//0 degrees
		hist1->SetLineColor(kBlue);
		hist1->Draw();
		//1 degree
		hist2->SetLineColor(kGreen);
		hist2->Draw("same"); //add "same" when also plotting 0 degrees
		//2 degrees
		hist3->SetLineColor(kRed);
		hist3->Draw("same");
		//3 degrees
		//hist4->SetLineColor(kGreen + 5);
		//hist4->Draw("same");
		//4 degrees
		hist5->SetLineColor(kGreen + 3); //gives a darker shade of green
		hist5->Draw("same");
		//5 degrees
		//hist6->SetLineColor(kRed);
		//hist6->Draw("same");
		//6 degrees
		hist7->SetLineColor(kRed);
		hist7->Draw("same");
		//7 degrees
		//hist8->SetLineColor(kRed);
		//hist8->Draw("same");
		//8 degrees
		//hist9->SetLineColor(kRed);
		//hist9->Draw("same");

   // Legends for the data
   leg = new TLegend(0.9,0.7,1,1);  //coordinates are fractions
   leg->SetHeader("Angles");
   leg->AddEntry(hist1,"0","l");
   leg->AddEntry(hist2,"1","l");
   leg->AddEntry(hist3,"2","l");
   leg->AddEntry(hist5,"4","l");
   leg->AddEntry(hist7,"6","l");
   leg->Draw();

   c3->SaveAs("AEDistrib" + Znum + ".png");
}

   in.close();

   f->Write();
}
