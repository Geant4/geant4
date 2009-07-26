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
   /*
   ntuple->SetMarkerStyle(5);
   ntuple->Draw("He:Energy","","l");
   ntuple->Draw("B:Energy","","l,Same");
   ntuple->Draw("H:Energy","","l,Same");
   ntuple->Draw("Li:Energy","","l,Same");
   ntuple->Draw("Be:Energy","","l,Same");
   printf(" found %d points\n",nlines);
*/ 
   //Let's pull in the simulation-data
   //TCanvas *mc = new TCanvas("mc", "Simulation");
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
/*
 //This is perhaps deprecated?
   Double_t ScaleHelium = 1/(MC_helium->Integral());
   Double_t ScaleHydrogen = 1/(MC_hydrogen->Integral()); 
   //x should also be scaled to per nucleon
   
   MC_helium->Scale(ScaleHelium);
//   MC_helium->Draw("");
   printf("Scaled helium by %.9f\n",ScaleHelium);

   MC_hydrogen->Scale(ScaleHydrogen);
   MC_hydrogen->SetLineColor(kRed);
//   MC_hydrogen->Draw("Same");
   printf("Scaled hydrogen by %.9f\n",ScaleHydrogen);

   
   TH1F *histH = new TH1F("histH", "Hydrogen", 60, 0.0, 450.0);
   TH1F *histHe = new TH1F("histHe", "Helium", 60, 0.0, 450.0);
   histHe->SetLineColor(kRed);
   TH1F *histLi = new TH1F("histLi", "Lithium", 60, 0.0, 450.0);
   histLi->SetLineColor(kBlue);
   TH1F *histBe = new TH1F("histBe", "Beryllium", 60, 0.0, 450.0);
   histBe->SetLineColor(kGreen);
   TH1F *histB = new TH1F("histB", "Boron", 60, 0.0, 450.0);
   histB->SetLineColor(kYellow);
   TH1F *histC = new TH1F("histC", "Carbon", 60, 0.0, 450.0);
*/

for(int k = 1; k <= 6; k++){
		TString Znum = Form("%i", k);
//the uglyugly hardcoded histograms
   TH1F *hist1 = new TH1F("hist1", "Z=" + Znum, 60, 0.0, 450.0);
   TH1F *hist2 = new TH1F("hist2", "2 degrees", 60, 0.0, 450.0);   
   TH1F *hist3 = new TH1F("hist3", "3 degrees", 60, 0.0, 450.0);
   TH1F *hist4 = new TH1F("hist4", "4 degrees", 60, 0.0, 450.0);
   TH1F *hist5 = new TH1F("hist5", "5 degrees", 60, 0.0, 450.0);
   TH1F *hist6 = new TH1F("hist6", "6 degrees", 60, 0.0, 450.0);
   TH1F *hist7 = new TH1F("hist7", "7 degrees", 60, 0.0, 450.0);
   TH1F *hist8 = new TH1F("hist8", "8 degrees", 60, 0.0, 450.0);
   TH1F *hist9 = new TH1F("hist9", "9 degrees", 60, 0.0, 450.0);

	//TH2F

	//TH1F* histPos = new TH1F("histPos", "check position",100,-2000,2000);

   fragments->SetLineColor(kRed);
   fragments->SetMarkerStyle(22);
   
   //fragments->Draw("posY:posZ", "abs(posZ) < 2000 && abs(posY) < 2000");


   fragments->SetLineColor(kGreen);

	//ALL UNITS ARE FOR NOW cm
	Double_t detectorSideLength = 4; //40mm
	Double_t scatteringDistance = detectorDistance - phantomCenterDistance; //temporarily hard-coded, should be distance from target-center to detector
	Double_t degrees = 1.0 * TMath::DegToRad(); //< Well okay, radians :)
	
	//Double_t r = TMath::Sqrt(posY*posY + posZ*posZ); //just fed through cut, but here for clarity
	
	//Distance from straight beam at the requested angle
	Double_t r = scatteringDistance * TMath::Tan(degrees);
	//now the "detector is rotated around all possible perpendicularlynangle values to beamline"
	Double_t rMin = TMath::Max(0.0,r - (detectorSideLength/2));
	Double_t rMax = r + (detectorSideLength/2);
	
	TString rMinString(Form("%f", rMin));
	TString rMaxString(Form("%f", rMax));
	
	//normalization.
	//fixme, oh my god there will be ridiculously much rounding error here
	Double_t deltaPhi = degrees - TMath::ATan(TMath::Tan(degrees) - detectorSideLength/scatteringDistance); // this should be around arctan(detectorsidelength/sd)
	Double_t deltaOmega = 2*TMath::Pi()*(TMath::Cos(degrees-deltaPhi) - TMath::Cos(degrees+deltaPhi));
	Double_t normFloat = deltaOmega * events;
	
	//std::cout << "deltaPhi " << deltaPhi << "deltaOmega" << deltaOmega << "normFloat" << normFloat;
	
	TString normString(Form("/%f", normFloat));
	
	//"(Z == 1 && energy > 45 && sqrt(posY^2 + posZ^2) < " + rMaxString + "&& sqrt(posY*posY + posZ*posZ) > " + rMinString + ")"
	TString same = "";
	TString histName;
////////////////////////////////////////
////////////////////////////////////////
////////////////////////////////////////
	TCanvas *c3 = new TCanvas("histograms", "Distribution (at different angles)");
	THStack hs("hs","test stacked histograms");
	int i = 0; //so that the degree steps can be varied a separate counter is needed
	for(Double_t j = 1.0; j <= 9.0; j=j+1.0){
	i++;
	degrees = j * TMath::DegToRad();
	std::cout << "plotting for " << j << " " << degrees << "\n";
	//Distance from straight beam at the requested angle
	r = scatteringDistance * TMath::Tan(degrees);
	//now the "detector is rotated around all possible perpendicularlynangle values to beamline"
	rMin = TMath::Max(0.0,r - (detectorSideLength/2));
	rMax = r + (detectorSideLength/2);
	
	rMinString = Form("%f", rMin);
	rMaxString = Form("%f", rMax);
	
	//normalization.
	//fixme, oh my god there will be ridiculously much rounding error here
	deltaPhi = degrees - TMath::ATan(TMath::Tan(degrees) - detectorSideLength/scatteringDistance); // this should be around arctan(detectorsidelength/sd)
	deltaOmega = 2*TMath::Pi()*(TMath::Cos(degrees-deltaPhi) - TMath::Cos(degrees+deltaPhi));
	normFloat = deltaOmega * events;
	normString = Form("/%f", normFloat);
	
	//This does not work, just redraws the same, DrawClone did "work" but it also redrew the labels every time making the graph undreadable
	//fragments->Draw("energy >> histH", "(Z == 1 && energy > 45 && sqrt(posY^2 + posZ^2) < " + rMaxString + "&& sqrt(posY*posY + posZ*posZ) > " + rMinString + ")" + normString, same);
/*
 * //because the code below and many of its variations did not work i will have to do this with an incredible cludgeball
	fragments->Project("histH","energy", "(Z == 1 && energy > 45 && sqrt(posY^2 + posZ^2) < " + rMaxString + "&& sqrt(posY*posY + posZ*posZ) > " + rMinString + ")" + normString);
	hs.Add(histH->Clone());
*/

// The following is veryvery ugly relies on a bunch of hardcoded histograms	
	histName = Form("hist%i", i);
	std::cout << histName << "\n";
	fragments->Project(histName,"energy", "(Z == " + Znum + " && energy > 0 && sqrt(posY^2 + posZ^2) < " + rMaxString + "&& sqrt(posY*posY + posZ*posZ) > " + rMinString + ")" + normString);

		}

//the uglyugly hardcoded histograms being filled
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

////////////////////////////////////////
////////////////////////////////////////
////////////////////////////////////////	
	//fragments->Draw("energy >> histH", "(Z == 1 && energy > 45 && sqrt(posY^2 + posZ^2) < " + rMaxString + "&& sqrt(posY*posY + posZ*posZ) > " + rMinString + ")" + normString);
	//fragments->Draw("energy >> histHe", "(Z == 2 && energy > 45 && sqrt(posY^2 + posZ^2) < " + rMaxString + "&& sqrt(posY*posY + posZ*posZ) > " + rMinString + ")" + normString, "same");
	
	
	//fragments->Draw("energy >> histHe", "(Z == 2 && energy > 45 && abs(posY) < 200 && abs(posZ) < 200)");
	//fragments->Scan("posY:posZ","(Z == 5 && energy > 45 && sqrt(posY^2 + posZ^2) < " + rMaxString + "&& sqrt(posY*posY + posZ*posZ) > " + rMinString + ")");
//   std::cout << "rmin and rmax are:" << rMinString << "/" << rMaxString << "\n";
	//std::cout << "energy >> histB", "(Z == 5 && energy > 45 && Sqrt(posY*posY + posZ*posZ) < " + rMaxString + "&& Sqrt(posY*posY + posZ*posZ) > " + rMinString + ")" + normalization;
    ///fragments->Draw("energy >> histB", "(Z == 5 && energy > 45 && sqrt(posY*posY + posZ*posZ) < " + rMaxString + "&& sqrt(posY*posY + posZ*posZ) > " + rMinString + ")");
    ///fragments->Draw("energy >> histH", "(Z == 1 && energy > 45 && sqrt(posY*posY + posZ*posZ) < " + rMaxString + "&& sqrt(posY*posY + posZ*posZ) > " + rMinString + ")" + normalization, "same");
   //fragments->Draw("energy >> histLi", "(Z == 3 && energy > 45 && abs(posY) < 200 && abs(posZ) < 200)" + normalization, "same");
   //fragments->Draw("energy >> histBe", "(Z == 4 && energy > 45 && abs(posY) < 200 && abs(posZ) < 200)" + normalization, "same");
   //fragments->Draw("energy >> histB", "(Z == 5 && energy > 45 && abs(posY) < 200 && abs(posZ) < 200)" + normalization, "same");
//   fragments->Draw("energy >> histC", "(Z == 6 && energy > 45 && abs(posY) < 200 && abs(posZ) < 200)" + normalization, "same");

/*
   //TCanvas *c3 = new TCanvas("histograms", "Histograms");
   Int_t nEve = events;



   cout <<"He : " << histHe->GetEntries() / nEve << endl;
   histHe->Draw("");


   cout <<"H : " << histH->GetEntries() / nEve << endl;
   histH->Draw("same");


    
   cout <<"Li : " << histLi->GetEntries() / nEve << endl;
   histLi->Draw("same");

   
   cout <<"Be : " << histBe->GetEntries() / nEve << endl;
   histBe->Draw("same");
   
   cout <<"B : " << histB->GetEntries() / nEve << endl;
   histB->Draw("same");


   ntuple->SetMarkerStyle(22);
   ntuple->Draw("H:Energy","","p,same");
   ntuple->SetMarkerColor(kRed);
   ntuple->Draw("He:Energy","","p,same");
   ntuple->SetMarkerColor(kBlue);
   ntuple->Draw("Li:Energy","","p,same");
   ntuple->SetMarkerColor(kGreen);
   ntuple->Draw("Be:Energy","","p,same");
   ntuple->SetMarkerColor(kYellow);
   ntuple->Draw("B:Energy","","p,same");

   c3->SaveAs("fig520.png");
*/
   in.close();

   f->Write();
}
