#include "Riostream.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TApplication.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TCut.h"
#include "TString.h"
#include "TMath.h"

/**
 * Macro for plotting the fragment energy distributions.
 *
 * Usage:
 * root -l RootScripts/iaeaBenchmark/fragmentEnergy.C++
 */
void fragmentEnergy() {

//   gROOT->SetStyle("clearRetro"); //For stylesheet
////////////////////////////////////////
//////     Importing data      /////////
////////////////////////////////////////
   TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
   dir.ReplaceAll("basic.C","");
   dir.ReplaceAll("/./","/");
   ifstream in;
   in.open(Form("experimentalData/iaeaBenchmark/fragmentEnergySpctra279mmWater0deg.dat",dir.Data()));
   Float_t f1,f2,f3, f4,f5,f6;
   Int_t nlines = 0;
   TFile *f = new TFile("fragmentEnergy.root","RECREATE");
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
   ntuple->SetMarkerStyle(5);
   ntuple->Draw("He:Energy","","l");
   ntuple->Draw("B:Energy","","l,Same");
   ntuple->Draw("H:Energy","","l,Same");
   ntuple->Draw("Li:Energy","","l,Same");
   ntuple->Draw("Be:Energy","","l,Same");
   printf(" found %d points\n",nlines);
   
   //Let's pull in the monte carlo simulation results
   TCanvas *mc = new TCanvas("mc", "Simulation");
   TFile *MCData = TFile::Open("IAEA.root");
   TH1F* MC_helium = (TH1F*)MCData->Get("heliumEnergyAfterPhantom");
   TH1F* MC_hydrogen = (TH1F*)MCData->Get("hydrogenEnergyAfterPhantom");
	//scale and plot
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

	//analysis numbers based on metadata
	Double_t scatteringDistance = detectorDistance - phantomCenterDistance;
    Double_t detectorSideLength = 4; //hardcoded, we have a zero angle square detector
	//good to keep for ref. G4 might give weird units due to change.
	metadata->Scan();

   
   Double_t binAmount = 50.0; //casting from int failed somehow, so in float temporarily, fixme
   Double_t maxEnergy = 450.0;
   Double_t binWidth = maxEnergy / binAmount;
   TH1F *histH = new TH1F("histH", "Hydrogen", binAmount, 0.0, maxEnergy);
   //HistH will be black
   TH1F *histHe = new TH1F("histHe", "Helium", binAmount, 0.0, maxEnergy);
   histHe->SetLineColor(kRed);
   TH1F *histLi = new TH1F("histLi", "Lithium", binAmount, 0.0, maxEnergy);
   histLi->SetLineColor(kBlue);
   TH1F *histBe = new TH1F("histBe", "Beryllium", binAmount, 0.0, maxEnergy);
   histBe->SetLineColor(kGreen);
   TH1F *histB = new TH1F("histB", "Boron", binAmount, 0.0, maxEnergy);
   histB->SetLineColor(kYellow);
   TH1F *histC = new TH1F("histC", "Carbon", binAmount, 0.0, maxEnergy);

	TH1F* histPos = new TH1F("histPos", "check position",100,-2000,2000);
	//Solid angle according to \Omega = 4 \arcsin \frac {\alpha\beta} {\sqrt{(4d^2+\alpha^2)(4d^2+\beta^2)}}
	Double_t steradians = 4 * TMath::ASin(pow(detectorSideLength,2.0) / (4*pow(scatteringDistance,2) + pow(detectorSideLength,2)) );
   std::cout << "Detector seen at solid angle: " << steradians << endl;
   TString normalization(Form("/%f", steradians*events*binWidth));

   fragments->SetLineColor(kRed);
   fragments->SetMarkerStyle(22);
   
   fragments->Draw("posY:posZ", "abs(posZ) < 2 && abs(posY) < 2");

	TString halfSideLengthString(Form("%f", detectorSideLength/2));

   
   fragments->SetLineColor(kGreen);
   fragments->Project("histHe", "energy", "(Z == 2 && energy > 45 && abs(posY) < " + halfSideLengthString + " && abs(posZ) < " + halfSideLengthString + " )" + normalization);
   fragments->Project("histB", "energy", "(Z == 5 && energy > 45 && abs(posY) < " + halfSideLengthString + " && abs(posZ) < " + halfSideLengthString + " )" + normalization);
   fragments->Project("histH", "energy", "(Z == 1 && energy > 45 && abs(posY) < " + halfSideLengthString + " && abs(posZ) < " + halfSideLengthString + " )" + normalization);
   fragments->Project("histLi", "energy", "(Z == 3 && energy > 45 && abs(posY) < " + halfSideLengthString + " && abs(posZ) < " + halfSideLengthString + " )" + normalization);
   fragments->Project("histBe", "energy", "(Z == 4 && energy > 45 && abs(posY) < " + halfSideLengthString + " && abs(posZ) < " + halfSideLengthString + " )" + normalization);
   fragments->Project("histB", "energy", "(Z == 5 && energy > 45 && abs(posY) < " + halfSideLengthString + " && abs(posZ) < " + halfSideLengthString + " )" + normalization);

   histH->Draw("");
   //TAxis *Yaxis = histH->GetYaxis();
   //Yaxis->set
   histH->SetMaximum(0.3);
   histHe->Draw("same");
   histLi->Draw("same");
   histBe->Draw("same");
   histB->Draw("same");

	/*
	 * This is cludgy but the previous plots can't contain <45MeV stuff but hte following calcualtions need them
	 * Thus there is another projection done with different selections, this is error-prone, 
	 * se to that changes are made in both places.
	 * */
   fragments->Project("histHe", "energy", "(Z == 2 && energy > 0 && abs(posY) < " + halfSideLengthString + " && abs(posZ) < " + halfSideLengthString + " )" + normalization);
   fragments->Project("histB", "energy", "(Z == 5 && energy > 0 && abs(posY) < " + halfSideLengthString + " && abs(posZ) < " + halfSideLengthString + " )" + normalization, "same");
   fragments->Project("histH", "energy", "(Z == 1 && energy > 0 && abs(posY) < " + halfSideLengthString + " && abs(posZ) < " + halfSideLengthString + " )" + normalization, "same");
   fragments->Project("histLi", "energy", "(Z == 3 && energy > 0 && abs(posY) < " + halfSideLengthString + " && abs(posZ) < " + halfSideLengthString + " )" + normalization, "same");
   fragments->Project("histBe", "energy", "(Z == 4 && energy > 0 && abs(posY) < " + halfSideLengthString + " && abs(posZ) < " + halfSideLengthString + " )" + normalization, "same");
   fragments->Project("histB", "energy", "(Z == 5 && energy > 0 && abs(posY) < " + halfSideLengthString + " && abs(posZ) < " + halfSideLengthString + " )" + normalization, "same");


	
   TCanvas *c3 = new TCanvas("histograms", "Histograms");
   Int_t nEve = events; //redundant

   cout <<"H : " << histH->GetEntries() / nEve << endl;
   histH->Draw("");

   cout <<"He : " << histHe->GetEntries() / nEve << endl;
   histHe->Draw("same");
    
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

   in.close();

   f->Write();
}
