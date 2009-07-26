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

//   gROOT->SetStyle("clearRetro");

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
   //Let's pull in the monte carlo analysis results

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

   //fragments->Scan();

   Double_t ScaleHelium = 1/(MC_helium->Integral());
   Double_t ScaleHydrogen = 1/(MC_hydrogen->Integral()); 
   //x should also be scaled to per nucleon
   
   MC_helium->Scale(ScaleHelium);
//   MC_helium->Draw("");
   //printf("Scaled helium by %.9f\n",ScaleHelium);

   MC_hydrogen->Scale(ScaleHydrogen);
   MC_hydrogen->SetLineColor(kRed);
//   MC_hydrogen->Draw("Same");
  // printf("Scaled hydrogen by %.9f\n",ScaleHydrogen);
   
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

	TH1F* histPos = new TH1F("histPos", "check position",100,-2000,2000);
	//The steradian calculation is a bit of an approx as we have a square detector and the covering area is approximated to be the are of the square which it isn't
	Double_t steradians = (detectorSideLength * detectorSideLength) / (4*TMath::Pi()*scatteringDistance);
   TString normalization(Form("/%f", steradians*events));

   fragments->SetLineColor(kRed);
   fragments->SetMarkerStyle(22);
   
   fragments->Draw("posY:posZ", "abs(posZ) < 20 && abs(posY) < 20");

   fragments->Draw("energy >> histHe", "(Z == 2 && energy > 45 && abs(posY) < 2 && abs(posZ) < 2)" + normalization);
   
   fragments->SetLineColor(kGreen);

   fragments->Draw("energy >> histB", "(Z == 5 && energy > 45 && abs(posY) < 200 && abs(posZ) < 200)" + normalization, "same");
   fragments->Draw("energy >> histH", "(Z == 1 && energy > 45 && abs(posY) < 200 && abs(posZ) < 200)" + normalization, "same");
   fragments->Draw("energy >> histLi", "(Z == 3 && energy > 45 && abs(posY) < 200 && abs(posZ) < 200)" + normalization, "same");
   fragments->Draw("energy >> histBe", "(Z == 4 && energy > 45 && abs(posY) < 200 && abs(posZ) < 200)" + normalization, "same");
   fragments->Draw("energy >> histB", "(Z == 5 && energy > 45 && abs(posY) < 200 && abs(posZ) < 200)" + normalization, "same");
//   fragments->Draw("energy >> histC", "(Z == 6 && energy > 45 && abs(posY) < 200 && abs(posZ) < 200)" + normalization, "same");


   TCanvas *c3 = new TCanvas("histograms", "Histograms");
   Int_t nEve = events; //redundant



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

   in.close();

   f->Write();
}
