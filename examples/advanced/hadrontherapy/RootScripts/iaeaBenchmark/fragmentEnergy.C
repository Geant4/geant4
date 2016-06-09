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
////////////////////////////////////////
////// Importing data /////////
////////////////////////////////////////
   TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
   dir.ReplaceAll("fragmentEnergy.C","");
   dir.ReplaceAll("/./","/");

gStyle->SetOptStat(0000000000); //remove the for this graphs totally redundant statbox


TString macroPath(gROOT->GetMacroPath());
gROOT->SetMacroPath(macroPath + ":RootScripts/iaeaBenchmark");
//gROOT->LoadMacro("rootlogon.C");
//gROOT->SetStyle("clearRetro"); //For stylesheet

   TString pDepth;
   cout << "Enter phantom depth (eg. 27.9, see experimentalData directory for choices): " << endl;
   cout << "Entering 27.9 will make script look for IAEA_27.9.root:";
   cin >> pDepth;

   TString simulationDataPath = "IAEA_" + pDepth + ".root";
   TNtuple *ntuple = new TNtuple("ntuple","Data from ascii file","Energy:He:B:H:Li:Be");
	if(pDepth == "27.9"){
	//there is only hand.read data for this depth, fixme, is ugly
	   ifstream in;
	   in.open(Form("experimentalData/iaeaBenchmark/fragmentEnergySpctra279mmWater0deg.dat",dir.Data()));
	   Float_t f1,f2,f3, f4,f5,f6;
	   Int_t nlines = 0;
	   TFile *f = new TFile("fragmentEnergy.root","RECREATE");
	 
	   Char_t DATAFLAG[4];
	   Int_t NDATA;
	   Char_t n1[6], n2[2], n3[2], n4[2], n5[2], n6[2];
	   in >> DATAFLAG >> NDATA ; // Read EXFOR line: 'DATA 6'
	   in >> n1 >> n2 >> n3 >> n4 >> n5 >> n6; // Read column titles: 'Energy He B [...]'
	 
	   cout <<n1<<" "<<n2<<" "<<n3<<" "<<n4<<" "<<n5<<" "<<n6<<"\n";
	   while (1) {
		  in >> f1 >> f2 >> f3 >>f4 >> f5 >> f6;
		  if (!in.good()) break;
		  if (nlines < 500 ) printf("%f %0.2f %0.2f %0.2f %0.2f %0.2f \n",f1,f2,f3,f4,f5,f6);
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
	}

   //Let's pull in the monte carlo simulation results
   TCanvas *mc = new TCanvas("mc", "Simulation");
   TFile *MCData = TFile::Open(simulationDataPath);
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

 
TString halfSideLengthString(Form("%f", detectorSideLength/2));
 
   
   fragments->SetLineColor(kGreen);
   fragments->Project("histHe", "energy", "(Z == 2 && energy > 45 && abs(posY) < " + halfSideLengthString + " && abs(posZ) < " + halfSideLengthString + " )" + normalization);
   fragments->Project("histB", "energy", "(Z == 5 && energy > 45 && abs(posY) < " + halfSideLengthString + " && abs(posZ) < " + halfSideLengthString + " )" + normalization);
   fragments->Project("histH", "energy", "(Z == 1 && energy > 45 && abs(posY) < " + halfSideLengthString + " && abs(posZ) < " + halfSideLengthString + " )" + normalization);
   fragments->Project("histLi", "energy", "(Z == 3 && energy > 45 && abs(posY) < " + halfSideLengthString + " && abs(posZ) < " + halfSideLengthString + " )" + normalization);
   fragments->Project("histBe", "energy", "(Z == 4 && energy > 45 && abs(posY) < " + halfSideLengthString + " && abs(posZ) < " + halfSideLengthString + " )" + normalization);
   fragments->Project("histB", "energy", "(Z == 5 && energy > 45 && abs(posY) < " + halfSideLengthString + " && abs(posZ) < " + halfSideLengthString + " )" + normalization);
 
   histH->SetMaximum(0.27);

 
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
 
    histH->SetXTitle("Energy per nucleon (MeV/u)");
    histH->SetYTitle("(N/N0) [1/(MeV*sr)]");
  
   TCanvas *c3 = new TCanvas("histograms", "Energy distribution of charged fragments");
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
  
   // Legends for the data
   leg = new TLegend(0.8,0.5,1,1);  //coordinates are fractions
   leg->SetHeader("Fragments");
   leg->AddEntry(histH,"Hydrogen","l");
   leg->AddEntry(histHe,"Helium","l");
   leg->AddEntry(histLi,"Lithium","l");
   leg->AddEntry(histBe,"Beryllium","l");
   leg->AddEntry(histB,"Boron","l");
   leg->Draw();
 
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
 
   c3->SaveAs("fragmentEnergyDistr.png");
 
   in.close();
 
   f->Write();
}
 
