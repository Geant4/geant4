// G.A.P.Cirrone (2009)

// Comparison beetwen experimental and simulated non modulated Bragg Peak
// Remember to run this file only if simulation has been run with the 
// proton_therapy.mac file

{  
gROOT->Reset();
#include "Riostream.h"

 ifstream in;


 // LOAD THE EXPERIMENTAL DATA FILE
 // CONTAINED IN THE DIRECTORY
 // hadrontherapy/experimentalData/proton/BraggPeak
 TFile *experimentalFile = new TFile("../../../experimentalData/proton/BraggPeak/62MeVInWater.root","READ");

 // HERE THE ROOT FILE IS INTERPRETED AS A TREE
 TTree *experimentalTree = (TTree*)experimentalFile -> Get("Experimental62MeVInWater");
 
 Float_t depthExp, EdepExp;
 experimentalTree -> SetBranchAddress("EdepExp", &EdepExp);
 experimentalTree -> SetBranchAddress("depthExp", &depthExp);

 
 // CREATION AND NORMALISATION TO THE FIRST POINT  OF AN NTUPLE CONTAINING THE EXPERIMENTAL DATA
 TNtuple *ntupleExperimental = new TNtuple("ntupleExperimental","Protons, exp. data", "depthExp:EdepExp");
 Int_t nentries = (Int_t)experimentalTree -> GetEntries();   
 for (Int_t i = 0; i<nentries; i++)
   {
     experimentalTree -> GetEntry(0);
     Float_t normFactor = EdepExp;
     experimentalTree -> GetEntry(i);
     ntupleExperimental -> Fill(depthExp, EdepExp/normFactor);
     
   }
 
 // LOAD THE SIMULATION RESULT FILE
 // CONTAINED IN THE DIRECTORY
 // hadrontherapy/simulationResults/proton/BraggPeak 
 TFile *simulationFile = new TFile("../../../SimulationOutputs/proton/BraggPeak/protonBraggPeak.root","READ");

 // EXTRACTION, FROM THE SIMULATION FILE OF THE INTERESTING HISTOGRAMS
 TH1D *simulatedPeak = (TH1D*) simulationFile -> Get("braggPeak");

 Float_t simulationNormalisationFactor =  simulatedPeak -> GetBinContent(1);
 simulatedPeak -> Scale(1/simulationNormalisationFactor);


 TCanvas *c1 = new TCanvas ("c1","c1",200,10,600,400);
 
 // PLOT
 ntupleExperimental -> SetMarkerStyle(4);
 simulatedPeak -> SetMarkerSize(2);

 ntupleExperimental -> Draw("EdepExp:depthExp");
 simulatedPeak-> Draw("same");
 
 // LEGEND
 leg = new TLegend(0.50,0.60,0.20,0.70); 
 leg -> SetTextSize(0.035);
 leg -> SetFillColor(0);
 leg -> AddEntry(ntupleExperimental, "Experiment","P");
 leg -> AddEntry(simulatedPeak, "Simulation");
 leg -> Draw();




};
