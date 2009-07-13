#include "Riostream.h"

/**
 * Example script for reading datafile.
 * Usage:
 * root -l RootScripts/iaeaBenchmark/readExfor.C
 */
void readExfor() {
   TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
   dir.ReplaceAll("basic.C","");
   dir.ReplaceAll("/./","/");
   ifstream in;
   //in.open(Form("../data/fragmentEnergySpctra279mmWater0deg.dat",dir.Data()));
   in.open(Form("experimentalData/fragmentEnergySpctra279mmWater0deg.dat",dir.Data()));
   Float_t f1,f2,f3, f4,f5,f6;
   Int_t nlines = 0;
   TFile *f = new TFile("basic.root","RECREATE");
   TH1F *h1 = new TH1F("h1","x distribution",100,-4,4);
   TNtuple *ntuple = new TNtuple("ntuple","data from ascii file","Energy:He:B:H:Li:Be");
	  
   Char_t DATAFLAG[4];
   Int_t NDATA;
   Char_t n1[6], n2[2], n3[2], n4[2], n5[2], n6[2],;
   in >> DATAFLAG >> NDATA ; // Read EXFOR line: 'DATA 6'
   in >> n1 >> n2 >> n3 >> n4 >> n5 >> n6; // Read  column titles: 'Energy He B [...]'

   cout <<n1<<" "<<n2<<" "<<n3<<"    "<<n4<<"    "<<n5<<"   "<<n6<<"\n";
   while (1) {
      in >> f1 >> f2 >> f3 >>f4 >> f5 >> f6;
      if (!in.good()) break;
      if (nlines < 500 ) printf("%i  %0.2f %0.2f %0.2f %0.2f %0.2f \n",f1,f2,f3,f4,f5,f6);
      //h1->Fill(f1);
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

   TCanvas *mc = new TCanvas("c2", "Simulation");
   TFile *MCData = TFile::Open("IAEA.root");
   TH1F* MC_helium = (TH1F*)MCData->Get("heliumEnergyAfterPhantom");
   TH1F* MC_hydrogen = (TH1F*)MCData->Get("hydrogenEnergyAfterPhantom");
		//scale and plot
   TNtuple *fragments = (TNtuple*) MCData->Get("fragmentNtuple");

   ScaleHelium = 1/(MC_helium->Integral());
   ScaleHydrogen = 1/(MC_hydrogen->Integral()); 
   //x should also be scaled to per nucleon
   
   MC_helium->Scale(ScaleHelium);
   MC_helium->Draw("");
   printf("Scaled helium by %.9f\n",ScaleHelium);

   MC_hydrogen->Scale(ScaleHydrogen);
   MC_hydrogen->SetLineColor(kRed);
   MC_hydrogen->Draw("Same");
   printf("Scaled hydrogen by %.9f\n",ScaleHydrogen);
   
   TCanvas *fc = new TCanvas("fc", "Fragments");
   fragments->Draw("energy");

   in.close();

   f->Write();
}
