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

/***************************
 * Root script that produces comparison between hadrontherapy
 * mc-simulation data and E.Haettner's experimental data.
 * 
 * prompts user if to normalize to zero depth=1
 * 
 * @author Gillis Danielsen
 * **************************/


void braggPeak() {
	TString normToZeroPos;
	cout << "Normalize to first bin? (Y/N):";
    cin >> normToZeroPos;
	
	TCanvas *c1 = new TCanvas("Bragg curve", "Bragg curve comparison");

	//TCanvas *c1 = new TCanvas("BraggPeaks", "Energy depositions along x-axis in phantom");
	gStyle->SetOptStat(0000000000); //remove the for this graphs totally redundant statbox

   gROOT->SetStyle("clearRetro");
 //this will be used as base for pulling the experimental data
   TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
   dir.ReplaceAll("basic.C","");
   dir.ReplaceAll("/./","/");
   ifstream in;
   in.open(Form("experimentalData/iaeaBenchmark/400tdk.dat",dir.Data())); //also set Znum!
   TString Znum = "2"; //set here what will be put in the selection for z, does not automatically change imoprted data.
   Float_t f1,f2,f3,f4;
   Int_t nlines = 0;
   TFile *f = new TFile("fragmentAngularDistribution.root","RECREATE");
   TNtuple *ntuple = new TNtuple("ntuple","Data from ascii file","d:i:err1:err2");
	  
   Char_t DATAFLAG[4];
   Int_t NDATA;
   Char_t n1[15], n2[15], n3[15], n4[15];
   in >> DATAFLAG >> NDATA ; // Read EXFOR line: 'DATA 4'
   in >> n1 >> n2 >> n3 >> n4; // Read  column titles: 'Energy He B [...]'

   cout <<n1<<"   "<<n2<<"   "<<n3<<"   "<<n4<<"\n";
   while (1) {
      in >> f1 >> f2 >> f3 >> f4;
      if (!in.good()){ 
		  break;
	  }
      if (nlines < 500 ) printf("%f %f %f %f\n",f1,f2,f3,f4);
      ntuple->Fill(f1,f2,f3,f4);
      nlines++;
   }   
   
   //Let's pull in the simulation data
   TFile* simulation = TFile::Open("IAEA_braggPeak.root");
   TH1F* simBragg = (TH1F*) simulation->Get("braggPeak");

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

	simBragg->GetXaxis()->SetLimits(0.0, waterThickness);
	simBragg->SetLineColor(kBlue);
	simBragg->SetXTitle("Depth in phantom (cm)");
	//simBragg->GetXaxis()->SetTitleOffset(1);
	simBragg->GetYaxis()->SetTitleOffset(1.5);
	std::cout << "Maximum (Bragg peak) for simulation data is at: " << simBragg->GetBinCenter(simBragg->GetMaximumBin()) + 0.478/2 + 0.027 + 0.073 << endl;
	std::cout << "Bin width is " << simBragg->GetBinWidth(simBragg->GetMaximumBin()) << endl;
	
	TString scaleTuple;
	if(normToZeroPos == "Y"){
	Float_t normElement;
	ntuple->SetBranchAddress("i",&normElement);
	ntuple->GetEntry(0);
	scaleTuple = Form("/%f", normElement);
	simBragg->Scale(1.0/simBragg->GetBinContent(0));
	simBragg->SetYTitle("Relative ionization");
	}else{
	simBragg->Scale(1.0/(100*events*simBragg->GetBinWidth(0))); // 100 for converting to MeV/m
	simBragg->SetYTitle("Ionization (MeV/m)");
	scaleTuple = "";
	}

	//simBragg->ShowPeaks(); //Can be used to mark out bragg peaks specificly.
	simBragg->Draw();
	ntuple->SetMarkerColor(kRed);
	ntuple->SetMarkerStyle(22);
	std::cout << ntuple->GetEntries() << endl;
	ntuple->Draw("i" + scaleTuple + ":d-(0.478+0.027+0.073)","","p,same"); // the minuses are the WE's, here only real water depth is plotted.
	c1->SaveAs("braggPeakComparisonToData_norm_" + normToZeroPos + ".png");
}
