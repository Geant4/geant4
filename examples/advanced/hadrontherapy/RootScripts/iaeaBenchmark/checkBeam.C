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
 * Macro for checking beam FWHM
 *
 * Usage:
 * root -l RootScripts/iaeaBenchmark/fragmentEnergy.C++
 */
void checkBeam() {
////////////////////////////////////////
////// Importing data /////////
////////////////////////////////////////
   TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
   dir.ReplaceAll("fragmentEnergy.C","");
   dir.ReplaceAll("/./","/");

TString macroPath(gROOT->GetMacroPath());
gROOT->SetMacroPath(macroPath + ":RootScripts/iaeaBenchmark");
//gROOT->LoadMacro("rootlogon.C");
//gROOT->SetStyle("clearRetro"); //For stylesheet

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
   in >> n1 >> n2 >> n3 >> n4 >> n5 >> n6; // Read column titles: 'Energy He B [...]'
 
   cout <<n1<<" "<<n2<<" "<<n3<<" "<<n4<<" "<<n5<<" "<<n6<<"\n";
   while (1) {
      in >> f1 >> f2 >> f3 >>f4 >> f5 >> f6;
      if (!in.good()) break;
      if (nlines < 500 ) printf("%f %0.2f %0.2f %0.2f %0.2f %0.2f \n",f1,f2,f3,f4,f5,f6);
      ntuple->Fill(f1,f2,f3,f4,f5,f6);
      nlines++;
   }

   printf(" found %d points\n",nlines);
   
   //Let's pull in the monte carlo simulation results
   TFile *MCData = TFile::Open("IAEA_15.9.root");
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

   //TCanvas* c3 = new TCanvas();
   fragments->SetLineColor(kRed);
   fragments->SetMarkerStyle(22);
   
    //fragments->Draw("posY:posZ", "abs(posZ) < 2 && abs(posY) < 2");
	TH1F* posYHisto = new TH1F("vertical","Vertical distribution of beam",200,-2.,2.);
	TH1F* posZHisto = new TH1F("horizontal","horizontal distribution of beam",200,-2.,2.);
	TH2F* posHisto = new TH2F("posHisto","Distribution of beam",200,-2.,2.,200,-2.,2.);
	//fragments->Project("posHisto","posY:posZ","abs(posZ) < 2 && abs(posY) < 2");
	fragments->Project("vertical","posZ","abs(posZ) < 2 && abs(posY) < 2 && Z == 6");
	fragments->Project("horizontal","posY","abs(posZ) < 2 && abs(posY) < 2 && Z == 6");
	Float_t maxVal = posYHisto->GetMaximum();
	std::cout << "maximum: " << maxVal << endl;
	//fragments->Scan("posY:posZ","abs(posZ) < 2 && abs(posY) < 2");
	//posYHisto->Draw();
//	posZHisto->Draw("same");

 	Float_t fwhm = 0.0, middle = 0.0, curVal;
	int fwhmBin;
	//So these dots have been binned and a fwhm is calculated.
	for(int i = 0; i < posYHisto->GetNbinsX(); i++){
		
		curVal = posYHisto->GetBinContent(i);
		if(pow(maxVal/2 - middle, 2.0) > pow(maxVal/2 - curVal, 2.0)){
				fwhm = 2*TMath::Abs(posYHisto->GetBinCenter(i));
				middle = curVal;
				fwhmBin = i;
			}else{
				}
		}
		//posYHisto->SetBinContent(fwhmBin, 0);
		/*
		TNtuple* where = new TNtuple("where","where","x:y");
		where->Fill(fwhmBin*posYHisto->GetBinWidth(0), 0);
		where->Fill((200-fwhmBin)*posYHisto->GetBinWidth(0), 0);
		where->->SetMarkerStyle(22); //triangle
		where->SetMarkerColor(kRed);
		*/
		TF1* fitgaus = new TF1("fitgaus","gaus");
		fitgaus->SetLineColor(2);
		posYHisto->Fit(fitgaus,"");
		posYHisto->Draw();
		std::cout << "fitted FWHM from normal distribution is: " << fitgaus->GetParameter(2)*10*2.35482 << endl;
		//where->Draw("x:y","same");
		//posYHisto->Smooth(150);
		//posYHisto->Draw();
		//posHisto->Draw();
		std::cout << "Calculated (closest point) FWHM of Monte-Carlo simulation to be: " << fwhm*10 << " mm" << endl; 
		std::cout << "beam contained " << fragments->GetEntries("Z == 6") << " carbon nuclei." << endl;
		
		
   //c3->SaveAs("checkBeam.png");
 
}
 
