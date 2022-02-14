
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1I.h"
#include "TF1.h"
#include "TProfile2D.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "string.h"
#include <sstream>

#define pi 3.1415926535

using namespace std;

void scattAnalysis()
{
	gROOT->Reset();
		
	//----------------------------------input--------------------------------  	  	 	 	
  	//open a simulation result file 
    TFile* file = TFile::Open("output.root");
    
    //scattering histogram parameters
    Double_t thetalimit = 0.5;		//deg
    Double_t Nbin = 120.;
    Double_t thetaMin = 0.;			//deg	
    Double_t thetaMax = 30.;		//deg
    
    //options
    Bool_t IwantTotalScatterig = false;
    Bool_t IWantPlotComptonScattering = false;
	Bool_t IWantPlotPlotComparison = false;
	Bool_t IWantPlotProcessDistrib = false;
	
	//export
    Bool_t IwantExportScattDistrib = false;
  	Char_t ExportFileName[128];
  	sprintf(ExportFileName,"SimResults.txt");
    //-----------------------------------------------------------------------

    //definitions
	Char_t cutproc[256];
	Double_t val, angle;

	//cut definition
	if (IwantTotalScatterig) {
		sprintf(cutproc,"(processID == 1 || processID == 2)");
	} else {
		sprintf(cutproc,"processID == 1");	
	}
	
	//define a tree from the ntuple contained in the input file
	TTree* scatt = (TTree*)file->Get("scatt");
	Int_t N = scatt->GetEntries();

	//Rayleigh/Total scattering
  	TH1D* thetaDistr = new TH1D("thetaDistr", "", Nbin, thetaMin, thetaMax);
  	scatt->Project("thetaDistr", "theta", cutproc); 
  	
  	Double_t Norig = thetaDistr->Integral();
  	
	for (Int_t i=1; i<=Nbin; i++) {
   		angle = thetaDistr->GetBinCenter(i);
   		if (angle >= thetalimit) {
			val = thetaDistr->GetBinContent(i);
			thetaDistr->SetBinContent(i,val/TMath::Sin(angle*pi/180));
		}
   	}
   	
	Double_t Ndiv = thetaDistr->Integral();
	thetaDistr->Scale(Norig/Ndiv);	
	
	cout << endl << "Scattering events: " << thetaDistr->Integral() << endl << endl; 
   	
  	thetaDistr->SetLineColor(kRed); 
  	thetaDistr->SetLineWidth(2);	
  	thetaDistr->SetXTitle("#theta (deg)");
  	thetaDistr->GetXaxis()->CenterTitle();
  	thetaDistr->GetXaxis()->SetTitleOffset(1.1); 
  	thetaDistr->GetXaxis()->SetRangeUser(0., thetaMax);
  	thetaDistr->SetYTitle("Counts (a.u.)");
  	thetaDistr->GetYaxis()->CenterTitle();
  	thetaDistr->GetYaxis()->SetTitleOffset(1.5); 
  	thetaDistr->SetTitle(""); 
  	TCanvas* c1 = new TCanvas("c1","",1200,800); 	
  	c1->cd(); 	
  	gStyle->SetOptStat(kFALSE);								
  	thetaDistr->Draw("HIST"); 
  	

	//Compton scattering
	if (IWantPlotComptonScattering) {
	  	TH1D* thetaDistrCompt = new TH1D("thetaDistrCompt", "", Nbin, thetaMin, thetaMax);
	  	scatt->Project("thetaDistrCompt", "theta", "processID == 2"); 
	  	
	  	Double_t NorigC = thetaDistrCompt->Integral();
	  		
		for (Int_t i=1; i<=Nbin; i++) {
	   		angle = thetaDistrCompt->GetBinCenter(i);
	   		if (angle >= thetalimit) {
				val = thetaDistrCompt->GetBinContent(i);
				thetaDistrCompt->SetBinContent(i,val/TMath::Sin(angle*pi/180));
			}
	   	}
	   	
	   	Double_t NdivC = thetaDistrCompt->Integral();
		thetaDistrCompt->Scale(NorigC/NdivC);
	   	
	  	thetaDistrCompt->SetLineColor(kBlue); 
	  	thetaDistrCompt->SetLineWidth(2);
	
		if (IWantPlotPlotComparison) {
	  		thetaDistrCompt->Draw("HIST SAME");  
		   	TLegend* leg = new TLegend(0.7,0.75,0.85,0.85);
		   	leg->SetBorderSize(0); 
			leg->AddEntry(thetaDistr,"Rayleigh","lp");
		   	leg->AddEntry(thetaDistrCompt,"Compton","lp");
		   	leg->Draw(); 
		} else {
			thetaDistrCompt->SetXTitle("#theta (deg)");
	  		thetaDistrCompt->GetXaxis()->CenterTitle();
	  		thetaDistrCompt->GetXaxis()->SetTitleOffset(1.1); 
	  		thetaDistrCompt->GetXaxis()->SetRangeUser(0., thetaMax);
	  		thetaDistr->SetYTitle("Counts (a.u.)");
	  		thetaDistrCompt->GetYaxis()->CenterTitle();
	  		thetaDistrCompt->GetYaxis()->SetTitleOffset(1.5); 
	  		thetaDistrCompt->SetTitle(""); 
		  	TCanvas* c2 = new TCanvas("c2","",1200,800); 	
		  	c2->cd(); 	
			thetaDistrCompt->Draw("HIST SAME"); 
		}				
	  	
	}


	//process distribution (0->transport, 1->Rayleigh, 2->Compton, 3->Photoel, 
	//                      4->Pair prod, 5->Nuclear, 6->Diffraction)
	if (IWantPlotProcessDistrib) {
		TH1D* procDistr = new TH1D("procDistr", "process distribution", 70, 0, 7);
	  	scatt->Project("procDistr", "processID"); 
	  	procDistr->SetLineColor(kBlack); 
	  	procDistr->SetLineWidth(2);	
	  	procDistr->SetXTitle("process");
	  	procDistr->GetXaxis()->CenterTitle();
	  	procDistr->GetXaxis()->SetTitleOffset(1.1); 
	  	procDistr->SetYTitle("Counts");
	  	procDistr->GetYaxis()->CenterTitle();
	  	procDistr->GetYaxis()->SetTitleOffset(1.2); 
		TCanvas* c3 = new TCanvas("c3","",1200,800); 
	  	c3->cd(); 								
	  	procDistr->Draw(); 
	}
	
	
	//export scattering data
  	if (IwantExportScattDistrib) {
		ofstream f(ExportFileName); 
		if (!f) {
    		cout << "Error opening the file!";
    		return;
		}
		for (Int_t i=1; i<=Nbin; i++) {
			angle = thetaDistr->GetBinCenter(i);
			if (angle >= thetalimit) {
				val = thetaDistr->GetBinContent(i);
				f << angle << " " << val << endl;
			}
		} 
		f.close(); 
		cout << "writing successful!" << endl << endl;
	} 

}

