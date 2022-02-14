
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

void ADXRD()
{	
	gROOT->Reset();

	//----------------------------------input--------------------------------  	  	 	 	
  	//open a simulation result file 
	TFile* file = TFile::Open("output.root");
    
  	//histogram parameters
    Double_t binXsize = 0.5; 		//mm
  	Double_t binXStart = -200.;		//mm
  	Double_t binXEnd = -binXStart; 
  	Int_t nbinsX = binXEnd*2./binXsize; 
  	
  	Double_t binYsize = 0.5; 		//mm   
  	Double_t binYStart = -200.;		//mm
  	Double_t binYEnd = binXEnd; 
  	Int_t nbinsY = binYEnd*2./binYsize;  
  	
    Double_t Nbin = 92;
    Double_t thetaMin = 1;			//deg
    Double_t thetaMax = 27;			//deg
  	
  	//cuts	
  	Double_t Emin = 0.;			    //keV
  	Double_t Emax = 140.;			//keV
  	Double_t angCut = 0.; 			//deg  	
  	Bool_t IWantOnlyScatt = true;		
  	Bool_t IWantOnlyRayleighScatt = false;
  	Bool_t IWantOnlyComptonScatt = false;
	Bool_t IWantOnlyDiffraction = false;
  	
  	//plot options
  	Bool_t IWantThetaDistribution = true;
  	Bool_t IWantEdistrib = true;  	
  	
  	Bool_t IWantPlot2D = true; 	
  	Bool_t IWantProfiles = false;
  	
  	Bool_t IWantBoxAnalysis = false;
  	Bool_t IWantPlotEtheta = true;  	
  	
  	Int_t nbp = 2;
  	gStyle->SetOptStat(kFALSE);	
  	//gStyle->SetPalette(52);   //52->gray, 53->hot
  	
  	//Energy Distribution Analysis
  	Double_t NbinE = 140;
    Double_t EMin = 0;				//keV
    Double_t EMax = 140;			//keV        
    Double_t Angle = 2.58;			//deg 
    Double_t DeltaAngle = 2.;		//deg   
    Bool_t ApplyThetaCut = false;
  	
  	//export
  	Bool_t IwantExportScattering = false;
  	Char_t ExportScattFileName[128];
  	sprintf(ExportScattFileName,"scatt.txt");
  	
  	Bool_t IwantExportVariables = false;
  	Char_t ExportVarFileName[128];
  	sprintf(ExportVarFileName,"var.txt");
  	
  	Bool_t IwantExportImage = false;
  	Char_t ExportImageFileName[128];
  	sprintf(ExportImageFileName,"image.txt");
    //-----------------------------------------------------------------------
  
  	//tree definition  
  	Double_t x,y,energy;
  	Int_t kind,ID,nri,nci,ndi;
	TTree* t1 = (TTree*)file->Get("part");
	t1->SetBranchAddress("e",&energy);
	t1->SetBranchAddress("posx",&x);
	t1->SetBranchAddress("posy",&y);
	t1->SetBranchAddress("type",&kind);
	t1->SetBranchAddress("trackID",&ID);
	t1->SetBranchAddress("NRi",&nri);
	t1->SetBranchAddress("NCi",&nci);
	t1->SetBranchAddress("NDi",&ndi);
	Int_t N = t1->GetEntries();
	const Int_t Nval = N;
    
    //scatt cut
    Char_t cutScatt[128];
    if (IWantOnlyScatt) {
    	if (IWantOnlyRayleighScatt) {
    		sprintf(cutScatt,"type==0 & trackID==1 & e>=%f & e<=%f & NRi>=1 & NCi==0 & NDi==0",Emin,Emax);
    	} else if (IWantOnlyComptonScatt) {
    		sprintf(cutScatt,"type==0 & trackID==1 & e>=%f & e<=%f & NRi==0 & NCi>=1 & NDi==0",Emin,Emax);
    	} else if (IWantOnlyDiffraction) {
    		sprintf(cutScatt,"type==0 & trackID==1 & e>=%f & e<=%f & NRi==0 & NCi==0 & NDi>1",Emin,Emax);
    	} else {
    		sprintf(cutScatt,"type==0 & trackID==1 & e>=%f & e<=%f & (NRi+NCi+NDi)>=1",Emin,Emax);
    	}
    } else {
    	sprintf(cutScatt,"type==0 & trackID==1 & e>=%f & e<=%f",Emin,Emax);
    }


    //spatial distribution
	TH2D* SpatialDistribution = new TH2D("SpatialDistribution", "Spatial Distribution", nbinsX, binXStart, binXEnd, nbinsY, binYStart, binYEnd); 
	t1->Project("SpatialDistribution","posy:posx",cutScatt); 
	if (IWantPlot2D) {
		SpatialDistribution->SetXTitle("X (mm)");
		SpatialDistribution->GetXaxis()->CenterTitle();
		SpatialDistribution->GetXaxis()->SetTitleOffset(1.2); 
		SpatialDistribution->SetYTitle("Y (mm)");
		SpatialDistribution->GetYaxis()->CenterTitle();
		SpatialDistribution->GetYaxis()->SetTitleOffset(1.3); 
		TCanvas* c1 = new TCanvas("c1","",1000,1100); 
		c1->cd();
		SpatialDistribution->SetContour(50);
		SpatialDistribution->Draw("colz");
	}
	
	//profiles (projections transformed in profiles by myself)
  	if (IWantProfiles) {
  		TCanvas* c2 = new TCanvas("c2","",1200,900); 
  		c2->Divide(2,1); 
  		
  		Int_t bcx = Int_t(nbinsX/2.);
  		Int_t bcy = Int_t(nbinsY/2.);
  		Double_t sf = 1./(2.*nbp+1.);
  		  		
  		TH1D* profileX = SpatialDistribution->ProjectionX("px", bcx-nbp, bcx+nbp);
  		profileX->Scale(sf);
  		profileX->SetTitle("X profile");
  		c2->cd(1);
  		//profileX->Fit("gaus");
  		profileX->Draw(); 
  		
  		TH1D* profileY = SpatialDistribution->ProjectionY("py", bcy-nbp, bcy+nbp);
  		profileY->Scale(sf);
  		profileY->SetTitle("Y profile");
  		c2->cd(2);
  		//profileY->Fit("gaus");
  		profileY->Draw();
  	
  		cout << endl;
  		cout << "profileX FWHM: " << profileX->GetRMS()*2.35 << " mm" << endl;
  		cout << "profileY FWHM: " << profileY->GetRMS()*2.35 << " mm" << endl;	
  		cout << endl;
  	}  

  	
	//theta distribution
  	TH1D* thetaDistr = new TH1D("thetaDistr", "", Nbin, thetaMin, thetaMax);
  	t1->Project("thetaDistr","acos(momz)*180/acos(-1)",cutScatt); 
  	
  	Double_t Norig = thetaDistr->Integral();  	
  	  	
	Double_t val, angle;
	for (Int_t i=1; i<=Nbin; i++) {
   		angle = thetaDistr->GetBinCenter(i)*acos(-1)/180;
		val = thetaDistr->GetBinContent(i);
		thetaDistr->SetBinContent(i,val/TMath::Sin(angle));
   	}	
   	
   	Double_t Ndiv = thetaDistr->Integral();
	thetaDistr->Scale(Norig/Ndiv);
	
	cout << endl << "total counts of thetaDistr: " << thetaDistr->Integral() << endl;
   	
	if (IWantThetaDistribution) {	   	
	  	thetaDistr->SetLineColor(kBlack); 
	  	thetaDistr->SetLineWidth(2);	
	  	thetaDistr->SetXTitle("#theta (degree)");
	  	thetaDistr->GetXaxis()->CenterTitle();
	  	thetaDistr->GetXaxis()->SetTitleOffset(1.1); 
	  	//thetaDistr->GetXaxis()->SetRangeUser(2., 12.);
	  	thetaDistr->SetYTitle("Counts");
	  	thetaDistr->GetYaxis()->CenterTitle();
	  	thetaDistr->GetYaxis()->SetTitleOffset(1.2); 
	  	TCanvas* c3 = new TCanvas("c3","",1100,1100); 
	  	c3->cd(); 								
	  	thetaDistr->Draw("HIST"); 
  	}

  	
  	//Energy Distribution Analysis
	if (IWantEdistrib) {
		//Theta cut
		Char_t cutTheta[256];
			
		if (ApplyThetaCut) {
			sprintf(cutTheta,"type==0 & TMath::Abs(acos(momz)-%f) <= %f",Angle*0.017453293,DeltaAngle*0.017453293);
		} else {
			sprintf(cutTheta,"");
		}
		  
		//Energy distribution of scattered photons
		TH1D* edistrib = new TH1D("edistrib", "", NbinE, EMin, EMax);
		t1->Project("edistrib", "e", cutTheta);
	  	edistrib->SetLineColor(kRed); 
	  	edistrib->SetLineWidth(2);	
	  	edistrib->SetXTitle("E (keV)");
	  	edistrib->GetXaxis()->CenterTitle();
	  	edistrib->GetXaxis()->SetTitleOffset(1.1); 
	  	edistrib->SetYTitle("Counts (a. u.)");
	  	edistrib->GetYaxis()->CenterTitle();
	  	edistrib->GetYaxis()->SetTitleOffset(1.3); 
	  	edistrib->SetTitle("Energy spectrum of the particles impinging on the detector"); 
	  	TCanvas* c4 = new TCanvas("c4","",1200,800); 	
	  	c4->cd(); 	
	  	gStyle->SetOptStat(kFALSE);								
	  	edistrib->Draw(); 
	  	
	  	Int_t Ntot = edistrib->Integral();
	  	cout << "total counts of edistrib: " << Ntot << endl;
	}

	
	//Box Score Analysis
	if (IWantBoxAnalysis) {
		Int_t xBins = 40; 	
	  	Int_t yBins = 40; 
	  	Int_t zBins = 20; 	
	
		Double_t xlow = -25; 			//mm
	  	Double_t xup = 25; 
	  	Double_t ylow = -25; 
	  	Double_t yup = 25; 
	  	Double_t zlow = -10; 
	  	Double_t zup = 10;  
	  	
	  	Double_t rngadd = 1.; 			//Box plot range margin (mm) 
	
		TH3D* XYZ = new TH3D("XYZ","Hot spots", xBins, xlow-rngadd, xup+rngadd, zBins, zlow-rngadd, zup+rngadd, yBins, ylow-rngadd, yup+rngadd);
		t1->Project("XYZ", "posy:posz:posx",cutScatt);
		XYZ->SetFillColor(2);
		XYZ->SetXTitle("x (mm)");
		XYZ->GetXaxis()->CenterTitle();
		XYZ->GetXaxis()->SetTitleOffset(1.8); 
		XYZ->SetYTitle("z (mm)");
		XYZ->GetYaxis()->CenterTitle();
		XYZ->GetYaxis()->SetTitleOffset(2.4);  	
		XYZ->SetZTitle("y (mm)");
		XYZ->GetZaxis()->CenterTitle();
		XYZ->GetZaxis()->SetTitleOffset(1.8); 
		gStyle->SetCanvasPreferGL(kTRUE);
		TCanvas* c5 = new TCanvas("c5","Box Score Analysis",1200,800); 
		c5->cd();
		XYZ->Draw("glbox");
	}

	
	//energy-angle correlation of impinging photons
	if (IWantPlotEtheta) {
		TH2D* Etheta = new TH2D("Etheta", "", NbinE, EMin, EMax, Nbin, thetaMin, thetaMax); 
		t1->Project("Etheta","acos(momz)*180/acos(-1):e"); 
		Etheta->SetXTitle("E (keV)");
		Etheta->GetXaxis()->CenterTitle();
		Etheta->GetXaxis()->SetTitleOffset(1.2); 
		Etheta->SetYTitle("#theta (degree)");
		Etheta->GetYaxis()->CenterTitle();
		Etheta->GetYaxis()->SetTitleOffset(1.3); 
		TCanvas* c6 = new TCanvas("c6","",1000,1100); 
		c6->cd();
		Etheta->SetContour(50);
		Etheta->Draw("colz");
	}

  	 	
	//export scattering data	
  	if (IwantExportScattering) {
  		//open a txt file
  		ofstream f(ExportScattFileName); 
    	if(!f) {
        	cout << "Error opening the file!";
        	return;
    	}
    	//variables extraction from tree
    	for (Int_t i=1; i<=Nbin; i++) {
    		Double_t ang = thetaDistr->GetBinCenter(i);
    		Double_t scatt = thetaDistr->GetBinContent(i);
    		if (ang >= angCut) {
				f << ang << " " << scatt << endl;
			}
    	} 
  		//close the txt file
    	f.close(); 
   		cout << "writing the file with the scattering data successful!" << endl << endl;
   	}	

  	
  	//export (x,y) variables	
  	if (IwantExportVariables) {
  		//open a txt file
  		ofstream f(ExportVarFileName); 
    	if(!f) {
        	cout << "Error opening the file!";
        	return;
    	}
    	for (Int_t i=0; i<Nval; i++) {
    		t1->GetEntry(i); 
			if (IWantOnlyScatt) {
				if (kind==0 && ID==1 && energy>=Emin && energy<=Emax && nri+nci+ndi>=1) {
					f << x << " " << y << endl;	
				}
			} else {
				f << x << " " << y << endl;
			}
    	} 
  		//close the txt file
    	f.close(); 
   		cout << "writing the file with the (x,y) variables successful!" << endl << endl;
   	}	

   	
   	//export Image	
  	if (IwantExportImage) {
  		//open a txt file
  		ofstream fout(ExportImageFileName); 
    	if (!fout) {
        	cout << "Error opening the file!";
        	return;
    	}
    	//variables extraction from histogram
    	Double_t counts = 0.;
    	for (Int_t i=1; i<=nbinsY; i++) {
    		for (Int_t j=1; j<=nbinsX; j++) {
    			counts = SpatialDistribution->GetBinContent(j,i);
    			fout << counts << " ";
    		}
    		fout << endl;
    	} 
  		//close the txt file
    	fout.close(); 
   		cout << "writing text image file successful!" << endl << endl;
   	}	

}

