//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// --------------------------------------------------------------------
// Class implementation
// --------------------------------------------------------------------

#include "G4TTool.h"


ClassImp(G4TTool)


using namespace std;
using namespace ROOT;
using namespace TMath;


//______________________________________________________________________________
void G4TTool::Initialize()
{
	fHxmin = 3;
	fHxmax = fPublication->fHeader.fTypeValue + 3;
	fHymin = 1.0001e-6;
	fHymax = 0;
	fHfbin	=  1;
	fHnbin	=  50;
	fDangle =  7.5;
	fPi = 3.141593;
	fDanrad = fPi * fDangle / 180;
	fHlxmin = log10(fHxmin);
	fHlxmax = log10(fHxmax);

	PrepareHistograms(fHnbin, fHlxmin, fHlxmax);
}

//______________________________________________________________________________
void G4TTool::PrepareHistograms(Double_t hnbin, Double_t hlxmin, Double_t hlxmax)
{
	// Check if deleted
	fHZL = (TH1F*)gDirectory->Get("hZL");
	fHDT = (TH1F*)gDirectory->Get("hDT");

	if(fHZL != 0) delete fHZL;
	if(fHDT != 0) delete fHDT;

	// create zero level and dT
	fHZL = new TH1F("hZL", "Zero Level", hnbin, hlxmin, hlxmax);
	fHDT = new TH1F("hDT", "dT in MeV", hnbin, hlxmin, hlxmax);

	// adjust zero levels, create dT histogram
	Double_t hzero =  0.1e-12;
	Double_t halfbn = 0.5 * (hlxmax - hlxmin) / hnbin;

	//cout << "hlxmax=" <<  hlxmax << ", hlxmin=" << hlxmin << ", hnbin=" << hnbin << ", halfbn=" << halfbn << endl;

	Double_t* vZeroLevels = new Double_t[(Int_t)hnbin];
	Double_t* vErrors = new Double_t[(Int_t)hnbin];

	for(Int_t i = 0; i< hnbin; ++i){
		 vZeroLevels[i] = hzero * hnbin;
		 vErrors[i] = vZeroLevels[i] * 0.8;
	}

	// set the content and error
	fHZL->SetContent(&vZeroLevels[0]);
	fHZL->SetError(&vErrors[0]);

	// get values of center of bins abscissa into vector
	TAxis* axis = fHDT->GetXaxis();
	for(Int_t i = 0; i< hnbin; i++){
		vZeroLevels[i] = axis->GetBinCenter(i);
	}

	for(Int_t i = 0; i< hnbin; i++){
		 vZeroLevels[i] = TMath::Power(10, vZeroLevels[i] + halfbn ) - TMath::Power(10, vZeroLevels[i] - halfbn);
	}
	fHDT->SetContent(&vZeroLevels[0]);

	// clear
	delete vZeroLevels;
	delete vErrors;

}

//______________________________________________________________________________
void G4TTool::RenderHSolid(TH1F* hist, Int_t hf, Int_t hn, Double_t m, Color_t color, Bool_t noDots )
{
	if(hist)
	{
		Int_t hnbin = hn - hf;
		Double_t* vWx = new Double_t[hnbin];
		Double_t* vWy = new Double_t[hnbin];
		Double_t* vWd = new Double_t[hnbin];
		Double_t* vWt = new Double_t[hnbin];
		Double_t* vWp = new Double_t[hnbin];
		Double_t* vWf = new Double_t[hnbin];

		// get values of center of bins abscissa into vector
		TAxis* axis = hist->GetXaxis();
		if(axis == 0){
			cout << "No axis found!" << endl;
		}

		for(Int_t i = hf; i< hn; ++i){
			vWx[i - hf] = axis->GetBinCenter(i);
		}

		for(Int_t i = hf; i< hn; ++i){
			vWy[i - hf] = hist->GetBinContent(i);
		}

		for(Int_t i = hf; i< hn; ++i){
			vWd[i - hf] = hist->GetBinError(i);
		}

		for(Int_t i = 0; i< hnbin; ++i){
			vWt[i] = TMath::Power(10.0, vWx[i]);
		}

		for(Int_t i = 0; i< hnbin; ++i){
			vWp[i] = sqrt(vWt[i] * ((2 * m) + vWt[i]));
		}

		for(Int_t i = 0; i< hnbin; ++i){
			vWf[i] = vWy[i] / vWp[i];
		}

		TGraph* graph = new TGraph(hnbin,vWt, vWf);
		graph->SetLineColor(color); // red
		graph->SetLineWidth(2);
		graph->Draw();

		if(!noDots)
		{
			TGraph* dots = new TGraph(hnbin,vWt, vWf);
			dots->SetLineStyle(4); // dot dot
			dots->SetLineWidth(2);
			dots->Draw();
		}

		// clean up
		delete[] vWx;
		delete[] vWy;
		delete[] vWd;
		delete[] vWp;
		delete[] vWt;
		delete[] vWf;

	}else{
		cout << "Error: no histogram object found!" << endl;
	}
}





