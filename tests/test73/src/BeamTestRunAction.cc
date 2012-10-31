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
#include "BeamTestRunAction.hh"
#include "BeamTestRun.hh"
#include "BeamTestEventAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include <assert.h>
#ifdef G4ANALYSIS_USEROOT
// Allow output of the std deviation
#include "TStyle.h"
#include "TChain.h"
#include "TApplication.h"
#endif
//#include "TPaveStats.h"
//#include "Fit/FitResult.h"
//#include "Fit/SparseData.h"

BeamTestRunAction::BeamTestRunAction() 
: mergedFileName("merged"), pT(0)
#ifdef G4ANALYSIS_USEROOT
, f(0), pTString(0),ipPerTrackb_x(0),ipPerTrackb_y(0), ipPerTrackb_z(0), ipPerTrackb_mag(0), EDepPerTrack(0), sigma(0)
, msigma_x(0), msigma_y(0), msigma_z(0), edep(0), fitresultx(0), fitresulty(0), fitresultz(0)
#endif
,messager(new BeamTestRunActionMessenger(this))
{}

void BeamTestRunAction::Initialize() {
#ifdef G4ANALYSIS_USEROOT
    gROOT->Reset();
    //gStyle->SetOptStat(11111111111);

    // Set pT to sting for naming of files
    //pT = parameter->pTransverse;
    pTString = fNameSuffix.str();
    //gStyle->SetOptStat(111111111111);
    TString Titlebx = "b_xPerTrack"+pTString;
    TString Titleby = "b_yPerTrack"+pTString;
    TString Titlebz = "b_zPerTrack"+pTString;
    TString Titlebmag = "b_magPerTrack"+pTString;
    TString Titleedep = "EDepPerTrack"+pTString;

    TString Name = "IPbinpT"+pTString+".root";
    f = new TFile(Name, "RECREATE");

    
    EDepPerTrack = new TH1F(Titleedep,Titleedep,1000,0,7);
    // Binning chosen to be 250 could try lower still
    if (pT/GeV < 0.57)
    {
        ipPerTrackb_x = new TH1F(Titlebx, Titlebx, 450, -0.35, 0.35);
        ipPerTrackb_y = new TH1F(Titleby, Titleby, 450, -0.35, 0.35);
        ipPerTrackb_z = new TH1F(Titlebz, Titlebz, 400, -0.05, 0.05);
        ipPerTrackb_mag = new TH1F(Titlebmag, Titlebmag, 450, 0.0, 0.6);
    }
    else if(pT/GeV >= 0.57 && pT/GeV < 1.33)
    {
        ipPerTrackb_x = new TH1F(Titlebx, Titlebx, 500, -0.25, 0.25);
        ipPerTrackb_y = new TH1F(Titleby, Titleby, 500, -0.25, 0.25);
        ipPerTrackb_z = new TH1F(Titlebz, Titlebz, 450, -0.04, 0.04);
        ipPerTrackb_mag = new TH1F(Titlebmag, Titlebmag, 450, 0.0, 0.5);
    }
    else if(pT/GeV >= 1.33)
    {
        ipPerTrackb_x = new TH1F(Titlebx, Titlebx, 500, -0.155, 0.155);
        ipPerTrackb_y = new TH1F(Titleby, Titleby, 500, -0.155, 0.155);
        ipPerTrackb_z = new TH1F(Titlebz, Titlebz, 450, -0.04, 0.04);
        ipPerTrackb_mag = new TH1F(Titlebmag, Titlebmag, 450, 0.0, 0.6);
    }
    sigma = new TTree("sigma", "sigma");

    fitresultx = new TF1("fitresultx","gaus");
    fitresulty = new TF1("fitresulty","gaus");
    fitresultz = new TF1("fitresultz","gaus");
#endif
}


void BeamTestRunAction::Write() {
#ifdef G4ANALYSIS_USEROOT
    ipPerTrackb_x->Write();
	ipPerTrackb_y->Write();
	ipPerTrackb_z->Write();
	ipPerTrackb_mag->Write();
    
    
	f->Write();
	f->Close();

    //While TF1s are  owned...
    delete fitresultx;
	delete fitresulty;
	delete fitresultz;
    fitresultx=0;
    fitresulty=0;
    fitresultz=0;
#endif
}

BeamTestRunAction::~BeamTestRunAction() 
{
    delete messager;
}

void BeamTestRunAction::BeginOfRunAction(const G4Run* aRun) 
{
	G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
	//inform the runManager to save random number seed
	//G4RunManager::GetRunManager()->SetRandomNumberStore(false);
}

G4Run*  BeamTestRunAction::GenerateRun()
{
	G4cout<<"Creating user define run class BeamTestRun"<<G4endl;
	return new BeamTestRun("BeamTest/TrackerChamberSD");
}

void BeamTestRunAction::EndOfRunAction(const G4Run* aRun)
{
	// Print interesting data
	G4cout <<"Number of Events Processed:" <<aRun->GetNumberOfEvent() << " events. " <<G4endl;
    G4cout <<"### Run " << aRun->GetRunID() << " end."<<G4endl;
	const BeamTestRun* theRun = dynamic_cast<const BeamTestRun*>(aRun);
	assert (0 != theRun);
    fNameSuffix.str("");    
#ifdef G4ANALYSIS_USEROOT
	// Now we wish to return the standard deviation of the plot so...
	gStyle->SetOptStat(11111111);
	
	if (pT/GeV < 0.57)
	{
		fitresultx->SetRange(-0.35,0.35);
		fitresulty->SetRange(-0.35,0.35);
		fitresultz->SetRange(-0.05,0.05);
	}
	else if(pT/GeV >= 0.57 && pT/GeV < 1.33)
	{
		fitresultx->SetRange(-0.25,0.25);
		fitresulty->SetRange(-0.25,0.25);
		fitresultz->SetRange(-0.04,0.04);
	}
	else if(pT/GeV >= 1.33)
	{
		fitresultx->SetRange(-0.155,0.155);
		fitresulty->SetRange(-0.15,0.155);
		fitresultz->SetRange(-0.04,0.04);
	}


	ipPerTrackb_x->Fit(fitresultx,"QRM0");
	ipPerTrackb_y->Fit(fitresulty,"QRM0");
	ipPerTrackb_z->Fit(fitresultz,"QRM0");
	//h1->DrawCopy();

	msigma_x=fitresultx->GetParameter(2);
	msigma_y=fitresulty->GetParameter(2);
	msigma_z=fitresultz->GetParameter(2);

	G4cout << "The std dev in x is: sigma_x = " << msigma_x << G4endl;
	//Double_t inversePt(0);

    Double_t pTGeV = pT / GeV;
	sigma->Branch("sigma_x", &msigma_x, "msigma_x/D");
	sigma->Branch("sigma_y", &msigma_y, "msigma_y/D");
	sigma->Branch("sigma_z", &msigma_z, "msigma_z/D");
	sigma->Branch("pT", &pTGeV, "pT/D");

	//invertPt(pT,inversePt);
	//sigma->Branch("1/pT", &inversePt, "inversePt/D");

	// Fill Tuple
	sigma->Fill();
#endif
    Write();
	theRun->DumpData();
}

void BeamTestRunAction::FillEvents(G4double bx, G4double by, G4double bz, G4double bmag)
{
	// We want output tofile for bx,by,bz,bmag and initial pT
	/*Double_t bx1(0), by1(0), bz1(0), bmag1(0), pTi(0);
	bx1 = bx;
	by1 = by;
	bz1 = bz;
	bmag1 = bmag;
	pTi = pT;*/
#ifdef G4ANALYSIS_USEROOT
	// Fill histograms
	ipPerTrackb_x->Fill(bx);
	ipPerTrackb_y->Fill(by);
	ipPerTrackb_z->Fill(bz);
	ipPerTrackb_mag->Fill(bmag);
#endif
	/*newTree->Branch("pTi", &pTi, "pTi/D");
	  newTree->Branch("bx1", &bx1, "bx1/D");
	  newTree->Branch("by1", &by1, "by1/D");
	  newTree->Branch("bz1", &bz1, "bz1/D");
	  newTree->Branch("bmag1", &bmag1, "bmag1/D");*/

	//newTree->Fill();


}

void BeamTestRunAction::MergeFiles(G4String pattern)
{
#ifdef G4ANALYSIS_USEROOT
    //Merge all ROOT files according to pattern
    G4cout<<"Merging files matching pattern: "<<pattern<<" to: "<<mergedFileName<<G4endl;
    TChain chain("sigma");
    chain.Add(pattern);
    TObjArray* files = chain.GetListOfFiles();
    for ( Int_t idx = 0 ; idx<files->GetEntries() ; ++idx )
    {
        G4cout<<"File "<<idx<<" "<<files->At(idx)->GetTitle()<<G4endl;
    }
    chain.Merge(mergedFileName);
#endif
}

void BeamTestRunAction::Finalize(G4String macroname)
{
#ifdef G4ANALYSIS_USEROOT
    G4cout<<"Executing analysis macro: "<<macroname<<G4endl;
    if ( macroname.isNull() || macroname.size()==0 ) return;
    TString cmd = ".x "+macroname;
    gROOT->ProcessLine(cmd);
#endif
}

void BeamTestRunAction::SetEnergyDeposit(G4double val ) 
{ 
#ifdef G4ANALYSIS_USEROOT
  EDepPerTrack->Fill(val); 
#endif
}
