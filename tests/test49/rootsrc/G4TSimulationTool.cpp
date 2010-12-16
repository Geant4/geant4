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

#include "G4TSimulationTool.h"

G4TSimulationTool *gSimulationTool = new G4TSimulationTool();


ClassImp(G4TSimulationTool)


using namespace std;
using namespace ROOT;
using namespace TMath;



//______________________________________________________________________________
int G4TSimulationTool::Run(TString const& publicationFile, Int_t crossSection, const TString& modelName,
		Int_t runsNumber, Bool_t useExistingData, TString const& dbPath, TString const& printQue)
{
	G4TSimHelper::LoadLibraries();

	gBenchmark->Start("Overall Benchmark");
	gROOT->Reset();

	// parse the header and load the publication

	cout << "Publication = " << publicationFile << endl;

	if(fPublication != 0) delete fPublication;
	fPublication = gTestDB->LoadData(publicationFile);

	// prepare arguments
	fProjectilePDG 	= fPublication->fHeader.fProjectilePDG;
	fTargetPDG 		= fPublication->fHeader.fTargetPDG;
	fModelName		= modelName;
	fRunsNumber		= runsNumber;
	fEventsNumber	= 2000 * runsNumber;
	fCrossSection	= crossSection;
	fUseExistingData= useExistingData;
	gTestDB->SetDirectory(dbPath);

	// get nn and np
	Int_t nn = gParticlesDAL->GetN(fTargetPDG);
	Int_t np = gParticlesDAL->GetZ(fTargetPDG);

	InternalExecute( np, nn, (Int_t)fPublication->fHeader.fTypeValue, printQue, 2, np );

	gBenchmark->Show("Overall Benchmark");
	return 0;
}

//______________________________________________________________________________
void G4TSimulationTool::InternalExecute(Int_t np, Int_t nn, Int_t e, const TString& pq, Int_t nzone , Int_t nvex, const TString& dir)
{
	// Prepare  the output simulation
	if(fSimulation != 0)  delete fSimulation;
	fSimulation  = new G4TData(fProjectilePDG, fTargetPDG, false, fModelName, E_Kin, e, MeV, 6);

	fSimulation->SetDirectory(gTestDB->GetDirectory());
	fSimulation->SetCrossSection(fCrossSection);
	fSimulation->SetNumberOfEvents(fEventsNumber);

	Initialize();

	if(fPublication == 0 || !fPublication->IsLoaded()){
		cout << "Publication is not defined! Aborting..." << endl;
	}
	else
	{
		Int_t runsNumber = fRunsNumber;

		// Execute Chips
		TString file = gParticlesDAL->GetFileName(fTargetPDG);
		if(fUseExistingData) runsNumber = 0;
		TFile* chipsFile = fSimHelper.ExecuteTest( np, nn, e, runsNumber, fEventsNumber, dir);

		Plot(np, (TTree*)chipsFile->Get("TuplInclTree"), e,  fCrossSection, file);
	}
}


//______________________________________________________________________________
void G4TSimulationTool::Plot(Int_t np, TTree* inclTree, Int_t e, Int_t cs, const TString& file)
{
	if(fPublication == 0 || !fPublication->IsLoaded()){
		cout << "Publication is not defined! Aborting..." << endl;
		return;
	}

	// Get number of events
	Int_t 	nevts = fSimHelper.GetEventsNumber();

	// Prepare the histograms
	cout << "Processing cuts and preparing histograms..." << endl;
	fPublication->PrepareHistograms(fHnbin, fHlxmin, fHlxmax, 0);
	gPlotHelper->PrepareCanvas();

	vector<G4TDataItem*> items = fPublication->GetItems();
	for(UInt_t x = 0; x < items.size(); ++x)
	{
		G4TDataItem* item = items[x];
		Int_t secondaryPDG = item->GetSecondaryParticlePDG();
		Double_t arad = item->GetAngleInRadians();

		// Prepare Cut
		TCut cutAtan = TString::Format("abs(atan2(sqrt(Px*Px+Py*Py),Pz)-%g)<%g", arad, fDanrad).Data();
		TCut cutForSec(gParticlesDAL->GetCut(secondaryPDG));
		TCut cutK = cutForSec && cutAtan;

		TString Expression = TString::Format("(log10(E-m+0.001)) >> %s", item->GetHistogram()->GetName());
		inclTree->Draw(Expression.Data(), cutK, "goff");

		//Save Data
		G4TDataItem* simulationItem = fSimulation->AddItem(item->GetHeader());
		// Copy the tree to the simulation item
		gROOT->cd(); // memory
		TTree* simData = new TTree(item->GetHeader().Data(), item->GetHeader().Data());
		while(simData == 0)
			gSystem->Sleep(10);

		Double_t* ptrV1 = inclTree->GetV1();
		Double_t  valV1 = (*ptrV1);
		simData->Branch("x", &valV1,"x/D");
		simData->Fill();

		for(UInt_t x = 0; x < inclTree->GetSelectedRows(); ++x)
		{
			ptrV1++;
			valV1 = (*ptrV1);
			simData->Fill();
		}
		simulationItem->SetData(simData);


		// scale
		Double_t dOmega = 2 * fPi * (cos(arad - fDanrad) - cos(arad + fDanrad));
		Double_t w		= cs / dOmega / nevts;


		TH1F* hK = item->GetHistogram();
		hK->Divide(fHDT);
		hK->Add(fHZL);
		hK->Scale(w);

	}
	vector<Int_t> fragments = fPublication->GetSecondaryPDGs();

	gPlotHelper->ClearCanvas();
	gPlotHelper->DrawBigTitle(TString::Format(" %s(p,f)X,E%d=MeV, #theta=%s#circ (%s)",
			gParticlesDAL->GetParticleName(fTargetPDG, true).Data(), e, gPlotHelper->PrepareCutEnumeration(fPublication).Data(), fModelName.Data()));
	Int_t padsPerRow = gPlotHelper->DivideForNumber(fragments.size());
	gPlotHelper->SetTitlePosition(33, 0.96, 0.96);


	// LOOP over secondary fragments
	for(UInt_t j = 0; j< fragments.size(); ++j)
	{
		gPlotHelper->PreparePad(j+1);


		// auto size feature
		DataItemLimit_t limits = fPublication->GetLimits(j,padsPerRow);
		if(limits.fMin != 0 && limits.fMax != 0)
		{
			// avoid negative or 0 log scale
			if(limits.fMin <= 0) limits.fMin = 0.00000001;
			if(limits.fMax <= 0) limits.fMax = 0.00000001;

			fHymin = limits.fMin / exp(  0.10 * log(limits.fMax / limits.fMin) );
			fHymax = limits.fMax * exp(  0.10 * log(limits.fMax / limits.fMin) );
		}

		Int_t pdgCode	= fragments[j];
		Double_t m		= gParticlesDAL->GetParticleMass(pdgCode);
		TString Title	= gParticlesDAL->GetParticleName(pdgCode);


		TH1F* frame = gPlotHelper->PrepareFrame(0,fHymin,fHxmax,fHymax,Title);
		if(j == 1){
			// Draw the legend on 2nd pad
			gPlotHelper->DrawAnglesLegend(fPublication);
		}

		gPlotHelper->DrawRightAxis(frame, fHxmax, fHymin, fHymax);

		if(j == 0) 						frame->SetYTitle("d#sigma/pdEd#Omega(mb MeV^{-2} sr^{-1})");
		if(j == fragments.size() - 1)	frame->SetXTitle("T (MeV)");
		gPad->Update();

		// Loop over angles
		vector<G4TDataItem*> sitems = fPublication->GetItemsForSecondary(pdgCode);
		for(UInt_t i = 0; i < sitems.size(); ++i)
		{
			gPad->Update();
			G4TDataItem* item = sitems[i];
			Int_t mk	= gPlotHelper->GetMarker(i);

			TTree* currentTree 	= item->GetData();
			currentTree->Draw(item->GetFormula().Data(),"","goff");

			TGraphErrors *gr = new TGraphErrors(currentTree->GetSelectedRows(), currentTree->GetV2(),
					currentTree->GetV1(), 0, currentTree->GetV3());
			gr->SetMarkerStyle(mk);
			gr->Draw("p");
			RenderHSolid(item->GetHistogram(), fHfbin,fHnbin, m);
		}

	}

	// Save to the DB
	gTestDB->SaveData(fSimulation);
}



