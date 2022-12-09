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
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayTelCalorimeterSD  ------
//           by  R.Giannitrapani, F.Longo & G.Santin (13 nov 2000)
//
// ************************************************************

#include "G4RunManager.hh"
#include "GammaRayTelCalorimeterSD.hh"
#include "GammaRayTelCalorimeterHit.hh"
#include "GammaRayTelDetectorConstruction.hh"

#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelCalorimeterSD::GammaRayTelCalorimeterSD(G4String name) :	G4VSensitiveDetector(name) {
	auto *runManager = G4RunManager::GetRunManager();
	detector = (GammaRayTelDetectorConstruction*) (runManager->GetUserDetectorConstruction());

	numberOfCALBars = detector->GetNbOfCALBars();
	numberOfCALLayers = detector->GetNbOfCALLayers();

	// G4cout << NbOfCALBars << " bars " << G4endl;
	// G4cout << NbOfCALLayers << " layers " << G4endl;

	numberOfCALChannels = numberOfCALBars * numberOfCALLayers;

	calHitXID = new G4int[numberOfCALChannels];
	calHitYID = new G4int[numberOfCALChannels];
	collectionName.insert("CalorimeterCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelCalorimeterSD::~GammaRayTelCalorimeterSD() {
	delete[] calHitXID;
	delete[] calHitYID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelCalorimeterSD::Initialize(G4HCofThisEvent*) {
    calorimeterCollection = new GammaRayTelCalorimeterHitsCollection(SensitiveDetectorName, collectionName[0]);
	for (auto i = 0; i < numberOfCALChannels; i++) {
		calHitXID[i] = -1;
		calHitYID[i] = -1;
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

auto GammaRayTelCalorimeterSD::ProcessHits(G4Step *step, G4TouchableHistory*) -> G4bool {
	G4double depositedEnergy = 0.;
	depositedEnergy = step->GetTotalEnergyDeposit();
	if (depositedEnergy == 0.) {
		return false;
	}

	// This TouchableHistory is used to obtain the physical volume of the hit
	auto *theTouchable = (G4TouchableHistory*) (step->GetPreStepPoint()->GetTouchable());
	auto *calorimeterBar = theTouchable->GetVolume();
	auto *calorimeterPlane = theTouchable->GetVolume(1);
	auto calorimeterBarNumber = calorimeterBar->GetCopyNo();
	auto calorimeterBarName = calorimeterBar->GetName();

	G4int planeNumber{0};
	planeNumber = calorimeterPlane->GetCopyNo();

	auto planeName = calorimeterPlane->GetName();

	G4int NChannel = 0;
	NChannel = planeNumber * numberOfCALBars + calorimeterBarNumber;

	if (planeName == "CALLayerX") { // The hit is on a X CsI (cesium iodide) plane
		if (calHitXID[NChannel] == -1) { // This is a new hit
			auto *calorimeterHit = new GammaRayTelCalorimeterHit;
			calorimeterHit->SetCALType(1);
			calorimeterHit->AddEnergy(depositedEnergy);
			calorimeterHit->SetPosition(step->GetPreStepPoint()->GetPosition());
			calorimeterHit->SetCALPlaneNumber(planeNumber);
			calorimeterHit->SetCALBarNumber(calorimeterBarNumber);
			calHitXID[NChannel] = calorimeterCollection->insert(calorimeterHit) - 1;
		} else { // This is not new
			(*calorimeterCollection)[calHitXID[NChannel]]->AddEnergy(depositedEnergy);
		}
	}

	if (planeName == "CALLayerY") { // The hit is on an Y CsI (cesium iodide) plane
		if (calHitYID[NChannel] == -1) { // This is a new hit
			auto *calorimeterHit = new GammaRayTelCalorimeterHit;
			calorimeterHit->SetCALType(0);
			calorimeterHit->AddEnergy(depositedEnergy);
			calorimeterHit->SetPosition(step->GetPreStepPoint()->GetPosition());
			calorimeterHit->SetCALPlaneNumber(planeNumber);
			calorimeterHit->SetCALBarNumber(calorimeterBarNumber);
			calHitYID[NChannel] = calorimeterCollection->insert(calorimeterHit) - 1;
		} else { // This is not new
			(*calorimeterCollection)[calHitYID[NChannel]]->AddEnergy(depositedEnergy);
		}
	}

	return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelCalorimeterSD::EndOfEvent(G4HCofThisEvent *HCE) {
	static G4int collectionIdentifier = -1;
	if (collectionIdentifier < 0) {
		collectionIdentifier = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
	}
	HCE->AddHitsCollection(collectionIdentifier, calorimeterCollection);

	for (auto i = 0; i < numberOfCALChannels; i++) {
		calHitXID[i] = -1;
		calHitYID[i] = -1;
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelCalorimeterSD::clear() {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelCalorimeterSD::DrawAll() {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelCalorimeterSD::PrintAll() {
}
