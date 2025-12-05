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
//---------------------------------------------------------------------------
//
// ClassName:   MicroElecSdSey
//
// Description: The process to kill e- to save CPU
//
// Author:      C. Inguimbert 16/02/2022
//				ONERA
//----------------------------------------------------------------------------
//
// Class description:
//
// SEY : Secondary Electron Emission Yield
// detecteor to be used to count the number of secondary electrons
// emitted by a surface of irradiated material
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#include	<iostream>
#include	<fstream>
#include	<sstream>


#include "G4RunManager.hh"

//#include "DepthSampling.hh"

#include	"SteppingAction.hh"
#include	"MicroElecSdSey.hh"
#include	"MicroElecRun.hh"
#include	"G4StepPoint.hh"
#include	"G4HCofThisEvent.hh"
#include	"G4Step.hh"
#include	"G4ThreeVector.hh"
#include	"G4SDManager.hh"
#include	"G4ios.hh"
#include	"G4VProcess.hh"
#include	"G4EmCalculator.hh"
#include	"G4Material.hh"
#include	"G4UnitsTable.hh"
#include	"G4EnergyLossTables.hh"  
#include	"G4VPhysicalVolume.hh"
#include	"PrimaryGeneratorAction.hh"
#include	"G4RunManager.hh"
#include	"DetectorConstruction.hh"
#include    "G4ParticleDefinition.hh"
#include	"G4SystemOfUnits.hh"
#include	"G4VHit.hh"
#include	"G4AttDef.hh"
//#include "Analysis.hh"



using namespace std;

#include	<vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
MicroElecSdSey::MicroElecSdSey(const G4String& name,
	const G4String& hitsCollectionName)
	: G4VSensitiveDetector(name), fHitsCollection(NULL)
{
	compteurPrimaire = 0.0, compteurSec= 0.0, compteurTot= 0.0, compteur50= 0.0, nbPrim= 0.0, nbSec= 0.0 , nbSup50= 0.0;
	ResetCounters();

	//G4String HCname;
	//collectionName.insert(HCname = "DetecteurMicroElecSdSey");
	G4String HCname = hitsCollectionName;
	collectionName.insert(HCname);
}

MicroElecSdSey::~MicroElecSdSey() { }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//Called at beginning of event
void MicroElecSdSey::Initialize(G4HCofThisEvent* HCE)
{

	static G4int HCID = -1;

	if (HCID < 0)
	{
		HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
	}

	fHitsCollection = new MicroElecHitSeyCollection();

	HCE->AddHitsCollection(HCID, fHitsCollection);
	

	nbPrim = 0; nbSec = 0; nbSup50 = 0;

	}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//called for each step in sensitive volume
G4bool MicroElecSdSey::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
	G4Track* aTrack = aStep->GetTrack();
	G4int trackID = aStep->GetTrack()->GetTrackID();
	G4int ParentID = aStep->GetTrack()->GetParentID();
	
	
	G4String processName;
	const G4VProcess* creatorProcess = aStep->GetTrack()->GetCreatorProcess();
	if (creatorProcess) { processName = creatorProcess->GetProcessName(); }
	if (aStep->GetPostStepPoint()) { G4String processNamePSP = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName(); }
	
	
	const G4DynamicParticle* particleDyn = aStep->GetTrack()->GetDynamicParticle();
	G4ParticleDefinition* particleDef = particleDyn->GetDefinition();
	G4int PDGEncoding = particleDef->GetPDGEncoding();
	G4String ParticleTypeName = particleDef->GetParticleType();
	G4String ParticleName = particleDef->GetParticleName();

	G4int Z = particleDef->GetAtomicNumber();
	G4int A = particleDef->GetAtomicMass();
	G4double Vertex_kinetic_energy = aStep->GetTrack()->GetVertexKineticEnergy();
	G4ThreeVector Vertex_Position = aStep->GetTrack()->GetVertexPosition();
	G4ThreeVector Vertex_Moment = aStep->GetTrack()->GetVertexMomentumDirection();

	G4double	Particle_kinetic_energy=999.;
	G4double	Post_kinetic_energy=999.;
	G4ThreeVector PreImpulsion;
	G4ThreeVector PostImpulsion;
	G4ThreeVector PrePosition;
	G4ThreeVector PostPosition;

	if (aStep->GetPreStepPoint()) { 
		Particle_kinetic_energy = aStep->GetPreStepPoint()->GetKineticEnergy(); 
		PreImpulsion = aStep->GetPreStepPoint()->GetMomentumDirection();
		PrePosition = aStep->GetPreStepPoint()->GetPosition();
	}
	if (aStep->GetPostStepPoint()) { 
		Post_kinetic_energy = aStep->GetPostStepPoint()->GetKineticEnergy(); 
		PostImpulsion = aStep->GetPostStepPoint()->GetMomentumDirection();
		PostPosition = aStep->GetPostStepPoint()->GetPosition();
	}


	G4double edep = aStep->GetTotalEnergyDeposit();
	G4double NIedep = aStep->GetNonIonizingEnergyDeposit();
	G4double stepLength = aStep->GetStepLength();

	G4String			originalVolume_;
	G4String			originalVolume_materialName;
	originalVolume_ = aStep->GetTrack()->GetTouchableHandle()->GetVolume()->GetName();
	originalVolume_materialName = aStep->GetTrack()->GetMaterial()->GetName();

	
	if (ParentID == 0)	{nbPrim++;}
	if (ParentID > 0)	{nbSec++;}
	if (Particle_kinetic_energy/eV >= 50.0) { nbSup50++; }



	MicroElecHitSey* newHit = new MicroElecHitSey();
	
	newHit->SetNbPrim(nbPrim);
	newHit->SetNbSec(nbSec);
	newHit->SetNbSup50(nbSup50); 
	newHit->SetPDGEncoding(PDGEncoding);
	newHit->SetParentID(ParentID);
	newHit->SetTrackID(trackID);
	newHit->SetParticleType(ParticleTypeName);
	newHit->SetParticleName(ParticleName);
	newHit->SetVolumeName(originalVolume_);
	newHit->SetZ(Z);
	newHit->SetA(A);
	newHit->SetVertexKineticEnergy(Vertex_kinetic_energy);
	newHit->SetVertexPos(Vertex_Position);
	newHit->SetVertexMomentum(Vertex_Moment);
	newHit->SetPreStepKineticEnergy(Particle_kinetic_energy);
	newHit->SetPostStepKineticEnergy(Post_kinetic_energy);
	newHit->SetEdep(edep);
	newHit->SetNi_Edep(NIedep);
	newHit->SetStepLength(stepLength);
	newHit->SetPreStepMomentum(PreImpulsion);
	newHit->SetPostStepMomentum(PostImpulsion);
	newHit->SetPrePos(PrePosition);
	newHit->SetPostPos(PostPosition);
	
	fHitsCollection->insert(newHit);
	aTrack->SetTrackStatus(fStopAndKill);
		
		return true;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//called at end of event
void MicroElecSdSey::EndOfEvent(G4HCofThisEvent*)
{

	MicroElecRun* run = static_cast<MicroElecRun*>(
		G4RunManager::GetRunManager()->GetNonConstCurrentRun());

	std::stringstream ss;
	ss << 1 + run->GetRunID() << G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

	compteurPrimaire = compteurPrimaire + nbPrim;
	compteurSec = compteurSec + nbSec;
	compteur50 = compteur50 + nbSup50;
	compteurTot = compteurTot+ nbPrim+ nbSec;

	run->AddElecPrimScorer(nbPrim);
	run->AddElecSecoScorer(nbSec);
	run->AddElecSup50Scorer(nbSup50);
	run->AddElecTotaScorer(nbPrim+ nbSec);
	

	nbPrim = 0.0;
	nbSec = 0.0;
	nbSup50 = 0.0;


}




