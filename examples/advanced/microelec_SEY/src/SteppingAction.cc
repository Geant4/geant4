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
// -------------------------------------------------------------------
// -------------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "SteppingAction.hh"
#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4SteppingManager.hh"
#include "G4VTouchable.hh"
#include "G4VPhysicalVolume.hh"

#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SteppingAction::UserSteppingAction(const G4Step* step)
{
	G4String ParticleName, PostProcessName, PreProcessName, VolumeNamePre, MaterialNamePre, VolumeNamePost, MaterialNamePost;
	PreProcessName = "oooo";

	if (step->GetPostStepPoint() && step->GetPreStepPoint() && step->GetPreStepPoint()->GetPhysicalVolume() &&
		step->GetPostStepPoint()->GetPhysicalVolume() && step->GetPostStepPoint()->GetProcessDefinedStep() &&
		step->GetTrack()) {

		PostProcessName = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
		if (step->GetPreStepPoint()->GetProcessDefinedStep()) PreProcessName = step->GetPreStepPoint()->GetProcessDefinedStep()->GetProcessName();
		MaterialNamePre = step->GetPreStepPoint()->GetMaterial()->GetName();
		VolumeNamePre = step->GetPreStepPoint()->GetPhysicalVolume()->GetName();
		MaterialNamePost = step->GetPostStepPoint()->GetMaterial()->GetName();
		VolumeNamePost = step->GetPostStepPoint()->GetPhysicalVolume()->GetName();

		ParticleName = step->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName();
		 
		
	}

 }
