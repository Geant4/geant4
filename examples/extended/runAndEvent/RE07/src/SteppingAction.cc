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
/// \file src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"

#include "G4PhysicalConstants.hh"
#include "G4Positron.hh"
#include "G4RunManager.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "Run.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det, EventAction* evt)
  : G4UserSteppingAction(), fDetector(det), fEventAct(evt)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  // track informations
  const G4StepPoint* prePoint = aStep->GetPreStepPoint();

  // if World, return
  //
  G4VPhysicalVolume* volume = prePoint->GetTouchableHandle()->GetVolume();
  // if sum of absorbers do not fill exactly a layer: check material, not
  // volume.
  const G4Material* mat = volume->GetLogicalVolume()->GetMaterial();
  if (mat == fDetector->GetWorldMaterial()) return;

  const G4ParticleDefinition* particle = aStep->GetTrack()->GetDefinition();

  // here we are in an absorber. Locate it
  //
  G4int absorNum = prePoint->GetTouchableHandle()->GetCopyNumber(0);
  // G4int layerNum  = prePoint->GetTouchableHandle()->GetCopyNumber(1);

  // get Run
  Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  // collect energy deposit taking into account track weight
  G4double edep = aStep->GetTotalEnergyDeposit() * aStep->GetTrack()->GetWeight();

  // collect step length of charged particles
  G4double stepl = 0.;
  if (particle->GetPDGCharge() != 0.) {
    stepl = aStep->GetStepLength();
    run->AddChargedStep();
  }
  else {
    run->AddNeutralStep();
  }

  // sum up per event
  fEventAct->SumEnergy(absorNum, edep, stepl);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
