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
/// \file ChemistryTrackingManager.cc
/// \brief Implementation of the ChemistryTrackingManager class

// The `molcounters` example(s) are provided as part of Geant4-DNA
// and any report or published result obtained using it shall cite
// the respective Geant4-DNA collaboration publications.
//
// Reports or results obtained using the spatially-aware `MoleculeCounter`
// provided in this example, shall further cite:
//
// Velten & Tomé, Radiation Physics and Chemistry, 2023 (10.1016/j.radphyschem.2023.111194)
//
//
// Author: Christian Velten (2025)
//

#include "ChemistryTrackingManager.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4MoleculeCounterManager.hh"
#include "G4UserSteppingAction.hh"
#include "G4VSensitiveDetector.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ChemistryTrackingManager::~ChemistryTrackingManager()
{
  // Ensure this manager's stepping action is not handled by the event manager
  auto eventManager = G4EventManager::GetEventManager();
  if (!(!eventManager || fUserSteppingAction == eventManager->GetUserSteppingAction()))
    delete fUserSteppingAction;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ChemistryTrackingManager::AppendStep(G4Track* /*track*/, G4Step* step)
{
  if (step->GetPreStepPoint()->GetPhysicalVolume() != nullptr
      && step->GetControlFlag() != AvoidHitInvocation)
  {
    auto sensitiveDetector = step->GetPreStepPoint()->GetSensitiveDetector();
    if (sensitiveDetector != nullptr) {
      sensitiveDetector->Hit(step);
    }
  }

  if (fUserSteppingAction) fUserSteppingAction->UserSteppingAction(step);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ChemistryTrackingManager::Finalize()
{
  if (G4MoleculeCounterManager::Instance()->GetIsActive())
    G4MoleculeCounterManager::Instance()->NotifyOfFinalize();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......