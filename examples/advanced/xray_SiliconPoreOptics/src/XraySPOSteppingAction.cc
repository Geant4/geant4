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
/// \file XraySPOSteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
// Authors: P.Dondero (paolo.dondero@cern.ch), R.Stanzani (ronny.stanzani@cern.ch)
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#include "XraySPOSteppingAction.hh"
#include "XraySPODetectorConstruction.hh"
#include "XraySPOHistoManager.hh"
#include "G4SteppingManager.hh"
#include "G4RunManager.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void XraySPOSteppingAction::UserSteppingAction(const G4Step* step)
{
  G4int parentID = step->GetTrack()->GetParentID();

  // Only primaries
  if (parentID == 0)
  {
    G4RunManager *rm = G4RunManager::GetRunManager();
    G4int eventID = rm->GetCurrentEvent()->GetEventID();

    G4double len_fraction = 0.1; // approx 1/20 of the lenght of the pore. Particles scatters an avg of 5 times, so by selecting all events with a distance >10mm we have all the scatters inside.
    G4bool proceed = false;

    if (fPrevEventID != eventID)
    {
      fNumReflections = -2;
    }

    G4ThreeVector direction = step->GetTrack()->GetMomentumDirection();
    auto theta = (G4double)direction.theta();
    auto phi = (G4double)direction.phi();

    G4double x = step->GetPreStepPoint()->GetPosition().x();
    G4double y = step->GetPreStepPoint()->GetPosition().y();
    G4double z = step->GetPreStepPoint()->GetPosition().z();

    const G4String& vol_name = step->GetTrack()->GetVolume()->GetName();
    G4int trackid = step->GetTrack()->GetTrackID();

    G4StepPoint * postStep = step->GetPostStepPoint();
    G4StepPoint* post_point = step->GetPostStepPoint();

    const G4String& proc_name = postStep->GetProcessDefinedStep()->GetProcessName();

    // Count the number of reflections
    G4String s = "";
    G4String& post_phys_name = s;
    if (fPrevEventID == eventID && step->GetStepLength() >= len_fraction)
    {
      if (post_point != nullptr)
      {
        G4VPhysicalVolume* post_phys = post_point->GetPhysicalVolume();
        if (post_phys != nullptr)
        {
          post_phys_name = post_phys->GetName();
        }

        if (post_phys_name == "cone_1" || post_phys_name == "cone_2" || post_phys_name == "cone_3" || post_phys_name == "cone_4" || post_phys_name == "cone_5" || post_phys_name == "cone_6" || post_phys_name == "cone_7" || post_phys_name == "cone_8")
        {
          if (fNumReflections < 0)
          {
            fNumReflections = 1;
          }
          else
          {
            fNumReflections += 1;
          }
          proceed = true;
        }
      }
    }

    if (vol_name == "pDummyEntrance" || vol_name == "pDummyExit" || vol_name == "pDummySphere")
    {
      proceed = true;
    }

    if (proceed)
    {
      // G4cout << "Save with " << fNumReflections << " reflections" << G4endl;
      G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
      analysisManager->FillNtupleIColumn(0, eventID);
      analysisManager->FillNtupleSColumn(1, vol_name);
      analysisManager->FillNtupleIColumn(2, trackid);
      analysisManager->FillNtupleDColumn(3, x);
      analysisManager->FillNtupleDColumn(4, y);
      analysisManager->FillNtupleDColumn(5, z);
      analysisManager->FillNtupleDColumn(6, theta);
      analysisManager->FillNtupleDColumn(7, phi);
      analysisManager->FillNtupleSColumn(8, proc_name);
      analysisManager->FillNtupleIColumn(9, parentID);
      analysisManager->FillNtupleIColumn(10, fNumReflections);
      analysisManager->AddNtupleRow();
    }
    fPrevx = x;
    fPrevy = y;
    fPrevz = z;
    fPrevEventID = eventID;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
