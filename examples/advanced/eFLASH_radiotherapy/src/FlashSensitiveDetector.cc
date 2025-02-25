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
#include "FlashSensitiveDetector.hh"

#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"

FlashSensitiveDetector::FlashSensitiveDetector(const G4String& name) : G4VSensitiveDetector(name)
{
    
}

FlashSensitiveDetector::~FlashSensitiveDetector()
{

}


G4bool FlashSensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*)

{
    G4double edep = step -> GetTotalEnergyDeposit();

    if (edep == 0) return false;

    const G4VTouchable *touchable = step -> GetPreStepPoint() -> GetTouchable();
    G4VPhysicalVolume *physVol = touchable -> GetVolume();
    G4ThreeVector posDetector = physVol -> GetTranslation();  // get the position

    G4Track *track = step -> GetTrack();

    auto material = track -> GetMaterial();
    auto density = material -> GetDensity();
    G4double cavityVolume = (track->GetVolume()->GetLogicalVolume()->GetSolid()->GetCubicVolume());
    G4double massOfCavity = cavityVolume * density;
    G4double dose = edep / massOfCavity / gray;
    edep = edep / MeV;



    G4AnalysisManager *fAnalysisManager = G4AnalysisManager::Instance();
    fAnalysisManager -> FillNtupleDColumn(0, 0, posDetector[0]);
    fAnalysisManager -> FillNtupleDColumn(0, 1, posDetector[1]);
    fAnalysisManager -> FillNtupleDColumn(0, 2, posDetector[2]);
    fAnalysisManager -> FillNtupleDColumn(0, 3, dose);
    fAnalysisManager -> FillNtupleDColumn(0, 4, edep);
    fAnalysisManager -> FillNtupleIColumn(0, 5, G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID());
    fAnalysisManager -> FillNtupleIColumn(0, 6, track->GetParentID());
    fAnalysisManager -> FillNtupleSColumn(0, 7, track->GetDynamicParticle()->GetDefinition()->GetParticleName());
    fAnalysisManager -> AddNtupleRow(0);

    return true;
}
