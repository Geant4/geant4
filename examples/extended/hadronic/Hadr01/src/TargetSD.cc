//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// -------------------------------------------------------------
//
//
//      ---------- TargetSD -------------
//              
//  Modified:
//
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "TargetSD.hh"

#include "G4RunManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4Positron.hh"
#include "globals.hh"
#include "HistoManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Step.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TargetSD::TargetSD(const G4String& name)
 :G4VSensitiveDetector(name),
  theHisto(HistoManager::GetPointer()),
  evno(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TargetSD::~TargetSD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TargetSD::Initialize(G4HCofThisEvent*)
{
  evno++;
  if(0 < theHisto->GetVerbose())
    G4cout << "TargetSD: Begin Of Event # " << evno << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool TargetSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  theHisto->AddStep(aStep->GetTrack()->GetDefinition(), 
		    aStep->GetPreStepPoint()->GetKineticEnergy() );
  G4double edep = aStep->GetTotalEnergyDeposit();
  if(0.0 == edep) return true;

  G4double length = aStep->GetStepLength();

  if(aStep->GetTrack()->GetTrackID() == 1) theHisto->AddTrackLength(length);

  G4ThreeVector p1 = aStep->GetPreStepPoint()->GetPosition();
  G4ThreeVector p2 = aStep->GetPostStepPoint()->GetPosition();
  p1 += p2;
  p1 *= 0.5;

  theHisto->AddEnergy(edep,length,p1);

  if(1 < theHisto->GetVerbose()) {
      G4cout << "TargetSD: energy = " << edep/MeV
             << " MeV is deposited at the step from " << p1
             << " to " << p2 << G4endl;
  }

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TargetSD::EndOfEvent(G4HCofThisEvent*)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TargetSD::clear()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void TargetSD::PrintAll()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

