// -------------------------------------------------------------
//
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------
//      GEANT4 hTest
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- hTestCalorimeterSD -------------
//              
//  Modified: 05.04.01 Vladimir Ivanchenko new design of hTest 
// 
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "hTestCalorimeterSD.hh"

#include "G4RunManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestCalorimeterSD::hTestCalorimeterSD(G4String name, 
                                       hTestDetectorConstruction* det)
 :G4VSensitiveDetector(name),
  theDet(det),
  theEvent(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestCalorimeterSD::~hTestCalorimeterSD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestCalorimeterSD::Initialize(G4HCofThisEvent*)
{
  if(!theEvent) theEvent = theDet->GetEventAction();
  verbose = theEvent->GetVerbose();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool hTestCalorimeterSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();

  if(0.0 < edep) {
    G4int j = aStep->GetTrack()->GetVolume()->GetCopyNo();
    theEvent->AddEnergy(edep, j);
    if(1 < verbose) {
      G4cout << "hTestCalorimeterSD: energy = " << edep/MeV
             << " MeV is deposited at " << j
             << "-th absorber slice " << G4endl;
    }
  }

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestCalorimeterSD::EndOfEvent(G4HCofThisEvent* HCE)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestCalorimeterSD::clear()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void hTestCalorimeterSD::PrintAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....







