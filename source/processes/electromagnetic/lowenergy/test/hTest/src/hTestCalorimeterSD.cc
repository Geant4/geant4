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

hTestCalorimeterSD::hTestCalorimeterSD(G4String name)
 :G4VSensitiveDetector(name)
{
  theRun = (G4RunManager::G4RunManager())->GetUserRunAction();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestCalorimeterSD::~hTestCalorimeterSD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestCalorimeterSD::Initialize(G4HCofThisEvent*HCE)
{
  verbose = run->GetVerbose();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool hTestCalorimeterSD::ProcessHits(G4Step* aStep,G4TouchableHistory* h)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  G4double z = 0.0;

  if(0.0 < edep) {
    G4double z1 = (aStep->GetPreStepPoint()->GetPosition()).z();
    G4double z2 = (aStep->GetPostStepPoint()->GetPosition()).z();
    z  = (z1 + z2)*0.5;
    run->AddEnergy(edep, z);
  }

  if(1 < verbose) {
    G4cout << "hTestCalorimeterSD: energy = " << edep/MeV
           << " MeV is deposited at Z = " << z/mm
           << " mm " << G4endl;

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







