// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Test17CalorimeterSD.hh"

#include "Test17DetectorConstruction.hh"

#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
  
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test17CalorimeterSD::Test17CalorimeterSD(G4String name,
                     Test17DetectorConstruction* det): 
  G4VSensitiveDetector(name),
  detector(det)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test17CalorimeterSD::~Test17CalorimeterSD()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17CalorimeterSD::Initialize(G4HCofThisEvent*HCE)
{
  if(!HCE)G4cout << "Sensitive detector is inicialized" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool Test17CalorimeterSD::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  
  G4double stepl = 0.;
  if (aStep->GetTrack()->GetDefinition()->GetPDGCharge() != 0.)
      stepl = aStep->GetStepLength();
      
  if ((edep==0.)&&(stepl==0.)) return false;      

  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
    
  G4VPhysicalVolume* physVol = theTouchable->GetVolume(); 
  G4LogicalVolume* logVol = physVol->GetLogicalVolume();
  //theTouchable->MoveUpHistory();
    
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17CalorimeterSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  if(!HCE)G4cout << "End of event" << G4endl;
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17CalorimeterSD::clear()
{;} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void Test17CalorimeterSD::PrintAll()
{;} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....









