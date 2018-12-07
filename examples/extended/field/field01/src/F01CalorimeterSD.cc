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
/// \file field/field01/src/F01CalorimeterSD.cc
/// \brief Implementation of the F01CalorimeterSD class
//
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "F01CalorimeterSD.hh"

#include "F01CalorHit.hh"
#include "F01DetectorConstruction.hh"

#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F01CalorimeterSD::F01CalorimeterSD(G4String name,
                                   F01DetectorConstruction* det)
: G4VSensitiveDetector(name),
  fCalCollection(0),
  fDetector(det),
  fHitID(new G4int[500])
{
  collectionName.insert("CalCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F01CalorimeterSD::~F01CalorimeterSD()
{
  delete [] fHitID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F01CalorimeterSD::Initialize(G4HCofThisEvent*)
{
  fCalCollection = new F01CalorHitsCollection
                       (SensitiveDetectorName,collectionName[0]);
  for (G4int j=0;j<1; j++) {fHitID[j] = -1;};
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool F01CalorimeterSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  G4double edep = step->GetTotalEnergyDeposit();
  G4double stepl = 0.;

  stepl = step->GetStepLength();

  if ((edep == 0.) && (stepl == 0.) ) return false;

  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(step->GetPreStepPoint()->GetTouchable());

  G4VPhysicalVolume* physVol = theTouchable->GetVolume();

  G4int number = 0;
  if (fHitID[number]==-1)
    {
      F01CalorHit* calHit = new F01CalorHit();
      if (physVol == fDetector->GetAbsorber()) calHit->AddAbs(edep,stepl);
      fHitID[number] = fCalCollection->insert(calHit) - 1;
      if (verboseLevel>0)
        G4cout << " New Calorimeter Hit on F01: " << number << G4endl;
    }
  else
    {
      if (physVol == fDetector->GetAbsorber())
         (*fCalCollection)[fHitID[number]]->AddAbs(edep,stepl);
      if (verboseLevel>0)
        G4cout << " Energy added to F01: " << number << G4endl; 
    }
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F01CalorimeterSD::EndOfEvent(G4HCofThisEvent* hce)
{
  static G4int hcID = -1;
  if (hcID<0)
  { hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  hce->AddHitsCollection(hcID,fCalCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
