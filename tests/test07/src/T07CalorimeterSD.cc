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
//
// $Id: T07CalorimeterSD.cc,v 1.5 2006-06-29 21:37:16 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "T07CalorimeterSD.hh"

#include "T07CalorHit.hh"
#include "T07DetectorConstruction.hh"

#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"

#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

T07CalorimeterSD::T07CalorimeterSD(G4String name,
                                   T07DetectorConstruction* det)
:G4VSensitiveDetector(name),Detector(det)
{
  collectionName.insert("CalCollection");
  HitID = new G4int[500];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

T07CalorimeterSD::~T07CalorimeterSD()
{
  delete [] HitID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void T07CalorimeterSD::Initialize(G4HCofThisEvent*)
{
  CalCollection = new T07CalorHitsCollection
                      (SensitiveDetectorName,collectionName[0]); 
  for (G4int j=0;j<Detector->GetNbOfLayers();j++) {HitID[j] = -1;};
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool T07CalorimeterSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  
  G4double stepl = 0.;
  if (aStep->GetTrack()->GetDefinition()->GetPDGCharge() != 0.)
      stepl = aStep->GetStepLength();
      
  if ((edep==0.)&&(stepl==0.)) return false;      

  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
    
  G4VPhysicalVolume* physVol = theTouchable->GetVolume(); 
  //theTouchable->MoveUpHistory();
  G4int LayerNumber = theTouchable->GetReplicaNumber(1);

  if (HitID[LayerNumber]==-1)
    { 
      T07CalorHit* calHit = new T07CalorHit();
      if (physVol == Detector->GetAbsorber()) calHit->AddAbs(edep,stepl);
      if (physVol == Detector->GetGap     ()) calHit->AddGap(edep,stepl);
      HitID[LayerNumber] = CalCollection->insert(calHit) - 1;
      if (verboseLevel>0)
        G4cout << " New Calorimeter Hit on layer: " << LayerNumber << G4endl;
    }
  else
    { 
      if (physVol == Detector->GetAbsorber())
         (*CalCollection)[HitID[LayerNumber]]->AddAbs(edep,stepl);
      if (physVol == Detector->GetGap())
         (*CalCollection)[HitID[LayerNumber]]->AddGap(edep,stepl);
      if (verboseLevel>0)
        G4cout << " Energy added to Layer: " << LayerNumber << G4endl; 
    }
    
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void T07CalorimeterSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  static G4int HCID = -1;
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection(HCID,CalCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void T07CalorimeterSD::clear()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void T07CalorimeterSD::DrawAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void T07CalorimeterSD::PrintAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

