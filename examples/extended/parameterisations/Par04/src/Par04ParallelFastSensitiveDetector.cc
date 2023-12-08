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
#include "Par04ParallelFastSensitiveDetector.hh"
#include <CLHEP/Vector/Rotation.h>     // for HepRotation
#include <CLHEP/Vector/ThreeVector.h>  // for Hep3Vector
#include <cmath>                       // for floor
#include <G4CollectionNameVector.hh>   // for G4CollectionNameVector
#include <G4FastHit.hh>                // for G4FastHit
#include <G4FastTrack.hh>              // for G4FastTrack
#include <G4RotationMatrix.hh>         // for G4RotationMatrix
#include <G4StepPoint.hh>              // for G4StepPoint
#include <G4THitsCollection.hh>        // for G4THitsCollection
#include <G4ThreeVector.hh>            // for G4ThreeVector
#include <G4VSensitiveDetector.hh>     // for G4VSensitiveDetector
#include <G4VUserEventInformation.hh>  // for G4VUserEventInformation
#include <cstddef>                     // for size_t
#include <vector>                      // for vector
#include "G4Event.hh"                  // for G4Event
#include "G4EventManager.hh"           // for G4EventManager
#include "G4HCofThisEvent.hh"          // for G4HCofThisEvent
#include "G4SDManager.hh"              // for G4SDManager
#include "G4Step.hh"                   // for G4Step
#include "G4Track.hh"                  // for G4Track
#include "Par04EventInformation.hh"    // for Par04EventInformation
#include "Par04Hit.hh"                 // for Par04Hit, Par04HitsCollection

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04ParallelFastSensitiveDetector::Par04ParallelFastSensitiveDetector(G4String aName)
 : G4VSensitiveDetector(aName)
{
  collectionName.insert("physicalCellsFastSim");
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04ParallelFastSensitiveDetector::Par04ParallelFastSensitiveDetector(G4String aName,
                                                                       G4int aNbOfLayers,
                                                                       G4int aNbOfSlices)
: G4VSensitiveDetector(aName)
  , fNbOfLayers(aNbOfLayers)
  , fNbOfSlices(aNbOfSlices)
{
  collectionName.insert("physicalCellsFastSim");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04ParallelFastSensitiveDetector::~Par04ParallelFastSensitiveDetector() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04ParallelFastSensitiveDetector::Initialize(G4HCofThisEvent* aHCE)
{
  fHitsCollection = new Par04HitsCollection(SensitiveDetectorName, collectionName[0]);
  if(fHitCollectionID < 0)
  {
    fHitCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
  }
  aHCE->AddHitsCollection(fHitCollectionID, fHitsCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool Par04ParallelFastSensitiveDetector::ProcessHits(G4Step*, G4TouchableHistory* )
{
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool Par04ParallelFastSensitiveDetector::ProcessHits(const G4FastHit* aHit,
                                                       const G4FastTrack*,
                                                       G4TouchableHistory* aTouchable)
{
  G4double edep = aHit->GetEnergy();
  if(edep == 0.)
    return true;

  G4int layerNo = aTouchable->GetCopyNumber(0);
  G4int sliceNo = aTouchable->GetCopyNumber(1);
  G4int rowNo = aTouchable->GetCopyNumber(2);
  
  G4int hitID = fNbOfLayers*fNbOfSlices*rowNo+fNbOfLayers*sliceNo+layerNo;
  auto hit = fHitsMap[hitID].get();
  if (hit==nullptr) {
    fHitsMap[hitID] = std::unique_ptr<Par04Hit>(new Par04Hit());
    hit = fHitsMap[hitID].get();
    hit->SetPhiId(sliceNo);
    hit->SetRhoId(layerNo);
    hit->SetZid(rowNo);
  }
  if (layerNo != hit->GetRhoId()) {
    G4cout << "ERROR, problem with  Layer IDs: " << layerNo << " != " << hit->GetRhoId() << G4endl;
    return false;
  }
  if (sliceNo != hit->GetPhiId()) {
    G4cout << "ERROR, problem with Slice IDs: " << sliceNo << " != " << hit->GetPhiId() << G4endl;
    return false;
  }
  if (rowNo != hit->GetZid()) {
    G4cout << "ERROR, problem with Row  IDs: " << rowNo << " != " << hit->GetZid() << G4endl;
    return false;
  }

  // Add energy deposit from G4Step
  hit->AddEdep(edep);
  // Increment the counter
  hit->AddNdep(1);
  // Set hit type to parallel world fast simulation (even if hit was already marked as full
  // sim, overwrite it)
  hit->SetType(3);

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Par04ParallelFastSensitiveDetector::EndOfEvent(G4HCofThisEvent*)
{
  for(const auto& hits: fHitsMap){
    fHitsCollection->insert(new Par04Hit(*hits.second.get()));
  }
  fHitsMap.clear();
}
