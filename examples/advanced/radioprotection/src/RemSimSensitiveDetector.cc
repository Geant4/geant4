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

#include "RemSimSensitiveDetector.hh"
#include "RemSimHit.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#ifdef G4ANALYSIS_USE
#include "RemSimAnalysisManager.hh"
#endif

RemSimSensitiveDetector::RemSimSensitiveDetector(G4String name)
:G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="trackerCollection");
}

RemSimSensitiveDetector::~RemSimSensitiveDetector(){;}

void RemSimSensitiveDetector::Initialize(G4HCofThisEvent* HCE)
{
  static G4int HCID = -1;
  trackerCollection = new RemSimHitsCollection
                      (SensitiveDetectorName,collectionName[0]); 
  if(HCID<0)
  { HCID = GetCollectionID(0); }
  HCE->AddHitsCollection(HCID,trackerCollection);
}

G4bool RemSimSensitiveDetector::ProcessHits(G4Step* aStep, 
                                            G4TouchableHistory* ROhist)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep==0.) return false;
 
  G4int j = ROhist->GetReplicaNumber();
  G4int k = ROhist->GetReplicaNumber(1);
  G4int i = ROhist->GetReplicaNumber(2);
  
  G4double x = aStep->GetPreStepPoint()->GetPosition().x();
  G4double y = aStep->GetPreStepPoint()->GetPosition().y();
  G4double z = aStep->GetPreStepPoint()->GetPosition().z();

  RemSimHit* newHit = new RemSimHit();
  newHit->SetEdep(edep);
  newHit->SetIndexes(i,j,k);
  newHit->SetPosition( aStep->GetPreStepPoint()->GetPosition());
  trackerCollection->insert( newHit );
  newHit->Draw();
 
  G4int numberOfVoxelX = 100;
  G4int numberOfVoxelY = 100;
  G4int numberOfVoxelZ = 100;
  G4double voxelWidthX = 0.1 *cm;
  G4double voxelWidthY = 0.1 *cm;
  G4double voxelWidthZ = 0.1 *cm;

  G4double xId = (-numberOfVoxelX+1+2*i)*voxelWidthX/2; 
  G4double yId = (- numberOfVoxelY+1+2*j)*voxelWidthY/2;
  G4double zId = (- numberOfVoxelZ+1+2*k)*voxelWidthZ/2;
 
  
  G4cout <<edep/keV <<":in astronaut in position"
    //       <<x/cm<<" "<<y/cm<<" "<<z/cm<<" "
    //   <<"indexes"<<i<<":i "<<j<<":j "<<k<<":k "
         <<"voxels: "<<xId/cm<<" "<<yId/cm<<" "<<zId/cm<<" "
         <<G4endl;
#ifdef G4ANALYSIS_USE
  RemSimAnalysisManager* analysis = RemSimAnalysisManager::getInstance();
  analysis->FillNtupleWithEnergy(xId/cm,yId/cm,zId/cm,edep/keV);
#endif

 return true;
}

void RemSimSensitiveDetector::EndOfEvent(G4HCofThisEvent*)
{
  /*
     G4int NbHits = trackerCollection->entries();
     G4cout << "\n-------->Hits Collection: in this event they are " << NbHits 
            << " hits in the astronaut: " << G4endl;
     for (G4int i=0;i<NbHits;i++) (*trackerCollection)[i]->Print();
  */
}
