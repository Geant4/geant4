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
//
//    **************************************
//    *                                    *
//    *           CellTrackerSD.cc         *
//    *                                    *
//    **************************************
//
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//	   Barbara Mascialino (Barbara.Mascialino@ge.infn.it)
//
// History:
// -----------
// 20 September 2006 S. Guatelli, B. Mascialino      first implementation
// -------------------------------------------------------------------

#include "CellTrackerSD.hh"
#include "CellRunAction.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "CellDetectorConstruction.hh"
#include "CellAnalysisManager.hh"
#include "CellTrackerHit.hh"
#include "G4RunManager.hh"

CellTrackerSD::CellTrackerSD(G4String name,
			       G4double XSize,
			       G4double YSize,
			       G4double ZSize,
			       G4int voxelX,
                               G4int voxelY,
			       G4int voxelZ,
			       CellDetectorConstruction* det)
  :G4VSensitiveDetector(name),Detector(det)
{   
  numberOfVoxelX = voxelX;
  numberOfVoxelY = voxelY; 
  numberOfVoxelZ = voxelZ;
       
  targetX = XSize; // Size along the Z axis
  targetY = YSize; 
  targetZ = ZSize;       

  G4String HCname;
  collectionName.insert(HCname="trackerCollection");
  totalEnergyDeposit = 0;
}

CellTrackerSD::~CellTrackerSD()
{
}

void CellTrackerSD::Initialize(G4HCofThisEvent* HCE)
{
  trackerCollection = new CellTrackerHitsCollection
   (SensitiveDetectorName,collectionName[0]); 

 static G4int HCID = -1;
  if(HCID<0)
  {HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE -> AddHitsCollection( HCID, trackerCollection ); 
 
  totalEnergyDeposit = 0;
}

G4bool CellTrackerSD::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist)
{
if(!ROhist)
    return false;

  // Retrieve the energy deposit in the target
  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep==0.) return false;

  G4double xx = aStep -> GetPreStepPoint() -> GetPosition().x();
  G4double yy = aStep -> GetPreStepPoint() -> GetPosition().y();
  G4double zz = aStep -> GetPreStepPoint() -> GetPosition().z();

  G4int k = ROhist -> GetReplicaNumber(1);
  G4int i = ROhist -> GetReplicaNumber(2);
  G4int j = ROhist -> GetReplicaNumber();

  G4double voxelWidthX = targetX/numberOfVoxelX;
  G4double voxelWidthY = targetY/numberOfVoxelY;
  G4double voxelWidthZ = targetZ/numberOfVoxelZ;
 
  G4double x = (- numberOfVoxelX+1+2*i)* voxelWidthX/2;
  G4double y = (- numberOfVoxelY+1+2*j)* voxelWidthY/2;
  G4double z = (- numberOfVoxelZ+1+2*k)* voxelWidthZ/2;

  
    G4cout << "Energy deposit (MeV): "<< edep/MeV
    << " in (" <<  xx/mm <<", " << yy/mm <<", "<< zz/mm << ")"
    << " mm"<<G4endl;
   
 /* 
  G4cout<< "ReadOut:" << edep/MeV <<" in ( " <<  x/mm <<", " << y/mm <<", "<< z/mm << " )"
	  << " mm"<<G4endl;          
  */

  CellAnalysisManager* analysis = CellAnalysisManager::getInstance();
  analysis -> FillEnergyDeposit(x/mm, y/mm, edep/MeV);
  analysis -> FillProfile(z/mm, edep/MeV);	 

  // Store the information in a hit collection
  CellTrackerHit* newHit = new CellTrackerHit();
  newHit -> SetEdep(edep); // Edep in the voxel
  trackerCollection -> insert( newHit );
 
  return true;
}

void CellTrackerSD::EndOfEvent(G4HCofThisEvent*)
{
G4int NbHits = trackerCollection->entries();

 for (G4int i=0;i<NbHits;i++)
     {
       totalEnergyDeposit =  totalEnergyDeposit + (*trackerCollection)[i]-> GetEdep();  
     };
 
CellRunAction* runAction = (CellRunAction*) G4RunManager::GetRunManager()->GetUserRunAction();
  
runAction-> IntegrateEnergyDeposit(totalEnergyDeposit);

 // G4RunManager::GetRunManager()-> GetUserRunAction() -> IntegrateEnergyDeposit(totalEnergyDeposit);

 G4cout << "Energy deposit fo the event " << totalEnergyDeposit/MeV << " MeV"<<G4endl;
}


