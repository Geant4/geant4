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
// $Id: TargetSD.cc,v 1.1 2003-10-08 12:32:11 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
#include "TargetSD.hh"

#include "TargetHit.hh"
#include "TargetConstruction.hh"

#include "G4Step.hh"
#include "G4VProcess.hh"
#include "G4StepPoint.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"

// #include "G4ios.hh"


TargetSD::TargetSD(G4String name, TargetConstruction* tgt)
 :G4VSensitiveDetector(name),Target(tgt)
{
  collectionName.insert("TgtCollection");
}


TargetSD::~TargetSD()
{
}


void TargetSD::Initialize(G4HCofThisEvent* HCE)
{
  TgtCollection = new TargetHitsCollection
                      (SensitiveDetectorName,collectionName[0]); 
  HitID = -1;
}


G4bool TargetSD::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist)
{
  G4Track* aTrack = aStep->GetTrack();
  const G4String& partName = aTrack->GetDefinition()->GetParticleName();

  // Kill original proton track after it undergoes inelastic scattering
  //  
  const G4VProcess* proc = aStep->GetPostStepPoint()->GetProcessDefinedStep();
  const G4String& procName = proc->GetProcessName();
  if (partName == "proton" && aTrack->GetTrackID() == 1) {
    if (procName == "ProtonInelastic") {
      aTrack->SetTrackStatus(fStopAndKill);
    }
  }

  // Store the pion information
  //
  if(partName == "pi+" || partName == "pi-") {

    G4int stepNo = aTrack->GetCurrentStepNumber();
    G4int parentID = aTrack->GetParentID();
    
    if(parentID == 1 && stepNo == 1) {
      // parent is track from generator AND daughter has taken only 1 step

      G4double costheta = aTrack->GetVertexMomentumDirection().z();
      G4double ke = aTrack->GetVertexKineticEnergy();
      G4double Q = aTrack->GetDefinition()->GetPDGCharge();
      //      G4cout << "Parent ID = " << parentID << G4endl;
      //      G4cout << "Step number = " << stepNo << G4endl;
      //      G4cout << " costheta, KE, Q = " << costheta << " , " 
      //                                      << ke << " , "
      //	                              << Q << G4endl;

      TargetHit* tgtHit = new TargetHit();
      tgtHit->StoreTheta(costheta);
      tgtHit->StoreEnergy(ke);
      tgtHit->StoreCharge(Q);
      HitID = TgtCollection->insert(tgtHit) - 1;

      return true;
    }
  }
  return false; 
}


void TargetSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  static G4int HCID = -1;
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection(HCID,TgtCollection);
}


void TargetSD::clear()
{} 


void TargetSD::DrawAll()
{} 


void TargetSD::PrintAll()
{} 


