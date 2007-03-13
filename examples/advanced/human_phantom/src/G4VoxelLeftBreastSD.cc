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
// Authors: S. Guatelli and M. G. Pia, INFN Genova, Italy
// 
// Based on code developed by the undergraduate student G. Guerrieri 
// Note: this is a preliminary beta-version of the code; an improved 
// version will be distributed in the next Geant4 public release, compliant
// with the design in a forthcoming publication, and subject to a 
// design and code review.
//
#include "G4VoxelLeftBreastSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#ifdef G4ANALYSIS_USE	
#include "G4HumanPhantomAnalysisManager.hh"
#endif

G4VoxelLeftBreastSD::G4VoxelLeftBreastSD(G4String name)
:G4VSensitiveDetector(name)
{
}

G4VoxelLeftBreastSD::~G4VoxelLeftBreastSD()
{;}

void G4VoxelLeftBreastSD::Initialize(G4HCofThisEvent*)
{}

G4bool G4VoxelLeftBreastSD::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist)
{
  if(!ROhist) return false;
   
    // Check the volume
  if(aStep -> GetPreStepPoint() -> GetPhysicalVolume() -> 
     GetName() != "LeftBreast") 
    return false;

  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep==0.) return false;

 if(edep != 0)                       
	    { 
             
#ifdef G4ANALYSIS_USE 
G4HumanPhantomAnalysisManager* analysis = 
  G4HumanPhantomAnalysisManager::getInstance();   

 if (aStep -> GetPreStepPoint() -> GetPhysicalVolume() -> 
     GetName() == "LeftBreast")
   {
     G4int sector = ROhist -> GetReplicaNumber();
     G4int slice = ROhist -> GetReplicaNumber(1);	             
     analysis -> voxelLeftBreastEnergyDeposit(slice,sector,edep/MeV);    
     //     G4cout << "LeftBreast:" << "slice: " << slice << ",sector: "<< sector << " "<< edep/MeV << G4endl;           
   }
#endif	
    }
  return true;
}

void G4VoxelLeftBreastSD::EndOfEvent(G4HCofThisEvent*)
{}
