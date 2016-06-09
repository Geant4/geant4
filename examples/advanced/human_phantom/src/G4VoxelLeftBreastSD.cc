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
#include "G4HumanPhantomAnalysis.hh"

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
 
// Get analysis manager
 /* G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

 if (aStep -> GetPreStepPoint() -> GetPhysicalVolume() -> 
     GetName() == "LeftBreast")
   {
     G4int sector = ROhist -> GetReplicaNumber();
     G4int slice = ROhist -> GetReplicaNumber(1);	             
     analysisManager -> FillH2(1,slice, sector,edep/MeV);    
     //     G4cout << "LeftBreast:" << "slice: " << slice << ",sector: "<< sector << " "<< edep/MeV << G4endl;           
   }
*/
    }
  return true;
}

void G4VoxelLeftBreastSD::EndOfEvent(G4HCofThisEvent*)
{}
