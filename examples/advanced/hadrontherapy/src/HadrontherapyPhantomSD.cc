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
// $Id: HadrontherapyPhantomSD.cc,v 2.0
// --------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// --------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone, G. Russo
// Laboratori Nazionali del Sud - INFN, Catania, Italy
//
// --------------------------------------------------------------


#include "HadrontherapyPhantomSD.hh"
#include "HadrontherapyPhantomHit.hh"
#include "HadrontherapyAnalysisManager.hh"
#include "HadrontherapyDetectorConstruction.hh"
#include "G4Track.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ios.hh"
#include "HadrontherapyRunAction.hh"


//-----------------------------------------------------------------------
HadrontherapyPhantomSD::HadrontherapyPhantomSD(G4String name):G4VSensitiveDetector(name)
{
 
  
}
// ----------------------------------------------------------------------
HadrontherapyPhantomSD::~HadrontherapyPhantomSD()
{
 
}
// ----------------------------------------------------------------------
void HadrontherapyPhantomSD::Initialize(G4HCofThisEvent*)
{

}
// ----------------------------------------------------------------------
G4bool HadrontherapyPhantomSD::ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist)
{
  if(!ROhist)
    return false;

  if(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() != "PhantomPhys")
    return false;

  G4double energyDeposit;
  energyDeposit = aStep->GetTotalEnergyDeposit();
  if(energyDeposit == 0.)
    return false;

// Read Voxel indexes: i is the x index, k is the z index
  
extern G4int k;
k  = ROhist->GetReplicaNumber(1);
extern G4int i;
i  = ROhist->GetReplicaNumber(2);
extern G4int j;
j  = ROhist->GetReplicaNumber(0);

// Definition of a 3D matrix for the collection of energy deposit
extern G4double matrix[40][40][40];

 

G4int numberOfVoxelZ = 40;
G4double voxelWidthZ = 1 *mm;

//  G4double x = (-numberOfVoxelZ+1+2*i)*voxelWidthZ/2; 
// G4double y = (- numberOfVoxelZ+1+2*j)*voxelWidthZ/2;
// G4double z = (- numberOfVoxelZ+1+2*k)*voxelWidthZ/2;

  G4Track * theTrack = aStep->GetTrack();
//G4double TrackID = theTrack->GetTrackID();
      
  G4String name = aStep->GetTrack()->GetDefinition()->GetParticleName();
  //  G4cout << "$$$$$$$$" << name << G4endl;
 
   
  if(energyDeposit != 0)                       
    { 
      matrix[i][j][k] = matrix[i][j][k] + energyDeposit;
    }
  
 return true;
}
// --------------------------------------------------------------
void HadrontherapyPhantomSD::EndOfEvent(G4HCofThisEvent*)
{
}
// --------------------------------------------------------------
void HadrontherapyPhantomSD::clear()
{
} 
// --------------------------------------------------------------
void HadrontherapyPhantomSD::DrawAll()
{
}
// --------------------------------------------------------------
void HadrontherapyPhantomSD::PrintAll()
{
}



