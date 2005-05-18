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
// $Id: HadrontherapyPhantomSD.cc,v 3.0, September 2004
// --------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// --------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone, F. Di Rosa, G. Russo
// Laboratori Nazionali del Sud - INFN, Catania, Italy
//
// --------------------------------------------------------------


#include "HadrontherapyPhantomSD.hh"
#include "HadrontherapyPhantomHit.hh"
#include "HadrontherapyAnalysisManager.hh"
#include "HadrontherapyDetectorConstruction.hh"
#include "HadrontherapyPhantomROGeometry.hh"
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


HadrontherapyPhantomSD::HadrontherapyPhantomSD(G4String name):G4VSensitiveDetector(name)
{ 
 G4String HCname;
 collectionName.insert(HCname="HadrontherapyPhantomHitsCollection");
 
 HitsCollection = NULL; 
 G4String sensitiveDetectorName = name;
}

HadrontherapyPhantomSD::~HadrontherapyPhantomSD()
{ 
}


void HadrontherapyPhantomSD::Initialize(G4HCofThisEvent*)
{ 
 HitsCollection = new HadrontherapyPhantomHitsCollection(sensitiveDetectorName,
                                                         collectionName[0]);
}


G4bool HadrontherapyPhantomSD::ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist)
{
  if(!ROhist)
    return false;
 
  if(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() != "PhantomPhys")
    return false;

  G4double energyDeposit = aStep -> GetTotalEnergyDeposit();
 
  if(energyDeposit == 0.) return false;

 // Read Voxel indexes: i is the x index, k is the z index

 G4int k  = ROhist->GetReplicaNumber(0);
 G4int i  = ROhist->GetReplicaNumber(2);
 G4int j  = ROhist->GetReplicaNumber(1);

  if(energyDeposit != 0)                       
    {       
      HadrontherapyPhantomHit* phantomHit = new HadrontherapyPhantomHit();
            
      phantomHit->SetEdepAndPosition(i, j, k, energyDeposit); 
      HitsCollection->insert(phantomHit);
    }
 
  //G4cout<< "Energy deposit in the sensitive detector:" << energyDeposit/MeV
  //	<< "in"<< i <<" " << j << " "<< "" << k<< G4endl;
  
 return true;
}

void HadrontherapyPhantomSD::EndOfEvent(G4HCofThisEvent* HCE)
{
 static G4int HCID = -1;
  if(HCID<0)
  	{ 
	HCID = GetCollectionID(0); 
	}
  HCE->AddHitsCollection(HCID,HitsCollection);
}

void HadrontherapyPhantomSD::clear()
{
}
 
void HadrontherapyPhantomSD::DrawAll()
{
}

void HadrontherapyPhantomSD::PrintAll()
{
}



