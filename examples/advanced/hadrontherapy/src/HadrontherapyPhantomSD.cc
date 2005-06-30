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
// $Id: HadrontherapyPhantomSD.cc; May 2005
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the National Institute for Nuclear Physics, Catania, Italy
// (b) National Institute for Nuclear Physics Section of Genova, genova, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------

#include "HadrontherapyPhantomSD.hh"
#include "HadrontherapyAnalysisManager.hh"
#include "HadrontherapyPhantomHit.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"

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
 
  if(aStep -> GetPreStepPoint() -> GetPhysicalVolume() -> GetName() != "PhantomPhys")
    return false;

  G4double energyDeposit = aStep -> GetTotalEnergyDeposit();
 
  if(energyDeposit == 0.) return false;

  // Read voxel indexes: i is the x index, k is the z index

  G4int k  = ROhist -> GetReplicaNumber(0);
  G4int i  = ROhist -> GetReplicaNumber(2);
  G4int j  = ROhist -> GetReplicaNumber(1);
 
  G4String particleName = aStep -> GetTrack() -> GetDynamicParticle() -> 
                           GetDefinition() -> GetParticleName();

  if(energyDeposit != 0)                       
    {  
      // Create a hit with the information of position in the phantom and energy deposit     
      HadrontherapyPhantomHit* phantomHit = new HadrontherapyPhantomHit();       
      phantomHit -> SetEdepAndPosition(i, j, k, energyDeposit); 
      HitsCollection -> insert(phantomHit);
    }

  // Energy deposit of secondary particles along X (integrated on Y and Z)
#ifdef G4ANALYSIS_USE 	

 HadrontherapyAnalysisManager* analysis = 
			HadrontherapyAnalysisManager::getInstance();

 if(energyDeposit != 0)                       
    {  
   
 if(aStep -> GetTrack() -> GetTrackID()!= 1)
   {
     if (particleName == "proton")
           analysis -> SecondaryProtonEnergyDeposit(i, energyDeposit/MeV);
  
     if (particleName == "neutron")
     analysis -> SecondaryNeutronEnergyDeposit(i, energyDeposit/MeV);

     if (particleName == "alpha")
       analysis -> SecondaryAlphaEnergyDeposit(i, energyDeposit/MeV);

     if (particleName == "gamma")
       analysis -> SecondaryGammaEnergyDeposit(i, energyDeposit/MeV);
       
     if (particleName == "e-")
       analysis -> SecondaryElectronEnergyDeposit(i, energyDeposit/MeV);
       
     if (particleName == "triton")
       analysis -> SecondaryTritonEnergyDeposit(i, energyDeposit/MeV);
  
     if (particleName == "deuteron")
       analysis -> SecondaryDeuteronEnergyDeposit(i, energyDeposit/MeV);
       
    if (particleName == "pi+" || particleName == "pi-" ||  particleName == "pi0")
       analysis -> SecondaryPionEnergyDeposit(i, energyDeposit/MeV);   	
   }
    }
#endif

  return true;
}

void HadrontherapyPhantomSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  static G4int HCID = -1;
  if(HCID < 0)
    { 
      HCID = GetCollectionID(0); 
    }
  HCE -> AddHitsCollection(HCID,HitsCollection);
}

