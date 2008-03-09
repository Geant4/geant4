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
// $Id: HadrontherapyDetectorSD.cc; 
// Last modified: G.A.P.Cirrone March 2008;
// 
// See more at: http://geant4infn.wikispaces.com/HadrontherapyExample
//
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

#include "HadrontherapyDetectorSD.hh"
#include "HadrontherapyAnalysisManager.hh"
#include "HadrontherapyDetectorHit.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"

HadrontherapyDetectorSD::HadrontherapyDetectorSD(G4String name):G4VSensitiveDetector(name)
{ 
  G4String HCname;
  collectionName.insert(HCname="HadrontherapyDetectorHitsCollection");
 
  HitsCollection = NULL; 
  G4String sensitiveDetectorName = name;
}

HadrontherapyDetectorSD::~HadrontherapyDetectorSD()
{ 
}

void HadrontherapyDetectorSD::Initialize(G4HCofThisEvent*)
{ 
  HitsCollection = new HadrontherapyDetectorHitsCollection(sensitiveDetectorName,
							  collectionName[0]);
}


G4bool HadrontherapyDetectorSD::ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist)
{
  if(!ROhist)
    return false;
 
  if(aStep -> GetPreStepPoint() -> GetPhysicalVolume() -> GetName() != "DetectorPhys")
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
      // Create a hit with the information of position is in the detector     
      HadrontherapyDetectorHit* detectorHit = new HadrontherapyDetectorHit();       
      detectorHit -> SetEdepAndPosition(i, j, k, energyDeposit); 
      HitsCollection -> insert(detectorHit);
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

void HadrontherapyDetectorSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  static G4int HCID = -1;
  if(HCID < 0)
    { 
      HCID = GetCollectionID(0); 
    }
  HCE -> AddHitsCollection(HCID,HitsCollection);
}

