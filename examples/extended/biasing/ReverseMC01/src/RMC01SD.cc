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
// $Id: RMC01SD.cc,v 1.1 2009-11-19 22:41:18 ldesorgh Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//////////////////////////////////////////////////////////////
//      Class Name:	RMC01SD
//	Author:       	L. Desorgher
// 	Organisation: 	SpaceIT GmbH
//	Contract:	ESA contract 21435/08/NL/AT
// 	Customer:     	ESA/ESTEC
//////////////////////////////////////////////////////////////
#include "RMC01SD.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4Track.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4VProcess.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"


//#include "G4AdjointAnalysisManager.hh"
#include "G4THitsCollection.hh"
#include "G4RunManager.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4Proton.hh"




class G4Step;

RMC01SD::RMC01SD(G4String name)
:G4VSensitiveDetector(name)
{  
   collectionName.insert("edep"); 
   collectionName.insert("current_electron"); 
   collectionName.insert("current_proton"); 
   collectionName.insert("current_gamma");
;
}
////////////////////////////////////////////////////////////////////////////////
//
RMC01SD::~RMC01SD()
{; 
}
////////////////////////////////////////////////////////////////////////////////
//
void RMC01SD::Initialize(G4HCofThisEvent* HCE)
{  TotalEventEdep=0.;
   static G4int HCID = -1;
   
   EventEdepCollection = new G4THitsCollection<RMC01DoubleWithWeightHit>(SensitiveDetectorName,collectionName[0]);
   HCID = GetCollectionID(0);
   HCE->AddHitsCollection(HCID,EventEdepCollection);
  
   ElectronCurrentCollection = new G4THitsCollection<RMC01DoubleWithWeightHit>(SensitiveDetectorName,collectionName[1]);
   HCID = GetCollectionID(1);
   HCE->AddHitsCollection(HCID,ElectronCurrentCollection);
 
   ProtonCurrentCollection = new G4THitsCollection<RMC01DoubleWithWeightHit>(SensitiveDetectorName,collectionName[2]);
   HCID = GetCollectionID(2);
   HCE->AddHitsCollection(HCID,ProtonCurrentCollection);
  
   GammaCurrentCollection = new G4THitsCollection<RMC01DoubleWithWeightHit>(SensitiveDetectorName,collectionName[3]);
   HCID = GetCollectionID(3);
   HCE->AddHitsCollection(HCID,GammaCurrentCollection); 
  
   
   
}
////////////////////////////////////////////////////////////////////////////////
//
G4bool RMC01SD::ProcessHits(G4Step*aStep,G4TouchableHistory* )
{
  G4double weight =aStep->GetTrack()->GetWeight();
  G4double edep= aStep->GetTotalEnergyDeposit();
  if (edep >0) EventEdepCollection->insert(new RMC01DoubleWithWeightHit(edep,weight)); 
  
  G4StepPoint* preStepPoint =aStep->GetPreStepPoint();
  
  if (preStepPoint->GetStepStatus() == fGeomBoundary ){
  
     // Entering the sensitive volume
     weight=preStepPoint->GetWeight();
     G4double eKin = preStepPoint->GetKineticEnergy();
     G4ParticleDefinition* thePartDef = aStep->GetTrack()->GetDefinition(); 
     if (thePartDef == G4Electron::Electron())  ElectronCurrentCollection->insert(new RMC01DoubleWithWeightHit(eKin,weight));
     else if (thePartDef == G4Gamma::Gamma())   GammaCurrentCollection->insert(new RMC01DoubleWithWeightHit(eKin,weight));
     else if (thePartDef == G4Proton::Proton()) ProtonCurrentCollection->insert(new RMC01DoubleWithWeightHit(eKin,weight));
  
  }   
  return true;   
}  
////////////////////////////////////////////////////////////////////////////////
//
void RMC01SD::EndOfEvent(G4HCofThisEvent*)
{ EventEdepCollection->insert(new RMC01DoubleWithWeightHit(TotalEventEdep,1.)); 
}
////////////////////////////////////////////////////////////////////////////////
//
void RMC01SD::clear()
{;
}
////////////////////////////////////////////////////////////////////////////////
//
void RMC01SD::DrawAll()
{;
} 
////////////////////////////////////////////////////////////////////////////////
//
void RMC01SD::PrintAll()
{;
}

