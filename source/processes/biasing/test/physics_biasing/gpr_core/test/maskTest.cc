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
// $Id: maskTest.cc,v 1.4 2007-09-11 03:01:44 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// J. Tinslay, August 2007. 
//
#include "G4GPRProcessWrappers.hh"
#include "G4VDiscreteProcess.hh"
#include "G4VParticleChange.hh"

#include "G4GPRSeedT.hh"
#include "G4GPRElementSuperStore.hh"
#include "G4GPRSimpleGenerator.hh"

#include "G4GPRMask.hh"

#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4GRSVolume.hh"
#include "G4GPRTriggerSuperStore.hh"
#include "TestSetup.hh"
#include "G4Gamma.hh"

#include "G4GPRPhysicsListManagerSuperStore.hh"
#include "G4GPRNode.hh"
#include "G4GPRManager.hh"

using namespace G4GPRProcessWrappers;

// Process function
G4double DiscreteGPIL1(const G4Track& track,
		       G4double previousStepSize,
		       G4ForceCondition* condition) 
{
  
  return 1.0;
}

G4double DiscreteGPIL2(const G4Track& track,
		       G4double previousStepSize,
		       G4ForceCondition* condition) 
{
  
  return 2.0;
}

G4double DiscreteGPIL3(const G4Track& track,
		       G4double previousStepSize,
		       G4ForceCondition* condition) 
{
  
  return 3.0;
}

G4double DiscreteGPIL4(const G4Track& track,
		       G4double previousStepSize,
		       G4ForceCondition* condition) 
{
  
  return 4.0;
}

G4bool MaskTrigger(const G4Track& track, const G4Step& step) 
{
  G4cout<<"jane vol "<<track.GetVolume()->GetName()<<G4endl;
  //  G4cout<<"jane executing MyTrigger "<<track->GetTrackID()<<G4endl;
  return (track.GetVolume()->GetName() == "VolB_Phys" ? true : false);
}

int main(int argc, char** argv) {

  // Construct world
  G4LogicalVolume* logicalWorld = new G4LogicalVolume(new G4Box("World", 2.0*m, 2.0*m, 2.0*m),
                                      G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic"),
                                      "World");
  
  G4VPhysicalVolume* world = new G4PVPlacement(0,
                                               G4ThreeVector(),
                                               logicalWorld,
                                               "World",
                                               0,   
                                               false,
                                               0);

  G4LogicalVolume* volA_log = new G4LogicalVolume(new G4Box("World", 1.0*m, 1.0*m, 1.0*m),
                                      G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic"),
                                      "VolA_Log");
  
  G4VPhysicalVolume* volA_phys = new G4PVPlacement(0,
                                               G4ThreeVector(0, 0, -1.0*m),
                                               volA_log,
                                               "VolA_Phys",
                                               0,   
                                               false,
                                               0);

  G4LogicalVolume* volB_log = new G4LogicalVolume(new G4Box("World", 1.0*m, 1.0*m, 1.0*m),
                                      G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic"),
                                      "VolB_Log");
  
  G4VPhysicalVolume* volB_phys = new G4PVPlacement(0,
                                               G4ThreeVector(0, 0, 1.0*m),
                                               volB_log,
                                               "VolB_Phys",
                                               0,   
                                               false,
                                               0);

  
  G4GRSVolume* touchable_A = new G4GRSVolume(volA_phys,NULL,G4ThreeVector(0,0,0));      
  G4GRSVolume* touchable_B = new G4GRSVolume(volB_phys,NULL,G4ThreeVector(0,0,0)); 

  typedef G4GPRSeedT<G4GPRProcessLists::DiscreteGPIL> Seed;

  Seed* seed1 = new Seed("Seed1", &DiscreteGPIL1, G4GPRPlacement::First);
  Seed* seed2 = new Seed("Seed2", &DiscreteGPIL2, G4GPRPlacement::Second);
  Seed* seed3 = new Seed("Seed3", &DiscreteGPIL3, G4GPRPlacement::Third);
  Seed* seed4 = new Seed("Seed4", &DiscreteGPIL4, G4GPRPlacement::Fourth);

  std::vector<unsigned> processes;
  processes.push_back(G4GPRPlacement::First);
  processes.push_back(G4GPRPlacement::Third);

  G4GPRMask* mask = new G4GPRMask("Mask",  processes);
  
  G4ParticleDefinition* def = G4Gamma::Definition();

  G4GPRPhysicsListManager* physicsListManager = &(*G4GPRPhysicsListManagerSuperStore::Instance())[def];
  G4GPRPhysicsList* physicsList = physicsListManager->GetDefaultList();

  G4GPRElementStore* elementStore = &(*G4GPRElementSuperStore::Instance())[def][physicsList];

  elementStore->G4GPRManagerT<Seed>::Register(seed1);
  elementStore->G4GPRManagerT<Seed>::Register(seed2);
  elementStore->G4GPRManagerT<Seed>::Register(seed3);
  elementStore->G4GPRManagerT<Seed>::Register(seed4);

  G4GPRElementStoreT<G4GPRProcessLists::DiscreteGPIL>* store(elementStore);
  store->G4GPRManagerT<G4GPRMask>::Register(mask);

  G4GPRTriggerStore* triggerStore = &(*G4GPRTriggerSuperStore::Instance())[def][physicsList];

  triggerStore->G4GPRTriggerManagerT<G4GPRTriggerTypes::Geometry::StartBoundary>::Register(&MaskTrigger, mask, 
										     &G4GPRMask::ChangeState);
  
    // Create and register key nodes with trigger manager so that know when an element has been activated or deactivated
  G4GPRNode* node1 = new G4GPRNode;

  triggerStore->G4GPRTriggerManagerT<G4GPRTriggerTypes::Geometry::StartBoundary>::Register(&MaskTrigger, node1, &G4GPRNode::FlipState);

  G4GPRKeyStore* keyStore = &(*G4GPRKeySuperStore::Instance())[def][physicsList];

  keyStore->G4GPRKeyManagerT<Seed::List>::AddNode(node1);

  // Generate process list
  typedef std::vector< G4GPRDiscreteGPIL > ProcessList;

  G4Track* track = new G4Track; 
  G4Step* step = new G4Step;

  track->SetTouchableHandle(touchable_A);

  // Each G4ParticleDefinition will have its own G4GPRManager to make processing quicker
  G4GPRManager gprManager(def);
  gprManager.Fire<G4GPRTriggerTypes::Geometry::StartBoundary>(*track, *step);  

  ProcessList* result(0);
  gprManager.GetList<G4GPRProcessLists::DiscreteGPIL>(result);
  G4cout<<"jane generated size should be 4 and is: "<<result->size()<<G4endl;
  // Iterate over process list
  for (ProcessList::iterator iter = result->begin(); iter != result->end(); ++iter) {
    G4cout<<"Executing functor :" <<iter->GetIdentifier()<<" : ";
    G4cout<< (*iter)(*track, 0, 0);
    G4cout<<G4endl;
  }

  track->SetTouchableHandle(touchable_B);

  gprManager.Fire<G4GPRTriggerTypes::Geometry::StartBoundary>(*track, *step);  
  gprManager.GetList<G4GPRProcessLists::DiscreteGPIL>(result);

  G4cout<<"jane generated size should be 2 and is: "<<result->size()<<G4endl;

  // Iterate over process list
  for (ProcessList::iterator iter = result->begin(); iter != result->end(); ++iter) {
    G4cout<<"Executing functor :" <<iter->GetIdentifier()<<" : ";
    G4cout<< (*iter)(*track, 0, 0);
    G4cout<<G4endl;
  }

  return 0;
}
