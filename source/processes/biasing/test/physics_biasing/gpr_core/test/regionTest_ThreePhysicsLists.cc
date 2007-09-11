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
// $Id: regionTest_ThreePhysicsLists.cc,v 1.3 2007-09-11 03:01:44 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// J. Tinslay, May 2007. 
//
#include "G4GPRProcessWrappers.hh"
#include "G4VDiscreteProcess.hh"
#include "G4VParticleChange.hh"

#include "G4GPRSeedT.hh"
#include "G4GPRElementSuperStore.hh"
#include "G4GPRSimpleGenerator.hh"
#include "G4Gamma.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4Region.hh"
#include "G4GRSVolume.hh"
#include "G4GPRTriggerSuperStore.hh"
#include "TestSetup.hh"
#include "G4GPRPhysicsList.hh"
#include "G4GPRPhysicsListManagerSuperStore.hh"
#include "G4GPRPhysicsListTriggerSuperStore.hh"
#include "G4GPRManager.hh"

using namespace TestSetup;

namespace Default {

  // Process functions
  G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&) 
  {
    G4cout<<"Execute AtRestDoIt_Default"<<G4endl;
    return 0;
  }

}
namespace NewDefault {

  // Process functions
  G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&) 
  {
    G4cout<<"Execute AtRestDoItNew_Default"<<G4endl;
    return 0;
  }

}

namespace RegionA {

  G4bool Trigger(const G4Track& track) 
  {
    return (track.GetVolume()->GetLogicalVolume()->GetRegion()->GetName() == "RegionA");
  }

  // Process functions
  G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&) 
  {
    G4cout<<"Execute AtRestDoIt_RegionA"<<G4endl;
    return 0;
  }
  
  G4VParticleChange* ContinuousDoIt(const G4Track&, const G4Step&) 
  {
    G4cout<<"Execute ContinuousDoIt_RegionA"<<G4endl;
    return 0;
  }
  
  G4VParticleChange* DiscreteDoIt(const G4Track&, const G4Step&) 
  {
    G4cout<<"Execute DiscreteDoIt_RegionA"<<G4endl;
    return 0;
  }
  
  G4double AtRestGPIL(const G4Track& track, G4ForceCondition* condition) 
  {
    G4cout<<"Execute AtRestGPIL_RegionA"<<G4endl;
    return 0;
  }
  
  G4double ContinuousGPIL(const G4Track& track,
			  G4double  previousStepSize,
			  G4double  currentMinimumStep,
			  G4double& proposedSafety,
			  G4GPILSelection* selection)
  {
    G4cout<<"Execute ContinuousGPIL_RegionA"<<G4endl;
    return 0;
  }
  
  G4double DiscreteGPIL(const G4Track& track,
			G4double   previousStepSize,
			G4ForceCondition* condition) {
    G4cout<<"Execute DiscreteGPIL_RegionA"<<G4endl;
    return 0;
  }
}

namespace RegionB {

  G4bool Trigger(const G4Track& track) 
  {
    return (track.GetVolume()->GetLogicalVolume()->GetRegion()->GetName() == "RegionB");
  }
  
  // Process functions
  G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&) 
  {
    G4cout<<"Execute AtRestDoIt_RegionB"<<G4endl;
    return 0;
  }
  
  G4VParticleChange* ContinuousDoIt(const G4Track&, const G4Step&) 
  {
    G4cout<<"Execute ContinuousDoIt_RegionB"<<G4endl;
    return 0;
  }
  
  G4VParticleChange* DiscreteDoIt(const G4Track&, const G4Step&) 
  {
    G4cout<<"Execute DiscreteDoIt_RegionB"<<G4endl;
    return 0;
  }
  
  G4double AtRestGPIL(const G4Track& track, G4ForceCondition* condition) 
  {
  G4cout<<"Execute AtRestGPIL_RegionB"<<G4endl;
  return 0;
  }
  
  G4double ContinuousGPIL(const G4Track& track,
			  G4double  previousStepSize,
			  G4double  currentMinimumStep,
			     G4double& proposedSafety,
			  G4GPILSelection* selection)
  {
    G4cout<<"Execute ContinuousGPIL_RegionB"<<G4endl;
    return 0;
  }
  
  G4double DiscreteGPIL(const G4Track& track,
			G4double   previousStepSize,
			G4ForceCondition* condition) 
  {
    G4cout<<"Execute DiscreteGPIL_RegionB"<<G4endl;
    return 0;
  }
}
///////////////////////////////

int main(int argc, char** argv) {

  // Construct world
  G4LogicalVolume* logicalWorld = new G4LogicalVolume(new G4Box("World", 2.0*m, 2.0*m, 2.0*m),
                                      G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic"),
                                      "World");
  
  new G4PVPlacement(0,
		    G4ThreeVector(),
		    logicalWorld,
		    "World",
		    0,   
		    false,
		    0);

  G4LogicalVolume* volA_log = new G4LogicalVolume(new G4Box("World", 1.0*m, 1.0*m, 1.0*m),
                                      G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic"),
                                      "volA_log");
  
  G4VPhysicalVolume* volA_phys = new G4PVPlacement(0,
                                               G4ThreeVector(0, 0, -1.0*m),
                                               volA_log,
                                               "volA_phys",
                                               0,   
                                               false,
                                               0);

  G4LogicalVolume* volB_log = new G4LogicalVolume(new G4Box("World", 1.0*m, 1.0*m, 1.0*m),
                                      G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic"),
                                      "volB_log");
  
  G4VPhysicalVolume* volB_phys = new G4PVPlacement(0,
                                               G4ThreeVector(0, 0, 1.0*m),
                                               volB_log,
                                               "volB_phys",
                                               0,   
                                               false,
                                               0);

  G4LogicalVolume* volC_log = new G4LogicalVolume(new G4Box("World", 1.0*m, 1.0*m, 1.0*m),
                                      G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic"),
                                      "volC_log");
  
  G4VPhysicalVolume* volC_phys = new G4PVPlacement(0,
                                               G4ThreeVector(0, 0, 1.0*m),
                                               volC_log,
                                               "volC_phys",
                                               0,   
                                               false,
                                               0);

  G4Region* regionA = new G4Region("RegionA");
  volA_log->SetRegion(regionA);

  G4Region* regionB = new G4Region("RegionB");
  volB_log->SetRegion(regionB);

  G4Region* regionC = new G4Region("RegionC");
  volC_log->SetRegion(regionC);

  G4GRSVolume* touchable_A = new G4GRSVolume(volA_phys,NULL,G4ThreeVector(0,0,0));      
  G4GRSVolume* touchable_B = new G4GRSVolume(volB_phys,NULL,G4ThreeVector(0,0,0));      
  G4GRSVolume* touchable_C = new G4GRSVolume(volC_phys,NULL,G4ThreeVector(0,0,0));      

  G4GPRSeedT<G4GPRProcessLists::AtRestDoIt>* atRestDoIt_Default = 
    new G4GPRSeedT<G4GPRProcessLists::AtRestDoIt>("AtRestDoIt_Default", &Default::AtRestDoIt, G4GPRPlacement::First);

  G4GPRSeedT<G4GPRProcessLists::AtRestDoIt>* atRestDoIt_A = 
    new G4GPRSeedT<G4GPRProcessLists::AtRestDoIt>("AtRestDoIt_A", &RegionA::AtRestDoIt, G4GPRPlacement::First);

  G4GPRSeedT<G4GPRProcessLists::AtRestDoIt>* atRestDoIt_B = 
    new G4GPRSeedT<G4GPRProcessLists::AtRestDoIt>("AtRestDoIt_B", &RegionB::AtRestDoIt, G4GPRPlacement::First);
  
  // Two physics lists in addition to the default one
  G4GPRPhysicsList* listA = new G4GPRPhysicsList("ListA", false);
  G4GPRPhysicsList* listB = new G4GPRPhysicsList("ListB", false);

  G4ParticleDefinition* def = G4Gamma::Definition();

  // Register physics lists with physics list manager
  G4GPRPhysicsListManager* physicsListManager = &(*G4GPRPhysicsListManagerSuperStore::Instance())[def];

  physicsListManager->Register(listA);
  physicsListManager->Register(listB);

  // Register physics list triggers with physics list trigger manager
  G4GPRTriggerStore* physicsListTrigger = &(*G4GPRPhysicsListTriggerSuperStore::Instance())[def];

  physicsListTrigger->G4GPRTriggerManagerT<G4GPRTriggerTypes::Geometry::NewRegion>::Register(&RegionA::Trigger, listA, 
										       &G4GPRPhysicsList::FlipState, listA->GetState());
  
  physicsListTrigger->G4GPRTriggerManagerT<G4GPRTriggerTypes::Geometry::NewRegion>::Register(&RegionB::Trigger, listB, 
										       &G4GPRPhysicsList::FlipState, listA->GetState());
   
  // Register elements with appropriate stores
  G4GPRElementStore* defaultStore = &(*G4GPRElementSuperStore::Instance())[def][physicsListManager->GetDefaultList()];
  G4GPRElementStore* regionAStore = &(*G4GPRElementSuperStore::Instance())[def][listA];
  G4GPRElementStore* regionBStore = &(*G4GPRElementSuperStore::Instance())[def][listB];

  defaultStore->G4GPRManagerT< G4GPRSeedT<G4GPRProcessLists::AtRestDoIt> >::Register(atRestDoIt_Default);
  regionAStore->G4GPRManagerT< G4GPRSeedT<G4GPRProcessLists::AtRestDoIt> >::Register(atRestDoIt_A);
  regionBStore->G4GPRManagerT< G4GPRSeedT<G4GPRProcessLists::AtRestDoIt> >::Register(atRestDoIt_B);

  // Process list

  // Dummy track and step
  G4DynamicParticle* dynamic= new G4DynamicParticle(def, 0, G4ThreeVector(0, 0, 0));
  G4Track* track = new G4Track(dynamic, 0, G4ThreeVector(0,0,0)); 

  track->SetTouchableHandle(touchable_A);

  // Each G4ParticleDefinition will have its own G4GPRManager to make processing quicker
  G4GPRManager gprManager(def);
  
  gprManager.Fire<G4GPRTriggerTypes::Geometry::NewRegion>(*track);
  
  // Generate list in region A
  typedef std::vector< G4GPRProcessWrappers::Wrappers<G4GPRProcessLists::AtRestDoIt>::SeedWrapper > ProcessList;

  G4cout<<"jane generating physics list "<<gprManager.GetActivePhysicsList()->GetName()<<G4endl;
  ProcessList* result1(0);
  gprManager.GetList<G4GPRProcessLists::AtRestDoIt>(result1);
  G4cout<<"jane process list for region A "<<result1<<G4endl;
  PrintList(result1);

  // Test caching
  ProcessList* result2(0);
  gprManager.GetList<G4GPRProcessLists::AtRestDoIt>(result2);
  G4cout<<"jane process list for region A "<<result2<<G4endl;
  PrintList(result2);

  assert(result1 == result2);

  // Switch to region B
  track->SetTouchableHandle(touchable_B);
  gprManager.Fire<G4GPRTriggerTypes::Geometry::NewRegion>(*track);

  gprManager.GetList<G4GPRProcessLists::AtRestDoIt>(result1);
  G4cout<<"jane process list for region B "<<result1<<G4endl;
  PrintList(result1);

  gprManager.GetList<G4GPRProcessLists::AtRestDoIt>(result2);
  G4cout<<"jane process list for region B "<<result2<<G4endl;
  PrintList(result2);

  assert(result1 == result2);

  // Switch to region C, default physics list should be used
  track->SetTouchableHandle(touchable_C);
  gprManager.Fire<G4GPRTriggerTypes::Geometry::NewRegion>(*track);

  gprManager.GetList<G4GPRProcessLists::AtRestDoIt>(result1);
  G4cout<<"jane process list for region C"<<G4endl;
  PrintList(result1);

  gprManager.GetList<G4GPRProcessLists::AtRestDoIt>(result2);
  G4cout<<"jane process list for region C"<<G4endl;
  PrintList(result2);

  assert(result1 == result2);

  // Now create a new default physics list, overriding the previous one. May be useful when want to create new default
  // physics list in place of prepackaged one 
  physicsListManager->SetDefaultList(new G4GPRPhysicsList("NewDefault", true));

  G4GPRSeedT<G4GPRProcessLists::AtRestDoIt>* atRestDoIt_NewDefault = 
    new G4GPRSeedT<G4GPRProcessLists::AtRestDoIt>("AtRestDoIt_NewDefault", &NewDefault::AtRestDoIt, G4GPRPlacement::First);

  defaultStore = &(*G4GPRElementSuperStore::Instance())[def][physicsListManager->GetDefaultList()];
  defaultStore->G4GPRManagerT< G4GPRSeedT<G4GPRProcessLists::AtRestDoIt> >::Register(atRestDoIt_NewDefault);
  gprManager.GetList<G4GPRProcessLists::AtRestDoIt>(result1);
  G4cout<<"jane process list for region C, new default"<<G4endl;
  PrintList(result1);  
  return 0;

}
